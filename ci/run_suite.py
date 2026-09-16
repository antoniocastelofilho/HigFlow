#!/usr/bin/env python3
"""Run the HiGFlow regression suite and check every case against its reference.

Each case is run at several process counts and compared against a single
stored reference.  That is deliberate: the metrics in the reference are
global over the domain, so they must not depend on how the mesh was
partitioned.  A case that passes at one process count and fails at another
is reporting a real parallel bug, not a tolerance problem.

Inputs are never modified in place.  Each run gets a private copy of the
case's input directory, with numsteps and dtp overridden there, so the
tracked YAML files stay untouched.

Usage:
  ./ci/run_suite.py                      # run the default suite
  ./ci/run_suite.py --generate           # (re)create the references
  ./ci/run_suite.py --case example2d_VOF # a single case
  ./ci/run_suite.py --np 1,2,3           # override process counts
  ./ci/run_suite.py --include-slow       # add the long-running cases

Exit status is 0 only if every selected case passed at every process count.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HIGFLOW = os.path.join(ROOT, "higflow")
REFDIR = os.path.join(HERE, "reference")

# Solver options used for every case.  Kept identical across cases so that a
# difference in results is a difference in the physics, not in the solver.
#
# The tolerance is tighter than the one the examples use by default.  With
# bjacobi the preconditioner depends on how the mesh was partitioned, so at
# 1e-10 the residual left behind is enough to move vel.w.max by half a percent
# between one and two processes in the 3D lid-driven case.  At 1e-12 the two
# agree, which is what makes a single reference valid at every process count.
KSP_OPTS = ["-ksp_type", "bcgs", "-pc_type", "bjacobi",
            "-ksp_atol", "1e-12", "-ksp_rtol", "1e-12"]


class Case:
    def __init__(self, name, binary, load, dim, multiphase=False,
                 numsteps=20, dtp=0.005, slow=False, known_broken=None,
                 max_np=None, max_np_reason="", overrides=None, example=None):
        self.name = name
        # Um mesmo exemplo pode render mais de um caso, quando caminhos
        # diferentes do solver sao escolhidos por configuracao e nao por
        # codigo.  `example` diz de qual diretorio o caso sai; `name` continua
        # sendo a identidade do caso e o nome do arquivo de referencia.
        self.example = example or name
        self.overrides = overrides or {}
        self.binary = binary          # link target of the example Makefile
        self.load = load              # input prefix, relative to the case dir
        self.dim = dim
        self.multiphase = multiphase  # mass conservation applies
        self.numsteps = numsteps
        self.dtp = dtp                # VTK write interval
        self.slow = slow
        self.known_broken = known_broken
        self.max_np = max_np            # malha fina demais para dividir mais
        self.max_np_reason = max_np_reason

    @property
    def dir(self):
        return os.path.join(HIGFLOW, self.example)


# The suite covers one case per physics path that the solver can take.
# numsteps is small on purpose: the reference guards against a change in
# behaviour, and a long run costs time without widening coverage.
CASES = [
    Case("example2d_Newt",         "ns-example",     "example-2d.load", 2),
    Case("example2d_Oldroyd",      "ns-example",     "example-2d.load", 2),
    Case("example2d_Gptt",         "ns-example",     "example-2d.load", 2),
    Case("example2d_Newt_contraction", "ns-example", "example-2d.load", 2),
    Case("example2d_VOF",          "ns-example",     "example-2d.load", 2, multiphase=True),
    Case("example2d_VOF_Gptt",     "ns-example",     "example-2d.load", 2, multiphase=True),
    Case("example2d_VOF_Oldroyd",  "ns-example",     "example-2d.load", 2, multiphase=True),
    # Como distribuido, este exemplo roda com eoflow desligado: apesar do nome,
    # ele cobre o caminho multifasico e nao o eletroosmotico.
    Case("example2d_ElectroOsmotic", "ns-example-2d",  "load", 2, multiphase=True),
    # Mesmo exemplo, com o acoplamento eletroosmotico ligado por configuracao.
    # E' o unico caso da suite que entra em step-multiphase-electroosmotic --
    # onde estava o Runge-Kutta que chamava o euler monofasico, um bug que
    # sobreviveu justamente por nada exercitar esse caminho.
    Case("example2d_ElectroOsmotic_eo", "ns-example-2d", "load", 2,
         multiphase=True, example="example2d_ElectroOsmotic",
         overrides={"eoflow": "true", "eoflow0": "true", "eoflow1": "true"}),
    # 10x10x10 = 1000 celulas.  Dividida em dois, a franja passa a ser uma
    # fracao grande de cada subdominio e vel.w.max muda 0,6% entre np=1 e np=2.
    # Nao e' bug de paralelismo: refinando para 20x20x20 a diferenca vai a zero
    # (medido).  Enquanto o exemplo usar esta malha, so np=1 e' comparavel.
    Case("example3d_lid_driven",   "ns-exemple-3d",  "example-3d.load", 3,
         max_np=1, max_np_reason="malha 10x10x10 grosseira demais para dividir"),
    # 162000 cells over 33 blocks; a single step takes minutes.  Opt-in.
    Case("example3d_complex",      "ns-complex-3d",  "example-3d.load", 3,
         numsteps=2, dtp=0.001, slow=True),
    # Diverges to NaN at step 3 with the integral viscoelastic solver.  Listed
    # so the suite reports it as known-broken instead of silently omitting it.
    Case("example2d_KBKZ",         "ns-example",     "example-2d.load", 2,
         known_broken="tensor diverges to NaN at step 3"),
    # Apurado.  O crescimento de FracVol de 0,72% por passo NAO e' defeito do
    # solver: a tampa (bc1) tem componente de velocidade normal ao contorno de
    # 0,1 enquanto as outras tres paredes tem normal zero, entao o dominio nao
    # e' fechado.  Zerando so' essa componente, a massa conserva (variacao
    # maxima 0,057%).  A checagem de sistema fechado de check_mass nao se
    # aplica a este exemplo, e por isso multiphase fica em False.
    #
    # De quebra, como so' a tampa tem fluxo normal, o fluxo liquido pelo
    # contorno e' nao nulo -- incompativel com incompressibilidade num dominio
    # fixo.  Era dai que vinha o "Time step is large!!!" com Tmin = -7,2e6 que
    # este caso exibia em np=2: com a normal zerada as velocidades ficam sas
    # (|V| < 0,35) e a explosao desaparece.
    #
    # O que sobra, e por isso o caso continua quebrado, e' anterior e mais
    # grave: em np=2 a corrida aborta com "free(): invalid next size", ou seja
    # corrupcao de heap.  O valgrind localiza um process_vm_readv dentro de
    # PMPI_Waitall lendo 0 bytes alem do fim de um bloco de 50.744 bytes
    # alocado em dp_create (higtree/src/pdomain.c:901), por
    # psfd_create_property <- higflow_create_distributed_properties <-
    # higflow_rebuild_with_amr, isto e', durante a reconstrucao por adaptacao
    # de malha.  Hipotese a confirmar: dp_sync faz
    # MPI_Irecv(dp->values, 1, psync[...].recv, ...) com tipos derivados que
    # guardam deslocamentos absolutos dentro de values, e esses tipos sao
    # criados uma unica vez sob `if(!psfd->dp_data.psync)` (pdomain.c:1016),
    # sendo liberados apenas em _dp_shared_destroy -- nada fora de pdomain.c
    # os invalida quando a malha muda.  Em np=1 nao ha troca e o caso passa.
    Case("example2d_DynamicMeshAdapt", "ns-example-2d", "load", 2,
         known_broken="corrupcao de heap em np=2: MPI le alem do buffer de dp_create"),
    # Entrada legada: .contr/.par posicionais, sem .par.contr.yaml, entao a
    # suite nao consegue encurtar a corrida -- ele roda o que o .par mandar.
    # dtp menor que o padrao de proposito: o dt deste caso e' 1e-4, entao com
    # 20 passos o tempo final e' 0,002 -- com dtp=0,005 so' o quadro inicial
    # seria escrito e a referencia nao cobriria o avanco no tempo.
    Case("example2d_BMP",          "ns-example",     "example-3d.load", 2,
         dtp=0.0005),
]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def build_for_dim(dim, cases, timeout):
    """Build the libraries and example binaries for one dimension.

    The object files in higflow/src carry no dimension in their names and the
    examples link them directly, so a tree can only hold one DIM at a time.
    The suite therefore builds per dimension and runs that group before moving
    on; it cannot have the 2D and 3D examples built at once.
    """
    higtree = os.path.join(ROOT, "higtree")
    env = dict(os.environ)
    env.setdefault("HIGTREE_DIR", higtree)
    env.setdefault("HIGFLOW_DIR", HIGFLOW)

    # "make clean" is not optional here.  The object files carry no dimension
    # in their names, so after a DIM=3 build the timestamps make them look
    # current to a DIM=2 build: make rebuilds nothing and the 2D examples link
    # 3D objects.  That fails as segfaults and wrong numbers, not as a build
    # error, which is the worst way for it to fail.
    steps = [(higtree, ["make", "-C", higtree, "clean"]),
             (higtree, ["make", "-C", higtree, "DIM=%d" % dim]),
             (HIGFLOW, ["make", "-C", HIGFLOW, "clean"]),
             (HIGFLOW, ["make", "-C", HIGFLOW, "DIM=%d" % dim])]
    for d in dict.fromkeys(c.dir for c in cases):     # ordem estavel, sem repetir
        steps.append((d, ["make", "-C", d, "clean"]))
        steps.append((d, ["make", "-C", d]))

    for where, cmd in steps:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, env=env)
        if r.returncode != 0 and cmd[-1] == "clean":
            continue          # nada a limpar nao e' erro
        if r.returncode != 0:
            err = next((l for l in (r.stdout + r.stderr).splitlines()
                        if "error:" in l), "")
            return False, "%s: %s" % (os.path.basename(where), err.strip()[:60])
    return True, ""


def _e_numero(texto):
    try:
        float(texto.strip())
        return True
    except ValueError:
        return False


def prepare_inputs(case, tmp, numsteps, dtp):
    """Copy the case inputs into tmp and override the run length there."""
    indir = os.path.join(tmp, "input")
    os.makedirs(indir, exist_ok=True)
    src = os.path.join(case.dir, "input")
    for f in os.listdir(src):
        if f.startswith(case.load):
            shutil.copy(os.path.join(src, f), indir)
    # Entrada legada: um .par de numeros soltos, sem rotulo, lido por
    # higflow_load_parameters (higflow/src/hig-flow-io.c) na ordem
    #   1 step | 2 numsteps | 3 t | 4 dt | 5 Re | 6 dts | 7 dtp | 8 frame ...
    # O example2d_BMP usa esse formato e traz numsteps=80000, o que faz a
    # corrida passar de quinze minutos.  Sem encurtar aqui, o caso nao cabe na
    # suite; com isso ele fica no mesmo pe' dos casos em yaml.
    legado = os.path.join(indir, case.load + ".par")
    if os.path.exists(legado) and not os.path.exists(legado + ".contr.yaml"):
        with open(legado) as fp:
            linhas = fp.read().split("\n")
        numericas = [l for l in linhas if l.strip()]
        if len(numericas) >= 7 and all(_e_numero(l) for l in numericas[:7]):
            linhas[1] = str(numsteps)
            linhas[6] = repr(float(dtp))
            with open(legado, "w") as fp:
                fp.write("\n".join(linhas))

    par = os.path.join(indir, case.load + ".par.contr.yaml")
    if os.path.exists(par):
        with open(par) as fp:
            text = fp.read()
        text = re.sub(r"numsteps:\s*\d+", "numsteps: %d" % numsteps, text)
        text = re.sub(r"dtp:\s*[0-9.eE+-]+", "dtp: %g" % dtp, text)
        # Sobrescritas de configuracao do caso.  Trocam so' o valor, deixando
        # o comentario da linha intacto, para que o diff contra o arquivo
        # original continue legivel.
        for key, value in case.overrides.items():
            text, n = re.subn(r"(?m)^(\s*%s:\s*)\S+" % re.escape(key),
                              lambda m: m.group(1) + value, text)
            if n == 0:
                raise SystemExit("caso %s: chave '%s' nao existe em %s"
                                 % (case.name, key, os.path.basename(par)))
        with open(par, "w") as fp:
            fp.write(text)
    return os.path.join(indir, case.load)


def run_case(case, np, numsteps, dtp, timeout):
    """Run one case at np processes.  Returns (ok, vtk_dir_or_None, detail)."""
    exe = os.path.join(case.dir, case.binary)
    if not os.path.exists(exe):
        return False, None, "binary not built: %s" % case.binary

    tmp = tempfile.mkdtemp(prefix="higflow-%s-np%d-" % (case.name, np))
    load = prepare_inputs(case, tmp, numsteps, dtp)
    vtk = os.path.join(tmp, "VTKS")
    os.makedirs(vtk, exist_ok=True)

    cmd = (["mpirun", "-use-hwthread-cpus", "-n", str(np), exe,
            load, os.path.join(tmp, "out.save"), os.path.join(vtk, "p.print")]
           + KSP_OPTS)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, cwd=case.dir)
    except subprocess.TimeoutExpired:
        return False, None, "timed out after %ds" % timeout

    log = r.stdout + r.stderr
    if r.returncode != 0:
        first = next((l for l in log.splitlines()
                      if re.search(r"SEGV|Segmentation|ERROR|Abort", l)), "")
        return False, vtk, "exit %d %s" % (r.returncode, first.strip()[:60])
    if not os.listdir(vtk):
        return False, vtk, "no VTK output produced"

    # NaN e' procurado nos DADOS, nao no log.  O example2d_BMP imprime um
    # diagnostico proprio de convergencia, fabs((old-novo)/novo)*100, que da
    # nan quando a grandeza comparada e' identicamente zero -- Txy no comeco da
    # corrida.  Isso nao diz nada sobre a solucao: os VTKs daquele caso nao tem
    # um unico NaN.  Procurar a substring no stdout reprovava o caso por causa
    # de uma divisao por zero na impressao.
    for f in sorted(os.listdir(vtk)):
        if not f.endswith(".vtk"):
            continue
        with open(os.path.join(vtk, f), errors="replace") as fp:
            if "nan" in fp.read().lower():
                return False, vtk, "NaN in output (%s)" % f
    return True, vtk, ""


def check(vtk, case, tolerance, generate):
    """Run the reference (and mass) checks over a VTK directory."""
    ref = os.path.join(REFDIR, case.name + ".yaml")
    script = os.path.join(HERE, "check_reference.py")

    if generate:
        os.makedirs(REFDIR, exist_ok=True)
        r = run([sys.executable, script, "--generate", vtk, "-o", ref])
        return r.returncode == 0, "wrote %s" % os.path.basename(ref)

    if not os.path.exists(ref):
        return False, "no reference (run with --generate)"

    r = run([sys.executable, script, vtk, ref, "--tolerance", str(tolerance)])
    if r.returncode != 0:
        diff = [l.strip() for l in r.stdout.splitlines() if "DIFF" in l]
        return False, diff[0][:70] if diff else "reference mismatch"

    if case.multiphase:
        m = run([sys.executable, os.path.join(HERE, "check_mass.py"), vtk])
        if m.returncode != 0:
            return False, "mass not conserved"
    return True, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--generate", action="store_true",
                    help="create references instead of checking against them")
    ap.add_argument("--case", action="append", default=[],
                    help="run only this case (repeatable)")
    ap.add_argument("--np", default="1,2,3",
                    help="process counts, comma separated (default 1,2,3)")
    ap.add_argument("--tolerance", type=float, default=0.001,
                    help="relative tolerance (default 0.001)")
    ap.add_argument("--timeout", type=int, default=900,
                    help="per-run timeout in seconds (default 900)")
    ap.add_argument("--include-slow", action="store_true",
                    help="also run the long cases (example3d_complex)")
    ap.add_argument("--include-broken", action="store_true",
                    help="also run cases known to be broken")
    ap.add_argument("--no-build", action="store_true",
                    help="use the binaries already built, do not run make")
    args = ap.parse_args()

    nps = [int(x) for x in args.np.split(",") if x.strip()]
    # References are generated from a single rank: with one block the weighted
    # and unweighted means coincide, so the file is the least partition-biased.
    if args.generate:
        nps = [1]

    cases = [c for c in CASES if not args.case or c.name in args.case]
    selected, skipped = [], []
    for c in cases:
        if c.known_broken and not args.include_broken:
            skipped.append((c, "known broken: " + c.known_broken))
        elif c.slow and not args.include_slow:
            skipped.append((c, "slow, use --include-slow"))
        else:
            selected.append(c)

    if not selected:
        print("No cases selected.")
        return 1

    verb = "Generating references" if args.generate else "Checking"
    print("%s for %d case(s) at np=%s\n" % (verb, len(selected),
                                            ",".join(map(str, nps))))
    print("%-24s %-5s  %-7s  %s" % ("CASE", "NP", "RESULT", "DETAIL"))
    print("-" * 78)

    failures = 0
    # Grouped by dimension: see build_for_dim for why they cannot share a tree.
    for dim in sorted({c.dim for c in selected}):
        group = [c for c in selected if c.dim == dim]
        if not args.no_build:
            ok, why = build_for_dim(dim, group, args.timeout)
            if not ok:
                for c in group:
                    print("%-24s %-5s  %-7s  build failed: %s"
                          % (c.name, "-", "FAIL", why))
                failures += len(group) * len(nps)
                continue
        for c in group:
            for np in nps:
                if c.max_np and np > c.max_np:
                    print("%-24s %-5d  %-7s  %s"
                          % (c.name, np, "skip", c.max_np_reason))
                    continue
                ok, vtk, detail = run_case(c, np, c.numsteps, c.dtp, args.timeout)
                if ok:
                    ok, detail = check(vtk, c, args.tolerance, args.generate)
                status = "ok" if ok else "FAIL"
                if not ok:
                    failures += 1
                print("%-24s %-5d  %-7s  %s" % (c.name, np, status, detail))

                # Apagar o temporario do caso que passou.  Sem isto cada
                # execucao deixa para tras a saida inteira: uma tarde de
                # rodadas acumulou 24 GB em 1011 diretorios e encheu o tmpfs,
                # a ponto de o compilador falhar com "Disk quota exceeded" e o
                # proprio shell nao iniciar.  O de um caso que falhou fica,
                # porque e' a evidencia para diagnosticar a falha.
                if ok and vtk:
                    shutil.rmtree(os.path.dirname(vtk), ignore_errors=True)

    for c, why in skipped:
        print("%-24s %-5s  %-7s  %s" % (c.name, "-", "skip", why))

    print("-" * 78)
    total = len(selected) * len(nps)
    print("%d of %d run(s) passed" % (total - failures, total))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
