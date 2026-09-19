#!/usr/bin/env python3
"""Run the HiGTree unit tests.

The HiGFlow suite (run_suite.py) BUILDS higtree and never tests it: the library
is treated as a prerequisite of the examples, and a defect inside it only shows
up as a wrong number several layers away -- if it shows up at all.  Every defect
fixed in domain.c during 2026-09-15..18 lived in higtree and was found by a 3D
simulation aborting, not by a test.

This driver runs the tests in higtree/tests, which assert VALUES rather than
shapes.  There is an older suite in higtree/atf-tests (71 cases, 2020) that needs
`atf-c`; it is not installed anywhere and those tests have not run in five years.
They are kept as reference material, not as a suite.

Each test binary prints one line per case:

    caso <name> PASS
    caso <name> FAIL <what failed, with the numbers>
    resumo <n> casos, <k> falharam
"""

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
HIGTREE = os.path.join(ROOT, "higtree")
TESTDIR = os.path.join(HIGTREE, "tests")


class Test:
    def __init__(self, name, dims=(2, 3)):
        self.name = name
        self.dims = dims


TESTS = [
    # O valor que o estencil produz fora do dominio, onde responde o fechamento
    # por condicao de contorno.  Roda em 2D e 3D: o defeito de compactacao do
    # ponto de consulta some quando proj_dir == DIM-1, entao uma dimensao so'
    # deixaria metade do espaco sem cobertura.
    Test("test-stencil-value", dims=(2, 3)),
    # Qual parede fecha o estencil, medido pelo VALOR: paredes com valores
    # constantes e DIFERENTES fazem o resultado identificar a escolha sozinho,
    # sem gancho de depuracao e sem observar o interior da biblioteca.
    Test("test-stencil-selection", dims=(2, 3)),
    # A mesma geometria como uma arvore e como duas tem de dar os mesmos valores.
    Test("test-partition-independence", dims=(2, 3)),
]


def build_for_dim(dim, timeout):
    """Build libhig<dim>d and the test binaries for one dimension.

    `make clean` is not optional, for the same reason it is not optional in
    run_suite.py: the object files carry no dimension in their names, so after a
    DIM=3 build their timestamps make them look current to a DIM=2 build.  Make
    then rebuilds nothing and the 2D tests link 3D objects, which fails as
    segfaults and wrong numbers rather than as a build error.
    """
    env = dict(os.environ)
    env.setdefault("HIGTREE_DIR", HIGTREE)

    if "PETSC_DIR" not in env:
        return False, ("PETSC_DIR nao esta' no ambiente: carregue o varsrc de "
                       "DENTRO desta arvore (ele usa $(pwd))")

    steps = [["make", "-C", HIGTREE, "clean"],
             ["make", "-C", HIGTREE, "DIM=%d" % dim],
             ["make", "-C", TESTDIR, "clean"],
             ["make", "-C", TESTDIR, "DIM=%d" % dim]]

    for cmd in steps:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, env=env)
        if r.returncode != 0:
            if cmd[-1] == "clean":
                continue                      # nada a limpar nao e' erro
            err = next((l for l in (r.stdout + r.stderr).splitlines()
                        if "error:" in l or "Error" in l), "")
            return False, err.strip()[:70]
    return True, ""


CASO = re.compile(r"^caso (\S+) (PASS|FAIL) ?(.*)$")
RESUMO = re.compile(r"^resumo (\d+) casos, (\d+) falharam$")


def run_test(test, dim, timeout):
    """Run one binary.  Returns (cases, erro) with cases a list of
    (name, ok, detail).  `erro` is set when the binary did not report at all."""
    exe = os.path.join(TESTDIR, "%s-%dd" % (test.name, dim))
    if not os.path.exists(exe):
        return [], "binario nao construido"
    try:
        r = subprocess.run([exe], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return [], "timeout apos %ds" % timeout

    cases, viu_resumo = [], False
    for line in r.stdout.splitlines():
        m = CASO.match(line)
        if m:
            cases.append((m.group(1), m.group(2) == "PASS", m.group(3)))
            continue
        if RESUMO.match(line):
            viu_resumo = True

    # A linha de resumo e' o controle de que o binario chegou ao fim.  Sem ela,
    # "nenhuma falha" pode significar "nenhum caso rodou" -- que e' como um
    # binario obsoleto, um segfault no meio ou um build parcial se apresentam.
    if not viu_resumo:
        detalhe = "sem linha de resumo (rc=%d)" % r.returncode
        if r.returncode < 0:
            detalhe += ", morto por sinal %d" % -r.returncode
        return cases, detalhe
    return cases, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dim", default="2,3",
                    help="dimensoes a exercitar (padrao: 2,3)")
    ap.add_argument("--test", action="append", default=[],
                    help="rodar apenas este teste (pode repetir)")
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--no-build", action="store_true",
                    help="nao reconstruir; so' faz sentido se a arvore ja' esta' "
                         "na dimensao pedida")
    args = ap.parse_args()

    dims = [int(d) for d in args.dim.split(",") if d.strip()]
    tests = [t for t in TESTS if not args.test or t.name in args.test]
    if not tests:
        print("nenhum teste corresponde a %s" % args.test)
        return 2

    print("%-26s %-4s %-34s %-7s %s" % ("TESTE", "DIM", "CASO", "RESULT", "DETALHE"))
    print("-" * 100)

    total = falhas = 0
    for dim in dims:
        alvos = [t for t in tests if dim in t.dims]
        if not alvos:
            continue
        if not args.no_build:
            ok, err = build_for_dim(dim, args.timeout)
            if not ok:
                print("%-26s %-4d %-34s %-7s %s"
                      % ("(build)", dim, "", "ERRO", err))
                falhas += 1
                continue
        for t in alvos:
            cases, erro = run_test(t, dim, args.timeout)
            # Um caso pode emitir VARIAS linhas de falha, uma por assercao.  A
            # contagem e' por CASO, nao por linha: "2 de 11" quando ha' 6 casos
            # mistura as duas unidades e o total deixa de significar algo.  As
            # linhas de detalhe continuam todas visiveis.
            for nome, ok, detalhe in cases:
                print("%-26s %-4d %-34s %-7s %s"
                      % (t.name, dim, nome, "ok" if ok else "FAIL", detalhe))
            por_caso = {}
            for nome, ok, _ in cases:
                por_caso[nome] = por_caso.get(nome, True) and ok
            total += len(por_caso)
            falhas += sum(1 for ok in por_caso.values() if not ok)
            if erro:
                falhas += 1
                print("%-26s %-4d %-34s %-7s %s" % (t.name, dim, "", "ERRO", erro))

    print("-" * 100)
    print("%d de %d caso(s) passaram" % (total - falhas, total))
    return 0 if falhas == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
