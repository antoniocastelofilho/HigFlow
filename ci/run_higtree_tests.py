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
    def __init__(self, name, dims=(2, 3), nps=(1,), mpi=False, so_t8code=False):
        self.name = name
        self.dims = dims
        # `mpi` diz que o binario chama higtree_initialize.  Ele entao SEMPRE vai
        # por mpirun, inclusive em np=1: rodado direto, o MPI_Init trava sem
        # imprimir nada e o teste aparece como timeout, nao como falha.
        self.mpi = mpi
        # Teste que so' existe com o t8code -- nao tem lado MTree.  Sem --t8code
        # ele nao e' construido nem cobrado, e a guarda de registro o ignora.
        self.so_t8code = so_t8code
        # Numeros de processos a exercitar.  O padrao e' (1,): a maioria dos testes
        # nao chama MPI_Init e, lancada sob mpirun, travaria.  Quem declara np>1 TEM
        # de chamar higtree_initialize e reduzir os veredictos entre os ranks antes
        # de imprimir -- senao cada rank imprime a sua versao e falha de um rank
        # some no meio das linhas dos outros.
        self.nps = nps


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
    # Estencil atravessando salto de nivel 4:1 em malha NAO graduada -- a
    # capacidade que qualquer segunda implementacao de malha tera' de preservar.
    Test("test-level-jump", dims=(2, 3), mpi=True),
    # Localizacao por ponto: a celula contem o ponto, o empate nao depende da
    # divisao do dominio, e a convencao de empate esta' fixada.
    Test("test-point-location", dims=(2, 3), mpi=True),
    # O ramo ON_BOUNDARY do despacho, que os outros cinco nao alcancam -- medido
    # com contador: zero chamadas aos fechamentos *_boundary neles, 56 neste.
    # Foi nesse ramo que sobreviveram tres dos sete sitios do defeito de
    # compactacao, verdes por ausencia de teste e nao por estarem certos.
    Test("test-boundary-path", dims=(2, 3)),

    # Montagem SERIAL de um sim_facet_domain.  Ate' 2026-09-19 era impossivel: o
    # sfbi[] so' era preenchido pelo caminho PARTICIONADO, entao um sfd montado
    # sem MPI estourava na primeira consulta.  A metade local passou a viver no
    # domain.c como sfd_compute_sfbi; este teste impede que ela volte a depender
    # do particionamento sem ninguem perceber.
    Test("test-facet-domain-serial", dims=(2, 3)),

    # As quatro consultas de celula que o contrato de Mesh vai congelar --
    # hig_get_center, hig_get_delta, hig_get_cid e o iterador -- somam mais de 850
    # usos em higflow/src e ate' aqui so' tinham cobertura INDIRETA: passavam
    # porque os outros testes dependem delas.  Oraculos analiticos, dois niveis de
    # refino.
    Test("test-cell-queries", dims=(2, 3)),

    # A franja, pelo criterio do que o estencil PEDE e nao por tamanho fixado.
    # build-fringe.cpp sao 815 linhas sem cobertura, e quando o t8code entrar ele
    # traz a propria camada de ghost -- o contrato precisa estar escrito antes.
    Test("test-fringe-support", dims=(2, 3)),

    # A franja sob particionamento REAL.  O teste serial acima monta a franja a
    # mao e afirma o contrato dela; este exercita quem a PRODUZ -- o
    # build-fringe.cpp e o balanceador --, que so' rodam sob MPI.  Unico teste
    # da suite com np>1: chama higtree_initialize e reduz os dados entre ranks
    # antes de concluir, no rank 0.
    Test("test-fringe-parallel", dims=(2, 3), nps=(1, 2, 3), mpi=True),

    # A FRONTEIRA DAS CONSULTAS.  86% das chamadas em laco quente sao leitura de
    # centro, tamanho e indice, e nenhuma precisa da arvore: o instantaneo as
    # atende com arranjo plano.  Com --t8code, o caso extra mostra os DOIS
    # backends preenchendo o mesmo conjunto por caminhos que nao se parecem --
    # o t8code sem materializar octree nenhum.
    Test("test-mesh-snapshot", dims=(2, 3), mpi=True),

    # O MESMO instantaneo, sob particionamento real.  O caso acima roda em np=1 e
    # monta a franja a mao -- ele confirma a convencao de numeracao que ele proprio
    # impos.  Quem a impoe de verdade e' o `psd_synced_mapper`, que so' roda sob
    # MPI e e' quem o solver usa.  O modo de falhar e' silencioso: franja com id em
    # [0, n) sobrescreve a linha de uma celula local e o instantaneo sai com o
    # tamanho certo e uma celula errada dentro.
    Test("test-mesh-snapshot-parallel", dims=(2, 3), nps=(1, 2, 3), mpi=True),

    # O instantaneo pendurado no dominio: quando nasce, e o que acontece se a
    # malha mudar depois.  Estado derivado num tipo de malha so' se sustenta
    # porque a malha e' imutavel apos a montagem -- e "e' imutavel" e' observacao
    # sobre o codigo de hoje, nao garantia.  Quem a transforma em garantia sao a
    # guarda (aborta) e o detector, e os dois sao afirmados aqui.  O caso da
    # guarda roda em processo filho e mede o SIGABRT: o processo que ela mata e'
    # justamente o que faria a afirmacao.
    # mpi=True mesmo em np=1: ele chama higtree_initialize, e standalone o
    # MPI_Init nao volta -- o sintoma e' timeout de 300s sem UMA linha de
    # saida, que nao se parece nada com "faltou mpirun".
    Test("test-domain-snapshot", dims=(2, 3), mpi=True),

    # C10 pela separacao certa: a malha fornece o SUPORTE, o ajuste de minimos
    # quadrados e' o mesmo para as duas.  Comparar sd_get_stencil contra uma
    # montagem propria do t8code mediria malha e discretizacao juntas e nao
    # diria qual falhou.
    Test("test-stencil-support", dims=(2, 3), mpi=True),

    # Onde o dominio termina, dito por dois mecanismos OPOSTOS: na HiGTree o
    # contorno e' explicito (arvores sim_boundary registradas), no t8code e'
    # implicito (face sem vizinho).  E' a metade de Mesh em C12 e C14; a
    # projecao e a interpolacao em DIM-1 sao de Discretization e ficam nos
    # testes que ja' existem.
    Test("test-boundary-faces", dims=(2, 3), mpi=True),

    # O despacho que escolhe o ramo ON_BOUNDARY (C14).  Os dois backends
    # decidem por criterios DIFERENTES -- caixa da arvore contra face sem
    # vizinho -- e o ultimo caso mede onde eles se separam, em vez de supor.
    Test("test-point-class", dims=(2, 3), mpi=True),

    # O t8code sob particionamento REAL: C13 e P1 a P4.  As sete clausulas ja'
    # verificadas rodam em um processo; esta familia exige a floresta
    # distribuida com ghost, e e' onde a substituicao do particionador de fato
    # acontece.  Espelha o test-fringe-parallel, que faz o mesmo para o MTree.
    Test("test-partition-t8code", dims=(2, 3), nps=(1, 2, 3), mpi=True,
         so_t8code=True),
]


def confere_registro():
    """O driver conhece todo teste que o Makefile constroi, e vice-versa?

    Um teste que existe, compila e passa mas nao esta' na lista do driver nunca
    roda -- e a suite reporta verde sem ele.  Aconteceu no dia em que os testes de
    salto de nivel e localizacao por ponto foram escritos: entraram no Makefile e
    nao aqui, e a tabela seguiu dizendo "10 de 10" como se nada faltasse.  E' a
    mesma familia de falha silenciosa que a suite existe para pegar, entao ela tem
    de se aplicar a si mesma.
    """
    mk = os.path.join(TESTDIR, "Makefile")
    if not os.path.exists(mk):
        return []
    no_make = set(re.findall(r"TESTS\s*\+?=\s*(test-[\w-]+)", open(mk).read()))
    no_driver = {t.name for t in TESTS}
    # Testes que so' existem com o t8code saem dos DOIS lados da comparacao.
    # Tirar de um lado so' faz a guarda disparar na direcao contraria -- foi o
    # que aconteceu na primeira tentativa.
    so_t8 = {t.name for t in TESTS if t.so_t8code}
    no_make -= so_t8
    no_driver -= so_t8
    faltam = []
    for nome in sorted(no_make - no_driver):
        faltam.append("%s esta' no Makefile e nao no driver: nunca roda" % nome)
    for nome in sorted(no_driver - no_make):
        faltam.append("%s esta' no driver e nao no Makefile: nao sera' construido" % nome)
    return faltam


def build_for_dim(dim, timeout, t8code=""):
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

    extra = ["T8CODE=%s" % t8code] if t8code else []
    steps = [["make", "-C", HIGTREE, "clean"],
             ["make", "-C", HIGTREE, "DIM=%d" % dim],
             ["make", "-C", TESTDIR, "clean"],
             ["make", "-C", TESTDIR, "DIM=%d" % dim] + extra]

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


CONTRATO = os.path.join(HIGTREE, "src", "hig-mesh-contract.h")

CLAUSULA = re.compile(r"^//\s+([PC]\d+)\s+\S")
# Uma prova marcada com [t8code] so' existe quando a suite roda com --t8code.
# Sem a marca, o driver nao distinguiria "clausula perdeu o teste" de "clausula
# tem prova extra que esta configuracao nao constroi", e a configuracao padrao
# acusaria falso alarme -- foi o que aconteceu quando a C8 ganhou a segunda prova.
APLICA   = re.compile(r"^//\s+(test-[\w-]+)\s*/\s*(\w+)\s*(\[t8code\])?\s*$")


def le_contrato():
    """Le as clausulas do hig-mesh-contract.h e os casos que as verificam.

    O CABECALHO E' A FONTE.  Manter a lista aqui, em paralelo, seria garantir que
    as duas divergissem -- e clausula que perdeu o teste deixa de ser clausula e
    vira comentario, que e' exatamente o que este arquivo existe para impedir.
    """
    if not os.path.exists(CONTRATO):
        return {}
    clausulas, atual = {}, None
    for linha in open(CONTRATO, encoding="utf-8"):
        linha = linha.rstrip("\n")
        m = CLAUSULA.match(linha)
        if m:
            atual = m.group(1)
            clausulas.setdefault(atual, [])
            continue
        m = APLICA.match(linha)
        if m and atual:
            opcional = m.group(3) is not None
            clausulas[atual].append((m.group(1), m.group(2), opcional))
            continue
        # linha em branco de comentario encerra a clausula corrente
        if linha.strip() in ("//", ""):
            atual = None
    return clausulas


CASO = re.compile(r"^caso (\S+) (PASS|FAIL) ?(.*)$")
RESUMO = re.compile(r"^resumo (\d+) casos, (\d+) falharam$")


def run_test(test, dim, np, timeout):
    """Run one binary.  Returns (cases, erro) with cases a list of
    (name, ok, detail).  `erro` is set when the binary did not report at all."""
    exe = os.path.join(TESTDIR, "%s-%dd" % (test.name, dim))
    if not os.path.exists(exe):
        return [], "binario nao construido"
    if np == 1 and not test.mpi:
        cmd = [exe]
    else:
        cmd = ["mpirun", "-use-hwthread-cpus", "-n", str(np), exe]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
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
    ap.add_argument("--t8code", default=os.environ.get("T8CODE", ""),
                    help="prefixo de instalacao do t8code.  Com ele, o "
                         "test-level-jump ganha o segundo produtor de malha e a "
                         "clausula C11 passa a exigir que os DOIS atravessem o "
                         "salto 4:1.  Sem ele nada muda.")
    ap.add_argument("--no-build", action="store_true",
                    help="nao reconstruir; so' faz sentido se a arvore ja' esta' "
                         "na dimensao pedida")
    args = ap.parse_args()

    # Caminho relativo e' o que se digita, mas ele e' repassado ao
    # `make -C higtree/tests` e resolvido A PARTIR DALI, nao da raiz -- o
    # resultado e' o t8code silenciosamente nao construido e clausulas
    # reportadas como sem teste.  Absolutiza-se aqui, uma vez.
    if args.t8code:
        args.t8code = os.path.abspath(args.t8code)

    dims = [int(d) for d in args.dim.split(",") if d.strip()]
    tests = [t for t in TESTS if not args.test or t.name in args.test]
    if not tests:
        print("nenhum teste corresponde a %s" % args.test)
        return 2

    descasados = confere_registro()
    for aviso in descasados:
        print("REGISTRO  %s" % aviso)
    if descasados:
        print()

    print("%-26s %-4s %-3s %-34s %-7s %s"
          % ("TESTE", "DIM", "NP", "CASO", "RESULT", "DETALHE"))
    print("-" * 100)

    total = falhas = 0
    passou = {}          # (teste, caso) -> passou em TODAS as combinacoes
    for dim in dims:
        alvos = [t for t in tests if dim in t.dims]
        if not alvos:
            continue
        if not args.no_build:
            ok, err = build_for_dim(dim, args.timeout, args.t8code)
            if not ok:
                print("%-26s %-4d %-3s %-34s %-7s %s"
                      % ("(build)", dim, "-", "", "ERRO", err))
                falhas += 1
                continue
        for t in alvos:
          if t.so_t8code and not args.t8code:
            continue
          for np in t.nps:
            cases, erro = run_test(t, dim, np, args.timeout)
            # Um caso pode emitir VARIAS linhas de falha, uma por assercao.  A
            # contagem e' por CASO, nao por linha: "2 de 11" quando ha' 6 casos
            # mistura as duas unidades e o total deixa de significar algo.  As
            # linhas de detalhe continuam todas visiveis.
            for nome, ok, detalhe in cases:
                chave = (t.name, nome)
                passou[chave] = passou.get(chave, True) and ok
                print("%-26s %-4d %-3d %-34s %-7s %s"
                      % (t.name, dim, np, nome, "ok" if ok else "FAIL", detalhe))
            por_caso = {}
            for nome, ok, _ in cases:
                por_caso[nome] = por_caso.get(nome, True) and ok
            total += len(por_caso)
            falhas += sum(1 for ok in por_caso.values() if not ok)
            if erro:
                falhas += 1
                print("%-26s %-4d %-3d %-34s %-7s %s"
                      % (t.name, dim, np, "", "ERRO", erro))

    print("-" * 100)
    print("%d de %d caso(s) passaram" % (total - falhas, total))

    # ------------------------------------------------ o contrato de Mesh
    clausulas = le_contrato()
    com_t8code = bool(args.t8code)
    if clausulas:
        sem_teste, reprovadas = [], []
        for cid, casos in sorted(clausulas.items(),
                                 key=lambda kv: (kv[0][0], int(kv[0][1:]))):
            if not casos:
                sem_teste.append(cid)
                continue
            # a clausula vale se TODOS os casos que a verificam passaram e
            # todos de fato rodaram
            # prova opcional que nao foi construida nao conta contra a clausula
            exigidos = [(t, c) for (t, c, opc) in casos
                        if not opc or com_t8code]
            estados = [passou.get(k) for k in exigidos]
            if not exigidos or any(e is None for e in estados):
                sem_teste.append(cid)
            elif not all(estados):
                reprovadas.append(cid)
        ok = len(clausulas) - len(sem_teste) - len(reprovadas)
        print()
        print("CONTRATO DE MESH: %d de %d clausula(s) verificada(s)"
              % (ok, len(clausulas)))
        for cid in sem_teste:
            print("  %s  SEM TESTE QUE RODE -- clausula sem teste e' comentario,"
                  " nao clausula" % cid)
            falhas += 1
        for cid in reprovadas:
            print("  %s  REPROVADA" % cid)
    if descasados:
        print("%d teste(s) fora de registro -- ver as linhas REGISTRO acima"
              % len(descasados))
    return 0 if falhas == 0 and not descasados else 1


if __name__ == "__main__":
    sys.exit(main())
