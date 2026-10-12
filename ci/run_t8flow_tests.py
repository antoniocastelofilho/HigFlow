#!/usr/bin/env python3
"""Roda a suite do t8Flow e COBRA o contrato.

Mesmo padrao do ci/run_higtree_tests.py, e as guardas sao as mesmas porque as
armadilhas sao as mesmas -- cada uma delas custou um ciclo de diagnostico naquela
suite:

  1. exigir os prefixos de dependencia e DIZER qual falta, em vez de reportar
     "0 de N construidos".  Caminho RELATIVO repassado a um `make -C` e'
     resolvido a partir do diretorio errado, e o sintoma e' "nao construido" --
     falha silenciosa, nao erro.  Por isso tudo e' absolutizado aqui;
  2. `make clean` antes de cada lado: os objetos do nucleo herdado sao por
     dimensao (`Point` e' `real[DIM]`) e nao carregam o DIM no nome, entao um
     build 2D sobre objetos 3D liga os errados -- falha como segfault e numero
     errado, nao como erro de build;
  3. exigir a linha `resumo` de cada binario: sem ela, "nenhuma falha" e "nenhum
     caso rodou" sao indistinguiveis, que e' exatamente como se apresentam
     binario obsoleto e segfault no meio da execucao.

E a quarta, que e' especifica deste projeto: todo binario vai por `mpirun`,
inclusive em np=1.  O arranque chama `higtree_initialize`, logo `MPI_Init`:
rodado direto, trava sem imprimir nada e aparece como timeout, nao como falha.
"""

import argparse
import os
import re
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
T8FLOW = os.path.join(RAIZ, "t8flow")
CONTRATO = os.path.join(T8FLOW, "src", "comum", "t8flow-contract.h")

# As dependencias vivem em `bibliotecas/`, que o git ignora -- logo nao estao
# neste worktree.  Os padroes apontam para a arvore principal; sobrescreva por
# argumento quando o prefixo for outro.
PADROES = {
    "T8CODE": "/home/castelo/HiGFlow-System/HigFlow/bibliotecas/t8code/install",
    "HIGTREE_SRC": "/home/castelo/HiGFlow-System/HigFlow/higtree/src",
    "PETSC_ROOT": "/home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4",
}

# CODIGO SEPARADO POR DIMENSAO (decisao 7): cada lado tem o seu alvo de make e os
# seus binarios, em vez de um fonte compilado duas vezes.  O lado 3D ainda nao
# existe, e esta' aqui como lista vazia em vez de ausente -- assim o driver o
# reporta como "sem testes" em vez de omiti-lo.
#
# Cada entrada: (nome do binario, numeros de processo a exercitar).
LADOS = {
    "2d": [
        ("test-mls2d", (1,)),
        ("test-solver2d", (1,)),
        ("test-floresta2d", (1,)),
        ("test-malha2d", (1,)),
        ("test-refino2d", (1,)),
        ("test-numeracao2d", (1,)),
    ],
    "3d": [],
}


def le_contrato(caminho):
    """Devolve [(clausula, [casos declarados])], na ordem do arquivo."""
    clausulas = []
    atual = None
    with open(caminho, encoding="utf-8") as f:
        for linha in f:
            m = re.match(r"\s*//!\s*([A-Z]\d+)\s+(.*)", linha)
            if m:
                atual = (m.group(1), [])
                clausulas.append(atual)
                continue
            m = re.match(r"\s*//!\s*COBERTA-POR:\s*(\S+):(\S+)", linha)
            if m and atual is not None:
                atual[1].append((m.group(1), m.group(2)))
    return clausulas


def constroi(lado, prefixos, verboso):
    """Limpa e constroi um lado (`2d` ou `3d`), devolvendo True se deu.

    O `make clean` vem ANTES e nao e' opcional: os objetos do nucleo herdado sao
    por dimensao (`Point` e' `real[DIM]`) e nao carregam o DIM no nome, entao
    objetos da dimensao anterior seriam ligados em silencio.  Um clean que FALHA e'
    tao ruim quanto nenhum, por isso o aviso -- foi assim que a guarda nasceu
    inoperante na primeira versao, com o Makefile exigindo DIM no alvo clean.
    """
    cmd = ["make", "-C", T8FLOW, lado]
    cmd += ["%s=%s" % (k, v) for k, v in prefixos.items()]
    limpa = subprocess.run(["make", "-C", T8FLOW, "clean"],
                           capture_output=True, text=True)
    if limpa.returncode != 0:
        print("  AVISO: `make clean` falhou:\n%s" % limpa.stderr.strip())
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("  BUILD FALHOU no lado %s:" % lado)
        print("\n".join(r.stderr.strip().splitlines()[-15:]))
        return False
    if verboso:
        print("  build %s ok" % lado)
    return True


def roda(binario, np, verboso):
    """Devolve (casos_ok, casos_falhos, teve_resumo, saida)."""
    cmd = ["mpirun", "-n", str(np), binario]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return set(), set(), False, "TIMEOUT apos 600 s"

    ok, falhos, resumo = set(), set(), False
    for linha in r.stdout.splitlines():
        m = re.match(r"caso (\S+) PASS", linha)
        if m:
            ok.add(m.group(1))
            continue
        m = re.match(r"caso (\S+) FAIL (.*)", linha)
        if m:
            falhos.add(m.group(1))
            print("    FAIL %s: %s" % (m.group(1), m.group(2)))
            continue
        if re.match(r"resumo \d+ casos, \d+ falharam", linha):
            resumo = True
    if not resumo:
        print("    SEM LINHA DE RESUMO -- binario obsoleto, segfault no meio, ou "
              "nenhum caso rodou.  Saida bruta:")
        print("\n".join(("      " + l) for l in (r.stdout + r.stderr).splitlines()[-12:]))
    return ok, falhos, resumo, r.stdout


def main():
    """Constroi, roda e COBRA o contrato; devolve 0 so' se tudo fechar.

    Falha se algum caso reprovar, se algum binario nao for construido, se algum
    binario nao imprimir a linha de resumo, ou se alguma clausula do contrato ficar
    sem caso que a afirme.  As quatro coisas contam igual: clausula sem teste e'
    falha, nao omissao.
    """
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--t8code", help="prefixo de instalacao do t8code")
    p.add_argument("--higtree-src", help="diretorio higtree/src com o nucleo herdado")
    p.add_argument("--petsc", help="raiz do PETSc")
    p.add_argument("--lado", choices=("2d", "3d"), help="rodar so' este lado (padrao: 2d)")
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args()

    prefixos = dict(PADROES)
    if a.t8code:      prefixos["T8CODE"] = os.path.abspath(a.t8code)
    if a.higtree_src: prefixos["HIGTREE_SRC"] = os.path.abspath(a.higtree_src)
    if a.petsc:       prefixos["PETSC_ROOT"] = os.path.abspath(a.petsc)

    # GUARDA 1: dizer QUAL falta, com o caminho que foi tentado.
    faltando = []
    for chave, marcador in (("T8CODE", "include/t8.h"),
                            ("HIGTREE_SRC", "wls.h"),
                            ("PETSC_ROOT", "include/petsc.h")):
        if not os.path.exists(os.path.join(prefixos[chave], marcador)):
            faltando.append("%s: nao achei %s em %s"
                            % (chave, marcador, prefixos[chave]))
    if faltando:
        print("Dependencias ausentes -- a suite NAO rodou:")
        for f in faltando:
            print("  " + f)
        return 2

    lados = (a.lado,) if a.lado else ("2d",)
    clausulas = le_contrato(CONTRATO)
    print("Contrato: %d clausulas em %s" % (len(clausulas), os.path.relpath(CONTRATO, RAIZ)))

    # casos que passaram, por (teste, caso), unindo todas as dimensoes
    passaram, falharam, problemas = set(), set(), []
    total_casos = 0

    for lado in lados:
        print("\n=== lado %s ===" % lado)
        testes = LADOS.get(lado, [])
        if not testes:
            print("  nenhum teste neste lado ainda")
            continue
        if not constroi(lado, prefixos, a.verbose):
            problemas.append("build do lado %s" % lado)
            continue
        for nome, nps in testes:
            binario = os.path.join(T8FLOW, "tests", lado, nome)
            if not os.path.exists(binario):
                problemas.append("%s nao foi construido" % nome)
                print("  %s NAO CONSTRUIDO" % nome)
                continue
            for np in nps:
                ok, bad, resumo, _ = roda(binario, np, a.verbose)
                total_casos += len(ok) + len(bad)
                print("  %-18s np=%d: %d passaram, %d falharam" % (nome, np, len(ok), len(bad)))
                if not resumo:
                    problemas.append("%s np=%d sem linha de resumo" % (nome, np))
                for c in ok:
                    passaram.add((nome, c))
                for c in bad:
                    falharam.add((nome, c))

    # GUARDA DO CONTRATO: cada clausula tem de ter caso declarado E passando.
    print("\n=== CONTRATO ===")
    reprovadas = 0
    for clausula, casos in clausulas:
        if not casos:
            print("  %-4s SEM TESTE QUE RODE" % clausula)
            reprovadas += 1
            continue
        estado = "OK"
        for ref in casos:
            if ref in falharam:
                estado = "REPROVADA (%s)" % ref[1]
            elif ref not in passaram:
                estado = "CASO AUSENTE (%s:%s)" % ref
                break
        if estado != "OK":
            reprovadas += 1
        print("  %-4s %s" % (clausula, estado))

    print("\nResumo: %d casos executados, %d falharam; contrato %d/%d"
          % (total_casos, len(falharam), len(clausulas) - reprovadas, len(clausulas)))
    if problemas:
        print("Problemas de execucao:")
        for pr in problemas:
            print("  " + pr)

    return 0 if (not falharam and not problemas and reprovadas == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
