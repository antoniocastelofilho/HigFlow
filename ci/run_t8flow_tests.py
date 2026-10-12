#!/usr/bin/env python3
"""Roda a suite do t8Flow e COBRA o contrato.

Mesmo padrao do ci/run_higtree_tests.py, e as guardas sao as mesmas porque as
armadilhas sao as mesmas -- cada uma delas custou um ciclo de diagnostico naquela
suite:

  1. exigir os prefixos de dependencia no ambiente e DIZER qual falta, em vez de
     reportar "0 de N construidos".  Caminho RELATIVO repassado a um `make -C` e'
     resolvido a partir do diretorio errado, e o sintoma e' "nao construido" --
     falha silenciosa, nao erro.  Por isso tudo e' absolutizado aqui;
  2. `make clean` por dimensao: os objetos nao carregam o DIM no nome, e um build
     2D sobre arvore 3D liga objetos 3D -- falha como segfault e numero errado,
     nao como erro de build;
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
CONTRATO = os.path.join(T8FLOW, "src", "t8flow-contract.h")

# As dependencias vivem em `bibliotecas/`, que o git ignora -- logo nao estao
# neste worktree.  Os padroes apontam para a arvore principal; sobrescreva por
# argumento quando o prefixo for outro.
PADROES = {
    "T8CODE": "/home/castelo/HiGFlow-System/HigFlow/bibliotecas/t8code/install",
    "HIGTREE_SRC": "/home/castelo/HiGFlow-System/HigFlow/higtree/src",
    "PETSC_ROOT": "/home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4",
}

# nome do teste -> dimensoes e numeros de processo a exercitar
TESTES = [
    ("test-nucleo-herdado", (2, 3), (1,)),
]


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


def constroi(dim, prefixos, verboso):
    cmd = ["make", "-C", T8FLOW, "DIM=%d" % dim]
    cmd += ["%s=%s" % (k, v) for k, v in prefixos.items()]
    limpa = subprocess.run(["make", "-C", T8FLOW, "clean"],
                           capture_output=True, text=True)
    if limpa.returncode != 0:
        print("  AVISO: `make clean` falhou:\n%s" % limpa.stderr.strip())
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("  BUILD FALHOU em DIM=%d:" % dim)
        print("\n".join(r.stderr.strip().splitlines()[-15:]))
        return False
    if verboso:
        print("  build DIM=%d ok" % dim)
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
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--t8code", help="prefixo de instalacao do t8code")
    p.add_argument("--higtree-src", help="diretorio higtree/src com o nucleo herdado")
    p.add_argument("--petsc", help="raiz do PETSc")
    p.add_argument("--dim", type=int, choices=(2, 3), help="rodar so' esta dimensao")
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

    dims = (a.dim,) if a.dim else (2, 3)
    clausulas = le_contrato(CONTRATO)
    print("Contrato: %d clausulas em %s" % (len(clausulas), os.path.relpath(CONTRATO, RAIZ)))

    # casos que passaram, por (teste, caso), unindo todas as dimensoes
    passaram, falharam, problemas = set(), set(), []
    total_casos = 0

    for dim in dims:
        print("\n=== DIM=%d ===" % dim)
        if not constroi(dim, prefixos, a.verbose):
            problemas.append("build DIM=%d" % dim)
            continue
        for nome, dims_do_teste, nps in TESTES:
            if dim not in dims_do_teste:
                continue
            binario = os.path.join(T8FLOW, "tests", "%s-%dd" % (nome, dim))
            if not os.path.exists(binario):
                problemas.append("%s-%dd nao foi construido" % (nome, dim))
                print("  %s-%dd NAO CONSTRUIDO" % (nome, dim))
                continue
            for np in nps:
                ok, bad, resumo, _ = roda(binario, np, a.verbose)
                total_casos += len(ok) + len(bad)
                print("  %s np=%d: %d passaram, %d falharam" % (nome, np, len(ok), len(bad)))
                if not resumo:
                    problemas.append("%s-%dd np=%d sem linha de resumo" % (nome, dim, np))
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
