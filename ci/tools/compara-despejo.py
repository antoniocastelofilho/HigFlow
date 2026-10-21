#!/usr/bin/env python3
"""Compara dois despejos de HIGFLOW_DESPEJA, com a coordenada como chave.

Existe porque o VTK deste caso não serve de oráculo entre decomposições: a
escrita interpola nos cantos das células e, em domínio com vão, fecha o
estêncil por contorno distante nos mesmos lugares em que a solução é difícil.
Aqui são os valores brutos de `dpp` e `dpu`, e a chave é o centro da célula ou
da faceta -- dois arranjos de blocos cobrem as mesmas células, então as chaves
casam uma a uma.

A diferença é relatada CONTRA A ESCALA DO CAMPO (max - min), não em absoluto:
`p` deste caso chega a 2e4 e `u` a 1,6, e um mesmo erro absoluto significa
coisas opostas nos dois.

Uso:  compara-despejo.py <prefixo-A> <prefixo-B>
      (lê todos os <prefixo>-<rank>.txt)
"""

import glob
import sys
from collections import defaultdict


def le(prefixo):
    campos = defaultdict(dict)
    arquivos = sorted(glob.glob(prefixo + "-*.txt"))
    if not arquivos:
        raise SystemExit("nenhum arquivo em %s-*.txt" % prefixo)
    for a in arquivos:
        for linha in open(a):
            p = linha.split()
            if len(p) != 5:
                continue
            campo = p[0]
            chave = tuple(round(float(v), 9) for v in p[1:4])
            valor = float(p[4])
            anterior = campos[campo].get(chave)
            if anterior is not None and anterior != valor:
                # Mesma coordenada com dois valores: franja repetida entre
                # ranks com conteúdo diferente.  É achado, não ruído.
                print("  AVISO %s %s: dois valores no mesmo ponto, %g e %g"
                      % (campo, chave, anterior, valor))
            campos[campo][chave] = valor
    return campos, arquivos


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    A, nA = le(sys.argv[1])
    B, nB = le(sys.argv[2])
    print("A: %d arquivo(s)   B: %d arquivo(s)" % (len(nA), len(nB)))

    pior = 0.0
    for campo in sorted(set(A) | set(B)):
        a, b = A.get(campo, {}), B.get(campo, {})
        so_a, so_b = set(a) - set(b), set(b) - set(a)
        comuns = set(a) & set(b)
        if not comuns:
            print("%-3s  SEM PONTO EM COMUM (%d x %d)" % (campo, len(a), len(b)))
            pior = float("inf")
            continue
        vals = [a[k] for k in comuns] + [b[k] for k in comuns]
        escala = max(vals) - min(vals)
        dmax, onde = 0.0, None
        for k in comuns:
            d = abs(a[k] - b[k])
            if d > dmax:
                dmax, onde = d, k
        rel = dmax / escala if escala > 0 else 0.0
        pior = max(pior, rel)
        extra = ""
        if so_a or so_b:
            extra = "   SO' EM A: %d   SO' EM B: %d" % (len(so_a), len(so_b))
        print("%-3s  %7d pontos  escala %10.4g  maior dif %10.4g  = %8.3g%% "
              "da escala  em %s%s"
              % (campo, len(comuns), escala, dmax, 100 * rel, onde, extra))

    print("\nPIOR: %.3g%% da escala" % (100 * pior))
    return 0


if __name__ == "__main__":
    sys.exit(main())
