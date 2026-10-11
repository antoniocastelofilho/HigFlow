#!/usr/bin/env python3
# Converte o despejo da frente do front-tracking 2D para o formato `x y` que as
# tabelas e os graficos leem, e mede area e centroide.
#
# O FORMATO tem duas linhas de cabecalho que NAO sao coordenadas:
#
#     # frente front-tracking v1          <- comentario
#     252 0.00623...                      <- numero de marcadores e ds_alvo
#     0.75 0.5                            <- dai' em diante, x y
#
# A segunda e' a armadilha: comeca com digito, entao um filtro ingenuo de
# "linha que nao comeca com letra" a deixa passar e ela entra como ponto.
import sys, math

def le(caminho):
    v = []
    with open(caminho) as f:
        linhas = [l.strip() for l in f if l.strip()]
    if linhas and linhas[0].startswith('#'):
        linhas = linhas[1:]
    if linhas:
        c = linhas[0].split()
        if len(c) == 2 and c[0].isdigit():       # contagem + ds_alvo
            n = int(c[0]); linhas = linhas[1:]
        else:
            n = None
    for l in linhas:
        c = l.split()
        if len(c) >= 2:
            try: v.append((float(c[0]), float(c[1])))
            except ValueError: pass
    if n is not None and len(v) != n:
        print("AVISO: cabecalho diz %d marcadores, li %d" % (n, len(v)), file=sys.stderr)
    return v

def medidas(v):
    n = len(v)
    A = 0.5*abs(sum(v[i][0]*v[(i+1)%n][1] - v[(i+1)%n][0]*v[i][1] for i in range(n)))
    P = sum(math.hypot(v[(i+1)%n][0]-v[i][0], v[(i+1)%n][1]-v[i][1]) for i in range(n))
    cx = sum(x for x,_ in v)/n; cy = sum(y for _,y in v)/n
    return n, A, P, cx, cy

if __name__ == "__main__":
    v = le(sys.argv[1])
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as g:
            g.write("x y\n")
            for x, y in v: g.write("%.8f %.8f\n" % (x, y))
    n, A, P, cx, cy = medidas(v)
    print("%d %.8f %.8f %.6f %.6f" % (n, A, P, cx, cy))
