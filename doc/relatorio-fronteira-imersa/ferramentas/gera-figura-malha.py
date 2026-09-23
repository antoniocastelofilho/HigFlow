#!/usr/bin/env python3
"""Gera figuras/malha-zoom.tex a partir da MESMA construcao que o codigo faz.

Os marcadores sao os centros dos subsegmentos do poligono, como em
`fi_cria_circulo` + `_percorre_curva` de higflow/src/hig-flow-fronteira-imersa.c;
as linhas de malha saem do h de fundo e do h do nivel fino dentro da caixa de
refino.

Mudou o experimento, rode isto de novo -- figura desenhada a mao envelhece
calada, e uma figura que discorda do texto e' pior que figura nenhuma.

    cd doc/relatorio-fronteira-imersa && python3 ferramentas/gera-figura-malha.py
"""
import io, math, os

# --- os parametros do experimento, em um lugar so' ---------------------------
CENTRO   = (2.0, 2.0)
RAIO     = 0.5
DS       = 0.025                      # espacamento alvo dos marcadores
H_FUNDO  = 0.05
H_FINO   = 0.025
CAIXA    = (1.0, 1.0, 3.5, 3.0)       # x0,y0,x1,y1 da regiao refinada
JANELA   = (1.2, 1.2, 2.9, 2.9)       # o zoom da figura
LARGURA  = 7.0                        # cm

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
x0, y0, x1, y1 = JANELA
cxlo, cylo, cxhi, cyhi = CAIXA
cx, cy = CENTRO
esc = LARGURA / (x1 - x0)
X = lambda v: (v - x0) * esc
Y = lambda v: (v - y0) * esc

# --- marcadores: a construcao do codigo --------------------------------------
n = int(round(math.pi / DS))
V = [(cx + RAIO*math.cos(2*math.pi*i/n), cy + RAIO*math.sin(2*math.pi*i/n))
     for i in range(n)]
marc = []
for i in range(n):
    a, b = V[i], V[(i+1) % n]
    comp = math.hypot(b[0]-a[0], b[1]-a[1])
    nsub = max(1, int(comp/DS + 0.5))
    for s in range(nsub):
        t = (s + 0.5)/nsub
        marc.append((a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1])))

io.open('dados/marcadores.dat', 'w').write(
    "# x y   -- %d marcadores, centros dos subsegmentos de %d lados\n" % (len(marc), n)
    + "".join("%.6f %.6f\n" % p for p in marc))

# --- a figura ----------------------------------------------------------------
L = ["%% GERADO por ferramentas/gera-figura-malha.py -- nao editar a mao.",
     "\\begin{tikzpicture}[x=1cm,y=1cm]",
     "\\clip (0.000,0.000) rectangle (%.3f,%.3f);" % (X(x1), Y(y1))]

def grade(passo, cor, larg, xa, xb, ya, yb):
    k = math.ceil(xa/passo)
    while k*passo <= xb + 1e-12:
        c = k*passo
        if x0 <= c <= x1:
            L.append("\\draw[%s,line width=%s] (%.3f,%.3f)--(%.3f,%.3f);"
                     % (cor, larg, X(c), Y(max(ya, y0)), X(c), Y(min(yb, y1))))
        k += 1
    k = math.ceil(ya/passo)
    while k*passo <= yb + 1e-12:
        c = k*passo
        if y0 <= c <= y1:
            L.append("\\draw[%s,line width=%s] (%.3f,%.3f)--(%.3f,%.3f);"
                     % (cor, larg, X(max(xa, x0)), Y(c), X(min(xb, x1)), Y(c)))
        k += 1

grade(H_FUNDO, "gray!45", "0.15pt", x0, x1, y0, y1)
grade(H_FINO,  "blue!30", "0.10pt", max(x0, cxlo), min(x1, cxhi),
                                    max(y0, cylo), min(y1, cyhi))
L.append("\\draw[blue!70!black,line width=0.9pt,dashed] (%.3f,%.3f) rectangle (%.3f,%.3f);"
         % (X(max(cxlo, x0)), Y(max(cylo, y0)), X(min(cxhi, x1)), Y(min(cyhi, y1))))
L.append("\\draw[black,line width=1.0pt] (%.3f,%.3f) circle (%.3fcm);"
         % (X(cx), Y(cy), RAIO*esc))
for p in marc:
    L.append("\\fill[red!80!black] (%.3f,%.3f) circle (0.9pt);" % (X(p[0]), Y(p[1])))
L.append("\\end{tikzpicture}")

io.open('figuras/malha-zoom.tex', 'w').write("\n".join(L) + "\n")
print("malha-zoom.tex: %d marcadores, %d comandos" % (len(marc), len(L)))
