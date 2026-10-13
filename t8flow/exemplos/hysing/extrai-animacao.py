#!/usr/bin/env python3
"""Extrai, por quadro, a interface e as linhas de corrente dos VTKs do t8Flow.

EMITE O MESMO ESQUEMA JSON do extrator da secao companheira --
`{"h", "em", "quadros": [{"q", "c": contornos, "l": linhas}]}` -- de proposito:
assim o montador de HTML dela serve sem alteracao, e as duas animacoes ficam
COMPARAVEIS lado a lado.  Um esquema proprio exigiria um montador proprio, e
dois montadores divergiriam -- e a diferenca entre as animacoes passaria a
mostrar a divergencia dos desenhistas em vez da dos metodos.

O QUE MUDA E' SO' A LEITURA.  La' o VTK e' `POLYDATA` tridimensional com a
velocidade nos PONTOS; aqui e' `UNSTRUCTURED_GRID` bidimensional com os campos
nas CELULAS, e a interface vem num arquivo proprio, ja' como poligonal.  O
calculo das linhas e' o mesmo: rede regular, interpolacao bilinear, Runge-Kutta
de segunda ordem integrado NOS DOIS SENTIDOS a partir da semente.

Uso:
    extrai-animacao.py <prefixo> <saida.json> [--cada N] [--sementes NxM]
"""

import glob
import json
import math
import os
import re
import sys


def le_malha(caminho):
    """O `UNSTRUCTURED_GRID` que o `SaidaVtk2d` escreve: quatro pontos por
    celula, campos em `CELL_DATA`."""
    L = open(caminho).read().split("\n")
    i, pontos, celulas, dados = 0, [], [], {}
    while i < len(L):
        t = L[i].split()
        if not t:
            i += 1
            continue
        if t[0] == "POINTS":
            n = int(t[1])
            for k in range(n):
                v = L[i + 1 + k].split()
                pontos.append((float(v[0]), float(v[1])))
            i += n
        elif t[0] == "CELLS":
            n = int(t[1])
            for k in range(n):
                v = L[i + 1 + k].split()
                celulas.append([int(x) for x in v[1:]])
            i += n
        elif t[0] in ("SCALARS", "VECTORS"):
            nome = t[1]
            ncomp = 1 if t[0] == "SCALARS" else 3
            salto = 2 if t[0] == "SCALARS" else 1
            vals = []
            for k in range(len(celulas)):
                v = L[i + salto + k].split()
                vals.append([float(x) for x in v[:ncomp]])
            dados[nome] = vals
            i += salto + len(celulas) - 1
        i += 1
    centros = [(sum(pontos[q][0] for q in c) / len(c),
                sum(pontos[q][1] for q in c) / len(c)) for c in celulas]
    return centros, dados


def le_interface(caminho):
    """A poligonal da frente, devolvida como UM contorno fechado em ordem.

    OS SEGMENTOS JA' VEM EM ORDEM no arquivo -- o escritor percorre a frente --,
    mas isso e' reconstruido a partir da conectividade em vez de suposto: se um
    dia a ordem mudar, o contorno sairia embaralhado e o desenho viraria um
    emaranhado, que e' um defeito dificil de atribuir olhando a animacao."""
    if not os.path.exists(caminho):
        return []
    L = open(caminho).read().split("\n")
    i, pts, segs = 0, [], []
    while i < len(L):
        t = L[i].split()
        if t and t[0] == "POINTS":
            n = int(t[1])
            for k in range(n):
                v = L[i + 1 + k].split()
                pts.append((float(v[0]), float(v[1])))
            i += n
        elif t and t[0] == "CELLS":
            n = int(t[1])
            for k in range(n):
                v = L[i + 1 + k].split()
                segs.append((int(v[1]), int(v[2])))
            i += n
        i += 1
    if not segs:
        return []
    prox = {}
    for a, b in segs:
        prox.setdefault(a, []).append(b)
    contornos, usados = [], set()
    for ini in list(prox):
        if ini in usados:
            continue
        cur, cad = ini, []
        while cur is not None and cur not in usados:
            usados.add(cur)
            cad.append(pts[cur])
            viz = [q for q in prox.get(cur, []) if q not in usados]
            cur = viz[0] if viz else None
        if len(cad) > 2:
            contornos.append(cad)
    return contornos


class Campo:
    def __init__(self, centros, dados):
        xs = sorted({round(c[0], 9) for c in centros})
        ys = sorted({round(c[1], 9) for c in centros})
        self.x0, self.y0 = xs[0], ys[0]
        self.hx = (xs[-1] - xs[0]) / (len(xs) - 1) if len(xs) > 1 else 1.0
        self.hy = (ys[-1] - ys[0]) / (len(ys) - 1) if len(ys) > 1 else 1.0
        self.nx, self.ny = len(xs), len(ys)
        self.lo = (xs[0] - self.hx / 2, ys[0] - self.hy / 2)
        self.hi = (xs[-1] + self.hx / 2, ys[-1] + self.hy / 2)
        self.u = [[0.0] * self.nx for _ in range(self.ny)]
        self.v = [[0.0] * self.nx for _ in range(self.ny)]
        vel = dados.get("velocidade")
        for k, (cx, cy) in enumerate(centros):
            i = int(round((cx - self.x0) / self.hx))
            j = int(round((cy - self.y0) / self.hy))
            if vel and 0 <= i < self.nx and 0 <= j < self.ny:
                self.u[j][i], self.v[j][i] = vel[k][0], vel[k][1]

    def em(self, x, y):
        fx, fy = (x - self.x0) / self.hx, (y - self.y0) / self.hy
        i, j = int(math.floor(fx)), int(math.floor(fy))
        if i < 0 or j < 0 or i + 1 >= self.nx or j + 1 >= self.ny:
            return (0.0, 0.0)
        a, b = fx - i, fy - j
        def mix(m):
            return ((1 - a) * (1 - b) * m[j][i] + a * (1 - b) * m[j][i + 1] +
                    (1 - a) * b * m[j + 1][i] + a * b * m[j + 1][i + 1])
        return (mix(self.u), mix(self.v))


def linha(campo, x0, y0, passo, n):
    """RK2, integrada NOS DOIS SENTIDOS: a linha passa PELA semente em vez de
    comecar nela.  So' para a frente, as linhas se acumulam onde o campo
    converge e deixam vazio onde ele diverge -- a densidade do desenho passa a
    mostrar a posicao dos semeadores, nao a estrutura do escoamento."""
    def integra(s):
        x, y, out = x0, y0, []
        for _ in range(n):
            ux, uy = campo.em(x, y)
            m = math.hypot(ux, uy)
            if m < 1e-12:
                break
            xm, ym = x + .5 * s * passo * ux / m, y + .5 * s * passo * uy / m
            vx, vy = campo.em(xm, ym)
            mm = math.hypot(vx, vy)
            if mm < 1e-12:
                break
            x, y = x + s * passo * vx / mm, y + s * passo * vy / mm
            if not (campo.lo[0] <= x <= campo.hi[0] and
                    campo.lo[1] <= y <= campo.hi[1]):
                break
            out.append((x, y))
        return out
    tras = integra(-1)
    tras.reverse()
    return tras + [(x0, y0)] + integra(+1)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    prefixo, saida = sys.argv[1], sys.argv[2]
    cada = 1
    nsx, nsy = 5, 17
    for a in sys.argv[3:]:
        if a.startswith("--cada"):
            cada = int(a.split("=")[-1] if "=" in a else a.split()[-1])
        elif a.startswith("--sementes"):
            v = a.split("=")[-1].split("x")
            nsx, nsy = int(v[0]), int(v[1])

    arqs = sorted(a for a in glob.glob(prefixo + "-[0-9]*-0000.vtk")
                  if "-interface-" not in a)
    if not arqs:
        print("nenhum quadro %s-NNNNNN-0000.vtk -- rode o exemplo com --vtk=N"
              % prefixo)
        return 1

    quadros, h, caixa, sem_campo = [], None, None, 0
    for n, a in enumerate(arqs):
        if n % cada:
            continue
        centros, dados = le_malha(a)
        campo = Campo(centros, dados)
        if h is None:
            h = campo.hx
            caixa = campo.hi[1] - campo.lo[1]
        lx = campo.hi[0] - campo.lo[0]
        ly = campo.hi[1] - campo.lo[1]

        q = int(re.search(r"-(\d{6})-", os.path.basename(a)).group(1))
        face = a.replace(prefixo + "-", prefixo + "-interface-", 1)
        reg = {"q": q,
               "c": [[[round(x, 5), round(y, 5)] for x, y in c]
                     for c in le_interface(face)]}
        if not reg["c"]:
            # SEM INTERFACE o quadro nao tem o que animar, e isso e' DITO em
            # vez de o quadro sair vazio parecendo que a bolha sumiu.
            sem_campo += 1

        if "velocidade" in dados:
            passo = 0.5 * min(campo.hx, campo.hy)
            ls = []
            for jj in range(1, nsy):
                for ii in range(1, nsx):
                    sx = campo.lo[0] + lx * ii / float(nsx)
                    sy = campo.lo[1] + ly * jj / float(nsy)
                    p = linha(campo, sx, sy, passo, 30)
                    if len(p) > 3:
                        ls.append([[round(x, 5), round(y, 5)] for x, y in p])
            reg["l"] = ls
        quadros.append(reg)

    json.dump({"h": h, "em": caixa, "quadros": quadros}, open(saida, "w"))
    # O QUE FALTOU E' DITO, e nao preenchido por interpolacao: uma animacao
    # montada entre tres instantes parece suave e nao e' o escoamento.
    print("%d quadros -> %s (%.2f MB)%s" %
          (len(quadros), saida, os.path.getsize(saida) / 1e6,
           "" if sem_campo == 0 else
           "   ATENCAO: %d quadros SEM interface" % sem_campo))
    return 0


if __name__ == "__main__":
    sys.exit(main())
