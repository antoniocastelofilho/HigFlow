#!/usr/bin/env python3
"""Anima a bolha: campo, LINHAS DE CORRENTE e a interface, em GIF.

POR QUE EM PYTHON PURO, com Pillow e nada mais: esta maquina nao tem numpy nem
matplotlib, e uma animacao que so' roda onde ha' a pilha cientifica instalada
nao serve para olhar o resultado no lugar onde ele e' produzido.  O custo e'
algumas dezenas de linhas de integracao e desenho; o beneficio e' nao haver
dependencia nenhuma entre gerar os dados e ver os dados.

AS LINHAS DE CORRENTE SAO INTEGRADAS AQUI, e nao aproximadas por setas.  Setas
num campo com recirculacao mostram a direcao e escondem a ESTRUTURA -- o par
de vortices na esteira da bolha, que e' o que se quer ver, so' aparece quando
a trajetoria e' seguida.

E ELAS SAO DO CAMPO DO INSTANTE, nao trajetorias de particula: num escoamento
que muda no tempo as duas coisas sao diferentes, e chamar uma pela outra e'
erro de leitura, nao de desenho.
"""

import glob
import math
import os
import re
import sys

from PIL import Image, ImageDraw


def le_vtk(caminho):
    """Le o VTK legado ASCII que o `SaidaVtk2d` escreve."""
    with open(caminho) as f:
        linhas = f.read().split("\n")
    i, n_pontos = 0, 0
    pontos, celulas, dados = [], [], {}
    while i < len(linhas):
        t = linhas[i].split()
        if not t:
            i += 1
            continue
        if t[0] == "POINTS":
            n_pontos = int(t[1])
            for k in range(n_pontos):
                v = linhas[i + 1 + k].split()
                pontos.append((float(v[0]), float(v[1])))
            i += n_pontos
        elif t[0] == "CELLS":
            nc = int(t[1])
            for k in range(nc):
                v = linhas[i + 1 + k].split()
                celulas.append([int(x) for x in v[1:]])
            i += nc
        elif t[0] in ("SCALARS", "VECTORS", "TENSORS"):
            nome, tipo = t[1], t[0]
            ncomp = {"SCALARS": 1, "VECTORS": 3, "TENSORS": 9}[tipo]
            salto = 2 if tipo == "SCALARS" else 1   # LOOKUP_TABLE
            nc = len(celulas)
            vals = []
            for k in range(nc):
                v = linhas[i + salto + k].split()
                vals.append([float(x) for x in v[:ncomp]])
            dados[nome] = vals
            i += salto + nc - 1
        i += 1
    centros = []
    for c in celulas:
        xs = [pontos[q][0] for q in c]
        ys = [pontos[q][1] for q in c]
        centros.append((sum(xs) / len(xs), sum(ys) / len(ys)))
    return centros, celulas, pontos, dados


def le_interface(caminho):
    """Le o arquivo de fronteira: pontos e segmentos."""
    if not os.path.exists(caminho):
        return []
    with open(caminho) as f:
        linhas = f.read().split("\n")
    pontos, segs = [], []
    i = 0
    while i < len(linhas):
        t = linhas[i].split()
        if t and t[0] == "POINTS":
            n = int(t[1])
            for k in range(n):
                v = linhas[i + 1 + k].split()
                pontos.append((float(v[0]), float(v[1])))
            i += n
        elif t and t[0] == "CELLS":
            n = int(t[1])
            for k in range(n):
                v = linhas[i + 1 + k].split()
                segs.append((int(v[1]), int(v[2])))
            i += n
        i += 1
    return [(pontos[a], pontos[b]) for a, b in segs]


class Campo:
    """Velocidade e escalares numa grade uniforme, com busca O(1)."""

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
        self.c = [[0.0] * self.nx for _ in range(self.ny)]
        vel = dados.get("velocidade")
        fra = dados.get("fracao")
        for k, (cx, cy) in enumerate(centros):
            i = int(round((cx - self.x0) / self.hx))
            j = int(round((cy - self.y0) / self.hy))
            if not (0 <= i < self.nx and 0 <= j < self.ny):
                continue
            if vel:
                self.u[j][i] = vel[k][0]
                self.v[j][i] = vel[k][1]
            if fra:
                self.c[j][i] = fra[k][0]

    def em(self, x, y):
        """Velocidade por interpolacao BILINEAR.  Fora do dominio devolve zero,
        o que faz a integracao da linha de corrente parar sozinha."""
        fx = (x - self.x0) / self.hx
        fy = (y - self.y0) / self.hy
        i, j = int(math.floor(fx)), int(math.floor(fy))
        if i < 0 or j < 0 or i + 1 >= self.nx or j + 1 >= self.ny:
            return (0.0, 0.0)
        a, b = fx - i, fy - j
        def mix(m):
            return ((1 - a) * (1 - b) * m[j][i] + a * (1 - b) * m[j][i + 1] +
                    (1 - a) * b * m[j + 1][i] + a * b * m[j + 1][i + 1])
        return (mix(self.u), mix(self.v))


def linha_de_corrente(campo, x0, y0, passo, n):
    """Integra pelo PONTO MEDIO (RK2), para TRAS e para a FRENTE do semeador.

    Euler encurva a trajetoria para fora num campo rotacional, e numa esteira
    com vortice isso vira espiral falsa -- a linha parece sair do vortice
    quando so' o integrador esta' errando.

    E INTEGRA NOS DOIS SENTIDOS porque a linha fica CENTRADA no semeador.  Só
    para a frente, as linhas se acumulam onde o campo converge e deixam vazio
    onde ele diverge: a densidade do desenho passa a mostrar a posicao dos
    semeadores em vez da estrutura do escoamento."""
    def integra(sinal):
        x, y, pts = x0, y0, []
        for _ in range(n):
            ux, uy = campo.em(x, y)
            m = math.hypot(ux, uy)
            if m < 1e-12:
                break
            xm, ym = x + 0.5 * sinal * passo * ux / m, y + 0.5 * sinal * passo * uy / m
            vx, vy = campo.em(xm, ym)
            mm = math.hypot(vx, vy)
            if mm < 1e-12:
                break
            x, y = x + sinal * passo * vx / mm, y + sinal * passo * vy / mm
            if not (campo.lo[0] <= x <= campo.hi[0] and
                    campo.lo[1] <= y <= campo.hi[1]):
                break
            pts.append((x, y))
        return pts
    tras = integra(-1)
    tras.reverse()
    return tras + [(x0, y0)] + integra(+1)


def cor_fracao(c):
    """Azul (fluido externo) a laranja (bolha).  Divergente de proposito: um
    mapa sequencial esconderia a interface no meio da faixa."""
    c = max(0.0, min(1.0, c))
    return (int(28 + 227 * c), int(61 + 100 * c), int(140 - 110 * c))


def desenha(arq_malha, arq_face, larg=360):
    centros, celulas, pontos, dados = le_vtk(arq_malha)
    campo = Campo(centros, dados)
    lx = campo.hi[0] - campo.lo[0]
    ly = campo.hi[1] - campo.lo[1]
    alt = int(larg * ly / lx)
    img = Image.new("RGB", (larg, alt), (255, 255, 255))
    dr = ImageDraw.Draw(img)

    def px(x, y):
        return ((x - campo.lo[0]) / lx * larg, (1 - (y - campo.lo[1]) / ly) * alt)

    # O CAMPO DE FRACAO, celula a celula.
    for k, c in enumerate(celulas):
        xs = [pontos[q][0] for q in c]
        ys = [pontos[q][1] for q in c]
        p0 = px(min(xs), max(ys))
        p1 = px(max(xs), min(ys))
        cx, cy = centros[k]
        i = int(round((cx - campo.x0) / campo.hx))
        j = int(round((cy - campo.y0) / campo.hy))
        v = campo.c[j][i] if (0 <= i < campo.nx and 0 <= j < campo.ny) else 0.0
        dr.rectangle([p0[0], p0[1], p1[0] + 1, p1[1] + 1], fill=cor_fracao(v))

    # AS LINHAS DE CORRENTE, semeadas numa grade regular.  Semear so' perto da
    # bolha mostraria o que se espera ver e esconderia o resto do dominio.
    #
    # O COMPRIMENTO E' LIMITADO de proposito.  Linha longa demais da' a volta
    # no vortice varias vezes e preenche a imagem: o quadro fica saturado e a
    # estrutura some justamente na regiao de maior interesse.  Trinta celulas
    # de cada lado mostram a curvatura local sem fechar o laco.
    passo = 0.5 * min(campo.hx, campo.hy)
    nsub = 30
    for jj in range(1, 17):
        for ii in range(1, 5):
            sx = campo.lo[0] + lx * ii / 5.0
            sy = campo.lo[1] + ly * jj / 17.0
            pts = linha_de_corrente(campo, sx, sy, passo, nsub)
            if len(pts) < 4:
                continue
            dr.line([px(*p) for p in pts], fill=(235, 235, 235), width=1)
            # A SETA no fim diz o SENTIDO, que a linha sozinha nao diz.
            a1, a2 = pts[-2], pts[-1]
            ang = math.atan2(a2[1] - a1[1], a2[0] - a1[0])
            for s2 in (2.4, -2.4):
                q = (a2[0] - 0.9 * passo * math.cos(ang + s2),
                     a2[1] - 0.9 * passo * math.sin(ang + s2))
                dr.line([px(*a2), px(*q)], fill=(235, 235, 235), width=1)

    # A INTERFACE POR CIMA de tudo, grossa: e' o que o metodo de fato
    # transporta, e o resto do quadro e' contexto.
    for a, b in le_interface(arq_face):
        dr.line([px(*a), px(*b)], fill=(20, 20, 20), width=3)

    m = re.search(r"-(\d{6})-", os.path.basename(arq_malha))
    if m:
        dr.text((6, 6), "quadro " + m.group(1), fill=(255, 255, 255))
    return img


def main():
    prefixo = sys.argv[1] if len(sys.argv) > 1 else "hysing"
    saida = sys.argv[2] if len(sys.argv) > 2 else prefixo + ".gif"
    malhas = sorted(glob.glob(prefixo + "-[0-9]*-0000.vtk"))
    malhas = [a for a in malhas if "-interface-" not in a]
    if not malhas:
        print("nenhum quadro %s-NNNNNN-0000.vtk encontrado -- rode o exemplo "
              "com --vtk=N primeiro" % prefixo)
        return 1
    quadros = []
    for a in malhas:
        face = a.replace(prefixo + "-", prefixo + "-interface-", 1)
        quadros.append(desenha(a, face))
        print("  %s" % a)
    quadros[0].save(saida, save_all=True, append_images=quadros[1:],
                    duration=120, loop=0)
    print("%d quadros -> %s" % (len(quadros), saida))
    return 0


if __name__ == "__main__":
    sys.exit(main())
