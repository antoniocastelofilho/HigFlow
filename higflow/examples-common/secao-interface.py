#!/usr/bin/env python3
"""Secao transversal da interface: corta a superficie triangulada por um plano
e devolve o contorno fechado, ordenado, pronto para o pgfplots.

POR QUE SECAO E NAO VISTA 3D.  E' o que as publicacoes de referencia usam para
comparar formas (a Figura 35 do artigo sobrepoe os contornos de dois softwares),
e e' o unico jeito de pôr FT e VOF na MESMA figura sem ambiguidade: duas vistas
tridimensionais lado a lado deixam a diferenca para o olho julgar, dois
contornos sobrepostos a mostram.

E nao precisa de renderizador.  O contorno sai da geometria: para cada
triangulo que o plano corta, a intersecao e' um SEGMENTO; os segmentos de uma
superficie fechada formam um ou mais lacos fechados, e o script os encadeia.

Serve aos dois metodos sem distincao, porque os dois escrevem VTK POLYDATA de
triangulos -- o do marching cubes no VOF e o da triangulacao no FT.

ORACULO.  No quadro 0 a superficie e' uma esfera de raio R centrada em c.  O
corte por um plano que passa pelo centro tem de dar um CIRCULO de raio R, e o
script mede |r - R| para que o contorno nao seja aceito sem conferencia.

Uso:
    secao-interface.py <arquivo.vtk> [--eixo y] [--em 0.5] [--dat saida.dat]
    secao-interface.py <dir> --quadro 120 [...]        (pega o quadro do dir)
"""
import sys, os, re, glob, math

# --------------------------------------------------------------------------
def le_polydata(caminho):
    """(vertices, triangulos) de um VTK POLYDATA ASCII.  Mesmo leitor do
    metricas-superficie.py -- se um dia divergirem, a secao e as grandezas
    passariam a descrever superficies diferentes."""
    L = open(caminho).read().split("\n")
    i, pts, tris = 0, [], []
    while i < len(L):
        l = L[i].strip()
        if l.startswith("POINTS"):
            n = int(l.split()[1]); i += 1
            vals = []
            while len(vals) < 3*n:
                vals.extend(float(x) for x in L[i].split()); i += 1
            pts = [vals[3*k:3*k+3] for k in range(n)]
            continue
        if l.startswith("POLYGONS"):
            nt = int(l.split()[1]); i += 1; lidos = 0
            while lidos < nt:
                c = L[i].split(); i += 1
                if not c: continue
                k = int(c[0]); idx = [int(x) for x in c[1:1+k]]
                for q in range(1, k-1):
                    tris.append([idx[0], idx[q], idx[q+1]])
                lidos += 1
            continue
        i += 1
    return pts, tris

# --------------------------------------------------------------------------
def segmentos(pts, tris, eixo, em):
    """Segmentos da intersecao do plano (coord[eixo] = em) com cada triangulo.
    Devolve pares de pontos 2-D nos DOIS eixos restantes, na ordem natural
    (x,y,z) menos `eixo` -- para eixo=1 (y) sai (x,z), que e' o plano vertical
    pelo centro da bolha."""
    outros = [d for d in range(3) if d != eixo]
    segs = []
    for t in tris:
        p = [pts[t[0]], pts[t[1]], pts[t[2]]]
        s = [q[eixo] - em for q in p]
        # o plano corta a aresta (a,b) quando os sinais diferem
        cortes = []
        for a, b in ((0,1), (1,2), (2,0)):
            # Convencao de sinal consistente em ZERO: com `>=`, um extremo
            # exatamente no plano fica de um lado so'.  Importa porque a
            # icosfera centrada em 0,5 tem vertices exatamente em y = 0,5 por
            # simetria, e nesse caso o triangulo e' apenas TOCADO num vertice:
            # as duas arestas incidentes dao o MESMO ponto de corte, o
            # segmento sai degenerado e `encadeia` o descarta.  Isso e'
            # correto, nao perda -- ver a contagem de degenerados que o script
            # imprime.
            if (s[a] >= 0) == (s[b] >= 0):
                continue
            if s[b] == s[a]:
                continue
            u = s[a] / (s[a] - s[b])
            cortes.append(tuple(p[a][d] + u*(p[b][d] - p[a][d]) for d in outros))
        # um triangulo generico da' 0 ou 2 cortes; 3 significa um vertice
        # exatamente no plano, e ai' dois dos cortes coincidem
        if len(cortes) >= 2:
            segs.append((cortes[0], cortes[1]))
    return segs

def encadeia(segs, tol=1e-9):
    """Encadeia os segmentos em lacos fechados.  Os extremos sao arredondados
    numa rede de passo `tol` para que o casamento seja exato em inteiros em vez
    de por comparacao de ponto flutuante -- comparar coordenadas por
    proximidade deixa laco aberto quando dois extremos caem a 2*tol."""
    def chave(p): return (round(p[0]/tol), round(p[1]/tol))
    viz = {}
    for a, b in segs:
        ka, kb = chave(a), chave(b)
        if ka == kb: continue
        viz.setdefault(ka, []).append((kb, b))
        viz.setdefault(kb, []).append((ka, a))
    # ponto de cada chave
    pos = {}
    for a, b in segs:
        pos[chave(a)] = a; pos[chave(b)] = b

    lacos, vistos = [], set()
    for k0 in viz:
        if k0 in vistos: continue
        if len(viz[k0]) < 2: continue
        laco, k, ant = [pos[k0]], k0, None
        vistos.add(k0)
        while True:
            prox = None
            for kk, pp in viz.get(k, []):
                if kk != ant and kk not in vistos:
                    prox = (kk, pp); break
            if prox is None:
                # fecha no inicio?
                break
            ant, k = k, prox[0]
            vistos.add(k); laco.append(prox[1])
            if len(laco) > 10*len(segs) + 10: break
        if len(laco) > 3:
            lacos.append(laco)
    lacos.sort(key=len, reverse=True)
    return lacos

# --------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    alvo = sys.argv[1]
    eixo, em, dat, quadro = 1, None, None, None
    for q, a in enumerate(sys.argv):
        if a == "--eixo" and q+1 < len(sys.argv):
            eixo = {"x":0, "y":1, "z":2}[sys.argv[q+1]]
        if a == "--em"     and q+1 < len(sys.argv): em = float(sys.argv[q+1])
        if a == "--dat"    and q+1 < len(sys.argv): dat = sys.argv[q+1]
        if a == "--quadro" and q+1 < len(sys.argv): quadro = int(sys.argv[q+1])

    if os.path.isdir(alvo):
        padroes = ["*_mult_mc_*-%d.vtk", "*_frente_*-%d.vtk"]
        arq = None
        for p in padroes:
            g = glob.glob(os.path.join(alvo, p % (quadro if quadro is not None else 0)))
            if g: arq = g[0]; break
        if arq is None:
            print("  nenhum quadro %s em %s" % (quadro, alvo)); return 1
    else:
        arq = alvo

    pts, tris = le_polydata(arq)
    if not tris:
        print("  %s: superficie vazia" % os.path.basename(arq)); return 1

    if em is None:   # por omissao, o plano pelo centroide da superficie
        em = sum(p[eixo] for p in pts) / len(pts)

    segs = segmentos(pts, tris, eixo, em)
    lacos = encadeia(segs)
    if not lacos:
        print("  %s: o plano nao corta a superficie em %g" % (arq, em)); return 1

    outros = [d for d in range(3) if d != eixo]
    nomes = "xyz"
    # Os DEGENERADOS entram na contagem, e sem separa-los o par
    # "190 segmentos, laco de 96" parece metade perdida.  Nao e': o plano
    # apenas toca um vertice em 94 dos triangulos, e 96 + 94 = 190.  Imprimir
    # os tres numeros evita a leitura errada -- que eu fiz.
    deg = sum(1 for a, b in segs
              if abs(a[0]-b[0]) < 1e-12 and abs(a[1]-b[1]) < 1e-12)
    print("  %-34s plano %s=%.4f: %d segmento(s) (%d so' tocam vertice), "
          "%d laco(s), maior com %d pontos"
          % (os.path.basename(arq), nomes[eixo], em, len(segs), deg,
             len(lacos), len(lacos[0])))

    if dat:
        with open(dat, "w") as f:
            f.write("# secao da interface: %s\n" % os.path.basename(arq))
            f.write("# plano %s = %.6f;  colunas: %s %s\n"
                    % (nomes[eixo], em, nomes[outros[0]], nomes[outros[1]]))
            f.write("# %d laco(s); lacos separados por linha em branco, que o\n"
                    "# pgfplots entende como quebra de curva\n" % len(lacos))
            f.write("%s %s\n" % (nomes[outros[0]], nomes[outros[1]]))
            for L in lacos:
                for p in L: f.write("%.8f %.8f\n" % (p[0], p[1]))
                p = L[0];   f.write("%.8f %.8f\n\n" % (p[0], p[1]))   # fecha
        print("    -> %s" % dat)
    return 0

if __name__ == "__main__":
    sys.exit(main())
