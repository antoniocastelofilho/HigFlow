# Extrai a INTERFACE do VOF de um VTK do HiGFlow como CURVA, para poder ser
# comparada com a polilinha do front-tracking -- mesmo tipo de objeto dos dois
# lados, e nao uma curva contra uma mancha de celulas.
#
# O VTK traz FracVol por CELULA (CELL_DATA).  Marching squares precisa de valor
# nos VERTICES, entao os valores de celula sao medianizados nos vertices que elas
# compartilham.  Isso SUAVIZA: o contorno resultante e' uma aproximacao da
# interface reconstruida, e nao a reconstrucao PLIC do solver -- que nao esta'
# exposta.  Serve para desenhar a forma; NAO serve para medir perimetro, e e' por
# isso que a circularidade do lado VOF continua nao sendo reportada.
import sys
from collections import defaultdict

def le_vtk(caminho):
    with open(caminho) as f:
        tok = f.read().split('\n')
    i = 0
    pontos, celulas, frac = [], [], []
    while i < len(tok):
        l = tok[i]
        if l.startswith('POINTS'):
            n = int(l.split()[1]); i += 1
            while len(pontos) < n:
                c = tok[i].split()
                for k in range(0, len(c), 3):
                    pontos.append((float(c[k]), float(c[k+1])))
                i += 1
            continue
        if l.startswith('CELLS'):
            n = int(l.split()[1]); i += 1
            while len(celulas) < n:
                c = [int(v) for v in tok[i].split()]
                if c: celulas.append(c[1:1+c[0]])
                i += 1
            continue
        if l.startswith('SCALARS FracVol'):
            i += 2                      # pula LOOKUP_TABLE
            while len(frac) < len(celulas):
                c = tok[i].split()
                frac.extend(float(v) for v in c)
                i += 1
            continue
        i += 1
    return pontos, celulas, frac

def valores_nos_vertices(pontos, celulas, frac):
    # Os pontos vem duplicados por celula; a chave e' a coordenada arredondada.
    soma, cont = defaultdict(float), defaultdict(int)
    def ch(p): return (round(p[0], 9), round(p[1], 9))
    for cel, v in zip(celulas, frac):
        for idp in cel:
            k = ch(pontos[idp]); soma[k] += v; cont[k] += 1
    return {k: soma[k]/cont[k] for k in soma}, ch

def marching(pontos, celulas, vert, ch, nivel=0.5):
    segs = []
    def interp(a, b, va, vb):
        if abs(vb - va) < 1e-30: return a
        t = (nivel - va) / (vb - va)
        t = min(1.0, max(0.0, t))
        return (a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1]))
    for cel in celulas:
        if len(cel) != 4: continue
        P = [pontos[i] for i in cel]
        V = [vert[ch(p)] for p in P]
        cortes = []
        for k in range(4):
            a, b = P[k], P[(k+1) % 4]
            va, vb = V[k], V[(k+1) % 4]
            if (va - nivel) * (vb - nivel) < 0:
                cortes.append(interp(a, b, va, vb))
        if len(cortes) == 2:
            segs.append((cortes[0], cortes[1]))
        elif len(cortes) == 4:          # celula ambigua: liga em duas partes
            segs.append((cortes[0], cortes[1])); segs.append((cortes[2], cortes[3]))
    return segs

def encadeia(segs, tol=1e-9):
    # Junta os segmentos em polilinhas, casando extremidades.
    adj = defaultdict(list)
    def k(p): return (round(p[0]/tol), round(p[1]/tol))
    for s in segs:
        adj[k(s[0])].append(s); adj[k(s[1])].append(s)
    usado, curvas = set(), []
    for s0 in segs:
        if id(s0) in usado: continue
        curva = [s0[0], s0[1]]; usado.add(id(s0))
        for _ in range(2):              # cresce pelas duas pontas
            while True:
                fim = curva[-1]
                prox = None
                for s in adj[k(fim)]:
                    if id(s) in usado: continue
                    prox = s; break
                if prox is None: break
                usado.add(id(prox))
                curva.append(prox[1] if k(prox[0]) == k(fim) else prox[0])
            curva.reverse()
        curvas.append(curva)
    curvas.sort(key=len, reverse=True)
    return curvas

if __name__ == '__main__':
    vtk, saida = sys.argv[1], sys.argv[2]
    pontos, celulas, frac = le_vtk(vtk)
    vert, ch = valores_nos_vertices(pontos, celulas, frac)
    curvas = encadeia(marching(pontos, celulas, vert, ch))
    with open(saida, 'w') as f:
        f.write('x y\n')
        for c in curvas:
            for p in c: f.write(f"{p[0]:.8f} {p[1]:.8f}\n")
            f.write('\n')               # linha vazia: pgfplots quebra o caminho
    tot = sum(len(c) for c in curvas)
    print(f"{vtk}: {len(celulas)} celulas, {len(curvas)} curva(s), "
          f"{tot} pontos, maior com {len(curvas[0]) if curvas else 0}")
