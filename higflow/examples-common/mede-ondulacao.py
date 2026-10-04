#!/usr/bin/env python3
"""Mede a ONDULACAO SUBMALHA de uma superficie triangulada, pelo mesmo
indicador que o TSUR usa como porta.

POR QUE ESTE INDICADOR.  Nos testes sinteticos a ondulacao foi medida como
desvio do raio em torno da media -- o que so' serve numa esfera.  Numa bolha
deformada nao ha' raio de referencia, e usar um mediria a DEFORMACAO, nao a
ondulacao.  O indicador do TSUR nao tem esse problema: para cada aresta, com
`ne` a normal da aresta (media das normais dos dois aneis) e `meio` o ponto
medio,

    h1 = (v1 - meio) . ne

mede o quanto a aresta esta' INCLINADA em relacao a' superficie local.  Numa
superficie lisa h1 e' residual, qualquer que seja a forma; num zigue-zague de
escala de malha ele e' da ordem da meia-aresta.  Normalizado pela meia-aresta,
o indicador e' adimensional e comparavel entre malhas.

E' o MESMO numero que `ft3_tsur` compara com o limiar, de modo que a medida
responde diretamente "o filtro teria agido aqui?".

Uso:  mede-ondulacao.py <arquivo.vtk | dir> [--quadro N] [--limiar 0.05]
"""
import sys, os, re, glob, math

def le_polydata(caminho):
    L = open(caminho).read().split("\n")
    i, pts, tris = 0, [], []
    while i < len(L):
        l = L[i].strip()
        if l.startswith("POINTS"):
            n = int(l.split()[1]); i += 1
            v = []
            while len(v) < 3*n:
                v.extend(float(x) for x in L[i].split()); i += 1
            pts = [v[3*k:3*k+3] for k in range(n)]
            continue
        if l.startswith("POLYGONS"):
            nt = int(l.split()[1]); i += 1; lid = 0
            while lid < nt:
                c = L[i].split(); i += 1
                if not c: continue
                k = int(c[0]); idx = [int(x) for x in c[1:1+k]]
                for q in range(1, k-1): tris.append([idx[0], idx[q], idx[q+1]])
                lid += 1
            continue
        i += 1
    return pts, tris

def sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def cruz(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def pto(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def nor(a):
    m = math.sqrt(pto(a,a))
    return (a[0]/m, a[1]/m, a[2]/m) if m > 1e-300 else (0.0,0.0,0.0)

def aneis(pts, tris):
    """Anel ordenado de cada vertice, da incidencia -- a mesma construcao do
    _anel do modulo em C."""
    inc = {}
    for t in tris:
        for k in range(3): inc.setdefault(t[k], []).append(t)
    anel = {}
    for v, ts in inc.items():
        de, pa = {}, set()
        for t in ts:
            c = t.index(v)
            a, b = t[(c+1)%3], t[(c+2)%3]
            de[a] = b; pa.add(a)
        if len(de) < 3: continue
        comeco = next(iter(de)); r = [comeco]; cur = de[comeco]
        while cur != comeco and cur in de and len(r) <= len(de):
            r.append(cur); cur = de[cur]
        if cur == comeco and len(r) == len(de): anel[v] = r
    return anel

def normal_anel(pts, anel_v):
    p = [0.0,0.0,0.0]
    for q in anel_v:
        for d in range(3): p[d] += pts[q][d]
    n = len(anel_v)
    p = [c/n for c in p]
    s = [0.0,0.0,0.0]
    for i in range(n):
        a = sub(pts[anel_v[i]], p); b = sub(pts[anel_v[(i+1)%n]], p)
        c = cruz(a,b)
        for d in range(3): s[d] += c[d]
    return p, nor(tuple(s))

def ondulacao(pts, tris):
    """Devolve a lista (indicador, z do ponto medio) por aresta."""
    an = aneis(pts, tris)
    arestas = set()
    for t in tris:
        for k in range(3):
            a, b = t[k], t[(k+1)%3]
            arestas.add((a,b) if a < b else (b,a))
    out = []
    for v1, v2 in arestas:
        if v1 not in an or v2 not in an: continue
        _, n1 = normal_anel(pts, an[v1])
        _, n2 = normal_anel(pts, an[v2])
        ne = nor((n1[0]+n2[0], n1[1]+n2[1], n1[2]+n2[2]))
        if ne == (0.0,0.0,0.0): continue
        meio = [0.5*(pts[v1][d]+pts[v2][d]) for d in range(3)]
        h1 = pto(sub(pts[v1], meio), ne)
        meia = 0.5*math.sqrt(pto(sub(pts[v2],pts[v1]), sub(pts[v2],pts[v1])))
        if meia <= 0: continue
        out.append((abs(h1)/meia, meio[2]))
    return out

def quadro_de(n):
    m = re.search(r'-(\d+)\.vtk$', n); return int(m.group(1)) if m else -1

def main():
    if len(sys.argv) < 2: print(__doc__); return 1
    alvo = sys.argv[1]; lim = 0.05; quadro = None
    for q,a in enumerate(sys.argv):
        if a == "--limiar" and q+1 < len(sys.argv): lim = float(sys.argv[q+1])
        if a == "--quadro" and q+1 < len(sys.argv): quadro = int(sys.argv[q+1])
    if os.path.isdir(alvo):
        g = sorted(glob.glob(os.path.join(alvo,"*_frente_*.vtk")) +
                   glob.glob(os.path.join(alvo,"*_mult_mc_*.vtk")), key=quadro_de)
        if not g: print("  sem superficie em %s"%alvo); return 1
        arq = g[-1] if quadro is None else [p for p in g if quadro_de(p)==quadro][0]
    else: arq = alvo

    pts, tris = le_polydata(arq)
    if not tris: print("  %s vazio"%arq); return 1
    ind = ondulacao(pts, tris)
    if not ind: print("  %s: nenhum anel fechado"%arq); return 1

    v = sorted(x for x,_ in ind)
    n = len(v)
    acima = sum(1 for x in v if x > lim)
    zs = [z for _,z in ind]
    zmin, zmax = min(zs), max(zs)
    # a METADE INFERIOR da bolha, que e' onde o usuario espera a ondulacao
    corte = 0.5*(zmin+zmax)
    baixo = [x for x,z in ind if z < corte]
    ab_baixo = sum(1 for x in baixo if x > lim)

    print("  %-30s %d arestas" % (os.path.basename(arq), n))
    print("    indicador |h1|/meia-aresta:  mediana %.4f  p90 %.4f  p99 %.4f  max %.4f"
          % (v[n//2], v[int(.90*n)], v[int(.99*n)], v[-1]))
    print("    acima do limiar %.2f:  %d de %d  (%.2f%%)" % (lim, acima, n, 100*acima/n))
    print("    na METADE INFERIOR:    %d de %d  (%.2f%%)   z em [%.3f, %.3f], corte %.3f"
          % (ab_baixo, len(baixo), 100*ab_baixo/max(1,len(baixo)), zmin, zmax, corte))
    return 0

if __name__ == "__main__":
    sys.exit(main())
