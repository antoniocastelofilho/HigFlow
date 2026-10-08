#!/usr/bin/env python3
"""Mede a NITIDEZ DA CUSPIDE na traseira de uma bolha triangulada.

POR QUE ISTO E NAO A ONDULACAO.  O indicador do TSUR responde "o filtro agiria
aqui?".  Ele NAO distingue zigue-zague de curvatura real -- os dois inclinam a
aresta --, e e' justamente por isso que aumentar a cadencia do filtro pode
apagar a cuspide junto com a ondulacao.  Para decidir a cadencia e' preciso uma
medida da FEATURE, nao do defeito.

AS DUAS GRANDEZAS.

1. ASSIMETRIA.  A = (z_c - z_min) / (z_max - z_c).  Numa esfera A = 1; uma
   cuspide na traseira estica a parte de baixo e A cresce.  E' integrada e
   portanto cega a' escala de malha -- nao confunde cuspide com ondulacao.

2. RAIO DE CURVATURA no ponto mais baixo, por ajuste de esfera ao vertice e ao
   seu anel (minimos quadrados).  R pequeno = ponta.  Normalizado pelo raio
   equivalente da bolha, R/R_eq, e' adimensional.  Esta e' LOCAL e portanto
   sensivel a' malha, e por isso vem acompanhada da assimetria: se as duas
   andarem juntas e' cuspide, se so' a local mudar e' malha.

ORACULOS (--oraculo).  Esfera de raio R: ajuste tem de dar R em todo vertice, e
A = 1.  Esferoide x^2/a^2 + y^2/a^2 + z^2/c^2 = 1: no polo as duas curvaturas
principais valem c/a^2, logo o ajuste tem de dar R = a^2/c.  Dois casos com
forma fechada, um deles com R != raio, que e' o que pega um ajuste que so'
devolve a distancia ao centroide.

Uso:  mede-cuspide.py <arquivo.vtk | dir> [--quadro N]
      mede-cuspide.py --oraculo
"""
import sys, os, re, glob, math

# ------------------------------------------------------------------ leitura
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

# ------------------------------------------------------- ajuste de esfera
def ajusta_esfera(P):
    """Centro e raio por minimos quadrados em |x-c|^2 = R^2, que lineariza em
    2*c.x + (R^2 - |c|^2) = |x|^2.  Quatro incognitas, sistema normal 4x4
    resolvido por eliminacao -- sem dependencia externa."""
    n = len(P)
    if n < 4: return None, None
    A = [[2*p[0], 2*p[1], 2*p[2], 1.0] for p in P]
    b = [p[0]**2 + p[1]**2 + p[2]**2 for p in P]
    # sistema normal A^T A x = A^T b
    M = [[sum(A[k][i]*A[k][j] for k in range(n)) for j in range(4)] + 
         [sum(A[k][i]*b[k] for k in range(n))] for i in range(4)]
    # eliminacao de Gauss com pivoteamento parcial
    for c in range(4):
        p = max(range(c,4), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300: return None, None
        M[c], M[p] = M[p], M[c]
        for r in range(4):
            if r == c: continue
            f = M[r][c]/M[c][c]
            for j in range(c,5): M[r][j] -= f*M[c][j]
    x = [M[i][4]/M[i][i] for i in range(4)]
    cen = x[:3]
    r2 = x[3] + cen[0]**2 + cen[1]**2 + cen[2]**2
    return (cen, math.sqrt(r2)) if r2 > 0 else (None, None)

# --------------------------------------------------------------- geometria
def volume_area(pts, tris):
    V = A = 0.0
    for t in tris:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        V += (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
              + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0
        u = [b[i]-a[i] for i in range(3)]; v = [c[i]-a[i] for i in range(3)]
        cr = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
        A += 0.5*math.sqrt(sum(x*x for x in cr))
    return abs(V), A

def aneis(tris):
    inc = {}
    for t in tris:
        for k in range(3):
            inc.setdefault(t[k], set()).update(t[j] for j in range(3) if j != k)
    return inc

def mede(pts, tris):
    V, A = volume_area(pts, tris)
    Req = (3*V/(4*math.pi))**(1/3) if V > 0 else 0.0
    # centroide do VOLUME (nao dos vertices): soma dos tetraedros com a origem
    cz = 0.0; vt = 0.0
    for t in tris:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        dv = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
              + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0
        cz += dv * (a[2]+b[2]+c[2]) / 4.0; vt += dv
    zc = cz/vt if abs(vt) > 0 else 0.0
    zs = [p[2] for p in pts]; zmin, zmax = min(zs), max(zs)
    assim = (zc - zmin)/(zmax - zc) if zmax > zc else float('nan')
    # curvatura no ponto mais baixo
    ib = min(range(len(pts)), key=lambda i: pts[i][2])
    an = aneis(tris)
    viz = sorted(an.get(ib, []))
    cen, R = ajusta_esfera([pts[ib]] + [pts[j] for j in viz]) if len(viz) >= 3 else (None, None)
    return dict(V=V, A=A, Req=Req, zc=zc, zmin=zmin, zmax=zmax, assim=assim,
                R=R, RRel=(R/Req if (R and Req) else float('nan')), nviz=len(viz))

# ----------------------------------------------------------------- oraculos
def icosfera(nsub, raio=1.0, centro=(0,0,0), escala=(1,1,1)):
    t = (1+math.sqrt(5))/2
    V = [(-1,t,0),(1,t,0),(-1,-t,0),(1,-t,0),(0,-1,t),(0,1,t),(0,-1,-t),(0,1,-t),
         (t,0,-1),(t,0,1),(-t,0,-1),(-t,0,1)]
    F = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),
         (10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),
         (2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    V = [list(v) for v in V]
    def nrm(p):
        m = math.sqrt(sum(x*x for x in p)); return [x/m for x in p]
    V = [nrm(v) for v in V]
    cache = {}
    def meio(a, b):
        k = (min(a,b), max(a,b))
        if k not in cache:
            V.append(nrm([(V[a][i]+V[b][i])/2 for i in range(3)])); cache[k] = len(V)-1
        return cache[k]
    for _ in range(nsub):
        F2 = []
        for a,b,c in F:
            ab, bc, ca = meio(a,b), meio(b,c), meio(c,a)
            F2 += [(a,ab,ca),(b,bc,ab),(c,ca,bc),(ab,bc,ca)]
        F = F2
    P = [[centro[i] + raio*escala[i]*v[i] for i in range(3)] for v in V]
    return P, [list(f) for f in F]

def oraculos():
    falhas = 0
    def ok(c, msg, det=""):
        nonlocal falhas
        print("  [%s] %-46s %s" % ("ok" if c else "FALHA", msg, det))
        if not c: falhas += 1
    print("\n-- esfera de raio 0,25: ajuste tem de devolver 0,25 e A = 1 --")
    P, F = icosfera(3, 0.25, (0.5,0.5,0.5))
    m = mede(P, F)
    ok(abs(m['R']-0.25) < 2e-3, "raio de curvatura", "R=%.6f (exato 0.25, erro %.1e)" % (m['R'], abs(m['R']-0.25)))
    ok(abs(m['assim']-1.0) < 1e-6, "assimetria", "A=%.8f" % m['assim'])
    ok(abs(m['Req']-0.25) < 2e-3, "raio equivalente do volume", "Req=%.6f" % m['Req'])
    print("\n-- esferoide a=0,25 c=0,40: no polo R = a^2/c = 0,15625 --")
    a, c = 0.25, 0.40
    P, F = icosfera(3, 1.0, (0.5,0.5,0.5), (a,a,c))
    m = mede(P, F)
    exato = a*a/c
    ok(abs(m['R']-exato) < 0.04*exato, "raio de curvatura no polo",
       "R=%.6f (exato %.6f, erro %.2f%%)" % (m['R'], exato, 100*abs(m['R']-exato)/exato))
    ok(abs(m['assim']-1.0) < 1e-6, "assimetria num esferoide simetrico", "A=%.8f" % m['assim'])
    ok(abs(m['R']-m['Req']) > 0.2*m['Req'], "R difere do raio equivalente",
       "R=%.4f Req=%.4f -- um ajuste que devolvesse a distancia ao centro falharia aqui" % (m['R'], m['Req']))
    print("\n-- esferoide ACHATADO a=0,40 c=0,25: no polo R = a^2/c = 0,64 --")
    a, c = 0.40, 0.25
    P, F = icosfera(3, 1.0, (0.5,0.5,0.5), (a,a,c))
    m = mede(P, F)
    exato = a*a/c
    ok(abs(m['R']-exato) < 0.06*exato, "raio de curvatura no polo",
       "R=%.6f (exato %.6f, erro %.2f%%)" % (m['R'], exato, 100*abs(m['R']-exato)/exato))
    print("\n%s  (%d falha(s))" % ("PASSOU" if falhas == 0 else "FALHOU", falhas))
    return falhas

# --------------------------------------------------------------------- main
def quadro_de(n):
    m = re.search(r'-(\d+)\.vtk$', n); return int(m.group(1)) if m else -1

def main():
    if "--oraculo" in sys.argv:
        return 1 if oraculos() else 0
    if len(sys.argv) < 2: print(__doc__); return 1
    alvo = sys.argv[1]; quadro = None
    for q, a in enumerate(sys.argv):
        if a == "--quadro" and q+1 < len(sys.argv): quadro = int(sys.argv[q+1])
    arqs = []
    if os.path.isdir(alvo):
        g = sorted(glob.glob(os.path.join(alvo, "*_frente_*.vtk")) +
                   glob.glob(os.path.join(alvo, "*_mult_mc_*.vtk")), key=quadro_de)
        arqs = g if quadro is None else [p for p in g if quadro_de(p) == quadro]
    else: arqs = [alvo]
    if not arqs: print("  nada em %s" % alvo); return 1
    print("  %-7s %9s %9s %9s %9s %9s %7s" %
          ("quadro", "z_c", "z_min", "assim.", "R", "R/Req", "viz"))
    for p in arqs:
        pts, tris = le_polydata(p)
        if not tris: continue
        m = mede(pts, tris)
        print("  %-7d %9.5f %9.5f %9.5f %9.5f %9.5f %7d" %
              (quadro_de(p), m['zc'], m['zmin'], m['assim'],
               m['R'] or float('nan'), m['RRel'], m['nviz']))
    return 0

if __name__ == "__main__":
    sys.exit(main())
