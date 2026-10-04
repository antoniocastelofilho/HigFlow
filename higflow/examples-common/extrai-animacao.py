#!/usr/bin/env python3
"""Extrai, por quadro, o contorno da interface e as linhas de corrente no plano
vertical pelo centro da bolha -- para animacao.

O QUE E' POSSIVEL, E O QUE NAO.  A interface sai de uma serie densa (121 a 164
quadros por corrida).  As LINHAS DE CORRENTE precisam do campo de velocidade no
VOLUME, e esse so' existe densamente nas corridas newtonianas, que foram rodadas
antes de a cadencia de campo ser reduzida; nas viscoelasticas ha' 3 quadros de
campo.  O script extrai o que cada corrida tem e diz o que faltou, em vez de
interpolar entre tres instantes e apresentar isso como animacao.

COMO AS LINHAS SAO CALCULADAS.  O VTK de campo e' um UNSTRUCTURED_GRID com oito
cantos por celula e `vel` nos PONTOS.  Os pontos sao reduzidos a uma rede
regular de passo h (a malha e' uniforme), por media dos que caem em cada no'; a
velocidade e' depois amostrada por interpolacao bilinear no plano.  As linhas
saem de sementes numa grade e sao integradas por Runge-Kutta de segunda ordem,
nos DOIS sentidos, para a linha passar pela semente em vez de comecar nela.

A componente FORA do plano e' descartada: o que se desenha e' a projecao do
escoamento no plano de corte, nao linha de corrente tridimensional.  Dito porque
a diferenca importa -- num escoamento com rotacao fora do plano as duas nao
coincidem, e aqui a simetria axial do caso faz com que coincidam quase.

Uso:
    extrai-animacao.py <dir VTKS> --saida <arq.json> [--cada 4] [--em 0.5]
"""
import sys, os, re, glob, json, math

# ----------------------------------------------------------------- leitura
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

def le_campo_vel(caminho):
    """Pontos e velocidade de um UNSTRUCTURED_GRID.  Le' so' POINTS e o
    VECTORS vel -- o resto do arquivo (CELLS, p, FracVol) e' pulado, porque sao
    dezenas de MB que nao entram no calculo das linhas."""
    pts, vel = [], []
    with open(caminho) as f:
        linha = f.readline()
        while linha:
            l = linha.strip()
            if l.startswith("POINTS"):
                n = int(l.split()[1]); v = []
                while len(v) < 3*n:
                    v.extend(float(x) for x in f.readline().split())
                pts = [v[3*k:3*k+3] for k in range(n)]
            elif l.startswith("VECTORS vel"):
                v = []
                while len(v) < 3*len(pts):
                    s = f.readline()
                    if not s: break
                    v.extend(float(x) for x in s.split())
                vel = [v[3*k:3*k+3] for k in range(len(pts))]
                break
            linha = f.readline()
    return pts, vel

# ------------------------------------------------------- contorno no plano
def contorno(pts, tris, eixo, em):
    outros = [d for d in range(3) if d != eixo]
    segs = []
    for t in tris:
        p = [pts[t[0]], pts[t[1]], pts[t[2]]]
        s = [q[eixo] - em for q in p]
        cor = []
        for a, b in ((0,1),(1,2),(2,0)):
            if (s[a] >= 0) == (s[b] >= 0): continue
            if s[b] == s[a]: continue
            u = s[a]/(s[a]-s[b])
            cor.append(tuple(p[a][d] + u*(p[b][d]-p[a][d]) for d in outros))
        if len(cor) >= 2 and cor[0] != cor[1]:
            segs.append((cor[0], cor[1]))
    # encadeia (mesma logica do secao-interface.py)
    tol = 1e-9
    ch = lambda p: (round(p[0]/tol), round(p[1]/tol))
    viz, pos = {}, {}
    for a, b in segs:
        ka, kb = ch(a), ch(b)
        if ka == kb: continue
        viz.setdefault(ka, []).append((kb, b)); viz.setdefault(kb, []).append((ka, a))
        pos[ka] = a; pos[kb] = b
    lac, vis = [], set()
    for k0 in viz:
        if k0 in vis or len(viz[k0]) < 2: continue
        L, k, ant = [pos[k0]], k0, None; vis.add(k0)
        while True:
            nx = None
            for kk, pp in viz.get(k, []):
                if kk != ant and kk not in vis: nx = (kk, pp); break
            if nx is None: break
            ant, k = k, nx[0]; vis.add(k); L.append(nx[1])
            if len(L) > 4*len(segs)+8: break
        if len(L) > 3: lac.append(L)
    lac.sort(key=len, reverse=True)
    return lac

# ------------------------------------------------- rede regular e as linhas
class Campo:
    """Velocidade do plano numa rede regular, para amostrar barato."""
    def __init__(self, pts, vel, eixo, em, h):
        o = [d for d in range(3) if d != eixo]
        # so' os pontos perto do plano entram, com tolerancia de meia celula
        acc, cnt = {}, {}
        self.lo = [1e30, 1e30]; self.hi = [-1e30, -1e30]
        for p, v in zip(pts, vel):
            if abs(p[eixo] - em) > 0.51*h: continue
            k = (round(p[o[0]]/h), round(p[o[1]]/h))
            a = acc.get(k, (0.0, 0.0)); c = cnt.get(k, 0)
            acc[k] = (a[0] + v[o[0]], a[1] + v[o[1]]); cnt[k] = c + 1
            for q in (0, 1):
                self.lo[q] = min(self.lo[q], p[o[q]]); self.hi[q] = max(self.hi[q], p[o[q]])
        self.h = h
        self.g = {k: (a[0]/cnt[k], a[1]/cnt[k]) for k, a in acc.items()}

    def em(self, x, y):
        """Bilinear na rede; (0,0) fora do dominio amostrado."""
        h = self.h
        i, j = math.floor(x/h), math.floor(y/h)
        fx, fy = x/h - i, y/h - j
        s = [0.0, 0.0]
        peso = 0.0
        for di, wx in ((0, 1-fx), (1, fx)):
            for dj, wy in ((0, 1-fy), (1, fy)):
                v = self.g.get((i+di, j+dj))
                if v is None: continue
                w = wx*wy
                s[0] += w*v[0]; s[1] += w*v[1]; peso += w
        if peso < 1e-9: return (0.0, 0.0)
        return (s[0]/peso, s[1]/peso)

def linha_corrente(campo, x0, y0, passo, n, sentido):
    """Runge-Kutta 2 no campo projetado.  Para quando sai do dominio amostrado
    ou quando a velocidade e' residual -- continuar ali desenharia uma linha
    que o campo nao sustenta."""
    out = [(x0, y0)]
    x, y = x0, y0
    for _ in range(n):
        u, v = campo.em(x, y)
        m = math.hypot(u, v)
        if m < 1e-8: break
        u, v = sentido*u/m, sentido*v/m
        xm, ym = x + 0.5*passo*u, y + 0.5*passo*v
        u2, v2 = campo.em(xm, ym)
        m2 = math.hypot(u2, v2)
        if m2 < 1e-8: break
        x += passo*u2/m2; y += passo*v2/m2
        if not (campo.lo[0] <= x <= campo.hi[0] and campo.lo[1] <= y <= campo.hi[1]): break
        out.append((x, y))
    return out

# -------------------------------------------------------------------- main
def quadro_de(n):
    m = re.search(r'-(\d+)\.vtk$', n); return int(m.group(1)) if m else -1

def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    d = sys.argv[1]
    saida, cada, em, eixo = None, 4, 0.5, 1
    for q, a in enumerate(sys.argv):
        if a == "--saida" and q+1 < len(sys.argv): saida = sys.argv[q+1]
        if a == "--cada"  and q+1 < len(sys.argv): cada = int(sys.argv[q+1])
        if a == "--em"    and q+1 < len(sys.argv): em = float(sys.argv[q+1])

    inter = sorted(glob.glob(os.path.join(d, "*_frente_*.vtk")) +
                   glob.glob(os.path.join(d, "*_mult_mc_*.vtk")), key=quadro_de)
    campos = {quadro_de(p): p for p in glob.glob(os.path.join(d, "*print_0-*.vtk"))
              if "_frente_" not in p and "_mult_mc_" not in p}
    if not inter:
        print("  sem serie de interface em %s" % d); return 1

    # h da malha: do nome do diretorio nao se sabe; estima-se do campo
    h = None
    if campos:
        p0 = campos[min(campos)]
        pts, _ = le_campo_vel(p0)
        xs = sorted(set(round(p[0], 6) for p in pts))
        if len(xs) > 2:
            h = min(b-a for a, b in zip(xs, xs[1:]) if b-a > 1e-9)

    quadros, com_linhas = [], 0
    for p in inter[::cada]:
        q = quadro_de(p)
        pts, tris = le_polydata(p)
        lac = contorno(pts, tris, eixo, em)
        reg = {"q": q, "c": [[[round(x,5), round(z,5)] for x, z in L] for L in lac]}
        if q in campos and h:
            cpts, cvel = le_campo_vel(campos[q])
            if cvel:
                cam = Campo(cpts, cvel, eixo, em, h)
                ls = []
                nx, nz = 7, 11
                for ix in range(nx):
                    for iz in range(nz):
                        x0 = cam.lo[0] + (ix+0.5)*(cam.hi[0]-cam.lo[0])/nx
                        z0 = cam.lo[1] + (iz+0.5)*(cam.hi[1]-cam.lo[1])/nz
                        a = linha_corrente(cam, x0, z0, 0.6*h, 90, +1)
                        b = linha_corrente(cam, x0, z0, 0.6*h, 90, -1)
                        L = b[::-1][:-1] + a
                        if len(L) > 8:
                            ls.append([[round(x,5), round(z,5)] for x, z in L[::2]])
                reg["l"] = ls
                com_linhas += 1
        quadros.append(reg)
        if len(quadros) % 10 == 0:
            print("    %d quadros..." % len(quadros), flush=True)

    print("  %s: %d quadros, %d com linhas de corrente%s"
          % (os.path.basename(os.path.dirname(d) or d), len(quadros), com_linhas,
             "" if com_linhas == len(quadros) else "  (os demais SEM campo)"))
    if saida:
        json.dump({"h": h, "em": em, "quadros": quadros}, open(saida, "w"))
        print("    -> %s (%.1f KB)" % (saida, os.path.getsize(saida)/1024))
    return 0

if __name__ == "__main__":
    sys.exit(main())
