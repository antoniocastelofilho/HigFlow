#!/usr/bin/env python3
"""Metricas de forma da bolha a partir da SUPERFICIE, para os dois metodos.

POR QUE UM SCRIPT SO' SERVE AOS DOIS.  O VOF escreve a interface como
isosuperficie de marching cubes (_mult_mc_) e o FT escreve a propria
triangulacao (_frente_), e os dois saem em VTK POLYDATA de triangulos.  Logo a
MESMA conta de volume, area e centroide se aplica aos dois, e a comparacao nao
herda diferenca de instrumento -- que seria o jeito mais facil de inventar uma
diferenca entre os metodos.

GRANDEZAS, as que o artigo de referencia e a tese publicam:

  centro de massa    z do centroide, ponderado por volume (nao media de
                     vertices: a triangulacao e' irregular e a media de
                     vertices pesa mais onde ha' mais triangulo)
  velocidade         d z_c / dt, por diferenca centrada na serie
  esfericidade       Psi = pi^(1/3) (6V)^(2/3) / A,  igual a 1 na esfera
  volume, area       pelo teorema da divergencia sobre os triangulos
  extensao           z minimo e maximo -- e' o que mostra a CUSPIDE: ela
                     aparece como z_min descolando do corpo da bolha

A esfericidade usa V e A da MESMA superficie, logo ela mede a forma do objeto
extraido.  Para o VOF isso inclui o erro da extracao (o marching cubes encerra
~4% menos volume que a integral de fracvol em N=16, 2.3% em N=25 -- cai como
h^2); para o FT inclui o erro da triangulacao.  Sao erros de REPRESENTACAO, que
e' exatamente o que a comparacao quer medir, mas por isso a esfericidade nao
deve ser lida como grandeza do escoamento em malha grossa.

Uso:
    metricas-superficie.py <dir> [--dt-quadro DT] [--csv saida.csv]

Varre <dir> por *_mult_mc_*.vtk e *_frente_*.vtk, ordena por numero de quadro e
imprime a serie.  Com --dt-quadro calcula a velocidade; sem ele, so' as
grandezas de forma.
"""
import sys, os, re, glob, math

def le_polydata(caminho):
    """Devolve (vertices, triangulos) de um VTK POLYDATA ASCII de triangulos."""
    with open(caminho) as f:
        linhas = f.read().split("\n")
    i = 0
    pts, tris = [], []
    while i < len(linhas):
        l = linhas[i].strip()
        if l.startswith("POINTS"):
            n = int(l.split()[1])
            i += 1
            # as coordenadas podem vir varias por linha
            vals = []
            while len(vals) < 3*n:
                vals.extend(float(x) for x in linhas[i].split())
                i += 1
            pts = [vals[3*k:3*k+3] for k in range(n)]
            continue
        if l.startswith("POLYGONS"):
            nt = int(l.split()[1])
            i += 1
            lidos = 0
            while lidos < nt:
                c = linhas[i].split()
                i += 1
                if not c:
                    continue
                k = int(c[0])
                if k == 3:
                    tris.append([int(c[1]), int(c[2]), int(c[3])])
                else:
                    # poligono de k lados -> leque de triangulos
                    idx = [int(x) for x in c[1:1+k]]
                    for q in range(1, k-1):
                        tris.append([idx[0], idx[q], idx[q+1]])
                lidos += 1
            continue
        i += 1
    return pts, tris

def grandezas(pts, tris):
    """Volume, area e centroide de um poliedro FECHADO.

    Decompoe-se o solido nos tetraedros que cada triangulo forma com a ORIGEM.
    O tetraedro (0,a,b,c) tem volume assinado  v = a.(b x c)/6  e centroide
    (0+a+b+c)/4.  Somando,

        V  = sum v                       (teorema da divergencia)
        C  = sum( v * (a+b+c)/4 ) / V

    Os sinais cuidam das partes que ficam do lado de fora quando a origem nao
    esta' dentro do solido, logo a identidade vale sem exigir nada da posicao da
    origem -- so' que a superficie seja fechada e consistentemente orientada.
    E' por isso que `arestas_orfas` e' verificada antes de o numero ser
    impresso: sem fechamento, estas somas nao significam nada.

    NAO se usa a media dos vertices como centroide: a triangulacao e' irregular
    e a media pesaria mais onde ha' mais triangulo, nao onde ha' mais volume.
    """
    V = 0.0
    A = 0.0
    cx = cy = cz = 0.0
    zmin, zmax = 1e300, -1e300
    for t in tris:
        a, b, c = pts[t[0]], pts[t[1]], pts[t[2]]
        # produto misto a.(b x c)
        bxc = (b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0])
        v6 = a[0]*bxc[0] + a[1]*bxc[1] + a[2]*bxc[2]
        v = v6/6.0
        V += v
        cx += v*(a[0]+b[0]+c[0])/4.0
        cy += v*(a[1]+b[1]+c[1])/4.0
        cz += v*(a[2]+b[2]+c[2])/4.0
        # area
        u = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
        w = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
        n = (u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0])
        A += 0.5*math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
        for p in (a, b, c):
            if p[2] < zmin: zmin = p[2]
            if p[2] > zmax: zmax = p[2]
    if V != 0.0:
        cx, cy, cz = cx/V, cy/V, cz/V
    return V, A, (cx, cy, cz), (zmin, zmax)

def esfericidade(V, A):
    """Psi = pi^(1/3) (6V)^(2/3) / A.  Igual a 1 na esfera, menor senao."""
    if A <= 0.0 or V <= 0.0:
        return float("nan")
    return (math.pi**(1.0/3.0)) * ((6.0*V)**(2.0/3.0)) / A

def arestas_orfas(tris):
    """Arestas usadas por um numero de triangulos != 2.  Zero numa superficie
    fechada.  Serve de guarda: se a superficie nao fechar, volume e esfericidade
    nao significam nada, e e' melhor dizer isso que imprimir o numero."""
    cont = {}
    for t in tris:
        for r in range(3):
            a, b = t[r], t[(r+1) % 3]
            k = (a, b) if a < b else (b, a)
            cont[k] = cont.get(k, 0) + 1
    return sum(1 for v in cont.values() if v != 2)

def quadro_de(nome):
    m = re.search(r'-(\d+)\.vtk$', nome)
    return int(m.group(1)) if m else -1

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    d = sys.argv[1]
    dtq = None
    csv = None
    for q, a in enumerate(sys.argv):
        if a == "--dt-quadro" and q+1 < len(sys.argv): dtq = float(sys.argv[q+1])
        if a == "--csv"       and q+1 < len(sys.argv): csv = sys.argv[q+1]

    arqs = sorted(glob.glob(os.path.join(d, "*_mult_mc_*.vtk")) +
                  glob.glob(os.path.join(d, "*_frente_*.vtk")),
                  key=quadro_de)
    if not arqs:
        print("nenhum *_mult_mc_*.vtk nem *_frente_*.vtk em %s" % d)
        return 1

    tipo = "marching cubes (VOF)" if "_mult_mc_" in arqs[0] else "frente (FT)"
    print("\n  %d quadros, %s\n" % (len(arqs), tipo))
    cab = "  %-6s %-8s %12s %12s %10s %10s %10s %8s" % (
          "quadro", "t", "volume", "area", "z_c", "z_min", "z_max", "Psi")
    if dtq: cab += " %10s" % "v_subida"
    print(cab)
    print("  " + "-"*(len(cab)-2))

    serie = []
    for a in arqs:
        pts, tris = le_polydata(a)
        if not tris:
            print("  %-6d %-8s  (vazio)" % (quadro_de(a), "-"))
            continue
        V, Ar, cc, ext = grandezas(pts, tris)
        orf = arestas_orfas(tris)
        t = quadro_de(a)*dtq if dtq else None
        serie.append((quadro_de(a), t, V, Ar, cc[2], ext[0], ext[1],
                      esfericidade(V, Ar), orf))

    for q, (nq, t, V, Ar, zc, z0, z1, psi, orf) in enumerate(serie):
        v = ""
        if dtq:
            if 0 < q < len(serie)-1:
                v = " %10.6f" % ((serie[q+1][4]-serie[q-1][4])/(2.0*dtq))
            elif q == 0 and len(serie) > 1:
                v = " %10.6f" % ((serie[1][4]-serie[0][4])/dtq)
            else:
                v = " %10s" % "-"
        marca = "" if orf == 0 else "   <-- NAO FECHADA (%d arestas orfas)" % orf
        print("  %-6d %-8s %12.8f %12.8f %10.6f %10.6f %10.6f %8.5f%s%s" % (
              nq, ("%.4f" % t) if t is not None else "-", V, Ar, zc, z0, z1, psi, v, marca))

    if serie:
        V0 = serie[0][2]
        Vf = serie[-1][2]
        print("\n  volume: %.8f -> %.8f   deriva %.3e (%.4f%%)" % (
              V0, Vf, abs(Vf-V0)/abs(V0) if V0 else 0.0,
              100*abs(Vf-V0)/abs(V0) if V0 else 0.0))
        nf = sum(1 for s in serie if s[8] != 0)
        print("  superficies nao fechadas: %d de %d" % (nf, len(serie)))

    if csv:
        with open(csv, "w") as f:
            f.write("quadro,t,volume,area,z_c,z_min,z_max,esfericidade,arestas_orfas\n")
            for s in serie:
                f.write("%d,%s,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%d\n" % (
                        s[0], ("%.6f" % s[1]) if s[1] is not None else "", *s[2:]))
        print("  csv: %s" % csv)
    print()
    return 0

if __name__ == "__main__":
    sys.exit(main())
