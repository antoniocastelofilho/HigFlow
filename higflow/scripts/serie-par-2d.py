#!/usr/bin/env python3
# Serie CASADA dos dois metodos no caso 2D: para cada instante, a area e o
# CENTROIDE DE AREA de cada interface.
#
#   FT  -> despejo lagrangeano da frente (a representacao nativa do metodo)
#   VOF -> contorno FracVol=0,5 por marching squares sobre o quadtree
#
# O CENTROIDE E' O DE AREA, nao a media dos vertices.  O benchmark define
# y_c = int y dA / int dA, e a media de vertices e' outra coisa: ela pondera
# por DENSIDADE DE MARCADOR, nao por area.  Medido na serie do VOF em t=3, a
# media de vertices erra 2,4e-3 contra 3e-6 do centroide -- 800 vezes.
#
# O PERIMETRO sai aqui para o FT e NAO para o VOF, de proposito: o contorno do
# VOF vem de valores de celula medianizados nos vertices, e essa suavizacao
# enviesa o comprimento.  Area e centroide sao integrais da regiao fechada e
# nao sofrem o mesmo; perimetro sofre.  E' por isso que a circularidade do lado
# VOF continua sem ser reportada.
import sys, os, glob, re, math

def le_xy(p):
    v=[]
    for l in open(p):
        l=l.strip()
        if not l or l[0].isalpha() or l.startswith('#'): continue
        c=l.split()
        if len(c)==2 and c[0].isdigit() and '.' not in c[0]:   # contagem + ds_alvo
            continue
        if len(c)>=2:
            try: v.append((float(c[0]),float(c[1])))
            except ValueError: pass
    return v

def medidas(v):
    n=len(v); A2=0.0; cx=0.0; cy=0.0; P=0.0
    for i in range(n):
        x0,y0=v[i]; x1,y1=v[(i+1)%n]
        cr=x0*y1-x1*y0
        A2+=cr; cx+=(x0+x1)*cr; cy+=(y0+y1)*cr
        P+=math.hypot(x1-x0,y1-y0)
    A=0.5*A2
    if abs(A)<1e-15: return None
    return abs(A), cx/(6*A), cy/(6*A), P, n

if __name__=="__main__":
    dir_ft, dir_vof, dt, dtp = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    fts={}
    for f in glob.glob(os.path.join(dir_ft,"frente_*.dat")):
        passo=int(re.search(r"frente_(\d+)",f).group(1)); fts[round(passo*dt,6)]=f
    vofs={}
    for f in glob.glob(os.path.join(dir_vof,"c_*.dat")):
        q=int(re.search(r"c_(\d+)",f).group(1)); vofs[round(q*dtp,6)]=f
    comuns=sorted(set(fts)&set(vofs))
    print("t areaft areavof difarea ycft ycvof difyc perimft nft nvof")
    for t in comuns:
        a=medidas(le_xy(fts[t])); b=medidas(le_xy(vofs[t]))
        if a is None or b is None: continue
        print("%.4f %.8f %.8f %.2e %.6f %.6f %.2e %.6f %d %d" %
              (t, a[0], b[0], abs(a[0]-b[0]), a[2], b[2], abs(a[2]-b[2]), a[3], a[4], b[4]))
