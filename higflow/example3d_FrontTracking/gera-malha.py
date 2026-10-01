#!/usr/bin/env python3
"""Malha cubica uniforme [0,L]^3 com N celulas por direcao, mais as seis
fronteiras.  Formato do .amr, inferido dos arquivos do example3d_lid_driven:
  linha 1: xlo xhi ylo yhi zlo zhi
  linha 2: numero de blocos
  linha 3: hx hy hz 1
  linha 4: 1 1 1 nx ny nz
"""
import sys, os
# N e' o numero de celulas POR UNIDADE DE COMPRIMENTO, nao por direcao: a caixa
# do artigo de referencia e' 1x1x2, e as celulas tem de ser CUBICAS nela.  Os
# roteiros antigos chamam com (N, L) e um cubo L^3 -- esses continuam valendo,
# porque Ly e Lz caem para Lx quando nao sao dados.
N  = int(sys.argv[1])   if len(sys.argv) > 1 else 20
Lx = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
Ly = float(sys.argv[3]) if len(sys.argv) > 3 else Lx
Lz = float(sys.argv[4]) if len(sys.argv) > 4 else Lx
h  = Lx / N
nx, ny, nz = N, int(round(Ly/h)), int(round(Lz/h))
os.makedirs("amrs/gota_d", exist_ok=True)
os.makedirs("amrs/gota_bc", exist_ok=True)

def escreve(caminho, caixa, hh, nn):
    with open(caminho, "w") as f:
        f.write(" ".join(f"{v:.10f}" for v in caixa) + " \n1\n")
        f.write(" ".join(f"{v:.10f}" for v in hh) + " 1\n")
        f.write("1 1 1 " + " ".join(str(v) for v in nn) + " \n")

escreve("amrs/gota_d/dominio-0.amr", (0,Lx,0,Ly,0,Lz), (h,h,h), (nx,ny,nz))
# seis faces: para cada eixo, o plano lo e o plano hi
faces = [
    ((0,0,   0,Ly,  0,Lz), (0,h,h), (1,ny,nz)),   # x = 0
    ((Lx,Lx, 0,Ly,  0,Lz), (0,h,h), (1,ny,nz)),   # x = Lx
    ((0,Lx,  0,0,   0,Lz), (h,0,h), (nx,1,nz)),   # y = 0
    ((0,Lx,  Ly,Ly, 0,Lz), (h,0,h), (nx,1,nz)),   # y = Ly
    ((0,Lx,  0,Ly,  0,0),  (h,h,0), (nx,ny,1)),   # z = 0
    ((0,Lx,  0,Ly,  Lz,Lz),(h,h,0), (nx,ny,1)),   # z = Lz
]
for i, (caixa, hh, nn) in enumerate(faces):
    escreve(f"amrs/gota_bc/bc-{i}.amr", caixa, hh, nn)
print(f"caixa {Lx}x{Ly}x{Lz}, h={h:.6f}, {nx}x{ny}x{nz} = {nx*ny*nz} celulas")
