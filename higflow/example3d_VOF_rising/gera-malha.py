#!/usr/bin/env python3
"""Malha de CAIXA Lx x Ly x Lz com h uniforme, mais as seis fronteiras.

O caso de bolha ascendente 3D usa 1 x 1 x 3 -- a direcao longa e' z, que e' onde
a gravidade age (dim == DIM-1).  E' o analogo tridimensional do Hysing 2D, que
usa 1 x 2 com a bolha subindo em y.
"""
import sys, os
N  = int(sys.argv[1]) if len(sys.argv) > 1 else 16     # celulas por unidade
Lx = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
Ly = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
Lz = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
h  = Lx / N
nx, ny, nz = round(Lx/h), round(Ly/h), round(Lz/h)
os.makedirs("amrs/d",  exist_ok=True)
os.makedirs("amrs/bc", exist_ok=True)

def escreve(caminho, caixa, hh, nn):
    with open(caminho, "w") as f:
        f.write(" ".join(f"{v:.10f}" for v in caixa) + " \n1\n")
        f.write(" ".join(f"{v:.10f}" for v in hh) + " 1\n")
        f.write("1 1 1 " + " ".join(str(v) for v in nn) + " \n")

escreve("amrs/d/dominio-0.amr", (0,Lx,0,Ly,0,Lz), (h,h,h), (nx,ny,nz))
faces = [((0,0, 0,Ly, 0,Lz),(0,h,h),(1,ny,nz)), ((Lx,Lx, 0,Ly, 0,Lz),(0,h,h),(1,ny,nz)),
         ((0,Lx, 0,0, 0,Lz),(h,0,h),(nx,1,nz)), ((0,Lx, Ly,Ly, 0,Lz),(h,0,h),(nx,1,nz)),
         ((0,Lx, 0,Ly, 0,0),(h,h,0),(nx,ny,1)), ((0,Lx, 0,Ly, Lz,Lz),(h,h,0),(nx,ny,1))]
for i,(c,hh,nn) in enumerate(faces):
    escreve(f"amrs/bc/bc-{i}.amr", c, hh, nn)
print(f"caixa {Lx}x{Ly}x{Lz}, h={h:.6f}, {nx}x{ny}x{nz} = {nx*ny*nz} celulas")
