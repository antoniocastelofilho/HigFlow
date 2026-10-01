#!/usr/bin/env python3
"""Malha cubica uniforme [0,L]^3 com N celulas por direcao, mais as seis
fronteiras.  Formato do .amr, inferido dos arquivos do example3d_lid_driven:
  linha 1: xlo xhi ylo yhi zlo zhi
  linha 2: numero de blocos
  linha 3: hx hy hz 1
  linha 4: 1 1 1 nx ny nz
"""
import sys, os
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
L = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
h = L / N
os.makedirs("amrs/gota_d", exist_ok=True)
os.makedirs("amrs/gota_bc", exist_ok=True)

def escreve(caminho, caixa, hh, nn):
    with open(caminho, "w") as f:
        f.write(" ".join(f"{v:.10f}" for v in caixa) + " \n1\n")
        f.write(" ".join(f"{v:.10f}" for v in hh) + " 1\n")
        f.write("1 1 1 " + " ".join(str(v) for v in nn) + " \n")

escreve("amrs/gota_d/dominio-0.amr", (0,L,0,L,0,L), (h,h,h), (N,N,N))
# seis faces: para cada eixo, o plano lo e o plano hi
faces = [
    ((0,0, 0,L, 0,L), (0,h,h), (1,N,N)),   # x = 0
    ((L,L, 0,L, 0,L), (0,h,h), (1,N,N)),   # x = L
    ((0,L, 0,0, 0,L), (h,0,h), (N,1,N)),   # y = 0
    ((0,L, L,L, 0,L), (h,0,h), (N,1,N)),   # y = L
    ((0,L, 0,L, 0,0), (h,h,0), (N,N,1)),   # z = 0
    ((0,L, 0,L, L,L), (h,h,0), (N,N,1)),   # z = L
]
for i, (caixa, hh, nn) in enumerate(faces):
    escreve(f"amrs/gota_bc/bc-{i}.amr", caixa, hh, nn)
print(f"malha {N}^3, h={h:.6f}, dominio [0,{L}]^3")
