#!/bin/bash
# BOLHA ASCENDENTE 3D -- os dois casos de Adelsberger et al., como reproduzidos
# em Silva, Fernandes, Organista, Souza & Castelo, "Validation of HiG-Flow
# Software for Simulating Two-Phase Flows with a 3D Geometric Volume of Fluid
# Algorithm", Mathematics 2023, 11(18), 3900, Tabela 3.
#
#   $1 = celulas por unidade de comprimento (50 reproduz a malha 50x50x100)
#   $2 = numsteps (1501 da' t = 3.0 com dt = 2e-3, o tempo final do artigo)
#   $3 = caso, 1 ou 2
#
# PARAMETROS DO ARTIGO (Tabela 3):
#   caso   rho1  rho2  mu1   mu2   g     sigma  Re  Bo   regime
#   1      1000  100   10    1     0.98  24.5   35  10   elipsoidal
#   2      1000  1     10    0.1   0.98  1.96   35  125  saiote/calota
#
# Dominio [0,1]x[0,1]x[0,2], bolha R = 0.25 em (0.5, 0.5, 0.5), nao-deslizamento
# nas paredes, t final 3.0, malha 50x50x100, dt = 2e-3.
#
# A CONVENCAO DOS ADIMENSIONAIS NESTA ARVORE, lida do example2d_VOF (o caso
# validado) e nao suposta: trabalha-se em variaveis DIMENSIONAIS, com
#
#   Re = 1          as densidades e viscosidades entram com os valores de tabela
#   Fr = 1/sqrt(g)  porque higflow_gravity_term devolve 1/Fr^2
#   Ca = 1/sigma    porque a tensao interfacial entra como IF/(Ca*Re*dens)
#
# Conferido: Fr = 1.010153 da' 1/Fr^2 = 0.979999, e Ca = 0.0408163265 da'
# 1/(Ca*Re) = 24.500000.  O comentario "(not used yet)" ao lado do Fr nos yamls
# e' FALSO -- o termo e' aplicado em 5 sitios do caminho multifasico.
#
# Grandezas de comparacao do artigo: centro de massa e velocidade de subida da
# bolha, contra OpenFOAM, DROPS, NaSt3D e FEEL (Metivet et al.).
M=/home/castelo/HiGFlow-System/HigFlow
cd $M && set -a && . ./varsrc && set +a
cd $M/higflow/example3d_VOF_rising

N=${1:-50}; PASSOS=${2:-1501}; CASO=${3:-1}

python3 gera-malha.py $N 1.0 1.0 2.0
sed -i "s/^  numsteps: .*/  numsteps: ${PASSOS}/" input/sobe.load.par.contr.yaml

if [ "$CASO" = "2" ]; then
    export SB_RHO0=1000.0 SB_RHO1=1.0 SB_MU0=10.0 SB_MU1=0.1
    CA=0.5102040816                      # 1/1.96
else
    export SB_RHO0=1000.0 SB_RHO1=100.0 SB_MU0=10.0 SB_MU1=1.0
    CA=0.0408163265                      # 1/24.5
fi
sed -i "s/^    Ca: [0-9.]*/    Ca: ${CA}/" input/sobe.load.par.contr.yaml

export SB_R=0.25 SB_CX=0.5 SB_CY=0.5 SB_CZ=0.5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/usr/bin/time -f "PAREDE %e s" mpirun -use-hwthread-cpus -n 1 ./ns-example-3d \
  input/sobe.load output/s.save VTKS/s.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  > /tmp/sobe-c${CASO}-N${N}.log 2>&1
echo "caso=$CASO N=$N exit=$? em $(date +%H:%M)" >> /tmp/sobe-status
