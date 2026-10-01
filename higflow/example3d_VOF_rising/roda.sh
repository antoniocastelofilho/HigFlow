#!/bin/bash
# LINHA DE BASE EM MALHA UNIFORME.  $1 = celulas por unidade  $2 = numsteps
#
# Existe para duas medidas, nao uma:
#   1. o comportamento -- a bolha sobe, o volume se conserva, nao explode
#   2. o CUSTO por celula e por passo, que e' a referencia contra a qual o
#      ganho do refino adaptativo sera' medido.  Sem ela, "o refino acelerou"
#      nao tem denominador.
M=/home/castelo/HiGFlow-System/HigFlow
cd $M && set -a && . ./varsrc && set +a
cd $M/higflow/example3d_VOF_rising
python3 gera-malha.py ${1:-16} 1.0 1.0 3.0
sed -i "s/^  numsteps: .*/  numsteps: ${2:-2001}/" input/sobe.load.par.contr.yaml
export SB_R=0.25 SB_CX=0.5 SB_CY=0.5 SB_CZ=0.5
export SB_RHO0=1000.0 SB_RHO1=100.0 SB_MU0=10.0 SB_MU1=1.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/usr/bin/time -f "PAREDE %e s" mpirun -use-hwthread-cpus -n 1 ./ns-example-3d \
  input/sobe.load output/s.save VTKS/s.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/sobe-N$1.log 2>&1
echo "N=$1 exit=$? em $(date +%H:%M)" >> /tmp/sobe-status
