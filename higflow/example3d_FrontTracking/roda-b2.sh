#!/bin/bash
# $1 = N (celulas por direcao)   $2 = nsub (subdivisoes da icosfera)
# O par (N,nsub) tem de manter ds/h ~ 0,75: front-tracking exige espacamento de
# marcador da ordem de h, e refinar so' a malha deixa a frente grossa demais.
cd /home/castelo/HiGFlow-System/HigFlow && set -a && . ./varsrc && set +a
cd higflow/example3d_FrontTracking
python3 gera-malha.py $1 1.0 >/dev/null
export FT3_R=0.25 FT3_SIGMA=1.0 FT3_NSUB=$2 FT3_POUT_Y=0.9
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mpirun -use-hwthread-cpus -n 1 ./ns-example-3d input/gota.load output/g.save VTKS/g.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/b2-3d-N$1-s$2.log 2>&1
echo "N=$1 nsub=$2 exit=$? em $(date +%H:%M)" >> /tmp/b2-3d-status
