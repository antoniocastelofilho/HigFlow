#!/bin/bash
# $1 = velocidade da tampa   $2 = sigma   $3 = numsteps   $4 = rotulo
cd /home/castelo/HiGFlow-System/HigFlow && set -a && . ./varsrc && set +a
cd higflow/example3d_FrontTracking
cp input/gota.load.par.contr.yaml /tmp/pc-$4.bak
sed -i "s/^  numsteps: .*/  numsteps: $3/" input/gota.load.par.contr.yaml
export FT3_R=${FT3_R:-0.25} FT3_SIGMA=$2 FT3_NSUB=3 FT3_POUT_Y=0.9 FT3_ADVECTA=1 FT3_TAMPA=$1
export FT3_RHO0=1.0 FT3_RHO1=1.0 FT3_MU0=1.0 FT3_MU1=1.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mpirun -use-hwthread-cpus -n 1 ./ns-example-3d input/gota.load output/g.save VTKS/g.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/adv-$4.log 2>&1
echo "$4 exit=$? em $(date +%H:%M)" >> /tmp/adv-status
cp /tmp/pc-$4.bak input/gota.load.par.contr.yaml
