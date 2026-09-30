#!/bin/bash
# $1 = velocidade da parede (taxa de cisalhamento, pois H=1)   $2 = numsteps
W=/home/castelo/HiGFlow-System/worktrees/vof-3d
cd $W && set -a && . ./varsrc && set +a
export PETSC_DIR=/home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4/x86_64
cd $W/higflow/example3d_VOF_VE
sed -i "s/^  numsteps: .*/  numsteps: ${2:-4001}/" input/cis.load.par.contr.yaml
export VE_CISALHA=$1 VE_RHO0=1.0 VE_RHO1=1.0 VE_MU0=1.0 VE_MU1=1.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mpirun -use-hwthread-cpus -n 1 ./ns-example-3d input/cis.load output/c.save VTKS/c.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/ve-cis.log 2>&1
echo "cisalhamento U=$1 exit=$? em $(date +%H:%M)" >> /tmp/cis-status
