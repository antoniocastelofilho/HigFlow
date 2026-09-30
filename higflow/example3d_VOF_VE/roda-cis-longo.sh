#!/bin/bash
# CONFIRMACAO: o residuo de 1,83% e' transiente ou discretizacao?
#
# exp(-t/De) com t=2 e De=0,5 da' exatamente 1,83% -- o erro medido em A_xy.
# Aqui dt dobra (5e-4 -> 1e-3) com o MESMO numero de passos, entao t vai a 4 e
# sao OITO tempos de relaxacao: exp(-8) = 0,034%.
#
# O TESTE E' DECISIVO NUM SENTIDO: se o erro CAIR para ~0,03%, e' transiente --
# um passo MAIOR daria MAIS erro de discretizacao, nao menos.  Se ficar em 1,8%,
# a causa e' outra e a hipotese morre.
W=/home/castelo/HiGFlow-System/worktrees/vof-3d
cd $W && set -a && . ./varsrc && set +a
export PETSC_DIR=/home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4/x86_64
cd $W/higflow/example3d_VOF_VE
sed -i "s/^  numsteps: .*/  numsteps: 4001/" input/cis.load.par.contr.yaml
sed -i "s/^  dt: .*/  dt: 0.001/"            input/cis.load.par.contr.yaml
export VE_CISALHA=1.0 VE_RHO0=1.0 VE_RHO1=1.0 VE_MU0=1.0 VE_MU1=1.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mpirun -use-hwthread-cpus -n 1 ./ns-example-3d input/cis.load output/c.save VTKS/c.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/ve-cis-8tau.log 2>&1
echo "cisalhamento 8 tau exit=$? em $(date +%H:%M)" >> /tmp/cis-status
sed -i "s/^  dt: .*/  dt: 0.0005/" input/cis.load.par.contr.yaml
