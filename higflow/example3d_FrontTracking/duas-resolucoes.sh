#!/bin/bash
# EXPERIMENTO DE DUAS RESOLUCOES: separa vies de interpolacao de erro de
# discretizacao.  So' UMA coisa varia -- a malha -- e a resolucao da SUPERFICIE
# acompanha (ds/h = 0,75 nos dois pares), senao seriam duas variaveis.
#
#   N=20 h=0,050 nsub=3 ds=0,0377   ds/h = 0,75
#   N=40 h=0,025 nsub=4 ds=0,0188   ds/h = 0,75
#
# Fisica fixa: R=0,25  sigma=0,2  tampa=4,0  dt=5e-4  ate' t=0,25 (500 passos).
# O corte em 500 e' anterior ao passo 660, onde o guarda da particao da unidade
# disparou no caso longo -- a gota ainda nao alcancou a tampa.
#
# SEQUENCIAL de proposito: as duas corridas compartilham amrs/, input/ e output/,
# e roda-las juntas faria uma sobrescrever a malha da outra.
cd /home/castelo/HiGFlow-System/HigFlow && set -a && . ./varsrc && set +a
cd higflow/example3d_FrontTracking
cp input/gota.load.par.contr.yaml /tmp/pc-2res.bak
sed -i "s/^  numsteps: .*/  numsteps: 501/" input/gota.load.par.contr.yaml
export FT3_R=0.25 FT3_SIGMA=0.2 FT3_TAMPA=4.0 FT3_ADVECTA=1 FT3_POUT_Y=0.9
export FT3_RHO0=1.0 FT3_RHO1=1.0 FT3_MU0=1.0 FT3_MU1=1.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

for par in "20 3" "40 4"; do
    set -- $par
    N=$1; NSUB=$2
    python3 gera-malha.py $N 1.0 >/dev/null
    export FT3_NSUB=$NSUB
    echo "== N=$N nsub=$NSUB inicio $(date +%H:%M)" >> /tmp/2res-status
    mpirun -use-hwthread-cpus -n 1 ./ns-example-3d input/gota.load output/g.save VTKS/g.print \
      -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 > /tmp/2res-N$N.log 2>&1
    echo "== N=$N exit=$? em $(date +%H:%M)" >> /tmp/2res-status
done
cp /tmp/pc-2res.bak input/gota.load.par.contr.yaml
python3 gera-malha.py 20 1.0 >/dev/null
echo "== EXPERIMENTO COMPLETO em $(date +%H:%M)" >> /tmp/2res-status
