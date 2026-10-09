#!/bin/bash
# Roda UM TRECHO de UMA cadencia.  Retoma de ponto valido se houver, pelo mesmo
# protocolo da producao: `.ok` presente, sem `.gravando`, e passo igual ao do
# init.yaml; se invalido, recua para a copia de seguranca; sem copia, recomeca.
set -u
CADA=${1:?uso: roda-um.sh <cadencia> <nucleo>}
NUC=${2:?uso: roda-um.sh <cadencia> <nucleo>}
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
D=$AQUI/cada-$CADA
M=/home/castelo/HiGFlow-System/HigFlow
TOTAL=2001
ST=$AQUI/registro/status.log
LOG=$AQUI/registro/cada-$CADA.log
cd "$M" && set -a && . ./varsrc && set +a
cd "$D"
diz() { echo "$(date '+%m-%d %H:%M') [cada=$CADA] $*" >> $ST; }

PASSO=0
if [ -f output/o.save.init.yaml ] && [ -f output/o.save.ok ]; then
    PI=$(awk '/step:/{print $2}' output/o.save.init.yaml)
    PO=$(awk '/^passo/{print $2}' output/o.save.ok)
    [ -f output/o.save.gravando ] && { diz "gravacao INTERROMPIDA"; PI=-1; }
    if [ "$PI" = "$PO" ] && [ "${PI:-0}" -gt 0 ]; then
        PASSO=$PI
    elif [ -f seguranca/o.save.ok ]; then
        PS=$(awk '/^passo/{print $2}' seguranca/o.save.ok)
        diz "ponto invalido (init=$PI ok=$PO); recuando para $PS"
        cp -f seguranca/o.save.* output/ 2>/dev/null; rm -f output/o.save.gravando
        PASSO=$PS
    else
        diz "ponto invalido e sem copia -- RECOMECANDO"; rm -f output/o.save.gravando; PASSO=0
    fi
fi

export FT3_RHO0=1000.0 FT3_RHO1=100.0 FT3_MU0=10.0 FT3_MU1=1.0 FT3_SIGMA=24.5
export FT3_R=0.25 FT3_CX=0.5 FT3_CY=0.5 FT3_CZ=0.5 FT3_NSUB=4 FT3_ADVECTA=1
export FT3_SALVA=1 FT3_CAMPO_CADA=0        # sem campo: a varredura mede a SUPERFICIE
export FT3_TSUR_CADA=$CADA FT3_TSUR_LIM=0.05 FT3_TSUR_RLX=0.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [ "$PASSO" -gt 0 ]; then
    REST=$(( TOTAL - PASSO ))
    [ "$REST" -le 0 ] && { diz "ja' completou $PASSO"; exit 0; }
    sed -i "s/^  numsteps: .*/  numsteps: $REST/" output/o.save.par.contr.yaml
    diz "RETOMANDO do passo $PASSO, faltam $REST"
    CARGA=output/o.save
else
    diz "comecando do zero, $TOTAL passos (filtro $([ "$CADA" = 0 ] && echo DESLIGADO || echo "a cada $CADA"))"
    CARGA=input/gota.load
fi

taskset -c $NUC stdbuf -oL -eL /usr/bin/time -f "PAREDE %e s" \
  mpirun --bind-to none -use-hwthread-cpus -n 1 \
  ./ns-example-3d "$CARGA" output/o.save VTKS/p.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 >> $LOG 2>&1
R=$?
diz "trecho fechou exit=$R, passo salvo=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null)"
exit $R
