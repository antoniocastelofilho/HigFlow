#!/bin/bash
# VOF viscoelastico, lambda=0.4 beta=0.2, N=50 -- a resolucao de referencia.
#
# RETOMAVEL.  Este roteiro roda UM TRECHO: se o diretorio ja' tem um ponto de
# salvamento valido, retoma dele; senao comeca do zero.  Chamado de novo, segue
# de onde parou.  Quem decide quando chamar e' o `lanca.sh`.
#
# O que faz um ponto ser VALIDO: o `.ok` existe e o passo dele e' o mesmo do
# init.yaml.  O `.ok` e' gravado por rename depois de todos os `._dp*`, logo
# divergencia entre os dois significa gravacao interrompida no meio -- conjunto
# metade novo e metade velho, que nao serve.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$AQUI/registro" 2>/dev/null
M=/home/castelo/HiGFlow-System/HigFlow
TOTAL=2001                       # passos da corrida inteira
cd "$M" && set -a && . ./varsrc && set +a
cd "$AQUI"

LOG=$AQUI/registro/corrida.log
# A impressao digital do binario no registro.  O diretorio da corrida tem uma
# COPIA do executavel, e ela ja' ficou atras de duas correcoes do ponto de
# salvamento -- a corrida teria rodado sem elas e o registro nao diria nada.
# Registrar faz a divergencia aparecer no mesmo lugar onde se olha o progresso.
echo "$(date +%H:%M) binario $(md5sum ns-example-3d | cut -c1-12) (fonte: $(md5sum "$M/higflow/example3d_VOF_rising/ns-example-3d" 2>/dev/null | cut -c1-12))" >> $AQUI/registro/status.log
PASSO=0
if [ -f output/o.save.init.yaml ] && [ -f output/o.save.ok ]; then
    PI=$(awk '/step:/{print $2}' output/o.save.init.yaml)
    PO=$(awk '/^passo/{print $2}' output/o.save.ok)
    # um `.gravando` sobrando e' informacao certa: aquela gravacao nao
    # terminou, e o conjunto corrente nao serve mesmo que os dois numeros
    # batam (o `.ok` pode ser o da geracao anterior).
    if [ -f output/o.save.gravando ]; then
        echo "$(date +%H:%M) gravacao INTERROMPIDA (sinal .gravando presente)" >> $AQUI/registro/status.log
        PI=-1
    fi
    if [ "$PI" = "$PO" ] && [ "${PI:-0}" -gt 0 ]; then
        PASSO=$PI
    else
        echo "$(date +%H:%M) ponto de salvamento INCONSISTENTE (init=$PI ok=$PO)" >> $AQUI/registro/status.log
        # a copia de seguranca, se houver, e' o conjunto anterior COMPLETO
        if [ -f seguranca/o.save.ok ]; then
            PS=$(awk '/^passo/{print $2}' seguranca/o.save.ok)
            echo "$(date +%H:%M) voltando a' copia de seguranca, passo $PS" >> $AQUI/registro/status.log
            cp -f seguranca/o.save.* output/ 2>/dev/null
            rm -f output/o.save.gravando
            PASSO=$PS
        else
            echo "$(date +%H:%M) sem copia de seguranca -- RECOMECANDO do zero" >> $AQUI/registro/status.log
            rm -f output/o.save.gravando
            PASSO=0
        fi
    fi
fi

export SB_RHO0=1000.0 SB_RHO1=100.0 SB_MU0=10.0 SB_MU1=1.0
export SB_R=0.25 SB_CX=0.5 SB_CY=0.5 SB_CZ=0.5
export SB_SALVA=1
# CADENCIA DO CAMPO.  Em N=50 um quadro de campo sai com ~300 MB e leva ~20 min
# -- mais que um passo.  A animacao da bolha NAO precisa dele: a interface ja'
# sai com p, vel, S e T interpolados nos vertices.  O campo serve as' linhas de
# corrente, e para isso tres amostras bastam: 80 quadros de dtp=0.025 e' t=2.0.
export SB_CAMPO_CADA=80
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [ "$PASSO" -gt 0 ]; then
    REST=$(( TOTAL - PASSO ))
    if [ "$REST" -le 0 ]; then
        echo "$(date +%H:%M) ja' completou $PASSO de $TOTAL -- nada a fazer" >> $AQUI/registro/status.log
        exit 0
    fi
    sed -i "s/^  numsteps: .*/  numsteps: $REST/" output/o.save.par.contr.yaml
    echo "$(date +%H:%M) RETOMANDO do passo $PASSO, faltam $REST" >> $AQUI/registro/status.log
    CARGA=output/o.save
else
    echo "$(date +%H:%M) comecando do zero, $TOTAL passos" >> $AQUI/registro/status.log
    CARGA=input/sobe.load
fi

# nucleo fisico 1 (o 0 fica para o sistema); os pares SMT deste i9 sao (k,k+8),
# logo 1 e 9 sao o MESMO nucleo -- usar numeros contiguos nao da' nucleos
# distintos, e e' por isso que ha' so' um numero aqui.
taskset -c 1 stdbuf -oL -eL /usr/bin/time -f "PAREDE %e s" mpirun --bind-to none -use-hwthread-cpus -n 1 \
  ./ns-example-3d "$CARGA" output/o.save VTKS/p.print \
  -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 \
  >> $LOG 2>&1
R=$?
echo "$(date +%H:%M) trecho terminou exit=$R, passo salvo=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null)" >> $AQUI/registro/status.log
exit $R
