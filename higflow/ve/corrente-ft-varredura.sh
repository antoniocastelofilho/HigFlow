#!/bin/bash
# Liga o fim do FT ao inicio da varredura de cadencia.
#
# Espera o FT fechar os 2001 passos, remonta as quatro corridas do zero (o
# `prepara.sh` e' idempotente e apaga o que houver), e passa a vez ao
# orquestrador, que espera a carga pela conta dele.
#
# Se o FT NAO fechar -- lancador morto, trechos sem avanco --, a varredura NAO
# comeca: ela compete pelos mesmos nucleos, e comecar com o FT incompleto
# trocaria o resultado principal por um resultado auxiliar.
set -u
VE=/home/castelo/HiGFlow-System/HigFlow/higflow/ve
FT=$VE/ft-N50-salva
VA=$VE/varredura-cadencia
TOTAL=2001
REG=/tmp/corrente-ft-varredura.log
diz() { echo "$(date '+%m-%d %H:%M') $*" | tee -a $REG; }

diz "=== corrente: espera o FT fechar os $TOTAL passos ==="
while :; do
    P=$(awk '/step:/{print $2}' $FT/output/o.save.init.yaml 2>/dev/null || echo 0)
    [ "${P:-0}" -ge "$TOTAL" ] && { diz "o FT fechou no passo $P"; break; }
    if ! pgrep -f "$FT/lanca\.sh$" >/dev/null; then
        diz "o lancador do FT morreu no passo $P de $TOTAL -- NAO inicio a varredura"
        diz "    (ela compete pelos mesmos nucleos; o FT e' o resultado principal)"
        exit 1
    fi
    sleep 900
done

diz "remontando as quatro corridas da varredura"
if $VA/prepara.sh >> $REG 2>&1; then diz "prepara.sh: ok"
else diz "*** prepara.sh FALHOU -- nao inicio ***"; exit 2; fi

diz "passando a vez ao orquestrador da varredura"
cd $VA
exec ./lanca.sh
