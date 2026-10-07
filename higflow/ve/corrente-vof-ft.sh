#!/bin/bash
# Liga o fim do VOF ao inicio do FT com filtro, com a janela de recarga entre os
# dois.  Tres etapas, e a do meio e' a que importa:
#
#   1. espera o VOF FECHAR os 2001 passos
#   2. verifica a FIACAO DO FILTRO antes de comprometer dias de maquina --
#      um filtro que silenciosamente nao aplica invalidaria a razao da corrida
#      existir, e isso so' apareceria no fim.  A janela de recarga e' de graca
#      para isso: a maquina esta' ociosa de todo jeito.
#   3. passa a vez ao lancador do FT, que espera a carga pela conta dele
#
# Nao duplica a espera de bateria: o lanca.sh do FT ja' tem as tres portas.
set -u
VE=/home/castelo/HiGFlow-System/HigFlow/higflow/ve
VOF=$VE/vof-N50-salva
FT=$VE/ft-N50-salva
TOTAL=2001
REG=/tmp/corrente-vof-ft.log
diz() { echo "$(date '+%m-%d %H:%M') $*" | tee -a $REG; }
B=/sys/class/power_supply/BAT0
carga() { printf "%s%% (razao %s%%) %s" "$(cat $B/capacity)" \
  "$(( $(cat $B/charge_now)*100/$(cat $B/charge_full) ))" "$(cat $B/status)"; }

diz "=== corrente: espera o VOF fechar os $TOTAL passos ==="
while :; do
    P=$(awk '/step:/{print $2}' $VOF/output/o.save.init.yaml 2>/dev/null || echo 0)
    [ "${P:-0}" -ge "$TOTAL" ] && { diz "o VOF fechou no passo $P"; break; }
    # se o lancador do VOF morreu E o passo nao chegou, nao e' hora de seguir
    if ! pgrep -f "$VOF/lanca\.sh$" >/dev/null; then
        diz "o lancador do VOF morreu no passo $P de $TOTAL -- NAO inicio o FT"
        diz "    (o FT competiria com um VOF incompleto; decida antes)"
        exit 1
    fi
    sleep 600
done

diz "o VOF terminou; carga $(carga).  Verificando a fiacao do filtro antes do FT"
if $FT/fumaca-filtro.sh > /tmp/fumaca-filtro-corrente.log 2>&1; then
    diz "fiacao do filtro: PASSOU"
    LC_ALL=C grep -a '\[ok\]' /tmp/fumaca-filtro-corrente.log | sed 's/^/    /' | tee -a $REG
else
    diz "*** fiacao do filtro FALHOU -- NAO inicio o FT ***"
    LC_ALL=C grep -a '\[FALHA\]' /tmp/fumaca-filtro-corrente.log | sed 's/^/    /' | tee -a $REG
    exit 2
fi

diz "passando a vez ao lancador do FT (ele espera a carga pela conta dele)"
cd $FT
exec ./lanca.sh
