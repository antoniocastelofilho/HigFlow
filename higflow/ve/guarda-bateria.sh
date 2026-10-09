#!/bin/bash
# GUARDA DE BATERIA INDEPENDENTE.  Rede de seguranca, nao o mecanismo normal.
#
# O mecanismo normal e' o piso do lancador (capacity <= 15%), que para o trecho
# num instante escolhido e espera recarga.  Este guarda existe para o caso em
# que esse mecanismo NAO age -- o lancador morreu, travou, ou o kill dele falhou
# e ele desistiu deixando o solver vivo.  Foi assim que a carga chegou a 4% em
# 05/10: o lancador achou que tinha parado o trecho e o solver seguiu 24 h.
#
# Por isso os limiares ficam ABAIXO do piso do lancador: se o lancador esta'
# fazendo o trabalho dele, este guarda nunca dispara.
#
#   <= GUARDA_AVISO (12%)  so' registra, para a divergencia aparecer
#   <= GUARDA_MATA  (9%)   mata a arvore do solver e o lancador
#   <= GUARDA_JA    (5%)   mata SEM esperar gravacao, que e' emergencia
#
# Limiares por ambiente para poderem ser TESTADOS -- um guarda nao exercitado
# da' falsa confianca, que e' pior que nao ter guarda.
set -u
B=${GUARDA_BAT:-/sys/class/power_supply/BAT0}
AVISO=${GUARDA_AVISO:-12}
MATA=${GUARDA_MATA:-9}
JA=${GUARDA_JA:-5}
ESPERA=${GUARDA_ESPERA:-60}
D=${GUARDA_DIR:-/home/castelo/HiGFlow-System/HigFlow/higflow/ve}
REG=${GUARDA_REG:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/registro-guarda.log}
PADRAO=${GUARDA_PADRAO:-'ns-example|prterun|mpirun|orted'}

diz() { echo "$(date '+%m-%d %H:%M:%S') $*" | tee -a $REG; }
cap() { cat $B/capacity; }
raz() { echo $(( $(cat $B/charge_now) * 100 / $(cat $B/charge_full) )); }
alvos() { ps -eo pid,comm | awk -v p="$PADRAO" '$2 ~ p {print $1}'; }

descendentes() {
    local p=$1 f
    for f in $(ps -eo pid=,ppid= | awk -v P="$p" '$2==P {print $1}'); do
        descendentes "$f"; echo "$f"
    done
}

# Mata os alvos e os descendentes deles.  NAO sobe a arvore: a primeira versao
# subia enquanto o pai fosse `bash`, e isso ia do roda.sh para o lanca.sh e dai'
# para o shell da sessao -- um guarda que mata o supervisor e quem o chamou.
# Nao precisa subir: `time`/`taskset` terminam quando o filho morre, e roda.sh
# segue para o fim sozinho.  O lanca.sh e' morto a parte, pelo caminho exato do
# roteiro, para nao lancar outro trecho.
#
# A lista e' colhida ANTES de matar -- matar o pai reparenta os filhos para o
# init e a caminhada os perde, que foi o defeito original do lancador.
mata_tudo() {
    local lista="" a p resto=""
    for a in $(alvos); do
        lista="$lista $(descendentes "$a") $a"
    done
    lista=$(echo $lista | tr ' ' '\n' | awk 'NF' | awk '!v[$0]++')
    [ -z "$lista" ] && { diz "nenhum alvo vivo"; return 0; }
    diz "matando: $(echo $lista | tr '\n' ' ')"
    for p in $lista; do kill -INT  "$p" 2>/dev/null; done; sleep 10
    for p in $lista; do kill -TERM "$p" 2>/dev/null; done; sleep 8
    for p in $lista; do kill -KILL "$p" 2>/dev/null; done; sleep 3
    for p in $lista; do kill -0 "$p" 2>/dev/null && resto="$resto $p"; done
    if [ -n "$resto" ]; then
        diz "ATENCAO: sobraram$resto -- segunda rodada de KILL"
        for p in $resto; do kill -KILL "$p" 2>/dev/null; done; sleep 3
        resto=""
        for p in $lista; do kill -0 "$p" 2>/dev/null && resto="$resto $p"; done
    fi
    [ -n "$resto" ] && { diz "FALHEI: ainda vivos$resto"; return 1; }
    diz "alvos encerrados e verificados por PID"
    return 0
}

diz "=== guarda de bateria: aviso<=${AVISO}% mata<=${MATA}% emergencia<=${JA}% ==="
avisou=0
while :; do
    C=$(cap); R=$(raz)
    if [ "$C" -le "$JA" ]; then
        diz "EMERGENCIA: ${C}% (razao ${R}%) -- matando SEM esperar gravacao"
        for p in $(pgrep -f "$D/.*/lanca\.sh$"); do kill -KILL $p 2>/dev/null; done
        mata_tudo
        diz "*** maquina salva do desligamento; retomar do ultimo salvamento valido ***"
        exit 2
    fi
    if [ "$C" -le "$MATA" ]; then
        diz "PISO DE SEGURANCA: ${C}% (razao ${R}%) -- o lancador NAO parou o solver"
        # espera a gravacao em curso fechar, para nao rasgar o conjunto
        for d in $D/*/output; do
            if [ -f "$d/o.save.gravando" ]; then
                diz "gravacao em curso em $d -- esperando ate' ${ESPERA}s"
                for i in $(seq 1 $ESPERA); do [ -f "$d/o.save.gravando" ] || break; sleep 1; done
            fi
        done
        for p in $(pgrep -f "$D/.*/lanca\.sh$"); do kill $p 2>/dev/null; done
        sleep 1
        mata_tudo
        diz "*** parei a corrida para a maquina nao cair; carga ${C}% ***"
        exit 1
    fi
    if [ "$C" -le "$AVISO" ] && [ "$avisou" -eq 0 ]; then
        diz "aviso: ${C}% (razao ${R}%) -- abaixo do piso do lancador, que deveria ter parado"
        avisou=1
    fi
    [ "$C" -gt "$AVISO" ] && avisou=0
    sleep 60
done
