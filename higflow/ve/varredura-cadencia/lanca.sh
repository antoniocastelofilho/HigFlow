#!/bin/bash
# Orquestra as QUATRO cadencias em paralelo, com as precaucoes da producao.
#
# Quatro em paralelo e nao em fila porque as quatro tem de ver a MESMA maquina:
# em fila, uma mudanca de carga termica ou de contencao entre a primeira e a
# quarta entraria na comparacao como se fosse efeito da cadencia.  Nucleos 1 a 4,
# que sao FISICOS distintos -- os pares SMT deste i9 sao (k,k+8), logo numeros
# contiguos abaixo de 8 dao nucleos distintos.
#
# A bateria para e retoma as QUATRO juntas, pela mesma razao: trechos
# desalinhados dariam a cada cadencia um historico de contencao diferente.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$AQUI"
B=/sys/class/power_supply/BAT0
ST=$AQUI/registro/status.log
TOTAL=2001
ALVOPCT=95; PISOPCT=15
CADS="0 20 10 5"
NUCS="1 2 3 4"
diz() { echo "$(date '+%m-%d %H:%M') $*" | tee -a $ST; }
pct() { printf "%s%% (razao %s%%)" "$(cat $B/capacity)" \
        "$(( $(cat $B/charge_now)*100/$(cat $B/charge_full) ))"; }
cheia() { local c=$(cat $B/capacity) s=$(cat $B/status)
  local r=$(( $(cat $B/charge_now)*100/$(cat $B/charge_full) ))
  [ "$c" -ge "$ALVOPCT" ] && return 0; [ "$r" -ge 99 ] && return 0
  [ "$s" = "Full" ] && [ "$r" -ge 90 ] && return 0; return 1; }
no_piso() { [ "$(cat $B/capacity)" -le "$PISOPCT" ]; }
passo_de() { awk '/step:/{print $2}' cada-$1/output/o.save.init.yaml 2>/dev/null || echo 0; }
faltam() { local f=0 c; for c in $CADS; do [ "$(passo_de $c)" -ge "$TOTAL" ] || f=$((f+1)); done; echo $f; }

descendentes() { local p=$1 f
  for f in $(ps -eo pid=,ppid= | awk -v P="$p" '$2==P {print $1}'); do descendentes "$f"; echo "$f"; done; }
para_todas() {   # recebe a lista de PIDs raiz
    local lista="" r p resto=""
    for r in "$@"; do lista="$lista $(descendentes "$r") $r"; done
    lista=$(echo $lista | tr ' ' '\n' | awk 'NF' | awk '!v[$0]++')
    diz "encerrando as quatro: $(echo $lista | tr '\n' ' ')"
    for p in $lista; do kill -INT  "$p" 2>/dev/null; done; sleep 15
    for p in $lista; do kill -TERM "$p" 2>/dev/null; done; sleep 10
    for p in $lista; do kill -KILL "$p" 2>/dev/null; done; sleep 3
    # NAO mata grupo: o grupo e' o do proprio lancador, e matar o grupo foi o
    # defeito que derrubou o lancador da producao em 06/10.
    for p in $lista; do kill -0 "$p" 2>/dev/null && resto="$resto $p"; done
    [ -n "$resto" ] && { diz "ATENCAO: sobraram$resto"; return 1; }
    diz "as quatro encerradas e verificadas por PID"; return 0
}

# copia de seguranca verificada das quatro, com o mesmo protocolo da producao
guarda() {
    while :; do
        sleep 300
        for c in $CADS; do
            d=cada-$c
            [ -f $d/output/o.save.ok ] || continue
            [ -f $d/output/o.save.gravando ] && continue
            A=$(cat $d/output/o.save.ok 2>/dev/null)
            PI=$(awk '/step:/{print $2}' $d/output/o.save.init.yaml 2>/dev/null)
            PO=$(echo "$A" | awk '/^passo/{print $2}')
            [ -n "${PI:-}" ] && [ "$PI" = "${PO:-x}" ] || continue
            [ -f $d/seguranca/o.save.ok ] && \
              [ "$(awk '/^passo/{print $2}' $d/seguranca/o.save.ok)" = "$PO" ] && continue
            rm -rf $d/seguranca.tmp && mkdir -p $d/seguranca.tmp
            cp -f $d/output/o.save.* $d/seguranca.tmp/ 2>/dev/null
            if [ "$(cat $d/output/o.save.ok 2>/dev/null)" = "$A" ] && \
               [ ! -f $d/output/o.save.gravando ]; then
                rm -rf $d/seguranca && mv $d/seguranca.tmp $d/seguranca
            else rm -rf $d/seguranca.tmp; fi
        done
    done
}

diz "=== varredura de cadencia: 0, 20, 10, 5 em N=25, 2001 passos, nucleos 1-4 ==="
diz "espera carga acima de ${ALVOPCT}%"
while ! cheia; do diz "aguardando: $(pct)"; sleep 900; done

guarda & GUARDA=$!
trap 'kill $GUARDA 2>/dev/null' EXIT
falhas=0
while [ "$(faltam)" -gt 0 ]; do
    ANT=""; PIDS=""
    set -- $NUCS
    for c in $CADS; do
        nuc=$1; shift
        if [ "$(passo_de $c)" -ge "$TOTAL" ]; then diz "[cada=$c] ja' completa"; continue; fi
        ANT="$ANT $c:$(passo_de $c)"
        ./roda-um.sh "$c" "$nuc" & PIDS="$PIDS $!"
    done
    [ -z "$PIDS" ] && break
    diz "quatro trechos no ar; partida de$ANT; carga $(pct)"
    # vigia a carga enquanto rodam
    while :; do
        vivo=0; for p in $PIDS; do kill -0 $p 2>/dev/null && vivo=1; done
        [ $vivo -eq 0 ] && break
        if no_piso; then
            diz "carga no piso ($(pct)) -- parando as quatro num instante escolhido"
            para_todas $PIDS || { diz "*** sobrou solver vivo; encerrando o lancador ***"; exit 1; }
            break
        fi
        sleep 120
    done
    for p in $PIDS; do wait $p 2>/dev/null; done
    DEP=""; for c in $CADS; do DEP="$DEP $c:$(passo_de $c)"; done
    diz "trechos fecharam; de$ANT para$DEP; carga $(pct)"
    # avancou alguma?
    av=0
    for par in $ANT; do c=${par%%:*}; a=${par##*:}; [ "$(passo_de $c)" -gt "${a:-0}" ] && av=1; done
    if [ $av -eq 0 ]; then
        falhas=$((falhas+1)); diz "NENHUMA avancou (falha $falhas de 3)"
        [ $falhas -ge 3 ] && { diz "*** tres rodadas sem avanco -- parando ***"; exit 1; }
    else falhas=0; fi
    [ "$(faltam)" -gt 0 ] && while ! cheia; do diz "recarregando: $(pct)"; sleep 900; done
done
diz "*** VARREDURA COMPLETA: $(for c in $CADS; do printf "cada-%s:%s " $c $(passo_de $c); done)***"
diz "rode ./compara.sh para a analise"
