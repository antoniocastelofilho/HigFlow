#!/bin/bash
# Avisa quando o FT FECHAR, ou quando algo exigir decisao.  Quatro saidas:
#   0 terminou   1 lancador morreu/encerrou   2 ORFAO   3 travou
#
# Vive no REPOSITORIO e nao no scratchpad: a versao anterior morava em /tmp e
# desapareceu no reboot de 09/10 junto com todos os registros.  Vigia que nao
# sobrevive ao que ele deveria vigiar nao e' vigia.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$AQUI/registro" 2>/dev/null
D=$AQUI
B=/sys/class/power_supply/BAT0
ST=$AQUI/registro/status.log
L=$AQUI/registro/corrida.log
TOTAL=2001

passo()  { awk '/step:/{print $2}' "$D/output/o.save.init.yaml" 2>/dev/null || echo 0; }
solver() { ps -eo comm | grep -qE 'ns-example|prterun'; }
carga()  { printf "%s%% (razao %s%%) %s" "$(cat $B/capacity)" \
           "$(( $(cat $B/charge_now)*100/$(cat $B/charge_full) ))" "$(cat $B/status)"; }

ant=$(passo); parado=0
while :; do
    P=$(passo)

    if [ "${P:-0}" -ge "$TOTAL" ]; then
        # espera o solver sair: o passo chega a 2001 DENTRO do laco, e os
        # oraculos de fim do FT -- salto de Laplace, corrente espuria --
        # imprimem depois dele.  Com teto, para solver travado nao engolir o aviso.
        esp=0
        while solver && [ $esp -lt 60 ]; do sleep 30; esp=$((esp+1)); done
        [ $esp -ge 60 ] && echo "(o solver nao saiu em 30 min; relato o que ha')"
        echo "*** O FT N=50 (com filtro) TERMINOU: passo $P de $TOTAL ***"
        echo "  quadros de superficie: $(ls "$D"/VTKS/*frente_* 2>/dev/null | wc -l)"
        echo "  quadros de campo:      $(ls "$D"/VTKS/p.print_0-* 2>/dev/null | wc -l)"
        echo "  carga: $(carga)"
        LC_ALL=C grep -a 'BOLHA 3D\|PAREDE\|Laplace\|espuria' "$L" 2>/dev/null | tail -5 | sed 's/^/  /'
        n=$(LC_ALL=C grep -ac '===> TSUR passo' "$L" 2>/dev/null); n=${n:-0}
        m=$(LC_ALL=C grep -a '===> TSUR passo' "$L" 2>/dev/null | LC_ALL=C grep -avc 'moveu 0, TSUR 0'); m=${m:-0}
        echo "  aplicacoes do filtro: $n, das quais com movimento: $m"
        LC_ALL=C grep -a 'RETOMANDO\|comecando\|trecho fechou\|CONCLUIDO' "$ST" 2>/dev/null | tail -6 | sed 's/^/  /'
        exit 0
    fi

    if LC_ALL=C grep -aq 'tres trechos sem avancar\|lancador encerrado\|sobrou solver vivo' "$ST" 2>/dev/null; then
        echo "*** O LANCADOR ENCERROU SEM COMPLETAR: passo $P de $TOTAL ***"
        tail -10 "$ST" 2>/dev/null; exit 1
    fi

    if ! pgrep -f "$D/lanca\.sh$" >/dev/null; then
        echo "*** O LANCADOR MORREU: passo $P de $TOTAL ***"
        tail -8 "$ST" 2>/dev/null
        LC_ALL=C grep -aiE 'abort|error|segmenta' "$L" 2>/dev/null | tail -4
        exit 1
    fi

    # ORFAO: solver vivo com a ultima linha do lancador dizendo que espera.
    # E' a assinatura do defeito de 05/10, que custou 24 h de bateria.
    if solver && tail -1 "$ST" 2>/dev/null | LC_ALL=C grep -aq 'recarregando'; then
        echo "*** ORFAO: solver rodando com o lancador esperando recarga ***"
        echo "  passo $P, carga $(carga)"
        ps -eo pid,ppid,etimes,comm | awk '$4 ~ /ns-example|prterun/ {print "  "$0}'
        exit 2
    fi

    # TRAVOU: solver vivo e o passo nao anda ha' 4 tiques (1 h).  Um quadro de
    # campo custa ~20 min, nao uma hora, logo 1 h parado nao e' impressao.
    if solver && [ "$P" = "$ant" ]; then
        parado=$((parado+1))
        if [ $parado -ge 4 ]; then
            echo "*** NAO AVANCA ha' ~1 h com solver vivo: passo $P ***"
            LC_ALL=C grep -a '===> Step' "$L" 2>/dev/null | tail -2 | sed 's/^/  /'
            exit 3
        fi
    else
        parado=0
    fi
    ant=$P
    sleep 900
done
