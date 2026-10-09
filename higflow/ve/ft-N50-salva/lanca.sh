#!/bin/bash
# Espera a carga passar de 95%, roda o VOF em trechos retomaveis, e mantem uma
# copia de seguranca VERIFICADA do ponto de salvamento.
#
# POR QUE EM TRECHOS.  Medido nesta maquina: com um solver rodando, a bateria
# DRENA mesmo ligada na tomada -- a margem de potencia do adaptador nao cobre a
# conta.  Logo a corrida nao cabe necessariamente numa carga, e o laco aqui e'
# o que a faz atravessar varias: roda ate' a carga chegar ao piso, para num
# instante ESCOLHIDO (logo depois de um salvamento, nao no meio de uma
# gravacao), e espera recarregar.
#
# POR QUE A COPIA DE SEGURANCA.  O `.ok` DETECTA gravacao interrompida mas nao
# recupera dela.  A copia e' feita de fora, e so' e' aceita se o `.ok` nao
# mudou entre o inicio e o fim da copia -- se mudou, o solver estava gravando e
# a copia e' descartada e repetida.  Com isso ha' sempre um conjunto completo:
# o corrente, ou o anterior.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$AQUI/registro" 2>/dev/null
cd "$AQUI"
B=/sys/class/power_supply/BAT0
TOTAL=2001
ST=$AQUI/registro/status.log

# OS DOIS NUMEROS DA BATERIA DISCORDAM, e escolher o errado quebra o lancador
# em silencio.  Medido agora: `capacity` = 21 enquanto
# charge_now/charge_full = 25,95%.  A referencia implicita de `capacity` e'
# ~8,6 Ah, nem charge_full (6,955) nem charge_full_design (8,790) -- o driver
# revisou charge_full para baixo depois da descarga funda e `capacity` seguiu
# outra referencia.
#
# A instrucao foi "acima de 95%", e o 95% que o usuario le' e' `capacity` --
# ele observou 97%-98% nesta maquina, valor que charge_now/charge_full ja' nao
# alcancaria se eu tomasse charge_full como teto.  Logo a PORTA e' `capacity`.
#
# Mas gate unico em `capacity` tem o risco oposto: se a referencia dele estiver
# velha, ele pode nao chegar a 95 nem com a bateria cheia, e o lancador
# esperaria para sempre sem nada rodando.  Por isso ha' a segunda porta --
# bateria efetivamente cheia pela razao propria -- e os dois numeros entram no
# registro a cada tique, para a divergencia ficar visivel em vez de decidir
# sozinha.
ALVOPCT=95
PISOPCT=15
cheia() {
    local c=$(cat $B/capacity)
    local s=$(cat $B/status)
    local r=$(( $(cat $B/charge_now) * 100 / $(cat $B/charge_full) ))
    [ "$c" -ge "$ALVOPCT" ] && return 0       # o numero que o usuario le'
    [ "$r" -ge 99 ] && return 0               # cheia pela razao propria
    # TERCEIRA PORTA, porque as duas primeiras podem passar RASPANDO.  Medido:
    # a referencia de `capacity` e' 8,36 Ah e charge_full esta' em 6,955 -- se
    # charge_full nao subir, `capacity` para em ~83% e a 1a porta nunca dispara;
    # se subir para ~8,3 enquanto a carga termina em ~7,9, a razao da' 95% e a
    # 2a tambem nao.  Nao ha' limite de carga configurado nesta maquina (sem
    # charge_control_end_threshold), logo o carregador vai ate' o fim, e o sinal
    # inequivoco disso e' o `status`.
    #
    # Mas `status` SOZINHO nao serve: ele ja' foi pego dizendo "Full" com a
    # bateria descarregando, e acreditar nele ali lancaria a corrida com 20% de
    # carga.  Por isso o par -- `Full` E carga de verdade.
    [ "$s" = "Full" ] && [ "$r" -ge 90 ] && return 0
    return 1
}
no_piso() { [ "$(cat $B/capacity)" -le "$PISOPCT" ]; }
pct() { printf "%s%% (razao %s%%)" "$(cat $B/capacity)" \
        "$(( $(cat $B/charge_now) * 100 / $(cat $B/charge_full) ))"; }
diz() { echo "$(date '+%m-%d %H:%M') $*" | tee -a $ST; }

# --------------------------------------------------- parar a ARVORE
# `pkill -P $SOLVER` mata so' os FILHOS DIRETOS de roda.sh.  A arvore real e'
#   lanca.sh -> roda.sh -> taskset -> stdbuf -> time -> prterun -> ns-example-3d
# logo aquilo matava o `taskset` e deixava prterun e ns-example-3d ORFAOS a 99%
# de CPU.  Foi o que aconteceu: o lancador passou a esperar recarga achando que
# o trecho tinha acabado, enquanto o solver seguia drenando -- a carga desceu a
# 4% em vez de parar no piso.  A corrida nao se perdeu (os salvamentos
# continuaram), mas o lancador poderia ter lancado um SEGUNDO trecho sobre o
# primeiro, dois solvers gravando no mesmo output.
#
# `setsid` + `kill -$grupo` tambem NAO resolve, e isto foi MEDIDO: no teste o
# grupo nao continha a arvore (setsid forkou, o PGID nao ficou igual ao $!) e
# quem matou foi a varredura.  Logo a varredura e' o mecanismo, nao o reserva.
#
# A lista de descendentes e' colhida ANTES de matar: se o pai morre primeiro, os
# filhos sao reparentados para o init e a caminhada os perde.  E a morte vai do
# mais fundo para o mais raso, pela mesma razao.
descendentes() {        # todos os descendentes de $1, do mais fundo para o raso
    local p=$1 f
    for f in $(ps -eo pid=,ppid= | awk -v P="$p" '$2==P {print $1}'); do
        descendentes "$f"
        echo "$f"
    done
}
para_arvore() {
    local raiz=$1 lista g p resto=""
    lista=$(descendentes "$raiz")           # ANTES de matar
    g=$(ps -o pgid= -p "$raiz" 2>/dev/null | tr -d ' ')
    diz "encerrando raiz $raiz (grupo ${g:-?}, so' informativo); descendentes: $(echo $lista | tr '\n' ' ')"
    for p in $lista;       do kill -INT  "$p" 2>/dev/null; done; sleep 15
    for p in $lista $raiz; do kill -TERM "$p" 2>/dev/null; done; sleep 10
    for p in $lista $raiz; do kill -KILL "$p" 2>/dev/null; done
    # NAO matar o grupo.  Eu havia deixado `kill -KILL -$g` como reserva depois
    # de MEDIR que o grupo nao continha a arvore, e conclui' que era inutil.
    # Era pior que inutil: sem controle de tarefa os filhos ficam no grupo do
    # SHELL, de modo que $g era o PID do PROPRIO LANCADOR -- e ele se matou ao
    # parar o trecho no piso de 15%.  O trecho parou certo, os tres descendentes
    # foram achados e mortos, e ninguem sobrou para retomar.
    #
    # A licao e' a mesma das outras vezes: reserva que eu nao exercitei nao e'
    # reserva, e' superficie de falha.  O kill por PID, que E' exercitado, da'
    # conta.  O $g fica so' no registro.
    sleep 3
    # VERIFICA por PID -- so' isso distingue "matei" de "achei que matei"
    for p in $lista $raiz; do kill -0 "$p" 2>/dev/null && resto="$resto $p"; done
    if [ -n "$resto" ]; then
        diz "ATENCAO: AINDA VIVOS da nossa arvore:$resto -- NAO lanco outro trecho"
        return 1
    fi
    # e avisa de solver ALHEIO, sem mata-lo: pode ser a outra corrida
    local outros
    outros=$(ps -eo pid,comm | awk '$2 ~ /ns-example|prterun|mpirun|orted/ {print $1}')
    [ -n "$outros" ] && diz "nota: ha' solver de OUTRA corrida vivo ($(echo $outros|tr '\n' ' ')) -- nao mexi"
    diz "arvore encerrada e verificada por PID"
    return 0
}

# --------------------------------------------------- copia de seguranca
guarda() {
    while :; do
        sleep 240
        [ -f output/o.save.ok ] || continue
        # RECUSA enquanto a gravacao esta' em curso.  Sem isto, comparar so' o
        # `.ok` antes e depois nao basta: os `._dp*` sao escritos ANTES de o
        # `.ok` ser renomeado, logo ha' uma janela em que os dados ja' sao da
        # geracao nova e a marca ainda diz a anterior -- copiar ali guardaria um
        # conjunto rasgado rotulado como bom, destruindo a copia de que se
        # precisaria.
        [ -f output/o.save.gravando ] && continue
        A=$(cat output/o.save.ok 2>/dev/null)
        PI=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null)
        PO=$(echo "$A" | awk '/^passo/{print $2}')
        [ -n "${PI:-}" ] && [ "$PI" = "${PO:-x}" ] || continue      # inconsistente: nao copia
        [ -f seguranca/o.save.ok ] && \
          [ "$(awk '/^passo/{print $2}' seguranca/o.save.ok)" = "$PO" ] && continue  # ja' e' esta
        mkdir -p seguranca.tmp && rm -f seguranca.tmp/*
        cp -f output/o.save.* seguranca.tmp/ 2>/dev/null
        # ACEITA so' se o solver nao gravou durante a copia
        if [ "$(cat output/o.save.ok 2>/dev/null)" = "$A" ] && \
           [ ! -f output/o.save.gravando ]; then
            rm -rf seguranca && mv seguranca.tmp seguranca
            diz "copia de seguranca do passo $PO ($(du -sh seguranca 2>/dev/null|cut -f1))"
        else
            rm -rf seguranca.tmp     # o solver gravou no meio: descarta e repete
        fi
    done
}

# --------------------------------------------------- espera a carga
diz "=== lancador do FT N=50 (com filtro): espera carga acima de ${ALVOPCT}% ==="
while ! cheia; do
    diz "aguardando: $(pct) de ${ALVOPCT}%"
    sleep 600
done
diz "carga em $(pct) -- LANCANDO"

guarda & GUARDA=$!
trap 'kill $GUARDA 2>/dev/null' EXIT

falhas=0
while :; do
    PASSO=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null || echo 0)
    if [ "${PASSO:-0}" -ge "$TOTAL" ]; then diz "*** CONCLUIDO: $PASSO de $TOTAL passos ***"; break; fi

    ./roda.sh & SOLVER=$!
    # vigia a carga enquanto o trecho roda
    while kill -0 $SOLVER 2>/dev/null; do
        sleep 120
        if no_piso; then
            diz "carga no piso ($(pct)) -- parando o trecho num instante escolhido"
            para_arvore $SOLVER || { diz "*** parando o lancador: sobrou solver vivo ***"; exit 1; }
            break
        fi
    done
    wait $SOLVER 2>/dev/null; R=$?

    NOVO=$(awk '/step:/{print $2}' output/o.save.init.yaml 2>/dev/null || echo 0)
    # CONSERTO DE UMA JANELA ESTREITA.  `save_all` grava sempre step+1, logo o
    # salvamento de partida de um trecho RETOMADO grava P+1 enquanto os campos
    # em disco ainda sao o estado de P.  Mesmo instante t, contador um a'
    # frente: cair entre a partida e o primeiro salvamento do laco faria a
    # retomada seguinte pular UM dt.
    #
    # Com dts = dtp = 0,025 o primeiro salvamento do laco so' pode vir no 13o
    # passo, logo NOVO == PASSO+1 significa INEQUIVOCAMENTE que o trecho morreu
    # nessa janela -- nenhum trecho que avancou de verdade para' em +1.  Ai' o
    # contador volta para onde os campos estao.
    #
    # So' para PASSO > 0: na primeira partida, step 1 com campos de t=0 custa um
    # dt, enquanto voltar para 0 custaria os 13 min de partida inteiros.
    if [ "${PASSO:-0}" -gt 0 ] && [ "${NOVO:-0}" -eq "$(( PASSO + 1 ))" ]; then
        diz "trecho morreu na partida -- contador de volta de $NOVO para $PASSO"
        sed -i "s/step: $NOVO/step: $PASSO/" output/o.save.init.yaml
        NOVO=$PASSO
    fi
    diz "trecho fechou exit=$R, passo $PASSO -> $NOVO, carga $(pct)"
    if [ "${NOVO:-0}" -ge "$TOTAL" ]; then diz "*** CONCLUIDO: $NOVO de $TOTAL ***"; break; fi
    if [ "${NOVO:-0}" -le "${PASSO:-0}" ]; then
        falhas=$((falhas+1))
        diz "o trecho NAO avancou (falha $falhas de 3)"
        [ $falhas -ge 3 ] && { diz "*** tres trechos sem avancar -- parando para nao girar em falso ***"; break; }
    else
        falhas=0
    fi
    # espera recarregar antes do proximo trecho
    while ! cheia; do
        diz "recarregando: $(pct) de ${ALVOPCT}%"
        sleep 900
    done
done
diz "lancador encerrado"
