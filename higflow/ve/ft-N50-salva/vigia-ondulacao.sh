#!/bin/bash
# Acompanha a ONDULACAO RETIDA na corrida filtrada, contra o GEMEO SEM FILTRO
# medido na MESMA janela de quadros.
#
# A primeira versao citava um numero fixo -- "sem filtro na mesma fase: 4,10%" --
# colhido na janela t=1,00-1,25.  Esse numero ENVELHECEU: na janela t=1,275-1,40
# o gemeo ja' estava em 7,49%.  Comparar o presente com uma referencia velha
# subestima o beneficio do filtro e, pior, nao detectaria o filtro perdendo
# terreno se a ondulacao de fundo crescesse junto.
#
# Agora o gemeo e' medido nos MESMOS quadros, e o criterio e' PAREADO enquanto
# ele existir.  Ele vai ate' o quadro 64 (t=1,6), onde a corrida sem filtro foi
# parada; dali em diante o criterio passa a ser absoluto, e o registro diz qual
# esta' em vigor -- um criterio que troca sem avisar e' pior que um criterio pior.
#
# O gemeo e' ESTATICO (corrida encerrada), logo cada quadro dele e' medido uma
# vez e guardado em cache: remedir a cada tique tomaria CPU do solver.
set -u
D=/home/castelo/HiGFlow-System/HigFlow/higflow/ve/ft-N50-salva
mkdir -p "$D/registro" 2>/dev/null
G=/home/castelo/HiGFlow-System/HigFlow/higflow/ve/ft-l04-b02-N50/VTKS-parcial-sem-filtro
M=/home/castelo/HiGFlow-System/HigFlow/higflow/examples-common
REG=${ONDUL_REG:-$D/registro/ondulacao.log}
CACHE=${ONDUL_CACHE:-$D/registro/gemeo.cache}
JAN=${ONDUL_JAN:-6}
RAZAO_MIN=${ONDUL_RAZAO:-3.0}     # criterio PAREADO: beneficio minimo aceitavel
ABS_MAX=${ONDUL_ABS:-4.0}         # criterio ABSOLUTO, quando o gemeo acaba
NUC=${ONDUL_NUC:-3}
touch $CACHE

mede() {        # $1 = arquivo vtk -> imprime a fracao acima do limiar
    taskset -c $NUC python3 $M/mede-ondulacao.py "$1" --limiar 0.05 2>/dev/null \
      | sed -n 's/.*acima do limiar [0-9.]*: *[0-9]* de [0-9]* *(\([0-9.]*\)%).*/\1/p'
}
mede_gemeo() {  # com cache, porque o gemeo nao muda
    local q=$1 v
    v=$(awk -F= -v q="$q" '$1==q {print $2; exit}' $CACHE)
    if [ -z "$v" ]; then
        local f=$(ls $G/*frente_*-$q.vtk 2>/dev/null | head -1)
        [ -n "$f" ] || { echo ""; return; }
        v=$(mede "$f")
        [ -n "$v" ] && echo "$q=$v" >> $CACHE
    fi
    echo "$v"
}

while :; do
    if ! pgrep -f "$D/lanca\.sh$" >/dev/null; then
        echo "$(date '+%m-%d %H:%M') o lancador saiu; encerrando" >> $REG; exit 0
    fi
    mapfile -t qs < <(ls $D/VTKS/*frente_*.vtk 2>/dev/null \
        | sed 's/.*-\([0-9]*\)\.vtk/\1/' | sort -n | tail -$JAN)
    if [ ${#qs[@]} -lt $JAN ]; then sleep 1800; continue; fi

    sc=0; nc=0; pc=0; sg=0; ng=0
    for q in "${qs[@]}"; do
        f=$(ls $D/VTKS/*frente_*-$q.vtk 2>/dev/null | head -1)
        [ -n "$f" ] || continue
        v=$(mede "$f"); [ -n "$v" ] || continue
        sc=$(python3 -c "print($sc+$v)"); nc=$((nc+1))
        pc=$(python3 -c "print(max($pc,$v))")
        g=$(mede_gemeo "$q")
        if [ -n "$g" ]; then sg=$(python3 -c "print($sg+$g)"); ng=$((ng+1)); fi
    done
    [ $nc -eq 0 ] && { sleep 1800; continue; }
    mc=$(python3 -c "print('%.2f'%($sc/$nc))")
    passo=$(LC_ALL=C grep -a '===> Step' $D/registro/corrida.log | tail -1 | grep -oE '[0-9]+' | head -1)
    t=$(python3 -c "print('%.2f'%(${qs[-1]}*0.025))")

    if [ $ng -eq $nc ]; then
        # PAREADO: o gemeo cobre a janela inteira
        mg=$(python3 -c "print('%.2f'%($sg/$ng))")
        raz=$(python3 -c "print('%.1f'%(($sg/$ng)/max(1e-9,$sc/$nc)))")
        echo "$(date '+%m-%d %H:%M')  passo $passo  t~$t  quadros ${qs[0]}-${qs[-1]}  COM ${mc}% (pior ${pc}%)  SEM ${mg}%  razao ${raz}x  [criterio PAREADO, minimo ${RAZAO_MIN}x]" >> $REG
        if python3 -c "import sys; sys.exit(0 if $raz < $RAZAO_MIN else 1)"; then
            echo "  *** O BENEFICIO CAIU A ${raz}x (minimo ${RAZAO_MIN}x) -- o filtro perdeu terreno ***" >> $REG
            echo "  *** com ${mc}% contra ${mg}% do gemeo na mesma janela; revisar a cadencia" >> $REG
            exit 3
        fi
    else
        # ABSOLUTO: o gemeo acabou (ele vai ate' o quadro 64, t=1,6)
        echo "$(date '+%m-%d %H:%M')  passo $passo  t~$t  quadros ${qs[0]}-${qs[-1]}  COM ${mc}% (pior ${pc}%)  SEM: fora do alcance do gemeo ($ng de $nc quadros)  [criterio ABSOLUTO, maximo ${ABS_MAX}%]" >> $REG
        if python3 -c "import sys; sys.exit(0 if $mc >= $ABS_MAX else 1)"; then
            echo "  *** A ONDULACAO RETIDA CHEGOU A ${mc}% (maximo ${ABS_MAX}%) ***" >> $REG
            echo "  *** sem gemeo para parear daqui em diante; revisar a cadencia" >> $REG
            exit 3
        fi
    fi
    sleep 1800
done
