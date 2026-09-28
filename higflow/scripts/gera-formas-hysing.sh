#!/bin/bash
# Gera os .dat de FORMA da bolha do Hysing, nos dois metodos e nos MESMOS
# instantes, para dados/formas/.  Re-executavel: quando as corridas terminarem,
# rodar de novo com outros instantes e o relatorio acompanha.
#
#   $@ = instantes desejados (ex: 0.0 0.5 1.0 2.0 3.0)
#
# FT   : despejo lagrangeano, DATA/frentes/frente_<passo>.dat  (passo = t/dt)
# VOF  : contorno FracVol=0,5 extraido do VTK do quadro (quadro = t/dtp)
set -e
H=/home/castelo/HiGFlow-System/HigFlow/higflow
D=/home/castelo/HiGFlow-System/HigFlow/doc/relatorio-front-tracking/dados/formas
S=/tmp/claude-1001/-home-castelo-Freeflow/7d349f3e-226e-48b6-8f28-a1c8ed0e93f3/scratchpad
DT=0.00025; DTP=0.01
k=0
for t in "$@"; do
    passo=$(python3 -c "print(round($t/$DT/40)*40)")
    quadro=$(python3 -c "print(round($t/$DTP))")
    ft=$H/example2d_FrontTracking/DATA/frentes/frente_$(printf %07d $passo).dat
    vtk=$H/example2d_VOF/VTKS/h.print_0-$quadro.vtk
    if [ -f "$ft" ]; then
        cp "$ft" $D/hy_ft_$k.dat
        echo "t=$t  FT  passo $passo  -> hy_ft_$k.dat ($(( $(wc -l < $D/hy_ft_$k.dat) - 1 )) pontos)"
    else
        echo "t=$t  FT  passo $passo  AUSENTE ($ft)"
    fi
    if [ -f "$vtk" ]; then
        python3 $S/contorno-vof.py "$vtk" $D/hy_vof_$k.dat
        mv_msg=$(printf "t=%s  VOF quadro %s -> hy_vof_%s.dat" $t $quadro $k); echo "$mv_msg"
    else
        echo "t=$t  VOF quadro $quadro AUSENTE ($vtk)"
    fi
    k=$((k+1))
done
