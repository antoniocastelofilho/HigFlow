#!/bin/bash
# Monta as quatro corridas da varredura de cadencia do filtro.
#
# A PERGUNTA: aumentar a frequencia do filtro protege mais a simulacao, ou apaga
# a cuspide junto com a ondulacao?  O indicador do TSUR nao distingue os dois --
# zigue-zague e curvatura real inclinam a aresta igual --, logo so' uma medida da
# FEATURE decide.
#
# OS PARAMETROS NAO SAO PALPITE, sairam de medida nos dados que ja' existem:
#
#  * 2001 PASSOS, nao truncado.  A assimetria (z_c-z_min)/(z_max-z_c) em t=1,2
#    vale 0,955, ABAIXO de 1: a cuspide ainda nao existe.  Em t=4 o caso da
#    cuspide da' 1,213 contra 0,88-0,91 nos outros.  Truncar em t=1,5 nao mediria
#    cuspide nenhuma.
#  * N=25, porque em N=50 cada corrida custaria 58 h e sao quatro.  O efeito a
#    medir -- ondulacao e cuspide -- existe em N=25: a corrida N=25 sem filtro
#    chegou a 21% das arestas acima do limiar e R/Req colapsado a 0,047.
#  * CADENCIAS 0, 20, 10, 5.  O 0 e' o controle sem filtro, sem o qual nenhuma
#    das outras tres significa algo.  O 20 e' a producao.
#  * lambda=0,4 beta=0,2, que e' o caso onde a tese poe a cuspide.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M=/home/castelo/HiGFlow-System/HigFlow
MOLDE=$M/higflow/ve/ft-l04-b02          # a corrida N=25 da matriz, como molde
cd "$AQUI"
for cada in 0 20 10 5; do
    d=cada-$cada
    rm -rf $d && mkdir -p $d/{input,output,VTKS,DATA,amrs}
    cp $MOLDE/gera-malha.py $d/
    cp $MOLDE/input/* $d/input/
    cp $M/higflow/example3d_FrontTracking/ns-example-3d $d/
    ( cd $d && python3 gera-malha.py 25 1.0 1.0 3.0 >/dev/null )
    # salvamento a cada QUADRO, as mesmas precaucoes da producao
    sed -i 's/^  dts: .*/  dts: 0.025/;s/^  numsteps: .*/  numsteps: 2001/' $d/input/*.par.contr.yaml
    B=$(ls $d/input/*.par.contr.yaml | sed 's#.*/##;s/\.par\.contr\.yaml//')
    ln -sfn ../input/$B.domain.yaml $d/output/o.save.domain.yaml
    ln -sfn ../input/$B.bc.yaml     $d/output/o.save.bc.yaml
    echo "  cada-$cada: malha $(find $d/amrs -type f | wc -l) arquivos, base $B, $(grep -c . $d/input/*.par.contr.yaml | head -1 | cut -d: -f2) linhas de par"
done
echo "  --- os caminhos dos yamls resolvem? ---"
for cada in 0 20 10 5; do
    n=0; mal=0
    for p in $(LC_ALL=C grep -aho 'path: *[^ ]*' cada-$cada/input/*.yaml | sed 's/path: *//'); do
        n=$((n+1)); [ -e "cada-$cada/$p" ] || mal=$((mal+1))
    done
    printf "  cada-%-3s %d caminhos, %d faltando\n" "$cada" "$n" "$mal"
done
