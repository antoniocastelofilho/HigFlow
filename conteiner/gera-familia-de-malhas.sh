#!/bin/bash
# As malhas do teste de escalonamento, geradas e NAO versionadas.
#
# Sao dado DERIVADO: 4 318 arquivos de ~150 bytes que saem exatos deste
# roteiro em segundos.  Versionar isso encheria o repositorio de arquivo
# minusculo para guardar o que uma linha reproduz.
#
# ---------------------------------------------------------------- o desenho
#
# FRACO: np = k^3 e m = 5k da' 20 250 celulas por rank EXATAMENTE, em todos os
# pontos.  Nao e' coincidencia -- celulas = 162*m^3, entao cubo com cubo
# cancela.  A alternativa, np em potencias de 2, obrigaria m nao inteiro e a
# carga por rank variaria ate' 12% entre pontos; a eficiencia fraca e' uma
# razao de tempos, e 12% de diferenca de carga vira 12% de erro direto nela.
# O ponto np=8 e' m=10: a malha que ja' tem linha de base medida.
#
# FORTE: uma malha so', m=32 (5,3 M celulas, 20 736 por rank em 256 -- de
# proposito igual a' carga do fraco, para os dois testes falarem da mesma
# quantidade de trabalho por processo).
#
# ------------------------------------------------- por que --para-np, e o custo
#
# A higtree reparte ARVORES INTEIRAS, entao o maior bloco e' um piso para o
# rank mais carregado.  Na malha original o maior bloco tem 14,8% das celulas
# e o ganho nao passa de 6,75x por mais nucleos que se de'.
#
# MAS LEVANTAR O TETO NAO ACELERA POR SI.  Medido em np=8, mesma malha, tres
# repeticoes: 33 blocos 153,1 s, 66 blocos 153,1 s, 264 blocos 157,4 s -- o teto
# subiu 26 vezes e o tempo ficou onde estava.  Era previsivel e eu nao previ: o
# teto era 6,75x e a medida ja' era 2,64x, entao o desequilibrio de particao
# nunca foi o termo que mandava.
#
# ENTAO POR QUE CORTAR.  Porque 6,75x e' teto ABSOLUTO: em 256 nucleos ele ata,
# e nenhuma curva ate' la' existe sem isto.  E' pre-condicao do regime de muitos
# processos, nao ganho hoje.
#
# E POR QUE --para-np E NAO O MAXIMO.  Cortar mais do que o alvo pede so' pode
# custar (mais franja, mais arvores por rank) sem nada a ganhar, ja' que o teto
# deixa de atar assim que cabe na fatia.  A folga 2 e' margem de empacotamento,
# nao perfeccionismo.
#
# A consequencia fica para quem for ler a curva forte: a malha e' fina o
# bastante para 256, e o tempo de referencia em np=1 pode sair inflado por isso.
# O ganho contra a propria malha e' a definicao padrao e e' o que se publica,
# mas o tempo serial da malha grossa tem de aparecer ao lado.
set -e
EX=${EX:-$PWD/higflow/example3d_complex}
DST=${DST:-$EX/malhas}
GERA=${GERA:-$PWD/ci/tools/gera-malha-complexa.py}

for k in 1 2 3 4 5 6; do
    np=$((k*k*k)); m=$((5*k))
    python3 "$GERA" --origem "$EX" --destino "$DST/fraco-np$np" \
        --m $m --para-np $np --raiz-relativa "malhas/fraco-np$np" | head -3
done

python3 "$GERA" --origem "$EX" --destino "$DST/forte-m32" \
    --m 32 --para-np 256 --raiz-relativa "malhas/forte-m32" | head -3

echo
echo "Para rodar uma delas, o prefixo e' argv[1] do executavel:"
echo "   ./ns-complex-3d malhas/fraco-np27/input/example-3d.load output/s VTKS/v ..."
