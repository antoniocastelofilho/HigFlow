#!/bin/bash
# A ANALISE da varredura.  Tres grandezas, e as tres juntas respondem a pergunta;
# nenhuma delas sozinha responde.
#
#  BENEFICIO   ondulacao retida -- quanto do zigue-zague o filtro tira.
#              Mais frequencia tem de reduzir.  Se nao reduzir, nao ha' motivo.
#  CUSTO       nitidez da cuspide -- assimetria (z_c-z_min)/(z_max-z_c) e raio de
#              curvatura na traseira.  Mais frequencia NAO pode reduzir a cuspide;
#              se reduzir, o filtro esta' apagando a feature junto com o defeito.
#  CONTROLE    cada-0, sem filtro.  Sem ele as outras tres nao tem referencia, e
#              a direcao de "melhor" seria escolhida e nao medida.
#
# O criterio de decisao esta' escrito ao fim, para nao ser inventado depois de
# ver o resultado.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$AQUI"
M=/home/castelo/HiGFlow-System/HigFlow/higflow/examples-common
CADS="0 20 10 5"
JAN=${JAN:-8}                 # quadros no fim, para descontar a fase do dente de serra

echo "=== estado das quatro ==="
for c in $CADS; do
    printf "  cada-%-3s passo %-5s de 2001   %3s quadros\n" "$c" \
      "$(awk '/step:/{print $2}' cada-$c/output/o.save.init.yaml 2>/dev/null)" \
      "$(ls cada-$c/VTKS/*frente_* 2>/dev/null | wc -l)"
done

echo
echo "=== 1. BENEFICIO: ondulacao retida, media nos $JAN ultimos quadros ==="
for c in $CADS; do
    python3 - "$c" "$JAN" "$M" <<'PY'
import sys, glob, re, subprocess, statistics as st
c, jan, M = sys.argv[1], int(sys.argv[2]), sys.argv[3]
g = sorted(glob.glob("cada-%s/VTKS/*frente_*.vtk"%c),
           key=lambda p: int(re.search(r'-(\d+)\.vtk$',p).group(1)))[-jan:]
v=[]
for f in g:
    o=subprocess.run(["python3",M+"/mede-ondulacao.py",f,"--limiar","0.05"],
                     capture_output=True,text=True).stdout
    m=re.search(r'acima do limiar [\d.]+:\s+\d+ de \d+\s+\(([\d.]+)%\)',o)
    if m: v.append(float(m.group(1)))
if v:
    print("  cada-%-3s  media %6.2f%%  min %6.2f%%  max %6.2f%%  (%d quadros, t~%.2f)"
          % (c, st.mean(v), min(v), max(v), len(v),
             int(re.search(r'-(\d+)\.vtk$',g[-1]).group(1))*0.025))
else: print("  cada-%-3s  sem medida" % c)
PY
done

echo
echo "=== 2. CUSTO: a cuspide no ultimo quadro ==="
printf "  %-8s %9s %9s %9s %9s\n" "cadencia" "z_c" "assim." "R" "R/Req"
for c in $CADS; do
    f=$(ls cada-$c/VTKS/*frente_*.vtk 2>/dev/null | sed 's/.*-\([0-9]*\)\.vtk/\1 &/' | sort -n | tail -1 | cut -d' ' -f2-)
    [ -n "$f" ] || { echo "  cada-$c: sem quadro"; continue; }
    python3 $M/mede-cuspide.py "$f" 2>/dev/null | tail -1 | \
      awk -v c="$c" '{printf "  cada-%-3s %9s %9s %9s %9s\n", c, $2, $4, $5, $6}'
done

echo
echo "=== 3. integridade: o filtro nao pode romper nada ==="
for c in $CADS; do
    L=$AQUI/registro/cada-$c.log
    # `grep -c || echo 0` imprime DUAS contagens quando grep acha zero: o proprio
    # 0 do grep e o do echo.  Sai da' um "0\n0" que o printf mostra em duas linhas.
    n=$(LC_ALL=C grep -ac '===> TSUR passo' $L 2>/dev/null); n=${n:-0}
    mal=$(LC_ALL=C grep -a '===> TSUR passo' $L 2>/dev/null | LC_ALL=C grep -avc 'chi=2'); mal=${mal:-0}
    dv=$(LC_ALL=C grep -ao 'volume [-+0-9.e]*' $L 2>/dev/null | awk '{s+=$2} END {printf "%.2e", s+0}')
    da=$(LC_ALL=C grep -ao 'AREA [-+0-9.e]*' $L 2>/dev/null | awk '{s+=$2} END {printf "%.2e", s+0}')
    pior=$(LC_ALL=C grep -ao 'moveu [0-9]*, TSUR [0-9]*' $L 2>/dev/null | \
           awk '{s=$2+$4; if(s>m) m=s} END {print m+0}')
    printf "  cada-%-3s %4s aplicacoes, %s com chi!=2, soma dV/V %9s, soma dA/A %9s, pior %s de 2562\n" \
      "$c" "$n" "$mal" "$dv" "$da" "$pior"
done

cat <<'TXT'

=== O CRITERIO DE DECISAO, escrito antes do resultado ===

  AUMENTAR a cadencia se, ao passar de 20 para 10 (ou 5):
    (a) a ondulacao retida cair de forma clara -- fora do envelope min-max da
        cadencia 20, nao so' na media;  E
    (b) a assimetria da cuspide NAO cair mais que 2% em relacao a cada-20;  E
    (c) o R/Req nao SUBIR (subir = ponta arredondada = cuspide apagada).

  NAO aumentar se (b) ou (c) falhar: o filtro estaria comprando ondulacao com
  cuspide, e a cuspide e' o resultado.  Preferir 20 tambem em caso de empate em
  (a): menos intervencao e' menos risco sobre a feature.

  E se cada-0 (sem filtro) ja' tiver R/Req colapsado (<0,1, como os 0,047
  medidos em t=4 na corrida N=25 da matriz) enquanto as filtradas nao, isso
  por si e' o resultado principal: sem filtro a cuspide NAO E' MENSURAVEL no
  front-tracking, porque a coisa mais pontuda da superficie e' ruido.
TXT
