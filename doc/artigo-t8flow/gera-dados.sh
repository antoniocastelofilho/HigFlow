#!/bin/bash
# Gera os .dat do artigo RODANDO o codigo, nunca digitando numeros.
#
# MESMA CONVENCAO DOS OUTROS ARTIGOS DESTA ARVORE: os numeros vivem em dados/*.dat
# e o texto os le' por pgfplotstable.  Regenerar uma medida troca o .dat e o texto
# acompanha, sem a chance de a tabela e a prosa discordarem.
#
# Uso:  cd doc/artigo-t8flow && ./gera-dados.sh [diretorio-do-t8flow]
set -u
T8=${1:-../../t8flow}
EX=$T8/exemplos/poisson-mms
D=$(dirname "$0")/dados

if [ ! -x "$EX" ]; then
  echo "ERRO: $EX nao existe.  Compilar com 'make 2d' primeiro." >&2
  exit 1
fi

# extrai (h, erro_inf, erro_l2, ordem_inf) da tabela do exemplo
corre() {   # corre <arquivo> <args...>
  local out=$1; shift
  echo "# gerado por gera-dados.sh -- NAO EDITAR" > "$out"
  echo "nivel celulas h einf eL2 ordinf ordL2" >> "$out"
  mpirun -np 1 "$EX" "$@" 2>/dev/null \
    | awk '/^ *[0-9]+ +[0-9]+ +[0-9.e+-]+ /{
        oi = ($5 == "exato" || $5 == "-") ? "---" : $5;
        ol = ($7 == "exato" || $7 == "-") ? "---" : $7;
        print $1, $2, $3, $4, $6, oi, ol }' >> "$out"
  printf '  %-34s %s\n' "$(basename "$out")" "$(($(wc -l < "$out") - 2)) niveis"
}

NIV=${NIV:-3,4,5}
echo "== malha uniforme, esquema de fluxo de segunda ordem"
for s in grau1 grau2 grau3 grau4 trig; do
  corre "$D/unif-fluxo2-$s.dat" --solucao=$s --niveis=$NIV --tol=1e-14
done
echo "== malha uniforme, compacto de quarta ordem"
for s in grau3 grau4 trig; do
  corre "$D/unif-lele4-$s.dat" --esquema=lele4 --solucao=$s --niveis=$NIV --tol=1e-14
done
echo "== malha uniforme, MLS + Gauss (ordem 3 e 4)"
for s in grau3 grau4 trig; do
  corre "$D/unif-mls3-$s.dat" --esquema=mlsgauss --mls-ordem=3 --solucao=$s --niveis=$NIV --tol=1e-14
  corre "$D/unif-mls4-$s.dat" --esquema=mlsgauss --mls-ordem=4 --mls-pontos=60 --solucao=$s --niveis=$NIV --tol=1e-14
done
echo "== o contraste da quadratura da fonte"
corre "$D/unif-mls4-pm-grau4.dat" --esquema=mlsgauss --mls-ordem=4 --mls-pontos=60 \
      --fonte-ponto-medio --solucao=grau4 --niveis=$NIV --tol=1e-14
echo "== malha adaptada (quadtree com refino aleatorio)"
for s in grau1 grau2 grau3 grau4 trig; do
  corre "$D/adap-fluxo2-$s.dat" --aleatorio=0.4 --semente=7 --solucao=$s --niveis=$NIV --tol=1e-14
  corre "$D/adap-mls3-$s.dat" --esquema=mlsgauss --mls-ordem=3 --aleatorio=0.4 --semente=7 \
        --solucao=$s --niveis=$NIV --tol=1e-14
done
echo "== o fechamento de dois pontos, como contraste"
corre "$D/unif-fluxo2-2pt-grau2.dat" --fechamento=2 --solucao=grau2 --niveis=$NIV --tol=1e-14
echo "pronto."
