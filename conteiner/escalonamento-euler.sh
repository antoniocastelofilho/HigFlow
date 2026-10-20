#!/usr/bin/env bash
#==============================================================================
# escalonamento-euler.sh -- o protocolo de medida de escalonamento.
#
# Uso, DENTRO da imagem (ou fora dela, numa arvore construida):
#   conteiner/escalonamento-euler.sh forte  --np 1,2,4,8,16,32,64
#   conteiner/escalonamento-euler.sh fraco  --np 1,8,64 --base 4
#
# PRINCIPIO: tempo so' vale depois de o resultado valer.  O portao 0 compara o
# paralelo com o serial ANTES de cronometrar qualquer coisa, e num caso NAO
# exato -- no caso exato tudo e' zero de maquina e nada discrimina.  Foi assim
# que tres defeitos de franja sobreviveram a uma suite verde.
#
# TRES DECISOES DE PROTOCOLO, e cada uma muda o que o numero significa:
#
#  1. AS ETAPAS SEPARADAS, e nao o tempo total.  O `poisson3d-mms` ja' imprime
#     franja, montagem e resolucao com `--tempos`; o total mistura arranque com
#     calculo, e cada um escala de um jeito.  (A primeira versao deste script
#     subtraia dois comprimentos de corrida para separar os dois -- truque
#     necessario no caso do HiGFlow, desnecessario aqui, e eu o escrevi antes
#     de ler as opcoes do exemplo.)
#  2. REPETICOES, com mediana e dispersao.  Uma corrida por ponto nao distingue
#     ganho de ruido de maquina compartilhada.  O padrao e' 3.
#  3. SEQUENCIAL.  Duas corridas disputando nucleo dao um tempo que nao e' de
#     nenhuma das duas.
#
# E AS CORRIDAS NAO GRAVAM VTK: escrita por rank e' E/S que nao escala, e num
# sistema de arquivos paralelo ela pode dominar.  Ver `[E4]`.
#==============================================================================
set -u

MODO="${1:-}"; shift || true
NPS="1,2,4,8"; REPS=3; BASE=3; CURTO=2; LONGO=12
ALVO="poisson"            # poisson | complexo
SAIDA="escalonamento"

while [ $# -gt 0 ]; do
  case "$1" in
    --np)    NPS="$2"; shift 2;;
    --reps)  REPS="$2"; shift 2;;
    --base)  BASE="$2"; shift 2;;
    --alvo)  ALVO="$2"; shift 2;;
    --saida) SAIDA="$2"; shift 2;;
    *) echo "opcao desconhecida: $1" >&2; exit 2;;
  esac
done

case "$MODO" in
  forte|fraco) ;;
  *) echo "uso: $0 {forte|fraco} [--np L] [--reps N] [--base N] [--alvo poisson|complexo]" >&2
     exit 2;;
esac

LANCA="${LANCA:-mpirun -np}"     # no Euler: `srun --mpi=pmix -n` -- ver [E2]
mkdir -p "$SAIDA"

#--------------------------------------------------------------- leitura
# A linha de resultado do `poisson3d-mms` e':
#   nivel celulas h erro_inf iters franja montagem resolucao total ...
erro_de() {   # $1 = np, $2 = nivel -> o erro, como texto
  $LANCA "$1" ./exemplos/poisson3d/poisson3d-mms --niveis="$2" \
      --solucao=trig --tempos 2>/dev/null \
    | awk -v n="$2" '$1 == n { print $4; exit }'
}

tempos_de() {  # $1 = np, $2 = nivel -> "montagem resolucao total"
  $LANCA "$1" ./exemplos/poisson3d/poisson3d-mms --niveis="$2" \
      --solucao=trig --tempos 2>/dev/null \
    | awk -v n="$2" '$1 == n { print $7, $8, $9; exit }'
}

#------------------------------------------------------------------ portao 0
# O PARALELO TEM DE DAR O MESMO QUE O SERIAL, em caso NAO exato.
portao_zero() {
  echo "== portao 0: paralelo identico ao serial (caso NAO exato) =="
  local maxnp
  maxnp=$(echo "$NPS" | tr ',' '\n' | sort -n | tail -1)
  if [ "$maxnp" -lt 2 ]; then
    echo "   PULADO: so' ha' np=1 na lista -- nada a comparar."
    echo "   ATENCAO: isto NAO e' 'passou'.  Sem np>=2 o portao nao mede nada."
    return 0
  fi
  # O caso trigonometrico NAO e' reproduzido exatamente pelo esquema, entao o
  # erro e' nao nulo e discrimina.  Num caso exato os dois dariam zero.
  local e1 e2
  e1=$(erro_de 1 4)
  e2=$(erro_de 2 4)
  echo "   np=1: $e1"
  echo "   np=2: $e2"
  if [ -z "$e1" ] || [ -z "$e2" ]; then
    echo "   FALHOU: nao consegui ler o erro -- NAO prosseguir para o tempo."
    echo "   (um portao que nao consegue olhar nao e' um portao que aprovou)"
    return 1
  fi
  [ "$e1" = "$e2" ] && { echo "   OK"; return 0; }
  echo "   FALHOU: os dois diferem.  Tempo medido agora seria tempo de"
  echo "   resposta errada.  Parar aqui."
  return 1
}

#------------------------------------------------------------------ varredura
portao_zero || exit 1

#------------------------------------------------------- portao 1: tamanho
# UM CASO PEQUENO DEMAIS PRODUZ UMA CURVA QUE NAO FALA DO METODO.
#
# Medido no nivel 3 (512 celulas): resolucao em 0,0013 s em serie e 0,0026 s em
# np=2 -- "eficiencia 25%", que e' so' a comunicacao contra trabalho nenhum.
# Uma tabela dessas tem a FORMA de um resultado de escalonamento e nao e' um.
#
# O limiar e' sobre o TEMPO DE RESOLUCAO EM SERIE, e nao sobre a contagem de
# celulas: o que importa e' haver trabalho que baste para a divisao vencer a
# comunicacao, e isso depende da maquina tanto quanto da malha.
LIMIAR_S="${LIMIAR_S:-1.0}"
portao_um() {
  echo "== portao 1: o caso e' grande o bastante? =="
  local t
  t=$(tempos_de 1 "$BASE" | awk '{print $2}')
  if [ -z "$t" ]; then
    echo "   FALHOU: nao li o tempo de resolucao em serie."
    return 1
  fi
  echo "   resolucao em serie, nivel $BASE: $t s   (limiar $LIMIAR_S s)"
  if python3 -c "import sys; sys.exit(0 if $t >= $LIMIAR_S else 1)"; then
    echo "   OK"
    return 0
  fi
  echo "   RECUSADO: com $t s de resolucao, a curva mede comunicacao contra"
  echo "   trabalho nenhum.  Suba --base ate' passar do limiar, ou ajuste"
  echo "   LIMIAR_S de propo'sito -- e, se ajustar, diga isso ao ler a tabela."
  return 1
}
portao_um || exit 1

ARQ="$SAIDA/$MODO-$ALVO.dat"
echo "# escalonamento $MODO, alvo $ALVO, $REPS repeticoes por ponto" > "$ARQ"
echo "# np nivel resol_mediana resol_min resol_max montagem ganho eficiencia" >> "$ARQ"

PRIMEIRO=""
for NP in $(echo "$NPS" | tr ',' ' '); do
  # ESCALONAMENTO FRACO: o nivel cresce com `np` para o trabalho POR PROCESSO
  # ficar constante.  Em 3D cada nivel multiplica as celulas por 8, logo `np`
  # tem de octuplicar junto -- e' por isso que a lista util e' 1,8,64.
  if [ "$MODO" = "fraco" ]; then
    NIVEL=$(python3 -c "
import math
print($BASE + round(math.log($NP, 8)))")
  else
    NIVEL=$BASE
  fi

  RESOL=(); MONT=()
  for R in $(seq 1 "$REPS"); do
    LIDO=$(tempos_de "$NP" "$NIVEL")
    [ -z "$LIDO" ] && { echo "np=$NP nivel=$NIVEL: NAO LI os tempos -- parando"; exit 1; }
    MONT+=("$(echo "$LIDO" | awk '{print $1}')")
    RESOL+=("$(echo "$LIDO" | awk '{print $2}')")
  done

  LINHA=$(python3 - "$NP" "$NIVEL" "${RESOL[*]}" "${MONT[*]}" "$PRIMEIRO" <<'PY'
import sys, statistics
np_, niv = int(sys.argv[1]), int(sys.argv[2])
resol = [float(x) for x in sys.argv[3].split()]
mont = [float(x) for x in sys.argv[4].split()]
prim = sys.argv[5]
med = statistics.median(resol)
base = float(prim) if prim else med
# ESCALONAMENTO FRACO: o ideal e' tempo CONSTANTE, nao tempo caindo -- por isso
# a coluna de eficiencia divide por `np` so' no forte.  Aqui ela sai como
# `base/med`, que no fraco e' a eficiencia direta.
print("%d %d %.4f %.4f %.4f %.4f %.3f %.1f" %
      (np_, niv, med, min(resol), max(resol), statistics.median(mont),
       base / med, 100.0 * base / med / np_))
PY
)
  [ -z "$PRIMEIRO" ] && PRIMEIRO=$(echo "$LINHA" | awk '{print $3}')
  echo "$LINHA" | tee -a "$ARQ"
done

echo
echo "gravado em $ARQ"
echo
echo "LEITURA: a eficiencia do ULTIMO ponto costuma ser pessimista quando o"
echo "numero de processos iguala os nucleos do no' -- nao sobra folga para o"
echo "sistema, e parte da perda e' contencao de maquina, nao do metodo."
