#!/usr/bin/env bash
#==============================================================================
# verificar-euler.sh -- roteiro de verificacao da imagem HiGFlow.
#
# Uso:
#   conteiner/verificar-euler.sh higflow-euler.sif            # local, sem Slurm
#   conteiner/verificar-euler.sh higflow-euler.sif --euler    # inclui dois nos
#
# PRINCIPIO: "compilou" nao e' verificacao.  Cada portao abaixo compara contra um
# valor de FORMA FECHADA ou contra uma invariante, e falha com o numero a vista.
# Os portoes estao em ordem crescente de custo: o que falha barato falha primeiro.
#
# GUARDA DE VACUIDADE: todo portao conta quantos casos ele de fato exercitou e
# FALHA se for zero.  Um arreio ausente, um binario que nao roda ou um grep que
# nao casa produzem "0 falhas" -- que e' indistinguivel de aprovacao, e foi
# exatamente assim que medidas vazias ja' passaram por boas neste projeto.
#==============================================================================
set -uo pipefail
SIF="${1:?uso: verificar-euler.sh <imagem.sif> [--euler]}"
EULER=0; [ "${2:-}" = "--euler" ] && EULER=1
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
falhas=0; portoes=0

# COMO O ROTEIRO SE EXERCITA.  Com VERIF_RAIZ=<arvore>, `ap` roda na arvore
# local em vez de no conteiner, mapeando /opt/HigFlow para ela.  Serve para
# testar o PROPRIO roteiro antes de haver imagem -- um roteiro de verificacao
# nao exercitado e' a ultima coisa que se deve entregar.  Sem a variavel, roda
# no conteiner, que e' o caso de uso.
if [ -n "${VERIF_RAIZ:-}" ]; then
  ap() { local c="$1"; shift
         case "$c" in /opt/HigFlow/*) c="$VERIF_RAIZ/${c#/opt/HigFlow/}" ;; esac
         "$c" "$@"; }
  echo ">> MODO LOCAL: raiz $VERIF_RAIZ (nao e' verificacao de imagem)"
else
  ap() { apptainer exec "$SIF" "$@"; }
fi
titulo() { echo; echo "=== $* ==="; }
ok()  { portoes=$((portoes+1)); printf "  [ok]    %-42s %s\n" "$1" "${2:-}"; }
nao() { portoes=$((portoes+1)); falhas=$((falhas+1)); printf "  [FALHA] %-42s %s\n" "$1" "${2:-}"; }

#------------------------------------------------------------------ 1. inventario
titulo "1. a imagem e' coerente consigo mesma"
if [ -n "${VERIF_RAIZ:-}" ]; then
  INV=$(PETSC_DIR="$VERIF_RAIZ/bibliotecas/petsc-3.25.4/x86_64" \
        T8CODE="$VERIF_RAIZ/bibliotecas/t8code/install" \
        bash "$VERIF_RAIZ/conteiner/inventario.sh" 2>&1) || true
else
  INV=$(ap /opt/HigFlow/conteiner/inventario.sh 2>&1) || true
fi
echo "$INV" | sed 's/^/     /'
# O inventario marca com *** tudo que esta' errado; nenhum *** e' o portao.
n=$(echo "$INV" | grep -c '\*\*\*' || true)
[ "$n" -eq 0 ] && ok "inventario sem alertas" || nao "inventario com $n alerta(s)" "ver acima"
# E o que NAO pode faltar, conferido por fora para o portao nao depender do proprio script
echo "$INV" | grep -q 'otimizado (sem PETSC_USE_DEBUG)' \
  && ok "PETSc otimizado" "o build de depuracao custa 2-5x" \
  || nao "PETSc em build de DEPURACAO" "reconfigurar com --with-debugging=0"

#--------------------------------------------- 2. oraculos de forma fechada
titulo "2. oraculos de forma fechada (os arreios em arvore)"
ARREIOS="teste-esfera3d teste-fracao3d teste-curvatura teste-tensao teste-indicadora
         teste-b1-3d teste-tsur teste-amr3d teste-persistencia"
rodados=0
for t in $ARREIOS; do
  b=/opt/HigFlow/higflow/tests/front-tracking/$t
  bb=$b; [ -n "${VERIF_RAIZ:-}" ] && bb="$VERIF_RAIZ/${b#/opt/HigFlow/}"
  if [ ! -x "$bb" ] && ! apptainer exec "$SIF" test -x "$b" 2>/dev/null; then
    nao "$t" "AUSENTE na imagem -- portao vazio, nao aprovacao"
    continue
  fi
  saida=$(ap "$b" 2>&1); rc=$?
  rodados=$((rodados+1))
  ult=$(echo "$saida" | grep -aE 'PASSOU|FALHOU|FALHA' | tail -1)
  if [ $rc -eq 0 ] && echo "$saida" | grep -aq 'PASSOU'; then
    ok "$t" "${ult:0:72}"
  else
    nao "$t" "exit=$rc ${ult:0:60}"
    echo "$saida" | tail -4 | sed 's/^/            /'
  fi
done
# GUARDA DE VACUIDADE do portao 2
[ "$rodados" -ge 7 ] && ok "cobertura dos arreios" "$rodados de 9 rodaram" \
                      || nao "cobertura dos arreios" "so' $rodados de 9 rodaram -- verificacao vazia"

#------------------------------------------- 3. um numero conhecido, ponta a ponta
titulo "3. um numero conhecido: o salto de Laplace (B2)"
# O teste-tensao ja' confere conservacao e equilibrio contra forma fechada; aqui
# se extrai o NUMERO e compara com o teto, para o portao nao depender do veredito
# impresso pelo proprio arreio.
S=$(ap /opt/HigFlow/higflow/tests/front-tracking/teste-tensao 2>&1 | grep -a 'PORTAO tensao' | tail -1)
if [ -n "$S" ]; then
  cons=$(echo "$S" | sed -n 's/.*conserv=\([0-9.e+-]*\).*/\1/p')
  equi=$(echo "$S" | sed -n 's/.*equilib=\([0-9.e+-]*\).*/\1/p')
  python3 - "$cons" "$equi" <<'PY' && ok "salto de Laplace" "conserv=$cons equilib=$equi" \
                                    || nao "salto de Laplace" "conserv=$cons equilib=$equi"
import sys
c, e = float(sys.argv[1]), float(sys.argv[2])
sys.exit(0 if (c < 1e-12 and e < 1e-3) else 1)
PY
else
  nao "salto de Laplace" "nao extrai o numero -- portao vazio"
fi

#------------------------------------------------------- 4. MPI dentro da imagem
titulo "4. MPI dentro da imagem: um no', dois ranks"
D=$TMP/caso1
mkdir -p "$D"
# Caso minimo: o proprio exemplo 3D com malha pequena.  O que se verifica NAO e'
# a fisica -- os oraculos acima fazem isso -- e' que dois ranks conversam e que a
# PARTICAO DE VOLUME fecha, que e' onde um MPI mal inicializado se denuncia.
cat > "$D/roda.sh" <<'EOS'
set -eu
cd /opt/HigFlow/higflow/example3d_VOF_rising
test -x ns-example-3d
mpirun --allow-run-as-root -n 2 ./ns-example-3d \
   input/sobe.load /tmp/o.save /tmp/p.print \
   -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 2>&1 | tail -40
EOS
SA=$(ap bash "$D/roda.sh" 2>&1) || true
if echo "$SA" | grep -aq 'Step'; then
  ok "dois ranks avancaram" "$(echo "$SA" | grep -ac 'Step') passo(s)"
else
  nao "dois ranks nao avancaram" "ver saida abaixo"
  echo "$SA" | tail -8 | sed 's/^/            /'
fi
# a invariante: volume da bolha reportado e finito
LB=$(echo "$SA" | grep -a 'BOLHA 3D' | tail -1)
if [ -n "$LB" ]; then
  echo "$LB" | grep -aqiE 'nan|inf' && nao "volume com nao-finito" "$LB" \
                                    || ok "volume finito em 2 ranks" "${LB:0:70}"
else
  nao "diagnostico da bolha ausente" "portao vazio: o exemplo nao chegou ao fim"
fi
# e nenhum aviso de inicializacao de MPI
echo "$SA" | grep -aiqE 'ompi_mpi_init|MPI_ABORT|orte_|pmix.*error' \
  && nao "MPI reclamou na inicializacao" "$(echo "$SA" | grep -aiE 'ompi_mpi_init|MPI_ABORT|pmix.*error' | head -1 | cut -c1-70)" \
  || ok "MPI inicializou sem queixa"

#------------------------------------------- 5. o PETSc da imagem e' o que roda
titulo "5. o binario liga o PETSc da imagem, nao outro"
L=$(ap ldd /opt/HigFlow/higflow/example3d_VOF_rising/ns-example-3d 2>&1)
echo "$L" | grep -q 'not found' \
  && nao "biblioteca nao resolvida" "$(echo "$L" | grep 'not found' | head -1)" \
  || ok "todas as bibliotecas resolvem"
echo "$L" | grep -aq '/opt/petsc-3.25.4' \
  && ok "libpetsc vem de /opt/petsc-3.25.4" \
  || nao "libpetsc vem de outro lugar" "$(echo "$L" | grep -a petsc | head -1)"
# um unico libmpi -- dois (OpenMPI e MPICH) abortam em tempo de execucao
nm=$(echo "$L" | grep -ac 'libmpi' || true)
[ "$nm" -le 2 ] && ok "um unico MPI ligado" "$nm entrada(s) libmpi*" \
                || nao "mais de um MPI ligado" "$nm entradas -- aborta em execucao"

#------------------------------------------------ 6. dois nos, so' no Euler  [E2]
if [ "$EULER" -eq 1 ]; then
  titulo "6. dois nos pelo Slurm (lancamento hibrido, PMIx)"
  command -v srun >/dev/null || { nao "srun ausente" "nao e' no' de submissao?"; }
  if command -v srun >/dev/null; then
    cat > "$TMP/dois.sh" <<'EOS'
set -eu
cd /opt/HigFlow/higflow/example3d_VOF_rising
./ns-example-3d input/sobe.load /tmp/o.save /tmp/p.print \
   -ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10 2>&1 | tail -30
EOS
    SB=$(srun --nodes=2 --ntasks=2 --mpi=pmix \
           apptainer exec "$SIF" bash "$TMP/dois.sh" 2>&1) || true
    if echo "$SB" | grep -aq 'Step'; then
      ok "dois NOS avancaram" "$(echo "$SB" | grep -ac 'Step') passo(s)"
    else
      nao "dois nos nao avancaram" "e' aqui que PMIx/interconexao aparecem"
      echo "$SB" | tail -10 | sed 's/^/            /'
    fi
    # a interconexao: se caiu para TCP, o OpenMPI costuma avisar
    echo "$SB" | grep -aiqE 'btl_tcp|falling back|no.*(openib|verbs|ucx)' \
      && nao "possivel queda para TCP" "conferir [E3]: verbs/UCX da imagem x host" \
      || ok "nenhum aviso de queda de transporte"
  fi
else
  echo; echo "=== 6. dois nos: NAO EXECUTADO (rode com --euler no no' de submissao) ==="
  echo "     E' o portao que descobre PMIx [E2] e interconexao [E3]."
  echo "     Sem ele, esta verificacao NAO cobre execucao entre nos."
fi

#---------------------------------------------------------------------- balanco
echo
echo "=============================================================="
printf " %d portao(oes) avaliado(s), %d falha(s)\n" "$portoes" "$falhas"
[ "$EULER" -eq 1 ] || echo " COBERTURA PARCIAL: faltou o portao de dois nos."
echo "=============================================================="
exit $([ "$falhas" -eq 0 ] && echo 0 || echo 1)
