#!/usr/bin/env bash
#
# Instalador do HiGFlow com o t8code.
#
# O QUE ELE FAZ:  constroi o t8code a partir do fonte em bibliotecas/, depois a
#                 HiGTree, a HiGFlow e os exemplos, e opcionalmente roda as suites.
#
# O QUE ELE NAO FAZ, de proposito:  instalar PETSc ou pacotes do sistema.  Ele
#                 CONFERE cada dependencia e diz exatamente o que falta e por que.
#                 Compilar PETSc dentro de um script falha de formas dificeis de
#                 diagnosticar -- depende de BLAS/LAPACK e MPI do sistema --, e um
#                 instalador que falha no meio deixa a arvore pior do que achou.
#
# ARMADILHAS QUE ELE EVITA, todas custaram tempo nesta arvore:
#   - build que "passa" sem produzir nada: aqui todo passo e' conferido por STATUS
#     DE SAIDA e, quando ha' artefato, pela EXISTENCIA do arquivo.  A guarda do
#     Makefile da HiGTree para com `*** Biblioteca ausente`, que nao casa com
#     `grep error:` -- quem confere por grep roda o binario velho e nao percebe.
#   - trocar de dimensao sem `make clean`: os objetos nao carregam a dimensao no
#     nome, entao o make os considera atuais e o binario 2D liga objetos 3D.  Isso
#     falha como segfault e numeros errados, nao como erro de compilacao.
#   - clonar o t8code com `--depth 1`: sem as tags o CMake para com
#     `VERSION ".." format invalid`, mensagem que nao menciona tag nenhuma.
#
# Uso:  ./instalar.sh [opcoes]
#       ./instalar.sh --ajuda

set -euo pipefail

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$RAIZ"

T8_VERSAO="v4.0.0-26.09"
T8_PREFIXO="$RAIZ/bibliotecas/t8code/install"

DIM=2
TRABALHOS="$(nproc 2>/dev/null || echo 4)"
COM_T8=1
REFAZER_T8=0
VERIFICAR=0
SO_CONFERIR=0

vermelho() { printf '\033[31m%s\033[0m\n' "$*"; }
verde()    { printf '\033[32m%s\033[0m\n' "$*"; }
amarelo()  { printf '\033[33m%s\033[0m\n' "$*"; }
titulo()   { printf '\n\033[1m== %s\033[0m\n' "$*"; }

ajuda() {
  cat <<'FIM'
Instalador do HiGFlow com o t8code.

  --dim N            dimensao a construir (2 ou 3).  Padrao: 2
  --jobs N           paralelismo do make/cmake.  Padrao: nproc
  --sem-t8code       nao constroi o t8code; o HiGFlow fica no caminho AMR
  --refazer-t8code   reconstroi o t8code mesmo que ja' exista instalado
  --verificar        ao final, roda as duas suites (demora ~20 min)
  --so-conferir      so' confere as dependencias e sai, sem construir nada
  --ajuda            esta mensagem

Depois de instalar, carregue o ambiente ANTES de rodar qualquer coisa:

    set -a; . ./varsrc; set +a

O varsrc usa $(pwd): carregue-o SEMPRE da raiz do repositorio, ou o PETSC_DIR
aponta para um diretorio que nao existe e a ligacao falha com
"undefined reference to PetscFinalize" -- que acusa o PETSc, e nao o diretorio.
FIM
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dim)          DIM="${2:?--dim exige um valor}"; shift 2 ;;
    --jobs)         TRABALHOS="${2:?--jobs exige um valor}"; shift 2 ;;
    --sem-t8code)   COM_T8=0; shift ;;
    --refazer-t8code) REFAZER_T8=1; shift ;;
    --verificar)    VERIFICAR=1; shift ;;
    --so-conferir)  SO_CONFERIR=1; shift ;;
    --ajuda|-h)     ajuda; exit 0 ;;
    *) vermelho "opcao desconhecida: $1"; echo; ajuda; exit 2 ;;
  esac
done

case "$DIM" in 2|3) ;; *) vermelho "--dim tem de ser 2 ou 3 (veio '$DIM')"; exit 2 ;; esac

# ---------------------------------------------------------------- dependencias
titulo "Dependencias"
FALTAM=0

exige_programa() {
  if command -v "$1" >/dev/null 2>&1; then
    printf '  %-22s %s\n' "$1" "$(command -v "$1")"
  else
    vermelho "  $1  AUSENTE -- $2"; FALTAM=1
  fi
}

exige_pkgconfig() {
  if pkg-config --exists "$1" 2>/dev/null; then
    printf '  %-22s %s\n' "$1" "$(pkg-config --modversion "$1" 2>/dev/null)"
  else
    vermelho "  $1  AUSENTE -- $2"; FALTAM=1
  fi
}

exige_programa mpic++   "OpenMPI.  Debian/Ubuntu: apt install libopenmpi-dev"
exige_programa make     "build-essential"
exige_programa git      "git"
exige_programa pkg-config "pkg-config"
[ "$COM_T8" = 1 ] && exige_programa cmake "CMake >= 3.16.  Debian/Ubuntu: apt install cmake"

exige_pkgconfig glib-2.0 "Debian/Ubuntu: apt install libglib2.0-dev"
exige_pkgconfig hdf5     "Debian/Ubuntu: apt install libhdf5-openmpi-dev"
exige_pkgconfig libfyaml "Debian/Ubuntu: apt install libfyaml-dev"

# Zoltan vem do Trilinos e nao tem .pc.  PERGUNTA-SE AO LIGADOR, que e' quem
# decide -- e nao ao `ldconfig -p`.
#
# Por que nao o ldconfig: sob `set -o pipefail`, `ldconfig -p | grep -q ...` falha
# mesmo com o grep casando, porque o proprio ldconfig sai com status nao-zero
# nesta maquina.  O verificador dizia AUSENTE para uma biblioteca presente, que e'
# o pior defeito possivel num verificador: manda instalar o que ja' existe.
if echo 'int main(){return 0;}' | mpic++ -x c++ - -ltrilinos_zoltan -o /tmp/.higflow-zoltan-teste >/dev/null 2>&1; then
  printf '  %-22s %s\n' "trilinos_zoltan" "o ligador encontra"
  rm -f /tmp/.higflow-zoltan-teste
else
  vermelho "  trilinos_zoltan  o ligador NAO encontra -- Debian/Ubuntu: apt install libtrilinos-zoltan-dev"
  FALTAM=1
fi

# PETSc: o varsrc e' quem define PETSC_DIR/PETSC_ARCH.
if [ -f "$RAIZ/varsrc" ]; then
  # shellcheck disable=SC1091
  ( set -a; . "$RAIZ/varsrc"; set +a
    # O varsrc deste repositorio poe o arch DENTRO do PETSC_DIR e deixa
    # PETSC_ARCH vazio -- entao nao se imprime o arch, que confundiria.
    if [ -n "${PETSC_DIR:-}" ] && [ -d "${PETSC_DIR}/${PETSC_ARCH:-}/lib" ]; then
      printf '  %-22s %s\n' "PETSc" "$PETSC_DIR"
    else
      printf '\033[31m  %-22s %s\033[0m\n' "PETSc" "NAO ENCONTRADO em ${PETSC_DIR:-<vazio>}/${PETSC_ARCH:-<vazio>}"
      echo "     O varsrc define PETSC_DIR a partir de \$(pwd).  Este instalador o carrega"
      echo "     da raiz, entao se falhou aqui e' porque o PETSc nao esta' construido."
      exit 1
    fi )
else
  vermelho "  varsrc  AUSENTE na raiz do repositorio"; FALTAM=1
fi

if [ "$FALTAM" = 1 ]; then
  echo
  vermelho "Faltam dependencias.  Instale-as e rode de novo."
  echo "Nada foi construido -- um instalador que para no meio deixa a arvore pior"
  echo "do que a achou."
  exit 1
fi
verde "  todas presentes"

[ "$SO_CONFERIR" = 1 ] && { echo; verde "Conferencia so'; nada construido."; exit 0; }

# shellcheck disable=SC1091
set -a; . "$RAIZ/varsrc"; set +a
export HIGTREE_DIR="$RAIZ/higtree"
export HIGFLOW_DIR="$RAIZ/higflow"

# ---------------------------------------------------------------------- t8code
if [ "$COM_T8" = 1 ]; then
  titulo "t8code $T8_VERSAO"
  if [ -f "$T8_PREFIXO/lib/libt8.so" ] && [ "$REFAZER_T8" = 0 ]; then
    echo "  ja' instalado em $T8_PREFIXO"
    echo "  (use --refazer-t8code para reconstruir)"
  else
    mkdir -p "$RAIZ/bibliotecas"
    if [ ! -d "$RAIZ/bibliotecas/t8code/.git" ]; then
      echo "  clonando (COM tags: sem elas o CMake para com VERSION \"..\" format invalid)"
      git clone --quiet https://github.com/DLR-AMR/t8code.git "$RAIZ/bibliotecas/t8code"
    fi
    ( cd "$RAIZ/bibliotecas/t8code"
      git fetch --quiet --tags
      git checkout --quiet "$T8_VERSAO"
      echo "  configurando"
      cmake -B build -S . \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$T8_PREFIXO" \
        -DT8CODE_ENABLE_MPI=ON \
        -DT8CODE_BUILD_AS_SHARED_LIBRARY=ON \
        -DT8CODE_BUILD_TESTS=OFF -DT8CODE_BUILD_EXAMPLES=OFF \
        -DT8CODE_BUILD_BENCHMARKS=OFF -DT8CODE_BUILD_TUTORIALS=OFF >/dev/null
      echo "  compilando (-j$TRABALHOS); leva alguns minutos"
      cmake --build build -j"$TRABALHOS" >/dev/null
      cmake --install build >/dev/null )
    [ -f "$T8_PREFIXO/lib/libt8.so" ] || { vermelho "  t8code: a instalacao nao produziu libt8.so"; exit 1; }
  fi
  verde "  ok"
fi

# ------------------------------------------------------------------- construir
#
# `passo` confere o STATUS DE SAIDA e, quando ha' artefato, a EXISTENCIA dele.
# Conferir por `grep error:` na saida deixa passar a guarda do Makefile, que usa
# `***`, e ai' o build "passa" sem produzir binario.
passo() {
  local descricao="$1" artefato="$2"; shift 2
  printf '  %-46s' "$descricao"
  if "$@" >/tmp/higflow-instalar.log 2>&1; then
    if [ -n "$artefato" ] && [ ! -e "$artefato" ]; then
      vermelho "FALHOU (sem artefato)"
      echo "     esperado: $artefato"
      tail -15 /tmp/higflow-instalar.log | sed 's/^/     /'
      exit 1
    fi
    verde "ok"
  else
    vermelho "FALHOU"
    tail -20 /tmp/higflow-instalar.log | sed 's/^/     /'
    exit 1
  fi
}

titulo "HiGTree e HiGFlow (DIM=$DIM)"
# O clean NAO e' opcional: ver o cabecalho.
passo "higtree clean"            ""  make -C "$RAIZ/higtree" clean
passo "higtree DIM=$DIM"         "$RAIZ/higtree/lib/libhig${DIM}d.a" \
                                 make -C "$RAIZ/higtree" -j"$TRABALHOS" "DIM=$DIM"
passo "higflow clean"            ""  make -C "$RAIZ/higflow" clean
passo "higflow DIM=$DIM"         "$RAIZ/higflow/src/hig-flow-kernel.o" \
                                 make -C "$RAIZ/higflow" -j"$TRABALHOS" "DIM=$DIM"

titulo "Exemplos"
# Sao os que a suite exercita nesta dimensao.  Com t8code, os que declaram fontes
# alternativas ganham tambem o caminho do t8code.
if [ "$DIM" = 2 ]; then
  EXEMPLOS=(example2d_Newt example2d_Oldroyd example2d_Gptt example2d_Newt_contraction
            example2d_VOF example2d_VOF_Gptt example2d_VOF_Oldroyd)
else
  EXEMPLOS=(example3d_lid_driven)
fi

ARGS_T8=()
[ "$COM_T8" = 1 ] && ARGS_T8=("T8CODE=$T8_PREFIXO")

for e in "${EXEMPLOS[@]}"; do
  [ -d "$RAIZ/higflow/$e" ] || { amarelo "  $e  (nao existe nesta arvore; pulado)"; continue; }
  bin="$RAIZ/higflow/$e/ns-example"
  [ -f "$RAIZ/higflow/$e/ns-exemple-3d.c" ] && bin="$RAIZ/higflow/$e/ns-exemple-3d"
  passo "$e clean" "" make -C "$RAIZ/higflow/$e" clean
  passo "$e"       "$bin" make -C "$RAIZ/higflow/$e" -j"$TRABALHOS" "${ARGS_T8[@]}"
done

# ------------------------------------------------------------------- verificar
if [ "$VERIFICAR" = 1 ]; then
  titulo "Verificacao"
  echo "  suite da biblioteca (contrato de malha)"
  if [ "$COM_T8" = 1 ]; then
    python3 "$RAIZ/ci/run_higtree_tests.py" --t8code "$T8_PREFIXO" | tail -3 | sed 's/^/     /'
  else
    python3 "$RAIZ/ci/run_higtree_tests.py" | tail -3 | sed 's/^/     /'
  fi
  echo "  suite de exemplos (fisica contra referencia)"
  if [ "$COM_T8" = 1 ]; then
    python3 "$RAIZ/ci/run_suite.py" --t8code "$T8_PREFIXO" | tail -2 | sed 's/^/     /'
  else
    python3 "$RAIZ/ci/run_suite.py" | tail -2 | sed 's/^/     /'
  fi
fi

# ----------------------------------------------------------------------- fecho
titulo "Pronto"
cat <<FIM
  Antes de rodar qualquer coisa, carregue o ambiente DA RAIZ do repositorio:

      cd $RAIZ
      set -a; . ./varsrc; set +a

  Rodar um exemplo pelo caminho de sempre (malha do arquivo AMR):

      python3 ci/run_suite.py --case example2d_Newt

FIM
if [ "$COM_T8" = 1 ]; then
cat <<FIM
  Rodar tambem pelas fontes do t8code, contra a MESMA referencia:

      python3 ci/run_suite.py --t8code $T8_PREFIXO

  E, num exemplo isolado, escolhendo a fonte em tempo de execucao:

      HIGFLOW_MALHA=t8code-particao mpirun -n 2 ./ns-example ...

  As quatro fontes e o limite de cada uma estao em doc/tutorial-mtree-t8code.pdf.
FIM
else
cat <<'FIM'
  Construido SEM t8code: o HiGFlow segue lendo a malha do arquivo AMR, que e' o
  caminho de sempre.  Para ganhar as fontes alternativas, rode de novo sem
  --sem-t8code.
FIM
fi
