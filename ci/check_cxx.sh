#!/usr/bin/env bash
#===============================================================================
# Verifica que os modulos do higflow continuam compilando como C++.
#
# A build ainda usa gcc.  Este teste nao troca isso: ele roda o g++ apenas como
# analisador (-fsyntax-only), porque o compilador C++ recusa coisas que o C
# aceita calado.  Foi assim que apareceu o bug em que quatro dos seis tipos de
# equacao eram atribuidos ao campo errado do struct e ficavam inalcancaveis.
#
# Sem esta checagem, nada impede que a proxima conversao implicita entre tipos
# volte a passar despercebida, e o trabalho da Fase 1 se perde em silencio.
#
# Uso:  ./ci/check_cxx.sh [DIM]      (DIM padrao: 2)
#===============================================================================
set -uo pipefail

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIM="${1:-2}"
PETSC="${PETSC_DIR:-$RAIZ/bibliotecas/petsc-3.25.4/x86_64}"

INCLUDES=(
  $(pkg-config --cflags glib-2.0 hdf5 libfyaml 2>/dev/null)
  -I"$RAIZ/higtree/src" -I"$RAIZ/higflow/src"
  -I/usr/include/trilinos
  $(mpicc -show 2>/dev/null | tr ' ' '\n' | grep '^-I' || true)
  -I"${PETSC%/x86_64}/include" -I"$PETSC/include"
)

# Mesma lista que o Makefile compila, lida dele para nao divergir com o tempo.
mapfile -t MODULOS < <(
  awk '/^MODULES/,0' "$RAIZ/higflow/Makefile" \
    | awk '/^$/{exit} {print}' | tr -d '\\' | tr -s ' \t' '\n' \
    | grep -vE '^$|MODULES|='
)

# Modulos do higtree, lidos do Makefile dele pelo mesmo motivo.
mapfile -t HT < <(
  awk '/^MODULES/,0' "$RAIZ/higtree/Makefile" \
    | awk '/^$/{exit} {print}' | tr -d '\\' | tr -s ' \t' '\n' \
    | grep -vE '^$|MODULES|='
)

falhas=0
for m in "${HT[@]}"; do
  f="$RAIZ/higtree/src/$m.c"
  [ -f "$f" ] || continue
  saida=$(g++ -DDIM="$DIM" -std=gnu++17 -fsyntax-only -w -x c++ \
              "${INCLUDES[@]}" "$f" 2>&1)
  n=$(grep -c 'error:' <<<"$saida")
  if [ "$n" -ne 0 ]; then
    printf '%-46s %d erro(s)\n' "$m.c" "$n"
    grep 'error:' <<<"$saida" | head -3 | sed 's/^/    /'
    falhas=$((falhas + n))
  fi
done

for m in "${MODULOS[@]}"; do
  f="$RAIZ/higflow/src/hig-flow-$m.c"
  [ -f "$f" ] || continue
  saida=$(g++ -DDIM="$DIM" -std=gnu++17 -fsyntax-only -w -x c++ \
              "${INCLUDES[@]}" "$f" 2>&1)
  n=$(grep -c 'error:' <<<"$saida")
  if [ "$n" -ne 0 ]; then
    printf '%-46s %d erro(s)\n' "hig-flow-$m.c" "$n"
    grep 'error:' <<<"$saida" | head -3 | sed 's/^/    /'
    falhas=$((falhas + n))
  fi
done

if [ "$falhas" -eq 0 ]; then
  echo "Os ${#HT[@]} modulos do higtree e ${#MODULOS[@]} do higflow compilam como C++ (DIM=$DIM)."
  exit 0
fi
echo "---"
echo "$falhas erro(s) ao compilar como C++ (DIM=$DIM)."
exit 1
