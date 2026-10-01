#!/bin/bash
# TESTE DE LINK do caminho multifasico VISCOELASTICO em DIM=3.
#
# A pergunta e' binaria: todo simbolo que esse caminho precisa resolve em tres
# dimensoes?  Se sim, o que falta e' CASO e VERIFICACAO, nao codigo.  Se nao, os
# simbolos ausentes nomeiam exatamente o que nunca foi portado.
#
# O QUE ELE PROVA: que nenhuma funcao falta.
# O QUE ELE NAO PROVA: que o resultado esta' certo, nem que roda sem quebrar.
#   Isso exige um caso, e e' o proximo degrau.
#
# Uso:  ./testa-link-3d.sh          (da raiz da arvore de trabalho)
set -e
RAIZ=$(cd "$(dirname "$0")/../../.." && pwd)
# GUARDA O PETSC_DIR DO CHAMADOR ANTES do varsrc: varsrc o redefine para
# $(pwd)/bibliotecas/..., e o worktree nasce SEM bibliotecas/ (nao e'
# versionado).  Sem guardar, o valor passado por quem chama era descartado em
# silencio e o ligador nao achava -lpetsc.
PETSC_PEDIDO="${PETSC_DIR:-}"
cd "$RAIZ" && set -a && . ./varsrc && set +a
[ -n "$PETSC_PEDIDO" ] && export PETSC_DIR="$PETSC_PEDIDO"
[ -d "$PETSC_DIR/lib" ] || { echo "PETSC_DIR=$PETSC_DIR sem lib/; aponte para um PETSc compilado" >&2; exit 2; }

S=$RAIZ/higflow/src
T=$RAIZ/higflow/tests/vof-3d
INC="$(pkg-config --cflags glib-2.0 libfyaml hdf5) -I$PETSC_DIR/include -I$RAIZ/higtree/src -I$S"
CXX="mpic++ -march=native -mtune=native -O2 -DDIM=3 -x c++ -std=gnu++17"

$CXX $INC -c $T/stub-link.c -o $T/stub-link.o

# COMPILA OS PROPRIOS OBJETOS, em diretorio separado.  Reusar os .o de
# higflow/src torna o teste refem da dimensao com que a arvore foi construida
# por ultimo -- e eles CARREGAM a dimensao.  Uma passagem do teste poderia
# significar apenas que alguem acabou de construir um exemplo 3D, e uma falha,
# que acabou de construir um 2D.  Nenhum dos dois e' o que se quer medir.
#
# hig-flow-vof-9-cells e hig-flow-vof-elvira ainda precisam de mencao: eles NAO
# estao na lista MODULES do Makefile de higflow, e sao compilados por exemplo.
O=$T/obj3d
mkdir -p $O
MODS="hig-flow-kernel hig-flow-io hig-flow-eval hig-flow-discret
 hig-flow-bc hig-flow-ic hig-flow-terms hig-flow-step
 hig-flow-step-multiphase hig-flow-vof-adap-hf hig-flow-vof-plic
 hig-flow-vof-9-cells hig-flow-vof-finite-difference-normal-curvature
 hig-flow-vof-elvira hig-flow-step-generalized-newtonian
 hig-flow-step-viscoelastic hig-flow-mittag-leffler
 hig-flow-step-multiphase-viscoelastic hig-flow-res
 hig-flow-linear-algebra hig-flow-viscoelastic-kernel
 hig-flow-vof-mehta hig-flow-vof-plic-3D hig-flow-vof-advection-3D
 hig-flow-vof-HF-3D hig-flow-vof-finite-difference-normal-curvature_3D
 hig-flow-timestep hig-flow-remalha"
for f in $MODS; do
    [ -f $O/$f.o ] || $CXX $INC -c $S/$f.c -o $O/$f.o
done

OBJ="$T/stub-link.o"
for f in $MODS; do OBJ="$OBJ $O/$f.o"; done

mpic++ -flto=8 -fno-fat-lto-objects $OBJ -o $T/teste-link \
  -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu/openmpi/lib \
  -L/usr/local/lib -L$RAIZ/higtree/lib \
  -lhig3d -lglib-2.0 -lhdf5 -lfyaml -lmpi \
  -Wl,-rpath,$PETSC_DIR/lib -L$PETSC_DIR/lib \
  -lpetsc -lHYPRE -lopenblas -lm -lX11 -lstdc++ -lquadmath -ltrilinos_zoltan -lrt

$T/teste-link
