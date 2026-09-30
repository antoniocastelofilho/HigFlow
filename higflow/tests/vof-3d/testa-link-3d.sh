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

# hig-flow-vof-9-cells e hig-flow-vof-elvira NAO estao na lista MODULES do
# Makefile de higflow: sao compilados por exemplo.  Aqui sao compilados a mao.
for f in hig-flow-vof-9-cells hig-flow-vof-elvira; do
    [ -f $S/$f.o ] || $CXX $INC -c $S/$f.c -o $S/$f.o
done

OBJ="$T/stub-link.o
 $S/hig-flow-kernel.o $S/hig-flow-io.o $S/hig-flow-eval.o $S/hig-flow-discret.o
 $S/hig-flow-bc.o $S/hig-flow-ic.o $S/hig-flow-terms.o $S/hig-flow-step.o
 $S/hig-flow-step-multiphase.o $S/hig-flow-vof-adap-hf.o $S/hig-flow-vof-plic.o
 $S/hig-flow-vof-9-cells.o $S/hig-flow-vof-finite-difference-normal-curvature.o
 $S/hig-flow-vof-elvira.o $S/hig-flow-step-generalized-newtonian.o
 $S/hig-flow-step-viscoelastic.o $S/hig-flow-mittag-leffler.o
 $S/hig-flow-step-multiphase-viscoelastic.o $S/hig-flow-res.o
 $S/hig-flow-linear-algebra.o $S/hig-flow-viscoelastic-kernel.o
 $S/hig-flow-vof-mehta.o $S/hig-flow-vof-plic-3D.o $S/hig-flow-vof-advection-3D.o
 $S/hig-flow-vof-HF-3D.o $S/hig-flow-vof-finite-difference-normal-curvature_3D.o
 $S/hig-flow-timestep.o"

mpic++ -flto=8 -fno-fat-lto-objects $OBJ -o $T/teste-link \
  -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu/openmpi/lib \
  -L/usr/local/lib -L$RAIZ/higtree/lib \
  -lhig3d -lglib-2.0 -lhdf5 -lfyaml -lmpi \
  -Wl,-rpath,$PETSC_DIR/lib -L$PETSC_DIR/lib \
  -lpetsc -lHYPRE -lopenblas -lm -lX11 -lstdc++ -lquadmath -ltrilinos_zoltan -lrt

$T/teste-link
