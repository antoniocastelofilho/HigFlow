# Teste dos operadores de transferencia.  Fora de qualquer build do projeto:
#
#     set -a; . ../../../varsrc; set +a
#     export HIGTREE_DIR=$PWD/../../../higtree
#     make DIM=2 && mpirun -n 3 ./teste-transferencia
DIM ?= 2
FI  = ../../src/hig-flow-fronteira-imersa.c

teste-transferencia: teste-transferencia.c $(FI)
	mpic++ -x c++ -std=gnu++17 -Wall -DDIM=$(DIM) $^ -o $@ \
	  -I$(HIGTREE_DIR)/src -I../../src \
	  -I$(PETSC_DIR)/include -I$(PETSC_DIR)/../include \
	  $(shell pkg-config --cflags glib-2.0 hdf5 libfyaml) \
	  -L$(HIGTREE_DIR)/lib -L$(PETSC_DIR)/lib -Wl,-rpath,$(PETSC_DIR)/lib \
	  -lhig$(DIM)d $(shell pkg-config --libs glib-2.0 hdf5 libfyaml) \
	  -lpetsc -ltrilinos_zoltan -lm

clean:
	$(RM) -f teste-transferencia

.PHONY: clean
