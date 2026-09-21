# Sondas de viabilidade da fronteira imersa.  NAO entram em build nenhum do
# projeto: sao programas de uma pergunta so', rodados a mao.
#
#     set -a; . ../../../varsrc; set +a
#     make && make rodar
#
# Elas nao dependem da HiGTree nem da HiGFlow de proposito -- o que se esta'
# julgando e' a biblioteca, nao a integracao.

SONDAS = sonda-dmswarm sonda-dmplex sonda-petscsf sonda-plex-sobreposicao

# A sonda da franja e' a unica que precisa da HiGTree -- ela mede a geometria
# entregue pela particao, nao uma capacidade de biblioteca.  Por isso tem regra
# propria, com DIM, e fica FORA de $(SONDAS): `make` nao a constroi.
#
#     make sonda-franja DIM=2      (ou DIM=3, com a libhig da dimensao pronta)
DIM ?= 2

CXX      = mpic++
CXXFLAGS = -x c++ -std=gnu++17 -Wall -I$(PETSC_DIR)/include -I$(PETSC_DIR)/../include
LDFLAGS  = -L$(PETSC_DIR)/lib -Wl,-rpath,$(PETSC_DIR)/lib
LIBS     = -lpetsc

all: $(SONDAS)

%: %.c
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

# np=2 e np=3: um caso par e um impar.  O impar ja' pegou erro de ordenacao na
# integracao do t8code que o par escondia.
rodar: $(SONDAS)
	@for s in $(SONDAS); do \
	  for np in 2 3; do \
	    echo "== $$s np=$$np"; \
	    mpirun -use-hwthread-cpus -n $$np ./$$s | grep "rank"; \
	  done; \
	done

sonda-franja: sonda-franja-euleriana.c
	$(CXX) -x c++ -std=gnu++17 -Wall -DDIM=$(DIM) $< -o $@ \
	  -I$(HIGTREE_DIR)/src -I$(PETSC_DIR)/include -I$(PETSC_DIR)/../include \
	  $(shell pkg-config --cflags glib-2.0 hdf5 libfyaml) \
	  -L$(HIGTREE_DIR)/lib -L$(PETSC_DIR)/lib -Wl,-rpath,$(PETSC_DIR)/lib \
	  -lhig$(DIM)d $(shell pkg-config --libs glib-2.0 hdf5 libfyaml) -lpetsc -ltrilinos_zoltan

clean:
	$(RM) -f $(SONDAS) sonda-franja

.PHONY: all rodar clean
