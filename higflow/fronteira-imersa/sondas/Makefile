# Sondas de viabilidade da fronteira imersa.  NAO entram em build nenhum do
# projeto: sao programas de uma pergunta so', rodados a mao.
#
#     set -a; . ../../../varsrc; set +a
#     make && make rodar
#
# Elas nao dependem da HiGTree nem da HiGFlow de proposito -- o que se esta'
# julgando e' a biblioteca, nao a integracao.

SONDAS = sonda-dmswarm sonda-dmplex sonda-petscsf sonda-plex-sobreposicao

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

clean:
	$(RM) -f $(SONDAS)

.PHONY: all rodar clean
