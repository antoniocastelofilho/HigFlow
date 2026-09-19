# Sonda de viabilidade do t8code, FORA da suite normal.
#
# Nao entra no ci/run_higtree_tests.py porque exige o t8code instalado, que nao
# e' dependencia da HiGTree.  Roda a mao, com T8CODE apontando para o prefixo de
# instalacao:
#
#     make T8CODE=<raiz>/bibliotecas/t8code/install
#     mpirun -n 1 ./probe-c11
#
# As tres definicoes abaixo nao sao opcionais -- ver o cabecalho do probe-c11.c.

T8CODE ?= ../../../bibliotecas/t8code/install

CXX      = mpic++
CXXFLAGS = -x c++ -std=gnu++20 -O1 -Wall
CXXFLAGS+= -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1 -DT8_TEST_LEVEL_INT=0
CXXFLAGS+= -I$(T8CODE)/include
LDFLAGS  = -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LIBS     = -lt8 -lp4est -lsc

probe-c11: probe-c11.c
	@test -d $(T8CODE)/include || { \
	  echo "t8code nao encontrado em $(T8CODE).  Passe T8CODE=<prefixo de instalacao>."; \
	  exit 1; }
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) -f probe-c11

.PHONY: clean
