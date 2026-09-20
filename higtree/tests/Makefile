# Suite de testes da HiGTree, sem dependencia externa.
#
# A suite ATF em ../atf-tests precisa de `atf-c`, que nao esta' instalado em lugar
# nenhum -- 71 casos de 2020 que nao rodam.  Aqui o requisito e' um compilador C, e
# o driver e' o ci/run_higtree_tests.py, no mesmo padrao do ci/run_suite.py.
#
# Uso: make DIM=2   ou   make DIM=3   (a partir de higtree/tests)

# Tres variantes porque o caminho mudou entre versoes do PETSc.  So' a terceira
# resolve no 3.25; a segunda, sozinha, deixava PETSC_LIB VAZIO -- e isso so'
# aparecia ao ligar um teste que referencia PETSc (o paralelo, via
# higtree_initialize), como `undefined reference to PetscFinalize`.
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ../Makefile.common

ifndef DIM
$(error Defina DIM=2 ou DIM=3)
endif

# GUARDA: sem a biblioteca, o make responde `No rule to make target 'test-X-2d'`,
# que se le como FONTE faltando -- e a fonte esta' la'.  Sao duas causas, nenhuma
# delas obvia pela mensagem:
#
#   1. a higtree esta' construida noutra dimensao.  Os objetos nao carregam a
#      dimensao no nome, entao apos um build DIM=3 o libhig2d.a simplesmente nao
#      existe;
#   2. HIGTREE_DIR nao esta' exportado.  Ai' HIGPATH e INSPATH vem vazios,
#      LIBPATH vira `/lib`, e a dependencia nunca resolve.  O `varsrc` NAO define
#      HIGTREE_DIR -- quem define e' o ci/run_higtree_tests.py.  Compilar um teste
#      a mao exige `export HIGTREE_DIR=<raiz>/higtree`.
#
# Custou uma rodada de diagnostico em 2026-09-19, com as duas causas somadas.
ifeq ($(filter clean,$(MAKECMDGOALS)),)
ifndef HIGTREE_DIR
$(error HIGTREE_DIR nao esta' exportado, entao LIBPATH fica vazio e a biblioteca \
nunca resolve. Use `export HIGTREE_DIR=<raiz>/higtree`, ou rode pelo driver \
ci/run_higtree_tests.py)
endif
ifeq ($(wildcard $(LIBPATH)/libhig$(DIM)d.a),)
$(error Biblioteca ausente: $(LIBPATH)/libhig$(DIM)d.a -- a higtree nao esta' \
construida em DIM=$(DIM). Construa com `make -C $(HIGPATH) clean && \
make -C $(HIGPATH) DIM=$(DIM)`. O clean nao e' opcional: os objetos nao carregam \
a dimensao no nome)
endif
endif

# A HiGTree e' compilada como C++ (os simbolos da .a vem mangled), entao o teste
# tem de usar o MESMO compilador -- com gcc o link falha em undefined reference.
CC = mpic++

CFLAGS  = -std=gnu++17 -g -O0 -Wall
# O hdf5 entra porque `lbal.h` inclui `higtree-io.h`.  So' o teste paralelo precisa
# do lbal.h, mas a flag e' inofensiva para os demais e evita um segundo conjunto de
# regras.  Sem ela o erro e' `hdf5.h: No such file or directory`, que aponta para um
# cabecalho que ninguem escreveu no teste.
CFLAGS += $(shell pkg-config --cflags glib-2.0 hdf5)
CFLAGS += $(PETSC_CC_INCLUDES)
CFLAGS += -I$(SRCPATH) -I.
CFLAGS += -DDIM=$(DIM)

LDFLAGS = $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5 libfyaml)
LDFLAGS+= -L$(LIBPATH)

LIBS  = -lhig$(DIM)d
LIBS += $(shell pkg-config --libs-only-l glib-2.0 hdf5 libfyaml)
LIBS += -lhdf5
LIBS += $(PETSC_LIB)
LIBS += $(PETSC_EXTRA_LIB)
# O balanceador (lbal.c) chama o Zoltan; so' o teste paralelo o alcança, mas o
# linker resolve por binario e nao por teste.  O nome da biblioteca no sistema e'
# `trilinos_zoltan` -- nao `zoltan`, e nao `trinilos_zoltan`, que e' um erro de
# digitacao que sobrevive no higflow/Makefile por ali so' se arquivar um .a.
LIBS += -ltrilinos_zoltan
LIBS += -lm

TESTS  = test-stencil-value
TESTS += test-stencil-selection
TESTS += test-partition-independence
TESTS += test-level-jump
TESTS += test-point-location
TESTS += test-boundary-path
TESTS += test-facet-domain-serial
TESTS += test-cell-queries
TESTS += test-fringe-support
TESTS += test-fringe-parallel
TESTS += test-mesh-snapshot

BIN = $(addsuffix -$(DIM)d, $(TESTS))

all: $(BIN)

# A dependencia da BIBLIOTECA e' obrigatoria.  Sem ela o make considera o teste
# atualizado quando so' a higtree mudou, e o binario obsoleto roda em silencio
# dando os numeros da versao anterior -- exatamente o modo de falha que esta suite
# existe para pegar.  Custou uma rodada de diagnostico no dia em que foi escrita.
# ---------------------------------------------------------------------------
# t8code, OPCIONAL.  Passe T8CODE=<prefixo de instalacao> para que o
# test-level-jump ganhe o segundo produtor de malha.  Sem isso nada muda e a
# suite nao depende da biblioteca.
#
# A unidade do produtor compila em C++20 porque o t8code exige, enquanto o resto
# da suite fica em gnu++17 -- so' a assinatura em extern "C" atravessa.
ifdef T8CODE
T8OBJ     = t8code/t8-mesh-producer-$(DIM)d.o t8code/t8-point-locator-$(DIM)d.o
CFLAGS   += -DHIGTREE_COM_T8CODE -I$(T8CODE)/include
LDFLAGS  += -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LIBS     += -lt8 -lp4est -lsc
T8FLAGS   = -x c++ -std=gnu++20 -O1 -DDIM=$(DIM) -I$(SRCPATH) -I. \
            -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1 -DT8_TEST_LEVEL_INT=0 \
            -I$(T8CODE)/include $(shell pkg-config --cflags glib-2.0 hdf5) \
            $(PETSC_CC_INCLUDES)

t8code/%-$(DIM)d.o: t8code/%.cxx t8code/%.h
	$(CC) $(T8FLAGS) -c $< -o $@
endif

%-$(DIM)d: %.c testing.h $(T8OBJ) $(LIBPATH)/libhig$(DIM)d.a
	$(CC) $(CFLAGS) $< $(T8OBJ) -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) -f *.o t8code/*.o test-*-2d test-*-3d

.PHONY: all clean
