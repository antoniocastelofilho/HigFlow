# Suite de testes da HiGTree, sem dependencia externa.
#
# A suite ATF em ../atf-tests precisa de `atf-c`, que nao esta' instalado em lugar
# nenhum -- 71 casos de 2020 que nao rodam.  Aqui o requisito e' um compilador C, e
# o driver e' o ci/run_higtree_tests.py, no mesmo padrao do ci/run_suite.py.
#
# Uso: make DIM=2   ou   make DIM=3   (a partir de higtree/tests)

-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
include ../Makefile.common

ifndef DIM
$(error Defina DIM=2 ou DIM=3)
endif

# A HiGTree e' compilada como C++ (os simbolos da .a vem mangled), entao o teste
# tem de usar o MESMO compilador -- com gcc o link falha em undefined reference.
CC = mpic++

CFLAGS  = -std=gnu++17 -g -O0 -Wall
CFLAGS += $(shell pkg-config --cflags glib-2.0)
CFLAGS += $(PETSC_CC_INCLUDES)
CFLAGS += -I$(SRCPATH) -I.
CFLAGS += -DDIM=$(DIM)

LDFLAGS = $(shell pkg-config --libs-only-L --libs-only-other glib-2.0)
LDFLAGS+= -L$(LIBPATH)

LIBS  = -lhig$(DIM)d
LIBS += $(shell pkg-config --libs-only-l glib-2.0)
LIBS += $(PETSC_LIB)
LIBS += $(PETSC_EXTRA_LIB)
LIBS += -lm

TESTS = test-stencil-value

BIN = $(addsuffix -$(DIM)d, $(TESTS))

all: $(BIN)

# A dependencia da BIBLIOTECA e' obrigatoria.  Sem ela o make considera o teste
# atualizado quando so' a higtree mudou, e o binario obsoleto roda em silencio
# dando os numeros da versao anterior -- exatamente o modo de falha que esta suite
# existe para pegar.  Custou uma rodada de diagnostico no dia em que foi escrita.
%-$(DIM)d: %.c testing.h $(LIBPATH)/libhig$(DIM)d.a
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) -f *.o test-*-2d test-*-3d

.PHONY: all clean
