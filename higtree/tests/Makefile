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

TESTS  = test-stencil-value
TESTS += test-stencil-selection
TESTS += test-partition-independence
TESTS += test-level-jump
TESTS += test-point-location
TESTS += test-boundary-path

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
