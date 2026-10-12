# Construcao do t8Flow.
#
# DECISAO 5 (doc/projeto-t8flow.md): C++ na norma do t8code.  O t8code impoe
# `cxx_std_20` pela INTERFACE do alvo (target_compile_features ... PUBLIC), e com
# as extensoes GNU no padrao isso e' `-std=gnu++20`.  Uma norma so', dos dois
# lados: nao ha' `extern "C"` a escrever entre o nucleo herdado e o codigo novo.
#
# E NAO HA' `-fpermissive`, de proposito.  A HiGTree precisa dele porque compila C
# como C++ (conversoes implicitas void* -> T*).  Medido: sem a flag, em gnu++20,
# 15 das 16 unidades do nucleo herdado ficam limpas -- so' o `solver-sor.c` falha,
# com 11 erros, e ele e' backend opcional que nao roda em producao.  Herdar a flag
# seria desligar uma checagem do compilador para sempre por causa desse arquivo.
#
# Uso, a partir de t8flow/:
#
#     make DIM=2 T8CODE=<prefixo> HIGTREE_SRC=<raiz>/higtree/src PETSC_DIR=<...>
#
# O driver ci/run_t8flow_tests.py monta tudo isso.  Os caminhos TEM de ser
# absolutos: caminho relativo repassado a um -C de make e' resolvido a partir do
# diretorio errado, e o sintoma e' "nao construido", que e' falha silenciosa.

# O `clean` nao exige DIM, e isto nao e' detalhe: o driver limpa ANTES de saber a
# dimensao, e um clean que falha deixa objetos da dimensao anterior no lugar --
# que e' justamente a falha silenciosa que a limpeza por dimensao previne (os
# objetos nao carregam o DIM no nome).  A primeira versao exigia DIM aqui, o
# clean do driver falhava com aviso, e a guarda ficava inoperante.
ifeq ($(filter clean,$(MAKECMDGOALS)),)
ifndef DIM
$(error Defina DIM=2 ou DIM=3)
endif
endif

# O t8code e o PETSc vivem em `bibliotecas/`, que o git IGNORA -- logo nao estao
# neste worktree.  Por isso os prefixos sao parametros, e nao caminhos relativos.
T8CODE      ?= /home/castelo/HiGFlow-System/HigFlow/bibliotecas/t8code/install
HIGTREE_SRC ?= /home/castelo/HiGFlow-System/HigFlow/higtree/src
PETSC_ROOT  ?= /home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4
PETSC_ARCH_DIR ?= $(PETSC_ROOT)/x86_64

OBJ = obj/$(DIM)d

CXX = mpic++

# A GUARDA ANTES DE QUALQUER REGRA: sem o t8code o erro sai como cabecalho
# faltando tres niveis abaixo, apontando para o libsc.  Dito aqui, aponta para a
# causa.
ifeq ($(filter clean,$(MAKECMDGOALS)),)
ifeq ($(wildcard $(T8CODE)/include/t8.h),)
$(error t8code nao encontrado em $(T8CODE) -- passe T8CODE=<prefixo de instalacao>)
endif
ifeq ($(wildcard $(HIGTREE_SRC)/wls.h),)
$(error nucleo herdado nao encontrado em $(HIGTREE_SRC) -- passe HIGTREE_SRC=<raiz>/higtree/src)
endif
ifeq ($(wildcard $(PETSC_ROOT)/include/petsc.h),)
$(error PETSc nao encontrado em $(PETSC_ROOT) -- passe PETSC_ROOT=<raiz do petsc>)
endif
endif

CXXFLAGS  = -std=gnu++20 -g -O1 -Wall -DDIM=$(DIM)
# As duas definicoes abaixo NAO sao opcionais: sem elas o `t8.h` para com "MPI
# configured differently in t8code and libsc", que aponta para o libsc e nao para
# a definicao que falta.
CXXFLAGS += -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1
CXXFLAGS += -I$(T8CODE)/include
CXXFLAGS += -I$(HIGTREE_SRC) -Isrc -Itests
# Os DOIS diretorios do PETSc: o `petsc.h` fica na raiz do pacote e o
# `petscconf.h` no diretorio do PETSC_ARCH.  Com um so', o erro e' cabecalho
# ausente no meio da cadeia de includes do PETSc.
CXXFLAGS += -I$(PETSC_ROOT)/include -I$(PETSC_ARCH_DIR)/include
CXXFLAGS += -I/usr/include/hypre
CXXFLAGS += $(shell pkg-config --cflags glib-2.0)

LDFLAGS = -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LDFLAGS+= -L$(PETSC_ARCH_DIR)/lib -Wl,-rpath,$(PETSC_ARCH_DIR)/lib
LIBS    = -lt8 -lp4est -lsc -lpetsc -lHYPRE
LIBS   += $(shell pkg-config --libs glib-2.0) -lm

# O NUCLEO HERDADO: as unidades independentes de malha, compiladas AQUI e nao
# linkadas da libhig*.a -- a biblioteca traz a estrutura, que e' justamente o que
# o t8Flow substitui.  Lista derivada da secao 1.3 do doc/projeto-t8flow.md.
#
# NOTA DE PROJETO, achada no primeiro link: o `solver.c` referencia
# `_slv_create_hypre` e `_slv_create_debug_write` SEM guarda de macro, enquanto o
# SOR e o ViennaCL ficam atras de `#ifdef USE_SOR` / `USE_VIENNACL`.  Entao nao da'
# para levar so' o PETSc: o link exige os tres motores.  Por isso o hypre entra
# aqui apesar de nao ser o motor de producao.  Fechar isso -- uma guarda por motor,
# como o SOR ja' tem -- e' uma limpeza barata e esta' na lista do t8Flow.
NUCLEO_SRC = wls.c coord.c utils.c allocator.c rect.c point-cloud.c \
             solver.c solver-petsc.c solver-hypre.c solver-debug-write.c solver-centralized.c
NUCLEO_CXX = rng.cpp
NUCLEO_OBJ = $(addprefix $(OBJ)/,$(NUCLEO_SRC:.c=.o) $(NUCLEO_CXX:.cpp=.o))

TESTS = tests/test-nucleo-herdado
BINS  = $(addsuffix -$(DIM)d,$(TESTS))

all: $(BINS)

# As unidades herdadas sao C compilado como C++ -- `-x c++` e' o que a HiGTree ja'
# faz, e e' o que mantem UMA norma no sistema todo.
$(OBJ)/%.o: $(HIGTREE_SRC)/%.c | $(OBJ)
	$(CXX) $(CXXFLAGS) -x c++ -c $< -o $@

$(OBJ)/%.o: $(HIGTREE_SRC)/%.cpp | $(OBJ)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ):
	mkdir -p $(OBJ)

# O binario DEPENDE dos objetos do nucleo, nao so' do .cxx: sem isso o make o
# considera atualizado quando so' o nucleo mudou, o teste roda obsoleto e da' os
# numeros da versao anterior.  Aconteceu na HiGTree e quase enterrou uma correcao
# como "nao funciona".
tests/%-$(DIM)d: tests/%.cxx tests/testing.h $(NUCLEO_OBJ)
	$(CXX) $(CXXFLAGS) $< $(NUCLEO_OBJ) -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) -rf obj
	$(RM) -f $(addsuffix -2d,$(TESTS)) $(addsuffix -3d,$(TESTS))
	@echo "limpo: objetos das DUAS dimensoes e binarios"

.PHONY: all clean
