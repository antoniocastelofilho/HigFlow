# Construcao do t8Flow.
#
# CODIGO SEPARADO POR DIMENSAO (decisao 7): o t8code tem malhas 2D e 3D nativas,
# com classes de elemento distintas, e o codigo NOVO do t8Flow e' escrito duas
# vezes em vez de uma vez com `#if DIM`.  Daqui sai um alvo por dimensao:
#
#     make 2d          constroi o lado 2D
#     make 3d          ainda nao existe, e o alvo DIZ isso
#     make clean       limpa as duas
#
# O que NAO muda: o nucleo herdado da HiGTree e' parametrizado por macro, porque
# `Point` E' `real[DIM]` (coord.h:31).  Logo as unidades herdadas sao compiladas
# com `-DDIM=2` para o lado 2D e `-DDIM=3` para o 3D, em objetos separados.  Isso
# nao e' custo novo: a HiGTree sempre construiu duas bibliotecas.
#
# NORMA (decisao 5): `-std=gnu++20`, a do t8code, que ele IMPOE pela interface do
# alvo (`target_compile_features(T8 PUBLIC c_std_11 cxx_std_20)`).  Uma norma so',
# dos dois lados: nao ha' `extern "C"` entre o nucleo herdado e o codigo novo.
#
# E NAO HA' `-fpermissive`, de proposito.  Medido: sem a flag, em gnu++20, 15 das
# 16 unidades do nucleo herdado ficam limpas -- so' o `solver-sor.c` falha, e ele
# e' backend opcional que nao roda em producao.
#
# Os caminhos de dependencia TEM de ser absolutos: caminho relativo repassado a um
# `make -C` resolve a partir do diretorio errado, e o sintoma e' "nao
# construido" -- falha silenciosa, nao erro.

# O t8code e o PETSc vivem em `bibliotecas/`, que o git IGNORA -- logo nao estao
# neste worktree.  Por isso sao parametros.
T8CODE      ?= /home/castelo/HiGFlow-System/HigFlow/bibliotecas/t8code/install
HIGTREE_SRC ?= /home/castelo/HiGFlow-System/HigFlow/higtree/src
PETSC_ROOT  ?= /home/castelo/HiGFlow-System/HigFlow/bibliotecas/petsc-3.25.4
PETSC_ARCH_DIR ?= $(PETSC_ROOT)/x86_64

CXX = mpic++

# GUARDAS ANTES DE QUALQUER REGRA: sem o t8code o erro sai como cabecalho
# faltando tres niveis abaixo, apontando para o libsc.  Dito aqui, aponta para a
# causa.  O `clean` nao as sofre.
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

BASE  = -std=gnu++20 -g -O1 -Wall
# As duas definicoes abaixo NAO sao opcionais: sem elas o `t8.h` para com "MPI
# configured differently in t8code and libsc", que aponta para o libsc e nao para
# a definicao que falta.
BASE += -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1
BASE += -I$(T8CODE)/include -I$(HIGTREE_SRC) -Isrc/comum
# Os DOIS diretorios do PETSc: `petsc.h` na raiz do pacote, `petscconf.h` no
# diretorio do PETSC_ARCH.  Com um so', o erro e' cabecalho ausente no meio da
# cadeia de includes do PETSc.
BASE += -I$(PETSC_ROOT)/include -I$(PETSC_ARCH_DIR)/include
BASE += -I/usr/include/hypre
BASE += $(shell pkg-config --cflags glib-2.0)

LDFLAGS = -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LDFLAGS+= -L$(PETSC_ARCH_DIR)/lib -Wl,-rpath,$(PETSC_ARCH_DIR)/lib
LIBS    = -lt8 -lp4est -lsc -lpetsc -lHYPRE
LIBS   += $(shell pkg-config --libs glib-2.0) -lm

# O NUCLEO HERDADO: as unidades independentes de malha, compiladas AQUI e nao
# linkadas da libhig*.a -- a biblioteca traz a estrutura, que e' o que o t8Flow
# substitui.  Lista derivada da secao 1.3 do doc/projeto-t8flow.md.
#
# NOTA, achada no primeiro link: o `solver.c` referencia `_slv_create_hypre` e
# `_slv_create_debug_write` SEM guarda de macro, enquanto o SOR e o ViennaCL ficam
# atras de `#ifdef`.  Entao nao da' para levar so' o PETSc -- o link exige hypre,
# debug-write e centralized.  Uma guarda por motor e' limpeza barata e esta' na
# lista do projeto.
NUCLEO_C   = wls.c coord.c utils.c allocator.c rect.c point-cloud.c \
             solver.c solver-petsc.c solver-hypre.c solver-debug-write.c \
             solver-centralized.c
NUCLEO_CPP = rng.cpp

COMUM_SRC = src/comum/arranque.cxx

# ---------------------------------------------------------------------------
# 2D
# ---------------------------------------------------------------------------
OBJ2 = obj/2d
CXXFLAGS2 = $(BASE) -DDIM=2 -Isrc/2d
SRC2   = src/2d/mls2d.cxx src/2d/solver2d.cxx
OBJ2_N = $(addprefix $(OBJ2)/,$(NUCLEO_C:.c=.o) $(NUCLEO_CPP:.cpp=.o))
OBJ2_T = $(addprefix $(OBJ2)/,$(notdir $(COMUM_SRC:.cxx=.o)) $(notdir $(SRC2:.cxx=.o)))
TESTES2 = $(basename $(wildcard tests/2d/test-*.cxx))
BINS2   = $(TESTES2)

2d: $(BINS2)

$(OBJ2):
	mkdir -p $(OBJ2)

# As unidades herdadas sao C compilado como C++ -- `-x c++` e' o que a HiGTree ja'
# faz, e e' o que mantem UMA norma no sistema todo.
$(OBJ2)/%.o: $(HIGTREE_SRC)/%.c | $(OBJ2)
	$(CXX) $(CXXFLAGS2) -x c++ -c $< -o $@

$(OBJ2)/%.o: $(HIGTREE_SRC)/%.cpp | $(OBJ2)
	$(CXX) $(CXXFLAGS2) -c $< -o $@

$(OBJ2)/%.o: src/comum/%.cxx | $(OBJ2)
	$(CXX) $(CXXFLAGS2) -c $< -o $@

$(OBJ2)/%.o: src/2d/%.cxx | $(OBJ2)
	$(CXX) $(CXXFLAGS2) -c $< -o $@

# O binario DEPENDE dos objetos: sem isso o make o considera atualizado quando so'
# a biblioteca mudou, o teste roda obsoleto e da' os numeros da versao anterior.
# Aconteceu na HiGTree e quase enterrou uma correcao como "nao funciona".
tests/2d/%: tests/2d/%.cxx tests/testing.h $(OBJ2_T) $(OBJ2_N)
	$(CXX) $(CXXFLAGS2) $< $(OBJ2_T) $(OBJ2_N) -o $@ $(LDFLAGS) $(LIBS)

# ---------------------------------------------------------------------------
# 3D -- ainda nao escrito, e o alvo diz isso em vez de falhar obscuro.
# ---------------------------------------------------------------------------
3d:
	@echo "O lado 3D ainda nao existe.  A decisao 7 e' comecar pelo 2D:"
	@echo "o codigo 3D sera' src/3d/ e tests/3d/, escrito e nao gerado por macro."
	@exit 1

clean:
	$(RM) -rf obj
	$(RM) -f $(wildcard tests/2d/test-*[!x]) 
	@echo "limpo: objetos das duas dimensoes e binarios de teste"

.PHONY: 2d 3d clean all
all: 2d
