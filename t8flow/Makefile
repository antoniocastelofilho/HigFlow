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

# A ESCOLHA DO BACK-END E' DE QUEM COMPILA, pela variavel `GPU`:
#
#     make 2d            AUTO -- detecta OpenCL; se houver, a GPU fica
#                        DISPONIVEL mas o padrao continua CPU
#     make 2d GPU=sim    EXIGE OpenCL e faz a GPU ser o padrao do lote.  Sem
#                        OpenCL, o build PARA com a causa -- quem pediu placa
#                        quer saber que nao ha', e nao receber CPU em silencio
#     make 2d GPU=nao    CPU so', mesmo que haja OpenCL instalado.  O codigo de
#                        placa sai do binario
#
# Tres estados e nao dois: "detectei e deixei disponivel" e' diferente de "o
# usuario pediu", e tratar os dois igual tira dele a escolha.  Com GPU=sim o lote
# comeca em GPU; nos outros dois comeca em CPU, e `usa_back_end()` troca em
# execucao quando houver placa.
GPU ?=
ifeq ($(GPU),nao)
  T8FLOW_GPU_MODO = desligada por GPU=nao
else
  OPENCL_H := $(firstword $(wildcard /usr/include/CL/cl.h))
  OPENCL_L := $(firstword $(wildcard /usr/lib/x86_64-linux-gnu/libOpenCL.so /usr/lib64/libOpenCL.so))
endif

ifeq ($(GPU),sim)
ifeq ($(and $(OPENCL_H),$(OPENCL_L)),)
$(error GPU=sim pedido, mas falta OpenCL: header /usr/include/CL/cl.h ou \
libOpenCL.so ausentes.  Instale `opencl-headers ocl-icd-opencl-dev` -- e o \
RUNTIME tambem (um ICD em /etc/OpenCL/vendors/), porque o loader sozinho nao \
calcula.  Ou use `make 2d` para detectar, ou `make 2d GPU=nao` para CPU so')
endif
endif

# OPENCL E' OPCIONAL, E A DETECCAO E' AUTOMATICA.  A versao de CPU do lote de
# estenceis funciona sem nada disto -- quem nao tem placa nao perde
# funcionalidade, so' velocidade.  Por isso a ausencia de OpenCL NAO e' erro de
# build: a macro simplesmente nao e' definida e o codigo de placa sai do binario.
#
# Exige o HEADER e a biblioteca de LINK (`libOpenCL.so`, nao so' o `.so.1`): o
# runtime sozinho nao compila nada, e detectar pelo `.so.1` daria erro de link no
# fim do build em vez de ausencia detectada no comeco.
ifneq ($(and $(OPENCL_H),$(OPENCL_L)),)
BASE += -DT8FLOW_COM_OPENCL
LIBS_CL = -lOpenCL
ifeq ($(GPU),sim)
BASE += -DT8FLOW_GPU_PADRAO
$(info t8flow: OpenCL encontrado, e GPU=sim -- a GPU e' o PADRAO do lote)
else
$(info t8flow: OpenCL encontrado -- GPU disponivel, padrao continua CPU (use GPU=sim))
endif
else
LIBS_CL =
$(info t8flow: lote so' em CPU -- $(if $(T8FLOW_GPU_MODO),$(T8FLOW_GPU_MODO),sem OpenCL: header ou libOpenCL.so ausentes))
endif

LDFLAGS = -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LDFLAGS+= -L$(PETSC_ARCH_DIR)/lib -Wl,-rpath,$(PETSC_ARCH_DIR)/lib
LIBS    = -lt8 -lp4est -lsc -lpetsc -lHYPRE
LIBS   += $(shell pkg-config --libs glib-2.0) -lm $(LIBS_CL)

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

COMUM_SRC = src/comum/arranque.cxx src/comum/solverlin.cxx

# ---------------------------------------------------------------------------
# A CONFIGURACAO ENTRA NA DEPENDENCIA, e isto nao e' refinamento: sem isso,
# trocar `GPU=sim` por `GPU=nao` NAO reconstroi nada -- os objetos e os binarios
# do build anterior continuam atualizados aos olhos do make, e o programa roda
# com a configuracao ERRADA sem um aviso.
#
# Aconteceu aqui: depois de `make GPU=nao`, um `make GPU=sim` nao recompilou, e o
# teste anunciou "compilado SEM OpenCL" num build pedido COM.  E' a mesma classe
# de "build sem o que se quer medir" que neste projeto ja' custou duas corridas
# medidas como t8code que eram malha uniforme.
#
# O remedio e' um SELO: a configuracao e' gravada em `obj/.config` e conferida no
# parse do Makefile.  Mudou, limpa.  Nao ha' como esquecer o `clean` porque nao
# ha' `clean` a lembrar.
CFG_GPU  := $(if $(filter sim,$(GPU)),sim,$(if $(filter nao,$(GPU)),nao,auto))
CFG_SELO := gpu=$(CFG_GPU) opencl=$(if $(and $(OPENCL_H),$(OPENCL_L)),sim,nao)
CFG_ARQ  := obj/.config
ifeq ($(filter clean,$(MAKECMDGOALS)),)
CFG_ANTES := $(shell cat $(CFG_ARQ) 2>/dev/null)
ifneq ($(CFG_SELO),$(CFG_ANTES))
$(info t8flow: configuracao mudou [$(if $(CFG_ANTES),$(CFG_ANTES),nenhuma) -> $(CFG_SELO)] -- limpando antes)
$(shell rm -rf obj; find tests -maxdepth 2 -type f ! -name '*.cxx' ! -name '*.h' -delete 2>/dev/null; mkdir -p obj; printf '%s' '$(CFG_SELO)' > $(CFG_ARQ))
endif
endif

# ---------------------------------------------------------------------------
# 2D
# ---------------------------------------------------------------------------
OBJ2 = obj/2d
CXXFLAGS2 = $(BASE) -DDIM=2 -Isrc/2d -Iexemplos/poisson

# O QUE E' BIBLIOTECA fica em `src/`: malha, estencil, contorno, campos, o ajuste
# por minimos quadrados, o lote, o invólucro de solver e a escolha de
# discretizacao.  Nada aqui sabe que equacao se esta' resolvendo.
SRC2   = src/2d/mls2d.cxx src/2d/discretizacao2d.cxx \
         src/2d/malha2d.cxx src/2d/campos2d.cxx \
         src/2d/lote2d.cxx src/2d/estencil2d.cxx \
         src/2d/contorno2d.cxx

# O QUE E' EXEMPLO fica em `exemplos/`, e o operador de Poisson com as solucoes
# manufaturadas e' exemplo: ele resolve UMA equacao escolhida, e o catalogo de
# solucoes com derivadas analiticas so' serve para verificar.
#
# CONSEQUENCIA DECLARADA, porque ela vai cobrar: o passo de projecao de
# Navier-Stokes E' uma equacao de Poisson, logo o t8SolverNS vai precisar de um
# operador de Poisson em `src/`.  Ele nascera' da extracao deste -- a parte que
# monta e resolve, sem as funcoes de erro e sem a fonte manufaturada, que sao
# verificacao.
SRC2_EX = exemplos/poisson/poisson2d.cxx exemplos/poisson/manufaturada2d.cxx
OBJ2_N = $(addprefix $(OBJ2)/,$(NUCLEO_C:.c=.o) $(NUCLEO_CPP:.cpp=.o))
OBJ2_T = $(addprefix $(OBJ2)/,$(notdir $(COMUM_SRC:.cxx=.o)) $(notdir $(SRC2:.cxx=.o)))
OBJ2_X = $(addprefix $(OBJ2)/,$(notdir $(SRC2_EX:.cxx=.o)))
TESTES2 = $(basename $(wildcard tests/2d/test-*.cxx))
BINS2   = $(TESTES2)

# OS EXEMPLOS SAO DECLARADOS, e nao varridos.
#
# O wildcard `exemplos/*/*.cxx` pegava as PECAS junto com os programas, e a lista
# de excecoes crescia a cada peca nova -- ate' que o diretorio do lado 3D
# apareceu e o alvo `2d` tentou construir `manufaturada3d` como exemplo seu.
# Um exemplo e' um programa com `main`; as pecas nao sao, e a diferenca nao se
# descobre por padrao de nome.
EXEMPLOS2 = exemplos/poisson/poisson-mms

2d: $(BINS2) $(EXEMPLOS2)

# OS EXEMPLOS entram no alvo `2d` e NAO num alvo separado, de proposito: exemplo
# que nao e' construido junto apodrece.  Eles usam os mesmos objetos dos testes.
exemplos/%: exemplos/%.cxx $(OBJ2_T) $(OBJ2_X) $(OBJ2_N)
	$(CXX) $(CXXFLAGS2) -Itests $< $(OBJ2_T) $(OBJ2_X) $(OBJ2_N) -o $@ $(LDFLAGS) $(LIBS)

# A MEDIDA DE ESCALABILIDADE e' ferramenta, nao teste: usa o t8code direto, roda
# em COMM_WORLD com np>1, e NAO entra na suite -- o que ela produz e' numero para
# decidir, nao veredicto.  Fica fora do alvo `2d` de proposito, para nao gastar
# tempo de suite, e e' construida com `make bench`.
#
# Compila com -DDIM=2 so' porque o cabecalho herdado exige a macro; o programa
# escolhe a dimensao em TEMPO DE EXECUCAO (`--dim=`), porque ele nao usa nada do
# nucleo herdado -- so' o t8code, que tem as duas dimensoes na mesma biblioteca.
bench: bench/escalabilidade bench/franja bench/densidade bench/franja3d

bench/escalabilidade: bench/escalabilidade.cxx
	$(CXX) $(CXXFLAGS2) $< -o $@ $(LDFLAGS) $(LIBS)

# O ARREIO DA FRANJA 3D fala direto com o t8code, como o de escalabilidade: a
# medida vem ANTES do codigo 3D, e nao pode depender de decisoes de arquitetura
# que ela mesma deve informar.
bench/franja3d: bench/franja3d.cxx
	$(CXX) $(CXXFLAGS2) $< -o $@ $(LDFLAGS) $(LIBS)

# O ARREIO DA FRANJA usa as classes do t8Flow -- `Malha`, `Poisson` --, e por isso
# liga os objetos de 2D, ao contrario do de escalabilidade, que fala direto com o
# t8code.  Medir a franja pelo t8code cru nao seria possivel: as camadas acima da
# primeira sao deste codigo, nao da biblioteca.
bench/densidade: bench/densidade.cxx $(OBJ2_T) $(OBJ2_X) $(OBJ2_N)
	$(CXX) $(CXXFLAGS2) $< $(OBJ2_T) $(OBJ2_X) $(OBJ2_N) -o $@ $(LDFLAGS) $(LIBS)

bench/franja: bench/franja.cxx $(OBJ2_T) $(OBJ2_X) $(OBJ2_N)
	$(CXX) $(CXXFLAGS2) $< $(OBJ2_T) $(OBJ2_X) $(OBJ2_N) -o $@ $(LDFLAGS) $(LIBS)

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

$(OBJ2)/%.o: exemplos/poisson/%.cxx | $(OBJ2)
	$(CXX) $(CXXFLAGS2) -c $< -o $@

# O binario DEPENDE dos objetos: sem isso o make o considera atualizado quando so'
# a biblioteca mudou, o teste roda obsoleto e da' os numeros da versao anterior.
# Aconteceu na HiGTree e quase enterrou uma correcao como "nao funciona".
tests/2d/%: tests/2d/%.cxx tests/testing.h $(OBJ2_T) $(OBJ2_X) $(OBJ2_N)
	$(CXX) $(CXXFLAGS2) $< $(OBJ2_T) $(OBJ2_X) $(OBJ2_N) -o $@ $(LDFLAGS) $(LIBS)

# ---------------------------------------------------------------------------
# 3D -- hexaedros.  Base SEPARADA do 2D, nao gerada por macro (decisao 7).
#
# O ESCOPO E' "SO' O QUE O SOLVER VAI USAR", decidido pelo usuario: malha, MLS,
# estencil, contorno, campos e Poisson com os dois esquemas de volumes finitos.
# Ficam de fora os dois compactos (que recusam malha adaptada -- sao contraste,
# nao ferramenta), o laplaciano pontual (instavel) e o lote em GPU.
# ---------------------------------------------------------------------------
OBJ3 = obj/3d
CXXFLAGS3 = $(BASE) -DDIM=3 -Isrc/3d -Iexemplos/poisson3d
SRC3   = src/3d/mls3d.cxx src/3d/malha3d.cxx src/3d/contorno3d.cxx src/3d/estencil3d.cxx src/3d/campos3d.cxx src/3d/discretizacao3d.cxx
OBJ3_N = $(addprefix $(OBJ3)/,$(NUCLEO_C:.c=.o) $(NUCLEO_CPP:.cpp=.o))
OBJ3_T = $(addprefix $(OBJ3)/,$(notdir $(COMUM_SRC:.cxx=.o)) $(notdir $(SRC3:.cxx=.o)))
TESTES3 = $(basename $(wildcard tests/3d/test-*.cxx))

# O EXEMPLO do Poisson 3D e o catalogo de solucoes manufaturadas: exemplo, e nao
# biblioteca -- nada em `src/` sabe que equacao se esta' resolvendo.
SRC3_EX = exemplos/poisson3d/poisson3d.cxx exemplos/poisson3d/manufaturada3d.cxx
OBJ3_X  = $(addprefix $(OBJ3)/,$(notdir $(SRC3_EX:.cxx=.o)))
EXEMPLOS3 = exemplos/poisson3d/poisson3d-mms

3d: $(TESTES3) $(EXEMPLOS3)

$(OBJ3):
	mkdir -p $(OBJ3)

$(OBJ3)/%.o: $(HIGTREE_SRC)/%.c | $(OBJ3)
	$(CXX) $(CXXFLAGS3) -x c++ -c $< -o $@

$(OBJ3)/%.o: $(HIGTREE_SRC)/%.cpp | $(OBJ3)
	$(CXX) $(CXXFLAGS3) -c $< -o $@

$(OBJ3)/%.o: src/comum/%.cxx | $(OBJ3)
	$(CXX) $(CXXFLAGS3) -c $< -o $@

$(OBJ3)/%.o: src/3d/%.cxx | $(OBJ3)
	$(CXX) $(CXXFLAGS3) -c $< -o $@

$(OBJ3)/%.o: exemplos/poisson3d/%.cxx | $(OBJ3)
	$(CXX) $(CXXFLAGS3) -c $< -o $@

tests/3d/%: tests/3d/%.cxx tests/testing.h $(OBJ3_T) $(OBJ3_X) $(OBJ3_N)
	$(CXX) $(CXXFLAGS3) $< $(OBJ3_T) $(OBJ3_X) $(OBJ3_N) -o $@ $(LDFLAGS) $(LIBS)

exemplos/poisson3d/%: exemplos/poisson3d/%.cxx $(OBJ3_T) $(OBJ3_X) $(OBJ3_N)
	$(CXX) $(CXXFLAGS3) -Itests $< $(OBJ3_T) $(OBJ3_X) $(OBJ3_N) -o $@ $(LDFLAGS) $(LIBS)

clean:
	$(RM) -rf obj
	$(RM) -f $(wildcard tests/2d/test-*[!x]) $(wildcard tests/3d/test-*[!x]) bench/escalabilidade bench/franja $(EXEMPLOS2) 
	@echo "limpo: objetos das duas dimensoes e binarios de teste"

.PHONY: 2d 3d bench clean all
all: 2d
