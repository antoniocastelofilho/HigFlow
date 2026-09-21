# Must be set outwards
# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# HigTree common
include ${HIGTREE_DIR}/Makefile.common

CC = mpic++

DIM=2

STANDARDS = -x c++ -std=gnu++17

CFILE = ns-example-2d

ifdef DEBUG
	CFLAGS+= -g -ggdb -Og -Wall
	GDB_FLAG = gdb --args
	DEBUG_DESC = "gdb debug"
else
	OPTSOPT = -flto=8 -fno-fat-lto-objects
	CFLAGS+= -march=native -mtune=native -Ofast -Werror
	DEBUG_DESC = "No Debugging/Optimized"
endif

CFLAGS+= $(shell pkg-config --cflags glib-2.0 libfyaml hdf5)
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(OPTSOPT) $(STANDARDS) -DDIM=$(DIM)

LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5 libfyaml)
LDFLAGS+= -L$(LIBPATH) $(OPTSOPT)

LIBS+= -lhig$(DIM)d
LIBS+= $(shell pkg-config --libs-only-l glib-2.0 hdf5 libfyaml)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -ltrilinos_zoltan -lHYPRE

# t8code, OPCIONAL.  Passe T8CODE=<prefixo> para que o exemplo possa escolher a
# malha produzida pelo t8code (em tempo de execucao, com HIGFLOW_MALHA=t8code).
# Sem isso nada muda e o exemplo nao depende da biblioteca.
#
# A unidade do produtor compila em C++20 porque o t8code exige, enquanto o resto
# fica em gnu++17 -- so' a assinatura em extern "C" atravessa.
ifdef T8CODE
CFLAGS  += -DHIGFLOW_COM_T8CODE -I$(T8CODE)/include
LDFLAGS += -L$(T8CODE)/lib -Wl,-rpath,$(T8CODE)/lib
LIBS    += -lt8 -lp4est -lsc
T8OBJ    = malha-t8.o $(HIGTREE_DIR)/tests/t8code/t8-mesh-rank.o $(HIGTREE_DIR)/tests/t8code/t8-forest-cache.o $(HIGTREE_DIR)/tests/t8code/t8-particao-grafo.o
T8FLAGS  = -x c++ -std=gnu++20 -O1 -DDIM=$(DIM) -I$(SRCPATH) -I../src -I$(HIGTREE_DIR)/tests/t8code \
           -DT8_ENABLE_MPI=1 -DT8_ENABLE_MPIIO=1 -I$(T8CODE)/include \
           $(shell pkg-config --cflags glib-2.0 hdf5) $(PETSC_CC_INCLUDES)
endif


ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
	LIBS+= -lrt
endif

ifeq ($(NP),)
	NP = 1
endif


ifdef RESTART
	LOAD_PATH = output/example-2d.save
else
	LOAD_PATH = input/example-2d.load
endif


OBJ1 = ns-example-2d.o 
OBJ1 += ../src/hig-flow-kernel.o ../src/hig-flow-io.o ../src/hig-flow-eval.o ../src/hig-flow-discret.o
OBJ1 += ../src/hig-flow-bc.o ../src/hig-flow-ic.o ../src/hig-flow-terms.o ../src/hig-flow-step.o
OBJ1 += ../src/hig-flow-step-multiphase.o ../src/hig-flow-vof-adap-hf.o ../src/hig-flow-vof-plic.o
OBJ1 += ../src/hig-flow-vof-9-cells.o ../src/hig-flow-vof-finite-difference-normal-curvature.o 
OBJ1 += ../src/hig-flow-vof-elvira.o
OBJ1 += ../src/hig-flow-step-generalized-newtonian.o ../src/hig-flow-step-viscoelastic.o ../src/hig-flow-mittag-leffler.o   
OBJ1 += ../src/hig-flow-step-multiphase-viscoelastic.o
OBJ1 += ../src/hig-flow-res.o
OBJ1 += ../src/hig-flow-linear-algebra.o
OBJ1 += ../src/hig-flow-viscoelastic-kernel.o
OBJ1 += $(T8OBJ)

ns-example: $(OBJ1)
	$(CC) $(LDFLAGS) $(OBJ1) -o $@ $(LIBS)

ifdef T8CODE
# A fonte vive em ../examples-common: ela le' a especificacao do proprio arquivo
# AMR, entao serve a qualquer exemplo de malha uniforme -- nao ha' nada decorado
# nela que a prenda a este.
malha-t8.o: ../examples-common/malha-t8.cxx
	$(CC) $(T8FLAGS) -c $< -o $@

# O adaptador do t8code vive sob higtree/tests porque e' opcional e e' la' que
# ele e' verificado, clausula por clausula.  O exemplo o compila de la' em vez de
# duplicar a decomposicao em caixas -- duplicar um algoritmo sutil e' pior que a
# dependencia estranha.  Quando o adaptador sair de `tests`, esta regra muda de
# caminho e mais nada.
$(HIGTREE_DIR)/tests/t8code/%.o: $(HIGTREE_DIR)/tests/t8code/%.cxx
	$(CC) $(T8FLAGS) -c $< -o $@
endif

%.o: %.c ${HIGTREE_DIR}/lib/libhig$(DIM)d.a
	$(CC) $(CFLAGS) -c $< -o $@

run:
	mkdir -p VTKS || true
	mkdir -p DATA || true
	mkdir -p output || true
	mpirun -use-hwthread-cpus -n $(NP) $(GDB_FLAG) ./ns-example \
		$(LOAD_PATH) \
		output/vof_oldroyd.save \
		VTKS/vof_oldroyd.print \
		-ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10

compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)

clean :
		$(RM) -f ns-example-2d.o ns-example

rebuild_higflow:
	$(MAKE) -C .. clean
	$(MAKE) -C .. DIM=$(DIM) DEBUG=$(DEBUG) RES=$(RES)

rebuild_higtree:
	$(MAKE) -C $(HIGTREE_DIR) clean
	$(MAKE) -C $(HIGTREE_DIR) DIM=$(DIM) DEBUG=$(DEBUG)

rebuild_hig:
	$(MAKE) rebuild_higtree DEBUG=$(DEBUG)
	$(MAKE) rebuild_higflow DEBUG=$(DEBUG) RES=$(RES)
rebuild_all:
	$(MAKE) rebuild_hig DEBUG=$(DEBUG) RES=$(RES)
	$(MAKE) clean
	$(MAKE)

.PHONY: run clean compilerflags
