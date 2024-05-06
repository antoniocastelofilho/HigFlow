# Must be set outwards
# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# HigTree common
include ${HIGTREE_DIR}/Makefile.common

CC = gcc

DIM=2

ifeq ($(DEBUG),)
	OPTSOPT = -Ofast -march=native -mtune=native -flto=8 -fno-fat-lto-objects
	DEBUG_DESC = "No Debugging/Optimized"
endif
ifeq ($(DEBUG),1)
	OPTSOPT = -Og -g -ggdb
	DEBUG_DESC = "Simple Debug/gdb"
	GDB_FLAG = gdb --args
endif
STANDARDS = -std=gnu99

CFLAGS+= $(shell pkg-config --cflags glib-2.0) -Werror
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(OPTSOPT) $(STANDARDS) -DDIM=$(DIM)

LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5 libfyaml)
LDFLAGS+= -L$(LIBPATH) $(OPTSOPT)

LIBS+= -lhig$(DIM)d
LIBS+= $(shell pkg-config --libs-only-l glib-2.0 hdf5 libfyaml)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -ltrilinos_zoltan -lHYPRE_krylov

ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
	LIBS+= -lrt
endif

ifeq ($(NP),)
	NP = 1
endif

OBJ1 = ns-example-2d.o 
OBJ1 += ../src/hig-flow-kernel.o ../src/hig-flow-io.o ../src/hig-flow-eval.o ../src/hig-flow-discret.o
OBJ1 += ../src/hig-flow-bc.o ../src/hig-flow-ic.o ../src/hig-flow-terms.o ../src/hig-flow-step.o
OBJ1 += ../src/hig-flow-step-multiphase.o ../src/hig-flow-vof-adap-hf.o ../src/hig-flow-vof-plic.o
OBJ1 += ../src/hig-flow-vof-9-cells.o ../src/hig-flow-vof-finite-difference-normal-curvature.o 
OBJ1 += ../src/hig-flow-vof-elvira.o

ns-example: $(OBJ1)
	$(CC) $(LDFLAGS) $(OBJ1) -o $@ $(LIBS)

%.o: %.c ${HIGTREE_DIR}/lib/libhig$(DIM)d.a
	$(CC) $(CFLAGS) -c $< -o $@

run:
	mkdir -p VTKS || true
	mkdir -p DATA || true
	mkdir -p output || true
	mpirun -use-hwthread-cpus -n $(NP) $(GDB_FLAG) ./ns-example \
		input/example-3d.load \
		output/example-3d.save \
		VTKS/example-3d.print \
		-ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-10 -ksp_rtol 1e-10

compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)

clean :
		$(RM) -f ns-example.o ns-example

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
