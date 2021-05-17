# Must be set outwards
# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# HigTree common
include ${HIGTREE_DIR}/Makefile.common

CC = gcc

DIM=2

OPTSOPT = -Ofast -march=native -mtune=native -flto=8 -fno-fat-lto-objects
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

OBJ1 = ns-example-2d.o ../src/hig-flow-io.o ../src/hig-flow-kernel.o ../src/hig-flow-step.o ../src/hig-flow-eval.o ../src/hig-flow-discret.o ../src/hig-flow-terms.o ../src/hig-flow-ic.o ../src/hig-flow-bc.o ../src/hig-flow-step-generalized-newtonian.o ../src/hig-flow-step-multiphase.o ../src/hig-flow-step-viscoelastic.o

OBJ2 = ns-example-2d.o ../src/hig-flow-io.o ../src/hig-flow-kernel.o ../src/hig-flow-step.o ../src/hig-flow-eval.o ../src/hig-flow-discret.o ../src/hig-flow-terms.o ../src/hig-flow-ic.o ../src/hig-flow-bc.o ../src/hig-flow-step-generalized-newtonian.o ../src/hig-flow-step-multiphase.o ../src/hig-flow-step-viscoelastic.o

ns-example: $(OBJ1)
	$(CC) $(LDFLAGS) $(OBJ1) -o $@ $(LIBS)

ns-ve: $(OBJ2)
	$(CC) $(LDFLAGS) $(OBJ2) -o $@ $(LIBS)

%.o: %.c ${HIGTREE_DIR}/lib/libhig$(DIM)d.a
	$(CC) $(CFLAGS) -c $< -o $@

run:
	mkdir -p VTKS || true
	mkdir -p DATA || true
	mkdir -p output || true
	mpirun -n 1 ./ns-example \
		input/example-3d.load \
		output/example-3d.save \
		VTKS/example-3d.print \
		-ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-6 -ksp_rtol 1e-6

run-ve:
	mkdir -p VTKS || true
	mkdir -p DATA || true
	mkdir -p output || true
	mpirun -n 1 ./ns-ve \
		input/example-3d.load \
		output/example-3d.save \
		VTKS/example-3d.print \
		-ksp_type bcgs -pc_type bjacobi -ksp_atol 1e-6 -ksp_rtol 1e-6

compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)

clean :
		$(RM) -f *.o VTKS/* DATA/* output/* ns-example

.PHONY: run clean compilerflags
