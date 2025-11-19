# Must be set outwards
# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# HigTree common
include ${HIGTREE_DIR}/Makefile.common

CC = gcc

DIM=2

STANDARDS = -std=gnu99

CFILE = ns-example-2d

ifdef DEBUG
	CFLAGS+= -g -ggdb -Og -Wall
	GDB_FLAG = gdb --args
	DEBUG_DESC = "gdb debug"
else
	OPTSOPT = -flto=8 -fno-fat-lto-objects
	CFLAGS+= -march=native -mtune=native -O3 -Werror
	DEBUG_DESC = "No Debugging/Optimized"
endif

CFLAGS+= $(shell pkg-config --cflags glib-2.0 libfyaml)
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(STANDARDS) -DDIM=$(DIM) $(OPTSOPT)

LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5 libfyaml)
LDFLAGS+= -L$(LIBPATH) $(OPTSOPT)

LIBS+= -lhig$(DIM)d
LIBS+= $(shell pkg-config --libs-only-l glib-2.0 hdf5 libfyaml)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -ltrilinos_zoltan -lHYPRE

ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
	LIBS+= -lrt
endif

MAKEFLAGS += --no-print-directory


ifdef MESH
	MESHNAME = $(MESH)
else
	MESHNAME = $(shell cat mesh/__using/name.txt)
endif

#### type of flow
INPUT = $(IN)
OUTNAME = $(shell python3 setup.py $(INPUT))__$(MESHNAME)

ifeq ($(OUTNAME),)
	OUTNAME = noname
endif

ifneq ($(ENAME),)
	EXTRANAME := __$(ENAME)
	OUTNAME := $(OUTNAME)$(EXTRANAME)
endif


### number of processes
ifeq ($(NP),)
	NP = 1
endif

ifeq ($(THREAD),)
	THREAD_FLAG = -use-hwthread-cpus
endif
ifeq ($(THREAD),0)
	THREAD_FLAG =
endif
ifeq ($(THREAD),1)
	THREAD_FLAG = --oversubscribe
endif


ifdef RESTART
	LOAD_PATH = output/$(OUTNAME)/save/save
else
	LOAD_PATH = input/load
endif


OBJ1 += $(CFILE).o 
OBJ1 += ../src/hig-flow-kernel.o ../src/hig-flow-io.o ../src/hig-flow-eval.o ../src/hig-flow-discret.o
OBJ1 += ../src/hig-flow-bc.o ../src/hig-flow-ic.o ../src/hig-flow-terms.o ../src/hig-flow-step.o
OBJ1 += ../src/hig-flow-res.o ../src/hig-flow-timestep.o
OBJ1 += ../src/hig-flow-step-generalized-newtonian.o  ../src/hig-flow-step-viscoelastic.o
OBJ1 += ../src/hig-flow-step-viscoelastic-integral.o ../src/hig-flow-mittag-leffler.o
OBJ1 += ../src/hig-flow-step-electroosmotic-viscoelastic.o ../src/hig-flow-step-electroosmotic.o
OBJ1 += ../src/hig-flow-step-multiphase.o ../src/hig-flow-vof-adap-hf.o ../src/hig-flow-vof-plic.o
OBJ1 += ../src/hig-flow-vof-9-cells.o ../src/hig-flow-vof-finite-difference-normal-curvature.o 
OBJ1 += ../src/hig-flow-vof-elvira.o
OBJ1 += ../src/hig-flow-step-multiphase-viscoelastic.o 
OBJ1 += ../src/hig-flow-step-multiphase-electroosmotic.o ../src/hig-flow-step-multiphase-electroosmotic-viscoelastic.o
OBJ1 += utilities-2d.o

$(CFILE): $(OBJ1)
	@$(CC) $(LDFLAGS) $(OBJ1) -o $@ $(LIBS)

%.o: %.c ${HIGTREE_DIR}/lib/libhig$(DIM)d.a
	@$(CC) $(CFLAGS) -c $< -o $@

setup_mesh:
	@find mesh/__using/ -type f -name "*.amr" -delete; \
	cp -r mesh/$(MESH)/. mesh/__using/; \
	mv mesh/__using/domain/$(MESH)-d.amr mesh/__using/domain/ch-d.amr; \
	mv mesh/__using/bc/$(MESH)-bc-0.amr mesh/__using/bc/ch-bc-0.amr; \
	mv mesh/__using/bc/$(MESH)-bc-1.amr mesh/__using/bc/ch-bc-1.amr; \
	mv mesh/__using/bc/$(MESH)-bc-2.amr mesh/__using/bc/ch-bc-2.amr; \
	mv mesh/__using/bc/$(MESH)-bc-3.amr mesh/__using/bc/ch-bc-3.amr; \
	touch mesh/__using/name.txt; \
	sed -i '1s/.*/$(MESH)/' mesh/__using/name.txt;

run:
	@mkdir -p "output/$(OUTNAME)/vtk" || true
	@mkdir -p "output/$(OUTNAME)/save" || true
	@mkdir -p "output/$(OUTNAME)/res" || true


	@echo "Running $(CFILE) with $(NP) processes - ${DEBUG_DESC}"
	@echo "Inputs from input/"
ifdef MESH
	@make setup_mesh MESH=$(MESH);
	@echo "Mesh copied from mesh/$(MESH)/ to mesh/__using/"
endif
	@echo "Setting checkpoints in output/$(OUTNAME)/save/save"
	@echo "VTK output in output/$(OUTNAME)/vtk"
	
	@mpirun $(THREAD_FLAG) -np $(NP) $(GDB_FLAG) ./$(CFILE) \
		$(LOAD_PATH) \
		output/$(OUTNAME)/save/save \
		output/$(OUTNAME)/vtk/ns$(EXTRANAME).print \
		-ksp_type fbcgsr -pc_type hypre -pc_hypre_type boomeramg -ksp_atol 1e-10 -ksp_rtol 1e-14 -ksp_initial_guess_nonzero true \
		| tee output/$(OUTNAME)/save/ns$(EXTRANAME).log; 
#		python3 scripts/compress_vtk.py output/$(OUTNAME)/vtk; \

compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)

clean:
	@$(RM) $(CFILE).o $(CFILE) utilities-2d.o

build_run:
	@make clean; make;
	@make run;

premake:
	$(CC) -E $(CFLAGS) ns-example-2d.c -o ns-example-2d.i

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

clean_all:
	$(MAKE) -C .. clean
	$(MAKE) -C $(HIGTREE_DIR) clean
	$(MAKE) clean

.PHONY: run clean_o compilerflags

outname:
	@echo "$(OUTNAME)"
