# Must be set outwards
# Try variables path of either PETSc versions (3.5.2 or 3.5.3):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ../Makefile.common

CC = mpicc

#OPTSOPT += -g

ifdef DEBUG
OPTSOPT += -DDEBUG -g
else
OPTSOPT += -Ofast -flto=4 -fno-fat-lto-objects -march=native -mtune=native -DNDEBUG
endif

LDFLAGS += $(OPTSOPT)

ifdef PROFILE
CFLAGS += -fprofile-arcs -ftest-coverage
endif

STANDARDS = -std=gnu99

CFLAGS+= $(shell pkg-config --cflags glib-2.0)
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(OPTSOPT) $(STANDARDS) #-fopenmp
CFLAGS+=-I/usr/include/trilinos/

LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5)
LDFLAGS+= -L$(LIBPATH)
LDFLAGS+= -L$(ZOLTAN_BUILD_DIR)src/

LIBS+= -lhig$(DIM)d
LIBS+= $(shell pkg-config --libs-only-l glib-2.0 hdf5)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -lhdf5 -ltrilinos_zoltan -lHYPRE_krylov

ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
  LIBS+= -lrt
endif

ifdef VIENNACL
LIBS += -lOpenCL -lstdc++
endif

ifdef SOR
LIBS += -lnuma -fopenmp
endif

DIM_EXFILES  = $(wildcard *-example-$(DIM)d.c)
DIMLESS_EXFILES = $(wildcard *-example.c)
EXAMPLES = $(DIM_EXFILES:%.c=%) $(DIMLESS_EXFILES:%.c=%-$(DIM)d)

.PHONY: all ex$(DIM)d clean

all :
	make DIM=2 ex2d
	make DIM=3 ex3d

ex$(DIM)d : $(EXAMPLES)

%-example-$(DIM)d : %-example-$(DIM)d.c ../lib/libhig$(DIM)d.a
ifndef DIM
	$(error DIM no set)
endif
	$(CC) $(CFLAGS) -DDIM=$(DIM) -c $< -o $@.o
	+$(CC) $(LDFLAGS) $@.o -o $@ $(LIBS)

# Support for DIM-agnostic examples.
%-example-$(DIM)d : %-example.c ../lib/libhig$(DIM)d.a
ifndef DIM
	$(error DIM no set)
endif
	$(CC) $(CFLAGS) -DDIM=$(DIM) -c $< -o $@.o
	+$(CC) $(LDFLAGS) $@.o -o $@ $(LIBS)


test : test-fast test-slow

test-fast :
	./test-termination-detection-example-2d
	./calc-bc-example-2d
	./dim-iterator-example-2d 1
	./face-iterator-example-2d
	./read-amr-calc-partition-example-2d mesh00000001-2d.amr
	./read-amr-interpolate-example-2d mesh00000001-2d.amr
	./read-amr-map-ids-print-adj-matrix-example-2d mesh2d-2.amr > /dev/null
	mpirun -n 2 ./read-amr-send-to-other-node-example-2d mesh2d-1.amr 1
	./read-amr-solve-derivatives-example-2d mesh2d-2.amr
	./read-amr-solve-poisson-force-example-2d mesh2d-2.amr
	./read-amr-solver-domain-example-2d mesh00000001-2d.amr
	./read-amr-traverse-all-leaves-print-all-cell-neighbours-example-2d mesh2d-2.amr
	./read-amr-traverse-all-leaves-print-center-example-2d mesh2d-2.amr
	./read-amr-write-vtk-example-2d mesh2d-2.amr /dev/null
	./read-amr-write-vtk-with-properties-example-2d mesh2d-2.amr /dev/null
	./read-amr-iterate-facets-example-2d mesh2d-2.amr
	./interpolate-example-3d mesh-test-3d.amr 2 0.01
	./read-amr-calc-partition-example-3d mesh-test-3d.amr
	./read-amr-calc-point-cloud-example-3d mesh-test-3d.amr
	./read-amr-distribute-example-3d mesh-test-3d.amr
	./read-amr-interpolate-add-example-3d mesh-test-3d.amr 1
	./read-amr-interpolate-example-3d mesh-test-3d.amr
	./read-amr-write-vtk-example-3d mesh-test-3d.amr /dev/null
	./test-solver-simple-example-2d
	./read-amr-write-xdmf-example-3d mesh-test-3d.amr /tmp/ttt

test-slow :
	mpirun -n 2 ./read-amr-solve-derivatives-example-2d mesh00000001-2d.amr
	mpirun -n 2 ./read-amr-solve-channel-ns-example-2d mesh2d-2.amr 0.01 10.0 3 1
	./poisson-variable-rho-example-2d
	mpirun -n 4 ./poisson-parallel-example-2d
	./poisson-example-2d
	./poisson-composed-higtree-example-2d
	./read-amr-solve-cavity-ns-incremental-example-3d mesh-test-3d.amr 0.01 10.0 3 1
	./read-amr-solve-cavity-ns-example-3d mesh-test-3d.amr 0.01 10.0 3 0
	mpirun -n 4 ./distributed-properties-example-2d mesh00000001-2d.amr
	mpirun -n 4 ./read-amr-solve-cavity-ns-example-2d mesh2d-2.amr 0.01 10.0 3 1 VTK
	./read-amr-solve-cavity-ns-example-2d mesh2d-2.amr 0.01 10.0 3 1 VTK
#	./read-amr-solve-ns-restart-example-2d \
		0 \
		1 /tmp/v.txt \
		2 \
			amrs/bfs-test/bfs-test-1-2d.amr \
			amrs/bfs-test/bfs-test-2-2d.amr \
		6 \
			amrs/bfs-test/bfs-test-bc-1-2d.amr \
			1 0 0 2 0 0   \
			amrs/bfs-test/bfs-test-bc-2-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-3-2d.amr \
			0 0 1 0 1 0   \
			amrs/bfs-test/bfs-test-bc-4-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-5-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-6-2d.amr \
			1 0 0 0 0 0   \
		0.0001 10.0 2 1  \
		4 \
			amrs/bfs-test/bfs-test-probe-3-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-1-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-2-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-4-2d.amr 0 \
		1.0 \
		/tmp \
		-ksp_type bcgs -pc_type none
	./read-amr-solve-ns-load-control-example-2d \
		2 \
			amrs/bfs-test/bfs-test-1-2d.amr \
			amrs/bfs-test/bfs-test-2-2d.amr \
		6 \
			amrs/bfs-test/bfs-test-bc-1-2d.amr \
			1 0 0 2 0 0   \
			amrs/bfs-test/bfs-test-bc-2-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-3-2d.amr \
			0 0 1 0 1 0   \
			amrs/bfs-test/bfs-test-bc-4-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-5-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-6-2d.amr \
			1 0 0 0 0 0   \
		0.0001 10.0 2 1  \
		4 \
			amrs/bfs-test/bfs-test-probe-3-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-1-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-2-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-4-2d.amr 0 \
		1.0 \
		/tmp \
		-ksp_type bcgs -pc_type none
	./read-amr-solve-ns-example-2d \
		2 \
			amrs/bfs-test/bfs-test-1-2d.amr \
			amrs/bfs-test/bfs-test-2-2d.amr \
		6 \
			amrs/bfs-test/bfs-test-bc-1-2d.amr \
			1 0 0 2 0 0   \
			amrs/bfs-test/bfs-test-bc-2-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-3-2d.amr \
			0 0 1 0 1 0   \
			amrs/bfs-test/bfs-test-bc-4-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-5-2d.amr \
			1 0 0 0 0 0   \
			amrs/bfs-test/bfs-test-bc-6-2d.amr \
			1 0 0 0 0 0   \
		0.0001 10.0 2 1  \
		4 \
			amrs/bfs-test/bfs-test-probe-3-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-1-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-2-2d.amr 0 \
			amrs/bfs-test/bfs-test-probe-4-2d.amr 0 \
		1.0 \
		/tmp VTK \
		-ksp_type bcgs -pc_type none
	./read-amr-solve-poisson-on-facets-example-3d mesh-test-3d.amr
	./read-amr-solve-poisson-on-facets-example-2d mesh00000001-2d.amr
	mpirun -n 4 ./read-amr-solve-channel-ns-example-2d mesh00000001-2d.amr 0.01 10.0 2 1
	./read-amr-solve-channel-ns-example-2d mesh00000001-2d.amr 0.01 10.0 2 1
	./read-amr-solve-poisson-parallel-example-2d mesh00000001-2d.amr
	mpirun -n 4 ./read-amr-solve-poisson-parallel-example-2d mesh00000001-2d.amr
	./read-amr-solve-poisson-utilde-example-2d mesh00000001-2d.amr
	./read-amr-solve-poisson-parallel-example-3d mesh-test-3d.amr


clean :
	$(RM) -f *.o *-example-?d
