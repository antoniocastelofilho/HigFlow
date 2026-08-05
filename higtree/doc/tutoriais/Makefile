# Must be set outwards
# Try variables path of either PETSc versions (3.5.2 or 3.5.3):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ${HIGTREE_DIR}/Makefile.common

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
LIBS+= -lm -lhdf5 -lzoltan -lHYPRE

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

#ifdef NN_WEIGHTS
# CFLAGS += -DNN_WEIGHTS
# CFLAGS+= -I$(LIBTORCH_DIR)/include -I$(LIBTORCH_DIR)/include/torch/csrc/api/include
# LIBS+=-L$(LIBTORCH_DIR)/lib -ltorch -ltorch_cpu -lc10 #Torch libs
#endif

TUTFILES  = $(wildcard tutorial*.c)

EXEC:= $(filter-out %.c, $(wildcard tutorial*))


TUTORIALS = $(TUTFILES:%.c=%)

TUTORIALS2D = $(filter %-2d, $(TUTORIALS))
TUTORIALS3D = $(filter %-3d, $(TUTORIALS))

.PHONY: all ex$(DIM)d clean

all :
	make DIM=2 tut2d
	make DIM=3 tut3d

tut2d : $(TUTORIALS2D)
tut3d : $(TUTORIALS3D)

tutorial%-$(DIM)d : tutorial%-$(DIM)d.c ${HIGTREE_DIR}/lib/libhig$(DIM)d.a
ifndef DIM
	$(error DIM no set)
endif
	$(CC) $(CFLAGS) -DDIM=$(DIM) -c $< -o $@.o
	+$(CC) $(LDFLAGS) $@.o -o $@ $(LIBS)

clean :
	$(RM) -f *.o $(EXEC)
