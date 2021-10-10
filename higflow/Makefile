# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ${HIGTREE_DIR}/Makefile.common
include Makefile.higflow.common

DEFINES=-DDIM=$(DIM)
TARGET=libhigflow$(DIM)d.a

OPTSOPT = -Ofast -march=native -mtune=native -flto=8 -fno-fat-lto-objects
#OPTSOPT = -g
STANDARDS = -std=gnu99

CFLAGS+= $(shell pkg-config --cflags glib-2.0) -Werror
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(OPTSOPT) $(STANDARDS) -DDIM=$(DIM)
CFLAGS+= -I$(INCPATH) 
CFLAGS+= -I$(HIGFLOW_SRCPATH) 

LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5 libfyaml)
LDFLAGS+= -L$(LIBPATH) $(OPTSOPT)

LIBS+= -lhig$(DIM)d
LIBS+= $(shell pkg-config --libs-only-l glib-2.0 hdf5 libfyaml)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -ltrinilos_zoltan -lHYPRE_krylov
#LIBS+= -lm -lzoltan -lHYPRE

CC = mpicc
AR = gcc-ar
ANLIB = gcc-ranlib

CXXSTANDARD = -std=gnu++11
CSTANDARD = -std=gnu99

ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
  LIBS+= -lrt
endif

MODULES = \
	kernel \
	io \
	ic \
	bc \
	eval \
	discret \
	terms \
	step \
	step-generalized-newtonian \
	step-multiphase \
	step-multiphase-viscoelastic \
	vof-plic \
	vof-adap-hf \
	vof-9-cells\
	vof-finite-difference-normal-curvature \
	step-viscoelastic \
	step-viscoelastic-integral \
	step-electroosmotic \

OBJS = $(MODULES:%=$(HIGFLOW_OBJPATH)/hig-flow-%.o)

.PHONY: all compilerflags clean

all: $(HIGFLOW_LIBPATH)/$(TARGET)
ifeq ($(wildcard $(HIGFLOW_INCPATH)),)
	mkdir $(HIGFLOW_INCPATH)
endif
	cp $(HIGFLOW_SRCPATH)/*.h $(HIGFLOW_INCPATH)

$(HIGFLOW_LIBPATH)/$(TARGET) : $(OBJS)
ifeq ($(wildcard $(HIGFLOW_LIBPATH)),)
	mkdir $(HIGFLOW_LIBPATH)
endif
	$(AR) -cr $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)
	mv $(TARGET) $(HIGFLOW_LIBPATH)

compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)
	@echo COMPILER_INCPATH=$(INCPATH)
	@echo COMPILER_SRCPATH=$(SRCPATH)
	@echo COMPILER_HIGFLOWSRCPATH=$(HIGFLOW_SRCPATH)
	@echo COMPILER_HIGTREEDIR=$(HIGTREE_DIR)

%.o : %.c
ifndef DIM
	$(error Please specify the dimension to be compiled with DIM=d)
endif
	$(CC) $(CSTANDARD) $(CFLAGS) -c $< -o $@

clean :
	$(RM) -f $(HIGFLOW_OBJPATH)/*.o $(HIGFLOW_OBJPATH)/*.gcda $(HIGFLOW_OBJPATH)/*.gcno $(HIGFLOW_OBJPATH)/*.gcov
	$(RM) -f $(HIGFLOW_LIBPATH)/*.a
