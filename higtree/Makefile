# Try variables path of either PETSc versions (3.5.2, 3.5.3, 3.6.2):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
-include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include Makefile.common

TARGET=libhig$(DIM)d.a
DEFINES=-DDIM=$(DIM)

CC = mpicc
CXX = mpic++
AR = gcc-ar
RANLIB = gcc-ranlib

ifdef DEBUG
OPTSOPT += -DDEBUG
else
OPTSOPT += -Ofast -flto -fno-fat-lto-objects -march=native -mtune=native -DNDEBUG
endif

OPTSOPT += -g

CXXSTANDARD = -std=gnu++14
CSTANDARD = -std=gnu99

CFLAGS+= -Wall -Wno-unused-result -Wno-unused-variable
CFLAGS+=$(OPTSOPT) $(DEFINES)
CFLAGS+=$(PETSC_CC_INCLUDES)
CFLAGS+=$(shell pkg-config --cflags glib-2.0 hdf5)
CFLAGS+=-I/usr/include/trilinos
CFLAGS+=-I/usr/include/hypre
CFLAGS+=-I/usr/include/petsc

LDFLAGS =$(shell pkg-config --libs-only-L --libs-only-other glib-2.0 hdf5)

LIBS+=$(PETSC_LIB)
LIBS+=$(PETSC_EXTRA_LIB)
LIBS+=$(shell pkg-config --libs-only-l glib-2.0 hdf5)
LIBS+=-lm
LIBS+=-lhdf5

ARCH = $(shell uname)
ifeq ($(ARCH), Linux)
  LIBS+= -lrt
endif

# Name of the files to be compiled
MODULES = \
	allocator \
	build-fringe \
	higtree \
	uniqueid \
	higtree-iterator \
	higtree-io \
	higtree-serialize \
	isend-pool \
	coord \
	rect \
	mapper \
	wls \
	point-cloud \
	solver \
	solver-centralized \
	solver-hypre \
	solver-petsc \
	solver-debug-write \
	domain \
	pdomain \
	point-mapper \
	rng \
	higtree-parallel \
	lbal \
	utils \
	term-det \
	mapper-syncer \
	build-fringe \
	global-stencil

ifdef VIENNACL
MODULES+= solver-viennacl
CFLAGS+= -DUSE_VIENNACL=1
endif

ifdef SOR
MODULES+= solver-sor
CFLAGS+= -DUSE_SOR=1 -fopenmp
endif

# Prepend $(OBJPATH) to all objects
OBJS = $(MODULES:%=$(OBJPATH)/%-$(DIM)d.o)

HEADERS = $(notdir $(wildcard src/*.h))

.PHONY: compilerflags clean all

all: $(addprefix $(INCPATH)/, $(HEADERS)) $(LIBPATH)/libhig$(DIM)d.a

$(INCPATH)/%.h: src/%.h | $(INCPATH)
	cp src/$(notdir $@) $@

$(INCPATH):
	mkdir $(INCPATH)

$(LIBPATH)/libhig$(DIM)d.a: $(OBJS) | $(LIBPATH)
	$(AR) -cr $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)
	mv $(TARGET) $(LIBPATH)

$(LIBPATH):
	mkdir $(LIBPATH)

.PHONY: compilerflags clean
compilerflags:
	@echo COMPILER_FLAGS=$(CFLAGS)

%-$(DIM)d.o : %.c
ifndef DIM
	$(error Please specify the dimension to be compiled with DIM=d)
endif
	$(CC) $(CSTANDARD) $(CFLAGS) -c $< -o $@

%-$(DIM)d.o : %.cpp
ifndef DIM
	$(error Please specify the dimension to be compiled with DIM=d)
endif
	$(CXX) $(CXXSTANDARD) $(CFLAGS) -c $< -o $@


test2d :
	make DIM=2

test3d :
	make DIM=3

test-examples : examples
	cd examples ; make test

test : test2d test3d test-examples
	cd $(ATFPATH) ; make RUN=1

test-all : test-examples

gcov : test2d
	cd $(ATFPATH) ; make RUN=1 rungcov
	cd $(SRCPATH) ; gcov higtree-2d
	cd $(SRCPATH) ; gcov higtree-io-2d
	cd $(SRCPATH) ; gcov higtree-iterator-2d
	cd $(SRCPATH) ; gcov mapper-2d
	cd $(SRCPATH) ; gcov wls-2d
	cd $(SRCPATH) ; gcov solver-2d
	cd $(SRCPATH) ; gcov domain-2d

examples2d :
	make DIM=2
	cd examples ; make DIM=2 RUN=1

examples3d :
	make DIM=3
	cd examples ; make DIM=3 RUN=1

examples : examples2d examples3d

clean :
	$(RM) -f $(OBJPATH)/*.o $(OBJPATH)/*.gcda $(OBJPATH)/*.gcno $(OBJPATH)/*.gcov
	$(RM) -f $(LIBPATH)/*.a
	$(RM) $(INCPATH)/*.h

clean-objs :
	$(RM) -f $(OBJPATH)/*.o

dist-clean : clean
	cd $(ATFPATH) ; make clean
	cd $(EXMPATH) ; make clean
	rmdir $(LIBPATH)
	rmdir $(INCPATH)
	cd $(DOXPATH) ; make clean
