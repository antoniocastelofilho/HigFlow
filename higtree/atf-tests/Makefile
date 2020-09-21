# Must be set outwards
# Try variables path of either PETSc versions (3.5.2 or 3.5.3):
-include ${PETSC_DIR}/conf/variables
-include ${PETSC_DIR}/lib/petsc-conf/petscvariables
include ../Makefile.common

CC = gcc

OPTSOPT = -g
STANDARDS = -std=c99

CFLAGS = $(shell pkg-config --cflags atf-c) 
CFLAGS+= $(shell pkg-config --cflags glib-2.0)
CFLAGS+= $(PETSC_CC_INCLUDES)
CFLAGS+= -I$(SRCPATH) $(OPTSOPT) $(STANDARDS)

LDFLAGS = $(shell pkg-config --libs-only-L --libs-only-other atf-c) 
LDFLAGS+= $(shell pkg-config --libs-only-L --libs-only-other glib-2.0)
LDFLAGS+= -L$(LIBPATH)

LIBS = $(shell pkg-config --libs-only-l atf-c) 
LIBS+= $(shell pkg-config --libs-only-l glib-2.0)
LIBS+= $(PETSC_LIB)
LIBS+= $(PETSC_EXTRA_LIB)
LIBS+= -lm -lgcov
LIBS+= -lhdf5

LIB2D = -lhig2d
LIB3D = -lhig3d

TESTS2D+= test-basic-operation-2d
TESTS2D+= test-higtree-iterators-2d
TESTS2D+= test-higtree-io-2d
TESTS2D+= test-mapper-2d
TESTS2D+= test-wls-2d
TESTS2D+= test-point-cloud-2d
TESTS2D+= test-domain-2d
TESTS2D+= test-example-2d
TESTS2D+= test-solver-2d

TESTS3D+= test-basic-operation-3d 
TESTS3D+= test-higtree-io-3d


TESTS = $(TESTS2D) $(TESTS3D) 

all : run

rungcov : test2d test3d
	atf-run | atf-report

run : test2d test3d
	atf-run | atf-report

run-verbose : test2d test3d
	atf-run

%-2d : %-2d.c
ifndef RUN
	$(error Please execute make test from $(HIGPATH))
endif
	$(CC) $(CFLAGS) -DDIM=2 -c $< -o $@.o
	$(CC) $(LDFLAGS) $@.o -o $@ $(LIB2D) $(LIBS)

%-3d : %-3d.c
ifndef RUN
	$(error Please execute make test from $(HIGPATH))
endif
	$(CC) $(CFLAGS) -DDIM=3 -c $< -o $@.o
	$(CC) $(LDFLAGS) $@.o -o $@ $(LIB3D) $(LIBS)

test2d : $(TESTS2D)

test3d : $(TESTS3D)

clean :
	$(RM) -f *.o *.gcda *.gcno $(TESTS)

