# Arreio do marching cubes, standalone -- o nucleo nao depende do HiGFlow,
# nem de PETSc, nem de MPI, nem da dimensao: a grade e' dada por argumento.
# Por isso este teste roda contra campos ANALITICOS, com oraculos exatos de
# fechamento e genero, antes de o modulo ver dado de solver.
CC     = gcc
CFLAGS = -O2 -Wall -Wextra -I../../src
MOD    = ../../src/hig-flow-mc.c
DEP    = $(MOD) ../../src/hig-flow-mc.h

all: teste-mc

teste-mc: teste-mc.c $(DEP)
	$(CC) $(CFLAGS) teste-mc.c $(MOD) -lm -o $@

check: teste-mc
	@./teste-mc

clean:
	rm -f teste-mc

.PHONY: all check clean
