# Arreio B1 do front-tracking: adveccao cinematica + cirurgia, campo prescrito.
# Standalone -- compila o modulo direto, sem a lib/PETSc/MPI (B1 nao acopla o
# solver).  DIM=2 (a frente do B1 e' polilinha; o 3D e' a fase B5).
CC      = gcc
CFLAGS  = -O2 -Wall -DDIM=2 -I../../../higtree/src -I../../src
MOD     = ../../src/hig-flow-front-tracking.c
teste-b1: teste-b1.c $(MOD) ../../src/hig-flow-front-tracking.h
	$(CC) $(CFLAGS) teste-b1.c $(MOD) -lm -o $@
run: teste-b1
	./teste-b1
clean:
	rm -f teste-b1
