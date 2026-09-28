# Arreios do front-tracking, standalone -- compilam o modulo direto, sem a
# lib/PETSc/MPI (as fases pre-acoplamento nao tocam o solver).  DIM=2 (a frente e'
# polilinha; o 3D e' a fase B5).
#
#   teste-b1         adveccao cinematica + cirurgia (Rider-Kothe reversivel)
#   teste-curvatura  curvatura contra oraculos analiticos (circulo, elipse)
#   teste-tensao     forca de tensao superficial espalhada (conservacao, equilibrio)
#   teste-indicadora fracao de area por celula (particao, celula cheia/vazia)
#
# Cada arreio devolve codigo de saida (0 passou, 1 falhou) -- sao oraculos, nao
# so' impressoes.  `make check` roda todos e falha se algum falhar.
CC      = gcc
CFLAGS  = -O2 -Wall -DDIM=2 -I../../../higtree/src -I../../src
MOD     = ../../src/hig-flow-front-tracking.c
DEP     = $(MOD) ../../src/hig-flow-front-tracking.h

ARREIOS = teste-b1 teste-curvatura teste-tensao teste-indicadora

all: $(ARREIOS)

teste-b1: teste-b1.c $(DEP)
	$(CC) $(CFLAGS) teste-b1.c $(MOD) -lm -o $@

teste-curvatura: teste-curvatura.c $(DEP)
	$(CC) $(CFLAGS) teste-curvatura.c $(MOD) -lm -o $@

teste-tensao: teste-tensao.c $(DEP)
	$(CC) $(CFLAGS) teste-tensao.c $(MOD) -lm -o $@

teste-indicadora: teste-indicadora.c $(DEP)
	$(CC) $(CFLAGS) teste-indicadora.c $(MOD) -lm -o $@

run: teste-b1
	./teste-b1

check: $(ARREIOS)
	@ok=1; for t in $(ARREIOS); do \
	  echo "=== $$t ==="; ./$$t || ok=0; echo; \
	done; \
	if [ $$ok -eq 1 ]; then echo "TODOS OS ARREIOS PASSARAM"; else echo "ALGUM ARREIO FALHOU"; exit 1; fi

clean:
	rm -f $(ARREIOS)
