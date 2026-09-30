# Arreios do front-tracking, standalone -- compilam o modulo direto, sem a
# lib/PETSc/MPI (as fases pre-acoplamento nao tocam o solver).  DIM=2 (a frente e'
# polilinha; o 3D e' a fase B5).
#
#   teste-b1         adveccao cinematica + cirurgia (Rider-Kothe reversivel)
#   teste-curvatura  curvatura contra oraculos analiticos (circulo, elipse)
#   teste-tensao     forca de tensao superficial espalhada (conservacao, equilibrio)
#   teste-indicadora fracao de area por celula (particao, celula cheia/vazia)
#   teste-persistencia  gravacao/leitura da frente (retomada de corrida)
#   teste-esfera3d   FASE B5: geometria e forca da superficie triangulada (DIM=3)
#
# Cada arreio devolve codigo de saida (0 passou, 1 falhou) -- sao oraculos, nao
# so' impressoes.  `make check` roda todos e falha se algum falhar.
CC      = gcc
CFLAGS  = -O2 -Wall -DDIM=2 -I../../../higtree/src -I../../src
MOD     = ../../src/hig-flow-front-tracking.c
DEP     = $(MOD) ../../src/hig-flow-front-tracking.h

ARREIOS = teste-b1 teste-curvatura teste-tensao teste-indicadora teste-persistencia \
          teste-esfera3d

# O 3D compila com DIM=3 e contra OUTRO modulo: a frente 2D e' polilinha, a 3D e'
# superficie triangulada, e os dois nucleos se recusam a compilar na dimensao do
# outro (#error nos dois .c).  Mesmo padrao dos testes da higtree, que ja'
# constroem -2d e -3d do mesmo fonte.
CFLAGS3 = -O2 -Wall -DDIM=3 -I../../../higtree/src -I../../src
MOD3    = ../../src/hig-flow-front-tracking-3d.c
DEP3    = $(MOD3) ../../src/hig-flow-front-tracking-3d.h

all: $(ARREIOS)

teste-b1: teste-b1.c $(DEP)
	$(CC) $(CFLAGS) teste-b1.c $(MOD) -lm -o $@

teste-esfera3d: teste-esfera3d.c $(DEP3)
	$(CC) $(CFLAGS3) teste-esfera3d.c $(MOD3) -lm -o $@

teste-curvatura: teste-curvatura.c $(DEP)
	$(CC) $(CFLAGS) teste-curvatura.c $(MOD) -lm -o $@

teste-tensao: teste-tensao.c $(DEP)
	$(CC) $(CFLAGS) teste-tensao.c $(MOD) -lm -o $@

teste-indicadora: teste-indicadora.c $(DEP)
	$(CC) $(CFLAGS) teste-indicadora.c $(MOD) -lm -o $@

teste-persistencia: teste-persistencia.c $(DEP)
	$(CC) $(CFLAGS) teste-persistencia.c $(MOD) -lm -o $@

run: teste-b1
	./teste-b1

check: $(ARREIOS)
	@ok=1; for t in $(ARREIOS); do \
	  echo "=== $$t ==="; ./$$t || ok=0; echo; \
	done; \
	if [ $$ok -eq 1 ]; then echo "TODOS OS ARREIOS PASSARAM"; else echo "ALGUM ARREIO FALHOU"; exit 1; fi

clean:
	rm -f $(ARREIOS)
