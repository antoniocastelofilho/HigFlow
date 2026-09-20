// Floresta de teste comum.  Ver t8-forest-cache.hxx.

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-forest-cache.hxx"

#include <cstring>

#if DIM == 2
#define ECLASSE T8_ECLASS_QUAD
#else
#define ECLASSE T8_ECLASS_HEX
#endif

#define LADO 8                  // celulas por direcao, como o teste do MTree

// Uma floresta por numero de arvores de cmesh.  Construir a cada consulta seria
// dominar o tempo do teste com montagem, nao com localizacao.
static t8_forest_t g_cache[T8_MAX_TREES + 1];
static int g_iniciado = 0;


void
t8_inicializa_uma_vez (void)
{
  if (g_iniciado) return;
  g_iniciado = 1;
  std::memset (g_cache, 0, sizeof g_cache);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
  t8_init (SC_LP_ERROR);
}

// [0,1]^DIM com LADO celulas por direcao, repartido em `ntrees_x` arvores de
// cmesh na direcao x.  Uma arvore -> nivel log2(LADO); duas -> nivel
// log2(LADO/2) em cada, e assim por diante.  A malha resultante e' a MESMA; o
// que muda e' a divisao, que e' exatamente o que a C9 poe a prova.
t8_forest_t
t8_floresta_de_teste (int ntrees_x)
{
  if (g_cache[ntrees_x] != NULL) return g_cache[ntrees_x];

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  if (ntrees_x == 1) {
    t8_cmesh_new_hypercube (&cmesh, ECLASSE, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else {
    const double caixa[24] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  1, 1, 0,
                               0, 0, 1,  1, 0, 1,  0, 1, 1,  1, 1, 1 };
    // Divide em TODAS as direcoes, nao so' em x.  Dividir apenas em x deixa cada
    // arvore com 0,5 de largura e 1,0 de altura, e o refino uniforme divide por
    // 2^nivel em cada direcao DA ARVORE -- o resultado sao celulas retangulares e
    // uma malha diferente da de uma arvore so'.  Medido: com ntrees_x=2, os
    // centros em y saiam em 0,125 e 0,375, isto e' 4 celulas em y em vez de 8.
    t8_cmesh_new_hypercube_pad (cmesh, ECLASSE, sc_MPI_COMM_WORLD, caixa,
                                ntrees_x, ntrees_x,
                                (DIM == 3 ? ntrees_x : 0), 0);
  }

  int nivel = 0;
  for (int k = LADO / ntrees_x; k > 1; k /= 2) nivel++;

  const t8_scheme_c *scheme = t8_scheme_new_default ();
  g_cache[ntrees_x] =
    t8_forest_new_uniform (cmesh, scheme, nivel, 0, sc_MPI_COMM_WORLD);
  return g_cache[ntrees_x];
}


void
t8_libera_florestas (void)
{
  for (int i = 0; i <= T8_MAX_TREES; i++) {
    if (g_cache[i] != NULL) { t8_forest_unref (&g_cache[i]); g_cache[i] = NULL; }
  }
}
