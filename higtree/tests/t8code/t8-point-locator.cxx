// Localizacao por ponto no t8code.  Ver t8-point-locator.h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-point-locator.h"

#include <cmath>
#include <cstring>

#if DIM == 2
#define ECLASSE T8_ECLASS_QUAD
#else
#define ECLASSE T8_ECLASS_HEX
#endif

#define LADO 8                  // celulas por direcao, como o teste do MTree
#define MAX_TREES 4

// Uma floresta por numero de arvores de cmesh.  Construir a cada consulta seria
// dominar o tempo do teste com montagem, nao com localizacao.
static t8_forest_t g_cache[MAX_TREES + 1];
static int g_iniciado = 0;

// Percorrer as folhas ao contrario.  Existe SO' para o teste: a ordem natural da
// curva de preenchimento visita a celula de menor coordenada primeiro, entao um
// "desempate" que apenas guardasse o primeiro candidato passaria por SORTE DE
// ORDEM, nao por regra.  MEDIDO: com a regra removida, os tres casos continuavam
// verdes.  Invertendo o percurso, a regra aparece -- ou a falta dela.
static int g_inverte = 0;

static void
inicializa_uma_vez (void)
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
static t8_forest_t
floresta (int ntrees_x)
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

extern "C" int
t8_localiza_ponto (int ntrees_x, const Point p, Point centro)
{
  if (ntrees_x < 1 || ntrees_x > MAX_TREES) return 0;
  inicializa_uma_vez ();

  t8_forest_t f = floresta (ntrees_x);
  if (f == NULL) return 0;

  double ponto[3] = { 0.0, 0.0, 0.0 };
  for (int d = 0; d < DIM; d++) ponto[d] = p[d];

  // --------------------------------------------------------------------
  // Varre as folhas e junta TODAS as que contem o ponto.  Sao varias quando
  // ele cai numa face: e' a propria documentacao do t8code que avisa.
  // --------------------------------------------------------------------
  int achou = 0;
  double melhor[3] = { 0.0, 0.0, 0.0 };

  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  for (t8_locidx_t k = 0; k < ntrees; k++) {
    const t8_locidx_t it = g_inverte ? (ntrees - 1 - k) : k;
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t q = 0; q < ne; q++) {
      const t8_locidx_t ie = g_inverte ? (ne - 1 - q) : q;
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      int dentro = 0;
      t8_forest_element_points_inside (f, it, e, ponto, 1, &dentro, 1e-12);
      if (!dentro) continue;

      double c[3];
      t8_forest_element_centroid (f, it, e, c);

      // ---------------- A REGRA DE DESEMPATE (clausula C8) ----------------
      // A primitiva do t8code nao decide entre vizinhos quando o ponto esta'
      // na face comum.  O contrato decide: fica a celula de MENOR coordenada,
      // porque o intervalo e' fechado em cima e aberto embaixo.  Entre os
      // candidatos, portanto, o de menor centro -- comparado direcao a
      // direcao, da primeira para a ultima.
      if (!achou) {
        for (int d = 0; d < 3; d++) melhor[d] = c[d];
        achou = 1;
        continue;
      }
      for (int d = 0; d < DIM; d++) {
        if (c[d] < melhor[d] - 1e-12) {                 // este e' menor
          for (int k = 0; k < 3; k++) melhor[k] = c[k];
          break;
        }
        if (c[d] > melhor[d] + 1e-12) break;            // o guardado e' menor
      }                                                  // empate: segue direcao
    }
  }

  if (!achou) return 0;
  for (int d = 0; d < DIM; d++) centro[d] = melhor[d];
  return 1;
}

extern "C" void
t8_localizador_inverte_percurso (int inverte)
{
  g_inverte = inverte;
}

extern "C" void
t8_localizador_encerra (void)
{
  for (int i = 0; i <= MAX_TREES; i++) {
    if (g_cache[i] != NULL) { t8_forest_unref (&g_cache[i]); g_cache[i] = NULL; }
  }
}
