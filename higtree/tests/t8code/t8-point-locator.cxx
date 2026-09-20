// Localizacao por ponto no t8code.  Ver t8-point-locator.h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-point-locator.h"
#include "t8-forest-cache.hxx"

#include <cmath>

// Percorrer as folhas ao contrario.  Existe SO' para o teste: a ordem natural da
// curva de preenchimento visita a celula de menor coordenada primeiro, entao um
// "desempate" que apenas guardasse o primeiro candidato passaria por SORTE DE
// ORDEM, nao por regra.  MEDIDO: com a regra removida, os tres casos continuavam
// verdes.  Invertendo o percurso, a regra aparece -- ou a falta dela.
static int g_inverte = 0;

extern "C" int
t8_localiza_ponto (int ntrees_x, const Point p, Point centro)
{
  if (ntrees_x < 1 || ntrees_x > T8_MAX_TREES) return 0;
  t8_inicializa_uma_vez ();

  t8_forest_t f = t8_floresta_de_teste (ntrees_x);
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
  t8_libera_florestas ();
}
