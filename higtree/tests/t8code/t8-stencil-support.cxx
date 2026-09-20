// Coleta do suporte de estencil na floresta do t8code.  Ver o .h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-stencil-support.h"
#include "t8-forest-cache.hxx"

#include <cmath>
#include <cstdlib>
#include <vector>

// Um elemento pelo seu indice na FLORESTA LOCAL -- que e' como
// `t8_forest_leaf_face_neighbors` devolve os vizinhos.  Indexar por (arvore,
// posicao na arvore) obrigaria a converter a cada vizinho, e converter errado e'
// silencioso: daria um elemento existente, so' que outro.
typedef t8_locidx_t Folha;

// Acha a folha que contem `x`.  Varredura simples: e' o ponto de PARTIDA da busca
// em largura, feito uma vez por consulta, e nao o mecanismo que se quer exercitar.
static bool
folha_do_ponto (t8_forest_t f, const double ponto[3], Folha *saida)
{
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_locidx_t base = t8_forest_get_tree_element_offset (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      int dentro = 0;
      t8_forest_element_points_inside (f, it, e, ponto, 1, &dentro, 1e-12);
      if (dentro) { *saida = base + ie; return true; }
    }
  }
  return false;
}

extern "C" int
t8_monta_suporte (int ntrees_x, const Point x, int minpts, int maxpts, Point pts[])
{
  t8_inicializa_uma_vez ();
  t8_forest_t f = t8_floresta_de_teste (ntrees_x);
  if (f == NULL) return 0;

  double ponto[3] = { 0.0, 0.0, 0.0 };
  for (int d = 0; d < DIM; d++) ponto[d] = x[d];

  Folha inicio;
  if (!folha_do_ponto (f, ponto, &inicio)) return 0;

  const t8_scheme_c *scheme = t8_forest_get_scheme (f);

  // ------------------------------------------------------------------------
  // Busca em largura por VIZINHANCA DE FACE.  E' a primitiva topologica real do
  // t8code -- a mesma que a C11 mediu atravessar salto 4:1 -- e nao um filtro de
  // caixa sobre todas as folhas, que serviria igual para qualquer backend e nao
  // exercitaria malha nenhuma.
  // ------------------------------------------------------------------------
  std::vector<Folha> fila;
  std::vector<Folha> vistos;
  fila.push_back (inicio);
  vistos.push_back (inicio);

  // Percorre em ANEIS e para ao FECHAR o anel em que juntou pontos bastante.
  // Parar no meio de um anel deixaria o suporte torto de um lado; nao parar deixa
  // o suporte virar a malha inteira -- medido: 64 pontos numa malha de 64, e ai'
  // a topologia de vizinhanca deixa de importar, porque tudo acaba entrando.
  int n = 0;
  size_t fim_do_anel = 1;
  for (size_t cab = 0; cab < fila.size () && n < maxpts; cab++) {
    if (cab == fim_do_anel) {            // fechou um anel
      if (n >= minpts) break;
      fim_do_anel = fila.size ();
    }
    const Folha atual = fila[cab];
    t8_locidx_t arvore = 0;
    const t8_element_t *e = t8_forest_get_leaf_element (f, atual, &arvore);
    const t8_eclass_t ec = t8_forest_get_tree_class (f, arvore);

    // o proprio elemento entra no suporte
    double c[3];
    t8_forest_element_centroid (f, arvore, e, c);
    for (int d = 0; d < DIM; d++) pts[n][d] = c[d];
    n++;

    const int nfaces = scheme->element_get_num_faces (ec, e);
    for (int fa = 0; fa < nfaces; fa++) {
      t8_element_t **vz = NULL;
      int *dual = NULL, nviz = 0;
      t8_locidx_t *idx = NULL;
      t8_eclass_t vec;
      t8_forest_leaf_face_neighbors (f, arvore, e,
                                     (const t8_element_t ***) &vz, fa, &dual,
                                     &nviz, &idx, &vec);
      for (int k = 0; k < nviz; k++) {
        const Folha v = idx[k];          // ja' vem no indice da floresta local
        bool ja = false;
        for (size_t q = 0; q < vistos.size () && !ja; q++) ja = (vistos[q] == v);
        if (!ja) { vistos.push_back (v); fila.push_back (v); }
      }
      T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
    }
  }

  return n;
}
