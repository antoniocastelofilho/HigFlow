// Classificacao de ponto no t8code.  Ver o .h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-point-class.h"
#include "t8-forest-cache.hxx"

#include <cmath>

extern "C" t8_classe_de_ponto
t8_classifica_ponto (int ntrees_x, const Point p)
{
  t8_inicializa_uma_vez ();
  t8_forest_t f = t8_floresta_de_teste (ntrees_x);
  if (f == NULL) return T8_FORA;

  double ponto[3] = { 0.0, 0.0, 0.0 };
  for (int d = 0; d < DIM; d++) ponto[d] = p[d];

  const t8_scheme_c *scheme = t8_forest_get_scheme (f);
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);

  bool contido = false;
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      int dentro = 0;
      t8_forest_element_points_inside (f, it, e, ponto, 1, &dentro, 1e-12);
      if (!dentro) continue;
      contido = true;

      // Esta' no elemento.  Sobre qual face, se alguma?  Face SEM VIZINHO e'
      // contorno do dominio; face com vizinho e' interface interna, e o ponto
      // segue DENTRO -- e' aqui que este criterio se separa do da caixa.
      const int nfaces = scheme->element_get_num_faces (ec, e);
      for (int fa = 0; fa < nfaces; fa++) {
        double c[3], nrm[3];
        t8_forest_element_face_centroid (f, it, e, fa, c);
        t8_forest_element_face_normal (f, it, e, fa, nrm);

        // ponto no plano da face?  (x - c) . n == 0
        double proj = 0.0;
        for (int d = 0; d < DIM; d++) proj += (ponto[d] - c[d]) * nrm[d];
        if (std::fabs (proj) > 1e-12) continue;

        t8_element_t **vz = NULL;
        int *dual = NULL, nviz = 0;
        t8_locidx_t *idx = NULL;
        t8_eclass_t vec;
        t8_forest_leaf_face_neighbors (f, it, e, (const t8_element_t ***) &vz,
                                       fa, &dual, &nviz, &idx, &vec);
        const bool sem_vizinho = (nviz == 0);
        T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
        if (sem_vizinho) return T8_NO_CONTORNO;
      }
    }
  }
  return contido ? T8_DENTRO : T8_FORA;
}
