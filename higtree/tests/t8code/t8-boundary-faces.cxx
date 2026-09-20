// Faces de contorno na floresta do t8code.  Ver o .h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-boundary-faces.h"
#include "t8-forest-cache.hxx"

extern "C" int
t8_faces_de_contorno (int ntrees_x, int maxn, Point centros[], Point normais[])
{
  t8_inicializa_uma_vez ();
  t8_forest_t f = t8_floresta_de_teste (ntrees_x);
  if (f == NULL) return 0;

  const t8_scheme_c *scheme = t8_forest_get_scheme (f);
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);

  int n = 0;
  for (t8_locidx_t it = 0; it < ntrees && n < maxn; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne && n < maxn; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      const int nfaces = scheme->element_get_num_faces (ec, e);
      for (int fa = 0; fa < nfaces && n < maxn; fa++) {
        t8_element_t **vz = NULL;
        int *dual = NULL, nviz = 0;
        t8_locidx_t *idx = NULL;
        t8_eclass_t vec;
        t8_forest_leaf_face_neighbors (f, it, e, (const t8_element_t ***) &vz,
                                       fa, &dual, &nviz, &idx, &vec);
        const bool contorno = (nviz == 0);
        T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
        if (!contorno) continue;

        // Face sem vizinho: e' contorno do dominio.  Note que isso vale porque a
        // malha de teste e' de UM processo; sob particionamento, face sem vizinho
        // LOCAL pode ser face de franja, e a distincao passa a exigir o ghost.
        double c[3], nrm[3];
        t8_forest_element_face_centroid (f, it, e, fa, c);
        t8_forest_element_face_normal (f, it, e, fa, nrm);
        for (int d = 0; d < DIM; d++) {
          centros[n][d] = c[d];
          normais[n][d] = nrm[d];
        }
        n++;
      }
    }
  }
  return n;
}
