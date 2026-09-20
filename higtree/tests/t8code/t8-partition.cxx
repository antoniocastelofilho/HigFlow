// O t8code sob particionamento real.  Ver o .h para o porque.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-partition.h"
#include "t8-forest-cache.hxx"

#include <cmath>

#if DIM == 2
#define ECLASSE T8_ECLASS_QUAD
#else
#define ECLASSE T8_ECLASS_HEX
#endif

#define NIVEL 3                 // 2^3 = 8 celulas por direcao

static t8_forest_t g_f = NULL;

static double
campo (const double p[3])
{
  return 1.0 + 2.0 * p[0] + 3.0 * p[1];
}

extern "C" void
t8_particao_monta (void)
{
  if (g_f != NULL) return;
  t8_inicializa_uma_vez ();

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, ECLASSE, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_scheme_c *scheme = t8_scheme_new_default ();

  // do_face_ghost = 1.  Sem a camada de ghost, "face sem vizinho" deixaria de
  // distinguir contorno do dominio de fronteira de particao -- que e' o limite
  // anotado na clausula C12.
  g_f = t8_forest_new_uniform (cmesh, scheme, NIVEL, 1, sc_MPI_COMM_WORLD);
}

extern "C" long
t8_particao_num_locais (void)
{
  return (long) t8_forest_get_local_num_leaf_elements (g_f);
}

extern "C" double
t8_particao_volume_local (void)
{
  double v = 0.0;
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (g_f);
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (g_f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (g_f, it, ie);
      v += t8_forest_element_volume (g_f, it, e);
    }
  }
  return v;
}

extern "C" long
t8_particao_num_ghosts (void)
{
  return (long) t8_forest_get_num_ghosts (g_f);
}

extern "C" long
t8_particao_suporte_alcanca_ghost (void)
{
  // Um vizinho de face e' de GHOST quando o indice que a API devolve cai depois
  // do ultimo elemento local.  E' o mesmo criterio de `lid >= n_local` que o
  // teste da franja do MTree usa, e nao precisa de API nova.
  const t8_locidx_t nlocal = t8_forest_get_local_num_leaf_elements (g_f);
  const t8_scheme_c *scheme = t8_forest_get_scheme (g_f);
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (g_f);

  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (g_f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (g_f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (g_f, it, ie);
      const int nfaces = scheme->element_get_num_faces (ec, e);
      for (int fa = 0; fa < nfaces; fa++) {
        t8_element_t **vz = NULL;
        int *dual = NULL, nviz = 0;
        t8_locidx_t *idx = NULL;
        t8_eclass_t vec;
        t8_forest_leaf_face_neighbors (g_f, it, e, (const t8_element_t ***) &vz,
                                       fa, &dual, &nviz, &idx, &vec);
        bool achou = false;
        for (int k = 0; k < nviz; k++) if (idx[k] >= nlocal) achou = true;
        T8_FREE (vz); T8_FREE (dual); T8_FREE (idx);
        if (achou) return 1;
      }
    }
  }
  return 0;
}

extern "C" int
t8_particao_interpola (const Point p, double *valor)
{
  double ponto[3] = { 0.0, 0.0, 0.0 };
  for (int d = 0; d < DIM; d++) ponto[d] = p[d];

  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (g_f);
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (g_f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (g_f, it, ie);
      int dentro = 0;
      t8_forest_element_points_inside (g_f, it, e, ponto, 1, &dentro, 1e-12);
      if (!dentro) continue;
      // O ponto cai estritamente dentro de uma celula (o teste o escolhe assim),
      // entao vale o valor do campo no centroide dela -- suficiente para a P4,
      // que pergunta se o valor MUDA com a particao, nao qual e' o valor.
      double c[3];
      t8_forest_element_centroid (g_f, it, e, c);
      *valor = campo (c);
      return 1;
    }
  }
  return 0;
}

extern "C" void
t8_particao_encerra (void)
{
  if (g_f != NULL) { t8_forest_unref (&g_f); g_f = NULL; }
}
