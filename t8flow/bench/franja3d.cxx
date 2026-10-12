// Quanto custa a franja de K camadas em TRES dimensoes.
//
// POR QUE ESTA MEDIDA VEM ANTES DO CODIGO 3D, e nao depois: a carga da troca
// cresce com a POTENCIA DA DIMENSAO, e o numero que sai daqui decide a
// arquitetura do resto.
//
//     vizinhanca de raio r      2D: (2r+1)^2 = 49 para r=2
//                               3D: (2r+1)^3 = 125 para r=2
//
//     vizinhos de face          2D: 4 faces x 2 (transicao 2:1) =  8
//                               3D: 6 faces x 4 (transicao 2:1) = 24
//
// Os dois se multiplicam na carga por folha.  Se o resultado for proibitivo, a
// TROCA precisa ser redesenhada -- e redesenha-la com o estencil e o Poisson ja'
// construidos em cima custa muito mais que descobrir agora.
//
// O QUE ESTE ARREIO NAO E': ele nao e' a `Malha3`.  Fala direto com o t8code, de
// proposito, para que a medida nao dependa de decisoes de arquitetura que ela
// mesma deve informar.  Mede tres coisas:
//
//   1. quantos vizinhos de face uma celula REALMENTE tem, em malha uniforme e
//      adaptada -- o 24 e' o teto teorico, e o que importa e' a media;
//   2. quantas folhas ha' na vizinhanca de raio r, que e' o que uma rodada de
//      troca publica -- o `(2r+1)^3` e' o teto de malha uniforme, e em adaptada
//      passa dele;
//   3. o tempo e o volume de uma troca com essa carga.
//
// Uso:  mpirun -np N ./bench/franja3d [--nivel=L] [--camadas=K] [--caixa]

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <vector>

#include <mpi.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

namespace {

int g_np = 1;
int g_rank = 0;
int g_nivel_caixa = 0;

//! Refina onde o centro cai na caixa `[0, 0.25]^3` -- o padrao real de AMR,
//! porque refina-se onde ha' interface, e isso e' sempre localizado.
int adapt_caixa(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree,
                const t8_eclass_t tclass, t8_locidx_t lelement,
                const t8_scheme_c *scheme, const int is_family,
                const int num_elements, t8_element_t *elements[])
{
  if (scheme->element_get_level(tclass, elements[0]) >= g_nivel_caixa) return 0;
  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(forest_from, ltree, elements[0], c);
  return (c[0] < 0.25 && c[1] < 0.25 && c[2] < 0.25) ? 1 : 0;
}

void relata(const char *etapa, int nivel, int camadas, long globais, long locais,
            long franja, double valor, const char *unidade)
{
  long lmin = 0, lmax = 0, fmax = 0;
  double vmax = 0.0;
  MPI_Allreduce(&locais, &lmin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&locais, &lmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&franja, &fmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&valor, &vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (g_rank == 0) {
    printf("medida %s nivel=%d np=%d camadas=%d globais=%ld locmin=%ld locmax=%ld "
           "franjamax=%ld %s=%.6f\n", etapa, nivel, g_np, camadas, globais, lmin,
           lmax, fmax, unidade, vmax);
    fflush(stdout);
  }
}

}  // namespace

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ERROR);
  t8_init(SC_LP_ERROR);
  MPI_Comm_size(MPI_COMM_WORLD, &g_np);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

  int nivel = 4;
  int camadas = 3;
  bool caixa = false;
  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--camadas=", 10)) camadas = std::atoi(a + 10);
    else if (!std::strcmp(a, "--caixa")) caixa = true;
  }
  g_nivel_caixa = nivel + 1;

  if (g_rank == 0) {
    printf("config dim=3 eclass=HEX nivel=%d np=%d camadas=%d refino=%s\n",
           nivel, g_np, camadas, caixa ? "CAIXA" : "nenhum");
    fflush(stdout);
  }

  // ---- A FLORESTA, com franja nativa de uma camada --------------------------
  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);
  t8_cmesh_new_hypercube(&cmesh, T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_scheme_c *esquema = t8_scheme_new_default();
  t8_forest_t f = t8_forest_new_uniform(cmesh, esquema, nivel, 1, sc_MPI_COMM_WORLD);

  if (caixa) {
    t8_forest_t ref;
    t8_forest_init(&ref);
    t8_forest_set_adapt(ref, f, adapt_caixa, 1);
    t8_forest_commit(ref);
    t8_forest_t eq;
    t8_forest_init(&eq);
    t8_forest_set_balance(eq, ref, 0);
    t8_forest_set_ghost(eq, 1, T8_GHOST_FACES);
    t8_forest_commit(eq);
    f = eq;
  }

  const long nloc = static_cast<long>(t8_forest_get_local_num_leaf_elements(f));
  const long nfan = static_cast<long>(t8_forest_get_num_ghosts(f));
  const long nglob = static_cast<long>(t8_forest_get_global_num_leaf_elements(f));

  // ---- 1. QUANTOS VIZINHOS DE FACE, de verdade -----------------------------
  //
  // O TETO TEORICO E' 24 (seis faces vezes quatro), e o que importa para a carga
  // e' a MEDIA -- dimensionar pelo teto desperdica, e o desperdicio multiplica
  // pelo numero de folhas publicadas.
  long soma_viz = 0, max_viz = 0;
  {
    const t8_locidx_t narv = t8_forest_get_num_local_trees(f);
    for (t8_locidx_t it = 0; it < narv; it++) {
      const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f, it);
      for (t8_locidx_t ie = 0; ie < ne; ie++) {
        const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f, it, ie);
        long deste = 0;
        for (int fa = 0; fa < 6; fa++) {
          t8_element_t **vz = nullptr;
          int *dual = nullptr, nviz = 0;
          t8_locidx_t *idx = nullptr;
          t8_eclass_t vec;
          t8_forest_leaf_face_neighbors(f, it, e,
                                        const_cast<const t8_element_t ***>(&vz), fa,
                                        &dual, &nviz, &idx, &vec);
          deste += nviz;
          T8_FREE(vz); T8_FREE(dual); T8_FREE(idx);
        }
        soma_viz += deste;
        if (deste > max_viz) max_viz = deste;
      }
    }
  }
  if (g_rank == 0) {
    printf("vizinhos nivel=%d np=%d media=%.2f maximo=%ld (teto teorico 24)\n",
           nivel, g_np, nloc > 0 ? double(soma_viz) / nloc : 0.0, max_viz);
    fflush(stdout);
  }

  // ---- 2. A CARGA DE CADA RODADA -------------------------------------------
  //
  // Publicar a vizinhanca de raio `r` custa `folhas(r) * campos` por folha
  // local.  `folhas(r)` e' `(2r+1)^3` em malha uniforme e mais em adaptada, e o
  // numero de campos e' 7 proprios mais dois por vizinho.
  const int CAMPOS = 7 + static_cast<int>(max_viz) * 2;
  for (int r = 1; r < camadas; r++) {
    const long teto_unif = 1;
    long folhas_r = teto_unif;
    for (int k = 0; k < 3; k++) folhas_r *= (2 * r + 1);
    const double kb = folhas_r * CAMPOS * 8.0 / 1024.0;
    if (g_rank == 0) {
      printf("carga rodada=%d raio=%d folhas=%ld campos=%d por_folha=%.1f KB "
             "total_local=%.1f MB\n", r, r, folhas_r, CAMPOS, kb,
             kb * nloc / 1024.0);
      fflush(stdout);
    }
  }

  // ---- 3. O TEMPO DE UMA TROCA com a carga da ultima rodada ----------------
  {
    long folhas_r = 1;
    for (int k = 0; k < 3; k++) folhas_r *= (2 * (camadas - 1) + 1);
    const size_t passo = static_cast<size_t>(folhas_r) * static_cast<size_t>(CAMPOS);

    double melhor = 1e300;
    for (int rep = 0; rep < 3; rep++) {
      sc_array_t *arr = sc_array_new_count(passo * sizeof(double),
                                           static_cast<size_t>(nloc + nfan));
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      if (nfan > 0) t8_forest_ghost_exchange_data(f, arr);
      const double dt = MPI_Wtime() - t0;
      sc_array_destroy(arr);
      if (rep > 0 && dt < melhor) melhor = dt;
    }
    relata("troca", nivel, camadas, nglob, nloc, nfan, melhor, "seg");
    if (g_rank == 0) {
      printf("volume total publicado: %.1f MB\n",
             passo * 8.0 * nloc / (1024.0 * 1024.0));
      fflush(stdout);
    }
  }

  t8_forest_unref(&f);
  sc_finalize();
  MPI_Finalize();
  return 0;
}
