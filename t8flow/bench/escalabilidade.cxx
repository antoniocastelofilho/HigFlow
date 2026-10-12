// MEDIDA DE ESCALABILIDADE DA FLORESTA DO T8CODE.
//
// Existe para responder a UMA pergunta, que o doc/projeto-t8flow.md levanta e
// nao pode responder por argumento: **a floresta escala no regime 3D?**  O 2D do
// HiGFlow rendeu so' 1,62x em np=2 e 1,68x em np=4, e a explicacao medida foi
// TAMANHO DO CASO -- 9377 celulas divididas por quatro dao ~2300 por processo, e
// a comunicacao come o que a divisao da'.  O dominio 3D tem 375 mil celulas, que
// e' o regime onde a divisao deveria vencer.  Sem medir, o escalonamento do
// t8Flow e' promessa.
//
// O QUE ISTO MEDE, e so' isto:
//
//   1. CONSTRUCAO  da floresta uniforme de nivel L em COMM_WORLD.  O t8code ja'
//      reparte na criacao, entao esta e' a etapa que a promessa de escalonamento
//      cobre: nenhum processo materializa a malha inteira;
//   2. PERCURSO    de todas as folhas locais calculando o centroide.  E' o proxy
//      do laco de celulas que o solver e o discretizador executam, e e' trabalho
//      puramente local -- o teto do que a divisao pode dar;
//   3. REFINO+EQUILIBRIO  de uma caixa, que e' comunicacao de verdade: o balance
//      precisa dos vizinhos, e e' aqui que a escalabilidade costuma morrer;
//   4. REPARTICAO depois do refino, e o percurso na malha ja' adaptada.  O refino
//      local CONCENTRA celulas onde a caixa cai, e a particao de antes dele nao
//      sabe disso -- esta etapa mede se reparticionar devolve o equilibrio, e
//      quanto custa.
//
// O QUE ISTO NAO MEDE, e dizer importa mais que medir: nao ha' fisica, nao ha'
// solver, nao ha' franja, e nao ha' estencil.  O numero daqui e' um TETO para a
// escalabilidade do t8Flow, nao uma previsao dela -- a franja e o solver linear
// entram depois e sao justamente as partes que comunicam mais.
//
// SAIDA parseavel, uma linha por etapa:
//
//     medida <etapa> dim=<2|3> nivel=<L> np=<n> globais=<N> locmin=<a> locmax=<b> seg=<t>
//
// `locmin`/`locmax` sao as folhas locais minima e maxima entre os ranks: e' o
// balanceamento, e um desequilibrio ali limita o ganho independentemente do
// resto.  O tempo e' o MAXIMO entre os ranks, que e' o que o relogio de parede ve'.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <mpi.h>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

namespace {

int g_dim = 3;
double g_caixa_lo[3] = {0.60, 0.60, 0.60};
double g_caixa_hi[3] = {1.00, 1.00, 1.00};
int g_nivel_alvo = 0;

//! Refina dentro da caixa, ate' `g_nivel_alvo`.  Mesmo criterio do malha2d: o
//! centro da celula, nao a interseccao -- ver o comentario la'.
int adapt_caixa(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                const t8_eclass_t tree_class, t8_locidx_t /*lelement_id*/,
                const t8_scheme_c *scheme, const int /*is_family*/,
                const int /*num_elements*/, t8_element_t *elements[])
{
  if (scheme->element_get_level(tree_class, elements[0]) >= g_nivel_alvo) return 0;
  double c[3];
  t8_forest_element_centroid(forest_from, which_tree, elements[0], c);
  for (int d = 0; d < g_dim; d++)
    if (c[d] < g_caixa_lo[d] || c[d] > g_caixa_hi[d]) return 0;
  return 1;
}

//! Imprime uma linha de medida, com os extremos e o tempo JA' reduzidos.
void relata(const char *etapa, int nivel, int np, long globais, long locais, double seg)
{
  long lmin = 0, lmax = 0;
  double tmax = 0.0;
  MPI_Allreduce(&locais, &lmin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&locais, &lmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&seg, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    printf("medida %s dim=%d nivel=%d np=%d globais=%ld locmin=%ld locmax=%ld seg=%.6f\n",
           etapa, g_dim, nivel, np, globais, lmin, lmax, tmax);
    fflush(stdout);
  }
}

}  // namespace

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  sc_init(sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ERROR);
  t8_init(SC_LP_ERROR);

  int np = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nivel = 5;
  for (int i = 1; i < argc; i++) {
    if (std::strncmp(argv[i], "--dim=", 6) == 0) g_dim = std::atoi(argv[i] + 6);
    else if (std::strncmp(argv[i], "--nivel=", 8) == 0) nivel = std::atoi(argv[i] + 8);
  }
  if (g_dim != 2 && g_dim != 3) {
    if (rank == 0) fprintf(stderr, "dim tem de ser 2 ou 3\n");
    MPI_Finalize();
    return 2;
  }
  g_nivel_alvo = nivel + 1;

  // A CONFIGURACAO VAI PARA A SAIDA, sempre.  Build sem o que se quer medir e'
  // o primeiro dos modos de medicao vazia catalogados neste projeto: duas
  // corridas "estaveis" minhas foram malha uniforme medida como se fosse t8code.
  if (rank == 0) {
    printf("config dim=%d nivel=%d np=%d eclass=%s\n", g_dim, nivel, np,
           g_dim == 2 ? "QUAD" : "HEX");
    fflush(stdout);
  }

  const t8_eclass_t ECLASSE = (g_dim == 2) ? T8_ECLASS_QUAD : T8_ECLASS_HEX;

  // ---- 1. CONSTRUCAO -------------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();

  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);
  t8_cmesh_new_hypercube(&cmesh, ECLASSE, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_scheme_c *esquema = t8_scheme_new_default();
  t8_forest_t f = t8_forest_new_uniform(cmesh, esquema, nivel, 0, sc_MPI_COMM_WORLD);

  double t_constroi = MPI_Wtime() - t0;
  long loc = static_cast<long>(t8_forest_get_local_num_leaf_elements(f));
  long glob = static_cast<long>(t8_forest_get_global_num_leaf_elements(f));
  relata("construcao", nivel, np, glob, loc, t_constroi);

  // ---- 2. PERCURSO (trabalho local) ---------------------------------------
  // Tres repeticoes, e o MENOR tempo vale: a maquina tem vizinhos (outras
  // corridas), e o minimo e' menos sensivel a eles que a media.
  double melhor = 1e300;
  for (int rep = 0; rep < 3; rep++) {
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    double soma = 0.0;
    const t8_locidx_t narv = t8_forest_get_num_local_trees(f);
    for (t8_locidx_t it = 0; it < narv; it++) {
      const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f, it);
      for (t8_locidx_t ie = 0; ie < ne; ie++) {
        const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f, it, ie);
        double c[3] = {0.0, 0.0, 0.0};
        t8_forest_element_centroid(f, it, e, c);
        soma += c[0] + c[1] + c[2];
      }
    }
    const double dt = MPI_Wtime() - t0;
    if (dt < melhor) melhor = dt;
    // `soma` existe para o compilador nao descartar o laco -- com -O1 ele
    // descartaria, e a medida seria de um laco vazio.
    if (soma == 1e300 && rank < 0) printf("impossivel\n");
  }
  relata("percurso", nivel, np, glob, loc, melhor);

  // ---- 3. REFINO + EQUILIBRIO (comunicacao) -------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();

  t8_forest_t refinada;
  t8_forest_init(&refinada);
  t8_forest_set_adapt(refinada, f, adapt_caixa, 1);
  t8_forest_commit(refinada);

  t8_forest_t equilibrada;
  t8_forest_init(&equilibrada);
  t8_forest_set_balance(equilibrada, refinada, 0);
  t8_forest_commit(equilibrada);

  const double t_refina = MPI_Wtime() - t0;
  loc = static_cast<long>(t8_forest_get_local_num_leaf_elements(equilibrada));
  glob = static_cast<long>(t8_forest_get_global_num_leaf_elements(equilibrada));
  relata("refino", nivel, np, glob, loc, t_refina);

  // ---- 4. REPARTICAO depois do refino -------------------------------------
  //
  // ESTA ETAPA ENTROU POR CAUSA DE UMA MEDIDA, nao por completude.  O refino
  // local CONCENTRA celulas onde a caixa cai, e a particao de antes dele nao sabe
  // disso: medido em np=4, nivel 7, o rank mais carregado ficou com 1,45 MILHAO de
  // folhas contra 524 mil do menos carregado -- 2,8x de desequilibrio.  Como o
  // relogio de parede e' o do rank mais lento, esse desequilibrio sozinho explica
  // a perda de escalabilidade do refino.
  //
  // A pergunta que esta etapa responde: reparticionar resolve, e a que custo?
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();

  t8_forest_t repartida;
  t8_forest_init(&repartida);
  t8_forest_set_partition(repartida, equilibrada, 0);
  t8_forest_commit(repartida);

  const double t_reparte = MPI_Wtime() - t0;
  loc = static_cast<long>(t8_forest_get_local_num_leaf_elements(repartida));
  glob = static_cast<long>(t8_forest_get_global_num_leaf_elements(repartida));
  relata("reparticao", nivel, np, glob, loc, t_reparte);

  // E o percurso DEPOIS de reparticionar: e' o que o solver e o discretizador
  // veriam numa malha adaptada de verdade.
  double melhor2 = 1e300;
  for (int rep = 0; rep < 3; rep++) {
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    double soma = 0.0;
    const t8_locidx_t narv2 = t8_forest_get_num_local_trees(repartida);
    for (t8_locidx_t it = 0; it < narv2; it++) {
      const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(repartida, it);
      for (t8_locidx_t ie = 0; ie < ne; ie++) {
        const t8_element_t *e = t8_forest_get_leaf_element_in_tree(repartida, it, ie);
        double c[3] = {0.0, 0.0, 0.0};
        t8_forest_element_centroid(repartida, it, e, c);
        soma += c[0] + c[1] + c[2];
      }
    }
    const double dt = MPI_Wtime() - t0;
    if (dt < melhor2) melhor2 = dt;
    if (soma == 1e300 && rank < 0) printf("impossivel\n");
  }
  relata("percurso-adaptado", nivel, np, glob, loc, melhor2);

  t8_forest_unref(&repartida);

  if (rank == 0) { printf("fim\n"); fflush(stdout); }

  sc_finalize();
  MPI_Finalize();
  return 0;
}
