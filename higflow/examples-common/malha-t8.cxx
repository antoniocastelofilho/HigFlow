// FONTES DE MALHA DO T8CODE, compartilhadas pelos exemplos.
//
// A malha NAO esta' decorada aqui: a especificacao vem da informacao AMR que o
// solver leria de qualquer jeito.  E' o que permite a mesma fonte servir a varios
// exemplos -- ela reproduz o arquivo daquele exemplo, seja ele qual for.
//
// POR QUE ISSO PODE SER VERIFICADO SEM REFERENCIA NOVA.  Para uma malha UNIFORME,
// um cmesh de BRICK do t8code com a mesma grade tem exatamente as mesmas celulas
// que o arquivo descreve.  Mesma malha, produtor diferente: a saida tem de bater
// com a referencia ja' gravada.  Gerar referencia a partir do codigo recem-escrito
// so' guardaria contra regressao futura.
//
// LIMITE, e ele decide quais exemplos podem usar isto: malha UNIFORME, um nivel e
// um patch por bloco.  Levantado nos exemplos: oito tem um bloco uniforme, o
// Newt_contraction tem DOIS (que esta fonte cobre, uma arvore por bloco), e cinco
// nao tem arquivo de dominio.  Com refino a fonte RECUSA em vez de aproximar --
// uma malha diferente daria numeros diferentes, e o silencio ai' seria pior.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "hig-flow-kernel.h"
#include "higtree.h"
#include "higtree-io.h"
#include "coord.h"
#include "t8-mesh-rank.h"
#include "t8-particao-grafo.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <mpi.h>

static int g_iniciado = 0;

static void
inicializa (void)
{
  if (g_iniciado) return;
  g_iniciado = 1;
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
  t8_init (SC_LP_ERROR);
}

// A especificacao de UM bloco: caixa e grade.  Recusa o que nao for uniforme.
static int
espec (higio_amr_info *mi, Point lo, Point hi, int nb[DIM], const char *quem)
{
  if (mi->numlevels != 1) {
    fprintf (stderr, "%s: a malha tem %d niveis; esta fonte so' reproduz malha "
                     "uniforme.  Recusar e' melhor que aproximar: malha diferente "
                     "da' numeros diferentes.\n", quem, mi->numlevels);
    return 0;
  }
  if (mi->levels[0].numpatches != 1) {
    fprintf (stderr, "%s: o nivel 0 tem %d patches; esperado 1\n", quem,
             mi->levels[0].numpatches);
    return 0;
  }
  POINT_ASSIGN (lo, mi->l);
  POINT_ASSIGN (hi, mi->h);
  POINT_ASSIGN (nb, mi->levels[0].patches[0].patchsize);
  for (int d = 0; d < DIM; d++) if (nb[d] < 1 || !(hi[d] > lo[d])) return 0;
  return 1;
}

#if DIM == 2
#define BRICK(cm, nb, comm) t8_cmesh_new_brick_2d (cm, (nb)[0], (nb)[1], 0, 0, comm)
#else
#define BRICK(cm, nb, comm) t8_cmesh_new_brick_3d (cm, (nb)[0], (nb)[1], (nb)[2], 0, 0, 0, comm)
#endif

// Um bloco -> uma arvore hig, materializada a partir das FOLHAS da floresta.
// Nada aqui refina por conta propria: se o t8code produzisse outra malha, a
// arvore sairia diferente -- e' o que impede este modulo de "passar" reproduzindo
// o arquivo por fora do t8code.
static hig_cell *
bloco (const Point lo, const Point hi, const int nb[DIM], const char *quem)
{
  inicializa ();
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  BRICK (cmesh, nb, sc_MPI_COMM_SELF);
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_SELF);

  hig_cell *raiz = hig_create_root ((real *) lo, (real *) hi);
  int nc[DIM];
  for (int d = 0; d < DIM; d++) nc[d] = nb[d];
  hig_refine_uniform (raiz, nc);

  long folhas = 0, ruins = 0;
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      Point p;
      for (int d = 0; d < DIM; d++)
        p[d] = lo[d] + (c[d] / (double) nb[d]) * (hi[d] - lo[d]);
      hig_cell *cel = hig_get_cell_with_point (raiz, p);
      if (cel == NULL) { ruins++; continue; }
      Point ce;
      hig_get_center (cel, ce);
      for (int d = 0; d < DIM; d++)
        if (fabs (ce[d] - p[d]) > 1.0e-12) { ruins++; break; }
      folhas++;
    }
  }
  t8_forest_unref (&f);

  long esperado = 1;
  for (int d = 0; d < DIM; d++) esperado *= nb[d];
  if (folhas != esperado || ruins != 0) {
    fprintf (stderr, "%s: a floresta deu %ld folha(s) e %ld divergencia(s) contra "
                     "a arvore (esperado %ld e 0)\n", quem, folhas, ruins, esperado);
    hig_destroy (raiz);
    return NULL;
  }
  return raiz;
}

// ---------------------------------------------------------------- em serie
extern "C" int
malha_t8_uniforme (void *ctx, higio_amr_info **mi, int numhigs,
                   hig_cell **arvores, int max)
{
  (void) ctx;
  int rank = 0, np = 1;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np);

  // CADA BLOCO E' CONTRIBUIDO POR UM RANK SO', pela mesma regra do caminho AMR
  // (`for i = myrank; i < numhigs; i += ntasks`).  Sem isso todo rank entrega a
  // malha inteira e o `lbal` recebe `np` copias dela.
  //
  // O defeito ficou latente ate' agora porque a fonte em serie so' tinha sido
  // exercitada em np=1, onde a regra e' trivialmente satisfeita.  Em np=2 a
  // referencia reprovou.
  int n = 0;
  for (int i = rank; i < numhigs; i += np) {
    if (n >= max) return 0;
    Point lo, hi;
    int nb[DIM];
    if (!espec (mi[i], lo, hi, nb, "malha_t8_uniforme")) return 0;
    hig_cell *t = bloco (lo, hi, nb, "malha_t8_uniforme");
    if (t == NULL) { for (int k = 0; k < n; k++) hig_destroy (arvores[k]); return 0; }
    arvores[n++] = t;
  }
  return n;
}

// -------------------------------------------------------------- por rank
// FONTE REFINADA EM TORNO DO CORPO.  Aqui a fonte deixa de REPRODUZIR o arquivo
// AMR e passa a PRODUZIR uma malha: o criterio de refino vem da caixa dada, nao
// do arquivo.  E' a diferenca que justifica existir uma fonte separada em vez de
// afrouxar o `espec`.
//
// A caixa e o numero de refinos vem do ambiente, para o exemplo nao precisar
// saber de t8code:
//
//   HIGFLOW_REFINO_CAIXA="x0,y0,x1,y1"   (ou com z0,z1 em 3D)
//   HIGFLOW_REFINO_NIVEIS=1
//
// Sem a caixa, recusa -- refinar "em algum lugar" nao e' comportamento util.
extern "C" int
malha_t8_refinada (void *ctx, higio_amr_info **mi, int numhigs,
                   hig_cell **arvores, int max)
{
  (void) ctx;
  if (numhigs != 1) {
    fprintf (stderr, "malha_t8_refinada: %d blocos; esta fonte cobre um so'\n",
             numhigs);
    return 0;
  }

  const char *scx = getenv ("HIGFLOW_REFINO_CAIXA");
  if (scx == NULL) {
    fprintf (stderr, "malha_t8_refinada: defina HIGFLOW_REFINO_CAIXA="
                     "\"x0,y0,x1,y1\" (2D) -- refinar sem caixa nao e' util\n");
    return 0;
  }
  Point cx_lo, cx_hi;
  {
    double v[2*DIM];
    const char *p = scx;
    for (int i = 0; i < 2*DIM; i++) {
      char *fim;
      v[i] = strtod (p, &fim);
      if (fim == p) {
        fprintf (stderr, "malha_t8_refinada: HIGFLOW_REFINO_CAIXA precisa de %d "
                         "numeros separados por virgula\n", 2*DIM);
        return 0;
      }
      p = (*fim == ',') ? fim + 1 : fim;
    }
    for (int d = 0; d < DIM; d++) { cx_lo[d] = v[d]; cx_hi[d] = v[DIM + d]; }
  }
  const char *sn = getenv ("HIGFLOW_REFINO_NIVEIS");
  const int niveis = (sn != NULL) ? atoi (sn) : 1;

  Point lo, hi; int nb[DIM];
  if (!espec (mi[0], lo, hi, nb, "malha_t8_refinada")) return 0;

  t8_producao_rank prod;
  if (!t8_produz_por_rank_brick_refinado (lo, hi, nb, cx_lo, cx_hi, niveis, &prod))
    return 0;

  // NAO TRUNCAR EM SILENCIO.  O solver passa max = 64 arvores; uma caixa de
  // refino grande produz MUITO mais caixas completas que isso, e descartar as
  // excedentes deixa buracos na malha.
  //
  // MEDIDO: com a caixa [1;8]x[1;3] o corpo ficou numa regiao descartada, a
  // interpolacao nao achou suporte nenhum, e a corrida deu Cd = 0, Cl = 0 e
  // residuo 0 -- tres zeros perfeitos, que parecem "ainda nao comecou" e nao
  // "a malha esta' furada".  Truncamento calado e' pior que falha.
  if (prod.n_locais > max) {
    int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
      fprintf (stderr,
        "malha_t8_refinada: a malha refinada precisa de %d caixas e o solver "
        "aceita %d.\n"
        "  A caixa de refino e' grande demais para esta representacao.\n"
        "  Ou diminua HIGFLOW_REFINO_CAIXA, ou aumente `raizes[]` em\n"
        "  higflow_partition_domain (hig-flow-kernel.c).\n",
        prod.n_locais, max);
    t8_producao_rank_destroi (&prod);
    return 0;
  }

  // Os limites de CADA arvore local, nao so' a contagem.  Uma regiao de refino
  // vira varias arvores por rank, e a costura entre duas arvores adjacentes e'
  // um sitio que "2 caixas" sozinho nao localiza.
  if (getenv ("HIGFLOW_REFINO_LISTA") != NULL) {
    int rk; MPI_Comm_rank (MPI_COMM_WORLD, &rk);
    // A contagem do T8CODE e' a autoridade.  A soma das folhas materializadas
    // tem de bater com ela: a mais e' sobreposicao de arvore, a menos e' buraco.
    if (rk == 0)
      printf ("GLOBAL t8code = %ld folhas\n", prod.n_global);
    for (int i = 0; i < prod.n_locais; i++) {
      Point tl, th;
      hig_get_lowpoint (prod.locais[i], tl);
      hig_get_highpoint (prod.locais[i], th);
      // A CONTAGEM DE FOLHAS E' A CARGA.  Estas arvores sao as locais DESTE
      // rank, entao a soma delas e' o que este rank vai resolver.  Sem este
      // numero, "desbalanceamento" fica em suposicao.
      long folhas = 0;
      {
        higcit_celliterator *it;
        for (it = higcit_create_all_leaves (prod.locais[i]);
             !higcit_isfinished (it); higcit_nextcell (it)) folhas++;
        higcit_destroy (it);
      }
      printf ("ARVORE rank %d [%d] x[%.4f,%.4f] y[%.4f,%.4f] folhas=%ld\n",
              rk, i, tl[0], th[0], tl[1], th[1], folhas);
      fflush (stdout);
    }
  }

  int n = 0;
  for (int i = 0; i < prod.n_locais; i++) arvores[n++] = prod.locais[i];
  // As arvores passam para o chamador; nao destruir aqui.
  prod.n_locais = 0;
  t8_producao_rank_destroi (&prod);

  int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (rank == 0)
    printf ("=+=+=+= malha do t8code REFINADA: %d niveis, %d caixas, na regiao "
            "[%g,%g]x[%g,%g] =+=+=+=\n", niveis, n,
            (double) cx_lo[0], (double) cx_hi[0],
            (double) cx_lo[1], (double) cx_hi[1]);
  return n;
}

extern "C" int
malha_t8_por_rank (void *ctx, higio_amr_info **mi, int numhigs,
                   hig_cell **arvores, int max)
{
  (void) ctx;
  if (numhigs != 1) {
    fprintf (stderr, "malha_t8_por_rank: %d blocos; esta fonte cobre um so'\n",
             numhigs);
    return 0;
  }
  Point lo, hi;
  int nb[DIM];
  if (!espec (mi[0], lo, hi, nb, "malha_t8_por_rank")) return 0;

  t8_producao_rank p;
  if (!t8_produz_por_rank_brick (lo, hi, nb, &p)) return 0;
  if (p.base_dividida != 0 || p.n_locais > max) {
    t8_producao_rank_destroi (&p);
    return 0;
  }
  for (int i = 0; i < p.n_locais; i++) arvores[i] = p.locais[i];
  const int n = p.n_locais;
  p.n_locais = 0;                    // as arvores passam a ser do `lbal`
  t8_producao_rank_destroi (&p);
  return n;
}

// ------------------------------------------------------------- particao
extern "C" int
particao_t8 (void *ctx, higio_amr_info **mi, int numhigs, sim_domain *sd,
             partition_graph *pg)
{
  (void) ctx;
  if (numhigs != 1) {
    fprintf (stderr, "particao_t8: %d blocos; esta fonte cobre um so'\n", numhigs);
    return 0;
  }
  Point lo, hi;
  int nb[DIM];
  if (!espec (mi[0], lo, hi, nb, "particao_t8")) return 0;
  return t8_monta_dominio_particionado (lo, hi, nb, sd, pg);
}

// ---------------------------------------------------------------------------
// INSTALADOR: uma linha por exemplo.
//
// -----------------------------------------------------------------------------
// A fonte guiada por CRITERIO (F3 do AMR dinamico).
//
// O exemplo avalia o criterio na malha corrente, preenche a tabela de
// nivel-alvo por celula base e chama `malha_t8_criterio_define` ANTES de
// reconstruir; a reconstrucao chega ao gancho de fonte de malha, que produz a
// floresta guiada pela tabela.  A tabela e' reduzida com MAX entre os ranks
// aqui dentro -- o chamador so' preenche o que e' dele.
// -----------------------------------------------------------------------------
static signed char *g_crit_tabela = NULL;
static long         g_crit_n = 0;
static int          g_crit_max_nivel = 0;

extern "C" void
malha_t8_criterio_define (const signed char *tabela, long n, int max_nivel)
{
  if (g_crit_tabela == NULL || g_crit_n != n) {
    free (g_crit_tabela);
    g_crit_tabela = (signed char *) malloc ((size_t) n);
    g_crit_n = n;
  }
  MPI_Allreduce (tabela, g_crit_tabela, (int) n, MPI_SIGNED_CHAR, MPI_MAX,
                 MPI_COMM_WORLD);
  g_crit_max_nivel = max_nivel;
}

//! A tabela corrente, para o chamador decidir se algo mudou (custo -> 0 do
//! portao da F3: tabela igual = remalhamento pulado).
extern "C" const signed char *
malha_t8_criterio_tabela (long *n)
{
  if (n) *n = g_crit_n;
  return g_crit_tabela;
}

static int
malha_t8_criterio (void *ctx, higio_amr_info **mi, int numhigs,
                   hig_cell **arvores, int max)
{
  (void) ctx;
  if (numhigs != 1) {
    fprintf (stderr, "malha_t8_criterio: so' um dominio (ha' %d)\n", numhigs);
    return 0;
  }
  if (g_crit_tabela == NULL) {
    fprintf (stderr, "malha_t8_criterio: chame malha_t8_criterio_define antes\n");
    return 0;
  }
  Point lo, hi; int nb[DIM];
  if (!espec (mi[0], lo, hi, nb, "malha_t8_criterio")) return 0;
  long esperado = 1;
  for (int d = 0; d < DIM; d++) esperado *= nb[d];
  if (esperado != g_crit_n) {
    fprintf (stderr, "malha_t8_criterio: tabela tem %ld celulas, dominio tem %ld\n",
             g_crit_n, esperado);
    return 0;
  }

  t8_producao_rank prod;
  if (!t8_produz_por_rank_brick_criterio (lo, hi, nb, g_crit_tabela,
                                          g_crit_max_nivel, &prod))
    return 0;
  if (prod.n_locais > max) {
    int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
      fprintf (stderr, "malha_t8_criterio: %d caixas e o solver aceita %d\n",
               prod.n_locais, max);
    t8_producao_rank_destroi (&prod);
    return 0;
  }
  if (getenv ("HIGFLOW_REFINO_LISTA") != NULL) {
    int rk; MPI_Comm_rank (MPI_COMM_WORLD, &rk);
    if (rk == 0)
      printf ("GLOBAL t8code = %ld folhas\n", prod.n_global);
  }
  int n = 0;
  for (int i = 0; i < prod.n_locais; i++) arvores[n++] = prod.locais[i];
  return n;
}

// Le' HIGFLOW_MALHA e instala o gancho correspondente.  Sem a variavel, nada e'
// instalado e o exemplo segue lendo o arquivo AMR -- que e' o que a suite padrao
// exercita.
//
// Existe para que aderir custe UMA chamada: sem ele, cada exemplo repetiria o
// `strcmp` de quatro fontes, e repetir isso em oito exemplos e' como os campos do
// `ns->cc` se perderam.
extern "C" void
malha_t8_instala (higflow_solver *ns, int myrank)
{
  const char *fonte = getenv ("HIGFLOW_MALHA");
  if (fonte == NULL) return;

  if (strcmp (fonte, "t8code") == 0) {
    if (myrank == 0) printf ("=+=+=+= Malha do t8code (serie) =+=+=+=\n");
    higflow_set_fonte_de_malha (ns, malha_t8_uniforme, NULL);
  } else if (strcmp (fonte, "t8code-rank") == 0) {
    if (myrank == 0) printf ("=+=+=+= Malha do t8code (por rank) =+=+=+=\n");
    higflow_set_fonte_de_malha (ns, malha_t8_por_rank, NULL);
  } else if (strcmp (fonte, "t8code-criterio") == 0) {
    if (myrank == 0) printf ("=+=+=+= Malha do t8code por CRITERIO =+=+=+=\n");
    higflow_set_fonte_de_malha (ns, malha_t8_criterio, NULL);
  } else if (strcmp (fonte, "t8code-refinada") == 0) {
    if (myrank == 0) printf ("=+=+=+= Malha do t8code REFINADA em caixa =+=+=+=\n");
    higflow_set_fonte_de_malha (ns, malha_t8_refinada, NULL);
  } else if (strcmp (fonte, "t8code-particao") == 0) {
    if (myrank == 0) printf ("=+=+=+= Malha e particao do t8code =+=+=+=\n");
    higflow_set_fonte_de_particao (ns, particao_t8, NULL);
  } else {
    if (myrank == 0)
      fprintf (stderr, "HIGFLOW_MALHA=%s nao e' uma fonte conhecida.  Use "
                       "t8code, t8code-rank, t8code-refinada, t8code-criterio ou "
                       "t8code-particao.\n", fonte);
  }
}
