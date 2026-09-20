// Adaptador: floresta do t8code -> arvore hig_cell.
//
// COMO ELE FUNCIONA, e por que nao e' uma reconstrucao em paralelo.
//
// A estrutura vem TODA do t8code.  A funcao constroi a floresta la', percorre as
// folhas dela, e para cada folha pede a' arvore hig que exista uma celula naquela
// posicao com aquele tamanho.  Se o t8code produzisse outra malha -- graduada, por
// exemplo -- a arvore hig sairia diferente, e e' isso que faz do teste uma
// comparacao e nao dois testes.  Nada aqui repete as decisoes de refino: elas sao
// lidas de volta do t8code, nivel a nivel.
//
// O CUSTO ESTA' A' VISTA, e e' o que a analise de estrategia antecipou: materializar
// a floresta linear como octree de ponteiros e' uma traducao, O(folhas) a cada
// producao.  Para responder ao contrato isso basta e e' honesto; para PRODUCAO, e'
// exatamente o custo que se queria evitar, e a saida sera' mover a fronteira das
// consultas, nao traduzir a cada passo.
//
// O t8code exige C++20 e a HiGTree constroi com gnu++17, entao esta unidade compila
// separada e so' a assinatura em `extern "C"` atravessa.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-mesh-producer.h"
#include "hig-mesh-snapshot.h"
#include "higtree.h"
#include "higtree-iterator.h"

#if DIM == 2
#define ECLASSE T8_ECLASS_QUAD
#else
#define ECLASSE T8_ECLASS_HEX
#endif

#define NIVEL_BASE 2            // 2^2 = 4 celulas por direcao, como o MTree

// Caixas de refino, nas mesmas posicoes do produtor de referencia.
static const double ALVO1 = 0.375;     // primeira refinada  -> h = 0,125
static const double ALVO2 = 0.3125;    // segunda, numa neta -> h = 0,0625

static int g_passo = 0;

// Refina a folha que contem o alvo do passo corrente.  Nenhum balanceamento e'
// pedido: `t8_forest_set_balance` NAO e' chamado, que e' o que permite o 4:1.
static int
adapt_alvo (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
            const t8_eclass_t tree_class, t8_locidx_t lelement_id,
            const t8_scheme_c *scheme, const int is_family, const int num_elements,
            t8_element_t *elements[])
{
  double c[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], c);

  const double alvo = (g_passo == 0) ? ALVO1 : ALVO2;
  const double meia = (g_passo == 0) ? 0.125 : 0.0625;   // meia largura da celula
  for (int d = 0; d < DIM; d++) {
    if (c[d] < alvo - meia || c[d] > alvo + meia) return 0;
  }
  return 1;
}

// Garante que a arvore hig tenha, na posicao `p`, uma celula de lado <= `h`.
// Refina 2 por direcao, um nivel de cada vez, como o t8code faz.
static void
refina_ate (hig_cell *raiz, const Point p, const double h)
{
  for (int guarda = 0; guarda < 32; guarda++) {
    hig_cell *c = hig_get_cell_with_point (raiz, p);
    if (c == NULL) return;                  // fora do dominio
    Point d;
    hig_get_delta (c, d);
    if (d[0] <= h * 1.0001) return;         // ja' esta' fino o bastante
    int nc[DIM];
    for (int k = 0; k < DIM; k++) nc[k] = 2;
    hig_refine_uniform (c, nc);
  }
}

// sc_init/t8_init sao por processo, nao por chamada.  O MPI ja' vem inicializado
// pelo higtree_initialize do teste -- ver o comentario no main dele.
static void
inicializa_uma_vez (void)
{
  static int feito = 0;
  if (feito) return;
  feito = 1;
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
  t8_init (SC_LP_ERROR);
}

// A floresta e' a mesma para os dois caminhos -- o que materializa em arvore e o
// que preenche o instantaneo direto.  Comum de proposito: se os dois montassem a
// floresta cada um a seu modo, compara-los nao diria nada.
static t8_forest_t
constroi_floresta (void)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);                    // obrigatorio ANTES do gerador
  t8_cmesh_new_hypercube (&cmesh, ECLASSE, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_scheme_c *scheme = t8_scheme_new_default ();

  t8_forest_t f =
    t8_forest_new_uniform (cmesh, scheme, NIVEL_BASE, 0, sc_MPI_COMM_WORLD);

  for (g_passo = 0; g_passo < 2; g_passo++) {
    t8_forest_t novo;
    t8_forest_init (&novo);
    t8_forest_set_adapt (novo, f, adapt_alvo, 0);
    t8_forest_set_ghost (novo, 1, T8_GHOST_FACES);
    t8_forest_commit (novo);
    f = novo;
  }
  return f;
}

// ---------------------------------------------------------------------------
// PRODUCAO: a malha que o dominio do solver vai usar sai daqui.
//
// A diferenca para `t8_produz_malha_nao_graduada` nao e' de mecanismo -- e' de
// PAPEL.  Aquela existe para responder a clausula C11 com uma malha decorada; esta
// recebe o nivel base, o alvo e quantas vezes refinar, e devolve a arvore que o
// `lb_add_input_tree` vai distribuir.  Quem decide a estrutura e' o t8code; a
// arvore hig so' a reproduz.
//
// A CAIXA CONTINUA [0,1]^DIM, e isso e' limite declarado: o hipercubo do t8code e'
// unitario, e mapear para uma caixa qualquer exige escalar tambem o criterio de
// parada do refinamento, que hoje compara `delta[0]`.  Enquanto a producao for
// verificada contra o contrato -- e nao contra a geometria de um exemplo -- a
// caixa unitaria basta e evita uma conversao sem teste.
static double g_alvo_prod = 0.375;
static double g_meia_prod = 0.125;
static int    g_refinos_prod = 0;

static int
adapt_prod (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
            const t8_eclass_t tree_class, t8_locidx_t lelement_id,
            const t8_scheme_c *scheme, const int is_family, const int num_elements,
            t8_element_t *elements[])
{
  double c[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], c);
  for (int d = 0; d < DIM; d++) {
    if (c[d] < g_alvo_prod - g_meia_prod || c[d] > g_alvo_prod + g_meia_prod)
      return 0;
  }
  return 1;
}

extern "C" hig_cell *
t8_produz_malha_para_dominio (int nivel_base, double alvo, int refinos,
                              long *folhas_out)
{
  if (folhas_out != NULL) *folhas_out = 0;
  if (nivel_base < 1 || refinos < 0) return NULL;
  inicializa_uma_vez ();

  // COMM_SELF NOS DOIS, cmesh e floresta.  Com o cmesh em COMM_WORLD ele sai
  // DISTRIBUIDO, e o rank 0 -- unico que produz -- enxerga so' a parte dele: a
  // malha saia com 47 folhas em np=2 e 36 em np=3, em vez das 79.  Aqui o t8code
  // produz a malha INTEIRA em serie; quem a reparte depois e' o `lbal`.
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, ECLASSE, sc_MPI_COMM_SELF, 0, 0, 0);
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, nivel_base, 0,
                                         sc_MPI_COMM_SELF);

  g_alvo_prod = alvo;
  g_meia_prod = 0.5 / (double) (1 << nivel_base);
  for (int r = 0; r < refinos; r++) {
    t8_forest_t novo;
    t8_forest_init (&novo);
    t8_forest_set_adapt (novo, f, adapt_prod, 0);
    t8_forest_commit (novo);
    f = novo;
    g_meia_prod *= 0.5;
  }

  Point lo, hi;
  for (int d = 0; d < DIM; d++) { lo[d] = 0.0; hi[d] = 1.0; }
  hig_cell *raiz = hig_create_root (lo, hi);
  int nc[DIM];
  for (int d = 0; d < DIM; d++) nc[d] = 1 << nivel_base;
  hig_refine_uniform (raiz, nc);

  const t8_scheme_c *sch = t8_forest_get_scheme (f);
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  long folhas = 0;
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      const int nivel = sch->element_get_level (ec, e);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      Point p;
      for (int d = 0; d < DIM; d++) p[d] = c[d];
      refina_ate (raiz, p, 1.0 / (double) (1 << nivel));
      folhas++;
    }
  }
  t8_forest_unref (&f);

  if (folhas == 0) { hig_destroy (raiz); return NULL; }
  if (folhas_out != NULL) *folhas_out = folhas;
  return raiz;
}

extern "C" hig_cell *
t8_produz_malha_nao_graduada (void)
{
  inicializa_uma_vez ();

  t8_forest_t f = constroi_floresta ();
  if (f == NULL) return NULL;
  const t8_scheme_c *scheme = t8_forest_get_scheme (f);

  // ------------------------------------------------ materializacao
  Point lo, hi;
  for (int d = 0; d < DIM; d++) { lo[d] = 0.0; hi[d] = 1.0; }
  hig_cell *raiz = hig_create_root (lo, hi);
  int nc[DIM];
  for (int d = 0; d < DIM; d++) nc[d] = 1 << NIVEL_BASE;
  hig_refine_uniform (raiz, nc);             // o nivel base, de uma vez

  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  long folhas = 0;
  int nivel_max = NIVEL_BASE;
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      const int nivel = scheme->element_get_level (ec, e);
      if (nivel > nivel_max) nivel_max = nivel;

      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      Point p;
      for (int d = 0; d < DIM; d++) p[d] = c[d];

      // lado da folha do t8code: o hipercubo tem lado 1 e cada nivel divide por 2
      refina_ate (raiz, p, 1.0 / (double) (1 << nivel));
      folhas++;
    }
  }

  t8_forest_unref (&f);

  // Guarda: se a floresta nao passou do nivel base, nao ha' salto nenhum e
  // devolver a malha seria entregar um uniforme fingindo de nao graduado.
  if (nivel_max < NIVEL_BASE + 2 || folhas == 0) {
    hig_destroy (raiz);
    return NULL;
  }
  return raiz;
}

// ---------------------------------------------------------------------------
// A FRONTEIRA MOVIDA: instantaneo preenchido DIRETO da floresta.
//
// Compare com t8_produz_malha_nao_graduada acima.  Aquela constroi a floresta e
// depois refina uma arvore hig folha a folha para reproduzi-la -- traducao, a
// cada producao.  Esta le as mesmas folhas e escreve centro e tamanho num
// arranjo plano.  Nao ha' octree de ponteiros em lugar nenhum, e e' esse o
// ponto: o backend passa a entregar o que as consultas precisam, em vez de
// entregar uma estrutura que elas saibam navegar.
// ---------------------------------------------------------------------------
extern "C" struct hig_mesh_snapshot *
t8_preenche_instantaneo (void)
{
  inicializa_uma_vez ();

  t8_forest_t f = constroi_floresta ();
  if (f == NULL) return NULL;

  const t8_scheme_c *scheme = t8_forest_get_scheme (f);
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);

  long n = 0;
  int nivel_max = NIVEL_BASE;
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    n += t8_forest_get_tree_num_leaf_elements (f, it);
  }

  hig_mesh_snapshot *s = hms_create ((int) n);
  if (s == NULL) { t8_forest_unref (&f); return NULL; }

  long i = 0;
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++, i++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      const int nivel = scheme->element_get_level (ec, e);
      if (nivel > nivel_max) nivel_max = nivel;

      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      const double lado = 1.0 / (double) (1 << nivel);
      for (int d = 0; d < DIM; d++) {
          // O instantaneo guarda A CAIXA.  O t8code entrega centroide e nivel,
          // entao a caixa sai dai' -- o backend se adapta ao que o consumidor
          // le', que e' o sentido certo da dependencia.
          s->low[i * DIM + d]  = c[d] - lado / 2.0;
          s->high[i * DIM + d] = c[d] + lado / 2.0;
      }
    }
  }

  t8_forest_unref (&f);

  if (nivel_max < NIVEL_BASE + 2) { hms_destroy (s); return NULL; }
  return s;
}
