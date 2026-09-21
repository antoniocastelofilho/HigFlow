// FONTE DE MALHA DO T8CODE, para este exemplo.
//
// O exemplo escolhe de onde vem a malha: sem `HIGFLOW_MALHA=t8code` no ambiente,
// nada muda e o solver le' o arquivo AMR como sempre.  Com a variavel, a malha e'
// produzida por uma floresta do t8code.
//
// POR QUE ISTO PODE SER VERIFICADO SEM REFERENCIA NOVA.  O `amrs/domain/ch-d-0.amr`
// descreve uma malha UNIFORME de 160x40 sobre [0,8]x[-1,1].  Um cmesh de BRICK do
// t8code com 160x40 arvores, no nivel 0, tem exatamente essas celulas -- entao a
// arvore produzida aqui e' a MESMA, e a saida do exemplo tem de ser identica a'
// do caminho AMR.  A verificacao e' diferenca de arquivo, nao tolerancia contra
// numero gravado: gerar referencia a partir do codigo que se acabou de escrever
// so' guardaria contra regressao futura, nao provaria nada hoje.
//
// A ESTRUTURA VEM DA FLORESTA, nao de um `hig_refine_uniform` escrito aqui.  Cada
// folha do t8code e' lida e pedida a' arvore hig; se o t8code produzisse outra
// malha, a arvore sairia diferente.  E' o que impede que este modulo "passe" por
// reproduzir a malha do arquivo por conta propria.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "higtree.h"
#include "coord.h"
#include "t8-mesh-rank.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>

// A malha deste exemplo, a mesma do amrs/domain/ch-d-0.amr.
#define LOX 0.0
#define HIX 8.0
#define LOY (-1.0)
#define HIY 1.0
#define NCX 160
#define NCY 40

static int g_iniciado = 0;

// Produz a malha INTEIRA, em serie, em todo rank.  Simples e suficiente para
// mostrar que a fonte funciona, mas mantem o gargalo: cada processo materializa
// 6400 celulas antes de o `lbal` repartir.  Quem tira o gargalo e'
// `malha_t8_por_rank`, abaixo.
extern "C" int
malha_t8_do_exemplo (void *ctx, hig_cell **arvores, int max)
{
  (void) ctx;
  if (max < 1) return 0;
#if DIM != 2
  fprintf (stderr, "malha_t8_do_exemplo: so' existe em DIM=2\n");
  return 0;
#else
  if (!g_iniciado) {
    g_iniciado = 1;
    sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ERROR);
    t8_init (SC_LP_ERROR);
  }

  // O brick cobre [0,NCX] x [0,NCY], uma arvore por celula.  Em serie (COMM_SELF)
  // porque cada rank precisa da malha inteira para o `lb_calc_partition`, que e'
  // quem reparte -- a producao por rank e' outra frente.
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_brick_2d (cmesh, NCX, NCY, 0, 0, sc_MPI_COMM_SELF);
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_SELF);

  Point lo, hi;
  lo[0] = LOX; lo[1] = LOY;
  hi[0] = HIX; hi[1] = HIY;
  hig_cell *raiz = hig_create_root (lo, hi);
  int nc[DIM];
  nc[0] = NCX; nc[1] = NCY;
  hig_refine_uniform (raiz, nc);

  // Confere celula a celula contra a floresta: cada centroide do t8code tem de
  // cair numa celula da arvore cujo centro seja o mesmo.  Se o brick nao
  // reproduzir a malha do arquivo, isto acusa aqui e nao tres horas depois.
  long folhas = 0, ruins = 0;
  const t8_locidx_t ntrees = t8_forest_get_num_local_trees (f);
  for (t8_locidx_t it = 0; it < ntrees; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      Point p;
      p[0] = lo[0] + (c[0] / (double) NCX) * (hi[0] - lo[0]);
      p[1] = lo[1] + (c[1] / (double) NCY) * (hi[1] - lo[1]);
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

  if (folhas != (long) NCX * NCY || ruins != 0) {
    fprintf (stderr, "malha_t8_do_exemplo: a floresta tem %ld folha(s) e %ld "
                     "divergencia(s) contra a arvore (esperado %ld e 0)\n",
             folhas, ruins, (long) NCX * NCY);
    hig_destroy (raiz);
    return 0;
  }
  arvores[0] = raiz;
  return 1;
#endif
}

// ---------------------------------------------------------------------------
// PRODUCAO POR RANK: nenhum processo chega a ter a malha inteira.
//
// A floresta nasce em COMM_WORLD, o t8code a reparte, e cada rank materializa so'
// a sua parte -- em CAIXAS COMPLETAS, porque arvore com buraco nao e' navegavel.
// As caixas entram no `lb_add_input_tree` como entrada distribuida, que e' o que
// o `lbal` ja' espera; ele reparte de novo e monta o grafo de vizinhanca e a
// franja.  Por isso esta fonte nao produz franja: nao e' dela.
//
// A particao final e' a do `lbal`, nao a do t8code.  O que muda em relacao a'
// fonte em serie e' de onde vem a malha e quanto cada processo precisa segurar.
extern "C" int
malha_t8_por_rank (void *ctx, hig_cell **arvores, int max)
{
  (void) ctx;
#if DIM != 2
  fprintf (stderr, "malha_t8_por_rank: so' existe em DIM=2\n");
  return 0;
#else
  Point lo, hi;
  lo[0] = LOX; lo[1] = LOY;
  hi[0] = HIX; hi[1] = HIY;
  int nb[DIM];
  nb[0] = NCX; nb[1] = NCY;

  t8_producao_rank p;
  if (!t8_produz_por_rank_brick (lo, hi, nb, &p)) {
    fprintf (stderr, "malha_t8_por_rank: o produtor falhou\n");
    return 0;
  }
  if (p.base_dividida != 0) {
    fprintf (stderr, "malha_t8_por_rank: %ld celula(s) com familia dividida\n",
             p.base_dividida);
    t8_producao_rank_destroi (&p);
    return 0;
  }
  if (p.n_locais > max) {
    fprintf (stderr, "malha_t8_por_rank: %d caixas, cabe %d\n", p.n_locais, max);
    t8_producao_rank_destroi (&p);
    return 0;
  }
  // As arvores passam a pertencer ao `lbal` (`managed = true`), entao aqui so' se
  // solta o vetor -- destruir as arvores seria destruir a malha que se acabou de
  // entregar.
  for (int i = 0; i < p.n_locais; i++) arvores[i] = p.locais[i];
  const int n = p.n_locais;
  p.n_locais = 0;
  t8_producao_rank_destroi (&p);
  return n;
#endif
}
