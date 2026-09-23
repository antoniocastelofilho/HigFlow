// Ver o .h para o porque das caixas completas.

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "t8-mesh-rank.h"
#include "t8-forest-cache.hxx"

#include <cmath>
#include <cstdlib>
#include <cstring>

#if DIM == 2
#define ECLASSE T8_ECLASS_QUAD
#else
#define ECLASSE T8_ECLASS_HEX
#endif

static double g_alvo[DIM];
static double g_meia = 0.0;

static int
adapt_rank (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
            const t8_eclass_t tree_class, t8_locidx_t lelement_id,
            const t8_scheme_c *scheme, const int is_family, const int num_elements,
            t8_element_t *elements[])
{
  double c[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], c);
  for (int d = 0; d < DIM; d++)
    if (c[d] < g_alvo[d] - g_meia || c[d] > g_alvo[d] + g_meia) return 0;
  return 1;
}

struct Folha { Point x; int nivel; };

// Materializa o interior de UMA celula, completo.  `fs` sao as folhas da floresta
// que caem dentro dela; `nivel` e' o nivel dessa celula.
//
// Se a familia estiver dividida entre ranks, algum neto ficara' sem folha e a
// contagem no fim acusa -- por isso nada aqui tenta adivinhar o que falta.
static void
materializa_completo (hig_cell *c, Folha *fs, long n, int nivel, long *dividida)
{
  if (n <= 0) { if (dividida) (*dividida)++; return; }
  int so_neste_nivel = 1;
  for (long k = 0; k < n; k++) if (fs[k].nivel > nivel) { so_neste_nivel = 0; break; }
  if (so_neste_nivel) return;                 // a celula E' a folha

  int dois[DIM];
  for (int d = 0; d < DIM; d++) dois[d] = 2;
  hig_refine_uniform (c, dois);               // COMPLETO: cria todos os filhos

  const int nf = hig_get_number_of_children (c);
  Folha *buf = (Folha *) malloc ((size_t) n * sizeof *buf);
  for (int i = 0; i < nf; i++) {
    hig_cell *f = hig_get_child (c, i);
    Point flo, fhi;
    hig_get_lowpoint (f, flo);
    hig_get_highpoint (f, fhi);
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < flo[d] || fs[k].x[d] > fhi[d]) dentro = 0;
      if (dentro) buf[m++] = fs[k];
    }
    materializa_completo (f, buf, m, nivel + 1, dividida);
  }
  free (buf);
}

// Agrupa celulas base possuidas em retangulos maximais.  Varredura por linhas: a
// primeira direcao cresce enquanto houver celula possuida, e a linha e' estendida
// nas demais direcoes enquanto a faixa inteira estiver possuida.  Nao e' a
// decomposicao minima -- e' simples e correta, e o numero de caixas so' importa
// para a contagem de arvores do dominio.
struct Caixa { int i0[DIM], i1[DIM]; };

static int
agrupa (const char *posse, const int nb[DIM], Caixa *cx, int max)
{
  long total = 1;
  for (int d = 0; d < DIM; d++) total *= nb[d];
  char *livre = (char *) malloc ((size_t) total);
  memcpy (livre, posse, (size_t) total);

  // indice linear a partir das coordenadas, na mesma ordem usada para preencher
  // `posse`: a direcao 0 varia mais rapido.
  #define LIN(p) ({ long _s = 0, _m = 1; \
                    for (int _d = 0; _d < DIM; _d++) { _s += (p)[_d] * _m; _m *= nb[_d]; } _s; })

  int n = 0;
  for (long lin = 0; lin < total; lin++) {
    if (!livre[lin]) continue;
    int idx[DIM];
    { long r = lin; for (int d = 0; d < DIM; d++) { idx[d] = (int) (r % nb[d]); r /= nb[d]; } }

    Caixa c;
    for (int d = 0; d < DIM; d++) { c.i0[d] = idx[d]; c.i1[d] = idx[d]; }

    // Cresce em cada direcao enquanto a FACE INTEIRA estiver livre.  A versao
    // anterior conferia so' uma linha da face, e em 3D isso fazia a caixa engolir
    // celulas que nao sao deste rank -- a materializacao as criava vazias, e o
    // dominio saia com celula a mais.
    int cresceu = 1;
    while (cresceu) {
      cresceu = 0;
      for (int d = 0; d < DIM; d++) {
        if (c.i1[d] + 1 >= nb[d]) continue;
        int p[DIM];
        for (int q = 0; q < DIM; q++) p[q] = c.i0[q];
        p[d] = c.i1[d] + 1;
        int ok = 1;
        // percorre a face inteira: todas as direcoes menos `d`
        while (ok) {
          if (!livre[LIN(p)]) { ok = 0; break; }
          int q = 0;
          while (q < DIM) {
            if (q == d) { q++; continue; }
            if (p[q] < c.i1[q]) { p[q]++; break; }
            p[q] = c.i0[q]; q++;
          }
          if (q == DIM) break;            // face inteira percorrida
        }
        if (ok) { c.i1[d]++; cresceu = 1; }
      }
    }

    // marca a caixa como consumida
    {
      int p[DIM];
      for (int q = 0; q < DIM; q++) p[q] = c.i0[q];
      while (1) {
        livre[LIN(p)] = 0;
        int q = 0;
        while (q < DIM) {
          if (p[q] < c.i1[q]) { p[q]++; break; }
          p[q] = c.i0[q]; q++;
        }
        if (q == DIM) break;
      }
    }
    if (n < max) cx[n++] = c;
  }
  #undef LIN
  free (livre);
  return n;
}

// A celula `box` e' COMPLETAMENTE possuida pelas folhas `fs`?
//
// MEDIDO, e por isso este predicado existe: `set_for_coarsening=1` NAO garante
// posse completa por celula base.  Ele preserva a familia do nivel mais fino, e
// uma celula base refinada duas vezes contem familia de familias -- em 3D com
// np=2, oito celulas base sairam divididas.  Decompor no nivel base e' portanto
// invalido; a caixa tem de ser emitida no nivel em que a posse E' completa.
static int
completo (const Point blo, const Point bhi, int nivel, Folha *fs, long n)
{
  if (n <= 0) return 0;
  int so_neste = 1;
  for (long k = 0; k < n; k++) if (fs[k].nivel > nivel) { so_neste = 0; break; }
  if (so_neste) return (n == 1);

  Point meio;
  for (int d = 0; d < DIM; d++) meio[d] = 0.5 * (blo[d] + bhi[d]);
  const int nf = 1 << DIM;
  Folha *buf = (Folha *) malloc ((size_t) n * sizeof *buf);
  int ok = 1;
  for (int i = 0; i < nf && ok; i++) {
    Point clo, chi;
    for (int d = 0; d < DIM; d++) {
      const int alto = (i >> d) & 1;
      clo[d] = alto ? meio[d] : blo[d];
      chi[d] = alto ? bhi[d]  : meio[d];
    }
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < clo[d] || fs[k].x[d] > chi[d]) dentro = 0;
      if (dentro) buf[m++] = fs[k];
    }
    if (!completo (clo, chi, nivel + 1, buf, m)) ok = 0;
  }
  free (buf);
  return ok;
}

// Emite caixas completas a partir de `box`: se ela ja' e' completa, sai uma
// arvore; se nao, desce e emite as sub-caixas que forem.
static void
emite_completas (const Point blo, const Point bhi, int nivel, Folha *fs, long n,
                 hig_cell **saida, int max, int *out, long *dividida)
{
  if (n <= 0 || *out >= max) return;
  if (completo (blo, bhi, nivel, fs, n)) {
    hig_cell *raiz = hig_create_root ((real *) blo, (real *) bhi);
    materializa_completo (raiz, fs, n, nivel, dividida);
    saida[(*out)++] = raiz;
    return;
  }
  Point meio;
  for (int d = 0; d < DIM; d++) meio[d] = 0.5 * (blo[d] + bhi[d]);
  const int nf = 1 << DIM;
  Folha *buf = (Folha *) malloc ((size_t) n * sizeof *buf);
  for (int i = 0; i < nf; i++) {
    Point clo, chi;
    for (int d = 0; d < DIM; d++) {
      const int alto = (i >> d) & 1;
      clo[d] = alto ? meio[d] : blo[d];
      chi[d] = alto ? bhi[d]  : meio[d];
    }
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < clo[d] || fs[k].x[d] > chi[d]) dentro = 0;
      if (dentro) buf[m++] = fs[k];
    }
    emite_completas (clo, chi, nivel + 1, buf, m, saida, max, out, dividida);
  }
  free (buf);
}

// Uma caixa -> uma arvore COMPLETA.
static hig_cell *
monta_caixa (const Point lo, const Point h, const Caixa *c, int nivel_base,
             Folha *fs, long n, long *dividida)
{
  Point rlo, rhi;
  int ext[DIM];
  for (int d = 0; d < DIM; d++) {
    rlo[d] = lo[d] + c->i0[d] * h[d];
    rhi[d] = lo[d] + (c->i1[d] + 1) * h[d];
    ext[d] = c->i1[d] - c->i0[d] + 1;
  }
  hig_cell *raiz = hig_create_root (rlo, rhi);
  hig_refine_uniform (raiz, ext);            // COMPLETO: sem buraco

  const int nf = hig_get_number_of_children (raiz);
  Folha *buf = (Folha *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *buf);
  for (int i = 0; i < nf; i++) {
    hig_cell *f = hig_get_child (raiz, i);
    Point flo, fhi;
    hig_get_lowpoint (f, flo);
    hig_get_highpoint (f, fhi);
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < flo[d] || fs[k].x[d] > fhi[d]) dentro = 0;
      if (dentro) buf[m++] = fs[k];
    }
    materializa_completo (f, buf, m, nivel_base, dividida);
  }
  free (buf);
  return raiz;
}

// `nb` e' a grade BASE (quantas celulas por direcao no nivel mais grosso) e
// `nivel_base` e' o nivel que os elementos dessa grade tem na floresta.  Para o
// hipercubo sao `1<<L` e `L`; para um BRICK de nb[0] x nb[1] arvores sao a
// propria grade e o nivel 0.  Separar os dois e' o que permite a mesma
// decomposicao servir aos dois produtores.
static int
monta_conjunto (const Point lo, const Point hi, const int nb[DIM], int nivel_base,
                Folha *fs, long n, hig_cell **saida, int max, long *dividida)
{
  if (n <= 0) return 0;
  Point h;
  for (int d = 0; d < DIM; d++) h[d] = (hi[d] - lo[d]) / (double) nb[d];

  long total = 1;
  for (int d = 0; d < DIM; d++) total *= nb[d];
  char *posse = (char *) calloc ((size_t) total, 1);
  for (long k = 0; k < n; k++) {
    long pos = 0, mul = 1;
    for (int d = 0; d < DIM; d++) {
      int i = (int) floor ((fs[k].x[d] - lo[d]) / h[d]);
      if (i < 0) i = 0;
      if (i >= nb[d]) i = nb[d] - 1;
      pos += i * mul; mul *= nb[d];
    }
    posse[pos] = 1;
  }

  // So' entram no agrupamento as celulas base COMPLETAS.  As incompletas descem
  // e viram caixas menores, no nivel em que a posse fecha.
  Folha *sel = (Folha *) malloc ((size_t) n * sizeof *sel);
  for (long pos = 0; pos < total; pos++) {
    if (!posse[pos]) continue;
    long r = pos; int ix[DIM];
    for (int d = 0; d < DIM; d++) { ix[d] = (int) (r % nb[d]); r /= nb[d]; }
    Point clo, chi;
    for (int d = 0; d < DIM; d++) {
      clo[d] = lo[d] + ix[d] * h[d];
      chi[d] = lo[d] + (ix[d] + 1) * h[d];
    }
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < clo[d] || fs[k].x[d] > chi[d]) dentro = 0;
      if (dentro) sel[m++] = fs[k];
    }
    if (!completo (clo, chi, nivel_base, sel, m)) posse[pos] = 2;   // incompleta
  }

  Caixa cx[T8_MAX_CAIXAS];
  char *so_completas = (char *) malloc ((size_t) total);
  for (long q = 0; q < total; q++) so_completas[q] = (posse[q] == 1) ? 1 : 0;
  const int ncx = agrupa (so_completas, nb, cx, T8_MAX_CAIXAS);
  free (so_completas);

  int out = 0;
  Folha *buf = (Folha *) malloc ((size_t) n * sizeof *buf);
  for (int c = 0; c < ncx && out < max; c++) {
    Point clo, chi;
    for (int d = 0; d < DIM; d++) {
      clo[d] = lo[d] + cx[c].i0[d] * h[d];
      chi[d] = lo[d] + (cx[c].i1[d] + 1) * h[d];
    }
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < clo[d] || fs[k].x[d] > chi[d]) dentro = 0;
      if (dentro) buf[m++] = fs[k];
    }
    saida[out++] = monta_caixa (lo, h, &cx[c], nivel_base, buf, m, dividida);
  }
  // As celulas base incompletas, uma a uma, descendo ate' onde a posse fecha.
  for (long pos = 0; pos < total && out < max; pos++) {
    if (posse[pos] != 2) continue;
    long r = pos; int ix[DIM];
    for (int d = 0; d < DIM; d++) { ix[d] = (int) (r % nb[d]); r /= nb[d]; }
    Point clo, chi;
    for (int d = 0; d < DIM; d++) {
      clo[d] = lo[d] + ix[d] * h[d];
      chi[d] = lo[d] + (ix[d] + 1) * h[d];
    }
    long m = 0;
    for (long k = 0; k < n; k++) {
      int dentro = 1;
      for (int d = 0; d < DIM && dentro; d++)
        if (fs[k].x[d] < clo[d] || fs[k].x[d] > chi[d]) dentro = 0;
      if (dentro) sel[m++] = fs[k];
    }
    emite_completas (clo, chi, nivel_base, sel, m, saida, max, &out, dividida);
  }
  free (sel);
  free (posse);
  free (buf);
  return out;
}

extern "C" int
t8_produz_por_rank (const Point lo, const Point hi, int nivel_base,
                    const Point alvo, int refinos, t8_producao_rank *out)
{
  if (out == NULL) return 0;
  memset (out, 0, sizeof *out);
  if (nivel_base < 1 || refinos < 0) return 0;
  for (int d = 0; d < DIM; d++) if (!(hi[d] > lo[d])) return 0;

  t8_inicializa_uma_vez ();

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, ECLASSE, sc_MPI_COMM_WORLD, 0, 0, 0);
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, nivel_base, 0,
                                         sc_MPI_COMM_WORLD);

  for (int d = 0; d < DIM; d++) g_alvo[d] = (alvo[d] - lo[d]) / (hi[d] - lo[d]);
  g_meia = 0.5 / (double) (1 << nivel_base);

  // Refino e distribuicao em COMMITS SEPARADOS: adapt, ghost e partition no mesmo
  // commit derruba com np>1 e funciona em np=1.
  for (int r = 0; r < refinos; r++) {
    t8_forest_t novo;
    t8_forest_init (&novo);
    t8_forest_set_adapt (novo, f, adapt_rank, 0);
    t8_forest_commit (novo);
    f = novo;
    g_meia *= 0.5;
  }
  {
    t8_forest_t fp;
    t8_forest_init (&fp);
    // set_for_coarsening = 1: pede ao t8code que NAO divida familia entre ranks.
    // E' o que torna a caixa completa uma representacao fiel; o teste confere.
    t8_forest_set_partition (fp, f, 1);
    t8_forest_set_ghost (fp, 1, T8_GHOST_FACES);
    t8_forest_commit (fp);
    f = fp;
  }

  out->n_global = (long) t8_forest_get_global_num_leaf_elements (f);
  const t8_scheme_c *sch = t8_forest_get_scheme (f);
  const t8_locidx_t nloc_trees = t8_forest_get_num_local_trees (f);

  // ---------------------------------------------------------------- locais
  long n = (long) t8_forest_get_local_num_leaf_elements (f);
  Folha *loc = (Folha *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *loc);
  long k = 0;
  for (t8_locidx_t it = 0; it < nloc_trees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      for (int d = 0; d < DIM; d++) loc[k].x[d] = lo[d] + c[d] * (hi[d] - lo[d]);
      loc[k].nivel = sch->element_get_level (ec, e);
      k++;
    }
  }
  out->n_local = k;
  int nb[DIM];
  for (int d = 0; d < DIM; d++) nb[d] = 1 << nivel_base;
  out->n_locais = monta_conjunto (lo, hi, nb, nivel_base, loc, k, out->locais,
                                  T8_MAX_CAIXAS, &out->base_dividida);
  free (loc);

  // ---------------------------------------------------------------- franja
  const t8_locidx_t ngh = t8_forest_ghost_num_trees (f);
  long ng = (long) t8_forest_get_num_ghosts (f);
  Folha *gh = (Folha *) malloc ((size_t) (ng > 0 ? ng : 1) * sizeof *gh);
  long kg = 0;
  for (t8_locidx_t gt = 0; gt < ngh; gt++) {
    const t8_eclass_t ec = t8_forest_ghost_get_tree_class (f, gt);
    const t8_locidx_t ne = t8_forest_ghost_tree_num_leaf_elements (f, gt);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_ghost_get_leaf_element (f, gt, ie);
      double c[3];
      t8_forest_element_centroid (f, nloc_trees + gt, e, c);
      for (int d = 0; d < DIM; d++) gh[kg].x[d] = lo[d] + c[d] * (hi[d] - lo[d]);
      gh[kg].nivel = sch->element_get_level (ec, e);
      kg++;
    }
  }
  out->n_franja = kg;
  // A franja NAO conta para `base_dividida`: ela e' um recorte da malha do
  // vizinho por construcao, e familia dividida ali e' esperada.
  long descartado = 0;
  out->n_franjas = monta_conjunto (lo, hi, nb, nivel_base, gh, kg, out->franjas,
                                   T8_MAX_CAIXAS, &descartado);
  free (gh);

  t8_forest_unref (&f);
  return 1;
}

// A MESMA decomposicao, alimentada por um cmesh de BRICK.  Serve para malha
// uniforme de dimensoes quaisquer -- 160x40, por exemplo --, que o hipercubo nao
// representa: ele e' uma arvore so', e o refino uniforme dele da' potencia de
// dois por direcao.
extern "C" int
t8_produz_por_rank_brick (const Point lo, const Point hi, const int nb[DIM],
                          t8_producao_rank *out)
{
  if (out == NULL) return 0;
  memset (out, 0, sizeof *out);
  for (int d = 0; d < DIM; d++) if (nb[d] < 1 || !(hi[d] > lo[d])) return 0;

  t8_inicializa_uma_vez ();

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
#if DIM == 2
  t8_cmesh_new_brick_2d (cmesh, nb[0], nb[1], 0, 0, sc_MPI_COMM_WORLD);
#else
  t8_cmesh_new_brick_3d (cmesh, nb[0], nb[1], nb[2], 0, 0, 0, sc_MPI_COMM_WORLD);
#endif
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, 0, 1, sc_MPI_COMM_WORLD);

  out->n_global = (long) t8_forest_get_global_num_leaf_elements (f);
  const t8_scheme_c *sch = t8_forest_get_scheme (f);
  const t8_locidx_t nloc_trees = t8_forest_get_num_local_trees (f);

  // O brick cobre [0,nb[0]] x [0,nb[1]] (x [0,nb[2]]); a volta para o dominio
  // divide pela extensao de cada direcao, e nao por um lado unitario.
  long n = (long) t8_forest_get_local_num_leaf_elements (f);
  Folha *loc = (Folha *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *loc);
  long k = 0;
  for (t8_locidx_t it = 0; it < nloc_trees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      for (int d = 0; d < DIM; d++)
        loc[k].x[d] = lo[d] + (c[d] / (double) nb[d]) * (hi[d] - lo[d]);
      loc[k].nivel = sch->element_get_level (ec, e);
      k++;
    }
  }
  out->n_local = k;
  out->n_locais = monta_conjunto (lo, hi, nb, 0, loc, k, out->locais,
                                  T8_MAX_CAIXAS, &out->base_dividida);
  free (loc);
  t8_forest_unref (&f);
  return 1;
}

// Caixa de refino, em coordenadas do BRICK.  O brick cobre [0,nb[d]] e nao o
// dominio, entao a caixa do usuario e' convertida uma vez, aqui, em vez de o
// callback converter a cada elemento.
static double g_cx_lo[DIM], g_cx_hi[DIM];

static int
adapt_caixa (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
             const t8_eclass_t tree_class, t8_locidx_t lelement_id,
             const t8_scheme_c *scheme, const int is_family, const int num_elements,
             t8_element_t *elements[])
{
  double c[3];
  t8_forest_element_centroid (forest_from, which_tree, elements[0], c);
  for (int d = 0; d < DIM; d++)
    if (c[d] < g_cx_lo[d] || c[d] > g_cx_hi[d]) return 0;
  return 1;
}

// Brick com REFINO numa caixa.  Duplica o corpo de `t8_produz_por_rank_brick`
// em vez de parametriza-lo: o caminho uniforme e' usado por exemplos que ja'
// tem referencia gravada, e nao vale arriscar mudanca de comportamento neles
// por economia de vinte linhas.
//
// NAO pede balanceamento (`t8_forest_set_balance` nao e' chamado), entao a razao
// entre celulas vizinhas pode chegar a 4:1.  Para a fronteira imersa isso
// importa: o corpo tem de ficar INTEIRAMENTE dentro do nivel fino, com folga
// maior que o suporte do nucleo, e a guarda do modulo verifica isso.
extern "C" int
t8_produz_por_rank_brick_refinado (const Point lo, const Point hi, const int nb[DIM],
                                   const Point caixa_lo, const Point caixa_hi,
                                   int refinos, t8_producao_rank *out)
{
  if (out == NULL) return 0;
  memset (out, 0, sizeof *out);
  for (int d = 0; d < DIM; d++) if (nb[d] < 1 || !(hi[d] > lo[d])) return 0;
  if (refinos < 0) return 0;

  // dominio -> brick
  for (int d = 0; d < DIM; d++) {
    g_cx_lo[d] = (caixa_lo[d] - lo[d]) / (hi[d] - lo[d]) * (double) nb[d];
    g_cx_hi[d] = (caixa_hi[d] - lo[d]) / (hi[d] - lo[d]) * (double) nb[d];
    if (!(g_cx_hi[d] > g_cx_lo[d])) return 0;
  }

  t8_inicializa_uma_vez ();

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
#if DIM == 2
  t8_cmesh_new_brick_2d (cmesh, nb[0], nb[1], 0, 0, sc_MPI_COMM_WORLD);
#else
  t8_cmesh_new_brick_3d (cmesh, nb[0], nb[1], nb[2], 0, 0, 0, sc_MPI_COMM_WORLD);
#endif
  const t8_scheme_c *scheme = t8_scheme_new_default ();
  t8_forest_t f = t8_forest_new_uniform (cmesh, scheme, 0, 1, sc_MPI_COMM_WORLD);

  // Um nivel de cada vez: o t8code refina UMA vez por passada de adapt.
  //
  // BALANCEAMENTO PEDIDO, ao contrario do resto deste arquivo.  Os produtores de
  // teste omitem `t8_forest_set_balance` de proposito -- e' o que permite exibir
  // o 4:1 que o contrato de Mesh tolera.  Para um SOLVER isso nao serve: com
  // salto 4:1 na fronteira de refino, o escoamento diverge.
  //
  // ATENCAO: pedir balanceamento NAO consertou a divergencia da caixa grande --
  // com UM nivel de refino o salto ja' e' no maximo 2:1 por construcao, e os
  // numeros da corrida saem IDENTICOS ate' o ultimo digito com e sem a chamada.
  // Fica porque para `refinos > 1` ela passa a importar, e porque um solver nao
  // deve depender de o usuario lembrar de pedir; mas nao foi a cura de nada.
  for (int r = 0; r < refinos; r++) {
    t8_forest_t novo;
    t8_forest_init (&novo);
    t8_forest_set_adapt (novo, f, adapt_caixa, 0);
    t8_forest_set_balance (novo, NULL, 0);
    t8_forest_set_partition (novo, NULL, 0);
    t8_forest_commit (novo);
    f = novo;
  }

  out->n_global = (long) t8_forest_get_global_num_leaf_elements (f);
  const t8_scheme_c *sch = t8_forest_get_scheme (f);
  const t8_locidx_t nloc_trees = t8_forest_get_num_local_trees (f);

  long n = (long) t8_forest_get_local_num_leaf_elements (f);
  Folha *loc = (Folha *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *loc);
  long k = 0;
  for (t8_locidx_t it = 0; it < nloc_trees; it++) {
    const t8_eclass_t ec = t8_forest_get_tree_class (f, it);
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements (f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree (f, it, ie);
      double c[3];
      t8_forest_element_centroid (f, it, e, c);
      for (int d = 0; d < DIM; d++)
        loc[k].x[d] = lo[d] + (c[d] / (double) nb[d]) * (hi[d] - lo[d]);
      loc[k].nivel = sch->element_get_level (ec, e);
      k++;
    }
  }
  out->n_local = k;
  out->n_locais = monta_conjunto (lo, hi, nb, 0, loc, k, out->locais,
                                  T8_MAX_CAIXAS, &out->base_dividida);
  free (loc);
  t8_forest_unref (&f);
  return 1;
}

extern "C" void
t8_producao_rank_destroi (t8_producao_rank *p)
{
  if (p == NULL) return;
  for (int i = 0; i < p->n_locais;  i++) if (p->locais[i])  hig_destroy (p->locais[i]);
  for (int i = 0; i < p->n_franjas; i++) if (p->franjas[i]) hig_destroy (p->franjas[i]);
  p->n_locais = p->n_franjas = 0;
}
