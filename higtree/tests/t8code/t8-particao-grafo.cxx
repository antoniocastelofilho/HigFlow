// Ver o .h para o contrato.  Nada aqui e' negociado por mensagem: os dois lados
// calculam o mesmo objeto a partir dos mesmos dados globais.

#include <mpi.h>
#include <glib.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "t8-particao-grafo.h"
#include "t8-mesh-rank.h"
#include "higtree-iterator.h"
#include "utils.h"

struct Caixa { int lo[DIM], hi[DIM]; };   // [lo, hi) em indices da grade base

static int
vazia (const Caixa *c)
{
  for (int d = 0; d < DIM; d++) if (c->hi[d] <= c->lo[d]) return 1;
  return 0;
}

// Interseccao de `a` com `b` dilatada de `f` celulas.  E' o que o dono de `b`
// precisa ver de `a`.
static Caixa
fatia (const Caixa *a, const Caixa *b, int f, const int nb[DIM])
{
  Caixa r;
  for (int d = 0; d < DIM; d++) {
    int lo = b->lo[d] - f; if (lo < 0) lo = 0;
    int hi = b->hi[d] + f; if (hi > nb[d]) hi = nb[d];
    r.lo[d] = (a->lo[d] > lo) ? a->lo[d] : lo;
    r.hi[d] = (a->hi[d] < hi) ? a->hi[d] : hi;
  }
  return r;
}

extern "C" int
t8_monta_dominio_particionado (const Point lo, const Point hi, const int nb[DIM],
                               sim_domain *sd, partition_graph *pg)
{
  int rank, np;
  MPI_Comm comm = pg_get_MPI_comm (pg);
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &np);
  const int f = (int) pg_get_fringe_size (pg);

  Point h;
  for (int d = 0; d < DIM; d++) h[d] = (hi[d] - lo[d]) / (double) nb[d];

  // ------------------------------------------------ as caixas deste rank
  t8_producao_rank p;
  if (!t8_produz_por_rank_brick (lo, hi, nb, &p)) {
    fprintf (stderr, "t8_monta_dominio_particionado: o produtor falhou\n");
    return 0;
  }
  if (p.base_dividida != 0) {
    fprintf (stderr, "t8_monta_dominio_particionado: %ld celula(s) com familia "
                     "dividida\n", p.base_dividida);
    t8_producao_rank_destroi (&p);
    return 0;
  }

  const int nloc = p.n_locais;
  Caixa *minhas = (Caixa *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *minhas);
  for (int i = 0; i < nloc; i++) {
    Point tlo, thi;
    hig_get_lowpoint (p.locais[i], tlo);
    hig_get_highpoint (p.locais[i], thi);
    for (int d = 0; d < DIM; d++) {
      minhas[i].lo[d] = (int) llround ((tlo[d] - lo[d]) / h[d]);
      minhas[i].hi[d] = (int) llround ((thi[d] - lo[d]) / h[d]);
    }
    sd_add_higtree (sd, p.locais[i]);
    struct tree_properties *tp =
        (struct tree_properties *) malloc (sizeof *tp);
    tp->partition_group = 0;
    tp->is_fringe = false;
    g_hash_table_insert (pg->tree_props, p.locais[i], tp);
  }
  p.n_locais = 0;                    // as arvores passaram a ser do dominio
  t8_producao_rank_destroi (&p);

  // ------------------------------------------- as caixas de TODOS os ranks
  // Allgather: cada rank precisa das caixas dos outros para calcular, sozinho, a
  // mesma faixa que o outro vai calcular.  E' isso que dispensa negociar ordem.
  int *contas = (int *) malloc ((size_t) np * sizeof *contas);
  MPI_Allgather (&nloc, 1, MPI_INT, contas, 1, MPI_INT, comm);
  int *desl = (int *) malloc ((size_t) np * sizeof *desl);
  int total = 0;
  for (int r = 0; r < np; r++) { desl[r] = total; total += contas[r]; }

  int *envio = (int *) malloc ((size_t) (nloc > 0 ? nloc : 1) * 2 * DIM * sizeof *envio);
  for (int i = 0; i < nloc; i++)
    for (int d = 0; d < DIM; d++) {
      envio[i * 2 * DIM + d]       = minhas[i].lo[d];
      envio[i * 2 * DIM + DIM + d] = minhas[i].hi[d];
    }
  int *contas_i = (int *) malloc ((size_t) np * sizeof *contas_i);
  int *desl_i   = (int *) malloc ((size_t) np * sizeof *desl_i);
  for (int r = 0; r < np; r++) { contas_i[r] = contas[r] * 2 * DIM; desl_i[r] = desl[r] * 2 * DIM; }
  int *todas_i = (int *) malloc ((size_t) (total > 0 ? total : 1) * 2 * DIM * sizeof *todas_i);
  MPI_Allgatherv (envio, nloc * 2 * DIM, MPI_INT, todas_i, contas_i, desl_i, MPI_INT, comm);
  free (envio); free (contas_i); free (desl_i);

  Caixa *todas = (Caixa *) malloc ((size_t) (total > 0 ? total : 1) * sizeof *todas);
  for (int k = 0; k < total; k++)
    for (int d = 0; d < DIM; d++) {
      todas[k].lo[d] = todas_i[k * 2 * DIM + d];
      todas[k].hi[d] = todas_i[k * 2 * DIM + DIM + d];
    }
  free (todas_i);

  // ------------------------------------------------------------ vizinhos
  unsigned idx = 0;
  for (int r = 0; r < np; r++) {
    if (r == rank) continue;

    // Aninhamento FIXO: minhas caixas por fora, as dele por dentro.  O outro lado
    // percorre o mesmo par (minha, dele) na mesma ordem, entao a i-esima faixa
    // que eu envio e' a i-esima arvore que ele recebe.  Trocar o aninhamento aqui
    // e nao la' entregaria valor de outra celula, sem erro nenhum.
    unsigned ns = 0, nr = 0;
    for (int a = 0; a < nloc; a++)
      for (int b = 0; b < contas[r]; b++) {
        Caixa s = fatia (&minhas[a], &todas[desl[r] + b], f, nb);
        if (!vazia (&s)) ns++;
      }
    for (int b = 0; b < contas[r]; b++)
      for (int a = 0; a < nloc; a++) {
        Caixa t = fatia (&todas[desl[r] + b], &minhas[a], f, nb);
        if (!vazia (&t)) nr++;
      }
    if (ns == 0 && nr == 0) continue;

    struct neighbor_proc *nbp =
        (struct neighbor_proc *) malloc (sizeof *nbp);
    nbp->idx = idx++;
    nbp->to_send_count = ns;
    nbp->to_recv_count = nr;
    nbp->to_send = (struct to_send_fringe *)
        malloc ((size_t) (ns > 0 ? ns : 1) * sizeof *nbp->to_send);
    nbp->to_recv_trees = (hig_cell **)
        malloc ((size_t) (nr > 0 ? nr : 1) * sizeof *nbp->to_recv_trees);

    // A REGRA DE ORDEM, e ela e' a coisa mais delicada deste arquivo: em ambas as
    // listas, a caixa DONA DO DADO fica no laco de fora e a que PRECISA dele no de
    // dentro.  Assim os dois lados percorrem o mesmo par na mesma ordem.
    //
    // MEDIDO: transpondo o aninhamento do recebimento, np=3 reprova com 228 de
    // 304 celulas de franja com valor de outra posicao -- e np=2 PASSA, porque
    // com uma caixa por rank a transposicao e' 1x1 e invisivel.  Um teste so' em
    // np=2 nao pegaria este defeito.
    //
    // Antes o envio aninhava (minhas, dele) e o recebimento tambem -- mas o
    // remetente emite com as caixas DELE por fora, entao as listas saiam
    // TRANSPOSTAS.  Com uma caixa de cada lado (np=2) 1x1 esconde o defeito; com
    // np=3, onde um rank tem duas ou tres caixas, 228 de 304 celulas de franja
    // recebiam o valor de outra posicao.
    unsigned is = 0;
    for (int a = 0; a < nloc; a++)            // dona do dado: minha caixa
      for (int b = 0; b < contas[r]; b++) {   // quem precisa: caixa dele
        Caixa s = fatia (&minhas[a], &todas[desl[r] + b], f, nb);
        if (!vazia (&s)) {
          struct to_send_fringe *sf = &nbp->to_send[is++];
          sf->local_tree = sd_get_higtree (sd, a);
          for (int d = 0; d < DIM; d++) {
            sf->lo_idx[d] = s.lo[d] - minhas[a].lo[d];   // coordenadas DA ARVORE
            sf->hi_idx[d] = s.hi[d] - minhas[a].lo[d];
          }
        }
      }

    unsigned ir = 0;
    for (int b = 0; b < contas[r]; b++)       // dona do dado: caixa dele
      for (int a = 0; a < nloc; a++) {        // quem precisa: minha caixa
        const Caixa *dele = &todas[desl[r] + b];
        Caixa t = fatia (dele, &minhas[a], f, nb);
        if (!vazia (&t)) {
          Point flo, fhi;
          int ext[DIM];
          for (int d = 0; d < DIM; d++) {
            flo[d] = lo[d] + t.lo[d] * h[d];
            fhi[d] = lo[d] + t.hi[d] * h[d];
            ext[d] = t.hi[d] - t.lo[d];
          }
          hig_cell *ft = hig_create_root (flo, fhi);
          hig_refine_uniform (ft, ext);
          sd_add_fringe_higtree (sd, ft);
          struct tree_properties *tp =
              (struct tree_properties *) malloc (sizeof *tp);
          tp->partition_group = 0;
          tp->is_fringe = true;
          g_hash_table_insert (pg->tree_props, ft, tp);
          nbp->to_recv_trees[ir++] = ft;
        }
      }

    g_hash_table_insert (pg->neighbors, GINT_TO_POINTER (r), nbp);
  }

  free (minhas); free (todas); free (contas); free (desl);
  return 1;
}
