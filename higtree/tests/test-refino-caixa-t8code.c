// O PRODUTOR DE MALHA REFINADA EM CAIXA -- que ate' agora nao tinha teste nenhum.
//
// `t8_produz_por_rank_brick_refinado` e' o que alimenta o refinamento local do
// HiGFlow, e dois defeitos dele so' apareceram como escoamento divergindo numa
// simulacao de meia hora.  Os dois sao geometricos e baratos de afirmar aqui.
//
//   nenhuma_arvore_de_uma_celula   uma arvore de UMA celula nao tem interior:
//       toda faceta dela e' de fronteira, e a franja declarada do HiGFlow e' de
//       cinco celulas.  Apareciam por dois caminhos distintos -- a familia
//       repartida DENTRO de um rank (curada fundindo irmaos) e a familia
//       repartida ENTRE ranks (curada por set_for_coarsening).  Este caso pega
//       os dois, e pega qualquer terceiro caminho que venha a existir.
//
//   razao_no_maximo_2_para_1       celulas vizinhas diferem de no maximo um
//       nivel.  E' o que `t8_forest_set_balance` promete; aqui se CONFERE, em
//       vez de confiar.
//
//   duas_camadas_entre_niveis      a interpolacao por minimos quadrados moveis
//       do HiGFlow quer DUAS camadas de um nivel antes de encontrar o proximo.
//       Refinar a mesma caixa em todas as passadas faz as interfaces de nivel 1
//       e 2 coincidirem na borda, e o balanceamento entao insere UMA camada --
//       uma nao basta.  O predicato: nenhuma celula de nivel L tem, a menos de
//       duas celulas, vizinha de nivel L-1 E vizinha de nivel L+1.
//
//       So' e' COMPLETO em np=1, onde a malha inteira e' local: o produtor
//       refinado nao materializa franja, entao com np>1 a vizinhanca que cruza
//       a fronteira de particao nao esta' visivel.  Com np>1 o caso confere o
//       que e' local -- que ja' e' a maior parte -- e a contagem global.
//
//   cobertura_global_bate          a soma das folhas de todos os ranks e' o
//       total que o t8code reporta.  Pega celula perdida e celula em dois ranks.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "t8code/t8-mesh-rank.h"

#define NB      16      // celulas base por direcao
#define REFINOS  2      // DOIS niveis: e' o caso em que a escada importa

static int falhou = 0;

static void caso(int rank, const char *nome, int ok)
{
    int todos = 0;
    MPI_Allreduce(&ok, &todos, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!todos) falhou = 1;
    if (rank == 0) printf("  %-34s %s\n", nome, todos ? "ok" : "FALHOU");
}

// O nivel de uma folha, a partir do tamanho dela contra a celula base.
static int nivel_de(hig_cell *c, real h_base)
{
    Point lo, hi;
    hig_get_lowpoint(c, lo);
    hig_get_highpoint(c, hi);
    const real h = hi[0] - lo[0];
    return (int) llround(log2((double) (h_base / h)));
}

// A folha que contem `x`, em qualquer arvore local.  NULL se nao for local --
// o que com np>1 acontece perto da fronteira de particao, e e' por isso que o
// caso das duas camadas so' e' completo em np=1.
static hig_cell *acha(t8_producao_rank *p, const Point x)
{
    for (int i = 0; i < p->n_locais; i++) {
        hig_cell *c = hig_get_cell_with_point(p->locais[i], x);
        if (c != NULL) return c;
    }
    return NULL;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    Point lo, hi, cx_lo, cx_hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    // Caixa interna, com folga das bordas do dominio, para que a escada tenha
    // espaco e o teste nao meca o recorte contra a fronteira.
    POINT_ASSIGN_SCALAR(cx_lo, 0.30);
    POINT_ASSIGN_SCALAR(cx_hi, 0.70);
    int nb[DIM];
    for (int d = 0; d < DIM; d++) nb[d] = NB;

    t8_producao_rank p;
    if (!t8_produz_por_rank_brick_refinado(lo, hi, nb, cx_lo, cx_hi, REFINOS, &p)) {
        if (rank == 0) printf("  producao FALHOU\n");
        MPI_Finalize();
        return 1;
    }

    const real h_base = (hi[0] - lo[0]) / (real) NB;

    // ---- nenhuma arvore de uma celula ---------------------------------------
    int uma_celula = 0;
    for (int i = 0; i < p.n_locais; i++) {
        long n = 0;
        higcit_celliterator *it;
        for (it = higcit_create_all_leaves(p.locais[i]);
             !higcit_isfinished(it); higcit_nextcell(it)) n++;
        higcit_destroy(it);
        if (n <= 1) uma_celula++;
    }
    caso(rank, "nenhuma_arvore_de_uma_celula", uma_celula == 0);

    // ---- razao 2:1 e duas camadas -------------------------------------------
    long viola_21 = 0, viola_2camadas = 0, conferidas = 0, vizinhos_fora = 0;
    for (int i = 0; i < p.n_locais; i++) {
        higcit_celliterator *it;
        for (it = higcit_create_all_leaves(p.locais[i]);
             !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            const int L = nivel_de(c, h_base);
            Point cc, cl, ch;
            hig_get_center(c, cc);
            hig_get_lowpoint(c, cl);
            hig_get_highpoint(c, ch);
            conferidas++;

            int viu_menor = 0, viu_maior = 0;     // dentro de DUAS celulas
            for (int d = 0; d < DIM; d++) {
                const real h = ch[d] - cl[d];
                for (int s = -1; s <= 1; s += 2) {
                    // 1 celula: afirma a razao 2:1.  2 celulas: a espessura.
                    for (int k = 1; k <= 2; k++) {
                        Point x;
                        POINT_ASSIGN(x, cc);
                        x[d] += s * (k - 0.5) * h + s * 0.5 * h;
                        if (x[d] <= lo[d] || x[d] >= hi[d]) continue;
                        hig_cell *v = acha(&p, x);
                        if (v == NULL) { vizinhos_fora++; continue; }
                        const int Lv = nivel_de(v, h_base);
                        if (k == 1 && abs(Lv - L) > 1) viola_21++;
                        // So' os vizinhos IMEDIATOS entram aqui.  Com a faixa de
                        // exatamente duas camadas, uma celula da faixa ve o nivel
                        // menor de um lado e o maior a DUAS celulas do outro --
                        // usar k=2 exigiria mais de duas camadas, que nao e' a
                        // regra.  Adjacente aos dois e' faixa de UMA celula.
                        if (k == 1 && Lv < L) viu_menor = 1;
                        if (k == 1 && Lv > L) viu_maior = 1;
                    }
                }
            }
            // Adjacente a um nivel mais grosso E a um mais fino: a faixa de
            // nivel L tem UMA celula de espessura, e o estencil do MLS
            // atravessa duas transicoes de uma vez.
            if (viu_menor && viu_maior) viola_2camadas++;
        }
        higcit_destroy(it);
    }

    long g[5] = { viola_21, viola_2camadas, conferidas, p.n_local, vizinhos_fora };
    long gs[5];
    MPI_Allreduce(g, gs, 5, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    caso(rank, "razao_no_maximo_2_para_1", gs[0] == 0);
    caso(rank, "duas_camadas_entre_niveis", gs[1] == 0);
    caso(rank, "cobertura_global_bate", gs[3] == p.n_global);

    if (rank == 0) {
        printf("     (%ld folhas conferidas, %ld de %ld globais; "
               "%ld sondas cairam fora deste rank)\n",
               gs[2], gs[3], p.n_global, gs[4]);
        if (ntasks > 1)
            printf("     np>1: o caso das duas camadas e' PARCIAL -- o produtor\n"
                   "     refinado nao materializa franja, entao a vizinhanca que\n"
                   "     cruza a particao nao esta' visivel.\n");
    }

    t8_producao_rank_destroi(&p);
    MPI_Finalize();
    return falhou;
}
