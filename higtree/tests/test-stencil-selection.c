// QUAL parede fecha o estencil -- medido pelo VALOR, sem espiar o interior.
//
// A classe de defeito: o criterio antigo escolhia a BC pela menor distancia
// NORMAL ate' o plano do retalho.  Em dominio nao convexo isso fecha uma derivada
// em x contra uma parede normal a y, cujo plano por acaso passa mais perto.
// Medido em 2026-09-17 sobre 510 trocas de selecao em duas geometrias: o criterio
// antigo escolhia parede PERPENDICULAR a' direcao do estencil em 510 de 510.
//
// COMO TESTAR ISSO SEM EXPOR O INTERIOR DA BIBLIOTECA: dando valores CONSTANTES e
// DIFERENTES a cada parede.  Ai' o valor interpolado identifica sozinho qual
// parede fechou -- nao e' preciso gancho de depuracao nem observar `proj_dir`.
// O teste fica sobre comportamento, nao sobre implementacao.
//
// GEOMETRIA (degrau em L, o caso mais simples com canto reentrante):
//
//     y=1  +--------+--------+
//          |   A    |   B    |      A = [0,1] x [0,1]
//          |        |        |      B = [1,2] x [0.5,1]
//     y=.5 +--------+--------+
//          |   A    |  FORA  |      o ponto de consulta fica aqui
//     y=0  +--------+
//          x=0      x=1      x=2
//
// O ponto x=(1.05, 0.48) esta' fora.  O estencil vem de origem=(0.95, 0.48),
// dentro de A, andando em +x: ele ATRAVESSA a parede x=1 (a face leste de A
// abaixo de y=0.5).  Mas o plano da parede y=0.5 (a face sul de B) passa a 0.02
// do ponto, enquanto o plano de x=1 passa a 0.05 -- ou seja, A PAREDE ERRADA E'
// DUAS VEZES E MEIA MAIS PROXIMA.  E' exatamente a armadilha em que o criterio
// antigo caia.

#include <stdlib.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

#define VAL_LESTE  7.0     // parede x=1, abaixo de y=0.5   <- a ATRAVESSADA
#define VAL_SUL   -3.0     // parede y=0.5, a leste de x=1  <- a mais PROXIMA

// Acrescenta um retalho Dirichlet de valor constante.
static void parede(sim_domain *sd, const Point lo, const Point hi,
                   const int nc[DIM], real valor) {
    hig_cell *bc = hig_create_root((real *) lo, (real *) hi);
    hig_refine_uniform(bc, (int *) nc);
    mp_mapper *bm = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(bc);
    mp_assign_from_celliterator(bm, it, 0);
    higcit_destroy(it);
    sim_boundary *sb = sb_create(bc, DIRICHLET, bm);
    for(it = higcit_create_all_leaves(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
        sb_set_value(sb, mp_lookup(bm, hig_get_id(higcit_getcell(it), 0)), valor);
    }
    higcit_destroy(it);
    sd_add_boundary(sd, sb);
}

int main(void) {
    const real h = 0.125;          // 8 celulas em [0,1]
    sim_stencil *stn = stn_create();

    // ---- o dominio em L, como duas arvores
    Point la, ha, lb, hb;
    la[0] = 0.0; la[1] = 0.0;  ha[0] = 1.0; ha[1] = 1.0;
    lb[0] = 1.0; lb[1] = 0.5;  hb[0] = 2.0; hb[1] = 1.0;
#if DIM == 3
    la[2] = 0.0; ha[2] = 1.0;  lb[2] = 0.0; hb[2] = 1.0;
#endif
    int nca[DIM], ncb[DIM];
    nca[0] = 8; nca[1] = 8;
    ncb[0] = 8; ncb[1] = 4;
#if DIM == 3
    nca[2] = 8; ncb[2] = 8;
#endif
    hig_cell *A = hig_create_root(la, ha);  hig_refine_uniform(A, nca);
    hig_cell *B = hig_create_root(lb, hb);  hig_refine_uniform(B, ncb);

    mp_mapper *m = mp_create();
    higcit_celliterator *its[2] = { higcit_create_all_leaves(A),
                                    higcit_create_all_leaves(B) };
    higcit_celliterator *todas = higcit_create_concat(its, 2);
    const int n = mp_assign_from_celliterator(m, todas, 0);
    higcit_destroy(todas);

    sim_domain *sd = sd_create(m);
    sd_add_higtree(sd, A);
    sd_add_higtree(sd, B);

    // ---- as duas paredes em disputa, com valores distintos
    {   // face LESTE de A, so' abaixo de y=0.5: e' a que o estencil atravessa
        Point lo, hi; int nc[DIM];
        lo[0] = 1.0 - EPSDELTA; hi[0] = 1.0 + EPSDELTA; nc[0] = 1;
        lo[1] = 0.0;            hi[1] = 0.5;            nc[1] = 4;
#if DIM == 3
        lo[2] = 0.0;            hi[2] = 1.0;            nc[2] = 8;
#endif
        parede(sd, lo, hi, nc, VAL_LESTE);
    }
    {   // face SUL de B: o plano dela passa MAIS PERTO do ponto de consulta
        Point lo, hi; int nc[DIM];
        lo[0] = 1.0;            hi[0] = 2.0;            nc[0] = 8;
        lo[1] = 0.5 - EPSDELTA; hi[1] = 0.5 + EPSDELTA; nc[1] = 1;
#if DIM == 3
        lo[2] = 0.0;            hi[2] = 1.0;            nc[2] = 8;
#endif
        parede(sd, lo, hi, nc, VAL_SUL);
    }

    // ---- a consulta
    Point org, x;
    org[0] = 1.0 - 0.4*h;  org[1] = 0.48;
    x[0]   = 1.0 + 0.4*h;  x[1]   = 0.48;
#if DIM == 3
    org[2] = 0.5;          x[2]   = 0.5;
#endif

    real *fval = (real *) calloc(n, sizeof *fval);   // campo nulo: so' a BC responde

    t_case("fecha_pela_parede_atravessada_e_nao_pela_mais_proxima");
    stn_reset(stn);
    sd_get_stencil(sd, org, x, 1.0, stn);
    const real v = stn_mult_vector(stn, fval) - stn_get_rhs(stn);

    // Com o campo interno nulo, o valor vem inteiro da parede que fechou.  A
    // extrapolacao de Lagrange pode escalar o valor, entao o criterio e' de QUAL
    // parede ele esta' mais perto, e nao igualdade exata: os dois valores tem
    // sinais opostos, o que torna a discriminacao robusta a esse fator.
    const real d_atravessada = fabs(v - VAL_LESTE);
    const real d_mais_perto  = fabs(v - VAL_SUL);
    T_CHECK_MSG(d_atravessada < d_mais_perto,
        "valor %.6f: esta' mais perto da parede SUL (%.1f, a mais proxima em plano) "
        "do que da parede LESTE (%.1f, a que o estencil atravessa) -- "
        "selecao por distancia normal em vez de travessia",
        v, VAL_SUL, VAL_LESTE);

    free(fval);
    stn_destroy(stn);
    sd_destroy(sd);
    return t_end();
}
