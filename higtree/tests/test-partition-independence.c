// A MESMA geometria, montada de duas formas, tem de dar os MESMOS valores.
//
// Um dominio nao muda de fisica por ser descrito como uma arvore ou como varias.
// Se muda, quem escolheu a divisao escolheu a resposta -- e nenhuma comparacao
// entre decomposicoes significa mais nada.
//
// Esta propriedade custou tres dias para ser estabelecida por via indireta, em
// 2026-09-15..17, porque so' era observavel atraves de uma simulacao 3D inteira.
// Aqui ela e' uma consulta de biblioteca, sem solver e sem MPI: as duas variantes
// cobrem a MESMA regiao com as MESMAS celulas e os MESMOS contornos, mudando so'
// o numero de arvores.
//
//     variante UMA      [0,1]x[0,1] em uma arvore de 8x8
//     variante DUAS     [0,0.5]x[0,1] e [0.5,1]x[0,1], de 4x8 cada
//
// Os pontos de consulta se concentram ao redor de x=0.5, que na variante DUAS e'
// interface entre blocos e na variante UMA nao e' nada -- e' la' que a diferenca
// aparece, se existir.
//
// NOTA DE ESTADO (2026-09-18): a saida VTK do example3d_complex FOI medida como
// dependente de particao, com 0,56% de divergencia em nos de fronteira entre
// blocos, e o usuario decidiu NAO consertar -- o defeito e' do caminho de
// visualizacao, e o desempate no corte do suporte da interpolacao e' compartilhado
// com o solver.  Se este teste falhar, e' esse mesmo mecanismo aparecendo na
// biblioteca, e nao regressao nova.  Ver a entrada `higtree-suite` na memoria.

#include <stdlib.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

// Malha fina de proposito.  O mecanismo que torna a saida dependente de particao
// e' o EMPATE na distancia de corte do suporte: a interpolacao guarda no maximo
// `maxpts` pontos (60 para ordem 3 em 2D), e quando ha' mais candidatos do que
// isso na mesma faixa de distancia, quais entram depende da ordem de visita das
// arvores.  Com malha grosseira o dominio inteiro cabe no suporte, nao ha' corte,
// e o teste passaria sem exercitar nada.
#define LADO 24

static real campo(const Point p) {
    real v = 0.5;
    const real coef[3] = { 2.0, -3.0, 1.5 };
    for(int d = 0; d < DIM; d++) v += coef[d] * p[d];
    return v;
}

typedef struct { sim_domain *sd; mp_mapper *m; real *fval; int n; } variante;

static void parede_dirichlet(sim_domain *sd, const Point lo, const Point hi,
                             const int nc[DIM]) {
    hig_cell *bc = hig_create_root((real *) lo, (real *) hi);
    hig_refine_uniform(bc, (int *) nc);
    mp_mapper *bm = mp_create();
    higcit_celliterator *it = higcit_create_all_leaves(bc);
    mp_assign_from_celliterator(bm, it, 0);
    higcit_destroy(it);
    sim_boundary *sb = sb_create(bc, DIRICHLET, bm);
    for(it = higcit_create_all_leaves(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ct; hig_get_center(c, ct);
        sb_set_value(sb, mp_lookup(bm, hig_get_id(c, 0)), campo(ct));
    }
    higcit_destroy(it);
    sd_add_boundary(sd, sb);
}

// `nblocos` arvores lado a lado ao longo de x, cobrindo [0,1]^DIM ao todo.
static variante monta(int nblocos) {
    variante v;
    hig_cell *raiz[4];
    higcit_celliterator *its[4];
    for(int b = 0; b < nblocos; b++) {
        Point lo, hi; int nc[DIM];
        for(int d = 0; d < DIM; d++) { lo[d] = 0.0; hi[d] = 1.0; nc[d] = LADO; }
        lo[0] = (real) b / nblocos;
        hi[0] = (real) (b + 1) / nblocos;
        nc[0] = LADO / nblocos;
        raiz[b] = hig_create_root(lo, hi);
        hig_refine_uniform(raiz[b], nc);
        its[b] = higcit_create_all_leaves(raiz[b]);
    }
    v.m = mp_create();
    higcit_celliterator *todas = higcit_create_concat(its, nblocos);
    v.n = mp_assign_from_celliterator(v.m, todas, 0);
    higcit_destroy(todas);

    v.sd = sd_create(v.m);
    for(int b = 0; b < nblocos; b++) sd_add_higtree(v.sd, raiz[b]);

    // O campo e' preenchido pela COORDENADA do centro, nao pelo id: os ids diferem
    // entre as variantes, o valor fisico nao.
    v.fval = (real *) calloc(v.n, sizeof *v.fval);
    for(int b = 0; b < nblocos; b++) {
        higcit_celliterator *it;
        for(it = higcit_create_all_leaves(raiz[b]); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            Point ct; hig_get_center(c, ct);
            v.fval[mp_lookup(v.m, hig_get_id(c, 0))] = campo(ct);
        }
        higcit_destroy(it);
    }

    // Contornos do dominio INTEIRO, identicos nas duas variantes.
    for(int d = 0; d < DIM; d++) {
        for(int lado = 0; lado < 2; lado++) {
            Point lo, hi; int nc[DIM];
            for(int k = 0; k < DIM; k++) { lo[k] = 0.0; hi[k] = 1.0; nc[k] = LADO; }
            const real plano = lado ? 1.0 : 0.0;
            lo[d] = plano - EPSDELTA; hi[d] = plano + EPSDELTA; nc[d] = 1;
            parede_dirichlet(v.sd, lo, hi, nc);
        }
    }

    return v;
}

static int maior_suporte = 0;

static real interpola(variante *v, const Point org, const Point x, sim_stencil *stn) {
    stn_reset(stn);
    sd_get_stencil(v->sd, org, x, 1.0, stn);
    const int ne = stn_get_numelems(stn);
    if(ne > maior_suporte) maior_suporte = ne;
    return stn_mult_vector(stn, v->fval) - stn_get_rhs(stn);
}

int main(void) {
    sim_stencil *stn = stn_create();
    variante uma  = monta(1);
    variante duas = monta(2);
    const real h = 1.0 / LADO;

    t_case("uma_arvore_contra_duas");
    int piores = 0;
    real maior = 0.0;
    Point pior_x;  POINT_ASSIGN_SCALAR(pior_x, 0.0);

    // Varre pontos ao redor de x=0.5, a interface da variante DUAS.
    for(int i = -3; i <= 3; i++) {
        for(int j = 1; j < LADO; j++) {
            Point org, x;
            for(int d = 0; d < DIM; d++) { org[d] = 0.5; x[d] = 0.5; }
            x[1]   = j * h;
            org[1] = j * h;
            x[0]   = 0.5 + i * h * 0.5;
            org[0] = x[0] - 0.5 * h;         // estencil vindo de oeste
            if(x[0] <= 0.0 || x[0] >= 1.0) continue;

            const real a = interpola(&uma,  org, x, stn);
            const real b = interpola(&duas, org, x, stn);
            if(fabs(a - b) > maior) { maior = fabs(a - b); POINT_ASSIGN(pior_x, x); }
            if(fabs(a - b) > 1e-10) piores++;
        }
    }

    // PRE-CONDICAO: se o suporte nunca chega perto de maxpts, nao houve corte e
    // portanto nao houve empate para desempatar -- o teste passaria sem exercitar
    // o mecanismo, e um verde assim nao significa nada.  Melhor falhar dizendo isso
    // do que reportar sucesso vazio.
    T_CHECK_MSG(maior_suporte >= 40,
        "suporte maximo foi %d pontos: a malha e' grossa demais para haver corte, "
        "entao este teste nao exercita o desempate e o resultado nao vale",
        maior_suporte);

    T_CHECK_MSG(piores == 0,
        "%d pontos diferem entre 1 arvore e 2; maior diferenca %.3e em "
        "(%.4f, %.4f) -- a mesma geometria muda de valor conforme a divisao",
        piores, maior, pior_x[0], pior_x[1]);

    // ------------------------------------------------------ O QUE FALTA
    //
    // A MESMA pergunta no dominio de FACETAS, consultado nos CANTOS das celulas a
    // partir do centro delas -- que e' o padrao do escritor VTK e o unico lugar
    // onde a dependencia de particao FOI medida (2026-09-18: 0,56% de divergencia
    // entre np=1 e np=2 em nos de fronteira entre blocos, contra 1e-10 medido
    // direto nas facetas do solver).
    //
    // Nao esta' aqui porque a montagem serial de um `sim_facet_domain` ainda nao
    // esta' entendida: `sfd_create(NULL, 0)` + `sfd_copy_higtrees_from_center_domain`
    // + `sfd_adjust_facet_ids` da' segfault na primeira consulta, com e sem
    // `sfd_set_interpolator_order` e `sfd_create_boundary`.  O higflow monta isso
    // em hig-flow-kernel.c por volta de 886, mas passando por `psfd_create` e pelo
    // dominio particionado, e a sequencia serial equivalente nao foi reproduzida.
    //
    // Quem for continuar: comece por hig-flow-kernel.c:886 e pela montagem de
    // contorno de faceta nos exemplos, nao por esta funcao.

    stn_destroy(stn);
    return t_end();
}
