// A franja: o que ela entrega ao estencil, e o que ela NAO entrega ao iterador.
//
// POR QUE ESTE TESTE EXISTE.  `build-fringe.cpp` sao 815 linhas sem cobertura, e a
// franja e' uma das quatro areas com acoplamento estrutural ao octree que a costura
// do Mesh vai ter de tocar.  Mais importante: quando o t8code entrar como segunda
// implementacao de particionamento, ele traz a PROPRIA camada de ghost
// (`t8_forest_ghost`).  Ou uma substitui a outra, ou as duas coexistem -- e em
// qualquer dos casos e' preciso ter escrito antes o que a franja tem de cumprir.
//
// O CRITERIO E' O QUE O ESTENCIL PEDE, nao um tamanho pre-fixado.  Afirmar "a franja
// tem N celulas" amarraria a implementacao; afirmar "o suporte do estencil atravessa
// a interface e carrega peso" amarra o CONTRATO, e qualquer implementacao que o
// cumpra passa.
//
// MONTAGEM (serial, sem MPI): o dominio [0,1]^DIM partido em dois blocos por x.  O
// esquerdo e' local; o direito entra como FRANJA, pelo mesmo padrao que o
// `_psd_setmapper` do pdomain.c usa -- ids locais pelo iterador de dominio, e cada
// arvore de franja por `higcit_create_all_leaves`.  Sem esse segundo laco as celulas
// da franja ficam sem id e o estencil as descarta em silencio: foi o primeiro
// resultado que esta sonda deu, e ele parecia dizer que a franja nao servia para
// nada.
//
// O QUE FOI MEDIDO, e virou asserção:
//
//                            sem franja   com franja
//     elementos do estencil       32           64
//     pontos alem da interface     0           32
//     peso carregado por eles     0,0        0,299
//     iterador: celulas / volume  32 / 0,5   32 / 0,5   (inalterado)
//
// UM ORACULO QUE NAO SERVE AQUI, e fica registrado para ninguem tentar de novo:
// reproducao polinomial NAO discrimina.  O campo linear e' reproduzido exato
// (4e-16) nos DOIS casos, porque minimos quadrados com suporte so' de um lado ainda
// acerta um linear por extrapolacao.  Ele entra abaixo como propriedade positiva --
// os dados da franja sao usados CERTO, nao apenas estao presentes -- e nao como
// guarda.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "testing.h"

#define NC_X     2      // celulas em x por bloco
#define NC_OUT   4      // celulas nas demais direcoes
#define INTERF   0.5    // onde os dois blocos se encontram

typedef struct {
    int    ne;          // elementos do estencil
    int    n_alem;      // quantos deles ficam alem da interface
    real   peso_alem;   // peso total que eles carregam
    real   soma_pesos;
    real   erro_linear; // |interpolado - exato| para f = 1 + 2x + 3y
    int    celulas;     // visitadas pelo iterador de dominio
    real   volume;      // somado pelo iterador de dominio
    int    n_arvores, n_locais, n_franja;
} Medida;

static hig_cell *bloco(real x0, real x1) {
    Point lo, hi;
    POINT_ASSIGN_SCALAR(lo, 0.0);
    POINT_ASSIGN_SCALAR(hi, 1.0);
    lo[0] = x0; hi[0] = x1;
    hig_cell *r = hig_create_root(lo, hi);
    int nc[DIM];
    for (int d = 0; d < DIM; d++) nc[d] = NC_OUT;
    nc[0] = NC_X;
    hig_refine_uniform(r, nc);
    return r;
}

static Medida mede(bool com_franja) {
    Medida md;
    sim_domain *sd = sd_create(NULL);
    sd_set_interpolator_order(sd, 2);
    sd_add_higtree(sd, bloco(0.0, INTERF));
    if (com_franja) sd_add_fringe_higtree(sd, bloco(INTERF, 1.0));

    // Ids: locais pelo iterador de dominio, franja por arvore.  Ver o
    // _psd_setmapper em pdomain.c -- sem o segundo laco a franja fica invisivel.
    mp_mapper *m = sd_get_domain_mapper(sd);
    higcit_celliterator *cit = sd_get_domain_celliterator(sd);
    unsigned livre = mp_assign_from_celliterator(m, cit, 0);
    higcit_destroy(cit);
    for (unsigned i = 0; i < sd_get_num_fringe_higtrees(sd); ++i) {
        cit = higcit_create_all_leaves(sd_get_fringe_higtree(sd, i));
        livre = mp_assign_from_celliterator(m, cit, livre);
        higcit_destroy(cit);
    }

    md.n_arvores = sd_get_num_higtrees(sd);
    md.n_locais  = (int) sd_get_num_local_higtrees(sd);
    md.n_franja  = (int) sd_get_num_fringe_higtrees(sd);

    md.volume = 0.0; md.celulas = 0;
    for (cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit);
         higcit_nextcell(cit)) {
        Point de;
        hig_get_delta(higcit_getcell(cit), de);
        real v = 1.0;
        for (int d = 0; d < DIM; d++) v *= de[d];
        md.volume += v; md.celulas++;
    }
    higcit_destroy(cit);

    // Consulta junto a' interface, do lado local.
    Point x;
    POINT_ASSIGN_SCALAR(x, INTERF);
    x[0] = INTERF - 0.5 * (INTERF / NC_X);      // centro da ultima celula local
    sim_stencil *stn = stn_create();
    sd_get_stencil(sd, x, x, 1.0, stn);
    md.ne = stn_get_numelems(stn);
    int  *ids = stn_get_ids(stn);
    real *w   = stn_get_vals(stn);

    md.soma_pesos = 0.0;
    for (int i = 0; i < md.ne; i++) md.soma_pesos += w[i];

    md.n_alem = 0; md.peso_alem = 0.0;
    real interp = 0.0;
    for (int t = 0; t < md.n_arvores; t++) {
        higcit_celliterator *ct = higcit_create_all_leaves(sd_get_higtree(sd, t));
        for (; !higcit_isfinished(ct); higcit_nextcell(ct)) {
            hig_cell *c = higcit_getcell(ct);
            int lid = mp_lookup(m, hig_get_cid(c));
            Point ce;
            hig_get_center(c, ce);
            real f = 1.0 + 2.0 * ce[0] + 3.0 * ce[1];
            for (int i = 0; i < md.ne; i++) {
                if (ids[i] != lid) continue;
                interp += w[i] * f;
                if (ce[0] > INTERF) { md.n_alem++; md.peso_alem += w[i]; }
            }
        }
        higcit_destroy(ct);
    }
    md.erro_linear = fabs(interp - (1.0 + 2.0 * x[0] + 3.0 * x[1]));

    stn_destroy(stn);
    return md;
}

int main(void) {
    Medida sem = mede(false);
    Medida com = mede(true);

    int cel_bloco = NC_X;
    for (int d = 1; d < DIM; d++) cel_bloco *= NC_OUT;

    t_case("contagens_de_arvore_separam_local_de_franja");
    T_CHECK_MSG(sem.n_arvores == 1 && sem.n_locais == 1 && sem.n_franja == 0,
        "sem franja: arvores %d, locais %d, franja %d (esperado 1, 1, 0)",
        sem.n_arvores, sem.n_locais, sem.n_franja);
    T_CHECK_MSG(com.n_arvores == 2 && com.n_locais == 1 && com.n_franja == 1,
        "com franja: arvores %d, locais %d, franja %d (esperado 2, 1, 1)",
        com.n_arvores, com.n_locais, com.n_franja);

    t_case("iterador_de_dominio_ignora_a_franja");
    // O contrato esta' escrito no proprio domain.h: "trees after this index will
    // not take part in the domain cell iterator".  Afirmado por VOLUME, que pega
    // celula a mais e celula a menos, e nao so' por contagem.
    T_CHECK_MSG(com.celulas == cel_bloco && sem.celulas == cel_bloco,
        "o iterador deveria ver so' o bloco local (%d celulas): sem franja %d, "
        "com franja %d", cel_bloco, sem.celulas, com.celulas);
    T_NEAR(sem.volume, INTERF, 1.0e-12, "volume visto pelo iterador, sem franja");
    T_NEAR(com.volume, INTERF, 1.0e-12, "volume visto pelo iterador, com franja");

    t_case("sem_franja_o_estencil_nao_atravessa");
    T_CHECK_MSG(sem.n_alem == 0 && sem.peso_alem == 0.0,
        "sem franja o suporte deveria parar na interface: %d pontos alem de "
        "x=%.3f, carregando peso %.6f", sem.n_alem, INTERF, sem.peso_alem);

    t_case("com_franja_o_estencil_atravessa_e_carrega_peso");
    // Este e' O criterio: a franja contem o que o estencil PEDE.  Nao se afirma
    // quantas celulas ela tem -- isso amarraria a implementacao -- e sim que o
    // suporte a alcanca e que ela pesa de verdade na interpolacao.
    T_CHECK_MSG(com.n_alem > 0,
        "com franja o suporte deveria atravessar a interface, e nao atravessou "
        "(0 pontos alem de x=%.3f)", INTERF);
    T_CHECK_MSG(com.ne > sem.ne,
        "o suporte deveria crescer com a franja: %d elementos sem, %d com",
        sem.ne, com.ne);
    T_CHECK_MSG(com.peso_alem > 0.05,
        "a franja deveria carregar peso nao desprezivel: %.6f (medido 0,299 na "
        "montagem de referencia)", com.peso_alem);

    t_case("pesos_somam_um_e_o_linear_e_exato");
    // Propriedade POSITIVA, nao guarda: os dados da franja sao usados certo.  Nao
    // discrimina -- o linear e' exato tambem sem a franja, porque minimos
    // quadrados de um lado so' ainda acerta um linear por extrapolacao.
    T_NEAR(com.soma_pesos, 1.0, 1.0e-12, "soma dos pesos com franja");
    T_NEAR(sem.soma_pesos, 1.0, 1.0e-12, "soma dos pesos sem franja");
    T_CHECK_MSG(com.erro_linear < 1.0e-12,
        "campo linear atraves da franja: erro %.3e", com.erro_linear);

    return t_end();
}
