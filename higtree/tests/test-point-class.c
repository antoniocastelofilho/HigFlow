// Dentro, sobre o contorno, ou fora -- o despacho que escolhe o ramo ON_BOUNDARY.
//
// C14 e' sobre o ramo ON_BOUNDARY do fechamento, e esse ramo e' escolhido por um
// DESPACHO que acontece antes de qualquer interpolacao (`cell_find_in_center`, em
// domain.c).  Metade do despacho de fechamento passa por ali, e foi nesse ramo que
// sobreviveram tres dos sete sitios do defeito do `bc_inter` -- verdes por
// ausencia de teste, nao por estarem certos.
//
// A DECISAO E' DE MALHA: "onde este ponto esta' em relacao ao dominio" nao depende
// de discretizacao nenhuma.  Quem interpola depois e' Discretization, e isso ja'
// esta' coberto em test-boundary-path.
//
// OS DOIS CRITERIOS SAO DIFERENTES, e e' o achado deste teste:
//
//   HiGTree  compara o ponto com a CAIXA de cada arvore do dominio.  Coordenada
//            igual a um limite da caixa => ON_BOUNDARY.
//   t8code   pergunta se o ponto cai sobre uma face SEM VIZINHO.  Face com
//            vizinho e' interface interna, e o ponto segue DENTRO.
//
// Onde o dominio e' uma arvore so', os dois concordam.  Onde ele e' dividido, a
// interface INTERNA entre arvores e' limite de caixa mas nao e' contorno do
// dominio -- e ai' os criterios se separam.  O caso final mede isso em vez de
// supor.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "rect.h"
#include "utils.h"
#include "testing.h"

#ifdef HIGTREE_COM_T8CODE
#include "t8code/t8-point-class.h"
#endif

#define LADO 8

typedef enum { DENTRO = 0, NO_CONTORNO = 1, FORA = 2 } Classe;
static const char *NOME[] = { "DENTRO", "NO_CONTORNO", "FORA" };

// ------------------------------------------------------------------ MTree
// CHAMA A FUNCAO REAL, `sd_classify_point`.  Ate' 2026-09-20 este teste
// REIMPLEMENTAVA o criterio, porque o classificador era estatico em domain.c --
// ou seja, ele afirmava coisas sobre uma copia.  Foi por isso que a divergencia
// da C14 precisou ser reconferida lendo a funcao antes de decidir o D3.
static Classe mtree_classe(int nblocos, const Point p) {
    sim_domain *sd = sd_create(NULL);
    for(int b = 0; b < nblocos; b++) {
        Point lo, hi;
        POINT_ASSIGN_SCALAR(lo, 0.0);
        POINT_ASSIGN_SCALAR(hi, 1.0);
        lo[0] = (real) b / nblocos;
        hi[0] = (real) (b + 1) / nblocos;
        hig_cell *raiz = hig_create_root(lo, hi);
        int nc[DIM];
        for(int d = 0; d < DIM; d++) nc[d] = LADO;
        nc[0] = LADO / nblocos;
        hig_refine_uniform(raiz, nc);
        sd_add_higtree(sd, raiz);
    }
    const Classe r = (Classe) sd_classify_point(sd, (CPPoint) p);
    sd_destroy(sd);
    return r;
}

typedef struct {
    const char *nome;
    Classe (*classe)(int nblocos, const Point p);
} Classificador;

#ifdef HIGTREE_COM_T8CODE
static Classe t8_classe(int nblocos, const Point p) {
    return (Classe) t8_classifica_ponto(nblocos, p);
}
#endif

static const Classificador CLASSIFICADORES[] = {
    { "mtree", mtree_classe },
#ifdef HIGTREE_COM_T8CODE
    { "t8code", t8_classe },
#endif
};

// Casos em que o contrato exige a MESMA resposta dos dois, com UMA arvore.
static void verifica(const Classificador *c) {
    const real h = 1.0 / LADO;

    t_case("interior_e_dentro");
    {
        int erradas = 0;
        for(int i = 1; i < LADO - 1; i++) {
            Point p; for(int d = 0; d < DIM; d++) p[d] = 0.5;
            p[0] = (i + 0.5) * h;
            if(c->classe(1, p) != DENTRO) erradas++;
        }
        T_CHECK_MSG(erradas == 0,
            "[%s] %d ponto(s) interior(es) nao classificados como DENTRO",
            c->nome, erradas);
    }

    t_case("parede_externa_e_contorno");
    {
        int erradas = 0;
        char primeira[256]; primeira[0] = '\0';
        for(int d = 0; d < DIM; d++) {
            for(int lado = 0; lado < 2; lado++) {
                Point p; for(int k = 0; k < DIM; k++) p[k] = 0.5;
                p[d] = lado ? 1.0 : 0.0;
                const Classe r = c->classe(1, p);
                if(r != NO_CONTORNO) {
                    if(!erradas) snprintf(primeira, sizeof primeira,
                        "parede do eixo %d em %.1f: %s", d, (double) p[d], NOME[r]);
                    erradas++;
                }
            }
        }
        T_CHECK_MSG(erradas == 0,
            "[%s] %d parede(s) externa(s) nao classificadas como NO_CONTORNO.  %s",
            c->nome, erradas, primeira);
    }

    t_case("fora_do_dominio_e_fora");
    {
        int erradas = 0;
        const real longe[] = { -0.1, 1.1 };
        for(int d = 0; d < DIM; d++) {
            for(int k = 0; k < 2; k++) {
                Point p; for(int q = 0; q < DIM; q++) p[q] = 0.5;
                p[d] = longe[k];
                if(c->classe(1, p) != FORA) erradas++;
            }
        }
        T_CHECK_MSG(erradas == 0,
            "[%s] %d ponto(s) fora do dominio nao classificados como FORA",
            c->nome, erradas);
    }

    t_case("interface_interna_de_celula_e_dentro");
    {
        // Face entre duas CELULAS, no meio do dominio: nao e' contorno por
        // criterio nenhum.  Serve de controle para o caso seguinte.
        Point p; for(int d = 0; d < DIM; d++) p[d] = 0.5625;
        p[0] = 4 * h;                    // face interna, longe das paredes
        const Classe r = c->classe(1, p);
        T_CHECK_MSG(r == DENTRO,
            "[%s] face entre duas celulas classificada como %s", c->nome, NOME[r]);
    }
}

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    for(unsigned i = 0; i < sizeof CLASSIFICADORES / sizeof *CLASSIFICADORES; i++) {
        verifica(&CLASSIFICADORES[i]);
    }

#ifdef HIGTREE_COM_T8CODE
    // ------------------------------------------------------------------
    // ONDE OS CRITERIOS SE SEPARAM, medido e nao suposto.
    //
    // Com o dominio dividido em duas arvores em x, o plano x=0,5 e' limite de
    // CAIXA para as duas, mas nao e' contorno do dominio -- ha' malha dos dois
    // lados.  O criterio da caixa diz NO_CONTORNO; o da face sem vizinho diz
    // DENTRO.
    //
    // DECIDIDO EM 2026-09-20 (C14): vale a semantica do t8code -- "sobre o
    // contorno" significa que o dominio TERMINA ali.  E O MTREE FOI CORRIGIDO no
    // mesmo dia, entao os dois CONCORDAM e este caso afirma a concordancia.
    //
    // O desvio antigo nascia no `cell_find_in_center`: ele parava no PRIMEIRO
    // higtree cuja caixa contem o ponto (`break`) e nunca perguntava se outro
    // bloco continuava o dominio -- a resposta chegava a depender da ORDEM das
    // arvores.  O classificador passou a ser `sd_classify_point`, que decide por
    // DIRECAO: so' e' contorno se o dominio nao continuar de um dos lados.
    //
    // ESTE CASO PEGOU A PROPRIA CORRECAO.  Ele estava escrito para falhar se o
    // MTree passasse a dizer DENTRO, com a mensagem de que isso seria boa
    // noticia -- e foi exatamente o que aconteceu quando a correcao entrou.
    t_case("os_dois_concordam_na_interface_entre_arvores");
    {
        Point p; for(int d = 0; d < DIM; d++) p[d] = 0.5625;
        p[0] = 0.5;                      // interface entre as duas arvores
        const Classe cx = mtree_classe(2, p);
        const Classe t8 = t8_classe(2, p);
        T_CHECK_MSG(t8 == DENTRO,
            "O CONTRATO: na interface entre arvores ha' malha dos dois lados, "
            "entao o ponto esta' DENTRO.  O t8code disse %s", NOME[t8]);
        T_CHECK_MSG(cx == DENTRO,
            "O MTree voltou a dizer %s na interface entre arvores.  Ele foi "
            "CORRIGIDO em 2026-09-20 para decidir por direcao; se regrediu, o "
            "suspeito e' o `sd_classify_point` ter voltado a parar no primeiro "
            "bloco que contem o ponto", NOME[cx]);
        T_CHECK_MSG(cx == t8,
            "os dois backends discordam na interface (%s contra %s) -- e' a C14 "
            "que esta' sendo violada por um dos dois", NOME[cx], NOME[t8]);
    }
#endif
    return t_end();
}
