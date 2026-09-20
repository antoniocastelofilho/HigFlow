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
// Espelha o criterio de `cell_find_in_center`: caixa da arvore, com igualdade de
// coordenada contando como contorno.  E' re-implementacao do criterio
// DOCUMENTADO, nao chamada da funcao -- ela e' estatica em domain.c.
static Classe mtree_classe(int nblocos, const Point p) {
    Classe r = FORA;
    for(int b = 0; b < nblocos; b++) {
        Rect bb;
        for(int d = 0; d < DIM; d++) { bb.lo[d] = 0.0; bb.hi[d] = 1.0; }
        bb.lo[0] = (real) b / nblocos;
        bb.hi[0] = (real) (b + 1) / nblocos;
        if(!rect_contains_point(&bb, (real *) p)) continue;
        r = DENTRO;
        for(int d = 0; d < DIM; d++) {
            if(fabs(p[d] - bb.lo[d]) < 1e-12 || fabs(p[d] - bb.hi[d]) < 1e-12) {
                r = NO_CONTORNO;
                break;
            }
        }
        break;                       // dominios nao se sobrepoem
    }
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
    // O caso registra a divergencia em vez de escolher um vencedor: qual dos
    // dois e' o desejado e' decisao de projeto, e muda o que o fechamento faz
    // numa interface interna.  O que NAO pode acontecer e' a divergencia passar
    // despercebida quando o t8code entrar.
    t_case("criterios_divergem_na_interface_entre_arvores");
    {
        Point p; for(int d = 0; d < DIM; d++) p[d] = 0.5625;
        p[0] = 0.5;                      // interface entre as duas arvores
        const Classe cx = mtree_classe(2, p);
        const Classe t8 = t8_classe(2, p);
        T_CHECK_MSG(cx == NO_CONTORNO,
            "o criterio da CAIXA deveria dizer NO_CONTORNO na interface entre "
            "arvores, e disse %s -- se isso mudou, o registro abaixo esta' velho",
            NOME[cx]);
        T_CHECK_MSG(t8 == DENTRO,
            "o criterio da FACE SEM VIZINHO deveria dizer DENTRO na interface "
            "entre arvores (ha' malha dos dois lados), e disse %s", NOME[t8]);
        T_CHECK_MSG(cx != t8,
            "os dois criterios passaram a CONCORDAR na interface entre arvores "
            "(%s).  Nao e' erro -- e' mudanca de comportamento, e este caso "
            "existe para que ela nao passe despercebida", NOME[cx]);
    }
#endif
    return t_end();
}
