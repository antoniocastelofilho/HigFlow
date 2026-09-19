// Localizacao por ponto: o desempate tem de ser DEFINIDO, nao emergente.
//
// `sd_get_cell_with_point` / `hig_get_cell_with_point` sao usadas em oito arquivos
// da HiGFlow, e todo estencil depende delas.  Para a etapa de incorporar uma
// segunda implementacao de malha (t8code), este e' o ponto mais fragil: a floresta
// dele e' arranjo linear sobre curva de preenchimento, sem ponteiro para pai, e a
// descida em arvore que a HiGTree faz hoje nao tem equivalente direto.
//
// O QUE ESTE TESTE FIXA.  Um ponto no INTERIOR de uma celula tem resposta obvia.
// Um ponto exatamente sobre uma FACE ou um CANTO pertence, geometricamente, a duas
// ou mais celulas -- e alguma tem de ganhar.  Se hoje quem ganha e' acidente da
// ordem de descida, isso vira comportamento observavel do qual o resto depende, e
// precisa estar escrito ANTES de existir uma segunda implementacao que possa
// desempatar de outro jeito.
//
// As tres afirmacoes, em ordem de forca:
//   1. a celula devolvida CONTEM o ponto (sanidade);
//   2. a resposta nao depende de como o dominio foi dividido em arvores;
//   3. a convencao de empate e' a que esta' registrada aqui.

#include <stdlib.h>
#include "testing.h"
#include "domain.h"
#include "utils.h"

#define LADO 8

// Monta [0,1]^DIM com `nblocos` arvores lado a lado em x.
static sim_domain *monta(int nblocos, hig_cell *raiz[]) {
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
    mp_mapper *m = mp_create();
    higcit_celliterator *todas = higcit_create_concat(its, nblocos);
    mp_assign_from_celliterator(m, todas, 0);
    higcit_destroy(todas);
    sim_domain *sd = sd_create(m);
    for(int b = 0; b < nblocos; b++) sd_add_higtree(sd, raiz[b]);
    return sd;
}

// Devolve o centro da celula que contem p, ou NULL.
static int centro_de(sim_domain *sd, const Point p, Point saida) {
    hig_cell *c = sd_get_cell_with_point(sd, (real *) p);
    if(c == NULL) return 0;
    hig_get_center(c, saida);
    return 1;
}

int main(void) {
    const real h = 1.0 / LADO;
    hig_cell *r1[4], *r2[4];
    sim_domain *uma  = monta(1, r1);
    sim_domain *duas = monta(2, r2);

    // ------------------------------------------------------------ sanidade
    t_case("celula_devolvida_contem_o_ponto");
    for(int i = 0; i < LADO; i++) {
        for(int j = 0; j < LADO; j++) {
            Point p;
            for(int d = 0; d < DIM; d++) p[d] = 0.5;
            p[0] = (i + 0.5) * h;      // interior, sem ambiguidade
            p[1] = (j + 0.5) * h;
            hig_cell *c = sd_get_cell_with_point(uma, p);
            T_CHECK_MSG(c != NULL, "ponto interior (%.4f, %.4f) nao achou celula",
                        p[0], p[1]);
            if(!c) continue;
            Rect bb; hig_get_bounding_box(c, &bb);
            for(int d = 0; d < DIM; d++) {
                T_CHECK_MSG(p[d] >= bb.lo[d] - 1e-12 && p[d] <= bb.hi[d] + 1e-12,
                    "celula devolvida nao contem o ponto na direcao %d: "
                    "p=%.6f fora de [%.6f, %.6f]", d, p[d], bb.lo[d], bb.hi[d]);
            }
        }
    }

    // --------------------------------------------- independencia da divisao
    // Faces e cantos INCLUSIVE.  Se o desempate mudar com o numero de arvores,
    // qualquer comparacao entre decomposicoes deixa de valer -- e uma segunda
    // implementacao de malha nao teria como ser conferida contra esta.
    t_case("empate_nao_depende_de_como_o_dominio_foi_dividido");
    int divergentes = 0;
    Point pior;  POINT_ASSIGN_SCALAR(pior, 0.0);
    for(int i = 0; i <= 2 * LADO; i++) {
        for(int j = 0; j <= 2 * LADO; j++) {
            Point p, c1, c2;
            for(int d = 0; d < DIM; d++) p[d] = 0.5;
            p[0] = i * h * 0.5;        // passo de meia celula: cai em centro,
            p[1] = j * h * 0.5;        // em face e em canto, alternadamente
            if(p[0] <= 0.0 || p[0] >= 1.0) continue;
            const int a = centro_de(uma,  p, c1);
            const int b = centro_de(duas, p, c2);
            T_CHECK_MSG(a == b,
                "ponto (%.4f, %.4f): uma arvore %s, duas arvores %s",
                p[0], p[1], a ? "achou" : "nao achou", b ? "achou" : "nao achou");
            if(!a || !b) continue;
            for(int d = 0; d < DIM; d++) {
                if(fabs(c1[d] - c2[d]) > 1e-12) {
                    if(divergentes == 0) POINT_ASSIGN(pior, p);
                    divergentes++;
                    break;
                }
            }
        }
    }
    T_CHECK_MSG(divergentes == 0,
        "%d pontos caem em celulas DIFERENTES conforme o dominio seja uma arvore "
        "ou duas; o primeiro e' (%.4f, %.4f) -- o desempate depende da divisao",
        divergentes, pior[0], pior[1]);

    // --------------------------------------------------- a convencao, fixada
    // Um ponto exatamente sobre uma face interna pertence a duas celulas.  Qual
    // ganha e' o que se registra aqui.  MEDIDO em 2026-09-19: ganha a celula do
    // lado de MENOR coordenada -- o intervalo da celula e' fechado em cima e
    // aberto embaixo.  (Eu esperava o contrario; a sonda desmentiu, e o valor
    // medido e' que vale.)
    //
    // Isto nao e' preferencia, e' contrato: oito arquivos da HiGFlow chamam
    // localizacao por ponto e o estencil inteiro depende dela.  Uma segunda
    // implementacao de malha que desempate ao contrario mudaria resultado em
    // silencio, e e' para isso que este caso existe.
    t_case("convencao_de_empate_na_face_esta_fixada");
    {
        Point p, c;
        for(int d = 0; d < DIM; d++) p[d] = 0.5 * h + 0.5 * h;   // centro em y,z
        p[0] = 4 * h;                 // EXATAMENTE sobre uma face interna
        for(int d = 1; d < DIM; d++) p[d] = 4.5 * h;
        if(centro_de(uma, p, c)) {
            const real esperado_maior = (4 + 0.5) * h;   // celula a' direita
            const real esperado_menor = (3 + 0.5) * h;   // celula a' esquerda
            const int ganhou_maior = fabs(c[0] - esperado_maior) < 1e-12;
            const int ganhou_menor = fabs(c[0] - esperado_menor) < 1e-12;
            T_CHECK_MSG(ganhou_maior || ganhou_menor,
                "centro devolvido %.6f nao e' nenhuma das duas celulas vizinhas "
                "(%.6f ou %.6f)", c[0], esperado_menor, esperado_maior);
            T_CHECK_MSG(ganhou_menor,
                "convencao de empate MUDOU: o ponto sobre a face foi para a celula "
                "de MAIOR coordenada (centro %.6f); a convencao registrada e' a de "
                "MENOR (%.6f). Se a mudanca e' intencional, atualize este caso -- "
                "ele existe para que uma segunda implementacao de malha nao "
                "desempate diferente em silencio", c[0], esperado_menor);
        } else {
            T_CHECK_MSG(0, "ponto sobre a face interna nao achou celula nenhuma");
        }
    }

    sd_destroy(uma);
    sd_destroy(duas);
    return t_end();
}
