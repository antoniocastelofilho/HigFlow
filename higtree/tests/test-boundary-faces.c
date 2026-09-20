// Onde o dominio termina, dito por dois mecanismos opostos -- C12 e C14, metade
// de Mesh.
//
// A HiGTree representa contorno EXPLICITAMENTE: arvores `sim_boundary`
// registradas com `sd_add_boundary`.  A malha nao sabe onde termina; alguem lhe
// conta.  O t8code representa IMPLICITAMENTE: face de elemento sem vizinho.
//
// Os dois mecanismos sao opostos, e e' por isso que vale afirmar que produzem a
// mesma geometria.  Um backend que errasse aqui faria o fechamento de contorno
// projetar sobre parede que nao existe, ou ignorar parede que existe -- e o
// fechamento e' onde viveu, por anos, o defeito dos sete sitios do `bc_inter`,
// invisivel porque os exemplos usam contorno constante.
//
// O QUE ESTE TESTE NAO FAZ, e o limite e' deliberado: nao fecha estencil.  A
// projecao sobre a parede, a escolha da parede ATRAVESSADA (C12) e a
// interpolacao em DIM-1 (C14) pertencem a Discretization/Boundary e sao as mesmas
// para qualquer malha -- ja' verificadas contra o MTree em test-stencil-selection
// e test-boundary-path.  Aqui se verifica o que e' da MALHA: onde o contorno
// esta'.
//
// O que se afirma:
//   contagem_de_faces_de_contorno   2*DIM paredes, LADO^(DIM-1) faces em cada
//   faces_caem_nas_paredes          toda face tem uma coordenada em 0 ou 1, e as
//                                   demais no centro de uma celula
//   normais_apontam_para_fora       a normal e' o eixo da parede, e aponta para
//                                   fora do dominio
//   os_dois_mecanismos_concordam    [t8code] o conjunto implicito e' igual ao
//                                   explicito, face a face

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "higtree.h"
#include "higtree-iterator.h"
#include "domain.h"
#include "utils.h"
#include "testing.h"

#ifdef HIGTREE_COM_T8CODE
#include "t8code/t8-boundary-faces.h"
#endif

#define LADO  8
#define MAXF  2048

// ----------------------------------------------------------------- MTree
// O contorno explicito: uma arvore de contorno por parede, com LADO^(DIM-1)
// celulas, exatamente como os exemplos e os outros testes montam.
static int mtree_faces(int maxn, Point centros[], Point normais[]) {
    int n = 0;
    for(int d = 0; d < DIM && n < maxn; d++) {
        for(int lado = 0; lado < 2 && n < maxn; lado++) {
            Point lo, hi; int nc[DIM];
            for(int k = 0; k < DIM; k++) { lo[k] = 0.0; hi[k] = 1.0; nc[k] = LADO; }
            const real plano = lado ? 1.0 : 0.0;
            lo[d] = plano; hi[d] = plano; nc[d] = 1;
            hig_cell *bc = hig_create_root(lo, hi);
            hig_refine_uniform(bc, nc);
            higcit_celliterator *it;
            for(it = higcit_create_all_leaves(bc);
                !higcit_isfinished(it) && n < maxn; higcit_nextcell(it)) {
                Point c; hig_get_center(higcit_getcell(it), c);
                for(int k = 0; k < DIM; k++) {
                    centros[n][k] = c[k];
                    normais[n][k] = 0.0;
                }
                normais[n][d] = lado ? 1.0 : -1.0;
                n++;
            }
            higcit_destroy(it);
        }
    }
    return n;
}

typedef struct {
    const char *nome;
    int (*faces)(int maxn, Point centros[], Point normais[]);
} FonteDeContorno;

#ifdef HIGTREE_COM_T8CODE
static int t8_faces(int maxn, Point centros[], Point normais[]) {
    return t8_faces_de_contorno(1, maxn, centros, normais);
}
#endif

static const FonteDeContorno FONTES[] = {
    { "mtree", mtree_faces },
#ifdef HIGTREE_COM_T8CODE
    { "t8code", t8_faces },
#endif
};

static void verifica(const FonteDeContorno *f) {
    const real h = 1.0 / LADO;
    static Point centros[MAXF], normais[MAXF];
    const int n = f->faces(MAXF, centros, normais);

    t_case("contagem_de_faces_de_contorno");
    {
        int esperado = 1;
        for(int d = 1; d < DIM; d++) esperado *= LADO;
        esperado *= 2 * DIM;               // 2*DIM paredes
        T_CHECK_MSG(n == esperado,
            "[%s] %d faces de contorno, esperado %d (%d paredes de %d faces)",
            f->nome, n, esperado, 2 * DIM, esperado / (2 * DIM));
    }

    t_case("faces_caem_nas_paredes");
    {
        int fora = 0;
        char primeira[256]; primeira[0] = '\0';
        for(int i = 0; i < n; i++) {
            int eixos_no_plano = 0, eixo = -1;
            for(int d = 0; d < DIM; d++) {
                if(fabs(centros[i][d]) < 1e-12 || fabs(centros[i][d] - 1.0) < 1e-12) {
                    eixos_no_plano++; eixo = d;
                }
            }
            // exatamente um eixo no plano da parede (canto contaria dois, e canto
            // nao e' centro de face)
            int ok = (eixos_no_plano == 1);
            if(ok) {
                for(int d = 0; d < DIM && ok; d++) {
                    if(d == eixo) continue;
                    const real k = centros[i][d] / h - 0.5;
                    if(fabs(k - round(k)) > 1e-9) ok = 0;
                }
            }
            if(!ok) {
                if(!fora) snprintf(primeira, sizeof primeira,
                    "face %d em (%.6f, %.6f): %d eixo(s) no plano da parede",
                    i, centros[i][0], centros[i][DIM > 1 ? 1 : 0], eixos_no_plano);
                fora++;
            }
        }
        T_CHECK_MSG(fora == 0,
            "[%s] %d de %d faces nao caem numa parede com as demais coordenadas "
            "em centro de celula.  %s", f->nome, fora, n, primeira);
    }

    t_case("normais_apontam_para_fora");
    {
        int erradas = 0;
        char primeira[256]; primeira[0] = '\0';
        for(int i = 0; i < n; i++) {
            int eixo = -1; real esperado = 0.0;
            for(int d = 0; d < DIM; d++) {
                if(fabs(centros[i][d]) < 1e-12)       { eixo = d; esperado = -1.0; }
                else if(fabs(centros[i][d] - 1.0) < 1e-12) { eixo = d; esperado = 1.0; }
            }
            if(eixo < 0) continue;
            int ok = fabs(normais[i][eixo] - esperado) < 1e-9;
            for(int d = 0; d < DIM && ok; d++) {
                if(d != eixo && fabs(normais[i][d]) > 1e-9) ok = 0;
            }
            if(!ok) {
                if(!erradas) snprintf(primeira, sizeof primeira,
                    "face em (%.4f, %.4f) na parede do eixo %d: normal "
                    "(%.4f, %.4f), esperada %+.0f naquele eixo e zero nos demais",
                    centros[i][0], centros[i][DIM > 1 ? 1 : 0], eixo,
                    normais[i][0], normais[i][DIM > 1 ? 1 : 0], esperado);
                erradas++;
            }
        }
        T_CHECK_MSG(erradas == 0,
            "[%s] %d de %d normais nao apontam para fora.  %s",
            f->nome, erradas, n, primeira);
    }
}

#ifdef HIGTREE_COM_T8CODE
// O conjunto implicito do t8code contra o explicito da HiGTree, face a face.
// Comparacao por CONJUNTO: a ordem e' consequencia de como cada um percorre.
static void verifica_concordancia(void) {
    t_case("os_dois_mecanismos_concordam");
    static Point ca[MAXF], na[MAXF], cb[MAXF], nb[MAXF];
    const int a = mtree_faces(MAXF, ca, na);
    const int b = t8_faces(MAXF, cb, nb);
    if(a != b) {
        T_CHECK_MSG(0, "contagens diferentes: contorno explicito %d faces, "
                       "implicito %d", a, b);
        return;
    }
    static char usado[MAXF];
    for(int j = 0; j < b; j++) usado[j] = 0;
    int sem_par = 0;
    char primeira[256]; primeira[0] = '\0';
    for(int i = 0; i < a; i++) {
        int achou = -1;
        for(int j = 0; j < b && achou < 0; j++) {
            if(usado[j]) continue;
            int bate = 1;
            for(int d = 0; d < DIM && bate; d++) {
                if(fabs(ca[i][d] - cb[j][d]) > 1e-12) bate = 0;
                if(fabs(na[i][d] - nb[j][d]) > 1e-9)  bate = 0;
            }
            if(bate) achou = j;
        }
        if(achou < 0) {
            if(!sem_par) snprintf(primeira, sizeof primeira,
                "face explicita em (%.6f, %.6f) com normal (%.1f, %.1f) nao tem "
                "par no conjunto implicito", ca[i][0], ca[i][DIM > 1 ? 1 : 0],
                na[i][0], na[i][DIM > 1 ? 1 : 0]);
            sem_par++;
        } else usado[achou] = 1;
    }
    T_CHECK_MSG(sem_par == 0,
        "%d de %d faces do contorno explicito nao tem par no implicito.  %s",
        sem_par, a, primeira);
}
#endif

int main(int argc, char *argv[]) {
    higtree_initialize(&argc, &argv);
    for(unsigned i = 0; i < sizeof FONTES / sizeof *FONTES; i++) {
        verifica(&FONTES[i]);
    }
#ifdef HIGTREE_COM_T8CODE
    verifica_concordancia();
#endif
    return t_end();
}
