// Arreio da fracao volumetrica INICIAL em 3D (examples-common/vof-geometry-3d.c).
//
// Standalone: declara `real` e `Point` por conta propria, porque o modulo e'
// incluido por #include e nao depende de PETSc.
//
// OS ORACULOS SAO EXATOS, e cada um exercita uma parte diferente do codigo:
//
//   1. PLANO PELO CENTRO, normal qualquer -> fracao = 1/2 por simetria,
//      na precisao da maquina.  Exercita o ajuste do plano e a REFLEXAO dos
//      eixos de normal negativa.  Se o sinal ou a reflexao estiverem trocados,
//      este teste cai na hora, para qualquer normal.
//   2. PLANO ALINHADO A UM EIXO, deslocado -> fracao = a fatia, exata.
//      Exercita a DEGENERACAO (n_i = 0 no denominador da formula fechada) e o
//      aparo em 1e-9.
//   3. PLANO CORTANDO UM CANTO -> a parte de dentro e' um tetraedro de volume
//      conhecido.  Deixa UM termo da inclusao-exclusao ativo; e' o unico teste
//      que distingue a formula certa de uma com termo faltando.
//   4. ESFERA, volume total -> 4/3 pi R^3 em segunda ordem.  O valor em N=16
//      tem de bater com 6.11e-05, medido fora daqui por outra implementacao.
//   5. PLATOS: celula inteiramente de um lado da' 0 ou 1 EXATOS -- e' isto que
//      permite rodar o marching cubes so' no nivel mais fino de uma malha
//      adaptativa.
//
// Nota sobre a formula fechada de Scardovelli & Zaleski: ela existe duas vezes
// na arvore, aqui e (estatica) no modulo de front-tracking 3D, porque a
// restricao do projeto e' manter os dois caminhos separados no codigo.  A
// protecao NAO e' comparar as duas copias -- duas copias da mesma formula
// errada concordariam.  Cada uma e' presa a oraculos independentes: esta aos
// de 1-3 acima, a do front-tracking a' fracao amostrada do seu proprio arreio.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double real;
typedef double Point[3];
#define DIM 3

// A forma inicial, trocada pelos testes.
static int    _forma = 0;
static real   _n[3], _off, _R = 0.25, _c[3] = {0.5,0.5,0.5};

real func3(Point p) {
    if (_forma == 0)       // plano: n.(p) - off,  >0 de um lado
        return _n[0]*p[0] + _n[1]*p[1] + _n[2]*p[2] - _off;
    else {                 // esfera: >0 DENTRO
        real d0=p[0]-_c[0], d1=p[1]-_c[1], d2=p[2]-_c[2];
        return _R - sqrt(d0*d0+d1*d1+d2*d2);
    }
}

#include "../../examples-common/vof-geometry-3d.c"

static int falhas = 0;
static void ok(const char *nome, int passou, const char *det) {
    printf("  [%s] %-44s %s\n", passou ? "ok" : "FALHOU", nome, det ? det : "");
    if (!passou) falhas++;
}

int main(void) {
    printf("\n=== fracao volumetrica inicial 3D ===\n\n");

    Point ctr = {0.5, 0.5, 0.5};
    Point del = {0.1, 0.1, 0.1};

    // -- 1. plano pelo centro: 1/2 por simetria, normal qualquer -----------
    printf("-- plano pelo centro: fracao = 1/2 exata, para toda normal --\n");
    {
        // normais escolhidas para cobrir os oito octantes e casos oblicuos
        real dirs[10][3] = {
            {1,1,1},{-1,1,1},{1,-1,1},{1,1,-1},{-1,-1,1},
            {-1,1,-1},{1,-1,-1},{-1,-1,-1},{0.3,-0.9,0.4},{-0.7,0.2,-0.6}
        };
        real pior = 0.0; int piorq = -1;
        for (int q = 0; q < 10; q++) {
            real m = sqrt(dirs[q][0]*dirs[q][0]+dirs[q][1]*dirs[q][1]+dirs[q][2]*dirs[q][2]);
            for (int i = 0; i < 3; i++) _n[i] = dirs[q][i]/m;
            _forma = 0;
            _off = _n[0]*ctr[0] + _n[1]*ctr[1] + _n[2]*ctr[2];   // passa pelo centro
            real f = vof3_fracao_celula(ctr, del, 0.0);
            real e = fabs(f - 0.5);
            if (e > pior) { pior = e; piorq = q; }
        }
        char d[140];
        snprintf(d, sizeof d, "pior erro %.3e (normal %d)", pior, piorq);
        ok("1/2 na precisao da maquina", pior < 1e-14, d);
    }

    // -- 2. plano alinhado a um eixo: a fatia, exata -----------------------
    printf("\n-- plano alinhado a um eixo: exercita n_i = 0 --\n");
    {
        real pior = 0.0; char d[140];
        for (int eixo = 0; eixo < 3; eixo++) {
            for (int q = 1; q <= 9; q++) {
                real s = q/10.0;                  // fracao da celula de dentro
                for (int i = 0; i < 3; i++) _n[i] = 0.0;
                _n[eixo] = -1.0;                  // dentro = lado inferior
                _forma = 0;
                // dentro onde -x_e + off > 0  <=>  x_e < off
                // func3 = n.p - off com n = -e_eixo, logo dentro e'
                // -x - off > 0  <=>  x < -off:  off = -(limiar).
                _off = -((ctr[eixo] - 0.5*del[eixo]) + s*del[eixo]);
                real f = vof3_fracao_celula(ctr, del, 0.0);
                real e = fabs(f - s);
                if (e > pior) pior = e;
            }
        }
        snprintf(d, sizeof d, "pior erro %.3e nas 27 fatias", pior);
        ok("fatia exata nos tres eixos", pior < 1e-8, d);
    }

    // -- 3. plano cortando um canto: tetraedro ----------------------------
    printf("\n-- plano cortando um canto: um termo da inclusao-exclusao --\n");
    {
        // plano x+y+z = off com off perto do canto inferior; dentro = o canto.
        // Intercepta os eixos a (off - canto) de cada lado; o tetraedro tem
        // volume a^3/6 com a = off medido do canto.
        real pior = 0.0; char d[140];
        for (int q = 1; q <= 5; q++) {
            real a = q*0.015;                      // < del=0.1, fica no canto
            real lo[3] = {ctr[0]-0.5*del[0], ctr[1]-0.5*del[1], ctr[2]-0.5*del[2]};
            for (int i = 0; i < 3; i++) _n[i] = -1.0/sqrt(3.0);
            _forma = 0;
            // dentro: -(x+y+z)/sqrt3 + off > 0  <=> x+y+z < off*sqrt3
            _off = -(lo[0]+lo[1]+lo[2] + a)/sqrt(3.0);
            real f = vof3_fracao_celula(ctr, del, 0.0);
            real fe = (a*a*a/6.0) / (del[0]*del[1]*del[2]);
            real e = fabs(f - fe);
            if (e > pior) pior = e;
        }
        snprintf(d, sizeof d, "pior erro %.3e em 5 tetraedros", pior);
        ok("tetraedro de canto exato", pior < 1e-13, d);
    }

    // -- 4. esfera: volume total em segunda ordem -------------------------
    printf("\n-- esfera R=0.25: volume total, subdivisao 16^3 --\n");
    {
        _forma = 1;
        real Ve = 4.0/3.0*M_PI*_R*_R*_R;
        int NS[3] = {16, 32, 64}; real er[3];
        for (int q = 0; q < 3; q++) {
            int N = NS[q]; real h = 1.0/N, V = 0.0;
            Point dd = {h, h, h};
            for (int i=0;i<N;i++) for (int j=0;j<N;j++) for (int k=0;k<N;k++) {
                Point cc = {(i+0.5)*h, (j+0.5)*h, (k+0.5)*h};
                V += vof3_fracao(cc, dd, 0.0, 16) * h*h*h;
            }
            er[q] = fabs(V-Ve)/Ve;
            printf("   N=%-4d  V=%.10f  err=%.3e%s\n", N, V, er[q],
                   q ? "" : "");
        }
        for (int q = 1; q < 3; q++)
            printf("   ordem %d->%d: %.2f\n", NS[q-1], NS[q],
                   log(er[q-1]/er[q])/log(2.0));
        char d[160];
        real p = log(er[1]/er[2])/log(2.0);
        snprintf(d, sizeof d, "ordem %.2f", p);
        ok("volume converge em segunda ordem", p > 1.8, d);
        // VALIDACAO CRUZADA.  Uma segunda implementacao, independente --
        // plano tangente com a distancia EXATA a' esfera, em vez do ajuste
        // pelos cantos -- deu 6.109e-05 em N=16.  As duas concordam em quatro
        // digitos, o que nao e' consequencia de nada compartilhado: elas nao
        // compartilham linha de codigo.  A rampa que isto substitui errava
        // 1.23e-02, duas ordens de grandeza acima.
        snprintf(d, sizeof d, "N=16: %.4e vs 6.109e-05 independente (rampa: 1.23e-02)", er[0]);
        ok("bate com implementacao independente", fabs(er[0]-6.109e-5)/6.109e-5 < 2e-3, d);
    }

    // -- 5. platos exatos -------------------------------------------------
    printf("\n-- platos: 0 e 1 EXATOS fora da interface --\n");
    {
        _forma = 1;
        Point dd = {0.01, 0.01, 0.01};
        Point dentro = {0.5, 0.5, 0.5};            // centro da esfera
        Point fora   = {0.1, 0.1, 0.1};            // longe
        real fd = vof3_fracao(dentro, dd, 0.0, 16);
        real ff = vof3_fracao(fora,   dd, 0.0, 16);
        char d[140];
        snprintf(d, sizeof d, "dentro=%.17g  fora=%.17g", fd, ff);
        ok("1.0 e 0.0 sem residuo", fd == 1.0 && ff == 0.0, d);
    }

    printf("\n%s  (%d falha(s))\n\n", falhas ? "FALHOU" : "PASSOU", falhas);
    return falhas ? 1 : 0;
}
