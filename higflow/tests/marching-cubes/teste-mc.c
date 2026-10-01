// Arreio do marching cubes.  Standalone: o nucleo nao depende do HiGFlow.
//
// Os oraculos sao de dois tipos, e a distincao importa:
//
//   EXATOS (passa ou falha, sem tolerancia):
//     - a tabela dos 256 casos e' consistente com a geometria do cubo
//     - a superficie e' FECHADA: toda aresta usada por exatamente 2 triangulos
//     - a caracteristica de Euler e' a do genero certo: 2 na esfera, 4 em
//       duas esferas, 0 no toro
//     - dentro de um plato (campo constante 0 ou 1) nao sai triangulo nenhum
//
//   CONVERGENTES (valor com ordem medida):
//     - volume -> 4/3 pi R^3,  area -> 4 pi R^2,  |r_vertice| -> R
//
// O toro esta' aqui de proposito: ele exercita muito mais casos da tabela que
// a esfera, e Euler = 0 e' um numero afiado -- um buraco na superficie o
// desloca na hora.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hig-flow-mc.h"

static int falhas = 0;
static void ok(const char *nome, int passou, const char *det) {
    printf("  [%s] %-46s %s\n", passou ? "ok" : "FALHOU", nome, det ? det : "");
    if (!passou) falhas++;
}

// ---- campos analiticos -------------------------------------------------
typedef struct { int n; double L, R, c[3]; int tipo; double Rt, rt; } Campo;

static double _campo(int i, int j, int k, void *ctx) {
    Campo *C = (Campo*) ctx;
    double h = C->L / C->n;
    double x = (i+0.5)*h, y = (j+0.5)*h, z = (k+0.5)*h;
    double dx = x-C->c[0], dy = y-C->c[1], dz = z-C->c[2];
    if (C->tipo == 0) {                       // esfera, campo = R - r
        return C->R - sqrt(dx*dx+dy*dy+dz*dz);
    } else if (C->tipo == 1) {                // duas esferas separadas
        double r1 = sqrt((x-C->L*0.3)*(x-C->L*0.3)+dy*dy+dz*dz);
        double r2 = sqrt((x-C->L*0.7)*(x-C->L*0.7)+dy*dy+dz*dz);
        double a = C->R - r1, b = C->R - r2;
        return a > b ? a : b;                 // uniao
    } else if (C->tipo == 2) {                // toro em torno do eixo z
        double q = sqrt(dx*dx+dy*dy) - C->Rt;
        return C->rt - sqrt(q*q + dz*dz);
    } else {                                  // esfera "tipo VOF": 0/1 com
        double r = sqrt(dx*dx+dy*dy+dz*dz);   // transicao de uma celula
        double t = 0.5 + (C->R - r)/h;
        return t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
    }
}

int main(void) {
    printf("\n=== marching cubes: isosuperficie da fracao de volume ===\n\n");

    // -- 1. a tabela ------------------------------------------------------
    printf("-- tabela dos 256 casos --\n");
    {
        int ruins = mc_verifica_tabela();
        char d[80]; snprintf(d, sizeof d, "%d casos inconsistentes", ruins);
        ok("conjunto de arestas = o que a geometria exige", ruins == 0, d);
    }

    // -- 2. esfera: fechamento, genero, convergencia ----------------------
    printf("\n-- esfera R=0.25 em caixa 1^3 --\n");
    double eV[3], eA[3], eR[3]; int ns[3] = {32, 64, 128};
    for (int q = 0; q < 3; q++) {
        Campo C = { ns[q], 1.0, 0.25, {0.5,0.5,0.5}, 0, 0, 0 };
        double h = C.L/C.n, x0[3] = {0.5*h,0.5*h,0.5*h}, hh[3] = {h,h,h};
        mc_malha m; mc_malha_zera(&m);
        if (mc_extrai(C.n, C.n, C.n, x0, hh, 0.0, _campo, &C, &m)) {
            ok("extracao", 0, "faltou memoria"); return 1;
        }
        double V = mc_volume(&m), A = mc_area(&m);
        double Ve = 4.0/3.0*M_PI*0.25*0.25*0.25, Ae = 4.0*M_PI*0.25*0.25;
        double rmax = 0.0;
        for (int p = 0; p < m.nv; p++) {
            double dx=m.v[3*p]-0.5, dy=m.v[3*p+1]-0.5, dz=m.v[3*p+2]-0.5;
            double e = fabs(sqrt(dx*dx+dy*dy+dz*dz) - 0.25);
            if (e > rmax) rmax = e;
        }
        eV[q] = fabs(V-Ve)/Ve; eA[q] = fabs(A-Ae)/Ae; eR[q] = rmax/0.25;

        char nome[80], det[160];
        snprintf(nome, sizeof nome, "N=%-4d fechada (arestas orfas = 0)", ns[q]);
        int orf = mc_arestas_orfas(&m), eu = mc_euler(&m);
        snprintf(det, sizeof det, "V=%d F=%d  orfas=%d", m.nv, m.nt, orf);
        ok(nome, orf == 0, det);
        snprintf(nome, sizeof nome, "N=%-4d Euler = 2 (esfera)", ns[q]);
        snprintf(det, sizeof det, "chi=%d", eu);
        ok(nome, eu == 2, det);
        mc_malha_libera(&m);
    }
    printf("\n   %-6s %12s %12s %12s\n", "N", "err V", "err A", "err raio");
    for (int q = 0; q < 3; q++)
        printf("   %-6d %12.2e %12.2e %12.2e\n", ns[q], eV[q], eA[q], eR[q]);
    printf("   %-6s %12s %12s %12s\n", "ordem", "", "", "");
    for (int q = 1; q < 3; q++)
        printf("   %d->%d  %12.2f %12.2f %12.2f\n", ns[q-1], ns[q],
               log(eV[q-1]/eV[q])/log(2.0), log(eA[q-1]/eA[q])/log(2.0),
               log(eR[q-1]/eR[q])/log(2.0));
    {
        double pV = log(eV[1]/eV[2])/log(2.0);
        double pR = log(eR[1]/eR[2])/log(2.0);
        char d[120];
        snprintf(d, sizeof d, "ordem do volume = %.2f (esperado ~2)", pV);
        ok("volume converge em segunda ordem", pV > 1.6, d);
        snprintf(d, sizeof d, "ordem do raio = %.2f (esperado ~2)", pR);
        ok("raio converge em segunda ordem", pR > 1.6, d);
    }

    // -- 3. genero: duas esferas e um toro --------------------------------
    printf("\n-- genero: a tabela exercitada fora da esfera --\n");
    {
        Campo C = { 96, 1.0, 0.15, {0.5,0.5,0.5}, 1, 0, 0 };
        double h = C.L/C.n, x0[3]={0.5*h,0.5*h,0.5*h}, hh[3]={h,h,h};
        mc_malha m; mc_malha_zera(&m);
        mc_extrai(C.n,C.n,C.n,x0,hh,0.0,_campo,&C,&m);
        char d[120];
        snprintf(d, sizeof d, "chi=%d  orfas=%d  F=%d",
                 mc_euler(&m), mc_arestas_orfas(&m), m.nt);
        ok("duas esferas: Euler = 4", mc_euler(&m) == 4 && mc_arestas_orfas(&m) == 0, d);
        mc_malha_libera(&m);
    }
    {
        // O toro a N=96 tem so' ~9.6 celulas no raio menor -- bem mais
        // apertado que a esfera.  Por isso o erro de volume aqui NAO se
        // compara ao da esfera; o que se exige e' que ele caia como h^2.
        int nt[2] = {96, 192}; double eVt[2]; int bom = 1;
        char d[200];
        for (int q = 0; q < 2; q++) {
            Campo C = { nt[q], 1.0, 0, {0.5,0.5,0.5}, 2, 0.30, 0.10 };
            double h = C.L/C.n, x0[3]={0.5*h,0.5*h,0.5*h}, hh[3]={h,h,h};
            mc_malha m; mc_malha_zera(&m);
            mc_extrai(C.n,C.n,C.n,x0,hh,0.0,_campo,&C,&m);
            double V = mc_volume(&m), Ve = 2.0*M_PI*M_PI*0.30*0.10*0.10;
            eVt[q] = fabs(V-Ve)/Ve;
            snprintf(d, sizeof d, "N=%d  chi=%d  orfas=%d  F=%d  err V=%.2e",
                     nt[q], mc_euler(&m), mc_arestas_orfas(&m), m.nt, eVt[q]);
            char nm[80];
            snprintf(nm, sizeof nm, "toro N=%-4d Euler = 0 (genero 1)", nt[q]);
            if (!(mc_euler(&m) == 0 && mc_arestas_orfas(&m) == 0)) bom = 0;
            ok(nm, mc_euler(&m) == 0 && mc_arestas_orfas(&m) == 0, d);
            mc_malha_libera(&m);
        }
        double pt = log(eVt[0]/eVt[1])/log(2.0);
        snprintf(d, sizeof d, "%.2e -> %.2e, ordem %.2f", eVt[0], eVt[1], pt);
        ok("toro: volume converge em segunda ordem", bom && pt > 1.6, d);
    }

    // -- 4. o caso que de fato chega do solver ----------------------------
    printf("\n-- campo tipo VOF: 0/1 exatos com transicao de uma celula --\n");
    {
        Campo C = { 64, 1.0, 0.25, {0.5,0.5,0.5}, 3, 0, 0 };
        double h = C.L/C.n, x0[3]={0.5*h,0.5*h,0.5*h}, hh[3]={h,h,h};
        mc_malha m; mc_malha_zera(&m);
        mc_extrai(C.n,C.n,C.n,x0,hh,0.5,_campo,&C,&m);
        int orf = mc_arestas_orfas(&m);
        double V = mc_volume(&m), Ve = 4.0/3.0*M_PI*0.25*0.25*0.25;
        char d[160];
        snprintf(d, sizeof d, "chi=%d orfas=%d F=%d  err V=%.2e",
                 mc_euler(&m), orf, m.nt, fabs(V-Ve)/Ve);
        ok("fechada mesmo com platos 0/1", orf == 0 && mc_euler(&m) == 2, d);

        // O PLATO.  Fora da transicao o campo e' constante; nenhum vertice
        // pode cair ali.  E' isto que torna legitimo rodar o marching cubes
        // so' no nivel mais fino de uma malha adaptativa.
        double rmin = 1e9, rmax = -1e9;
        for (int p = 0; p < m.nv; p++) {
            double dx=m.v[3*p]-0.5, dy=m.v[3*p+1]-0.5, dz=m.v[3*p+2]-0.5;
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            if (r < rmin) rmin = r; if (r > rmax) rmax = r;
        }
        snprintf(d, sizeof d, "raios em [%.4f, %.4f], faixa=%.4f = %.2f h",
                 rmin, rmax, rmax-rmin, (rmax-rmin)/h);
        ok("todo vertice dentro da transicao de 1 celula", (rmax-rmin) < 2.2*h, d);
        mc_malha_libera(&m);
    }

    printf("\n%s  (%d falha(s))\n\n", falhas ? "FALHOU" : "PASSOU", falhas);
    return falhas ? 1 : 0;
}
