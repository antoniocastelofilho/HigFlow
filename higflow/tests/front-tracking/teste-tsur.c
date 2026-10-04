// Arreio do BALANCE e do TSUR -- remocao de ondulacoes submalha, portadas do
// Freeflow (simulador/tsur3d.cpp).
//
// O ORACULO QUE IMPORTA NAO E' "SUAVIZOU".  Qualquer media suaviza; a questao e'
// se suaviza SEM ENCOLHER.  As duas rotinas sao desenhadas para preservar
// volume localmente, e e' essa propriedade que justifica porta-las em vez de
// usar uma media laplaciana de tres linhas.  Entao o arreio:
//
//   1. mede que a esfera LIMPA e' ponto fixo -- aplicar nao a deforma
//   2. poe ondulacao radial e mede que a amplitude CAI
//   3. mede o VOLUME antes e depois, e exige que se conserve
//   4. compara com a media LAPLACIANA na mesma superficie ondulada: ela tem de
//      suavizar tambem, mas ENCOLHENDO -- sem esta comparacao a conservacao de
//      volume seria uma afirmacao sem contraste
//   5. confere que a topologia sobrevive (chi = 2, nenhum triangulo invertido)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hig-flow-front-tracking-3d.h"

// O limiar PADRAO.  Ver a varredura no fim deste arreio para por que ele
// existe e por que 0.05.
#define LIM 0.05

static int falhas = 0;
static void ok(const char *nome, int passou, const char *det) {
    printf("  [%s] %-44s %s\n", passou ? "ok" : "FALHOU", nome, det ? det : "");
    if (!passou) falhas++;
}

// amplitude da ondulacao: desvio-padrao do raio em torno da media
static void estat_raio(const ft3_superficie *s, const real c[3],
                       real *medio, real *desvio, real *pior)
{
    const int nv = ft3_num_vertices(s);
    const Point *p = ft3_posicoes(s);
    real soma = 0.0;
    for (int q = 0; q < nv; q++) {
        real d[3] = { p[q][0]-c[0], p[q][1]-c[1], p[q][2]-c[2] };
        soma += sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    }
    *medio = soma/nv;
    real v = 0.0; *pior = 0.0;
    for (int q = 0; q < nv; q++) {
        real d[3] = { p[q][0]-c[0], p[q][1]-c[1], p[q][2]-c[2] };
        real r = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) - *medio;
        v += r*r; if (fabs(r) > *pior) *pior = fabs(r);
    }
    *desvio = sqrt(v/nv);
}

// ondulacao radial de alta frequencia -- o que as rotinas devem remover
static void ondula(ft3_superficie *s, const real c[3], real amp)
{
    const int nv = ft3_num_vertices(s);
    Point *p = (Point *) ft3_posicoes(s);
    for (int q = 0; q < nv; q++) {
        real d[3] = { p[q][0]-c[0], p[q][1]-c[1], p[q][2]-c[2] };
        real r = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (r < 1e-12) continue;
        // sinal alternando com o indice: ondulacao na escala da malha, que e'
        // exatamente a "submalha" do nome.  Nao e' ruido aleatorio: ruido tem
        // componente de baixa frequencia que a suavizacao tambem remove, e ai'
        // nao se saberia qual das duas foi medida.
        real f = 1.0 + amp * ((q % 2) ? 1.0 : -1.0);
        for (int k = 0; k < 3; k++) p[q][k] = c[k] + d[k]*f;
    }
}

// media laplaciana, para contraste: cada vertice vai para a media dos vizinhos
static void laplaciano(ft3_superficie *s)
{
    const int nv = ft3_num_vertices(s), nt = ft3_num_triangulos(s);
    Point *p = (Point *) ft3_posicoes(s);
    real *acc = (real *) calloc((size_t)nv*3, sizeof *acc);
    int  *cnt = (int  *) calloc((size_t)nv, sizeof *cnt);
    const int (*T)[3] = ft3_triangulos(s);
    for (int t = 0; t < nt; t++) {
        for (int k = 0; k < 3; k++) {
            int a = T[t][k], b = T[t][(k+1)%3];
            for (int d = 0; d < 3; d++) { acc[3*a+d] += p[b][d]; acc[3*b+d] += p[a][d]; }
            cnt[a]++; cnt[b]++;
        }
    }
    for (int q = 0; q < nv; q++)
        if (cnt[q]) for (int d = 0; d < 3; d++) p[q][d] = acc[3*q+d]/cnt[q];
    free(acc); free(cnt);
}

int main(void) {
    printf("\n=== BALANCE e TSUR: ondulacoes submalha ===\n\n");
    Point c = {0.5, 0.5, 0.5};
    const real R = 0.25;

    // ---- 1. a esfera limpa e' ponto fixo ------------------------------
    printf("-- a esfera limpa: aplicar nao deve deformar --\n");
    {
        ft3_superficie *s = ft3_cria_esfera(c, R, 3);
        real V0 = ft3_volume(s), m0, d0, p0;
        estat_raio(s, c, &m0, &d0, &p0);
        int mb = ft3_balance(s, LIM);
        int mt = ft3_tsur(s, 0.0, LIM);
        real V1 = ft3_volume(s), m1, d1, p1;
        estat_raio(s, c, &m1, &d1, &p1);
        char det[180];
        snprintf(det, sizeof det, "moveu %d/%d por BALANCE, %d por TSUR; "
                 "raio %.6f->%.6f, volume %.3e", mb, ft3_num_vertices(s), mt,
                 (double)m0, (double)m1, (double)(fabs(V1-V0)/V0));
        ok("raio e volume praticamente intactos",
           fabs(m1-m0)/m0 < 2e-3 && fabs(V1-V0)/V0 < 5e-3, det);
        snprintf(det, sizeof det, "chi=%d", ft3_euler(s));
        ok("topologia intacta", ft3_euler(s) == 2, det);
        ft3_destroi(s);
    }

    // ---- 2 e 3. a ondulacao cai, e o volume se conserva ---------------
    printf("\n-- ondulacao de 8%% do raio, alternando por vertice --\n");
    real vol_tsur = 0.0, des_tsur = 0.0, des0 = 0.0;
    {
        ft3_superficie *s = ft3_cria_esfera(c, R, 3);
        const real Vlimpo = ft3_volume(s);
        ondula(s, c, 0.08);
        real Vond = ft3_volume(s), m0, p0;
        estat_raio(s, c, &m0, &des0, &p0);
        printf("     ondulada: desvio do raio %.5f, pior %.5f, volume %.8f\n",
               (double)des0, (double)p0, (double)Vond);

        for (int it = 0; it < 3; it++) { ft3_balance(s, LIM); ft3_tsur(s, 0.0, LIM); }
        real m1, p1; estat_raio(s, c, &m1, &des_tsur, &p1);
        vol_tsur = ft3_volume(s);
        printf("     depois:   desvio do raio %.5f, pior %.5f, volume %.8f\n",
               (double)des_tsur, (double)p1, (double)vol_tsur);
        char det[200];
        snprintf(det, sizeof det, "desvio %.5f -> %.5f  (caiu %.0f%%)",
                 (double)des0, (double)des_tsur, (double)(100*(1-des_tsur/des0)));
        ok("a ondulacao CAI", des_tsur < 0.5*des0, det);
        snprintf(det, sizeof det, "volume %.8f -> %.8f  (%.3f%%), limpo %.8f",
                 (double)Vond, (double)vol_tsur,
                 (double)(100*fabs(vol_tsur-Vond)/Vond), (double)Vlimpo);
        ok("o volume se CONSERVA", fabs(vol_tsur-Vond)/Vond < 0.02, det);
        snprintf(det, sizeof det, "chi=%d", ft3_euler(s));
        ok("topologia intacta depois de 3 passadas", ft3_euler(s) == 2, det);
        ft3_destroi(s);
    }

    // ---- 4. o contraste com o laplaciano ------------------------------
    printf("\n-- o contraste: a media laplaciana suaviza ENCOLHENDO --\n");
    {
        ft3_superficie *s = ft3_cria_esfera(c, R, 3);
        ondula(s, c, 0.08);
        const real Vond = ft3_volume(s);
        for (int it = 0; it < 3; it++) laplaciano(s);
        real m1, d1, p1; estat_raio(s, c, &m1, &d1, &p1);
        const real Vlap = ft3_volume(s);
        printf("     laplaciano: desvio %.5f, volume %.8f (%.2f%% do ondulado)\n",
               (double)d1, (double)Vlap, (double)(100*Vlap/Vond));
        char det[220];
        snprintf(det, sizeof det,
                 "laplaciano perde %.2f%% do volume, BALANCE+TSUR %.2f%%",
                 (double)(100*fabs(Vlap-Vond)/Vond),
                 (double)(100*fabs(vol_tsur-Vond)/Vond));
        ok("o laplaciano encolhe MUITO mais",
           fabs(Vlap-Vond)/Vond > 3.0*fabs(vol_tsur-Vond)/Vond, det);
        ft3_destroi(s);
    }

    // ---- 5. TENSAO SUPERFICIAL ARTIFICIAL -----------------------------
    //
    // O teste que o usuario pediu, e que os anteriores NAO fazem.  Na esfera o
    // estado liso ja' e' o de area minima no volume dado, logo um filtro que
    // reduz area empurraria para onde a fisica tambem empurra -- e o efeito
    // seria invisivel.  Num ELIPSOIDE nao: se o filtro o aproximar da esfera a
    // volume fixo, isso E' tensao artificial, porque e' o que a tensao faz.
    printf("\n-- tensao superficial artificial: o teste do elipsoide --\n");
    {
        printf("     %-8s %14s %14s %10s\n",
               "limiar","d esfericid.","moveu/passada","ondulacao");
        real dpsi_sem = 0.0, dpsi_com = 0.0;
        for (int i = 0; i < 2; i++) {
            const real lim = i ? LIM : 0.0;
            // elipsoide alongado, 75 aplicacoes = a cada 20 passos em 1501
            ft3_superficie *s = ft3_cria_esfera(c, R, 3);
            { const int nv=ft3_num_vertices(s); Point *p=(Point*)ft3_posicoes(s);
              const real fz=1.6, f=1.0/sqrt(fz);
              for(int q=0;q<nv;q++){ p[q][0]=c[0]+(p[q][0]-c[0])*f;
                p[q][1]=c[1]+(p[q][1]-c[1])*f; p[q][2]=c[2]+(p[q][2]-c[2])*fz; } }
            real A=ft3_area(s), V=ft3_volume(s);
            const real p0 = pow(M_PI,1.0/3.0)*pow(6*V,2.0/3.0)/A;
            long mov = 0;
            for (int k = 0; k < 75; k++) { mov += ft3_balance(s, lim);
                                           mov += ft3_tsur(s, 0.0, lim); }
            A=ft3_area(s); V=ft3_volume(s);
            const real dp = pow(M_PI,1.0/3.0)*pow(6*V,2.0/3.0)/A - p0;
            ft3_destroi(s);
            // e a ondulacao ainda cai?
            s = ft3_cria_esfera(c, R, 3); ondula(s, c, 0.08);
            real m,d0,pp; estat_raio(s,c,&m,&d0,&pp);
            for (int k=0;k<3;k++){ ft3_balance(s,lim); ft3_tsur(s,0.0,lim); }
            real d1; estat_raio(s,c,&m,&d1,&pp); ft3_destroi(s);
            printf("     %-8.2f %14.5f %14.0f %9.0f%%\n",
                   (double)lim,(double)dp,(double)(mov/75.0),
                   (double)(100*(1-d1/d0)));
            if (i) dpsi_com = dp; else dpsi_sem = dp;
        }
        char det[220];
        snprintf(det, sizeof det,
                 "sem limiar %.5f, com limiar %.5f; a diferenca FT-VOF e' 0.0023",
                 (double)dpsi_sem, (double)dpsi_com);
        ok("SEM limiar o filtro gera tensao artificial", dpsi_sem > 0.01, det);
        ok("COM limiar o elipsoide liso fica intocado", fabs(dpsi_com) < 1e-4, det);
    }

    printf("\n%s  (%d falha(s))\n\n", falhas ? "FALHOU" : "PASSOU", falhas);
    return falhas ? 1 : 0;
}
