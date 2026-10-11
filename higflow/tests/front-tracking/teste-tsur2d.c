// Arreio do TSUR 2D -- remocao de ondulacao submalha na polilinha.
//
// MESMO DESENHO DO teste-tsur.c (3D), e pela mesma razao: o oraculo que importa
// NAO e' "suavizou".  Qualquer media suaviza; a questao e' se suaviza SEM
// ENCOLHER.  Entao o arreio:
//
//   1. a circunferencia LIMPA e' ponto fixo -- aplicar nao a deforma
//   2. com ondulacao radial, a amplitude CAI
//   3. a AREA se conserva enquanto a amplitude cai
//   4. a media LAPLACIANA na MESMA curva ondulada suaviza tambem, mas
//      ENCOLHENDO -- sem este contraste, "conserva area" nao diz nada
//   5. o LIMIAR deixa a ELIPSE LISA intocada.  E' o teste que separa filtro de
//      tensao superficial artificial: na circunferencia o efeito e' invisivel
//      (o estado liso ja' e' perimetro minimo), e so' a elipse expoe.
//
// GUARDA DE VACUIDADE: cada bloco conta quantos casos exercitou e falha se for
// zero.  Arreio que nao toca o caso em teste imprime "0 falhas", que e'
// indistinguivel de aprovacao.

#include "hig-flow-front-tracking.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int falhas = 0, casos = 0;

static void ok(const char *nome, int cond, const char *detalhe)
{
    casos++;
    if (!cond) { falhas++; printf("  FALHOU  %-44s %s\n", nome, detalhe); }
    else       { printf("  ok      %-44s %s\n", nome, detalhe); }
}

// amplitude radial em torno do raio medio
static real amplitude(const ft_frente *f, const real c[2])
{
    const int n = ft_num(f);
    Point *p = (Point *) malloc((size_t) n * sizeof(Point));
    ft_posicoes(f, p);
    real med = 0.0;
    for (int i = 0; i < n; i++)
        med += sqrt((p[i][0]-c[0])*(p[i][0]-c[0]) + (p[i][1]-c[1])*(p[i][1]-c[1]));
    med /= n;
    real amp = 0.0;
    for (int i = 0; i < n; i++) {
        const real r = sqrt((p[i][0]-c[0])*(p[i][0]-c[0]) + (p[i][1]-c[1])*(p[i][1]-c[1]));
        const real d = fabs(r - med);
        if (d > amp) amp = d;
    }
    free(p);
    return amp;
}

// circunferencia com ondulacao radial de modo `modo` e amplitude `eps`
static ft_frente *ondulada(const real c[2], real R, int n, int modo, real eps)
{
    Point *v = (Point *) malloc((size_t) n * sizeof(Point));
    for (int i = 0; i < n; i++) {
        const real th = 2.0 * M_PI * i / n;
        const real r  = R * (1.0 + eps * cos(modo * th));
        v[i][0] = c[0] + r * cos(th);
        v[i][1] = c[1] + r * sin(th);
    }
    ft_frente *f = ft_cria_curva(v, n, 2.0 * M_PI * R / n * 1.0001);
    free(v);
    return f;
}

int main(void)
{
    const real c[2] = { 0.0, 0.0 };
    printf("=== TSUR 2D ===\n");

    // --- 1. a circunferencia limpa e' ponto fixo ---------------------------
    {
        ft_frente *f = ft_cria_circulo(c, 1.0, 120);
        const real A0 = ft_area(f);
        const int mov = ft_tsur(f, 0.0, 0.02);
        const real A1 = ft_area(f);
        char d[128];
        snprintf(d, sizeof d, "movidos=%d  dA/A=%.2e", mov, fabs(A1-A0)/A0);
        ok("circunferencia limpa: nenhum marcador move", mov == 0, d);
        ft_destroi(f);
    }

    // --- 2 e 3. o modo ALTERNADO cai, e a area se conserva -----------------
    //
    // O modo e' n/2: dois marcadores por onda, o zigue-zague que alterna marcador
    // a marcador.  E' o caso para o qual o TSUR existe -- o 3D testa o mesmo
    // ("alternando por vertice").  Testar um modo RESOLVIDO aqui seria pedir ao
    // filtro que apagasse forma legitima, e foi o engano da primeira versao
    // deste arreio: no modo 40 de 160 marcadores (quatro por onda) a queda e' de
    // 44% em 30 passadas, e deve mesmo ser -- ver o caso 6.
    real dA_tsur = 0.0;
    {
        ft_frente *f = ondulada(c, 1.0, 160, 80, 0.01);
        const real a0 = amplitude(f, c), A0 = ft_area(f);
        int mov = 0;
        for (int k = 0; k < 10; k++) mov += ft_tsur(f, 0.0, 0.0);
        const real a1 = amplitude(f, c), A1 = ft_area(f);
        dA_tsur = fabs(A1 - A0) / A0;
        char d[180];
        snprintf(d, sizeof d, "amp %.3e -> %.3e (queda %.1f%%), movidos=%d",
                 a0, a1, 100.0*(1.0 - a1/a0), mov);
        ok("modo alternado cai mais de 90% em 10 passadas", a1 < 0.10 * a0, d);
        snprintf(d, sizeof d, "dA/A = %.3e", dA_tsur);
        ok("area conservada enquanto a ondulacao cai", dA_tsur < 1.0e-3, d);
        ft_destroi(f);
    }

    // --- 4. o contraste: a laplaciana suaviza ENCOLHENDO -------------------
    {
        ft_frente *f = ondulada(c, 1.0, 160, 80, 0.01);
        const real A0 = ft_area(f);
        const int n = ft_num(f);
        Point *p = (Point *) malloc((size_t) n * sizeof(Point));
        Point *q = (Point *) malloc((size_t) n * sizeof(Point));
        ft_posicoes(f, p);
        for (int k = 0; k < 10; k++) {
            for (int i = 0; i < n; i++) {
                const real *a = p[(i-1+n)%n], *b = p[(i+1)%n];
                q[i][0] = 0.5*(a[0]+b[0]); q[i][1] = 0.5*(a[1]+b[1]);
            }
            for (int i = 0; i < n; i++) { p[i][0]=q[i][0]; p[i][1]=q[i][1]; }
        }
        real soma = 0.0;
        for (int i = 0; i < n; i++) {
            const real *a = p[i], *b = p[(i+1)%n];
            soma += a[0]*b[1] - b[0]*a[1];
        }
        const real A1 = 0.5 * fabs(soma);
        const real enc = (A0 - A1) / A0;
        char d[160];
        snprintf(d, sizeof d, "laplaciana encolhe %.3e; TSUR %.3e -- razao %.0fx",
                 enc, dA_tsur, dA_tsur > 0 ? enc/dA_tsur : 0.0);
        ok("contraste: a laplaciana encolhe mais que o TSUR", enc > 10.0 * dA_tsur, d);
        free(p); free(q); ft_destroi(f);
    }

    // --- 5. o limiar deixa a ELIPSE LISA intocada --------------------------
    {
        const int n = 160;
        Point *v = (Point *) malloc((size_t) n * sizeof(Point));
        const real A = 1.3, B = 0.7;
        for (int i = 0; i < n; i++) {
            const real th = 2.0 * M_PI * i / n;
            v[i][0] = A * cos(th); v[i][1] = B * sin(th);
        }
        ft_frente *f = ft_cria_curva(v, n, 0.05);
        free(v);
        real a0, b0, a1, b1;
        ft_semieixos(f, &a0, &b0);
        const real A0 = ft_area(f);
        int mov = 0;
        for (int k = 0; k < 20; k++) mov += ft_tsur(f, 0.0, 0.02);
        ft_semieixos(f, &a1, &b1);
        const real A1 = ft_area(f);
        char d[180];
        snprintf(d, sizeof d, "movidos=%d em 20 passadas, dA/A=%.2e, semieixos %.4f/%.4f -> %.4f/%.4f",
                 mov, fabs(A1-A0)/A0, a0, b0, a1, b1);
        ok("limiar 0,02: elipse lisa intocada", mov == 0, d);
        ft_destroi(f);
    }

    // --- 6. SELETIVIDADE: o modo RESOLVIDO sobrevive ----------------------
    //
    // Sem este caso, "o filtro remove ondulacao" nao distingue um filtro de uma
    // borracha.  Quatro marcadores por onda e' forma que a discretizacao
    // REPRESENTA, e apaga-la seria perder fisica, nao ruido.
    {
        ft_frente *f = ondulada(c, 1.0, 160, 40, 0.01);
        const real a0 = amplitude(f, c);
        for (int k = 0; k < 10; k++) ft_tsur(f, 0.0, 0.0);
        const real a1 = amplitude(f, c);
        char d[180];
        snprintf(d, sizeof d, "modo resolvido (4 marcadores/onda): amp %.3e -> %.3e, sobrou %.0f%%",
                 a0, a1, 100.0*a1/a0);
        ok("seletivo: o modo resolvido NAO e' apagado", a1 > 0.5 * a0, d);
        ft_destroi(f);
    }

    // --- guarda de vacuidade ----------------------------------------------
    printf("--- %d casos, %d falhas ---\n", casos, falhas);
    if (casos < 6) { printf("GUARDA DE VACUIDADE: so' %d casos exercitados\n", casos); return 2; }
    return falhas ? 1 : 0;
}
