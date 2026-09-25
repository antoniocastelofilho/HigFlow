// Passo 4: forca de tensao superficial sigma*kappa*n espalhada pelo kernel de
// Roma numa malha uniforme.  Dois oraculos:
//   (1) CONSERVACAO  SUM_grade F*h^2 == sigma * integral do vetor curvatura,
//       porque o nucleo soma 1.  Precisao de maquina (suporte nao vaza).
//   (2) EQUILIBRIO   sigma * integral do vetor curvatura -> 0 para curva
//       fechada, e cai com o refino da frente.
//
// Um circulo em repouso: a forca liquida tem de ser nula (Laplace), e cada
// pedaco puxa para o centro com sigma/R -- e' o que o B2 acoplado vai impor
// como salto de pressao.

#include "hig-flow-front-tracking.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const real sigma = 0.7;
    const real R = 0.15;
    Point centro = {0.5, 0.5};

    // Malha uniforme cobrindo [0,1]^2.
    const int N = 64;
    const real h = 1.0 / N, ox = 0.0, oy = 0.0;
    real *fx = calloc((size_t) N * N, sizeof(real));
    real *fy = calloc((size_t) N * N, sizeof(real));

    printf("=== Tensao superficial: circulo R=%.2f, sigma=%.2f, malha %dx%d ===\n",
           R, sigma, N, N);
    printf("%6s | %14s %14s | %14s %14s\n",
           "nmarc", "conserv_x", "conserv_y", "equilib_|F|", "esperado~0");

    real eq_fino = 0, cons_fino = 0;
    int Ns[] = {64, 128, 256, 512};
    for (int t = 0; t < 4; t++) {
        int nmarc = Ns[t];
        ft_frente *f = ft_cria_circulo(centro, R, nmarc);

        // Lado dos marcadores: integral do vetor curvatura (sem sigma).
        real integ[DIM];
        ft_integral_curvatura(f, integ);
        real marker_x = sigma * integ[0], marker_y = sigma * integ[1];

        // Espalha na grade (zera antes).
        for (int i = 0; i < N * N; i++) { fx[i] = 0; fy[i] = 0; }
        ft_espalha_tensao(f, sigma, ox, oy, h, N, N, fx, fy);

        // Soma da grade * h^2.
        real gx = 0, gy = 0;
        for (int i = 0; i < N * N; i++) { gx += fx[i]; gy += fy[i]; }
        gx *= h * h; gy *= h * h;

        real cons_x = fabs(gx - marker_x);
        real cons_y = fabs(gy - marker_y);
        real equilib = sqrt(marker_x * marker_x + marker_y * marker_y);

        printf("%6d | %14.4e %14.4e | %14.4e %14.4e\n",
               nmarc, cons_x, cons_y, equilib, 0.0);

        cons_fino = cons_x > cons_y ? cons_x : cons_y;
        eq_fino = equilib;
        ft_destroi(f);
    }

    free(fx); free(fy);

    // Portao: conservacao no nivel de maquina (< 1e-12) e equilibrio pequeno e
    // caindo (< 1e-3 no mais fino).  A conservacao e' o teste do nucleo; o
    // equilibrio e' o teste da curvatura fechar em zero ao redor do circulo.
    const real TETO_CONS = 1.0e-12, TETO_EQ = 1.0e-3;
    int falhou = (cons_fino > TETO_CONS) || (eq_fino > TETO_EQ);
    printf("\nPORTAO tensao: conserv=%.3e (teto %.0e)  equilib=%.3e (teto %.0e)"
           "  ==> %s\n", cons_fino, TETO_CONS, eq_fino, TETO_EQ,
           falhou ? "FALHOU" : "PASSOU");
    return falhou ? 1 : 0;
}
