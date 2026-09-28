// Passo 8: a funcao indicadora (fracao de area da frente por celula), que e' o
// insumo da formulacao BALANCEADA da tensao superficial.
//
// Oraculos:
//   (1) PARTICAO  somar ft_area_na_caixa sobre uma grade que cobre a frente tem
//       de dar ft_area -- exatamente, porque as caixas particionam o plano.
//   (2) CELULA CHEIA / VAZIA  caixa dentro da gota da' a area da caixa; caixa
//       fora da' zero.
//   (3) CONVERGENCIA  a soma nao depende da resolucao da grade (e' exata em
//       qualquer uma), e a area do poligono converge para pi R^2 com N.

#include "hig-flow-front-tracking.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Soma da indicadora sobre uma grade N x N cobrindo [0,1]^2.
static real soma_grade(const ft_frente *f, int N)
{
    real h = 1.0 / N, total = 0.0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            real lo[2] = { i * h, j * h };
            real hi[2] = { (i + 1) * h, (j + 1) * h };
            total += ft_area_na_caixa(f, lo, hi);
        }
    }
    return total;
}

int main(void)
{
    Point centro = {0.5, 0.5};
    const real R = 1.0 / 6.0;

    printf("=== Indicadora: particao da grade reproduz a area da frente ===\n");
    printf("%6s %6s | %14s %14s %12s\n",
           "nmarc", "N", "soma_grade", "ft_area", "erro_rel");
    real pior = 0.0;
    int Ns[]  = {16, 32, 61, 128};
    int Nm[]  = {32, 64, 128, 256};
    for (int a = 0; a < 4; a++) {
        ft_frente *f = ft_cria_circulo(centro, R, Nm[a]);
        real area = ft_area(f);
        for (int b = 0; b < 4; b++) {
            real s = soma_grade(f, Ns[b]);
            real e = fabs(s - area) / area;
            if (e > pior) pior = e;
            if (b == 2)  // imprime uma grade por frente, para nao poluir
                printf("%6d %6d | %14.10f %14.10f %12.3e\n",
                       Nm[a], Ns[b], (double) s, (double) area, (double) e);
        }
        ft_destroi(f);
    }
    printf("  pior erro sobre 4 frentes x 4 grades: %.3e\n", (double) pior);

    // (2) celula cheia e celula vazia
    printf("\n=== Celula cheia / vazia ===\n");
    ft_frente *f = ft_cria_circulo(centro, R, 256);
    real h = 1.0 / 61.0;
    real lo_in[2]  = {0.5 - 0.5*h, 0.5 - 0.5*h};      // no centro da gota
    real hi_in[2]  = {0.5 + 0.5*h, 0.5 + 0.5*h};
    real lo_out[2] = {0.05, 0.05};                     // canto do dominio
    real hi_out[2] = {0.05 + h, 0.05 + h};
    real a_in  = ft_area_na_caixa(f, lo_in,  hi_in);
    real a_out = ft_area_na_caixa(f, lo_out, hi_out);
    real a_cel = h * h;
    printf("  dentro: %.12e  (area da celula %.12e)  erro=%.3e\n",
           (double) a_in, (double) a_cel,
           (double) fabs(a_in - a_cel) / a_cel);
    printf("  fora:   %.12e  (esperado 0)\n", (double) a_out);

    // (3) a area do poligono converge para o circulo
    printf("\n=== Area do poligono -> pi R^2 ===\n");
    real exato = M_PI * R * R;
    for (int a = 0; a < 4; a++) {
        ft_frente *g = ft_cria_circulo(centro, R, Nm[a]);
        printf("  N=%4d  area=%.10f  erro_rel=%.3e\n",
               Nm[a], (double) ft_area(g),
               (double) fabs(ft_area(g) - exato) / exato);
        ft_destroi(g);
    }

    const real TETO_PART = 1.0e-12, TETO_CEL = 1.0e-12;
    int falhou = (pior > TETO_PART)
              || (fabs(a_in - a_cel) / a_cel > TETO_CEL)
              || (a_out != 0.0);
    printf("\nPORTAO indicadora: particao=%.3e (teto %.0e)  cheia=%.3e  "
           "vazia=%.1e  ==> %s\n",
           (double) pior, TETO_PART,
           (double) (fabs(a_in - a_cel) / a_cel), (double) a_out,
           falhou ? "FALHOU" : "PASSOU");
    ft_destroi(f);
    return falhou ? 1 : 0;
}
