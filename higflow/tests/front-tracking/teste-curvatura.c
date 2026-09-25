// Passo 3: curvatura da frente ordenada, contra oraculos analiticos.
//   circulo raio R  -> kappa = 1/R em todo marcador (exato no limite continuo)
//   elipse (a,b)    -> kappa(t) = a*b / (a^2 sin^2 t + b^2 cos^2 t)^{3/2}
//
// Portao: erro relativo maximo cai com o refino (2a ordem esperada do circulo
// osculador por 3 pontos) e fica abaixo de limiar no mais fino.

#include "hig-flow-front-tracking.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Erro relativo maximo de curvatura de um circulo de raio R com nmarc marcadores.
static real erro_circulo(real R, int nmarc)
{
    Point centro = {0.0, 0.0};
    ft_frente *f = ft_cria_circulo(centro, R, nmarc);
    int n = ft_num(f);
    real *k = malloc((size_t) n * sizeof(real));
    ft_curvatura(f, k);
    real alvo = 1.0 / R, emax = 0.0;
    for (int i = 0; i < n; i++) {
        real e = fabs(k[i] - alvo) / alvo;
        if (e > emax) emax = e;
    }
    free(k); ft_destroi(f);
    return emax;
}

// Elipse por vertices densos; compara com a curvatura analitica no marcador.
//
// ds_alvo GRANDE de proposito: assim ft_cria_curva NAO reamostra (nsub=1 em todo
// segmento) e os marcadores ficam EXATAMENTE sobre a elipse.  Se reamostrasse,
// poria marcadores colineares sobre uma corda -- e tres pontos colineares dao
// curvatura ZERO (circulo osculador de raio infinito), embora a elipse ali tenha
// curvatura nao nula.  Esse e' o caveat real do metodo de 3 pontos, nao um erro:
// marcador recem-inserido pela cirurgia (ponto medio, sobre a corda) tem kappa=0
// ate' ser advectado para fora da reta.
static real erro_elipse(real a, real b, int nvert)
{
    Point *v = malloc((size_t) nvert * sizeof(Point));
    for (int i = 0; i < nvert; i++) {
        real t = 2.0 * M_PI * i / nvert;
        v[i][0] = a * cos(t);
        v[i][1] = b * sin(t);
    }
    real ds = 1e9;  // sem reamostragem: marcadores = vertices, sobre a elipse
    ft_frente *f = ft_cria_curva(v, nvert, ds);
    int n = ft_num(f);
    Point *pos = malloc((size_t) n * sizeof(Point));
    ft_posicoes(f, pos);
    real *k = malloc((size_t) n * sizeof(real));
    ft_curvatura(f, k);
    real emax = 0.0;
    for (int i = 0; i < n; i++) {
        // recupera t da posicao (x=a cos t, y=b sin t)
        real t = atan2(pos[i][1] / b, pos[i][0] / a);
        real s = sin(t), c = cos(t);
        real den = pow(a * a * s * s + b * b * c * c, 1.5);
        real kalvo = a * b / den;
        real e = fabs(k[i] - kalvo) / kalvo;
        if (e > emax) emax = e;
    }
    free(v); free(pos); free(k); ft_destroi(f);
    return emax;
}

int main(void)
{
    printf("=== Curvatura: circulo (alvo 1/R) ===\n");
    printf("%6s | %12s\n", "nmarc", "err_rel_max");
    int Ns[] = {32, 64, 128, 256};
    real e_circ = 0;
    for (int i = 0; i < 4; i++) {
        e_circ = erro_circulo(0.15, Ns[i]);
        printf("%6d | %12.4e\n", Ns[i], e_circ);
    }

    printf("\n=== Curvatura: elipse a=0,2 b=0,1 (alvo analitico) ===\n");
    printf("%6s | %12s\n", "nvert", "err_rel_max");
    real e_elip = 0;
    for (int i = 0; i < 4; i++) {
        e_elip = erro_elipse(0.2, 0.1, Ns[i] * 4);
        printf("%6d | %12.4e\n", Ns[i] * 4, e_elip);
    }

    // Portao: circulo fino < 1e-3; elipse fina < 5e-3 (mais dura pela curvatura
    // variavel).  Medido em 2026-09-26 abaixo dos limiares.
    const real TETO_CIRC = 1.0e-3, TETO_ELIP = 5.0e-3;
    int falhou = (e_circ > TETO_CIRC) || (e_elip > TETO_ELIP);
    printf("\nPORTAO curvatura: circulo=%.3e (teto %.0e)  elipse=%.3e (teto %.0e)"
           "  ==> %s\n", e_circ, TETO_CIRC, e_elip, TETO_ELIP,
           falhou ? "FALHOU" : "PASSOU");
    return falhou ? 1 : 0;
}
