// Arreio B1: teste reversivel de Rider-Kothe (vortice unico).
//
// Campo (funcao de corrente psi = (1/pi) sin^2(pi x) sin^2(pi y)), divergencia
// nula, multiplicado por cos(pi t / T) para reverter em t=T/2 e VOLTAR a' forma
// inicial em t=T:
//     u =  sin^2(pi x) sin(2 pi y) cos(pi t/T)
//     v = -sin(2 pi x) sin^2(pi y) cos(pi t/T)
//
// Escoamento incompressivel PRESERVA a area fechada em todo instante; a
// reversibilidade acrescenta que a FORMA volta em t=T.  Dois oraculos:
//   (1) conservacao: max |A(t)-A0|/A0 ao longo da corrida  (deve -> 0 com refino)
//   (2) reversao:    |A(T)-A0|/A0                            (idem)

#include "hig-flow-front-tracking.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { real T; } vortex_ctx;

static void campo_rider_kothe(const Point x, real t, void *vctx, real u[DIM])
{
    vortex_ctx *c = (vortex_ctx *) vctx;
    real g = cos(M_PI * t / c->T);
    real sx = sin(M_PI * x[0]), sy = sin(M_PI * x[1]);
    u[0] =  sx * sx * sin(2.0 * M_PI * x[1]) * g;
    u[1] = -sin(2.0 * M_PI * x[0]) * sy * sy * g;
}

// Roda uma corrida e devolve (err_max_conservacao, err_reversao, n_final).
static void roda(int nmarc, real dt, real T, int cirurgia, int escreve_vtk,
                 real *err_cons, real *err_rev, int *n_final)
{
    Point centro = {0.5, 0.75};
    ft_frente *f = ft_cria_circulo(centro, 0.15, nmarc);
    vortex_ctx ctx = { T };

    real A0 = ft_area(f);
    real emax = 0.0;
    int nsteps = (int) llround(T / dt);
    int quadro = 0;
    if (escreve_vtk) ft_escreve_vtk(f, "scratchpad/rk", quadro++);

    for (int s = 0; s < nsteps; s++) {
        real t = s * dt;
        ft_advecta(f, campo_rider_kothe, &ctx, t, dt);
        if (cirurgia) ft_cirurgia(f);
        real e = fabs(ft_area(f) - A0) / A0;
        if (e > emax) emax = e;
        if (escreve_vtk && (s % (nsteps / 20 + 1) == 0))
            ft_escreve_vtk(f, "scratchpad/rk", quadro++);
    }
    if (escreve_vtk) ft_escreve_vtk(f, "scratchpad/rk", quadro++);

    *err_cons = emax;
    *err_rev  = fabs(ft_area(f) - A0) / A0;
    *n_final  = ft_num(f);
    ft_destroi(f);
}

int main(void)
{
    const real T = 4.0;

    printf("=== Rider-Kothe reversivel, T=%.1f ===\n\n", T);

    printf("--- SEM cirurgia (mostra a perda de resolucao no filamento) ---\n");
    printf("%6s %8s | %12s %12s %8s\n",
           "nmarc", "dt", "cons_max", "reversao", "n_final");
    int Ns[]   = {64, 128, 256};
    real dts[] = {0.02, 0.01, 0.005};
    for (int i = 0; i < 3; i++) {
        real ec, er; int nf;
        roda(Ns[i], dts[i], T, 0, 0, &ec, &er, &nf);
        printf("%6d %8.4f | %12.4e %12.4e %8d\n", Ns[i], dts[i], ec, er, nf);
    }

    printf("\n--- COM cirurgia (Delta s ~ ds_alvo mantido) ---\n");
    printf("%6s %8s | %12s %12s %8s\n",
           "nmarc", "dt", "cons_max", "reversao", "n_final");
    real ec_fino = 0, er_fino = 0; int nf_fino = 0;
    for (int i = 0; i < 3; i++) {
        real ec, er; int nf;
        roda(Ns[i], dts[i], T, 1, 0, &ec, &er, &nf);
        printf("%6d %8.4f | %12.4e %12.4e %8d\n", Ns[i], dts[i], ec, er, nf);
        ec_fino = ec; er_fino = er; nf_fino = nf;  // guarda o mais fino
    }
    (void) nf_fino;

    // PORTAO: no mais fino (256, dt=0,005) com cirurgia, a conservacao de area
    // ao longo da corrida e a reversao em t=T tem de ficar abaixo dos limiares.
    // Medido em 2026-09-26: cons_max=9,9e-5, reversao=3,6e-5 -- os limiares dao
    // margem folgada e pegam regressao de sinal, escala ou ordem.
    const real TETO_CONS = 2.0e-4;
    const real TETO_REV  = 1.0e-4;
    int falhou = (ec_fino > TETO_CONS) || (er_fino > TETO_REV);
    printf("\nPORTAO B1 (256, dt=0,005): cons_max=%.3e (teto %.0e)  "
           "reversao=%.3e (teto %.0e)  ==> %s\n",
           ec_fino, TETO_CONS, er_fino, TETO_REV, falhou ? "FALHOU" : "PASSOU");
    return falhou ? 1 : 0;
}
