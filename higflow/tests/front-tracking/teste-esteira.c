// Oraculo da PARTE PURA do sensor de esteira: |rot u| do gradiente.
//
// A funcao testada e' a MESMA que roda em producao -- vem de
// examples-common/criterio-esteira-pura.h, incluida pelos dois lados.  Nao ha'
// copia: a primeira versao disto copiava a funcao e trazia um roteiro para
// comparar os textos, e a guarda acusava diferenca de FORMATACAO como se fosse
// divergencia.  Guardar uma duplicacao e' pior que nao tê-la.
//
// Compila sem solver e sem MPI, que e' a razao de a conta estar separada da
// amostragem: a amostragem das facetas se exercita na corrida, a conta e' o que
// erra sinal ou troca eixo e se pega aqui contra forma fechada.
#include <stdio.h>
#include <string.h>
typedef double real;
#include "../../examples-common/criterio-esteira-pura.h"

static int falhas = 0;
static void ok(int c, const char *msg, double v, double e)
{
    printf("  [%s] %-44s |omega|=%.6f  exato %.6f\n", c ? "ok" : "FALHA", msg, v, e);
    if (!c) falhas++;
}

int main(void)
{
    real du[3][3];
    printf("\n-- |rot u| contra forma fechada --\n");

    // u = (g y, 0, 0)  ->  du[0][1] = g,  rot u = (0,0,-g),  |rot u| = |g|
    const real g1 = 3.7;
    memset(du, 0, sizeof du); du[0][1] = g1;
    ok(fabs(esteira_omega_grad(du)-g1) < 1e-12, "cisalhamento u=(g y,0,0)",
       esteira_omega_grad(du), g1);

    // v = (0, g z, 0)  ->  du[1][2] = g,  rot u = (-g,0,0)
    const real g2 = -2.3;
    memset(du, 0, sizeof du); du[1][2] = g2;
    ok(fabs(esteira_omega_grad(du)-fabs(g2)) < 1e-12, "cisalhamento v=(0,g z,0)",
       esteira_omega_grad(du), fabs(g2));

    // w = (0, 0, g x)  ->  du[2][0] = g,  rot u = (0,-g,0)
    const real g3 = 1.4;
    memset(du, 0, sizeof du); du[2][0] = g3;
    ok(fabs(esteira_omega_grad(du)-g3) < 1e-12, "cisalhamento w=(0,0,g x)",
       esteira_omega_grad(du), g3);

    // ROTACAO DE CORPO RIGIDO em torno de z: u=(-W y, W x, 0)
    // du[0][1]=-W, du[1][0]=+W  ->  wz = W-(-W) = 2W.  Os dois termos SOMAM, e
    // e' este caso que pega sinal trocado: nos cisalhamentos simples so' um
    // termo e' nao nulo e o modulo esconderia o erro.
    const real W = 1.9;
    memset(du, 0, sizeof du); du[0][1] = -W; du[1][0] = W;
    ok(fabs(esteira_omega_grad(du)-2.0*W) < 1e-12, "rotacao de corpo rigido (termos SOMAM)",
       esteira_omega_grad(du), 2.0*W);

    // DEFORMACAO PURA irrotacional: u=(g x, -g y, 0).  du diagonal, rot = 0.
    // Pega um sensor que respondesse a GRADIENTE em vez de rotacional -- e' o
    // erro mais plausivel de todos, porque |grad u| e' grande aqui.
    const real g4 = 2.6;
    memset(du, 0, sizeof du); du[0][0] = g4; du[1][1] = -g4;
    ok(fabs(esteira_omega_grad(du)) < 1e-14, "deformacao pura (irrotacional)",
       esteira_omega_grad(du), 0.0);

    // SIMETRICO puro: du[a][b] = du[b][a].  O rotacional ve' so' a parte
    // ANTISSIMETRICA, logo tem de dar zero por maior que seja o cisalhamento.
    memset(du, 0, sizeof du);
    du[0][1] = du[1][0] = 7.3;  du[0][2] = du[2][0] = -4.1;  du[1][2] = du[2][1] = 2.2;
    ok(fabs(esteira_omega_grad(du)) < 1e-13, "gradiente SIMETRICO (parte antissim. nula)",
       esteira_omega_grad(du), 0.0);

    // CAMPO UNIFORME: gradiente nulo.  Guarda de vacuidade.
    memset(du, 0, sizeof du);
    ok(fabs(esteira_omega_grad(du)) < 1e-14, "campo uniforme", esteira_omega_grad(du), 0.0);

    printf("\n%s  (%d falha(s))\n\n", falhas == 0 ? "PASSOU" : "FALHOU", falhas);
    return falhas != 0;
}
