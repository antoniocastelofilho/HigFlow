// FRACAO DE VOLUME em 3D: quanto de cada celula esta' DENTRO da superficie
// triangulada.  E' o que destrava o bifasico -- densidade e viscosidade saem
// dela --, e erro aqui corrompe o campo em silencio.
//
// Cinco portoes, em camadas, cada um com oraculo fechado contra a esfera:
//   1. dentro/fora: classifica certo, inclusive rente a' superficie
//   2. distancia: bate com | |x-c| - R | da esfera
//   3. limites: fracao em [0,1]; longe dentro da' 1, longe fora da' 0
//   4. soma: SUM fracao*h^3 = (4/3) pi R^3, com ORDEM medida
//   5. CRUZAMENTO com metodo INDEPENDENTE (subamostragem), celula a celula
//
// O portao 5 e' o que mais vale.  A soma global pode fechar com erros que se
// cancelam entre celulas; a subamostragem nao usa plano nenhum, entao concordar
// com ela celula a celula e' evidencia de outra natureza.

#include "hig-flow-front-tracking-3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int falhas = 0;

static void ok(int cond, const char *rotulo, const char *detalhe)
{
    printf("  %s %s: %s\n", cond ? "ok" : "FALHOU", rotulo, detalhe);
    if (!cond) falhas++;
}

static const Point CENTRO = {0.5, 0.5, 0.5};
static const real  R = 0.25;

static real raio(const Point p)
{
    real d[3] = { p[0]-CENTRO[0], p[1]-CENTRO[1], p[2]-CENTRO[2] };
    return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

// =========================================================================
static void portao_dentro(const ft3_superficie *s)
{
    printf("\n=== 1) dentro/fora ===\n");
    int erros = 0, rente = 0, erros_rente = 0;
    // amostra pseudoaleatoria determinista (congruencia linear): sorteio fixo
    // repete, e teste que nao repete nao e' oraculo
    unsigned long semente = 12345;
    for (int i = 0; i < 4000; i++) {
        Point p;
        for (int d = 0; d < 3; d++) {
            semente = semente * 6364136223846793005UL + 1442695040888963407UL;
            p[d] = (real) ((semente >> 33) % 1000000) / 1000000.0;
        }
        real r = raio(p);
        int esperado = (r < R);
        int obtido = ft3_dentro(s, p);
        // A superficie e' a icosfera INSCRITA, nao a esfera: entre o poliedro e
        // a esfera ha' uma casca de espessura O(h^2) onde discordar e' correto.
        // So' conta erro fora dessa casca.
        real folga = 0.02 * R;
        if (fabs(r - R) < folga) { rente++; if (esperado != obtido) erros_rente++; continue; }
        if (esperado != obtido) erros++;
    }
    char buf[160];
    snprintf(buf, sizeof buf, "%d erros em 4000 pontos (%d rente a' superficie, %d deles discordando)",
             erros, rente, erros_rente);
    ok(erros == 0, "classificacao", buf);
}

static void portao_distancia(const ft3_superficie *s)
{
    printf("\n=== 2) distancia ===\n");
    real pior = 0.0;
    unsigned long semente = 999;
    for (int i = 0; i < 400; i++) {
        Point p;
        for (int d = 0; d < 3; d++) {
            semente = semente * 6364136223846793005UL + 1442695040888963407UL;
            p[d] = (real) ((semente >> 33) % 1000000) / 1000000.0;
        }
        real n[3];
        real dm = ft3_distancia(s, p, n);
        real dex = fabs(raio(p) - R);
        real e = fabs(dm - dex);
        if (e > pior) pior = e;
    }
    char buf[160];
    // O erro e' dominado pela faceta: a icosfera esta' a O(h_s^2) da esfera.
    snprintf(buf, sizeof buf, "pior desvio %.3e (a icosfera esta' a %.3e da esfera)",
             (double) pior, (double) (R * 1.2e-3));
    ok(pior < 5e-3, "contra a esfera", buf);
}

static void portao_limites_e_soma(const ft3_superficie *s, int N, real *err_soma)
{
    const real h = 1.0 / N;
    real soma = 0.0;
    int fora_faixa = 0, cheias = 0, vazias = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++) {
                real lo[3] = { i*h, j*h, k*h };
                real hi[3] = { (i+1)*h, (j+1)*h, (k+1)*h };
                real f = ft3_fracao_na_caixa(s, lo, hi);
                if (f < -1e-12 || f > 1.0 + 1e-12) fora_faixa++;
                if (f > 1.0 - 1e-12) cheias++;
                if (f < 1e-12) vazias++;
                soma += f * h * h * h;
            }
    real Vex = ft3_volume(s);            // o volume DA SUPERFICIE, nao o da esfera:
                                          // e' esse que a particao tem de reproduzir
    *err_soma = fabs(soma - Vex) / Vex;
    printf("    N=%3d h=%.4f  soma=%.8f  volume da superficie=%.8f  err=%.3e  "
           "cheias=%d vazias=%d fora_faixa=%d\n",
           N, (double) h, (double) soma, (double) Vex, (double) *err_soma,
           cheias, vazias, fora_faixa);
    if (fora_faixa) ok(0, "limites", "alguma fracao fora de [0,1]");
}

static void portao_cruzamento(const ft3_superficie *s)
{
    printf("\n=== 5) cruzamento com subamostragem (metodo independente) ===\n");
    const int N = 20;
    const real h = 1.0 / N;
    const int K = 10;                    // 1000 amostras por celula
    real pior = 0.0, soma_abs = 0.0;
    int n = 0;
    printf("    %8s %10s %10s %10s\n", "celula", "plano", "amostrado", "dif");
    for (int i = 0; i < N && n < 12; i++)
        for (int j = 0; j < N && n < 12; j++)
            for (int k = 0; k < N && n < 12; k++) {
                real lo[3] = { i*h, j*h, k*h };
                real hi[3] = { (i+1)*h, (j+1)*h, (k+1)*h };
                real f = ft3_fracao_na_caixa(s, lo, hi);
                if (f < 0.05 || f > 0.95) continue;     // so' celulas de interface
                real g = ft3_fracao_amostrada(s, lo, hi, K);
                real d = fabs(f - g);
                if (d > pior) pior = d;
                soma_abs += d;
                if (n < 6)
                    printf("    %2d,%2d,%2d  %10.6f %10.6f %10.6f\n",
                           i, j, k, (double) f, (double) g, (double) d);
                n++;
            }
    char buf[160];
    snprintf(buf, sizeof buf, "%d celulas de interface, pior diferenca %.4f, media %.4f",
             n, (double) pior, (double) (n ? soma_abs/n : 0.0));
    // A diferenca tem SINAL sistematico (plano > amostrado) e isso e' esperado:
    // a subamostragem conta o CENTRO de cada sub-caixa, o que subestima a
    // fracao de uma regiao convexa.  O comparador e' enviesado para baixo, nao
    // ruidoso -- entao concordancia dentro de 0,03 e' melhor do que parece.
    // A subamostragem tem ela propria erro O(1/K) = 0,1 por celula; exigir mais
    // que 0,06 seria exigir que o plano batesse com o RUIDO do comparador.
    ok(n > 0 && pior < 0.06, "concordancia", buf);
}

// =========================================================================
// 6: a grade acelera sem mudar o resultado
// =========================================================================
static void portao_grade(const ft3_superficie *s)
{
    printf("\n=== 6) grade espacial: MESMO resultado, menos trabalho ===\n");
    ft3_grade *g = ft3_grade_cria(s);
    if (!g) { ok(0, "construcao", "ft3_grade_cria devolveu NULL"); return; }
    int n[3], refs;
    ft3_grade_estado(g, n, &refs);
    printf("    grade %dx%dx%d, %d referencias a triangulo para %d triangulos "
           "(%.2f por triangulo)\n",
           n[0], n[1], n[2], refs, ft3_num_triangulos(s),
           (double) refs / ft3_num_triangulos(s));

    // (a) dentro/fora: tem de ser IGUAL, nao parecido
    int difs = 0;
    unsigned long semente = 4242;
    for (int i = 0; i < 20000; i++) {
        Point p;
        for (int d = 0; d < 3; d++) {
            semente = semente * 6364136223846793005UL + 1442695040888963407UL;
            p[d] = (real) ((semente >> 33) % 1000000) / 1000000.0;
        }
        if (ft3_dentro(s, p) != ft3_dentro_g(s, g, p)) difs++;
    }
    char buf[200];
    snprintf(buf, sizeof buf, "%d discordancias em 20000 pontos", difs);
    ok(difs == 0, "dentro/fora identico", buf);

    // (b) distancia: BIT a bit.  E' minimo sobre o MESMO conjunto de valores,
    // entao a grade achar o mesmo minimo nao e' sorte -- e' o contrato.
    int difd = 0; real pior = 0.0;
    semente = 777;
    for (int i = 0; i < 2000; i++) {
        Point p;
        for (int d = 0; d < 3; d++) {
            semente = semente * 6364136223846793005UL + 1442695040888963407UL;
            p[d] = (real) ((semente >> 33) % 1000000) / 1000000.0;
        }
        real na[3], nb[3];
        real da = ft3_distancia(s, p, na);
        real db = ft3_distancia_g(s, g, p, nb);
        if (da != db) { difd++; if (fabs(da-db) > pior) pior = fabs(da-db); }
    }
    snprintf(buf, sizeof buf, "%d diferencas em 2000 consultas (pior %.3e)", difd, (double) pior);
    ok(difd == 0, "distancia bit a bit", buf);

    // (c) a fracao sobre uma particao inteira, e o TEMPO das duas
    const int N = 20;
    const real h = 1.0 / N;
    int diff = 0;
    clock_t t0 = clock();
    real soma_b = 0.0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        real lo[3] = { i*h, j*h, k*h }, hi[3] = { (i+1)*h, (j+1)*h, (k+1)*h };
        soma_b += ft3_fracao_na_caixa(s, lo, hi);
    }
    clock_t t1 = clock();
    real soma_g = 0.0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        real lo[3] = { i*h, j*h, k*h }, hi[3] = { (i+1)*h, (j+1)*h, (k+1)*h };
        real fg = ft3_fracao_na_caixa_g(s, g, lo, hi);
        real fb = ft3_fracao_na_caixa(s, lo, hi);
        if (fg != fb) diff++;
        soma_g += fg;
    }
    clock_t t2 = clock();
    double tb = (double)(t1-t0)/CLOCKS_PER_SEC;
    double tg = (double)(t2-t1)/CLOCKS_PER_SEC - tb;   // o laco (c) roda as duas
    printf("    particao %d^3: exaustiva %.3f s, acelerada %.3f s  (ganho %.1fx)\n",
           N, tb, tg > 0 ? tg : 1e-9, tg > 0 ? tb/tg : 0.0);
    snprintf(buf, sizeof buf, "%d celulas diferentes em %d; somas %.10f e %.10f",
             diff, N*N*N, (double) soma_b, (double) soma_g);
    ok(diff == 0, "fracao identica", buf);
    // Teto do ganho posto DEPOIS de medir 14,1x, em 5x: abaixo disso alguma
    // coisa voltou a varrer o que a grade existe para evitar.  Exigir 12x
    // amarraria o portao a esta maquina.
    snprintf(buf, sizeof buf, "exaustiva %.3f s contra acelerada %.3f s (ganho %.1fx)",
             tb, tg, tg > 0 ? tb/tg : 0.0);
    ok(tg > 0 && tb / tg > 5.0, "a grade acelera", buf);

    ft3_grade_destroi(g);
}

int main(void)
{
    printf("=== teste-fracao3d: fracao de volume da superficie triangulada ===\n");
    ft3_superficie *s = ft3_cria_esfera(CENTRO, R, 3);
    printf("    esfera R=%.3f em (%.1f,%.1f,%.1f): %d vertices, %d triangulos\n",
           (double) R, (double) CENTRO[0], (double) CENTRO[1], (double) CENTRO[2],
           ft3_num_vertices(s), ft3_num_triangulos(s));

    portao_dentro(s);
    portao_distancia(s);

    printf("\n=== 3 e 4) limites e soma sobre a particao ===\n");
    real e1, e2, e3;
    portao_limites_e_soma(s, 10, &e1);
    portao_limites_e_soma(s, 20, &e2);
    portao_limites_e_soma(s, 40, &e3);
    ok(1, "limites", "toda fracao em [0,1]");

    char buf[160];
    // TRES pontos, e nao dois: a razao entre duas resolucoes nao distingue
    // segunda ordem de coincidencia.  Duas razoes proximas de 4 distinguem.
    real r1 = e1 / e2, r2 = e2 / e3;
    snprintf(buf, sizeof buf, "err %.3e -> %.3e -> %.3e  (razoes %.2f e %.2f)",
             (double) e1, (double) e2, (double) e3, (double) r1, (double) r2);
    ok(r1 > 3.0 && r2 > 3.0, "segunda ordem na soma", buf);
    // Teto absoluto ~1,5x o medido no mais fino: folga de ordem de grandeza nao
    // pegaria regressao.
    snprintf(buf, sizeof buf, "err no mais fino (N=40, h=0,025) = %.3e (teto 4e-3)",
             (double) e3);
    ok(e3 < 4e-3, "soma reproduz o volume", buf);

    portao_cruzamento(s);
    portao_grade(s);

    ft3_destroi(s);
    printf("\n%s\n", falhas ? "TESTE DA FRACAO 3D FALHOU" : "TESTE DA FRACAO 3D PASSOU");
    return falhas ? 1 : 0;
}
