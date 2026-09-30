// FASE B5, primeiro degrau: geometria e forca da superficie triangulada,
// verificadas contra a ESFERA ANALITICA, sem solver.
//
// Cinco portoes, todos com oraculo fechado:
//   1. area -> 4 pi R^2 e volume -> (4/3) pi R^3, com ORDEM medida
//   2. topologia: V - A + F = 2 em toda resolucao
//   3. forca: soma exatamente zero (estrutural) e |F|/area -> 2 sigma / R
//   4. cirurgia: preserva topologia, nao inverte triangulo, equaliza aresta
//   5. persistencia: ida-e-volta bit a bit COM conectividade, e recusa lixo
//
// O 2 SOBRE R DO PORTAO 3 NAO E' ERRO DE DIGITACAO.  Em 2D o salto de Laplace
// e' sigma/R; em 3D e' 2*sigma/R.  O caso 2D deste repositorio fechou a 0,02%
// contra sigma/R, e copiar aquele numero para ca' validaria o errado.

#include "hig-flow-front-tracking-3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int falhas = 0;

static void ok(int cond, const char *rotulo, const char *detalhe)
{
    printf("  %s %s: %s\n", cond ? "ok" : "FALHOU", rotulo, detalhe);
    if (!cond) falhas++;
}

// =========================================================================
// 1 e 2: geometria e topologia
// =========================================================================
static void portao_geometria(void)
{
    const Point centro = {0.3, -0.2, 0.7};       // longe da origem de proposito:
    const real  R = 0.45;                        // pega erro de origem no volume
    const real  A_ex = 4.0 * M_PI * R * R;
    const real  V_ex = 4.0 / 3.0 * M_PI * R * R * R;

    printf("\n=== 1) area e volume da esfera ===\n");
    printf("    exatos: A = %.10f   V = %.10f\n", (double) A_ex, (double) V_ex);

    real eA[5], eV[5];
    int  niveis = 5;
    for (int n = 0; n < niveis; n++) {
        ft3_superficie *s = ft3_cria_esfera(centro, R, n);
        real A = ft3_area(s), V = ft3_volume(s);
        eA[n] = fabs(A - A_ex) / A_ex;
        eV[n] = fabs(V - V_ex) / V_ex;
        int chi = ft3_euler(s);
        printf("    nsub=%d  nv=%5d nt=%5d  A=%.8f (err %.3e)  V=%.8f (err %.3e)  chi=%d\n",
               n, ft3_num_vertices(s), ft3_num_triangulos(s),
               (double) A, (double) eA[n], (double) V, (double) eV[n], chi);
        if (chi != 2) { ok(0, "topologia", "V - A + F != 2"); }
        if (V <= 0.0) { ok(0, "orientacao", "volume nao-positivo: normais invertidas"); }
        ft3_destroi(s);
    }
    ok(1, "topologia", "V - A + F = 2 em todos os niveis");
    ok(1, "orientacao", "volume positivo: normais para fora");

    // Ordem: a icosfera INSCRITA converge como h^2, e h cai pela metade a cada
    // subdivisao, entao o erro tem de cair por ~4.
    real pA = log(eA[niveis-2] / eA[niveis-1]) / log(2.0);
    real pV = log(eV[niveis-2] / eV[niveis-1]) / log(2.0);
    char buf[160];
    snprintf(buf, sizeof buf, "ordem da area p=%.2f (esperado 2)", (double) pA);
    ok(pA > 1.8 && pA < 2.2, "convergencia A", buf);
    snprintf(buf, sizeof buf, "ordem do volume p=%.2f (esperado 2)", (double) pV);
    ok(pV > 1.8 && pV < 2.2, "convergencia V", buf);
}

// =========================================================================
// 3: forca de tensao superficial
// =========================================================================
static void portao_forca(void)
{
    const Point centro = {0.0, 0.0, 0.0};
    const real  R = 0.5, sigma = 24.5;
    const real  alvo = 2.0 * sigma / R;          // DOIS sobre R

    printf("\n=== 3) forca: soma zero, e |F|/area -> 2 sigma / R ===\n");
    printf("    alvo analitico: 2*sigma/R = %.8f\n", (double) alvo);

    real err[5];
    int niveis = 5;
    for (int n = 0; n < niveis; n++) {
        ft3_superficie *s = ft3_cria_esfera(centro, R, n);
        int nv = ft3_num_vertices(s);
        Point *pos = malloc((size_t) nv * sizeof *pos);
        Point *f   = malloc((size_t) nv * sizeof *f);
        real  *w   = malloc((size_t) nv * sizeof *w);
        ft3_forcas_tensao(s, sigma, pos, f, w);

        // (a) soma das forcas: zero por construcao (integral de linha fechada)
        real soma[3] = {0,0,0}, escala = 0.0;
        for (int i = 0; i < nv; i++)
            for (int d = 0; d < 3; d++) {
                soma[d] += f[i][d];
                escala += fabs(f[i][d]);
            }
        real res = fabs(soma[0]) + fabs(soma[1]) + fabs(soma[2]);
        real rel = escala > 0 ? res / escala : 0.0;

        // (b) magnitude por area, e direcao (tem de apontar para DENTRO)
        real med = 0.0, wtot = 0.0;
        int  fora = 0;
        for (int i = 0; i < nv; i++) {
            real m = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2]);
            if (w[i] > 0.0) { med += m; wtot += w[i]; }
            // normal para fora na esfera centrada na origem e' a propria posicao
            real dot = f[i][0]*pos[i][0] + f[i][1]*pos[i][1] + f[i][2]*pos[i][2];
            if (dot >= 0.0) fora++;
        }
        real razao = wtot > 0 ? med / wtot : 0.0;
        err[n] = fabs(razao - alvo) / alvo;

        printf("    nsub=%d  nv=%5d  |sum F|/escala=%.3e  |F|/area=%.6f (err %.3e)  apontando p/ fora: %d\n",
               n, nv, (double) rel, (double) razao, (double) err[n], fora);

        if (rel > 1e-12) ok(0, "soma da forca", "residuo acima de 1e-12 da escala");
        if (fora != 0)   ok(0, "direcao da forca", "algum vertice com forca para fora");

        free(pos); free(f); free(w);
        ft3_destroi(s);
    }
    ok(1, "soma da forca", "residuo <= 1e-12 da escala em todos os niveis");
    ok(1, "direcao da forca", "todos os vertices com forca para dentro");

    real p = log(err[niveis-2] / err[niveis-1]) / log(2.0);
    char buf[160];
    snprintf(buf, sizeof buf, "err final %.3e, ordem p=%.2f", (double) err[niveis-1], (double) p);
    ok(err[niveis-1] < 1e-2 && p > 1.5, "convergencia da forca", buf);
}

// =========================================================================
// 4: cirurgia
// =========================================================================
static int inverteu_algum(const ft3_superficie *s)
{
    // Numa superficie estrelada em relacao ao centroide, todo triangulo bem
    // orientado tem normal com componente positiva na direcao radial.
    const Point *x = ft3_posicoes(s);
    const int (*tri)[3] = ft3_triangulos(s);
    int nt = ft3_num_triangulos(s), nv = ft3_num_vertices(s);
    real c[3] = {0,0,0};
    for (int i = 0; i < nv; i++) for (int d = 0; d < 3; d++) c[d] += x[i][d] / nv;
    int mau = 0;
    for (int t = 0; t < nt; t++) {
        real u[3], v[3], n[3], r[3];
        for (int d = 0; d < 3; d++) {
            u[d] = x[tri[t][1]][d] - x[tri[t][0]][d];
            v[d] = x[tri[t][2]][d] - x[tri[t][0]][d];
            r[d] = x[tri[t][0]][d] - c[d];
        }
        n[0] = u[1]*v[2]-u[2]*v[1]; n[1] = u[2]*v[0]-u[0]*v[2]; n[2] = u[0]*v[1]-u[1]*v[0];
        if (n[0]*r[0] + n[1]*r[1] + n[2]*r[2] <= 0.0) mau++;
    }
    return mau;
}

static void portao_cirurgia(void)
{
    printf("\n=== 4) cirurgia: elipsoide achatado, arestas desiguais ===\n");
    const Point centro = {0,0,0};
    ft3_superficie *s = ft3_cria_esfera(centro, 0.5, 3);

    // Achata em z e estica em x: cria aresta curta e longa na mesma malha.
    // Mexer nas posicoes por fora e' legitimo -- e' o que a adveccao fara'.
    {
        Point *x = (Point *) ft3_posicoes(s);
        for (int i = 0; i < ft3_num_vertices(s); i++) { x[i][0] *= 2.2; x[i][2] *= 0.35; }
    }
    real amin, amax, pior;
    ft3_qualidade(s, &amin, &amax, &pior);
    real ds = ft3_ds_alvo(s);
    printf("    antes:  nv=%5d nt=%5d  aresta [%.5f , %.5f]  pior aspecto %.2f  chi=%d\n",
           ft3_num_vertices(s), ft3_num_triangulos(s),
           (double) amin, (double) amax, (double) pior, ft3_euler(s));
    printf("    ds_alvo=%.5f  ==>  faixa aceita [%.5f , %.5f]\n",
           (double) ds, (double) (0.5*ds), (double) (2.0*ds));

    int inv0 = inverteu_algum(s);
    printf("    triangulos invertidos ANTES da cirurgia: %d\n", inv0);
    real A0 = ft3_area(s), V0 = ft3_volume(s);
    int n = ft3_cirurgia(s);
    real A1 = ft3_area(s), V1 = ft3_volume(s);
    ft3_qualidade(s, &amin, &amax, &pior);
    int chi = ft3_euler(s), inv = inverteu_algum(s);
    printf("    depois: nv=%5d nt=%5d  aresta [%.5f , %.5f]  pior aspecto %.2f  chi=%d\n",
           ft3_num_vertices(s), ft3_num_triangulos(s),
           (double) amin, (double) amax, (double) pior, chi);
    printf("    %d operacoes; area %.6f -> %.6f (%.2f%%)  volume %.6f -> %.6f (%.2f%%)\n",
           n, (double) A0, (double) A1, (double) ((A1-A0)/A0*100),
           (double) V0, (double) V1, (double) ((V1-V0)/V0*100));

    char buf[160];
    ok(n > 0, "cirurgia agiu", "houve operacao de divisao/colapso/giro");
    snprintf(buf, sizeof buf, "chi=%d (tem de ser 2)", chi);
    ok(chi == 2, "topologia preservada", buf);
    snprintf(buf, sizeof buf, "%d invertidos antes, %d depois", inv0, inv);
    ok(inv <= inv0, "sem inversao criada", buf);
    snprintf(buf, sizeof buf, "maior aresta %.5f <= 2*ds = %.5f", (double) amax, (double) (2.0*ds));
    ok(amax <= 2.0 * ds * 1.02, "aresta longa tratada", buf);
    // O volume nao pode saltar: cirurgia e' remalhamento, nao fisica.
    snprintf(buf, sizeof buf, "variacao de volume %.3f%%", (double) ((V1-V0)/V0*100));
    ok(fabs((V1-V0)/V0) < 0.02, "volume preservado", buf);

    ft3_destroi(s);
}

// =========================================================================
// 5: persistencia
// =========================================================================
static void portao_persistencia(void)
{
    printf("\n=== 5) persistencia com conectividade ===\n");
    const char *arq = "/tmp/ft3-teste.superficie";
    const Point centro = {0.1, 0.2, 0.3};
    ft3_superficie *a = ft3_cria_esfera(centro, 0.37, 2);
    ft3_cirurgia(a);                     // grava algo que ja' passou por cirurgia

    ok(ft3_grava(a, arq) == 0, "gravacao", "arquivo escrito");
    ft3_superficie *b = ft3_le(arq);
    if (!b) { ok(0, "leitura", "devolveu NULL"); ft3_destroi(a); return; }

    int igual = (ft3_num_vertices(a) == ft3_num_vertices(b))
             && (ft3_num_triangulos(a) == ft3_num_triangulos(b));
    const Point *xa = ft3_posicoes(a), *xb = ft3_posicoes(b);
    for (int i = 0; igual && i < ft3_num_vertices(a); i++)
        for (int d = 0; d < 3; d++)
            if (xa[i][d] != xb[i][d]) igual = 0;
    const int (*ta)[3] = ft3_triangulos(a), (*tb)[3] = ft3_triangulos(b);
    for (int t = 0; igual && t < ft3_num_triangulos(a); t++)
        for (int k = 0; k < 3; k++)
            if (ta[t][k] != tb[t][k]) igual = 0;

    char buf[160];
    snprintf(buf, sizeof buf, "%d vertices e %d triangulos identicos bit a bit",
             ft3_num_vertices(a), ft3_num_triangulos(a));
    ok(igual, "ida-e-volta", buf);
    ft3_destroi(b);

    // Recusa: ausente, cabecalho estranho, truncado, indice fora de faixa.
    ok(ft3_le("/tmp/ft3-nao-existe-mesmo") == NULL, "ausencia", "arquivo inexistente -> NULL");

    FILE *fp = fopen(arq, "w"); fprintf(fp, "isto nao e' ft3\n"); fclose(fp);
    ok(ft3_le(arq) == NULL, "cabecalho", "formato estranho -> NULL");

    fp = fopen(arq, "w"); fprintf(fp, "ft3 1 100 200 0.1\n0 0 0\n"); fclose(fp);
    ok(ft3_le(arq) == NULL, "truncado", "menos vertices do que promete -> NULL");

    fp = fopen(arq, "w");
    fprintf(fp, "ft3 1 4 4 0.1\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n0 1 2\n0 1 3\n0 2 3\n1 2 99\n");
    fclose(fp);
    ok(ft3_le(arq) == NULL, "indice invalido", "triangulo fora de faixa -> NULL");

    remove(arq);
    ft3_destroi(a);
}

int main(void)
{
    printf("=== teste-esfera3d: geometria e forca da frente 3D ===\n");
    portao_geometria();
    portao_forca();
    portao_cirurgia();
    portao_persistencia();
    printf("\n%s\n", falhas ? "TESTE DA ESFERA 3D FALHOU" : "TESTE DA ESFERA 3D PASSOU");
    return falhas ? 1 : 0;
}
