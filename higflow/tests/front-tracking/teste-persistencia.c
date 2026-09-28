// Arreio da PERSISTENCIA da frente: ft_grava / ft_le.
//
// O que esta' em jogo.  O `h.save` do solver guarda campos na malha; a frente
// lagrangeana nao e' campo, e nenhum campo a determina.  Retomar uma corrida sem
// gravar a frente poe os campos no instante salvo e os marcadores em t=0 -- e
// roda, em silencio.  Por isso a gravacao existe, e por isso ela precisa de
// oraculo: uma retomada que "quase" reproduz a corrida original nao serve, pois
// o erro nao se manifesta como falha, mas como resultado diferente.
//
// Quatro portoes, do mais fraco ao mais forte:
//   (1) ida-e-volta EXATA numa frente regular (circulo)
//   (2) ida-e-volta EXATA numa frente deformada pela adveccao + cirurgia,
//       onde n mudou e as posicoes nao tem simetria nenhuma
//   (3) ausencia e corrupcao do arquivo devolvem NULL -- nao lixo, nao crash
//   (4) CONTINUACAO: advectar 100 passos direto e advectar 50, gravar, ler e
//       advectar os outros 50 dao o MESMO resultado, bit a bit.  E' a unica
//       propriedade que de fato importa -- que retomar seja a mesma corrida.

#include "hig-flow-front-tracking.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *ARQ = "/tmp/ft-teste-persistencia.frente";

typedef struct { real T; } vortex_ctx;

static void campo_rider_kothe(const Point x, real t, void *vctx, real u[DIM])
{
    vortex_ctx *c = (vortex_ctx *) vctx;
    real g = cos(M_PI * t / c->T);
    real sx = sin(M_PI * x[0]), sy = sin(M_PI * x[1]);
    u[0] =  sx * sx * sin(2.0 * M_PI * x[1]) * g;
    u[1] = -sin(2.0 * M_PI * x[0]) * sy * sy * g;
}

// Compara duas frentes BIT A BIT.  Devolve 0 se identicas, 1 caso contrario.
static int difere(const ft_frente *a, const ft_frente *b, const char *rotulo)
{
    if (ft_num(a) != ft_num(b)) {
        printf("  FALHOU %s: n %d != %d\n", rotulo, ft_num(a), ft_num(b));
        return 1;
    }
    if (ft_ds_alvo(a) != ft_ds_alvo(b)) {
        printf("  FALHOU %s: ds_alvo %.17g != %.17g\n", rotulo,
               (double) ft_ds_alvo(a), (double) ft_ds_alvo(b));
        return 1;
    }
    int n = ft_num(a);
    Point *pa = malloc((size_t) n * sizeof(Point));
    Point *pb = malloc((size_t) n * sizeof(Point));
    ft_posicoes(a, pa); ft_posicoes(b, pb);
    int ruim = 0;
    real pior = 0.0;
    for (int i = 0; i < n && !ruim; i++)
        for (int d = 0; d < DIM; d++) {
            real dif = fabs((double) pa[i][d] - (double) pb[i][d]);
            if (dif > pior) pior = dif;
            if (pa[i][d] != pb[i][d]) {
                printf("  FALHOU %s: marcador %d coord %d: "
                       "%.17g != %.17g (dif %.3e)\n", rotulo, i, d,
                       (double) pa[i][d], (double) pb[i][d], (double) dif);
                ruim = 1;
            }
        }
    free(pa); free(pb);
    if (!ruim) printf("  ok %s: %d marcadores identicos bit a bit\n", rotulo, n);
    return ruim;
}

// (1) e (2): ida-e-volta
static int porta_ida_e_volta(ft_frente *f, const char *rotulo)
{
    if (ft_grava(f, ARQ) != ft_num(f)) {
        printf("  FALHOU %s: ft_grava nao devolveu %d\n", rotulo, ft_num(f));
        return 1;
    }
    ft_frente *g = ft_le(ARQ);
    if (g == NULL) { printf("  FALHOU %s: ft_le devolveu NULL\n", rotulo); return 1; }
    int ruim = difere(f, g, rotulo);
    // A area e' funcao das posicoes: se elas sao identicas, ela tambem e'.
    // Conferir e' redundante de proposito -- pega um erro de leitura que
    // trocasse a ORDEM dos marcadores sem mudar o conjunto.
    if (!ruim && ft_area(f) != ft_area(g)) {
        printf("  FALHOU %s: area %.17g != %.17g (ordem trocada?)\n", rotulo,
               (double) ft_area(f), (double) ft_area(g));
        ruim = 1;
    }
    ft_destroi(g);
    return ruim;
}

// (3) ausencia e corrupcao
static int porta_arquivo_ruim(void)
{
    int ruim = 0;
    remove("/tmp/ft-teste-nao-existe.frente");
    if (ft_le("/tmp/ft-teste-nao-existe.frente") != NULL) {
        printf("  FALHOU ausencia: ft_le de arquivo inexistente nao deu NULL\n");
        ruim = 1;
    } else {
        printf("  ok ausencia: arquivo inexistente -> NULL\n");
    }

    // Cabecalho errado.
    FILE *fp = fopen(ARQ, "w");
    fprintf(fp, "posicoes da gota\n4 0.1\n0 0\n1 0\n1 1\n0 1\n");
    fclose(fp);
    if (ft_le(ARQ) != NULL) {
        printf("  FALHOU cabecalho: aceitou arquivo sem o cabecalho\n");
        ruim = 1;
    } else {
        printf("  ok cabecalho: formato estranho -> NULL\n");
    }

    // TRUNCADO: promete 64 marcadores e entrega 3.  E' o arquivo que uma queda
    // de energia no meio da escrita deixa, e o modo de falha mais perigoso --
    // aceitar isso daria uma frente com marcadores em memoria nao inicializada.
    fp = fopen(ARQ, "w");
    fprintf(fp, "# frente front-tracking v1\n64 0.1\n0 0\n1 0\n1 1\n");
    fclose(fp);
    if (ft_le(ARQ) != NULL) {
        printf("  FALHOU truncado: aceitou arquivo cortado no meio\n");
        ruim = 1;
    } else {
        printf("  ok truncado: n prometido maior que o entregue -> NULL\n");
    }

    // n absurdo.
    fp = fopen(ARQ, "w");
    fprintf(fp, "# frente front-tracking v1\n2 0.1\n0 0\n1 0\n");
    fclose(fp);
    if (ft_le(ARQ) != NULL) {
        printf("  FALHOU n<3: aceitou frente de 2 marcadores\n");
        ruim = 1;
    } else {
        printf("  ok n<3: poligono impossivel -> NULL\n");
    }
    return ruim;
}

// (4) o portao forte: a retomada e' a MESMA corrida.
static int porta_continuacao(void)
{
    const real T = 2.0, dt = 0.005;
    const int  nsteps = 100, corte = 50;
    vortex_ctx ctx = { T };
    Point centro = {0.5, 0.75};

    // Corrida A: direto, 100 passos.
    ft_frente *A = ft_cria_circulo(centro, 0.15, 64);
    for (int s = 0; s < nsteps; s++) {
        ft_advecta(A, campo_rider_kothe, &ctx, (real) s * dt, dt);
        ft_cirurgia(A);
    }

    // Corrida B: 50 passos, grava, le, e os outros 50 -- exatamente o que a
    // retomada faz.
    ft_frente *B = ft_cria_circulo(centro, 0.15, 64);
    for (int s = 0; s < corte; s++) {
        ft_advecta(B, campo_rider_kothe, &ctx, (real) s * dt, dt);
        ft_cirurgia(B);
    }
    ft_grava(B, ARQ);
    ft_destroi(B);
    B = ft_le(ARQ);
    if (B == NULL) { printf("  FALHOU continuacao: ft_le deu NULL\n"); ft_destroi(A); return 1; }
    for (int s = corte; s < nsteps; s++) {
        ft_advecta(B, campo_rider_kothe, &ctx, (real) s * dt, dt);
        ft_cirurgia(B);
    }

    printf("  continuacao: A tem %d marcadores (area %.10f), "
           "B %d (area %.10f)\n", ft_num(A), (double) ft_area(A),
           ft_num(B), (double) ft_area(B));
    int ruim = difere(A, B, "continuacao 100 = 50 + retomada + 50");
    ft_destroi(A); ft_destroi(B);
    return ruim;
}

int main(void)
{
    int ruim = 0;

    printf("=== (1) ida-e-volta, frente regular ===\n");
    Point centro = {0.5, 0.5};
    ft_frente *c = ft_cria_circulo(centro, 1.0 / 6.0, 137);
    ruim |= porta_ida_e_volta(c, "circulo de 137 marcadores");
    ft_destroi(c);

    printf("=== (2) ida-e-volta, frente deformada ===\n");
    // Deforma de verdade: o vortice estica a frente, a cirurgia insere e remove
    // marcadores, e as posicoes ficam sem simetria nenhuma.  Uma gravacao que
    // arredondasse na 15a casa passaria em (1) e falharia aqui.
    vortex_ctx ctx = { 2.0 };
    Point c0 = {0.5, 0.75};
    ft_frente *d = ft_cria_circulo(c0, 0.15, 64);
    for (int s = 0; s < 200; s++) {
        ft_advecta(d, campo_rider_kothe, &ctx, (real) s * 0.005, 0.005);
        ft_cirurgia(d);
    }
    printf("  frente deformada: %d marcadores, area %.10f\n",
           ft_num(d), (double) ft_area(d));
    ruim |= porta_ida_e_volta(d, "frente deformada");
    ft_destroi(d);

    printf("=== (3) arquivo ausente e corrompido ===\n");
    ruim |= porta_arquivo_ruim();

    printf("=== (4) a retomada e' a mesma corrida ===\n");
    ruim |= porta_continuacao();

    remove(ARQ);
    printf("\n%s\n", ruim ? "TESTE DE PERSISTENCIA FALHOU"
                          : "TESTE DE PERSISTENCIA PASSOU");
    return ruim;
}
