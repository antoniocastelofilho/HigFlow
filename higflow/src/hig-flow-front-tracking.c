// Ver o .h para o contrato.  FASE B1: adveccao cinematica + cirurgia, campo
// prescrito, frente ordenada e replicada.  Sem solver, sem MPI -- a frente vive
// inteira aqui.  O paralelismo (posicoes distribuidas, topologia replicada) e' a
// decisao #1 do projeto e entra quando o B4 acoplar a malha.

#include "hig-flow-front-tracking.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if DIM != 2
// O B1 e' 2D (a frente e' uma polilinha).  O 3D (superficie triangulada,
// Laplace-Beltrami) e' a fase B5 e usara' outra representacao de topologia.
// Compilar este arquivo em 3D e' engano de build, nao caso a tratar.
#error "hig-flow-front-tracking.c: B1 e' 2D; ver a fase B5 para 3D"
#endif

struct ft_frente {
    Point *x;        // posicoes, EM ORDEM; x[i] liga em x[i+1], x[n-1] em x[0]
    int    n;        // numero de marcadores
    int    cap;      // capacidade alocada de x[]
    real   ds_alvo;  // espacamento alvo (criterio da cirurgia)
};

// Garante capacidade para pelo menos `nec` marcadores, preservando o conteudo.
static void _garante_cap(ft_frente *f, int nec)
{
    if (f->cap >= nec) return;
    int nova = f->cap ? f->cap : 8;
    while (nova < nec) nova *= 2;
    f->x = (Point *) realloc(f->x, (size_t) nova * sizeof(Point));
    if (f->x == NULL) {
        fprintf(stderr, "ft: realloc de %d marcadores falhou\n", nova);
        exit(1);
    }
    f->cap = nova;
}

// --- criacao / destruicao -------------------------------------------------

ft_frente *ft_cria_curva(const Point *vertices, int nvert, real ds_alvo)
{
    if (nvert < 3 || ds_alvo <= 0.0) {
        fprintf(stderr, "ft_cria_curva: nvert>=3 e ds_alvo>0 exigidos "
                        "(recebi nvert=%d ds_alvo=%g)\n", nvert, ds_alvo);
        exit(1);
    }
    ft_frente *f = (ft_frente *) calloc(1, sizeof *f);
    f->ds_alvo = ds_alvo;

    // Reamostra cada segmento do poligono de entrada a um espacamento ~ds_alvo,
    // pondo marcadores nas EXTREMIDADES dos subsegmentos (o marcador i e' o
    // inicio do segmento i->i+1).  O primeiro vertice nao se repete no fim: o
    // fechamento e' implicito.
    for (int v = 0; v < nvert; v++) {
        const real *a = vertices[v];
        const real *b = vertices[(v + 1) % nvert];
        real dx = b[0] - a[0], dy = b[1] - a[1];
        real L = sqrt(dx * dx + dy * dy);
        int nsub = (int) ceil(L / ds_alvo);
        if (nsub < 1) nsub = 1;
        // Poe nsub marcadores no segmento [a,b): a, a+L/nsub, ..., sem incluir b
        // (b entra como primeiro do proximo segmento).
        for (int k = 0; k < nsub; k++) {
            _garante_cap(f, f->n + 1);
            real s = (real) k / (real) nsub;
            f->x[f->n][0] = a[0] + s * dx;
            f->x[f->n][1] = a[1] + s * dy;
            f->n++;
        }
    }
    return f;
}

ft_frente *ft_cria_circulo(const Point centro, real raio, int nmarc)
{
    if (nmarc < 3 || raio <= 0.0) {
        fprintf(stderr, "ft_cria_circulo: nmarc>=3 e raio>0 exigidos\n");
        exit(1);
    }
    ft_frente *f = (ft_frente *) calloc(1, sizeof *f);
    _garante_cap(f, nmarc);
    for (int k = 0; k < nmarc; k++) {
        real th = 2.0 * M_PI * (real) k / (real) nmarc;
        f->x[k][0] = centro[0] + raio * cos(th);
        f->x[k][1] = centro[1] + raio * sin(th);
    }
    f->n = nmarc;
    // Espacamento alvo = corda entre marcadores vizinhos do circulo inicial.
    f->ds_alvo = 2.0 * M_PI * raio / (real) nmarc;
    return f;
}

void ft_destroi(ft_frente *f)
{
    if (f == NULL) return;
    free(f->x);
    free(f);
}

// --- consulta -------------------------------------------------------------

int  ft_num(const ft_frente *f)     { return f->n; }
real ft_ds_alvo(const ft_frente *f) { return f->ds_alvo; }

real ft_area(const ft_frente *f)
{
    // Formula do laco (shoelace): A = 1/2 |sum_i (x_i * y_{i+1} - x_{i+1} * y_i)|,
    // com i+1 modulo n para fechar a curva.
    real soma = 0.0;
    for (int i = 0; i < f->n; i++) {
        const real *a = f->x[i];
        const real *b = f->x[(i + 1) % f->n];
        soma += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * fabs(soma);
}

real ft_perimetro(const ft_frente *f)
{
    real p = 0.0;
    for (int i = 0; i < f->n; i++) {
        const real *a = f->x[i];
        const real *b = f->x[(i + 1) % f->n];
        real dx = b[0] - a[0], dy = b[1] - a[1];
        p += sqrt(dx * dx + dy * dy);
    }
    return p;
}

void ft_posicoes(const ft_frente *f, Point *dst)
{
    memcpy(dst, f->x, (size_t) f->n * sizeof(Point));
}

// --- geometria: curvatura e normal ----------------------------------------

// Centro do circulo que passa por tres pontos (circuncentro).  Devolve 0 e nao
// escreve `centro` se os tres forem quase colineares (D ~ 0).
static int _circuncentro(const real *A, const real *B, const real *C,
                         real centro[2])
{
    real D = 2.0 * (A[0] * (B[1] - C[1])
                  + B[0] * (C[1] - A[1])
                  + C[0] * (A[1] - B[1]));
    if (fabs(D) < 1e-300) return 0;
    real A2 = A[0] * A[0] + A[1] * A[1];
    real B2 = B[0] * B[0] + B[1] * B[1];
    real C2 = C[0] * C[0] + C[1] * C[1];
    centro[0] = (A2 * (B[1] - C[1]) + B2 * (C[1] - A[1]) + C2 * (A[1] - B[1])) / D;
    centro[1] = (A2 * (C[0] - B[0]) + B2 * (A[0] - C[0]) + C2 * (B[0] - A[0])) / D;
    return 1;
}

void ft_vetor_curvatura(const ft_frente *f, Point *kv)
{
    for (int i = 0; i < f->n; i++) {
        const real *p0 = f->x[(i - 1 + f->n) % f->n];
        const real *p1 = f->x[i];
        const real *p2 = f->x[(i + 1) % f->n];
        real c[2];
        if (!_circuncentro(p0, p1, p2, c)) {
            kv[i][0] = 0.0; kv[i][1] = 0.0;   // trecho reto: curvatura nula
            continue;
        }
        // Vetor de p1 ao centro; |dif| = R.  kappa*n = dif / R^2, modulo 1/R,
        // apontando para o centro (lado concavo).
        real dx = c[0] - p1[0], dy = c[1] - p1[1];
        real R2 = dx * dx + dy * dy;
        kv[i][0] = dx / R2;
        kv[i][1] = dy / R2;
    }
}

void ft_curvatura(const ft_frente *f, real *kappa)
{
    Point *kv = (Point *) malloc((size_t) f->n * sizeof(Point));
    ft_vetor_curvatura(f, kv);
    for (int i = 0; i < f->n; i++)
        kappa[i] = sqrt(kv[i][0] * kv[i][0] + kv[i][1] * kv[i][1]);
    free(kv);
}

// --- adveccao -------------------------------------------------------------

void ft_advecta(ft_frente *f, ft_campo_u u, void *ctx, real t, real dt)
{
    // RK2 do ponto medio, por marcador:
    //   k1 = u(x, t);   xm = x + (dt/2) k1;   k2 = u(xm, t + dt/2)
    //   x += dt k2
    // Segunda ordem no tempo -- Euler explicito nao fecha a area de volta no
    // teste reversivel.  A frente e' replicada, entao cada marcador anda sozinho.
    for (int i = 0; i < f->n; i++) {
        real k1[DIM], k2[DIM];
        Point xm;
        u(f->x[i], t, ctx, k1);
        for (int d = 0; d < DIM; d++) xm[d] = f->x[i][d] + 0.5 * dt * k1[d];
        u(xm, t + 0.5 * dt, ctx, k2);
        for (int d = 0; d < DIM; d++) f->x[i][d] += dt * k2[d];
    }
}

// --- cirurgia -------------------------------------------------------------

// Area (com sinal) do triangulo (p, a, b) -- o que a remocao de `a` cortaria da
// polilinha ao ligar p direto em b.  E' a medida de curvatura que decide se
// remover `a` e' seguro: num trecho quase reto o triangulo e' minusculo.
static real _tri_area(const real *p, const real *a, const real *b)
{
    return 0.5 * fabs((a[0] - p[0]) * (b[1] - p[1])
                    - (b[0] - p[0]) * (a[1] - p[1]));
}

int ft_cirurgia(ft_frente *f)
{
    const real l_max = 2.0 * f->ds_alvo;   // segmento acima disto: inserir
    const real l_min = 0.5 * f->ds_alvo;   // segmento abaixo disto: candidato a remover
    // A remocao so' e' SEGURA onde a curva e' quase reta -- senao corta canto e
    // perde area (medido: a remocao cega dominava o erro do B1).  Teto do
    // triangulo cortado, calibrado em fracao de ds_alvo^2.
    const real area_tol = 1e-3 * f->ds_alvo * f->ds_alvo;

    Point *novo = (Point *) malloc((size_t) (2 * f->n + 8) * sizeof(Point));
    int m = 0;

    // PASSE 1 -- REMOCAO guardada por curvatura.  Percorre em ordem mantendo o
    // ultimo marcador COPIADO como `p`; remove x[i] so' se estiver perto de p
    // (segmento curto) E quase sobre a corda p->b (triangulo minusculo).
    for (int i = 0; i < f->n; i++) {
        const real *a = f->x[i];
        const real *b = f->x[(i + 1) % f->n];
        if (m > 0 && f->n - (i - m) > 3) {
            const real *p = novo[m - 1];
            real dx = a[0] - p[0], dy = a[1] - p[1];
            real Lpa = sqrt(dx * dx + dy * dy);
            if (Lpa < l_min && _tri_area(p, a, b) < area_tol)
                continue;  // seguro remover `a`: funde p->b
        }
        novo[m][0] = a[0];
        novo[m][1] = a[1];
        m++;
    }

    // PASSE 2 -- INSERCAO do ponto medio (sobre a corda: neutra em area) onde o
    // segmento ficou longo demais.  Reconstroi in-place em `f->x`.
    _garante_cap(f, 2 * m + 8);
    int q = 0;
    for (int i = 0; i < m; i++) {
        const real *a = novo[i];
        const real *b = novo[(i + 1) % m];
        f->x[q][0] = a[0];
        f->x[q][1] = a[1];
        q++;
        real dx = b[0] - a[0], dy = b[1] - a[1];
        if (sqrt(dx * dx + dy * dy) > l_max) {
            f->x[q][0] = a[0] + 0.5 * dx;
            f->x[q][1] = a[1] + 0.5 * dy;
            q++;
        }
    }
    f->n = q;
    free(novo);
    return f->n;
}

// --- saida ----------------------------------------------------------------

void ft_escreve_vtk(const ft_frente *f, const char *prefixo, int quadro)
{
    char nome[512];
    snprintf(nome, sizeof nome, "%s_ft_%d.vtk", prefixo, quadro);
    FILE *fp = fopen(nome, "w");
    if (fp == NULL) {
        fprintf(stderr, "ft_escreve_vtk: nao abriu %s\n", nome);
        return;
    }
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "front-tracking quadro %d\n", quadro);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %d double\n", f->n);
    for (int i = 0; i < f->n; i++)
        fprintf(fp, "%.10g %.10g 0\n", f->x[i][0], f->x[i][1]);
    // Uma unica polilinha fechada: n+1 indices, repetindo o 0 no fim.
    fprintf(fp, "LINES 1 %d\n", f->n + 2);
    fprintf(fp, "%d", f->n + 1);
    for (int i = 0; i < f->n; i++) fprintf(fp, " %d", i);
    fprintf(fp, " 0\n");
    fclose(fp);
}
