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

// --- funcao indicadora (para a forca balanceada) ---------------------------

// Recorte de Sutherland-Hodgman por UMA aresta do retangulo.  `lado` escolhe a
// aresta (0=x>=lo, 1=x<=hi, 2=y>=lo, 3=y<=hi) e `val` o seu valor.  Escreve o
// resultado em `saida` e devolve quantos vertices ele tem.
//
// A caixa e' CONVEXA, que e' a condicao para o recorte sucessivo pelas quatro
// arestas dar a intersecao correta.
static int _recorta_aresta(const Point *ent, int n, int lado, real val,
                           Point *saida)
{
    int m = 0;
    if (n == 0) return 0;
    for (int i = 0; i < n; i++) {
        const real *A = ent[i];
        const real *B = ent[(i + 1) % n];
        // "dentro" para cada uma das quatro arestas
        real da, db;
        switch (lado) {
            case 0: da = A[0] - val; db = B[0] - val; break;  // x >= lo
            case 1: da = val - A[0]; db = val - B[0]; break;  // x <= hi
            case 2: da = A[1] - val; db = B[1] - val; break;  // y >= lo
            default: da = val - A[1]; db = val - B[1]; break; // y <= hi
        }
        const int A_dentro = (da >= 0.0), B_dentro = (db >= 0.0);

        if (A_dentro) { saida[m][0] = A[0]; saida[m][1] = A[1]; m++; }
        if (A_dentro != B_dentro) {
            // Interseccao com a reta da aresta: parametro pela razao das
            // distancias com sinal, que e' exata para segmento contra reta.
            real t = da / (da - db);
            saida[m][0] = A[0] + t * (B[0] - A[0]);
            saida[m][1] = A[1] + t * (B[1] - A[1]);
            m++;
        }
    }
    return m;
}

real ft_area_na_caixa(const ft_frente *f, const real lo[DIM], const real hi[DIM])
{
    // Caminho rapido: caixa fora da envoltoria da frente nao pode conter area
    // nenhuma dela -- tudo que a frente fecha esta' dentro da propria
    // envoltoria.
    real bx0 = f->x[0][0], bx1 = f->x[0][0];
    real by0 = f->x[0][1], by1 = f->x[0][1];
    for (int i = 1; i < f->n; i++) {
        if (f->x[i][0] < bx0) bx0 = f->x[i][0];
        if (f->x[i][0] > bx1) bx1 = f->x[i][0];
        if (f->x[i][1] < by0) by0 = f->x[i][1];
        if (f->x[i][1] > by1) by1 = f->x[i][1];
    }
    if (hi[0] <= bx0 || lo[0] >= bx1 || hi[1] <= by0 || lo[1] >= by1) return 0.0;

    // Recorte sucessivo pelas quatro arestas.  Cada recorte pode acrescentar no
    // maximo um vertice por aresta do poligono, dai' a folga na alocacao.
    const int cap = 2 * f->n + 8;
    Point *a = (Point *) malloc((size_t) cap * sizeof(Point));
    Point *b = (Point *) malloc((size_t) cap * sizeof(Point));
    memcpy(a, f->x, (size_t) f->n * sizeof(Point));
    int n = f->n;

    n = _recorta_aresta(a, n, 0, lo[0], b);  memcpy(a, b, (size_t) n * sizeof(Point));
    n = _recorta_aresta(a, n, 1, hi[0], b);  memcpy(a, b, (size_t) n * sizeof(Point));
    n = _recorta_aresta(a, n, 2, lo[1], b);  memcpy(a, b, (size_t) n * sizeof(Point));
    n = _recorta_aresta(a, n, 3, hi[1], b);  memcpy(a, b, (size_t) n * sizeof(Point));

    real soma = 0.0;
    for (int i = 0; i < n; i++) {
        const real *p = b[i], *q = b[(i + 1) % n];
        soma += p[0] * q[1] - q[0] * p[1];
    }
    free(a); free(b);
    return 0.5 * fabs(soma);
}

// --- forca de tensao superficial ------------------------------------------

real ft_delta_roma(real r)
{
    r = fabs(r);
    if (r <= 0.5) return (1.0 / 3.0) * (1.0 + sqrt(1.0 - 3.0 * r * r));
    if (r <= 1.5) {
        real s = 1.0 - r;
        return (1.0 / 6.0) * (5.0 - 3.0 * r - sqrt(1.0 - 3.0 * s * s));
    }
    return 0.0;
}

// Comprimento de arco que o marcador `i` representa: metade de cada segmento
// vizinho (o do lado anterior e o do proximo).  A soma sobre i da' o perimetro.
static real _ds_marcador(const ft_frente *f, int i)
{
    const real *p = f->x[(i - 1 + f->n) % f->n];
    const real *a = f->x[i];
    const real *q = f->x[(i + 1) % f->n];
    real d0 = sqrt((a[0]-p[0])*(a[0]-p[0]) + (a[1]-p[1])*(a[1]-p[1]));
    real d1 = sqrt((q[0]-a[0])*(q[0]-a[0]) + (q[1]-a[1])*(q[1]-a[1]));
    return 0.5 * (d0 + d1);
}

void ft_integral_curvatura(const ft_frente *f, real integral[DIM])
{
    Point *kv = (Point *) malloc((size_t) f->n * sizeof(Point));
    ft_vetor_curvatura(f, kv);
    integral[0] = 0.0; integral[1] = 0.0;
    for (int i = 0; i < f->n; i++) {
        real ds = _ds_marcador(f, i);
        integral[0] += kv[i][0] * ds;
        integral[1] += kv[i][1] * ds;
    }
    free(kv);
}

void ft_forcas_tensao(const ft_frente *f, real sigma,
                      Point *pos, Point *forca, real *ds)
{
    Point *kv = (Point *) malloc((size_t) f->n * sizeof(Point));
    ft_vetor_curvatura(f, kv);
    for (int i = 0; i < f->n; i++) {
        pos[i][0] = f->x[i][0];
        pos[i][1] = f->x[i][1];
        forca[i][0] = sigma * kv[i][0];
        forca[i][1] = sigma * kv[i][1];
        ds[i] = _ds_marcador(f, i);
    }
    free(kv);
}

void ft_espalha_tensao(const ft_frente *f, real sigma,
                       real ox, real oy, real h, int nx, int ny,
                       real *fx, real *fy)
{
    Point *kv = (Point *) malloc((size_t) f->n * sizeof(Point));
    ft_vetor_curvatura(f, kv);

    for (int k = 0; k < f->n; k++) {
        real ds = _ds_marcador(f, k);
        real Fx = sigma * kv[k][0] * ds;
        real Fy = sigma * kv[k][1] * ds;
        const real *X = f->x[k];

        // Celula que contem o marcador, e janela de +-2 celulas (suporte 1,5h).
        int ic = (int) floor((X[0] - ox) / h - 0.5);
        int jc = (int) floor((X[1] - oy) / h - 0.5);
        for (int j = jc - 2; j <= jc + 2; j++) {
            if (j < 0 || j >= ny) continue;
            real yc = oy + (j + 0.5) * h;
            real wy = ft_delta_roma((yc - X[1]) / h);
            if (wy == 0.0) continue;
            for (int i = ic - 2; i <= ic + 2; i++) {
                if (i < 0 || i >= nx) continue;
                real xc = ox + (i + 0.5) * h;
                real wx = ft_delta_roma((xc - X[0]) / h);
                if (wx == 0.0) continue;
                // delta_h = (1/h^2) phi(x) phi(y); a forca DENSIDADE recebe
                // F_k * delta_h.  (Conservacao: SUM delta_h * h^2 = 1.)
                real d = wx * wy / (h * h);
                fx[j * nx + i] += Fx * d;
                fy[j * nx + i] += Fy * d;
            }
        }
    }
    free(kv);
}

void ft_semieixos(const ft_frente *f, real *a, real *b)
{
    // Formulas de poligono para area, centroide e momentos de segunda ordem.
    real A2 = 0.0, cx = 0.0, cy = 0.0, Mxx = 0.0, Myy = 0.0;
    for (int i = 0; i < f->n; i++) {
        const real *p = f->x[i], *q = f->x[(i + 1) % f->n];
        real cr = p[0] * q[1] - q[0] * p[1];          // produto cruzado
        A2  += cr;
        cx  += (p[0] + q[0]) * cr;
        cy  += (p[1] + q[1]) * cr;
        Mxx += (p[0]*p[0] + p[0]*q[0] + q[0]*q[0]) * cr;
        Myy += (p[1]*p[1] + p[1]*q[1] + q[1]*q[1]) * cr;
    }
    real A = 0.5 * A2;
    if (fabs(A) < 1e-300) { *a = *b = 0.0; return; }
    cx /= (6.0 * A);  cy /= (6.0 * A);
    Mxx = Mxx / 12.0 - A * cx * cx;                   // transporte ao centroide
    Myy = Myy / 12.0 - A * cy * cy;
    // Elipse de semieixos a,b tem Mxx = A a^2/4 e Myy = A b^2/4.
    *a = 2.0 * sqrt(fabs(Mxx / A));
    *b = 2.0 * sqrt(fabs(Myy / A));
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

// --- persistencia ---------------------------------------------------------
//
// POR QUE ISTO EXISTE.  O `h.save` do solver guarda CAMPOS na malha; a frente
// nao e' campo -- e' estado lagrangeano que nenhum campo determina.  Sem
// gravacao, retomar uma corrida punha os campos no instante salvo e os
// marcadores em t=0, o que e' pior do que nao retomar: roda sem reclamar.
//
// MEDIDO: uma corrida de Hysing de 12000 passos leva ~6 h nesta maquina, e ja'
// morreu uma vez no meio.  A retomada nao e' conveniencia.
//
// Formato de TEXTO, com %.17g -- que faz ida-e-volta exata em double.  Texto
// porque o arquivo tambem serve de despejo para as figuras, e porque um estado
// que nao se pode ler com `cat` e' um estado que nao se audita.

int ft_grava(const ft_frente *f, const char *arquivo)
{
    FILE *fp = fopen(arquivo, "w");
    if (fp == NULL) {
        fprintf(stderr, "ft_grava: nao abriu %s para escrita\n", arquivo);
        return -1;
    }
    fprintf(fp, "# frente front-tracking v1\n");
    fprintf(fp, "%d %.17g\n", f->n, (double) f->ds_alvo);
    for (int i = 0; i < f->n; i++)
        fprintf(fp, "%.17g %.17g\n", (double) f->x[i][0], (double) f->x[i][1]);
    int erro = ferror(fp);
    if (fclose(fp) != 0 || erro) {
        fprintf(stderr, "ft_grava: escrita de %s falhou\n", arquivo);
        return -1;
    }
    return f->n;
}

ft_frente *ft_le(const char *arquivo)
{
    FILE *fp = fopen(arquivo, "r");
    if (fp == NULL) return NULL;        // ausencia NAO e' erro: quem chama decide

    char linha[256];
    if (fgets(linha, sizeof linha, fp) == NULL) { fclose(fp); return NULL; }
    if (strncmp(linha, "# frente front-tracking v1", 26) != 0) {
        fprintf(stderr, "ft_le: %s nao tem o cabecalho esperado\n", arquivo);
        fclose(fp);
        return NULL;
    }
    int n = 0; double ds = 0.0;
    if (fscanf(fp, "%d %lf", &n, &ds) != 2 || n < 3 || !(ds > 0.0)) {
        fprintf(stderr, "ft_le: %s com n=%d ds_alvo=%g invalidos\n",
                arquivo, n, ds);
        fclose(fp);
        return NULL;
    }
    // Constroi a struct DIRETO, sem passar por ft_cria_curva: aquela reamostra,
    // e reamostrar ao retomar deslocaria todos os marcadores -- a retomada
    // deixaria de ser a continuacao da mesma corrida.
    ft_frente *f = (ft_frente *) calloc(1, sizeof *f);
    f->ds_alvo = (real) ds;
    _garante_cap(f, n);
    for (int i = 0; i < n; i++) {
        double x = 0.0, y = 0.0;
        if (fscanf(fp, "%lf %lf", &x, &y) != 2) {
            fprintf(stderr, "ft_le: %s truncado no marcador %d de %d\n",
                    arquivo, i, n);
            fclose(fp); ft_destroi(f);
            return NULL;
        }
        f->x[i][0] = (real) x;
        f->x[i][1] = (real) y;
    }
    f->n = n;
    fclose(fp);
    return f;
}
