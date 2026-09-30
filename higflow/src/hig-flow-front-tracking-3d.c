// Ver o .h para o contrato.  FASE B5, primeiro degrau: GEOMETRIA E ORACULOS,
// sem solver e sem MPI.  A superficie vive inteira aqui, replicada, como a
// frente 2D vivia no B1.

#include "hig-flow-front-tracking-3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ft3_superficie {
    Point *x;          // vertices
    int    nv, capv;
    int  (*tri)[3];    // triangulos, orientados para FORA
    int    nt, capt;
    real   ds_alvo;    // espacamento alvo de ARESTA
};

// --- vetores --------------------------------------------------------------

static void _sub(const real a[3], const real b[3], real r[3])
{ r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2]; }

static void _cruz(const real a[3], const real b[3], real r[3])
{
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

static real _ponto(const real a[3], const real b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

static real _norma(const real a[3]) { return sqrt(_ponto(a,a)); }

// Vetor normal do triangulo, com MODULO igual a DUAS vezes a area.  Sai assim
// de proposito: area e normal unitaria vem dos dois do mesmo calculo, sem
// repetir o produto vetorial.
static void _normal2(const ft3_superficie *s, int t, real n[3])
{
    real u[3], v[3];
    _sub(s->x[s->tri[t][1]], s->x[s->tri[t][0]], u);
    _sub(s->x[s->tri[t][2]], s->x[s->tri[t][0]], v);
    _cruz(u, v, n);
}

static real _area_tri(const ft3_superficie *s, int t)
{ real n[3]; _normal2(s, t, n); return 0.5 * _norma(n); }

// --- capacidade -----------------------------------------------------------

static void _cap_v(ft3_superficie *s, int nec)
{
    if (s->capv >= nec) return;
    int c = s->capv ? s->capv : 16;
    while (c < nec) c *= 2;
    s->x = (Point *) realloc(s->x, (size_t) c * sizeof(Point));
    if (!s->x) { fprintf(stderr, "ft3: sem memoria para %d vertices\n", c); abort(); }
    s->capv = c;
}

static void _cap_t(ft3_superficie *s, int nec)
{
    if (s->capt >= nec) return;
    int c = s->capt ? s->capt : 16;
    while (c < nec) c *= 2;
    s->tri = (int (*)[3]) realloc(s->tri, (size_t) c * sizeof *s->tri);
    if (!s->tri) { fprintf(stderr, "ft3: sem memoria para %d triangulos\n", c); abort(); }
    s->capt = c;
}

// --- arestas --------------------------------------------------------------

// Uma aresta da malha, com os (ate' dois) triangulos que a compartilham.  `k`
// e' o lado dentro do triangulo: a aresta vai de tri[t][k] a tri[t][(k+1)%3].
typedef struct { int v0, v1, t0, k0, t1, k1; } ft3_aresta;

static int _cmp_aresta(const void *A, const void *B)
{
    const int *a = (const int *) A, *b = (const int *) B;
    if (a[0] != b[0]) return a[0] < b[0] ? -1 : 1;
    if (a[1] != b[1]) return a[1] < b[1] ? -1 : 1;
    return 0;
}

// Constroi a tabela de arestas.  Devolve o vetor (o chamador libera) e escreve
// `*na`.  Aresta com mais de dois triangulos e' superficie nao-variedade: a
// funcao NAO tenta consertar, marca t1 = -2 e quem chamar decide.
static ft3_aresta *_arestas(const ft3_superficie *s, int *na)
{
    const int m = 3 * s->nt;
    int (*bruto)[4] = (int (*)[4]) malloc((size_t) m * sizeof *bruto);
    for (int t = 0; t < s->nt; t++)
        for (int k = 0; k < 3; k++) {
            int a = s->tri[t][k], b = s->tri[t][(k+1)%3];
            int i = 3*t + k;
            bruto[i][0] = a < b ? a : b;
            bruto[i][1] = a < b ? b : a;
            bruto[i][2] = t;
            bruto[i][3] = k;
        }
    qsort(bruto, (size_t) m, sizeof *bruto, _cmp_aresta);

    ft3_aresta *ar = (ft3_aresta *) malloc((size_t) m * sizeof *ar);
    int n = 0;
    for (int i = 0; i < m; ) {
        int j = i;
        while (j < m && bruto[j][0] == bruto[i][0] && bruto[j][1] == bruto[i][1]) j++;
        ar[n].v0 = bruto[i][0]; ar[n].v1 = bruto[i][1];
        ar[n].t0 = bruto[i][2]; ar[n].k0 = bruto[i][3];
        if (j - i == 1)      { ar[n].t1 = -1; ar[n].k1 = -1; }   // borda
        else if (j - i == 2) { ar[n].t1 = bruto[i+1][2]; ar[n].k1 = bruto[i+1][3]; }
        else                 { ar[n].t1 = -2; ar[n].k1 = -1; }   // nao-variedade
        n++;
        i = j;
    }
    free(bruto);
    *na = n;
    return ar;
}

// --- criacao --------------------------------------------------------------

static ft3_superficie *_vazia(void)
{
    ft3_superficie *s = (ft3_superficie *) calloc(1, sizeof *s);
    if (!s) { fprintf(stderr, "ft3: sem memoria\n"); abort(); }
    return s;
}

ft3_superficie *ft3_cria_malha(const Point *vert, int nv,
                               const int (*tri)[3], int nt, real ds_alvo)
{
    if (nv < 4 || nt < 4 || !vert || !tri) {
        fprintf(stderr, "ft3_cria_malha: %d vertices, %d triangulos\n", nv, nt);
        return NULL;
    }
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++)
            if (tri[t][k] < 0 || tri[t][k] >= nv) {
                fprintf(stderr, "ft3_cria_malha: triangulo %d refere vertice %d de %d\n",
                        t, tri[t][k], nv);
                return NULL;
            }
    ft3_superficie *s = _vazia();
    _cap_v(s, nv); _cap_t(s, nt);
    memcpy(s->x, vert, (size_t) nv * sizeof(Point));
    memcpy(s->tri, tri, (size_t) nt * sizeof *tri);
    s->nv = nv; s->nt = nt;

    if (ds_alvo > 0.0) {
        s->ds_alvo = ds_alvo;
    } else {
        int na; ft3_aresta *ar = _arestas(s, &na);
        real soma = 0.0;
        for (int e = 0; e < na; e++) {
            real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
            soma += _norma(d);
        }
        free(ar);
        s->ds_alvo = na ? soma / na : 0.0;
    }
    return s;
}

// Cache de ponto medio para a subdivisao: cada aresta gera UM vertice, e os
// dois triangulos vizinhos tem de reusar o mesmo -- senao a malha racha em
// vertices duplicados que parecem colados mas nao sao.
typedef struct { int a, b, m; } _meio;

static int _meio_de(ft3_superficie *s, _meio *cache, int *ncache,
                    int a, int b, const Point centro, real raio)
{
    int lo = a < b ? a : b, hi = a < b ? b : a;
    for (int i = 0; i < *ncache; i++)
        if (cache[i].a == lo && cache[i].b == hi) return cache[i].m;

    _cap_v(s, s->nv + 1);
    real p[3];
    for (int d = 0; d < 3; d++) p[d] = 0.5 * (s->x[a][d] + s->x[b][d]);
    // projeta na esfera: sem isto a subdivisao daria um poliedro cada vez mais
    // facetado, e nao uma esfera cada vez melhor
    real r[3]; _sub(p, centro, r);
    real nr = _norma(r);
    for (int d = 0; d < 3; d++) s->x[s->nv][d] = centro[d] + raio * r[d] / nr;

    cache[*ncache].a = lo; cache[*ncache].b = hi; cache[*ncache].m = s->nv;
    (*ncache)++;
    return s->nv++;
}

ft3_superficie *ft3_cria_esfera(const Point centro, real raio, int nsub)
{
    if (raio <= 0.0 || nsub < 0 || nsub > 7) {
        fprintf(stderr, "ft3_cria_esfera: raio=%g nsub=%d\n", (double) raio, nsub);
        return NULL;
    }
    const real t = (1.0 + sqrt(5.0)) / 2.0;
    real v0[12][3] = {
        {-1, t, 0}, { 1, t, 0}, {-1,-t, 0}, { 1,-t, 0},
        { 0,-1, t}, { 0, 1, t}, { 0,-1,-t}, { 0, 1,-t},
        { t, 0,-1}, { t, 0, 1}, {-t, 0,-1}, {-t, 0, 1}
    };
    // Faces do icosaedro, orientadas para FORA.
    const int f0[20][3] = {
        {0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},
        {1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
        {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},
        {4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}
    };
    ft3_superficie *s = _vazia();
    _cap_v(s, 12); _cap_t(s, 20);
    for (int i = 0; i < 12; i++) {
        real n = _norma(v0[i]);
        for (int d = 0; d < 3; d++) s->x[i][d] = centro[d] + raio * v0[i][d] / n;
    }
    s->nv = 12;
    memcpy(s->tri, f0, sizeof f0);
    s->nt = 20;

    for (int it = 0; it < nsub; it++) {
        const int nt0 = s->nt;
        // cada aresta da um vertice; arestas = 3*F/2
        const int maxc = 3 * nt0 / 2 + 8;
        _meio *cache = (_meio *) malloc((size_t) maxc * sizeof *cache);
        int ncache = 0;
        _cap_t(s, 4 * nt0);
        int (*novo)[3] = (int (*)[3]) malloc((size_t) (4*nt0) * sizeof *novo);
        int nn = 0;
        for (int f = 0; f < nt0; f++) {
            int a = s->tri[f][0], b = s->tri[f][1], c = s->tri[f][2];
            int ab = _meio_de(s, cache, &ncache, a, b, centro, raio);
            int bc = _meio_de(s, cache, &ncache, b, c, centro, raio);
            int ca = _meio_de(s, cache, &ncache, c, a, centro, raio);
            // orientacao preservada em cada um dos quatro
            novo[nn][0]=a;  novo[nn][1]=ab; novo[nn][2]=ca; nn++;
            novo[nn][0]=b;  novo[nn][1]=bc; novo[nn][2]=ab; nn++;
            novo[nn][0]=c;  novo[nn][1]=ca; novo[nn][2]=bc; nn++;
            novo[nn][0]=ab; novo[nn][1]=bc; novo[nn][2]=ca; nn++;
        }
        _cap_t(s, nn);
        memcpy(s->tri, novo, (size_t) nn * sizeof *novo);
        s->nt = nn;
        free(novo); free(cache);
    }

    // ds_alvo = media dos comprimentos de aresta da malha construida
    int na; ft3_aresta *ar = _arestas(s, &na);
    real soma = 0.0;
    for (int e = 0; e < na; e++) {
        real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
        soma += _norma(d);
    }
    free(ar);
    s->ds_alvo = na ? soma / na : 0.0;
    return s;
}

void ft3_destroi(ft3_superficie *s)
{
    if (!s) return;
    free(s->x); free(s->tri); free(s);
}

// --- consulta -------------------------------------------------------------

int  ft3_num_vertices(const ft3_superficie *s)   { return s ? s->nv : 0; }
int  ft3_num_triangulos(const ft3_superficie *s) { return s ? s->nt : 0; }
real ft3_ds_alvo(const ft3_superficie *s)        { return s ? s->ds_alvo : 0.0; }
const Point *ft3_posicoes(const ft3_superficie *s) { return s ? (const Point *) s->x : NULL; }
const int (*ft3_triangulos(const ft3_superficie *s))[3]
{ return s ? (const int (*)[3]) s->tri : NULL; }

real ft3_area(const ft3_superficie *s)
{
    if (!s) return 0.0;
    real a = 0.0;
    for (int t = 0; t < s->nt; t++) a += _area_tri(s, t);
    return a;
}

real ft3_volume(const ft3_superficie *s)
{
    if (!s) return 0.0;
    // Teorema do divergente: V = (1/3) int x.n dA = (1/6) sum det[x0,x1,x2].
    // O resultado NAO depende da origem numa superficie fechada; usar o
    // primeiro vertice como origem reduz o cancelamento catastrofico quando a
    // superficie esta' longe da origem.
    real v = 0.0;
    const real *o = s->x[0];
    for (int t = 0; t < s->nt; t++) {
        real a[3], b[3], c[3], cr[3];
        _sub(s->x[s->tri[t][0]], o, a);
        _sub(s->x[s->tri[t][1]], o, b);
        _sub(s->x[s->tri[t][2]], o, c);
        _cruz(b, c, cr);
        v += _ponto(a, cr);
    }
    return v / 6.0;
}

int ft3_euler(const ft3_superficie *s)
{
    if (!s) return 0;
    int na; ft3_aresta *ar = _arestas(s, &na);
    free(ar);
    return s->nv - na + s->nt;
}

void ft3_qualidade(const ft3_superficie *s, real *amin, real *amax, real *pior)
{
    if (!s || s->nt == 0) {
        if (amin) *amin = 0.0;
        if (amax) *amax = 0.0;
        if (pior) *pior = 0.0;
        return;
    }
    int na; ft3_aresta *ar = _arestas(s, &na);
    real lo = 1e300, hi = 0.0;
    for (int e = 0; e < na; e++) {
        real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
        real L = _norma(d);
        if (L < lo) lo = L;
        if (L > hi) hi = L;
    }
    free(ar);
    real pa = 1.0;
    for (int t = 0; t < s->nt; t++) {
        real L[3];
        for (int k = 0; k < 3; k++) {
            real d[3];
            _sub(s->x[s->tri[t][(k+1)%3]], s->x[s->tri[t][k]], d);
            L[k] = _norma(d);
        }
        real mn = L[0] < L[1] ? (L[0] < L[2] ? L[0] : L[2]) : (L[1] < L[2] ? L[1] : L[2]);
        real mx = L[0] > L[1] ? (L[0] > L[2] ? L[0] : L[2]) : (L[1] > L[2] ? L[1] : L[2]);
        if (mn > 0.0 && mx / mn > pa) pa = mx / mn;
    }
    if (amin) *amin = lo;
    if (amax) *amax = hi;
    if (pior) *pior = pa;
}

// --- tensao superficial ---------------------------------------------------

void ft3_forcas_tensao(const ft3_superficie *s, real sigma,
                       Point *pos, Point *forca, real *peso)
{
    if (!s) return;
    for (int i = 0; i < s->nv; i++) {
        for (int d = 0; d < 3; d++) { pos[i][d] = s->x[i][d]; forca[i][d] = 0.0; }
        peso[i] = 0.0;
    }
    for (int t = 0; t < s->nt; t++) {
        real n2[3];
        _normal2(s, t, n2);                 // modulo = 2*area
        real a2 = _norma(n2);
        if (a2 <= 0.0) continue;            // triangulo degenerado: nao contribui
        real n[3] = { n2[0]/a2, n2[1]/a2, n2[2]/a2 };
        const real area = 0.5 * a2;

        for (int k = 0; k < 3; k++) {
            int ia = s->tri[t][k], ib = s->tri[t][(k+1)%3];
            real e[3], f[3];
            _sub(s->x[ib], s->x[ia], e);
            // f = sigma * (b - a) x n : o vetor tangente a' superficie e
            // NORMAL a' aresta, apontando para FORA do triangulo.  A soma dos
            // tres da zero num triangulo plano -- a forca vem da diferenca de
            // normais entre triangulos vizinhos da mesma aresta.
            _cruz(e, n, f);
            // SINAL, derivado e nao adivinhado.  `f` e' o puxao que a vizinhanca
            // exerce ATRAVES da aresta, e somar as duas arestas incidentes a um
            // vertice da' exatamente +sigma * grad_a(Area) -- conferido contra o
            // gradiente numerico, que bate casa a casa.  A forca de tensao e'
            // MENOS esse gradiente, porque tensao MINIMIZA area:
            //
            //     F_a = -sigma * grad_a A = -(sigma/2) (b - c) x n .
            //
            // Sem o sinal a forca aponta para FORA e a gota explodiria em vez de
            // contrair.  O oraculo pega (todos os vertices para fora), mas a
            // razao de estar certo e' a variacional, nao o oraculo.
            for (int d = 0; d < 3; d++) {
                forca[ia][d] -= 0.5 * sigma * f[d];
                forca[ib][d] -= 0.5 * sigma * f[d];
            }
        }
        // peso do vertice: um terco da area de cada triangulo incidente.  E' o
        // peso que a fronteira imersa espera para superficie em 3D (codim 1).
        for (int k = 0; k < 3; k++) peso[s->tri[t][k]] += area / 3.0;
    }
}

// --- cirurgia -------------------------------------------------------------

// Vertice oposto a' aresta (v0,v1) dentro do triangulo t.
static int _oposto(const ft3_superficie *s, int t, int v0, int v1)
{
    for (int k = 0; k < 3; k++) {
        int v = s->tri[t][k];
        if (v != v0 && v != v1) return v;
    }
    return -1;
}

// Normal (nao unitaria) que o triangulo TERIA se `de` fosse trocado por `para`.
static void _normal2_com(const ft3_superficie *s, int t, int de, const real para[3],
                         real n[3])
{
    const real *p[3];
    for (int k = 0; k < 3; k++)
        p[k] = (s->tri[t][k] == de) ? para : s->x[s->tri[t][k]];
    real u[3], v[3];
    _sub(p[1], p[0], u);
    _sub(p[2], p[0], v);
    _cruz(u, v, n);
}

// Divide a aresta e no ponto medio.  Os dois triangulos viram quatro.
static void _divide(ft3_superficie *s, const ft3_aresta *e)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0) return;

    _cap_v(s, s->nv + 1);
    int m = s->nv++;
    for (int k = 0; k < 3; k++) s->x[m][k] = 0.5 * (s->x[a][k] + s->x[b][k]);

    // Reconstroi os quatro respeitando a ORIENTACAO de cada original: em t0 a
    // aresta vai de tri[t0][k0] para o seguinte, e e' isso que diz quem e' "a"
    // no sentido de percurso.
    int a0 = s->tri[e->t0][e->k0];                 // primeiro da aresta em t0
    int b0 = s->tri[e->t0][(e->k0+1)%3];
    _cap_t(s, s->nt + 2);
    s->tri[e->t0][0] = a0; s->tri[e->t0][1] = m;  s->tri[e->t0][2] = c;
    s->tri[s->nt][0] = m;  s->tri[s->nt][1] = b0; s->tri[s->nt][2] = c; s->nt++;

    int a1 = s->tri[e->t1][e->k1];
    int b1 = s->tri[e->t1][(e->k1+1)%3];
    s->tri[e->t1][0] = a1; s->tri[e->t1][1] = m;  s->tri[e->t1][2] = d;
    s->tri[s->nt][0] = m;  s->tri[s->nt][1] = b1; s->tri[s->nt][2] = d; s->nt++;
}

// Colapsa a aresta (a,b) no ponto medio.  RECUSA se violar a condicao de elo ou
// se inverter algum triangulo -- malha pior e' recuperavel, nao-variedade nao.
// Devolve 1 se colapsou.
static int _colapsa(ft3_superficie *s, const ft3_aresta *e, char *morto)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0 || c == d) return 0;

    // CONDICAO DE ELO: os vizinhos comuns de a e b tem de ser EXATAMENTE {c,d}.
    // Se houver um terceiro, o colapso cria aresta dupla e a superficie deixa
    // de ser variedade.
    int viz_a[256], nva = 0, viz_b[256], nvb = 0;
    for (int t = 0; t < s->nt; t++) {
        if (morto[t]) continue;
        int tem_a = 0, tem_b = 0;
        for (int k = 0; k < 3; k++) {
            if (s->tri[t][k] == a) tem_a = 1;
            if (s->tri[t][k] == b) tem_b = 1;
        }
        for (int k = 0; k < 3; k++) {
            int v = s->tri[t][k];
            if (v == a || v == b) continue;
            if (tem_a && nva < 256) { int j=0; while (j<nva && viz_a[j]!=v) j++; if (j==nva) viz_a[nva++]=v; }
            if (tem_b && nvb < 256) { int j=0; while (j<nvb && viz_b[j]!=v) j++; if (j==nvb) viz_b[nvb++]=v; }
        }
    }
    int comuns = 0;
    for (int i = 0; i < nva; i++)
        for (int j = 0; j < nvb; j++)
            if (viz_a[i] == viz_b[j]) comuns++;
    if (comuns != 2) return 0;

    real novo[3];
    for (int k = 0; k < 3; k++) novo[k] = 0.5 * (s->x[a][k] + s->x[b][k]);

    // Nenhum triangulo sobrevivente pode INVERTER.  Comparar o sinal do produto
    // escalar entre a normal antiga e a nova pega inversao sem depender de
    // escala; area indo a zero tambem e' recusa.
    for (int t = 0; t < s->nt; t++) {
        if (morto[t] || t == e->t0 || t == e->t1) continue;
        int tem = 0, qual = -1;
        for (int k = 0; k < 3; k++)
            if (s->tri[t][k] == a || s->tri[t][k] == b) { tem = 1; qual = s->tri[t][k]; }
        if (!tem) continue;
        real n_v[3], n_n[3];
        _normal2(s, t, n_v);
        _normal2_com(s, t, qual, novo, n_n);
        if (_ponto(n_v, n_n) <= 0.0) return 0;
        if (_norma(n_n) < 1e-12 * _norma(n_v)) return 0;
    }

    // Aplica: a recebe o ponto medio, b some, os dois triangulos da aresta morrem.
    for (int k = 0; k < 3; k++) s->x[a][k] = novo[k];
    morto[e->t0] = 1;
    morto[e->t1] = 1;
    for (int t = 0; t < s->nt; t++) {
        if (morto[t]) continue;
        for (int k = 0; k < 3; k++) if (s->tri[t][k] == b) s->tri[t][k] = a;
    }
    return 1;
}

// Menor angulo (em cosseno: MAIOR cosseno = menor angulo) de um triangulo dado
// por tres pontos.  Devolve o cosseno do maior angulo interno -- quanto MENOR,
// melhor o triangulo.  Nao chama acos: comparar cossenos basta e e' monotono.
static real _pior_cos(const real *p0, const real *p1, const real *p2)
{
    const real *P[3] = {p0, p1, p2};
    real pior = -1.0;
    for (int k = 0; k < 3; k++) {
        real u[3], v[3];
        _sub(P[(k+1)%3], P[k], u);
        _sub(P[(k+2)%3], P[k], v);
        real nu = _norma(u), nv = _norma(v);
        if (nu <= 0.0 || nv <= 0.0) return 1.0;          // degenerado: o pior possivel
        real c = _ponto(u, v) / (nu * nv);
        if (c > pior) pior = c;
    }
    return pior;
}

// Gira a aresta (a,b) para a diagonal (c,d), se isso melhorar o PIOR ANGULO e
// os dois triangulos forem quase coplanares.
//
// A COPLANARIDADE NAO E' ZELO: girar atraves de uma quina DOBRA o triangulo
// sobre a superficie e muda a forma, em vez de melhorar a malha.  Medido: com o
// criterio anterior (girar se a diagonal nova for mais curta) o giro era a UNICA
// das tres operacoes que produzia triangulo invertido -- divisao e colapso
// davam zero, e qualquer combinacao contendo o giro dava de 1 a 15.
static int _gira(ft3_superficie *s, const ft3_aresta *e, const char *morto)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0 || c == d) return 0;

    // (c,d) ja' existir torna o giro uma aresta dupla.
    for (int t = 0; t < s->nt; t++) {
        if (morto[t]) continue;
        int tc = 0, td = 0;
        for (int k = 0; k < 3; k++) {
            if (s->tri[t][k] == c) tc = 1;
            if (s->tri[t][k] == d) td = 1;
        }
        if (tc && td) return 0;
    }

    real nv0[3], nv1[3];
    _normal2(s, e->t0, nv0);
    _normal2(s, e->t1, nv1);
    real m0 = _norma(nv0), m1 = _norma(nv1);
    if (m0 <= 0.0 || m1 <= 0.0) return 0;
    real u0[3] = {nv0[0]/m0, nv0[1]/m0, nv0[2]/m0};
    real u1[3] = {nv1[0]/m1, nv1[1]/m1, nv1[2]/m1};

    // COPLANARIDADE: normais dentro de ~25 graus uma da outra.
    const real COS_LIM = 0.90;
    if (_ponto(u0, u1) < COS_LIM) return 0;

    // Melhora do pior angulo.  Cosseno MAIOR = angulo pior, entao o novo pior
    // cosseno tem de ser MENOR que o velho.
    const real *pa = s->x[a], *pb = s->x[b], *pc = s->x[c], *pd = s->x[d];
    real velho = _pior_cos(pa, pb, pc);
    real w = _pior_cos(pb, pa, pd); if (w > velho) velho = w;
    real novo  = _pior_cos(pa, pd, pc);
    w = _pior_cos(pb, pc, pd); if (w > novo) novo = w;
    if (novo >= velho) return 0;

    // Normais novas alinhadas com a media das velhas -- pega dobra que o teste
    // de coplanaridade sozinho deixaria passar.
    real med[3] = {u0[0]+u1[0], u0[1]+u1[1], u0[2]+u1[2]};
    real mm = _norma(med);
    if (mm <= 0.0) return 0;
    for (int k = 0; k < 3; k++) med[k] /= mm;

    real uu[3], vv[3], n0[3], n1[3];
    _sub(pd, pa, uu); _sub(pc, pa, vv); _cruz(uu, vv, n0);
    _sub(pc, pb, uu); _sub(pd, pb, vv); _cruz(uu, vv, n1);
    real q0 = _norma(n0), q1 = _norma(n1);
    if (q0 <= 0.0 || q1 <= 0.0) return 0;
    for (int k = 0; k < 3; k++) { n0[k] /= q0; n1[k] /= q1; }
    if (_ponto(n0, med) < COS_LIM || _ponto(n1, med) < COS_LIM) return 0;

    s->tri[e->t0][0] = a; s->tri[e->t0][1] = d; s->tri[e->t0][2] = c;
    s->tri[e->t1][0] = b; s->tri[e->t1][1] = c; s->tri[e->t1][2] = d;
    return 1;
}

// Remove triangulos mortos e vertices que ficaram sem uso.
static void _compacta(ft3_superficie *s, const char *morto)
{
    int nt = 0;
    for (int t = 0; t < s->nt; t++) {
        if (morto[t]) continue;
        if (nt != t) memcpy(s->tri[nt], s->tri[t], sizeof s->tri[t]);
        nt++;
    }
    s->nt = nt;

    int *mapa = (int *) malloc((size_t) s->nv * sizeof *mapa);
    for (int i = 0; i < s->nv; i++) mapa[i] = -1;
    for (int t = 0; t < s->nt; t++)
        for (int k = 0; k < 3; k++) mapa[s->tri[t][k]] = 0;
    int nv = 0;
    for (int i = 0; i < s->nv; i++)
        if (mapa[i] == 0) {
            if (nv != i) memcpy(s->x[nv], s->x[i], sizeof s->x[i]);
            mapa[i] = nv++;
        }
    for (int t = 0; t < s->nt; t++)
        for (int k = 0; k < 3; k++) s->tri[t][k] = mapa[s->tri[t][k]];
    s->nv = nv;
    free(mapa);
}

// Quais etapas rodar.  FT3_CIRURGIA_ETAPAS seleciona um subconjunto de "dcg"
// (divide, colapsa, gira) -- gancho de DIAGNOSTICO, para atribuir um defeito a
// uma etapa em vez de adivinhar qual das tres o causou.  Sem a variavel, roda
// as tres.
static int _etapa(char c)
{
    const char *e = getenv("FT3_CIRURGIA_ETAPAS");
    return e ? (strchr(e, c) != NULL) : 1;
}

// Uma rodada de divide/colapsa/gira.  Devolve quantas operacoes fez.
static int _rodada(ft3_superficie *s)
{
    const real longa = 2.0 * s->ds_alvo;
    const real curta = 0.5 * s->ds_alvo;
    int feitas = 0;

    if (_etapa('d'))
    // 1) DIVIDE.  Uma passada por vez: dividir muda a tabela de arestas, entao
    // ela e' reconstruida.  O laco externo para quando nada mais se divide.
    for (int passe = 0; passe < 8; passe++) {
        int na; ft3_aresta *ar = _arestas(s, &na);
        int fez = 0;
        char *usado = (char *) calloc((size_t) s->nt, 1);
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;                    // borda ou nao-variedade
            if (usado[ar[e].t0] || usado[ar[e].t1]) continue;
            real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
            if (_norma(d) <= longa) continue;
            usado[ar[e].t0] = usado[ar[e].t1] = 1;
            _divide(s, &ar[e]);
            fez++; feitas++;
        }
        free(usado); free(ar);
        if (!fez) break;
    }

    if (_etapa('c'))
    // 2) COLAPSA.  Uma aresta por vez, com a tabela reconstruida a cada passada
    // -- colapso muda vizinhanca longe da aresta tratada.
    for (int passe = 0; passe < 8; passe++) {
        int na; ft3_aresta *ar = _arestas(s, &na);
        char *morto = (char *) calloc((size_t) s->nt, 1);
        // A TABELA DE ARESTAS ENVELHECE A CADA COLAPSO.  Colapsar (a,b) faz `b`
        // sumir de todos os triangulos, mas as arestas seguintes da tabela ainda
        // o referenciam -- e `_oposto` num triangulo que ja' nao contem `b`
        // devolve vertice errado, que e' como triangulo invertido nascia aqui.
        // Marcar os vertices tocados e pular quem os cita mantem a passada
        // valida sem reconstruir a tabela a cada operacao.
        char *tocado = (char *) calloc((size_t) s->nv, 1);
        int fez = 0;
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;
            if (morto[ar[e].t0] || morto[ar[e].t1]) continue;
            if (tocado[ar[e].v0] || tocado[ar[e].v1]) continue;
            real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
            if (_norma(d) >= curta) continue;
            int a = ar[e].v0, b = ar[e].v1;
            if (_colapsa(s, &ar[e], morto)) {
                fez++; feitas++;
                // tudo que compartilha triangulo com a ou b sai desta passada
                for (int t = 0; t < s->nt; t++) {
                    if (morto[t]) continue;
                    int toca = 0;
                    for (int k = 0; k < 3; k++)
                        if (s->tri[t][k] == a || s->tri[t][k] == b) toca = 1;
                    if (toca) for (int k = 0; k < 3; k++) tocado[s->tri[t][k]] = 1;
                }
                tocado[a] = tocado[b] = 1;
            }
        }
        free(tocado);
        if (fez) _compacta(s, morto);
        free(morto); free(ar);
        if (!fez) break;
    }

    if (_etapa('g'))
    // 3) GIRA para qualidade.
    {
        int na; ft3_aresta *ar = _arestas(s, &na);
        char *morto = (char *) calloc((size_t) s->nt, 1);
        char *usado = (char *) calloc((size_t) s->nt, 1);
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;
            if (usado[ar[e].t0] || usado[ar[e].t1]) continue;
            if (_gira(s, &ar[e], morto)) {
                usado[ar[e].t0] = usado[ar[e].t1] = 1;
                feitas++;
            }
        }
        free(usado); free(morto); free(ar);
    }
    return feitas;
}

int ft3_cirurgia(ft3_superficie *s)
{
    if (!s || s->ds_alvo <= 0.0) return 0;
    // O CICLO SE REPETE porque as tres operacoes interferem: colapsar move um
    // vertice para o ponto medio e ALONGA as arestas vizinhas, e girar troca
    // qual aresta e' a longa.  Uma passada so' deixa aresta acima do teto --
    // medido: 0,157 contra teto 0,151.  Repete ate' nada mudar.
    int total = 0;
    for (int r = 0; r < 12; r++) {
        int n = _rodada(s);
        total += n;
        if (n == 0) break;
    }
    return total;
}

// --- persistencia e saida -------------------------------------------------

#define FT3_CABECALHO "ft3 1"

int ft3_grava(const ft3_superficie *s, const char *caminho)
{
    if (!s || !caminho) return -1;
    FILE *fp = fopen(caminho, "w");
    if (!fp) return -1;
    fprintf(fp, "%s %d %d %.17g\n", FT3_CABECALHO, s->nv, s->nt, (double) s->ds_alvo);
    for (int i = 0; i < s->nv; i++)
        fprintf(fp, "%.17g %.17g %.17g\n",
                (double) s->x[i][0], (double) s->x[i][1], (double) s->x[i][2]);
    for (int t = 0; t < s->nt; t++)
        fprintf(fp, "%d %d %d\n", s->tri[t][0], s->tri[t][1], s->tri[t][2]);
    int ok = (fflush(fp) == 0);
    fclose(fp);
    return ok ? 0 : -1;
}

ft3_superficie *ft3_le(const char *caminho)
{
    if (!caminho) return NULL;
    FILE *fp = fopen(caminho, "r");
    if (!fp) return NULL;
    char c0[8], c1[8];
    int nv, nt; double ds;
    if (fscanf(fp, "%7s %7s %d %d %lf", c0, c1, &nv, &nt, &ds) != 5
        || strcmp(c0, "ft3") != 0 || strcmp(c1, "1") != 0) {
        fprintf(stderr, "ft3_le: %s nao tem o cabecalho esperado\n", caminho);
        fclose(fp); return NULL;
    }
    if (nv < 4 || nt < 4 || !(ds > 0.0)) {
        fprintf(stderr, "ft3_le: %s com nv=%d nt=%d ds=%g invalidos\n", caminho, nv, nt, ds);
        fclose(fp); return NULL;
    }
    ft3_superficie *s = _vazia();
    _cap_v(s, nv); _cap_t(s, nt);
    for (int i = 0; i < nv; i++) {
        double a, b, c;
        if (fscanf(fp, "%lf %lf %lf", &a, &b, &c) != 3) {
            fprintf(stderr, "ft3_le: %s truncado no vertice %d de %d\n", caminho, i, nv);
            ft3_destroi(s); fclose(fp); return NULL;
        }
        s->x[i][0] = a; s->x[i][1] = b; s->x[i][2] = c;
    }
    for (int t = 0; t < nt; t++) {
        int a, b, c;
        if (fscanf(fp, "%d %d %d", &a, &b, &c) != 3) {
            fprintf(stderr, "ft3_le: %s truncado no triangulo %d de %d\n", caminho, t, nt);
            ft3_destroi(s); fclose(fp); return NULL;
        }
        if (a < 0 || a >= nv || b < 0 || b >= nv || c < 0 || c >= nv) {
            fprintf(stderr, "ft3_le: %s, triangulo %d refere vertice fora de faixa\n", caminho, t);
            ft3_destroi(s); fclose(fp); return NULL;
        }
        s->tri[t][0] = a; s->tri[t][1] = b; s->tri[t][2] = c;
    }
    fclose(fp);
    s->nv = nv; s->nt = nt; s->ds_alvo = ds;
    return s;
}

int ft3_escreve_vtk(const ft3_superficie *s, const char *caminho)
{
    if (!s || !caminho) return -1;
    FILE *fp = fopen(caminho, "w");
    if (!fp) return -1;
    fprintf(fp, "# vtk DataFile Version 3.0\nfrente 3d\nASCII\nDATASET POLYDATA\n");
    fprintf(fp, "POINTS %d double\n", s->nv);
    for (int i = 0; i < s->nv; i++)
        fprintf(fp, "%.8f %.8f %.8f\n",
                (double) s->x[i][0], (double) s->x[i][1], (double) s->x[i][2]);
    fprintf(fp, "POLYGONS %d %d\n", s->nt, 4 * s->nt);
    for (int t = 0; t < s->nt; t++)
        fprintf(fp, "3 %d %d %d\n", s->tri[t][0], s->tri[t][1], s->tri[t][2]);
    int ok = (fflush(fp) == 0);
    fclose(fp);
    return ok ? 0 : -1;
}
