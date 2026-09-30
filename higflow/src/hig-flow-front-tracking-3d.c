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

// Incidencia vertice -> triangulos, em formato comprimido (ini/lst).  SEM isto
// a condicao de elo e a checagem de aresta duplicada varrem todos os triangulos
// por candidato, e a cirurgia fica O(nt^2) -- inviavel quando ela roda a cada
// passo de tempo, que e' o caso do B1.
static void _incidencia(const ft3_superficie *s, int **ini, int **lst)
{
    int *c = (int *) calloc((size_t) s->nv + 1, sizeof *c);
    for (int t = 0; t < s->nt; t++)
        for (int k = 0; k < 3; k++) c[s->tri[t][k] + 1]++;
    for (int i = 0; i < s->nv; i++) c[i+1] += c[i];
    int *L = (int *) malloc((size_t) (3 * s->nt) * sizeof *L);
    int *pos = (int *) malloc((size_t) s->nv * sizeof *pos);
    memcpy(pos, c, (size_t) s->nv * sizeof *pos);
    for (int t = 0; t < s->nt; t++)
        for (int k = 0; k < 3; k++) L[pos[s->tri[t][k]]++] = t;
    free(pos);
    *ini = c; *lst = L;
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

// --- adveccao -------------------------------------------------------------

void ft3_advecta(ft3_superficie *s, ft3_campo_u u, void *ctx, real t, real dt)
{
    // RK2 do ponto medio, por vertice:
    //   k1 = u(x, t);   xm = x + (dt/2) k1;   k2 = u(xm, t + dt/2);   x += dt k2
    // A superficie e' replicada, entao cada vertice anda sozinho -- a
    // conectividade nao entra aqui, so' na cirurgia.
    if (!s) return;
    for (int i = 0; i < s->nv; i++) {
        real k1[3], k2[3];
        Point xm;
        u(s->x[i], t, ctx, k1);
        for (int d = 0; d < 3; d++) xm[d] = s->x[i][d] + 0.5 * dt * k1[d];
        u(xm, t + 0.5 * dt, ctx, k2);
        for (int d = 0; d < 3; d++) s->x[i][d] += dt * k2[d];
    }
}

// --- cirurgia -------------------------------------------------------------

// Normais nos VERTICES, media ponderada por area dos triangulos incidentes.
// Precisas o bastante para posicionar ponto novo; nao entram na forca, que sai
// da integral de linha e nao usa normal de vertice.
static void _normais_vertices(const ft3_superficie *s, Point *nv)
{
    for (int i = 0; i < s->nv; i++) { nv[i][0] = nv[i][1] = nv[i][2] = 0.0; }
    for (int t = 0; t < s->nt; t++) {
        real n[3]; _normal2(s, t, n);       // modulo = 2*area: ja' e' o peso
        for (int k = 0; k < 3; k++)
            for (int d = 0; d < 3; d++) nv[s->tri[t][k]][d] += n[d];
    }
    for (int i = 0; i < s->nv; i++) {
        real m = _norma(nv[i]);
        if (m > 0.0) for (int d = 0; d < 3; d++) nv[i][d] /= m;
    }
}

// Ponto medio SOBRE A SUPERFICIE, e nao sobre a corda.
//
// POR QUE ISTO EXISTE.  Inserir no meio da corda puxa o ponto para DENTRO numa
// superficie convexa, e cada divisao perde um pouco de volume.  Medido no B1
// 3D: com ponto medio de corda, ~1500 operacoes ao longo da corrida custavam
// 5,2% do volume, e refinar o passo de tempo pela metade nao mudava nada --
// prova de que o erro era da cirurgia e nao do integrador.
//
// A colocacao e' o meio da curva de Hermite cubica que interpola os dois
// extremos COM AS NORMAIS, a mesma dos "PN triangles":
//
//     m = (a+b)/2 - (1/8) [ ((b-a).n_a) n_a + ((a-b).n_b) n_b ]
//
// Na esfera isso empurra para fora exatamente na direcao que a corda tinha
// cortado.  NAO e' correcao de volume imposta a posteriori: o volume continua
// livre para errar, e o oraculo continua valendo.
static void _meio_curvo(const real *a, const real *b,
                        const real *na, const real *nb, real m[3])
{
    real ab[3], ba[3];
    _sub(b, a, ab);
    _sub(a, b, ba);
    real wa = _ponto(ab, na), wb = _ponto(ba, nb);
    for (int d = 0; d < 3; d++)
        m[d] = 0.5 * (a[d] + b[d]) - 0.125 * (wa * na[d] + wb * nb[d]);
}

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
static void _divide(ft3_superficie *s, const ft3_aresta *e, const Point *nvert)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0) return;

    real pm[3];
    _meio_curvo(s->x[a], s->x[b], nvert[a], nvert[b], pm);
    _cap_v(s, s->nv + 1);
    int m = s->nv++;
    for (int k = 0; k < 3; k++) s->x[m][k] = pm[k];

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
static int _colapsa(ft3_superficie *s, const ft3_aresta *e, char *morto,
                    const int *ini, const int *lst, const Point *nvert)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0 || c == d) return 0;

    // CONDICAO DE ELO: os vizinhos comuns de a e b tem de ser EXATAMENTE {c,d}.
    // Se houver um terceiro, o colapso cria aresta dupla e a superficie deixa
    // de ser variedade.  Percorre so' os triangulos INCIDENTES, nao todos.
    int viz_a[256], nva = 0, viz_b[256], nvb = 0;
    for (int p = ini[a]; p < ini[a+1]; p++) {
        int t = lst[p];
        if (morto[t]) continue;
        for (int k = 0; k < 3; k++) {
            int v = s->tri[t][k];
            if (v == a || v == b) continue;
            if (nva < 256) { int j=0; while (j<nva && viz_a[j]!=v) j++; if (j==nva) viz_a[nva++]=v; }
        }
    }
    for (int p = ini[b]; p < ini[b+1]; p++) {
        int t = lst[p];
        if (morto[t]) continue;
        for (int k = 0; k < 3; k++) {
            int v = s->tri[t][k];
            if (v == a || v == b) continue;
            if (nvb < 256) { int j=0; while (j<nvb && viz_b[j]!=v) j++; if (j==nvb) viz_b[nvb++]=v; }
        }
    }
    int comuns = 0;
    for (int i = 0; i < nva; i++)
        for (int j = 0; j < nvb; j++)
            if (viz_a[i] == viz_b[j]) comuns++;
    if (comuns != 2) return 0;

    real novo[3];
    _meio_curvo(s->x[a], s->x[b], nvert[a], nvert[b], novo);

    // Nenhum triangulo sobrevivente pode INVERTER.  Comparar o sinal do produto
    // escalar entre a normal antiga e a nova pega inversao sem depender de
    // escala; area indo a zero tambem e' recusa.
    for (int pp = ini[a]; pp < ini[b+1]; pp++) {
        if (pp >= ini[a+1] && pp < ini[b]) continue;   // so' as duas faixas
        int t = lst[pp];
        if (morto[t] || t == e->t0 || t == e->t1) continue;
        int qual = -1;
        for (int k = 0; k < 3; k++)
            if (s->tri[t][k] == a || s->tri[t][k] == b) qual = s->tri[t][k];
        if (qual < 0) continue;
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
    for (int p = ini[b]; p < ini[b+1]; p++) {
        int t = lst[p];
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
static int _gira(ft3_superficie *s, const ft3_aresta *e, const char *morto,
                 const int *ini, const int *lst)
{
    int a = e->v0, b = e->v1;
    int c = _oposto(s, e->t0, a, b);
    int d = _oposto(s, e->t1, a, b);
    if (c < 0 || d < 0 || c == d) return 0;

    // (c,d) ja' existir torna o giro uma aresta dupla.  Basta olhar os
    // triangulos incidentes a c.
    for (int p = ini[c]; p < ini[c+1]; p++) {
        int t = lst[p];
        if (morto[t]) continue;
        for (int k = 0; k < 3; k++) if (s->tri[t][k] == d) return 0;
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
        Point *nvert = (Point *) malloc((size_t) s->nv * sizeof *nvert);
        _normais_vertices(s, nvert);
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;                    // borda ou nao-variedade
            if (usado[ar[e].t0] || usado[ar[e].t1]) continue;
            real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
            if (_norma(d) <= longa) continue;
            usado[ar[e].t0] = usado[ar[e].t1] = 1;
            _divide(s, &ar[e], nvert);
            fez++; feitas++;
        }
        free(nvert); free(usado); free(ar);
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
        int *ini, *lst; _incidencia(s, &ini, &lst);
        Point *nvert = (Point *) malloc((size_t) s->nv * sizeof *nvert);
        _normais_vertices(s, nvert);
        int fez = 0;
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;
            if (morto[ar[e].t0] || morto[ar[e].t1]) continue;
            if (tocado[ar[e].v0] || tocado[ar[e].v1]) continue;
            real d[3]; _sub(s->x[ar[e].v1], s->x[ar[e].v0], d);
            if (_norma(d) >= curta) continue;
            int a = ar[e].v0, b = ar[e].v1;
            if (_colapsa(s, &ar[e], morto, ini, lst, nvert)) {
                fez++; feitas++;
                // tudo que compartilha triangulo com a ou b sai desta passada
                int viz[2] = {a, b};
                for (int j = 0; j < 2; j++)
                    for (int pp = ini[viz[j]]; pp < ini[viz[j]+1]; pp++) {
                        int t = lst[pp];
                        if (morto[t]) continue;
                        for (int k = 0; k < 3; k++) tocado[s->tri[t][k]] = 1;
                    }
                tocado[a] = tocado[b] = 1;
            }
        }
        free(nvert); free(ini); free(lst); free(tocado);
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
        // O giro tambem ENVELHECE a incidencia: apos girar (a,b) para (c,d), os
        // quatro vertices mudam de vizinhanca, e a checagem de aresta duplicada
        // de um giro seguinte leria estado velho e poderia criar aresta dupla.
        char *tocado = (char *) calloc((size_t) s->nv, 1);
        int *ini, *lst; _incidencia(s, &ini, &lst);
        for (int e = 0; e < na; e++) {
            if (ar[e].t1 < 0) continue;
            if (usado[ar[e].t0] || usado[ar[e].t1]) continue;
            if (tocado[ar[e].v0] || tocado[ar[e].v1]) continue;
            int c = _oposto(s, ar[e].t0, ar[e].v0, ar[e].v1);
            int d = _oposto(s, ar[e].t1, ar[e].v0, ar[e].v1);
            if (c < 0 || d < 0) continue;
            if (tocado[c] || tocado[d]) continue;
            if (_gira(s, &ar[e], morto, ini, lst)) {
                usado[ar[e].t0] = usado[ar[e].t1] = 1;
                tocado[ar[e].v0] = tocado[ar[e].v1] = tocado[c] = tocado[d] = 1;
                feitas++;
            }
        }
        free(ini); free(lst); free(tocado); free(usado); free(morto); free(ar);
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


// --- fracao de volume -----------------------------------------------------

// Caixa envolvente da superficie, com folga zero.  Caminho rapido de tudo que
// se segue: ponto fora dela esta' fora da superficie, sem consultar triangulo.
static void _envoltoria(const ft3_superficie *s, real lo[3], real hi[3])
{
    for (int d = 0; d < 3; d++) { lo[d] = s->x[0][d]; hi[d] = s->x[0][d]; }
    for (int i = 1; i < s->nv; i++)
        for (int d = 0; d < 3; d++) {
            if (s->x[i][d] < lo[d]) lo[d] = s->x[i][d];
            if (s->x[i][d] > hi[d]) hi[d] = s->x[i][d];
        }
}

// Cruzamento do raio +x a partir de `x` com o triangulo t.  Mesma conta da
// versao sem grade -- extraida para as duas usarem UMA implementacao, porque
// duas copias divergem com o tempo e o arreio deixaria de valer.
static int _cruza_raio(const ft3_superficie *s, int t, const Point x)
{
    const real *a = s->x[s->tri[t][0]];
    const real *b = s->x[s->tri[t][1]];
    const real *c = s->x[s->tri[t][2]];
    real d1 = (b[1]-a[1])*(c[2]-a[2]) - (c[1]-a[1])*(b[2]-a[2]);
    if (d1 == 0.0) return 0;
    real u = ((x[1]-a[1])*(c[2]-a[2]) - (c[1]-a[1])*(x[2]-a[2])) / d1;
    real v = ((b[1]-a[1])*(x[2]-a[2]) - (x[1]-a[1])*(b[2]-a[2])) / d1;
    if (u < 0.0 || v < 0.0 || u + v > 1.0) return 0;
    real xi = a[0] + u*(b[0]-a[0]) + v*(c[0]-a[0]);
    return (xi > x[0]);
}

int ft3_dentro(const ft3_superficie *s, const Point x)
{
    if (!s) return 0;
    real lo[3], hi[3]; _envoltoria(s, lo, hi);
    for (int d = 0; d < 3; d++) if (x[d] < lo[d] || x[d] > hi[d]) return 0;
    int cruz = 0;
    for (int t = 0; t < s->nt; t++) cruz += _cruza_raio(s, t, x);
    return (cruz & 1);
}

// Distancia do ponto ao triangulo t (algoritmo padrao de regiao: projeta no
// plano e, se cair fora, cai para a aresta ou o vertice mais proximo).
static real _dist_tri(const ft3_superficie *s, int t, const real p[3])
{
    const real *a = s->x[s->tri[t][0]];
    const real *b = s->x[s->tri[t][1]];
    const real *c = s->x[s->tri[t][2]];
    real ab[3], ac[3], ap[3];
    _sub(b, a, ab); _sub(c, a, ac); _sub(p, a, ap);
    real d1 = _ponto(ab, ap), d2 = _ponto(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return _norma(ap);

    real bp[3]; _sub(p, b, bp);
    real d3 = _ponto(ab, bp), d4 = _ponto(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return _norma(bp);

    real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        real w = d1 / (d1 - d3), q[3];
        for (int k = 0; k < 3; k++) q[k] = a[k] + w*ab[k] - p[k];
        return _norma(q);
    }
    real cp[3]; _sub(p, c, cp);
    real d5 = _ponto(ab, cp), d6 = _ponto(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return _norma(cp);

    real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        real w = d2 / (d2 - d6), q[3];
        for (int k = 0; k < 3; k++) q[k] = a[k] + w*ac[k] - p[k];
        return _norma(q);
    }
    real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
        real w = (d4-d3) / ((d4-d3) + (d5-d6)), q[3];
        for (int k = 0; k < 3; k++) q[k] = b[k] + w*(c[k]-b[k]) - p[k];
        return _norma(q);
    }
    // interior: distancia ao plano
    real den = va + vb + vc;
    real w1 = vb / den, w2 = vc / den, q[3];
    for (int k = 0; k < 3; k++) q[k] = a[k] + w1*ab[k] + w2*ac[k] - p[k];
    return _norma(q);
}

real ft3_distancia(const ft3_superficie *s, const Point x, real normal[DIM])
{
    if (!s || s->nt == 0) return 0.0;
    real melhor = 1e300; int tm = -1;
    for (int t = 0; t < s->nt; t++) {
        real d = _dist_tri(s, t, x);
        // DESEMPATE POR INDICE.  Distancia igual com triangulo diferente da'
        // NORMAL diferente, plano diferente e fracao diferente.  Sem esta regra
        // o caminho exaustivo e o acelerado divergem em celula rente a um plano
        // de simetria -- medido, 6 celulas em 8000 -- e a divergencia nao e' de
        // precisao, e' de qual empate venceu.
        if (d < melhor || (d == melhor && t < tm)) { melhor = d; tm = t; }
    }
    if (tm < 0) tm = 0;
    if (tm == 0x7fffffff) tm = 0;
    if (normal) {
        real n[3]; _normal2(s, tm, n);
        real m = _norma(n);
        for (int d = 0; d < 3; d++) normal[d] = (m > 0.0) ? n[d]/m : 0.0;
    }
    return melhor;
}

// VOLUME DA CAIXA DO LADO NEGATIVO DO PLANO n.x = alfa, com a caixa em
// [0,L0]x[0,L1]x[0,L2] e a origem no canto.  Formula fechada de
// Scardovelli & Zaleski: inclusao-exclusao sobre os cantos que o plano corta.
//
//   V = ( alfa^3 - SUM max(0, alfa - n_i L_i)^3
//                + SUM_{i<j} max(0, alfa - n_i L_i - n_j L_j)^3
//                - max(0, alfa - SUM n_i L_i)^3 ) / (6 n_0 n_1 n_2)
//
// Exige n_i >= 0; o chamador reflete os eixos de normal negativa, o que troca
// o canto de referencia mas nao o volume.
static real _volume_plano_caixa(const real n[3], real alfa, const real L[3])
{
    real den = 6.0 * n[0] * n[1] * n[2];
    real soma = alfa*alfa*alfa;
    real m[3];
    for (int i = 0; i < 3; i++) {
        m[i] = alfa - n[i]*L[i];
        if (m[i] > 0.0) soma -= m[i]*m[i]*m[i];
    }
    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++) {
            real w = alfa - n[i]*L[i] - n[j]*L[j];
            if (w > 0.0) soma += w*w*w;
        }
    real w = alfa - n[0]*L[0] - n[1]*L[1] - n[2]*L[2];
    if (w > 0.0) soma -= w*w*w;
    return soma / den;
}

// Definidas adiante, junto da grade; declaradas aqui porque `_fracao` as usa.
static real _dist_g_teto(const ft3_superficie *s, ft3_grade *g, const Point x,
                         real teto, real normal[3]);
static int  _ha_triangulo_perto(const ft3_grade *g, const real p[3], real raio);
static void _grade_envoltoria(const ft3_grade *g, real lo[3], real hi[3]);

// Corpo UNICO.  `g` nulo cai na busca exaustiva; nao nulo, na acelerada.  Duas
// copias desta funcao divergiriam com o tempo e o arreio que as compara
// deixaria de significar alguma coisa.
static real _fracao(const ft3_superficie *s, ft3_grade *g,
                    const real lo[DIM], const real hi[DIM])
{
    if (!s) return 0.0;
    real L[3], vol = 1.0, centro[3];
    for (int d = 0; d < 3; d++) {
        L[d] = hi[d] - lo[d];
        if (!(L[d] > 0.0)) return 0.0;
        vol *= L[d];
        centro[d] = 0.5 * (lo[d] + hi[d]);
    }
    const real meia_diag = 0.5 * sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);

    // Caminho rapido por envoltoria: caixa disjunta da envoltoria nao tem nada
    // dentro.  (O caso "envoltoria dentro da caixa" nao e' atalho: a caixa pode
    // conter a superficie inteira e a fracao nao ser 1.)
    //
    // COM GRADE, a envoltoria vem DELA.  Recalcula-la aqui varre os nv vertices
    // a cada celula -- medido: era o que segurava o ganho em 2,2x, porque a
    // grade acelera a busca por triangulo e nao a varredura de vertices.
    real blo[3], bhi[3];
    if (g) _grade_envoltoria(g, blo, bhi); else _envoltoria(s, blo, bhi);
    int disjunta = 0;
    for (int d = 0; d < 3; d++) if (hi[d] <= blo[d] || lo[d] >= bhi[d]) disjunta = 1;
    if (disjunta) return 0.0;

    real n[3];

    // TESTE DE OCUPACAO, antes de qualquer distancia.  Se nenhuma celula da
    // grade dentro da meia-diagonal guarda triangulo, a caixa esta' inteira de
    // um lado e a distancia exata nao e' usada -- basta dentro/fora.  Medido:
    // 94% das celulas de uma particao 20^3 caem aqui (7242 vazias e 304 cheias
    // de 8000), e sem este atalho a busca em cascas roda para todas elas.
    if (g && !_ha_triangulo_perto(g, centro, meia_diag))
        return ft3_dentro_g(s, g, centro) ? 1.0 : 0.0;

    // O teto e' a meia-diagonal: acima dela a caixa esta' inteira de um lado e o
    // valor exato nao e' usado.  Os dois caminhos tomam o MESMO ramo.
    real dist = g ? _dist_g_teto(s, g, centro, meia_diag, n)
                  : ft3_distancia(s, centro, n);
    const int dentro = g ? ft3_dentro_g(s, g, centro) : ft3_dentro(s, centro);

    // Longe da interface, a caixa inteira esta' de um lado so'.
    if (dist >= meia_diag) return dentro ? 1.0 : 0.0;

    // A normal do nucleo aponta para FORA.  A distancia COM SINAL do centro,
    // positiva fora, sai da paridade -- e nao do produto escalar com a normal
    // do triangulo mais proximo, que erra o lado perto de aresta.
    const real phi = dentro ? -dist : dist;

    // Plano n.(x - centro) = -phi  <=>  n.x = n.centro - phi.  O lado DENTRO e'
    // n.x < alfa.  Reflete os eixos de normal negativa para a formula fechada.
    real nn[3], alfa = -phi;
    for (int d = 0; d < 3; d++) {
        nn[d] = n[d];
        // coordenada local com origem no canto `lo`; o centro fica em L/2
        alfa += n[d] * 0.5 * L[d];
    }
    for (int d = 0; d < 3; d++)
        if (nn[d] < 0.0) { alfa -= nn[d] * L[d]; nn[d] = -nn[d]; }

    // Normal degenerada (alinhada a um eixo) zeraria o denominador: trata pelo
    // limite, que e' o corte por um plano perpendicular a um eixo.
    const real EPS = 1e-9;
    int nz = 0;
    for (int d = 0; d < 3; d++) if (nn[d] < EPS) nz++;
    if (nz > 0) {
        for (int d = 0; d < 3; d++) if (nn[d] < EPS) nn[d] = EPS;
        real m = sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2]);
        for (int d = 0; d < 3; d++) nn[d] /= m;
    }

    if (alfa <= 0.0) return 0.0;
    real total = nn[0]*L[0] + nn[1]*L[1] + nn[2]*L[2];
    if (alfa >= total) return 1.0;

    real f = _volume_plano_caixa(nn, alfa, L) / vol;
    return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
}

real ft3_distancia_g(const ft3_superficie *s, ft3_grade *g, const Point x,
                     real normal[DIM])
{
    if (!g) return ft3_distancia(s, x, normal);
    return _dist_g_teto(s, g, x, 1e300, normal);
}

real ft3_fracao_na_caixa(const ft3_superficie *s, const real lo[DIM], const real hi[DIM])
{ return _fracao(s, NULL, lo, hi); }

real ft3_fracao_na_caixa_g(const ft3_superficie *s, ft3_grade *g,
                           const real lo[DIM], const real hi[DIM])
{ return _fracao(s, g, lo, hi); }

real ft3_fracao_amostrada(const ft3_superficie *s, const real lo[DIM],
                          const real hi[DIM], int k)
{
    if (!s || k < 1) return 0.0;
    int dentro = 0;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            for (int l = 0; l < k; l++) {
                Point p;
                p[0] = lo[0] + (i + 0.5) * (hi[0]-lo[0]) / k;
                p[1] = lo[1] + (j + 0.5) * (hi[1]-lo[1]) / k;
                p[2] = lo[2] + (l + 0.5) * (hi[2]-lo[2]) / k;
                if (ft3_dentro(s, p)) dentro++;
            }
    return (real) dentro / (real) (k*k*k);
}


// --- grade espacial -------------------------------------------------------

struct ft3_grade {
    real lo[3], hi[3], h[3];   // extensao da GRADE (com folga)
    real blo[3], bhi[3];       // envoltoria EXATA da superficie, sem folga
    int  n[3];
    int *ini, *lst;        // CSR: triangulos da celula c em lst[ini[c]..ini[c+1])
    int *marca;            // carimbo por triangulo, contra contagem dupla
    int  carimbo;
};

static void _celula_de(const ft3_grade *g, const real p[3], int c[3])
{
    for (int d = 0; d < 3; d++) {
        int i = (int) floor((p[d] - g->lo[d]) / g->h[d]);
        if (i < 0) i = 0;
        if (i >= g->n[d]) i = g->n[d] - 1;
        c[d] = i;
    }
}

static int _lin(const ft3_grade *g, int i, int j, int k)
{ return (i * g->n[1] + j) * g->n[2] + k; }

ft3_grade *ft3_grade_cria(const ft3_superficie *s)
{
    if (!s || s->nt < 1) return NULL;
    ft3_grade *g = (ft3_grade *) calloc(1, sizeof *g);
    _envoltoria(s, g->lo, g->hi);
    // Guarda a envoltoria EXATA antes da folga: o teste de disjuncao tem de
    // usar a mesma caixa nos dois caminhos.  Com a folga, celula que encosta na
    // envoltoria era descartada pelo exaustivo e processada pelo acelerado, e as
    // duas fracoes diferiam por um residuo -- 12 celulas em 8000, medido.
    for (int d = 0; d < 3; d++) { g->blo[d] = g->lo[d]; g->bhi[d] = g->hi[d]; }
    // Folga: superficie plana em algum eixo daria celula de tamanho zero.
    for (int d = 0; d < 3; d++) {
        real ext = g->hi[d] - g->lo[d];
        real folga = (ext > 0.0) ? 1e-6 * ext : 1e-9;
        g->lo[d] -= folga; g->hi[d] += folga;
    }
    // Uma superficie ocupa um subconjunto BIDIMENSIONAL da grade, entao com n^3
    // celulas as ocupadas sao ~n^2.  Escolher n ~ sqrt(nt) da' da ordem de um
    // triangulo por celula ocupada -- que e' o ponto em que a busca deixa de ser
    // dominada pela lista e passa a ser dominada pela vizinhanca.
    // Celula do tamanho do LADO do triangulo.  Menor que isso, cada triangulo
    // se registra em dezenas de celulas (medido: 26,4) e a lista incha sem
    // ganho.  Maior, a celula devolve triangulos demais por consulta.
    real lado = s->ds_alvo;
    if (!(lado > 0.0)) lado = (g->hi[0] - g->lo[0]) / 16.0;
    real ext = g->hi[0] - g->lo[0];
    for (int d = 1; d < 3; d++) if (g->hi[d] - g->lo[d] > ext) ext = g->hi[d] - g->lo[d];
    int n = (int) (ext / lado + 0.5);
    if (n < 4)   n = 4;
    if (n > 128) n = 128;
    for (int d = 0; d < 3; d++) {
        g->n[d] = n;
        g->h[d] = (g->hi[d] - g->lo[d]) / n;
        if (!(g->h[d] > 0.0)) g->h[d] = 1e-12;
    }
    const int ncel = g->n[0] * g->n[1] * g->n[2];

    // Duas passadas: conta e preenche, para nao realocar por celula.
    int *cont = (int *) calloc((size_t) ncel + 1, sizeof *cont);
    for (int t = 0; t < s->nt; t++) {
        real blo[3], bhi[3];
        for (int d = 0; d < 3; d++) { blo[d] = 1e300; bhi[d] = -1e300; }
        for (int k = 0; k < 3; k++) {
            const real *v = s->x[s->tri[t][k]];
            for (int d = 0; d < 3; d++) {
                if (v[d] < blo[d]) blo[d] = v[d];
                if (v[d] > bhi[d]) bhi[d] = v[d];
            }
        }
        int a[3], b[3];
        _celula_de(g, blo, a); _celula_de(g, bhi, b);
        for (int i = a[0]; i <= b[0]; i++)
            for (int j = a[1]; j <= b[1]; j++)
                for (int k = a[2]; k <= b[2]; k++) cont[_lin(g,i,j,k) + 1]++;
    }
    for (int c = 0; c < ncel; c++) cont[c+1] += cont[c];
    const int nref = cont[ncel];
    g->ini = cont;
    g->lst = (int *) malloc((size_t) (nref > 0 ? nref : 1) * sizeof *g->lst);
    int *pos = (int *) malloc((size_t) ncel * sizeof *pos);
    memcpy(pos, cont, (size_t) ncel * sizeof *pos);
    for (int t = 0; t < s->nt; t++) {
        real blo[3], bhi[3];
        for (int d = 0; d < 3; d++) { blo[d] = 1e300; bhi[d] = -1e300; }
        for (int k = 0; k < 3; k++) {
            const real *v = s->x[s->tri[t][k]];
            for (int d = 0; d < 3; d++) {
                if (v[d] < blo[d]) blo[d] = v[d];
                if (v[d] > bhi[d]) bhi[d] = v[d];
            }
        }
        int a[3], b[3];
        _celula_de(g, blo, a); _celula_de(g, bhi, b);
        for (int i = a[0]; i <= b[0]; i++)
            for (int j = a[1]; j <= b[1]; j++)
                for (int k = a[2]; k <= b[2]; k++) g->lst[pos[_lin(g,i,j,k)]++] = t;
    }
    free(pos);
    g->marca = (int *) calloc((size_t) s->nt, sizeof *g->marca);
    g->carimbo = 0;
    return g;
}

void ft3_grade_destroi(ft3_grade *g)
{
    if (!g) return;
    free(g->ini); free(g->lst); free(g->marca); free(g);
}

void ft3_grade_estado(const ft3_grade *g, int n[DIM], int *refs)
{
    if (!g) return;
    for (int d = 0; d < 3; d++) n[d] = g->n[d];
    if (refs) *refs = g->ini[g->n[0]*g->n[1]*g->n[2]];
}

int ft3_dentro_g(const ft3_superficie *s, ft3_grade *g, const Point x)
{
    if (!s) return 0;
    if (!g)  return ft3_dentro(s, x);
    for (int d = 0; d < 3; d++) if (x[d] < g->lo[d] || x[d] > g->hi[d]) return 0;

    // O raio segue +x com (y,z) fixos: so' a COLUNA de celulas (i, j, k) com j e
    // k do ponto pode conter triangulo que ele cruza.  Um triangulo que cruza o
    // raio tem o ponto de cruzamento na propria envoltoria, logo esta' registrado
    // em alguma celula dessa coluna.
    int c[3]; _celula_de(g, x, c);
    g->carimbo++;
    int cruz = 0;
    for (int i = c[0]; i < g->n[0]; i++) {
        int cl = _lin(g, i, c[1], c[2]);
        for (int p = g->ini[cl]; p < g->ini[cl+1]; p++) {
            int t = g->lst[p];
            if (g->marca[t] == g->carimbo) continue;   // ja' contado nesta consulta
            g->marca[t] = g->carimbo;
            cruz += _cruza_raio(s, t, x);
        }
    }
    return (cruz & 1);
}

// Ha' algum triangulo registrado a menos de `raio` de `p`?  Varre APENAS a
// ocupacao das celulas, sem tocar em triangulo -- e' o atalho das celulas longe
// da interface.  Conservador: usa a caixa que contem a bola de raio `raio`,
// entao pode responder "sim" onde a resposta exata seria "nao".  Responder sim a
// mais so' custa trabalho; responder nao a mais mudaria o resultado, e isso nao
// pode acontecer.
static void _grade_envoltoria(const ft3_grade *g, real lo[3], real hi[3])
{
    for (int d = 0; d < 3; d++) { lo[d] = g->blo[d]; hi[d] = g->bhi[d]; }
}

static int _ha_triangulo_perto(const ft3_grade *g, const real p[3], real raio)
{
    int a[3], b[3];
    for (int d = 0; d < 3; d++) {
        a[d] = (int) floor((p[d] - raio - g->lo[d]) / g->h[d]);
        b[d] = (int) floor((p[d] + raio - g->lo[d]) / g->h[d]);
        if (a[d] < 0) a[d] = 0;
        if (b[d] >= g->n[d]) b[d] = g->n[d] - 1;
        if (a[d] > b[d]) return 0;               // caixa fora da grade
    }
    for (int i = a[0]; i <= b[0]; i++)
        for (int j = a[1]; j <= b[1]; j++)
            for (int k = a[2]; k <= b[2]; k++) {
                int cl = _lin(g, i, j, k);
                if (g->ini[cl+1] > g->ini[cl]) return 1;
            }
    return 0;
}

// `teto`: distancia acima da qual o chamador nao se importa com o valor exato.
// Sem ela, uma celula no CENTRO da gota expande cascas ate' o raio inteiro --
// (2r+1)^3 celulas -- e a grade fica MAIS LENTA que a varredura.  Medida antes
// do teto: ganho de 1,2x, ou seja nenhum.
static real _dist_g_teto(const ft3_superficie *s, ft3_grade *g, const Point x,
                         real teto, real normal[3])
{
    if (!s || s->nt == 0) return 0.0;

    int c[3]; _celula_de(g, x, c);
    real melhor = 1e300; int tm = 0x7fffffff;
    g->carimbo++;
    const int rmax = g->n[0] + g->n[1] + g->n[2];

    for (int r = 0; r <= rmax; r++) {
        int a[3], b[3];
        for (int d = 0; d < 3; d++) {
            a[d] = c[d] - r; if (a[d] < 0) a[d] = 0;
            b[d] = c[d] + r; if (b[d] >= g->n[d]) b[d] = g->n[d] - 1;
        }
        // So' a CASCA nova: o miolo ja' foi visto nas voltas anteriores.
        for (int i = a[0]; i <= b[0]; i++)
            for (int j = a[1]; j <= b[1]; j++)
                for (int k = a[2]; k <= b[2]; k++) {
                    if (r > 0 && i > c[0]-r && i < c[0]+r
                             && j > c[1]-r && j < c[1]+r
                             && k > c[2]-r && k < c[2]+r) continue;
                    int cl = _lin(g, i, j, k);
                    for (int p = g->ini[cl]; p < g->ini[cl+1]; p++) {
                        int t = g->lst[p];
                        if (g->marca[t] == g->carimbo) continue;
                        g->marca[t] = g->carimbo;
                        real d = _dist_tri(s, t, x);
                        if (d < melhor || (d == melhor && t < tm)) { melhor = d; tm = t; }
                    }
                }

        // PARADA SEGURA: distancia de `x` ate' sair do bloco ja' varrido.  Nada
        // fora dele pode estar mais perto que isso, entao achar algo melhor e'
        // impossivel e a busca termina.  Calculada da extensao FISICA do bloco,
        // e nao de r*h -- o ponto pode estar rente a' borda da propria celula.
        real fora = 1e300;
        int cobre_tudo = 1;
        for (int d = 0; d < 3; d++) {
            real blo = g->lo[d] + a[d] * g->h[d];
            real bhi = g->lo[d] + (b[d] + 1) * g->h[d];
            if (a[d] > 0)            { real q = x[d] - blo; if (q < fora) fora = q; cobre_tudo = 0; }
            if (b[d] < g->n[d] - 1)  { real q = bhi - x[d]; if (q < fora) fora = q; cobre_tudo = 0; }
        }
        if (cobre_tudo) break;
        if (melhor <= fora) break;
        // Nada mais perto que `fora` existe fora do bloco; se `fora` ja' passou
        // do teto, o chamador nao precisa do valor exato e a busca para.
        if (fora >= teto) { if (melhor > teto) melhor = fora; break; }
    }

    // Parar pela cota SEM ter achado triangulo e' o caso comum -- celula longe
    // da interface --, e ai' `tm` continua invalido.  Usa-lo indexaria fora do
    // vetor de triangulos.  O chamador nao le' a normal nesse ramo, mas calcular
    // uma normal invalida segmenta antes de ele decidir nao usa-la.
    if (tm == 0x7fffffff) {
        if (normal) normal[0] = normal[1] = normal[2] = 0.0;
        return melhor;
    }
    if (normal) {
        real n[3]; _normal2(s, tm, n);
        real m = _norma(n);
        for (int d = 0; d < 3; d++) normal[d] = (m > 0.0) ? n[d]/m : 0.0;
    }
    return melhor;
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
