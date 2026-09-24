// Ver o .h para o contrato e para o par adjunto, que e' o que importa aqui.

#include "hig-flow-fronteira-imersa.h"

#include <petsc.h>
#include <petscdmswarm.h>
#include <petscsf.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct fi_corpo {
    DM       enxame;      // DMSwarm: posicao, peso, velocidade, forca, h
    real     ds;          // espacamento ALVO entre marcadores (geometrico)
    //
    // O h do NUCLEO nao e' global: e' o tamanho da celula onde cada marcador
    // esta', guardado no campo "h" do enxame.  Em malha uniforme todos sao
    // iguais e a distincao nao aparece; numa malha graduada ela e' a diferenca
    // entre funcionar e nao funcionar.
    //
    // `ds` e o h do nucleo SAO COISAS DIFERENTES que em malha uniforme calham
    // de ser iguais: um e' espacamento de marcador (geometrico, escolhido), o
    // outro e' tamanho de celula (da malha, lido).
    int      codim;       // codimensao do corpo -- ver _volume_marcador
    real     forca_passo[DIM];  // forca do PASSO, somada sobre as iteracoes
    MPI_Comm comm;
};

// O VOLUME DE UM MARCADOR e' `peso * h^codim`, e o expoente e' a CODIMENSAO do
// corpo, nao a dimensao do espaco.
//
//   curva em 2D      peso = comprimento, codim 1  ->  dV = ds * h
//   superficie em 3D peso = area,        codim 1  ->  dV = dA * h
//   curva em 3D      peso = comprimento, codim 2  ->  dV = ds * h^2
//
// Em 2D codimensao e DIM-1 coincidem, e por isso `h^(DIM-1)` passou nos testes.
// Em 3D daria h^2 para uma superficie: forca 1/h vezes pequena demais, e o
// sintoma seria um corpo POROSO -- o fluido atravessando devagar -- que se le'
// como "malha grosseira" em vez de como erro.
static real _volume_marcador(const fi_corpo *c, real peso, real h_mar)
{
    real v = peso;
    for (int i = 0; i < c->codim; i++) v *= h_mar;
    return v;
}

//! O tamanho da celula que contem `x`, na direcao 0.  E' o h do nucleo para um
//! marcador ali.  Devolve 0 se nao houver celula.
static real _h_da_celula(sim_facet_domain *sfd, const Point x)
{
    hig_cell *cel = sfd_get_cell_with_point(sfd, (real *) x);
    if (cel == NULL) {
        // Ponto sobre fronteira de celula: cutuca.  Mesmo motivo do _suporte.
        const real eps = 1e-9;
        const int combos = 1 << DIM;
        for (int m = 0; m < combos && cel == NULL; m++) {
            Point q;
            for (int d = 0; d < DIM; d++)
                q[d] = x[d] + (((m >> d) & 1) ? eps : -eps);
            cel = sfd_get_cell_with_point(sfd, q);
        }
    }
    if (cel == NULL) return 0.0;
    Point delta;
    hig_get_delta(cel, delta);
    return delta[0];
}

// ------------------------------------------------------------------ o nucleo

real fi_delta_roma(real r)
{
    const real a = fabs(r);
    if (a <= 0.5)
        return (1.0 + sqrt(1.0 - 3.0*r*r)) / 3.0;
    if (a <= 1.5) {
        const real t = 1.0 - a;
        return (5.0 - 3.0*a - sqrt(1.0 - 3.0*t*t)) / 6.0;
    }
    return 0.0;
}

// ------------------------------------------------------- posse de um ponto

// Verdadeiro se a celula que contem `x` pertence a uma arvore PROPRIA deste
// rank.  A copia de franja da mesma celula tem raiz numa arvore de franja, cujo
// indice e' >= o numero de arvores locais -- entao so' um rank responde sim.
static int _possui(sim_facet_domain *sfd, const Point x)
{
    hig_cell *c = sfd_get_cell_with_point(sfd, (real *) x);
    if (c == NULL) return 0;

    // CONVENCAO MEIO-ABERTA [lo, hi): marcador exatamente sobre a face ALTA da
    // celula pertence a' vizinha.  Sem isto, um marcador sobre a face entre
    // arvores de ranks diferentes e' achado como "meu" pelos DOIS lados -- a
    // busca de cada rank e' inclusiva na borda da sua arvore -- e a guarda de
    // reivindicacao aborta.  MEDIDO com a malha do criterio (F3): 252
    // reivindicacoes para 251 marcadores.  A malha da caixa nunca poe fronteira
    // de particao sobre marcador; a do criterio pos.
    //
    // A re-consulta com o ponto empurrado escolhe a vizinha; se nao houver
    // (borda externa do dominio), a celula original fica -- ali nao ha' outro
    // rank para disputar.
    {
        Point cl, ch, q;
        hig_get_lowpoint(c, cl);
        hig_get_highpoint(c, ch);
        POINT_ASSIGN(q, x);
        int na_face = 0;
        for (int d = 0; d < DIM; d++)
            if (fabs(x[d] - ch[d]) <= 1e-9 * (ch[d] - cl[d])) {
                q[d] = ch[d] + 1e-6 * (ch[d] - cl[d]);
                na_face = 1;
            }
        if (na_face) {
            hig_cell *v = sfd_get_cell_with_point(sfd, q);
            if (v != NULL) c = v;
        }
    }

    hig_cell *raiz = c;
    while (hig_get_parent(raiz) != NULL) raiz = hig_get_parent(raiz);

    sim_domain *sd = sfd->cdom;
    const unsigned n = sd_get_num_local_higtrees(sd);
    for (unsigned i = 0; i < n; i++)
        if (sd_get_higtree(sd, i) == raiz) return 1;
    return 0;
}

// ------------------------------------------------------------ construcao

// Comum aos geradores: todo rank produz a MESMA lista de candidatos -- a
// geometria e' analitica e pequena, entao isso evita comunicacao -- e cada um
// fica com os que possui.
static fi_corpo *_de_candidatos(sim_facet_domain *sfd, real h, int codim,
                                const Point *pts, const real *pesos, int ncand)
{
    // O PETSc e' inicializado PREGUICOSAMENTE pela HiGTree, ao criar um solver
    // (utils.c).  Este modulo usa PETSc sem solver nenhum, entao a inicializacao
    // pode nao ter acontecido -- e o sintoma e' um erro em MPI_Comm_get_attr,
    // que nao menciona inicializacao.  Declarado no sitio como em
    // higtree/src/solver-petsc.c:362, que e' a convencao da casa.
    void _try_initialize_petsc(void);
    _try_initialize_petsc();

    fi_corpo *c = (fi_corpo *) malloc(sizeof *c);
    c->ds    = h;          // o argumento e' o ESPACAMENTO alvo, nao o h do nucleo
    c->codim = codim;
    c->comm  = MPI_COMM_WORLD;
    for (int d = 0; d < DIM; d++) c->forca_passo[d] = 0.0;

    PetscCallAbort(c->comm, DMCreate(c->comm, &c->enxame));
    PetscCallAbort(c->comm, DMSetType(c->enxame, DMSWARM));
    PetscCallAbort(c->comm, DMSetDimension(c->enxame, DIM));
    PetscCallAbort(c->comm, DMSwarmSetType(c->enxame, DMSWARM_BASIC));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "posicao",    DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "peso",         1, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "h",            1, PETSC_REAL));
    // INDICE AO LONGO DA CURVA.  A distribuicao por posse euleriana destroi a
    // ordem -- cada rank fica com um subconjunto arbitrario --, entao sem este
    // campo nao ha' como reconstruir a curva, nem para desenhar nem para
    // calcular qualquer coisa que dependa de vizinhanca.
    //
    // E' tambem o que o caso de INTERFACE entre fluidos vai exigir: curvatura
    // precisa de vizinhos.  Guardar agora custa um campo.
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "indice",       1, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "velocidade", DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "forca",      DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmFinalizeFieldRegister(c->enxame));

    int meus = 0;
    for (int i = 0; i < ncand; i++) if (_possui(sfd, pts[i])) meus++;

    PetscCallAbort(c->comm, DMSwarmSetLocalSizes(c->enxame, meus, 4));
    PetscReal *pos = NULL, *peso = NULL, *hmar = NULL, *idx = NULL;
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "peso",    NULL, NULL, (void **) &peso));
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "h",       NULL, NULL, (void **) &hmar));
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "indice",  NULL, NULL, (void **) &idx));
    int k = 0;
    for (int i = 0; i < ncand; i++) {
        if (!_possui(sfd, pts[i])) continue;
        for (int d = 0; d < DIM; d++) pos[DIM*k + d] = pts[i][d];
        peso[k] = pesos[i];
        // O h DO NUCLEO vem da malha, nao do argumento.  E' isto que permite
        // malha graduada: cada marcador usa o tamanho da celula onde esta'.
        idx[k]  = (PetscReal) i;     // posicao GLOBAL na curva, antes de distribuir
        hmar[k] = _h_da_celula(sfd, pts[i]);
        if (!(hmar[k] > 0.0)) {
            fprintf(stderr, "fronteira imersa: marcador em (%g,%g) sem celula\n",
                    (double) pts[i][0], (double) pts[i][1]);
            abort();
        }
        k++;
    }
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "peso",    NULL, NULL, (void **) &peso));
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "h",       NULL, NULL, (void **) &hmar));
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "indice",  NULL, NULL, (void **) &idx));

    // Cada marcador tem de ser reivindicado por EXATAMENTE um rank.  Marcador
    // sobre a face entre celulas de ranks diferentes e' o caso que quebra isso,
    // e seguir daria forca errada sem sintoma visivel.
    int soma = 0;
    MPI_Allreduce(&meus, &soma, 1, MPI_INT, MPI_SUM, c->comm);
    if (soma != ncand) {
        int rank; MPI_Comm_rank(c->comm, &rank);
        if (rank == 0)
            fprintf(stderr,
                "fronteira imersa: %d marcadores no total, mas os ranks reivindicaram %d.  "
                "Marcador sobre fronteira de particao reivindicado por dois ranks (ou por "
                "nenhum).  Seguir daria forca errada sem sintoma visivel.\n",
                ncand, soma);
        MPI_Barrier(c->comm);
        abort();
    }
    return c;
}

// Percorre a curva fechada e emite (ponto, peso) nos CENTROS dos subsegmentos.
// Devolve quantos; com `pts`/`pesos` nulos so' conta.
static int _percorre_curva(const Point *v, int nvert, real h,
                           Point *pts, real *pesos)
{
    int n = 0;
    for (int i = 0; i < nvert; i++) {
        const real *a = v[i];
        const real *b = v[(i + 1) % nvert];        // fecha sozinha
        real comp = 0.0;
        for (int d = 0; d < DIM; d++) comp += (b[d] - a[d]) * (b[d] - a[d]);
        comp = sqrt(comp);
        if (comp <= 0.0) continue;
        int nsub = (int) (comp / h + 0.5);
        if (nsub < 1) nsub = 1;
        const real ds = comp / nsub;
        for (int s = 0; s < nsub; s++) {
            if (pts != NULL) {
                const real t = (s + 0.5) / nsub;
                for (int d = 0; d < DIM; d++) pts[n][d] = a[d] + t * (b[d] - a[d]);
                pesos[n] = ds;
            }
            n++;
        }
    }
    return n;
}

fi_corpo *fi_cria_curva(sim_facet_domain *sfd, const Point *vertices, int nvert,
                        real h)
{
    if (nvert < 3 || h <= 0.0) {
        fprintf(stderr, "fi_cria_curva: curva com %d vertices e h=%g\n", nvert, h);
        abort();
    }
    const int n = _percorre_curva(vertices, nvert, h, NULL, NULL);
    Point *pts  = (Point *) malloc(n * sizeof *pts);
    real  *peso = (real  *) malloc(n * sizeof *peso);
    _percorre_curva(vertices, nvert, h, pts, peso);
    fi_corpo *c = _de_candidatos(sfd, h, 1, pts, peso, n);   // curva em 2D: codim 1
    free(pts); free(peso);
    return c;
}

#if DIM == 3
fi_corpo *fi_cria_extrusao(sim_facet_domain *sfd, const Point *vertices, int nvert,
                           real z0, real z1, real h)
{
    if (nvert < 3 || h <= 0.0 || !(z1 > z0)) {
        fprintf(stderr, "fi_cria_extrusao: %d vertices, h=%g, z de %g a %g\n",
                nvert, h, z0, z1);
        abort();
    }
    const int nc = _percorre_curva(vertices, nvert, h, NULL, NULL);
    Point *curva = (Point *) malloc(nc * sizeof *curva);
    real  *dscur = (real  *) malloc(nc * sizeof *dscur);
    _percorre_curva(vertices, nvert, h, curva, dscur);

    // Em z, os marcadores ficam nos CENTROS dos intervalos, pelo mesmo motivo
    // que na curva: centro nao e' compartilhado entre intervalos vizinhos, e a
    // soma dos pesos da a area exata.
    int nz = (int) ((z1 - z0) / h + 0.5);
    if (nz < 1) nz = 1;
    const real dz = (z1 - z0) / nz;

    const int n = nc * nz;
    Point *pts  = (Point *) malloc(n * sizeof *pts);
    real  *peso = (real  *) malloc(n * sizeof *peso);
    int k = 0;
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nz; j++) {
            pts[k][0] = curva[i][0];
            pts[k][1] = curva[i][1];
            pts[k][2] = z0 + (j + 0.5) * dz;
            peso[k]   = dscur[i] * dz;          // AREA do retalho
            k++;
        }

    // Superficie em 3D: codimensao 1, igual a' curva em 2D.  NAO e' DIM-1.
    fi_corpo *c = _de_candidatos(sfd, h, 1, pts, peso, n);
    free(curva); free(dscur); free(pts); free(peso);
    return c;
}

fi_corpo *fi_cria_cilindro(sim_facet_domain *sfd, real cx, real cy, real raio,
                           int nlados, real z0, real z1, real h)
{
    Point *v = (Point *) malloc(nlados * sizeof *v);
    for (int i = 0; i < nlados; i++) {
        const real a = 2.0 * M_PI * i / nlados;
        v[i][0] = cx + raio * cos(a);
        v[i][1] = cy + raio * sin(a);
        v[i][2] = z0;
    }
    fi_corpo *c = fi_cria_extrusao(sfd, (const Point *) v, nlados, z0, z1, h);
    free(v);
    return c;
}
#endif

fi_corpo *fi_cria_circulo(sim_facet_domain *sfd, const Point centro, real raio,
                          int nlados, real h)
{
    Point *v = (Point *) malloc(nlados * sizeof *v);
    for (int i = 0; i < nlados; i++) {
        const real a = 2.0 * M_PI * i / nlados;
        v[i][0] = centro[0] + raio * cos(a);
        v[i][1] = centro[1] + raio * sin(a);
        for (int d = 2; d < DIM; d++) v[i][d] = centro[d];
    }
    fi_corpo *c = fi_cria_curva(sfd, (const Point *) v, nlados, h);
    free(v);
    return c;
}

void fi_destroi(fi_corpo *c)
{
    if (c == NULL) return;
    DMDestroy(&c->enxame);
    free(c);
}

int fi_num_locais(const fi_corpo *c)
{
    PetscInt n = 0;
    DMSwarmGetLocalSize(c->enxame, &n);
    return (int) n;
}

real fi_peso_total(const fi_corpo *c)
{
    PetscInt n = 0;
    PetscReal *peso = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "peso", NULL, NULL, (void **) &peso);
    real soma = 0.0;
    for (PetscInt i = 0; i < n; i++) soma += peso[i];
    DMSwarmRestoreField(c->enxame, "peso", NULL, NULL, (void **) &peso);
    real total = 0.0;
    MPI_Allreduce(&soma, &total, 1, MPI_DOUBLE, MPI_SUM, c->comm);
    return total;
}

void fi_escreve_vtk(const fi_corpo *c, const char *prefixo, int quadro)
{
    int rank;
    MPI_Comm_rank(c->comm, &rank);

    // UM ARQUIVO POR RANK, como o escritor euleriano ja' faz -- os marcadores
    // sao distribuidos por posse, entao nao ha' arquivo global sem comunicacao.
    char nome[512];
    snprintf(nome, sizeof nome, "%s_lag_%d-%d.vtk", prefixo, rank, quadro);
    FILE *fp = fopen(nome, "w");
    if (fp == NULL) {
        fprintf(stderr, "fi_escreve_vtk: nao abriu %s\n", nome);
        return;
    }

    PetscInt n = 0;
    PetscReal *pos = NULL, *vel = NULL, *f = NULL, *peso = NULL, *hmar = NULL, *idx = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmGetField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",       NULL, NULL, (void **) &peso);
    DMSwarmGetField(c->enxame, "h",          NULL, NULL, (void **) &hmar);
    DMSwarmGetField(c->enxame, "indice",     NULL, NULL, (void **) &idx);

    fprintf(fp, "# vtk DataFile Version 3.0\n"
                "malha lagrangeana da fronteira imersa\n"
                "ASCII\nDATASET POLYDATA\nPOINTS %ld float\n", (long) n);
    for (PetscInt k = 0; k < n; k++) {
        // VTK quer sempre tres coordenadas, mesmo em 2D.
        fprintf(fp, "%g %g %g\n", (double) pos[DIM*k],
                (double) pos[DIM*k + 1],
                (double) (DIM > 2 ? pos[DIM*k + 2] : 0.0));
    }
    // VERTICES: sem isto o ParaView abre o arquivo e nao mostra nada.
    fprintf(fp, "VERTICES %ld %ld\n", (long) n, (long) (2*n));
    for (PetscInt k = 0; k < n; k++) fprintf(fp, "1 %ld\n", (long) k);

    fprintf(fp, "POINT_DATA %ld\n", (long) n);
    fprintf(fp, "VECTORS forca float\n");
    for (PetscInt k = 0; k < n; k++)
        fprintf(fp, "%g %g %g\n", (double) f[DIM*k], (double) f[DIM*k + 1],
                (double) (DIM > 2 ? f[DIM*k + 2] : 0.0));
    fprintf(fp, "VECTORS velocidade float\n");
    for (PetscInt k = 0; k < n; k++)
        fprintf(fp, "%g %g %g\n", (double) vel[DIM*k], (double) vel[DIM*k + 1],
                (double) (DIM > 2 ? vel[DIM*k + 2] : 0.0));
    // O peso e' a medida geometrica (comprimento em 2D, area em 3D), nao o
    // volume -- ver `fi_forca_total` para a diferenca, que ja' custou caro.
    fprintf(fp, "SCALARS peso float 1\nLOOKUP_TABLE default\n");
    for (PetscInt k = 0; k < n; k++) fprintf(fp, "%g\n", (double) peso[k]);
    // O h DA CELULA onde o marcador esta'.  Num corpo sobre malha graduada ele
    // varia, e ver isso e' o jeito mais direto de conferir que o corpo ficou no
    // nivel que se pretendia.
    fprintf(fp, "SCALARS h_celula float 1\nLOOKUP_TABLE default\n");
    for (PetscInt k = 0; k < n; k++) fprintf(fp, "%g\n", (double) hmar[k]);
    // O INDICE ao longo da curva.  Sem ele a malha lagrangeana e' uma nuvem de
    // pontos: a posse euleriana embaralha a ordem entre os ranks.  No ParaView,
    // colorir ou ordenar por este campo reconstitui a curva.
    fprintf(fp, "SCALARS indice float 1\nLOOKUP_TABLE default\n");
    for (PetscInt k = 0; k < n; k++) fprintf(fp, "%g\n", (double) idx[k]);

    DMSwarmRestoreField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmRestoreField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",       NULL, NULL, (void **) &peso);
    DMSwarmRestoreField(c->enxame, "h",          NULL, NULL, (void **) &hmar);
    DMSwarmRestoreField(c->enxame, "indice",     NULL, NULL, (void **) &idx);
    fclose(fp);
}

void fi_zera_forca_passo(fi_corpo *c)
{
    for (int d = 0; d < DIM; d++) c->forca_passo[d] = 0.0;
}

void fi_acumula_forca_passo(fi_corpo *c)
{
    real parcial[DIM];
    fi_forca_total(c, parcial);
    for (int d = 0; d < DIM; d++) c->forca_passo[d] += parcial[d];
}

void fi_forca_passo(const fi_corpo *c, real forca[DIM])
{
    for (int d = 0; d < DIM; d++) forca[d] = c->forca_passo[d];
}

real fi_residuo_max(const fi_corpo *c)
{
    PetscInt n = 0;
    PetscReal *vel = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    // NaN TEM DE SOBREVIVER AO MAXIMO.  `m > pior` com NaN e' FALSO, entao um
    // maximo ingenuo devolve ZERO num campo que explodiu -- foi o que aconteceu
    // na primeira corrida do Uhlmann: residuo 0,000000e+00 com max|u| em 1e96.
    // Maximo que nao enxerga NaN e' maximo que mente.
    real pior = 0.0;
    for (PetscInt k = 0; k < n; k++) {
        real m = 0.0;
        for (int d = 0; d < DIM; d++) m += vel[DIM*k + d] * vel[DIM*k + d];
        m = sqrt(m);
        if (!isfinite(m)) { pior = m; break; }
        if (m > pior) pior = m;
    }
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    real global = 0.0;
    MPI_Allreduce(&pior, &global, 1, MPI_DOUBLE, MPI_MAX, c->comm);
    return global;
}

void fi_forca_total(const fi_corpo *c, real total[DIM])
{
    PetscInt n = 0;
    PetscReal *f = NULL, *peso = NULL, *hmar = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "forca", NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",  NULL, NULL, (void **) &peso);
    DMSwarmGetField(c->enxame, "h",     NULL, NULL, (void **) &hmar);
    real local[DIM];
    for (int d = 0; d < DIM; d++) local[d] = 0.0;
    for (PetscInt i = 0; i < n; i++)
        for (int d = 0; d < DIM; d++)
            local[d] += f[DIM*i + d] * _volume_marcador(c, peso[i], hmar[i]);
    DMSwarmRestoreField(c->enxame, "forca", NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",  NULL, NULL, (void **) &peso);
    DMSwarmRestoreField(c->enxame, "h",     NULL, NULL, (void **) &hmar);
    MPI_Allreduce(local, total, DIM, MPI_DOUBLE, MPI_SUM, c->comm);
}

// ------------------------------------------------- o suporte de um marcador

// O suporte do nucleo de Roma e' 1,5 celulas para cada lado, entao 5 posicoes
// por direcao bastam com folga.  Nao se varre faceta nenhuma: a malha e'
// uniforme, os centros do suporte estao em DESLOCAMENTOS CONHECIDOS, e cada um
// se acha com uma localizacao de ponto.  Varrer seria O(marcadores x facetas).
#define FI_LARGURA 5
#define FI_MEIO    2

// Quantos pontos de suporte foram descartados por id inutilizavel.  Nao e'
// contador de depuracao: descartar reduz a soma do nucleo abaixo de 1 naquele
// marcador, entao a forca sai menor ali -- e o sintoma seria corpo levemente
// poroso, que se le' como "malha grosseira".  Se este numero nao for zero, ha'
// o que investigar.
static long _suporte_perdidos = 0;
// Separados porque significam coisas DIFERENTES: id espelhado e' convencao e
// tem conserto; faceta fora do mapeador e' suporte que sai do dominio mapeado,
// e conserto seria outro.  Contar junto foi o que me fez "corrigir" tres vezes
// o ramo errado -- o numero repetia identico e eu lia como correcao que falhou.
static long _perdidos_espelho = 0;   // fid < 0
static long _perdidos_mapa    = 0;   // fid >= 0 mas nao esta' no mapeador
static long _perdidos_sem_faceta = 0;// nao ha' faceta naquele ponto
// Pontos de suporte que cairam numa celula de tamanho DIFERENTE do h do
// marcador -- isto e', o suporte atravessou uma fronteira de refinamento.
//
// DEVE SER ZERO.  Atravessando, o nucleo deixa de ser normalizado: a particao
// da unidade vale para um h so'.  O tratamento correto e' o do
// Roma-Peskin-Berger (a versao ADAPTATIVA do metodo), que nao esta'
// implementado -- entao a restricao e' que o corpo fique inteiramente dentro de
// um nivel, com folga maior que o suporte (1,5 celulas).
//
// Contar em vez de abortar porque o numero diz QUANTO do corpo esta' fora da
// regiao fina, o que orienta o ajuste da caixa de refino.
static long _suporte_nivel_trocado = 0;
// Quantas facetas de suporte a ULTIMA interpolacao achou, somadas sobre
// marcadores e direcoes.  Com `fi_num_locais`, separa duas causas que dao o
// mesmo sintoma (forca zero): marcador que sumiu do rank, e marcador que esta'
// la' mas nao acha faceta.
static long _suporte_achados = 0;
static int  _diag_suporte = 0;   // ligado por FI_DIAG_SUPORTE, um marcador so'

long fi_suporte_achados(void) { return _suporte_achados; }

//! Maior |u| e maior |f| entre os marcadores DESTE rank, SEM reducao.
//! Existe para separar "o campo e' zero" de "a reducao esta' quebrada": os dois
//! dao o mesmo sintoma no valor global.
void fi_locais_cru(const fi_corpo *c, real *umax, real *fmax)
{
    PetscInt n = 0;
    PetscReal *vel = NULL, *f = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmGetField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    real um = 0.0, fm = 0.0;
    for (PetscInt k = 0; k < n; k++) {
        for (int d = 0; d < DIM; d++) {
            const real a = fabs(vel[DIM*k + d]); if (a > um) um = a;
            const real b = fabs(f[DIM*k + d]);   if (b > fm) fm = b;
        }
    }
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmRestoreField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    *umax = um; *fmax = fm;
}

long fi_suporte_nivel_trocado(void) { return _suporte_nivel_trocado; }

long fi_perdidos_espelho(void)     { return _perdidos_espelho; }
long fi_perdidos_mapa(void)        { return _perdidos_mapa; }
long fi_perdidos_sem_faceta(void)  { return _perdidos_sem_faceta; }

long fi_suporte_perdidos(void) { return _suporte_perdidos; }

// Preenche `lids` e `pesos` com as facetas do suporte e o peso do nucleo.
// Devolve quantas.  `capac` deve ser >= FI_LARGURA^DIM.
static int _suporte(sim_facet_domain *sfd, int dim, const Point X, real h,
                    int *lids, real *pesos, int capac)
{
    // ANCORA, e ela tem de sair da CELULA, nao de uma faceta.
    //
    // `sfd_get_facet_with_point` acha a faceta cujo PLANO contem o ponto.  Um
    // marcador no meio de uma celula nao esta' sobre plano nenhum na direcao
    // `dim`, e a busca falha -- nao por degenerescencia, mas porque nao existe
    // tal faceta.  Ancorar numa faceta so' funcionava para a direcao em que o
    // marcador por acaso caisse sobre um plano de face.
    //
    // O sintoma foi forca EXATAMENTE pela metade numa curva alinhada com a
    // malha, e a conservacao NAO o pegou: os dois lados dela usam os mesmos
    // marcadores.  Quem pegou foi a particao da unidade, que compara com valor
    // absoluto conhecido.
    hig_cell *cel = sfd_get_cell_with_point(sfd, (real *) X);
    if (cel == NULL) {
        // Ponto exatamente sobre fronteira de celula: cutuca e tenta de novo.
        const real eps = 1e-6 * h;
        const int combos = 1 << DIM;
        for (int m = 0; m < combos && cel == NULL; m++) {
            Point q;
            for (int d = 0; d < DIM; d++)
                q[d] = X[d] + (((m >> d) & 1) ? eps : -eps);
            cel = sfd_get_cell_with_point(sfd, q);
        }
    }
    if (cel == NULL) return 0;

    // Centro da face BAIXA da celula na direcao `dim`: e' um centro de faceta
    // legitimo, e serve de origem para a grade de deslocamentos.
    Point cbaixo, ccentro, c0;
    hig_get_lowpoint(cel, cbaixo);
    hig_get_center(cel, ccentro);
    for (int d = 0; d < DIM; d++) c0[d] = (d == dim) ? cbaixo[d] : ccentro[d];

    int total = FI_LARGURA;
    for (int d = 1; d < DIM; d++) total *= FI_LARGURA;

    int n = 0;
    for (int idx = 0; idx < total; idx++) {
        Point p;
        int resto = idx;
        for (int d = 0; d < DIM; d++) {
            const int passo = (resto % FI_LARGURA) - FI_MEIO;
            resto /= FI_LARGURA;
            p[d] = c0[d] + passo * h;
        }

        hig_facet f;
        if (!sfd_get_facet_with_point(sfd, p, &f)) { _perdidos_sem_faceta++; continue; }
        Point cf;
        hig_get_facet_center(&f, cf);

        // O peso sai do centro ENCONTRADO, nao do ponto tentado: se a
        // localizacao devolver faceta com centro ligeiramente outro, o peso
        // continua certo.
        real w = 1.0;
        for (int d = 0; d < DIM; d++) w *= fi_delta_roma((cf[d] - X[d]) / h);
        if (w == 0.0) continue;

        if (n >= capac) {
            fprintf(stderr, "_suporte: mais de %d facetas no suporte\n", capac);
            abort();
        }
        // O SUPORTE ATRAVESSOU UM NIVEL DE REFINAMENTO?
        //
        // O nucleo so' e' normalizado para um h.  Se a celula desta faceta tem
        // outro tamanho, a soma dos pesos deixa de dar 1 e a forca sai errada
        // naquele marcador -- silenciosamente.
        {
            hig_cell *cf_cel = sfd_get_cell_with_point(sfd, cf);
            if (cf_cel != NULL) {
                Point dcel;
                hig_get_delta(cf_cel, dcel);
                if (dcel[0] < 0.9*h || dcel[0] > 1.1*h) _suporte_nivel_trocado++;
            }
        }

        // O ID PODE NAO SER USAVEL, e isso nao e' excepcional.
        //
        // `sfd_adjust_facet_ids` (domain.h) marca a faceta compartilhada entre
        // duas celulas guardando `-id` na de coordenada maior, e `sfd_get_local_id`
        // devolve -1 para ela (domain.c:3043).  Usar isso como indice e' escrita
        // em dp[-1]: corrupcao de heap, e o abort sai longe dali.
        //
        // O teste de `fronteira-imersa/testes` NAO pega este caso -- o dominio
        // dele nao passa por ajuste de ids.  Quem pegou foi o exemplo real.
        int lid = sfd_get_local_id(sfd, &f);
        if (lid < 0) {
            // O ID NEGATIVO NAO E' LIXO: E' A CONVENCAO.
            //
            // `sfd_adjust_facet_ids` (domain.h:464) diz que, numa faceta
            // compartilhada, a celula de coordenada MENOR guarda o id e a outra
            // guarda `-id`.  `sfd_get_local_id` devolve -1 para essa (domain.c),
            // mas o id canonico e' simplesmente o simetrico.
            //
            // A tentativa anterior -- deslocar o ponto e reconsultar -- nao
            // podia funcionar: o plano da faceta e' perpendicular a `dim`, entao
            // andar ao longo de `dim` tira o ponto do plano e nenhuma faceta o
            // contem.  O sintoma foi o contador repetir o MESMO numero,
            // 46080, que e' sinal de caminho novo sem efeito, nao de correcao
            // que errou por pouco.
            // MEDIDO: o ramo do espelho NUNCA e' tomado neste caso (contador
            // espelho = 0 em corrida completa).  Ele fica porque a convencao de
            // `sfd_adjust_facet_ids` existe e pode ocorrer noutra geometria --
            // mas nao era a causa das perdas, e eu "corrigi" esse ramo duas
            // vezes antes de separar os contadores.  Numero que repete
            // IDENTICO e' ramo nao tomado, nao correcao que errou por pouco.
            const uniqueid fid = hig_get_fid(&f);
            if (fid < 0) {
                _perdidos_espelho++;
                const mp_value_t v = mp_lookup(sfd_get_domain_mapper(sfd), -fid);
                if (v != MP_UNDEF) lid = (int) v;
            } else {
                // A faceta existe na geometria e NAO esta' no mapeador deste
                // dominio -- suporte saindo da regiao mapeada.
                _perdidos_mapa++;
            }
        }
        if (lid < 0) { _suporte_perdidos++; continue; }
        // COMPARACAO PEDIDA x ENCONTRADA.  A pergunta e' se a faceta devolvida
        // esta' ONDE se pediu.  Ligada por FI_DIAG_SUPORTE para nao poluir
        // corrida normal; imprime so' o primeiro marcador de cada chamada.
        if (_diag_suporte) {
            real dist = 0.0;
            for (int d = 0; d < DIM; d++) dist += (cf[d]-p[d])*(cf[d]-p[d]);
            dist = sqrt(dist);
            fprintf(stderr, "  SUP dim=%d pedido=(%.5f,%.5f) achado=(%.5f,%.5f) "
                            "dist=%.3e (h=%.4f)  lid=%d  peso=%.4f\n",
                    dim, (double) p[0], (double) p[1],
                    (double) cf[0], (double) cf[1], (double) dist, (double) h,
                    lid, (double) w);
        }
        lids[n]  = lid;
        pesos[n] = w;
        n++;
    }
    return n;
}

// ------------------------------------------------------------ interpolacao

void fi_interpola(fi_corpo *c, sim_facet_domain *sfd[DIM],
                  distributed_property *dpu[DIM])
{
    PetscInt n = 0;
    PetscReal *pos = NULL, *vel = NULL, *hmar = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmGetField(c->enxame, "h",          NULL, NULL, (void **) &hmar);
    _suporte_achados = 0;      // conta a chamada CORRENTE, nao o acumulado
    static int _ja = 0;
    _diag_suporte = (getenv("FI_DIAG_SUPORTE") != NULL) && (_ja++ % 200 == 0);

    int capac = FI_LARGURA;
    for (int d = 1; d < DIM; d++) capac *= FI_LARGURA;
    int  *lids  = (int  *) malloc(capac * sizeof *lids);
    real *pesos = (real *) malloc(capac * sizeof *pesos);

    for (PetscInt k = 0; k < n; k++) {
        Point X;
        for (int d = 0; d < DIM; d++) X[d] = pos[DIM*k + d];
        for (int dim = 0; dim < DIM; dim++) {
            const int dg = _diag_suporte && (k == 0);
            const int guarda = _diag_suporte; _diag_suporte = dg;
            const int m = _suporte(sfd[dim], dim, X, hmar[k], lids, pesos, capac);
            _diag_suporte = guarda;
            if (dg) fprintf(stderr, "  SUP marcador 0 em (%.5f,%.5f), h=%.4f, %d facetas\n",
                            (double) X[0], (double) X[1], (double) hmar[k], m);
            _suporte_achados += m;
            real u = 0.0;
            // u(X) = SUM u(x) d_h(x-X) h^DIM, e o h^DIM cancela com o 1/h^DIM
            // do nucleo -- sobra a soma ponderada pelo produto de phi.
            for (int i = 0; i < m; i++) u += dp_get_value(dpu[dim], lids[i]) * pesos[i];
            vel[DIM*k + dim] = u;
        }
    }

    free(lids); free(pesos);
    DMSwarmRestoreField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmRestoreField(c->enxame, "h",          NULL, NULL, (void **) &hmar);
}

void fi_forca_corpo_rigido(fi_corpo *c, real dt)
{
    PetscInt n = 0;
    PetscReal *vel = NULL, *f = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmGetField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    // Corpo rigido e FIXO: a velocidade desejada e' zero.
    for (PetscInt k = 0; k < n; k++)
        for (int d = 0; d < DIM; d++) f[DIM*k + d] = -vel[DIM*k + d] / dt;
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmRestoreField(c->enxame, "forca",      NULL, NULL, (void **) &f);
}

// --------------------------------------------- espalhamento, com acumulacao

// O espalhamento escreve em facetas que podem ser de FRANJA -- pertencer a
// outro rank.  O dp_sync manda dono->franja e SOBRESCREVE (pdomain.c), entao
// essas contribuicoes sumiriam caladas.  O PetscSF faz a direcao que falta.
//
// O grafo se monta do que a HiGTree ja' tem: gid_map (id local -> global),
// firstid e local_count de cada rank.  Nada na estrutura dela e' tocado.
static PetscSF _grafo_de(distributed_property *dp, MPI_Comm comm,
                         int **lids_franja, int *n_franja)
{
    const _dp_shared *sh = dp->pdata;
    const int nlocal = sh->local_count;
    const int ntotal = sh->total_count;
    const int nfr    = ntotal - nlocal;

    int np, eu;
    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &eu);

    int *primeiro = (int *) malloc(np * sizeof *primeiro);
    int *quantos  = (int *) malloc(np * sizeof *quantos);
    int meu_primeiro = sh->firstid, meu_quantos = nlocal;
    MPI_Allgather(&meu_primeiro, 1, MPI_INT, primeiro, 1, MPI_INT, comm);
    MPI_Allgather(&meu_quantos,  1, MPI_INT, quantos,  1, MPI_INT, comm);

    PetscSFNode *remoto = (PetscSFNode *) malloc((nfr > 0 ? nfr : 1) * sizeof *remoto);
    int *lids = (int *) malloc((nfr > 0 ? nfr : 1) * sizeof *lids);

    int n = 0;
    for (int lid = nlocal; lid < ntotal; lid++) {
        const int gid = sh->gid_map[lid];
        int dono = -1;
        for (int r = 0; r < np; r++)
            if (gid >= primeiro[r] && gid < primeiro[r] + quantos[r]) { dono = r; break; }
        if (dono < 0 || dono == eu) continue;    // sem dono conhecido: fica local
        remoto[n].rank  = dono;
        remoto[n].index = gid - primeiro[dono];
        lids[n] = lid;
        n++;
    }

    PetscSF sf;
    PetscCallAbort(comm, PetscSFCreate(comm, &sf));
    PetscCallAbort(comm, PetscSFSetGraph(sf, nlocal, n, NULL, PETSC_COPY_VALUES,
                                         remoto, PETSC_COPY_VALUES));
    PetscCallAbort(comm, PetscSFSetUp(sf));

    free(primeiro); free(quantos); free(remoto);
    *lids_franja = lids;
    *n_franja    = n;
    return sf;
}

void fi_espalha(fi_corpo *c, sim_facet_domain *sfd[DIM],
                distributed_property *dpF[DIM])
{
    fi_espalha_com_escala(c, sfd, dpF, 1.0);
}

void fi_espalha_com_escala(fi_corpo *c, sim_facet_domain *sfd[DIM],
                           distributed_property *dpF[DIM], real escala)
{
    PetscInt n = 0;
    PetscReal *pos = NULL, *f = NULL, *peso = NULL, *hmar = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao", NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "forca",   NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",    NULL, NULL, (void **) &peso);
    DMSwarmGetField(c->enxame, "h",       NULL, NULL, (void **) &hmar);

    int capac = FI_LARGURA;
    for (int d = 1; d < DIM; d++) capac *= FI_LARGURA;
    int  *lids  = (int  *) malloc(capac * sizeof *lids);
    real *pesos = (real *) malloc(capac * sizeof *pesos);

    // VOLUME DO MARCADOR, e e' aqui que a geometria vira fisica.
    //
    // O peso guardado e' `ds`, COMPRIMENTO de arco -- geometrico, e por isso
    // testavel contra o perimetro com exatidao.  Mas o que o espalhamento pede
    // e' o VOLUME que o marcador representa: em 2D, ds*h (a curva com espessura
    // de uma celula); em 3D, dA*1.  Dai' dV = peso * h^(DIM-1).
    //
    // Usar `ds` direto faz a forca sair 1/h vezes maior -- vinte vezes, nesta
    // malha -- e o resultado NAO e' erro visivel: e' realimentacao com ganho, e
    // a velocidade explode em algumas dezenas de passos.
    //
    // A CONSERVACAO NAO PEGA ISTO.  Ela afirma SUM F h^DIM = SUM f w, que e'
    // identidade em w seja qual for o significado dele.  Quem pega e' o
    // acoplamento com a equacao, ou uma conta de unidades.

    for (PetscInt k = 0; k < n; k++) {
        Point X;
        for (int d = 0; d < DIM; d++) X[d] = pos[DIM*k + d];
        // h^DIM do MARCADOR: numa malha graduada ele muda de marcador para
        // marcador, e usar um valor global erraria a normalizacao do nucleo.
        real hd = 1.0;
        for (int d = 0; d < DIM; d++) hd *= hmar[k];

        for (int dim = 0; dim < DIM; dim++) {
            const int m = _suporte(sfd[dim], dim, X, hmar[k], lids, pesos, capac);
            // F(x) = SUM f_k d_h(x-X_k) w_k, com d_h = (1/h^DIM) prod phi.
            const real esc = escala * f[DIM*k + dim]
                             * _volume_marcador(c, peso[k], hmar[k]) / hd;
            for (int i = 0; i < m; i++)
                dp_add_value(dpF[dim], lids[i], esc * pesos[i]);
        }
    }

    free(lids); free(pesos);
    DMSwarmRestoreField(c->enxame, "posicao", NULL, NULL, (void **) &pos);
    DMSwarmRestoreField(c->enxame, "forca",   NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",    NULL, NULL, (void **) &peso);
    DMSwarmRestoreField(c->enxame, "h",       NULL, NULL, (void **) &hmar);

    // Agora a parte que o dp_sync nao faz: somar a franja no dono.
    for (int dim = 0; dim < DIM; dim++) {
        int *fr = NULL, nfr = 0;
        PetscSF sf = _grafo_de(dpF[dim], c->comm, &fr, &nfr);

        PetscReal *folha = (PetscReal *) malloc((nfr > 0 ? nfr : 1) * sizeof *folha);
        for (int i = 0; i < nfr; i++) folha[i] = dp_get_value(dpF[dim], fr[i]);

        PetscCallAbort(c->comm, PetscSFReduceBegin(sf, MPIU_REAL, folha,
                                                   dpF[dim]->values, MPI_SUM));
        PetscCallAbort(c->comm, PetscSFReduceEnd  (sf, MPIU_REAL, folha,
                                                   dpF[dim]->values, MPI_SUM));

        // A franja ja' entregou o que tinha; zera para nao entregar duas vezes
        // se `fi_espalha` for chamada de novo antes de um sync.
        for (int i = 0; i < nfr; i++) dp_set_value(dpF[dim], fr[i], 0.0);

        free(folha); free(fr);
        PetscSFDestroy(&sf);

        // Agora o dono tem o total; o sync reparte de volta para as franjas,
        // que e' a direcao que ele sabe fazer.
        dp_sync(dpF[dim]);
    }
}
