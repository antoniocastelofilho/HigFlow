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
    DM       enxame;      // DMSwarm: posicao, peso, velocidade, forca
    real     h;           // espacamento euleriano (malha uniforme)
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
static real _volume_marcador(const fi_corpo *c, real peso)
{
    real v = peso;
    for (int i = 0; i < c->codim; i++) v *= c->h;
    return v;
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
    c->h     = h;
    c->codim = codim;
    c->comm  = MPI_COMM_WORLD;
    for (int d = 0; d < DIM; d++) c->forca_passo[d] = 0.0;

    PetscCallAbort(c->comm, DMCreate(c->comm, &c->enxame));
    PetscCallAbort(c->comm, DMSetType(c->enxame, DMSWARM));
    PetscCallAbort(c->comm, DMSetDimension(c->enxame, DIM));
    PetscCallAbort(c->comm, DMSwarmSetType(c->enxame, DMSWARM_BASIC));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "posicao",    DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "peso",         1, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "velocidade", DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "forca",      DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmFinalizeFieldRegister(c->enxame));

    int meus = 0;
    for (int i = 0; i < ncand; i++) if (_possui(sfd, pts[i])) meus++;

    PetscCallAbort(c->comm, DMSwarmSetLocalSizes(c->enxame, meus, 4));
    PetscReal *pos = NULL, *peso = NULL;
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
    PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "peso",    NULL, NULL, (void **) &peso));
    int k = 0;
    for (int i = 0; i < ncand; i++) {
        if (!_possui(sfd, pts[i])) continue;
        for (int d = 0; d < DIM; d++) pos[DIM*k + d] = pts[i][d];
        peso[k] = pesos[i];
        k++;
    }
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
    PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "peso",    NULL, NULL, (void **) &peso));

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
    PetscReal *pos = NULL, *vel = NULL, *f = NULL, *peso = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmGetField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",       NULL, NULL, (void **) &peso);

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

    DMSwarmRestoreField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmRestoreField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);
    DMSwarmRestoreField(c->enxame, "forca",      NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",       NULL, NULL, (void **) &peso);
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
    PetscReal *f = NULL, *peso = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "forca", NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",  NULL, NULL, (void **) &peso);
    real local[DIM];
    for (int d = 0; d < DIM; d++) local[d] = 0.0;
    for (PetscInt i = 0; i < n; i++)
        for (int d = 0; d < DIM; d++)
            local[d] += f[DIM*i + d] * _volume_marcador(c, peso[i]);
    DMSwarmRestoreField(c->enxame, "forca", NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",  NULL, NULL, (void **) &peso);
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
    PetscReal *pos = NULL, *vel = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao",    NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "velocidade", NULL, NULL, (void **) &vel);

    int capac = FI_LARGURA;
    for (int d = 1; d < DIM; d++) capac *= FI_LARGURA;
    int  *lids  = (int  *) malloc(capac * sizeof *lids);
    real *pesos = (real *) malloc(capac * sizeof *pesos);

    for (PetscInt k = 0; k < n; k++) {
        Point X;
        for (int d = 0; d < DIM; d++) X[d] = pos[DIM*k + d];
        for (int dim = 0; dim < DIM; dim++) {
            const int m = _suporte(sfd[dim], dim, X, c->h, lids, pesos, capac);
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
    PetscReal *pos = NULL, *f = NULL, *peso = NULL;
    DMSwarmGetLocalSize(c->enxame, &n);
    DMSwarmGetField(c->enxame, "posicao", NULL, NULL, (void **) &pos);
    DMSwarmGetField(c->enxame, "forca",   NULL, NULL, (void **) &f);
    DMSwarmGetField(c->enxame, "peso",    NULL, NULL, (void **) &peso);

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
    real hd = 1.0;
    for (int d = 0; d < DIM; d++) hd *= c->h;

    for (PetscInt k = 0; k < n; k++) {
        Point X;
        for (int d = 0; d < DIM; d++) X[d] = pos[DIM*k + d];
        for (int dim = 0; dim < DIM; dim++) {
            const int m = _suporte(sfd[dim], dim, X, c->h, lids, pesos, capac);
            // F(x) = SUM f_k d_h(x-X_k) w_k, com d_h = (1/h^DIM) prod phi.
            const real esc = escala * f[DIM*k + dim] * _volume_marcador(c, peso[k]) / hd;
            for (int i = 0; i < m; i++)
                dp_add_value(dpF[dim], lids[i], esc * pesos[i]);
        }
    }

    free(lids); free(pesos);
    DMSwarmRestoreField(c->enxame, "posicao", NULL, NULL, (void **) &pos);
    DMSwarmRestoreField(c->enxame, "forca",   NULL, NULL, (void **) &f);
    DMSwarmRestoreField(c->enxame, "peso",    NULL, NULL, (void **) &peso);

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
