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
    MPI_Comm comm;
};

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

fi_corpo *fi_cria_curva(sim_facet_domain *sfd, const Point *vertices, int nvert,
                        real h)
{
    if (nvert < 3 || h <= 0.0) {
        fprintf(stderr, "fi_cria_curva: curva com %d vertices e h=%g\n", nvert, h);
        abort();
    }

    // O PETSc e' inicializado PREGUICOSAMENTE pela HiGTree, ao criar um solver
    // (utils.c).  Este modulo usa PETSc sem solver nenhum, entao a inicializacao
    // pode nao ter acontecido -- e o sintoma e' um erro em MPI_Comm_get_attr,
    // que nao menciona inicializacao.  Declarado no sitio como em
    // higtree/src/solver-petsc.c:362, que e' a convencao da casa.
    void _try_initialize_petsc(void);
    _try_initialize_petsc();

    fi_corpo *c = (fi_corpo *) malloc(sizeof *c);
    c->h    = h;
    c->comm = MPI_COMM_WORLD;

    PetscCallAbort(c->comm, DMCreate(c->comm, &c->enxame));
    PetscCallAbort(c->comm, DMSetType(c->enxame, DMSWARM));
    PetscCallAbort(c->comm, DMSetDimension(c->enxame, DIM));
    PetscCallAbort(c->comm, DMSwarmSetType(c->enxame, DMSWARM_BASIC));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "posicao",    DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "peso",         1, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "velocidade", DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmRegisterPetscDatatypeField(c->enxame, "forca",      DIM, PETSC_REAL));
    PetscCallAbort(c->comm, DMSwarmFinalizeFieldRegister(c->enxame));

    // Primeira passada: contar os marcadores DESTE rank.  Todo rank percorre a
    // curva inteira -- ela e' pequena e analitica, e assim a subdivisao e'
    // identica em todos, sem comunicacao.
    int esperados = 0, meus = 0;
    for (int passada = 0; passada < 2; passada++) {
        PetscReal *pos = NULL, *peso = NULL;
        if (passada == 1) {
            PetscCallAbort(c->comm, DMSwarmSetLocalSizes(c->enxame, meus, 4));
            PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
            PetscCallAbort(c->comm, DMSwarmGetField(c->enxame, "peso",    NULL, NULL, (void **) &peso));
            meus = 0;
        }
        esperados = 0;

        for (int v = 0; v < nvert; v++) {
            const real *a = vertices[v];
            const real *b = vertices[(v + 1) % nvert];   // fecha sozinha
            real comp = 0.0;
            for (int d = 0; d < DIM; d++) comp += (b[d] - a[d]) * (b[d] - a[d]);
            comp = sqrt(comp);
            if (comp <= 0.0) continue;

            int nsub = (int) (comp / h + 0.5);
            if (nsub < 1) nsub = 1;
            const real ds = comp / nsub;

            for (int s = 0; s < nsub; s++) {
                // CENTRO do subsegmento, nao o vertice -- ver o .h.
                const real t = (s + 0.5) / nsub;
                Point x;
                for (int d = 0; d < DIM; d++) x[d] = a[d] + t * (b[d] - a[d]);
                esperados++;
                if (!_possui(sfd, x)) continue;
                if (passada == 1) {
                    for (int d = 0; d < DIM; d++) pos[DIM*meus + d] = x[d];
                    peso[meus] = ds;
                }
                meus++;
            }
        }

        if (passada == 1) {
            PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "posicao", NULL, NULL, (void **) &pos));
            PetscCallAbort(c->comm, DMSwarmRestoreField(c->enxame, "peso",    NULL, NULL, (void **) &peso));
        }
    }

    // Cada marcador tem de ser reivindicado por EXATAMENTE um rank.  Marcador
    // sobre a face entre celulas de ranks diferentes e' o caso que quebra isso,
    // e seguir daria forca errada sem sintoma visivel.
    int soma = 0;
    MPI_Allreduce(&meus, &soma, 1, MPI_INT, MPI_SUM, c->comm);
    if (soma != esperados) {
        int rank; MPI_Comm_rank(c->comm, &rank);
        if (rank == 0)
            fprintf(stderr,
                "fi_cria_curva: %d marcadores no total, mas os ranks reivindicaram %d.  "
                "Marcador sobre fronteira de particao reivindicado por dois ranks (ou por "
                "nenhum).  Seguir daria forca errada sem sintoma visivel.\n",
                esperados, soma);
        MPI_Barrier(c->comm);
        abort();
    }

    return c;
}

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
        for (int d = 0; d < DIM; d++) local[d] += f[DIM*i + d] * peso[i];
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
        if (!sfd_get_facet_with_point(sfd, p, &f)) continue;
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
        lids[n]  = sfd_get_local_id(sfd, &f);
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

void fi_calcula_forca(fi_corpo *c, real dt)
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

    real hd = 1.0;
    for (int d = 0; d < DIM; d++) hd *= c->h;

    for (PetscInt k = 0; k < n; k++) {
        Point X;
        for (int d = 0; d < DIM; d++) X[d] = pos[DIM*k + d];
        for (int dim = 0; dim < DIM; dim++) {
            const int m = _suporte(sfd[dim], dim, X, c->h, lids, pesos, capac);
            // F(x) = SUM f_k d_h(x-X_k) w_k, com d_h = (1/h^DIM) prod phi.
            const real esc = f[DIM*k + dim] * peso[k] / hd;
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
