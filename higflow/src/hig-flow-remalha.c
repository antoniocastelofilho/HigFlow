#include "hig-flow-remalha.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "higtree.h"
#include "higtree-iterator.h"

// A QUANTIZACAO.  A posicao vira um inteiro para poder ser comparada e
// espalhada.  O passo tem de ser MUITO mais fino que a menor celula -- para nao
// juntar duas celulas distintas -- e MUITO mais grosso que o ruido de ponto
// flutuante -- para nao separar a mesma celula vista pelos dois lados.
//
// Com coordenada da ordem de 20, o ulp e' ~4e-15; o passo de 1e-9 esta' seis
// ordens acima disso e seis ordens abaixo de qualquer celula util.  A margem e'
// larga dos dois lados de proposito: ela e' o que permite que a chave seja
// EXATA sem ser fragil.
//
// Nao e' tolerancia de comparacao: duas posicoes que caiam em lados opostos de
// um degrau da grade nao se encontram.  No caso para o qual isto foi escrito --
// a mesma construcao geometrica dos dois lados -- as coordenadas sao bit a bit
// iguais e a questao nao se poe.  Para malha que MUDOU, o caminho nao e'
// afrouxar a chave, e' interpolar; ver o cabecalho.
#define REM_PASSO 1.0e-9

typedef struct {
    long long k[DIM];
    real      v;
    int       origem;     // de que rank o valor veio; so' para o diagnostico
} Par;

typedef struct {
    long long k[DIM];
} Chave;

struct rem_colheita {
    Par     *hosp;        // o que ESTE rank hospeda, ordenado por chave
    long     n_hosp;
    MPI_Comm comm;
};

static long _vieram_de_fora = 0;
static long _posicoes       = 0;
static long _colisoes       = 0;

long rem_vieram_de_outro_rank(void) { return _vieram_de_fora; }
long rem_posicoes_hospedadas(void)  { return _posicoes; }
long rem_colisoes(void)             { return _colisoes; }

static void _chave(const Point x, long long k[DIM])
{
    for (int d = 0; d < DIM; d++) k[d] = (long long) llround(x[d] / REM_PASSO);
}

// O anfitriao da posicao.  Qualquer funcao serve, desde que TODOS os ranks
// computem a mesma; esta e' um FNV-1a sobre os bytes da chave, que espalha bem
// e nao tem periodicidade alinhada com a malha -- um `soma % P` poria planos
// inteiros no mesmo rank.
static int _anfitriao(const long long k[DIM], int ntasks)
{
    unsigned long long h = 1469598103934665603ULL;
    const unsigned char *b = (const unsigned char *) k;
    for (size_t i = 0; i < sizeof(long long) * DIM; i++) {
        h ^= b[i];
        h *= 1099511628211ULL;
    }
    return (int) (h % (unsigned long long) ntasks);
}

static int _cmp_chave(const long long a[DIM], const long long b[DIM])
{
    for (int d = 0; d < DIM; d++) {
        if (a[d] < b[d]) return -1;
        if (a[d] > b[d]) return  1;
    }
    return 0;
}

static int _cmp_par(const void *pa, const void *pb)
{
    return _cmp_chave(((const Par *) pa)->k, ((const Par *) pb)->k);
}

// Busca binaria na tabela do anfitriao.
static const Par *_acha(const Par *t, long n, const long long k[DIM])
{
    long lo = 0, hi = n - 1;
    while (lo <= hi) {
        const long m = lo + (hi - lo) / 2;
        const int c = _cmp_chave(t[m].k, k);
        if (c == 0) return &t[m];
        if (c < 0) lo = m + 1; else hi = m - 1;
    }
    return NULL;
}

// -----------------------------------------------------------------------------
// A troca todos-para-todos, feita uma vez e usada pelos dois lados.
//
// Recebe itens ja' etiquetados com o destino, devolve o que chegou.  `tam` e' o
// tamanho de um item em bytes, para servir a `Par` e a `Chave` sem duplicar o
// codigo de contagem e deslocamento -- que e' onde se erra.
// -----------------------------------------------------------------------------
static void *_troca(const void *envio, const int *destino, long n, size_t tam,
                    MPI_Comm comm, long *n_recebido,
                    int *env_cnt, int *env_desl, int *rec_cnt, int *rec_desl)
{
    int ntasks;
    MPI_Comm_size(comm, &ntasks);

    for (int p = 0; p < ntasks; p++) env_cnt[p] = 0;
    for (long i = 0; i < n; i++) env_cnt[destino[i]]++;

    MPI_Alltoall(env_cnt, 1, MPI_INT, rec_cnt, 1, MPI_INT, comm);

    env_desl[0] = rec_desl[0] = 0;
    for (int p = 1; p < ntasks; p++) {
        env_desl[p] = env_desl[p-1] + env_cnt[p-1];
        rec_desl[p] = rec_desl[p-1] + rec_cnt[p-1];
    }
    const long nrec = (long) rec_desl[ntasks-1] + rec_cnt[ntasks-1];

    // Ordenar o envio por destino, preservando a ordem dentro de cada destino --
    // e' ela que permite casar a resposta com o pedido, mais adiante.
    char *buf = (char *) malloc ((size_t) (n > 0 ? n : 1) * tam);
    int *cursor = (int *) malloc ((size_t) ntasks * sizeof *cursor);
    memcpy (cursor, env_desl, (size_t) ntasks * sizeof *cursor);
    for (long i = 0; i < n; i++)
        memcpy (buf + (size_t) cursor[destino[i]]++ * tam,
                (const char *) envio + (size_t) i * tam, tam);
    free (cursor);

    char *rec = (char *) malloc ((size_t) (nrec > 0 ? nrec : 1) * tam);

    // MPI_Alltoallv conta em ELEMENTOS de um tipo; aqui o tipo e' o item inteiro
    // em bytes, criado uma vez.
    MPI_Datatype item;
    MPI_Type_contiguous ((int) tam, MPI_BYTE, &item);
    MPI_Type_commit (&item);
    MPI_Alltoallv (buf, env_cnt, env_desl, item, rec, rec_cnt, rec_desl, item, comm);
    MPI_Type_free (&item);

    free (buf);
    *n_recebido = nrec;
    return rec;
}

// -----------------------------------------------------------------------------
// Colher
// -----------------------------------------------------------------------------
static rem_colheita *_colhe(const Par *local, long n, MPI_Comm comm)
{
    int ntasks;
    MPI_Comm_size (comm, &ntasks);

    int *destino = (int *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *destino);
    for (long i = 0; i < n; i++) destino[i] = _anfitriao (local[i].k, ntasks);

    int *a = (int *) malloc ((size_t) ntasks * sizeof *a);
    int *b = (int *) malloc ((size_t) ntasks * sizeof *b);
    int *c = (int *) malloc ((size_t) ntasks * sizeof *c);
    int *d = (int *) malloc ((size_t) ntasks * sizeof *d);
    long nrec = 0;
    Par *rec = (Par *) _troca (local, destino, n, sizeof (Par), comm, &nrec,
                               a, b, c, d);
    free (destino); free (a); free (b); free (c); free (d);

    qsort (rec, (size_t) nrec, sizeof *rec, _cmp_par);

    // Duas entidades distintas na mesma chave seria defeito da quantizacao, e
    // faria a busca devolver a errada em silencio.  Contado, nao suposto.
    _colisoes = 0;
    for (long i = 1; i < nrec; i++)
        if (_cmp_chave (rec[i-1].k, rec[i].k) == 0) _colisoes++;
    _posicoes = nrec;

    rem_colheita *h = (rem_colheita *) malloc (sizeof *h);
    h->hosp = rec;
    h->n_hosp = nrec;
    h->comm = comm;
    return h;
}

rem_colheita *rem_colhe_centro(sim_domain *sd, distributed_property *dp)
{
    mp_mapper *m = sd_get_domain_mapper (sd);
    const int nloc = dp->pdata->local_count;
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank (comm, &rank);

    Par *local = (Par *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *local);
    long n = 0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator (sd); !higcit_isfinished (it);
         higcit_nextcell (it)) {
        hig_cell *cel = higcit_getcell (it);
        const int lid = mp_lookup (m, hig_get_cid (cel));
        if (lid < 0 || lid >= nloc) continue;
        Point cc;
        hig_get_center (cel, cc);
        _chave (cc, local[n].k);
        local[n].v = dp_get_value (dp, lid);
        local[n].origem = rank;
        n++;
    }
    higcit_destroy (it);

    rem_colheita *h = _colhe (local, n, comm);
    free (local);
    return h;
}

rem_colheita *rem_colhe_faceta(sim_facet_domain *sfd, distributed_property *dp)
{
    mp_mapper *m = sfd->fm;
    const int nloc = dp->pdata->local_count;
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank (comm, &rank);

    Par *local = (Par *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *local);
    long n = 0;
    sim_domain *cd = sfd->cdom;
    for (int k = 0; k < sd_get_num_higtrees (cd); k++) {
        hig_cell *root = sd_get_higtree (cd, k);
        Point blo, bhi;
        POINT_ASSIGN_SCALAR (blo, -1.0e30);
        POINT_ASSIGN_SCALAR (bhi,  1.0e30);
        higfit_facetiterator *fit;
        for (fit = higfit_create_bounding_box_facets (root, sfd->dimofinterest,
                                                      blo, bhi);
             !higfit_isfinished (fit); higfit_nextfacet (fit)) {
            hig_facet *f = higfit_getfacet (fit);
            const int lid = mp_lookup (m, hig_get_fid (f));
            if (lid < 0 || lid >= nloc) continue;
            Point fc;
            hig_get_facet_center (f, fc);
            _chave (fc, local[n].k);
            local[n].v = dp_get_value (dp, lid);
            local[n].origem = rank;
            n++;
        }
        higfit_destroy (fit);
    }

    rem_colheita *h = _colhe (local, n, comm);
    free (local);
    return h;
}

// -----------------------------------------------------------------------------
// Plantar
// -----------------------------------------------------------------------------
//
// Tres trocas: o pedido vai ao anfitriao, a resposta volta, e o valor cai no id
// local certo.  O que amarra o pedido a' resposta e' a ORDEM: `_troca` preserva
// a ordem dentro de cada destino, entao a i-esima resposta vinda do rank p
// corresponde ao i-esimo pedido enviado a p.
static long _planta(rem_colheita *h, const Chave *pedido, long n,
                    real *saida_por_indice, char *achado_por_indice)
{
    int ntasks, rank;
    MPI_Comm_size (h->comm, &ntasks);
    MPI_Comm_rank (h->comm, &rank);

    int *destino = (int *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *destino);
    for (long i = 0; i < n; i++) destino[i] = _anfitriao (pedido[i].k, ntasks);

    int *env_cnt  = (int *) malloc ((size_t) ntasks * sizeof *env_cnt);
    int *env_desl = (int *) malloc ((size_t) ntasks * sizeof *env_desl);
    int *rec_cnt  = (int *) malloc ((size_t) ntasks * sizeof *rec_cnt);
    int *rec_desl = (int *) malloc ((size_t) ntasks * sizeof *rec_desl);

    long nped = 0;
    Chave *chegou = (Chave *) _troca (pedido, destino, n, sizeof (Chave),
                                      h->comm, &nped,
                                      env_cnt, env_desl, rec_cnt, rec_desl);

    // O anfitriao responde, na MESMA ordem em que os pedidos chegaram.
    Par *resp = (Par *) malloc ((size_t) (nped > 0 ? nped : 1) * sizeof *resp);
    for (long i = 0; i < nped; i++) {
        const Par *p = _acha (h->hosp, h->n_hosp, chegou[i].k);
        if (p != NULL) {
            resp[i].v = p->v;
            resp[i].origem = p->origem;
        } else {
            resp[i].v = 0.0;
            resp[i].origem = -1;          // nao achado, e dito
        }
        memcpy (resp[i].k, chegou[i].k, sizeof resp[i].k);
    }
    free (chegou);

    // A volta: os papeis de contagem e deslocamento se invertem.
    long nvolta = 0;
    {
        const long prev = (long) env_desl[ntasks-1] + env_cnt[ntasks-1];
        Par *devolvido = (Par *) malloc ((size_t) (prev > 0 ? prev : 1) * sizeof *devolvido);
        MPI_Datatype item;
        MPI_Type_contiguous ((int) sizeof (Par), MPI_BYTE, &item);
        MPI_Type_commit (&item);
        MPI_Alltoallv (resp, rec_cnt, rec_desl, item,
                       devolvido, env_cnt, env_desl, item, h->comm);
        MPI_Type_free (&item);
        nvolta = prev;

        // Desfazer a ordenacao por destino: percorrer os pedidos na ordem
        // original e consumir a resposta do bloco do destino correspondente.
        int *cursor = (int *) malloc ((size_t) ntasks * sizeof *cursor);
        memcpy (cursor, env_desl, (size_t) ntasks * sizeof *cursor);
        long perdidas = 0;
        _vieram_de_fora = 0;
        for (long i = 0; i < n; i++) {
            const Par *r = &devolvido[cursor[destino[i]]++];
            if (r->origem < 0) {
                perdidas++;
                saida_por_indice[i] = 0.0;
                if (achado_por_indice) achado_por_indice[i] = 0;
                continue;
            }
            saida_por_indice[i] = r->v;
            if (achado_por_indice) achado_por_indice[i] = 1;
            if (r->origem != rank) _vieram_de_fora++;
        }
        free (cursor);
        free (devolvido);
        free (resp); free (destino);
        free (env_cnt); free (env_desl); free (rec_cnt); free (rec_desl);
        (void) nvolta;
        return perdidas;
    }
}

long rem_planta_centro(rem_colheita *h, sim_domain *sd,
                       distributed_property *dp, char *achado)
{
    mp_mapper *m = sd_get_domain_mapper (sd);
    const int nloc = dp->pdata->local_count;

    Chave *pedido = (Chave *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *pedido);
    int   *lids   = (int *)   malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *lids);
    long n = 0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator (sd); !higcit_isfinished (it);
         higcit_nextcell (it)) {
        hig_cell *cel = higcit_getcell (it);
        const int lid = mp_lookup (m, hig_get_cid (cel));
        if (lid < 0 || lid >= nloc) continue;
        Point cc;
        hig_get_center (cel, cc);
        _chave (cc, pedido[n].k);
        lids[n] = lid;
        n++;
    }
    higcit_destroy (it);

    real *val = (real *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *val);
    char *ach = (char *) malloc ((size_t) (n > 0 ? n : 1));
    const long perdidas = _planta (h, pedido, n, val, ach);
    for (long i = 0; i < n; i++) {
        dp_set_value (dp, lids[i], val[i]);
        if (achado) achado[lids[i]] = ach[i];   // por ID LOCAL, nao por ordem
    }

    free (ach); free (val); free (lids); free (pedido);
    return perdidas;
}

long rem_planta_faceta(rem_colheita *h, sim_facet_domain *sfd,
                       distributed_property *dp, char *achado)
{
    mp_mapper *m = sfd->fm;
    const int nloc = dp->pdata->local_count;

    Chave *pedido = (Chave *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *pedido);
    int   *lids   = (int *)   malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *lids);
    long n = 0;
    sim_domain *cd = sfd->cdom;
    for (int k = 0; k < sd_get_num_higtrees (cd); k++) {
        hig_cell *root = sd_get_higtree (cd, k);
        Point blo, bhi;
        POINT_ASSIGN_SCALAR (blo, -1.0e30);
        POINT_ASSIGN_SCALAR (bhi,  1.0e30);
        higfit_facetiterator *fit;
        for (fit = higfit_create_bounding_box_facets (root, sfd->dimofinterest,
                                                      blo, bhi);
             !higfit_isfinished (fit); higfit_nextfacet (fit)) {
            hig_facet *f = higfit_getfacet (fit);
            const int lid = mp_lookup (m, hig_get_fid (f));
            if (lid < 0 || lid >= nloc) continue;
            Point fc;
            hig_get_facet_center (f, fc);
            _chave (fc, pedido[n].k);
            lids[n] = lid;
            n++;
        }
        higfit_destroy (fit);
    }

    real *val = (real *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *val);
    char *ach = (char *) malloc ((size_t) (n > 0 ? n : 1));
    const long perdidas = _planta (h, pedido, n, val, ach);
    for (long i = 0; i < n; i++) {
        dp_set_value (dp, lids[i], val[i]);
        if (achado) achado[lids[i]] = ach[i];   // por ID LOCAL, nao por ordem
    }

    free (ach); free (val); free (lids); free (pedido);
    return perdidas;
}

void rem_destroi(rem_colheita *h)
{
    if (h == NULL) return;
    free (h->hosp);
    free (h);
}
