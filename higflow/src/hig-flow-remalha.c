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

// -----------------------------------------------------------------------------
// Interpolacao conservativa
// -----------------------------------------------------------------------------

static long _por_refino = 0, _por_engrossamento = 0, _parciais = 0;

void rem_interpolou(long *r, long *e)
{
    if (r) *r = _por_refino;
    if (e) *e = _por_engrossamento;
}

long rem_engrossamento_parcial(void) { return _parciais; }

// Consulta em lote: `nk` chaves por item, `n` itens.  Devolve, por chave, se foi
// achada e o valor.  E' o `_planta` sem a escrita, e sem amarrar a id local.
static void _consulta(rem_colheita *h, const Chave *chaves, long total,
                      char *achou, real *valor)
{
    int ntasks, rank;
    MPI_Comm_size (h->comm, &ntasks);
    MPI_Comm_rank (h->comm, &rank);

    int *destino = (int *) malloc ((size_t) (total > 0 ? total : 1) * sizeof *destino);
    for (long i = 0; i < total; i++) destino[i] = _anfitriao (chaves[i].k, ntasks);

    int *env_cnt  = (int *) malloc ((size_t) ntasks * sizeof *env_cnt);
    int *env_desl = (int *) malloc ((size_t) ntasks * sizeof *env_desl);
    int *rec_cnt  = (int *) malloc ((size_t) ntasks * sizeof *rec_cnt);
    int *rec_desl = (int *) malloc ((size_t) ntasks * sizeof *rec_desl);

    long nped = 0;
    Chave *chegou = (Chave *) _troca (chaves, destino, total, sizeof (Chave),
                                      h->comm, &nped,
                                      env_cnt, env_desl, rec_cnt, rec_desl);

    Par *resp = (Par *) malloc ((size_t) (nped > 0 ? nped : 1) * sizeof *resp);
    for (long i = 0; i < nped; i++) {
        const Par *p = _acha (h->hosp, h->n_hosp, chegou[i].k);
        resp[i].v = (p != NULL) ? p->v : 0.0;
        resp[i].origem = (p != NULL) ? p->origem : -1;
        memcpy (resp[i].k, chegou[i].k, sizeof resp[i].k);
    }
    free (chegou);

    const long prev = (long) env_desl[ntasks-1] + env_cnt[ntasks-1];
    Par *devolvido = (Par *) malloc ((size_t) (prev > 0 ? prev : 1) * sizeof *devolvido);
    MPI_Datatype item;
    MPI_Type_contiguous ((int) sizeof (Par), MPI_BYTE, &item);
    MPI_Type_commit (&item);
    MPI_Alltoallv (resp, rec_cnt, rec_desl, item,
                   devolvido, env_cnt, env_desl, item, h->comm);
    MPI_Type_free (&item);

    int *cursor = (int *) malloc ((size_t) ntasks * sizeof *cursor);
    memcpy (cursor, env_desl, (size_t) ntasks * sizeof *cursor);
    for (long i = 0; i < total; i++) {
        const Par *r = &devolvido[cursor[destino[i]]++];
        achou[i] = (r->origem >= 0);
        valor[i] = r->v;
    }
    free (cursor); free (devolvido); free (resp); free (destino);
    free (env_cnt); free (env_desl); free (rec_cnt); free (rec_desl);
}

// O nucleo, comum a celula e faceta: dadas as posicoes que faltam e o tamanho
// de cada uma, tenta PAI (2^DIM candidatos, no maximo um existe) e depois
// FILHAS (2^DIM, todas tem de existir).
//
// `dir_livre[d]` diz em que direcoes a entidade se subdivide.  Para celula sao
// todas; para faceta, todas MENOS a normal -- uma faceta nao se parte na
// direcao em que ela e' plana.
// `exige_todas` distingue celula de faceta.  Para CELULA as 2^DIM filhas sempre
// existem, e aceitar menos mascararia defeito.  Para FACETA nao: numa interface
// de refino a malha de origem pode nao carregar as sub-facetas como graus de
// liberdade -- MEDIDO, a malha refinada do teste tem 16 facetas no plano da
// interface onde a grossa tem 24.
static long _interpola_nucleo(rem_colheita *h, const Point *pos, const Point *tam,
                              const int *dir_livre, long n, real *saida, char *ok,
                              int exige_todas)
{
    // SEM SAIDA ANTECIPADA POR n == 0.  `_consulta` e' COLETIVA: um rank que
    // nao tenha nada a interpolar e volte aqui deixa os outros esperando para
    // sempre no MPI_Alltoall.  Eu escrevi o aviso disso em `rem_interpola_centro`
    // e mesmo assim pus um `if (n <= 0) return` aqui -- e travou em np=3.
    // Com n == 0 os lacos nao iteram e a troca vai com contagem zero, que e'
    // legitima.

    int nlivres = 0;
    for (int d = 0; d < DIM; d++) if (dir_livre[d]) nlivres++;
    const int nc = 1 << nlivres;          // candidatos a pai, e filhas

    Chave *q = (Chave *) malloc ((size_t) n * nc * sizeof *q);
    char  *a = (char *)  malloc ((size_t) n * nc);
    real  *v = (real *)  malloc ((size_t) n * nc * sizeof *v);

    // --- PAI: centro +- tam/2 nas direcoes livres -----------------------------
    for (long i = 0; i < n; i++) {
        for (int c = 0; c < nc; c++) {
            Point x;
            POINT_ASSIGN (x, pos[i]);
            int bit = 0;
            for (int d = 0; d < DIM; d++) {
                if (!dir_livre[d]) continue;
                x[d] += ((c >> bit) & 1) ? 0.5 * tam[i][d] : -0.5 * tam[i][d];
                bit++;
            }
            _chave (x, q[i * nc + c].k);
        }
    }
    _consulta (h, q, (long) n * nc, a, v);

    long faltam = 0;
    for (long i = 0; i < n; i++) {
        if (ok[i]) continue;
        int achados = 0;
        real val = 0.0;
        for (int c = 0; c < nc; c++)
            if (a[i * nc + c]) { achados++; val = v[i * nc + c]; }
        if (achados == 1) {          // exatamente um pai: o caso sem ambiguidade
            saida[i] = val;
            ok[i] = 1;
            _por_refino++;
        }
    }

    // --- FILHAS: centro +- tam/4 nas direcoes livres, TODAS necessarias -------
    for (long i = 0; i < n; i++) {
        for (int c = 0; c < nc; c++) {
            Point x;
            POINT_ASSIGN (x, pos[i]);
            int bit = 0;
            for (int d = 0; d < DIM; d++) {
                if (!dir_livre[d]) continue;
                x[d] += ((c >> bit) & 1) ? 0.25 * tam[i][d] : -0.25 * tam[i][d];
                bit++;
            }
            _chave (x, q[i * nc + c].k);
        }
    }
    _consulta (h, q, (long) n * nc, a, v);

    for (long i = 0; i < n; i++) {
        if (ok[i]) continue;
        int achados = 0;
        real soma = 0.0;
        for (int c = 0; c < nc; c++)
            if (a[i * nc + c]) { achados++; soma += v[i * nc + c]; }
        if (achados == nc) {         // volumes iguais: media simples E' em volume
            saida[i] = soma / (real) nc;
            ok[i] = 1;
            _por_engrossamento++;
        } else if (!exige_todas && achados > 0) {
            // Interface de refino: a malha de origem pode nao carregar todas as
            // sub-facetas como graus de liberdade.  A media do que existe e' a
            // melhor informacao disponivel -- exata para campo constante -- e
            // fica contada a' parte, para ser escolha visivel e nao silencio.
            saida[i] = soma / (real) achados;
            ok[i] = 1;
            _por_engrossamento++;
            _parciais++;
        }
    }

    for (long i = 0; i < n; i++) if (!ok[i]) faltam++;
    free (v); free (a); free (q);
    return faltam;
}

long rem_interpola_centro(rem_colheita *h, sim_domain *sd,
                          distributed_property *dp, char *achado)
{
    mp_mapper *m = sd_get_domain_mapper (sd);
    const int nloc = dp->pdata->local_count;

    Point *pos = (Point *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *pos);
    Point *tam = (Point *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *tam);
    int   *lids = (int *)  malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *lids);
    long n = 0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator (sd); !higcit_isfinished (it);
         higcit_nextcell (it)) {
        hig_cell *cel = higcit_getcell (it);
        const int lid = mp_lookup (m, hig_get_cid (cel));
        if (lid < 0 || lid >= nloc) continue;
        if (achado != NULL && achado[lid]) continue;     // ja' veio exato
        Point cc, cl, ch;
        hig_get_center (cel, cc);
        hig_get_lowpoint (cel, cl);
        hig_get_highpoint (cel, ch);
        POINT_ASSIGN (pos[n], cc);
        for (int d = 0; d < DIM; d++) tam[n][d] = ch[d] - cl[d];
        lids[n] = lid;
        n++;
    }
    higcit_destroy (it);

    // TODOS os ranks entram no nucleo, inclusive os que nao tem nada a
    // interpolar: `_consulta` e' COLETIVA.  Um rank que saisse aqui travaria os
    // outros -- e' a mesma armadilha do MPI_Allreduce dentro de `if (rank==0)`.
    int dir[DIM];
    for (int d = 0; d < DIM; d++) dir[d] = 1;
    real *val = (real *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *val);
    char *ok  = (char *) calloc ((size_t) (n > 0 ? n : 1), 1);
    const long faltam = _interpola_nucleo (h, pos, tam, dir, n, val, ok, 1);
    for (long i = 0; i < n; i++)
        if (ok[i]) {
            dp_set_value (dp, lids[i], val[i]);
            if (achado) achado[lids[i]] = 1;
        }

    free (ok); free (val); free (lids); free (tam); free (pos);
    return faltam;
}

long rem_interpola_faceta(rem_colheita *h, sim_facet_domain *sfd,
                          distributed_property *dp, char *achado)
{
    mp_mapper *m = sfd->fm;
    const int nloc = dp->pdata->local_count;
    const int dim = sfd_get_dim (sfd);

    Point *pos = (Point *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *pos);
    Point *tam = (Point *) malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *tam);
    int   *lids = (int *)  malloc ((size_t) (nloc > 0 ? nloc : 1) * sizeof *lids);
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
            if (achado != NULL && achado[lid]) continue;
            Point fc, cl, ch;
            hig_get_facet_center (f, fc);
            hig_cell *cel = hig_get_facet_cell (f);
            hig_get_lowpoint (cel, cl);
            hig_get_highpoint (cel, ch);
            POINT_ASSIGN (pos[n], fc);
            for (int d = 0; d < DIM; d++) tam[n][d] = ch[d] - cl[d];
            lids[n] = lid;
            n++;
        }
        higfit_destroy (fit);
    }

    // A faceta nao se parte na direcao NORMAL a ela.
    int dir[DIM];
    for (int d = 0; d < DIM; d++) dir[d] = (d == dim) ? 0 : 1;

    real *val = (real *) malloc ((size_t) (n > 0 ? n : 1) * sizeof *val);
    char *ok  = (char *) calloc ((size_t) (n > 0 ? n : 1), 1);
    long faltam = _interpola_nucleo (h, pos, tam, dir, n, val, ok, 0);

    // O CASO QUE A CELULA NAO TEM: a faceta do plano do MEIO.  Refinar uma
    // celula cria, na direcao normal, uma faceta que no nivel grosso era
    // interior -- ela nao e' sub-faceta de faceta nenhuma, e o laco acima nao a
    // acha.  Ela recebe a media das duas facetas paralelas que limitavam a
    // celula antiga, a +- tam na direcao normal, o que preserva o balanco de
    // fluxo atraves dela.
    {
        // As duas facetas paralelas ficam a +- tam na direcao normal, E a
        // +- tam/2 nas TRANSVERSAIS: a faceta do meio pertence a uma celula
        // fina cujo centro transversal esta' deslocado do centro da mae.
        // Esquecer esse deslocamento foi o primeiro defeito aqui.
        // CANDIDATOS TRANSVERSAIS, em ordem de prioridade.  O primeiro e' o
        // deslocamento ZERO -- a vizinha na direcao normal, na MESMA posicao
        // transversal.  Depois vem os +- tam/2, que sao os pais quando a
        // faceta e' de fato do plano do meio.
        //
        // Deslocar SEMPRE era defeito: a faceta que simplesmente nao tem
        // contraparte na origem (interface de refino) tem vizinha na mesma
        // transversal, e a busca so' com deslocamento nunca a achava.  MEDIDO:
        // 8 facetas com ZERO lados achados, todas no plano da interface.
        const int nt = 1 + (1 << (DIM - 1));
        const long nq = (long) n * 2 * nt;
        Chave *q = (Chave *) malloc ((size_t) (nq > 0 ? nq : 1) * sizeof *q);
        char  *a = (char *)  malloc ((size_t) (nq > 0 ? nq : 1));
        real  *v = (real *)  malloc ((size_t) (nq > 0 ? nq : 1) * sizeof *v);
        for (long i = 0; i < n; i++) {
            for (int s = 0; s < 2; s++) {
                for (int c = 0; c < nt; c++) {
                    Point x;
                    POINT_ASSIGN (x, pos[i]);
                    x[dim] += (s ? 1.0 : -1.0) * tam[i][dim];
                    if (c > 0) {              // c == 0 e' o deslocamento zero
                        const int cc = c - 1;
                        int bit = 0;
                        for (int d = 0; d < DIM; d++) {
                            if (d == dim) continue;
                            x[d] += ((cc >> bit) & 1) ? 0.5 * tam[i][d] : -0.5 * tam[i][d];
                            bit++;
                        }
                    }
                    _chave (x, q[(i * 2 + s) * nt + c].k);
                }
            }
        }
        _consulta (h, q, nq, a, v);
        long c0 = 0, c1 = 0, c2 = 0, cmuitos = 0;
        for (long i = 0; i < n; i++) {
            if (ok[i]) continue;
            real lado[2] = {0.0, 0.0};
            int achou_lado[2] = {0, 0};
            for (int s = 0; s < 2; s++) {
                // O PRIMEIRO que casar, na ordem de prioridade: deslocamento
                // zero antes dos deslocados.
                for (int c = 0; c < nt; c++)
                    if (a[(i * 2 + s) * nt + c]) {
                        lado[s] = v[(i * 2 + s) * nt + c];
                        achou_lado[s] = 1;
                        break;
                    }
            }
            const int nl = achou_lado[0] + achou_lado[1];
            if (nl == 0) c0++; else if (nl == 1) c1++; else c2++;
            if (nl == 2) {
                val[i] = 0.5 * (lado[0] + lado[1]);
                ok[i] = 1; _por_refino++; faltam--;
            } else if (nl == 1) {
                // UM LADO SO': a faceta do meio encosta na borda EXTERNA do
                // dominio, e ali nao ha' faceta paralela do outro lado.
                //
                // MEDIDO, e nao suposto: a sonda contou ZERO facetas proprias no
                // plano x = 1,0 da malha de teste.  O iterador de facetas visita
                // uma faceta por celula -- a de baixo --, entao a de cima da
                // ultima celula nao e' visitada por ninguem.  A condicao de
                // contorno ali vive em outro dominio, nao neste.
                //
                // Usar o unico lado disponivel e' exato para campo constante e
                // de primeira ordem no geral.  Fica contado a' parte para que o
                // chamador saiba quantas foram assim.
                // O SINALIZADOR, nao o valor.  Escrito como `lado[0] ? ...`
                // isto testava se o valor achado e' nao nulo -- e um valor
                // legitimamente zero caia no outro lado, que esta' sem
                // inicializar.  Aparecia como 1.000 onde se esperava 107.25.
                val[i] = achou_lado[0] ? lado[0] : lado[1];
                ok[i] = 1; _por_refino++; faltam--;
            }
        }
        if (getenv ("REMALHA_ONDE") != NULL) {
            int rk; MPI_Comm_rank (h->comm, &rk);
            printf ("     [plano do meio] rank %d dim %d: %ld com 0 lados, "
                    "%ld com 1, %ld com 2 (%ld candidatos por lado)\n",
                    rk, dim, c0, c1, c2, (long) nt);
            fflush (stdout);
        }
        (void) cmuitos;
        free (v); free (a); free (q);
    }

    for (long i = 0; i < n; i++)
        if (ok[i]) {
            dp_set_value (dp, lids[i], val[i]);
            if (achado) achado[lids[i]] = 1;
        }

    free (ok); free (val); free (lids); free (tam); free (pos);
    return faltam;
}

void rem_destroi(rem_colheita *h)
{
    if (h == NULL) return;
    free (h->hosp);
    free (h);
}
