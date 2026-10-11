#include "vof-amr-3d.h"
#include "malha-adaptativa.h"
#include "../src/hig-flow-amr-criterio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ------------------------------------------------------------ banda -> celmin
//
// O piso de cinco vence e AVISA.  Obedecer em silencio a uma banda que pede
// menos de cinco celulas finas por lado quebraria o PLIC e o transporte da
// fracao, e o sintoma apareceria como perda de massa longe daqui.
int vof3_celmin_da_banda(const char *var, int celmin,
                         real lx, int nx, int niveis)
{
    const char *b = getenv(var);
    if (b == NULL) return celmin;
    real banda = atof(b);
    if (!(banda > 0.0)) return celmin;
    real h_fino = lx / (real) nx;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    int novo = (int) (banda / h_fino + 0.5);
    if (novo < 5) {
        print0f("===> %s=%g pede %d celulas finas em h=%.5f, e o PLIC exige 5. "
                "O piso vence: a banda efetiva sera' %.5f, nao %.5f\n",
                var, (double) banda, novo, (double) h_fino,
                (double) (5 * h_fino), (double) banda);
        novo = 5;
    } else {
        print0f("===> banda fisica %.5f em h_fino=%.5f  ==>  celmin = %d (era %d)\n",
                (double) banda, (double) h_fino, novo, celmin);
    }
    return novo;
}

// ------------------------------------------------- os limiares, como no 2D
//
// Mesma conta do `vof_escreve_amr_criterio` bidimensional: o nivel mais fino
// cobre `cel_min` celulas finas de cada lado, mais a margem que a discretizacao
// exige, e cada nivel acima soma a sua propria contagem.
static void _limiares(real h_base, int niveis, int cel_min, real *thr)
{
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    const real margem = h_base + h_fino;
    const real banda_fina = cel_min * h_fino + margem;
    for (int l = niveis - 1; l >= 0; l--) {
        thr[l] = banda_fina;
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        for (int k = 0; k < niveis - 1 - l; k++) thr[l] += cel_min * h_nivel;
    }
}

// ------------------------------------------------------- o criterio dinamico
void vof3_escreve_amr_criterio(higflow_solver *ns, real lx, real ly, real lz,
                               int nx, int ny, int nz,
                               int niveis, int cel_min, const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny, hz = lz / nz;
    real h_base = hx;
    if (hy < h_base) h_base = hy;
    if (hz < h_base) h_base = hz;

    real *thr = (real *) malloc((size_t) niveis * sizeof(real));
    _limiares(h_base, niveis, cel_min, thr);

    // SEMENTES: celulas de interface do CAMPO de fracao.  O laco percorre o
    // dominio LOCAL -- cada rank ve' so' o seu pedaco --, e e' por isso que a
    // reuniao abaixo existe.
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper  *mp  = sd_get_domain_mapper(sdm);
    long cap = 0;
    { // capacidade: o numero de celulas LOCAIS, que limita quantas sementes ha'
      higcit_celliterator *ct;
      for (ct = sd_get_domain_celliterator(sdm); !higcit_isfinished(ct); higcit_nextcell(ct)) cap++;
      higcit_destroy(ct);
    }
    if (cap < 1) cap = 1;
    Point *sem = (Point *) malloc((size_t) cap * sizeof(Point));
    long nsem = 0;

    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        real f;
        if (ns->par.step == 0) {
            // No passo 0 o campo ainda nao foi preenchido; a fracao vem da
            // forma analitica, a MESMA que inicializa o fracvol.  Sem isto o
            // primeiro remalhamento refinaria em torno de um campo nulo.
            Point xc, d; hig_get_center(c, xc); hig_get_delta(c, d);
            f = ns->ed.mult.problem->fracvol(xc, d, ns->par.t);
        } else {
            f = dp_get_value(ns->ed.mult.dpfracvol, clid);
        }
        if (f > 0.001 && f < 0.999 && nsem < cap) {
            hig_get_center(c, sem[nsem]); nsem++;
        }
    }
    higcit_destroy(it);

    // A REUNIAO ENTRE RANKS, e ela NAO e' opcional aqui.  As sementes do VOF
    // saem do campo DISTRIBUIDO, nao de uma geometria replicada: sem a reuniao o
    // .amr sairia refinado apenas em torno do pedaco de quem escreve, e a regra
    // das celulas minimas ficaria violada exatamente onde ninguem esta' olhando.
    // Em np=1 a reuniao e' identidade, logo o caminho serial nao muda.
    long nsem_g = 0;
    Point *sem_g = malha_adapt_reune_sementes(sem, nsem, &nsem_g);
    free(sem);

    // A tabela por DILATACAO, nao por forca bruta.  Em 2D a busca exaustiva
    // O(ncel x nsem) passa; em 3D, 375 mil celulas contra milhares de sementes
    // dariam ~1e9 operacoes POR REMALHAMENTO.  `amr_tabela` e' exata, nao
    // aproximada: o nivel e' monotono decrescente na distancia.
    const long ncel = (long) nx * ny * nz;
    signed char *tab = (signed char *) malloc((size_t) ncel);
    real *plano = (real *) malloc((size_t) (nsem_g > 0 ? nsem_g : 1) * 3 * sizeof(real));
    for (long k = 0; k < nsem_g; k++) {
        plano[3*k+0] = sem_g[k][0];
        plano[3*k+1] = sem_g[k][1];
        plano[3*k+2] = sem_g[k][2];
    }
    amr_tabela(plano, nsem_g, hx, hy, hz, nx, ny, nz, thr, niveis, tab);

    // SO' O RANK 0 ESCREVE, E TODOS ESPERAM.  Sem a guarda todos truncam o
    // MESMO arquivo ao mesmo tempo; sem a barreira alguem o rele' (em
    // higflow_reconstroi_dominio, logo adiante) antes de estar pronto.  Os dois
    // sao necessarios: a guarda sozinha nao faz ninguem esperar.
    int meurank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &meurank);
    if (meurank == 0) {
        if (amr_escreve_arquivo(caminho, lx, ly, lz, nx, ny, nz, niveis, tab) != 0)
            print0f("VOF3 amr: FALHEI ao escrever %s\n", caminho);
        else {
            long marc = 0;
            for (long q = 0; q < ncel; q++) if (tab[q] > 0) marc++;
            print0f("VOF3 amr: %ld sementes (%ld locais) de %ld celulas base, "
                    "%ld marcadas (%.2f%%), niveis=%d, >=%d celulas finas por lado\n",
                    nsem_g, nsem, ncel, marc, 100.0*marc/(double)ncel, niveis, cel_min);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(plano); free(sem_g); free(tab); free(thr);
}

// -------------------------------------------------------- o criterio inicial
void vof3_escreve_amr_inicial(real lx, real ly, real lz, int nx, int ny, int nz,
                              int niveis, int cel_min,
                              real (*fracvol)(Point, Point, real),
                              const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny, hz = lz / nz;
    real h_base = hx;
    if (hy < h_base) h_base = hy;
    if (hz < h_base) h_base = hz;

    real *thr = (real *) malloc((size_t) niveis * sizeof(real));
    _limiares(h_base, niveis, cel_min, thr);

    // Sem dominio ainda, logo sem rank nem particao: a fracao sai da forma
    // analitica e o laco cobre a malha base inteira.  E' O(ncel) uma vez so',
    // no arranque -- nao vale a envoltoria que o caminho dinamico usa.
    const long ncel = (long) nx * ny * nz;
    real *sem = (real *) malloc((size_t) ncel * 3 * sizeof(real));
    long nsem = 0;
    for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            Point xc, d;
            xc[0] = (i + 0.5) * hx; d[0] = hx;
            xc[1] = (j + 0.5) * hy; d[1] = hy;
            xc[2] = (k + 0.5) * hz; d[2] = hz;
            const real f = fracvol(xc, d, 0.0);
            if (f > 0.001 && f < 0.999) {
                sem[3*nsem+0] = xc[0];
                sem[3*nsem+1] = xc[1];
                sem[3*nsem+2] = xc[2];
                nsem++;
            }
        }

    signed char *tab = (signed char *) malloc((size_t) ncel);
    amr_tabela(sem, nsem, hx, hy, hz, nx, ny, nz, thr, niveis, tab);

    int meurank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &meurank);
    if (meurank == 0) {
        if (amr_escreve_arquivo(caminho, lx, ly, lz, nx, ny, nz, niveis, tab) != 0)
            print0f("VOF3 amr inicial: FALHEI ao escrever %s\n", caminho);
        else {
            long marc = 0;
            for (long q = 0; q < ncel; q++) if (tab[q] > 0) marc++;
            print0f("VOF3 amr inicial: %ld sementes de %ld celulas base, "
                    "%ld marcadas (%.2f%%), niveis=%d\n",
                    nsem, ncel, marc, 100.0*marc/(double)ncel, niveis);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(sem); free(tab); free(thr);
}
