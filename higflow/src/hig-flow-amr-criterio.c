// *******************************************************************
//  CRITERIO DE REFINO ADAPTATIVO: limiares, tabela de niveis, arquivo .amr
// *******************************************************************
//
// POR QUE ESTE MODULO E' SEPARADO DOS DOIS METODOS.
//
// Nada aqui sabe de front-tracking nem de VOF: as funcoes recebem uma LISTA DE
// SEMENTES e devolvem a tabela de niveis e o arquivo.  O que e' especifico de
// cada metodo e' apenas DE ONDE as sementes vem -- da geometria da superficie,
// no front-tracking, ou do campo de fracao advectado, no VOF.
//
// A separacao tem uma consequencia que importa para a comparacao entre os dois:
// como os dois usam a MESMA construcao de tabela e o MESMO escritor, qualquer
// diferenca entre as malhas que eles produzem vem das SEMENTES, isto e', da
// fracao de cada metodo -- e nao de dois criterios escritos em lugares
// diferentes, que divergiriam com o tempo.
//
// E nao depende do solver: nem PETSc, nem MPI, nem higflow_solver.  E' o que
// permite o arreio em tests/front-tracking/teste-amr3d.c compila-lo direto e
// exigir igualdade exata entre a tabela por dilatacao e a exaustiva.
//
// Este arquivo nasceu dentro do modulo de front-tracking 3D, e saiu de la' no
// momento em que o VOF passou a precisar do mesmo nucleo: deixa-lo onde estava
// faria o VOF depender do front-tracking, o que a arvore mantem separado de
// proposito.

#include "hig-flow-amr-criterio.h"
#include <math.h>
#include <string.h>

void amr_limiares(real h_base, int niveis, int cel_min, real *thr)
{
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    // Um h_base INTEIRO de margem, pelo mesmo motivo do 2D: o criterio decide
    // por celula base, entao a borda da banda e' incerta em h_base/2 pelo teste
    // no centro E a celula se estende outro h_base/2 alem.
    const real margem = h_base + h_fino;
    const real banda_fina = cel_min * h_fino + margem;
    for (int l = niveis - 1; l >= 0; l--) {
        thr[l] = banda_fina;
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        for (int k = 0; k < niveis - 1 - l; k++) thr[l] += cel_min * h_nivel;
    }
}

//! A tabela de niveis por DILATACAO das sementes.  O caminho de producao.
void amr_tabela(const real *sem, long nsem,
                               real hx, real hy, real hz,
                               int nx, int ny, int nz,
                               const real *thr, int niveis,
                               signed char *tab)
{
    const long ncel = (long) nx * ny * nz;
    for (long q = 0; q < ncel; q++) tab[q] = 0;
    if (nsem <= 0) return;

    // Raio de carimbo, por direcao, do MAIOR limiar.
    const int Rx = (int) (thr[0] / hx) + 1;
    const int Ry = (int) (thr[0] / hy) + 1;
    const int Rz = (int) (thr[0] / hz) + 1;

    for (long k = 0; k < nsem; k++) {
        const real sx = sem[3*k+0], sy = sem[3*k+1], sz = sem[3*k+2];
        const int ci = (int) (sx / hx), cj = (int) (sy / hy), ck = (int) (sz / hz);
        int i0 = ci - Rx, i1 = ci + Rx;
        int j0 = cj - Ry, j1 = cj + Ry;
        int k0 = ck - Rz, k1 = ck + Rz;
        if (i0 < 0) i0 = 0;
        if (i1 >= nx) i1 = nx - 1;
        if (j0 < 0) j0 = 0;
        if (j1 >= ny) j1 = ny - 1;
        if (k0 < 0) k0 = 0;
        if (k1 >= nz) k1 = nz - 1;
        for (int kk = k0; kk <= k1; kk++) {
            const real dz = (kk + 0.5) * hz - sz;
            for (int jj = j0; jj <= j1; jj++) {
                const real dy = (jj + 0.5) * hy - sy;
                for (int ii = i0; ii <= i1; ii++) {
                    const real dx = (ii + 0.5) * hx - sx;
                    const real d2 = dx*dx + dy*dy + dz*dz;
                    // o nivel mais FINO cuja banda contem o ponto
                    int nivel = 0;
                    for (int l = niveis - 1; l >= 0; l--)
                        if (d2 <= thr[l]*thr[l]) { nivel = l + 1; break; }
                    if (nivel == 0) continue;
                    const long q = ii + (long) nx * (jj + (long) ny * kk);
                    if ((signed char) nivel > tab[q]) tab[q] = (signed char) nivel;
                }
            }
        }
    }
}

//! A MESMA tabela por busca exaustiva.  Existe para o arreio: o contrato da
//! dilatacao e' igualdade EXATA com isto, e sem a segunda implementacao a
//! aceleracao viraria uma aproximacao silenciosa.  Nao usar em producao --
//! e' O(ncel x nsem).
void amr_tabela_exaustiva(const real *sem, long nsem,
                                         real hx, real hy, real hz,
                                         int nx, int ny, int nz,
                                         const real *thr, int niveis,
                                         signed char *tab)
{
    for (int kk = 0; kk < nz; kk++) {
        const real cz = (kk + 0.5) * hz;
        for (int jj = 0; jj < ny; jj++) {
            const real cy = (jj + 0.5) * hy;
            for (int ii = 0; ii < nx; ii++) {
                const real cx = (ii + 0.5) * hx;
                real d2 = 1e300;
                for (long k = 0; k < nsem; k++) {
                    const real a = cx - sem[3*k+0];
                    const real b = cy - sem[3*k+1];
                    const real c = cz - sem[3*k+2];
                    const real s = a*a + b*b + c*c;
                    if (s < d2) d2 = s;
                }
                int nivel = 0;
                for (int l = niveis - 1; l >= 0; l--)
                    if (d2 <= thr[l]*thr[l]) { nivel = l + 1; break; }
                tab[ii + (long) nx * (jj + (long) ny * kk)] = (signed char) nivel;
            }
        }
    }
}

// Escrita do .amr multinivel.  Separada e INDEPENDENTE DO SOLVER para o arreio
// poder reler o arquivo e conferir as POSICOES, nao so' as contagens -- em 2D a
// primeira versao escreveu indice base 0 e a contagem de celulas finas saiu
// CERTA com as posicoes ERRADAS, deslocadas de uma celula do nivel, e a banda
// medida caiu de 5 para 1,5 celulas.  Contagem certa com posicao errada e' o
// modo de falha que nao aparece em nenhum total.
//
// Formato conferido no parser higio_read_amr_info, que e' generico em DIM:
//   l[d] h[d] para cada dimensao
//   numero de niveis
//   por nivel: delta[d], numero de retalhos, e initcell[d] patchsize[d]
// INDICE BASE 1.
int amr_escreve_arquivo(const char *caminho, real lx, real ly, real lz,
                            int nx, int ny, int nz, int niveis,
                            const signed char *tab)
{
    const real hx = lx / nx, hy = ly / ny, hz = lz / nz;
    const long ncel = (long) nx * ny * nz;
    FILE *f = fopen(caminho, "w");
    if (f == NULL) return 1;

    long cont[16] = {0};
    for (long q = 0; q < ncel; q++)
        for (int l = 1; l <= niveis; l++) if (tab[q] >= l) cont[l]++;
    int niv_esc = 1;
    for (int l = 1; l <= niveis; l++) if (cont[l] > 0) niv_esc = l + 1;

    fprintf(f, "0.0 %.10g 0.0 %.10g 0.0 %.10g\n",
            (double) lx, (double) ly, (double) lz);
    fprintf(f, "%d\n", niv_esc);
    fprintf(f, "%.10g %.10g %.10g\n1\n1 1 1 %d %d %d\n",
            (double) hx, (double) hy, (double) hz, nx, ny, nz);
    for (int l = 1; l < niv_esc; l++) {
        const int fator = 1 << l;
        fprintf(f, "%.10g %.10g %.10g\n%ld\n",
                (double) (hx/fator), (double) (hy/fator), (double) (hz/fator),
                cont[l]);
        for (int kk = 0; kk < nz; kk++)
          for (int jj = 0; jj < ny; jj++)
            for (int ii = 0; ii < nx; ii++)
                if (tab[ii + (long) nx * (jj + (long) ny * kk)] >= l)
                    fprintf(f, "%d %d %d %d %d %d\n",
                            fator*ii + 1, fator*jj + 1, fator*kk + 1,
                            fator, fator, fator);
    }
    fclose(f);
    return 0;
}
