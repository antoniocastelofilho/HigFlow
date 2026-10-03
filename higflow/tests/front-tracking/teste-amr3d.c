// Arreio do CRITERIO DE REFINO 3D.  Standalone: as funcoes de criterio vivem no
// modulo FT justamente por nao tocarem o solver.
//
// Quatro oraculos, e a ordem importa -- cada um so' significa algo se o anterior
// passou:
//
//   1. LIMIARES conferidos contra a regra recalculada aqui, de forma
//      independente.  Se os limiares estiverem errados, nada adiante significa.
//   2. DILATACAO igual a BUSCA EXAUSTIVA, byte a byte, em varias configuracoes.
//      E' o contrato da aceleracao: ela restringe quais sementes sao testadas,
//      nao o que se calcula.  Sem este teste a dilatacao viraria uma
//      aproximacao silenciosa.
//   3. A BANDA de fato entrega `cel_min` celulas finas de cada lado da
//      interface.  E' o SENTIDO da regra, nao a sua aritmetica: mede-se, da
//      superficie, a distancia ate' a celula mais proxima que NAO esta' no nivel
//      mais fino.
//   4. O ARQUIVO .amr rele' com as POSICOES certas.  Em 2D a primeira versao
//      escreveu indice base 0: a contagem de celulas finas saiu CERTA e as
//      posicoes ERRADAS, e a banda medida caiu de 5 para 1,5 celulas.  Contagem
//      certa com posicao errada nao aparece em nenhum total, e por isso este
//      teste confere retalho por retalho contra a tabela.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hig-flow-amr-criterio.h"

static int falhas = 0;
static void ok(const char *nome, int passou, const char *det) {
    printf("  [%s] %-46s %s\n", passou ? "ok" : "FALHOU", nome, det ? det : "");
    if (!passou) falhas++;
}

int main(void) {
    printf("\n=== criterio de refino 3D ===\n\n");

    // ---- 1. os limiares ------------------------------------------------
    printf("-- limiares da regra das celulas minimas --\n");
    {
        int pior_n = -1; real pior = 0.0;
        for (int niveis = 1; niveis <= 3; niveis++) {
            const real h = 1.0/50.0; const int cm = 5;
            real thr[4];
            amr_limiares(h, niveis, cm, thr);
            // recalculo independente da regra
            real hf = h; for (int l = 0; l < niveis; l++) hf *= 0.5;
            real banda = cm*hf + (h + hf);
            for (int l = niveis-1; l >= 0; l--) {
                real esp = banda;
                real hn = h; for (int k = 0; k < l+1; k++) hn *= 0.5;
                for (int k = 0; k < niveis-1-l; k++) esp += cm*hn;
                real e = fabs(thr[l] - esp);
                if (e > pior) { pior = e; pior_n = niveis; }
            }
            // e sao DECRESCENTES
            for (int l = 1; l < niveis; l++)
                if (thr[l] >= thr[l-1]) { pior = 1.0; pior_n = niveis; }
        }
        char d[120];
        snprintf(d, sizeof d, "pior desvio %.3e (niveis=%d)", pior, pior_n);
        ok("limiares batem com a regra, e decrescem", pior < 1e-15, d);
    }

    // ---- 2. dilatacao == exaustiva -------------------------------------
    printf("\n-- dilatacao contra busca exaustiva: igualdade EXATA --\n");
    {
        // varias formas de nuvem de sementes, para nao testar so' a esfera
        struct { const char *nome; int nx, ny, nz, niveis, cm; } casos[] = {
            {"esfera centrada",        20, 20, 20, 1, 3},
            {"esfera centrada",        20, 20, 20, 2, 3},
            {"caixa achatada",         16, 16, 32, 2, 5},
            {"semente unica",          12, 12, 12, 3, 2},
            {"sementes na borda",      14, 14, 14, 2, 4},
        };
        int piores = 0; char det[160] = "";
        for (int c = 0; c < 5; c++) {
            const int nx = casos[c].nx, ny = casos[c].ny, nz = casos[c].nz;
            const real hx = 1.0/nx, hy = 1.0/ny, hz = (c==2 ? 2.0 : 1.0)/nz;
            real h = hx; if (hy<h) h=hy; if (hz<h) h=hz;
            real thr[4];
            amr_limiares(h, casos[c].niveis, casos[c].cm, thr);

            // monta as sementes
            long nsem = 0;
            real *sem = (real*) malloc((size_t)nx*ny*nz*3*sizeof(real));
            if (c == 3) {                         // uma semente so'
                sem[0]=0.5; sem[1]=0.5; sem[2]=0.5; nsem=1;
            } else if (c == 4) {                  // sementes rentes a borda
                for (int q = 0; q < 6; q++) {
                    sem[3*nsem+0] = (q==0?0.5*hx:(q==1?1.0-0.5*hx:0.5));
                    sem[3*nsem+1] = (q==2?0.5*hy:(q==3?1.0-0.5*hy:0.5));
                    sem[3*nsem+2] = (q==4?0.5*hz:(q==5?1.0-0.5*hz:0.5));
                    nsem++;
                }
            } else {                              // casca de esfera
                const real R = 0.25;
                for (int kk=0;kk<nz;kk++) for (int jj=0;jj<ny;jj++) for (int ii=0;ii<nx;ii++) {
                    real x=(ii+0.5)*hx-0.5, y=(jj+0.5)*hy-0.5;
                    real z=(kk+0.5)*hz-(c==2?1.0:0.5);
                    real r=sqrt(x*x+y*y+z*z);
                    if (fabs(r-R) < 0.6*h) {
                        sem[3*nsem+0]=(ii+0.5)*hx; sem[3*nsem+1]=(jj+0.5)*hy;
                        sem[3*nsem+2]=(kk+0.5)*hz; nsem++;
                    }
                }
            }

            const long ncel = (long)nx*ny*nz;
            signed char *a = (signed char*) malloc((size_t)ncel);
            signed char *b = (signed char*) malloc((size_t)ncel);
            amr_tabela          (sem,nsem,hx,hy,hz,nx,ny,nz,thr,casos[c].niveis,a);
            amr_tabela_exaustiva(sem,nsem,hx,hy,hz,nx,ny,nz,thr,casos[c].niveis,b);
            long dif = 0;
            for (long q=0;q<ncel;q++) if (a[q]!=b[q]) dif++;
            char nm[120];
            snprintf(nm, sizeof nm, "%s niveis=%d cm=%d (%ld sementes)",
                     casos[c].nome, casos[c].niveis, casos[c].cm, nsem);
            snprintf(det, sizeof det, "%ld de %ld celulas diferem", dif, ncel);
            ok(nm, dif == 0, det);
            if (dif) piores++;
            free(a); free(b); free(sem);
        }
        (void) piores;
    }

    // ---- 3. a banda entrega cel_min celulas finas por lado -------------
    printf("\n-- a banda entrega o piso de celulas finas --\n");
    {
        const int nx=40, ny=40, nz=40, niveis=2, cm=5;
        const real hx=1.0/nx, hy=1.0/ny, hz=1.0/nz, R=0.25;
        real h=hx; real thr[4];
        amr_limiares(h,niveis,cm,thr);
        long nsem=0; real *sem=(real*)malloc((size_t)nx*ny*nz*3*sizeof(real));
        for (int kk=0;kk<nz;kk++) for (int jj=0;jj<ny;jj++) for (int ii=0;ii<nx;ii++) {
            real x=(ii+0.5)*hx-0.5,y=(jj+0.5)*hy-0.5,z=(kk+0.5)*hz-0.5;
            real r=sqrt(x*x+y*y+z*z);
            if (fabs(r-R) < 0.6*h) {
                sem[3*nsem+0]=(ii+0.5)*hx;sem[3*nsem+1]=(jj+0.5)*hy;
                sem[3*nsem+2]=(kk+0.5)*hz;nsem++;
            }
        }
        signed char *tab=(signed char*)malloc((size_t)nx*ny*nz);
        amr_tabela(sem,nsem,hx,hy,hz,nx,ny,nz,thr,niveis,tab);
        // distancia da SUPERFICIE (r = R) ate' a celula mais proxima que NAO
        // esta' no nivel mais fino, em celulas FINAS
        real h_fino = h; for (int l=0;l<niveis;l++) h_fino*=0.5;
        real dmin = 1e300;
        for (int kk=0;kk<nz;kk++) for (int jj=0;jj<ny;jj++) for (int ii=0;ii<nx;ii++) {
            if (tab[ii+(long)nx*(jj+(long)ny*kk)] >= niveis) continue;
            real x=(ii+0.5)*hx-0.5,y=(jj+0.5)*hy-0.5,z=(kk+0.5)*hz-0.5;
            real d=fabs(sqrt(x*x+y*y+z*z)-R);
            if (d<dmin) dmin=d;
        }
        char det[160];
        snprintf(det,sizeof det,"%.2f celulas finas (piso %d), %ld sementes",
                 dmin/h_fino, cm, nsem);
        ok("piso de celulas finas por lado respeitado", dmin/h_fino >= (real)cm, det);
        free(tab); free(sem);
    }

    // ---- 4. o arquivo .amr: POSICOES, nao so' contagens ----------------
    printf("\n-- o .amr rele' com as posicoes certas --\n");
    {
        const int nx=12, ny=12, nz=20, niveis=2, cm=2;
        const real lx=1.0, ly=1.0, lz=2.0;
        const real hx=lx/nx, hy=ly/ny, hz=lz/nz;
        real h=hx; if(hy<h)h=hy; if(hz<h)h=hz;
        real thr[4]; amr_limiares(h,niveis,cm,thr);
        real sem[3]={0.5,0.5,1.0}; 
        const long ncel=(long)nx*ny*nz;
        signed char *tab=(signed char*)malloc((size_t)ncel);
        amr_tabela(sem,1,hx,hy,hz,nx,ny,nz,thr,niveis,tab);

        const char *cam="/tmp/teste-amr3d.amr";
        int e = amr_escreve_arquivo(cam,lx,ly,lz,nx,ny,nz,niveis,tab);
        ok("arquivo escrito", e==0, cam);

        // rele' com o MESMO formato que o parser da higtree documenta
        FILE *f=fopen(cam,"r");
        double l0,h0,l1,h1,l2,h2; int nniv;
        int lidos = fscanf(f,"%lf %lf %lf %lf %lf %lf %d",&l0,&h0,&l1,&h1,&l2,&h2,&nniv);
        char det[200];
        snprintf(det,sizeof det,"caixa [%g,%g]x[%g,%g]x[%g,%g], %d niveis",l0,h0,l1,h1,l2,h2,nniv);
        ok("cabecalho: caixa e numero de niveis", lidos==7 && h0==lx && h1==ly && h2==lz, det);

        int mau_pos=0, mau_niv=0; long total_ret=0;
        for (int niv=0; niv<nniv; niv++) {
            double d0,d1,d2; long npat;
            if (fscanf(f,"%lf %lf %lf %ld",&d0,&d1,&d2,&npat)!=4) { mau_pos=1; break; }
            const int fator = 1<<niv;
            for (long p=0;p<npat;p++) {
                int ci,cj,ck,si,sj,sk;
                if (fscanf(f,"%d %d %d %d %d %d",&ci,&cj,&ck,&si,&sj,&sk)!=6) { mau_pos=1; break; }
                total_ret++;
                if (niv==0) continue;            // o nivel base e' um retalho so'
                // BASE 1: (ci-1)/fator tem de ser o indice da celula base, e
                // aquela celula tem de ter nivel >= niv na tabela.
                if ((ci-1)%fator || (cj-1)%fator || (ck-1)%fator) { mau_pos++; continue; }
                int ii=(ci-1)/fator, jj=(cj-1)/fator, kk=(ck-1)/fator;
                if (ii<0||ii>=nx||jj<0||jj>=ny||kk<0||kk>=nz) { mau_pos++; continue; }
                if (tab[ii+(long)nx*(jj+(long)ny*kk)] < niv) mau_niv++;
                if (si!=fator||sj!=fator||sk!=fator) mau_pos++;
            }
        }
        fclose(f);
        snprintf(det,sizeof det,"%ld retalhos; %d com posicao fora da rede, %d em celula de nivel menor",
                 total_ret, mau_pos, mau_niv);
        ok("cada retalho cai na celula base que a tabela refina", mau_pos==0 && mau_niv==0, det);
        free(tab);
    }

    // ---- 5. o GANHO, medido: quantas celulas o criterio poupa -----------
    // Nao e' passa-ou-falha de precisao: e' o numerador do ganho do refino, e
    // sem ele "o refino acelerou" nao tem denominador.  O oraculo e' fraco de
    // proposito (ganho > 1 e banda >= piso); o valor e' o que interessa.
    printf("\n-- ganho de celulas, e a largura da banda --\n");
    printf("     %-6s %-7s %12s %12s %8s %10s\n",
           "base", "niveis", "criterio", "unif. fina", "ganho", "banda");
    {
        int ruim = 0; real banda_pior = 1e300, ganho_pior = 1e300;
        for (int nb = 25; nb <= 50; nb *= 2) {
          for (int niveis = 1; niveis <= 3; niveis++) {
            const int nx=nb, ny=nb, nz=2*nb;
            const real hx=1.0/nx, hy=1.0/ny, hz=2.0/nz;
            real h=hx; if(hy<h)h=hy; if(hz<h)h=hz;
            real thr[4]; amr_limiares(h,niveis,5,thr);
            long nsem=0;
            real *sem=(real*)malloc((size_t)nx*(size_t)ny*(size_t)nz*3*sizeof(real));
            for (int kk=0;kk<nz;kk++) for (int jj=0;jj<ny;jj++) for (int ii=0;ii<nx;ii++) {
                real x=(ii+0.5)*hx-0.5,y=(jj+0.5)*hy-0.5,z=(kk+0.5)*hz-0.5;
                real r=sqrt(x*x+y*y+z*z);
                // meia diagonal da celula: o criterio real marca a celula que a
                // superficie CORTA, nao a que tem o centro perto dela
                if (fabs(r-0.25) < 0.87*h) {
                    sem[3*nsem+0]=(ii+0.5)*hx;sem[3*nsem+1]=(jj+0.5)*hy;
                    sem[3*nsem+2]=(kk+0.5)*hz;nsem++;
                }
            }
            const long ncel=(long)nx*(long)ny*(long)nz;
            signed char *tab=(signed char*)malloc((size_t)ncel);
            amr_tabela(sem,nsem,hx,hy,hz,nx,ny,nz,thr,niveis,tab);
            long adap=0;
            for (long q=0;q<ncel;q++) { long fa=1; for(int l=0;l<tab[q];l++) fa*=8; adap+=fa; }
            long fina=ncel; for(int l=0;l<niveis;l++) fina*=8;
            real h_fino=h; for(int l=0;l<niveis;l++) h_fino*=0.5;
            real dmin=1e300;
            for (int kk=0;kk<nz;kk++) for (int jj=0;jj<ny;jj++) for (int ii=0;ii<nx;ii++) {
                if (tab[ii+(long)nx*(jj+(long)ny*kk)]>=niveis) continue;
                real x=(ii+0.5)*hx-0.5,y=(jj+0.5)*hy-0.5,z=(kk+0.5)*hz-0.5;
                real d=fabs(sqrt(x*x+y*y+z*z)-0.25);
                if(d<dmin)dmin=d;
            }
            real ganho=(real)fina/adap, banda=dmin/h_fino;
            if (banda < banda_pior) banda_pior = banda;
            if (ganho < ganho_pior) ganho_pior = ganho;
            if (ganho <= 1.0 || banda < 5.0) ruim++;
            printf("     %-6d %-7d %12ld %12ld %7.1fx %9.1f\n",
                   nb,niveis,adap,fina,(double)ganho,(double)banda);
            free(tab); free(sem);
          }
        }
        char d[200];
        snprintf(d,sizeof d,"menor ganho %.1fx, banda mais estreita %.1f celulas finas",
                 (double)ganho_pior,(double)banda_pior);
        ok("ganho > 1 e banda >= piso em toda configuracao", ruim==0, d);
        printf("\n     NOTA: a banda CRESCE com os niveis (8 -> 12 -> 16 contra piso 5)\n");
        printf("     porque a margem e' fixa em unidades da celula BASE e portanto\n");
        printf("     vale 2^niveis celulas finas.  Conservador, e custa celulas.\n");
    }

    printf("\n%s  (%d falha(s))\n\n", falhas ? "FALHOU" : "PASSOU", falhas);
    return falhas ? 1 : 0;
}
