// *******************************************************************
//  Marching cubes -- extracao da isosuperficie fracvol = 1/2
// *******************************************************************
//
// POR QUE ESTE ARQUIVO E' INDEPENDENTE DO HIGFLOW.
//
// O nucleo nao conhece higflow_solver, nem PETSc, nem distributed_property.
// Ele recebe uma GRADE REGULAR de nos e um AMOSTRADOR, e devolve triangulos
// com vertices SOLDADOS (cada aresta cortada gera um vertice so', partilhado
// por todos os triangulos que a usam).  Duas consequencias:
//
//   1. a superficie e' CONTINUA -- e' isso que o PLIC nao da', porque o PLIC
//      reconstroi um plano por celula sem casar os planos vizinhos;
//   2. o nucleo e' verificavel contra um campo ANALITICO, sem solver.
//      Ver tests/marching-cubes/.
//
// MALHA ADAPTATIVA.  A grade e' a do NIVEL MAIS FINO.  Isso basta porque
// a interface esta' sempre toda no nivel mais fino -- o criterio de banda
// garante cinco celulas finas por lado -- e fora dela a fracao e' 0 ou 1
// exatos, onde o marching cubes nao gera triangulo nenhum.  Logo o estencil
// nunca alcanca a costura entre niveis, que e' o que normalmente quebra
// extracao de isosuperficie em malha adaptativa.
//
// A GRADE E' DUAL.  A fracao de volume vive no CENTRO da celula, nao no
// vertice.  Os nos da grade do marching cubes sao os centros das celulas:
// o no (i,j,k) fica em x0 + (i,j,k)*h, com x0 = origem + h/2.  O cubo
// (i,j,k) usa os oito centros (i..i+1, j..j+1, k..k+1).  Daqui sai a
// superficie encolhida de meia celula em cada borda do dominio, o que e'
// correto: nao ha' dado de fracao fora da ultima celula para interpolar.

#ifndef HIG_FLOW_MC_H
#define HIG_FLOW_MC_H

#include <stdlib.h>

// Malha triangular de saida.  Cresce sozinha; liberar com mc_malha_libera.
typedef struct {
    int      nv, nt;            // vertices e triangulos usados
    int      cv, ct;            // capacidade alocada
    double  *v;                 // 3*cv  coordenadas
    int     *t;                 // 3*ct  indices de vertice
    double  *vn;                // 3*cv  normal por vertice (gradiente do campo)
} mc_malha;

// Amostrador: valor do campo no no' (i,j,k) da grade.  Chamado uma vez por
// no' por fatia; o nucleo nao presume continuidade nem suavidade.
typedef double (*mc_amostra)(int i, int j, int k, void *ctx);

void mc_malha_zera(mc_malha *m);
void mc_malha_libera(mc_malha *m);

// Extrai a isosuperficie campo = iso.
// ni,nj,nk: numero de NOS em cada direcao (nao de celulas).
// Devolve 0 em caso de sucesso, != 0 se faltou memoria.
int mc_extrai(int ni, int nj, int nk,
              const double x0[3], const double h[3],
              double iso, mc_amostra f, void *ctx, mc_malha *m);

// --- oraculos, usados pelos testes e pelo codigo de producao -------------
//
// Caracteristica de Euler V - A + F.  Para uma superficie fechada de genero
// zero vale 2.  Serve como teste de FECHAMENTO: se faltar um triangulo por
// erro de tabela, a aresta orfa derruba este numero.
int    mc_euler(const mc_malha *m);
// Numero de arestas usadas por um numero de triangulos != 2.  Zero numa
// superficie fechada.  Mais direto que Euler para localizar o defeito.
int    mc_arestas_orfas(const mc_malha *m);
// Volume encerrado, pelo teorema da divergencia sobre os triangulos.
double mc_volume(const mc_malha *m);
double mc_area(const mc_malha *m);

// Verifica a tabela de triangulos contra o que a geometria do cubo exige:
// o conjunto de arestas que a tabela usa no caso c tem de ser EXATAMENTE
// o conjunto de arestas cujos extremos trocam de sinal em c.  Devolve o
// numero de casos inconsistentes (zero se a tabela esta' boa).  Um numero
// de aresta trocado na tabela quase certamente viola isto.
int    mc_verifica_tabela(void);

#endif
