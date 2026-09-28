// Refino adaptativo compartilhado pelos exemplos bifasicos (VOF e front-tracking).
// Ver malha-adaptativa.c para o criterio e para a regra das celulas minimas.
#ifndef MALHA_ADAPTATIVA_H
#define MALHA_ADAPTATIVA_H

#include "../src/hig-flow-kernel.h"

#ifdef __cplusplus
extern "C" {
#endif

//! Limiares de distancia derivados da REGRA DAS CELULAS MINIMAS: o nivel mais
//! fino cobre ao menos `cel_min` celulas finas de cada lado da interface.
//! `thr` deve caber `niveis+1` reais.  Derivado, nao digitado -- ver o .c.
void malha_adapt_limiares(real h_base, int niveis, int cel_min, real *thr);

//! MEDE na malha de verdade a banda fina, em unidades de h_fino: a menor
//! distancia de uma semente de interface a uma celula que NAO esta' no nivel
//! mais fino, ja' descontada a incerteza da semente.  E' o ORACULO da regra --
//! tem de dar >= cel_min.  Devolve -1 se nao houver o que medir.
real malha_adapt_mede_banda(higflow_solver *ns, int niveis, real h_base);

//! Reconstroi o dominio com a malha adaptada, com o MESMO criterio nos dois
//! metodos (o criterio le' fracvol, que o VOF advecta e o front-tracking
//! escreve da geometria da frente).
void malha_adapt_reconstroi(higflow_solver *ns, higflow_solver *ns2,
                            int myrank, int ntasks, int cache, int order_center,
                            real *limiares, const char *amr_base);

//! Refina a arvore do solver no lugar, pelo criterio.
void higflow_refine_tree_inplace(higflow_solver *ns, real *thresholds);

//! Reune em TODOS os ranks as sementes de interface locais de todos eles, e
//! devolve o vetor global (que o chamador libera com `free`).  `n_total` recebe
//! o numero de sementes reunidas.
//!
//! POR QUE E' NECESSARIA.  O criterio de refino percorre o dominio LOCAL: cada
//! rank ve' apenas o pedaco de interface que possui.  Escrever o .amr dessas
//! sementes -- ainda que so' o rank 0 escreva -- produz uma malha refinada
//! apenas em torno daquele pedaco, e a regra das celulas minimas fica violada
//! exatamente onde ninguem esta' olhando.  MEDIDO em np=2: o contador do
//! remalhamento acusou 265 posicoes sem valor no primeiro remalhamento e 2709 no
//! segundo, com 77% de erro em Vc.
//!
//! Em np=1 e' identidade (devolve copia das locais), entao o caminho serial
//! continua bit a bit o mesmo.
Point *malha_adapt_reune_sementes(const Point *locais, long n_locais,
                                  long *n_total);

#ifdef __cplusplus
}
#endif
#endif
