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

#ifdef __cplusplus
}
#endif
#endif
