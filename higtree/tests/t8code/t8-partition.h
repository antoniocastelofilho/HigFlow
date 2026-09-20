#ifndef T8_PARTITION_H
#define T8_PARTITION_H

// O t8code sob particionamento REAL -- clausulas C13 e P1 a P4.
//
// As sete clausulas ja' verificadas para o t8code rodam em um processo.  Esta
// familia e' a que exige a floresta distribuida entre ranks, com camada de ghost,
// e e' onde a substituicao do particionador de fato acontece: o MTree particiona
// pelo `lbal.c` com Zoltan, o t8code por curva de preenchimento de espaco.
//
// UM LIMITE REGISTRADO NA C12 CAI AQUI.  "Face sem vizinho e' contorno do
// dominio" so' vale em um processo: sob particionamento, face sem vizinho LOCAL
// pode ser face de franja.  Com a camada de ghost pedida, o
// `t8_forest_leaf_face_neighbors` passa a enxergar o vizinho remoto, e a
// distincao volta a ser possivel.
//
// Todas as respostas aqui sao LOCAIS ao rank.  Quem reduz globalmente e' o teste:
// um teste paralelo em que cada rank conclui sozinho esconde a falha de um rank
// no meio das linhas dos outros.

#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Monta (uma vez) a floresta distribuida com ghost, nivel 3 -> 8 celulas
//! por direcao no total.  Deve ser chamada depois de MPI_Init.
void t8_particao_monta(void);

//! Numero de elementos LOCAIS a este rank.
long t8_particao_num_locais(void);

//! Soma dos volumes dos elementos locais.
double t8_particao_volume_local(void);

//! Numero de elementos de GHOST neste rank (0 se nao houver vizinho).
long t8_particao_num_ghosts(void);

//! 1 se existe elemento local cujo suporte de face alcanca um elemento de
//! ghost; 0 caso contrario.  E' o criterio da C13: o que a franja entrega.
long t8_particao_suporte_alcanca_ghost(void);

//! \brief Interpola o campo linear 1 + 2x + 3y no ponto `p`, se ele pertencer a
//! este rank.  Devolve 1 e escreve em `valor`; devolve 0 se o ponto nao e' deste
//! rank.  Usado pela P4: o valor nao pode depender da particao.
int t8_particao_interpola(const Point p, double *valor);

void t8_particao_encerra(void);

#ifdef __cplusplus
}
#endif

#endif
