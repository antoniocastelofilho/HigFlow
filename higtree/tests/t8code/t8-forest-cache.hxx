#ifndef T8_FOREST_CACHE_HXX
#define T8_FOREST_CACHE_HXX

// Floresta de teste, comum aos modulos do t8code.
//
// O localizador e o coletor de suporte trabalham sobre a MESMA floresta, de
// proposito: se cada um montasse a sua, compara-los nao diria nada sobre a
// malha, so' sobre duas construcoes independentes.  Cabecalho C++ -- o t8code
// exige C++20 e nenhum destes tipos atravessa para o lado C.

#include <t8_forest/t8_forest_general.h>

//! Maior numero de arvores de cmesh por direcao que o cache aceita.
#define T8_MAX_TREES 4

//! Constroi (ou reaproveita) [0,1]^DIM com 8 celulas por direcao, repartido em
//! `ntrees_x` arvores de cmesh em cada direcao.
t8_forest_t t8_floresta_de_teste(int ntrees_x);

//! sc_init/t8_init, uma vez por processo.
void t8_inicializa_uma_vez(void);

//! Libera as florestas em cache.
void t8_libera_florestas(void);

#endif
