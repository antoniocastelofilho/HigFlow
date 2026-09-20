#ifndef T8_POINT_CLASS_H
#define T8_POINT_CLASS_H

// Classificacao de ponto: dentro, sobre o contorno, ou fora -- metade de Mesh na
// clausula C14.
//
// C14 e' sobre o ramo ON_BOUNDARY do fechamento, e esse ramo e' escolhido por um
// DESPACHO: `cell_find_in_center`, em domain.c, decide entre IN_DOMAIN_PROPER,
// ON_BOUNDARY e OUTSIDE_DOMAIN antes de qualquer interpolacao.  Metade do
// despacho de fechamento passa por ali, e foi nesse ramo que sobreviveram tres
// dos sete sitios do defeito do bc_inter -- verdes por ausencia de teste.
//
// A DECISAO E' DE MALHA: "onde este ponto esta' em relacao ao dominio" nao depende
// de discretizacao nenhuma.  Quem interpola depois e' Discretization.
//
// Os dois backends decidem por caminhos diferentes, como em C12: a HiGTree compara
// o ponto com a CAIXA da arvore, o t8code pergunta se o ponto cai numa face sem
// vizinho.  Sao respostas que podem divergir, e onde divergem importa.

#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

//! Os mesmos tres estados do `point_location` de domain.c, na mesma ordem.
typedef enum {
    T8_DENTRO       = 0,
    T8_NO_CONTORNO  = 1,
    T8_FORA         = 2
} t8_classe_de_ponto;

//! \brief Onde `p` esta' em relacao ao dominio, decidido no t8code.
//!
//! Malha [0,1]^DIM com 8 celulas por direcao, em `ntrees_x` arvores de cmesh por
//! direcao.  FORA se nenhum elemento o contem; NO_CONTORNO se cai sobre uma face
//! sem vizinho; DENTRO caso contrario.
t8_classe_de_ponto t8_classifica_ponto(int ntrees_x, const Point p);

#ifdef __cplusplus
}
#endif

#endif
