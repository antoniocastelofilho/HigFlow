#ifndef T8_POINT_LOCATOR_H
#define T8_POINT_LOCATOR_H

// Localizacao por ponto apoiada no t8code, para as clausulas C7, C8 e C9.
//
// AQUI O T8CODE RESPONDE, nao e' traduzido.  O produtor de malha
// (t8-mesh-producer.h) materializa a floresta como arvore hig e deixa o MTree
// responder tudo -- util para a C11, inutil para julgar as consultas.  Este
// modulo nao constroi arvore hig nenhuma: recebe um ponto e devolve o centro da
// celula que o contem, decidido dentro do t8code.
//
// O OBSTACULO MEDIDO, e por que este modulo existe.  A documentacao do
// `t8_forest_element_points_inside` diz, textualmente, que o retorno e' true
// tambem quando o ponto esta' na FRONTEIRA do elemento, e portanto "may return
// true for different leaf elements, if they are neighbors and the point lies on
// the common boundary".  A primitiva NAO DESEMPATA.
//
// A clausula C8 exige resposta definida: ponto sobre face interna vai para a
// celula de MENOR coordenada -- intervalo fechado em cima, aberto embaixo.  Essa
// convencao nao se herda do t8code; e' imposta aqui, sobre a primitiva dele, e e'
// o unico ponto do modulo onde ha' decisao de projeto em vez de consulta.

#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Centro da celula que contem `p`, decidido no t8code.
//!
//! A malha e' [0,1]^DIM com 8 celulas por direcao, dividida em `ntrees_x`
//! arvores de cmesh na direcao x -- o analogo de "uma arvore contra duas" que a
//! clausula C9 exige.  A floresta e' construida uma vez por `ntrees_x` e
//! reaproveitada.
//!
//! Devolve 1 e preenche `centro` se achou; 0 se o ponto esta' fora.
int t8_localiza_ponto(int ntrees_x, const Point p, Point centro);

//! \brief Percorre as folhas ao contrario.  So' para teste.
//!
//! A ordem natural da curva de preenchimento visita a celula de menor coordenada
//! primeiro, entao um desempate que apenas guardasse o PRIMEIRO candidato
//! passaria por sorte de ordem e nao por regra -- medido: removendo a regra, os
//! tres casos continuavam verdes.  Invertendo o percurso, a regra aparece.
void t8_localizador_inverte_percurso(int inverte);

//! \brief Libera as florestas em cache.  Chamar ao fim do teste.
void t8_localizador_encerra(void);

#ifdef __cplusplus
}
#endif

#endif
