#ifndef T8_STENCIL_SUPPORT_H
#define T8_STENCIL_SUPPORT_H

// Coleta do SUPORTE de estencil na floresta do t8code -- clausula C10.
//
// A DIVISAO DE TRABALHO IMPORTA E E' A DO CONTRATO.  C10 diz que interpolacao de
// ordem k reproduz polinomio de grau <= k.  Quem reproduz e' o ajuste de minimos
// quadrados moveis, que pertence a Discretization e e' o MESMO para qualquer
// malha (`wls_set_points_and_calc`, em higtree/src/wls.h).  O que pertence a Mesh
// -- e o que uma segunda implementacao precisa entregar -- e' o SUPORTE: quais
// celulas ficam perto do ponto de consulta.
//
// Entao este modulo nao interpola.  Ele responde "quais centros usar", e o teste
// passa esses centros ao mesmo `wls` que o MTree usa.  Se a reproducao falhar com
// o suporte do t8code e acertar com o do MTree, o defeito esta' na malha, que e'
// exatamente o que se quer poder afirmar.
//
// COMO O SUPORTE E' COLETADO: busca em largura a partir do elemento que contem o
// ponto, por VIZINHANCA DE FACE (`t8_forest_leaf_face_neighbors`).  Nao e' varredura
// de todas as folhas com filtro de caixa -- isso seria O(folhas) por consulta e
// nao usaria topologia nenhuma, servindo igual para qualquer backend.  A busca em
// largura exercita a primitiva de vizinhanca do t8code, inclusive atraves de salto
// de nivel, que e' onde ela tem historia de ser fragil (ver C11).

#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Centros das celulas do suporte de estencil em torno de `x`.
//!
//! A malha e' a mesma do localizador: [0,1]^DIM com 8 celulas por direcao,
//! dividida em `ntrees_x` arvores de cmesh.  Percorre vizinhos de face em
//! largura ate' juntar ao menos `minpts` centros ou esgotar `maxpts`.
//!
//! Devolve o numero de centros escritos em `pts`, ou 0 se `x` estiver fora.
int t8_monta_suporte(int ntrees_x, const Point x, int minpts, int maxpts,
                     Point pts[]);

#ifdef __cplusplus
}
#endif

#endif
