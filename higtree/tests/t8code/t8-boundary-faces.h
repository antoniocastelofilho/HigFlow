#ifndef T8_BOUNDARY_FACES_H
#define T8_BOUNDARY_FACES_H

// Enumeracao das faces de CONTORNO na floresta do t8code -- metade de Mesh nas
// clausulas C12 e C14.
//
// A DIFERENCA DE MODELO E' O PONTO.  Na HiGTree o contorno e' EXPLICITO: arvores
// `sim_boundary` registradas com `sd_add_boundary`.  A malha nao sabe onde
// termina; alguem lhe conta.  No t8code o contorno e' IMPLICITO: e' toda face de
// elemento cujo `t8_forest_leaf_face_neighbors` devolve zero vizinhos.
//
// Os dois mecanismos sao opostos, e por isso vale afirmar que produzem a MESMA
// geometria.  Um backend de malha que errasse aqui faria o fechamento de contorno
// projetar sobre paredes que nao existem, ou ignorar paredes que existem -- e o
// fechamento e' onde viveu, por anos, o defeito dos sete sitios do bc_inter.
//
// O QUE ESTE MODULO NAO FAZ: fechar o estencil.  A projecao sobre a parede, a
// escolha da parede atravessada e a interpolacao em DIM-1 pertencem a
// Discretization/Boundary e sao as mesmas para qualquer malha.  Aqui so' se
// responde ONDE o contorno esta'.

#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Centros e normais das faces de contorno da malha de teste.
//!
//! [0,1]^DIM com 8 celulas por direcao, em `ntrees_x` arvores de cmesh por
//! direcao.  Escreve ate' `maxn` faces.  Devolve quantas escreveu.
int t8_faces_de_contorno(int ntrees_x, int maxn, Point centros[], Point normais[]);

#ifdef __cplusplus
}
#endif

#endif
