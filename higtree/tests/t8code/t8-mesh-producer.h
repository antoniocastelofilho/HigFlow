#ifndef T8_MESH_PRODUCER_H
#define T8_MESH_PRODUCER_H

// Produtor de malha apoiado no t8code, na forma que o contrato de Mesh pede.
//
// O contrato (higtree/src/hig-mesh-contract.h) diz que uma segunda implementacao
// nao herda de uma classe: ela PRODUZ as estruturas que as consultas leem.  E' o
// que esta funcao faz -- constroi uma floresta no t8code e materializa dela uma
// arvore `hig_cell`, que o test-level-jump consome sem saber de onde veio.
//
// A fronteira e' `extern "C"` de proposito: o t8code exige C++20 e a HiGTree
// constroi com gnu++17, entao a unidade do produtor compila separada, no padrao
// dela, e so' esta assinatura atravessa.

#include "higtree.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Constroi, VIA T8CODE, uma malha nao graduada com salto 4:1, e a
//! materializa como arvore hig_cell.
//!
//! A geometria e' a mesma do produtor de referencia (MTree) para que a
//! comparacao seja uma comparacao: raiz [0,1]^DIM em nivel 2 (4 por direcao),
//! a celula que contem 0,375 refinada uma vez e a que contem 0,3125 refinada de
//! novo, com a vizinha imediata permanecendo grossa.
//!
//! Devolve NULL se o t8code nao produzir a malha esperada.
hig_cell *t8_produz_malha_nao_graduada(void);

#ifdef __cplusplus
}
#endif

#endif
