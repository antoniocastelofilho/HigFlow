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

//! \brief A malha de PRODUCAO: a arvore que o `lb_add_input_tree` vai
//! distribuir, com a estrutura decidida pelo t8code.
//!
//! `nivel_base` da' 2^nivel_base celulas por direcao; `refinos` diz quantas vezes
//! refinar em torno de `alvo`.  Devolve em `folhas_out` quantas folhas o t8code
//! produziu -- o numero vem DELE, e e' contra ele que a contagem do dominio e'
//! conferida depois.
//!
//! A caixa e' [0,1]^DIM.  Limite declarado: o hipercubo do t8code e' unitario, e
//! mapear para outra caixa exige escalar o criterio de parada do refinamento.
hig_cell *t8_produz_malha_para_dominio(int nivel_base, double alvo, int refinos,
                                       long *folhas_out);

//! \brief Preenche um instantaneo DIRETO da floresta do t8code.
//!
//! E' aqui que a fronteira se move.  A funcao acima materializa a floresta como
//! octree de ponteiros -- O(folhas) por producao -- so' para que as consultas
//! antigas possam navega-la.  Esta nao constroi arvore nenhuma: le as folhas do
//! t8code e escreve centro e tamanho nos arranjos, que e' tudo o que 86% das
//! chamadas em laco quente precisam.
//!
//! Devolve NULL se a floresta nao tiver o salto 4:1 esperado.
struct hig_mesh_snapshot *t8_preenche_instantaneo(void);

#ifdef __cplusplus
}
#endif

#endif
