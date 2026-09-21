#ifndef T8_PARTICAO_GRAFO_H
#define T8_PARTICAO_GRAFO_H

// O `partition_graph` CONSTRUIDO A PARTIR DAS CAIXAS DO T8CODE.
//
// Ate' aqui o t8code produzia a malha e o `lbal` a reparticionava.  Aqui a
// particao FINAL e' a do t8code: cada rank fica com as celulas que a curva de
// preenchimento lhe deu, e o grafo de vizinhanca e' montado a partir das caixas,
// sem `lb_calc_partition`.
//
// O CONTRATO QUE ISTO TEM DE CUMPRIR, lido do pdomain.c e nao suposto:
//
//   quem envia   para cada `to_send[i]`, enumera os elementos de `local_tree` no
//                intervalo [lo_idx, hi_idx) e manda, concatenados na ordem de i.
//   quem recebe  para cada `to_recv_trees[k]`, enumera a arvore INTEIRA e consome
//                o buffer na ordem de k.
//
// Logo a i-esima FAIXA que eu envio tem de ter a mesma forma e a mesma ordem de
// enumeracao que a i-esima ARVORE que o vizinho recebe.  O cabecalho do `pg` diz
// "strictly on the same order as the corresponding list on the remote process", e
// ordem trocada nao quebra nada visivel -- entrega o valor de outra celula.  Por
// isso os dois lados calculam o MESMO objeto geometrico a partir dos MESMOS dados
// globais, no MESMO aninhamento de lacos: nada e' negociado por mensagem.
//
// LIMITE DECLARADO: malha UNIFORME.  Com refino, a faixa do remetente e a arvore
// do receptor precisariam ter a mesma subdivisao interna, e isso exige trocar a
// estrutura -- nao so' os indices.  O teste cobre o caso uniforme, que e' o do
// example2d_Newt.

#include "higtree.h"
#include "domain.h"
#include "pdomain.h"

#ifdef __cplusplus
extern "C" {
#endif

//! \brief Monta `sd` e `pg` com a particao do t8code.  Coletivo.
//!
//! `nb` e' a grade uniforme (160x40, por exemplo).  A franja usa
//! `pg_get_fringe_size(pg)`, que o chamador ja' deve ter ajustado.
//!
//! Devolve 0 se nao conseguir -- e o motivo vai para stderr, porque falhar aqui
//! em silencio significaria dominio com vizinhanca errada.
int t8_monta_dominio_particionado(const Point lo, const Point hi, const int nb[DIM],
                                  sim_domain *sd, partition_graph *pg);

#ifdef __cplusplus
}
#endif
#endif
