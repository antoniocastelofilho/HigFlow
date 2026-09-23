// =============================================================================
// TRANSFERENCIA DE SOLUCAO ENTRE MALHAS
//
// Quando a malha e' refeita no meio da corrida, a particao muda: uma celula que
// era do rank 2 passa a ser do rank 5.  O campo tem de ACOMPANHAR.  Este modulo
// faz so' isso, e nao conhece o solver -- recebe dominios e propriedades.
//
// A CHAVE E' A POSICAO, NAO O ID LOCAL.  Id local nao sobrevive a uma
// reparticao: e' um indice num mapeador que foi jogado fora.  O que sobrevive e'
// onde a celula esta'.
//
// COMO A BUSCA E' DISTRIBUIDA.  Nem o dono antigo sabe quem sera' o dono novo,
// nem o novo sabe quem era o antigo.  Descobrir isso com uma varredura global
// custaria O(N) de memoria por rank.  Em vez disso, cada posicao tem um rank
// ANFITRIAO, calculado da propria posicao por uma funcao que todos os ranks
// computam igual:
//
//     colher:  o dono antigo manda (posicao, valor) ao anfitriao
//     plantar: o dono novo pede a posicao ao anfitriao, que responde
//
// Duas trocas todos-para-todos, O(N/P) de memoria.  E' o encontro marcado
// ("rendezvous") classico, e a unica coisa que os dois lados precisam combinar
// e' a funcao que leva posicao em anfitriao.
//
// LIMITE, E ELE E' DURO: a busca e' EXATA.  So' acha a posicao se a entidade
// existe nas DUAS malhas, com a mesma coordenada.  Isso cobre o caso para o qual
// este modulo foi escrito -- malha refeita, particao nova, geometria igual -- e
// NAO cobre celula genuinamente nova, que nao existia antes.  Para essa e'
// preciso INTERPOLAR, e interpolar conservativamente nao e' o mesmo problema.
// Quem pedir uma posicao que nao existe recebe isso contado, nao mascarado: ver
// o retorno de `rem_planta_*`.
// =============================================================================

#ifndef HIG_FLOW_REMALHA_H
#define HIG_FLOW_REMALHA_H

#include "domain.h"
#include "pdomain.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rem_colheita rem_colheita;

//! \brief Colhe o valor das celulas PROPRIAS, indexado pela posicao.
//!
//! Coletivo.  Depois disto o dominio e a propriedade podem ser destruidos: a
//! colheita nao guarda ponteiro para nenhum dos dois.
rem_colheita *rem_colhe_centro(sim_domain *sd, distributed_property *dp);

//! \brief Idem, para as facetas proprias de um dominio de facetas.
rem_colheita *rem_colhe_faceta(sim_facet_domain *sfd, distributed_property *dp);

//! \brief Planta a colheita na malha NOVA, casando por posicao.
//!
//! Coletivo.  Escreve so' nas entidades PROPRIAS; a franja fica por conta do
//! `dp_sync` do chamador, como em qualquer outro preenchimento.
//!
//! \return quantas posicoes proprias NAO foram achadas na colheita.  Zero e' o
//!         que se espera quando a geometria nao mudou.  Diferente de zero
//!         significa celula nova, que este modulo nao sabe preencher -- e o
//!         chamador precisa saber disso em vez de receber zero calado.
//! \param achado opcional (pode ser NULL): arranjo de pelo menos `local_count`
//!        bytes, onde entra 1 para o id local que foi achado e 0 para o que nao
//!        foi.  QUEM VAI INTERPOLAR PRECISA DA LISTA, nao da contagem: uma
//!        celula nova sem valor e' um buraco no campo, e so' o chamador sabe com
//!        o que preenche-lo.
long rem_planta_centro(rem_colheita *c, sim_domain *sd, distributed_property *dp,
                       char *achado);

//! \brief Idem, para facetas.  Ver `rem_planta_centro`.
long rem_planta_faceta(rem_colheita *c, sim_facet_domain *sfd,
                       distributed_property *dp, char *achado);

void rem_destroi(rem_colheita *c);

//! \brief Quantos valores da ultima plantacao vieram de um rank DIFERENTE do que
//!        os possuia antes.
//!
//! Existe para guardar contra teste vazio.  Se a particao nova sair igual a'
//! antiga, nada se move, e uma transferencia quebrada passa no oraculo de
//! identidade sem nunca ter transferido nada.  Um teste que nao afirma este
//! numero nao afirma a transferencia.
long rem_vieram_de_outro_rank(void);

//! \brief Quantas posicoes distintas a ultima colheita registrou no anfitriao.
long rem_posicoes_hospedadas(void);

//! \brief Quantas colisoes de chave a ultima colheita viu -- duas entidades
//!        distintas na mesma posicao quantizada.  Deve ser zero; diferente de
//!        zero e' defeito da quantizacao, nao do escoamento.
long rem_colisoes(void);

#ifdef __cplusplus
}
#endif

#endif
