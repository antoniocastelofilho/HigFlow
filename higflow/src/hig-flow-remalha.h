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

//! \brief Preenche por INTERPOLACAO CONSERVATIVA o que o casamento exato nao
//!        achou.  Chamar DEPOIS de `rem_planta_*`, com o mesmo `achado`.
//!
//! A separacao entre as duas etapas e' deliberada: o transporte exato e'
//! verificavel bit a bit, e a interpolacao nao e'.  Misturar as duas tiraria o
//! unico oraculo forte que existe aqui.
//!
//! COBRE SALTO DE UM NIVEL, que e' o que uma adaptacao produz por passo:
//!
//!   celula nova mais FINA   recebe o valor da antiga que a continha.  Para
//!       campo constante por celula isso E' conservativo: a soma dos filhos
//!       vezes o volume de cada um da' o volume do pai vezes o valor do pai.
//!
//!   celula nova mais GROSSA recebe a media das antigas que ela contem.  Os
//!       2^DIM filhos tem volumes iguais, entao media simples E' media em
//!       volume.
//!
//! COMO O PAI E' ACHADO SEM CONHECER A MALHA ANTIGA.  O plantador nao tem a
//! arvore de origem, so' a colheita indexada por posicao.  Mas os 2^DIM
//! candidatos a pai sao calculaveis do proprio box -- centro +- h/2 em cada
//! direcao -- e NO MAXIMO UM existe: dois candidatos seriam celulas de mesmo
//! tamanho que se sobrepoem parcialmente, o que malha nenhuma admite.  Logo a
//! escolha e' unica sem precisar saber qual filho a celula e'.
//!
//! So' as posicoes que faltaram pagam as consultas extras.
//!
//! \return quantas posicoes continuaram sem valor -- salto de mais de um nivel,
//!         ou regiao que nao existia de forma alguma na malha antiga.
long rem_interpola_centro(rem_colheita *c, sim_domain *sd,
                          distributed_property *dp, char *achado);

//! \brief Idem, para facetas.  Ver `rem_interpola_centro`.
//!
//! A FACETA TEM UM CASO A MAIS QUE A CELULA.  Refinar uma celula cria facetas
//! que nao sao sub-facetas de faceta nenhuma: as do plano do MEIO, que no nivel
//! grosso era interior.  Essa recebe a media das duas facetas paralelas que
//! limitam a celula antiga, o que preserva o balanco de fluxo atraves dela.
//!
//! CONSERVA FLUXO, NAO DIVERGENCIA.  A soma de u*A e' preservada faceta a
//! faceta, mas o campo resultante NAO e' de divergencia nula na malha nova.
//! Depois de remalhar, chame `higflow_projecao_remalha` (hig-flow-step.h): ela
//! projeta com o MESMO par do passo de tempo.  Em malha uniforme remove a
//! divergencia ate' o solver linear; em malha refinada sobra o residuo de
//! interface do proprio par do solver -- ver o comentario dela.
long rem_interpola_faceta(rem_colheita *c, sim_facet_domain *sfd,
                          distributed_property *dp, char *achado);

//! \brief Quantas posicoes a ultima interpolacao preencheu como filha (refino)
//!        e como mae (engrossamento).  Para guardar contra teste vazio: se as
//!        duas derem zero, a interpolacao nao foi exercitada.
void rem_interpolou(long *por_refino, long *por_engrossamento);

//! \brief Quantas facetas foram preenchidas com media PARCIAL das filhas.
//!
//! Numa interface de refino a malha de origem pode nao carregar todas as
//! sub-facetas como graus de liberdade.  Ali a media do que existe e' a melhor
//! informacao disponivel: exata para campo constante, de primeira ordem no
//! geral.  Fica contado para que isso seja escolha visivel, e nao silencio.
long rem_engrossamento_parcial(void);

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
