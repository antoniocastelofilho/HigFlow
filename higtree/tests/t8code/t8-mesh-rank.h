#ifndef T8_MESH_RANK_H
#define T8_MESH_RANK_H

// PRODUCAO POR RANK, em CAIXAS COMPLETAS.
//
// A tentativa anterior representava a regiao de um rank como UMA arvore com
// buracos.  Nao funciona: `hig_get_cell_coords_of_point` desreferencia os filhos
// para ler as caixas deles, e `hig_get_cell_with_point` ainda faz
// `cell = cell->children[p]` sem testar nulo -- toda localizacao por ponto morre
// em SEGV.  O `lbal` nunca produz buraco: ele usa `hig_refine_empty` para ALOCAR
// o vetor e em seguida preenche TODOS os compartimentos.  Por isso um dominio tem
// VARIAS arvores por rank, e nao uma.
//
// Aqui a regiao vira um CONJUNTO DE CAIXAS COMPLETAS, no nivel base: as celulas
// base que pertencem ao rank sao agrupadas em retangulos maximais, e cada
// retangulo vira uma arvore completa -- todo compartimento preenchido, o refino
// interno materializado inteiro.
//
// O QUE PODE DAR ERRADO, e por isso e' MEDIDO em vez de suposto: o t8code
// reparte por elemento, entao ele pode dividir uma FAMILIA entre ranks -- parte
// dos filhos de uma celula base ficando com o vizinho.  Uma caixa completa nao
// representa isso.  O campo `base_dividida` conta os casos; o teste exige zero, e
// se algum dia for diferente de zero a saida e' pedir ao t8code que nao divida
// familias, nao remendar a materializacao.

#include "higtree.h"
#include "coord.h"

#ifdef __cplusplus
extern "C" {
#endif

#define T8_MAX_CAIXAS 4096

typedef struct t8_producao_rank {
    hig_cell *locais[T8_MAX_CAIXAS];   //!< caixas completas deste rank
    int       n_locais;
    hig_cell *franjas[T8_MAX_CAIXAS];  //!< caixas completas da camada de ghost
    int       n_franjas;
    long      n_local;    //!< folhas materializadas nas caixas locais
    long      n_franja;   //!< folhas materializadas nas caixas de franja
    long      n_global;   //!< folhas da floresta inteira, segundo o t8code
    long      base_dividida;  //!< celulas base cuja familia ficou dividida
} t8_producao_rank;

//! \brief Produz a parte deste rank, em caixas completas.  Coletivo.
int t8_produz_por_rank(const Point lo, const Point hi, int nivel_base,
                       const Point alvo, int refinos, t8_producao_rank *out);

//! \brief A mesma decomposicao, com cmesh de BRICK: malha uniforme de dimensoes
//! quaisquer (160x40, por exemplo), que o hipercubo nao representa.
//!
//! Nao preenche franja: quem reparte depois e' o `lbal`, e e' ele quem monta a
//! franja e o grafo de vizinhanca.
int t8_produz_por_rank_brick(const Point lo, const Point hi, const int nb[DIM],
                             t8_producao_rank *out);

//! \brief Brick com REFINO numa caixa.  `caixa_lo`/`caixa_hi` delimitam a
//! regiao a refinar (em coordenadas do DOMINIO), e `refinos` e' quantos niveis.
//!
//! Nao pede balanceamento: a razao entre celulas vizinhas pode chegar a 4:1.
//! Para fronteira imersa o corpo tem de ficar INTEIRAMENTE dentro do nivel
//! fino, com folga maior que o suporte do nucleo (1,5 celulas) -- a guarda
//! `fi_suporte_nivel_trocado` verifica isso.
int t8_produz_por_rank_brick_refinado(const Point lo, const Point hi, const int nb[DIM],
                                      const Point caixa_lo, const Point caixa_hi,
                                      int refinos, t8_producao_rank *out);

void t8_producao_rank_destroi(t8_producao_rank *p);

#ifdef __cplusplus
}
#endif
#endif
