// Fronteira imersa com forca direta -- o nucleo: os dois operadores de
// transferencia entre a malha euleriana (facetas) e a lagrangeana (marcadores).
//
// PROJETO em doc/projeto-fronteira-imersa.md.  Decisoes: corpo RIGIDO e FIXO,
// forca DEFASADA, so' o newtoniano, estruturas do PETSc.
//
// ESTE ARQUIVO NAO CONHECE O SOLVER de proposito.  Ele toma um sim_facet_domain
// e uma distributed_property por direcao, o que o torna testavel sem montar um
// higflow_solver -- e foi assim que a conservacao da forca foi verificada.
//
// O PAR ADJUNTO E' O QUE IMPORTA.  Espalhar e' o transposto de interpolar:
//
//     u(X_k) = SUM_x  u(x) d_h(x - X_k) h^DIM        interpolar
//     F(x)   = SUM_k  f_k d_h(x - X_k) w_k           espalhar
//
// Com o MESMO nucleo d_h nos dois, a forca total se conserva: o que sai dos
// marcadores chega as facetas.  Trocar um dos dois nucleos quebra a conservacao
// sem quebrar nada que se veja na tela -- e' o defeito que o teste persegue.

#ifndef HIG_FLOW_FRONTEIRA_IMERSA_H
#define HIG_FLOW_FRONTEIRA_IMERSA_H

#include "domain.h"
#include "pdomain.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fi_corpo fi_corpo;

//! Nucleo regularizado de Roma, 3 pontos.  Argumento r = (x - X)/h.
//! Suporte 1,5 celulas para cada lado; particao da unidade em r inteiro+fase.
real fi_delta_roma(real r);

//! Cria o corpo a partir de uma CURVA FECHADA POR SEGMENTOS, dada pelos
//! vertices em ordem.  A curva fecha sozinha: NAO repita o primeiro vertice no
//! fim.
//!
//! Cada segmento e' subdividido para que o espacamento fique proximo de `h`, e
//! os marcadores ficam nos CENTROS dos subsegmentos -- nao nos vertices.  Isso
//! nao e' detalhe: marcador em vertice duplicaria o vertice compartilhado entre
//! dois segmentos consecutivos, e o peso total deixaria de ser o perimetro.
//! Com centros, o peso de cada marcador e' o comprimento do seu subsegmento e a
//! soma da o perimetro exato -- que e' o que `fi_peso_total` afirma.
//!
//! Os marcadores sao distribuidos entre os ranks pela POSSE EULERIANA: cada um
//! fica no rank que possui a celula que o contem, o que e' a condicao para o
//! espalhamento ser local.  `sfd` serve so' para essa localizacao.
//!
//! ABORTA se a soma dos marcadores locais nao bater com o total esperado -- e'
//! marcador reivindicado por dois ranks ou por nenhum, e seguir daria forca
//! errada sem sintoma visivel.
fi_corpo *fi_cria_curva(sim_facet_domain *sfd, const Point *vertices, int nvert,
                        real h);

//! Conveniencia: poligono fechado de `nlados` lados inscrito no circulo, e
//! entao `fi_cria_curva`.  Com nlados grande e' um circulo.
fi_corpo *fi_cria_circulo(sim_facet_domain *sfd, const Point centro, real raio,
                          int nlados, real h);

void fi_destroi(fi_corpo *c);

//! Numero de marcadores DESTE rank.
int fi_num_locais(const fi_corpo *c);

//! Soma global dos pesos -- deve dar o perimetro do circulo.  Usado em teste.
real fi_peso_total(const fi_corpo *c);

//! Interpola a velocidade das facetas nos marcadores.
void fi_interpola(fi_corpo *c, sim_facet_domain *sfd[DIM],
                  distributed_property *dpu[DIM]);

//! Forca direta defasada, CORPO RIGIDO E FIXO: f_k = (0 - u_k) / dt.
//!
//! O nome diz a LEI, nao a operacao, porque ha' duas leis e elas nao sao
//! variantes uma da outra (ver secao 8 do projeto):
//!
//!   contorno rigido curvo   a forca e' MULTIPLICADOR DE LAGRANGE -- vale o que
//!                           for preciso para u = U_corpo.  Sem lei
//!                           constitutiva, magnitude crescendo como 1/dt.
//!   interface entre fluidos a forca e' CONSTITUTIVA -- tensao superficial
//!                           sigma*kappa*n.  Exige CURVATURA, e portanto
//!                           vizinhos: ordem da curva em 2D, triangulacao em 3D.
//!
//! A transferencia (`fi_interpola`, `fi_espalha`) serve aos dois sem diferenca;
//! a lei, nao.  E a distribuicao por posse euleriana que `fi_cria_curva` faz
//! DESTROI a ordem da curva -- correto para o caso rigido, inviavel para o
//! outro.
void fi_forca_corpo_rigido(fi_corpo *c, real dt);

//! Espalha a forca dos marcadores nas facetas, ACUMULANDO (dp_add_value).
//! Contribuicoes que caem em faceta de franja sao somadas no dono por
//! PetscSFReduce -- o dp_sync nao serve, ele sobrescreve.
void fi_espalha(fi_corpo *c, sim_facet_domain *sfd[DIM],
                distributed_property *dpF[DIM]);

//! Soma global de f_k * w_k, por direcao.  O que DEVE chegar as facetas.
void fi_forca_total(const fi_corpo *c, real total[DIM]);

#ifdef __cplusplus
}
#endif

#endif
