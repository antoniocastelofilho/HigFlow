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

#if DIM == 3
//! Em 3D o corpo e' uma SUPERFICIE, nao uma curva.  Esta funcao extruda a curva
//! fechada em z, colocando marcadores nos centros dos retalhos e o peso igual a'
//! AREA de cada um -- a soma da a area lateral exata.
//!
//! Nao exige topologia, e isso e' proposital: corpo RIGIDO tem peso fixo na
//! criacao.  Topologia (o `DMPlex`) so' faz falta quando for preciso CURVATURA,
//! isto e', no caso de interface entre fluidos.  Ver a secao 8 do projeto.
fi_corpo *fi_cria_extrusao(sim_facet_domain *sfd, const Point *vertices, int nvert,
                           real z0, real z1, real h);

//! Conveniencia: cilindro de eixo z, entre `z0` e `z1`.
fi_corpo *fi_cria_cilindro(sim_facet_domain *sfd, real cx, real cy, real raio,
                           int nlados, real z0, real z1, real h);
#endif

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

//! ATENCAO AO CONTRATO: `fi_espalha` E' DONO DO CAMPO QUE RECEBE.
//!
//! Ele acumula franja->dono com `PetscSFReduce`, o que so' e' correto se o
//! campo contiver APENAS contribuicoes espalhadas.  Passar um campo ja'
//! povoado -- `dpustar`, por exemplo -- faz o reduce somar os valores de franja
//! nos donos, isto e', somar um campo inteiro nas fronteiras de particao.
//!
//! Isso nao e' hipotese: eu mesmo violei o contrato um dia depois de escreve-lo,
//! passando `dpustar` na rota Uhlmann.  A corrida divergia para 1e96 e o
//! diagnostico so' fechou quando anular a forca (alfa = 0) NAO mudou nada --
//! sinal de que a corrupcao nao vinha da forca.
//!
//! Para somar num campo povoado: espalhe num campo ZERADO proprio e some depois.
//!
//! Espalha a forca dos marcadores nas facetas, ACUMULANDO (dp_add_value).
//! Contribuicoes que caem em faceta de franja sao somadas no dono por
//! PetscSFReduce -- o dp_sync nao serve, ele sobrescreve.
void fi_espalha(fi_corpo *c, sim_facet_domain *sfd[DIM],
                distributed_property *dpF[DIM]);

//! O mesmo, multiplicando a contribuicao por `escala`.
//!
//! Existe por causa do Uhlmann: ali a forca nao vai para o campo de FONTE (que a
//! equacao le' e multiplica por dt sozinha), vai CORRIGIR a velocidade
//! provisoria -- `u* += dt * F` -- entao o dt entra aqui.
void fi_espalha_com_escala(fi_corpo *c, sim_facet_domain *sfd[DIM],
                           distributed_property *dpF[DIM], real escala);

//! Grava a malha lagrangeana em VTK: posicoes, forca, velocidade e peso.
//! UM ARQUIVO POR RANK -- `<prefixo>_lag_<rank>-<quadro>.vtk` -- porque os
//! marcadores sao distribuidos por posse euleriana.
//!
//! Sem isto o corpo e' INVISIVEL: o VTK do solver grava so' a malha euleriana,
//! e nao ha' como conferir se os marcadores estao onde se pensa que estao.
void fi_escreve_vtk(const fi_corpo *c, const char *prefixo, int quadro);

//! Quantos pontos de suporte foram descartados por id de faceta inutilizavel
//! (`sfd_get_local_id` devolve -1 para a face marcada por `sfd_adjust_facet_ids`).
//! DEVE SER ZERO.  Diferente de zero significa nucleo somando menos que 1 em
//! algum marcador, e forca menor ali -- corpo levemente poroso, que se le' como
//! defeito de malha.  Acumulado no processo, nao reduzido entre ranks.
long fi_suporte_perdidos(void);

//! A decomposicao de `fi_suporte_perdidos`, que e' o que diz o QUE consertar:
//!   espelho     id negativo pela convencao de `sfd_adjust_facet_ids` -- o
//!               canonico e' o simetrico, e ha' conserto
//!   mapa        a faceta existe mas nao esta' no mapeador deste dominio
//!   sem_faceta  nao ha' faceta naquele ponto (suporte saindo do dominio)
//! Pontos de suporte caidos em celula de tamanho DIFERENTE do h do marcador --
//! o suporte atravessou uma fronteira de refinamento.  DEVE SER ZERO.
//!
//! O nucleo so' e' normalizado para um h; atravessando, a particao da unidade
//! deixa de valer e a forca sai errada naquele marcador, em silencio.  O
//! tratamento correto e' o da versao ADAPTATIVA do metodo (Roma-Peskin-Berger),
//! nao implementado -- entao a restricao e' que o corpo fique inteiramente
//! dentro de um nivel de refinamento, com folga maior que o suporte.
long fi_suporte_nivel_trocado(void);

long fi_perdidos_espelho(void);
long fi_perdidos_mapa(void);
long fi_perdidos_sem_faceta(void);

//! A FORCA DO PASSO, somada sobre as ITERACOES do forcamento.
//!
//! Com forcamento iterado, `fi_forca_total` devolve so' a parcela da ULTIMA
//! iteracao -- que e' justamente a menor.  Ler dali daria arrasto uma ordem de
//! grandeza pequeno demais.  Use: zerar no inicio do passo, acumular apos cada
//! `fi_forca_corpo_rigido`, e ler no fim.
//!
//! O sinal: `f` e' a forca que o corpo aplica AO FLUIDO.  A forca hidrodinamica
//! sobre o corpo e' a reacao, `-forca_passo`.  Para escoamento em +x sobre
//! corpo fixo, o arrasto sai POSITIVO: `f_x` e' negativo (freia o fluido) e o
//! sinal troca.
void fi_zera_forca_passo(fi_corpo *c);
void fi_acumula_forca_passo(fi_corpo *c);
void fi_forca_passo(const fi_corpo *c, real forca[DIM]);

//! Maior |u| entre TODOS os marcadores, de todos os ranks -- o RESIDUO DE NAO
//! ESCORREGAMENTO.  E' o primeiro oraculo do metodo, e o mais barato.
//!
//! Com forca defasada ele NAO vai a zero: vai a O(dt).  Entao o teste nao e' o
//! valor, e' a TAXA -- refinar dt e ver o residuo cair na mesma ordem.  Erro de
//! sinal, de escala ou de peso nao passa nisso; figura de esteira passa.
//!
//! Le' o campo `velocidade`, entao so' vale depois de um `fi_interpola`.
real fi_residuo_max(const fi_corpo *c);

//! Soma global de f_k * dV_k, por direcao -- a forca que DEVE chegar as
//! facetas.  dV_k = peso * h^(DIM-1) e' o VOLUME do marcador; o peso guardado e'
//! comprimento de arco (2D) ou area (3D), que e' geometrico e testavel contra o
//! perimetro, mas nao e' o que a equacao pede.
void fi_forca_total(const fi_corpo *c, real total[DIM]);

#ifdef __cplusplus
}
#endif

#endif
