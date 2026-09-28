// Front-tracking -- a frente lagrangeana ORDENADA que se move com o fluido.
//
// ESTE ARQUIVO E' O ALTERNATIVO AO VOF, NAO UM COMPLEMENTO.  Ver
// doc/projeto-bifasico-front-tracking.md.  VOF captura a interface por um campo
// de fracao advectado (euleriano); front-tracking a RASTREIA por marcadores
// conectados.  Os dois ficam separados no codigo, e este e' o comeco do lado
// front-tracking.
//
// A DIFERENCA PARA O CORPO RIGIDO (hig-flow-fronteira-imersa.h): la' os
// marcadores sao DESCONEXOS e distribuidos por posse euleriana, o que DESTROI a
// ordem da curva -- correto para forca de multiplicador de Lagrange, que nao
// precisa de vizinho.  Aqui a frente e' CONECTADA: a curvatura (B2+) e a medida
// de area (B1) exigem a ordem.  Por isso a frente e' REPLICADA e ORDENADA, nao
// um DMSwarm por posse.  (Decisao #1 do projeto: "comecar replicado".)
//
// FASE B1: adveccao cinematica em campo PRESCRITO.  Sem solver, sem tensao
// superficial.  Isola adveccao + cirurgia de toda a fisica.  O oraculo e' o
// teste reversivel de Rider-Kothe: a frente volta a' forma inicial em t=T, e o
// erro de area e' a medida.  Por isso o nucleo NAO conhece o solver -- recebe a
// velocidade por CALLBACK, analitica no B1 e interpolada da malha depois.

#ifndef HIG_FLOW_FRONT_TRACKING_H
#define HIG_FLOW_FRONT_TRACKING_H

#include "coord.h"   // Point, real, DIM, macros POINT_*

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ft_frente ft_frente;

//! Campo de velocidade avaliado numa posicao e num instante.  E' a UNICA porta
//! da fisica para dentro da adveccao: no B1 e' um vortice analitico; a partir do
//! B2 sera' a interpolacao da velocidade da malha nos marcadores (o mesmo papel
//! que `fi_interpola` faz para o corpo rigido).  Escreve `u[DIM]`.
typedef void (*ft_campo_u)(const Point x, real t, void *ctx, real u[DIM]);

// --- criacao / destruicao -------------------------------------------------

//! Cria a frente a partir de uma CURVA FECHADA POR SEGMENTOS, dada pelos
//! vertices em ordem.  A curva fecha sozinha: NAO repita o primeiro vertice.
//!
//! Ao contrario do corpo rigido, os marcadores ficam NOS VERTICES da polilinha
//! reamostrada, e a conectividade e' implicita na ordem do vetor: o marcador `i`
//! liga em `i+1`, e o ultimo liga no primeiro.  A reamostragem poe marcadores a
//! um espacamento proximo de `ds_alvo` ao longo da curva.
ft_frente *ft_cria_curva(const Point *vertices, int nvert, real ds_alvo);

//! Conveniencia: circulo de `nmarc` marcadores igualmente espacados.  Com
//! `nmarc` grande aproxima o circulo; `ds_alvo` sai de 2*pi*raio/nmarc.
ft_frente *ft_cria_circulo(const Point centro, real raio, int nmarc);

void ft_destroi(ft_frente *f);

// --- consulta -------------------------------------------------------------

//! Numero de marcadores da frente.
int ft_num(const ft_frente *f);

//! Espacamento alvo entre marcadores (o criterio da cirurgia).
real ft_ds_alvo(const ft_frente *f);

//! AREA fechada pela frente, pela formula do laco (shoelace).  E' o ORACULO do
//! B1: no teste reversivel deve voltar ao valor inicial em t=T.  Devolve o
//! modulo, entao independe da orientacao (horaria ou anti-horaria) da curva.
real ft_area(const ft_frente *f);

//! Perimetro da frente -- a soma dos comprimentos dos segmentos, fechando o
//! ultimo->primeiro.  Cresce muito na deformacao reversivel; e' o que a cirurgia
//! reamostra.
real ft_perimetro(const ft_frente *f);

//! Copia as posicoes dos marcadores para `dst` (que deve caber `ft_num`
//! pontos).  Usado pelo arreio para medir o desvio de forma sem abrir a struct.
void ft_posicoes(const ft_frente *f, Point *dst);

// --- geometria: curvatura e normal ----------------------------------------

//! VETOR CURVATURA de cada marcador -- kappa*n, apontando para o centro do
//! circulo osculador (o lado concavo), com modulo 1/R.  E' a quantidade que a
//! forca de tensao superficial pede: F = sigma * (vetor curvatura), espalhada
//! pelo nucleo.  Independe da orientacao (horaria/anti-horaria) da curva: o
//! centro osculador esta' sempre do lado concavo, entao numa gota convexa
//! aponta para dentro e a pressao sobe dentro -- Laplace.
//!
//! Calculado pelo circulo que passa por tres marcadores consecutivos
//! (i-1, i, i+1).  Onde os tres sao quase colineares (curva reta), o vetor e'
//! ~zero.  `kv` deve caber `ft_num` pontos.
void ft_vetor_curvatura(const ft_frente *f, Point *kv);

//! Curvatura ESCALAR (modulo 1/R) de cada marcador -- o oraculo: para um circulo
//! de raio R, todo marcador da' 1/R.  `kappa` deve caber `ft_num` reais.
void ft_curvatura(const ft_frente *f, real *kappa);

// --- funcao indicadora (para a forca balanceada) ---------------------------

//! AREA da regiao fechada pela frente que cai DENTRO da caixa alinhada aos
//! eixos [lo,hi].  Recorte de Sutherland-Hodgman do poligono pela caixa, e
//! entao a formula do laco.  Caixa inteiramente dentro da frente devolve a area
//! da caixa; inteiramente fora devolve 0; cortada devolve a fracao exata.
//!
//! E' a fracao volumetrica calculada DA GEOMETRIA da frente -- o analogo exato
//! do campo de cor do VOF, so' que derivado da interface em vez de advectado.
//! Serve a' formulacao BALANCEADA da tensao superficial: com ela, a forca vira
//! sigma*kappa*grad(H), e o gradiente pode ser o MESMO operador discreto que a
//! projecao inverte -- que e' a condicao para as correntes parasitas
//! cancelarem.  Ver Francois et al. (2006).
//!
//! ORACULO: somar sobre uma particao do plano que contenha a frente devolve
//! `ft_area`.
real ft_area_na_caixa(const ft_frente *f, const real lo[DIM], const real hi[DIM]);

// --- forca de tensao superficial ------------------------------------------

//! Nucleo regularizado de Roma, 3 pontos (o MESMO do corpo rigido -- repetido
//! aqui para o modulo ficar standalone).  Argumento r = (x - X)/h.  Suporte 1,5
//! celulas para cada lado; particao da unidade em r inteiro+fase.
real ft_delta_roma(real r);

//! Espalha a forca de tensao superficial F = sigma * kappa * n dos marcadores
//! numa MALHA UNIFORME CENTRADA NA CELULA (a versao standalone; o acoplamento
//! com a malha escalonada do solver vem na fase seguinte).  ACUMULA em fx/fy.
//!
//! Cada marcador leva f_k = sigma * (vetor curvatura)_k com peso ds_k (o
//! comprimento de arco que ele representa: metade de cada segmento vizinho).
//! O espalhamento usa o produto tensorial do nucleo de Roma nas duas direcoes.
//!
//! Grade: nx*ny celulas, origem (ox,oy), espacamento h; a celula (i,j) tem
//! centro em (ox+(i+0.5)h, oy+(j+0.5)h).  fx/fy tem nx*ny reais (indice j*nx+i).
//!
//! DOIS ORACULOS que este arranjo permite:
//!   conservacao  SUM_grade F * h^2 == sigma * SUM_k (vetor curvatura)_k * ds_k,
//!                porque o nucleo soma 1 (particao da unidade).  Precisao de
//!                maquina se o suporte nao vazar da grade.
//!   equilibrio   sigma * SUM_k (vetor curvatura)_k * ds_k -> 0 para curva
//!                fechada (a integral do vetor curvatura ao redor e' nula) --
//!                a forca liquida da tensao superficial numa gota em repouso.
void ft_espalha_tensao(const ft_frente *f, real sigma,
                       real ox, real oy, real h, int nx, int ny,
                       real *fx, real *fy);

//! A soma dos vetores curvatura ponderados pelo arco, por direcao:
//! SUM_k (vetor curvatura)_k * ds_k.  E' a forca total (sem sigma) que DEVE
//! chegar a' grade, e que tende a zero para curva fechada.  Usado no oraculo.
void ft_integral_curvatura(const ft_frente *f, real integral[DIM]);

//! Preenche, por marcador: posicao, forca de tensao superficial
//! F = sigma * (vetor curvatura), e peso ds (comprimento de arco -- metade de
//! cada segmento vizinho).  Os arrays devem caber `ft_num` elementos.  E' o que
//! o adaptador do solver espalha na malha escalonada: mantem a geometria (que e'
//! testada aqui) fora do adaptador.  A soma de F*ds tende a zero (Laplace).
void ft_forcas_tensao(const ft_frente *f, real sigma,
                      Point *pos, Point *forca, real *ds);

//! Semieixos EQUIVALENTES da frente, pelos momentos de area de segunda ordem:
//! a = 2*sqrt(Mxx/A), b = 2*sqrt(Myy/A), com Mxx = integral de (x-xc)^2 dA.
//! Para uma elipse verdadeira devolve os semieixos exatos.
//!
//! E' a medida de deformacao usada na comparacao com o VOF: la' os mesmos
//! momentos saem da integral da fracao volumetrica, entao os dois lados medem a
//! MESMA grandeza, e nao duas parecidas.  A deformacao com sinal e'
//! D = (a-b)/(a+b), que troca de sinal a cada meio periodo de oscilacao --
//! diferente da circularidade, que e' sempre <= 1 e oscila no DOBRO da
//! frequencia.
void ft_semieixos(const ft_frente *f, real *a, real *b);

// --- adveccao + cirurgia --------------------------------------------------

//! Avanca a frente um passo `dt` a partir do instante `t`, no campo `u`:
//!
//!     X^{n+1} = X^n + dt * u(X, t)
//!
//! integrado por RK2 do ponto medio (segunda ordem no tempo) -- o teste
//! reversivel de Rider-Kothe distingue a ordem, e Euler explicito nao fecha a
//! area de volta.  `ctx` e' repassado ao callback sem interpretacao.
void ft_advecta(ft_frente *f, ft_campo_u u, void *ctx, real t, real dt);

//! CIRURGIA: reamostra a frente para manter o espacamento proximo de
//! `ds_alvo`.  Insere um marcador (ponto medio) onde um segmento passa de
//! `2*ds_alvo`; remove um marcador onde o segmento cai abaixo de `ds_alvo/2`.
//! Preserva a ordem e o fechamento.  E' maquinaria LOCAL e barata -- o mesmo
//! criterio Delta s ~ h do corpo rigido -- e e' o que impede a frente esticada
//! de perder resolucao (e portanto area) no filamento fino do vortice.
//!
//! NAO trata mudanca de topologia (coalescencia/quebra): o B1 assume topologia
//! fixa.  Devolve o numero de marcadores DEPOIS da cirurgia.
int ft_cirurgia(ft_frente *f);

// --- saida ----------------------------------------------------------------

//! Grava a frente em VTK como POLILINHA (celulas de linha, nao pontos), para
//! que ela apareca COMO CURVA no visualizador -- fechando o ultimo no primeiro.
//! Nome: `<prefixo>_ft_<quadro>.vtk`.  Frente replicada: um arquivo so'.
void ft_escreve_vtk(const ft_frente *f, const char *prefixo, int quadro);

// --- persistencia ---------------------------------------------------------

//! Grava a frente em TEXTO (numero de marcadores, ds_alvo, e as posicoes), com
//! precisao de ida-e-volta exata.  Devolve o numero de marcadores gravados, ou
//! -1 em erro.
//!
//! POR QUE E' NECESSARIA.  O `h.save` do solver guarda campos na malha, e a
//! frente NAO e' campo: nenhum campo a determina (a fracao volumetrica que ela
//! gera nao devolve a ordem nem o espacamento dos marcadores).  Retomar uma
//! corrida sem este arquivo poe os campos no instante salvo e a frente em t=0 --
//! e roda, em silencio, sobre um estado incoerente.
int ft_grava(const ft_frente *f, const char *arquivo);

//! Le uma frente gravada por `ft_grava`.  Devolve NULL se o arquivo nao existe
//! (ausencia nao e' erro -- quem chama decide se cria a frente inicial) ou se o
//! conteudo nao confere.
//!
//! NAO reamostra: os marcadores voltam nas posicoes exatas em que estavam.
//! Passar pelo `ft_cria_curva` os deslocaria, e a corrida retomada deixaria de
//! ser a continuacao da mesma.
ft_frente *ft_le(const char *arquivo);

#ifdef __cplusplus
}
#endif

#endif
