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

#ifdef __cplusplus
}
#endif

#endif
