// Front-tracking 3D -- a frente e' uma SUPERFICIE TRIANGULADA fechada.
//
// POR QUE E' OUTRO ARQUIVO, E NAO UM #if NO NUCLEO 2D.  A frente 2D e' uma
// polilinha cuja conectividade e' IMPLICITA no indice: o marcador i liga em
// i+1.  Aqui a conectividade e' EXPLICITA -- lista de triangulos mais
// adjacencia por aresta -- e nenhuma das operacoes 2D sobrevive: area de laco
// vira soma de triangulos, curvatura por circuncirculo nao generaliza, e a
// cirurgia deixa de ser inserir/remover ponto para ser dividir, colapsar e
// GIRAR aresta preservando a variedade.  Os dois nucleos quase nao tem codigo
// em comum, e o 2D esta' verificado fechando um benchmark: nao se mexe nele.
//
// O QUE SE HERDA PRONTO.  A maquinaria de transferencia
// (hig-flow-fronteira-imersa.h) ja' e' 3D: trata superficie como CODIMENSAO 1
// com peso = area, e a lei do espalhamento e' dV = peso * h^(DIM-1), generica.
// Busca de suporte na malha escalonada, nucleo de Roma e interpolacao vem de
// la' sem alteracao.  O que faltava era a representacao da frente.
//
// A TENSAO SUPERFICIAL NAO PASSA POR CURVATURA.  A forca sai da integral de
// linha na borda de cada triangulo,
//
//     F = sigma * contorno (t x n) ds ,
//
// que para aresta reta de `a` a `b` num triangulo de normal `n` vale
// (b - a) x n.  Curvatura nunca e' calculada.  Duas consequencias que valem
// mais que a economia:
//
//   1. Num triangulo PLANO a soma das tres arestas e' (sum (b-a)) x n = 0.
//      A forca liquida por triangulo e' EXATAMENTE zero em ponto flutuante, e
//      portanto a forca total numa superficie fechada tambem e'.  Isso e'
//      estrutura, nao tolerancia -- e por isso NAO e' um oraculo forte.  Ver o
//      paragrafo seguinte.
//
//   2. A forca surge porque a aresta INTERIOR e' compartilhada por dois
//      triangulos de normais DIFERENTES, e as duas contribuicoes nao se
//      cancelam.  E' a curvatura entrando pela diferenca de normais, sem
//      nunca ser nomeada.
//
// O ORACULO DE VERDADE E' OUTRO: numa esfera de raio R a forca por unidade de
// area tem de tender a 2*sigma/R -- DOIS sobre R, nao um.  Em 2D o salto de
// Laplace e' sigma/R e em 3D e' 2*sigma/R; copiar o oraculo do B2 bidimensional
// para ca' valida o numero errado.  Esta advertencia esta' aqui porque o caso
// 2D fechou a 0,02% e a tentacao de reaproveitar o numero e' real.
//
// ALCANCE DESTA FASE: geometria e oraculos, SEM solver.  Area, volume, forca e
// cirurgia verificados contra esfera analitica antes de qualquer acoplamento --
// o mesmo degrau que, em 2D, pegou os defeitos antes de custarem horas de
// maquina.  Adveccao nao esta' aqui de proposito: entra com o B1 3D, junto com
// o teste que a exercita.

#ifndef HIG_FLOW_FRONT_TRACKING_3D_H
#define HIG_FLOW_FRONT_TRACKING_3D_H

#include "coord.h"
#include "types.h"

#if DIM != 3
#error "hig-flow-front-tracking-3d.c e' 3D; a frente 2D esta' em hig-flow-front-tracking.c"
#endif

typedef struct ft3_superficie ft3_superficie;

// --- criacao / destruicao -------------------------------------------------

//! Esfera por subdivisao de icosaedro (`icosfera`).  `nsub` subdivisoes dao
//! 20*4^nsub triangulos e 10*4^nsub+2 vertices, com triangulos quase uniformes
//! -- o que importa para a cirurgia e para a qualidade da forca.
//!
//! Os triangulos saem orientados para FORA (normal apontando para longe do
//! centro), e `ft3_volume` devolve valor positivo.  A orientacao e' contrato:
//! a forca e a normal dependem dela.
ft3_superficie *ft3_cria_esfera(const Point centro, real raio, int nsub);

//! Superficie a partir de vertices e triangulos dados.  Copia os dois vetores.
//! `ds_alvo` e' o espacamento alvo de ARESTA (criterio da cirurgia); se vier
//! <= 0, sai da media dos comprimentos de aresta.
ft3_superficie *ft3_cria_malha(const Point *vert, int nv,
                               const int (*tri)[3], int nt, real ds_alvo);

void ft3_destroi(ft3_superficie *s);

// --- consulta -------------------------------------------------------------

int  ft3_num_vertices(const ft3_superficie *s);
int  ft3_num_triangulos(const ft3_superficie *s);
real ft3_ds_alvo(const ft3_superficie *s);

//! Vetor de posicoes dos vertices, com `ft3_num_vertices` entradas.  Aponta
//! para dentro da estrutura: nao liberar, e invalido apos cirurgia.
const Point *ft3_posicoes(const ft3_superficie *s);

//! Triangulos, com `ft3_num_triangulos` entradas de 3 indices.  Mesmas regras.
const int (*ft3_triangulos(const ft3_superficie *s))[3];

//! AREA da superficie: soma das areas dos triangulos.  Oraculo da esfera:
//! 4*pi*R^2.
real ft3_area(const ft3_superficie *s);

//! VOLUME fechado, pelo teorema do divergente: V = (1/6) sum det[x1,x2,x3].
//! Exige orientacao para fora e superficie FECHADA; devolve negativo se a
//! orientacao estiver invertida, o que e' diagnostico e nao acidente.
//! Oraculo da esfera: (4/3)*pi*R^3.  E' o analogo 3D da area do B1: no teste
//! reversivel o volume tem de voltar.
real ft3_volume(const ft3_superficie *s);

//! Caracteristica de Euler V - A + F.  Para superficie fechada de genero 0
//! (esfera topologica) vale 2.  E' o oraculo de TOPOLOGIA da cirurgia: qualquer
//! divisao, colapso ou giro de aresta tem de preserva-la.
int ft3_euler(const ft3_superficie *s);

//! Comprimento minimo e maximo de aresta, e a razao de aspecto do PIOR
//! triangulo (maior lado sobre o menor).  Diagnostico da qualidade da malha,
//! que e' o que a cirurgia existe para manter.
void ft3_qualidade(const ft3_superficie *s, real *amin, real *amax, real *pior_aspecto);

// --- tensao superficial ---------------------------------------------------

//! Forca de tensao superficial nos VERTICES, pela integral de linha na borda de
//! cada triangulo.  Escreve `pos[i]`, `forca[i]` e `peso[i]` para i em
//! [0, ft3_num_vertices), todos alocados pelo chamador.
//!
//! `peso[i]` e' a AREA associada ao vertice (um terco da area de cada triangulo
//! incidente), que e' o peso que a fronteira imersa espera para superficie em
//! 3D -- codimensao 1, dV = peso * h^2.
//!
//! A forca aponta para DENTRO numa superficie convexa orientada para fora: a
//! tensao contrai.  A soma de `forca` sobre todos os vertices e' zero por
//! construcao (ver o cabecalho); o que se verifica de verdade e'
//! |forca[i]| / peso[i] -> 2*sigma/R na esfera.
void ft3_forcas_tensao(const ft3_superficie *s, real sigma,
                       Point *pos, Point *forca, real *peso);

// --- adveccao -------------------------------------------------------------

//! Campo de velocidade avaliado numa posicao e num instante.  E' a UNICA porta
//! da fisica para dentro da adveccao: no B1 3D e' o campo de deformacao de
//! LeVeque, analitico; no acoplamento sera' a interpolacao da malha.
typedef void (*ft3_campo_u)(const Point x, real t, void *ctx, real u[DIM]);

//! Avanca a superficie um passo `dt` a partir de `t`, por RK2 do ponto medio --
//! segunda ordem no tempo.  A conectividade nao muda: so' as posicoes andam.
//! Euler explicito nao fecha o volume de volta no teste reversivel, do mesmo
//! jeito que nao fechava a area em 2D.
void ft3_advecta(ft3_superficie *s, ft3_campo_u u, void *ctx, real t, real dt);

// --- cirurgia -------------------------------------------------------------

//! Mantem o espacamento de aresta proximo de `ds_alvo`: divide aresta acima de
//! 2*ds_alvo, colapsa abaixo de ds_alvo/2, e gira aresta para melhorar o pior
//! angulo.  Devolve quantas operacoes fez.
//!
//! PARA EM VEZ DE CORROMPER.  Colapso que violaria a condicao de elo (criando
//! superficie nao-variedade) ou que inverteria triangulo incidente e' RECUSADO,
//! nao forcado.  Uma aresta que nao pode ser tratada fica como esta': malha
//! pior e' recuperavel, malha nao-variedade nao e'.
int ft3_cirurgia(ft3_superficie *s);

// --- persistencia e saida -------------------------------------------------

//! Grava a superficie em texto, com `%.17g` -- ida-e-volta bit a bit.  Ao
//! contrario do 2D, a CONECTIVIDADE tambem vai: em 3D ela nao se reconstroi da
//! ordem dos vertices.  Devolve 0 em sucesso.
int ft3_grava(const ft3_superficie *s, const char *caminho);

//! Le o que `ft3_grava` escreveu.  Devolve NULL em arquivo ausente, cabecalho
//! estranho, arquivo truncado ou indice de triangulo fora de faixa -- nunca uma
//! superficie meio lida.
ft3_superficie *ft3_le(const char *caminho);

//! VTK polydata com os triangulos, para visualizacao.
int ft3_escreve_vtk(const ft3_superficie *s, const char *caminho);

#endif
