// *******************************************************************
//  |rot u| do gradiente de velocidade -- a parte PURA do sensor de esteira
// *******************************************************************
//
// POR QUE EXISTE SEPARADO.  A amostragem das facetas depende do solver e so' se
// exercita numa corrida; a conta de |rot u| e' o que erra SINAL ou TROCA EIXO, e
// isso se testa contra forma fechada em milissegundos.  Em `static inline` num
// cabecalho, a mesma definicao serve a producao e ao arreio -- sem copia a
// manter em sincronia.
//
// A forma do argumento NAO e' escolha minha: e' a que o sensor em
// `ft3_escreve_amr_criterio` ja' monta, `du[a][b] = d u_a / d x_b` por diferenca
// centrada nas facetas.  Fazer a funcao receber doze amostras obrigaria o
// chamador a reorganizar o que ele ja' tem, por nada.
//
// Nao depende de mais nada: <math.h> e `real` ja' definido.

#ifndef CRITERIO_ESTEIRA_PURA_H
#define CRITERIO_ESTEIRA_PURA_H

#include <math.h>

//! |rot u| de du[a][b] = d u_a / d x_b.
//!
//!   rot u = ( dw/dy - dv/dz ,  du/dz - dw/dx ,  dv/dx - du/dy )
//!         = ( du[2][1]-du[1][2], du[0][2]-du[2][0], du[1][0]-du[0][1] )
static inline real esteira_omega_grad(const real du[3][3])
{
    const real wx = du[2][1] - du[1][2];
    const real wy = du[0][2] - du[2][0];
    const real wz = du[1][0] - du[0][1];
    return sqrt(wx*wx + wy*wy + wz*wz);
}

#endif
