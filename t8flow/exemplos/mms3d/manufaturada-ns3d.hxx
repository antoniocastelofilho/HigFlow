// UMA SOLUCAO MANUFATURADA PARA NAVIER-STOKES EM TRES DIMENSOES.
//
// ===========================================================================
// COMO ELA FOI CONSTRUIDA, E POR QUE ASSIM
// ===========================================================================
//
// O CAMPO VEM DE UM POTENCIAL VETORIAL, `u = rot(psi)`, e por isso e'
// solenoidal POR CONSTRUCAO -- nao por conta conferida.  Com
//
//     psi = (0, chi, phi),   phi = chi = sin(pi x) sin(pi y) sin(pi z)
//
// sai
//
//     u_x = d_y phi - d_z chi = pi S_x (C_y S_z - S_y C_z)
//     u_y = -d_x phi          = -pi C_x S_y S_z
//     u_z =  d_x chi          =  pi C_x S_y S_z
//
// com `S_x = sin(pi x)` e `C_x = cos(pi x)`.  As tres componentes sao nao
// nulas, e a divergencia e' zero ate' o ultimo bit porque e' um rotacional.
//
// O LAPLACIANO SAI DE GRACA: cada termo e' um produto de tres funcoes
// trigonometricas de argumento `pi .`, logo `lap = -3 pi^2` vezes o proprio
// termo.  Entao `lap(u) = -3 pi^2 u`, componente a componente.  Isso elimina
// a parte mais arriscada da derivacao a mao.
//
// A PRESSAO E' `cos(pi x) cos(pi y) cos(pi z)`, e a escolha nao e' estetica: a
// derivada NORMAL dela se anula nas SEIS paredes do cubo, que e' o que a
// condicao `d(dp)/dn = 0` do metodo de projecao exige.  Com uma pressao que
// nao tivesse essa propriedade, o erro da projecao na parede contaminaria a
// medida de ordem e seria atribuido ao esquema.
//
// ===========================================================================
// O TERMO FONTE: DERIVADO A MAO, E CONFERIDO
// ===========================================================================
//
//     f_i = d_t u_i + sum_j u_j d_j u_i - (1/Re) lap(u_i) + d_i p
//
// As NOVE derivadas primeiras `d_j u_i` e as tres `d_i p` estao escritas a
// mao, uma a uma.  O convectivo e' MONTADO pelo codigo a partir delas -- isso
// e' algebra sobre formulas exatas, e nao diferenciacao numerica.
//
// A DISTINCAO IMPORTA.  Se `f` fosse calculado por diferencas finitas, o
// oraculo mediria a consistencia entre dois erros meus e nao o esquema.  Mas
// CONFERIR a formula a mao contra diferencas finitas e' outra coisa: e' uma
// guarda contra escorregao de derivacao, e o caso
// `as_derivadas_a_mao_batem_com_diferencas_finitas` faz exatamente isso.

#ifndef T8FLOW_EXEMPLOS_MMS3D_MANUFATURADA_NS3D_HXX
#define T8FLOW_EXEMPLOS_MMS3D_MANUFATURADA_NS3D_HXX

#include "malha3d.hxx"

namespace t8flow3d {

//! A solucao, as derivadas e a fonte.  `T(t) = exp(-t)` multiplica a
//! velocidade; a pressao nao depende do tempo.
struct ManufaturadaNS3 {
  double re = 1.0;

  //! A componente `i` da velocidade.
  double u(const Ponto3 &x, int i, double t) const;

  //! `d u_i / d x_j`.
  double dudx(const Ponto3 &x, int i, int j, double t) const;

  //! `d u_i / d t`.
  double dudt(const Ponto3 &x, int i, double t) const;

  //! `lap(u_i)` -- vale `-3 pi^2 u_i`, e essa e' a razao de o campo ser
  //! produto de tres trigonometricas de mesmo argumento.
  double lap(const Ponto3 &x, int i, double t) const;

  double p(const Ponto3 &x) const;
  double dpdx(const Ponto3 &x, int i) const;

  //! O termo fonte, MONTADO das derivadas acima.
  double f(const Ponto3 &x, int i, double t) const;
};

}  // namespace t8flow3d

#endif
