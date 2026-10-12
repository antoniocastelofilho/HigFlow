// A SOLUCAO MANUFATURADA DA CAIXA COM TAMPA DESLIZANTE.
//
// CONSTRUIDA POR FUNCAO DE CORRENTE, e nao escolhida componente a componente:
// com `u = d(psi)/dy` e `v = -d(psi)/dx` a divergencia e' nula POR CONSTRUCAO,
// exatamente, e nao por cancelamento de termos que alguem conferiu uma vez.
// Numa solucao manufaturada incompressivel essa e' a propriedade que nao pode
// falhar: se `div u != 0`, a projecao tenta corrigir o que e' a PROPRIA
// solucao, e o erro medido deixa de ser do esquema.
//
//     psi(x,y,t) = f(t) G(x) H(y)
//     G(x) = sin^2(pi x)
//     H(y) = sin(pi y) - (1/2) sin(2 pi y)  =  sin(pi y) (1 - cos(pi y))
//
// E E' TRIGONOMETRICA, E NAO POLINOMIAL, de proposito.  Um polinomio de grau
// baixo pode ser reproduzido EXATAMENTE por um esquema de segunda ordem, e
// nesse caso a ordem medida nao e' a do esquema -- e' a exatidao acidental.
// Com senos e cossenos nenhum esquema polinomial e' exato, e o que se mede e'
// o truncamento.
//
// AS CONDICOES SAEM DAS RAIZES DE `G` E `H`, e cada uma tem um proposito:
//
//     G(0) = G(1) = 0          -> u = 0 nas paredes LATERAIS
//     G'(0) = G'(1) = 0        -> v = 0 nas paredes laterais
//     H(0) = H(1) = 0          -> v = 0 no FUNDO e na TAMPA
//     H'(0) = 0                -> u = 0 no FUNDO
//     H'(1) = -2 pi  (nao nulo)-> u != 0 na TAMPA
//
// Logo `u(x,1,t) = -2 pi f(t) sin^2(pi x)` e todo o resto do contorno e' nulo
// em TODO instante: so' a tampa depende do tempo, que e' o que o caso pede.
//
// A PRESSAO E' `cos(pi x) cos(pi y)`, e a escolha tambem nao e' livre: o
// metodo de projecao impoe `d(dp)/dn = 0` nas paredes de Neumann, e uma
// pressao manufaturada com derivada normal NAO nula ali seria incompativel com
// o proprio esquema -- o erro medido misturaria o esquema com a
// inconsistencia que eu mesmo introduzi.  Esta tem derivada normal nula nos
// quatro lados do quadrado unitario (`sin(0) = sin(pi) = 0`) e media zero.

#ifndef T8FLOW_EXEMPLOS_CAIXA_MANUFATURADA_HXX
#define T8FLOW_EXEMPLOS_CAIXA_MANUFATURADA_HXX

#include <cmath>

#include "malha2d.hxx"

namespace caixa {

using t8flow2d::Ponto2;

constexpr double PI = 3.14159265358979323846;

//! O fator temporal.  `f(0) = 0`, logo o REPOUSO e' condicao inicial exata --
//! nao ha' erro de inicializacao somado ao do esquema.
inline double f(double t) { return std::sin(2.0 * PI * t); }
inline double fl(double t) { return 2.0 * PI * std::cos(2.0 * PI * t); }

inline double G(double x) { const double s = std::sin(PI * x); return s * s; }
inline double G1(double x) { return PI * std::sin(2 * PI * x); }
inline double G2(double x) { return 2 * PI * PI * std::cos(2 * PI * x); }
inline double G3(double x) { return -4 * PI * PI * PI * std::sin(2 * PI * x); }

inline double H(double y)
{
  return std::sin(PI * y) - 0.5 * std::sin(2 * PI * y);
}
inline double H1(double y)
{
  return PI * std::cos(PI * y) - PI * std::cos(2 * PI * y);
}
inline double H2(double y)
{
  return -PI * PI * std::sin(PI * y) + 2 * PI * PI * std::sin(2 * PI * y);
}
inline double H3(double y)
{
  return -PI * PI * PI * std::cos(PI * y)
         + 4 * PI * PI * PI * std::cos(2 * PI * y);
}

//! A componente `d` da velocidade exata.
inline double u_exata(const Ponto2 &x, int d, double t)
{
  return (d == 0) ? f(t) * G(x[0]) * H1(x[1])
                  : -f(t) * G1(x[0]) * H(x[1]);
}

inline double p_exata(const Ponto2 &x, double)
{
  return std::cos(PI * x[0]) * std::cos(PI * x[1]);
}

//! O TERMO FONTE: o que sobra ao substituir a solucao na equacao.
//
//     S_d = du_d/dt + (u . grad) u_d + dp/dx_d - (1/Re) lap(u_d)
//
// Escrito na forma ADVECTIVA porque a divergencia e' nula por construcao e as
// duas formas coincidem no continuo.  O esquema usa a forma divergente; usar a
// advectiva aqui nao e' inconsistencia, e' a mesma funcao.
inline double fonte(const Ponto2 &x, int d, double t, double re)
{
  const double X = x[0], Y = x[1];
  const double u = f(t) * G(X) * H1(Y);
  const double v = -f(t) * G1(X) * H(Y);

  if (d == 0) {
    const double ut = fl(t) * G(X) * H1(Y);
    const double ux = f(t) * G1(X) * H1(Y);
    const double uy = f(t) * G(X) * H2(Y);
    const double uxx = f(t) * G2(X) * H1(Y);
    const double uyy = f(t) * G(X) * H3(Y);
    const double px = -PI * std::sin(PI * X) * std::cos(PI * Y);
    return ut + u * ux + v * uy + px - (uxx + uyy) / re;
  }
  const double vt = -fl(t) * G1(X) * H(Y);
  const double vx = -f(t) * G2(X) * H(Y);
  const double vy = -f(t) * G1(X) * H1(Y);
  const double vxx = -f(t) * G3(X) * H(Y);
  const double vyy = -f(t) * G1(X) * H2(Y);
  const double py = -PI * std::cos(PI * X) * std::sin(PI * Y);
  return vt + u * vx + v * vy + py - (vxx + vyy) / re;
}

}  // namespace caixa

#endif  // T8FLOW_EXEMPLOS_CAIXA_MANUFATURADA_HXX
