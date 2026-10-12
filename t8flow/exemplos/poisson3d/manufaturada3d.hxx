// O CATALOGO DE SOLUCOES MANUFATURADAS EM TRES DIMENSOES.
//
// Escolhe-se `u`, DERIVA-SE `f = lap(u)` analiticamente, impoem-se as condicoes
// de contorno a partir do proprio `u`, e compara-se a solucao discreta com `u`.
// A vantagem decisiva e' que a solucao exata e' conhecida em todo ponto, logo o
// erro e' mensuravel sem refino de referencia.
//
// TODAS AS DERIVADAS SAO ESCRITAS A MAO, nunca obtidas por diferenciacao
// numerica -- ela introduziria exatamente o erro que se quer medir.
//
// OS COEFICIENTES SAO DISTINTOS EM CADA MONOMIO, e em tres dimensoes isso
// importa mais que em duas: com `x`, `y` e `z` ha' mais maneiras de trocar dois
// eixos, e um operador que os troque passaria em qualquer campo simetrico.
//
// ESTE ARQUIVO E' EXEMPLO E NAO BIBLIOTECA, e por isso vive em `exemplos/`: um
// catalogo de solucoes com derivadas analiticas so' serve para VERIFICAR, e
// nada em `src/` sabe que equacao se esta' resolvendo.

#ifndef T8FLOW_EX_MANUFATURADA3D_HXX
#define T8FLOW_EX_MANUFATURADA3D_HXX

#include <functional>
#include <string>

#include "malha3d.hxx"

namespace t8flow3d {

//! Uma solucao manufaturada: o campo, as derivadas, a fonte.
struct Manufaturada {
  std::string nome;
  std::string descricao;

  std::function<double(const Ponto3 &)> u;       //!< a solucao
  std::function<double(const Ponto3 &)> dudx;
  std::function<double(const Ponto3 &)> dudy;
  std::function<double(const Ponto3 &)> dudz;
  std::function<double(const Ponto3 &)> f;       //!< `lap(u)`

  //! Grau do polinomio, ou -1 se nao for polinomial.
  int grau = -1;

  //! A derivada NORMAL, para as condicoes de Neumann.
  double dudn(const Ponto3 &p, const double n[3]) const
  {
    return dudx(p) * n[0] + dudy(p) * n[1] + dudz(p) * n[2];
  }
};

//! Busca pelo nome; lanca com a LISTA dos validos quando nao reconhece.
const Manufaturada &manufaturada(const std::string &nome);

//! Os nomes do catalogo, separados por virgula.
std::string nomes_manufaturadas();

}  // namespace t8flow3d

#endif  // T8FLOW_EX_MANUFATURADA3D_HXX
