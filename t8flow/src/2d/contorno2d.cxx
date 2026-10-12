// Implementacao das condicoes de fronteira.  A DERIVACAO de (A, B) para cada tipo
// esta' no `contorno2d.hxx`, com as contas -- aqui ficam as decisoes de
// implementacao, que sao duas:
//
//   1. UMA conta para os tres tipos, parametrizada por (alfa, beta).  Tres contas
//      separadas divergiriam, e divergir aqui significa Robin certo e Neumann
//      errado -- ou pior, os dois plausiveis;
//   2. o caso DEGENERADO e' recusado, nao aproximado.  `alfa*h/2 + beta == 0`
//      anula o denominador: a condicao nao determina o fantasma naquela celula, e
//      devolver qualquer numero ali seria inventar dado.
#include "contorno2d.hxx"

#include <cmath>
#include <string>

namespace t8flow2d {

void Contorno::afim(const Ponto2 &x_parede, double h, double *A, double *B) const
{
  if (!(h > 0.0))
    throw std::invalid_argument("Contorno::afim: h = " + std::to_string(h) +
                                " -- a distancia ao fantasma tem de ser positiva");
  if (!g)
    throw std::invalid_argument("Contorno::afim: condicao sem dado `g` instalado");

  const double gv = g(x_parede);

  // Dirichlet tem caminho PROPRIO no construtor de estencil -- o ponto de suporte
  // e' a parede, com valor conhecido.  Aqui ele aparece como o limite beta=0, que
  // e' a reflexao antissimetrica (A=-1, B=2g); serve de verificacao da formula, e
  // e' o que o teste compara.
  const double den = alfa * h / 2.0 + beta;
  if (std::fabs(den) <= 1e-300)
    throw std::invalid_argument(
        "Contorno::afim: condicao degenerada -- alfa*h/2 + beta = 0 com alfa=" +
        std::to_string(alfa) + ", beta=" + std::to_string(beta) + ", h=" +
        std::to_string(h) + "; ela nao determina o valor do fantasma");

  *A = (beta - alfa * h / 2.0) / den;
  *B = (gv * h) / den;
}

}  // namespace t8flow2d
