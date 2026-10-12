#include "mls2d.hxx"

#include <stdexcept>
#include <string>

#include "coord.h"
#include "wls.h"

namespace t8flow2d {

static_assert(sizeof(Ponto2) == sizeof(double) * 2,
              "Ponto2 tem de ter o layout de double[2] para atravessar ao wls");

namespace {

//! A travessia para o herdado: conversao de TIPO, nao de dados.
//
// `std::array<double,2>` tem layout garantido de `double[2]`, logo um vetor dela
// e' contiguamente um `Point[]` valido quando DIM==2.  O `static_assert` abaixo
// e' o que impede isto de virar suposicao se alguem construir com DIM=3.
Point *como_pontos(std::vector<Ponto2> &v)
{
  static_assert(DIM == 2, "mls2d e' o codigo 2D: construa esta unidade com DIM=2");
  return reinterpret_cast<Point *>(v.data());
}

}  // namespace

Mls::Mls(int ordem, int max_pontos)
    : ordem_(ordem), max_pontos_(max_pontos)
{
  if (ordem < 1)
    throw std::invalid_argument("Mls: ordem tem de ser >= 1");
  if (max_pontos < min_pontos(ordem))
    throw std::invalid_argument(
        "Mls: max_pontos (" + std::to_string(max_pontos) + ") abaixo do minimo " +
        std::to_string(min_pontos(ordem)) + " para ordem " + std::to_string(ordem));

  wls_ = wls_create(2, ordem, max_pontos);
  if (wls_ == nullptr)
    throw std::runtime_error("Mls: wls_create devolveu nulo");
}

Mls::~Mls()
{
  if (wls_ != nullptr) wls_destroy(wls_);
}

int Mls::min_pontos(int ordem)
{
  return wls_num_min_points(2, ordem);
}

std::vector<double> Mls::pesos(const Ponto2 &alvo, const std::vector<Ponto2> &suporte) const
{
  const int n = static_cast<int>(suporte.size());
  if (n < min_pontos(ordem_))
    throw std::invalid_argument(
        "Mls::pesos: suporte de " + std::to_string(n) + " pontos, minimo " +
        std::to_string(min_pontos(ordem_)) + " para ordem " + std::to_string(ordem_));
  if (n > max_pontos_)
    throw std::invalid_argument(
        "Mls::pesos: suporte de " + std::to_string(n) + " pontos excede max_pontos " +
        std::to_string(max_pontos_));

  // O herdado recebe os pontos por ponteiro NAO-const e pode reordena-los
  // internamente, entao a copia nao e' desperdicio: e' o que preserva o suporte
  // do chamador.
  std::vector<Ponto2> copia = suporte;
  std::vector<double> w(static_cast<size_t>(n), 0.0);
  Point x = {alvo[0], alvo[1]};

  wls_set_points_and_calc(wls_, n, como_pontos(copia), x, w.data());
  return w;
}

double Mls::interpola(const Ponto2 &alvo, const std::vector<Ponto2> &suporte,
                      const std::vector<double> &valores) const
{
  if (valores.size() != suporte.size())
    throw std::invalid_argument(
        "Mls::interpola: " + std::to_string(valores.size()) + " valores para " +
        std::to_string(suporte.size()) + " pontos de suporte");

  const std::vector<double> w = pesos(alvo, suporte);
  double v = 0.0;
  for (size_t i = 0; i < w.size(); i++) v += w[i] * valores[i];
  return v;
}

}  // namespace t8flow2d
