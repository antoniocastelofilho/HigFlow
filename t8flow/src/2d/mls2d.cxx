// Implementacao do discretizador 2D.  O PORQUE de cada escolha de interface esta'
// no `mls2d.hxx`; aqui ficam as decisoes de implementacao, que sao tres:
//
//   1. a travessia para o nucleo herdado e' `reinterpret_cast<Point *>` sobre um
//      vetor de `std::array<double,2>`, sustentada por `static_assert` -- conversao
//      de TIPO, nao de dados;
//   2. o suporte do chamador e' COPIADO antes de ir ao herdado, porque a API dele
//      recebe ponteiro nao-const e pode reordenar os pontos;
//   3. as recusas lancam `std::invalid_argument` com os NUMEROS na mensagem --
//      "suporte de 5 pontos, minimo 6 para ordem 2" poupa a rodada de depuracao
//      que um "argumento invalido" obriga.
#include "mls2d.hxx"

#include <algorithm>
#include <cmath>
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

namespace {

//! Monomios da base de ordem `ord` em 2D, em coordenadas RELATIVAS ao alvo, e as
//! derivadas deles no alvo (que e' a origem relativa).
//
// Base relativa nao e' gosto: em coordenadas absolutas a matriz normal fica mal
// condicionada quando o dominio esta' longe da origem.  E com ela, o vetor de
// avaliacao no alvo e' trivial -- `e_0` para valor, `e_1` ou `e_2` para derivada --,
// que e' o que torna os pesos de derivada uma troca de lado direito.
int base_rel(double dx, double dy, int ord, double *b)
{
  int n = 0;
  b[n++] = 1.0;
  if (ord >= 1) { b[n++] = dx; b[n++] = dy; }
  if (ord >= 2) { b[n++] = dx * dx; b[n++] = dx * dy; b[n++] = dy * dy; }
  if (ord >= 3) {
    b[n++] = dx * dx * dx; b[n++] = dx * dx * dy;
    b[n++] = dx * dy * dy; b[n++] = dy * dy * dy;
  }
  return n;
}

constexpr int NMAX = 10;   // ordem 3 em 2D

//! Resolve `M y = rhs` por Gauss com pivotamento parcial.
//
// Com pivotamento porque `M[0][0]` pode ser pequeno quando o suporte e' simetrico
// em torno do alvo, e sem ele o resultado degrada em silencio.
bool resolve_pequeno(double M[NMAX][NMAX], double *rhs, int n)
{
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int i = c + 1; i < n; i++)
      if (std::fabs(M[i][c]) > std::fabs(M[piv][c])) piv = i;
    if (std::fabs(M[piv][c]) < 1e-300) return false;
    if (piv != c) {
      for (int j = 0; j < n; j++) std::swap(M[c][j], M[piv][j]);
      std::swap(rhs[c], rhs[piv]);
    }
    for (int i = c + 1; i < n; i++) {
      const double m = M[i][c] / M[c][c];
      for (int j = c; j < n; j++) M[i][j] -= m * M[c][j];
      rhs[i] -= m * rhs[c];
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    double sv = rhs[i];
    for (int j = i + 1; j < n; j++) sv -= M[i][j] * rhs[j];
    rhs[i] = sv / M[i][i];
  }
  return true;
}

//! O nucleo comum: pesos do funcional `alvo_base` (valor ou derivada).
std::vector<double> pesos_de(const Ponto2 &alvo, const std::vector<Ponto2> &sup,
                             int ord, int qual)
{
  const int n = static_cast<int>(sup.size());
  double M[NMAX][NMAX] = {};
  double b[NMAX];
  int nmon = 0;

  // O PESO E' O MESMO DO HERDADO -- WLSMOVING, `D_i = 1/sqrt(dist^2 + EPSMACH)`,
  // aplicado ao quadrado na matriz normal.  Usar outro peso daria um ajuste
  // diferente que reproduz os mesmos polinomios com pesos diferentes, e a
  // comparacao com o herdado perderia o sentido -- foi exatamente o defeito que o
  // lote de estenceis em GPU teve.
  std::vector<double> p2(static_cast<size_t>(n));
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    nmon = base_rel(dx, dy, ord, b);
    p2[static_cast<size_t>(k)] = 1.0 / (dx * dx + dy * dy + EPSMACH);
    for (int i = 0; i < nmon; i++)
      for (int j = 0; j < nmon; j++) M[i][j] += p2[static_cast<size_t>(k)] * b[i] * b[j];
  }

  // O LADO DIREITO E' O VETOR DE AVALIACAO: `e_0` para valor, `e_1` para `d/dx`,
  // `e_2` para `d/dy` -- porque a base e' relativa ao alvo.
  double y[NMAX] = {};
  if (qual < 0 || qual >= nmon) return std::vector<double>(static_cast<size_t>(n), 0.0);
  y[qual] = 1.0;
  if (!resolve_pequeno(M, y, nmon))
    return std::vector<double>(static_cast<size_t>(n), 0.0);

  std::vector<double> w(static_cast<size_t>(n), 0.0);
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    base_rel(dx, dy, ord, b);
    double sv = 0.0;
    for (int i = 0; i < nmon; i++) sv += y[i] * b[i];
    w[static_cast<size_t>(k)] = p2[static_cast<size_t>(k)] * sv;
  }
  return w;
}

}  // namespace

std::vector<double> Mls::pesos_proprios(const Ponto2 &alvo,
                                        const std::vector<Ponto2> &suporte) const
{
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls::pesos_proprios: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return pesos_de(alvo, suporte, ordem_, 0);
}

std::vector<double> Mls::pesos_derivada(const Ponto2 &alvo,
                                        const std::vector<Ponto2> &suporte, int dir) const
{
  if (dir < 0 || dir > 1)
    throw std::invalid_argument("Mls::pesos_derivada: direcao " + std::to_string(dir) +
                                " -- em 2D sao 0 (x) e 1 (y)");
  if (ordem_ < 1)
    throw std::invalid_argument("Mls::pesos_derivada: ordem 0 nao tem derivada");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls::pesos_derivada: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)) + " para ordem " +
                                std::to_string(ordem_));
  // `e_1` e' `d/dx`, `e_2` e' `d/dy` -- a ordem dos monomios em `base_rel`.
  return pesos_de(alvo, suporte, ordem_, 1 + dir);
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
