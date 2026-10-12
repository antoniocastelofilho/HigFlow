#include "mls3d.hxx"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "coord.h"
#include "wls.h"

namespace t8flow3d {

static_assert(sizeof(Ponto3) == sizeof(double) * 3,
              "Ponto3 tem de ter o layout de double[3] para atravessar ao wls");

namespace {

//! A travessia para o herdado: conversao de TIPO, nao de dados.
//
// `std::array<double,3>` tem layout garantido de `double[3]`, logo um vetor dela
// e' contiguamente um `Point[]` valido quando DIM==3.  O `static_assert` e' o
// que impede isto de virar suposicao se a unidade for construida com DIM=2.
Point *como_pontos(std::vector<Ponto3> &v)
{
  static_assert(DIM == 3, "mls3d e' o codigo 3D: construa esta unidade com DIM=3");
  return reinterpret_cast<Point *>(v.data());
}

//! O TETO DA BASE PROPRIA: 35 monomios, que e' a ordem 4 em tres dimensoes.
//
// Em duas dimensoes eram 15.  A matriz normal e' `NMAX x NMAX` e resolve-la e'
// cubico nisso, logo o ajuste de ordem 4 em 3D custa `(35/15)^3 = 13` vezes o de
// 2D -- por ponto.  E' um numero a ter em mente antes de escolher a ordem.
constexpr int NMAX = 35;

//! Abaixo desta razao de pivos o ajuste nao merece confianca e a ordem cai.
constexpr double RAZAO_MINIMA = 1e-11;

//! Codigo do funcional LAPLACIANO.  Fora da faixa dos `e_q` de proposito.
constexpr int LAPL = -100;

//! Ate' onde a base propria vai.
constexpr int ORDEM_PROPRIA_MAX = 4;

//! Os monomios de `(dx, dy, dz)` ate' a ordem `ord`, na ordem canonica.
//
// A ORDEM DOS MONOMIOS NAO E' ARBITRARIA: os funcionais sao vetores canonicos
// nesta base, e muda-la muda quais.  Com a ordem abaixo:
//
//     indice 0        o VALOR
//     indices 1,2,3   `d/dx`, `d/dy`, `d/dz`
//     indices 4,7,9   `x^2`, `y^2`, `z^2` -- o LAPLACIANO e' `2(c4 + c7 + c9)`
//
// Reordenar sem acertar esses tres lugares produz um operador que reproduz
// polinomios e calcula outra coisa.
int base_rel(double dx, double dy, double dz, int ord, double *b)
{
  int n = 0;
  b[n++] = 1.0;
  if (ord >= 1) { b[n++] = dx; b[n++] = dy; b[n++] = dz; }
  if (ord >= 2) {
    b[n++] = dx * dx;   // 4
    b[n++] = dx * dy;   // 5
    b[n++] = dx * dz;   // 6
    b[n++] = dy * dy;   // 7
    b[n++] = dy * dz;   // 8
    b[n++] = dz * dz;   // 9
  }
  if (ord >= 3) {
    b[n++] = dx * dx * dx; b[n++] = dx * dx * dy; b[n++] = dx * dx * dz;
    b[n++] = dx * dy * dy; b[n++] = dx * dy * dz; b[n++] = dx * dz * dz;
    b[n++] = dy * dy * dy; b[n++] = dy * dy * dz; b[n++] = dy * dz * dz;
    b[n++] = dz * dz * dz;
  }
  if (ord >= 4) {
    const double x2 = dx * dx, y2 = dy * dy, z2 = dz * dz;
    b[n++] = x2 * x2;      b[n++] = x2 * dx * dy; b[n++] = x2 * dx * dz;
    b[n++] = x2 * y2;      b[n++] = x2 * dy * dz; b[n++] = x2 * z2;
    b[n++] = dx * dy * y2; b[n++] = dx * y2 * dz; b[n++] = dx * dy * z2;
    b[n++] = dx * z2 * dz; b[n++] = y2 * y2;      b[n++] = y2 * dy * dz;
    b[n++] = y2 * z2;      b[n++] = dy * z2 * dz; b[n++] = z2 * z2;
  }
  return n;
}

//! Resolve `M y = rhs` por Gauss com pivotamento parcial.
//
// DEVOLVE A RAZAO `pivo_min / pivo_max`, e nao um `bool`: recusar apenas a
// singularidade exata deixa passar o ajuste que degradou, e devolver zeros
// nesse caso produz uma linha de matriz silenciosamente errada.
double resolve_pequeno(double M[NMAX][NMAX], double *rhs, int n)
{
  double pmax = 0.0, pmin = 1e300;
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int i = c + 1; i < n; i++)
      if (std::fabs(M[i][c]) > std::fabs(M[piv][c])) piv = i;
    const double ap = std::fabs(M[piv][c]);
    pmax = std::max(pmax, ap);
    pmin = std::min(pmin, ap);
    if (ap < 1e-300) return 0.0;
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
  return (pmax > 0.0) ? pmin / pmax : 0.0;
}

std::vector<double> pesos_de_1(const Ponto3 &alvo, const std::vector<Ponto3> &sup,
                               int ord, int qual);

//! Com rebaixamento de ordem registrado.
std::vector<double> pesos_de(const Ponto3 &alvo, const std::vector<Ponto3> &sup,
                             int ord, int qual, int *ord_usada)
{
  // O PISO DO REBAIXAMENTO DEPENDE DO FUNCIONAL: o valor sobrevive ate' a ordem
  // 0 (media ponderada), a derivada precisa de 1, e o laplaciano precisa de 2.
  const int piso = (qual == LAPL) ? 2 : 0;
  for (int o = ord; o >= piso; o--) {
    const std::vector<double> w = pesos_de_1(alvo, sup, o, qual);
    if (!w.empty()) {
      if (ord_usada) *ord_usada = o;
      return w;
    }
  }
  if (ord_usada) *ord_usada = 0;
  return std::vector<double>(sup.size(), 0.0);
}

//! Uma tentativa, numa ordem so'.  Devolve VAZIO quando o ajuste e' indigno.
std::vector<double> pesos_de_1(const Ponto3 &alvo, const std::vector<Ponto3> &sup,
                               int ord, int qual)
{
  const int n = static_cast<int>(sup.size());
  double M[NMAX][NMAX] = {};
  double b[NMAX];
  int nmon = 0;

  // A BASE E' ADIMENSIONALIZADA pelo raio do suporte.  Sem isso o condicionamento
  // da matriz normal depende de `h` e o ajuste de ordem alta degrada em silencio.
  double L = 0.0;
  for (int k = 0; k < n; k++) {
    double d2 = 0.0;
    for (int c = 0; c < 3; c++) {
      const double t = sup[static_cast<size_t>(k)][c] - alvo[c];
      d2 += t * t;
    }
    L = std::max(L, std::sqrt(d2));
  }
  if (L <= 0.0) L = 1.0;   // suporte degenerado num ponto: escala neutra

  // O REGULARIZADOR DO PESO, diferente para o LAPLACIANO.
  //
  // O peso herdado e' `1/(d^2 + EPSMACH)`.  Quando o alvo COINCIDE com um ponto
  // do suporte -- o que acontece sempre no laplaciano pontual, porque o alvo e'
  // o centro da celula -- esse ponto recebe peso da ordem de 1e32 e o ajuste
  // passa a INTERPOLAR: a matriz normal fica efetivamente de posto 1, os pivos
  // despencam, e o ajuste REBAIXA em silencio.
  //
  // Para o funcional VALOR esse comportamento e' o certo; para a segunda
  // derivada e' fatal.  O regularizador troca interpolacao por APROXIMACAO.
  const double TAU = 0.25;
  const double reg = (qual == LAPL) ? (TAU * L) * (TAU * L) : EPSMACH;

  std::vector<double> p2(static_cast<size_t>(n));
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    const double dz = sup[static_cast<size_t>(k)][2] - alvo[2];
    nmon = base_rel(dx / L, dy / L, dz / L, ord, b);
    // O PESO FICA NA ESCALA FISICA: e' o peso do herdado, e mexer nele mudaria o
    // ajuste em vez de so' recondiciona-lo -- e a comparacao com o nucleo
    // perderia o sentido.
    p2[static_cast<size_t>(k)] = 1.0 / (dx * dx + dy * dy + dz * dz + reg);
    for (int i = 0; i < nmon; i++)
      for (int j = 0; j < nmon; j++) M[i][j] += p2[static_cast<size_t>(k)] * b[i] * b[j];
  }

  // O LADO DIREITO E' O VETOR DE AVALIACAO.  O LAPLACIANO e' o unico que nao e'
  // um `e_q`: ele e' `2(c_xx + c_yy + c_zz)`, e em tres dimensoes esses
  // coeficientes estao nas posicoes 4, 7 e 9 da base -- nao 3 e 5 como em duas.
  double y[NMAX] = {};
  if (qual == LAPL) {
    if (nmon < 10) return {};        // ordem < 2 nao tem segunda derivada
    y[4] = 2.0;
    y[7] = 2.0;
    y[9] = 2.0;
  } else {
    if (qual < 0 || qual >= nmon) return {};
    y[qual] = 1.0;
  }
  if (resolve_pequeno(M, y, nmon) < RAZAO_MINIMA) return {};

  std::vector<double> w(static_cast<size_t>(n), 0.0);
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    const double dz = sup[static_cast<size_t>(k)][2] - alvo[2];
    base_rel(dx / L, dy / L, dz / L, ord, b);
    double sv = 0.0;
    for (int i = 0; i < nmon; i++) sv += y[i] * b[i];
    w[static_cast<size_t>(k)] = p2[static_cast<size_t>(k)] * sv;
  }

  // DESFAZ A ESCALA: a primeira derivada e' de ordem 1 em `1/L`, e o laplaciano,
  // que e' de segunda ordem, e' de `1/L^2`.  Esquecer o expoente daria um
  // laplaciano com fator errado que AINDA reproduz a constante e AINDA soma
  // zero -- as duas conferencias baratas nao o pegariam.
  if (qual == LAPL) {
    for (int k = 0; k < n; k++) w[static_cast<size_t>(k)] /= (L * L);
  } else if (qual >= 1) {
    for (int k = 0; k < n; k++) w[static_cast<size_t>(k)] /= L;
  }
  return w;
}

}  // namespace

Mls3::Mls3(int ordem, int max_pontos) : ordem_(ordem), max_pontos_(max_pontos)
{
  if (ordem < 1) throw std::invalid_argument("Mls3: ordem tem de ser >= 1");
  if (max_pontos < min_pontos(ordem))
    throw std::invalid_argument(
        "Mls3: max_pontos (" + std::to_string(max_pontos) + ") abaixo do minimo " +
        std::to_string(min_pontos(ordem)) + " para ordem " + std::to_string(ordem) +
        " em TRES dimensoes -- a base tem 4, 10, 20 e 35 monomios nas ordens 1 a 4");

  wls_ = wls_create(3, ordem, max_pontos);
  if (wls_ == nullptr) throw std::runtime_error("Mls3: wls_create devolveu nulo");
}

Mls3::~Mls3()
{
  if (wls_ != nullptr) wls_destroy(static_cast<wls_interpolator *>(wls_));
}

int Mls3::min_pontos(int ordem) { return wls_num_min_points(3, ordem); }

std::vector<double> Mls3::pesos(const Ponto3 &alvo,
                                const std::vector<Ponto3> &suporte) const
{
  const int n = static_cast<int>(suporte.size());
  if (n < min_pontos(ordem_))
    throw std::invalid_argument(
        "Mls3::pesos: suporte de " + std::to_string(n) + " pontos, minimo " +
        std::to_string(min_pontos(ordem_)) + " para ordem " + std::to_string(ordem_));
  if (n > max_pontos_)
    throw std::invalid_argument(
        "Mls3::pesos: suporte de " + std::to_string(n) + " pontos acima do teto " +
        std::to_string(max_pontos_) + " -- o herdado escreve fora do buffer nisso");

  std::vector<Ponto3> sup = suporte;   // o herdado recebe ponteiro nao-const
  Ponto3 x = alvo;
  std::vector<double> w(static_cast<size_t>(n), 0.0);
  wls_set_points_and_calc(static_cast<wls_interpolator *>(wls_), n, como_pontos(sup),
                          reinterpret_cast<real *>(x.data()), w.data());
  return w;
}

std::vector<double> Mls3::pesos_proprios(const Ponto3 &alvo,
                                         const std::vector<Ponto3> &suporte) const
{
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls3::pesos_proprios: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls3::pesos_proprios: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return conta(pesos_de(alvo, suporte, ordem_, 0, &ord_usada_));
}

std::vector<double> Mls3::pesos_derivada(const Ponto3 &alvo,
                                         const std::vector<Ponto3> &suporte,
                                         int dir) const
{
  if (dir < 0 || dir > 2)
    throw std::invalid_argument("Mls3::pesos_derivada: direcao " + std::to_string(dir) +
                                " -- em 3D ha' 0 (x), 1 (y) e 2 (z)");
  if (ordem_ < 1)
    throw std::invalid_argument("Mls3::pesos_derivada: ordem 0 nao tem derivada");
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls3::pesos_derivada: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls3::pesos_derivada: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return conta(pesos_de(alvo, suporte, ordem_, 1 + dir, &ord_usada_));
}

std::vector<double> Mls3::pesos_laplaciano(const Ponto3 &alvo,
                                           const std::vector<Ponto3> &suporte) const
{
  if (ordem_ < 2)
    throw std::invalid_argument(
        "Mls3::pesos_laplaciano: ordem " + std::to_string(ordem_) +
        " -- um ajuste de ordem menor que 2 nao tem segunda derivada");
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls3::pesos_laplaciano: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls3::pesos_laplaciano: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return conta(pesos_de(alvo, suporte, ordem_, LAPL, &ord_usada_));
}

double Mls3::interpola(const Ponto3 &alvo, const std::vector<Ponto3> &suporte,
                       const std::vector<double> &valores) const
{
  if (valores.size() != suporte.size())
    throw std::invalid_argument(
        "Mls3::interpola: " + std::to_string(valores.size()) + " valores para " +
        std::to_string(suporte.size()) + " pontos de suporte");
  const std::vector<double> w = pesos(alvo, suporte);
  double s = 0.0;
  for (size_t k = 0; k < w.size(); k++) s += w[k] * valores[k];
  return s;
}

const std::vector<double> &Mls3::conta(const std::vector<double> &w) const
{
  ajustes_++;
  if (ord_usada_ < ordem_) rebaixos_++;
  if (ord_usada_ == 0) nulos_++;
  return w;
}

}  // namespace t8flow3d
