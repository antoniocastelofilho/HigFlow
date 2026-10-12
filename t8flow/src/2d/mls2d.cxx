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
  if (ord >= 4) {
    b[n++] = dx * dx * dx * dx; b[n++] = dx * dx * dx * dy;
    b[n++] = dx * dx * dy * dy; b[n++] = dx * dy * dy * dy;
    b[n++] = dy * dy * dy * dy;
  }
  return n;
}

constexpr int NMAX = 15;   // ordem 4 em 2D: 15 monomios

//! Resolve `M y = rhs` por Gauss com pivotamento parcial.
//
// Com pivotamento porque `M[0][0]` pode ser pequeno quando o suporte e' simetrico
// em torno do alvo, e sem ele o resultado degrada em silencio.
//
// DEVOLVE A RAZAO `pivo_min / pivo_max`, nao um `bool`.  A versao anterior so'
// recusava pivo abaixo de `1e-300` -- ou seja, so' a singularidade exata -- e
// devolvia um vetor de ZEROS quando recusava.  Em malha uniforme isso nunca
// acontecia; em malha ADAPTADA com ordem 4 acontecia, e o efeito era uma linha da
// matriz silenciosamente errada: erro de 3,3e-01 num nivel cercado por 1e-10 nos
// vizinhos.  A razao de pivos nao e' o numero de condicao, mas e' barata, monotona
// nele, e basta para decidir se o ajuste merece confianca.
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

//! O nucleo comum: pesos do funcional `alvo_base` (valor ou derivada).
//! Abaixo desta razao de pivos o ajuste nao merece confianca e a ordem cai.
//
// `1e-11` foi escolhido por MEDIDA, nao por gosto: com a base adimensionalizada um
// suporte saudavel da' razao entre `1e-2` e `1e-5`, e os ajustes que produziam o
// erro de 3,3e-01 ficavam abaixo de `1e-13`.  O limiar fica na faixa vazia entre
// os dois regimes.
constexpr double RAZAO_MINIMA = 1e-11;

//! Codigo do funcional LAPLACIANO em `pesos_de`.  Fora da faixa dos `e_q` de
//! proposito, para que um `qual` valido nunca colida com ele.
constexpr int LAPL = -100;

std::vector<double> pesos_de_1(const Ponto2 &alvo, const std::vector<Ponto2> &sup,
                               int ord, int qual);

std::vector<double> pesos_de(const Ponto2 &alvo, const std::vector<Ponto2> &sup,
                             int ord, int qual, int *ord_usada)
{
  // REBAIXAMENTO REGISTRADO.  Quando o suporte nao sustenta a ordem pedida -- em
  // malha adaptada ele as vezes nao sustenta, porque os vizinhos podem cair quase
  // alinhados --, a saida certa e' ajustar com ordem MENOR, nao devolver
  // zeros.  Quem chama fica sabendo por `ord_usada`, e a `Mls` conta.
  // DESCE ATE' A ORDEM 0, que ainda e' um ajuste honesto para VALOR -- a media
  // ponderada do suporte.  Para DERIVADA a ordem 0 nao tem o funcional, logo o
  // loop termina sem resposta e `ajustes_nulos()` conta: e' o sinal certo, porque
  // derivada nula fingida de derivada verdadeira seria pior que a recusa.
  // O PISO DO REBAIXAMENTO DEPENDE DO FUNCIONAL: o valor sobrevive ate' a ordem
  // 0 (media ponderada), a derivada precisa de 1, e o laplaciano precisa de 2.
  // Descer abaixo do piso produziria um funcional que a base nao tem, e devolver
  // zero ali seria afirmar que a segunda derivada e' nula.
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
std::vector<double> pesos_de_1(const Ponto2 &alvo, const std::vector<Ponto2> &sup,
                               int ord, int qual)
{
  const int n = static_cast<int>(sup.size());
  double M[NMAX][NMAX] = {};
  double b[NMAX];
  int nmon = 0;

  // A BASE E' ADIMENSIONALIZADA pelo raio do suporte.  Sem isso a matriz normal da
  // ordem 4 tem condicionamento da ordem de `h^-8`: com `h = 3e-2` isso e' 1e12, e
  // o ajuste degrada em SILENCIO -- os pesos saem, a reproducao polinomial quase
  // passa, e so' a ordem de convergencia acusa.  Com a escala, o condicionamento
  // nao depende de `h`.
  //
  // O ajuste em variavel escalada `s = d/L` e' o MESMO ajuste; so' a derivada muda
  // de unidade, e por isso o peso de derivada e' dividido por `L` no fim.
  double L = 0.0;
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    L = std::max(L, std::sqrt(dx * dx + dy * dy));
  }
  if (L <= 0.0) L = 1.0;   // suporte degenerado num ponto: escala neutra

  // O PESO E' O MESMO DO HERDADO -- WLSMOVING, `D_i = 1/sqrt(dist^2 + EPSMACH)`,
  // aplicado ao quadrado na matriz normal.  Usar outro peso daria um ajuste
  // diferente que reproduz os mesmos polinomios com pesos diferentes, e a
  // comparacao com o herdado perderia o sentido -- foi exatamente o defeito que o
  // lote de estenceis em GPU teve.
  // O REGULARIZADOR DO PESO, e ele e' diferente para o LAPLACIANO.
  //
  // O peso herdado e' `1/(d^2 + EPSMACH)`, com `EPSMACH` da ordem de `1e-32`.
  // Quando o alvo COINCIDE com um ponto do suporte -- o que acontece sempre no
  // laplaciano pontual, porque o alvo e' o centro da celula e o centro esta' no
  // suporte -- esse ponto recebe peso da ordem de `1e32` e o ajuste passa a
  // INTERPOLAR: ele reproduz o valor ali exatamente, as custas de tudo o mais.
  //
  // Para o funcional VALOR isso e' o comportamento certo, e a clausula C5 o
  // registra: no centro da celula o ajuste e' degenerado, peso 1 na propria
  // celula.  Para a SEGUNDA DERIVADA e' fatal -- a matriz normal fica
  // efetivamente de posto 1, os pivos despencam, e o ajuste REBAIXA para a
  // ordem 2 em silencio.
  //
  // MEDIDO: com o peso herdado, o laplaciano de ordem 3 e 4 num suporte 5x5
  // centrado rebaixava sempre, e dava exatamente o resultado da ordem 2 -- o que
  // apareceu, no Poisson, como ordem 2 de convergencia num esquema que deveria
  // dar 3.
  //
  // O regularizador `(TAU*L)^2` limita o peso do ponto coincidente a
  // `1/(TAU*L)^2`, que e' da ordem do peso dos vizinhos mais proximos.  Troca
  // interpolacao por APROXIMACAO, que e' o que a segunda derivada precisa.
  const double TAU = 0.25;
  const double reg = (qual == LAPL) ? (TAU * L) * (TAU * L) : EPSMACH;

  std::vector<double> p2(static_cast<size_t>(n));
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    nmon = base_rel(dx / L, dy / L, ord, b);
    // O PESO FICA NA ESCALA FISICA, de proposito: e' o peso do herdado, e mexer
    // nele mudaria o ajuste em vez de so' recondiciona-lo.  O regularizador e' a
    // unica licenca, e so' para o laplaciano -- ver acima.
    p2[static_cast<size_t>(k)] = 1.0 / (dx * dx + dy * dy + reg);
    for (int i = 0; i < nmon; i++)
      for (int j = 0; j < nmon; j++) M[i][j] += p2[static_cast<size_t>(k)] * b[i] * b[j];
  }

  // O LADO DIREITO E' O VETOR DE AVALIACAO: `e_0` para valor, `e_1` para `d/dx`,
  // `e_2` para `d/dy` -- porque a base e' relativa ao alvo, e no alvo so'
  // sobrevive o coeficiente do monomio correspondente.
  //
  // O LAPLACIANO (`qual == LAPL`) e' o unico que nao e' um `e_q`: ele e'
  // `2*c_xx + 2*c_yy`, os coeficientes de `dx^2` e `dy^2` na base escalada, que
  // ocupam as posicoes 3 e 5.  O fator 2 vem de `d^2(x^2)/dx^2 = 2`.
  double y[NMAX] = {};
  if (qual == LAPL) {
    if (nmon < 6) return {};        // ordem < 2 nao tem segunda derivada
    y[3] = 2.0;
    y[5] = 2.0;
  } else {
    if (qual < 0 || qual >= nmon) return {};   // a ordem nao tem este funcional
    y[qual] = 1.0;
  }
  if (resolve_pequeno(M, y, nmon) < RAZAO_MINIMA) return {};

  std::vector<double> w(static_cast<size_t>(n), 0.0);
  for (int k = 0; k < n; k++) {
    const double dx = sup[static_cast<size_t>(k)][0] - alvo[0];
    const double dy = sup[static_cast<size_t>(k)][1] - alvo[1];
    base_rel(dx / L, dy / L, ord, b);
    double sv = 0.0;
    for (int i = 0; i < nmon; i++) sv += y[i] * b[i];
    w[static_cast<size_t>(k)] = p2[static_cast<size_t>(k)] * sv;
  }
  // DESFAZ A ESCALA DA VARIAVEL: a primeira derivada e' de ordem 1 em `1/L`, e o
  // laplaciano, que e' de segunda ordem, e' de `1/L^2`.  Esquecer o expoente aqui
  // daria um laplaciano com o fator errado -- que ainda reproduz a constante e
  // ainda soma zero, logo as duas conferencias baratas NAO o pegariam; so' a
  // reproducao de um campo quadratico pega.
  if (qual == LAPL) {
    for (int k = 0; k < n; k++) w[static_cast<size_t>(k)] /= (L * L);
  } else if (qual >= 1) {
    for (int k = 0; k < n; k++) w[static_cast<size_t>(k)] /= L;
  }
  return w;
}

}  // namespace

constexpr int ORDEM_PROPRIA_MAX = 4;   // ate' onde `base_rel` vai

std::vector<double> Mls::pesos_proprios(const Ponto2 &alvo,
                                        const std::vector<Ponto2> &suporte) const
{
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls::pesos_proprios: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls::pesos_proprios: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return conta(pesos_de(alvo, suporte, ordem_, 0, &ord_usada_));
}

std::vector<double> Mls::pesos_derivada(const Ponto2 &alvo,
                                        const std::vector<Ponto2> &suporte, int dir) const
{
  if (dir < 0 || dir > 1)
    throw std::invalid_argument("Mls::pesos_derivada: direcao " + std::to_string(dir) +
                                " -- em 2D sao 0 (x) e 1 (y)");
  if (ordem_ < 1)
    throw std::invalid_argument("Mls::pesos_derivada: ordem 0 nao tem derivada");
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls::pesos_derivada: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls::pesos_derivada: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)) + " para ordem " +
                                std::to_string(ordem_));
  // `e_1` e' `d/dx`, `e_2` e' `d/dy` -- a ordem dos monomios em `base_rel`.
  return conta(pesos_de(alvo, suporte, ordem_, 1 + dir, &ord_usada_));
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

const std::vector<double> &Mls::conta(const std::vector<double> &w) const
{
  ajustes_++;
  if (ord_usada_ < ordem_) rebaixos_++;
  if (ord_usada_ == 0) nulos_++;
  return w;
}

std::vector<double> Mls::pesos_laplaciano(const Ponto2 &alvo,
                                          const std::vector<Ponto2> &suporte) const
{
  if (ordem_ < 2)
    throw std::invalid_argument(
        "Mls::pesos_laplaciano: ordem " + std::to_string(ordem_) +
        " -- um ajuste de ordem menor que 2 nao tem segunda derivada");
  if (ordem_ > ORDEM_PROPRIA_MAX)
    throw std::invalid_argument(
        "Mls::pesos_laplaciano: ordem " + std::to_string(ordem_) + " acima de " +
        std::to_string(ORDEM_PROPRIA_MAX) + ", que e' o limite da base propria");
  if (static_cast<int>(suporte.size()) < min_pontos(ordem_))
    throw std::invalid_argument("Mls::pesos_laplaciano: suporte de " +
                                std::to_string(suporte.size()) + " pontos, minimo " +
                                std::to_string(min_pontos(ordem_)));
  return conta(pesos_de(alvo, suporte, ordem_, LAPL, &ord_usada_));
}

}  // namespace t8flow2d
