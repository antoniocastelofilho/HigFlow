#include "manufaturada3d.hxx"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace t8flow3d {

namespace {

// Os coeficientes sao distintos em cada monomio de proposito -- ver o cabecalho.
//
//   g1 = 1 + 2x - 3y + 1.5z
//   g2 = g1 + 0.5x^2 - 0.6y^2 + 0.9z^2 + 0.25xy - 0.4xz + 0.7yz
//   g3 = g2 + 4x^3 - 2.5y^3 + 1.3z^3 + 0.8xyz
//   g4 = g3 + 1.7x^4 - 2.2y^4 + 0.6z^4
//
// E os laplacianos, derivados a mao:
//
//   lap(g1) = 0
//   lap(g2) = 2(0.5) + 2(-0.6) + 2(0.9) = 1.6
//   lap(g3) = lap(g2) + 24x - 15y + 7.8z        (o termo xyz tem laplaciano 0)
//   lap(g4) = lap(g3) + 20.4x^2 - 26.4y^2 + 7.2z^2

double g1(const Ponto3 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 1.5 * p[2]; }

double g2(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return g1(p) + 0.5 * x * x - 0.6 * y * y + 0.9 * z * z
         + 0.25 * x * y - 0.4 * x * z + 0.7 * y * z;
}

double g3(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return g2(p) + 4.0 * x * x * x - 2.5 * y * y * y + 1.3 * z * z * z
         + 0.8 * x * y * z;
}

double g4(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return g3(p) + 1.7 * x * x * x * x - 2.2 * y * y * y * y + 0.6 * z * z * z * z;
}

//! TRIGONOMETRICA, nao polinomial: e' nela que ha' ordem de convergencia a
//! medir, porque num polinomio do grau da exatidao o erro e' de maquina e a
//! razao entre dois niveis divide ruido por ruido.
double trig(const Ponto3 &p)
{
  return std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]) * std::sin(M_PI * p[2]);
}

std::vector<Manufaturada> faz_catalogo()
{
  std::vector<Manufaturada> c;

  {
    Manufaturada m;
    m.nome = "grau1";
    m.descricao = "linear: todo esquema consistente e' exato";
    m.grau = 1;
    m.u = g1;
    m.dudx = [](const Ponto3 &) { return 2.0; };
    m.dudy = [](const Ponto3 &) { return -3.0; };
    m.dudz = [](const Ponto3 &) { return 1.5; };
    m.f = [](const Ponto3 &) { return 0.0; };
    c.push_back(m);
  }
  {
    Manufaturada m;
    m.nome = "grau2";
    m.descricao = "quadratica com todos os termos mistos";
    m.grau = 2;
    m.u = g2;
    m.dudx = [](const Ponto3 &p) { return 2.0 + p[0] + 0.25 * p[1] - 0.4 * p[2]; };
    m.dudy = [](const Ponto3 &p) { return -3.0 - 1.2 * p[1] + 0.25 * p[0] + 0.7 * p[2]; };
    m.dudz = [](const Ponto3 &p) { return 1.5 + 1.8 * p[2] - 0.4 * p[0] + 0.7 * p[1]; };
    m.f = [](const Ponto3 &) { return 1.6; };
    c.push_back(m);
  }
  {
    Manufaturada m;
    m.nome = "grau3";
    m.descricao = "cubica, com o termo xyz que so' existe em 3D";
    m.grau = 3;
    m.u = g3;
    m.dudx = [](const Ponto3 &p) {
      return 2.0 + p[0] + 0.25 * p[1] - 0.4 * p[2] + 12.0 * p[0] * p[0]
             + 0.8 * p[1] * p[2];
    };
    m.dudy = [](const Ponto3 &p) {
      return -3.0 - 1.2 * p[1] + 0.25 * p[0] + 0.7 * p[2] - 7.5 * p[1] * p[1]
             + 0.8 * p[0] * p[2];
    };
    m.dudz = [](const Ponto3 &p) {
      return 1.5 + 1.8 * p[2] - 0.4 * p[0] + 0.7 * p[1] + 3.9 * p[2] * p[2]
             + 0.8 * p[0] * p[1];
    };
    m.f = [](const Ponto3 &p) {
      return 1.6 + 24.0 * p[0] - 15.0 * p[1] + 7.8 * p[2];
    };
    c.push_back(m);
  }
  {
    Manufaturada m;
    m.nome = "grau4";
    m.descricao = "quartica: acima da exatidao da ordem 3";
    m.grau = 4;
    m.u = g4;
    m.dudx = [](const Ponto3 &p) {
      return 2.0 + p[0] + 0.25 * p[1] - 0.4 * p[2] + 12.0 * p[0] * p[0]
             + 0.8 * p[1] * p[2] + 6.8 * p[0] * p[0] * p[0];
    };
    m.dudy = [](const Ponto3 &p) {
      return -3.0 - 1.2 * p[1] + 0.25 * p[0] + 0.7 * p[2] - 7.5 * p[1] * p[1]
             + 0.8 * p[0] * p[2] - 8.8 * p[1] * p[1] * p[1];
    };
    m.dudz = [](const Ponto3 &p) {
      return 1.5 + 1.8 * p[2] - 0.4 * p[0] + 0.7 * p[1] + 3.9 * p[2] * p[2]
             + 0.8 * p[0] * p[1] + 2.4 * p[2] * p[2] * p[2];
    };
    m.f = [](const Ponto3 &p) {
      return 1.6 + 24.0 * p[0] - 15.0 * p[1] + 7.8 * p[2]
             + 20.4 * p[0] * p[0] - 26.4 * p[1] * p[1] + 7.2 * p[2] * p[2];
    };
    c.push_back(m);
  }
  {
    Manufaturada m;
    m.nome = "trig";
    m.descricao = "sin(pi x) cos(pi y) sin(pi z): nao polinomial";
    m.grau = -1;
    m.u = trig;
    m.dudx = [](const Ponto3 &p) {
      return M_PI * std::cos(M_PI * p[0]) * std::cos(M_PI * p[1]) * std::sin(M_PI * p[2]);
    };
    m.dudy = [](const Ponto3 &p) {
      return -M_PI * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]) * std::sin(M_PI * p[2]);
    };
    m.dudz = [](const Ponto3 &p) {
      return M_PI * std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]) * std::cos(M_PI * p[2]);
    };
    // lap(u) = -3 pi^2 u, porque cada uma das tres derivadas segundas traz um
    // fator `-pi^2`.  EM DUAS DIMENSOES O FATOR E' -2 pi^2, e traduzir aquele
    // numero daria uma fonte errada por 50% -- com a solucao ainda convergindo,
    // para a equacao errada.
    m.f = [](const Ponto3 &p) { return -3.0 * M_PI * M_PI * trig(p); };
    c.push_back(m);
  }

  return c;
}

const std::vector<Manufaturada> &catalogo()
{
  static const std::vector<Manufaturada> c = faz_catalogo();
  return c;
}

}  // namespace

std::string nomes_manufaturadas()
{
  std::string s;
  for (const Manufaturada &m : catalogo()) {
    if (!s.empty()) s += ", ";
    s += m.nome;
  }
  return s;
}

const Manufaturada &manufaturada(const std::string &nome)
{
  for (const Manufaturada &m : catalogo())
    if (m.nome == nome) return m;
  // A LISTA VAI NA MENSAGEM: um erro que diz apenas "nao existe" obriga quem o
  // recebe a ler o codigo.
  throw std::invalid_argument("manufaturada: \"" + nome +
                              "\" nao existe.  Validas: " + nomes_manufaturadas());
}

}  // namespace t8flow3d
