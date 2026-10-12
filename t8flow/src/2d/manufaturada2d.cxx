// O catalogo de solucoes manufaturadas.  A RAZAO de cada uma estar aqui esta' no
// `manufaturada2d.hxx`.  Aqui ficam as contas, e elas sao a parte que tem de estar
// certa: `f` ERRADO faz o oraculo medir a consistencia entre dois erros meus, e
// nao o esquema.
//
// Por isso cada `f` abaixo vem com a derivacao escrita ao lado.  Conferir a conta
// e' barato; descobrir depois que o oraculo estava errado nao e'.
#include "manufaturada2d.hxx"

#include <cmath>
#include <map>
#include <stdexcept>

namespace t8flow2d {

namespace {

const double PI = 3.14159265358979323846;

std::map<std::string, Manufaturada> monta_catalogo()
{
  std::map<std::string, Manufaturada> c;

  // ---- grau 1:  u = 1 + 2x - 3y ------------------------------------------
  //   du/dx = 2,  du/dy = -3,  d2u/dx2 = d2u/dy2 = 0  ->  lap = 0
  {
    Manufaturada m;
    m.nome = "grau1";
    m.descricao = "u = 1 + 2x - 3y   (lap u = 0)";
    m.grau = 1;
    m.u    = [](const Ponto2 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1]; };
    m.dudx = [](const Ponto2 &)  { return 2.0; };
    m.dudy = [](const Ponto2 &)  { return -3.0; };
    m.f    = [](const Ponto2 &)  { return 0.0; };
    c[m.nome] = m;
  }

  // ---- grau 2:  u = 1 + 2x - 3y + 4x^2 - 5y^2 + 6xy ----------------------
  //   du/dx = 2 + 8x + 6y,   du/dy = -3 - 10y + 6x
  //   d2u/dx2 = 8,           d2u/dy2 = -10
  //   lap = 8 - 10 = -2   (constante)
  //
  // O TERMO CRUZADO 6xy esta' aqui de proposito: sem ele, um esquema que trate as
  // duas direcoes de forma desacoplada passaria do mesmo jeito, e o caso nao
  // distinguiria.  E os coeficientes de x^2 e y^2 sao DIFERENTES, senao `lap` daria
  // zero e o termo de fonte nao seria exercitado.
  {
    Manufaturada m;
    m.nome = "grau2";
    m.descricao = "u = 1 + 2x - 3y + 4x^2 - 5y^2 + 6xy   (lap u = -2)";
    m.grau = 2;
    m.u    = [](const Ponto2 &p) {
      return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 4.0 * p[0] * p[0]
             - 5.0 * p[1] * p[1] + 6.0 * p[0] * p[1];
    };
    m.dudx = [](const Ponto2 &p) { return 2.0 + 8.0 * p[0] + 6.0 * p[1]; };
    m.dudy = [](const Ponto2 &p) { return -3.0 - 10.0 * p[1] + 6.0 * p[0]; };
    m.f    = [](const Ponto2 &)  { return -2.0; };
    c[m.nome] = m;
  }

  // ---- grau 3:  u = x^3 - 2y^3 + x^2 y ------------------------------------
  //   du/dx = 3x^2 + 2xy,    du/dy = -6y^2 + x^2
  //   d2u/dx2 = 6x + 2y,     d2u/dy2 = -12y
  //   lap = 6x + 2y - 12y = 6x - 10y   (linear)
  {
    Manufaturada m;
    m.nome = "grau3";
    m.descricao = "u = x^3 - 2y^3 + x^2 y   (lap u = 6x - 10y)";
    m.grau = 3;
    m.u    = [](const Ponto2 &p) {
      return p[0]*p[0]*p[0] - 2.0*p[1]*p[1]*p[1] + p[0]*p[0]*p[1];
    };
    m.dudx = [](const Ponto2 &p) { return 3.0*p[0]*p[0] + 2.0*p[0]*p[1]; };
    m.dudy = [](const Ponto2 &p) { return -6.0*p[1]*p[1] + p[0]*p[0]; };
    m.f    = [](const Ponto2 &p) { return 6.0*p[0] - 10.0*p[1]; };
    c[m.nome] = m;
  }

  // ---- trig:  u = sin(pi x) cos(2 pi y) -----------------------------------
  //   du/dx =  pi cos(pi x) cos(2 pi y)
  //   du/dy = -2pi sin(pi x) sin(2 pi y)
  //   d2u/dx2 = -pi^2 u,   d2u/dy2 = -4pi^2 u
  //   lap = -5 pi^2 u
  //
  // Frequencias DIFERENTES nas duas direcoes (1 e 2), senao a solucao fica
  // simetrica na diagonal e troca de eixo nao aparece.
  {
    Manufaturada m;
    m.nome = "trig";
    m.descricao = "u = sin(pi x) cos(2 pi y)   (lap u = -5 pi^2 u)";
    m.grau = -1;
    m.u    = [](const Ponto2 &p) { return std::sin(PI*p[0]) * std::cos(2.0*PI*p[1]); };
    m.dudx = [](const Ponto2 &p) { return PI * std::cos(PI*p[0]) * std::cos(2.0*PI*p[1]); };
    m.dudy = [](const Ponto2 &p) { return -2.0*PI * std::sin(PI*p[0]) * std::sin(2.0*PI*p[1]); };
    m.f    = [](const Ponto2 &p) {
      return -5.0 * PI * PI * std::sin(PI*p[0]) * std::cos(2.0*PI*p[1]);
    };
    c[m.nome] = m;
  }

  return c;
}

const std::map<std::string, Manufaturada> &catalogo()
{
  static const std::map<std::string, Manufaturada> c = monta_catalogo();
  return c;
}

}  // namespace

const Manufaturada &manufaturada(const std::string &nome)
{
  const auto &c = catalogo();
  const auto it = c.find(nome);
  if (it != c.end()) return it->second;

  // A LISTA VAI NA MENSAGEM.  Errar o nome e receber a solucao padrao em silencio
  // seria medir uma coisa e reportar outra -- e o relatorio sairia com o nome
  // errado, que e' pior que nao sair.
  std::string nomes;
  for (const auto &kv : c) nomes += (nomes.empty() ? "" : ", ") + kv.first;
  throw std::invalid_argument("manufaturada: \"" + nome + "\" nao existe.  Validas: " + nomes);
}

std::vector<std::string> manufaturadas_disponiveis()
{
  std::vector<std::string> v;
  for (const auto &kv : catalogo()) v.push_back(kv.first);
  return v;
}

}  // namespace t8flow2d
