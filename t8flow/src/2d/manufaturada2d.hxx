// SOLUCOES MANUFATURADAS 2D: o catalogo, para trocar a funcao sem recompilar.
//
// O metodo das solucoes manufaturadas e' escolher `u` conhecido, DERIVAR o termo
// de fonte `f = lap(u)` e os dados de contorno a partir dele, e medir o erro
// contra o `u` que se escolheu.  A forca do metodo esta' em `f` e nos dados serem
// EXATOS -- derivados a mao, nao calculados numericamente.  Por isso cada solucao
// aqui carrega as suas derivadas, e nao ha' diferenciacao automatica: derivada
// numerica de `f` reintroduziria no oraculo o erro que ele existe para medir.
//
// COMO AS SOLUCOES FORAM ESCOLHIDAS, e cada uma responde a uma pergunta:
//
//   grau1   lap = 0.  Qualquer esquema consistente acerta; serve de piso -- se ela
//           falhar, o defeito nao esta' na ordem do esquema, esta' na montagem.
//   grau2   lap = constante.  E' o caso do pedido: o esquema de balanco com
//           fechamento de TRES pontos e' exato, logo o erro tem de ser de MAQUINA.
//           Com fechamento de dois pontos NAO e', e esse contraste e' o oraculo do
//           fechamento.
//   grau3   lap linear.  O interior (diferenca central) ainda e' exato, o
//           fechamento de tres pontos NAO -- entao o erro deixa de ser de maquina
//           sem o esquema estar errado.  MEDIDO: ordem 2,90 -> 2,95, ou seja ~3,
//           e nao 2 -- o fechamento de tres pontos e' melhor do que o minimo que
//           se exigiria dele.  Serve para o contraste nao ser confundido com
//           defeito.
//   trig    nao-polinomial.  E' a solucao para MEDIR ORDEM DE CONVERGENCIA: com
//           polinomio exato nao ha' taxa para medir, porque o erro ja' e' zero.
//
// Os coeficientes sao DISTINTOS por direcao em todas elas, de proposito: com campo
// simetrico ou unidirecional, eixo trocado e sinal de normal errado ficam
// invisiveis.
#ifndef T8FLOW_MANUFATURADA2D_HXX
#define T8FLOW_MANUFATURADA2D_HXX

#include <functional>
#include <string>
#include <vector>

#include "mls2d.hxx"

namespace t8flow2d {

//! Uma solucao manufaturada: `u`, as derivadas primeiras, e `f = lap(u)`.
struct Manufaturada {
  std::string nome;
  std::string descricao;
  std::function<double(const Ponto2 &)> u;
  std::function<double(const Ponto2 &)> dudx;
  std::function<double(const Ponto2 &)> dudy;
  std::function<double(const Ponto2 &)> f;      //!< lap(u), derivado a mao
  int grau = -1;                                //!< -1 para nao-polinomial

  //! `du/dn` no ponto, dada a normal exterior -- o dado de Neumann.
  double dudn(const Ponto2 &p, const double n[2]) const
  {
    return dudx(p) * n[0] + dudy(p) * n[1];
  }
};

//! O catalogo, por nome: "grau1", "grau2", "grau3", "trig".
//
// Lanca `std::invalid_argument` com a LISTA dos nomes validos se o nome nao
// existir -- errar o nome e receber silenciosamente a solucao padrao seria medir
// uma coisa e reportar outra.
const Manufaturada &manufaturada(const std::string &nome);

//! Os nomes disponiveis, para o exemplo poder listar.
std::vector<std::string> manufaturadas_disponiveis();

}  // namespace t8flow2d

#endif
