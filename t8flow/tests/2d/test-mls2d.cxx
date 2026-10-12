// O DISCRETIZADOR 2D: o que ele reproduz, o que ele NAO reproduz, e o que ele
// recusa.
//
// ORACULO: reproducao polinomial COM CONTRASTE.  Um ajuste de ordem k reproduz
// exatamente polinomio de grau <= k.  Sozinho isso nao discrimina quase nada --
// nao pega franja ausente, nem parede de fechamento errada com campo linear, nem
// suporte deslocado --, entao cada caso de reproducao vem com um campo de grau
// ACIMA da ordem, cujo erro tem de ser ordens de grandeza maior.  Se os dois
// derem o mesmo erro, nao houve reproducao: houve coincidencia, e o caso reprova
// dizendo isso.
#include "../testing.h"

#include <array>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "mls2d.hxx"

using t8flow2d::Mls;
using t8flow2d::Ponto2;

// Coeficientes DISTINTOS por direcao de proposito: assim o valor identifica a
// posicao, e um suporte que troque os eixos nao passa por acaso.
static double grau2(const Ponto2 &p)
{
  return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 0.5 * p[0] * p[0] - 0.6 * p[1] * p[1]
         + 0.25 * p[0] * p[1];
}

static double grau3(const Ponto2 &p)
{
  return grau2(p) + 4.0 * p[0] * p[0] * p[0] - 2.5 * p[1] * p[1] * p[1];
}

static std::vector<Ponto2> grade(int lado, double h)
{
  std::vector<Ponto2> v;
  for (int i = 0; i < lado; i++)
    for (int j = 0; j < lado; j++) v.push_back({i * h, j * h});
  return v;
}

static std::vector<double> amostra(const std::vector<Ponto2> &pts, double (*f)(const Ponto2 &))
{
  std::vector<double> v;
  v.reserve(pts.size());
  for (const Ponto2 &p : pts) v.push_back(f(p));
  return v;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const std::vector<Ponto2> pts = grade(4, 0.25);

  // PONTO DE CONSULTA ASSIMETRICO, e isto o teste descobriu sozinho: com o alvo
  // no centro exato da nuvem, o ajuste de ordem 2 reproduziu o campo CUBICO com
  // erro 1,1e-16 -- os termos impares cancelam contra pesos simetricos --, e o
  // contraste deu 0,01x em vez de 1000x.  O caso reprovou dizendo "o oraculo nao
  // discrimina", que e' o trabalho dele.
  const Ponto2 alvo = {0.30, 0.44};

  Mls mls(2, static_cast<int>(pts.size()));

  t_case("mls_reproduz_polinomio_do_grau_da_ordem");
  const double obtido2 = mls.interpola(alvo, pts, amostra(pts, grau2));
  const double exato2 = grau2(alvo);
  T_NEAR(obtido2, exato2, 1e-9, "campo de grau 2 com ajuste de ordem 2");

  t_case("mls_erra_mais_no_grau_acima_da_ordem");
  const double obtido3 = mls.interpola(alvo, pts, amostra(pts, grau3));
  const double err2 = std::fabs(obtido2 - exato2);
  const double err3 = std::fabs(obtido3 - grau3(alvo));
  T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
              "grau 2 contra grau 3 no mesmo ajuste de ordem 2");

  // Os pesos de um ajuste que reproduz constante somam 1.  E' barato, e e'
  // independente do campo: pega escala errada nos pesos, que a reproducao de um
  // campo especifico pode absorver.
  t_case("pesos_somam_um");
  const std::vector<double> w = mls.pesos(alvo, pts);
  double soma = 0.0;
  for (double wi : w) soma += wi;
  T_NEAR(soma, 1.0, 1e-10, "soma dos pesos (reproducao da constante)");
  T_CHECK_MSG(w.size() == pts.size(), "%zu pesos para %zu pontos", w.size(), pts.size());

  // RECUSA EM VEZ DE APROXIMAR.  Com suporte abaixo do minimo da ordem o ajuste
  // herdado devolve numeros, e numeros ruins em silencio sao o modo de falha que
  // este projeto mais paga.  A fachada lanca.
  t_case("suporte_insuficiente_e_recusado");
  const int minimo = Mls::min_pontos(2);
  std::vector<Ponto2> poucos(pts.begin(), pts.begin() + (minimo - 1));
  bool lancou = false;
  try {
    mls.pesos(alvo, poucos);
  } catch (const std::invalid_argument &) {
    lancou = true;
  }
  T_CHECK_MSG(lancou, "suporte de %d pontos (minimo %d) passou sem reclamar",
              minimo - 1, minimo);

  t_case("valores_de_tamanho_errado_sao_recusados");
  lancou = false;
  try {
    mls.interpola(alvo, pts, std::vector<double>(pts.size() - 1, 0.0));
  } catch (const std::invalid_argument &) {
    lancou = true;
  }
  T_CHECK_MSG(lancou, "numero de valores diferente do suporte passou sem reclamar");

  return t_end();
}
