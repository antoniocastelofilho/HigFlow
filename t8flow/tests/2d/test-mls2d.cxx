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

//! Campo linear puro, para o contraste da derivada.
static double grau1_so(const Ponto2 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1]; }

static double grau3(const Ponto2 &p)
{
  return grau2(p) + 4.0 * p[0] * p[0] * p[0] - 2.5 * p[1] * p[1] * p[1];
}

//! Grau 4 e grau 5, para o teto da base propria e o contraste acima dele.
//
// Coeficientes distintos em cada monomio, pela mesma razao de sempre: um ajuste que
// troque `x` por `y` ou perca um termo nao passa por acaso.
static double grau4(const Ponto2 &p)
{
  const double x = p[0], y = p[1];
  return grau3(p) + 1.7 * x * x * x * x - 0.9 * x * x * x * y + 0.4 * x * x * y * y
         + 1.1 * x * y * y * y - 2.2 * y * y * y * y;
}

static double dgrau4_dx(const Ponto2 &p)
{
  const double x = p[0], y = p[1];
  // d/dx de grau2 + d/dx dos cubicos + d/dx dos quarticos, a mao -- nao por
  // diferenca numerica, que traria o erro que o caso mede.
  return 2.0 + x + 0.25 * y                                   // de grau2
         + 12.0 * x * x                                       // de 4x^3
         + 6.8 * x * x * x - 2.7 * x * x * y + 0.8 * x * y * y + 1.1 * y * y * y;
}

static double grau5(const Ponto2 &p)
{
  return grau4(p) + 0.8 * p[0] * p[0] * p[0] * p[0] * p[0];
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

  // ---- PESOS DE DERIVADA -----------------------------------------------------
  //
  // A peca que faltava para reconstruir fluxo de alta ordem em malha adaptada.  A
  // algebra e' a mesma do valor com outro vetor de avaliacao, e os casos abaixo
  // afirmam as duas coisas que isso exige: que a derivada e' EXATA ate' a ordem do
  // ajuste, e que a implementacao PROPRIA concorda com o nucleo herdado no caso de
  // valor -- sem a segunda, nao ha' oraculo para a primeira.

  t_case("os_pesos_proprios_de_valor_concordam_com_o_nucleo_herdado");
  {
    // DUAS IMPLEMENTACOES DO MESMO AJUSTE: o herdado resolve por QR, a propria por
    // equacoes normais.  Se divergirem, a derivada -- que so' existe na propria --
    // fica sem verificacao possivel.  A tolerancia e' o patamar de dois algoritmos
    // diferentes, o mesmo 1e-10 que o lote de estenceis em GPU usa.
    Mls m2(2, static_cast<int>(pts.size()));
    const std::vector<double> wh = m2.pesos(alvo, pts);
    const std::vector<double> wp = m2.pesos_proprios(alvo, pts);
    T_CHECK_MSG(wh.size() == wp.size(), "%zu pesos herdados, %zu proprios",
                wh.size(), wp.size());
    double pior = 0.0;
    for (size_t k = 0; k < wh.size() && k < wp.size(); k++)
      pior = std::fmax(pior, std::fabs(wh[k] - wp[k]));
    T_CHECK_MSG(pior < 1e-10,
                "pior diferenca de peso entre o herdado (QR) e o proprio (equacoes "
                "normais): %.3e", pior);
    T_CHECK_MSG(pior > 0.0,
                "diferenca EXATAMENTE zero entre dois algoritmos diferentes -- "
                "suspeito: um dos dois caminhos provavelmente nao rodou");
    printf("  (pesos de valor: herdado contra proprio, pior diferenca %.3e)\n", pior);
  }

  t_case("a_derivada_e_exata_para_polinomio_do_grau_da_ordem");
  {
    // Ordem 2: a derivada tem de ser exata para campo de grau 2.  As derivadas
    // analiticas de `grau2` sao `2 + x` e `-1,2y + 0,25x`... (ver as funcoes no
    // alto deste arquivo) -- escritas aqui explicitamente para o caso nao depender
    // de diferenciacao numerica, que introduziria o erro que ele mede.
    Mls m2(2, static_cast<int>(pts.size()));
    const std::vector<double> wx = m2.pesos_derivada(alvo, pts, 0);
    const std::vector<double> wy = m2.pesos_derivada(alvo, pts, 1);

    const std::vector<double> v = amostra(pts, grau2);
    double dx = 0.0, dy = 0.0;
    for (size_t k = 0; k < v.size(); k++) { dx += wx[k] * v[k]; dy += wy[k] * v[k]; }

    // grau2(p) = 1 + 2x - 3y + 0,5x^2 - 0,6y^2 + 0,25xy
    const double ex = 2.0 + 1.0 * alvo[0] + 0.25 * alvo[1];          // d/dx
    const double ey = -3.0 - 1.2 * alvo[1] + 0.25 * alvo[0];         // d/dy
    T_NEAR(dx, ex, 1e-9, "du/dx de um campo de grau 2 com ajuste de ordem 2");
    T_NEAR(dy, ey, 1e-9, "du/dy de um campo de grau 2 com ajuste de ordem 2");
  }

  t_case("a_derivada_erra_no_grau_acima_da_ordem");
  {
    // CONTRASTE: sem ele, "a derivada esta' certa" nao distingue ajuste de
    // coincidencia -- e a derivada de um campo linear sai certa em quase qualquer
    // esquema.
    Mls m1(1, static_cast<int>(pts.size()));
    const std::vector<double> w1 = m1.pesos_derivada(alvo, pts, 0);
    const std::vector<double> v1 = amostra(pts, grau1_so);
    const std::vector<double> v2 = amostra(pts, grau2);
    double d1 = 0.0, d2 = 0.0;
    for (size_t k = 0; k < v1.size(); k++) { d1 += w1[k] * v1[k]; d2 += w1[k] * v2[k]; }
    const double e1 = std::fabs(d1 - 2.0);                            // d/dx de grau1
    const double e2 = std::fabs(d2 - (2.0 + 1.0 * alvo[0] + 0.25 * alvo[1]));
    T_NEAR(e1, 0.0, 1e-10, "du/dx de campo linear com ajuste de ordem 1");
    T_CONTRASTE(e1 > 1e-14 ? e1 : 1e-14, e2, 1e3,
                "derivada de grau 1 contra grau 2, ajuste de ordem 1");
  }

  t_case("os_pesos_de_derivada_somam_zero");
  {
    // A derivada de uma CONSTANTE e' zero, logo os pesos somam zero -- a conferencia
    // barata que pega erro de escala e de sinal, independente de qualquer campo.
    // (Os pesos de VALOR somam 1; estes somam 0.)
    Mls m2(2, static_cast<int>(pts.size()));
    for (int dir = 0; dir < 2; dir++) {
      const std::vector<double> w = m2.pesos_derivada(alvo, pts, dir);
      double soma = 0.0;
      for (double x : w) soma += x;
      T_NEAR(soma, 0.0, 1e-9, dir == 0 ? "soma dos pesos de d/dx" : "soma dos pesos de d/dy");
    }
  }

  t_case("derivada_malformada_e_recusada");
  {
    Mls m2(2, static_cast<int>(pts.size()));
    bool lancou = false;
    try { m2.pesos_derivada(alvo, pts, 2); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "direcao 2 passou sem reclamar");

    lancou = false;
    try {
      std::vector<Ponto2> poucos(pts.begin(), pts.begin() + 3);
      m2.pesos_derivada(alvo, poucos, 0);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "suporte insuficiente para derivada passou sem reclamar");
  }

  // -------------------------------------------------------------------------
  // D1, D2, D3 -- a ordem 4, a escala da base e o rebaixamento registrado
  // -------------------------------------------------------------------------

  t_case("a_derivada_de_ordem_quatro_e_exata_no_grau_quatro");
  {
    // Ordem 4 e' o TETO da base propria (15 monomios).  Preciso de nuvem maior:
    // 6x6 = 36 pontos, com folga sobre os 15.
    const std::vector<Ponto2> p36 = grade(6, 0.2);
    Mls m4(4, static_cast<int>(p36.size()));
    const std::vector<double> w = m4.pesos_derivada(alvo, p36, 0);

    double d4 = 0.0, d5 = 0.0;
    const std::vector<double> v4 = amostra(p36, grau4);
    const std::vector<double> v5 = amostra(p36, grau5);
    for (size_t k = 0; k < w.size(); k++) { d4 += w[k] * v4[k]; d5 += w[k] * v5[k]; }

    const double e4 = std::fabs(d4 - dgrau4_dx(alvo));
    // o grau 5 acrescenta 0,8 x^5, cuja derivada em x e' 4 x^4
    const double e5 = std::fabs(d5 - (dgrau4_dx(alvo) + 4.0 * std::pow(alvo[0], 4)));
    T_NEAR(e4, 0.0, 1e-9, "du/dx de campo de grau 4 com ajuste de ordem 4");
    T_CONTRASTE(e4 > 1e-14 ? e4 : 1e-14, e5, 1e3,
                "grau 4 contra grau 5 no mesmo ajuste de ordem 4");
    T_CHECK_MSG(m4.rebaixamentos() == 0,
                "a nuvem 6x6 nao deveria exigir rebaixamento de ordem");
  }

  t_case("a_derivada_de_ordem_quatro_sobrevive_a_suporte_minusculo");
  {
    // D2.  A MESMA conta num suporte mil vezes menor.  Em base NAO adimensionalizada
    // a matriz normal da ordem 4 aqui tem condicionamento acima de 1e24 e o ajuste
    // perde todos os digitos -- devolvendo pesos, nao erro.  Com a escala, o
    // condicionamento nao depende de `h`, e o erro relativo e' o mesmo de cima.
    const double h = 1e-3;
    std::vector<Ponto2> p36;
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) p36.push_back({i * h, j * h});
    const Ponto2 alvo_p = {2.4 * h, 3.1 * h};

    Mls m4(4, static_cast<int>(p36.size()));
    const std::vector<double> w = m4.pesos_derivada(alvo_p, p36, 0);
    double d4 = 0.0;
    const std::vector<double> v4 = amostra(p36, grau4);
    for (size_t k = 0; k < w.size(); k++) d4 += w[k] * v4[k];

    // o erro ADMISSIVEL escala com a derivada, que aqui vale ~2; `1e-6` absoluto
    // ainda e' seis ordens abaixo dela, e um ajuste colapsado erra na primeira casa
    T_NEAR(d4, dgrau4_dx(alvo_p), 1e-6, "du/dx de grau 4 em suporte de h = 1e-3");
    T_CHECK_MSG(m4.rebaixamentos() == 0,
                "suporte pequeno nao deveria exigir rebaixamento: a base e' escalada");
  }

  t_case("o_rebaixamento_de_ordem_e_contado_e_os_pesos_continuam_validos");
  {
    // D3.  Suporte QUASE colinear: zigue-zague de `1e-6` em torno de uma reta.
    // A ordem 1 e' sustentada (a direcao transversal tem informacao de 1e-4); a
    // ordem 2 nao e' (os monomios em `y` caem para 1e-16 depois da escala, e a
    // razao de pivos fica abaixo do limiar).  O ajuste honesto aqui e' o linear.
    //
    // O desvio foi CALIBRADO: com `1e-6` nem a ordem 1 sobrevive -- a razao de
    // pivos da' 4e-12, abaixo do limiar de 1e-11 -- e o caso media a recusa.
    //
    // QUASE e nao exatamente colinear de proposito: numa reta EXATA nem a ordem 1
    // e' determinada em 2D, e o caso mediria a recusa em vez do rebaixamento.
    // Quase-degenerado e' tambem o que a malha adaptada produz de verdade.
    //
    // O que este caso proibe e' o comportamento ANTERIOR: devolver um vetor de
    // ZEROS e seguir.  Em malha adaptada com ordem 4 isso produziu uma linha da
    // matriz silenciosamente errada -- erro de 3,3e-01 num nivel cercado por 1e-10.
    std::vector<Ponto2> reta;
    for (int i = 0; i < 20; i++)
      reta.push_back({0.05 * i, 0.3 + ((i % 2) ? 1e-4 : -1e-4)});
    const Ponto2 alvo_r = {0.47, 0.3};

    Mls m3(3, static_cast<int>(reta.size()));
    const std::vector<double> w = m3.pesos_derivada(alvo_r, reta, 0);

    T_CHECK_MSG(m3.rebaixamentos() == 1,
                "suporte colinear tinha de contar UM rebaixamento de ordem");
    T_CHECK_MSG(m3.ajustes_nulos() == 0,
                "nenhuma ordem serviu -- o rebaixamento devia ter achado uma");

    // OS PESOS CONTINUAM VALIDOS: a derivada ao longo da reta sai certa, e e' o
    // unico que um suporte colinear pode dar.  Campo linear em x sobre a reta.
    double soma = 0.0, dx = 0.0;
    for (size_t k = 0; k < w.size(); k++) {
      soma += w[k];
      dx += w[k] * grau1_so(reta[k]);
    }
    T_NEAR(soma, 0.0, 1e-6, "pesos de derivada do ajuste rebaixado somam zero");
    T_NEAR(dx, 2.0, 1e-4, "du/dx de campo linear com ajuste rebaixado");
  }

  return t_end();
}
