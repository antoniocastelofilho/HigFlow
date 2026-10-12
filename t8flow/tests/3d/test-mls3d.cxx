// O ajuste por minimos quadrados moveis em TRES dimensoes.
//
// OS ORACULOS SAO OS QUE DUAS DIMENSOES ENSINOU A TER, mais o que so' tres pede.
// Cada um existe porque um defeito concreto passou por ele:
//
//   a concordancia com o nucleo HERDADO, que resolve o mesmo problema por
//   decomposicao QR -- e que rejeita a concordancia EXATAMENTE nula, porque
//   diferenca zero entre duas implementacoes significa que so' uma rodou;
//
//   o CONTRASTE no grau acima da ordem, sem o qual "reproduz o polinomio" nao
//   distingue ajuste de coincidencia;
//
//   o ALVO ASSIMETRICO, porque num alvo centrado os termos impares cancelam
//   contra pesos simetricos e um ajuste de ordem menor passa pelo de ordem
//   maior;
//
//   o regularizador do peso no LAPLACIANO, com o alvo DENTRO do suporte -- que
//   e' onde o laplaciano pontual sempre o coloca.
//
// E o que e' de tres dimensoes: a base tem 4, 10, 20 e 35 monomios, e o
// laplaciano soma TRES segundas derivadas e nao duas.

#include "../testing.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "mls3d.hxx"

using namespace t8flow3d;

namespace {

// COEFICIENTES DISTINTOS EM CADA MONOMIO, e em tres dimensoes isso importa
// mais: com `x`, `y` e `z` ha' mais maneiras de trocar dois eixos, e um ajuste
// que os troque passaria em qualquer campo simetrico.
double grau1(const Ponto3 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 1.5 * p[2]; }

double grau2(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return grau1(p) + 0.5 * x * x - 0.6 * y * y + 0.9 * z * z
         + 0.25 * x * y - 0.4 * x * z + 0.7 * y * z;
}

double grau3(const Ponto3 &p)
{
  const double x = p[0], y = p[1], z = p[2];
  return grau2(p) + 4.0 * x * x * x - 2.5 * y * y * y + 1.3 * z * z * z
         + 0.8 * x * y * z;
}

//! `lap(grau2)` e' CONSTANTE: `2*(0,5) + 2*(-0,6) + 2*(0,9)`.
double lap_grau2() { return 2.0 * (0.5 - 0.6 + 0.9); }

//! `lap(grau3)` acrescenta `24x - 15y + 7,8z`; o termo `xyz` tem laplaciano nulo.
double lap_grau3(const Ponto3 &p)
{
  return lap_grau2() + 24.0 * p[0] - 15.0 * p[1] + 7.8 * p[2];
}

std::vector<Ponto3> grade(int lado, double h)
{
  std::vector<Ponto3> v;
  for (int i = 0; i < lado; i++)
    for (int j = 0; j < lado; j++)
      for (int k = 0; k < lado; k++) v.push_back({i * h, j * h, k * h});
  return v;
}

std::vector<double> amostra(const std::vector<Ponto3> &pts, double (*f)(const Ponto3 &))
{
  std::vector<double> v;
  v.reserve(pts.size());
  for (const Ponto3 &p : pts) v.push_back(f(p));
  return v;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // 4x4x4 = 64 pontos, que cobre a ordem 3 (20 monomios) com folga.
  const std::vector<Ponto3> pts = grade(4, 0.25);

  // ALVO ASSIMETRICO.  Num alvo centrado na nuvem os termos impares cancelam
  // contra pesos simetricos, e o ajuste de ordem 2 reproduz o campo CUBICO --
  // em duas dimensoes isso fez o contraste dar 0,01x em vez de 1000x, e o caso
  // reprovou dizendo "o oraculo nao discrimina".
  const Ponto3 alvo = {0.30, 0.44, 0.19};

  t_case("a_base_de_tres_dimensoes_tem_4_10_20_e_35_monomios");
  {
    // A CONTAGEM E' A CONFERENCIA MAIS BARATA de que a base e' a de tres
    // dimensoes: a de duas daria 3, 6, 10 e 15.
    const int esperado[4] = {4, 10, 20, 35};
    for (int o = 1; o <= 4; o++)
      T_CHECK_MSG(Mls3::min_pontos(o) == esperado[o - 1],
                  "ordem %d: minimo %d pontos, esperado %d", o, Mls3::min_pontos(o),
                  esperado[o - 1]);

    // E SUPORTE ABAIXO DO MINIMO E' RECUSADO na construcao, com o numero na
    // mensagem -- o herdado escreve fora do buffer se isso passar.
    bool lancou = false;
    try { Mls3 m(3, 10); } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "ordem 3 com teto de 10 pontos passou; ela pede 20");
  }

  t_case("os_pesos_proprios_de_valor_concordam_com_o_nucleo_herdado");
  {
    // DOIS CAMINHOS PARA O MESMO PROBLEMA: o herdado resolve por decomposicao
    // QR, o proprio pela matriz normal.  A concordancia e' um oraculo que
    // nenhum dos dois produz sozinho.
    Mls3 m(2, static_cast<int>(pts.size()));
    const std::vector<double> wh = m.pesos(alvo, pts);
    const std::vector<double> wp = m.pesos_proprios(alvo, pts);

    double pior = 0.0;
    for (size_t k = 0; k < wh.size(); k++)
      pior = std::fmax(pior, std::fabs(wh[k] - wp[k]));

    T_CHECK_MSG(pior < 1e-12, "pior diferenca entre herdado e proprio: %.3e", pior);
    // E A CONCORDANCIA EXATAMENTE NULA E' SUSPEITA: ela significaria que os dois
    // caminhos sao o mesmo codigo, e o oraculo nao estaria comparando nada.
    T_CHECK_MSG(pior > 0.0,
                "a diferenca deu EXATAMENTE zero -- dois caminhos distintos nao "
                "coincidem bit a bit, e isto sugere que so' um deles rodou");
    printf("  (pesos de valor: herdado contra proprio, pior diferenca %.3e)\n", pior);
  }

  t_case("o_ajuste_reproduz_o_polinomio_do_grau_da_ordem_e_erra_acima");
  {
    Mls3 m2(2, static_cast<int>(pts.size()));
    const double v2 = m2.interpola(alvo, pts, amostra(pts, grau2));
    const double v3 = m2.interpola(alvo, pts, amostra(pts, grau3));
    const double e2 = std::fabs(v2 - grau2(alvo));
    const double e3 = std::fabs(v3 - grau3(alvo));

    T_NEAR(e2, 0.0, 1e-10, "campo de grau 2 com ajuste de ordem 2");
    // O CONTRASTE: sem ele, "reproduz" nao distingue ajuste de coincidencia.
    T_CONTRASTE(e2 > 1e-14 ? e2 : 1e-14, e3, 1e3,
                "grau 2 contra grau 3 no mesmo ajuste de ordem 2");

    // E a ordem 3 reproduz o cubico, que e' a outra metade da afirmacao.
    Mls3 m3(3, static_cast<int>(pts.size()));
    T_NEAR(m3.interpola(alvo, pts, amostra(pts, grau3)), grau3(alvo), 1e-9,
           "campo de grau 3 com ajuste de ordem 3");
  }

  t_case("as_tres_derivadas_sao_exatas_e_os_pesos_somam_zero");
  {
    // TRES DIRECOES, e pedir a quarta e' recusado.  Em duas dimensoes eram duas,
    // e um laco que nao acompanhe a dimensao passaria calado.
    Mls3 m(2, static_cast<int>(pts.size()));
    const std::vector<double> v = amostra(pts, grau2);

    // grau2 = grau1 + 0,5x^2 - 0,6y^2 + 0,9z^2 + 0,25xy - 0,4xz + 0,7yz
    const double x = alvo[0], y = alvo[1], z = alvo[2];
    const double esperado[3] = {
        2.0 + 1.0 * x + 0.25 * y - 0.4 * z,     // d/dx
        -3.0 - 1.2 * y + 0.25 * x + 0.7 * z,    // d/dy
        1.5 + 1.8 * z - 0.4 * x + 0.7 * y,      // d/dz
    };

    for (int dir = 0; dir < 3; dir++) {
      const std::vector<double> w = m.pesos_derivada(alvo, pts, dir);
      double soma = 0.0, der = 0.0;
      for (size_t k = 0; k < w.size(); k++) { soma += w[k]; der += w[k] * v[k]; }

      // OS PESOS DE DERIVADA SOMAM ZERO -- a derivada de uma constante.
      // Conferencia barata, independente do campo, que pega erro de escala e de
      // sinal que a reproducao de um campo especifico pode absorver.
      T_NEAR(soma, 0.0, 1e-9, "pesos de derivada somam zero");
      T_NEAR(der, esperado[dir], 1e-9, "derivada na direcao do campo de grau 2");
    }

    bool lancou = false;
    try { m.pesos_derivada(alvo, pts, 3); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "a direcao 3 passou num codigo de tres dimensoes");
  }

  t_case("o_laplaciano_soma_as_TRES_segundas_derivadas");
  {
    // EM TRES DIMENSOES O LAPLACIANO TEM TRES TERMOS, e o erro de somar so' dois
    // -- o que uma traducao descuidada do codigo 2D produziria -- daria um
    // operador que ainda soma zero nos pesos e ainda reproduz a constante.
    // Nenhuma das duas conferencias baratas o pegaria; o que o pega e' o valor.
    Mls3 m2(2, static_cast<int>(pts.size()));
    const std::vector<double> w = m2.pesos_laplaciano(alvo, pts);

    double soma = 0.0, l1 = 0.0, l2 = 0.0, l3 = 0.0;
    const std::vector<double> v1 = amostra(pts, grau1);
    const std::vector<double> v2 = amostra(pts, grau2);
    const std::vector<double> v3 = amostra(pts, grau3);
    for (size_t k = 0; k < w.size(); k++) {
      soma += w[k];
      l1 += w[k] * v1[k];
      l2 += w[k] * v2[k];
      l3 += w[k] * v3[k];
    }

    T_NEAR(soma, 0.0, 1e-8, "pesos do laplaciano somam zero");
    T_NEAR(l1, 0.0, 1e-8, "laplaciano de campo LINEAR e' zero");
    T_NEAR(l2, lap_grau2(), 1e-8, "laplaciano de campo de grau 2");

    // CONTRASTE: a ordem 2 nao reproduz o laplaciano de um cubico.
    const double erro3 = std::fabs(l3 - lap_grau3(alvo));
    T_CONTRASTE(std::fabs(l2 - lap_grau2()) > 1e-14 ? std::fabs(l2 - lap_grau2()) : 1e-14,
                erro3, 1e2, "laplaciano: grau 2 contra grau 3 no ajuste de ordem 2");

    // E a ordem 3 reproduz o cubico.
    Mls3 m3(3, static_cast<int>(pts.size()));
    const std::vector<double> w3 = m3.pesos_laplaciano(alvo, pts);
    double c3 = 0.0;
    for (size_t k = 0; k < w3.size(); k++) c3 += w3[k] * v3[k];
    T_NEAR(c3, lap_grau3(alvo), 1e-7, "laplaciano de grau 3 com ajuste de ordem 3");

    // ORDEM 1 E' RECUSADA: um ajuste linear nao tem segunda derivada.
    bool lancou = false;
    try {
      Mls3 m1(1, static_cast<int>(pts.size()));
      m1.pesos_laplaciano(alvo, pts);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "ordem 1 passou no laplaciano; nao ha' segunda derivada ali");
  }

  t_case("o_peso_do_laplaciano_e_regularizado_quando_o_alvo_esta_no_suporte");
  {
    // ONDE O LAPLACIANO PONTUAL SEMPRE COLOCA O ALVO: no centro da celula, que
    // ESTA' no suporte.  Com o peso herdado `1/(d^2 + EPSMACH)` esse ponto
    // recebe peso da ordem de 1e32 e o ajuste passa a INTERPOLAR -- a matriz
    // normal fica de posto 1, os pivos despencam, e o ajuste REBAIXA em
    // silencio, dando o resultado da ordem 2 onde se pediu 3.
    const double h = 0.1;
    const Ponto3 centro = {0.5, 0.5, 0.5};
    std::vector<Ponto3> sup;
    for (int i = -2; i <= 2; i++)
      for (int j = -2; j <= 2; j++)
        for (int k = -2; k <= 2; k++)
          sup.push_back({centro[0] + i * h, centro[1] + j * h, centro[2] + k * h});

    bool achou = false;
    for (const Ponto3 &p : sup)
      if (std::fabs(p[0] - centro[0]) < 1e-15 && std::fabs(p[1] - centro[1]) < 1e-15 &&
          std::fabs(p[2] - centro[2]) < 1e-15)
        achou = true;
    T_CHECK_MSG(achou, "o caso nao mede nada: o alvo nao esta' no suporte");

    for (int ord : {2, 3}) {
      Mls3 m(ord, static_cast<int>(sup.size()));
      const std::vector<double> w = m.pesos_laplaciano(centro, sup);

      // NENHUM REBAIXAMENTO e' a afirmacao central.  Com o peso singular, as
      // ordens acima de 2 rebaixavam aqui.
      T_CHECK_MSG(m.rebaixamentos() == 0,
                  "ordem %d rebaixou com o alvo dentro do suporte -- o peso nao "
                  "esta' regularizado", ord);

      double lap = 0.0;
      const std::vector<double> v2 = amostra(sup, grau2);
      for (size_t k = 0; k < w.size(); k++) lap += w[k] * v2[k];
      T_NEAR(lap, lap_grau2(), 1e-7, "laplaciano de grau 2 com o alvo no suporte");
    }

    // E AQUI A MEDIDA CORRIGIU UMA EXPECTATIVA MINHA.  Eu esperava que a ordem
    // 3 reproduzisse o laplaciano de um cubico e a ordem 2 nao -- e as duas
    // reproduzem, em 4e-13.
    //
    // A RAZAO E' REAL E VALE REGISTRAR: num suporte SIMETRICO em torno do alvo,
    // os monomios IMPARES nao contribuem para o laplaciano no centro.  O termo
    // `x^3` tem laplaciano `6x`, que no alvo vale `6*x_alvo` -- e o ajuste de
    // ordem 2, que nao tem `x^3`, acerta assim mesmo porque a contribuicao
    // impar se cancela contra pesos simetricos.
    //
    // ISSO E' SUPERCONVERGENCIA DO ESTENCIL SIMETRICO, e nao um defeito: em
    // malha uniforme o centro de uma celula ve' exatamente esse suporte, logo a
    // ordem 2 ja' entrega grau 3 no laplaciano.  Medido igual em duas dimensoes.
    //
    // O caso afirma a propriedade em vez de esperar o contraste que nao existe.
    Mls3 m2(2, static_cast<int>(sup.size()));
    Mls3 m3(3, static_cast<int>(sup.size()));
    const std::vector<double> w2 = m2.pesos_laplaciano(centro, sup);
    const std::vector<double> w3 = m3.pesos_laplaciano(centro, sup);
    const std::vector<double> v3 = amostra(sup, grau3);
    double a2 = 0.0, a3 = 0.0;
    for (size_t k = 0; k < w2.size(); k++) { a2 += w2[k] * v3[k]; a3 += w3[k] * v3[k]; }
    const double e3 = std::fabs(a3 - lap_grau3(centro));
    const double e2 = std::fabs(a2 - lap_grau3(centro));
    T_NEAR(e3, 0.0, 1e-6, "laplaciano de grau 3 com ajuste de ordem 3");
    T_NEAR(e2, 0.0, 1e-6,
           "laplaciano de grau 3 com ajuste de ordem 2 em suporte SIMETRICO -- os "
           "monomios impares nao contribuem no centro, logo a ordem 2 ja' basta");

    // E O CONTRASTE EXISTE, so' que fora da simetria: deslocando o alvo para
    // dentro do suporte mas FORA do centro, os termos impares passam a
    // contribuir e a ordem 2 erra.  Sem esta metade, o caso acima poderia passar
    // com um laplaciano que ignora a ordem.
    const Ponto3 desloc = {centro[0] + 0.5 * h, centro[1] + 0.3 * h, centro[2]};
    Mls3 d2(2, static_cast<int>(sup.size()));
    Mls3 d3(3, static_cast<int>(sup.size()));
    const std::vector<double> u2 = d2.pesos_laplaciano(desloc, sup);
    const std::vector<double> u3 = d3.pesos_laplaciano(desloc, sup);
    double b2 = 0.0, b3 = 0.0;
    for (size_t k = 0; k < u2.size(); k++) { b2 += u2[k] * v3[k]; b3 += u3[k] * v3[k]; }
    const double f3 = std::fabs(b3 - lap_grau3(desloc));
    const double f2 = std::fabs(b2 - lap_grau3(desloc));
    T_NEAR(f3, 0.0, 1e-6, "laplaciano de grau 3 com ordem 3, alvo deslocado");
    T_CONTRASTE(f3 > 1e-14 ? f3 : 1e-14, f2, 1e2,
                "laplaciano de grau 3 fora da simetria: ordem 3 contra ordem 2");

    printf("  (alvo no suporte: no centro ordem 2 e 3 empatam em %.1e e %.1e "
           "(simetria); deslocado, %.1e contra %.1e)\n", e3, e2, f3, f2);
  }

  t_case("o_rebaixamento_de_ordem_e_contado_e_os_pesos_continuam_validos");
  {
    // SUPORTE QUASE PLANO: todos os pontos quase no mesmo plano `z`.  A direcao
    // `z` tem informacao de apenas `1e-4`, logo os monomios em `z` de grau 2
    // caem para `1e-8` depois da escala e a ordem 2 nao e' sustentada; a ordem 1
    // e'.
    //
    // QUASE e nao exatamente plano de proposito: num plano EXATO nem a ordem 1 e'
    // determinada em tres dimensoes, e o caso mediria a recusa em vez do
    // rebaixamento.  Quase-degenerado e' tambem o que a malha adaptada produz.
    std::vector<Ponto3> plano;
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
        plano.push_back({0.1 * i, 0.1 * j, 0.3 + (((i + j) % 2) ? 1e-4 : -1e-4)});
    const Ponto3 alvo_p = {0.27, 0.33, 0.3};

    Mls3 m(2, static_cast<int>(plano.size()));
    const std::vector<double> w = m.pesos_derivada(alvo_p, plano, 0);

    T_CHECK_MSG(m.rebaixamentos() == 1,
                "suporte quase plano tinha de contar UM rebaixamento, contou %ld",
                m.rebaixamentos());
    T_CHECK_MSG(m.ajustes_nulos() == 0,
                "nenhuma ordem serviu -- o rebaixamento devia ter achado uma");

    // OS PESOS CONTINUAM VALIDOS: a derivada no plano sai certa, e e' o unico
    // que um suporte quase plano pode dar.  O que este caso PROIBE e' o
    // comportamento anterior -- devolver um vetor de zeros e seguir.
    double soma = 0.0, dx = 0.0;
    for (size_t k = 0; k < w.size(); k++) {
      soma += w[k];
      dx += w[k] * grau1(plano[k]);
    }
    T_NEAR(soma, 0.0, 1e-6, "pesos de derivada do ajuste rebaixado somam zero");
    T_NEAR(dx, 2.0, 1e-4, "du/dx de campo linear com ajuste rebaixado");
  }

  return t_end();
}
