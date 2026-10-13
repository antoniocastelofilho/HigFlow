// O NUCLEO GEOMETRICO DO PLIC: as duas vias, e a ida e volta.
//
// O ORACULO PRINCIPAL NAO E' DE ORDEM, E' DE EXATIDAO.  A area cortada de um
// retangulo por uma reta e' um poligono, e a area de um poligono sai EXATA da
// formula do sapateiro -- logo a forma fechada tem de concordar com o recorte
// em precisao de maquina, em qualquer normal e qualquer razao de aspecto.  Erro
// de ordem certa com constante errada nao existe aqui: ou bate, ou ha' defeito.
//
// E CADA AFIRMACAO DE CONCORDANCIA VEM COM OS REGIMES ONDE UMA DELAS ERRARIA:
// normal quase alinhada a um eixo (onde a reducao canonica divide por `a`
// pequeno), celula muito achatada, e area perto de 0 e de 1 (onde os ramos da
// quadratica se emendam).  Varrer so' o caso bonito -- normal a 45 graus em
// celula quadrada com area 1/2 -- passaria com metade dos ramos errados.

#include "../testing.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "arranque.hxx"
#include "vof2d.hxx"

using namespace t8flow2d;

namespace {

//! As celulas da varredura: quadrada, achatada nos dois sentidos, e uma razao
//! que nao e' potencia de dois (as potencias cancelam erros de arredondamento
//! que uma razao qualquer expoe).
const Ponto2 CELULAS[] = {{1.0, 1.0},
                          {0.0625, 0.0625},
                          {1.0, 0.125},
                          {0.125, 1.0},
                          {0.3, 0.7}};

//! As normais: os quatro eixos, as quatro diagonais, duas QUASE alinhadas (o
//! regime em que a reducao canonica tem `a` perto de zero) e duas obliquas
//! quaisquer.
const Ponto2 NORMAIS[] = {
    {1.0, 0.0},   {-1.0, 0.0},  {0.0, 1.0},   {0.0, -1.0},
    {1.0, 1.0},   {-1.0, 1.0},  {1.0, -1.0},  {-1.0, -1.0},
    {1.0, 1e-9},  {-1e-9, 1.0}, {1.0, 1e-14}, {1e-14, -1.0},
    {0.73, -0.41}, {-0.17, 0.92}};

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("a_area_concorda_pelas_duas_vias");
  {
    // As fracoes-alvo incluem 0 e 1 (os patamares) e valores PERTO deles, que
    // e' onde os ramos `t <= a` e `t >= b` da quadratica vivem.
    const double FRACS[] = {0.0,  1e-12, 1e-6, 0.01, 0.1,  0.25, 0.5,
                            0.75, 0.9,   0.99, 1.0 - 1e-6, 1.0};

    double pior = 0.0, pior_rel = 0.0;
    long conta = 0, cortadas = 0;
    Ponto2 pior_n{0.0, 0.0}, pior_h{0.0, 0.0};
    double pior_d = 0.0;

    for (const Ponto2 &h : CELULAS)
      for (const Ponto2 &n : NORMAIS) {
        const double vol = h[0] * h[1];
        const double nn = std::hypot(n[0], n[1]);
        const double S = (std::fabs(n[0]) / nn) * h[0]
                         + (std::fabs(n[1]) / nn) * h[1];
        for (const double f : FRACS) {
          // Varre `d` DIRETO, e nao pela inversa: assim este caso mede a ida
          // sozinha, sem depender da volta.
          const double d = -S / 2.0 + f * S;
          const double a1 = plic_area(n, h, d);
          const double a2 = plic_area_por_recorte(n, h, d);
          conta++;
          if (a1 > 1e-13 * vol && a1 < vol * (1.0 - 1e-13)) cortadas++;
          const double e = std::fabs(a1 - a2);
          if (e / vol > pior_rel) {
            pior_rel = e / vol;
            pior = e;
            pior_n = n; pior_h = h; pior_d = d;
          }
        }
      }

    // A GUARDA DE VACUIDADE.  Se nenhuma combinacao cortasse a celula de
    // verdade -- todas no patamar 0 ou no patamar cheio --, as duas vias
    // concordariam trivialmente em zero e o caso nao afirmaria nada.
    T_CHECK_MSG(cortadas > conta / 3,
                "so' %ld de %ld combinacoes cortam a celula de verdade -- a "
                "varredura ficou nos patamares e o caso nao mede nada",
                cortadas, conta);
    T_CHECK_MSG(pior_rel <= 1e-14,
                "as duas vias discordam em %.3e (relativo %.3e) com "
                "n = (%.3g, %.3g), h = (%.3g, %.3g), d = %.6g",
                pior, pior_rel, pior_n[0], pior_n[1], pior_h[0], pior_h[1],
                pior_d);
  }

  t_case("a_distancia_e_a_inversa_exata_da_area");
  {
    // A IDA E VOLTA.  Mais forte do que parece: a area e' estritamente
    // decrescente em `d` no intervalo ativo, logo a inversa e' UNICA.  Quem
    // reproduz a area do HiGFlow e inverte exatamente TEM de dar o `d` dele --
    // sem precisar linkar o codigo dele para comparar.
    const double FRACS[] = {1e-9, 1e-4, 0.01, 0.123, 0.25, 0.5, 0.5001,
                            0.75, 0.877, 0.99, 1.0 - 1e-4, 1.0 - 1e-9};
    double pior = 0.0;
    Ponto2 pior_n{0.0, 0.0}, pior_h{0.0, 0.0};
    double pior_f = 0.0;

    for (const Ponto2 &h : CELULAS)
      for (const Ponto2 &n : NORMAIS)
        for (const double f : FRACS) {
          const double alvo = f * h[0] * h[1];
          const double d = plic_distancia(n, h, alvo);
          const double volta = plic_area(n, h, d);
          const double e = std::fabs(volta - alvo) / (h[0] * h[1]);
          if (e > pior) { pior = e; pior_n = n; pior_h = h; pior_f = f; }
        }
    T_CHECK_MSG(pior <= 1e-14,
                "a ida e volta erra %.3e com n = (%.3g, %.3g), h = (%.3g, "
                "%.3g), fracao %.6g", pior, pior_n[0], pior_n[1], pior_h[0],
                pior_h[1], pior_f);
  }

  t_case("a_distancia_concorda_com_a_bisseccao");
  {
    // A VOLTA PELAS DUAS VIAS.  A ida e volta acima nao pega um erro que esteja
    // nas DUAS -- area e distancia partilham a reducao canonica.  A bisseccao
    // sobre o recorte nao partilha algebra nenhuma com elas.
    //
    // A TOLERANCIA AQUI E' DA BISSECCAO, nao da forma fechada: 200 bisseccoes
    // sobre um intervalo de tamanho `S` chegam ao ultimo bit, mas a comparacao
    // justa e' em unidades de `S`.
    const double FRACS[] = {0.02, 0.137, 0.3, 0.5, 0.64, 0.91, 0.998};
    double pior = 0.0;
    Ponto2 pior_n{0.0, 0.0}, pior_h{0.0, 0.0};

    for (const Ponto2 &h : CELULAS)
      for (const Ponto2 &n : NORMAIS)
        for (const double f : FRACS) {
          const double alvo = f * h[0] * h[1];
          const double nn = std::hypot(n[0], n[1]);
          const double S = (std::fabs(n[0]) / nn) * h[0]
                           + (std::fabs(n[1]) / nn) * h[1];
          const double d1 = plic_distancia(n, h, alvo);
          const double d2 = plic_distancia_por_bisseccao(n, h, alvo);
          const double e = std::fabs(d1 - d2) / S;
          if (e > pior) { pior = e; pior_n = n; pior_h = h; }
        }
    T_CHECK_MSG(pior <= 1e-13,
                "a forma fechada e a bisseccao discordam em %.3e de S, com "
                "n = (%.3g, %.3g), h = (%.3g, %.3g)",
                pior, pior_n[0], pior_n[1], pior_h[0], pior_h[1]);
  }

  t_case("os_casos_que_tem_area_conhecida_a_mao");
  {
    // CONTRA NUMERO CONHECIDO, e nao so' contra a outra via: duas vias podem
    // concordar na MESMA convencao errada -- o lado trocado, por exemplo, daria
    // concordancia perfeita e area complementar.
    const Ponto2 q{1.0, 1.0};
    // Reta pelo centro: metade, em qualquer normal.
    T_NEAR(plic_area({0.0, 1.0}, q, 0.0), 0.5, 1e-15, "reta pelo centro, n = +y");
    T_NEAR(plic_area({1.0, 1.0}, q, 0.0), 0.5, 1e-15, "reta pelo centro, n na diagonal");
    // A NORMAL APONTA PARA A REGIAO CONTADA: com `n = +y` e `d = +0.3`, a
    // regiao e' `y >= 0.3`, de area 0,2.  Com o lado trocado daria 0,8.
    T_NEAR(plic_area({0.0, 1.0}, q, 0.3), 0.2, 1e-15, "n = +y, d = 0,3: a regiao e' y >= 0,3");
    T_NEAR(plic_area({0.0, -1.0}, q, 0.3), 0.2, 1e-15, "n = -y, d = 0,3: a regiao e' y <= -0,3");
    // Triangulo de canto: `x + y >= t` tem area `(1-t)^2/2`; com `t = 0.5` da'
    // 0,125, e `d = t/sqrt(2)`.
    T_NEAR(plic_area({1.0, 1.0}, q, 0.5 / std::sqrt(2.0)), 0.125, 1e-15, "triangulo de canto");
    // Fora do suporte: patamares exatos, e nao "quase".
    T_NEAR(plic_area({1.0, 1.0}, q, 10.0), 0.0, 0.0, "patamar vazio, EXATO");
    T_NEAR(plic_area({1.0, 1.0}, q, -10.0), 1.0, 0.0, "patamar cheio, EXATO");
  }

  t_case("a_normal_nula_e_a_area_fora_do_intervalo_sao_RECUSADAS");
  {
    // NORMAL NULA NAO TEM LADO.  Aceitar calado daria area arbitraria numa
    // celula de interface -- e e' justamente o que o estimador de normal
    // devolve numa celula PURA, onde `grad f` e' zero.  Quem chamar com ela
    // tem de ser informado, nao servido.
    bool lancou = false;
    try { plic_area({0.0, 0.0}, {1.0, 1.0}, 0.0); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "normal nula foi aceita");

    lancou = false;
    try { plic_distancia({1.0, 1.0}, {1.0, 1.0}, 0.0); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(!lancou, "area ZERO e' legitima e foi recusada");

    // E A AREA FORA DE [0, vol] LANCA.  A fracao sai do intervalo em advecao
    // direcional partida -- e' conhecido --, e saturar AQUI apagaria a unica
    // evidencia de que saiu.  O ajuste e' de quem advecta, e contado.
    lancou = false;
    try { plic_distancia({1.0, 1.0}, {1.0, 1.0}, 1.0000001); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "area acima do volume da celula foi aceita");

    lancou = false;
    try { plic_distancia({1.0, 1.0}, {1.0, 1.0}, -1e-12); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "area negativa foi aceita");
  }

  // =========================================================================
  // A FRACAO VOLUMETRICA INICIAL
  // =========================================================================

  t_case("a_fracao_inicial_e_EXATA_para_interface_plana");
  {
    // O ORACULO MAIS FORTE DISPONIVEL, e ele e' cruzado: com `forma` AFIM a
    // interface dentro da celula E' uma reta, logo a marcha de quadrados nao
    // aproxima nada -- tem de dar o mesmo que `plic_area`, que chega la' por
    // caminho completamente diferente (reducao canonica, nenhum canto
    // amostrado).
    //
    // Duas implementacoes independentes concordando num caso onde as duas sao
    // exatas afirma bem mais que qualquer teste de ordem.
    double pior = 0.0;
    long conta = 0, cortadas = 0;
    for (const Ponto2 &h : CELULAS)
      for (const Ponto2 &n : NORMAIS) {
        const double nn = std::hypot(n[0], n[1]);
        const double ux = n[0] / nn, uy = n[1] / nn;
        const double S = std::fabs(ux) * h[0] + std::fabs(uy) * h[1];
        for (double f = 0.05; f < 0.96; f += 0.1) {
          const double d = -S / 2.0 + f * S;
          // `forma > 0` e' dentro, e `plic_area` conta `{n.x >= d}`: a forma
          // que corresponde e' `n.x - d`.
          auto forma = [&](const Ponto2 &x) {
            return ux * x[0] + uy * x[1] - d;
          };
          const double a1 = plic_area(n, h, d);
          const double a2 = area_da_forma(forma, {0.0, 0.0}, h);
          conta++;
          if (a1 > 1e-12 && a1 < h[0] * h[1] * (1.0 - 1e-12)) cortadas++;
          pior = std::max(pior, std::fabs(a1 - a2) / (h[0] * h[1]));
        }
      }
    T_CHECK_MSG(cortadas > conta / 2,
                "so' %ld de %ld combinacoes cortam a celula -- nada medido",
                cortadas, conta);
    T_CHECK_MSG(pior <= 1e-14,
                "a marcha de quadrados e o PLIC discordam em %.3e numa "
                "interface PLANA, onde as duas sao exatas", pior);
  }

  t_case("a_subdivisao_e_o_que_faz_o_volume_inicial_ficar_certo");
  {
    // O PONTO DESTE CASO E' O CONTRASTE, e nao o valor: `n = 16` da' area
    // melhor que `n = 1`, e por quanto.  Sem medir os dois, "a subdivisao
    // existe" seria decoracao -- e e' justamente ela que o envoltorio do
    // HiGFlow acrescenta por cima de `get_fracvolN`.
    //
    // O circulo de Hysing, na malha de Hysing.
    const double R = 0.25;
    const Ponto2 c{0.5, 0.5};
    auto circulo = [&](const Ponto2 &x) {
      return R - std::hypot(x[0] - c[0], x[1] - c[1]);   // > 0 dentro
    };
    const double exata = M_PI * R * R;

    for (const double hh : {0.025, 0.00625}) {
      // Varre uma caixa que contem o circulo inteiro, celula a celula.
      auto soma_com = [&](int n) {
        double s = 0.0;
        const int k = static_cast<int>(std::lround(1.0 / hh));
        for (int i = 0; i < k; i++)
          for (int j = 0; j < k; j++) {
            const Ponto2 cc{(i + 0.5) * hh, (j + 0.5) * hh};
            s += fracao_inicial(circulo, cc, {hh, hh}, n) * hh * hh;
          }
        return s;
      };
      const double e1 = std::fabs(soma_com(1) - exata) / exata;
      const double e16 = std::fabs(soma_com(16) - exata) / exata;

      T_CHECK_MSG(e16 < e1 / 20.0,
                  "h = %.5g: a subdivisao 16x16 melhorou so' de %.3e para "
                  "%.3e -- esperado mais de 20 vezes", hh, e1, e16);

      // O LIMIAR VEM DE UM NUMERO MEDIDO, e nao de um que eu escolhi.  O FT do
      // HiGFlow deriva 2,44e-03 de massa ate' `t = 3` neste caso; o erro de
      // INICIALIZACAO so' e' irrelevante se ficar muito abaixo disso.  Duas
      // ordens de margem e' o criterio.
      //
      // A primeira versao deste caso exigia 1e-9, numero que eu nao derivei de
      // nada -- e ele falhou.  Medido: 8,1e-06 em `h = 0,025` e 5,1e-07 em
      // `h = 0,00625`, caindo como `h^2`, que e' o que a conta preve (erro de
      // linearizacao `O((h/n)^3/R)` vezes `O(R/(h/n))` sub-celulas cortadas).
      const double DERIVA_DO_METODO = 2.44e-3;
      T_CHECK_MSG(e16 < DERIVA_DO_METODO / 100.0,
                  "h = %.5g: com 16x16 a area do circulo erra %.3e, que nao "
                  "esta' duas ordens abaixo da deriva do proprio metodo (%.3e)",
                  hh, e16, DERIVA_DO_METODO);
    }
  }

  t_case("o_caso_ambiguo_da_diagonal_e_RECUSADO");
  {
    // Dois cantos dentro NA DIAGONAL: a celula tanto pode ter duas quinas
    // separadas quanto uma garganta que as liga, e os quatro cantos nao
    // distinguem.  Dar um numero aqui seria ESCOLHER, e a escolha ficaria
    // invisivel.
    auto sela = [](const Ponto2 &x) { return x[0] * x[1]; };   // > 0 nos
    bool lancou = false;                                       // quadrantes 1 e 3
    try { area_da_forma(sela, {0.0, 0.0}, {1.0, 1.0}); }
    catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "o caso ambiguo da diagonal devolveu numero");

    // CONTRASTE: o par de ARESTA tem resposta, e ela e' a obvia.
    auto meio = [](const Ponto2 &x) { return -x[1]; };          // y < 0 dentro
    T_NEAR(area_da_forma(meio, {0.0, 0.0}, {1.0, 1.0}), 0.5, 1e-15,
           "metade de baixo");
  }

  return t_end();
}
