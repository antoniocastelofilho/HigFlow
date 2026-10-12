// O ACOPLAMENTO DA FRENTE COM A MALHA: fase B2.
//
// TRES ORACULOS, e os tres sao de forma fechada ou de conservacao exata --
// nenhum depende de uma tolerancia escolhida:
//
//   A PARTICAO DA UNIDADE.  A soma de `fracao x volume` da' a area do
//   poligono em precisao de MAQUINA, e so' e' afirmavel assim porque o recorte
//   e' exato.  Uma fracao por amostragem daria erro de amostragem, e o teste
//   teria de aceitar uma tolerancia -- que e' justamente onde um defeito de
//   conservacao se esconde;
//
//   A FORCA TOTAL DISTRIBUIDA E' NULA.  A tensao superficial nao exerce forca
//   liquida sobre uma curva FECHADA, e o nucleo de Peskin soma um -- logo a
//   soma sobre as facetas tem de dar zero de maquina.  Valor nao nulo
//   significa nucleo que nao soma um, ou marcador cuja forca caiu fora;
//
//   O SALTO DE LAPLACE, `delta p = sigma / R`.  Este e' o unico que envolve o
//   solver, e por isso vem por ultimo: se um dos dois acima falha, o salto
//   nao diz nada.

#include "../testing.h"

#include <cmath>
#include <vector>

#include "acoplamento2d.hxx"
#include "arranque.hxx"
#include "campos2d.hxx"
#include "contorno2d.hxx"
#include "frente2d.hxx"
#include "malha2d.hxx"
#include "navierstokes2d.hxx"

using namespace t8flow2d;

namespace {
constexpr double PI = 3.14159265358979323846;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("o_recorte_de_poligono_por_caixa_e_exato");
  {
    // UM QUADRADO recortado por caixas de posicao conhecida: as areas tem
    // forma fechada e sao conferidas uma a uma.
    const std::vector<Ponto2> quad = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    T_NEAR(area_recortada(quad, {0.0, 0.0}, {1.0, 1.0}), 1.0, 1e-15,
           "a caixa igual ao poligono");
    T_NEAR(area_recortada(quad, {0.25, 0.25}, {0.75, 0.75}), 0.25, 1e-15,
           "a caixa inteiramente DENTRO do poligono");
    T_NEAR(area_recortada(quad, {2.0, 2.0}, {3.0, 3.0}), 0.0, 1e-15,
           "a caixa inteiramente FORA");
    T_NEAR(area_recortada(quad, {0.5, -1.0}, {2.0, 2.0}), 0.5, 1e-15,
           "a caixa cortando pela metade");

    // UM POLIGONO NAO CONVEXO -- em L --, porque Sutherland-Hodgman produz
    // arestas degeneradas nele e a area tem de sair certa assim mesmo.
    const std::vector<Ponto2> ele = {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}};
    T_NEAR(area_recortada(ele, {-1, -1}, {3, 3}), 3.0, 1e-15,
           "a area do L inteiro (o poligono nao convexo)");
    T_NEAR(area_recortada(ele, {0.0, 0.0}, {1.0, 2.0}), 2.0, 1e-15,
           "a metade esquerda do L");
  }

  t_case("a_fracao_soma_a_area_do_poligono_em_precisao_de_maquina");
  {
    for (int nivel : {3, 5}) {
      Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
      Frente f = Frente::circulo({0.5, 0.5}, 0.3, 97);
      CampoEscalar c(m);
      escreve_fracao(f, m, &c);

      double soma = 0.0;
      long cortadas = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          const Ponto2 hh = m.lados_da_folha(b, i);
          const double v = c.em(m.id_global(b, i));
          soma += v * hh[0] * hh[1];
          if (v > 1e-12 && v < 1.0 - 1e-12) cortadas++;
        }
      T_NEAR(soma, f.area(), 1e-14,
             "a soma de fracao x volume contra a area do poligono");
      // O CASO TEM DE TER CELULA CORTADA.  Sem elas o recorte nunca roda, e a
      // soma bateria com qualquer implementacao de "0 ou 1".
      T_CHECK_MSG(cortadas > 0,
                  "nivel %d: nenhuma celula cortada -- o recorte nao foi "
                  "exercitado", nivel);
    }
  }

  t_case("a_forca_de_tensao_tem_resultante_nula_na_curva_fechada");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 5}}, sc_MPI_COMM_WORLD);
    DominioFacetas fd(m);
    Frente f = Frente::circulo({0.5, 0.5}, 0.25, 160);
    CampoVetorial forca(fd);
    forca_tensao(f, m, fd, 1.0, &forca);

    const Ponto2 tot = forca_total(fd, forca);
    T_BELOW(std::fabs(tot[0]), 1e-12, "a resultante em x da tensao superficial");
    T_BELOW(std::fabs(tot[1]), 1e-12, "a resultante em y da tensao superficial");

    // E A FORCA NAO E' NULA PONTO A PONTO -- sem isso, "resultante zero" nao
    // distingue a distribuicao correta de um campo todo nulo.
    double pico = 0.0;
    for (int d = 0; d < 2; d++)
      for (long q = 0; q < fd.num_facetas(d); q++)
        pico = std::max(pico, std::fabs(forca.em(d, q)));
    T_CHECK_MSG(pico > 1.0,
                "o pico da forca foi %.3e -- campo praticamente nulo, e a "
                "resultante nula nao afirma nada", pico);

    // A ASSIMETRIA e' o par: uma SEMICIRCUNFERENCIA nao e' fechada, e a
    // resultante dela NAO e' nula.  E' o contraste que mostra que o oraculo
    // mede a curva e nao devolve zero sempre.
    std::vector<Ponto2> meia;
    for (int k = 0; k < 80; k++) {
      const double a = PI * k / 79.0;
      meia.push_back({0.5 + 0.25 * std::cos(a), 0.5 + 0.25 * std::sin(a)});
    }
    meia.push_back({0.25, 0.5});           // fecha por uma corda reta
    Frente g(meia);
    CampoVetorial fg(fd);
    forca_tensao(g, m, fd, 1.0, &fg);
    const Ponto2 tg = forca_total(fd, fg);
    T_CHECK_MSG(std::hypot(tg[0], tg[1]) > 1e-3,
                "a meia-circunferencia deu resultante %.3e -- se tambem fosse "
                "nula, o caso do circulo nao afirmaria nada",
                std::hypot(tg[0], tg[1]));
  }

  t_case("o_salto_de_laplace_segue_sigma_sobre_r");
  {
    // O UNICO ORACULO QUE ENVOLVE O SOLVER, e por isso vem por ultimo.
    //
    // E O CONTRASTE E' DUPLO: o salto tem de seguir `sigma` LINEARMENTE e `R`
    // INVERSAMENTE.  Um unico par (sigma, R) nao distingue `sigma/R` de
    // qualquer constante que acerte aquele caso.
    // O RAIO TEM DE SATISFAZER `R > 4h`, e nao e' folga: a banda excluida da
    // medida tem `2,5h` de cada lado -- a largura do nucleo --, e sobra para o
    // interior um disco de raio `R - 2,5h`.  Com menos que `1,5h` ali nao ha'
    // centro de celula nenhum, e a medida fica sem o lado de dentro.
    //
    // A PRIMEIRA VERSAO USAVA `R = 0,20` no nivel 4, onde `4h = 0,25`: zero
    // celulas puras dentro.  A guarda de vacuidade pegou, que e' para o que
    // ela existe -- sem ela o caso teria comparado `0,0` com `5,0` e falhado
    // dizendo "erro relativo 1,0", mandando procurar defeito no metodo.
    struct Caso { double sigma, R; };
    const Caso casos[3] = {{1.0, 0.25}, {2.0, 0.25}, {1.0, 0.3125}};
    double dp[3];

    for (int q = 0; q < 3; q++) {
      const int nivel = 4;
      const double h = 1.0 / (1 << nivel);
      Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
      DominioFacetas fd(m);
      const Ponto2 centro{0.5, 0.5};
      const int nm = std::max(16, static_cast<int>(4.0 * 2.0 * PI * casos[q].R / h));
      Frente f = Frente::circulo(centro, casos[q].R, nm);

      const double E = 1e-9;
      auto zero = [](const Ponto2 &) { return 0.0; };
      ContornoPorRegiao b0, b1, bp;
      b0.padrao(Contorno::dirichlet(zero));
      b1.padrao(Contorno::dirichlet(zero));
      bp.registra("pino",
                  [E, h](const Ponto2 &x) { return std::fabs(x[1]) < E && x[0] < h; },
                  Contorno::dirichlet(zero));
      bp.padrao(Contorno::neumann(zero));

      NavierStokes ns(m, fd, b0, b1, bp);
      ns.reynolds(1.0);
      const double dt = h * h / 8.0;
      ns.passo_de_tempo(dt);

      CampoEscalar frac(m);
      CampoVetorial forca(fd);
      escreve_fracao(f, m, &frac);
      ns.forca_de_faceta(&forca);
      forca_tensao(f, m, fd, casos[q].sigma, &forca);

      for (long k = 0; k < static_cast<long>(0.02 / dt); k++) {
        ns.passo();
        T_CHECK_MSG(ns.pressao_convergiu(), "o Poisson nao convergiu no caso %d", q);
      }

      double pin = 0.0, pout = 0.0;
      long nin = 0, nout = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          const long id = m.id_global(b, i);
          const Ponto2 c = m.centroide(b, i);
          const double r = std::hypot(c[0] - centro[0], c[1] - centro[1]);
          if (std::fabs(r - casos[q].R) < 2.5 * h) continue;
          if (frac.em(id) > 0.99) { pin += ns.pressao().em(id); nin++; }
          else if (frac.em(id) < 0.01) { pout += ns.pressao().em(id); nout++; }
        }
      T_CHECK_MSG(nin > 0 && nout > 0,
                  "caso %d (sigma=%.1f R=%.4f, h=%.4f): %ld celulas puras "
                  "dentro, %ld fora -- o oraculo nao mediu nada.  O interior "
                  "fora da banda tem raio %.4f, e precisa de ao menos 1,5h = "
                  "%.4f: a condicao e' R > 4h",
                  q, casos[q].sigma, casos[q].R, h, nin, nout,
                  casos[q].R - 2.5 * h, 1.5 * h);
      dp[q] = (nin > 0 && nout > 0) ? (pin / nin - pout / nout) : 0.0;
    }

    for (int q = 0; q < 3; q++) {
      const double exato = casos[q].sigma / casos[q].R;
      const double rel = std::fabs(dp[q] - exato) / exato;
      T_CHECK_MSG(rel < 5e-3,
                  "caso sigma=%.1f R=%.2f: salto %.6f contra %.6f exato, erro "
                  "relativo %.3e", casos[q].sigma, casos[q].R, dp[q], exato, rel);
    }

    // O CONTRASTE: dobrar sigma dobra o salto; reduzir R aumenta-o na razao
    // inversa.  Sao as duas dependencias que `sigma/R` afirma.
    T_NEAR(dp[1] / dp[0], 2.0, 0.02, "a razao do salto ao dobrar sigma");
    T_NEAR(dp[2] / dp[0], 0.25 / 0.3125, 0.02, "a razao do salto ao aumentar R");
  }

  return t_end();
}
