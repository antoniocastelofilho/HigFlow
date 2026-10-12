// O SOLVER DE NAVIER-STOKES: os operadores que ja' tratam a transicao 2:1.
//
// O QUE ESTE ARQUIVO VERIFICA HOJE e' o operador de DIVERGENCIA sobre malha
// adaptada -- o primeiro pedaco da restricao 2:1.  O estencil da quantidade de
// movimento ainda nao a trata, e o `passo` recusa; quando tratar, os casos dele
// entram aqui.
//
// O ORACULO E' DE EXATIDAO, nao de ordem: a divergencia de volumes finitos e' a
// regra do ponto medio sobre a superficie, que e' EXATA para campo linear.
// Logo um campo linear com divergencia nula tem de dar ZERO DE MAQUINA em
// qualquer malha -- uniforme ou adaptada.  Ordem de convergencia certa convive
// com constante errada; exatidao nao.
//
// E CADA AFIRMACAO DE ZERO VEM COM O PAR QUE NAO E' ZERO.  Sem ele, "a
// divergencia deu 1e-16" nao distingue um operador correto de um que devolve
// zero sempre -- e um operador que descarta faces devolve justamente isso em
// metade dos casos.

#include "../testing.h"

#include <cmath>
#include <algorithm>
#include <functional>

#include "arranque.hxx"
#include "campos2d.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "navierstokes2d.hxx"

using namespace t8flow2d;

namespace {

//! Um contorno minimo: velocidade nula e pressao com UM ponto de Dirichlet.
//! O solver recusa pressao toda-Neumann, e com razao -- aqui so' se quer
//! construir o objeto para medir o operador.
void arma(ContornoPorRegiao *b0, ContornoPorRegiao *b1, ContornoPorRegiao *bp,
          double h)
{
  const double E = 1e-9;
  auto zero = [](const Ponto2 &) { return 0.0; };
  b0->padrao(Contorno::dirichlet(zero));
  b1->padrao(Contorno::dirichlet(zero));
  bp->registra("pino",
               [E, h](const Ponto2 &x) { return std::fabs(x[1]) < E && x[0] < h; },
               Contorno::dirichlet(zero));
  bp->padrao(Contorno::neumann(zero));
}

//! Escreve um campo analitico nas facetas.
void poe(const DominioFacetas &fd, CampoVetorial *u,
         const std::function<double(const Ponto2 &, int)> &f)
{
  for (int d = 0; d < 2; d++)
    for (long k = 0; k < fd.num_facetas(d); k++)
      u->em(d, k) = f(fd.centro(d, k), d);
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Os tres campos lineares do arquivo, e o que cada um afirma.
  auto solenoidal = [](const Ponto2 &x, int d) {
    return (d == 0) ? x[0] : -x[1];          // div = 1 - 1 = 0
  };
  auto cisalhante = [](const Ponto2 &x, int d) {
    return (d == 0) ? x[1] : 0.0;            // div = 0, mas VARIA ao longo da face
  };
  auto com_fonte = [](const Ponto2 &x, int d) {
    return (d == 0) ? x[0] : 0.0;            // div = 1, exatamente
  };

  t_case("a_divergencia_e_exata_em_malha_uniforme");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_WORLD);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, bp;
    arma(&b0, &b1, &bp, 1.0 / 8.0);
    NavierStokes ns(m, fd, b0, b1, bp);

    poe(fd, &ns.velocidade(), solenoidal);
    T_BELOW(ns.divergencia_maxima(), 1e-13,
            "a divergencia de um campo linear solenoidal em malha uniforme");

    // O PAR QUE NAO E' ZERO.  Sem ele o caso acima passa num operador que
    // devolve zero sempre.
    poe(fd, &ns.velocidade(), com_fonte);
    T_NEAR(ns.divergencia_maxima(), 1.0, 1e-13,
           "a divergencia de `u = x, v = 0`, que vale 1 exatamente");
  }

  t_case("a_divergencia_e_exata_em_malha_adaptada");
  {
    // O CASO QUE MOTIVA TUDO.  Numa transicao 2:1 a face da celula GROSSA nao
    // tem faceta propria: ela e' coberta por DUAS facetas das vizinhas finas, e
    // o centro da face grossa nao e' o centro de nenhuma das duas.
    //
    // A versao anterior buscava a faceta por esse centro, recebia -1 e
    // DESCARTAVA o fluxo da face inteira -- em silencio.  O erro resultante nao
    // e' de truncamento, e' de O(1).
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_WORLD);
    m.refina_caixa({0.0, 0.0}, {0.5, 0.5}, 4);
    T_CHECK_MSG(m.num_folhas_em_transicao() > 0,
                "a malha deste caso nao tem transicao 2:1 -- ele nao verifica "
                "nada do que diz verificar");

    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, bp;
    arma(&b0, &b1, &bp, 1.0 / 16.0);
    NavierStokes ns(m, fd, b0, b1, bp);

    poe(fd, &ns.velocidade(), solenoidal);
    T_BELOW(ns.divergencia_maxima(), 1e-13,
            "a divergencia de um campo linear solenoidal em malha ADAPTADA");

    // CISALHANTE: `u = y` varia AO LONGO da face, logo as duas sub-facetas de
    // uma face grossa tem valores DIFERENTES.  E' o campo que distingue somar
    // as duas de tomar uma delas.
    poe(fd, &ns.velocidade(), cisalhante);
    T_BELOW(ns.divergencia_maxima(), 1e-13,
            "a divergencia de `u = y, v = 0` em malha adaptada -- as duas "
            "sub-facetas de uma face grossa tem valores diferentes");

    poe(fd, &ns.velocidade(), com_fonte);
    T_NEAR(ns.divergencia_maxima(), 1.0, 1e-13,
           "a divergencia de `u = x, v = 0` em malha adaptada, que vale 1");
  }

  t_case("o_passo_fecha_a_divergencia_em_malha_adaptada");
  {
    // O ORACULO DA RESTRICAO 2:1, e o unico que nao precisa de solucao de
    // referencia: o metodo de projecao EXISTE para zerar a divergencia.  Se o
    // gradiente que a projecao aplica nao for o mesmo que a divergencia do
    // Poisson supos, o campo corrigido nao e' solenoidal -- e NADA diz, porque
    // o solver converge igual.  E' o par `div`/`grad` que fecha, nao cada um.
    //
    // Na transicao 2:1 ha' tres jeitos de o par deixar de fechar: a face da
    // grossa coberta por duas facetas finas (somatorio), a distancia entre
    // centros que nao e' `h` (coeficiente), e os centros desalinhados da
    // normal (correcao transversal).  Qualquer um deles quebra isto.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_WORLD);
    m.refina_caixa({0.25, 0.25}, {0.75, 0.75}, 4);
    T_CHECK_MSG(m.num_folhas_em_transicao() > 0,
                "a malha deste caso nao tem transicao 2:1");

    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, bp;
    arma(&b0, &b1, &bp, 1.0 / 8.0);
    NavierStokes ns(m, fd, b0, b1, bp);
    ns.reynolds(1.0);
    ns.passo_de_tempo(1.0 / 8.0 / 8.0 / 16.0 / 16.0);

    // Uma forca de corpo qualquer, para que haja o que projetar: sem fonte o
    // campo fica no repouso e "divergencia nula" nao afirma nada.
    ns.termo_fonte_faceta([](const Ponto2 &x, int d, double) {
      return (d == 0) ? std::sin(3.0 * x[1]) : std::cos(3.0 * x[0]);
    });

    for (int k = 0; k < 5; k++) {
      ns.passo();
      T_CHECK_MSG(ns.pressao_convergiu(),
                  "o Poisson da projecao nao convergiu no passo %d", k + 1);
    }
    T_BELOW(ns.divergencia_maxima(), 1e-12,
            "a divergencia depois de cinco passos em malha ADAPTADA");

    // A COMPATIBILIDADE: numa caixa fechada o que entra e' o que sai.
    T_BELOW(std::fabs(ns.soma_do_lado_direito()), 1e-10,
            "a soma dos fluxos do Poisson numa caixa fechada");

    // CONTRASTE: o campo NAO e' trivialmente nulo -- sem isto, "divergencia
    // 1e-14" nao distingue a projecao funcionando de nada ter acontecido.
    double maxu = 0.0;
    for (int d = 0; d < 2; d++)
      for (long f = 0; f < fd.num_facetas(d); f++)
        maxu = std::max(maxu, std::fabs(ns.velocidade().em(d, f)));
    T_CHECK_MSG(maxu > 1e-6,
                "a velocidade ficou em %.3e -- nada se moveu, e a divergencia "
                "nula nao afirma nada", maxu);
  }

  t_case("a_pressao_toda_neumann_e_recusada");
  {
    // A CAIXA FECHADA poe a constante no nucleo do operador, e o PETSc converge
    // para um elemento QUALQUER do subespaco -- isso da' numero.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 2}}, sc_MPI_COMM_WORLD);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, bp;
    auto zero = [](const Ponto2 &) { return 0.0; };
    b0.padrao(Contorno::dirichlet(zero));
    b1.padrao(Contorno::dirichlet(zero));
    bp.padrao(Contorno::neumann(zero));

    bool lancou = false;
    std::string msg;
    try { NavierStokes ns(m, fd, b0, b1, bp); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "a pressao toda-Neumann foi aceita");
    T_CHECK_MSG(msg.find("Dirichlet") != std::string::npos,
                "a mensagem nao diz a saida: \"%s\"", msg.c_str());

    // CONTRASTE: com UM ponto de Dirichlet, constroi.
    ContornoPorRegiao bp2;
    arma(&b0, &b1, &bp2, 1.0 / 4.0);
    bool ok = true;
    try { NavierStokes ns(m, fd, b0, b1, bp2); }
    catch (const std::exception &) { ok = false; }
    T_CHECK_MSG(ok, "um ponto de Dirichlet nao bastou para construir");
  }

  return t_end();
}
