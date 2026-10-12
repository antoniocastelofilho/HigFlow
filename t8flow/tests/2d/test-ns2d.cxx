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

  t_case("o_passo_recusa_a_malha_adaptada_dizendo_o_que_falta");
  {
    // O FLUXO JA' TRATA 2:1; o estencil da quantidade de movimento ainda nao.
    // A recusa esta' no PASSO e nao na construcao, de proposito: medir o fluxo
    // numa malha adaptada e' parte de chegar la'.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_WORLD);
    m.refina_caixa({0.0, 0.0}, {0.5, 0.5}, 4);
    DominioFacetas fd(m);
    ContornoPorRegiao b0, b1, bp;
    arma(&b0, &b1, &bp, 1.0 / 16.0);

    bool construiu = true;
    try { NavierStokes ns(m, fd, b0, b1, bp); }
    catch (const std::exception &) { construiu = false; }
    T_CHECK_MSG(construiu,
                "a construcao recusou malha adaptada -- os operadores de fluxo "
                "precisam ser mediveis nela");

    NavierStokes ns(m, fd, b0, b1, bp);
    ns.passo_de_tempo(1e-4);
    bool lancou = false;
    std::string msg;
    try { ns.passo(); }
    catch (const std::runtime_error &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "o passo aceitou malha com transicao 2:1");
    T_CHECK_MSG(msg.find("quantidade de movimento") != std::string::npos,
                "a mensagem nao diz O QUE falta: \"%s\"", msg.c_str());
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
