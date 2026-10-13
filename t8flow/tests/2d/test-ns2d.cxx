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

  // =========================================================================
  // O OPERADOR DE DIFUSAO PELOS DOIS CAMINHOS
  // =========================================================================
  //
  // O semi-implicito precisa do operador em forma de MATRIZ; o explicito o
  // avalia como NUMERO.  Sao duas implementacoes da mesma coisa, e o historico
  // deste codigo diz o que isso custa: o gradiente transversal tinha duas, elas
  // concordavam no interior e discordavam NA PAREDE, e o par div/grad deixou de
  // fechar sem nada acusar.
  //
  // Aqui as duas ficam de proposito, porque colapsa-las trocaria a ordem das
  // somas e os casos monofasicos deixariam de reproduzir bit a bit.  O que as
  // mantem honestas e' este caso.
  t_case("o_operador_de_difusao_concorda_pelos_dois_caminhos");
  {
    // `mu` NAO UNIFORME, e isso e' o ponto: com `mu` constante os pesos ficam
    // simetricos e um erro de media na quina nao aparece.
    auto mu_campo = [](const Ponto2 &x) {
      return 1.0 + 3.0 * x[0] + 7.0 * x[1] * x[1];
    };
    // Campo de velocidade SEM simetria: um campo simetrico cancela erros de
    // sinal entre os dois lados.
    auto vel = [](const Ponto2 &x, int d) {
      return (d == 0) ? 0.3 + 1.7 * x[0] - 0.9 * x[1] + 2.1 * x[0] * x[1]
                      : -0.4 + 0.6 * x[0] + 1.3 * x[1] - 1.1 * x[0] * x[0];
    };

    struct Caso { const char *nome; bool adaptada; };
    const Caso casos[] = {{"uniforme", false}, {"adaptada 2:1", true}};

    for (const Caso &caso : casos) {
      const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, 4}};
      Malha m(blocos, sc_MPI_COMM_WORLD);
      if (caso.adaptada) {
        // Refino que ENCOSTA NA PAREDE: foi exatamente ali que as duas
        // implementacoes do gradiente transversal discordaram.  Um bloco
        // refinado so' no interior nao discriminaria.
        m.refina_caixa({0.0, 0.0}, {0.4, 0.45}, 5);
      }
      DominioFacetas fd(m);

      ContornoPorRegiao b0, b1, bp;
      // OS TRES TIPOS DE CONTORNO NO MESMO CASO: Dirichlet numa componente,
      // Neumann na outra.  O fechamento afim tem de valer para os dois, e e'
      // o Neumann que um tratamento preguicoso (fantasma inteiro no lado
      // direito) erra sem falhar.
      auto zero = [](const Ponto2 &) { return 0.0; };
      auto rampa = [](const Ponto2 &x) { return 0.25 * x[0] - 0.1; };
      b0.padrao(Contorno::dirichlet(rampa));
      b1.padrao(Contorno::neumann(rampa));
      // O PINO TEM DE PEGAR EXATAMENTE UMA FACETA DE PAREDE, e quem decide
      // isso e' o `h` DA PAREDE -- nao o da base e nao o minimo.  Aqui a
      // parede de baixo e' nivel 4 no caso uniforme (primeiro centro em
      // `1/32`) e nivel 5 no adaptado (centros em `1/64` e `3/64`).  O limiar
      // `1/24` fica acima do primeiro centro nos DOIS e abaixo do segundo no
      // fino.
      //
      // Com `x < 1/32` -- que parecia o natural -- o caso UNIFORME nao pegava
      // nenhuma: o centro esta' EXATAMENTE em `1/32` e o teste e' estrito.  O
      // solver entao recusou a pressao toda-Neumann, com razao, e o binario
      // abortou.
      bp.registra("pino",
                  [](const Ponto2 &x) {
                    return std::fabs(x[1]) < 1e-9 && x[0] < 1.0 / 24.0;
                  },
                  Contorno::dirichlet(zero));
      bp.padrao(Contorno::neumann(zero));

      NavierStokes ns(m, fd, b0, b1, bp);
      CampoEscalar mu(m);
      for (long bb = 0; bb < m.num_blocos(); bb++)
        for (long ii = 0; ii < m.num_folhas_do_bloco(bb); ii++)
          mu.em(m.id_global(bb, ii)) = mu_campo(m.centroide(bb, ii));
      ns.viscosidade(&mu);
      poe(fd, &ns.velocidade(), vel);

      double pior = 0.0, escala = 0.0;
      long onde_d = -1, onde_f = -1;
      for (int d = 0; d < 2; d++)
        for (long f = 0; f < fd.num_facetas(d); f++) {
          const double a = ns.difusao_explicita(d, f);
          const double bb = ns.difusao_por_pesos(d, f);
          escala = std::max(escala, std::fabs(a));
          if (std::fabs(a - bb) > pior) {
            pior = std::fabs(a - bb);
            onde_d = d;
            onde_f = f;
          }
        }
      T_CHECK_MSG(escala > 1.0,
                  "%s: o operador deu quase zero em todo lado (%.3e) -- o caso "
                  "nao exercita nada", caso.nome, escala);
      T_CHECK_MSG(pior <= 1e-11 * escala,
                  "%s: os dois caminhos discordam em %.3e (escala %.3e), pior "
                  "na faceta %ld da direcao %ld",
                  caso.nome, pior, escala, onde_f, onde_d);
    }
  }

  // =========================================================================
  t_case("os_dois_avancos_temporais_convergem_para_o_mesmo_campo");
  {
    // O CONTRASTE QUE DA' SENTIDO AO ZERO: com `dt` grande os dois DIFEREM --
    // sao esquemas diferentes --, e a diferenca tem de CAIR com `dt`.  Sem
    // medir os dois `dt`, "a diferenca e' pequena" nao distingue um
    // semi-implicito correto de um que ignorou o termo viscoso.
    auto corre = [](double dt, Tempo tempo, std::vector<double> *fora) {
      const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, 4}};
      Malha m(blocos, sc_MPI_COMM_WORLD);
      DominioFacetas fd(m);
      ContornoPorRegiao b0, b1, bp;
      arma(&b0, &b1, &bp, 1.0 / 16.0);
      NavierStokes ns(m, fd, b0, b1, bp);
      ns.passo_de_tempo(dt);
      ns.tempo(tempo);
      // Um campo inicial com cisalhamento: e' o que o termo viscoso ataca.
      poe(fd, &ns.velocidade(), [](const Ponto2 &x, int d) {
        return (d == 0) ? std::sin(M_PI * x[1]) * std::sin(M_PI * x[0]) : 0.0;
      });
      const int n = static_cast<int>(std::lround(0.02 / dt));
      for (int k = 0; k < n; k++) ns.passo();
      fora->clear();
      for (int d = 0; d < 2; d++)
        for (long f = 0; f < fd.num_facetas(d); f++)
          fora->push_back(ns.velocidade().em(d, f));
    };

    auto difere = [](const std::vector<double> &a,
                     const std::vector<double> &b) {
      double pior = 0.0;
      for (size_t k = 0; k < a.size(); k++)
        pior = std::max(pior, std::fabs(a[k] - b[k]));
      return pior;
    };

    std::vector<double> ex1, im1, ex2, im2;
    corre(1.0e-3, Tempo::EulerExplicito, &ex1);
    corre(1.0e-3, Tempo::EulerSemiImplicito, &im1);
    corre(2.5e-4, Tempo::EulerExplicito, &ex2);
    corre(2.5e-4, Tempo::EulerSemiImplicito, &im2);

    const double d1 = difere(ex1, im1);
    const double d2 = difere(ex2, im2);
    T_CHECK_MSG(d1 > 1e-9,
                "com dt = 1e-3 os dois avancos deram o MESMO campo (%.3e) -- "
                "o semi-implicito provavelmente nao esta' fazendo nada", d1);
    T_CHECK_MSG(d2 < d1 / 2.0,
                "a diferenca entre os avancos nao caiu com dt: %.3e com 1e-3 "
                "contra %.3e com 2,5e-4", d1, d2);
  }

  // =========================================================================
  t_case("o_semi_implicito_e_estavel_acima_do_limite_viscoso_explicito");
  {
    // O GANHO DECLARADO, MEDIDO.  O limite viscoso explicito e'
    // `dt < rho h^2 / (4 mu)`; com `rho = 1`, `mu = 1` e `h = 1/16` isso da'
    // 9,77e-4.  Em `4x` esse valor o explicito TEM de estourar e o
    // semi-implicito TEM de sobreviver -- se os dois sobrevivessem, o limite
    // nao seria o que eu digo que e' e o ganho seria afirmacao vazia.
    auto corre = [](double dt, Tempo tempo) {
      const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, 4}};
      Malha m(blocos, sc_MPI_COMM_WORLD);
      DominioFacetas fd(m);
      ContornoPorRegiao b0, b1, bp;
      arma(&b0, &b1, &bp, 1.0 / 16.0);
      NavierStokes ns(m, fd, b0, b1, bp);
      ns.passo_de_tempo(dt);
      ns.tempo(tempo);
      poe(fd, &ns.velocidade(), [](const Ponto2 &x, int d) {
        return (d == 0) ? std::sin(M_PI * x[1]) * std::sin(M_PI * x[0]) : 0.0;
      });
      double pior = 0.0;
      for (int k = 0; k < 60; k++) {
        ns.passo();
        for (int d = 0; d < 2; d++)
          for (long f = 0; f < fd.num_facetas(d); f++)
            pior = std::max(pior, std::fabs(ns.velocidade().em(d, f)));
        if (!(pior < 1e4)) break;
      }
      return pior;
    };

    const double h = 1.0 / 16.0;
    const double dt_lim = h * h / 4.0;            // rho = mu = 1
    const double dt = 4.0 * dt_lim;

    const double ex = corre(dt, Tempo::EulerExplicito);
    const double im = corre(dt, Tempo::EulerSemiImplicito);
    T_CHECK_MSG(!(ex < 10.0),
                "com dt = %.3e (4x o limite viscoso %.3e) o EXPLICITO ficou "
                "limitado em %.3e -- o limite nao e' o que se diz",
                dt, dt_lim, ex);
    T_CHECK_MSG(im < 2.0,
                "com dt = %.3e o semi-implicito estourou (%.3e)", dt, im);
  }

  // =========================================================================
  // A ADJUNCAO `D = -G^T`
  // =========================================================================
  //
  // AFIRMA SOBRE O PAR, e nao sobre cada operador.  Divergencia e gradiente
  // tem casos proprios acima -- exatidao em campo linear, fechamento da
  // projecao --, e os dois podem passar com o PAR desencontrado: a adjuncao
  // diz se a troca discreta de energia entre velocidade e pressao se conserva,
  //
  //     <Du, p>_V = -<u, Gp>_W,
  //
  // e nenhum teste de um operador so' a ve'.
  //
  // O PESO NAO E' ESCOLHA.  Igualando os coeficientes de `u_f p_i` nas duas
  // formas bilineares, a igualdade exige `W_f = medida_f x dist_f` -- o volume
  // de controle dual.  Com outro peso ela nao sai nem no caso uniforme e a
  // pergunta perde sentido.
  //
  // A VELOCIDADE E' NULA NAS FACETAS DE CONTORNO, e isso tambem nao e'
  // conveniencia: ali ela NAO E' INCOGNITA -- vem da condicao --, e a forma
  // bilinear so' envolve incognitas.  `grad_campo` devolve zero na faceta de
  // contorno e o volume dual tambem, entao os dois lados as descartam de modo
  // consistente.
  //
  // MEDIDO AO DESENHAR UM SOLVER ACOPLADO, QUE FICOU DE LADO.  O numero
  // sobreviveu ao plano: no metodo segregado a nao adjuncao e' invisivel -- o
  // Poisson la' e' montado literalmente como `D G`, consistente por construcao
  // --, mas e' ela que caracteriza o par, e e' dela que sairia a escolha do
  // Krylov se houver acoplado algum dia.
  t_case("a_divergencia_e_o_gradiente_sao_adjuntos");
  {
    auto mede_gap = [](bool adaptada, double *rel, int niv = 4) {
      const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, niv}};
      Malha m(blocos, sc_MPI_COMM_WORLD);
      if (adaptada) m.refina_caixa({0.25, 0.25}, {0.6, 0.65}, niv + 1);
      DominioFacetas fd(m);
      ContornoPorRegiao b0, b1, bp;
      arma(&b0, &b1, &bp, 1.0 / (1 << niv));
      NavierStokes ns(m, fd, b0, b1, bp);

      // CAMPOS SEM SIMETRIA NEM DIVERGENCIA NULA: um campo solenoidal daria
      // `Du = 0` e os dois lados seriam zero -- concordancia vazia.
      auto fu = [](const Ponto2 &x, int d) {
        return (d == 0) ? 0.7 + 1.3 * x[0] - 0.4 * x[1] + 2.0 * x[0] * x[1]
                        : -0.2 + 0.9 * x[0] * x[0] + 1.1 * x[1];
      };
      poe(fd, &ns.velocidade(), fu);

      // ZERA O CONTORNO: ali a velocidade nao e' incognita.
      for (int d = 0; d < 2; d++)
        for (long f = 0; f < fd.num_facetas(d); f++)
          if (ns.volume_dual_da_faceta(d, f) == 0.0)
            ns.velocidade().em(d, f) = 0.0;

      CampoEscalar q(m);
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          const Ponto2 c = m.centroide(b, i);
          q.em(m.id_global(b, i)) =
              0.3 + 1.7 * c[0] - 0.8 * c[1] + 2.3 * c[0] * c[1]
              + 0.6 * c[1] * c[1];
        }

      double lado_d = 0.0, lado_g = 0.0, escala = 0.0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
          const Ponto2 hh = m.lados_da_folha(b, i);
          const double v = hh[0] * hh[1];
          const double termo =
              ns.divergencia_em(b, i) * q.em(m.id_global(b, i)) * v;
          lado_d += termo;
          escala += std::fabs(termo);
        }
      for (int d = 0; d < 2; d++)
        for (long f = 0; f < fd.num_facetas(d); f++) {
          const double w = ns.volume_dual_da_faceta(d, f);
          if (w == 0.0) continue;
          const double termo =
              ns.velocidade().em(d, f) * ns.gradiente_em(q, d, f) * w;
          lado_g += termo;
          escala += std::fabs(termo);
        }
      *rel = std::fabs(lado_d + lado_g) / escala;
      return escala;
    };

    double rel_u = 0.0, rel_a = 0.0;
    const double esc_u = mede_gap(false, &rel_u);
    const double esc_a = mede_gap(true, &rel_a);

    T_CHECK_MSG(esc_u > 1e-3 && esc_a > 1e-3,
                "as formas bilineares deram quase zero (%.3e, %.3e) -- os "
                "campos nao exercitam nada", esc_u, esc_a);
    T_CHECK_MSG(rel_u <= 1e-13,
                "MALHA UNIFORME: D e G nao sao adjuntos, gap relativo %.3e",
                rel_u);
    // NA MALHA ADAPTADA A ADJUNCAO NAO VALE EXATAMENTE, e isso era previsivel
    // antes de medir: o gradiente leva uma CORRECAO TRANSVERSAL na transicao
    // 2:1 -- o centro da celula grossa esta' deslocado da linha da faceta -- e
    // a divergencia nao tem termo que a case.  Sem a correcao o Poisson tinha
    // ordem NEGATIVA na malha adaptada, entao ela nao e' negociavel.
    //
    // O QUE IMPORTA NAO E' SE VALE, E' EM QUE ORDEM FALHA.  Medido:
    //
    //     nivel 4  8,47e-06      nivel 5  1,78e-06      nivel 6  4,18e-07
    //     razoes      4,75                  4,26
    //
    // SEGUNDA ORDEM.  A nao adjuncao e' perturbacao `O(h^2)` que some no
    // limite, e nao defeito estrutural.  Consequencias, as tres:
    //
    //   * a troca discreta de energia entre velocidade e pressao se conserva
    //     ate' `O(h^2)`, e nao exatamente -- e e' isso que o caso afirma;
    //   * num sistema ACOPLADO, se vier a haver um, a matriz seria nao
    //     simetrica por essa mesma ordem: Krylov de GMRES e nao de MINRES, com
    //     o complemento de Schur montado como se simetrico, valido ate'
    //     `O(h^2)`.
    //
    // A ASSERCAO E' A ORDEM, e nao o valor: um valor pequeno num nivel so' nao
    // distingue "perturbacao que some" de "defeito pequeno que fica".
    double r4 = 0.0, r5 = 0.0, r6 = 0.0;
    mede_gap(true, &r4, 4);
    mede_gap(true, &r5, 5);
    mede_gap(true, &r6, 6);
    T_CHECK_MSG(r4 > 1e4 * rel_u,
                "o gap da adaptada (%.3e) nao se distingue do da uniforme "
                "(%.3e) -- a medida de ordem abaixo seria ruido", r4, rel_u);
    T_CHECK_MSG(r4 / r5 > 3.5 && r5 / r6 > 3.5,
                "a nao adjuncao NAO cai em segunda ordem: %.3e, %.3e, %.3e -- "
                "razoes %.2f e %.2f.  Se nao cai, nao e' perturbacao: e' "
                "defeito estrutural do par div/grad na transicao 2:1",
                r4, r5, r6, r4 / r5, r5 / r6);
  }

  return t_end();
}
