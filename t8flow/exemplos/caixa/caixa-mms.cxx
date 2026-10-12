// A CAIXA COM TAMPA DESLIZANTE, verificada por SOLUCAO MANUFATURADA.
//
// POR QUE ESTE CASO DEPOIS DO CANAL: o canal nao exercita o termo convectivo.
// Nele `u` nao depende de `x` e `v` e' nulo, de modo que `(u.grad)u` se anula
// IDENTICAMENTE -- central, upwind de primeira e CUBISTA dao o mesmo
// resultado.  Aqui as duas componentes variam nas duas direcoes e o convectivo
// e' o maior termo da equacao em boa parte do dominio.
//
// E POR QUE MANUFATURADA e nao o problema classico de Ghia: a cavidade com
// tampa de velocidade constante tem SINGULARIDADE nas duas quinas de cima --
// a velocidade salta de `U` para zero num ponto --, e nenhum esquema converge
// na ordem nominal perto dela.  Medir ordem ali mede a singularidade.  A
// solucao manufaturada e' suave ate' a parede e a ordem medida e' a do
// esquema.
//
// SO' A TAMPA DEPENDE DO TEMPO: `u(x,1,t) = -f(t) x^2 (1-x)^2`, e os outros
// tres lados sao nulos em todo instante.  Ver `manufaturada_caixa.hxx`.
//
// O PASSO DE TEMPO ACOMPANHA `h^2`, e nao e' conveniencia: Euler explicito e'
// de PRIMEIRA ordem no tempo, logo o erro total e' `C1 dt + C2 h^2`.  Com `dt`
// fixo o termo temporal domina e a ordem medida no espaco seria 0; com
// `dt ~ h^2` os dois caem juntos e o que se le' e' a ordem do esquema.  O
// limite difusivo do Euler explicito e' `Re h^2 / 4`, de modo que `h^2 / 8`
// tambem e' estavel -- as duas exigencias apontam para o mesmo lugar.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "arranque.hxx"
#include "campos2d.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "manufaturada_caixa.hxx"
#include "navierstokes2d.hxx"
#include "saida_vtk2d.hxx"

using namespace t8flow2d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0, np = 1;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &np);

  double re = 1.0, T = 0.0625, cfl_dif = 8.0;
  int nivel = 4;
  long cada_vtk = 0;
  double refina = 0.0;
  std::string disc_conv = "cubista", disc_proj = "incremental";

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--re=", 5)) re = std::atof(a + 5);
    else if (!std::strncmp(a, "--T=", 4)) T = std::atof(a + 4);
    else if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--divisor-dt=", 13)) cfl_dif = std::atof(a + 13);
    else if (!std::strncmp(a, "--convectivo=", 13)) disc_conv = a + 13;
    else if (!std::strncmp(a, "--projecao=", 11)) disc_proj = a + 11;
    else if (!std::strncmp(a, "--vtk=", 6)) cada_vtk = std::atol(a + 6);
    else if (!std::strncmp(a, "--refina=", 9)) refina = std::atof(a + 9);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf(
            "uso: ./exemplos/caixa/caixa-mms [opcoes]\n\n"
            "  --re           numero de Reynolds\n"
            "  --T            instante final\n"
            "  --nivel        refino (quadrado unitario)\n"
            "  --divisor-dt   dt = h^2 / divisor (padrao 8; limite difusivo 4)\n"
            "  --convectivo   central, upwind1, cubista, quick\n"
            "  --projecao     incremental, nao_incremental\n"
            "  --vtk=N        grava VTK a cada N passos\n");
      return 0;
    }
  }

  if (np > 1) {
    if (rank == 0)
      std::printf("ESTE CASO AINDA E' SERIAL (%d processos pedidos).\n", np);
    return 2;
  }

  Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  // REFINO OPCIONAL numa caixa interior, para exercitar a transicao 2:1 com a
  // mesma solucao manufaturada -- o erro passa a ser comparavel entre a malha
  // uniforme e a adaptada, que e' o que torna a comparacao uma medida.
  if (refina > 0.0)
    m.refina_caixa({0.25, 0.25}, {0.25 + refina, 0.25 + refina}, nivel + 1);
  DominioFacetas fd(m);
  // O PASSO SAI DO `h` MINIMO, e nao do da base.
  //
  // O limite difusivo do Euler explicito e' `Re h^2 / 4`, com o `h` da MENOR
  // celula -- refinar um pedaco aperta o passo por 4 em todo o dominio, nao
  // so' ali.  Usar o `h` da base deu campo divergindo em silencio na primeira
  // tentativa: o solver nao reclama, o resultado sai NaN depois de algumas
  // centenas de passos.
  const double h_base = 1.0 / (1 << nivel);
  const double hh = (refina > 0.0) ? h_base / 2.0 : h_base;
  const double dt = hh * hh / cfl_dif;

  // O RELOGIO DO CONTORNO vive AQUI, e quem o move e' o solver.  Ver
  // `relogio_do_contorno` no cabecalho do modulo: a condicao precisa do
  // relogio e o solver precisa da condicao, logo o relogio tem de nascer
  // antes dos dois.
  double tbc = 0.0;

  const double EPSG = 1e-9;
  auto na_tampa = [EPSG](const Ponto2 &x) { return std::fabs(x[1] - 1.0) < EPSG; };
  auto no_fundo = [EPSG](const Ponto2 &x) { return std::fabs(x[1]) < EPSG; };
  auto nas_laterais = [EPSG](const Ponto2 &x) {
    return std::fabs(x[0]) < EPSG || std::fabs(x[0] - 1.0) < EPSG;
  };

  ContornoPorRegiao bcu0, bcu1, bcp;

  // SO' A TAMPA DEPENDE DO TEMPO.  As outras tres paredes sao Dirichlet NULO
  // em todo instante, e isso e' propriedade da solucao manufaturada (as
  // raizes de `G` e `H`), nao uma simplificacao do caso.
  bcu0.registra("tampa", na_tampa, Contorno::dirichlet([&tbc](const Ponto2 &x) {
    return caixa::u_exata(x, 0, tbc);
  }));
  bcu0.registra("fundo", no_fundo, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));
  bcu0.registra("laterais", nas_laterais, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));

  bcu1.registra("tampa", na_tampa, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));
  bcu1.registra("fundo", no_fundo, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));
  bcu1.registra("laterais", nas_laterais, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));

  // A PRESSAO: Neumann em toda a volta, MENOS num ponto.
  //
  // Com velocidade de Dirichlet nas quatro paredes, a condicao natural da
  // pressao e' Neumann em todas, e ai' a constante fica no nucleo do operador
  // -- a matriz e' singular e o solver devolve um elemento qualquer do
  // subespaco.  O modulo RECUSA esse caso em vez de resolver.
  //
  // A saida e' Dirichlet em UM ponto de parede, e numa solucao manufaturada o
  // valor ali e' conhecido: `p_exata`.  O ponto e' a primeira faceta do fundo,
  // e a regiao e' definida por `x < h` para que seja EXATAMENTE uma faceta em
  // qualquer nivel de refino -- uma faixa de largura fixa cobriria mais
  // facetas a cada refino, e o caso deixaria de ser o mesmo.
  // O PINO USA O `h` DA PAREDE, e nao o `h` minimo do dominio.
  //
  // Com refino no interior o minimo cai pela metade e `x < h_min` nao casa com
  // faceta nenhuma do fundo -- a pressao fica toda-Neumann e o solver recusa.
  // A regiao tem de ser a largura de UMA faceta DAQUELA parede.
  auto no_pino = [EPSG, h_base](const Ponto2 &x) {
    return std::fabs(x[1]) < EPSG && x[0] < h_base;
  };
  auto no_fundo_sem_pino = [EPSG, h_base](const Ponto2 &x) {
    return std::fabs(x[1]) < EPSG && x[0] >= h_base;
  };
  bcp.registra("pino", no_pino, Contorno::dirichlet([](const Ponto2 &x) {
    return caixa::p_exata(x, 0.0);
  }));
  bcp.registra("fundo", no_fundo_sem_pino, Contorno::neumann([](const Ponto2 &) { return 0.0; }));
  bcp.registra("tampa", na_tampa, Contorno::neumann([](const Ponto2 &) { return 0.0; }));
  bcp.registra("laterais", nas_laterais, Contorno::neumann([](const Ponto2 &) { return 0.0; }));

  NavierStokes ns(m, fd, bcu0, bcu1, bcp);
  ns.reynolds(re);
  ns.passo_de_tempo(dt);
  ns.projecao(projecao_do_nome(disc_proj));
  ns.convectivo(convectivo_do_nome(disc_conv));
  ns.relogio_do_contorno(&tbc);
  ns.termo_fonte_faceta([re](const Ponto2 &x, int d, double t) {
    return caixa::fonte(x, d, t, re);
  });

  // A PRESSAO COMECA NA EXATA, e a velocidade no repouso -- que E' a exata em
  // `t = 0`, porque `f(0) = 0`.  Comecar `p` em zero somaria ao erro do
  // esquema um transiente de inicializacao que nao e' do esquema.
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      ns.pressao().em(m.id_global(b, i)) = caixa::p_exata(m.centroide(b, i), 0.0);

  const long passos = static_cast<long>(std::llround(T / dt));

  // `T` TEM DE SER MULTIPLO INTEIRO DE `dt`, e isto e' recusa e nao aviso.
  //
  // Com `dt ~ h^2`, um `T` qualquer da' numero de passos fracionario e cada
  // nivel para num INSTANTE DIFERENTE.  A solucao exata depende do tempo, de
  // modo que a comparacao passa a misturar erro de discretizacao com
  // diferenca de instante -- e a ordem medida sai errada sem que nada acuse.
  // Medido: com `T = 0,05` o nivel 3 parava em 0,050781.
  //
  // `T = 1/16` serve para os niveis 3, 4 e 5 (32, 128 e 512 passos), porque
  // `dt = 1/(8 . 4^n)`.
  if (std::fabs(passos * dt - T) > 1e-14 * std::max(1.0, T)) {
    if (rank == 0)
      std::printf(
          "T = %.17g NAO e' multiplo inteiro de dt = %.17g (daria %.6f\n"
          "passos).  Cada nivel pararia num instante diferente e a ordem\n"
          "medida seria a mistura do erro com a diferenca de instante.\n"
          "Use T = k/(8*4^nivel); para os niveis 3 a 5, T = 0.0625 serve.\n",
          T, dt, T / dt);
    return 2;
  }

  if (rank == 0) {
    std::printf("caso        caixa com tampa deslizante, solucao MANUFATURADA\n");
    std::printf("dominio     [0,1]^2, %ld celulas, h_min = %.5g, em transicao 2:1: %ld\n",
                m.num_folhas_globais(), hh, m.num_folhas_em_transicao());
    std::printf("Re          %.4g\n", re);
    std::printf("dt          %.5g  (= h^2 / %.3g; limite difusivo h^2 Re / 4)\n",
                dt, cfl_dif);
    std::printf("T           %.4g  em %ld passos\n", T, passos);
    std::printf("projecao    %s\n", nome_de(ns.projecao()).c_str());
    std::printf("convectivo  %s  (EXERCITADO: as duas componentes variam nas duas direcoes)\n",
                nome_de(ns.convectivo()).c_str());
    std::printf("tempo       %s\n", nome_de(ns.tempo()).c_str());
    std::printf("pressao     Neumann em toda a volta, Dirichlet em 1 faceta do fundo\n\n");
    std::fflush(stdout);
  }

  t8flow::SaidaVtk2d vtk(m);
  auto grava = [&](int quadro) {
    vtk.limpa();
    vtk.vetorial("velocidade", fd, ns.velocidade());
    vtk.escalar("pressao", ns.pressao());
    vtk.escreve("caixa", quadro);
  };
  int quadro = 0;
  if (cada_vtk > 0) grava(quadro++);

  double pior_compat = 0.0;
  for (long k = 0; k < passos; k++) {
    ns.passo();
    pior_compat = std::max(pior_compat, std::fabs(ns.soma_do_lado_direito()));
    if (!ns.pressao_convergiu()) {
      if (rank == 0)
        std::printf("PARANDO no passo %ld: o Poisson da projecao NAO convergiu\n", k + 1);
      return 1;
    }
    if (cada_vtk > 0 && (k + 1) % cada_vtk == 0) grava(quadro++);
  }
  grava(quadro);

  // ---- o erro contra a manufaturada, no instante final --------------------
  const double tf = ns.instante();
  double eu = 0.0, ev = 0.0, ep = 0.0;
  for (int d = 0; d < 2; d++)
    for (long f2 = 0; f2 < fd.num_facetas(d); f2++) {
      const double e = std::fabs(ns.velocidade().em(d, f2) -
                                 caixa::u_exata(fd.centro(d, f2), d, tf));
      if (d == 0) eu = std::max(eu, e); else ev = std::max(ev, e);
    }
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m.centroide(b, i);
      ep = std::max(ep, std::fabs(ns.pressao().em(m.id_global(b, i)) -
                                  caixa::p_exata(c, tf)));
    }

  if (rank == 0) {
    std::printf("t final  %.6f   (f(t) = %.6f)\n", tf, caixa::f(tf));
    std::printf("\nerro contra a solucao manufaturada, norma do maximo:\n");
    std::printf("  u   %.4e\n  v   %.4e\n  p   %.4e\n", eu, ev, ep);
    std::printf("  divergencia maxima            %.4e\n", ns.divergencia_maxima());
    std::printf("  pior |soma do lado direito|   %.4e\n", pior_compat);
    std::printf("  rebaixamentos do convectivo   %ld\n", ns.rebaixamentos_convectivos());
    std::printf(
        "\nLEITURA: a solucao e' TRIGONOMETRICA, logo nenhum esquema polinomial\n"
        "e' exato nela e o que se mede e' truncamento.  Com dt ~ h^2 o erro de\n"
        "tempo (primeira ordem) e o de espaco (segunda) caem juntos, e a ordem\n"
        "que se le' e' a do esquema.\n");
  }
  return 0;
}
