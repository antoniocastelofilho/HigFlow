// NAVIER-STOKES 3D POR SOLUCAO MANUFATURADA: a ordem, medida.
//
// Este e' o exemplo que transforma o solver 3D de "implementado" em
// "verificado".  Os quatro operadores tem casos proprios no arreio -- cada um
// com oraculo de exatidao --, mas nenhum deles diz se o PASSO os combina
// certo: a ordem no tempo, o lugar do contorno, o nivel da pressao.
//
// ===========================================================================
// O QUE O CASO EXIGE, E POR QUE
// ===========================================================================
//
// `T` TEM DE SER MULTIPLO INTEIRO DE `dt`, e o exemplo RECUSA quando nao e'.
// Senao cada nivel para num instante diferente e a ordem medida mistura erro
// de discretizacao com diferenca de instante -- que e' de primeira ordem e
// contamina tudo.  O mesmo tropeco ja' aconteceu no lado 2D.
//
// O `dt` E' O MESMO EM TODOS OS NIVEIS, de proposito.  Querendo a ordem no
// ESPACO, deixar `dt` cair com `h` mediria a soma das duas ordens e atribuiria
// o resultado ao espaco.  O `dt` fica fixo e pequeno o bastante para o erro
// temporal nao dominar -- e o exemplo imprime os dois erros para que a
// suposicao fique a' vista.
//
// A PRESSAO TEM UM PONTO DE DIRICHLET, com o valor da manufaturada.  O Poisson
// todo-Neumann e' singular e o solver o recusa; e com o valor conhecido a
// pressao sai em valor ABSOLUTO, nao a menos de constante.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "arranque.hxx"
#include "campos3d.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"
#include "manufaturada-ns3d.hxx"
#include "navierstokes3d.hxx"

using namespace t8flow3d;

namespace {

struct Erros {
  double u[3] = {0.0, 0.0, 0.0};
  double ui[3] = {0.0, 0.0, 0.0};   //!< so' o MIOLO, longe da parede
  double p = 0.0;
  long celulas = 0;
  int its = 0;
  long rebaixos = 0;
};

Erros corre(int nivel, double dt, double T, double re, Convectivo conv,
            Projecao proj, int rank)
{
  ManufaturadaNS3 ms;
  ms.re = re;

  Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  DominioFacetas fd(m);
  const double h = 1.0 / (1 << nivel);

  // O RELOGIO DO CONTORNO.  As condicoes dependem do tempo, e e' o SOLVER que
  // move este valor entre o preditor e a aplicacao final -- dentro do passo
  // sao dois instantes, e de fora so' se veem as fronteiras dele.
  double t_bc = 0.0;

  ContornoPorRegiao bc[3], bcp;
  for (int d = 0; d < 3; d++)
    bc[d].padrao(Contorno::dirichlet(
        [&ms, &t_bc, d](const Ponto3 &x) { return ms.u(x, d, t_bc); }));

  // O PINO: uma face de parede, com o valor da manufaturada.  A regiao usa o
  // `h` DA PAREDE -- com `h` errado ela pega zero faces (e o solver recusa) ou
  // duas (e o Poisson fica sobredeterminado).
  bcp.registra("pino",
               [h](const Ponto3 &x) {
                 return std::fabs(x[2]) < 1e-9 && x[0] < h && x[1] < h;
               },
               Contorno::dirichlet([&ms](const Ponto3 &x) { return ms.p(x); }));
  bcp.padrao(Contorno::neumann([](const Ponto3 &) { return 0.0; }));

  NavierStokes ns(m, fd, bc[0], bc[1], bc[2], bcp);
  ns.reynolds(re);
  ns.convectivo(conv);
  ns.projecao(proj);
  ns.passo_de_tempo(dt);
  ns.termo_fonte_faceta(
      [&ms](const Ponto3 &x, int d, double t) { return ms.f(x, d, t); });
  ns.relogio_do_contorno(&t_bc);

  // DADO INICIAL EXATO: a manufaturada em `t = 0`.
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd.num_facetas(d); f++)
      ns.velocidade().em(d, f) = ms.u(fd.centro(d, f), d, 0.0);
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      ns.pressao().em(m.id_global(b, i)) = ms.p(m.centroide(b, i));

  const long passos = static_cast<long>(std::llround(T / dt));

  // SINAL DE VIDA DENTRO DO NIVEL.  Sem ele o exemplo so' fala quando o nivel
  // FECHA, e um nivel de duas horas fica indistinguivel de um travado.
  //
  // Nao e' zelo: nesta mesma arvore um processo ficou 39 h bloqueado a 0,2% de
  // CPU parecendo que calculava, e no nivel 6 deste exemplo eu tive de contar
  // linhas de log do PETSc para saber em que passo estava -- inferencia a
  // partir de ruido de biblioteca, que e' a resposta errada para a pergunta
  // certa.
  const long cada = (passos >= 10) ? passos / 10 : 1;
  for (long k = 0; k < passos; k++) {
    ns.passo();
    if (rank == 0 && ((k + 1) % cada == 0 || k + 1 == passos)) {
      // LINHA DE VERDADE, e nao `\r`: a saida deste exemplo quase sempre vai
      // para ARQUIVO, e ali `\r` nao apaga coisa nenhuma -- fica um caractere
      // de controle no meio do texto e a linha de resultado sai grudada na de
      // progresso.  Dez linhas por nivel nao atrapalham ninguem; uma linha
      // corrompida atrapalha quem for ler depois.
      std::printf("...  nivel %d: passo %ld/%ld, t = %.5g, its %d\n", nivel,
                  k + 1, passos, ns.instante(), ns.iteracoes_da_pressao());
      std::fflush(stdout);
    }
  }

  Erros e;
  e.celulas = m.num_folhas_globais();
  e.its = ns.iteracoes_da_pressao();
  e.rebaixos = ns.rebaixamentos_convectivos();
  const double tf = ns.instante();
  // DOIS ERROS, E A DIFERENCA ENTRE ELES E' INFORMACAO.
  //
  // O fechamento de parede -- a reflexao `u_f = 2 g - u_d` -- e' de ordem
  // MENOR que o esquema interior.  Medir so' o maximo global esconde isso:
  // ve-se a ordem cair e atribui-se ao esquema.  Com as duas colunas, a
  // diferenca diz onde o erro mora.
  for (int d = 0; d < 3; d++)
    for (long f = 0; f < fd.num_facetas(d); f++) {
      const Ponto3 x = fd.centro(d, f);
      const double er =
          std::fabs(ns.velocidade().em(d, f) - ms.u(x, d, tf));
      e.u[d] = std::max(e.u[d], er);
      // A MARGEM DO MIOLO E' FISICA, e nao em contagem de celulas.
      //
      // Com `2h` o miolo e' uma regiao DIFERENTE a cada nivel -- em `h = 1/8`
      // sobra `[0,25; 0,75]`, em `h = 1/32` sobra `[0,0625; 0,9375]` --, e a
      // ordem medida compara erros de dominios distintos.  Deu 0,169 entre
      // dois niveis e 1,500 entre os seguintes, oscilacao que nao tem nada a
      // ver com o esquema.
      //
      // E' A MESMA ARMADILHA DA BANDA DE REFINO DO AMR, registrada em
      // `doc/relatorio-bifasico-ft-2d.md`: largura em numero de celulas muda
      // de tamanho fisico com o refino, e o estudo deixa de convergir por
      // razoes alheias a' discretizacao.  Custou uma corrida de 15 h naquela
      // vez; aqui custou uma leitura de tabela.
      const double margem = 0.15;
      bool borda = false;
      for (int q = 0; q < 3; q++)
        if (x[q] < margem || x[q] > 1.0 - margem) borda = true;
      if (!borda) e.ui[d] = std::max(e.ui[d], er);
    }
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      e.p = std::max(e.p, std::fabs(ns.pressao().em(m.id_global(b, i)) -
                                    ms.p(m.centroide(b, i))));
  return e;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);

  int n0 = 2, n1 = 4;
  double dt = 2.5e-4, T = 0.01, re = 10.0;
  std::string nome_conv = "cubista", nome_proj = "incremental";

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--niveis=", 9))
      std::sscanf(a + 9, "%d,%d", &n0, &n1);
    else if (!std::strncmp(a, "--dt=", 5)) dt = std::atof(a + 5);
    else if (!std::strncmp(a, "--T=", 4)) T = std::atof(a + 4);
    else if (!std::strncmp(a, "--re=", 5)) re = std::atof(a + 5);
    else if (!std::strncmp(a, "--conv=", 7)) nome_conv = a + 7;
    else if (!std::strncmp(a, "--proj=", 7)) nome_proj = a + 7;
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf("uso: ./exemplos/mms3d/ns3d-mms [opcoes]\n\n"
                    "  --niveis=A,B  do nivel A ao B (padrao 2,4)\n"
                    "  --dt, --T     passo e instante final\n"
                    "  --re=R        numero de Reynolds\n"
                    "  --conv=NOME   central, first_order, cubista, quick\n"
                    "  --proj=NOME   incremental ou nao_incremental\n");
      return 0;
    }
  }

  // A RECUSA: `T` MULTIPLO INTEIRO DE `dt`.
  const double passos_reais = T / dt;
  const long passos = static_cast<long>(std::llround(passos_reais));
  if (std::fabs(passos_reais - passos) > 1e-9 * passos_reais || passos < 1) {
    if (rank == 0)
      std::printf("T = %.6g NAO e' multiplo inteiro de dt = %.6g (daria "
                  "%.6f passos).\nCada nivel pararia num instante diferente e "
                  "a ordem medida misturaria erro de discretizacao com\n"
                  "diferenca de instante.  Use T = %.6g (%ld passos) ou "
                  "T = %.6g (%ld passos).\n",
                  T, dt, passos_reais, dt * passos, passos,
                  dt * (passos + 1), passos + 1);
    return 2;
  }

  const Convectivo conv = convectivo_do_nome(nome_conv);
  const Projecao proj = (nome_proj == "nao_incremental")
                            ? Projecao::NaoIncremental
                            : Projecao::Incremental;

  if (rank == 0) {
    std::printf("caso        Navier-Stokes 3D por solucao manufaturada\n");
    std::printf("dominio     [0,1]^3, Re = %.4g\n", re);
    std::printf("tempo       dt = %.4g, T = %.4g (%ld passos), Euler "
                "explicito\n", dt, T, passos);
    std::printf("esquemas    convectivo %s, projecao %s\n\n",
                nome_conv.c_str(), nome_proj.c_str());
    std::printf("%6s %10s %12s %12s %12s %12s %12s\n", "nivel", "celulas",
                "erro u", "miolo u", "erro w", "miolo w", "erro p");
    std::fflush(stdout);
  }

  std::vector<Erros> es;
  for (int nv = n0; nv <= n1; nv++) {
    const Erros e = corre(nv, dt, T, re, conv, proj, rank);
    es.push_back(e);
    if (rank == 0) {
      std::printf("%6d %10ld %12.4e %12.4e %12.4e %12.4e %12.4e\n", nv,
                  e.celulas, e.u[0], e.ui[0], e.u[2], e.ui[2], e.p);
      std::fflush(stdout);
    }
  }

  if (rank == 0 && es.size() >= 2) {
    std::printf("\n%6s %12s %12s %12s %12s %12s\n", "niveis", "ordem u",
                "miolo u", "ordem w", "miolo w", "ordem p");
    for (size_t k = 1; k < es.size(); k++) {
      const double ou = std::log2(es[k - 1].u[0] / es[k].u[0]);
      const double oiu = std::log2(es[k - 1].ui[0] / es[k].ui[0]);
      const double ow = std::log2(es[k - 1].u[2] / es[k].u[2]);
      const double oiw = std::log2(es[k - 1].ui[2] / es[k].ui[2]);
      const double op = std::log2(es[k - 1].p / es[k].p);
      std::printf("%3zu->%-2zu %12.3f %12.3f %12.3f %12.3f %12.3f\n",
                  k - 1 + static_cast<size_t>(n0), k + static_cast<size_t>(n0),
                  ou, oiu, ow, oiw, op);
    }
    std::printf(
        "\nLEITURA: a ordem aqui e' a do ESPACO, com `dt` FIXO em todos os\n"
        "niveis.  Se ela cair ao refinar, o erro temporal passou a dominar --\n"
        "nao e' o esquema espacial degradando.  Baixe `dt` e remeca.\n");
  }
  return 0;
}
