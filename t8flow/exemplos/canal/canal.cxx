// O CANAL: escoamento de Poiseuille, que tem SOLUCAO ANALITICA.
//
// E' o primeiro caso do solver de Navier-Stokes, e a escolha nao e' de
// conveniencia -- e' porque ele da' um oraculo de EXATIDAO e nao so' de ordem:
//
//     u(y) = (3/2) (1 - y^2),   v = 0,   dp/dx = -2 u_max / Re
//
// O perfil e' QUADRATICO em `y`, o laplaciano de tres pontos e' exato para
// quadraticas, o gradiente de pressao e' constante e o termo convectivo se
// ANULA IDENTICAMENTE (`u` nao depende de `x`, `v` e' nulo).  Logo o
// estacionario discreto tem de bater com o analitico em PRECISAO DE MAQUINA.
// Ordem de convergencia certa convive com constante errada; exatidao nao.
//
// E JUSTAMENTE POR ISSO O CANAL NAO VERIFICA O TERMO CONVECTIVO.  Central,
// upwind de primeira e CUBISTA dao o MESMO resultado aqui, e um caso que nao
// distingue os tres nao esta' exercitando nenhum.  Isso esta' dito e nao
// escondido: o convectivo precisa de outro caso, e ate' ele existir "CUBISTA
// esta' certo" nao e' afirmacao verificada.
//
// A geometria e os numeros sao os do `example2d_Newt` do HiGFlow, para que a
// comparacao entre os dois seja o mesmo caso e nao dois parecidos:
// `x` em [0,8], `y` em [-1,1], 160x40 celulas, Re = 1.

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
#include "navierstokes2d.hxx"
#include "saida_vtk2d.hxx"

using namespace t8flow2d;

namespace {

constexpr double U_MAX = 1.5;     //!< media unitaria no canal de meia-altura 1

double u_exata(const Ponto2 &p) { return U_MAX * (1.0 - p[1] * p[1]); }

//! `dp/dx = -2 u_max / Re`, e a pressao e' linear com zero na saida.
double p_exata(const Ponto2 &p, double re, double L)
{
  return (2.0 * U_MAX / re) * (L - p[0]);
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0, np = 1;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &np);

  double re = 1.0, dt = 1e-3, tol_estac = 1e-10;
  long cada_vtk = 0;                 //!< 0 = so' o quadro final
  double refina_parede = 0.0;        //!< espessura da faixa refinada na parede
  long max_passos = 200000;
  int nivel = 5;
  std::string disc_conv = "cubista", disc_proj = "incremental";

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--re=", 5)) re = std::atof(a + 5);
    else if (!std::strncmp(a, "--dt=", 5)) dt = std::atof(a + 5);
    else if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--passos=", 9)) max_passos = std::atol(a + 9);
    else if (!std::strncmp(a, "--tol=", 6)) tol_estac = std::atof(a + 6);
    else if (!std::strncmp(a, "--convectivo=", 13)) disc_conv = a + 13;
    else if (!std::strncmp(a, "--projecao=", 11)) disc_proj = a + 11;
    else if (!std::strncmp(a, "--vtk=", 6)) cada_vtk = std::atol(a + 6);
    else if (!std::strncmp(a, "--refina-parede=", 16))
      refina_parede = std::atof(a + 16);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf(
            "uso: mpirun -np N ./exemplos/canal/canal [opcoes]\n\n"
            "  --re           numero de Reynolds\n"
            "  --dt           passo de tempo\n"
            "  --nivel        refino (o dominio e' quadrado por bloco)\n"
            "  --passos       teto de passos\n"
            "  --tol          criterio de estacionariedade\n"
            "  --convectivo   central, upwind1, cubista, quick\n"
            "  --projecao     incremental, nao_incremental\n"
            "  --vtk=N        grava VTK a cada N passos (0 = so' o final)\n"
            "  --refina-parede=E  refina uma faixa de espessura E nas paredes\n");
      return 0;
    }
  }

  // O DOMINIO E' 8 x 2, e e' coberto por QUATRO BLOCOS QUADRADOS de lado 2.
  // Blocos quadrados porque a folha do t8code e' quadrada: um bloco 8x2 com
  // refino uniforme daria celula de proporcao 4:1, e a comparacao com o
  // HiGFlow deixaria de ser do mesmo caso.
  const double L = 8.0;
  std::vector<Bloco> blocos;
  for (int k = 0; k < 4; k++)
    blocos.push_back({{2.0 * k, -1.0}, {2.0 * (k + 1), 1.0}, nivel});
  Malha m(blocos, sc_MPI_COMM_WORLD);

  // A RECUSA DO PARALELO VEM AQUI, antes de `DominioFacetas`, e nao so' dentro
  // do solver: o dominio de facetas tambem percorre todos os blocos, e em
  // paralelo quebra primeiro -- com "arvore 2 nao e' local", que nao diz nada
  // a quem pediu dois processos.  A mensagem util tem de chegar antes do
  // sintoma.
  if (np > 1) {
    if (rank == 0)
      std::printf(
          "ESTE CASO AINDA E' SERIAL (%d processos pedidos).  Faltam a troca\n"
          "de franja dos campos -- o estencil le' celula e faceta do vizinho,\n"
          "e os campos cobrem so' os ids locais -- e o percurso por blocos,\n"
          "que supoe todas as arvores locais.  Rodar assim leria lixo perto\n"
          "da fronteira de particao e daria numero.\n", np);
    return 2;
  }

  // REFINO JUNTO A' PAREDE, por pedido: as duas faixas de altura `espessura`
  // coladas nas paredes sobem um nivel.  O solver RECUSA malha com transicao
  // 2:1 -- ver o comentario do construtor dele --, entao isto serve hoje para
  // mostrar a recusa, e nao para medir.
  if (refina_parede > 0.0) {
    m.refina_caixa({0.0, -1.0}, {L, -1.0 + refina_parede}, nivel + 1);
    m.refina_caixa({0.0, 1.0 - refina_parede}, {L, 1.0}, nivel + 1);
    if (rank == 0)
      std::printf("malha REFINADA junto as paredes: faixa de %.4g, nivel %d\n"
                  "celulas: %ld, em transicao 2:1: %ld\n\n",
                  refina_parede, nivel + 1, m.num_folhas_globais(),
                  m.num_folhas_em_transicao());
  }

  DominioFacetas fd(m);

  // AS CONDICOES SAO FUNCOES, e por regiao.  Entrada com o perfil analitico,
  // paredes sem escorregamento, saida com gradiente nulo; a pressao e' o
  // inverso disso -- Neumann onde a velocidade e' dada, Dirichlet onde ela
  // nao e'.  Essa troca nao e' escolha do caso: e' a condicao de contorno do
  // metodo de projecao.
  const double EPSG = 1e-9;
  ContornoPorRegiao bcu0, bcu1, bcp;

  auto na_entrada = [EPSG](const Ponto2 &x) { return std::fabs(x[0]) < EPSG; };
  auto na_saida = [EPSG, L](const Ponto2 &x) { return std::fabs(x[0] - L) < EPSG; };
  auto na_parede = [EPSG](const Ponto2 &x) {
    return std::fabs(std::fabs(x[1]) - 1.0) < EPSG;
  };

  bcu0.registra("entrada", na_entrada, Contorno::dirichlet(u_exata));
  bcu0.registra("saida", na_saida, Contorno::neumann([](const Ponto2 &) { return 0.0; }));
  bcu0.registra("parede", na_parede, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));

  bcu1.registra("entrada", na_entrada, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));
  bcu1.registra("saida", na_saida, Contorno::neumann([](const Ponto2 &) { return 0.0; }));
  bcu1.registra("parede", na_parede, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));

  bcp.registra("entrada", na_entrada, Contorno::neumann([](const Ponto2 &) { return 0.0; }));
  bcp.registra("saida", na_saida, Contorno::dirichlet([](const Ponto2 &) { return 0.0; }));
  bcp.registra("parede", na_parede, Contorno::neumann([](const Ponto2 &) { return 0.0; }));

  NavierStokes ns(m, fd, bcu0, bcu1, bcp);
  ns.reynolds(re);
  ns.passo_de_tempo(dt);
  ns.projecao(projecao_do_nome(disc_proj));
  ns.convectivo(convectivo_do_nome(disc_conv));

  if (rank == 0) {
    // A CONFIGURACAO INTEIRA, no alto.  Relatorio sem ela e' numero sem
    // procedencia, e neste projeto um binario sem o que se queria medir ja'
    // produziu corrida "estavel" que era outra coisa.
    std::printf("caso        canal (Poiseuille), solucao analitica\n");
    std::printf("dominio     [0,%.0f] x [-1,1], %ld celulas, h = %.4g\n", L,
                m.num_folhas_globais(), 2.0 / (1 << nivel));
    std::printf("Re          %.4g\n", re);
    std::printf("dt          %.4g\n", dt);
    std::printf("projecao    %s\n", nome_de(ns.projecao()).c_str());
    std::printf("convectivo  %s  (NAO exercitado neste caso: v = 0)\n",
                nome_de(ns.convectivo()).c_str());
    std::printf("tempo       %s\n", nome_de(ns.tempo()).c_str());
    std::printf("processos   %d\n\n", np);
    std::printf("  passo        t      residuo      div max   iters\n");
    std::fflush(stdout);
  }

  // A SAIDA DE VISUALIZACAO.  Um escritor, com os campos REGISTRADOS -- e o
  // registro e' refeito a cada quadro porque as referencias podem mudar.
  //
  // E E' VISUALIZACAO, NAO RESULTADO: os numeros que este programa afirma no
  // fim saem dos campos em memoria, com precisao cheia, e nao do arquivo.  O
  // VTK interpola o vetor das facetas para o centro, e essa interpolacao nao
  // e' a mesma com e sem vizinho local -- um oraculo construido sobre ele
  // mediria o erro da saida.
  t8flow::SaidaVtk2d vtk(m);
  auto grava = [&](int quadro) {
    vtk.limpa();
    vtk.vetorial("velocidade", fd, ns.velocidade());
    vtk.escalar("pressao", ns.pressao());
    vtk.escreve("canal", quadro);
  };
  int quadro = 0;
  grava(quadro++);

  long k = 0;
  for (; k < max_passos; k++) {
    ns.passo();
    if (cada_vtk > 0 && (k + 1) % cada_vtk == 0) grava(quadro++);
    if (!ns.pressao_convergiu()) {
      if (rank == 0)
        std::printf("PARANDO no passo %ld: o Poisson da projecao NAO convergiu\n",
                    k + 1);
      return 1;
    }
    if (rank == 0 && ((k + 1) % 500 == 0 || k == 0)) {
      std::printf("%7ld %8.4f %12.4e %12.4e %7d\n", k + 1, ns.instante(),
                  ns.residuo_estacionario(), ns.divergencia_maxima(),
                  ns.iteracoes_da_pressao());
      std::fflush(stdout);
    }
    if (ns.residuo_estacionario() < tol_estac) break;
  }

  // ---- o erro contra a analitica ------------------------------------------
  double eu = 0.0, ev = 0.0;
  for (long f = 0; f < fd.num_facetas(0); f++)
    eu = std::max(eu, std::fabs(ns.velocidade().em(0, f) - u_exata(fd.centro(0, f))));
  for (long f = 0; f < fd.num_facetas(1); f++)
    ev = std::max(ev, std::fabs(ns.velocidade().em(1, f)));

  // A PRESSAO E' COMPARADA A MENOS DE CONSTANTE?  NAO: a saida tem Dirichlet
  // zero, logo a constante esta' fixada e o valor absoluto e' comparavel.
  double ep = 0.0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m.centroide(b, i);
      ep = std::max(ep, std::fabs(ns.pressao().em(m.id_global(b, i)) -
                                  p_exata(c, re, L)));
    }

  // AS FAIXAS, e nao so' o erro: um erro maior que a propria solucao nao diz
  // ONDE esta' o problema, e as faixas dizem na hora se o campo esta' deslocado
  // por constante, escalado, ou simplesmente errado.
  double pmin = 1e300, pmax = -1e300, exmin = 1e300, exmax = -1e300;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m.centroide(b, i);
      const double v = ns.pressao().em(m.id_global(b, i));
      const double e = p_exata(c, re, L);
      pmin = std::min(pmin, v); pmax = std::max(pmax, v);
      exmin = std::min(exmin, e); exmax = std::max(exmax, e);
    }

  grava(quadro);

  if (rank == 0) {
    std::printf("\n  p calculada  [%.4e, %.4e]\n  p exata      [%.4e, %.4e]\n",
                pmin, pmax, exmin, exmax);
    std::printf("\nparou em %ld passos, t = %.4f, residuo %.3e\n", k + 1,
                ns.instante(), ns.residuo_estacionario());
    std::printf("\nerro contra a solucao analitica, norma do maximo:\n");
    std::printf("  u   %.4e\n", eu);
    std::printf("  v   %.4e   (a analitica e' identicamente nula)\n", ev);
    std::printf("  p   %.4e\n", ep);
    std::printf("  divergencia maxima  %.4e\n", ns.divergencia_maxima());
    std::printf("  rebaixamentos do convectivo  %ld\n",
                ns.rebaixamentos_convectivos());
    std::printf(
        "\nLEITURA: o perfil e' QUADRATICO e o esquema e' de segunda ordem,\n"
        "logo o erro esperado e' de MAQUINA -- nao uma taxa de convergencia.\n"
        "O que sobra e' residuo do solver e do criterio de parada.\n");
  }
  return 0;
}
