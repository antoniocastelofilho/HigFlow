// A GOTA ESTATICA: o oraculo de Laplace, `delta p = sigma / R`.
//
// E' o degrau B2 da escada -- o primeiro em que a frente CONVERSA com o
// solver.  E e' um oraculo de FORMA FECHADA: uma gota circular em repouso, sem
// gravidade, tem salto de pressao exatamente `sigma/R` atraves da interface, e
// velocidade identicamente nula.
//
// O CASO E' DE PROPRIEDADES IGUAIS nos dois fluidos, de proposito.  Com
// densidade e viscosidade iguais, a UNICA coisa que age e' a tensao
// superficial -- logo o que o erro mede e' a cadeia frente -> curvatura ->
// nucleo -> facetas -> projecao, e nada mais.  Variar as propriedades junto
// faria um erro medido ficar sem dono.
//
// AS CORRENTES ESPURIAS SAO O SEGUNDO ORACULO, e sao o artefato classico do
// metodo: a forca discreta nao e' exatamente um gradiente, de modo que a
// projecao nao a anula por completo e sobra um escoamento parasita junto a'
// interface.  Ele NAO vai a zero com o refino da malha -- escala com
// `sigma/mu` --, e por isso o caso o imprime em vez de o esconder: um numero
// que nao converge tem de ser visto, nao tolerado em silencio.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "acoplamento2d.hxx"
#include "arranque.hxx"
#include "campos2d.hxx"
#include "contorno2d.hxx"
#include "frente2d.hxx"
#include "malha2d.hxx"
#include "navierstokes2d.hxx"
#include "saida_vtk2d.hxx"

using namespace t8flow2d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0, np = 1;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &np);

  double sigma = 1.0, R = 0.25, re = 1.0, T = 0.05;
  int nivel = 5, nmarc = 0;
  long cada_vtk = 0;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--sigma=", 8)) sigma = std::atof(a + 8);
    else if (!std::strncmp(a, "--raio=", 7)) R = std::atof(a + 7);
    else if (!std::strncmp(a, "--re=", 5)) re = std::atof(a + 5);
    else if (!std::strncmp(a, "--T=", 4)) T = std::atof(a + 4);
    else if (!std::strncmp(a, "--nivel=", 8)) nivel = std::atoi(a + 8);
    else if (!std::strncmp(a, "--marcadores=", 13)) nmarc = std::atoi(a + 13);
    else if (!std::strncmp(a, "--vtk=", 6)) cada_vtk = std::atol(a + 6);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf("uso: ./exemplos/gota/gota-estatica [opcoes]\n\n"
                    "  --sigma        tensao superficial\n"
                    "  --raio         raio da gota\n"
                    "  --re           numero de Reynolds\n"
                    "  --T            instante final\n"
                    "  --nivel        refino do quadrado unitario\n"
                    "  --marcadores   marcadores da frente (0 = 4 por celula)\n"
                    "  --vtk=N        grava VTK a cada N passos\n");
      return 0;
    }
  }
  if (np > 1) {
    if (rank == 0) std::printf("ESTE CASO AINDA E' SERIAL (%d pedidos).\n", np);
    return 2;
  }

  Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  DominioFacetas fd(m);
  const double h = 1.0 / (1 << nivel);
  const double dt = h * h / 8.0;

  // QUATRO MARCADORES POR CELULA de interface, por padrao.  Menos que isso e a
  // frente nao resolve a celula; muito mais nao melhora a curvatura e so'
  // aumenta o custo da distribuicao.
  if (nmarc <= 0)
    nmarc = std::max(16, static_cast<int>(4.0 * 2.0 * M_PI * R / h));

  const Ponto2 centro{0.5, 0.5};
  Frente frente = Frente::circulo(centro, R, nmarc);

  // CAIXA FECHADA: velocidade nula nas quatro paredes.  A pressao fica
  // toda-Neumann e o solver recusa -- o pino de Dirichlet vai numa faceta do
  // fundo, LONGE da gota, onde a pressao exata e' a do fluido externo.
  const double E = 1e-9;
  auto zero = [](const Ponto2 &) { return 0.0; };
  ContornoPorRegiao bcu0, bcu1, bcp;
  bcu0.padrao(Contorno::dirichlet(zero));
  bcu1.padrao(Contorno::dirichlet(zero));
  bcp.registra("pino",
               [E, h](const Ponto2 &x) { return std::fabs(x[1]) < E && x[0] < h; },
               Contorno::dirichlet(zero));
  bcp.padrao(Contorno::neumann(zero));

  NavierStokes ns(m, fd, bcu0, bcu1, bcp);
  ns.reynolds(re);
  ns.passo_de_tempo(dt);

  CampoEscalar frac(m);
  CampoVetorial forca(fd);
  escreve_fracao(frente, m, &frac);
  ns.forca_de_faceta(&forca);

  // A PARTICAO DA UNIDADE, verificada ANTES de rodar: a soma de `fracao x
  // volume` tem de dar a area do poligono em precisao de maquina.  Se nao
  // der, nao adianta medir salto de pressao nenhum.
  double soma = 0.0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 hh = m.lados_da_folha(b, i);
      soma += frac.em(m.id_global(b, i)) * hh[0] * hh[1];
    }
  const double area_poli = frente.area();

  const long passos = static_cast<long>(std::llround(T / dt));
  if (rank == 0) {
    std::printf("caso        gota estatica -- oraculo de Laplace\n");
    std::printf("dominio     [0,1]^2, %ld celulas, h = %.5g\n",
                m.num_folhas_globais(), h);
    std::printf("gota        R = %.4g em (%.2f, %.2f), %d marcadores "
                "(%.1f por celula de interface)\n",
                R, centro[0], centro[1], nmarc, nmarc * h / (2 * M_PI * R));
    std::printf("sigma       %.4g      Re  %.4g\n", sigma, re);
    std::printf("dt          %.5g   T  %.4g em %ld passos\n", dt, T, passos);
    std::printf("propriedades IGUAIS nos dois fluidos: so' a tensao age\n\n");
    std::printf("particao da unidade: soma(c.V) = %.17g, area do poligono = "
                "%.17g, desvio %.3e\n\n", soma, area_poli,
                std::fabs(soma - area_poli));
    std::fflush(stdout);
  }

  t8flow::SaidaVtk2d vtk(m);
  auto grava = [&](int quadro) {
    vtk.limpa();
    vtk.vetorial("velocidade", fd, ns.velocidade());
    vtk.escalar("pressao", ns.pressao());
    vtk.escalar("fracao", frac);
    vtk.fronteira("frente", {frente.marcadores()});
    vtk.escreve("gota", quadro);
  };
  int quadro = 0;
  if (cada_vtk > 0) grava(quadro++);

  Ponto2 ftot{0.0, 0.0};
  for (long k = 0; k < passos; k++) {
    // A FRENTE NAO SE MOVE NESTE CASO: a gota e' estatica, e advecta-la com as
    // correntes espurias faria o raio mudar e o oraculo deixar de ser `sigma/R`
    // com o `R` que se conhece.  E' simplificacao DECLARADA: o que se verifica
    // aqui e' a cadeia curvatura -> nucleo -> projecao, nao o transporte, que
    // e' B1 e ja' foi verificado.
    forca_tensao(frente, m, fd, sigma, &forca);
    if (k == 0) ftot = forca_total(fd, forca);
    ns.passo();
    if (!ns.pressao_convergiu()) {
      if (rank == 0)
        std::printf("PARANDO no passo %ld: o Poisson nao convergiu\n", k + 1);
      return 1;
    }
    if (cada_vtk > 0 && (k + 1) % cada_vtk == 0) grava(quadro++);
  }
  grava(quadro);

  // ---- o salto de pressao ------------------------------------------------
  //
  // MEDIDO EM CELULAS PURAS, longe da interface: dentro da banda a pressao
  // esta' em transicao e a media ali misturaria os dois lados.  "Pura" e'
  // fracao acima de 0,99 ou abaixo de 0,01, e o centro a mais de `2h` da
  // interface -- a largura do nucleo.
  double p_in = 0.0, p_out = 0.0;
  long n_in = 0, n_out = 0;
  double umax = 0.0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const long id = m.id_global(b, i);
      const Ponto2 c = m.centroide(b, i);
      const double r = std::hypot(c[0] - centro[0], c[1] - centro[1]);
      if (std::fabs(r - R) < 2.5 * h) continue;      // banda do nucleo
      if (frac.em(id) > 0.99) { p_in += ns.pressao().em(id); n_in++; }
      else if (frac.em(id) < 0.01) { p_out += ns.pressao().em(id); n_out++; }
    }
  for (int d = 0; d < 2; d++)
    for (long q = 0; q < fd.num_facetas(d); q++)
      umax = std::max(umax, std::fabs(ns.velocidade().em(d, q)));

  if (rank == 0) {
    if (n_in == 0 || n_out == 0) {
      std::printf("SEM CELULAS PURAS de um dos lados (dentro %ld, fora %ld) -- "
                  "o oraculo nao mediu nada\n", n_in, n_out);
      return 1;
    }
    p_in /= n_in;
    p_out /= n_out;
    const double dp = p_in - p_out;
    const double exato = sigma / R;
    std::printf("forca total distribuida (tem de ser nula numa curva fechada):"
                " (%.3e, %.3e)\n", ftot[0], ftot[1]);
    std::printf("\ncelulas puras: %ld dentro, %ld fora\n", n_in, n_out);
    std::printf("  delta p calculado  %.6f\n", dp);
    std::printf("  delta p exato      %.6f   (sigma/R)\n", exato);
    std::printf("  erro relativo      %.3e\n", std::fabs(dp - exato) / exato);
    std::printf("  correntes espurias %.3e   (|u| maximo)\n", umax);
    std::printf("  numero capilar das espurias  %.3e   (mu |u| / sigma)\n",
                umax / (re * sigma));
    std::printf(
        "\nLEITURA: o salto tem forma fechada e o erro dele mede a cadeia\n"
        "curvatura -> nucleo -> projecao.  As correntes espurias NAO vao a\n"
        "zero com o refino -- sao o artefato do metodo, e estao aqui para\n"
        "serem vistas e nao toleradas em silencio.\n");
  }
  return 0;
}
