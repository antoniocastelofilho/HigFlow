// O POISSON EM TRES DIMENSOES, verificado pelo metodo das solucoes manufaturadas.
//
// Parametrizado para trocar a solucao, a discretizacao e a malha sem recompilar
// -- e a configuracao inteira vai IMPRESSA no alto da saida, porque relatorio
// sem a configuracao e' numero sem procedencia.
//
// ONDE O ERRO E' DE MAQUINA NAO HA' ORDEM A MEDIR, e a tabela diz "exato" em vez
// de imprimir um numero: com erro de maquina em dois niveis consecutivos, a
// razao divide ruido por ruido e produz qualquer valor.
//
// Uso:  mpirun -np N ./exemplos/poisson3d/poisson3d-mms [opcoes]

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "arranque.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"
#include "manufaturada3d.hxx"
#include "poisson3d.hxx"

using namespace t8flow3d;

namespace {

struct Linha {
  int nivel = 0;
  long celulas = 0;
  double h = 0.0;
  double e_inf = 0.0;
  double e_l2 = 0.0;
  int iteracoes = 0;
  double residuo = 0.0;
  bool convergiu = false;
  long rebaixos = 0;
  long truncados = 0;
};

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);
  int rank = 0;
  sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);

  std::string nome_sol = "grau2";
  std::string disc = "vf-fluxo2";
  std::string niveis_txt = "2,3";
  double tol = 1e-12;
  double fracao = 0.0;
  unsigned long semente = 2026;
  int mls_pontos = 200;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--solucao=", 10)) nome_sol = a + 10;
    else if (!std::strncmp(a, "--discretizacao=", 16)) disc = a + 16;
    else if (!std::strncmp(a, "--niveis=", 9)) niveis_txt = a + 9;
    else if (!std::strncmp(a, "--tol=", 6)) tol = std::atof(a + 6);
    else if (!std::strncmp(a, "--aleatorio=", 12)) fracao = std::atof(a + 12);
    else if (!std::strncmp(a, "--semente=", 10)) semente = std::strtoul(a + 10, nullptr, 10);
    else if (!std::strncmp(a, "--mls-pontos=", 13)) mls_pontos = std::atoi(a + 13);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      if (rank == 0)
        std::printf(
            "uso: mpirun -np N ./exemplos/poisson3d/poisson3d-mms [opcoes]\n\n"
            "  --solucao        %s\n"
            "  --discretizacao  %s\n"
            "  --niveis         lista, p.ex. 2,3,4\n"
            "  --aleatorio=F    quadtree adaptada: fracao F de celulas refinadas\n"
            "  --semente        a do refino aleatorio\n"
            "  --tol            tolerancia do solver\n"
            "  --mls-pontos     teto do suporte do ajuste\n",
            nomes_manufaturadas().c_str(),
            Discretizacao::nomes_validos().c_str());
      return 0;
    }
  }

  std::vector<int> niveis;
  {
    std::string t = niveis_txt;
    size_t p = 0;
    while ((p = t.find(',')) != std::string::npos) {
      niveis.push_back(std::atoi(t.substr(0, p).c_str()));
      t = t.substr(p + 1);
    }
    if (!t.empty()) niveis.push_back(std::atoi(t.c_str()));
  }

  const Manufaturada &sol = manufaturada(nome_sol);
  const Discretizacao d = Discretizacao::do_nome(disc);

  int np = 1;
  sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &np);
  if (rank == 0) {
    // A CONFIGURACAO INTEIRA, no alto: relatorio sem ela e' numero sem
    // procedencia, e neste projeto um build sem o que se queria medir ja'
    // produziu corridas "estaveis" que eram outra coisa.
    std::printf("dimensao    3 (HEX)\n");
    std::printf("solucao     %s -- %s (grau %d)\n", sol.nome.c_str(),
                sol.descricao.c_str(), sol.grau);
    std::printf("discretiz.  %s (familia %s, ordem %d)\n", d.nome().c_str(),
                d.familia() == Familia::VolumesFinitos ? "volumes finitos"
                                                       : "diferencas finitas",
                d.ordem());
    std::printf("malha       %s\n", fracao > 0.0 ? "ADAPTADA (refino aleatorio)"
                                                 : "uniforme");
    if (fracao > 0.0) std::printf("            fracao %.2f, semente %lu\n", fracao, semente);
    std::printf("processos   %d\n", np);
    std::printf("tolerancia  %.1e\n\n", tol);
    std::printf("nivel  celulas          h      erro inf   ordem       erro L2"
                "   ordem   iters     residuo\n");
    std::fflush(stdout);
  }

  std::vector<Linha> tabela;
  for (int nivel : niveis) {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
    if (fracao > 0.0) m.refina_aleatorio(fracao, semente, nivel + 1);

    ContornoPorRegiao bc;
    bc.padrao(Contorno::dirichlet([&sol](const Ponto3 &p) { return sol.u(p); }));

    Poisson po(m, bc);
    po.fonte([&sol](const Ponto3 &p) { return sol.f(p); });
    po.discretizacao(d);
    po.mls(mls_pontos);
    po.resolve(Motor::Petsc, tol);

    Linha l;
    l.nivel = nivel;
    l.celulas = m.num_folhas_globais();
    l.h = 1.0 / (1 << nivel);
    l.e_inf = po.erro_inf([&sol](const Ponto3 &p) { return sol.u(p); });
    l.e_l2 = po.erro_l2([&sol](const Ponto3 &p) { return sol.u(p); });
    l.iteracoes = po.iteracoes();
    l.residuo = po.residuo();
    l.convergiu = po.convergiu();
    l.rebaixos = po.rebaixamentos_mls();
    l.truncados = po.truncados_pela_franja();
    tabela.push_back(l);

    if (rank == 0 && (l.rebaixos > 0 || l.truncados > 0))
      std::printf("AVISO nivel %d: %ld ajustes REBAIXARAM a ordem, %ld suportes "
                  "truncados pela franja -- o esquema nao entrega a ordem pedida "
                  "em todo o dominio\n", nivel, l.rebaixos, l.truncados);
  }

  if (rank != 0) return 0;

  for (size_t k = 0; k < tabela.size(); k++) {
    const Linha &l = tabela[k];
    char oi[16] = "      -", ol[16] = "      -";
    if (k > 0) {
      const Linha &a = tabela[k - 1];
      // ONDE O ESQUEMA E' EXATO NAO HA' ORDEM A MEDIR: com erro de maquina em
      // dois niveis, a razao divide ruido por ruido.  Imprimir o numero
      // sugeriria uma taxa onde ha' so' aritmetica de ponto flutuante.
      if (a.e_inf > 1e-11 && l.e_inf > 1e-11)
        std::snprintf(oi, sizeof oi, "%7.2f", std::log2(a.e_inf / l.e_inf));
      else
        std::snprintf(oi, sizeof oi, "  exato");
      if (a.e_l2 > 1e-11 && l.e_l2 > 1e-11)
        std::snprintf(ol, sizeof ol, "%7.2f", std::log2(a.e_l2 / l.e_l2));
      else
        std::snprintf(ol, sizeof ol, "  exato");
    }
    std::printf("%5d %8ld %10.3e %13.4e %7s %13.4e %7s %7d %11.2e%s\n",
                l.nivel, l.celulas, l.h, l.e_inf, oi, l.e_l2, ol,
                l.iteracoes, l.residuo, l.convergiu ? "" : "  NAO CONVERGIU");
  }

  // A LEITURA, dita em vez de deixada para quem le'.
  const int grau_ex = d.grau_exato(fracao > 0.0);
  std::printf("\n");
  if (sol.grau >= 0 && sol.grau <= grau_ex)
    std::printf("LEITURA: o esquema e' EXATO para esta solucao (grau %d <= %d, que e'\n"
                "o grau de exatidao de \"%s\" nesta malha), logo o erro esperado e' de\n"
                "MAQUINA -- e o que sobra e' residuo do solver, nao do esquema.\n",
                sol.grau, grau_ex, d.nome().c_str());
  else if (sol.grau < 0)
    std::printf("LEITURA: solucao nao-polinomial, logo ha' erro de truncamento e a\n"
                "ordem de convergencia e' mensuravel -- esperada ~%d.\n", d.ordem());
  else
    std::printf("LEITURA: o esquema NAO e' exato para esta combinacao (grau %d acima\n"
                "do grau de exatidao %d de \"%s\" nesta malha), logo o erro cai com h.\n",
                sol.grau, grau_ex, d.nome().c_str());
  return 0;
}
