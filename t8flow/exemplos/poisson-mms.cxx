// O PRIMEIRO EXEMPLO DO T8FLOW: Poisson pelo metodo das solucoes manufaturadas.
//
// Serial, parametrizado, e entra na suite -- a pedido do usuario.  Troca-se a
// solucao, as condicoes de contorno, o fechamento de parede e os niveis de malha
// SEM RECOMPILAR, porque um exemplo que exige recompilar para mudar o caso vira
// dois exemplos divergentes na primeira vez que alguem precisa do outro caso.
//
//     ./exemplos/poisson-mms --solucao=grau2 --niveis=3,4,5 --contorno=dirichlet
//
// O QUE ELE MEDE, e a diferenca entre as duas coisas e' o ponto do metodo:
//
//   ERRO DE MAQUINA quando o esquema e' EXATO para aquela solucao.  Com
//   `grau2` e fechamento de tres pontos, o balanco de fluxos reproduz a solucao
//   exatamente: o erro que sobra e' do SOLVER LINEAR, nao do esquema, e por isso o
//   residuo e' impresso ao lado -- ele e' o piso do que se pode exigir.
//
//   ORDEM DE CONVERGENCIA quando nao e'.  Com `trig`, o erro cai com `h` e a taxa
//   `log2(e_h / e_{h/2})` tem de dar ~2.  Com polinomio exato nao HA' taxa para
//   medir, porque o erro ja' e' zero -- pedir ordem de convergencia de um caso
//   exato produz divisao de ruido por ruido, e e' por isso que a tabela diz
//   "exato" em vez de imprimir um numero sem sentido.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "manufaturada2d.hxx"
#include "poisson2d.hxx"

using namespace t8flow2d;

namespace {

//! Normal EXTERIOR do quadrado unitario no ponto da parede.
void normal_em(const Ponto2 &p, double n[2])
{
  const double eps = 1e-9;
  n[0] = n[1] = 0.0;
  if (std::fabs(p[0]) < eps) n[0] = -1.0;
  else if (std::fabs(p[0] - 1.0) < eps) n[0] = 1.0;
  else if (std::fabs(p[1]) < eps) n[1] = -1.0;
  else if (std::fabs(p[1] - 1.0) < eps) n[1] = 1.0;
}

//! As quatro paredes do quadrado unitario, como regioes nomeadas.
//
// Nomeadas porque o nome aparece no diagnostico: "nenhuma regiao cobre (1.000,
// 0.375)" e' acionavel, "contorno sem condicao" nao e'.
ContornoPorRegiao monta_contorno(const Manufaturada &sol, const std::string &modo)
{
  const double eps = 1e-9;
  const Manufaturada *s = &sol;
  const auto dado_u = [s](const Ponto2 &p) { return s->u(p); };
  const auto dado_dn = [s](const Ponto2 &p) {
    double n[2];
    normal_em(p, n);
    return s->dudn(p, n);
  };

  ContornoPorRegiao r;
  const auto sul   = [eps](const Ponto2 &p) { return std::fabs(p[1]) < eps; };
  const auto norte = [eps](const Ponto2 &p) { return std::fabs(p[1] - 1.0) < eps; };
  const auto oeste = [eps](const Ponto2 &p) { return std::fabs(p[0]) < eps; };
  const auto leste = [eps](const Ponto2 &p) { return std::fabs(p[0] - 1.0) < eps; };

  if (modo == "dirichlet") {
    r.registra("sul", sul, Contorno::dirichlet(dado_u));
    r.registra("norte", norte, Contorno::dirichlet(dado_u));
    r.registra("oeste", oeste, Contorno::dirichlet(dado_u));
    r.registra("leste", leste, Contorno::dirichlet(dado_u));
  } else if (modo == "misto") {
    // Dirichlet em duas paredes e Neumann nas outras duas.  NAO so' Neumann: o
    // problema puramente de Neumann e' SINGULAR -- a solucao fica a menos de uma
    // constante --, e medir erro contra a exata ali mediria a constante, nao o
    // esquema.
    r.registra("sul", sul, Contorno::dirichlet(dado_u));
    r.registra("norte", norte, Contorno::dirichlet(dado_u));
    r.registra("oeste", oeste, Contorno::neumann(dado_dn));
    r.registra("leste", leste, Contorno::neumann(dado_dn));
  } else if (modo == "umdirichlet") {
    // UMA parede Dirichlet e as outras TRES Neumann.  E' o minimo para o problema
    // ser bem-posto: com Neumann em TODAS as paredes ele e' SINGULAR -- a solucao
    // fica a menos de uma constante, o solver nao tem o que fixar, e medir erro
    // contra a exata mediria a constante, nao o esquema.  Uma unica parede de
    // Dirichlet fixa o nivel.
    r.registra("sul", sul, Contorno::dirichlet(dado_u));
    r.registra("norte", norte, Contorno::neumann(dado_dn));
    r.registra("oeste", oeste, Contorno::neumann(dado_dn));
    r.registra("leste", leste, Contorno::neumann(dado_dn));
  } else {
    throw std::invalid_argument("contorno: \"" + modo +
                                "\" nao existe.  Validos: dirichlet, misto, umdirichlet");
  }
  return r;
}

struct Linha {
  int nivel = 0;
  long celulas = 0;
  double h = 0.0;
  double e_inf = 0.0;
  double e_l2 = 0.0;
  int iteracoes = 0;
  double residuo = 0.0;
  bool convergiu = false;
};

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  std::string nome_sol = "grau2";
  std::string modo_bc = "dirichlet";
  std::string niveis_txt = "3,4,5";
  std::string fech = "3";
  double tol = 1e-14;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--solucao=", 10)) nome_sol = a + 10;
    else if (!std::strncmp(a, "--contorno=", 11)) modo_bc = a + 11;
    else if (!std::strncmp(a, "--niveis=", 9)) niveis_txt = a + 9;
    else if (!std::strncmp(a, "--fechamento=", 13)) fech = a + 13;
    else if (!std::strncmp(a, "--tol=", 6)) tol = std::atof(a + 6);
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      std::printf("uso: poisson-mms [--solucao=N] [--contorno=dirichlet|misto]\n"
                  "                 [--niveis=3,4,5] [--fechamento=3|2] [--tol=1e-14]\n\n"
                  "solucoes disponiveis:\n");
      for (const std::string &n : manufaturadas_disponiveis())
        std::printf("  %-8s %s\n", n.c_str(), manufaturada(n).descricao.c_str());
      return 0;
    }
  }

  std::vector<int> niveis;
  {
    std::string t = niveis_txt;
    size_t p = 0;
    while (!t.empty()) {
      p = t.find(',');
      niveis.push_back(std::atoi(t.substr(0, p).c_str()));
      if (p == std::string::npos) break;
      t = t.substr(p + 1);
    }
  }

  const Manufaturada &sol = manufaturada(nome_sol);
  const ContornoPorRegiao bc = monta_contorno(sol, modo_bc);
  const FechamentoParede fp = (fech == "2") ? FechamentoParede::DoisPontos
                                            : FechamentoParede::TresPontos;

  // A CONFIGURACAO VAI PARA A SAIDA, sempre.  Relatorio sem a configuracao e'
  // numero sem procedencia, e neste projeto build sem o que se queria medir ja'
  // produziu duas corridas "estaveis" que eram outra coisa.
  std::printf("=== Poisson 2D, solucoes manufaturadas ===\n");
  std::printf("solucao     %s  --  %s\n", sol.nome.c_str(), sol.descricao.c_str());
  std::printf("contorno    %s (%zu regioes)\n", modo_bc.c_str(), bc.num_regioes());
  std::printf("fechamento  %s pontos na parede\n", fech == "2" ? "DOIS" : "tres");
  std::printf("solver      tolerancia %.1e\n", tol);
  std::printf("dominio     [0,1]x[0,1], malha UNIFORME\n\n");

  std::vector<Linha> tab;
  for (int nivel : niveis) {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_SELF);

    Poisson po(m, bc);
    po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
    po.fechamento(fp);
    po.resolve(Motor::Petsc, tol);

    Linha l;
    l.nivel = nivel;
    l.celulas = m.num_folhas();
    l.h = 1.0 / (1 << nivel);
    l.e_inf = po.erro_inf([&sol](const Ponto2 &p) { return sol.u(p); });
    l.e_l2 = po.erro_l2([&sol](const Ponto2 &p) { return sol.u(p); });
    l.iteracoes = po.iteracoes();
    l.residuo = po.residuo();
    l.convergiu = po.convergiu();
    tab.push_back(l);
  }

  std::printf("%5s %8s %10s %13s %7s %13s %7s %7s %11s\n",
              "nivel", "celulas", "h", "erro inf", "ordem", "erro L2", "ordem",
              "iters", "residuo");
  for (size_t k = 0; k < tab.size(); k++) {
    const Linha &l = tab[k];
    char oi[16] = "     -", ol[16] = "     -";
    if (k > 0) {
      // A ORDEM SO' E' IMPRESSA SE HOUVER ERRO PARA MEDIR.  Com esquema exato o
      // erro e' ruido do solver, e `log2(ruido/ruido)` e' um numero aleatorio que
      // pareceria uma taxa -- e' assim que se "mede" ordem de convergencia de um
      // caso que nao tem nenhuma.
      const Linha &a = tab[k - 1];
      if (a.e_inf > 1e-11 && l.e_inf > 0.0)
        std::snprintf(oi, sizeof oi, "%6.2f", std::log2(a.e_inf / l.e_inf));
      else
        std::snprintf(oi, sizeof oi, " exato");
      if (a.e_l2 > 1e-11 && l.e_l2 > 0.0)
        std::snprintf(ol, sizeof ol, "%6.2f", std::log2(a.e_l2 / l.e_l2));
      else
        std::snprintf(ol, sizeof ol, " exato");
    }
    std::printf("%5d %8ld %10.3e %13.4e %7s %13.4e %7s %7d %11.2e%s\n",
                l.nivel, l.celulas, l.h, l.e_inf, oi, l.e_l2, ol,
                l.iteracoes, l.residuo, l.convergiu ? "" : "  NAO CONVERGIU");
  }

  // A LEITURA, dita em vez de deixada para quem le'.
  std::printf("\n");
  if (sol.grau >= 0 && fp == FechamentoParede::TresPontos && sol.grau <= 2) {
    std::printf("LEITURA: o esquema e' EXATO para esta solucao (grau %d, fechamento de\n"
                "tres pontos), logo o erro esperado e' de MAQUINA -- e o que sobra e'\n"
                "residuo do solver linear, nao do esquema.  Nao ha' ordem de\n"
                "convergencia a medir aqui, e e' por isso que a coluna diz \"exato\".\n",
                sol.grau);
  } else if (sol.grau < 0) {
    std::printf("LEITURA: solucao nao-polinomial, logo ha' erro de truncamento e a\n"
                "ordem de convergencia e' mensuravel -- esperada ~2.\n");
  } else {
    std::printf("LEITURA: o esquema NAO e' exato para esta combinacao (grau %d,\n"
                "fechamento de %s pontos), logo o erro cai com h e a ordem e'\n"
                "mensuravel.  Isso nao e' defeito: e' o contraste que mostra que o\n"
                "fechamento esta' influenciando o resultado.\n",
                sol.grau, fech == "2" ? "DOIS" : "tres");
  }
  return 0;
}
