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
  double fracao = 0.0;          // 0 = malha uniforme
  std::string esq = "fluxo2";
  std::string disc;
  int mls_ordem = 3;
  int mls_maxpts = 40;
  bool quad_pm = false;
  unsigned long semente = 2026;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!std::strncmp(a, "--solucao=", 10)) nome_sol = a + 10;
    else if (!std::strncmp(a, "--contorno=", 11)) modo_bc = a + 11;
    else if (!std::strncmp(a, "--niveis=", 9)) niveis_txt = a + 9;
    else if (!std::strncmp(a, "--fechamento=", 13)) fech = a + 13;
    else if (!std::strncmp(a, "--tol=", 6)) tol = std::atof(a + 6);
    else if (!std::strncmp(a, "--aleatorio=", 12)) fracao = std::atof(a + 12);
    else if (!std::strncmp(a, "--semente=", 10)) semente = std::strtoul(a + 10, nullptr, 10);
    else if (!std::strncmp(a, "--esquema=", 10)) esq = a + 10;
    else if (!std::strncmp(a, "--discretizacao=", 16)) disc = a + 16;
    else if (!std::strncmp(a, "--mls-ordem=", 12)) mls_ordem = std::atoi(a + 12);
    else if (!std::strncmp(a, "--mls-pontos=", 13)) mls_maxpts = std::atoi(a + 13);
    else if (!std::strcmp(a, "--fonte-ponto-medio")) quad_pm = true;
    else if (!std::strcmp(a, "--ajuda") || !std::strcmp(a, "-h")) {
      std::printf("uso: poisson-mms [--solucao=N] [--contorno=dirichlet|misto]\n"
                  "                 [--niveis=3,4,5] [--fechamento=3|2] [--tol=1e-14]\n"
                  "                 [--aleatorio=FRACAO] [--semente=N]\n"
                  "                 [--discretizacao=NOME]  (preferido)\n"
                  "                 [--esquema=fluxo2|compacto4|lele4|mlsgauss]\n"
                  "                 [--mls-ordem=N] [--mls-pontos=N]\n"
                  "                 [--fonte-ponto-medio]\n\n"
                  "  --esquema      fluxo2    balanco de fluxos, 2a ordem, 5 pontos\n"
                  "                 mlsgauss  derivada MLS + Gauss de 2 pontos na face\n"
                  "                           (VOLUMES FINITOS: a incognita e' a media)\n"
                  "                 mlspontual  laplaciano MLS avaliado no centro\n"
                  "                           (DIFERENCAS FINITAS: a incognita e' o valor)\n"
                  "                           (o unico que passa de 2a ordem em malha\n"
                  "                            ADAPTADA; ordem e suporte por --mls-*)\n"
                  "                 compacto4 Mehrstellen, 4a ordem, 9 pontos (usa lap f)\n"
                  "                 lele4     compacto de Lele, 4a ordem, 9 pontos (usa f)\n\n"
                  "  --aleatorio=F  refina aleatoriamente uma fracao F das celulas (quadtree\n"
                  "                 adaptada, com transicao 2:1).  0 = malha uniforme.\n\n"
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
  std::printf("esquema     %s\n", esq.c_str());
  std::printf("fechamento  %s pontos na parede (so' no esquema de fluxo)\n",
              fech == "2" ? "DOIS" : "tres");
  std::printf("solver      tolerancia %.1e\n", tol);
  if (fracao > 0.0)
    std::printf("dominio     [0,1]x[0,1], quadtree ADAPTADA (fracao %.2f, semente %lu)\n\n",
                fracao, semente);
  else
    std::printf("dominio     [0,1]x[0,1], malha UNIFORME\n\n");

  std::vector<Linha> tab;
  for (int nivel : niveis) {
    // COMM_WORLD, e nao COMM_SELF: e' o que faz este exemplo ser o mesmo
    // programa em np=1 e em np>1.  Com COMM_SELF cada processo resolveria o
    // problema INTEIRO sozinho e o paralelo nao seria exercitado -- e os numeros
    // sairiam iguais, o que faria o teste de consistencia entre np passar sem
    // medir nada.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
    if (fracao > 0.0) {
      // UMA passada, e com a MESMA semente em todos os niveis: assim as malhas da
      // sequencia tem o mesmo PADRAO de adaptacao, e a taxa medida e' de refino de
      // malha e nao de mudanca de padrao.  Padroes diferentes por nivel dariam uma
      // "ordem" que mistura as duas coisas.
      m.refina_aleatorio(fracao, semente, nivel + 1);
    }

    Poisson po(m, bc);
    po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
    po.fechamento(fp);

    // `--discretizacao=` E' A FORMA PREFERIDA, e ela instala a familia e a ordem
    // de uma vez: `vf-mlsgauss4` e' "volumes finitos de quarta ordem", nomeado.
    // O `--esquema=` antigo continua, para quem ja' o usa.
    if (!disc.empty()) {
      const Discretizacao d = Discretizacao::do_nome(disc, mls_ordem);
      po.discretizacao(d);
      if (d.variante() == Variante::CompactoMehrstellen)
        po.fonte_laplaciano([&sol](const Ponto2 &p) { return sol.lap_f(p); });
      if (fracao > 0.0 && !d.estavel_em_malha_adaptada())
        std::printf("AVISO: \"%s\" nao e' estavel em malha ADAPTADA -- ver a "
                    "clausula W6 do contrato.  Seguindo porque foi pedido.\n",
                    d.nome().c_str());
    } else if (esq == "compacto4") {
      po.esquema(Esquema::Compacto4);
      po.fonte_laplaciano([&sol](const Ponto2 &p) { return sol.lap_f(p); });
    } else if (esq == "lele4") {
      po.esquema(Esquema::Compacto4Lele);
    } else if (esq == "mlspontual" || esq == "df-mlspontual") {
      po.esquema(Esquema::MlsPontual);
      po.mls(mls_ordem, mls_maxpts);
    } else if (esq == "mlsgauss") {
      po.esquema(Esquema::MlsGauss);
      po.mls(mls_ordem, mls_maxpts);
      if (quad_pm) po.quadratura_fonte(Poisson::QuadFonte::PontoMedio);
    } else if (esq != "fluxo2") {
      throw std::invalid_argument("esquema: \"" + esq +
                                  "\" nao existe.  Validos: fluxo2, compacto4, "
                                  "lele4, mlsgauss, mlspontual");
    }
    po.resolve(Motor::Petsc, tol);
    if (esq == "mlsgauss" && po.truncados_pela_franja() > 0)
      std::printf("AVISO: a FRANJA truncou %ld suportes -- os pesos perto das\n"
                  "       fronteiras de particao nao sao os do serial, logo o\n"
                  "       resultado DEPENDE do numero de processos.\n",
                  po.truncados_pela_franja());
    if (esq == "mlsgauss" && po.rebaixamentos_mls() > 0)
      std::printf("AVISO: %ld dos %ld ajustes MLS REBAIXARAM a ordem -- o esquema\n"
                  "       nao entrega a ordem pedida em todo o dominio.\n",
                  po.rebaixamentos_mls(), po.ajustes_mls());

    Linha l;
    l.nivel = nivel;
    // CONTAGEM GLOBAL, nao a da particao: em np=3 as locais sao 86, 85 e 85, e
    // imprimir a local faria a tabela parecer tres malhas diferentes.
    l.celulas = m.num_folhas_globais();
    l.h = 1.0 / (1 << nivel);      // h do nivel BASE; a malha adaptada tem h/2 em parte
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
  //
  // O GRAU DE EXATIDAO DEPENDE DO ESQUEMA, e dizer "grau 2" para todos faria o
  // texto MENTIR exatamente onde ele serve para algo.  Medido:
  //
  //   fluxo2              grau 2 em malha uniforme; grau 1 em malha adaptada
  //                       (a quadratura do ponto medio na faceta e' o limite)
  //   compacto4 / lele4   grau 3 (o fantasma de quatro pontos na parede e' o limite)
  //   mlsgauss            grau = ordem do ajuste MLS, em uniforme E adaptada
  const int grau_exato =
      (esq == "mlsgauss" || esq == "mlspontual") ? mls_ordem
      : (esq == "compacto4" || esq == "lele4")   ? 3
      : (fracao > 0.0)                           ? 1
                                                  : 2;

  std::printf("\n");
  if (sol.grau >= 0 && (fp == FechamentoParede::TresPontos || esq != "fluxo2") &&
      sol.grau <= grau_exato) {
    std::printf("LEITURA: o esquema e' EXATO para esta solucao (grau %d <= %d, que e' o\n"
                "grau de exatidao de \"%s\" aqui), logo o erro esperado e' de MAQUINA --\n"
                "e o que sobra e' residuo do solver linear, nao do esquema.  Nao ha'\n"
                "ordem de convergencia a medir, e e' por isso que a coluna diz \"exato\".\n",
                sol.grau, grau_exato, esq.c_str());
  } else if (sol.grau < 0) {
    std::printf("LEITURA: solucao nao-polinomial, logo ha' erro de truncamento e a\n"
                "ordem de convergencia e' mensuravel -- esperada ~%d.\n",
                (esq == "mlsgauss" || esq == "mlspontual") ? mls_ordem
                                                           : (esq == "fluxo2" ? 2 : 4));
  } else {
    std::printf("LEITURA: o esquema NAO e' exato para esta combinacao (grau %d,\n"
                "fechamento de %s pontos), logo o erro cai com h e a ordem e'\n"
                "mensuravel.  Isso nao e' defeito: e' o contraste que mostra que o\n"
                "fechamento esta' influenciando o resultado.  Grau de exatidao de\n"
                "\"%s\" aqui: %d.\n",
                sol.grau, fech == "2" ? "DOIS" : "tres", esq.c_str(), grau_exato);
  }
  return 0;
}
