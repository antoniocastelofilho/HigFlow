// O Poisson em tres dimensoes, verificado pelo metodo das solucoes manufaturadas.
//
// E' ONDE TUDO O QUE FOI CONSTRUIDO SE ENCONTRA: a malha de hexaedros, a franja
// geometrica, o suporte por caixa de raio 2, a medida de faceta como AREA, e a
// transicao 2:1 com faces partidas em quatro.
//
// O ORACULO PRINCIPAL E' A EXATIDAO POLINOMIAL, e nao a ordem de convergencia:
// onde o esquema e' exato o erro e' de maquina, e isso e' uma afirmacao muito
// mais forte -- a ordem pode estar certa com a constante errada, a exatidao nao.
//
// E CADA AFIRMACAO DE EXATIDAO VEM COM UM CONTRASTE: a mesma medida onde ela
// DEVE falhar.  Sem o par, "o erro e' de maquina" nao distingue um esquema
// correto de um oraculo cego.

#include "../testing.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "arranque.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"
#include "manufaturada3d.hxx"
#include "poisson3d.hxx"

using namespace t8flow3d;

namespace {

//! O piso do erro e' o SOLVER, e nao o esquema.  Com tolerancia 1e-12 o residuo
//! que sobra fica na faixa de 1e-13 a 1e-10; pedir menos seria exigir do solver
//! algo que ele nao promete.
constexpr double MAQUINA = 1e-9;

struct Resultado {
  double e_inf = 0.0;
  double e_l2 = 0.0;
  double resid_op = 0.0;
  long celulas = 0;
  long rebaixos = 0;
  long truncados = 0;
  bool convergiu = false;
  int iteracoes = 0;
};

Resultado roda(const std::string &sol_nome, int nivel, const std::string &disc,
               double fracao = 0.0)
{
  const Manufaturada &sol = manufaturada(sol_nome);
  Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nivel}}, sc_MPI_COMM_WORLD);
  if (fracao > 0.0) m.refina_aleatorio(fracao, 2026, nivel + 1);

  ContornoPorRegiao bc;
  bc.padrao(Contorno::dirichlet([&sol](const Ponto3 &p) { return sol.u(p); }));

  Poisson po(m, bc);
  po.fonte([&sol](const Ponto3 &p) { return sol.f(p); });
  po.discretizacao(Discretizacao::do_nome(disc));
  po.resolve(Motor::Petsc, 1e-12);

  Resultado r;
  r.e_inf = po.erro_inf([&sol](const Ponto3 &p) { return sol.u(p); });
  r.e_l2 = po.erro_l2([&sol](const Ponto3 &p) { return sol.u(p); });
  r.resid_op = po.residuo_da_exata([&sol](const Ponto3 &p) { return sol.u(p); });
  r.celulas = m.num_folhas_globais();
  r.rebaixos = po.rebaixamentos_mls();
  r.truncados = po.truncados_pela_franja();
  r.convergiu = po.convergiu();
  r.iteracoes = po.iteracoes();
  return r;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("o_fluxo_de_segunda_ordem_e_exato_para_campo_linear");
  {
    // GRAU 1 E' O MINIMO que se exige de um esquema consistente, e em volumes
    // finitos ele e' exato: o fluxo de um campo linear pela faceta e' calculado
    // sem erro pela diferenca entre os centros.
    const Resultado r = roda("grau1", 2, "vf-fluxo2");
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu: %d iteracoes", r.iteracoes);
    T_BELOW(r.e_inf, MAQUINA, "erro inf do grau1 com vf-fluxo2");
    T_BELOW(r.resid_op, MAQUINA, "residuo do operador do grau1");
    printf("  (grau1 fluxo2: %ld celulas, inf %.3e, residuo do operador %.3e)\n",
           r.celulas, r.e_inf, r.resid_op);
  }

  t_case("o_fluxo_de_segunda_ordem_converge_em_segunda_ordem");
  {
    // E O CONTRASTE do caso acima: no grau 2 o esquema deixa de ser exato em
    // tres dimensoes -- o fluxo pela faceta usa a diferenca entre centros, que
    // erra a curvatura -- e converge em ~2.
    const Resultado a = roda("grau2", 2, "vf-fluxo2");
    const Resultado b = roda("grau2", 3, "vf-fluxo2");
    T_CHECK_MSG(a.e_inf > 1e-5,
                "o grau2 deu %.3e com fluxo2 -- se fosse de maquina, o caso do "
                "grau1 nao estaria medindo nada", a.e_inf);
    const double ordem = std::log2(a.e_inf / b.e_inf);
    T_CHECK_MSG(ordem > 1.5 && ordem < 2.5,
                "ordem %.2f para grau2 com fluxo2 (esperava ~2)", ordem);
    printf("  (grau2 fluxo2: %.3e -> %.3e, ordem %.2f)\n", a.e_inf, b.e_inf, ordem);
  }

  t_case("o_mlsgauss_e_exato_ate_o_grau_da_ordem_do_ajuste");
  {
    // A PROPRIEDADE QUE DISTINGUE ESTE ESQUEMA: o grau de exatidao acompanha a
    // ordem do ajuste.  Com ordem 2, exato ate' grau 2 -- onde o fluxo de
    // segunda ordem ja' erra.
    const Resultado r = roda("grau2", 2, "vf-mlsgauss2");
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu: %d iteracoes", r.iteracoes);
    T_CHECK_MSG(r.rebaixos == 0,
                "%ld ajustes MLS rebaixaram a ordem -- a exatidao afirmada abaixo "
                "nao seria da ordem pedida em todo o dominio", r.rebaixos);
    T_BELOW(r.e_inf, MAQUINA, "erro inf do grau2 com vf-mlsgauss2");

    // O CONTRASTE E' O MESMO GRAU COM O OUTRO ESQUEMA, na mesma malha.
    const Resultado f2 = roda("grau2", 2, "vf-fluxo2");
    T_CONTRASTE(r.e_inf > 1e-16 ? r.e_inf : 1e-16, f2.e_inf, 1e4,
                "grau 2: mlsgauss contra o fluxo de segunda ordem");
    printf("  (grau2: mlsgauss2 %.3e, fluxo2 %.3e)\n", r.e_inf, f2.e_inf);
  }

  t_case("o_mlsgauss_de_ordem_tres_e_exato_no_grau_tres");
  {
    // ORDEM 3 EM TRES DIMENSOES pede 20 monomios e suporte por caixa de raio 2,
    // que por sua vez pede tres camadas de franja.  Este caso e' o que exercita
    // a cadeia inteira.
    const Resultado r = roda("grau3", 2, "vf-mlsgauss3");
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu: %d iteracoes", r.iteracoes);
    T_CHECK_MSG(r.rebaixos == 0, "%ld ajustes rebaixaram a ordem", r.rebaixos);
    T_CHECK_MSG(r.truncados == 0,
                "a franja truncou %ld suportes -- os pesos nao seriam os do "
                "serial", r.truncados);
    T_BELOW(r.e_inf, MAQUINA, "erro inf do grau3 com vf-mlsgauss3");
    T_BELOW(r.resid_op, MAQUINA, "residuo do operador do grau3");

    // CONTRASTE: o grau 4 esta' acima da ordem, e o esquema erra nele.
    const Resultado g4 = roda("grau4", 2, "vf-mlsgauss3");
    T_CONTRASTE(r.e_inf > 1e-16 ? r.e_inf : 1e-16, g4.e_inf, 1e3,
                "mlsgauss3: grau 3 contra grau 4");
    printf("  (mlsgauss3: grau3 %.3e, grau4 %.3e)\n", r.e_inf, g4.e_inf);
  }

  t_case("na_malha_adaptada_o_mlsgauss_mantem_a_exatidao");
  {
    // A RAZAO DE EXISTIR DO ESQUEMA.  Em malha adaptada a face de um hexaedro e'
    // partida em QUATRO, e a quadratura de sub-faceta tem de somar os quatro
    // fluxos -- integrar a face inteira de uma vez daria o fluxo medio, que nao
    // e' a soma dos quatro.
    const Resultado r = roda("grau2", 2, "vf-mlsgauss2", 0.4);
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu em malha adaptada: %d iteracoes",
                r.iteracoes);
    T_CHECK_MSG(r.rebaixos == 0, "%ld ajustes rebaixaram em malha adaptada",
                r.rebaixos);
    T_BELOW(r.e_inf, MAQUINA, "erro inf do grau2 em malha ADAPTADA");

    // O CONTRASTE E' O FLUXO NA MESMA MALHA: ele perde a exatidao do grau 2 na
    // transicao, como em duas dimensoes.
    const Resultado f2 = roda("grau2", 2, "vf-fluxo2", 0.4);
    T_CHECK_MSG(f2.e_inf > 1e-5,
                "o fluxo2 deu %.3e no grau2 adaptado -- se fosse de maquina, a "
                "malha nao estaria de fato adaptada", f2.e_inf);
    printf("  (adaptada, %ld celulas: mlsgauss2 %.3e, fluxo2 %.3e)\n",
           r.celulas, r.e_inf, f2.e_inf);
  }

  t_case("a_solucao_trigonometrica_converge_na_ordem_do_ajuste");
  {
    // NAO POLINOMIAL, logo ha' ordem a medir -- num polinomio do grau da
    // exatidao o erro e' de maquina e a razao divide ruido por ruido.
    const Resultado a = roda("trig", 2, "vf-mlsgauss2");
    const Resultado b = roda("trig", 3, "vf-mlsgauss2");
    const double ordem = std::log2(a.e_inf / b.e_inf);
    T_CHECK_MSG(ordem > 1.5,
                "ordem %.2f para a trigonometrica com mlsgauss2 (%.3e -> %.3e)",
                ordem, a.e_inf, b.e_inf);
    printf("  (trig mlsgauss2: %.3e -> %.3e, ordem %.2f)\n", a.e_inf, b.e_inf, ordem);
  }

  t_case("o_poisson_3d_recusa_discretizacao_que_nao_e_de_volumes_finitos");
  {
    // O ESCOPO E' DECLARADO E VERIFICADO: o lado 3D traz so' volumes finitos, e
    // pedir outra coisa tem de falhar dizendo por que -- nao produzir numero.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1}}, sc_MPI_COMM_WORLD);
    ContornoPorRegiao bc;
    bc.padrao(Contorno::dirichlet([](const Ponto3 &) { return 0.0; }));
    Poisson po(m, bc);

    // As variantes de diferencas finitas nem existem no `enum` do lado 3D, logo
    // a recusa acontece ja' no nome.
    bool lancou = false;
    std::string msg;
    try { Discretizacao::do_nome("df-lele4"); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "\"df-lele4\" foi aceito num lado que so' tem volumes finitos");
    T_CHECK_MSG(msg.find("vf-mlsgauss") != std::string::npos,
                "a mensagem nao traz a lista de validos: \"%s\"", msg.c_str());

    // E contorno vazio e' recusado na construcao: sem condicao em parede nenhuma
    // o sistema e' singular.
    bool sem_bc = false;
    try { Poisson p2(m, ContornoPorRegiao()); }
    catch (const std::invalid_argument &) { sem_bc = true; }
    T_CHECK_MSG(sem_bc, "contorno vazio foi aceito");
  }

  t_case("o_operador_nao_depende_da_particao");
  {
    // O ORACULO QUE FALTAVA, e cuja ausencia deixou a suite VERDE com o paralelo
    // quebrado: este arquivo so' rodava em np=1, onde nao ha' franja, de modo
    // que nenhum caso tocava o caminho que estava errado.
    //
    // O QUE SE AFIRMA AQUI NAO PRECISA DE REFERENCIA SERIAL -- e' propriedade
    // do operador, verificavel em qualquer numero de processos:
    //
    //   UM LAPLACIANO ANIQUILA CONSTANTES, logo a linha de uma celula INTERIOR
    //   soma zero.  Quando a franja fazia o fechamento tomar uma face de
    //   particao por parede do dominio, o peso saia da linha e ia para o lado
    //   direito: a soma da linha dava -1,907e-02 onde o serial dava 1e-17.  O
    //   erro contra a solucao exata mudava 11% com `np`, e nada acusava;
    //
    //   O SUPORTE PEDIDO TEM DE CHEGAR INTEIRO.  `truncados_pela_franja`
    //   diferente de zero significa pesos diferentes dos do serial perto das
    //   fronteiras de particao -- resultado que depende de `np`.
    //
    // O primeiro oraculo que tentei NAO servia: `residuo_da_exata` e' um
    // MAXIMO, dominado pelas linhas de parede, e ficou identico em 1,1591e+00
    // nos quatro numeros de processos enquanto o interior divergia.
    const Manufaturada &sol = manufaturada("grau2");
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_WORLD);
    ContornoPorRegiao bc;
    bc.padrao(Contorno::dirichlet([&sol](const Ponto3 &p) { return sol.u(p); }));
    Poisson po(m, bc);
    po.fonte([&sol](const Ponto3 &p) { return sol.f(p); });
    po.discretizacao(Discretizacao::do_nome("vf-mlsgauss2"));
    po.despeja_operador(std::string("/tmp/t8flow-op-teste"));
    po.resolve(Motor::Petsc, 1e-12);

    T_CHECK_MSG(po.truncados_pela_franja() == 0,
                "%ld suportes truncados pela franja: os pesos perto das "
                "fronteiras de particao nao sao os do serial",
                po.truncados_pela_franja());

    char nome[256];
    std::snprintf(nome, sizeof nome, "/tmp/t8flow-op-teste.%d", m.processo());
    std::FILE *f = std::fopen(nome, "r");
    T_CHECK_MSG(f != nullptr, "o despejo do operador nao foi escrito");
    if (f) {
      long pior_id = -1;
      double pior = 0.0;
      long linhas = 0;
      char buf[1 << 16];
      while (std::fgets(buf, sizeof buf, f)) {
        long id = 0, nnz = 0;
        double rhs = 0.0;
        char *q = buf;
        id = std::strtol(q, &q, 10);
        q += 4;                               // " rhs"
        rhs = std::strtod(q, &q);
        q += 4;                               // " nnz"
        nnz = std::strtol(q, &q, 10);
        (void) rhs;
        (void) nnz;

        // SO' AS CELULAS INTERIORES: na parede a condicao de Dirichlet move
        // peso para o lado direito e a soma da linha NAO e' zero -- e deve
        // mesmo nao ser.
        const long idx = m.indice_do_id(id);
        if (idx < 0) continue;
        bool interior = true;
        for (int fa = 0; fa < 6; fa++)
          if (m.vizinhos_do_indice(idx, fa).empty()) interior = false;
        if (!interior) continue;

        double soma = 0.0;
        for (;;) {
          while (*q == ' ') q++;
          // O FIM DE LINHA TERMINA O LACO, e nao so' o `\0`.  Sem isto o
          // `strtol` nao consome o `\n`, nao avanca, e o laco NUNCA SAI -- o
          // teste trava em vez de falhar, e travar e' pior que falhar.
          if (*q == '\0' || *q == '\n') break;
          char *antes = q;
          std::strtol(q, &q, 10);             // coluna
          if (*q == ':') q++;
          soma += std::strtod(q, &q);
          if (q == antes) break;              // nao avancou: linha malformada
        }
        linhas++;
        if (std::fabs(soma) > std::fabs(pior)) { pior = soma; pior_id = id; }
      }
      std::fclose(f);
      std::remove(nome);
      T_CHECK_MSG(linhas > 0, "nenhuma celula interior neste processo");
      T_CHECK_MSG(std::fabs(pior) < 1e-10,
                  "a linha interior %ld soma %.3e, e um laplaciano aniquila "
                  "constantes: ha' peso saindo da linha -- fechamento de "
                  "contorno aplicado numa face de PARTICAO",
                  pior_id, pior);
    }
  }

  return t_end();
}
