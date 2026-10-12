// POISSON COM SOLUCOES MANUFATURADAS: o primeiro exemplo do sistema, afirmado.
//
// O METODO: escolhe-se `u` conhecido, DERIVA-SE `f = lap(u)` e os dados de
// contorno a partir dele, e mede-se o erro contra o `u` escolhido.  A forca esta'
// em `f` e nos dados serem EXATOS -- se fossem calculados numericamente, o oraculo
// mediria a consistencia entre dois erros meus, nao o esquema.
//
// O QUE ESTE ARQUIVO AFIRMA, e por que cada caso existe:
//
//   EXATIDAO   grau 2 com fechamento de tres pontos: o esquema de balanco de
//              fluxos reproduz a solucao, logo o erro tem de ser de MAQUINA nas
//              duas normas.  E' o pedido do usuario.
//   CONTRASTE  o mesmo grau 2 com fechamento de DOIS pontos NAO e' exato -- erro
//              1e-2, dez ordens de grandeza acima.  Sem esse par, "erro de
//              maquina" nao distingue esquema exato de oraculo que nao mede nada:
//              um sistema montado como identidade tambem daria 1e-13.
//   CONVERGENCIA  solucao trigonometrica: ha' erro de truncamento, e a ordem tem de
//              dar ~2.  Com polinomio exato nao HA' taxa para medir.
//   NEUMANN    grau 2 com Neumann em duas paredes tambem e' exato, porque no
//              balanco de fluxos o dado de Neumann E' o fluxo -- nenhuma
//              aproximacao.
//
// A TOLERANCIA de "erro de maquina" e' 1e-10, e o numero nao e' arbitrario: o erro
// medido e' ~3e-13 e o RESIDUO do solver e' ~5e-15.  O piso e' o solver, nao o
// esquema -- exigir 1e-15 seria exigir do solver linear algo que ele nao promete.
#include "../testing.h"

#include <cmath>
#include <vector>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "manufaturada2d.hxx"
#include "poisson2d.hxx"

using namespace t8flow2d;

namespace {

void normal_em(const Ponto2 &p, double n[2])
{
  const double eps = 1e-9;
  n[0] = n[1] = 0.0;
  if (std::fabs(p[0]) < eps) n[0] = -1.0;
  else if (std::fabs(p[0] - 1.0) < eps) n[0] = 1.0;
  else if (std::fabs(p[1]) < eps) n[1] = -1.0;
  else if (std::fabs(p[1] - 1.0) < eps) n[1] = 1.0;
}

//! Quantas paredes de Dirichlet: 4 (todas), 2 (sul e norte), ou 1 (so' o sul).
//
// UMA e' o MINIMO para o problema ser bem-posto: com Neumann em TODAS as paredes
// ele e' SINGULAR -- a solucao fica a menos de uma constante, o solver nao tem o
// que fixar, e medir erro contra a exata mediria a constante e nao o esquema.
enum class Paredes { TodasDirichlet, DuasNeumann, TresNeumann };

ContornoPorRegiao paredes(const Manufaturada &sol, Paredes modo)
{
  const double eps = 1e-9;
  const Manufaturada *s = &sol;
  const auto du = [s](const Ponto2 &p) { return s->u(p); };
  const auto dn = [s](const Ponto2 &p) {
    double n[2]; normal_em(p, n); return s->dudn(p, n);
  };
  ContornoPorRegiao r;
  // O sul e' SEMPRE Dirichlet: e' a parede que fixa o nivel da solucao.
  r.registra("sul",   [eps](const Ponto2 &p) { return std::fabs(p[1]) < eps; },
             Contorno::dirichlet(du));
  const bool norte_dir = (modo == Paredes::TodasDirichlet || modo == Paredes::DuasNeumann);
  const bool lados_dir = (modo == Paredes::TodasDirichlet);
  r.registra("norte", [eps](const Ponto2 &p) { return std::fabs(p[1] - 1.0) < eps; },
             norte_dir ? Contorno::dirichlet(du) : Contorno::neumann(dn));
  r.registra("oeste", [eps](const Ponto2 &p) { return std::fabs(p[0]) < eps; },
             lados_dir ? Contorno::dirichlet(du) : Contorno::neumann(dn));
  r.registra("leste", [eps](const Ponto2 &p) { return std::fabs(p[0] - 1.0) < eps; },
             lados_dir ? Contorno::dirichlet(du) : Contorno::neumann(dn));
  return r;
}

struct Resultado {
  double e_inf = 0.0, e_l2 = 0.0, residuo = 0.0;
  bool convergiu = false;
  int iteracoes = 0;
};

Resultado roda(const std::string &nome_sol, int nivel, FechamentoParede fp, Paredes modo)
{
  const Manufaturada &sol = manufaturada(nome_sol);
  Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_SELF);
  Poisson po(m, paredes(sol, modo));
  po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
  po.fechamento(fp);
  po.resolve(Motor::Petsc, 1e-14);

  Resultado r;
  r.e_inf = po.erro_inf([&sol](const Ponto2 &p) { return sol.u(p); });
  r.e_l2 = po.erro_l2([&sol](const Ponto2 &p) { return sol.u(p); });
  r.residuo = po.residuo();
  r.convergiu = po.convergiu();
  r.iteracoes = po.iteracoes();
  return r;
}

}  // namespace

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const double MAQUINA = 1e-10;   // ver o cabecalho: o piso e' o residuo do solver

  // ---- O PEDIDO: grau 2, Dirichlet em tudo, erro de maquina nas duas normas --

  t_case("grau2_com_dirichlet_da_erro_de_maquina_nas_duas_normas");
  double g2_inf = 0.0, g2_l2 = 0.0;
  {
    const Resultado r = roda("grau2", 4, FechamentoParede::TresPontos, Paredes::TodasDirichlet);
    g2_inf = r.e_inf;
    g2_l2 = r.e_l2;
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu: %d iteracoes, residuo %.3g",
                r.iteracoes, r.residuo);
    T_CHECK_MSG(r.e_inf < MAQUINA, "erro na norma do maximo: %.4e (limite %.0e)",
                r.e_inf, MAQUINA);
    T_CHECK_MSG(r.e_l2 < MAQUINA, "erro na norma 2: %.4e (limite %.0e)",
                r.e_l2, MAQUINA);
    printf("  (grau2 Dirichlet, nivel 4: inf %.3e, L2 %.3e, residuo %.2e)\n",
           r.e_inf, r.e_l2, r.residuo);
  }

  t_case("a_exatidao_do_grau2_nao_depende_do_refino");
  {
    // Esquema EXATO nao melhora nem piora com h -- o erro continua sendo o do
    // solver.  Se caisse como h^2, o esquema nao seria exato e o caso de cima
    // estaria passando por a malha ser fina, nao por exatidao.
    const Resultado a = roda("grau2", 3, FechamentoParede::TresPontos, Paredes::TodasDirichlet);
    const Resultado b = roda("grau2", 5, FechamentoParede::TresPontos, Paredes::TodasDirichlet);
    T_CHECK_MSG(a.e_inf < MAQUINA && b.e_inf < MAQUINA,
                "nivel 3: %.3e, nivel 5: %.3e -- esperava maquina nos dois",
                a.e_inf, b.e_inf);
    // E a razao entre eles NAO pode ser de ordem 4 (= h^2 com h/4): seria
    // convergencia, nao exatidao.
    const double razao = (b.e_inf > 0.0) ? a.e_inf / b.e_inf : 1.0;
    T_CHECK_MSG(razao < 4.0 && razao > 0.25,
                "o erro mudou por fator %.1f entre os niveis 3 e 5 -- isso e' "
                "convergencia, nao exatidao", razao);
  }

  // ---- O CONTRASTE que da' sentido ao numero de cima --------------------------

  t_case("o_fechamento_de_dois_pontos_perde_a_exatidao");
  {
    const Resultado r = roda("grau2", 4, FechamentoParede::DoisPontos, Paredes::TodasDirichlet);
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu com dois pontos");
    T_CHECK_MSG(r.e_inf > 1e-4,
                "com fechamento de DOIS pontos o erro foi %.4e -- pequeno demais. "
                "Se ele tambem fosse de maquina, o caso de cima nao distinguiria "
                "esquema exato de sistema montado como identidade", r.e_inf);
    T_CONTRASTE(g2_inf > 1e-16 ? g2_inf : 1e-16, r.e_inf, 1e6,
                "tres pontos contra dois pontos, norma do maximo");
    T_CONTRASTE(g2_l2 > 1e-16 ? g2_l2 : 1e-16, r.e_l2, 1e6,
                "tres pontos contra dois pontos, norma 2");
    printf("  (grau2 com DOIS pontos: inf %.3e, L2 %.3e -- contra %.3e e %.3e)\n",
           r.e_inf, r.e_l2, g2_inf, g2_l2);
  }

  t_case("o_fechamento_de_dois_pontos_converge_em_segunda_ordem");
  {
    // Nao e' exato, mas tambem nao esta' errado: tem de convergir em ~2.  Sem este
    // caso, "o de dois pontos erra mais" nao distinguiria segunda ordem de um
    // esquema inconsistente.
    const Resultado a = roda("grau2", 3, FechamentoParede::DoisPontos, Paredes::TodasDirichlet);
    const Resultado b = roda("grau2", 5, FechamentoParede::DoisPontos, Paredes::TodasDirichlet);
    const double ordem = std::log2(a.e_inf / b.e_inf) / 2.0;   // dois refinos
    T_CHECK_MSG(ordem > 1.7 && ordem < 2.3,
                "ordem medida %.2f entre os niveis 3 e 5 (esperava ~2)", ordem);
    printf("  (dois pontos: ordem %.2f)\n", ordem);
  }

  // ---- NEUMANN: exato, porque o dado E' o fluxo ------------------------------

  t_case("grau2_com_neumann_em_duas_paredes_tambem_da_erro_de_maquina");
  {
    const Resultado r = roda("grau2", 4, FechamentoParede::TresPontos, Paredes::DuasNeumann);
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu no caso misto");
    T_CHECK_MSG(r.e_inf < MAQUINA && r.e_l2 < MAQUINA,
                "misto: inf %.4e, L2 %.4e (limite %.0e)", r.e_inf, r.e_l2, MAQUINA);
    printf("  (grau2 com Neumann em duas paredes: inf %.3e, L2 %.3e)\n",
           r.e_inf, r.e_l2);
  }

  // ---- UMA parede Dirichlet e TRES Neumann -----------------------------------
  //
  // Pedido do usuario, e e' o caso de contorno mais exigente que o problema admite:
  // com Neumann nas QUATRO o problema e' singular.

  t_case("grau2_com_uma_parede_dirichlet_e_tres_neumann_da_erro_de_maquina");
  {
    const Resultado r = roda("grau2", 4, FechamentoParede::TresPontos,
                             Paredes::TresNeumann);
    T_CHECK_MSG(r.convergiu, "o solver nao convergiu com tres paredes de Neumann: "
                "%d iteracoes, residuo %.3g", r.iteracoes, r.residuo);
    T_CHECK_MSG(r.e_inf < MAQUINA && r.e_l2 < MAQUINA,
                "tres Neumann: inf %.4e, L2 %.4e (limite %.0e)",
                r.e_inf, r.e_l2, MAQUINA);
    printf("  (grau2, uma parede Dirichlet e tres Neumann: inf %.3e, L2 %.3e)\n",
           r.e_inf, r.e_l2);
  }

  t_case("com_tres_neumann_o_grau3_converge_em_segunda_ordem");
  {
    // E AQUI ESTA' A EXPLICACAO DE POR QUE O GRAU 2 E' EXATO E O GRAU 3 NAO, e ela
    // vale guardar porque e' contra-intuitiva: o dado de Neumann e' o fluxo EXATO
    // no CENTRO da faceta, mas o balanco de volumes finitos precisa do fluxo
    // MEDIO sobre ela.  A quadratura do ponto medio e' exata para integrando
    // LINEAR:
    //
    //   grau 2  ->  du/dn varia LINEARMENTE ao longo da faceta  ->  exato
    //   grau 3  ->  du/dn varia quadraticamente                 ->  erro O(h^2)
    //
    // O mesmo vale para a fonte: `f(c)*area` e' exato enquanto `f` for linear.  Por
    // isso a exatidao para o grau 2 nao e' sorte nem coincidencia de coeficientes:
    // e' a quadratura do ponto medio sendo exata nos DOIS lugares.
    const Resultado a = roda("grau3", 3, FechamentoParede::TresPontos, Paredes::TresNeumann);
    const Resultado b = roda("grau3", 5, FechamentoParede::TresPontos, Paredes::TresNeumann);
    const double ordem = std::log2(a.e_inf / b.e_inf) / 2.0;
    T_CHECK_MSG(ordem > 1.7 && ordem < 2.3,
                "ordem medida %.2f para grau3 com tres Neumann (esperava ~2)", ordem);
    T_CHECK_MSG(b.e_inf > 1e-6,
                "o erro do nivel 5 foi %.3e, perto do ruido -- a ordem medida ali "
                "nao significaria nada", b.e_inf);
    printf("  (grau3 com tres Neumann: ordem %.2f, erro %.3e -> %.3e)\n",
           ordem, a.e_inf, b.e_inf);
  }

  t_case("com_tres_neumann_a_trigonometrica_tambem_converge_em_segunda_ordem");
  {
    const Resultado a = roda("trig", 4, FechamentoParede::TresPontos, Paredes::TresNeumann);
    const Resultado b = roda("trig", 5, FechamentoParede::TresPontos, Paredes::TresNeumann);
    const double oi = std::log2(a.e_inf / b.e_inf);
    const double ol = std::log2(a.e_l2 / b.e_l2);
    T_CHECK_MSG(oi > 1.8 && oi < 2.2, "ordem na norma do maximo: %.2f", oi);
    T_CHECK_MSG(ol > 1.8 && ol < 2.2, "ordem na norma 2: %.2f", ol);
    printf("  (trig com tres Neumann: ordem inf %.2f, L2 %.2f)\n", oi, ol);
  }

  // ---- CONVERGENCIA, onde ha' erro para medir --------------------------------

  t_case("a_solucao_trigonometrica_converge_em_segunda_ordem");
  {
    const Resultado a = roda("trig", 4, FechamentoParede::TresPontos, Paredes::TodasDirichlet);
    const Resultado b = roda("trig", 5, FechamentoParede::TresPontos, Paredes::TodasDirichlet);
    const double oi = std::log2(a.e_inf / b.e_inf);
    const double ol = std::log2(a.e_l2 / b.e_l2);
    T_CHECK_MSG(oi > 1.8 && oi < 2.2, "ordem na norma do maximo: %.2f", oi);
    T_CHECK_MSG(ol > 1.8 && ol < 2.2, "ordem na norma 2: %.2f", ol);
    // E o erro tem de ser GRANDE o suficiente para a taxa significar algo: medir
    // ordem em cima de ruido do solver da' numero aleatorio.
    T_CHECK_MSG(b.e_inf > 1e-6,
                "o erro do nivel 5 foi %.3e -- perto do ruido do solver, e a ordem "
                "medida ali nao significaria nada", b.e_inf);
    printf("  (trig: ordem inf %.2f, L2 %.2f; erro %.3e -> %.3e)\n",
           oi, ol, a.e_inf, b.e_inf);
  }

  // ---- as recusas ------------------------------------------------------------

  t_case("fonte_nao_instalada_e_recusada");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    Poisson po(m, paredes(manufaturada("grau1"), Paredes::TodasDirichlet));
    bool lancou = false;
    try { po.resolve(); } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "resolver sem termo de fonte passou sem reclamar");
  }

  t_case("contorno_vazio_e_recusado");
  {
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    bool lancou = false;
    try { Poisson po(m, ContornoPorRegiao{}); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "Poisson sem condicao de contorno passou sem reclamar");
  }

  t_case("malha_com_transicao_2_para_1_e_recusada");
  {
    // A montagem e' de malha UNIFORME: na face grossa o fluxo e' a soma dos dois
    // lados finos, e a linha deixa de ser a mesma.  Recusar e' o que impede um
    // resultado plausivel de sair de um caso nao suportado.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.3, 0.3}, {0.7, 0.7}, 4);
    const Manufaturada &sol = manufaturada("grau2");
    Poisson po(m, paredes(sol, Paredes::TodasDirichlet));
    po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
    bool lancou = false;
    try { po.resolve(); } catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "malha com transicao 2:1 passou sem reclamar");
  }

  t_case("nome_de_solucao_inexistente_e_recusado_com_a_lista");
  {
    bool lancou = false;
    std::string msg;
    try { manufaturada("grau7"); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "nome de solucao inexistente passou sem reclamar");
    T_CHECK_MSG(msg.find("grau2") != std::string::npos,
                "a mensagem nao lista as solucoes validas: \"%s\"", msg.c_str());
  }

  return t_end();
}
