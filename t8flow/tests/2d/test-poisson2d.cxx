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

#include <algorithm>
#include <cmath>
#include <vector>

#include "arranque.hxx"
#include "contorno2d.hxx"
#include "estencil2d.hxx"
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
  double e_inf = 0.0, e_l2 = 0.0, residuo = 0.0, resid_operador = 0.0;
  bool convergiu = false;
  int iteracoes = 0;
};

//! `fracao > 0` refina aleatoriamente: quadtree ADAPTADA, com transicao 2:1.
//
// A semente e' a MESMA em todos os niveis de proposito: assim as malhas da
// sequencia tem o mesmo PADRAO de adaptacao, e a taxa medida e' de refino e nao de
// mudanca de padrao.  Padroes diferentes por nivel dariam uma "ordem" que mistura
// as duas coisas.
Resultado roda(const std::string &nome_sol, int nivel, FechamentoParede fp, Paredes modo,
               double fracao = 0.0, long *transicao = nullptr,
               Esquema esq = Esquema::Fluxo2, int mls_ordem = 3,
               bool fonte_ponto_medio = false, long *rebaixos = nullptr)
{
  const Manufaturada &sol = manufaturada(nome_sol);
  Malha m({{{0.0, 0.0}, {1.0, 1.0}, nivel}}, sc_MPI_COMM_SELF);
  if (fracao > 0.0) m.refina_aleatorio(fracao, 2026, nivel + 1);
  if (transicao != nullptr) *transicao = m.num_folhas_em_transicao();
  Poisson po(m, paredes(sol, modo));
  po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
  po.fechamento(fp);
  po.esquema(esq);
  if (esq == Esquema::MlsGauss || esq == Esquema::MlsPontual) {
    po.mls(mls_ordem, 60);
    if (fonte_ponto_medio) po.quadratura_fonte(Poisson::QuadFonte::PontoMedio);
  }
  if (esq == Esquema::Compacto4)
    po.fonte_laplaciano([&sol](const Ponto2 &p) { return sol.lap_f(p); });
  po.resolve(Motor::Petsc, 1e-14);

  Resultado r;
  r.e_inf = po.erro_inf([&sol](const Ponto2 &p) { return sol.u(p); });
  r.e_l2 = po.erro_l2([&sol](const Ponto2 &p) { return sol.u(p); });
  r.residuo = po.residuo();
  r.convergiu = po.convergiu();
  r.iteracoes = po.iteracoes();
  r.resid_operador = po.residuo_da_exata([&sol](const Ponto2 &p) { return sol.u(p); });
  if (rebaixos != nullptr) *rebaixos = po.rebaixamentos_mls();
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

  // ---- QUADTREE ADAPTADA, com transicao 2:1 ----------------------------------
  //
  // O ORACULO PRINCIPAL AQUI NAO E' O ERRO DA SOLUCAO, e sim o RESIDUO DO
  // OPERADOR: `max |A*u_exata - b| / area`, sem passar pelo solver.  A razao e'
  // que o erro da solucao mistura consistencia, condicionamento e tolerancia do
  // solver, e quando ele e' grande nao se sabe qual dos tres -- foi exatamente o
  // que faltou quando a primeira versao do esquema de transicao deu ordem
  // NEGATIVA.  O residuo do operador isola a montagem.

  t_case("na_malha_adaptada_o_operador_e_exato_para_campo_linear");
  {
    long trans = 0;
    const Resultado r = roda("grau1", 3, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4, &trans);
    T_CHECK_MSG(trans > 0,
                "a malha adaptada nao tem transicao 2:1 -- o caso nao exercita o "
                "esquema que afirma");
    // CAMPO LINEAR E' O PISO: se o operador nao reproduz linear, ele nao e'
    // consistente, e nenhuma ordem de convergencia medida depois significa nada.
    // Medido: SEM a correcao transversal este numero era 53,3 -- erro de O(1), nao
    // de O(h).  Com ela, 4,5e-13.
    T_CHECK_MSG(r.resid_operador < 1e-9,
                "residuo do operador %.4e para campo LINEAR na malha adaptada "
                "(%ld folhas em transicao) -- o esquema de transicao nao e' "
                "consistente", r.resid_operador, trans);
    T_CHECK_MSG(r.e_inf < MAQUINA,
                "erro da solucao %.4e para campo linear na malha adaptada", r.e_inf);
    printf("  (adaptada, grau1: %ld em transicao, residuo do operador %.3e, "
           "erro %.3e)\n", trans, r.resid_operador, r.e_inf);
  }

  t_case("na_malha_adaptada_o_grau2_deixa_de_ser_exato_e_converge_em_segunda_ordem");
  {
    // A EXATIDAO DO GRAU 2 ERA PROPRIEDADE DA MALHA UNIFORME, e perde-se na
    // adaptada -- o que a analise da quadratura ja' previa: na transicao, a
    // distancia entre os centros (3h/4) nao e' a distancia do centro a' face
    // (h/2), e o termo de curvatura nao cancela.  Isso NAO e' defeito; afirmar que
    // deixa de ser exato E' o que impede alguem de "consertar" o que esta' certo.
    const Resultado a = roda("grau2", 3, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4);
    const Resultado b = roda("grau2", 5, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4);
    T_CHECK_MSG(a.e_inf > 1e-5,
                "na malha adaptada o grau2 deu erro %.3e -- se fosse de maquina, a "
                "distincao entre malha uniforme e adaptada nao estaria sendo "
                "exercitada", a.e_inf);
    const double ordem = std::log2(a.e_inf / b.e_inf) / 2.0;
    T_CHECK_MSG(ordem > 1.7 && ordem < 2.3,
                "ordem %.2f para grau2 na malha adaptada (esperava ~2)", ordem);
    printf("  (adaptada, grau2: erro %.3e -> %.3e, ordem %.2f)\n",
           a.e_inf, b.e_inf, ordem);
  }

  t_case("na_malha_adaptada_a_trigonometrica_converge_em_segunda_ordem");
  {
    const Resultado a = roda("trig", 4, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4);
    const Resultado b = roda("trig", 6, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4);
    const double oi = std::log2(a.e_inf / b.e_inf) / 2.0;
    const double ol = std::log2(a.e_l2 / b.e_l2) / 2.0;
    T_CHECK_MSG(oi > 1.6 && oi < 2.3, "ordem na norma do maximo: %.2f", oi);
    T_CHECK_MSG(ol > 1.7 && ol < 2.3, "ordem na norma 2: %.2f", ol);
    printf("  (adaptada, trig: ordem inf %.2f, L2 %.2f)\n", oi, ol);
  }

  t_case("na_malha_adaptada_com_neumann_o_linear_continua_exato");
  {
    // A transicao 2:1 e o fechamento de Neumann se encontram aqui: a celula de
    // parede pode ter o vizinho interior mais fino, e o fechamento de tres pontos
    // tem de usar a MEDIA deles.  A primeira versao exigia vizinho UNICO e lancou
    // na primeira malha aleatoria -- a guarda funcionando.
    long trans = 0;
    const Resultado r = roda("grau1", 4, FechamentoParede::TresPontos,
                             Paredes::TresNeumann, 0.4, &trans);
    T_CHECK_MSG(trans > 0, "sem transicao: o caso nao exercita o que afirma");
    T_CHECK_MSG(r.resid_operador < 1e-9,
                "residuo do operador %.4e para linear na malha adaptada com tres "
                "paredes de Neumann", r.resid_operador);
    printf("  (adaptada + 3 Neumann, grau1: residuo do operador %.3e)\n",
           r.resid_operador);
  }

  t_case("o_residuo_do_operador_distingue_malha_uniforme_de_adaptada_no_grau2");
  {
    // O par que da' sentido aos dois casos de cima: na uniforme o operador e'
    // exato para grau 2, na adaptada nao.  Se os dois dessem o mesmo numero, o
    // esquema de transicao nao estaria sendo exercitado em nenhum deles.
    const Resultado u = roda("grau2", 3, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.0);
    const Resultado a = roda("grau2", 3, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.4);
    T_CHECK_MSG(u.resid_operador < 1e-9,
                "na malha UNIFORME o residuo do operador para grau2 foi %.4e -- "
                "esperava maquina", u.resid_operador);
    T_CONTRASTE(u.resid_operador > 1e-16 ? u.resid_operador : 1e-16,
                a.resid_operador, 1e6,
                "residuo do operador, uniforme contra adaptada, grau 2");
    printf("  (residuo do operador, grau2: uniforme %.3e, adaptada %.3e)\n",
           u.resid_operador, a.resid_operador);
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

  // ---- O ESQUEMA COMPACTO DE QUARTA ORDEM ------------------------------------
  //
  // Duas variantes do MESMO esquema de nove pontos, com lados direitos diferentes:
  // Mehrstellen usa `f + (h^2/12) lap(f)` e exige `lap(f)` ANALITICO; Lele e'
  // implicito na segunda derivada e usa `f` nos nove pontos vizinhos.  A algebra
  // diz que sao equivalentes a O(h^4); os casos abaixo VERIFICAM em vez de confiar.

  t_case("o_compacto_e_exato_ate_grau_3_e_o_limite_e_o_fechamento");
  {
    // O INTERIOR do esquema e' exato ate' grau 5 (o erro e' proporcional a
    // derivadas de SEXTA ordem).  Mas o FECHAMENTO de parede usa quatro pontos --
    // parede e tres celulas --, que e' exato ate' grau 3.  Logo o corte esta'
    // entre grau 3 e grau 4, e nao entre 5 e 6.
    //
    // MEDIDO, e e' por isso que a clausula afirma o LIMITE e nao so' a exatidao:
    //   grau1 1,44e-14   grau2 1,12e-13   grau3 3,58e-13   grau4 3,99e-05
    for (const char *nome : {"grau1", "grau2", "grau3"}) {
      const Resultado r = roda(nome, 4, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::Compacto4Lele);
      T_CHECK_MSG(r.convergiu, "%s: o solver nao convergiu no compacto", nome);
      T_CHECK_MSG(r.e_inf < MAQUINA, "%s no compacto: erro %.4e (limite %.0e)",
                  nome, r.e_inf, MAQUINA);
    }
    // E O GRAU 4 NAO: se ele tambem fosse exato, a clausula estaria afirmando
    // menos do que parece -- nao se saberia onde o esquema deixa de ser exato.
    const Resultado r4 = roda("grau4", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr,
                              Esquema::Compacto4Lele);
    T_CHECK_MSG(r4.e_inf > 1e-7,
                "o grau 4 deu erro %.4e no compacto -- se fosse de maquina, o "
                "limite do fechamento nao estaria sendo exercitado", r4.e_inf);
    printf("  (compacto, nivel 4: grau3 exato, grau4 %.3e)\n", r4.e_inf);
  }

  t_case("o_compacto_converge_em_QUARTA_ordem");
  {
    const Resultado a4 = roda("trig", 3, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr,
                              Esquema::Compacto4Lele);
    const Resultado b4 = roda("trig", 5, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr,
                              Esquema::Compacto4Lele);
    const double oi = std::log2(a4.e_inf / b4.e_inf) / 2.0;
    const double ol = std::log2(a4.e_l2 / b4.e_l2) / 2.0;
    T_CHECK_MSG(oi > 3.6 && oi < 4.4, "ordem na norma do maximo: %.2f (esperava ~4)", oi);
    T_CHECK_MSG(ol > 3.6 && ol < 4.4, "ordem na norma 2: %.2f (esperava ~4)", ol);
    printf("  (compacto, trig: ordem inf %.2f, L2 %.2f; erro %.3e -> %.3e)\n",
           oi, ol, a4.e_inf, b4.e_inf);
  }

  t_case("o_compacto_e_muito_mais_preciso_que_o_de_segunda_ordem");
  {
    // O par que mostra que a quarta ordem nao e' so' uma taxa: no MESMO nivel, o
    // erro tem de ser ordens de grandeza menor.  Sem este caso, "ordem 4" poderia
    // conviver com um erro absoluto pior.
    const Resultado dois = roda("trig", 5, FechamentoParede::TresPontos,
                                Paredes::TodasDirichlet);
    const Resultado quatro = roda("trig", 5, FechamentoParede::TresPontos,
                                  Paredes::TodasDirichlet, 0.0, nullptr,
                                  Esquema::Compacto4Lele);
    T_CONTRASTE(quatro.e_inf > 1e-16 ? quatro.e_inf : 1e-16, dois.e_inf, 20.0,
                "erro do compacto contra o de segunda ordem, mesmo nivel");
    printf("  (nivel 5, trig: fluxo2 %.3e, compacto %.3e -- %.0fx)\n",
           dois.e_inf, quatro.e_inf, dois.e_inf / quatro.e_inf);
  }

  t_case("mehrstellen_e_lele_concordam_na_ordem_do_truncamento");
  {
    // A VERIFICACAO DA ALGEBRA, que e' o motivo de ter as duas variantes.
    // `A_x f / (6/5)` expande em `f + (h^2/12) f_xx + O(h^4)`, porque
    // `(1/10)/(6/5) = 1/12`.  Logo os dois lados direitos sao o mesmo a O(h^4) --
    // e os dois esquemas tem de dar resultados que diferem na ORDEM DO ERRO, nao
    // mais que isso.
    //
    // Para polinomio de grau <= 4 eles dao valores IDENTICOS, porque ali o termo
    // de O(h^4) do lado direito envolve derivadas de quarta ordem de `f`, que sao
    // nulas.  Para a trigonometrica diferem em ~1%, que e' a ordem do proprio
    // truncamento.
    const Resultado m = roda("grau4", 4, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.0, nullptr, Esquema::Compacto4);
    const Resultado l = roda("grau4", 4, FechamentoParede::TresPontos,
                             Paredes::TodasDirichlet, 0.0, nullptr, Esquema::Compacto4Lele);
    T_NEAR(m.e_inf, l.e_inf, 1e-12,
           "grau4: Mehrstellen e Lele tem de coincidir (lap f e' constante)");

    const Resultado mt = roda("trig", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr, Esquema::Compacto4);
    const Resultado lt = roda("trig", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr, Esquema::Compacto4Lele);
    const double dif = std::fabs(mt.e_inf - lt.e_inf) / mt.e_inf;
    T_CHECK_MSG(dif < 0.25,
                "trig: Mehrstellen %.4e e Lele %.4e diferem %.1f%% -- mais que a "
                "ordem do truncamento, logo nao sao o mesmo esquema",
                mt.e_inf, lt.e_inf, 100.0 * dif);
    printf("  (grau4 identicos; trig: Mehrstellen %.3e, Lele %.3e, dif %.1f%%)\n",
           mt.e_inf, lt.e_inf, 100.0 * dif);
  }

  t_case("o_compacto_recusa_o_que_nao_sabe");
  {
    const Manufaturada &sol = manufaturada("grau2");

    // MALHA ADAPTADA: o estencil usa vizinhos DIAGONAIS, que na transicao 2:1 nao
    // estao definidos sem reconstrucao propria.
    {
      Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
      m.refina_aleatorio(0.4, 1, 4);
      Poisson po(m, paredes(sol, Paredes::TodasDirichlet));
      po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
      po.esquema(Esquema::Compacto4Lele);
      bool lancou = false;
      try { po.resolve(); } catch (const std::runtime_error &) { lancou = true; }
      T_CHECK_MSG(lancou, "o compacto aceitou malha adaptada");
    }

    // MEHRSTELLEN SEM `lap(f)`: rodar com zero daria SEGUNDA ordem e pareceria
    // quarta falhando -- o pior resultado possivel, porque o numero sai plausivel.
    {
      Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
      Poisson po(m, paredes(sol, Paredes::TodasDirichlet));
      po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
      po.esquema(Esquema::Compacto4);          // sem fonte_laplaciano
      bool lancou = false;
      try { po.resolve(); } catch (const std::invalid_argument &) { lancou = true; }
      T_CHECK_MSG(lancou, "Mehrstellen rodou sem lap(f)");
    }

    // NEUMANN: o fechamento do compacto e' de Dirichlet; usa-lo para Neumann daria
    // quarta ordem no interior e segunda na parede, sem aviso.
    {
      Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
      Poisson po(m, paredes(sol, Paredes::TresNeumann));
      po.fonte([&sol](const Ponto2 &p) { return sol.f(p); });
      po.esquema(Esquema::Compacto4Lele);
      bool lancou = false;
      try { po.resolve(); } catch (const std::runtime_error &) { lancou = true; }
      T_CHECK_MSG(lancou, "o compacto aceitou Neumann");
    }
  }

  // -------------------------------------------------------------------------
  // D4 .. D7 -- derivada MLS na faceta e quadratura de Gauss
  // -------------------------------------------------------------------------

  t_case("o_mlsgauss_e_exato_no_grau_tres_onde_o_fluxo2_nao_e");
  {
    // D4.  O par de medidas na MESMA malha, mudando SO' o esquema.  O grau 3 e' o
    // grau onde a distincao aparece: o fluxo de segunda ordem ja' erra nele, e o
    // `MlsGauss` de ordem 3 e' exato.
    long reb = 0;
    const Resultado mg = roda("grau3", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr,
                              Esquema::MlsGauss, 3, false, &reb);
    const Resultado f2 = roda("grau3", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet);
    T_CHECK_MSG(mg.e_inf < MAQUINA,
                "MlsGauss no grau 3: erro %.3e (limite %.0e)", mg.e_inf, MAQUINA);
    T_CHECK_MSG(reb == 0, "%ld ajustes MLS rebaixaram a ordem em malha UNIFORME -- "
                "aqui o suporte e' regular e nao deveria haver nenhum", reb);
    T_CONTRASTE(mg.e_inf > 1e-16 ? mg.e_inf : 1e-16, f2.e_inf, 1e4,
                "grau 3: MlsGauss contra o fluxo de segunda ordem");
    printf("  (grau3 uniforme: mlsgauss %.3e, fluxo2 %.3e)\n", mg.e_inf, f2.e_inf);
  }

  t_case("na_malha_adaptada_o_mlsgauss_recupera_a_exatidao_do_grau2");
  {
    // D5.  A RAZAO DE EXISTIR DO ESQUEMA, medida contra a clausula T2: la' o grau 2
    // DEIXA de ser exato em malha adaptada; aqui ele volta a ser, na mesma malha --
    // mesma fracao, mesma semente -- e o grau 3 tambem.
    long reb2 = 0, reb3 = 0;
    const Resultado g2 = roda("grau2", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.4, nullptr,
                              Esquema::MlsGauss, 3, false, &reb2);
    const Resultado g3 = roda("grau3", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.4, nullptr,
                              Esquema::MlsGauss, 3, false, &reb3);
    const Resultado f2 = roda("grau2", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.4);

    T_CHECK_MSG(g2.e_inf < MAQUINA,
                "grau2 em malha ADAPTADA com MlsGauss: erro %.3e (limite %.0e)",
                g2.e_inf, MAQUINA);
    T_CHECK_MSG(g3.e_inf < MAQUINA,
                "grau3 em malha ADAPTADA com MlsGauss: erro %.3e (limite %.0e)",
                g3.e_inf, MAQUINA);
    T_CHECK_MSG(reb2 == 0 && reb3 == 0,
                "houve rebaixamento de ordem (%ld e %ld) -- a exatidao afirmada "
                "acima nao seria de ordem 3 em todo o dominio", reb2, reb3);
    // o CONTRASTE e' contra o esquema de fluxo na mesma malha: sem ele, "exato em
    // malha adaptada" nao distingue o esquema de uma malha que virou uniforme
    T_CHECK_MSG(f2.e_inf > 1e-5,
                "o fluxo2 deu %.3e no grau2 adaptado -- se fosse de maquina, a "
                "malha nao estaria de fato adaptada e o caso nao mediria nada",
                f2.e_inf);
    printf("  (grau2 adaptada: mlsgauss %.3e, fluxo2 %.3e; grau3 mlsgauss %.3e)\n",
           g2.e_inf, f2.e_inf, g3.e_inf);
  }

  t_case("a_fonte_no_ponto_medio_trava_o_esquema_em_segunda_ordem");
  {
    // D6.  O interruptor existe para esta medida.  Grau 4 com ajuste de ordem 4: o
    // fluxo da faceta e' EXATO nos dois casos, e a unica diferenca e' a quadratura
    // de `integral f dA`.  Com o ponto medio o esquema cai para segunda ordem.
    long reb = 0;
    const Resultado g4 = roda("grau4", 4, FechamentoParede::TresPontos,
                              Paredes::TodasDirichlet, 0.0, nullptr,
                              Esquema::MlsGauss, 4, false, &reb);
    const Resultado pm4 = roda("grau4", 4, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsGauss, 4, true);
    const Resultado pm5 = roda("grau4", 5, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsGauss, 4, true);
    T_CHECK_MSG(g4.e_inf < MAQUINA,
                "grau4 com ordem 4 e Gauss 2x2: erro %.3e (limite %.0e)",
                g4.e_inf, MAQUINA);
    T_CHECK_MSG(reb == 0, "%ld rebaixamentos em malha uniforme com ordem 4", reb);
    T_CONTRASTE(g4.e_inf > 1e-16 ? g4.e_inf : 1e-16, pm4.e_inf, 1e4,
                "grau4: Gauss 2x2 contra ponto medio na fonte");
    const double ordem = std::log2(pm4.e_inf / pm5.e_inf);
    T_CHECK_MSG(ordem > 1.7 && ordem < 2.3,
                "com a fonte no ponto medio a ordem deu %.2f (esperava ~2: e' o "
                "travamento que o interruptor existe para mostrar)", ordem);
    printf("  (grau4: gauss2x2 %.3e, ponto medio %.3e -> %.3e, ordem %.2f)\n",
           g4.e_inf, pm4.e_inf, pm5.e_inf, ordem);
  }

  t_case("o_fluxo_da_faceta_de_transicao_e_o_mesmo_dos_dois_lados");
  {
    // D7.  CONSERVACAO, afirmada pelo caminho em que ela pode falhar: a celula
    // GROSSA parte a faceta em duas sub-facetas e a FINA nao parte a sua.  Se os
    // pontos de quadratura dos dois lados nao coincidirem, os fluxos nao se
    // cancelam -- e um esquema nao conservativo AINDA CONVERGE, logo nenhum dos
    // casos de convergencia acima pegaria isso.
    //
    // O oraculo refaz a geometria por conta propria, de proposito: conferir o
    // codigo com o proprio codigo nao discrimina.
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.0, 0.0}, {0.5, 0.5}, 4);

    ConstrutorEstencil ce(m, 3, 60, 6);
    ce.dirichlet([](const Ponto2 &) { return 0.0; });

    const double no = 1.0 / std::sqrt(3.0);
    int pares = 0;
    double pior_ponto = 0.0, pior_peso = 0.0;

    for (long b = 0; b < m.num_blocos() && pares < 8; b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b) && pares < 8; i++) {
        const int nself = m.nivel(b, i);
        const Ponto2 c = m.centroide(b, i);
        const Ponto2 h = m.lados_da_folha(b, i);

        for (int fa = 0; fa < 4 && pares < 8; fa++) {
          const std::vector<int> nvz = m.niveis_vizinhos(b, i, fa);
          if (nvz.empty() || nvz[0] <= nself) continue;   // so' a GROSSA da transicao

          const int d = fa / 2, t = 1 - d;
          const double sinal = (fa % 2 == 0) ? -1.0 : 1.0;
          const double larg = h[t] / 2.0;                 // sub-faceta = face da fina

          for (int sub = 0; sub < 2; sub++) {
            Ponto2 cf = c;
            cf[d] += sinal * h[d] / 2.0;
            cf[t] += -h[t] / 2.0 + (sub + 0.5) * larg;

            // o CENTRO DA FACETA DA FINA, achado pelo outro caminho: a celula fina
            // e' a que contem o ponto logo DEPOIS da faceta
            Ponto2 dentro = cf;
            dentro[d] += sinal * larg / 4.0;
            const long id_fina = m.id_do_ponto(dentro);
            T_CHECK_MSG(id_fina >= 0, "nao achei a celula fina do outro lado");

            // reconstroi o centro da faceta a partir da celula fina
            long bf = 0, jf = 0;
            bool achou = false;
            for (long bb = 0; bb < m.num_blocos() && !achou; bb++)
              for (long jj = 0; jj < m.num_folhas_do_bloco(bb); jj++)
                if (m.id_global(bb, jj) == id_fina) { bf = bb; jf = jj; achou = true; break; }
            T_CHECK_MSG(achou, "id %ld da celula fina nao mapeou de volta", id_fina);

            const Ponto2 cfina = m.centroide(bf, jf);
            const Ponto2 hfina = m.lados_da_folha(bf, jf);
            Ponto2 cf_fina = cfina;
            cf_fina[d] -= sinal * hfina[d] / 2.0;   // a faceta OPOSTA, vista da fina

            pior_ponto = std::max(pior_ponto,
                                  std::max(std::fabs(cf[0] - cf_fina[0]),
                                           std::fabs(cf[1] - cf_fina[1])));
            T_NEAR(hfina[t], larg, 1e-14, "largura da faceta vista dos dois lados");

            // E OS ESTENCEIS NOS PONTOS DE GAUSS: identicos, porque o suporte
            // depende do PONTO.  Comparados peso a peso, nao por uma soma -- somas
            // iguais com pesos diferentes seriam um falso verde.
            for (int q = 0; q < 2; q++) {
              Ponto2 xg = cf, yg = cf_fina;
              xg[t] += (q == 0 ? -no : no) * larg / 2.0;
              yg[t] += (q == 0 ? -no : no) * hfina[t] / 2.0;
              const Estencil ex = ce.derivada_em(xg, d);
              const Estencil ey = ce.derivada_em(yg, d);
              T_CHECK_MSG(ex.ids == ey.ids,
                          "os dois lados usaram SUPORTES diferentes no mesmo ponto "
                          "de quadratura (%zu contra %zu pontos)",
                          ex.ids.size(), ey.ids.size());
              for (size_t k = 0; k < ex.pesos.size(); k++)
                pior_peso = std::max(pior_peso, std::fabs(ex.pesos[k] - ey.pesos[k]));
            }
            pares++;
          }
        }
      }

    T_CHECK_MSG(pares >= 4,
                "so' %d sub-facetas de transicao examinadas -- a malha nao tem a "
                "transicao 2:1 que o caso precisa exercitar", pares);
    T_NEAR(pior_ponto, 0.0, 1e-14, "centro da sub-faceta visto dos dois lados");
    T_NEAR(pior_peso, 0.0, 1e-12, "pesos do estencil de derivada vistos dos dois lados");
    printf("  (%d sub-facetas de transicao: ponto %.1e, peso %.1e)\n",
           pares, pior_ponto, pior_peso);
  }

  // -------------------------------------------------------------------------
  // W5 -- as duas familias, com tudo o mais igual
  // -------------------------------------------------------------------------

  t_case("as_duas_familias_diferem_com_tudo_o_mais_igual");
  {
    // MESMO AJUSTE, MESMA ORDEM, MESMO SUPORTE, MESMO FECHAMENTO DE CONTORNO.  A
    // unica diferenca entre as duas corridas e' integrar a equacao contra
    // impo-la no ponto -- e e' isso que o caso mede.
    //
    // Afirmar que DIFEREM e' o que impede alguem de tratar a escolha da familia
    // como detalhe de implementacao.  E o grau 3 e' onde as duas concordam
    // (ambas exatas), logo nao serve: a solucao trigonometrica serve.
    long r1 = 0, r2 = 0;
    const Resultado vf4 = roda("trig", 4, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsGauss, 3, false, &r1);
    const Resultado vf5 = roda("trig", 5, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsGauss, 3, false, &r1);
    const Resultado df4 = roda("trig", 4, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsPontual, 3, false, &r2);
    const Resultado df5 = roda("trig", 5, FechamentoParede::TresPontos,
                               Paredes::TodasDirichlet, 0.0, nullptr,
                               Esquema::MlsPontual, 3, false, &r2);

    const double ovf = std::log2(vf4.e_inf / vf5.e_inf);
    const double odf = std::log2(df4.e_inf / df5.e_inf);

    // VOLUMES FINITOS chega perto de 3 com ajuste de ordem 3; DIFERENCAS FINITAS
    // fica em 2, porque o truncamento do laplaciano pontual exato ate' grau 3 e'
    // `O(h^2)` -- vem da quarta derivada.
    T_CHECK_MSG(ovf > 2.5, "volumes finitos deu ordem %.2f, esperada ~2,8", ovf);
    T_CHECK_MSG(odf > 1.7 && odf < 2.4,
                "diferencas finitas deu ordem %.2f, esperada ~2", odf);
    // A ORDEM DOS ARGUMENTOS E' (melhor, pior), e inverte-la foi o primeiro
    // resultado deste caso: o contraste saiu 0,023x e o arreio reprovou dizendo
    // "o oraculo nao discrimina", que e' exatamente o trabalho dele.
    T_CONTRASTE(vf5.e_inf, df5.e_inf, 5.0,
                "erro das duas familias no nivel 5 com o MESMO ajuste de ordem 3");
    printf("  (W5: mesma ordem de ajuste; VF ordem %.2f erro %.3e, "
           "DF ordem %.2f erro %.3e)\n", ovf, vf5.e_inf, odf, df5.e_inf);
  }

  return t_end();
}
