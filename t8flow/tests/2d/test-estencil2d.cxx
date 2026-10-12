// O ESTENCIL COM A ORDEM DESEJADA, e a clausula F3.
//
// ORACULO PRINCIPAL -- reproducao polinomial INCLUSIVE PERTO DO CONTORNO.
// Preenchendo o campo das celulas E os valores de contorno com o MESMO polinomio
// analitico, o valor interpolado tem de bater em qualquer ponto do dominio,
// inclusive junto a' parede, onde parte do suporte vem do fechamento.  E' o
// oraculo da suite da HiGTree, e e' o unico que exercita o `rhs`.
//
// E O CAMPO DE CONTORNO VARIA EM TODAS AS DIRECOES de proposito.  Com valor de
// contorno CONSTANTE -- parede de nao-deslizamento, que e' quase todo exemplo --
// fechar contra a parede ERRADA da' o mesmo numero, e foi assim que nove defeitos
// de fechamento sobreviveram cinco anos na HiGTree.
//
// A CLAUSULA F3, que a reproducao polinomial NAO ve': "o suporte cai na grade e e'
// local".  Medido na suite deste projeto: deslocando todos os centroides meia
// celula, a reproducao continua exata -- ela verifica a DISCRETIZACAO, nao a
// malha.  Entao o suporte e' afirmado por GEOMETRIA: todo ponto e' centroide de
// uma celula existente (ou centro de faceta de contorno), e o diametro e' limitado
// -- nao a malha inteira.
#include "../testing.h"

#include <set>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "estencil2d.hxx"
#include "lote2d.hxx"
#include "malha2d.hxx"

using t8flow2d::Bloco;
using t8flow2d::ConstrutorEstencil;
using t8flow2d::Estencil;
using t8flow2d::LoteEstenceis;
using t8flow2d::Malha;
using t8flow2d::Mls;
using t8flow2d::Ponto2;

//! Campo de grau 1 que varia nas DUAS direcoes, com coeficientes distintos.
static double grau1(const Ponto2 &p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1]; }
static double grau2(const Ponto2 &p) { return grau1(p) + 4.0 * p[0] * p[0] - 2.5 * p[1] * p[1]; }

//! Valor do estencil aplicado a um campo: pesos nas celulas mais o `rhs`.
static double aplica(const Estencil &e, const Malha &m,
                     double (*f)(const Ponto2 &))
{
  double v = e.rhs;
  for (size_t k = 0; k < e.ids.size(); k++) {
    // o valor da celula vem da GEOMETRIA dela, pelo id -- e e' por isso que a
    // malha precisa saber achar o centroide de um id
    bool achou = false;
    for (long b = 0; b < m.num_blocos() && !achou; b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        if (m.id_global(b, i) == e.ids[k]) {
          v += e.pesos[k] * f(m.centroide(b, i));
          achou = true;
          break;
        }
  }
  return v;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Um bloco uniforme basta para o estencil, e nivel 3 da' 64 celulas -- malha
  // pequena o suficiente para a varredura do teste e grande o suficiente para
  // haver interior de verdade.
  const std::vector<Bloco> blocos = {{{0.0, 0.0}, {1.0, 1.0}, 3}};
  Malha m(blocos, sc_MPI_COMM_SELF);

  const int ORDEM = 1;
  ConstrutorEstencil ce(m, ORDEM, 4 * Mls::min_pontos(ORDEM));

  t_case("sem_condicao_de_contorno_o_suporte_perto_da_parede_e_recusado");
  {
    // Sem Dirichlet instalado, o contorno nao contribui.  Perto da parede o
    // suporte fica curto, e a chamada LANCA -- em vez de devolver um ajuste de
    // ordem menor que a pedida, que daria numeros plausiveis e errados.
    //
    // O caso usa max_aneis=1 para garantir que o suporte nao se salve crescendo
    // para dentro: com aneis suficientes, a celula do canto alcanca o interior.
    ConstrutorEstencil curto(m, 2, 40, 1);
    bool lancou = false;
    try { curto.em({0.02, 0.02}); }
    catch (const std::runtime_error &) { lancou = true; }
    T_CHECK_MSG(lancou, "suporte curto no canto passou sem reclamar");
  }

  // A CONDICAO DE DIRICHLET, com o MESMO campo do interior: e' o que faz a
  // reproducao valer ate' a parede.
  ce.dirichlet([](const Ponto2 &p) { return grau1(p); });

  t_case("o_estencil_reproduz_o_grau_da_ordem_no_interior");
  double pior_int = 0.0;
  {
    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const Ponto2 c = m.centroide(0, i);
      if (c[0] < 0.3 || c[0] > 0.7 || c[1] < 0.3 || c[1] > 0.7) continue;  // interior
      const Estencil e = ce.em(c);
      const double d = std::fabs(aplica(e, m, grau1) - grau1(c));
      if (d > pior_int) pior_int = d;
    }
    T_NEAR(pior_int, 0.0, 1e-10, "pior erro no interior, campo de grau 1");
  }

  // O CASO QUE EXERCITA O `rhs`: pontos JUNTO A' PAREDE, onde parte do suporte vem
  // do fechamento.  Se o `rhs` estiver errado -- ou se o contorno entrar como
  // celula -- aqui aparece, e so' aqui.
  t_case("o_estencil_reproduz_junto_a_parede_pelo_fechamento");
  double pior_par = 0.0;
  long com_rhs = 0;
  {
    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const Ponto2 c = m.centroide(0, i);
      const bool junto = (c[0] < 0.1 || c[0] > 0.9 || c[1] < 0.1 || c[1] > 0.9);
      if (!junto) continue;
      const Estencil e = ce.em(c);
      if (e.rhs != 0.0) com_rhs++;
      const double d = std::fabs(aplica(e, m, grau1) - grau1(c));
      if (d > pior_par) pior_par = d;
    }
    T_CHECK_MSG(com_rhs > 0,
                "nenhum estencil junto a' parede teve rhs -- o fechamento nao foi "
                "exercitado, e o caso nao afirma nada");
    T_NEAR(pior_par, 0.0, 1e-10, "pior erro junto a' parede, campo de grau 1");
  }

  t_case("o_estencil_erra_mais_no_grau_acima_da_ordem");
  {
    // CONTRASTE, com a condicao de contorno TAMBEM no campo de grau 2 -- senao o
    // erro mediria a incoerencia entre interior e parede, nao a ordem do ajuste.
    ConstrutorEstencil ce2(m, ORDEM, 4 * Mls::min_pontos(ORDEM));
    ce2.dirichlet([](const Ponto2 &p) { return grau2(p); });
    double pior = 0.0;
    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const Ponto2 c = m.centroide(0, i);
      if (c[0] < 0.3 || c[0] > 0.7 || c[1] < 0.3 || c[1] > 0.7) continue;
      const Estencil e = ce2.em(c);
      const double d = std::fabs(aplica(e, m, grau2) - grau2(c));
      if (d > pior) pior = d;
    }
    T_CONTRASTE(pior_int > 1e-16 ? pior_int : 1e-16, pior, 1e3,
                "grau 1 contra grau 2 no estencil de ordem 1");
  }

  // ---- a clausula F3: o suporte cai na grade e e' local ---------------------

  t_case("F3_todo_ponto_do_suporte_esta_na_grade");
  {
    // "Na grade" = cada ponto e' centroide de uma celula que EXISTE (id >= 0) ou
    // centro de uma faceta de contorno (id == -1).  Deslocar o suporte meia celula
    // -- o defeito que a reproducao polinomial NAO ve' -- quebra esta afirmacao.
    long conferidos = 0, fora = 0;
    for (long i = 0; i < m.num_folhas_do_bloco(0); i += 7) {
      const ConstrutorEstencil::Suporte s = ce.suporte_em(m.centroide(0, i));
      for (size_t k = 0; k < s.pontos.size(); k++) {
        conferidos++;
        if (s.ids[k] >= 0) {
          // ponto de celula: o id que a GEOMETRIA daquele ponto devolve tem de ser
          // o mesmo que o suporte declara
          if (m.id_do_ponto(s.pontos[k]) != s.ids[k]) fora++;
        } else {
          // ponto de contorno: tem de estar SOBRE a borda do dominio
          const double eps = 1e-12;
          const bool na_borda = std::fabs(s.pontos[k][0] - 0.0) < eps ||
                                std::fabs(s.pontos[k][0] - 1.0) < eps ||
                                std::fabs(s.pontos[k][1] - 0.0) < eps ||
                                std::fabs(s.pontos[k][1] - 1.0) < eps;
          if (!na_borda) fora++;
        }
      }
    }
    T_CHECK_MSG(conferidos > 0, "nenhum ponto conferido -- o caso nao mediu nada");
    T_CHECK_MSG(fora == 0, "%ld de %ld pontos do suporte fora da grade",
                fora, conferidos);
  }

  t_case("F3_o_suporte_e_local_e_nao_a_malha_inteira");
  {
    // "Local" = o diametro do suporte e' da ordem de poucas celulas, nao do
    // dominio.  Busca em largura SEM parada devolveria a malha toda (medido na
    // HiGTree: 64 pontos numa malha de 64), e o ajuste reproduziria polinomio
    // mesmo assim -- o verde nao acusaria.
    const double h = 1.0 / 8.0;          // lado da celula no nivel 3
    long pior_pontos = 0;
    double pior_diam = 0.0;
    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const ConstrutorEstencil::Suporte s = ce.suporte_em(m.centroide(0, i));
      pior_pontos = std::max<long>(pior_pontos, static_cast<long>(s.pontos.size()));
      pior_diam = std::fmax(pior_diam, s.diametro);
    }
    T_CHECK_MSG(pior_pontos < m.num_folhas(),
                "o maior suporte tem %ld pontos para uma malha de %ld celulas -- "
                "isso e' a malha inteira, nao um suporte local",
                pior_pontos, m.num_folhas());
    T_CHECK_MSG(pior_diam <= 6.0 * h,
                "o maior diametro de suporte e' %.4f, mais de seis celulas (%.4f) "
                "-- o suporte deixou de ser local", pior_diam, 6.0 * h);
    printf("  (maior suporte: %ld pontos, diametro %.4f = %.1f celulas)\n",
           pior_pontos, pior_diam, pior_diam / h);
  }

  t_case("o_anel_e_fechado_e_o_suporte_nao_depende_da_ordem_de_percurso");
  {
    // Duas construcoes do MESMO estencil tem de dar o MESMO conjunto de ids.  Se o
    // crescimento parasse no meio de um anel, o suporte dependeria da ordem de
    // visita -- e' como a clausula de desempate da HiGTree passava por sorte.
    const Ponto2 x = {0.41, 0.56};
    const ConstrutorEstencil::Suporte a = ce.suporte_em(x);
    const ConstrutorEstencil::Suporte b = ce.suporte_em(x);
    std::set<long> ia(a.ids.begin(), a.ids.end()), ib(b.ids.begin(), b.ids.end());
    T_CHECK_MSG(ia == ib, "duas construcoes do mesmo suporte deram conjuntos diferentes");
    // E o anel tem de ter fechado: com ordem 1 (minimo 3) e vizinhanca de face, um
    // anel da' 5 pontos (a celula e os quatro vizinhos) -- nunca 3 ou 4.
    T_CHECK_MSG(a.pontos.size() >= 5,
                "suporte de %zu pontos: menos que o anel fechado de 5, logo parou no "
                "meio do anel", a.pontos.size());
  }

  // ---- o lote, com o fechamento preservado ---------------------------------

  t_case("o_lote_de_estenceis_da_o_mesmo_que_um_por_um");
  {
    std::vector<Ponto2> alvos;
    for (long i = 0; i < m.num_folhas_do_bloco(0); i += 3)
      alvos.push_back(m.centroide(0, i));

    LoteEstenceis lote(ORDEM, 4 * Mls::min_pontos(ORDEM));
    lote.usa_back_end(t8flow2d::BackEnd::Cpu);
    std::vector<Estencil> em_lote;
    ce.em_lote(alvos, lote, em_lote);

    T_CHECK_MSG(em_lote.size() == alvos.size(), "%zu estenceis para %zu alvos",
                em_lote.size(), alvos.size());
    double pior = 0.0;
    for (size_t q = 0; q < alvos.size(); q++) {
      const Estencil um = ce.em(alvos[q]);
      // compara pelo VALOR aplicado, porque a ordem dos ids pode diferir
      const double a = aplica(um, m, grau1);
      const double b = aplica(em_lote[q], m, grau1);
      pior = std::fmax(pior, std::fabs(a - b));
      // e o rhs tem de estar presente nos dois onde ha' contorno
      if ((um.rhs != 0.0) != (em_lote[q].rhs != 0.0)) {
        T_CHECK_MSG(false, "alvo %zu: rhs presente em um caminho e ausente no outro", q);
        break;
      }
    }
    T_NEAR(pior, 0.0, 1e-10, "pior diferenca entre o lote e o um-por-um");
  }

  t_case("ponto_fora_do_dominio_e_recusado");
  {
    bool lancou = false;
    try { ce.em({5.0, 5.0}); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "ponto fora do dominio passou sem reclamar");
  }

  return t_end();
}
