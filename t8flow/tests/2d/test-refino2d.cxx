// REFINO LOCAL DENTRO DE UM BLOCO, e o que ele tem de preservar.
//
// ORACULOS, e a escolha de cada um:
//
//   contagem por regiao   diz que refinou ONDE se pediu, e so' ali;
//   AREA TOTAL            diz que nao se perdeu nem se duplicou celula.  Contagem
//                         nao pega isso: uma celula com geometria errada entra na
//                         contagem do mesmo jeito;
//   RAZAO 2:1             e' o oraculo que o `t8_forest_set_balance` existe para
//                         satisfazer, e sem ele a malha COMMITA, as contagens
//                         fecham, e o salto de nivel atravessa o estencil do
//                         discretizador sem nenhum aviso.  Nenhum dos outros dois
//                         o ve;
//   juncoes sobrevivem    refino perto da interface entre blocos nao pode
//                         transformar face interna em contorno -- que e' o modo de
//                         falha mais caro desta frente.
//
// A razao 2:1 e' afirmada sobre TODAS as folhas e TODAS as faces, com os niveis
// vindos dos elementos vizinhos (nao dos indices: o indice que
// `t8_forest_leaf_face_neighbors` devolve e' da floresta, nao da arvore, e
// confundi-los entrega um elemento existente, so' que outro).
#include "../testing.h"

#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "malha2d.hxx"
#include "mls2d.hxx"

using t8flow2d::Bloco;
using t8flow2d::Malha;
using t8flow2d::Mls;
using t8flow2d::Ponto2;

static double grau2(const Ponto2 &p)
{
  return 1.0 + 2.0 * p[0] - 3.0 * p[1] + 0.5 * p[0] * p[0] - 0.6 * p[1] * p[1]
         + 0.25 * p[0] * p[1];
}

static double grau3(const Ponto2 &p)
{
  return grau2(p) + 4.0 * p[0] * p[0] * p[0] - 2.5 * p[1] * p[1] * p[1];
}

//! Area total de um bloco, somada folha a folha.
static double area_do_bloco(const Malha &m, long b)
{
  double a = 0.0;
  for (long i = 0; i < m.num_folhas_do_bloco(b); i++) a += m.area(b, i);
  return a;
}

//! A pior violacao de 2:1 na malha: maior |nivel(folha) - nivel(vizinho)|.
static int pior_salto(const Malha &m, long *bloco_pior, long *folha_pior)
{
  int pior = 0;
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const int nv = m.nivel(b, i);
      for (int fa = 0; fa < 4; fa++)
        for (int nviz : m.niveis_vizinhos(b, i, fa)) {
          const int d = (nv > nviz) ? (nv - nviz) : (nviz - nv);
          if (d > pior) {
            pior = d;
            if (bloco_pior) *bloco_pior = b;
            if (folha_pior) *folha_pior = i;
          }
        }
    }
  return pior;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Dois blocos lado a lado, ambos em nivel 2, e o refino pedido com DOIS niveis
  // de salto (2 -> 4) numa caixa pequena perto do canto do bloco 0.  Dois niveis
  // de uma vez e' o que faz o equilibrio importar: um nivel so' nunca produz
  // salto maior que 2:1, e o teste passaria sem o balance.
  const std::vector<Bloco> blocos = {
      {{0.0, 0.0}, {1.0, 1.0}, 2},
      {{1.0, 0.0}, {2.0, 1.0}, 2},
  };

  Malha m(blocos, sc_MPI_COMM_SELF);

  const long folhas_antes_0 = m.num_folhas_do_bloco(0);
  const long folhas_antes_1 = m.num_folhas_do_bloco(1);
  const double area_antes_0 = area_do_bloco(m, 0);
  const double area_antes_1 = area_do_bloco(m, 1);

  t_case("a_malha_comeca_uniforme_em_cada_bloco");
  T_CHECK_MSG(folhas_antes_0 == 16 && folhas_antes_1 == 16,
              "esperava 16 e 16 folhas, obtive %ld e %ld", folhas_antes_0, folhas_antes_1);
  T_CHECK_MSG(m.nivel_minimo() == 2 && m.nivel_maximo() == 2,
              "niveis antes do refino: min %d, max %d (esperava 2 e 2)",
              m.nivel_minimo(), m.nivel_maximo());

  // A caixa cobre o canto superior direito do bloco 0, encostando na interface
  // com o bloco 1 -- de proposito: e' ali que o refino poderia estragar a juncao.
  const Ponto2 cx_lo = {0.80, 0.80};
  const Ponto2 cx_hi = {1.00, 1.00};
  m.refina_caixa(cx_lo, cx_hi, 4);

  t_case("o_refino_subiu_ao_nivel_pedido");
  T_CHECK_MSG(m.nivel_maximo() == 4, "nivel maximo %d, esperava 4", m.nivel_maximo());

  t_case("o_refino_nao_rebaixou_nada");
  T_CHECK_MSG(m.nivel_minimo() == 2, "nivel minimo %d, esperava 2 (nada rebaixa)",
              m.nivel_minimo());

  // REFINOU ONDE SE PEDIU, e a formulacao certa disto NAO e' a obvia.
  //
  // A primeira versao deste caso afirmava "toda folha de nivel 4 tem o centro
  // dentro da caixa", e REPROVOU com 7 folhas fora.  O teste estava errado, nao o
  // codigo: o criterio e' aplicado a' celula ANTES de dividir, e os filhos herdam
  // a decisao do pai.  Uma celula de nivel 3 com centro na caixa gera quatro
  // filhos de nivel 4, e os de fora ficam a meia celula do pai -- ou seja, a
  // regiao refinada e' a caixa DILATADA por ate' uma celula do nivel final.
  //
  // O que vale afirmar, entao, e' que o refino maximo nao VAZA: a distancia do
  // centro de cada folha fina ate' a caixa e' no maximo o lado dela.  Isso ainda
  // reprova refino aplicado na regiao errada, que e' o defeito que importa.
  t_case("o_refino_maximo_nao_vaza_alem_de_uma_celula_da_caixa");
  {
    long longe = 0;
    double pior = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        if (m.nivel(b, i) != 4) continue;
        const Ponto2 c = m.centroide(b, i);
        const double lado = std::sqrt(m.area(b, i));
        double d = 0.0;   // distancia de Chebyshev do centro ate' a caixa
        for (int k = 0; k < 2; k++) {
          const double lo = cx_lo[k], hi = cx_hi[k];
          const double dk = (c[k] < lo) ? (lo - c[k]) : ((c[k] > hi) ? (c[k] - hi) : 0.0);
          if (dk > d) d = dk;
        }
        if (d > pior) pior = d;
        if (d > lado + 1e-12) longe++;
      }
    T_CHECK_MSG(longe == 0,
                "%ld folhas de nivel 4 a mais de um lado de celula da caixa "
                "(pior distancia %.6f)", longe, pior);
  }

  t_case("o_bloco_que_a_caixa_nao_alcanca_ganhou_so_transicao");
  {
    // A caixa toca o bloco 1 apenas no plano x=1, que nao tem area dentro dele:
    // nenhuma folha dele chega a nivel 4, mas o equilibrio pode inserir nivel 3
    // junto da interface -- e isso e' correto, nao defeito.
    int maior = 0;
    for (long i = 0; i < m.num_folhas_do_bloco(1); i++)
      if (m.nivel(1, i) > maior) maior = m.nivel(1, i);
    T_CHECK_MSG(maior <= 3, "bloco 1 chegou a nivel %d; esperava no maximo 3 (transicao)",
                maior);
  }

  // AREA: o oraculo que a contagem nao da'.  Refino divide celula; a soma das
  // areas nao muda, e nao muda em NENHUM digito que importe.
  t_case("o_refino_conserva_a_area_de_cada_bloco");
  T_NEAR(area_do_bloco(m, 0), area_antes_0, 1e-13, "area do bloco 0 depois do refino");
  T_NEAR(area_do_bloco(m, 1), area_antes_1, 1e-13, "area do bloco 1 depois do refino");
  T_NEAR(area_do_bloco(m, 0), 1.0, 1e-13, "area do bloco 0 contra a caixa dele");
  T_NEAR(area_do_bloco(m, 1), 1.0, 1e-13, "area do bloco 1 contra a caixa dele");

  t_case("o_refino_aumentou_as_folhas_so_onde_pediu");
  T_CHECK_MSG(m.num_folhas_do_bloco(0) > folhas_antes_0,
              "bloco 0 continuou com %ld folhas", m.num_folhas_do_bloco(0));

  // A RAZAO 2:1, sobre todas as folhas e todas as faces.  E' o que o
  // `t8_forest_set_balance` existe para garantir, e e' o unico caso aqui que o
  // detecta.
  t_case("a_razao_entre_niveis_vizinhos_e_no_maximo_2_para_1");
  {
    long bp = -1, fp = -1;
    const int pior = pior_salto(m, &bp, &fp);
    T_CHECK_MSG(pior <= 1,
                "salto de %d niveis entre vizinhos (bloco %ld, folha %ld); 2:1 "
                "admite no maximo 1 -- sem o balance o salto e' o que o criterio pedir",
                pior, bp, fp);
  }

  // AS JUNCOES DEPOIS DO REFINO.  A interface entre os blocos tem agora celulas
  // de niveis diferentes dos dois lados; o que NAO pode e' aparecer face sem
  // vizinho ali, que seria contorno fisico no meio do dominio.
  t_case("a_interface_entre_blocos_continua_com_vizinho_depois_do_refino");
  {
    long sem_vizinho = 0, com_vizinho = 0;
    const Bloco &k0 = m.bloco(0);
    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const Ponto2 c = m.centroide(0, i);
      const double lado = std::sqrt(m.area(0, i));
      // folha encostada na face x+ do bloco 0, que e' a interface
      if (c[0] > k0.hi[0] - lado) {
        if (m.num_vizinhos(0, i, 1) > 0) com_vizinho++; else sem_vizinho++;
      }
    }
    T_CHECK_MSG(sem_vizinho == 0 && com_vizinho > 0,
                "na interface depois do refino: %ld com vizinho, %ld SEM "
                "(sem vizinho ali seria contorno no meio do dominio)",
                com_vizinho, sem_vizinho);
  }

  // O DISCRETIZADOR SOBRE A MALHA REFINADA, com o suporte atravessando a
  // transicao de nivel -- que e' onde um estencil mal formado apareceria.
  t_case("mls_reproduz_polinomio_sobre_suporte_com_transicao_de_nivel");
  double err2 = 0.0, err3 = 0.0;
  {
    const std::vector<Ponto2> todos = m.centroides();
    const int ORDEM = 2;
    const int MAXPTS = 4 * Mls::min_pontos(ORDEM);
    const Ponto2 alvo = {0.78, 0.83};   // junto da borda da caixa de refino

    std::vector<int> ordem(todos.size());
    for (size_t i = 0; i < todos.size(); i++) ordem[i] = static_cast<int>(i);
    const auto d2 = [&](int k) {
      const double dx = todos[k][0] - alvo[0], dy = todos[k][1] - alvo[1];
      return dx * dx + dy * dy;
    };
    for (int k = 0; k < MAXPTS && k < static_cast<int>(ordem.size()); k++) {
      int melhor = k;
      for (int i = k + 1; i < static_cast<int>(ordem.size()); i++)
        if (d2(ordem[i]) < d2(ordem[melhor])) melhor = i;
      const int t = ordem[k]; ordem[k] = ordem[melhor]; ordem[melhor] = t;
    }

    std::vector<Ponto2> suporte;
    std::vector<double> v2, v3;
    const int numpts = (MAXPTS < static_cast<int>(todos.size()))
                           ? MAXPTS : static_cast<int>(todos.size());
    for (int i = 0; i < numpts; i++) {
      suporte.push_back(todos[ordem[i]]);
      v2.push_back(grau2(suporte.back()));
      v3.push_back(grau3(suporte.back()));
    }

    Mls mls(ORDEM, numpts);
    const double obtido2 = mls.interpola(alvo, suporte, v2);
    err2 = std::fabs(obtido2 - grau2(alvo));
    T_NEAR(obtido2, grau2(alvo), 1e-9, "campo de grau 2 com suporte atravessando niveis");
    err3 = std::fabs(mls.interpola(alvo, suporte, v3) - grau3(alvo));
  }

  t_case("suporte_com_transicao_erra_mais_no_grau_acima");
  T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
              "grau 2 contra grau 3 com suporte atravessando niveis");

  // AS RECUSAS do refino.
  t_case("caixa_degenerada_e_recusada");
  {
    bool lancou = false;
    try { m.refina_caixa({0.5, 0.5}, {0.5, 0.9}, 5); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "caixa de largura zero passou sem reclamar");
  }

  t_case("caixa_fora_de_todos_os_blocos_e_recusada");
  {
    bool lancou = false;
    try { m.refina_caixa({5.0, 5.0}, {6.0, 6.0}, 5); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "caixa que nao toca bloco nenhum passou sem reclamar");
  }

  // Refino que NAO refinaria nada e' engano de quem chamou, e aceito calado
  // deixaria a malha PARECER refinada.
  t_case("nivel_alvo_que_nao_aprofunda_e_recusado");
  {
    bool lancou = false;
    try { m.refina_caixa({0.1, 0.1}, {0.3, 0.3}, 2); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "nivel_alvo igual ao do bloco passou sem reclamar");
  }

  // ---- REFINO ALEATORIO -----------------------------------------------------
  //
  // E' a malha que o teste de convergencia em quadtree precisa.  O que estes casos
  // afirmam nao e' "o refino aconteceu" -- e' que ele e' REPRODUTIVEL, que a
  // semente faz diferenca, e que a malha resultante tem transicao 2:1 DE VERDADE.
  // Sem a ultima, um teste de esquema na face de transicao poderia rodar numa
  // malha sem transicao nenhuma e passar sem testar nada.

  t_case("o_refino_aleatorio_e_reprodutivel_pela_semente");
  {
    const std::vector<Bloco> bl = {{{0.0, 0.0}, {1.0, 1.0}, 3}};
    Malha a(bl, sc_MPI_COMM_SELF);
    Malha b(bl, sc_MPI_COMM_SELF);
    a.refina_aleatorio(0.3, 12345, 5);
    b.refina_aleatorio(0.3, 12345, 5);

    T_CHECK_MSG(a.num_folhas() == b.num_folhas(),
                "a mesma semente deu %ld e %ld folhas", a.num_folhas(), b.num_folhas());

    // Contagem igual nao basta: as malhas tem de ser a MESMA.  Duas malhas
    // diferentes com o mesmo numero de folhas sao perfeitamente possiveis, e a
    // contagem nao as distingue.
    const std::vector<Ponto2> ca = a.centroides();
    const std::vector<Ponto2> cb = b.centroides();
    double pior = 0.0;
    for (size_t k = 0; k < ca.size() && k < cb.size(); k++)
      for (int d = 0; d < 2; d++)
        pior = std::fmax(pior, std::fabs(ca[k][d] - cb[k][d]));
    T_NEAR(pior, 0.0, 0.0, "pior diferenca entre os centroides das duas malhas");
  }

  t_case("sementes_diferentes_dao_malhas_diferentes");
  {
    // Se a semente nao mudasse nada, o caso de cima estaria afirmando apenas que o
    // codigo e' deterministico -- o que qualquer codigo sem aleatoriedade e'.
    const std::vector<Bloco> bl = {{{0.0, 0.0}, {1.0, 1.0}, 3}};
    Malha a(bl, sc_MPI_COMM_SELF);
    Malha b(bl, sc_MPI_COMM_SELF);
    a.refina_aleatorio(0.3, 12345, 5);
    b.refina_aleatorio(0.3, 999, 5);
    T_CHECK_MSG(a.num_folhas() != b.num_folhas(),
                "sementes diferentes deram o MESMO numero de folhas (%ld) -- a "
                "semente pode nao estar sendo usada", a.num_folhas());
    printf("  (semente 12345: %ld folhas; semente 999: %ld folhas)\n",
           a.num_folhas(), b.num_folhas());
  }

  t_case("o_refino_aleatorio_respeita_2_para_1_e_tem_transicao_de_verdade");
  {
    Malha m2({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    const double area_antes = area_do_bloco(m2, 0);
    m2.refina_aleatorio(0.4, 2026, 5);
    m2.refina_aleatorio(0.4, 777, 6);        // segunda passada, outra semente

    long bp = -1, fp = -1;
    const int pior = pior_salto(m2, &bp, &fp);
    T_CHECK_MSG(pior <= 1,
                "salto de %d niveis entre vizinhos (bloco %ld folha %ld) depois do "
                "refino aleatorio -- o equilibrio nao fechou", pior, bp, fp);

    // AREA: o oraculo que a contagem nao da'.
    T_NEAR(area_do_bloco(m2, 0), area_antes, 1e-13,
           "area do bloco depois de duas passadas de refino aleatorio");

    // E A MALHA TEM DE TER TRANSICAO: sem isso, qualquer teste de esquema na face
    // de transicao rodaria sem exercitar o que afirma.
    const long trans = m2.num_folhas_em_transicao();
    T_CHECK_MSG(trans > 0,
                "a malha aleatoria nao tem NENHUMA folha em transicao 2:1 -- ela nao"
                " serve para testar esquema na face de transicao");
    T_CHECK_MSG(m2.nivel_maximo() > m2.nivel_minimo(),
                "niveis min %d e max %d -- a malha ficou uniforme",
                m2.nivel_minimo(), m2.nivel_maximo());
    printf("  (malha aleatoria: %ld folhas, niveis %d a %d, %ld em transicao = %.0f%%)\n",
           m2.num_folhas(), m2.nivel_minimo(), m2.nivel_maximo(), trans,
           100.0 * trans / m2.num_folhas());
  }

  t_case("a_fracao_controla_quanto_se_refina");
  {
    // Fracao maior tem de refinar mais.  Sem este caso, a fracao poderia estar
    // sendo ignorada e os casos de cima passariam igual.
    const std::vector<Bloco> bl = {{{0.0, 0.0}, {1.0, 1.0}, 3}};
    Malha pouco(bl, sc_MPI_COMM_SELF), muito(bl, sc_MPI_COMM_SELF);
    pouco.refina_aleatorio(0.1, 55, 5);
    muito.refina_aleatorio(0.8, 55, 5);
    T_CHECK_MSG(muito.num_folhas() > pouco.num_folhas(),
                "fracao 0,8 deu %ld folhas e fracao 0,1 deu %ld -- a fracao nao "
                "esta' controlando nada", muito.num_folhas(), pouco.num_folhas());
    printf("  (fracao 0,1: %ld folhas; fracao 0,8: %ld folhas)\n",
           pouco.num_folhas(), muito.num_folhas());
  }

  t_case("refino_aleatorio_malformado_e_recusado");
  {
    Malha m3({{{0.0, 0.0}, {1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    bool lancou = false;
    try { m3.refina_aleatorio(0.0, 1, 5); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "fracao zero passou sem reclamar");

    lancou = false;
    try { m3.refina_aleatorio(1.5, 1, 5); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "fracao maior que 1 passou sem reclamar");

    lancou = false;
    try { m3.refina_aleatorio(0.5, 1, 3); }   // nivel_max == nivel atual
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "nivel_max que nao aprofunda passou sem reclamar");
  }

  return t_end();
}
