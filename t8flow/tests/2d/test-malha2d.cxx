// A MALHA MULTIBLOCOS 2D: uma arvore por bloco, nivel por bloco, geometria
// conferida pelo centroide, e as interfaces unidas.
//
// ORACULOS, e cada um existe porque um defeito especifico nao e' visto pelos
// outros:
//
//   contagem de folhas por bloco  pega nivel aplicado na arvore errada -- e so'
//                                 discrimina se os niveis forem DIFERENTES entre
//                                 blocos.  Com todos iguais, uma malha uniforme
//                                 global passaria;
//   centroide dentro da caixa     pega geometria que SAI da caixa -- vertice no
//                                 bloco errado, caixa trocada.  NAO pega
//                                 espelhamento: trocar v1 com v2 reflete o
//                                 mapeamento DENTRO da mesma caixa, e todo
//                                 centroide continua dentro dela.  Isso foi
//                                 MEDIDO, nao suposto: a primeira versao deste
//                                 arquivo afirmava que o centroide pegava a troca
//                                 de vertices, e a sabotagem mostrou que nao;
//   ordem de percurso em x        pega o espelhamento.  A curva do t8code varre x
//                                 primeiro, entao a folha 1 esta' deslocada em x
//                                 da folha 0 -- com os vertices trocados, em y;
//   vizinhos por face             pega join ausente.  Sem `t8_cmesh_set_join`, a
//                                 malha tambem commita, tambem tem as folhas
//                                 certas, e toda face de borda de bloco passa a
//                                 ser "sem vizinho" -- que e' o criterio de
//                                 CONTORNO.  Nenhuma contagem pega isso.
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

static long quatro_elevado(int n)
{
  long v = 1;
  for (int i = 0; i < n; i++) v *= 4;
  return v;
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Tres blocos em L, com NIVEIS DISTINTOS e TAMANHOS distintos.  Os niveis
  // diferentes sao o que faz a contagem discriminar; os tamanhos diferentes fazem
  // o centroide discriminar (com blocos de mesmo tamanho, trocar a geometria de
  // dois deles seria invisivel).
  //
  //   b0: [0,1]x[0,1]  nivel 2    b1 encosta em b0 pela face x+
  //   b1: [1,3]x[0,1]  nivel 1    e b2 encosta em b0 pela face y+
  //   b2: [0,1]x[1,2]  nivel 3
  const std::vector<Bloco> blocos = {
      {{0.0, 0.0}, {1.0, 1.0}, 2},
      {{1.0, 0.0}, {3.0, 1.0}, 1},
      {{0.0, 1.0}, {1.0, 2.0}, 3},
  };

  Malha m(blocos, sc_MPI_COMM_SELF);

  t_case("uma_arvore_por_bloco");
  T_CHECK_MSG(m.num_blocos() == 3, "esperava 3 blocos, obtive %ld", m.num_blocos());

  t_case("folhas_por_bloco_seguem_o_nivel_do_bloco");
  for (long b = 0; b < m.num_blocos(); b++) {
    const long esperado = quatro_elevado(m.bloco(b).nivel);
    const long obtido = m.num_folhas_do_bloco(b);
    T_CHECK_MSG(obtido == esperado,
                "bloco %ld (nivel %d): esperava %ld folhas, obtive %ld",
                b, m.bloco(b).nivel, esperado, obtido);
  }

  t_case("total_de_folhas_e_a_soma_dos_blocos");
  {
    long soma = 0;
    for (long b = 0; b < m.num_blocos(); b++) soma += m.num_folhas_do_bloco(b);
    T_CHECK_MSG(m.num_folhas() == soma, "total %ld, soma dos blocos %ld",
                m.num_folhas(), soma);
    // 16 + 4 + 64 = 84, escrito a mao para o caso nao se validar sozinho.
    T_CHECK_MSG(soma == 84, "esperava 84 folhas (16+4+64), obtive %ld", soma);
  }

  t_case("centroide_cai_na_caixa_do_seu_bloco");
  for (long b = 0; b < m.num_blocos(); b++) {
    const Bloco &k = m.bloco(b);
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m.centroide(b, i);
      const bool dentro = c[0] > k.lo[0] && c[0] < k.hi[0] &&
                          c[1] > k.lo[1] && c[1] < k.hi[1];
      if (!dentro) {
        T_CHECK_MSG(false,
                    "bloco %ld folha %ld: centroide (%.6f, %.6f) fora da caixa "
                    "[%.3f,%.3f]x[%.3f,%.3f]",
                    b, i, c[0], c[1], k.lo[0], k.hi[0], k.lo[1], k.hi[1]);
        break;
      }
    }
  }

  // O TAMANHO DA CELULA e' consequencia do nivel E da caixa, e conferi-lo
  // separado da contagem pega o caso em que o numero de folhas esta' certo e a
  // geometria nao: o bloco 1 tem o dobro da largura e metade do nivel, logo
  // celulas quatro vezes mais largas que as do bloco 0.
  t_case("celulas_do_bloco_largo_sao_maiores_na_proporcao_certa");
  {
    // Menor distancia entre centroides vizinhos aproxima o lado da celula; aqui
    // basta a extensao: com 4^n folhas num quadrado, o lado e' (hi-lo)/2^n.
    const double lado0 = (m.bloco(0).hi[0] - m.bloco(0).lo[0]) / (1 << m.bloco(0).nivel);
    const double lado1 = (m.bloco(1).hi[0] - m.bloco(1).lo[0]) / (1 << m.bloco(1).nivel);
    // O bloco 1 e' 2x mais largo e tem 1 nivel menos -> 4x o lado.
    T_NEAR(lado1 / lado0, 4.0, 1e-12, "razao entre os lados de celula (bloco 1 / bloco 0)");
    // E o centroide da PRIMEIRA folha de cada bloco tem de estar a meio lado do
    // canto inferior do bloco -- o que amarra geometria e nivel juntos.
    const Ponto2 c0 = m.centroide(0, 0);
    T_NEAR(c0[0] - m.bloco(0).lo[0], lado0 / 2, 1e-12, "x do centroide da 1a folha do bloco 0");
    T_NEAR(c0[1] - m.bloco(0).lo[1], lado0 / 2, 1e-12, "y do centroide da 1a folha do bloco 0");
  }

  // A ORIENTACAO DO MAPEAMENTO, que o centroide-na-caixa nao ve.  A curva de
  // preenchimento do t8code varre x primeiro, logo a folha 0 de uma arvore esta'
  // no canto (lo, lo) e a folha 1 esta' deslocada em X dela, com o MESMO y.
  // Trocar a ordem dos vertices do QUAD espelha o mapeamento na diagonal: a caixa
  // continua a mesma, o conjunto de centroides continua o mesmo, e so' este caso
  // acusa.
  t_case("a_ordem_de_percurso_varre_x_primeiro");
  for (long b = 0; b < m.num_blocos(); b++) {
    const Bloco &k = m.bloco(b);
    const int div = 1 << k.nivel;
    const double hx = (k.hi[0] - k.lo[0]) / div;
    const double hy = (k.hi[1] - k.lo[1]) / div;
    if (div < 2) continue;   // com uma celula por direcao nao ha' segunda folha

    const Ponto2 c0 = m.centroide(b, 0);
    const Ponto2 c1 = m.centroide(b, 1);
    T_NEAR(c0[0], k.lo[0] + hx / 2, 1e-12, "x da folha 0 (canto inferior do bloco)");
    T_NEAR(c0[1], k.lo[1] + hy / 2, 1e-12, "y da folha 0 (canto inferior do bloco)");
    T_NEAR(c1[0] - c0[0], hx, 1e-12, "a folha 1 esta' deslocada um lado em X");
    T_NEAR(c1[1] - c0[1], 0.0, 1e-12, "a folha 1 tem o MESMO y da folha 0");
  }

  // AS JUNCOES.  b0-b1 encostam em x (faces coincidentes, [0,1] em y) e b0-b2
  // encostam em y ([0,1] em x).  b1-b2 so' se tocam no canto, que nao e' face.
  t_case("as_interfaces_coincidentes_foram_unidas");
  T_CHECK_MSG(m.num_juncoes() == 2, "esperava 2 juncoes (b0-b1 e b0-b2), obtive %ld",
              m.num_juncoes());

  // E O QUE A JUNCAO SIGNIFICA, medido no lugar que importa: a face que da' para
  // outro bloco TEM vizinho, e a face que da' para fora NAO tem.  Sem o join a
  // primeira metade falha -- e nenhuma contagem de folhas acusaria.
  t_case("face_entre_blocos_tem_vizinho_e_face_externa_nao");
  {
    const Bloco &k0 = m.bloco(0);
    const double lado0 = (k0.hi[0] - k0.lo[0]) / (1 << k0.nivel);
    long com_vizinho_na_interface = 0, sem_vizinho_na_interface = 0;
    long com_vizinho_fora = 0, sem_vizinho_fora = 0;

    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      const Ponto2 c = m.centroide(0, i);
      // face x+ do bloco 0: interface com o bloco 1
      if (c[0] > k0.hi[0] - lado0) {
        (m.num_vizinhos(0, i, 1) > 0 ? com_vizinho_na_interface : sem_vizinho_na_interface)++;
      }
      // face x- do bloco 0: borda externa do dominio
      if (c[0] < k0.lo[0] + lado0) {
        (m.num_vizinhos(0, i, 0) > 0 ? com_vizinho_fora : sem_vizinho_fora)++;
      }
    }

    T_CHECK_MSG(com_vizinho_na_interface == 4 && sem_vizinho_na_interface == 0,
                "na interface com o bloco 1: %ld com vizinho, %ld sem "
                "(esperava 4 e 0 -- sem vizinho ali seria contorno no meio do dominio)",
                com_vizinho_na_interface, sem_vizinho_na_interface);
    T_CHECK_MSG(sem_vizinho_fora == 4 && com_vizinho_fora == 0,
                "na borda externa: %ld sem vizinho, %ld com (esperava 4 e 0)",
                sem_vizinho_fora, com_vizinho_fora);
  }

  // O CONTRATO MALHA/DISCRETIZACAO sobre a malha multiblocos: o suporte sai dos
  // centroides e os pesos do MLS.  Alvo DENTRO do bloco mais fino, assimetrico
  // pelo motivo registrado no test-mls2d.
  t_case("mls_reproduz_polinomio_sobre_suporte_da_malha");
  double err2 = 0.0, err3 = 0.0;
  {
    const std::vector<Ponto2> todos = m.centroides();
    T_CHECK_MSG(static_cast<long>(todos.size()) == m.num_folhas(),
                "centroides() deu %zu para %ld folhas", todos.size(), m.num_folhas());

    const int ORDEM = 2;
    const int MAXPTS = 4 * Mls::min_pontos(ORDEM);
    const Ponto2 alvo = {0.37, 1.41};   // interior do bloco 2, o de nivel 3

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
    T_NEAR(obtido2, grau2(alvo), 1e-9, "campo de grau 2 sobre a malha multiblocos");
    err3 = std::fabs(mls.interpola(alvo, suporte, v3) - grau3(alvo));
  }

  t_case("suporte_da_malha_erra_mais_no_grau_acima");
  T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
              "grau 2 contra grau 3 sobre o suporte da malha multiblocos");

  // AS RECUSAS.  Cada uma e' um caso que, aceito calado, daria malha plausivel e
  // errada.
  t_case("lista_vazia_e_recusada");
  {
    bool lancou = false;
    try { Malha vazia({}, sc_MPI_COMM_SELF); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "lista de blocos vazia passou sem reclamar");
  }

  t_case("bloco_degenerado_e_recusado");
  {
    bool lancou = false;
    try { Malha d({{{0.0, 0.0}, {0.0, 1.0}, 1}}, sc_MPI_COMM_SELF); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "bloco de largura zero passou sem reclamar");
  }

  t_case("blocos_sobrepostos_sao_recusados");
  {
    bool lancou = false;
    try {
      Malha s({{{0.0, 0.0}, {2.0, 2.0}, 1}, {{1.0, 1.0}, {3.0, 3.0}, 1}}, sc_MPI_COMM_SELF);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "blocos com area em comum passaram sem reclamar");
  }

  t_case("faces_parcialmente_coincidentes_sao_recusadas");
  {
    // b1 encosta em b0 por x, mas cobre so' METADE da face: o join do t8code nao
    // representa isso.
    bool lancou = false;
    try {
      Malha p({{{0.0, 0.0}, {1.0, 1.0}, 1}, {{1.0, 0.0}, {2.0, 0.5}, 1}}, sc_MPI_COMM_SELF);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "faces com sobreposicao parcial passaram sem reclamar");
  }

  return t_end();
}
