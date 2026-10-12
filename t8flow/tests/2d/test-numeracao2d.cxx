// A NUMERACAO GLOBAL: contigua por processo, determinada pela GEOMETRIA, e
// consistente entre quem escreve e quem le'.
//
// ESTE ARQUIVO EXISTE POR CAUSA DE UMA TAUTOLOGIA MEDIDA NO HIGFLOW.  La', a
// identidade entre o id do mapeador e o contador do percurso e' verdadeira por
// construcao, e trocar `mp_lookup` por `n_visitadas++` deixava a suite INTEIRA
// verde -- nos dois DIM e nos tres np.  Qualquer teste que verifique a numeracao
// percorrendo na mesma ordem que a gerou tem esse defeito.
//
// As tres vias independentes usadas aqui, em ordem de forca:
//
//   1. ARITMETICA DE MORTON.  Num bloco uniforme de nivel L, a folha de indice k
//      esta' na posicao obtida DESINTERCALANDO os bits de k.  E' a definicao da
//      curva, calculada aqui sem olhar a malha, e comparada com a geometria da
//      folha.  Pega numeracao que nao e' a curva de preenchimento;
//   2. BIJECAO.  Os ids locais cobrem [primeiro, primeiro + n_local) exatamente
//      uma vez cada.  Pega furo e duplicata, que e' o defeito mais comum;
//   3. ESCREVE POR PERCURSO, LE' POR GEOMETRIA.  Um sistema linear e' montado com
//      os ids vindos do percurso, e a solucao e' lida com os ids vindos de
//      `id_do_ponto` -- geometria, nao contagem.  Se as duas vias discordarem, o
//      valor lido na posicao e' o de OUTRA posicao.  E' o defeito do sdmult/sdp do
//      HiGFlow, que sobreviveu anos porque em malha uniforme os dois indices
//      coincidem.
//
// O campo de teste tem COEFICIENTES DISTINTOS POR DIRECAO de proposito: assim o
// valor identifica a posicao, e troca de indices nao passa por acaso.
#include "../testing.h"

#include <set>
#include <vector>

#include "arranque.hxx"
#include "malha2d.hxx"
#include "solver2d.hxx"

using t8flow2d::Bloco;
using t8flow2d::Malha;
using t8flow2d::Motor;
using t8flow2d::Ponto2;
using t8flow2d::Solver;

//! Valor que identifica a posicao: coeficientes distintos e nao comensuraveis.
static double marca(const Ponto2 &p)
{
  return 1000.0 * p[0] + 7.0 * p[1] + 0.125;
}

//! Desintercala os bits de `k`: a posicao (ix, iy) na curva de Morton 2D.
//
// Calculado aqui, da definicao da curva, sem consultar a malha -- e' isso que
// torna esta via independente.
static void morton_para_ij(long k, int nivel, int *ix, int *iy)
{
  *ix = 0;
  *iy = 0;
  for (int b = 0; b < nivel; b++) {
    // A convencao do t8code (p4est) e' x no bit MENOS significativo de cada par.
    *ix |= static_cast<int>((k >> (2 * b)) & 1) << b;
    *iy |= static_cast<int>((k >> (2 * b + 1)) & 1) << b;
  }
}

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // Dois blocos, niveis distintos: 16 + 4 = 20 folhas.  O bloco 0 e' QUADRADO e
  // uniforme, que e' o que permite a conta de Morton; o bloco 1 e' retangular, o
  // que faz os lados de celula diferirem por direcao e pega quem usar
  // `sqrt(area)`.
  const std::vector<Bloco> blocos = {
      {{0.0, 0.0}, {1.0, 1.0}, 2},
      {{1.0, 0.0}, {3.0, 1.0}, 1},
  };
  Malha m(blocos, sc_MPI_COMM_SELF);

  t_case("a_numeracao_local_e_contigua_e_comeca_no_primeiro_id");
  {
    const long p0 = m.primeiro_id();
    const long n = m.num_folhas();
    T_CHECK_MSG(p0 == 0, "em np=1 o primeiro id tem de ser 0, obtive %ld", p0);
    T_CHECK_MSG(n == 20, "esperava 20 folhas (16+4), obtive %ld", n);
    T_CHECK_MSG(m.num_folhas_globais() == n,
                "em np=1 o total global (%ld) tem de ser o local (%ld)",
                m.num_folhas_globais(), n);
  }

  // VIA 2: BIJECAO.  Cada id de [primeiro, primeiro + n) aparece exatamente uma
  // vez.  Furo e duplicata sao o defeito mais comum de numeracao, e nenhuma
  // contagem de folhas os ve.
  t_case("os_ids_cobrem_o_intervalo_local_sem_furo_nem_repeticao");
  {
    std::set<long> vistos;
    long repetidos = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const long id = m.id_global(b, i);
        if (!vistos.insert(id).second) repetidos++;
      }
    T_CHECK_MSG(repetidos == 0, "%ld ids repetidos", repetidos);
    T_CHECK_MSG(static_cast<long>(vistos.size()) == m.num_folhas(),
                "%zu ids distintos para %ld folhas", vistos.size(), m.num_folhas());
    const long p0 = m.primeiro_id();
    long fora = 0;
    for (long id : vistos)
      if (id < p0 || id >= p0 + m.num_folhas()) fora++;
    T_CHECK_MSG(fora == 0, "%ld ids fora de [%ld, %ld)", fora, p0, p0 + m.num_folhas());
  }

  // VIA 1: ARITMETICA DE MORTON.  O bloco 0 e' uniforme de nivel 2, logo a folha
  // de indice local k esta' na posicao (ix, iy) dada pelos bits de k
  // desintercalados.  A conta sai da definicao da curva; a posicao vem da
  // geometria da malha.  As duas tem de bater.
  t_case("a_ordem_local_e_a_curva_de_morton");
  {
    const Bloco &k0 = m.bloco(0);
    const int nivel = 2;
    const int div = 1 << nivel;
    const double hx = (k0.hi[0] - k0.lo[0]) / div;
    const double hy = (k0.hi[1] - k0.lo[1]) / div;

    for (long i = 0; i < m.num_folhas_do_bloco(0); i++) {
      int ix = 0, iy = 0;
      morton_para_ij(i, nivel, &ix, &iy);
      const Ponto2 esperado = {k0.lo[0] + (ix + 0.5) * hx, k0.lo[1] + (iy + 0.5) * hy};
      const Ponto2 obtido = m.centroide(0, i);
      if (std::fabs(obtido[0] - esperado[0]) > 1e-12 ||
          std::fabs(obtido[1] - esperado[1]) > 1e-12) {
        T_CHECK_MSG(false,
                    "folha %ld: Morton diz (%d,%d) -> (%.6f, %.6f), a malha da' "
                    "(%.6f, %.6f)", i, ix, iy, esperado[0], esperado[1],
                    obtido[0], obtido[1]);
        break;
      }
    }
  }

  // A BUSCA POR PONTO acha a folha certa, e o id dela e' o da folha onde o ponto
  // esta' -- conferido no centroide de cada folha, que e' o ponto menos ambiguo.
  t_case("a_busca_por_ponto_devolve_o_id_da_folha_que_contem_o_ponto");
  {
    long erradas = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const long esperado = m.id_global(b, i);
        const long obtido = m.id_do_ponto(m.centroide(b, i));
        if (obtido != esperado) {
          if (erradas == 0)
            T_CHECK_MSG(false, "bloco %ld folha %ld: id %ld, busca devolveu %ld",
                        b, i, esperado, obtido);
          erradas++;
        }
      }
    T_CHECK_MSG(erradas == 0, "%ld folhas com id divergente da busca", erradas);
  }

  t_case("ponto_fora_de_todos_os_blocos_nao_tem_id");
  T_CHECK_MSG(m.id_do_ponto({5.0, 5.0}) == -1,
              "ponto fora do dominio devolveu id %ld", m.id_do_ponto({5.0, 5.0}));

  // O DESEMPATE.  Um ponto exatamente sobre a face comum entre duas celulas
  // pertence a' de coordenada MAIOR, pelo intervalo semiaberto -- e a resposta nao
  // pode depender da ordem de visita.  Na HiGTree a clausula analoga passava por
  // SORTE: a curva visita a celula de menor coordenada primeiro, entao guardar o
  // primeiro candidato acertava por acidente.
  t_case("o_desempate_na_face_comum_e_a_celula_de_coordenada_maior");
  {
    // x = 0,25 e' a face entre a primeira e a segunda coluna de celulas do bloco
    // 0 (nivel 2, lado 0,25).  y = 0,125 esta' no meio da primeira linha.
    const long id_face = m.id_do_ponto({0.25, 0.125});
    const long id_direita = m.id_do_ponto({0.26, 0.125});
    const long id_esquerda = m.id_do_ponto({0.24, 0.125});
    T_CHECK_MSG(id_face == id_direita,
                "ponto na face deu id %ld; a celula de coordenada maior e' %ld "
                "(a de menor e' %ld)", id_face, id_direita, id_esquerda);
    T_CHECK_MSG(id_direita != id_esquerda,
                "as duas celulas vizinhas tem o mesmo id (%ld) -- o caso nao "
                "discrimina nada", id_direita);
  }

  t_case("o_limite_superior_do_dominio_pertence_a_ultima_celula");
  {
    // No limite de cima o semiaberto deixaria o ponto sem celula; a ultima o
    // inclui.  Sem esta excecao, o canto do dominio nao teria id.
    const long id = m.id_do_ponto({1.0, 1.0});
    T_CHECK_MSG(id >= 0, "o canto superior do bloco 0 nao tem id (%ld)", id);
  }

  // VIA 3: ESCREVE POR PERCURSO, LE' POR GEOMETRIA.
  //
  // Sistema identidade: x_i = marca(centro da folha i), montado com os ids do
  // PERCURSO.  Depois, para cada folha, o valor e' lido no id que a GEOMETRIA
  // daquele centro devolve.  Se as duas vias discordarem, o valor lido e' o de
  // outra posicao -- exatamente o defeito do sdmult/sdp, que sobreviveu anos
  // porque em malha uniforme os dois indices coincidem celula a celula.
  t_case("o_valor_escrito_por_id_de_percurso_e_lido_por_id_de_geometria");
  {
    const long n = m.num_folhas();
    Solver s(Motor::Petsc, static_cast<int>(m.primeiro_id()),
             static_cast<std::size_t>(n), 1);
    s.tolerancia(1e-14);
    s.max_iteracoes(100);

    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const int id = static_cast<int>(m.id_global(b, i));
        s.linha(id, {id}, {1.0});
        s.lado_direito(id, marca(m.centroide(b, i)));
      }
    s.monta();
    s.resolve();
    T_CHECK_MSG(s.convergiu(), "o sistema identidade nao convergiu: %d iteracoes",
                s.iteracoes());

    double pior = 0.0;
    long pior_id = -1;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto2 c = m.centroide(b, i);
        const long id_geo = m.id_do_ponto(c);
        T_CHECK_MSG(id_geo >= 0, "o centro da folha %ld do bloco %ld nao tem id", i, b);
        const double lido = s.x(static_cast<int>(id_geo));
        const double e = std::fabs(lido - marca(c));
        if (e > pior) { pior = e; pior_id = id_geo; }
      }
    T_CHECK_MSG(pior <= 1e-10,
                "pior divergencia %.3g no id %ld -- o valor lido pela geometria "
                "nao e' o escrito pelo percurso", pior, pior_id);
  }

  // E A MESMA COISA DEPOIS DO REFINO, que e' quando os dois caminhos tem a chance
  // real de divergir: a numeracao muda inteira e os niveis deixam de ser
  // uniformes.
  t_case("a_consistencia_sobrevive_ao_refino_local");
  {
    m.refina_caixa({0.60, 0.60}, {1.00, 1.00}, 4);

    std::set<long> vistos;
    long repetidos = 0, divergentes = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const long id = m.id_global(b, i);
        if (!vistos.insert(id).second) repetidos++;
        if (m.id_do_ponto(m.centroide(b, i)) != id) divergentes++;
      }
    T_CHECK_MSG(repetidos == 0, "%ld ids repetidos depois do refino", repetidos);
    T_CHECK_MSG(static_cast<long>(vistos.size()) == m.num_folhas(),
                "%zu ids distintos para %ld folhas depois do refino",
                vistos.size(), m.num_folhas());
    T_CHECK_MSG(divergentes == 0,
                "%ld folhas com id divergente da busca depois do refino", divergentes);
  }

  return t_end();
}
