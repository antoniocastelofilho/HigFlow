// O ARRANJO DESLOCADO: escalares e tensores nos centros, vetores nas facetas.
//
// ORACULOS, e cada um existe porque os outros nao o veem:
//
//   contagem de facetas    N(N+1) por direcao em malha uniforme N x N -- a conta
//                          classica do MAC, calculada no teste a partir de N e
//                          comparada com o que o dominio produziu.  Pega posse
//                          errada (contar a faceta nas duas pontas daria 2N^2);
//   POSSE UNICA            cada faceta e' dona de exatamente UMA folha.  A
//                          contagem certa nao garante isso: duas trocas que se
//                          compensam dariam o mesmo total;
//   geometria da faceta    o centro esta' a meia celula do centro da celula, na
//                          direcao da normal, e a MEDIDA e' o lado TRANSVERSAL --
//                          pega quem usar o lado errado em celula nao-quadrada;
//   fluxo de um campo      um campo vetorial constante tem divergencia ZERO em
//     constante            toda celula interna.  E' o oraculo que amarra posse,
//                          geometria e sinal de uma vez, e nenhuma contagem o da';
//   independencia das      escrever numa componente do tensor nao altera as
//     componentes          outras, e escrever numa direcao do vetor nao altera a
//                          outra.  Pega layout sobreposto, que passa em todo
//                          teste de valor unico.
#include "../testing.h"

#include <set>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "campos2d.hxx"
#include "malha2d.hxx"

using t8flow2d::Bloco;
using t8flow2d::CampoEscalar;
using t8flow2d::CampoTensor;
using t8flow2d::CampoVetorial;
using t8flow2d::DominioFacetas;
using t8flow2d::Malha;
using t8flow2d::Ponto2;

//! Marca que identifica a posicao: coeficientes distintos por direcao.
static double marca(const Ponto2 &p) { return 100.0 * p[0] + 3.0 * p[1] + 0.5; }

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  // UM bloco uniforme para as contagens fecharem com a conta do MAC, e um
  // segundo bloco RETANGULAR para a medida da faceta discriminar: nele a celula
  // nao e' quadrada, entao usar o lado errado como medida aparece.
  const std::vector<Bloco> blocos = {
      {{0.0, 0.0}, {1.0, 1.0}, 2},   // 4x4 celulas de 0,25 x 0,25
      {{1.0, 0.0}, {3.0, 1.0}, 1},   // 2x2 celulas de 1,00 x 0,50
  };
  Malha m(blocos, sc_MPI_COMM_SELF);
  DominioFacetas fd(m);

  // --- os campos de centro ---------------------------------------------------

  t_case("campo_escalar_tem_um_valor_por_celula");
  {
    CampoEscalar p(m);
    T_CHECK_MSG(p.tamanho() == m.num_folhas(), "campo com %ld valores para %ld folhas",
                p.tamanho(), m.num_folhas());
  }

  t_case("escalar_escrito_por_id_e_lido_por_geometria");
  {
    CampoEscalar p(m);
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        p.em(m.id_global(b, i)) = marca(m.centroide(b, i));

    double pior = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto2 c = m.centroide(b, i);
        const double e = std::fabs(p.em_ponto(c) - marca(c));
        if (e > pior) pior = e;
      }
    T_NEAR(pior, 0.0, 0.0, "pior divergencia entre escrita por id e leitura por ponto");
  }

  t_case("escalar_recusa_id_fora_da_faixa_local");
  {
    CampoEscalar p(m);
    bool lancou = false;
    try { p.em(m.primeiro_id() + m.num_folhas()); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "id fora da faixa passou sem reclamar");
  }

  // O TENSOR: quatro componentes independentes.  Escrever numa nao pode mexer nas
  // outras -- layout sobreposto passa em qualquer teste de valor unico.
  t_case("componentes_do_tensor_sao_independentes");
  {
    CampoTensor T(m);
    T.preenche(0.0);
    const long id = m.id_global(0, 0);
    T.em(id, 0, 1) = 7.0;
    T_NEAR(T.em(id, 0, 1), 7.0, 0.0, "a componente escrita");
    T_NEAR(T.em(id, 0, 0), 0.0, 0.0, "componente (0,0) depois de escrever em (0,1)");
    T_NEAR(T.em(id, 1, 0), 0.0, 0.0, "componente (1,0) depois de escrever em (0,1)");
    T_NEAR(T.em(id, 1, 1), 0.0, 0.0, "componente (1,1) depois de escrever em (0,1)");
    // e a componente vista como campo escalar e' a MESMA memoria
    T_NEAR(T.componente(0, 1).em(id), 7.0, 0.0, "a componente lida como campo escalar");
  }

  t_case("tensor_recusa_componente_inexistente");
  {
    CampoTensor T(m);
    bool lancou = false;
    try { T.em(m.id_global(0, 0), 2, 0); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "componente (2,0) passou sem reclamar");
  }

  // --- as facetas ------------------------------------------------------------

  // A CONTA DO MAC, calculada aqui a partir de N: num bloco uniforme N x N, as
  // facetas de normal x sao N^2 (a "menos" de cada celula) mais N (a borda
  // direita) = N(N+1).  Com dois blocos UNIDOS, a interface entre eles NAO conta
  // duas vezes: a faceta ali e' a "menos" das celulas do bloco 1, e as celulas do
  // bloco 0 nao sao donas da "mais" porque ela tem vizinho.
  t_case("a_contagem_de_facetas_segue_a_conta_do_mac");
  {
    // bloco 0: 4x4 (nivel 2); bloco 1: 2x2 (nivel 1).  A JUNCAO EM `x = 1` TEM
    // NIVEIS DIFERENTES, e e' ali que a conta discrimina.
    //
    // A FACETA PERTENCE AO LADO MAIS FINO.  Na junção o bloco 0 e' o fino, logo
    // as QUATRO celulas dele encostadas em `x = 1` possuem cada uma a sua
    // "mais" (medida 0,25), e as duas celulas grossas do bloco 1 NAO tem
    // "menos" ali.
    //
    // Facetas de normal x:
    //   bloco 0: 16 "menos" + 4 "mais" na junção = 20.
    //   bloco 1: 2 "menos" na junção? NAO -- o lado de la' e' mais fino.
    //            Sobram as 2 "menos" da coluna interna e as 2 da borda
    //            direita = 4.
    //
    // ESTE NUMERO JA' FOI 22, com a regra antiga ("a 'mais' so' existe no
    // contorno"): a junção era representada por DUAS facetas grossas em vez de
    // quatro finas.  A regra antiga fazia o conjunto de facetas de uma mesma
    // interface depender de PARA QUE LADO o refino aponta -- duas finas numa
    // orientacao, uma grossa na espelhada.
    const long esperado_x = 20 + 4;
    //   normal y: bloco 0: 16 "menos" + 4 do topo = 20;
    //             bloco 1: 4 "menos" + 2 do topo = 6.
    const long esperado_y = 20 + 6;
    T_CHECK_MSG(fd.num_facetas(0) == esperado_x,
                "facetas de normal x: %ld, esperava %ld", fd.num_facetas(0), esperado_x);
    T_CHECK_MSG(fd.num_facetas(1) == esperado_y,
                "facetas de normal y: %ld, esperava %ld", fd.num_facetas(1), esperado_y);
  }

  // POSSE UNICA: cada id de faceta aparece como propriedade de exatamente uma
  // folha, e os ids cobrem [0, num_facetas) sem furo.
  t_case("cada_faceta_tem_exatamente_um_dono");
  for (int d = 0; d < 2; d++) {
    std::set<long> vistos;
    long repetidos = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        for (int lado = 0; lado < 2; lado++) {
          const long id = fd.id_faceta(b, i, d, lado);
          if (id < 0) continue;
          if (!vistos.insert(id).second) repetidos++;
        }
    T_CHECK_MSG(repetidos == 0, "direcao %d: %ld facetas com mais de um dono", d, repetidos);
    T_CHECK_MSG(static_cast<long>(vistos.size()) == fd.num_facetas(d),
                "direcao %d: %zu facetas com dono para %ld facetas", d, vistos.size(),
                fd.num_facetas(d));
  }

  // GEOMETRIA: o centro da faceta esta' a meia celula do centro da celula, na
  // direcao da normal e so' nela; a MEDIDA e' o lado TRANSVERSAL.
  t_case("geometria_e_medida_da_faceta");
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m.centroide(b, i);
      const Ponto2 h = m.lados_da_folha(b, i);
      bool ok = true;
      for (int d = 0; d < 2 && ok; d++)
        for (int lado = 0; lado < 2 && ok; lado++) {
          const long id = fd.id_faceta(b, i, d, lado);
          if (id < 0) continue;
          const Ponto2 cf = fd.centro(d, id);
          const double esperado_n = c[d] + (lado == 0 ? -h[d] / 2.0 : h[d] / 2.0);
          if (std::fabs(cf[d] - esperado_n) > 1e-12 ||
              std::fabs(cf[1 - d] - c[1 - d]) > 1e-12) {
            T_CHECK_MSG(false,
                        "bloco %ld folha %ld dir %d lado %d: centro da faceta "
                        "(%.6f, %.6f), esperava (%.6f, %.6f)",
                        b, i, d, lado, cf[0], cf[1],
                        d == 0 ? esperado_n : c[0], d == 0 ? c[1] : esperado_n);
            ok = false;
          } else if (std::fabs(fd.medida(d, id) - h[1 - d]) > 1e-12) {
            T_CHECK_MSG(false,
                        "bloco %ld folha %ld dir %d: medida %.6f, esperava o lado "
                        "transversal %.6f (a celula e' %.3f x %.3f)",
                        b, i, d, fd.medida(d, id), h[1 - d], h[0], h[1]);
            ok = false;
          }
        }
      if (!ok) break;
    }

  // --- o campo vetorial nas facetas -----------------------------------------

  t_case("campo_vetorial_tem_um_valor_por_faceta_de_cada_direcao");
  {
    CampoVetorial u(fd);
    T_CHECK_MSG(u.tamanho(0) == fd.num_facetas(0) && u.tamanho(1) == fd.num_facetas(1),
                "campo com %ld/%ld valores para %ld/%ld facetas",
                u.tamanho(0), u.tamanho(1), fd.num_facetas(0), fd.num_facetas(1));
  }

  t_case("as_duas_direcoes_do_vetorial_sao_independentes");
  {
    CampoVetorial u(fd);
    u.preenche(0.0);
    u.em(0, 0) = 5.0;
    T_NEAR(u.em(0, 0), 5.0, 0.0, "a direcao escrita");
    T_NEAR(u.em(1, 0), 0.0, 0.0, "a direcao y depois de escrever na x");
  }

  // O ORACULO QUE AMARRA TUDO: um campo vetorial CONSTANTE tem divergencia zero
  // em toda celula.  Isso usa posse, geometria, medida e sinal de uma vez -- e a
  // celula so' sabe somar suas facetas se conseguir obter as duas pontas, inclusive
  // a "mais" que pertence ao vizinho.
  t_case("divergencia_de_campo_constante_e_zero");
  {
    CampoVetorial u(fd);
    for (int d = 0; d < 2; d++)
      for (long id = 0; id < fd.num_facetas(d); id++)
        u.em(d, id) = (d == 0) ? 2.0 : -3.0;   // constante, distinta por direcao

    double pior = 0.0;
    long pior_b = -1, pior_i = -1, internas = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        // So' celulas sem face de contorno: ali o fluxo entra e sai pelas quatro.
        bool interna = true;
        for (int fa = 0; fa < 4; fa++)
          if (m.num_vizinhos(b, i, fa) == 0) interna = false;
        if (!interna) continue;
        internas++;

        const Ponto2 h = m.lados_da_folha(b, i);
        double div = 0.0;
        for (int d = 0; d < 2; d++) {
          // faceta "menos": desta folha.  faceta "mais": do vizinho, e a folha
          // tem de alcanca-la -- e' por isso que o laco pergunta ao vizinho.
          const long id_menos = fd.id_faceta(b, i, d, 0);
          T_CHECK_MSG(id_menos >= 0, "a folha nao e' dona da propria faceta menos");

          // A FACETA DA FACE "MAIS" E' MINHA OU DO VIZINHO, conforme quem e'
          // mais fino -- nao e' sempre a "menos" do vizinho.
          //
          // Este laco ja' supos que era, e a suposicao valia enquanto a regra
          // de posse era "a 'mais' so' existe no contorno".  Na junção dos
          // dois blocos deste arquivo, que tem NIVEIS diferentes, ela e' falsa:
          // ali a celula fina possui a propria "mais" e a grossa nao tem
          // "menos".  O teste passava porque media so' a orientacao favoravel.
          long id_mais = fd.id_faceta(b, i, d, 1);
          if (id_mais < 0) {
            const Ponto2 c = m.centroide(b, i);
            Ponto2 alem = c;
            alem[d] += h[d];
            const long id_viz = m.id_do_ponto(alem);
            if (id_viz >= 0)
              for (long bb = 0; bb < m.num_blocos() && id_mais < 0; bb++)
                for (long ii = 0; ii < m.num_folhas_do_bloco(bb); ii++)
                  if (m.id_global(bb, ii) == id_viz) {
                    id_mais = fd.id_faceta(bb, ii, d, 0);
                    break;
                  }
          }
          T_CHECK_MSG(id_mais >= 0, "nao achei a faceta mais da folha %ld (dir %d)", i, d);
          if (id_mais < 0) continue;

          div += (u.em(d, id_mais) - u.em(d, id_menos)) * fd.medida(d, id_menos);
        }
        if (std::fabs(div) > pior) { pior = std::fabs(div); pior_b = b; pior_i = i; }
      }
    T_CHECK_MSG(internas > 0, "nenhuma celula interna no caso -- o oraculo nao mediu nada");
    T_CHECK_MSG(pior <= 1e-13,
                "pior divergencia %.3g (bloco %ld folha %ld) para campo constante",
                pior, pior_b, pior_i);
  }

  // DEPOIS DO REFINO 2:1: a face "mais" de uma celula grossa encostada em duas
  // finas NAO e' faceta -- as finas possuem as suas "menos", e o lado fino fica
  // com DOIS valores.  E' a consequencia da convencao, e e' o que o metodo pede.
  t_case("no_refino_2_para_1_o_lado_fino_tem_duas_facetas");
  {
    m.refina_caixa({0.30, 0.30}, {0.70, 0.70}, 4);
    DominioFacetas fd2(m);

    // Posse unica continua valendo, que e' o invariante que o refino poderia
    // quebrar.
    for (int d = 0; d < 2; d++) {
      std::set<long> vistos;
      long repetidos = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
          for (int lado = 0; lado < 2; lado++) {
            const long id = fd2.id_faceta(b, i, d, lado);
            if (id < 0) continue;
            if (!vistos.insert(id).second) repetidos++;
          }
      T_CHECK_MSG(repetidos == 0, "direcao %d depois do refino: %ld facetas com dois donos",
                  d, repetidos);
      T_CHECK_MSG(static_cast<long>(vistos.size()) == fd2.num_facetas(d),
                  "direcao %d depois do refino: %zu com dono para %ld facetas",
                  d, vistos.size(), fd2.num_facetas(d));
    }

    // E a assinatura do 2:1: existe pelo menos uma celula GROSSA cuja face "mais"
    // nao e' faceta dela (tem vizinho) e cujo vizinho fino tem DUAS facetas "menos"
    // naquele plano -- ou seja, facetas de medida MENOR que o lado da grossa.
    long grossas_com_vizinho_fino = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const int nv = m.nivel(b, i);
        for (int d = 0; d < 2; d++) {
          const std::vector<int> vz = m.niveis_vizinhos(b, i, 2 * d + 1);
          if (vz.size() >= 2) {
            // mais de um vizinho pela mesma face: o lado de la' e' mais fino
            bool todos_finos = true;
            for (int n : vz) if (n <= nv) todos_finos = false;
            if (todos_finos) grossas_com_vizinho_fino++;
          }
        }
      }
    T_CHECK_MSG(grossas_com_vizinho_fino > 0,
                "nenhuma transicao 2:1 no caso -- o caso nao exercita o que afirma");
  }

  t_case("toda_face_e_coberta_por_facetas_de_medida_exata");
  {
    // O INVARIANTE QUE FALTAVA, e e' o que a transicao 2:1 pode quebrar sem
    // quebrar a posse unica: para CADA celula e CADA lado, a soma das medidas
    // das facetas que cobrem aquela face tem de dar o lado TRANSVERSAL da
    // celula -- nem mais, nem menos.
    //
    // DOIS DEFEITOS DISTINTOS passam pela posse unica e morrem aqui:
    //
    //   FACE SEM FACETA -- soma zero.  O fluxo por ela seria descartado em
    //   silencio, e a divergencia da celula sairia errada;
    //
    //   FACETA DE MEDIDA GROSSA cobrindo duas faces finas -- soma o dobro.  A
    //   celula fina contaria o dobro da propria face.
    //
    // Era o segundo que acontecia: com a regra antiga, a celula grossa do lado
    // "mais" criava a sua "menos" com a medida dela, e as duas finas do outro
    // lado ficavam sem faceta propria.  Numa orientacao, duas facetas finas; na
    // espelhada, uma grossa.  A mesma geometria com graus de liberdade
    // diferentes conforme para que lado o refino aponta.
    m.refina_caixa({0.10, 0.10}, {0.45, 0.45}, 4);
    DominioFacetas fd3(m);

    long faces = 0, ruins = 0;
    double pior = 0.0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto2 c = m.centroide(b, i);
        const Ponto2 h = m.lados_da_folha(b, i);
        const int nivel_self = m.nivel(b, i);
        for (int d = 0; d < 2; d++) {
          const int t = 1 - d;
          for (int lado = 0; lado < 2; lado++) {
            const std::vector<int> nvz = m.niveis_vizinhos(b, i, 2 * d + lado);
            int nsub = 1;
            if (!nvz.empty() && nvz[0] > nivel_self)
              nsub = 1 << (nvz[0] - nivel_self);
            const double larg = h[t] / nsub;

            double soma = 0.0;
            for (int sub = 0; sub < nsub; sub++) {
              Ponto2 cf = c;
              cf[d] += (lado == 0 ? -1.0 : 1.0) * h[d] / 2.0;
              cf[t] += -h[t] / 2.0 + (sub + 0.5) * larg;
              // acha a faceta cujo CENTRO e' esse ponto
              for (long q = 0; q < fd3.num_facetas(d); q++) {
                const Ponto2 cq = fd3.centro(d, q);
                if (std::fabs(cq[0] - cf[0]) < 1e-12 &&
                    std::fabs(cq[1] - cf[1]) < 1e-12) {
                  soma += fd3.medida(d, q);
                  break;
                }
              }
            }
            faces++;
            const double erro = std::fabs(soma - h[t]);
            if (erro > pior) pior = erro;
            if (erro > 1e-12) ruins++;
          }
        }
      }
    T_CHECK_MSG(faces > 0, "nenhuma face percorrida -- o oraculo nao mediu nada");
    T_CHECK_MSG(ruins == 0,
                "%ld de %ld faces nao sao cobertas exatamente pelas facetas "
                "(pior desvio %.3g)", ruins, faces, pior);
  }

  return t_end();
}
