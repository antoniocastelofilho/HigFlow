// A malha em TRES dimensoes: floresta de hexaedros, multibloco, com franja.
//
// OS ORACULOS AQUI NAO SAO OS DE DUAS DIMENSOES TRADUZIDOS.  Eles sao os que o
// 2D ENSINOU a ter -- cada um existe porque um defeito concreto passou por ele
// --, mais os que so' tres dimensoes pede:
//
//   a ordem dos OITO vertices do hexaedro, que nao tem a intuicao de "sentido
//   anti-horario" que salva em duas dimensoes, e cujo erro COMMITA sem reclamar;
//
//   a transicao 2:1, onde uma face e' partida em QUATRO e nao em duas;
//
//   o maximo de vizinhos, que a franja usa para se dimensionar -- e que foi
//   MEDIDO antes deste arquivo existir, com o resultado mudando a arquitetura.

#include "../testing.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "arranque.hxx"
#include "malha3d.hxx"

using namespace t8flow3d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("um_bloco_uniforme_tem_8_elevado_ao_nivel_folhas");
  {
    // EM TRES DIMENSOES O REFINO E' OITO PARA UM, e a contagem e' a conferencia
    // mais barata de que a classe de elemento esta' certa: um QUAD daria 4^n.
    for (int nv : {1, 2, 3}) {
      Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
      long esperado = 1;
      for (int k = 0; k < nv; k++) esperado *= 8;
      T_CHECK_MSG(m.num_folhas_globais() == esperado,
                  "nivel %d devia dar %ld folhas, deu %ld", nv, esperado,
                  m.num_folhas_globais());
    }
  }

  t_case("o_centroide_de_cada_folha_cai_dentro_da_caixa_do_bloco");
  {
    // O ORACULO DA ORDEM DOS VERTICES, e em tres dimensoes ele importa mais.
    //
    // O `cmesh` espera os oito vertices na ordem de bits -- `x` mais rapido,
    // depois `y`, depois `z` --, e trocar dois produz um hexaedro TORCIDO que
    // commita sem erro nenhum.  Em duas dimensoes a intuicao de "sentido
    // anti-horario" ajuda a acertar; com oito vertices ela nao se aplica, e o
    // unico recurso e' verificar a geometria resultante.
    //
    // DOIS BLOCOS, e o segundo DESLOCADO: com um bloco so' na origem, uma troca
    // de vertices que espelhe o cubo em torno do centro passaria -- o centroide
    // continuaria dentro da caixa.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2},
              {{1.0, 0.0, 0.0}, {3.0, 1.0, 1.0}, 1}}, sc_MPI_COMM_SELF);

    long fora = 0;
    double soma_vol = 0.0;
    for (long b = 0; b < m.num_blocos(); b++) {
      const Bloco &k = m.bloco(b);
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        for (int d = 0; d < 3; d++)
          if (c[d] < k.lo[d] - 1e-12 || c[d] > k.hi[d] + 1e-12) fora++;
        soma_vol += m.volume(b, i);
      }
    }
    T_CHECK_MSG(fora == 0, "%ld coordenadas de centroide fora da caixa do bloco", fora);

    // E A CONSERVACAO DE VOLUME, que e' independente e pega o que o centroide
    // nao pega: um hexaedro torcido pode ter o centroide no lugar certo e volume
    // errado.  Os dois blocos somam 1 + 2 = 3.
    T_NEAR(soma_vol, 3.0, 1e-12, "volume total dos dois blocos");
    T_CHECK_MSG(m.num_juncoes() == 1, "os dois blocos encostam numa face, e a "
                "malha contou %ld juncoes", m.num_juncoes());
  }

  t_case("o_hex_tem_seis_faces_e_a_transicao_parte_cada_uma_em_quatro");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, 3);

    // SEIS FACES, e pedir a setima e' recusado.
    bool lancou = false;
    try { m.ids_vizinhos(0, 0, 6); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "a face 6 passou num hexaedro, que tem 6 faces (0 a 5)");

    // NA TRANSICAO 2:1 UMA FACE TEM QUATRO VIZINHOS, nao dois: a face de um
    // hexaedro e' um quadrado, e refinar o vizinho a parte em quatro.  Afirmar o
    // numero e' o que distingue HEX de QUAD tratado como HEX.
    long com_quatro = 0, maior = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
        for (int fa = 0; fa < 6; fa++) {
          const long n = static_cast<long>(m.ids_vizinhos(b, i, fa).size());
          maior = std::max(maior, n);
          if (n == 4) com_quatro++;
        }
    T_CHECK_MSG(maior == 4,
                "o maior numero de vizinhos numa face deu %ld; em HEX com "
                "transicao 2:1 tem de ser 4", maior);
    T_CHECK_MSG(com_quatro > 0, "nenhuma face com quatro vizinhos -- a malha nao "
                "tem a transicao que este caso precisa exercitar");

    // E O MAXIMO POR FOLHA, que a franja usa para se dimensionar.  O teto
    // teorico e' 24 (seis faces vezes quatro) e o medido e' muito menor, porque
    // uma celula tem no maximo uma face com vizinhos refinados em cada direcao.
    T_CHECK_MSG(m.max_vizinhos_por_folha() >= 6 && m.max_vizinhos_por_folha() < 24,
                "max vizinhos por folha deu %d; esperado entre 6 e 23 -- 24 seria "
                "o teto teorico, que a medida mostrou nunca ser atingido",
                m.max_vizinhos_por_folha());
    printf("  (transicao: %ld faces com 4 vizinhos, max por folha %d)\n",
           com_quatro, m.max_vizinhos_por_folha());
  }

  t_case("a_razao_entre_vizinhos_e_no_maximo_2_para_1");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.0, 0.0, 0.0}, {0.4, 0.4, 0.4}, 4);

    long violacoes = 0;
    int pior = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const int nv = m.nivel(b, i);
        for (int fa = 0; fa < 6; fa++)
          for (int viz : m.niveis_vizinhos(b, i, fa)) {
            const int dif = std::abs(viz - nv);
            pior = std::max(pior, dif);
            if (dif > 1) violacoes++;
          }
      }
    T_CHECK_MSG(violacoes == 0,
                "%ld pares de vizinhos com diferenca de nivel maior que 1 (pior: "
                "%d).  O balanceamento 2:1 e' o que os esquemas pressupoem",
                violacoes, pior);
    T_CHECK_MSG(pior == 1, "nenhuma transicao na malha (pior diferenca %d) -- o "
                "caso nao mediu nada", pior);
  }

  t_case("a_busca_por_ponto_devolve_a_folha_que_contem_o_ponto");
  {
    // O ORACULO DO INDICE ESPACIAL, e ele e' independente: `id_do_ponto` e
    // `id_global` sao dois caminhos de codigo distintos, e o centro de uma folha
    // tem de levar ao id dela pelos dois.
    //
    // Em duas dimensoes este caso pegou uma chave de indice que esquecia o
    // BLOCO, fazendo a celula (0,0) do bloco 1 colidir com a do bloco 0.  Com
    // tres coordenadas ha' mais bits a repartir, e mais chance de errar.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2},
              {{1.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    m.refina_caixa({0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, 3);

    long erros = 0, nao_achados = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const long id = m.id_do_ponto(m.centroide(b, i));
        if (id < 0) nao_achados++;
        else if (id != m.id_global(b, i)) erros++;
      }
    T_CHECK_MSG(nao_achados == 0, "%ld centroides nao foram achados pela busca",
                nao_achados);
    T_CHECK_MSG(erros == 0, "%ld folhas cujo id pela busca difere do id pelo "
                "percurso -- o indice espacial esta' colidindo", erros);

    // FORA DO DOMINIO devolve -1, e nao um id plausivel.
    T_CHECK_MSG(m.id_do_ponto({5.0, 0.5, 0.5}) < 0,
                "um ponto fora do dominio recebeu um id");
  }

  t_case("arvore_de_uma_celula_e_recusada");
  {
    // UMA CELULA NAO TEM INTERIOR: as SEIS facetas dela sao de fronteira.  O
    // levantamento do t8code registra o preco de deixar passar -- "so' apareceu
    // como escoamento divergindo numa simulacao de meia hora".
    bool lancou = false;
    std::string msg;
    try { Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 0}}, sc_MPI_COMM_SELF); }
    catch (const std::invalid_argument &e) { lancou = true; msg = e.what(); }
    T_CHECK_MSG(lancou, "bloco de nivel 0 -- uma celula so' -- passou");
    T_CHECK_MSG(msg.find("interior") != std::string::npos,
                "a mensagem nao explica por que uma celula nao serve: \"%s\"",
                msg.c_str());

    // CONTRASTE: nivel 1 passa, e da' oito celulas.  Sem ele, a guarda poderia
    // estar recusando tudo.
    Malha3 ok({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1}}, sc_MPI_COMM_SELF);
    T_CHECK_MSG(ok.num_folhas_globais() == 8, "nivel 1 devia dar 8 folhas, deu %ld",
                ok.num_folhas_globais());
  }

  t_case("blocos_que_se_sobrepoem_ou_encostam_parcialmente_sao_recusados");
  {
    // SOBREPOSICAO EM VOLUME: a busca por ponto dependeria da ordem dos blocos.
    bool lancou = false;
    try {
      Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1},
                {{0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, 1}}, sc_MPI_COMM_SELF);
    } catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "blocos sobrepostos em volume passaram");

    // ENCOSTE PARCIAL: em tres dimensoes ha' mais maneiras de acontecer que em
    // duas -- as DUAS direcoes transversais podem casar uma e nao a outra --, e
    // o `join` do t8code nao representa isso.  Aceitar daria contorno no meio do
    // dominio.
    bool lancou2 = false;
    try {
      Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1},
                {{1.0, 0.0, 0.0}, {2.0, 0.5, 1.0}, 1}}, sc_MPI_COMM_SELF);
    } catch (const std::invalid_argument &) { lancou2 = true; }
    T_CHECK_MSG(lancou2, "blocos que encostam com faces de tamanhos diferentes "
                "passaram -- o join do t8code nao representa isso");
  }

  t_case("o_refino_que_nao_refina_nada_e_recusado");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);

    // CAIXA QUE NAO ALCANCA bloco algum.
    bool fora = false;
    try { m.refina_caixa({5.0, 5.0, 5.0}, {6.0, 6.0, 6.0}, 3); }
    catch (const std::invalid_argument &) { fora = true; }
    T_CHECK_MSG(fora, "caixa fora do dominio foi aceita");

    // NIVEL QUE NAO APROFUNDA.  Aceito calado, deixaria a malha PARECER
    // refinada -- e quem olhasse a contagem depois acharia que o refino falhou
    // em vez de nao ter sido pedido.
    bool raso = false;
    try { m.refina_caixa({0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, 2); }
    catch (const std::invalid_argument &) { raso = true; }
    T_CHECK_MSG(raso, "refino para o nivel que a malha ja' tem foi aceito");
  }

  t_case("o_refino_aleatorio_nao_depende_da_ordem_de_percurso");
  {
    // A MESMA SEMENTE DA' A MESMA MALHA, e a decisao vem do HASH DO CENTROIDE e
    // nao de um gerador sequencial.
    //
    // A RAZAO E' O PARALELO: um gerador sequencial produz resultado que depende
    // da ordem de visita, logo da particao -- a mesma malha refinaria diferente
    // em np=2 e np=4, e nenhuma comparacao entre eles valeria.
    Malha3 a({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    Malha3 b({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    a.refina_aleatorio(0.4, 2026, 3);
    b.refina_aleatorio(0.4, 2026, 3);
    T_CHECK_MSG(a.num_folhas_globais() == b.num_folhas_globais(),
                "a mesma semente deu %ld e %ld folhas", a.num_folhas_globais(),
                b.num_folhas_globais());

    // E SEMENTE DIFERENTE DA' MALHA DIFERENTE -- o contraste, sem o qual o caso
    // acima passaria com um refino que ignora a semente.
    Malha3 c({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    c.refina_aleatorio(0.4, 777, 3);
    T_CHECK_MSG(c.num_folhas_globais() != a.num_folhas_globais(),
                "sementes diferentes deram a MESMA contagem (%ld) -- o refino "
                "pode estar ignorando a semente", a.num_folhas_globais());

    // A FRACAO E' RESPEITADA, dentro do ruido do sorteio: 40% de 64 celulas
    // refinadas em oito cada da' cerca de 64 + 0,4*64*7 = 243.
    const long base = 64;
    const double esperado = base + 0.4 * base * 7.0;
    const double obtido = static_cast<double>(a.num_folhas_globais());
    T_CHECK_MSG(obtido > 0.6 * esperado && obtido < 1.4 * esperado,
                "%ld folhas depois de refinar 40%%; esperava cerca de %.0f",
                a.num_folhas_globais(), esperado);
    printf("  (aleatorio: %ld folhas, esperado ~%.0f)\n", a.num_folhas_globais(),
           esperado);
  }

  return t_end();
}
