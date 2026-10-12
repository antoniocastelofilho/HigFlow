// Os campos em tres dimensoes, no arranjo deslocado (MAC).
//
// O ORACULO CENTRAL E' A CONTAGEM DE FACETAS, e ele prova a regra de
// propriedade: a faceta pertence a' celula para a qual e' a face "menos", e
// exatamente a uma.  Numa malha uniforme `N x N x N`, cada direcao tem
// `N^2 (N+1)` facetas -- `N^2` colunas de `N+1` planos.
//
// Contar e' barato e pega erro de propriedade ANTES que ele apareca como
// divergencia que nao fecha, que e' onde ele apareceria sem este caso.

#include "../testing.h"

#include <cmath>
#include <set>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "campos3d.hxx"
#include "malha3d.hxx"

using namespace t8flow3d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("a_contagem_de_facetas_segue_a_regra_de_propriedade");
  {
    // `N^2 (N+1)` POR DIRECAO.  Em duas dimensoes era `N(N+1)`; o expoente a
    // mais e' a terceira dimensao, e conferir o numero e' o que distingue a
    // regra certa de "cada celula tem suas seis facetas" (que daria `6N^3`).
    for (int nv : {1, 2, 3}) {
      const long N = 1L << nv;
      Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
      DominioFacetas fd(m);
      const long esperado = N * N * (N + 1);
      for (int d = 0; d < 3; d++)
        T_CHECK_MSG(fd.num_facetas(d) == esperado,
                    "nivel %d, direcao %d: %ld facetas, esperado N^2(N+1) = %ld",
                    nv, d, fd.num_facetas(d), esperado);
    }
  }

  t_case("cada_faceta_tem_exatamente_um_dono");
  {
    // A REGRA DITA DE OUTRO JEITO, e por isso vale como oraculo separado: se um
    // id de faceta aparecesse para duas folhas, o laco que soma fluxos contaria
    // a mesma faceta duas vezes -- e isso nao falha, da' um resultado.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);

    for (int d = 0; d < 3; d++) {
      std::set<long> vistos;
      long repetidos = 0, total = 0;
      for (long b = 0; b < m.num_blocos(); b++)
        for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
          for (int lado = 0; lado < 2; lado++) {
            const long id = fd.id_faceta(b, i, d, lado);
            if (id < 0) continue;     // nao e' dona: pertence ao vizinho
            total++;
            if (!vistos.insert(id).second) repetidos++;
          }
      T_CHECK_MSG(repetidos == 0, "direcao %d: %ld facetas com mais de um dono",
                  d, repetidos);
      T_CHECK_MSG(total == fd.num_facetas(d),
                  "direcao %d: %ld facetas com dono, mas %ld numeradas -- alguma "
                  "ficou orfa", d, total, fd.num_facetas(d));
    }
  }

  t_case("a_medida_da_faceta_e_uma_AREA_e_nao_um_comprimento");
  {
    // A TROCA QUE MAIS FACILMENTE PASSA DESPERCEBIDA numa traducao de 2D: la' a
    // medida e' `h[1-d]`, um lado so'; aqui e' o produto dos DOIS lados
    // transversais.  Um fluxo com a dimensao errada ainda converge -- para a
    // equacao errada.
    //
    // O ORACULO E' A SOMA: as facetas de uma direcao, somadas, cobrem o dominio
    // `N+1` vezes (os `N+1` planos perpendiculares aquela direcao), e cada plano
    // tem area 1.  Logo a soma das medidas e' exatamente `N+1`.
    for (int nv : {1, 2, 3}) {
      const long N = 1L << nv;
      Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, nv}}, sc_MPI_COMM_SELF);
      DominioFacetas fd(m);
      for (int d = 0; d < 3; d++) {
        double soma = 0.0;
        for (long id = 0; id < fd.num_facetas(d); id++) soma += fd.medida(d, id);
        T_NEAR(soma, static_cast<double>(N + 1), 1e-12,
               "soma das areas das facetas de uma direcao");
      }
    }

    // E O CONTRASTE COM A TRADUCAO ERRADA: se a medida fosse um COMPRIMENTO
    // (`h` em vez de `h*h`), a soma daria `(N+1)*N` em vez de `N+1` -- um fator
    // `N`, que cresce com o refino.  Afirmar a soma certa e' o que separa os
    // dois.
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 3}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    const double h = 1.0 / 8.0;
    T_NEAR(fd.medida(0, 0), h * h, 1e-15,
           "a medida de uma faceta e' o produto dos dois lados transversais");
  }

  t_case("o_centro_da_faceta_cai_na_face_da_celula_dona");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);

    long erros = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        const Ponto3 hh = m.lados_da_folha(b, i);
        for (int d = 0; d < 3; d++)
          for (int lado = 0; lado < 2; lado++) {
            const long id = fd.id_faceta(b, i, d, lado);
            if (id < 0) continue;
            const Ponto3 cf = fd.centro(d, id);
            // O centro da faceta difere do centro da celula por meio lado NA
            // DIRECAO `d`, e coincide nas outras duas.
            const double esperado_d = c[d] + (lado == 0 ? -0.5 : 0.5) * hh[d];
            if (std::fabs(cf[d] - esperado_d) > 1e-12) erros++;
            for (int t = 0; t < 3; t++)
              if (t != d && std::fabs(cf[t] - c[t]) > 1e-12) erros++;
          }
      }
    T_CHECK_MSG(erros == 0, "%ld centros de faceta fora da face da celula dona",
                erros);
  }

  t_case("o_tensor_tem_nove_componentes_e_o_indice_recusa_o_que_nao_existe");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1}}, sc_MPI_COMM_SELF);
    CampoTensor t(m);

    // NOVE COMPONENTES INDEPENDENTES: escrever numa nao pode mexer noutra.  O
    // teste escreve um valor distinto em cada e confere todos -- um indice
    // errado (por exemplo `2*lin+col`, que era o de duas dimensoes) faria duas
    // componentes colidirem.
    const long id = m.id_global(0, 0);
    for (int l = 0; l < 3; l++)
      for (int c = 0; c < 3; c++) t.em(id, l, c) = 10.0 * l + c;
    for (int l = 0; l < 3; l++)
      for (int c = 0; c < 3; c++)
        T_NEAR(t.em(id, l, c), 10.0 * l + c, 1e-15,
               "componente do tensor depois de escrever as nove");

    bool lancou = false;
    try { t.em(id, 3, 0); } catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "a linha 3 passou num tensor 3x3");
  }

  t_case("o_campo_vetorial_tem_uma_componente_por_direcao_de_faceta");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    DominioFacetas fd(m);
    CampoVetorial u(fd);

    for (int d = 0; d < 3; d++)
      T_CHECK_MSG(u.tamanho(d) == fd.num_facetas(d),
                  "direcao %d: o campo tem %ld valores para %ld facetas", d,
                  u.tamanho(d), fd.num_facetas(d));

    // AS TRES COMPONENTES SAO INDEPENDENTES: preencher uma nao pode tocar as
    // outras.  Em duas dimensoes eram duas, e um laco que nao acompanhe a
    // dimensao deixaria `z` sem inicializar -- em silencio.
    u.preenche(0.0);
    for (long id = 0; id < fd.num_facetas(1); id++) u.em(1, id) = 7.0;
    for (int d = 0; d < 3; d++)
      for (long id = 0; id < fd.num_facetas(d); id++)
        T_NEAR(u.em(d, id), (d == 1) ? 7.0 : 0.0, 1e-15,
               "componente do campo vetorial depois de escrever so' em y");
  }

  t_case("o_campo_escalar_le_pelo_ponto_e_pelo_id");
  {
    Malha3 m({{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 2}}, sc_MPI_COMM_SELF);
    CampoEscalar p(m);

    // DOIS CAMINHOS PARA O MESMO VALOR -- por id e pela geometria --, e e' isso
    // que faz deste um oraculo e nao uma tautologia.
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        p.em(m.id_global(b, i)) = c[0] + 10.0 * c[1] + 100.0 * c[2];
      }

    long erros = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto3 c = m.centroide(b, i);
        if (std::fabs(p.em_ponto(c) - p.em(m.id_global(b, i))) > 1e-15) erros++;
      }
    T_CHECK_MSG(erros == 0, "%ld celulas onde a leitura por ponto difere da "
                "leitura por id", erros);

    bool lancou = false;
    try { p.em_ponto({5.0, 0.5, 0.5}); }
    catch (const std::out_of_range &) { lancou = true; }
    T_CHECK_MSG(lancou, "um ponto fora do dominio foi lido");
  }

  return t_end();
}
