// O CRITERIO DE REFINO ADAPTATIVO.
//
// O ORACULO CENTRAL DESTE ARQUIVO E' A IGUALDADE BIT A BIT entre a tabela por
// DILATACAO -- o caminho de producao, `O(nsem . R^2)` -- e a por BUSCA
// EXAUSTIVA, `O(ncel . nsem)`.  Nao e' redundancia: a dilatacao so' e' exata
// porque o nivel e' monotono decrescente na distancia, e sem a segunda
// implementacao ao lado essa propriedade seria uma crenca.
//
// E IGUALDADE EXATA, nao "proxima": as duas calculam o MESMO numero por
// caminhos diferentes, e qualquer diferenca e' defeito de uma delas -- nao
// erro de arredondamento tolerado.

#include "../testing.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "arranque.hxx"
#include "criterio2d.hxx"
#include "malha2d.hxx"

using namespace t8flow2d;

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  t_case("os_limiares_saem_decrescentes_e_cobrem_as_celulas_minimas");
  {
    const double h = 0.125;
    const std::vector<double> thr = amr_limiares(h, 3, 4);
    T_CHECK_MSG(thr.size() == 3, "deviam sair 3 limiares, sairam %zu", thr.size());
    T_CHECK_MSG(thr[0] > thr[1] && thr[1] > thr[2],
                "os limiares nao sao decrescentes: %.4g %.4g %.4g",
                thr[0], thr[1], thr[2]);

    // O MAIS FINO cobre `cel_min` celulas finas mais a margem de uma celula
    // base.  Afirmar o numero aqui e' o que impede que alguem o "ajuste".
    const double h_fino = h / 8.0;
    T_NEAR(thr[2], 4 * h_fino + h + h_fino, 1e-15,
           "o limiar do nivel mais fino");

    // E NENHUM e' aceito com cel_min zero: banda de zero celulas nao refina
    // nada, e aceitar calado deixaria a malha PARECER refinada.
    bool lancou = false;
    try { amr_limiares(h, 3, 0); }
    catch (const std::invalid_argument &) { lancou = true; }
    T_CHECK_MSG(lancou, "cel_min = 0 foi aceito");
  }

  t_case("a_dilatacao_e_a_exaustiva_dao_a_mesma_tabela");
  {
    // O CASO E' DESENHADO PARA DISCRIMINAR: sementes espalhadas, grade nao
    // quadrada (hx != hy), e bandas que se sobrepoem -- e' na sobreposicao que
    // o `max` da dilatacao tem de dar o mesmo que o `min` de distancia da
    // exaustiva.
    const Ponto2 lo{0.0, 0.0};
    const int nx = 48, ny = 32;
    const double hx = 1.0 / nx, hy = 0.75 / ny;
    const std::vector<double> thr = amr_limiares(std::fmin(hx, hy), 3, 3);

    std::vector<Ponto2> sem;
    for (int k = 0; k < 17; k++) {
      const double t = k / 17.0;
      sem.push_back({0.1 + 0.8 * t, 0.15 + 0.45 * std::sin(6.0 * t)});
    }

    const std::vector<signed char> a = amr_tabela(sem, lo, hx, hy, nx, ny, thr);
    const std::vector<signed char> b =
        amr_tabela_exaustiva(sem, lo, hx, hy, nx, ny, thr);

    long difs = 0, refinadas = 0;
    for (size_t q = 0; q < a.size(); q++) {
      if (a[q] != b[q]) difs++;
      if (a[q] > 0) refinadas++;
    }
    T_CHECK_MSG(difs == 0,
                "%ld de %zu celulas diferem entre a dilatacao e a exaustiva",
                difs, a.size());

    // E O CASO TEM DE REFINAR ALGUMA COISA, sem refinar tudo.  Duas tabelas
    // todas zero tambem sao iguais, e nao verificam nada.
    T_CHECK_MSG(refinadas > 0, "nenhuma celula refinada -- as duas tabelas sao zero");
    T_CHECK_MSG(refinadas < static_cast<long>(a.size()),
                "TODAS as %zu celulas refinadas -- o criterio nao discrimina",
                a.size());

    // OS TRES NIVEIS TEM DE APARECER: com so' um nivel presente, a monotonia
    // da distancia nunca e' exercitada.
    bool visto[4] = {false, false, false, false};
    for (signed char v : a) if (v >= 0 && v <= 3) visto[v] = true;
    T_CHECK_MSG(visto[1] && visto[2] && visto[3],
                "nem todos os niveis aparecem (1:%d 2:%d 3:%d) -- a monotonia "
                "nao foi exercitada", visto[1], visto[2], visto[3]);
  }

  t_case("a_malha_segue_a_tabela");
  {
    // DA TABELA A' MALHA: o nivel que o criterio pede tem de ser o nivel que a
    // folha tem -- ou MAIS, porque o equilibrio 2:1 pode subir vizinhas.
    const int base = 3;
    Malha m({{{0.0, 0.0}, {1.0, 1.0}, base}}, sc_MPI_COMM_WORLD);
    const int nx = 1 << base, ny = 1 << base;
    const double h = 1.0 / nx;
    const std::vector<double> thr = amr_limiares(h, 1, 1);
    const std::vector<Ponto2> sem = {{0.5, 0.5}};
    const std::vector<signed char> tab =
        amr_tabela(sem, {0.0, 0.0}, h, h, nx, ny, thr);

    const long antes = m.num_folhas_globais();
    m.refina_por_criterio([&](const Ponto2 &x) {
      return base + amr_nivel_em(tab, {0.0, 0.0}, h, h, nx, ny, x);
    });
    const long depois = m.num_folhas_globais();

    T_CHECK_MSG(depois > antes,
                "a malha nao refinou: %ld folhas antes, %ld depois", antes, depois);
    T_CHECK_MSG(m.num_folhas_em_transicao() > 0,
                "refino sem transicao 2:1 -- ou refinou tudo, ou nada");

    // CADA FOLHA TEM AO MENOS o nivel que a tabela pede no centro dela.  "Ao
    // menos" porque o equilibrio 2:1 sobe vizinhas que a tabela nao pediu --
    // e exigir IGUALDADE aqui falharia por causa do equilibrio, nao do
    // criterio.
    long faltou = 0;
    for (long b = 0; b < m.num_blocos(); b++)
      for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
        const Ponto2 c = m.centroide(b, i);
        const int pedido = base + amr_nivel_em(tab, {0.0, 0.0}, h, h, nx, ny, c);
        if (m.nivel(b, i) < pedido) faltou++;
      }
    T_CHECK_MSG(faltou == 0, "%ld folhas com nivel ABAIXO do que a tabela pede",
                faltou);
  }

  return t_end();
}
