// A FLORESTA DO T8CODE COMO ESTRUTURA, em 2D, alimentando o discretizador 2D.
//
// Este caso afirma a conviviencia que nenhuma das metades prova sozinha: o
// nucleo herdado (C da HiGTree, compilado como C++20) e o t8code na MESMA
// unidade de traducao, com o SUPORTE saindo da floresta e os PESOS saindo do MLS.
// Se nao convivessem, isto nao linkaria.
//
// Em 2D o elemento e' `T8_ECLASS_QUAD`, 4 vertices e 4 faces -- escrito
// explicitamente, sem `#if DIM`, que e' a decisao 7.
#include "../testing.h"

#include <array>
#include <vector>

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include "arranque.hxx"
#include "mls2d.hxx"

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

int main(int argc, char **argv)
{
  t8flow::Arranque arranque(&argc, &argv);

  const int NIVEL = 3;   // 8 por direcao

  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);   // obrigatorio ANTES do gerador: sem ele o segfault
                           // sai tres quadros abaixo, em t8_cmesh_set_tree_class
  t8_cmesh_new_hypercube(&cmesh, T8_ECLASS_QUAD, sc_MPI_COMM_SELF, 0, 0, 0);
  const t8_scheme_c *esquema = t8_scheme_new_default();
  t8_forest_t f = t8_forest_new_uniform(cmesh, esquema, NIVEL, 0, sc_MPI_COMM_SELF);

  t_case("floresta_uniforme_tem_as_folhas_do_nivel");
  const long esperadas = (1L << NIVEL) * (1L << NIVEL);
  const long folhas = static_cast<long>(t8_forest_get_local_num_leaf_elements(f));
  T_CHECK_MSG(folhas == esperadas, "nivel %d em 2D: esperava %ld folhas, obtive %ld",
              NIVEL, esperadas, folhas);

  // Os centroides das folhas SAO o suporte: e' o contrato malha/discretizacao --
  // a malha diz quais celulas ficam perto, a discretizacao faz o ajuste.
  std::vector<Ponto2> centros;
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f, it, ie);
      double c[3] = {0.0, 0.0, 0.0};
      t8_forest_element_centroid(f, it, e, c);
      centros.push_back({c[0], c[1]});
    }
  }

  t_case("suporte_colhido_da_floresta_tem_uma_entrada_por_folha");
  T_CHECK_MSG(static_cast<long>(centros.size()) == folhas,
              "colhi %ld centroides para %ld folhas", static_cast<long>(centros.size()), folhas);

  // SUPORTE LOCAL, nao a floresta inteira: busca sem parada devolve a malha toda,
  // o ajuste reproduz mesmo assim, e o teste passa enquanto a topologia deixa de
  // importar.  Ponto de consulta assimetrico, pelo motivo do test-mls2d.
  const int ORDEM = 2;
  const int MAXPTS = 4 * Mls::min_pontos(ORDEM);
  const Ponto2 alvo = {0.41, 0.56};

  std::vector<int> ordem(centros.size());
  for (size_t i = 0; i < centros.size(); i++) ordem[i] = static_cast<int>(i);
  const auto d2 = [&](int k) {
    const double dx = centros[k][0] - alvo[0], dy = centros[k][1] - alvo[1];
    return dx * dx + dy * dy;
  };
  for (int k = 0; k < MAXPTS && k < static_cast<int>(ordem.size()); k++) {
    int melhor = k;
    for (int i = k + 1; i < static_cast<int>(ordem.size()); i++)
      if (d2(ordem[i]) < d2(ordem[melhor])) melhor = i;
    const int t = ordem[k]; ordem[k] = ordem[melhor]; ordem[melhor] = t;
  }

  const int numpts = (MAXPTS < static_cast<int>(centros.size()))
                         ? MAXPTS : static_cast<int>(centros.size());
  std::vector<Ponto2> suporte;
  std::vector<double> v2, v3;
  for (int i = 0; i < numpts; i++) {
    suporte.push_back(centros[ordem[i]]);
    v2.push_back(grau2(suporte.back()));
    v3.push_back(grau3(suporte.back()));
  }

  Mls mls(ORDEM, numpts);

  t_case("mls_reproduz_polinomio_sobre_suporte_da_floresta");
  const double obtido2 = mls.interpola(alvo, suporte, v2);
  const double exato2 = grau2(alvo);
  T_NEAR(obtido2, exato2, 1e-9, "campo de grau 2 nos centroides das folhas");

  t_case("suporte_da_floresta_erra_mais_no_grau_acima");
  const double obtido3 = mls.interpola(alvo, suporte, v3);
  const double err2 = std::fabs(obtido2 - exato2);
  const double err3 = std::fabs(obtido3 - grau3(alvo));
  T_CONTRASTE(err2 > 1e-14 ? err2 : 1e-14, err3, 1e3,
              "grau 2 contra grau 3 sobre o suporte da floresta");

  t8_forest_unref(&f);
  return t_end();
}
