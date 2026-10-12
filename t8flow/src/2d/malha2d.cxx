#include "malha2d.hxx"

#include <cmath>
#include <stdexcept>
#include <string>

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_schemes/t8_default/t8_default.hxx>

namespace t8flow2d {

namespace {

constexpr double EPS = 1e-12;

bool igual(double a, double b) { return std::fabs(a - b) <= EPS * (1.0 + std::fabs(a) + std::fabs(b)); }

//! Os dois intervalos se sobrepoem em MEDIDA (toque de borda nao conta).
bool sobrepoe(double a0, double a1, double b0, double b1)
{
  return (a0 < b1 - EPS) && (b0 < a1 - EPS);
}

//! O nivel alvo de cada arvore, visto pelo callback de refino.
//
// Vai pelo `user_data` da floresta porque o callback do t8code nao aceita
// contexto proprio -- e vai por ponteiro para vetor, nao por estatico de arquivo,
// que e' o que permitiria duas malhas vivas ao mesmo tempo.
struct Alvos {
  const std::vector<int> *nivel;
};

int adapt_por_arvore(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                     const t8_eclass_t tree_class, t8_locidx_t /*lelement_id*/,
                     const t8_scheme_c *scheme, const int /*is_family*/,
                     const int /*num_elements*/, t8_element_t *elements[])
{
  const Alvos *a = static_cast<const Alvos *>(t8_forest_get_user_data(forest));
  if (a == nullptr || a->nivel == nullptr) return 0;

  // `which_tree` e' o indice da arvore LOCAL, e em serie ele e' o indice do
  // bloco.  Sob particionamento isso deixa de valer e a traducao passa a ser
  // explicita -- anotado porque e' exatamente a classe de erro que este projeto
  // decidiu nao herdar.
  const long b = static_cast<long>(t8_forest_get_first_local_tree_id(forest_from)) + which_tree;
  if (b < 0 || b >= static_cast<long>(a->nivel->size())) return 0;

  const int nivel = scheme->element_get_level(tree_class, elements[0]);
  return (nivel < (*a->nivel)[static_cast<size_t>(b)]) ? 1 : 0;
}

}  // namespace

Malha::Malha(const std::vector<Bloco> &blocos, sc_MPI_Comm comm)
    : blocos_(blocos)
{
  const long n = static_cast<long>(blocos.size());
  if (n == 0) throw std::invalid_argument("Malha: lista de blocos vazia");

  for (long b = 0; b < n; b++) {
    const Bloco &k = blocos[static_cast<size_t>(b)];
    for (int d = 0; d < 2; d++)
      if (!(k.hi[d] > k.lo[d] + EPS))
        throw std::invalid_argument(
            "Malha: bloco " + std::to_string(b) + " degenerado na direcao " +
            std::to_string(d) + ": lo=" + std::to_string(k.lo[d]) + " hi=" + std::to_string(k.hi[d]));
    if (k.nivel < 0)
      throw std::invalid_argument("Malha: bloco " + std::to_string(b) +
                                  " com nivel negativo (" + std::to_string(k.nivel) + ")");
  }

  // SOBREPOSICAO EM AREA e' recusada.  O HiGFlow assume dominios que nao se
  // sobrepoem (o `cell_find_in_center` para no primeiro bloco que contem o ponto,
  // com o comentario "assume domains do not overlap"), e aceitar sobreposicao
  // aqui faria a resposta depender da ORDEM dos blocos.
  for (long b = 0; b < n; b++)
    for (long c = b + 1; c < n; c++) {
      const Bloco &p = blocos[static_cast<size_t>(b)];
      const Bloco &q = blocos[static_cast<size_t>(c)];
      if (sobrepoe(p.lo[0], p.hi[0], q.lo[0], q.hi[0]) &&
          sobrepoe(p.lo[1], p.hi[1], q.lo[1], q.hi[1]))
        throw std::invalid_argument("Malha: blocos " + std::to_string(b) + " e " +
                                    std::to_string(c) + " se sobrepoem em area");
    }

  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);   // obrigatorio ANTES de qualquer set_: sem ele o
                           // segfault sai dentro de t8_cmesh_set_tree_class
  // A GEOMETRIA LINEAR, pela API C -- e a escolha nao e' de estilo.  A forma C++
  // (`t8_cmesh_register_geometry<t8_geometry_linear>`) vive em
  // `t8_cmesh/t8_cmesh.hxx`, que o PREFIXO INSTALADO expoe mas NAO e' compilavel
  // de fora: ele inclui `t8_cmesh_internal/t8_cmesh_types.h`, que o install nao
  // leva.  A sobrecarga C (`t8_cmesh_register_geometry(cmesh, t8_geometry_c *)`
  // com `t8_geometry_linear_new()`) esta' no `t8_cmesh.h` publico e o cmesh toma
  // posse do objeto.
  t8_cmesh_register_geometry(cmesh, t8_geometry_linear_new());

  for (long b = 0; b < n; b++) {
    const Bloco &k = blocos[static_cast<size_t>(b)];
    t8_cmesh_set_tree_class(cmesh, b, T8_ECLASS_QUAD);

    // Z-ORDER, x variando mais rapido -- a convencao do t8code, lida do
    // `vertices_coords` dos exemplos.  Trocar v1 com v2 produz uma malha que
    // COMMITA e tem geometria torcida, sem erro nenhum: e' por isso que existe o
    // caso "centroide cai na caixa do seu bloco".
    const double v[12] = {
        k.lo[0], k.lo[1], 0.0,
        k.hi[0], k.lo[1], 0.0,
        k.lo[0], k.hi[1], 0.0,
        k.hi[0], k.hi[1], 0.0,
    };
    t8_cmesh_set_tree_vertices(cmesh, b, v, 4);
  }

  // AS JUNCOES.  Faces do QUAD: 0:x-, 1:x+, 2:y-, 3:y+.  Duas caixas vizinhas
  // compartilham face quando um limite coincide e o intervalo transversal e' o
  // MESMO; coincidencia parcial e' recusada, porque o join do t8code nao a
  // representa e o silencio ali daria contorno no meio do dominio.
  for (long b = 0; b < n; b++)
    for (long c = b + 1; c < n; c++) {
      const Bloco &p = blocos[static_cast<size_t>(b)];
      const Bloco &q = blocos[static_cast<size_t>(c)];

      for (int d = 0; d < 2; d++) {
        const int t = 1 - d;                 // direcao transversal
        const int face_mais = 2 * d + 1;     // d+ de quem esta' antes
        const int face_menos = 2 * d;        // d- de quem esta' depois

        // p antes de q, ou q antes de p
        for (int sentido = 0; sentido < 2; sentido++) {
          const Bloco &a = (sentido == 0) ? p : q;
          const Bloco &z = (sentido == 0) ? q : p;
          const long ia = (sentido == 0) ? b : c;
          const long iz = (sentido == 0) ? c : b;
          if (!igual(a.hi[d], z.lo[d])) continue;

          const bool coincide = igual(a.lo[t], z.lo[t]) && igual(a.hi[t], z.hi[t]);
          if (coincide) {
            // Orientacao 0: os blocos estao alinhados aos eixos e varrem a face
            // na mesma direcao, logo os vertices casam em ordem.
            t8_cmesh_set_join(cmesh, ia, iz, face_mais, face_menos, 0);
            juncoes_++;
          } else if (sobrepoe(a.lo[t], a.hi[t], z.lo[t], z.hi[t])) {
            t8_cmesh_destroy(&cmesh);
            throw std::invalid_argument(
                "Malha: blocos " + std::to_string(ia) + " e " + std::to_string(iz) +
                " encostam na direcao " + std::to_string(d) +
                " com faces que se sobrepoem PARCIALMENTE; o join do t8code nao"
                " representa isso, e aceitar daria contorno no meio do dominio");
          }
        }
      }
    }

  t8_cmesh_commit(cmesh, comm);

  // NIVEL POR BLOCO.  A floresta nasce no menor nivel pedido e sobe por adapt
  // RECURSIVO, uma arvore de cada vez conforme o alvo dela.
  //
  // E o refino vem em commit SEPARADO de qualquer ghost/particao: pedir
  // `set_adapt`, `set_ghost` e `set_partition` no mesmo commit do t8code derruba
  // com np>1 e FUNCIONA em np=1 -- armadilha medida, e das piores, porque o teste
  // serial nao a ve.
  int nivel_min = blocos[0].nivel;
  std::vector<int> alvo(static_cast<size_t>(n));
  for (long b = 0; b < n; b++) {
    alvo[static_cast<size_t>(b)] = blocos[static_cast<size_t>(b)].nivel;
    if (blocos[static_cast<size_t>(b)].nivel < nivel_min) nivel_min = blocos[static_cast<size_t>(b)].nivel;
  }

  const t8_scheme_c *esquema = t8_scheme_new_default();
  f_ = t8_forest_new_uniform(cmesh, esquema, nivel_min, 0, comm);

  int nivel_max = nivel_min;
  for (int v : alvo) if (v > nivel_max) nivel_max = v;

  if (nivel_max > nivel_min) {
    Alvos ctx{&alvo};
    t8_forest_t novo;
    t8_forest_init(&novo);
    t8_forest_set_user_data(novo, &ctx);
    t8_forest_set_adapt(novo, f_, adapt_por_arvore, 1 /* recursivo */);
    t8_forest_commit(novo);
    f_ = novo;
  }
}

Malha::~Malha()
{
  if (f_ != nullptr) t8_forest_unref(&f_);
}

const Bloco &Malha::bloco(long b) const
{
  if (b < 0 || b >= num_blocos())
    throw std::out_of_range("Malha::bloco: " + std::to_string(b) + " fora de [0, " +
                            std::to_string(num_blocos() - 1) + "]");
  return blocos_[static_cast<size_t>(b)];
}

long Malha::num_folhas() const
{
  return static_cast<long>(t8_forest_get_local_num_leaf_elements(f_));
}

long Malha::num_folhas_do_bloco(long b) const
{
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha::num_folhas_do_bloco: arvore " + std::to_string(b) +
                            " nao e' local");
  return static_cast<long>(t8_forest_get_tree_num_leaf_elements(f_, static_cast<t8_locidx_t>(b)));
}

Ponto2 Malha::centroide(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha::centroide: arvore " + std::to_string(b) + " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha::centroide: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));
  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(f_, it, e, c);
  return Ponto2{c[0], c[1]};
}

std::vector<Ponto2> Malha::centroides() const
{
  std::vector<Ponto2> v;
  v.reserve(static_cast<size_t>(num_folhas()));
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, ie);
      double c[3] = {0.0, 0.0, 0.0};
      t8_forest_element_centroid(f_, it, e, c);
      v.push_back(Ponto2{c[0], c[1]});
    }
  }
  return v;
}

int Malha::num_vizinhos(long b, long i, int fa) const
{
  if (fa < 0 || fa > 3)
    throw std::out_of_range("Malha::num_vizinhos: face " + std::to_string(fa) +
                            " -- o QUAD tem 4 (0:x-, 1:x+, 2:y-, 3:y+)");
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));

  t8_element_t **vz = nullptr;
  int *dual = nullptr, nviz = 0;
  t8_locidx_t *idx = nullptr;
  t8_eclass_t vec;
  t8_forest_leaf_face_neighbors(f_, it, e, const_cast<const t8_element_t ***>(&vz), fa,
                                &dual, &nviz, &idx, &vec);
  T8_FREE(vz);
  T8_FREE(dual);
  T8_FREE(idx);
  return nviz;
}

}  // namespace t8flow2d
