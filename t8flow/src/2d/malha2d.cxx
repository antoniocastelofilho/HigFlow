// Implementacao da malha multiblocos 2D.  O PORQUE do desenho -- multiarvore,
// juncoes, limites -- esta' no `malha2d.hxx`.  Aqui ficam as decisoes de
// implementacao, e cada uma tem uma razao que o cabecalho nao cobre:
//
//   1. `EPS` e as comparacoes.  `igual()` e' RELATIVA (escala com os valores) e
//      `sobrepoe()` exige sobreposicao em MEDIDA -- toque de borda nao conta.
//      Sem essa distincao, dois blocos encostados contariam como sobrepostos e a
//      malha seria recusada;
//   2. o nivel alvo por arvore vai pelo `user_data` da floresta, em `struct
//      Alvos`, porque o callback de adapt do t8code nao aceita contexto proprio.
//      Por PONTEIRO para vetor e nao por estatico de arquivo: estatico impediria
//      duas malhas vivas ao mesmo tempo;
//   3. o criterio do refino local julga o CENTRO da celula, e a consequencia
//      disso esta' comentada em `adapt_caixa`;
//   4. os commits do t8code sao SEPARADOS (refino num, equilibrio em outro), e o
//      `t8_forest_commit` CONSOME a floresta de origem -- as duas coisas estao
//      comentadas onde acontecem, porque errar qualquer uma falha em silencio.
#include "malha2d.hxx"

#include <cmath>
#include <cstring>
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

//! Contexto do refino ALEATORIO, visto pelo callback.
struct Sorteio {
  unsigned long semente;
  double fracao;
  int nivel_max;
};

//! Hash de 64 bits dos bits de um double, no estilo splitmix.
//
// Escrito aqui e nao tirado de `<random>` de proposito: o que se quer e' uma
// funcao DETERMINISTICA da posicao, igual em qualquer compilador e em qualquer
// ordem de chamada.  Um `std::mt19937` daria uma SEQUENCIA, e sequencia depende da
// ordem em que e' consumida -- que e' justamente o que nao pode acontecer aqui.
unsigned long mistura(unsigned long x)
{
  x += 0x9E3779B97F4A7C15UL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBUL;
  return x ^ (x >> 31);
}

//! Sorteio em [0,1) determinado pela POSICAO da celula e pela semente.
double sorteio_do_ponto(const double c[3], unsigned long semente)
{
  unsigned long h = mistura(semente);
  for (int d = 0; d < 2; d++) {
    unsigned long bits = 0;
    std::memcpy(&bits, &c[d], sizeof bits);
    h = mistura(h ^ bits);
  }
  // 53 bits para o mantissa do double: o resto nao acrescenta informacao.
  return static_cast<double>(h >> 11) / 9007199254740992.0;
}

int adapt_sorteio(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                  const t8_eclass_t tree_class, t8_locidx_t /*lelement_id*/,
                  const t8_scheme_c *scheme, const int /*is_family*/,
                  const int /*num_elements*/, t8_element_t *elements[])
{
  const Sorteio *so = static_cast<const Sorteio *>(t8_forest_get_user_data(forest));
  if (so == nullptr) return 0;
  if (scheme->element_get_level(tree_class, elements[0]) >= so->nivel_max) return 0;

  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(forest_from, which_tree, elements[0], c);
  return (sorteio_do_ponto(c, so->semente) < so->fracao) ? 1 : 0;
}

//! Contexto do refino local, visto pelo callback.
struct Caixa {
  double lo[2];
  double hi[2];
  int nivel_alvo;
};

int adapt_caixa(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                const t8_eclass_t tree_class, t8_locidx_t /*lelement_id*/,
                const t8_scheme_c *scheme, const int /*is_family*/,
                const int /*num_elements*/, t8_element_t *elements[])
{
  const Caixa *cx = static_cast<const Caixa *>(t8_forest_get_user_data(forest));
  if (cx == nullptr) return 0;

  if (scheme->element_get_level(tree_class, elements[0]) >= cx->nivel_alvo) return 0;

  // O CRITERIO E' O CENTRO DA CELULA, nao a interseccao com a caixa, e a
  // diferenca importa: por interseccao, uma celula grossa que toca a caixa de
  // raspao e' refinada inteira, e o refino vaza para fora da regiao pedida a cada
  // nivel.  Pelo centro, a regiao refinada converge para a caixa.
  //
  // CONSEQUENCIA A NAO ESQUECER, e ela derrubou a primeira versao do teste: o
  // criterio julga a celula ANTES de dividir, e os filhos HERDAM a decisao do
  // pai.  Logo a regiao efetivamente refinada e' a caixa DILATADA por ate' uma
  // celula do nivel final -- um filho de uma celula cujo centro estava na caixa
  // pode ter o seu proprio centro fora dela.  Quem precisa de refino contido na
  // caixa tem de pedir a caixa menor, ou julgar a celula pelos cantos.
  double c[3];
  t8_forest_element_centroid(forest_from, which_tree, elements[0], c);
  for (int d = 0; d < 2; d++)
    if (c[d] < cx->lo[d] || c[d] > cx->hi[d]) return 0;
  return 1;
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

int Malha::nivel(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha::nivel: arvore " + std::to_string(b) + " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha::nivel: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));
  const t8_scheme_c *esquema = t8_forest_get_scheme(f_);
  return esquema->element_get_level(t8_forest_get_tree_class(f_, it), e);
}

double Malha::area(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha::area: arvore " + std::to_string(b) + " nao e' local");
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));
  return t8_forest_element_volume(f_, it, e);
}

std::vector<int> Malha::niveis_vizinhos(long b, long i, int fa) const
{
  if (fa < 0 || fa > 3)
    throw std::out_of_range("Malha::niveis_vizinhos: face " + std::to_string(fa) +
                            " -- o QUAD tem 4 (0:x-, 1:x+, 2:y-, 3:y+)");
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));

  t8_element_t **vz = nullptr;
  int *dual = nullptr, nviz = 0;
  t8_locidx_t *idx = nullptr;
  t8_eclass_t vec;
  t8_forest_leaf_face_neighbors(f_, it, e, const_cast<const t8_element_t ***>(&vz), fa,
                                &dual, &nviz, &idx, &vec);

  // O NIVEL SAI DO ELEMENTO DEVOLVIDO, nao do indice.  O indice que esta funcao
  // devolve e' da FLORESTA local, nao da arvore -- trata-lo como indice de arvore
  // entrega um elemento existente, so' que OUTRO, e o erro e' silencioso.
  const t8_scheme_c *esquema = t8_forest_get_scheme(f_);
  std::vector<int> niveis;
  niveis.reserve(static_cast<size_t>(nviz));
  for (int k = 0; k < nviz; k++)
    niveis.push_back(esquema->element_get_level(vec, vz[k]));

  T8_FREE(vz);
  T8_FREE(dual);
  T8_FREE(idx);
  return niveis;
}

std::vector<long> Malha::ids_vizinhos(long b, long i, int fa) const
{
  if (fa < 0 || fa > 3)
    throw std::out_of_range("Malha::ids_vizinhos: face " + std::to_string(fa) +
                            " -- o QUAD tem 4 (0:x-, 1:x+, 2:y-, 3:y+)");
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));

  t8_element_t **vz = nullptr;
  int *dual = nullptr, nviz = 0;
  t8_locidx_t *idx = nullptr;
  t8_eclass_t vec;
  t8_forest_leaf_face_neighbors(f_, it, e, const_cast<const t8_element_t ***>(&vz), fa,
                               &dual, &nviz, &idx, &vec);

  // O INDICE QUE O T8CODE DEVOLVE E' DA FLORESTA LOCAL, nao da arvore.  Tratar
  // como indice de arvore entrega um elemento existente, so' que OUTRO -- erro
  // silencioso, e esta' na lista de armadilhas desta frente.  O id global e'
  // `primeiro_id() + indice_de_floresta`, direto, porque a numeracao local E' a
  // ordem de percurso da floresta.
  const long p0 = primeiro_id();
  std::vector<long> ids;
  ids.reserve(static_cast<size_t>(nviz));
  for (int k = 0; k < nviz; k++) ids.push_back(p0 + static_cast<long>(idx[k]));

  T8_FREE(vz);
  T8_FREE(dual);
  T8_FREE(idx);
  return ids;
}

int Malha::nivel_minimo() const
{
  int v = -1;
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const int n = nivel(it, ie);
      if (v < 0 || n < v) v = n;
    }
  }
  return v;
}

int Malha::nivel_maximo() const
{
  int v = -1;
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const int n = nivel(it, ie);
      if (n > v) v = n;
    }
  }
  return v;
}

long Malha::primeiro_id() const
{
  return static_cast<long>(t8_forest_get_first_local_leaf_element_id(f_));
}

long Malha::num_folhas_globais() const
{
  return static_cast<long>(t8_forest_get_global_num_leaf_elements(f_));
}

long Malha::id_global(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha::id_global: arvore " + std::to_string(b) + " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha::id_global: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));

  // O deslocamento da arvore dentro da floresta local vem do t8code, nao de uma
  // soma nossa: somar as contagens das arvores anteriores daria o mesmo numero
  // hoje e passaria a divergir no dia em que a floresta tiver arvore vazia.
  const long offset = static_cast<long>(t8_forest_get_tree_element_offset(f_, it));
  return primeiro_id() + offset + i;
}

Ponto2 Malha::lados_da_folha(long b, long i) const
{
  const Bloco &k = bloco(b);
  const int div = 1 << nivel(b, i);
  return Ponto2{(k.hi[0] - k.lo[0]) / div, (k.hi[1] - k.lo[1]) / div};
}

long Malha::id_do_ponto(const Ponto2 &p) const
{
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const Bloco &k = blocos_[static_cast<size_t>(it)];

    // O bloco e' semiaberto como as celulas, com a mesma excecao no limite
    // superior do dominio -- senao um ponto na borda de cima nao pertenceria a
    // bloco nenhum.
    bool no_bloco = true;
    for (int d = 0; d < 2; d++)
      if (p[d] < k.lo[d] - EPS || p[d] > k.hi[d] + EPS) no_bloco = false;
    if (!no_bloco) continue;

    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const Ponto2 c = centroide(it, ie);
      const Ponto2 h = lados_da_folha(it, ie);

      bool dentro = true;
      for (int d = 0; d < 2; d++) {
        const double lo = c[d] - h[d] / 2.0;
        const double hi = c[d] + h[d] / 2.0;
        // SEMIABERTO [lo, hi): ponto sobre a face comum cai na celula de
        // coordenada maior.  No limite superior do DOMINIO o fechado volta, senao
        // a borda de cima nao pertenceria a celula nenhuma.
        const bool ultima = (hi > k.hi[d] - EPS);
        const bool ok = (p[d] >= lo - EPS) && (ultima ? (p[d] <= hi + EPS) : (p[d] < hi - EPS));
        if (!ok) { dentro = false; break; }
      }
      if (dentro) return id_global(it, ie);
    }
  }
  return -1;
}

void Malha::refina_aleatorio(double fracao, unsigned long semente, int nivel_max)
{
  if (!(fracao > 0.0) || fracao > 1.0)
    throw std::invalid_argument("Malha::refina_aleatorio: fracao = " +
                                std::to_string(fracao) + " fora de (0, 1]");
  if (nivel_max <= nivel_maximo())
    throw std::invalid_argument(
        "Malha::refina_aleatorio: nivel_max " + std::to_string(nivel_max) +
        " nao e' maior que o nivel atual " + std::to_string(nivel_maximo()) +
        " -- nada seria refinado, e aceitar isso deixaria a malha PARECER adaptada");

  Sorteio so{semente, fracao, nivel_max};

  // COMMIT 1: o refino.  NAO recursivo, de proposito: recursivo faria cada celula
  // sorteada descer ate' `nivel_max` de uma vez, e o resultado seria um refino
  // profundo em poucos lugares em vez de espalhado.  Uma passada sobe UM nivel; o
  // chamador repete se quiser mais.
  t8_forest_t refinada;
  t8_forest_init(&refinada);
  t8_forest_set_user_data(refinada, &so);
  t8_forest_set_adapt(refinada, f_, adapt_sorteio, 0);
  t8_forest_commit(refinada);

  // COMMIT 2: o equilibrio 2:1, separado -- ver o `refina_caixa`.  Sem ele, duas
  // celulas sorteadas lado a lado com uma nao-sorteada no meio dariam salto de
  // dois niveis.
  t8_forest_t equilibrada;
  t8_forest_init(&equilibrada);
  t8_forest_set_balance(equilibrada, refinada, 0);
  t8_forest_commit(equilibrada);

  // Sem `unref`: o commit CONSOME a floresta de origem.  Liberar de novo deixa a
  // malha vazia sem erro nenhum -- ja' aconteceu aqui.
  f_ = equilibrada;
}

long Malha::num_folhas_em_transicao() const
{
  long n = 0;
  for (long b = 0; b < num_blocos(); b++)
    for (long i = 0; i < num_folhas_do_bloco(b); i++) {
      const int nv = nivel(b, i);
      bool transicao = false;
      for (int fa = 0; fa < 4 && !transicao; fa++)
        for (int viz : niveis_vizinhos(b, i, fa))
          if (viz != nv) { transicao = true; break; }
      if (transicao) n++;
    }
  return n;
}

void Malha::refina_caixa(const Ponto2 &lo, const Ponto2 &hi, int nivel_alvo)
{
  for (int d = 0; d < 2; d++)
    if (!(hi[d] > lo[d] + EPS))
      throw std::invalid_argument(
          "Malha::refina_caixa: caixa degenerada na direcao " + std::to_string(d));

  // A caixa tem de alcancar algum bloco, E tem de pedir mais nivel do que esse
  // bloco ja' tem.  Refino que nao refina nada e' quase sempre engano de quem
  // chamou -- e aceito calado deixaria a malha PARECER refinada.
  bool alcanca = false, aprofunda = false;
  for (const Bloco &k : blocos_) {
    if (sobrepoe(k.lo[0], k.hi[0], lo[0], hi[0]) && sobrepoe(k.lo[1], k.hi[1], lo[1], hi[1])) {
      alcanca = true;
      if (nivel_alvo > k.nivel) aprofunda = true;
    }
  }
  if (!alcanca)
    throw std::invalid_argument("Malha::refina_caixa: a caixa nao intersecta bloco nenhum");
  if (!aprofunda)
    throw std::invalid_argument(
        "Malha::refina_caixa: nivel_alvo " + std::to_string(nivel_alvo) +
        " nao e' maior que o nivel de nenhum bloco que a caixa alcanca -- nada seria refinado");

  Caixa cx{{lo[0], lo[1]}, {hi[0], hi[1]}, nivel_alvo};

  // COMMIT 1: o refino.  Recursivo, para subir mais de um nivel numa passada.
  t8_forest_t refinada;
  t8_forest_init(&refinada);
  t8_forest_set_user_data(refinada, &cx);
  t8_forest_set_adapt(refinada, f_, adapt_caixa, 1 /* recursivo */);
  t8_forest_commit(refinada);

  // COMMIT 2: o equilibrio 2:1.  Em commit SEPARADO -- ver o cabecalho.  Sem
  // ele, um refino de dois niveis numa caixa deixa 4:1 na borda dela, a malha
  // commita, as contagens fecham, e o salto atravessa o estencil do
  // discretizador sem aviso.
  t8_forest_t equilibrada;
  t8_forest_init(&equilibrada);
  t8_forest_set_balance(equilibrada, refinada, 0);
  t8_forest_commit(equilibrada);

  // NAO HA' `unref` AQUI, e isto custou uma malha vazia: o `t8_forest_commit`
  // CONSOME a floresta de origem passada em `set_adapt`/`set_balance` -- ele faz
  // o unref dela.  Liberar de novo nao da' erro nenhum: a floresta fica invalida,
  // as contagens vao a ZERO, e a malha parece simplesmente... vazia.  O caso
  // `o_refino_conserva_a_area_de_cada_bloco` pegou (area 0 contra 1), e nenhum
  // "funcionou" teria pegado.
  f_ = equilibrada;
}

}  // namespace t8flow2d
