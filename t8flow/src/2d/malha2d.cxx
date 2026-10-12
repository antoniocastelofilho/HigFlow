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

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
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
    : blocos_(blocos), comm_(comm)
{
  sc_MPI_Comm_rank(comm, &rank_);
  sc_MPI_Comm_size(comm, &nprocs_);

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
  // A FRANJA JA' NASCE AQUI (`do_face_ghost = 1`).  Pedi-la so' quando o
  // paralelo chegar seria deixar o caminho serial e o paralelo diferentes
  // exatamente onde o serial nao discrimina: sem franja, um vizinho de outro
  // processo some, a face passa por CONTORNO DO DOMINIO, e o sistema resolve
  // outro problema sem reclamar.
  f_ = t8_forest_new_uniform(cmesh, esquema, nivel_min, 1, comm);

  int nivel_max = nivel_min;
  for (int v : alvo) if (v > nivel_max) nivel_max = v;

  if (nivel_max > nivel_min) {
    Alvos ctx{&alvo};
    t8_forest_t novo;
    t8_forest_init(&novo);
    t8_forest_set_user_data(novo, &ctx);
    t8_forest_set_adapt(novo, f_, adapt_por_arvore, 1 /* recursivo */);
    t8_forest_set_ghost(novo, 1, T8_GHOST_FACES);
    t8_forest_commit(novo);
    f_ = novo;
  }

  monta_franja();
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

  // O INDICE QUE O T8CODE DEVOLVE E' DA FLORESTA, nao da arvore.  Tratar como
  // indice de arvore entrega um elemento existente, so' que OUTRO -- erro
  // silencioso, e esta' na lista de armadilhas desta frente.
  //
  // E A SOMA `primeiro_id() + indice` ESTA' ERRADA EM PARALELO, que e' o segundo
  // nivel da mesma armadilha.  Ela vale para as folhas LOCAIS, cuja numeracao e'
  // a ordem de percurso; nao vale para as de FRANJA, que o t8code indexa logo
  // depois das locais e cujo id pertence a OUTRO processo.  Somar o deslocamento
  // deste processo a um indice de franja da' o id de uma folha local que existe
  // -- so' que outra.  Em np=1 nao ha' franja e a soma acerta sempre, e e' por
  // isso que o caminho serial nunca acusaria isto.
  //
  // `id_global_do_indice` cobre as duas faixas, e os ids de franja vem da TROCA.
  std::vector<long> ids;
  ids.reserve(static_cast<size_t>(nviz));
  for (int k = 0; k < nviz; k++)
    ids.push_back(id_global_do_indice(static_cast<long>(idx[k])));

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
  // PROPRIEDADE DA MALHA INTEIRA, logo reduz.  Um processo que recebeu apenas a
  // regiao fina reportaria um minimo maior que o verdadeiro, e as guardas que
  // comparam nivel pedido com nivel atual decidiriam com o numero errado.
  //
  // `v < 0` significa particao VAZIA -- possivel em paralelo com poucas folhas --
  // e e' neutralizado para o maximo inteiro antes do minimo, para que a particao
  // vazia nao vire o minimo global.
  if (nprocs_ > 1) {
    int envia = (v < 0) ? std::numeric_limits<int>::max() : v;
    int r = envia;
    sc_MPI_Allreduce(&envia, &r, 1, sc_MPI_INT, sc_MPI_MIN, comm_);
    v = (r == std::numeric_limits<int>::max()) ? -1 : r;
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
  // Mesma razao do minimo: e' da malha inteira.  Aqui `-1` da particao vazia ja'
  // e' neutro para o maximo, e nao precisa de tratamento.
  if (nprocs_ > 1) {
    int r = v;
    sc_MPI_Allreduce(&v, &r, 1, sc_MPI_INT, sc_MPI_MAX, comm_);
    v = r;
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

// ---------------------------------------------------------------------------
// A FRANJA
// ---------------------------------------------------------------------------

void Malha::monta_franja() { monta_franja_camadas(camadas_); }

void Malha::exige_camadas(int k)
{
  if (k < 1)
    throw std::invalid_argument("Malha::exige_camadas: " + std::to_string(k) +
                                " -- a franja tem pelo menos uma camada");
  if (k <= camadas_) return;      // idempotente: pedir menos nao desfaz
  monta_franja_camadas(k);
}

namespace {

//! Indices de FLORESTA dos vizinhos de face -- crus, como o t8code os da'.
//
// Devolve INDICE e nao id de proposito: na primeira fase da montagem da franja o
// id global de um vizinho de franja AINDA NAO EXISTE neste processo -- ele chega
// pela troca.  Converter indice em id antes disso exigiria a aritmetica
// `primeiro_id() + indice`, que e' exatamente o erro que a franja existe para
// impedir: para um indice de franja ela da' o id de uma folha local que existe,
// so' que outra.
std::vector<long> vizinhos_crus(t8_forest_t f, t8_locidx_t it, t8_locidx_t ie, int fa)
{
  const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f, it, ie);
  t8_element_t **vz = nullptr;
  int *dual = nullptr, nviz = 0;
  t8_locidx_t *idx = nullptr;
  t8_eclass_t vec;
  t8_forest_leaf_face_neighbors(f, it, e, const_cast<const t8_element_t ***>(&vz), fa,
                                &dual, &nviz, &idx, &vec);
  std::vector<long> v;
  v.reserve(static_cast<size_t>(nviz));
  for (int k = 0; k < nviz; k++) v.push_back(static_cast<long>(idx[k]));
  T8_FREE(vz);
  T8_FREE(dual);
  T8_FREE(idx);
  return v;
}

}  // namespace

void Malha::monta_franja_camadas(int k)
{
  camadas_ = k;
  saturacoes_ = 0;
  n_locais_ = static_cast<long>(t8_forest_get_local_num_leaf_elements(f_));
  const long nf = static_cast<long>(t8_forest_get_num_ghosts(f_));
  const long p0 = primeiro_id();

  struct Folha {
    double cx = 0.0, cy = 0.0, hx = 0.0, hy = 0.0;
    int nivel = -1;
    std::array<std::vector<long>, 4> viz;

    // O NIVEL DE CADA VIZINHO VIAJA JUNTO COM O ID, e isso economiza uma camada
    // inteira de franja.  Sem ele, descobrir o nivel de um vizinho exigiria que
    // ele proprio estivesse na tabela com geometria -- uma camada mais fundo --,
    // e a medida mostrou que nem quatro camadas resolviam: a carga satura antes,
    // porque em malha adaptada a vizinhanca de raio 3 passa de `(2k+1)^2` folhas.
    //
    // Com o nivel junto, o fluxo precisa de DUAS camadas e os niveis vem de graca.
    std::array<std::vector<int>, 4> viz_nivel;

    // MASCARA DE FACES CONHECIDAS, um bit por face -- e nao um `bool` por folha.
    //
    // A distincao e' obrigatoria e custou um erro de O(1): uma face sem vizinhos
    // na lista pode significar DUAS coisas -- contorno do dominio, ou "o dono
    // desta folha me contou sobre as outras faces e nao sobre esta".  Com um
    // `bool` por folha as duas se confundiam, e a segunda era tratada como a
    // primeira: o esquema aplicava condicao de contorno no MEIO do dominio.
    //
    // Medido: `grau1` em malha adaptada, que o esquema de fluxo resolve
    // EXATAMENTE, deu 3,97e-13 em np=1 e 1,57e-01 em np=2.
    unsigned char faces = 0;

    bool completa() const { return faces == 0x0F; }
  };
  // A tabela e' por ID GLOBAL e nao por indice: acima da primeira camada a mesma
  // folha chega por dois caminhos, e a tabela por id a unifica sem duplicar.
  std::unordered_map<long, Folha> tab;

  // =========================================================================
  // FASE A -- identidade e geometria da camada 1
  // =========================================================================
  //
  // A troca minima: seis campos por folha.  Ela existe separada porque resolve o
  // problema do ovo e da galinha -- os vizinhos de franja de uma folha local sao
  // conhecidos por INDICE, e o id deles so' chega aqui.
  const t8_scheme_c *esq = t8_forest_get_scheme(f_);
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);

  std::vector<long> id_por_indice(static_cast<size_t>(n_locais_ + nf), -1);
  std::vector<std::array<std::vector<long>, 4>> viz_cru(static_cast<size_t>(n_locais_));
  {
    sc_array_t *arr = sc_array_new_count(6 * sizeof(double),
                                         static_cast<size_t>(n_locais_ + nf));
    long kk = 0;
    for (t8_locidx_t it = 0; it < narv; it++) {
      const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
      for (t8_locidx_t ie = 0; ie < ne; ie++) {
        const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, ie);
        double c[3] = {0.0, 0.0, 0.0};
        t8_forest_element_centroid(f_, it, e, c);
        const int nv = esq->element_get_level(t8_forest_get_tree_class(f_, it), e);

        // OS LADOS VEM DO BLOCO E DO NIVEL, como em `lados_da_folha`, e nao do
        // diametro do elemento: a mesma formula nos dois lugares e' o que impede
        // que a folha local e a MESMA folha vista como franja tenham medidas
        // ligeiramente diferentes.  A clausula de conservacao compara pesos de
        // estencil em `0,0e+00`, e uma diferenca de arredondamento ali apareceria
        // como esquema nao conservativo.
        const Bloco &kb = blocos_[static_cast<size_t>(it)];
        const int div = 1 << nv;

        const long id = p0 + kk;
        id_por_indice[static_cast<size_t>(kk)] = id;
        Folha &fo = tab[id];
        fo.cx = c[0];
        fo.cy = c[1];
        fo.hx = (kb.hi[0] - kb.lo[0]) / div;
        fo.hy = (kb.hi[1] - kb.lo[1]) / div;
        fo.nivel = nv;
        for (int fa = 0; fa < 4; fa++)
          viz_cru[static_cast<size_t>(kk)][static_cast<size_t>(fa)] =
              vizinhos_crus(f_, it, ie, fa);

        double *d = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(kk)));
        d[0] = static_cast<double>(id);
        d[1] = fo.cx; d[2] = fo.cy; d[3] = fo.hx; d[4] = fo.hy;
        d[5] = static_cast<double>(nv);
        kk++;
      }
    }

    if (nf > 0) {
      t8_forest_ghost_exchange_data(f_, arr);
      for (long i = n_locais_; i < n_locais_ + nf; i++) {
        const double *d =
            static_cast<const double *>(sc_array_index_int(arr, static_cast<int>(i)));
        const long id = static_cast<long>(d[0]);
        id_por_indice[static_cast<size_t>(i)] = id;
        Folha &fo = tab[id];
        fo.cx = d[1]; fo.cy = d[2]; fo.hx = d[3]; fo.hy = d[4];
        fo.nivel = static_cast<int>(d[5]);
      }
    }
    sc_array_destroy(arr);
  }

  // Agora os vizinhos das LOCAIS podem virar ids: o indice de franja tem id.
  for (long i = 0; i < n_locais_; i++) {
    Folha &fo = tab[id_por_indice[static_cast<size_t>(i)]];
    for (int fa = 0; fa < 4; fa++) {
      for (long x : viz_cru[static_cast<size_t>(i)][static_cast<size_t>(fa)]) {
        const long id = (x >= 0 && x < n_locais_ + nf)
                            ? id_por_indice[static_cast<size_t>(x)] : -1;
        if (id < 0) continue;
        fo.viz[static_cast<size_t>(fa)].push_back(id);
        const auto itv = tab.find(id);
        fo.viz_nivel[static_cast<size_t>(fa)].push_back(
            (itv != tab.end() && itv->second.nivel >= 0) ? itv->second.nivel : -1);
      }
    }
    // A FOLHA LOCAL CONHECE AS QUATRO FACES: o t8code as respondeu todas, e
    // lista vazia aqui significa contorno do dominio, sem ambiguidade.
    fo.faces = 0x0F;
  }

  // =========================================================================
  // FASE B -- as camadas acima da primeira
  // =========================================================================
  //
  // Cada rodada publica, por folha LOCAL, a vizinhanca que este processo ja'
  // conhece em torno dela -- com os vizinhos de cada folha publicada --, e
  // recebe o mesmo dos donos das folhas de franja.  Isso estende o alcance em um
  // anel por rodada: a rodada 1 entrega os vizinhos da camada 1 (logo a camada 2
  // existe como ponto), a rodada 2 entrega os vizinhos da camada 2, e assim por
  // diante.
  //
  // O CUSTO CRESCE COM O QUADRADO DE K, e e' por isso que K e' escolha explicita
  // de quem monta o esquema e nao um padrao generoso.
  if (nf > 0 && k > 1) {
    const int MAXV = (2 * k + 1) * (2 * k + 1);
    // 7 campos proprios mais 8 vizinhos de DOIS campos cada (codigo e nivel).
    const int CAMPOS = 7 + 8 * 2;
    const size_t passo = static_cast<size_t>(MAXV) * static_cast<size_t>(CAMPOS);

    for (int rodada = 1; rodada < k; rodada++) {
      sc_array_t *arr = sc_array_new_count(passo * sizeof(double),
                                           static_cast<size_t>(n_locais_ + nf));

      for (long i = 0; i < n_locais_ + nf; i++) {
        double *b = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(i)));
        for (size_t z = 0; z < passo; z++) b[z] = -1.0;
      }

      for (long i = 0; i < n_locais_; i++) {
        double *b = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(i)));

        // o conjunto a publicar: a propria e tudo a distancia <= `rodada`
        const long eu = id_por_indice[static_cast<size_t>(i)];
        std::vector<long> fila{eu};
        std::set<long> posto{eu};
        for (int r = 0; r < rodada; r++) {
          std::vector<long> prox;
          for (long id : fila) {
            const auto it = tab.find(id);
            if (it == tab.end() || !it->second.completa()) continue;
            for (int fa = 0; fa < 4; fa++)
              for (long v : it->second.viz[static_cast<size_t>(fa)])
                if (posto.insert(v).second) prox.push_back(v);
          }
          fila = prox;
        }

        int slot = 0;
        for (long id : posto) {
          // A CARGA SATUROU.  Em malha adaptada a vizinhanca de raio r passa de
          // `(2k+1)^2` folhas, porque uma celula grossa tem muitos vizinhos
          // finos.  Cortar aqui e' inevitavel -- o canal e' de tamanho fixo --,
          // mas cortar EM SILENCIO nao e': foi isso que fez quatro camadas ainda
          // deixarem niveis ausentes, e levou a' busca pela camada errada.
          if (slot >= MAXV) { saturacoes_++; break; }
          const auto it = tab.find(id);
          if (it == tab.end() || it->second.nivel < 0) continue;
          const Folha &fo = it->second;
          double *c = b + static_cast<size_t>(slot) * CAMPOS;
          c[0] = static_cast<double>(id);
          c[1] = fo.cx; c[2] = fo.cy; c[3] = fo.hx; c[4] = fo.hy;
          c[5] = static_cast<double>(fo.nivel);
          // A MASCARA VAI JUNTO, e e' o que permite a quem recebe distinguir
          // "face de contorno do dominio" de "face sobre a qual nao me contaram".
          c[6] = static_cast<double>(fo.faces);
          int z = 7;
          for (int fa = 0; fa < 4; fa++)
            for (size_t q = 0; q < fo.viz[static_cast<size_t>(fa)].size(); q++) {
              if (z + 1 >= CAMPOS) break;
              // A FACE vai codificada no id, para que quem recebe saiba por qual
              // face o vizinho entra -- sem isso a malha remota chegaria como uma
              // lista de ids sem topologia.  E o NIVEL vai no campo seguinte.
              c[z++] = static_cast<double>(fo.viz[static_cast<size_t>(fa)][q] * 4 + fa);
              c[z++] = static_cast<double>(fo.viz_nivel[static_cast<size_t>(fa)][q]);
            }
          slot++;
        }
      }

      t8_forest_ghost_exchange_data(f_, arr);

      for (long i = n_locais_; i < n_locais_ + nf; i++) {
        const double *b =
            static_cast<const double *>(sc_array_index_int(arr, static_cast<int>(i)));
        for (int slot = 0; slot < MAXV; slot++) {
          const double *c = b + static_cast<size_t>(slot) * CAMPOS;
          if (c[0] < 0.0) continue;
          const long id = static_cast<long>(c[0]);
          Folha &fo = tab[id];
          // A FOLHA LOCAL E' A VERSAO AUTORITATIVA DELA MESMA: o que chega pela
          // franja nao sobrescreve o que este processo calculou.  No melhor caso
          // seria trocar dados iguais por dados iguais; no pior, esconderia
          // divergencia entre os dois caminhos.
          if (fo.nivel < 0) {
            fo.cx = c[1]; fo.cy = c[2]; fo.hx = c[3]; fo.hy = c[4];
            fo.nivel = static_cast<int>(c[5]);
          }
          if (!fo.completa()) {
            std::array<std::vector<long>, 4> nv;
            std::array<std::vector<int>, 4> nn;
            for (int z = 7; z + 1 < CAMPOS; z += 2) {
              if (c[z] < 0.0) continue;
              const long cod = static_cast<long>(c[z]);
              const size_t fa = static_cast<size_t>(cod % 4);
              nv[fa].push_back(cod / 4);
              nn[fa].push_back(static_cast<int>(c[z + 1]));
            }
            // A MASCARA RECEBIDA MANDA, e nao a presenca de vizinhos na lista:
            // uma face de contorno do dominio chega com a mascara marcada e a
            // lista vazia, que e' exatamente o que se quer poder afirmar.
            const unsigned char mk = static_cast<unsigned char>(c[6]);
            if (mk != 0) { fo.viz = nv; fo.viz_nivel = nn; fo.faces = mk; }
          }
        }
      }

      sc_array_destroy(arr);
    }
  }

  // =========================================================================
  // ACHATA a tabela em arrays indexados por indice de floresta
  // =========================================================================
  //
  // A ORDEM AQUI NAO E' LIVRE, e supor que era custou um erro de O(1).
  //
  // O t8code numera as folhas que `leaf_face_neighbors` devolve assim, e a
  // documentacao dele e' explicita: `0 .. num_local-1` para as locais e
  // `num_local .. num_local+num_ghosts-1` para as da CAMADA 1.  Logo as duas
  // primeiras faixas tem de sair exatamente nessa ordem.
  //
  // A primeira versao ordenou TODA a franja por id, porque id e' determinista e
  // parecia a escolha limpa.  O efeito: o indice que o t8code devolvia para o
  // terceiro fantasma passava a designar outra folha, e `ids_vizinhos` entregava
  // um id EXISTENTE -- so' que de outra celula.  Medido em `grau1` adaptado, que
  // o esquema resolve exatamente: 3,97e-13 em np=1 e 1,70e-01 em np=2.
  //
  // So' as camadas 2 e acima podem ser ordenadas por id: elas nao tem indice do
  // t8code, e a ordenacao por id as torna reprodutiveis entre execucoes.
  id_g_.clear(); cen_.clear(); lad_.clear(); niv_.clear(); viz_.clear();
  viz_niv_.clear(); expansivel_.clear(); idx_por_id_.clear();

  auto empurra = [&](long id, const Folha &fo) {
    idx_por_id_[id] = static_cast<long>(id_g_.size());
    id_g_.push_back(id);
    cen_.push_back(fo.cx); cen_.push_back(fo.cy);
    lad_.push_back(fo.hx); lad_.push_back(fo.hy);
    niv_.push_back(fo.nivel);
    viz_.push_back(fo.viz);
    viz_niv_.push_back(fo.viz_nivel);
    expansivel_.push_back(fo.completa() ? 1 : 0);
  };

  for (long i = 0; i < n_locais_ + nf; i++) {
    const long id = id_por_indice[static_cast<size_t>(i)];
    if (id < 0) continue;
    empurra(id, tab[id]);
  }
  std::vector<long> resto;
  for (const auto &kv : tab)
    if (kv.second.nivel >= 0 && idx_por_id_.find(kv.first) == idx_por_id_.end())
      resto.push_back(kv.first);
  std::sort(resto.begin(), resto.end());
  for (long id : resto) empurra(id, tab[id]);

  n_fantasmas_ = static_cast<long>(id_g_.size()) - n_locais_;
}

std::vector<long> Malha::vizinhos_do_indice(long idx, int fa) const
{
  if (idx < 0 || idx >= static_cast<long>(id_g_.size()))
    throw std::out_of_range("Malha::vizinhos_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  if (fa < 0 || fa > 3)
    throw std::out_of_range("Malha::vizinhos_do_indice: face " + std::to_string(fa) +
                            " -- o QUAD tem 4 (0:x-, 1:x+, 2:y-, 3:y+)");
  return viz_[static_cast<size_t>(idx)][static_cast<size_t>(fa)];
}

std::vector<int> Malha::niveis_vizinhos_do_indice(long idx, int fa) const
{
  if (idx < 0 || idx >= static_cast<long>(viz_niv_.size()))
    throw std::out_of_range("Malha::niveis_vizinhos_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  if (fa < 0 || fa > 3)
    throw std::out_of_range("Malha::niveis_vizinhos_do_indice: face " +
                            std::to_string(fa) + " -- o QUAD tem 4");
  return viz_niv_[static_cast<size_t>(idx)][static_cast<size_t>(fa)];
}

bool Malha::expansivel_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(expansivel_.size()))
    throw std::out_of_range("Malha::expansivel_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  return expansivel_[static_cast<size_t>(idx)] != 0;
}

long Malha::id_global_do_indice(long idx) const
{
  if (idx < 0 || idx >= n_locais_ + n_fantasmas_)
    throw std::out_of_range("Malha::id_global_do_indice: " + std::to_string(idx) +
                            " fora de [0, " + std::to_string(n_locais_ + n_fantasmas_) +
                            ") -- " + std::to_string(n_locais_) + " locais mais " +
                            std::to_string(n_fantasmas_) + " de franja");
  return id_g_[static_cast<size_t>(idx)];
}

Ponto2 Malha::centroide_do_indice(long idx) const
{
  if (idx < 0 || idx >= n_locais_ + n_fantasmas_)
    throw std::out_of_range("Malha::centroide_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas locais e de franja");
  return Ponto2{cen_[static_cast<size_t>(2 * idx)], cen_[static_cast<size_t>(2 * idx + 1)]};
}

Ponto2 Malha::lados_do_indice(long idx) const
{
  if (idx < 0 || idx >= n_locais_ + n_fantasmas_)
    throw std::out_of_range("Malha::lados_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas locais e de franja");
  return Ponto2{lad_[static_cast<size_t>(2 * idx)], lad_[static_cast<size_t>(2 * idx + 1)]};
}

int Malha::nivel_do_indice(long idx) const
{
  if (idx < 0 || idx >= n_locais_ + n_fantasmas_)
    throw std::out_of_range("Malha::nivel_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas locais e de franja");
  return niv_[static_cast<size_t>(idx)];
}

long Malha::indice_do_id(long id) const
{
  const auto it = idx_por_id_.find(id);
  return (it == idx_por_id_.end()) ? -1 : it->second;
}

void Malha::reparticiona()
{
  // PARTICAO EM COMMIT PROPRIO, e a franja junto com ela.
  //
  // O comentario do construtor registra a armadilha: `set_adapt`, `set_ghost` e
  // `set_partition` no MESMO commit derruba com np>1 e funciona em np=1.  Aqui
  // ha' so' particao e franja, que o t8code combina, e nenhum adapt.
  t8_forest_t novo;
  t8_forest_init(&novo);
  t8_forest_set_partition(novo, f_, 0);
  t8_forest_set_ghost(novo, 1, T8_GHOST_FACES);
  t8_forest_commit(novo);
  f_ = novo;
  monta_franja();
}

long Malha::id_do_ponto(const Ponto2 &p) const
{
  // PERCORRE AS LOCAIS E AS DE FRANJA, pela mesma regra, SEM preferir as locais.
  //
  // Preferir a local seria o reflexo natural e estragaria a conservacao: o
  // ponto de quadratura sobre uma faceta que e' fronteira de PARTICAO seria
  // atribuido a uma celula diferente em cada um dos dois processos, os suportes
  // divergiriam, e os dois fluxos pela mesma faceta deixariam de se cancelar.
  //
  // A regra SEMIABERTA [lo, hi) e' global e nao depende de quem pergunta, e e'
  // isso que faz os dois processos chegarem a' MESMA celula.  Com ela, o
  // criterio aqui e' puramente geometrico, e a unica diferenca entre serial e
  // paralelo e' a faixa de folhas que este processo conhece.
  const long n = n_locais_ + n_fantasmas_;
  for (long x = 0; x < n; x++) {
    const Ponto2 c = centroide_do_indice(x);
    const Ponto2 h = lados_do_indice(x);

    // O BLOCO da folha e' identificado pela geometria, porque a folha de franja
    // nao traz a arvore do dono.  Com blocos disjuntos o centro decide sozinho.
    const Bloco *kb = nullptr;
    for (const Bloco &k : blocos_) {
      bool no_bloco = true;
      for (int d = 0; d < 2; d++)
        if (c[d] < k.lo[d] - EPS || c[d] > k.hi[d] + EPS) no_bloco = false;
      if (no_bloco) { kb = &k; break; }
    }
    if (kb == nullptr) continue;

    bool dentro = true;
    for (int d = 0; d < 2; d++) {
      const double lo = c[d] - h[d] / 2.0;
      const double hi = c[d] + h[d] / 2.0;
      // SEMIABERTO [lo, hi): ponto sobre a face comum cai na celula de
      // coordenada maior.  No limite superior do DOMINIO o fechado volta, senao
      // a borda de cima nao pertenceria a celula nenhuma.
      const bool ultima = (hi > kb->hi[d] - EPS);
      const bool ok = (p[d] >= lo - EPS) && (ultima ? (p[d] <= hi + EPS) : (p[d] < hi - EPS));
      if (!ok) { dentro = false; break; }
    }
    if (dentro) return id_global_do_indice(x);
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
  // A FRANJA ACOMPANHA O BALANCEAMENTO.  O t8code combina `set_balance` com
  // `set_ghost` no mesmo commit; o que ele nao aceita e' `set_adapt` junto --
  // e por isso o refino ficou no commit anterior.
  t8_forest_set_ghost(equilibrada, 1, T8_GHOST_FACES);
  t8_forest_commit(equilibrada);

  // Sem `unref`: o commit CONSOME a floresta de origem.  Liberar de novo deixa a
  // malha vazia sem erro nenhum -- ja' aconteceu aqui.
  f_ = equilibrada;
  monta_franja();
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
  // CONTAGEM DA MALHA INTEIRA.  Os testes a usam para afirmar que a malha de
  // fato TEM transicao 2:1 -- e um processo cuja particao caiu toda dentro da
  // regiao uniforme contaria zero, fazendo o caso passar por nao medir nada.
  if (nprocs_ > 1) {
    long r = n;
    sc_MPI_Allreduce(&n, &r, 1, sc_MPI_LONG, sc_MPI_SUM, comm_);
    n = r;
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
  // A FRANJA ACOMPANHA O BALANCEAMENTO.  O t8code combina `set_balance` com
  // `set_ghost` no mesmo commit; o que ele nao aceita e' `set_adapt` junto --
  // e por isso o refino ficou no commit anterior.
  t8_forest_set_ghost(equilibrada, 1, T8_GHOST_FACES);
  t8_forest_commit(equilibrada);

  // NAO HA' `unref` AQUI, e isto custou uma malha vazia: o `t8_forest_commit`
  // CONSOME a floresta de origem passada em `set_adapt`/`set_balance` -- ele faz
  // o unref dela.  Liberar de novo nao da' erro nenhum: a floresta fica invalida,
  // as contagens vao a ZERO, e a malha parece simplesmente... vazia.  O caso
  // `o_refino_conserva_a_area_de_cada_bloco` pegou (area 0 contra 1), e nenhum
  // "funcionou" teria pegado.
  f_ = equilibrada;
  monta_franja();
}

}  // namespace t8flow2d
