#include "malha3d.hxx"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>

#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_schemes/t8_default/t8_default.hxx>

namespace t8flow3d {

namespace {

//! Tolerancia geometrica.  A mesma de duas dimensoes, e pela mesma razao: as
//! coordenadas vem de divisoes por potencias de dois, logo os limites coincidem
//! em bits quando deveriam coincidir, e o `EPS` cobre apenas o acumulo.
constexpr double EPS = 1e-12;

bool coincide(double a, double b) { return std::fabs(a - b) < EPS; }

bool sobrepoe(double a0, double a1, double b0, double b1)
{
  return a0 < b1 - EPS && b0 < a1 - EPS;
}

//! A chave do indice espacial: `(bloco, nivel, i, j, k)` numa palavra.
//
// A REPARTICAO DOS BITS E' MAIS APERTADA QUE EM DUAS DIMENSOES, porque ha' uma
// coordenada a mais: 8 para o bloco (256 blocos), 5 para o nivel (32 niveis, e o
// t8code nao passa de 30 em hexaedro), 17 para cada coordenada (131 mil celulas
// por direcao, que em tres dimensoes ja' sao 2e15 celulas -- muito acima do que
// cabe em memoria).  Somam 59 dos 64.
//
// O BLOCO ENTRA NA CHAVE, e esquece-lo foi um defeito medido em duas dimensoes:
// as coordenadas `(i, j, k)` sao RELATIVAS ao bloco, logo a celula `(0,0,0)` do
// bloco 1 e a do bloco 0 colidem sem ele.  O arreio de numeracao o pegou na
// primeira execucao.
unsigned long long chave_celula(long bloco, int nivel, long i, long j, long k)
{
  return (static_cast<unsigned long long>(bloco & 0xFF) << 56) |
         (static_cast<unsigned long long>(nivel & 0x1F) << 51) |
         (static_cast<unsigned long long>(i & 0x1FFFF) << 34) |
         (static_cast<unsigned long long>(j & 0x1FFFF) << 17) |
         static_cast<unsigned long long>(k & 0x1FFFF);
}

//! Indices de FLORESTA dos vizinhos de face -- crus, como o t8code os da'.
//
// Devolve INDICE e nao id de proposito: na primeira fase da montagem da franja o
// id global de um vizinho de franja AINDA NAO EXISTE neste processo.  Converter
// antes exigiria a aritmetica `primeiro_id() + indice`, que para um indice de
// franja da' o id de uma folha local que existe -- so' que outra.
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

//! Contexto do refino por nivel-alvo por arvore.
struct Alvos {
  const std::vector<int> *alvo;
};

int adapt_por_arvore(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree,
                     const t8_eclass_t tclass, t8_locidx_t lelement,
                     const t8_scheme_c *scheme, const int is_family,
                     const int num_elements, t8_element_t *elements[])
{
  const Alvos *ctx = static_cast<const Alvos *>(t8_forest_get_user_data(forest));
  const t8_gloidx_t gt = t8_forest_global_tree_id(forest_from, ltree);
  const int alvo = (*ctx->alvo)[static_cast<size_t>(gt)];
  return (scheme->element_get_level(tclass, elements[0]) < alvo) ? 1 : 0;
}

//! Contexto do refino em caixa.
struct Caixa {
  Ponto3 lo, hi;
  int nivel_alvo;
};

int adapt_caixa(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree,
                const t8_eclass_t tclass, t8_locidx_t lelement,
                const t8_scheme_c *scheme, const int is_family,
                const int num_elements, t8_element_t *elements[])
{
  const Caixa *cx = static_cast<const Caixa *>(t8_forest_get_user_data(forest));
  if (scheme->element_get_level(tclass, elements[0]) >= cx->nivel_alvo) return 0;
  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(forest_from, ltree, elements[0], c);
  for (int d = 0; d < 3; d++)
    if (c[d] < cx->lo[d] - EPS || c[d] > cx->hi[d] + EPS) return 0;
  return 1;
}

//! Contexto do refino aleatorio.
struct Sorteio {
  double fracao;
  unsigned long semente;
  int nivel_max;
};

//! Decide pelo HASH do centroide, e nao por um gerador sequencial.
//
// A RAZAO E' O PARALELO: um gerador sequencial da' resultados que dependem da
// ORDEM DE PERCURSO, logo da particao -- a mesma malha refinaria diferente em
// np=2 e np=4, e nenhuma comparacao entre eles valeria.  O hash do centroide
// depende so' da geometria.
int adapt_aleatorio(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree,
                    const t8_eclass_t tclass, t8_locidx_t lelement,
                    const t8_scheme_c *scheme, const int is_family,
                    const int num_elements, t8_element_t *elements[])
{
  const Sorteio *so = static_cast<const Sorteio *>(t8_forest_get_user_data(forest));
  if (scheme->element_get_level(tclass, elements[0]) >= so->nivel_max) return 0;

  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(forest_from, ltree, elements[0], c);

  unsigned long long h = so->semente;
  for (int d = 0; d < 3; d++) {
    unsigned long long bits = 0;
    static_assert(sizeof(bits) == sizeof(c[d]), "double de 64 bits");
    std::memcpy(&bits, &c[d], sizeof bits);
    h ^= bits + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
  }
  // splitmix, para espalhar os bits do hash
  h += 0x9E3779B97F4A7C15ULL;
  h = (h ^ (h >> 30)) * 0xBF58476D1CE4E5B9ULL;
  h = (h ^ (h >> 27)) * 0x94D049BB133111EBULL;
  h ^= h >> 31;

  const double u = static_cast<double>(h >> 11) / 9007199254740992.0;  // [0,1)
  return (u < so->fracao) ? 1 : 0;
}

}  // namespace

// ---------------------------------------------------------------------------
// CONSTRUCAO
// ---------------------------------------------------------------------------

Malha3::Malha3(const std::vector<Bloco> &blocos, sc_MPI_Comm comm)
    : blocos_(blocos), comm_(comm)
{
  sc_MPI_Comm_rank(comm, &rank_);
  sc_MPI_Comm_size(comm, &nprocs_);

  const long n = static_cast<long>(blocos.size());
  if (n == 0) throw std::invalid_argument("Malha3: lista de blocos vazia");

  for (long b = 0; b < n; b++) {
    const Bloco &k = blocos[static_cast<size_t>(b)];
    for (int d = 0; d < 3; d++)
      if (!(k.hi[d] > k.lo[d] + EPS))
        throw std::invalid_argument(
            "Malha3: bloco " + std::to_string(b) + " degenerado na direcao " +
            std::to_string(d) + ": lo=" + std::to_string(k.lo[d]) +
            " hi=" + std::to_string(k.hi[d]));

    // ARVORE DE UMA CELULA E' RECUSADA.  Uma celula nao tem interior -- as SEIS
    // facetas dela sao de fronteira --, e nenhum estencil centrado nela encontra
    // vizinho dentro do proprio bloco.  O levantamento do t8code registra isto
    // como defeito medido no sistema antecessor, e com o preco: "so' apareceu
    // como escoamento divergindo numa simulacao de meia hora".
    if (k.nivel < 1)
      throw std::invalid_argument(
          "Malha3: bloco " + std::to_string(b) + " com nivel " +
          std::to_string(k.nivel) + " tem UMA celula so', e uma celula nao tem "
          "interior -- as seis facetas dela sao de fronteira.  Use nivel >= 1");
  }

  for (long b = 0; b < n; b++)
    for (long c = b + 1; c < n; c++) {
      const Bloco &p = blocos[static_cast<size_t>(b)];
      const Bloco &q = blocos[static_cast<size_t>(c)];
      if (sobrepoe(p.lo[0], p.hi[0], q.lo[0], q.hi[0]) &&
          sobrepoe(p.lo[1], p.hi[1], q.lo[1], q.hi[1]) &&
          sobrepoe(p.lo[2], p.hi[2], q.lo[2], q.hi[2]))
        throw std::invalid_argument("Malha3: blocos " + std::to_string(b) + " e " +
                                    std::to_string(c) + " se sobrepoem em volume");
    }

  // AS JUNCOES SAO VALIDADAS ANTES DE O `cmesh` EXISTIR, e isso nao e' estilo:
  // um `throw` depois do `t8_cmesh_init` VAZA o cmesh, e o `sc_finalize` detecta
  // -- "Memory balance (t8)" -- e ABORTA o programa.  O teste de recusa passava e
  // o processo morria depois, o que num arreio aparece como suite verde com
  // codigo de saida errado.
  //
  // Validar antes de alocar tambem e' o desenho melhor: a recusa nao precisa
  // desmontar nada.
  for (long b = 0; b < n; b++)
    for (long c = b + 1; c < n; c++) {
      const Bloco &p = blocos[static_cast<size_t>(b)];
      const Bloco &q = blocos[static_cast<size_t>(c)];
      for (int d = 0; d < 3; d++) {
        const int t1 = (d + 1) % 3;
        const int t2 = (d + 2) % 3;
        for (int sentido = 0; sentido < 2; sentido++) {
          const Bloco &a = (sentido == 0) ? p : q;
          const Bloco &z = (sentido == 0) ? q : p;
          const long ia = (sentido == 0) ? b : c;
          const long iz = (sentido == 0) ? c : b;
          if (!coincide(a.hi[d], z.lo[d])) continue;

          const bool casa_t1 = coincide(a.lo[t1], z.lo[t1]) && coincide(a.hi[t1], z.hi[t1]);
          const bool casa_t2 = coincide(a.lo[t2], z.lo[t2]) && coincide(a.hi[t2], z.hi[t2]);
          const bool toca = sobrepoe(a.lo[t1], a.hi[t1], z.lo[t1], z.hi[t1]) &&
                            sobrepoe(a.lo[t2], a.hi[t2], z.lo[t2], z.hi[t2]);
          if (toca && !(casa_t1 && casa_t2))
            throw std::invalid_argument(
                "Malha3: blocos " + std::to_string(ia) + " e " + std::to_string(iz) +
                " encostam na direcao " + std::to_string(d) +
                " com faces que se sobrepoem PARCIALMENTE; o join do t8code nao"
                " representa isso, e aceitar daria contorno no meio do dominio");
        }
      }
    }

  t8_cmesh_t cmesh;
  t8_cmesh_init(&cmesh);
  t8_cmesh_register_geometry(cmesh, t8_geometry_linear_new());

  for (long b = 0; b < n; b++) {
    const Bloco &k = blocos[static_cast<size_t>(b)];
    t8_cmesh_set_tree_class(cmesh, b, T8_ECLASS_HEX);

    // A ORDEM DOS OITO VERTICES E' A DE BITS: `x` mais rapido, depois `y`,
    // depois `z` -- isto e', o vertice `v` tem a coordenada `hi` na direcao `d`
    // quando o bit `d` de `v` esta' ligado.  E' a convencao do t8code.
    //
    // TROCAR DOIS VERTICES PRODUZ UM HEXAEDRO TORCIDO QUE COMMITA SEM ERRO, e
    // em tres dimensoes isso e' mais facil de errar que em duas: ha' oito
    // vertices e a intuicao de "sentido anti-horario" nao se aplica.  E' por
    // isso que o teste afirma que o centroide de cada folha cai dentro da caixa
    // do bloco dela.
    double v[24];
    for (int vi = 0; vi < 8; vi++) {
      v[3 * vi + 0] = (vi & 1) ? k.hi[0] : k.lo[0];
      v[3 * vi + 1] = (vi & 2) ? k.hi[1] : k.lo[1];
      v[3 * vi + 2] = (vi & 4) ? k.hi[2] : k.lo[2];
    }
    t8_cmesh_set_tree_vertices(cmesh, b, v, 8);
  }

  // AS JUNCOES.  Faces do HEX: 0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+.  Dois blocos
  // compartilham face quando um limite coincide e o RETANGULO transversal e' o
  // mesmo nas duas direcoes restantes.
  //
  // COINCIDENCIA PARCIAL E' RECUSADA, e em tres dimensoes ha' mais maneiras de
  // ela acontecer: as duas direcoes transversais podem casar uma e nao a outra.
  // O `join` do t8code nao representa isso, e aceitar em silencio daria contorno
  // no meio do dominio.
  for (long b = 0; b < n; b++)
    for (long c = b + 1; c < n; c++) {
      const Bloco &p = blocos[static_cast<size_t>(b)];
      const Bloco &q = blocos[static_cast<size_t>(c)];

      for (int d = 0; d < 3; d++) {
        const int t1 = (d + 1) % 3;
        const int t2 = (d + 2) % 3;
        const int face_mais = 2 * d + 1;
        const int face_menos = 2 * d;

        for (int sentido = 0; sentido < 2; sentido++) {
          const Bloco &a = (sentido == 0) ? p : q;
          const Bloco &z = (sentido == 0) ? q : p;
          const long ia = (sentido == 0) ? b : c;
          const long iz = (sentido == 0) ? c : b;

          if (!coincide(a.hi[d], z.lo[d])) continue;

          const bool casa_t1 = coincide(a.lo[t1], z.lo[t1]) && coincide(a.hi[t1], z.hi[t1]);
          const bool casa_t2 = coincide(a.lo[t2], z.lo[t2]) && coincide(a.hi[t2], z.hi[t2]);
          // A VALIDACAO JA' FOI FEITA acima, antes de o `cmesh` existir; aqui
          // so' se registra o que passou.
          if (casa_t1 && casa_t2) {
            t8_cmesh_set_join(cmesh, ia, iz, face_mais, face_menos, 0);
            juncoes_++;
          }
        }
      }
    }

  t8_cmesh_commit(cmesh, comm);

  int nivel_min = blocos[0].nivel;
  std::vector<int> alvo(static_cast<size_t>(n));
  for (long b = 0; b < n; b++) {
    alvo[static_cast<size_t>(b)] = blocos[static_cast<size_t>(b)].nivel;
    nivel_min = std::min(nivel_min, blocos[static_cast<size_t>(b)].nivel);
  }

  const t8_scheme_c *esquema = t8_scheme_new_default();

  // A FRANJA JA' NASCE AQUI (`do_face_ghost = 1`).  Sem ela, um vizinho de outro
  // processo some, a face passa por CONTORNO DO DOMINIO, e o sistema resolve
  // outro problema sem reclamar.
  f_ = t8_forest_new_uniform(cmesh, esquema, nivel_min, 1, comm);

  int nivel_max = nivel_min;
  for (int v : alvo) nivel_max = std::max(nivel_max, v);

  if (nivel_max > nivel_min) {
    // O REFINO VEM EM COMMIT SEPARADO de qualquer ghost ou particao: pedir
    // `set_adapt`, `set_ghost` e `set_partition` no mesmo commit do t8code
    // derruba com np>1 e FUNCIONA em np=1 -- armadilha medida em duas dimensoes,
    // e das piores, porque o teste serial nao a ve.
    Alvos ctx{&alvo};
    t8_forest_t novo;
    t8_forest_init(&novo);
    t8_forest_set_user_data(novo, &ctx);
    t8_forest_set_adapt(novo, f_, adapt_por_arvore, 1 /* recursivo */);
    t8_forest_set_ghost(novo, 1, T8_GHOST_FACES);
    t8_forest_commit(novo);
    f_ = novo;
  }

  mede_vizinhos();
  monta_franja();
}

Malha3::~Malha3()
{
  if (f_ != nullptr) t8_forest_unref(&f_);
}

void Malha3::mede_vizinhos()
{
  // MEDE, nao presume.  O teto teorico e' 24 -- seis faces vezes quatro --, e o
  // maximo real medido e' 9 em malha adaptada e 6 em uniforme, porque uma celula
  // tem no maximo uma face com vizinhos refinados em cada direcao.
  //
  // Dimensionar a carga da franja pelo teto custaria 2,7 vezes mais memoria: em
  // tres dimensoes, a diferenca entre 55 KB e 24 KB por folha.
  long maior = 0;
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      long deste = 0;
      for (int fa = 0; fa < 6; fa++)
        deste += static_cast<long>(vizinhos_crus(f_, it, ie, fa).size());
      maior = std::max(maior, deste);
    }
  }

  // O MAXIMO E' GLOBAL, e tem de ser: a carga publicada por cada processo precisa
  // caber no buffer de TODOS, senao o que um publica o outro nao consegue ler.
  if (nprocs_ > 1) {
    long r = maior;
    sc_MPI_Allreduce(&maior, &r, 1, sc_MPI_LONG, sc_MPI_MAX, comm_);
    maior = r;
  }
  max_viz_ = static_cast<int>(std::max(6L, maior));
}

// ---------------------------------------------------------------------------
// CONSULTAS
// ---------------------------------------------------------------------------

const Bloco &Malha3::bloco(long b) const
{
  if (b < 0 || b >= num_blocos())
    throw std::out_of_range("Malha3::bloco: " + std::to_string(b) + " fora de [0, " +
                            std::to_string(num_blocos()) + ")");
  return blocos_[static_cast<size_t>(b)];
}

long Malha3::num_folhas() const
{
  return static_cast<long>(t8_forest_get_local_num_leaf_elements(f_));
}

long Malha3::num_folhas_do_bloco(long b) const
{
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha3::num_folhas_do_bloco: arvore " +
                            std::to_string(b) + " nao e' local");
  return static_cast<long>(
      t8_forest_get_tree_num_leaf_elements(f_, static_cast<t8_locidx_t>(b)));
}

long Malha3::num_folhas_globais() const
{
  return static_cast<long>(t8_forest_get_global_num_leaf_elements(f_));
}

long Malha3::primeiro_id() const
{
  return static_cast<long>(t8_forest_get_first_local_leaf_element_id(f_));
}

long Malha3::id_global(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha3::id_global: arvore " + std::to_string(b) +
                            " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha3::id_global: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));
  // O deslocamento da arvore dentro da floresta local mais o indice na arvore.
  long antes = 0;
  for (t8_locidx_t t = 0; t < it; t++)
    antes += static_cast<long>(t8_forest_get_tree_num_leaf_elements(f_, t));
  return primeiro_id() + antes + i;
}

Ponto3 Malha3::centroide(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha3::centroide: arvore " + std::to_string(b) +
                            " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha3::centroide: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));
  const t8_element_t *e =
      t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));
  double c[3] = {0.0, 0.0, 0.0};
  t8_forest_element_centroid(f_, it, e, c);
  return Ponto3{c[0], c[1], c[2]};
}

std::vector<Ponto3> Malha3::centroides() const
{
  std::vector<Ponto3> v;
  v.reserve(static_cast<size_t>(num_folhas()));
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) {
      const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, ie);
      double c[3] = {0.0, 0.0, 0.0};
      t8_forest_element_centroid(f_, it, e, c);
      v.push_back(Ponto3{c[0], c[1], c[2]});
    }
  }
  return v;
}

int Malha3::nivel(long b, long i) const
{
  const t8_locidx_t it = static_cast<t8_locidx_t>(b);
  if (b < 0 || b >= static_cast<long>(t8_forest_get_num_local_trees(f_)))
    throw std::out_of_range("Malha3::nivel: arvore " + std::to_string(b) +
                            " nao e' local");
  if (i < 0 || i >= num_folhas_do_bloco(b))
    throw std::out_of_range("Malha3::nivel: folha " + std::to_string(i) +
                            " fora da arvore " + std::to_string(b));
  const t8_element_t *e =
      t8_forest_get_leaf_element_in_tree(f_, it, static_cast<t8_locidx_t>(i));
  const t8_scheme_c *esq = t8_forest_get_scheme(f_);
  return esq->element_get_level(t8_forest_get_tree_class(f_, it), e);
}

Ponto3 Malha3::lados_da_folha(long b, long i) const
{
  const Bloco &k = bloco(b);
  const int div = 1 << nivel(b, i);
  return Ponto3{(k.hi[0] - k.lo[0]) / div, (k.hi[1] - k.lo[1]) / div,
                (k.hi[2] - k.lo[2]) / div};
}

double Malha3::volume(long b, long i) const
{
  const Ponto3 h = lados_da_folha(b, i);
  return h[0] * h[1] * h[2];
}

int Malha3::nivel_minimo() const
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
  // PROPRIEDADE DA MALHA INTEIRA, logo reduz.  `v < 0` e' particao VAZIA, e e'
  // neutralizada para o maximo inteiro antes do minimo, senao a particao vazia
  // viraria o minimo global.
  if (nprocs_ > 1) {
    int envia = (v < 0) ? std::numeric_limits<int>::max() : v;
    int r = envia;
    sc_MPI_Allreduce(&envia, &r, 1, sc_MPI_INT, sc_MPI_MIN, comm_);
    v = (r == std::numeric_limits<int>::max()) ? -1 : r;
  }
  return v;
}

int Malha3::nivel_maximo() const
{
  int v = -1;
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);
  for (t8_locidx_t it = 0; it < narv; it++) {
    const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
    for (t8_locidx_t ie = 0; ie < ne; ie++) v = std::max(v, nivel(it, ie));
  }
  if (nprocs_ > 1) {
    int r = v;
    sc_MPI_Allreduce(&v, &r, 1, sc_MPI_INT, sc_MPI_MAX, comm_);
    v = r;
  }
  return v;
}

int Malha3::num_vizinhos(long b, long i, int fa) const
{
  return static_cast<int>(ids_vizinhos(b, i, fa).size());
}

std::vector<long> Malha3::ids_vizinhos(long b, long i, int fa) const
{
  if (fa < 0 || fa > 5)
    throw std::out_of_range("Malha3::ids_vizinhos: face " + std::to_string(fa) +
                            " -- o HEX tem 6 (0:x-, 1:x+, 2:y-, 3:y+, 4:z-, 5:z+)");

  // A SOMA `primeiro_id() + indice` ESTA' ERRADA EM PARALELO, e e' o defeito mais
  // caro que a versao de duas dimensoes teve: ela vale para as folhas LOCAIS e
  // nao para as de FRANJA, cujo id pertence a outro processo.  Somar o
  // deslocamento deste processo a um indice de franja da' o id de uma folha
  // local que EXISTE -- so' que outra.  Em np=1 nao ha' franja e a soma acerta
  // sempre, e e' por isso que o caminho serial nunca acusaria.
  const std::vector<long> crus =
      vizinhos_crus(f_, static_cast<t8_locidx_t>(b), static_cast<t8_locidx_t>(i), fa);
  std::vector<long> ids;
  ids.reserve(crus.size());
  for (long x : crus) ids.push_back(id_global_do_indice(x));
  return ids;
}

std::vector<int> Malha3::niveis_vizinhos(long b, long i, int fa) const
{
  if (fa < 0 || fa > 5)
    throw std::out_of_range("Malha3::niveis_vizinhos: face " + std::to_string(fa) +
                            " -- o HEX tem 6");
  const std::vector<long> crus =
      vizinhos_crus(f_, static_cast<t8_locidx_t>(b), static_cast<t8_locidx_t>(i), fa);
  std::vector<int> niveis;
  niveis.reserve(crus.size());
  for (long x : crus) niveis.push_back(nivel_do_indice(x));
  return niveis;
}

long Malha3::num_folhas_em_transicao() const
{
  long n = 0;
  for (long b = 0; b < num_blocos(); b++)
    for (long i = 0; i < num_folhas_do_bloco(b); i++) {
      const int nv = nivel(b, i);
      bool transicao = false;
      for (int fa = 0; fa < 6 && !transicao; fa++)
        for (int viz : niveis_vizinhos(b, i, fa))
          if (viz != nv) { transicao = true; break; }
      if (transicao) n++;
    }
  // CONTAGEM DA MALHA INTEIRA: um processo cuja particao caiu toda na regiao
  // uniforme contaria zero, e o teste que afirma "a malha TEM transicao"
  // passaria por nao medir nada.
  if (nprocs_ > 1) {
    long r = n;
    sc_MPI_Allreduce(&n, &r, 1, sc_MPI_LONG, sc_MPI_SUM, comm_);
    n = r;
  }
  return n;
}

// ---------------------------------------------------------------------------
// REFINO
// ---------------------------------------------------------------------------

void Malha3::refina_caixa(const Ponto3 &lo, const Ponto3 &hi, int nivel_alvo)
{
  for (int d = 0; d < 3; d++)
    if (!(hi[d] > lo[d] + EPS))
      throw std::invalid_argument("Malha3::refina_caixa: caixa degenerada na "
                                  "direcao " + std::to_string(d));

  // A caixa tem de ALCANCAR algum bloco E pedir mais nivel do que ele ja' tem.
  // Refino que nao refina nada e' quase sempre engano de quem chamou, e aceito
  // calado deixaria a malha PARECER refinada.
  bool alcanca = false, aprofunda = false;
  for (const Bloco &k : blocos_) {
    bool toca = true;
    for (int d = 0; d < 3; d++)
      if (!sobrepoe(k.lo[d], k.hi[d], lo[d], hi[d])) toca = false;
    if (toca) {
      alcanca = true;
      if (nivel_alvo > k.nivel) aprofunda = true;
    }
  }
  if (!alcanca)
    throw std::invalid_argument("Malha3::refina_caixa: a caixa nao alcanca bloco algum");
  if (!aprofunda)
    throw std::invalid_argument(
        "Malha3::refina_caixa: nivel alvo " + std::to_string(nivel_alvo) +
        " nao e' maior que o dos blocos que a caixa alcanca -- isto nao refinaria "
        "nada, e aceitar calado deixaria a malha parecer refinada");

  Caixa cx{lo, hi, nivel_alvo};
  t8_forest_t refinada;
  t8_forest_init(&refinada);
  t8_forest_set_user_data(refinada, &cx);
  t8_forest_set_adapt(refinada, f_, adapt_caixa, 1);
  t8_forest_commit(refinada);

  t8_forest_t equilibrada;
  t8_forest_init(&equilibrada);
  t8_forest_set_balance(equilibrada, refinada, 0);
  t8_forest_set_ghost(equilibrada, 1, T8_GHOST_FACES);
  t8_forest_commit(equilibrada);

  // Sem `unref`: o commit CONSOME a floresta de origem.  Liberar de novo deixa a
  // malha VAZIA em silencio -- medido em duas dimensoes, e o que pegou foi a
  // conservacao de volume, nao a contagem de folhas.
  f_ = equilibrada;
  mede_vizinhos();
  monta_franja();
}

void Malha3::refina_aleatorio(double fracao, unsigned long semente, int nivel_max)
{
  if (fracao <= 0.0 || fracao > 1.0)
    throw std::invalid_argument("Malha3::refina_aleatorio: fracao " +
                                std::to_string(fracao) + " fora de (0, 1]");
  if (nivel_max <= nivel_maximo())
    throw std::invalid_argument(
        "Malha3::refina_aleatorio: nivel_max " + std::to_string(nivel_max) +
        " nao e' maior que o nivel atual " + std::to_string(nivel_maximo()) +
        " -- isto nao refinaria nada");

  Sorteio so{fracao, semente, nivel_max};
  t8_forest_t refinada;
  t8_forest_init(&refinada);
  t8_forest_set_user_data(refinada, &so);
  t8_forest_set_adapt(refinada, f_, adapt_aleatorio, 0 /* NAO recursivo */);
  t8_forest_commit(refinada);

  t8_forest_t equilibrada;
  t8_forest_init(&equilibrada);
  t8_forest_set_balance(equilibrada, refinada, 0);
  t8_forest_set_ghost(equilibrada, 1, T8_GHOST_FACES);
  t8_forest_commit(equilibrada);

  f_ = equilibrada;
  mede_vizinhos();
  monta_franja();
}


// ---------------------------------------------------------------------------
// A FRANJA DE K CAMADAS
// ---------------------------------------------------------------------------

void Malha3::monta_franja() { monta_franja_camadas(camadas_); }

void Malha3::exige_camadas(int k)
{
  if (k < 1)
    throw std::invalid_argument("Malha3::exige_camadas: " + std::to_string(k) +
                                " -- a franja tem pelo menos uma camada");
  if (k <= camadas_) return;      // idempotente: pedir menos nao desfaz

  // SEM FANTASMAS NAO HA' NADA A TROCAR.  A franja de um processo sozinho e'
  // vazia em qualquer numero de camadas, e refazer a fase A custaria o percurso
  // inteiro para nao trocar dado nenhum -- medido em duas dimensoes: 22 ms numa
  // malha de 4096 celulas, 26% do tempo de resolver.
  if (t8_forest_get_num_ghosts(f_) == 0) {
    camadas_ = k;
    return;
  }
  monta_franja_camadas(k);
}

void Malha3::monta_franja_camadas(int k)
{
  camadas_ = k;
  saturacoes_ = 0;
  bytes_pub_ = 0;
  n_locais_ = static_cast<long>(t8_forest_get_local_num_leaf_elements(f_));
  const long nf = static_cast<long>(t8_forest_get_num_ghosts(f_));
  const long p0 = primeiro_id();

  struct Folha {
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double hx = 0.0, hy = 0.0, hz = 0.0;
    int nivel = -1;
    std::array<std::vector<long>, 6> viz;
    std::array<std::vector<int>, 6> viz_nivel;
    unsigned char faces = 0;     //!< bit por face conhecida; 0x3F = as seis
    bool completa() const { return faces == 0x3F; }
  };
  std::unordered_map<long, Folha> tab;

  // =========================================================================
  // FASE A -- identidade e geometria da camada 1
  // =========================================================================
  //
  // Resolve o problema do ovo e da galinha: os vizinhos de franja de uma folha
  // local sao conhecidos por INDICE, e o id deles so' chega aqui.
  const t8_scheme_c *esq = t8_forest_get_scheme(f_);
  const t8_locidx_t narv = t8_forest_get_num_local_trees(f_);

  std::vector<long> id_por_indice(static_cast<size_t>(n_locais_ + nf), -1);
  std::vector<std::array<std::vector<long>, 6>> viz_cru(static_cast<size_t>(n_locais_));
  {
    sc_array_t *arr = sc_array_new_count(8 * sizeof(double),
                                         static_cast<size_t>(n_locais_ + nf));
    bytes_pub_ += 8 * static_cast<long>(sizeof(double)) * n_locais_;

    long kk = 0;
    for (t8_locidx_t it = 0; it < narv; it++) {
      const t8_locidx_t ne = t8_forest_get_tree_num_leaf_elements(f_, it);
      for (t8_locidx_t ie = 0; ie < ne; ie++) {
        const t8_element_t *e = t8_forest_get_leaf_element_in_tree(f_, it, ie);
        double c[3] = {0.0, 0.0, 0.0};
        t8_forest_element_centroid(f_, it, e, c);
        const int nv = esq->element_get_level(t8_forest_get_tree_class(f_, it), e);

        // OS LADOS VEM DO BLOCO E DO NIVEL, e nao do diametro do elemento: a
        // mesma formula aqui e em `lados_da_folha` e' o que impede que a folha
        // local e a MESMA folha vista como franja tenham medidas ligeiramente
        // diferentes.  Uma diferenca de arredondamento ali apareceria como
        // esquema nao conservativo.
        const Bloco &kb = blocos_[static_cast<size_t>(it)];
        const int div = 1 << nv;

        const long id = p0 + kk;
        id_por_indice[static_cast<size_t>(kk)] = id;
        Folha &fo = tab[id];
        fo.cx = c[0]; fo.cy = c[1]; fo.cz = c[2];
        fo.hx = (kb.hi[0] - kb.lo[0]) / div;
        fo.hy = (kb.hi[1] - kb.lo[1]) / div;
        fo.hz = (kb.hi[2] - kb.lo[2]) / div;
        fo.nivel = nv;
        for (int fa = 0; fa < 6; fa++)
          viz_cru[static_cast<size_t>(kk)][static_cast<size_t>(fa)] =
              vizinhos_crus(f_, it, ie, fa);

        double *d = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(kk)));
        d[0] = static_cast<double>(id);
        d[1] = fo.cx; d[2] = fo.cy; d[3] = fo.cz;
        d[4] = fo.hx; d[5] = fo.hy; d[6] = fo.hz;
        d[7] = static_cast<double>(nv);
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
        fo.cx = d[1]; fo.cy = d[2]; fo.cz = d[3];
        fo.hx = d[4]; fo.hy = d[5]; fo.hz = d[6];
        fo.nivel = static_cast<int>(d[7]);
      }
    }
    sc_array_destroy(arr);
  }

  // Agora os vizinhos das LOCAIS podem virar ids, e as SEIS faces delas sao
  // conhecidas -- o t8code as respondeu todas, logo lista vazia aqui significa
  // contorno do dominio, sem ambiguidade.
  for (long i = 0; i < n_locais_; i++) {
    Folha &fo = tab[id_por_indice[static_cast<size_t>(i)]];
    for (int fa = 0; fa < 6; fa++)
      for (long x : viz_cru[static_cast<size_t>(i)][static_cast<size_t>(fa)]) {
        const long id = (x >= 0 && x < n_locais_ + nf)
                            ? id_por_indice[static_cast<size_t>(x)] : -1;
        if (id < 0) continue;
        fo.viz[static_cast<size_t>(fa)].push_back(id);
        const auto itv = tab.find(id);
        fo.viz_nivel[static_cast<size_t>(fa)].push_back(
            (itv != tab.end() && itv->second.nivel >= 0) ? itv->second.nivel : -1);
      }
    fo.faces = 0x3F;
  }

  // =========================================================================
  // FASE B -- a caixa GEOMETRICA, numa rodada so'
  // =========================================================================
  //
  // O QUE O SUPORTE PRECISA E' VIZINHANCA GEOMETRICA, e nao topologica, e
  // confundi-las custou caro: a celula diagonal `(2,2,2)` esta' a raio
  // geometrico 2 e a SEIS faces de distancia.  Propagando por faces, cobri-la
  // exigia cinco camadas e 200,7 KB por folha publicada -- medido.
  //
  // UMA RODADA BASTA, e a razao e' simples quando se olha de quem publica para
  // quem recebe: o dono de um fantasma meu tem como LOCAIS justamente as
  // celulas que me faltam.  Se ele publica a caixa de raio `R` em torno de cada
  // folha dele, eu recebo -- no lugar do fantasma `F`, que esta' a uma face de
  // mim -- a caixa de raio `R` em torno de `F`.
  //
  // E ela cobre a minha: a celula mais distante da minha caixa de raio `R` e'
  // `(R,R,R)`, e a distancia dela a `F = (1,0,0)` em norma do maximo e'
  // `max(R-1, R, R) = R`.  Esta' na caixa de `F`.
  //
  // PUBLICO O QUE CONHECO, nao o que preciso -- nao ha' circularidade.  Depois
  // da fase A conheco minhas locais e a camada 1; a caixa de raio `R` em torno
  // de uma local minha inclui celulas que conheco e outras que nao, e publico as
  // que conheco.  Quem recebe precisa justamente das que sao LOCAIS para mim.
  if (nf > 0 && k > 1) {
    // 8 PROPRIOS, 1 DE MASCARA, e (id+nivel) por vizinho.
    //
    // O `9` ja' foi `8`: a mascara de faces entrou no campo 8 sem que a conta
    // crescesse, e com isso cabiam `max_viz_ - 1` vizinhos em vez de
    // `max_viz_`.  Numa malha uniforme, onde cada folha tem exatamente seis
    // vizinhos, o SEXTO nao cabia -- e o sexto e' a face 5, a ultima do laco.
    // A folha chegava com a mascara dizendo "conheco as seis" e a lista sem a
    // face 5, de modo que o fechamento lia "face sem vizinho numa celula que
    // conhece todas as faces" e aplicava DIRICHLET NO MEIO DO DOMINIO.
    const int CAMPOS = 9 + max_viz_ * 2;

    // O RAIO E' `k - 1`: `exige_camadas(k)` pede cobertura de raio geometrico
    // `k - 1` em torno de cada folha local.  A conversao mantem a interface
    // antiga -- quem pedia 3 camadas topologicas pede 3 e recebe raio 2, que e'
    // o que o suporte por caixa de raio 2 precisa.
    const int R = k - 1;

    // INDICE ESPACIAL TEMPORARIO sobre o que ja' se conhece.  O definitivo so'
    // existe depois do achatamento, e a caixa precisa dele agora.
    std::unordered_map<unsigned long long, long> onde;
    auto chave_de_folha = [&](double cx, double cy, double cz, int nv) -> unsigned long long {
      long ib = -1;
      for (long bb = 0; bb < static_cast<long>(blocos_.size()); bb++) {
        const Bloco &kk = blocos_[static_cast<size_t>(bb)];
        bool dentro = true;
        for (int d = 0; d < 3; d++) {
          const double c[3] = {cx, cy, cz};
          if (c[d] < kk.lo[d] - EPS || c[d] > kk.hi[d] + EPS) dentro = false;
        }
        if (dentro) { ib = bb; break; }
      }
      if (ib < 0) return ~0ULL;
      const Bloco &kb = blocos_[static_cast<size_t>(ib)];
      const int div = 1 << nv;
      const long i = static_cast<long>((cx - kb.lo[0]) / ((kb.hi[0] - kb.lo[0]) / div));
      const long j = static_cast<long>((cy - kb.lo[1]) / ((kb.hi[1] - kb.lo[1]) / div));
      const long kq = static_cast<long>((cz - kb.lo[2]) / ((kb.hi[2] - kb.lo[2]) / div));
      return chave_celula(ib, nv, i, j, kq);
    };
    // RECONSTRUIDO A CADA RODADA: a tabela cresce, e um indice de uma rodada
    // anterior faria a rodada seguinte republicar o que ja' sabe e nao ver o
    // que acabou de aprender.
    auto reindexa = [&]() {
      onde.clear();
      nivel_local_max_ = -1;
      for (const auto &kv : tab)
        if (kv.second.nivel >= 0) {
          const unsigned long long ch =
              chave_de_folha(kv.second.cx, kv.second.cy, kv.second.cz, kv.second.nivel);
          if (ch != ~0ULL) onde[ch] = kv.first;
          if (kv.second.nivel > nivel_local_max_) nivel_local_max_ = kv.second.nivel;
        }
    };
    reindexa();

    // Busca geometrica sobre o indice temporario, com a MESMA regra semiaberta
    // da busca definitiva -- duas regras produziriam duas malhas.
    auto acha = [&](double px, double py, double pz) -> long {
      const double p[3] = {px, py, pz};
      for (int nv = nivel_local_max_ >= 0 ? nivel_local_max_ : 30; nv >= 0; nv--) {
        for (long ib = 0; ib < static_cast<long>(blocos_.size()); ib++) {
          const Bloco &kk = blocos_[static_cast<size_t>(ib)];
          bool no_bloco = true;
          for (int d = 0; d < 3; d++)
            if (p[d] < kk.lo[d] - EPS || p[d] > kk.hi[d] + EPS) no_bloco = false;
          if (!no_bloco) continue;
          const int div = 1 << nv;
          long c[3];
          bool fora = false;
          for (int d = 0; d < 3; d++) {
            const double hh = (kk.hi[d] - kk.lo[d]) / div;
            long q = static_cast<long>((p[d] - kk.lo[d] + EPS) / hh);
            if (q >= div) q = div - 1;
            if (q < 0) { fora = true; break; }
            c[d] = q;
          }
          if (fora) continue;
          const auto it = onde.find(chave_celula(ib, nv, c[0], c[1], c[2]));
          if (it != onde.end()) return it->second;
        }
      }
      return -1;
    };

    const int MAXV = (2 * R + 1) * (2 * R + 1) * (2 * R + 1);
    const size_t passo = static_cast<size_t>(MAXV) * static_cast<size_t>(CAMPOS);

    // QUANTOS PONTOS DA CAIXA DE RAIO `R` AINDA FALTAM, somado sobre as folhas
    // locais.  E' o criterio de parada das rodadas, e e' o MESMO criterio que o
    // estencil usa depois para contar truncamento -- medir uma coisa aqui e
    // cobrar outra la' foi o que deixou passar a falha de quina.
    auto falta_na_caixa = [&]() -> long {
      long falta = 0;
      for (long i = 0; i < n_locais_; i++) {
        const Folha &me = tab[id_por_indice[static_cast<size_t>(i)]];
        for (int ii = -R; ii <= R; ii++)
          for (int jj = -R; jj <= R; jj++)
            for (int kk2 = -R; kk2 <= R; kk2++) {
              const Ponto3 p{me.cx + ii * me.hx, me.cy + jj * me.hy,
                             me.cz + kk2 * me.hz};
              if (!dentro_do_dominio(p)) continue;   // fora e' legitimo
              if (acha(p[0], p[1], p[2]) < 0) falta++;
            }
      }
      return falta;
    };

    sc_array_t *arr = sc_array_new_count(passo * sizeof(double),
                                         static_cast<size_t>(n_locais_ + nf));

    // AS RODADAS SAO ITERADAS ATE' COBRIR, e nao fixadas em uma.
    //
    // UMA RODADA SO' ALCANCA QUEM E' VIZINHO DE FACE.  A franja nativa do
    // t8code e' `T8_GHOST_FACES` -- a documentacao diz que so' esse tipo e'
    // suportado --, de modo que o processo diagonalmente oposto numa quina nao
    // e' vizinho de face de ninguem que fale comigo, e uma rodada nunca traz a
    // celula dele.  Em 8 processos isso aparecia como 192 suportes truncados no
    // centro do dominio, onde as oito particoes se tocam, e AUMENTAR AS CAMADAS
    // NAO RESOLVIA: a informacao nao viajava mais longe, ela nao viajava por
    // aquele caminho.  Cada rodada avanca um salto de POSSE, e tres saltos
    // bastam para a quina de oito.
    //
    // A PARADA E' MEDIDA E COLETIVA: para quando nenhuma folha local tem ponto
    // de caixa faltando em processo nenhum.  Um numero fixo de rodadas estaria
    // certo para esta particao e errado para a proxima.
    const int TETO = 16;   // rede de seguranca; chegar nele e' defeito a apurar
    rodadas_franja_ = 0;
    faltando_franja_ = 0;
    for (int rodada = 0; rodada < TETO; rodada++) {
      reindexa();
      long falta = falta_na_caixa(), falta_g = 0;
      sc_MPI_Allreduce(&falta, &falta_g, 1, sc_MPI_LONG, sc_MPI_SUM, comm_);
      faltando_franja_ = falta_g;
      if (falta_g == 0) break;
      rodadas_franja_++;
      bytes_pub_ += static_cast<long>(passo * sizeof(double)) * n_locais_;

    for (long i = 0; i < n_locais_ + nf; i++) {
      double *b = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(i)));
      for (size_t z = 0; z < passo; z++) b[z] = -1.0;
    }

    for (long i = 0; i < n_locais_; i++) {
      double *b = static_cast<double *>(sc_array_index_int(arr, static_cast<int>(i)));
      const long eu = id_por_indice[static_cast<size_t>(i)];
      const Folha &me = tab[eu];

      int slot = 0;
      bool descartou = false;
      std::set<long> posto;
      for (int ii = -R; ii <= R; ii++)
        for (int jj = -R; jj <= R; jj++)
          for (int kk2 = -R; kk2 <= R; kk2++) {
            const long id = acha(me.cx + ii * me.hx, me.cy + jj * me.hy,
                                 me.cz + kk2 * me.hz);
            if (id < 0 || !posto.insert(id).second) continue;
            const auto it = tab.find(id);
            if (it == tab.end() || it->second.nivel < 0) continue;

            // SATURACAO E' PONTO DESCARTADO, e nao ultimo slot usado.
            //
            // A primeira versao contava `slot >= MAXV` depois do laco, o que
            // marcava saturacao sempre que a caixa preenchia exatamente os
            // slots -- 216 falsos positivos com raio 1, onde a caixa tem 27
            // pontos e `MAXV` e' 27.  Um contador que acusa quando nada se
            // perdeu desacredita o oraculo inteiro.
            if (slot >= MAXV) { descartou = true; continue; }
            const Folha &fo = it->second;

            double *c = b + static_cast<size_t>(slot) * CAMPOS;
            c[0] = static_cast<double>(id);
            c[1] = fo.cx; c[2] = fo.cy; c[3] = fo.cz;
            c[4] = fo.hx; c[5] = fo.hy; c[6] = fo.hz;
            c[7] = static_cast<double>(fo.nivel);
            // A MASCARA DE FACES CONHECIDAS, no primeiro campo do bloco de
            // vizinhos.  Sem ela, uma face sem vizinhos na lista e' ambigua:
            // contorno do dominio, ou "nao me contaram"?  Confundi-las aplica
            // condicao de contorno no MEIO do dominio.
            // A MASCARA SO' AFIRMA AS FACES QUE COUBERAM INTEIRAS.
            //
            // Antes ela era escrita de uma vez, antes da carga, e a carga era
            // truncada por capacidade -- de modo que a folha viajava afirmando
            // conhecer uma face cuja lista de vizinhos tinha ficado para tras.
            // Quem recebe nao tem como desconfiar: lista vazia mais bit ligado
            // e' exatamente a assinatura de uma parede do dominio.
            //
            // UMA FACE PARCIAL E' PIOR QUE UMA FACE DESCONHECIDA: desconhecida,
            // o fechamento se RECUSA a aplicar condicao de contorno e o esquema
            // para; parcial, ele aplica a condicao errada e o resultado sai
            // plausivel.  Por isso o bit cai junto com a lista.
            int z = 9;
            unsigned char mascara = 0;
            for (int fa = 0; fa < 6; fa++) {
              const size_t nv = fo.viz[static_cast<size_t>(fa)].size();
              if (!(fo.faces & (1u << fa))) continue;
              if (z + 2 * static_cast<int>(nv) > CAMPOS) { descartou = true; continue; }
              for (size_t q = 0; q < nv; q++) {
                c[z++] = static_cast<double>(fo.viz[static_cast<size_t>(fa)][q] * 6 + fa);
                c[z++] = static_cast<double>(fo.viz_nivel[static_cast<size_t>(fa)][q]);
              }
              mascara |= static_cast<unsigned char>(1u << fa);
            }
            c[8] = static_cast<double>(mascara);
            slot++;
          }
      if (descartou) saturacoes_++;
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
        // A FOLHA LOCAL E' A VERSAO AUTORITATIVA DELA MESMA.
        if (fo.nivel < 0) {
          fo.cx = c[1]; fo.cy = c[2]; fo.cz = c[3];
          fo.hx = c[4]; fo.hy = c[5]; fo.hz = c[6];
          fo.nivel = static_cast<int>(c[7]);
        }
        if (!fo.completa()) {
          const unsigned char mk = static_cast<unsigned char>(c[8]);
          if (mk != 0) {
            std::array<std::vector<long>, 6> nv;
            std::array<std::vector<int>, 6> nn;
            for (int z = 9; z + 1 < CAMPOS; z += 2) {  // pares (id*6+face, nivel)
              if (c[z] < 0.0) continue;
              const long cod = static_cast<long>(c[z]);
              const size_t fa = static_cast<size_t>(cod % 6);
              nv[fa].push_back(cod / 6);
              nn[fa].push_back(static_cast<int>(c[z + 1]));
            }
            fo.viz = nv;
            fo.viz_nivel = nn;
            fo.faces = mk;
          }
        }
      }
    }
    }   // fim da rodada
    sc_array_destroy(arr);
  }

  // =========================================================================
  // ACHATA a tabela em arrays indexados por indice de floresta
  // =========================================================================
  //
  // A ORDEM AQUI NAO E' LIVRE.  O t8code numera `0 .. num_local-1` para as
  // locais e `num_local .. num_local+num_ghosts-1` para as da CAMADA 1, e
  // ordenar TODA a franja por id -- que parece a escolha limpa -- faz o indice
  // do terceiro fantasma designar OUTRA folha.  Em duas dimensoes isso custou
  // 1,70e-01 de erro num esquema EXATO para campo linear.
  //
  // So' as camadas 2 e acima podem ser ordenadas por id: elas nao tem indice do
  // t8code, e a ordenacao as torna reprodutiveis entre execucoes.
  id_g_.clear(); cen_.clear(); lad_.clear(); niv_.clear();
  viz_.clear(); viz_niv_.clear(); expansivel_.clear();
  idx_por_id_.clear(); ponto_por_celula_.clear();

  auto chave_de = [&](double cx, double cy, double cz, int nv) -> unsigned long long {
    long ib = -1;
    const Ponto3 centro{cx, cy, cz};
    for (long bb = 0; bb < static_cast<long>(blocos_.size()); bb++) {
      const Bloco &kk = blocos_[static_cast<size_t>(bb)];
      bool dentro = true;
      for (int d = 0; d < 3; d++)
        if (centro[d] < kk.lo[d] - EPS || centro[d] > kk.hi[d] + EPS) dentro = false;
      if (dentro) { ib = bb; break; }
    }
    if (ib < 0) return ~0ULL;
    const Bloco &kb = blocos_[static_cast<size_t>(ib)];
    const int div = 1 << nv;
    const long i = static_cast<long>((cx - kb.lo[0]) / ((kb.hi[0] - kb.lo[0]) / div));
    const long j = static_cast<long>((cy - kb.lo[1]) / ((kb.hi[1] - kb.lo[1]) / div));
    const long kq = static_cast<long>((cz - kb.lo[2]) / ((kb.hi[2] - kb.lo[2]) / div));
    return chave_celula(ib, nv, i, j, kq);
  };

  auto empurra = [&](long id, const Folha &fo) {
    const long ix = static_cast<long>(id_g_.size());
    const unsigned long long ch = chave_de(fo.cx, fo.cy, fo.cz, fo.nivel);
    if (ch != ~0ULL) ponto_por_celula_[ch] = ix;
    idx_por_id_[id] = ix;
    id_g_.push_back(id);
    cen_.push_back(fo.cx); cen_.push_back(fo.cy); cen_.push_back(fo.cz);
    lad_.push_back(fo.hx); lad_.push_back(fo.hy); lad_.push_back(fo.hz);
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

  // OS NIVEIS LOCAIS, para a busca espacial.  Sem reducao, de proposito:
  // `id_do_ponto` e' chamada de caminhos que divergem entre os processos, e uma
  // coletiva ali da' `MPI_ERR_TRUNCATE` -- medido em duas dimensoes.
  nivel_local_min_ = -1;
  nivel_local_max_ = -1;
  for (int nv : niv_) {
    if (nivel_local_min_ < 0 || nv < nivel_local_min_) nivel_local_min_ = nv;
    nivel_local_max_ = std::max(nivel_local_max_, nv);
  }
}

// ---------------------------------------------------------------------------
// ACESSO POR INDICE
// ---------------------------------------------------------------------------

long Malha3::id_global_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(id_g_.size()))
    throw std::out_of_range("Malha3::id_global_do_indice: " + std::to_string(idx) +
                            " fora de [0, " + std::to_string(id_g_.size()) + ") -- " +
                            std::to_string(n_locais_) + " locais mais " +
                            std::to_string(n_fantasmas_) + " de franja");
  return id_g_[static_cast<size_t>(idx)];
}

Ponto3 Malha3::centroide_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(niv_.size()))
    throw std::out_of_range("Malha3::centroide_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  const size_t t = static_cast<size_t>(3 * idx);
  return Ponto3{cen_[t], cen_[t + 1], cen_[t + 2]};
}

Ponto3 Malha3::lados_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(niv_.size()))
    throw std::out_of_range("Malha3::lados_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  const size_t t = static_cast<size_t>(3 * idx);
  return Ponto3{lad_[t], lad_[t + 1], lad_[t + 2]};
}

int Malha3::nivel_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(niv_.size()))
    throw std::out_of_range("Malha3::nivel_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  return niv_[static_cast<size_t>(idx)];
}

long Malha3::indice_do_id(long id) const
{
  const auto it = idx_por_id_.find(id);
  return (it == idx_por_id_.end()) ? -1 : it->second;
}

std::vector<long> Malha3::vizinhos_do_indice(long idx, int fa) const
{
  if (idx < 0 || idx >= static_cast<long>(viz_.size()))
    throw std::out_of_range("Malha3::vizinhos_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  if (fa < 0 || fa > 5)
    throw std::out_of_range("Malha3::vizinhos_do_indice: face " + std::to_string(fa) +
                            " -- o HEX tem 6");
  return viz_[static_cast<size_t>(idx)][static_cast<size_t>(fa)];
}

std::vector<int> Malha3::niveis_vizinhos_do_indice(long idx, int fa) const
{
  if (idx < 0 || idx >= static_cast<long>(viz_niv_.size()))
    throw std::out_of_range("Malha3::niveis_vizinhos_do_indice: " +
                            std::to_string(idx) + " fora da faixa conhecida");
  if (fa < 0 || fa > 5)
    throw std::out_of_range("Malha3::niveis_vizinhos_do_indice: face " +
                            std::to_string(fa) + " -- o HEX tem 6");
  return viz_niv_[static_cast<size_t>(idx)][static_cast<size_t>(fa)];
}

bool Malha3::expansivel_do_indice(long idx) const
{
  if (idx < 0 || idx >= static_cast<long>(expansivel_.size()))
    throw std::out_of_range("Malha3::expansivel_do_indice: " + std::to_string(idx) +
                            " fora da faixa de folhas conhecidas");
  return expansivel_[static_cast<size_t>(idx)] != 0;
}

// ---------------------------------------------------------------------------
// BUSCA GEOMETRICA
// ---------------------------------------------------------------------------

bool Malha3::dentro_do_dominio(const Ponto3 &p) const
{
  for (const Bloco &k : blocos_) {
    bool dentro = true;
    for (int d = 0; d < 3; d++)
      if (p[d] < k.lo[d] - EPS || p[d] > k.hi[d] + EPS) dentro = false;
    if (dentro) return true;
  }
  return false;
}

long Malha3::id_do_ponto(const Ponto3 &p) const
{
  // CONSULTA PELO INDICE ESPACIAL: O(niveis), e nao O(folhas).  Uma varredura
  // linear aqui torna a montagem do operador QUADRATICA -- medido em duas
  // dimensoes, 8,74 s contra 0,679 s no mesmo problema, com o sintoma se
  // disfarcando de ganho de paralelismo.
  //
  // A SEMANTICA E' A SEMIABERTA `[lo, hi)`, e ela e' GLOBAL: nao depende de quem
  // pergunta.  E' dela que vem a conservacao do esquema de fluxo em paralelo.
  if (ponto_por_celula_.empty()) return -1;

  // OS NIVEIS SAO OS LOCAIS, e nao `nivel_minimo()`/`nivel_maximo()`: aqueles
  // sao COLETIVOS, e esta funcao e' chamada de lacos que divergem entre os
  // processos.  O sintoma foi `MPI_ERR_TRUNCATE` em duas dimensoes.
  const int nmin = nivel_local_min_;
  const int nmax = nivel_local_max_;
  if (nmin < 0 || nmax < 0) return -1;

  // DO MAIS FINO PARA O MAIS GROSSO: numa malha adaptada o ponto cai dentro de
  // celulas de varios niveis ao mesmo tempo, e a FOLHA e' a mais fina.
  for (int nv = nmax; nv >= nmin; nv--) {
    for (long ib = 0; ib < static_cast<long>(blocos_.size()); ib++) {
      const Bloco &k = blocos_[static_cast<size_t>(ib)];

      bool no_bloco = true;
      for (int d = 0; d < 3; d++)
        if (p[d] < k.lo[d] - EPS || p[d] > k.hi[d] + EPS) no_bloco = false;
      if (!no_bloco) continue;

      const int div = 1 << nv;
      long c[3];
      bool fora = false;
      for (int d = 0; d < 3; d++) {
        const double h = (k.hi[d] - k.lo[d]) / div;
        // A QUANTIZACAO REPRODUZ A REGRA SEMIABERTA: o `EPS` empurra o ponto que
        // esta' exatamente sobre uma faceta para a celula de coordenada MAIOR, e
        // no limite superior do dominio ele e' puxado de volta, para que a borda
        // pertenca a' ultima celula.
        long q = static_cast<long>((p[d] - k.lo[d] + EPS) / h);
        if (q >= div) q = div - 1;
        if (q < 0) { fora = true; break; }
        c[d] = q;
      }
      if (fora) continue;

      const auto it = ponto_por_celula_.find(chave_celula(ib, nv, c[0], c[1], c[2]));
      if (it != ponto_por_celula_.end()) return id_global_do_indice(it->second);
    }
  }
  return -1;
}

// ---------------------------------------------------------------------------
// PARTICAO
// ---------------------------------------------------------------------------

double Malha3::desbalanceamento() const
{
  if (nprocs_ <= 1) return 1.0;

  long n = num_folhas(), lmin = 0, lmax = 0;
  sc_MPI_Allreduce(&n, &lmin, 1, sc_MPI_LONG, sc_MPI_MIN, comm_);
  sc_MPI_Allreduce(&n, &lmax, 1, sc_MPI_LONG, sc_MPI_MAX, comm_);

  // PARTICAO VAZIA da' razao infinita, e e' a resposta certa: um processo sem
  // folha nenhuma nao esta' "bem balanceado", esta' ocioso.
  if (lmin <= 0) return (lmax > 0) ? std::numeric_limits<double>::infinity() : 1.0;
  return static_cast<double>(lmax) / static_cast<double>(lmin);
}

void Malha3::reparticiona()
{
  // PARTICAO EM COMMIT PROPRIO, com a franja junto.  O que o t8code nao aceita
  // e' `set_adapt` no mesmo commit -- e aqui nao ha' adapt.
  t8_forest_t novo;
  t8_forest_init(&novo);
  t8_forest_set_partition(novo, f_, 0);
  t8_forest_set_ghost(novo, 1, T8_GHOST_FACES);
  t8_forest_commit(novo);
  f_ = novo;
  mede_vizinhos();
  monta_franja();
}


}  // namespace t8flow3d
