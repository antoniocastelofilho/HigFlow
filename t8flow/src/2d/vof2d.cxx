#include "vof2d.hxx"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace t8flow2d {
namespace {

//! A REDUCAO CANONICA, e e' ela que dispensa a analise de casos por quadrante.
//
// Com `p` e `q` medidos a partir do canto de onde a normal SAI -- isto e'
// `p = 0` no lado para o qual `-n` aponta --, as duas componentes entram com o
// mesmo sinal e a regiao `{n.x >= d}` vira `{a p + b q >= t}` no quadrado
// unitario, com
//
//     a = |n0| h0,   b = |n1| h1,   t = d + (a + b)/2.
//
// O reflexo que leva `x` em `p` nao muda AREA, so' orientacao -- e por isso a
// formula nao depende dos sinais de `n`.  Quatro quadrantes, uma conta.
struct Canonico {
  double a;      //!< o MENOR dos dois, por construcao
  double b;      //!< o maior
  double vol;    //!< `h0 * h1`
  double soma;   //!< `a + b`
};

Canonico canoniza(const Ponto2 &n, const Ponto2 &h)
{
  const double nn = std::hypot(n[0], n[1]);
  if (!(nn > 0.0))
    throw std::invalid_argument(
        "plic: a normal e' nula -- nao ha' lado para a reta cortar");
  if (!(h[0] > 0.0) || !(h[1] > 0.0))
    throw std::invalid_argument("plic: lado de celula nao positivo");

  double a = std::fabs(n[0]) / nn * h[0];
  double b = std::fabs(n[1]) / nn * h[1];
  // `a <= b` SEMPRE.  A troca e' livre (a formula e' simetrica em `(a,b)`) e e'
  // o que faz o caso eixo-alinhado -- `a = 0` -- cair sozinho no ramo do meio,
  // em vez de precisar de um ramo proprio como no HiGFlow.
  if (a > b) std::swap(a, b);
  return {a, b, h[0] * h[1], a + b};
}

//! A fracao do quadrado unitario com `a p + b q <= t`.  Quadratica por pedacos,
//! crescente, e continua nas emendas `t = a` e `t = b`.
double fracao_abaixo(double t, double a, double b)
{
  const double S = a + b;
  if (t <= 0.0) return 0.0;
  if (t >= S) return 1.0;
  if (a <= 0.0) return t / b;            // reta paralela a um eixo
  if (t <= a) return t * t / (2.0 * a * b);
  if (t <= b) return (t - a / 2.0) / b;
  const double r = S - t;
  return 1.0 - r * r / (2.0 * a * b);
}

}  // namespace

double plic_area(const Ponto2 &n, const Ponto2 &h, double d)
{
  const Canonico c = canoniza(n, h);
  return c.vol * (1.0 - fracao_abaixo(d + c.soma / 2.0, c.a, c.b));
}

double plic_distancia(const Ponto2 &n, const Ponto2 &h, double area)
{
  const Canonico c = canoniza(n, h);
  if (!(area >= 0.0) || area > c.vol)
    throw std::invalid_argument(
        "plic_distancia: area fora de [0, h0*h1] -- a fracao volumetrica saiu "
        "do intervalo e isso tem de ser tratado (e CONTADO) por quem advecta, "
        "nao saturado aqui");

  // `g` e' a fracao ABAIXO da reta; inverte-se `fracao_abaixo`, ramo a ramo.
  const double g = 1.0 - area / c.vol;
  const double a = c.a, b = c.b;
  double t;
  if (a <= 0.0) {
    t = g * b;
  } else if (g <= a / (2.0 * b)) {        // `fracao_abaixo(a)` vale `a/(2b)`
    t = std::sqrt(2.0 * a * b * g);
  } else if (g <= 1.0 - a / (2.0 * b)) {  // e `fracao_abaixo(b)`, `1 - a/(2b)`
    t = b * g + a / 2.0;
  } else {
    t = c.soma - std::sqrt(2.0 * a * b * (1.0 - g));
  }
  return t - c.soma / 2.0;
}

// ---------------------------------------------------------------------------
// AS GEMEAS
// ---------------------------------------------------------------------------

double plic_area_por_recorte(const Ponto2 &n, const Ponto2 &h, double d)
{
  const double nn = std::hypot(n[0], n[1]);
  if (!(nn > 0.0))
    throw std::invalid_argument(
        "plic: a normal e' nula -- nao ha' lado para a reta cortar");
  if (!(h[0] > 0.0) || !(h[1] > 0.0))
    throw std::invalid_argument("plic: lado de celula nao positivo");
  const double nx = n[0] / nn, ny = n[1] / nn;

  // O retangulo da celula, em sentido ANTI-HORARIO para o sapateiro sair
  // positivo.
  std::vector<Ponto2> P = {{-h[0] / 2.0, -h[1] / 2.0},
                           {h[0] / 2.0, -h[1] / 2.0},
                           {h[0] / 2.0, h[1] / 2.0},
                           {-h[0] / 2.0, h[1] / 2.0}};

  // SUTHERLAND-HODGMAN CONTRA UM SEMIPLANO SO'.  `dentro` e' `n.x >= d`, que e'
  // a convencao declarada no cabecalho.
  auto val = [&](const Ponto2 &p) { return nx * p[0] + ny * p[1] - d; };
  std::vector<Ponto2> Q;
  Q.reserve(P.size() + 1);
  for (std::size_t k = 0; k < P.size(); k++) {
    const Ponto2 &A = P[k];
    const Ponto2 &B = P[(k + 1) % P.size()];
    const double va = val(A), vb = val(B);
    if (va >= 0.0) Q.push_back(A);
    if ((va > 0.0 && vb < 0.0) || (va < 0.0 && vb > 0.0)) {
      const double lam = va / (va - vb);
      Q.push_back({A[0] + lam * (B[0] - A[0]), A[1] + lam * (B[1] - A[1])});
    }
  }
  if (Q.size() < 3) return 0.0;

  double s = 0.0;
  for (std::size_t k = 0; k < Q.size(); k++) {
    const Ponto2 &A = Q[k];
    const Ponto2 &B = Q[(k + 1) % Q.size()];
    s += A[0] * B[1] - B[0] * A[1];
  }
  return std::fabs(s) / 2.0;
}

double plic_distancia_por_bisseccao(const Ponto2 &n, const Ponto2 &h,
                                    double area)
{
  const double nn = std::hypot(n[0], n[1]);
  if (!(nn > 0.0))
    throw std::invalid_argument(
        "plic: a normal e' nula -- nao ha' lado para a reta cortar");
  const double vol = h[0] * h[1];
  if (!(area >= 0.0) || area > vol)
    throw std::invalid_argument("plic_distancia_por_bisseccao: area fora de "
                                "[0, h0*h1]");

  // O SUPORTE EXATO: fora dele a area e' o volume inteiro ou zero, e bissecar
  // num intervalo maior gastaria iteracoes no patamar.
  const double S =
      (std::fabs(n[0]) / nn) * h[0] + (std::fabs(n[1]) / nn) * h[1];
  double lo = -S / 2.0, hi = S / 2.0;
  // A area e' DECRESCENTE em `d`: `lo` tem area `vol`, `hi` tem zero.
  for (int it = 0; it < 200; it++) {
    const double meio = 0.5 * (lo + hi);
    if (plic_area_por_recorte(n, h, meio) > area) lo = meio;
    else hi = meio;
  }
  return 0.5 * (lo + hi);
}

// ===========================================================================
// A FRACAO VOLUMETRICA INICIAL
// ===========================================================================

namespace {

//! A area assinalada de um poligono, NA FORMA DO HIGFLOW: soma de trapezios
//! `(y_i + y_{i+1})(x_{i+1} - x_i)/2`, aresta a aresta.
//
// Podia ser a forma do produto cruzado, que da' o mesmo numero a menos do
// ultimo bit.  E' esta porque o objetivo declarado e' o volume inicial
// COINCIDIR com o do HiGFlow, e usar a mesma aritmetica tira a diferenca de
// arredondamento da conversa -- custa nada e remove uma pergunta.
double area_trapezios(const std::vector<Ponto2> &P)
{
  double s = 0.0;
  for (std::size_t k = 0; k < P.size(); k++) {
    const Ponto2 &A = P[k];
    const Ponto2 &B = P[(k + 1) % P.size()];
    s += (A[1] + B[1]) * (B[0] - A[0]) / 2.0;
  }
  return std::fabs(s);
}

//! O ponto onde a forma cruza o segmento, por interpolacao LINEAR dos valores
//! nos extremos.  E' `Intersec` do HiGFlow.
Ponto2 cruzamento(const Ponto2 &a, const Ponto2 &b, double fa, double fb)
{
  const double lam = -fa / (fb - fa);
  return {a[0] + lam * (b[0] - a[0]), a[1] + lam * (b[1] - a[1])};
}

}  // namespace

double area_da_forma(const std::function<double(const Ponto2 &)> &forma,
                     const Ponto2 &centro, const Ponto2 &h)
{
  if (!(h[0] > 0.0) || !(h[1] > 0.0))
    throw std::invalid_argument("area_da_forma: lado de celula nao positivo");

  // OS CANTOS NA ORDEM DO HIGFLOW: 0 inferior-esquerdo, 1 inferior-direito,
  // 2 superior-esquerdo, 3 superior-DIREITO.  Note que NAO e' a ordem de
  // percurso do poligono -- 2 e 3 estao trocados em relacao a ela --, e e'
  // por isso que os casos abaixo listam os vertices um a um em vez de usar o
  // indice.
  const Ponto2 p0{centro[0] - h[0] / 2.0, centro[1] - h[1] / 2.0};
  const Ponto2 p1{centro[0] + h[0] / 2.0, centro[1] - h[1] / 2.0};
  const Ponto2 p2{centro[0] - h[0] / 2.0, centro[1] + h[1] / 2.0};
  const Ponto2 p3{centro[0] + h[0] / 2.0, centro[1] + h[1] / 2.0};
  const double f0 = forma(p0), f1 = forma(p1), f2 = forma(p2), f3 = forma(p3);
  const double vol = h[0] * h[1];

  int dentro = 0;
  if (f0 > 0.0) dentro++;
  if (f1 > 0.0) dentro++;
  if (f2 > 0.0) dentro++;
  if (f3 > 0.0) dentro++;

  if (dentro == 0) return 0.0;
  if (dentro == 4) return vol;

  if (dentro == 1 || dentro == 3) {
    // UM CANTO DE UM LADO, TRES DO OUTRO: o corte e' um triangulo nesse canto.
    // Com tres dentro, calcula-se o triangulo de FORA e subtrai.
    const bool um = (dentro == 1);
    auto tri = [&](const Ponto2 &c, const Ponto2 &a, double fc, double fa,
                   const Ponto2 &b, double fb) {
      return area_trapezios({c, cruzamento(c, a, fc, fa),
                             cruzamento(c, b, fc, fb)});
    };
    double t = 0.0;
    if (um ? (f0 > 0.0) : (f0 <= 0.0)) t = tri(p0, p1, f0, f1, p2, f2);
    else if (um ? (f1 > 0.0) : (f1 <= 0.0)) t = tri(p1, p0, f1, f0, p3, f3);
    else if (um ? (f2 > 0.0) : (f2 <= 0.0)) t = tri(p2, p0, f2, f0, p3, f3);
    else t = tri(p3, p1, f3, f1, p2, f2);
    return um ? t : vol - t;
  }

  // DOIS A DOIS.  So' os pares de ARESTA tem resposta; o par DIAGONAL nao.
  if (f0 > 0.0 && f1 > 0.0)
    return area_trapezios({p0, p1, cruzamento(p1, p3, f1, f3),
                           cruzamento(p0, p2, f0, f2)});
  if (f0 > 0.0 && f2 > 0.0)
    return area_trapezios({p2, p0, cruzamento(p0, p1, f0, f1),
                           cruzamento(p2, p3, f2, f3)});
  if (f1 > 0.0 && f3 > 0.0)
    return area_trapezios({p1, p3, cruzamento(p2, p3, f2, f3),
                           cruzamento(p0, p1, f0, f1)});
  if (f2 > 0.0 && f3 > 0.0)
    return area_trapezios({p3, p2, cruzamento(p0, p2, f0, f2),
                           cruzamento(p1, p3, f1, f3)});

  // O CASO AMBIGUO.  Dois cantos dentro na DIAGONAL: pela amostragem nos
  // cantos, a celula tanto pode conter duas quinas separadas quanto uma
  // garganta que as liga, e as duas tem areas diferentes.  Nao ha' resposta a
  // dar, e dar uma seria escolher.
  //
  // O HIGFLOW ABORTA AQUI com "malha pouco refinada para este problema", e a
  // leitura e' essa mesma: nao e' defeito do codigo, e' a malha nao resolvendo
  // a forma.  Lanca em vez de abortar porque quem chama pode estar numa
  // SUB-celula, e ai' a saida e' subdividir mais -- nao morrer.
  throw std::runtime_error(
      "area_da_forma: os dois cantos de dentro estao na DIAGONAL, e os quatro "
      "cantos nao distinguem duas quinas de uma garganta.  A malha nao resolve "
      "a forma neste ponto -- refine, ou aumente a subdivisao");
}

double fracao_inicial(const std::function<double(const Ponto2 &)> &forma,
                      const Ponto2 &centro, const Ponto2 &h, int n)
{
  if (n < 1) throw std::invalid_argument("fracao_inicial: subdivisao < 1");
  const double vol = h[0] * h[1];

  // O ATALHO DOS CANTOS, antes de subdividir: celula inteiramente dentro ou
  // inteiramente fora nao precisa de conta nenhuma, e e' a esmagadora maioria.
  const Ponto2 p0{centro[0] - h[0] / 2.0, centro[1] - h[1] / 2.0};
  const Ponto2 p1{centro[0] + h[0] / 2.0, centro[1] - h[1] / 2.0};
  const Ponto2 p2{centro[0] - h[0] / 2.0, centro[1] + h[1] / 2.0};
  const Ponto2 p3{centro[0] + h[0] / 2.0, centro[1] + h[1] / 2.0};
  int dentro = 0;
  if (forma(p0) > 0.0) dentro++;
  if (forma(p1) > 0.0) dentro++;
  if (forma(p2) > 0.0) dentro++;
  if (forma(p3) > 0.0) dentro++;
  if (dentro == 0) return 0.0;
  if (dentro == 4) return 1.0;
  if (n == 1) return area_da_forma(forma, centro, h) / vol;

  const Ponto2 hs{h[0] / n, h[1] / n};
  double soma = 0.0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      const Ponto2 c{p0[0] + (i + 0.5) * hs[0], p0[1] + (j + 0.5) * hs[1]};
      soma += area_da_forma(forma, c, hs);
    }
  return soma / vol;
}

void escreve_fracao_inicial(const Malha &m,
                            const std::function<double(const Ponto2 &)> &forma,
                            CampoEscalar *frac, int n)
{
  if (frac == nullptr)
    throw std::invalid_argument("escreve_fracao_inicial: campo nulo");
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      frac->em(m.id_global(b, i)) =
          fracao_inicial(forma, m.centroide(b, i), m.lados_da_folha(b, i), n);
}

}  // namespace t8flow2d
