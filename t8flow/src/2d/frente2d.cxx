#include "frente2d.hxx"

#include <cmath>
#include <stdexcept>

namespace t8flow2d {

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double EPS = 1e-300;
}  // namespace

Frente::Frente(std::vector<Ponto2> marcadores) : m_(std::move(marcadores))
{
  // TRES E' O MINIMO, e nao dois: a curvatura e' do circunraio de tres
  // consecutivos, e com dois marcadores `k-1` e `k+1` sao o mesmo ponto.
  if (m_.size() < 3)
    throw std::invalid_argument(
        "Frente: " + std::to_string(m_.size()) +
        " marcadores -- o minimo e' 3, porque a curvatura sai do circunraio de "
        "tres consecutivos e com menos eles se repetem");
}

Frente Frente::circulo(const Ponto2 &centro, double raio, int n)
{
  if (!(raio > 0.0))
    throw std::invalid_argument("Frente::circulo: raio tem de ser positivo");
  if (n < 3)
    throw std::invalid_argument("Frente::circulo: ao menos 3 marcadores");
  std::vector<Ponto2> m;
  m.reserve(static_cast<size_t>(n));
  for (int k = 0; k < n; k++) {
    const double a = 2.0 * PI * k / n;     // ANTI-HORARIO: area positiva
    m.push_back({centro[0] + raio * std::cos(a), centro[1] + raio * std::sin(a)});
  }
  return Frente(std::move(m));
}

Frente Frente::elipse(const Ponto2 &centro, double a, double b, int n)
{
  if (!(a > 0.0) || !(b > 0.0))
    throw std::invalid_argument("Frente::elipse: semieixos tem de ser positivos");
  if (n < 3)
    throw std::invalid_argument("Frente::elipse: ao menos 3 marcadores");
  std::vector<Ponto2> m;
  m.reserve(static_cast<size_t>(n));
  for (int k = 0; k < n; k++) {
    const double t = 2.0 * PI * k / n;
    m.push_back({centro[0] + a * std::cos(t), centro[1] + b * std::sin(t)});
  }
  return Frente(std::move(m));
}

const Ponto2 &Frente::marcador(long k) const
{
  const long n = static_cast<long>(m_.size());
  // INDICE CIRCULAR: a frente e' um anel, e `k-1` em `k = 0` e' o ultimo.
  // Fazer isso na chamada espalharia a aritmetica modular por todo o arquivo,
  // e e' exatamente onde o erro de um a mais mora.
  long q = k % n;
  if (q < 0) q += n;
  return m_[static_cast<size_t>(q)];
}

double Frente::area() const
{
  // FORMULA DO CADARCO, exata para o poligono.
  const long n = static_cast<long>(m_.size());
  double s = 0.0;
  for (long k = 0; k < n; k++) {
    const Ponto2 &a = m_[static_cast<size_t>(k)];
    const Ponto2 &b = marcador(k + 1);
    s += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * s;
}

double Frente::perimetro() const
{
  const long n = static_cast<long>(m_.size());
  double s = 0.0;
  for (long k = 0; k < n; k++) {
    const Ponto2 &a = m_[static_cast<size_t>(k)];
    const Ponto2 &b = marcador(k + 1);
    s += std::hypot(b[0] - a[0], b[1] - a[1]);
  }
  return s;
}

Ponto2 Frente::centroide() const
{
  // O CENTROIDE DA REGIAO, e nao a media dos marcadores.  A media so' coincide
  // com ele quando os marcadores sao equidistantes, e depois da cirurgia eles
  // nao sao -- usar a media faria a posicao da gota depender de onde a
  // cirurgia aconteceu.
  const long n = static_cast<long>(m_.size());
  double cx = 0.0, cy = 0.0, a2 = 0.0;
  for (long k = 0; k < n; k++) {
    const Ponto2 &a = m_[static_cast<size_t>(k)];
    const Ponto2 &b = marcador(k + 1);
    const double cruz = a[0] * b[1] - b[0] * a[1];
    a2 += cruz;
    cx += (a[0] + b[0]) * cruz;
    cy += (a[1] + b[1]) * cruz;
  }
  if (std::fabs(a2) < EPS) return m_[0];
  return {cx / (3.0 * a2), cy / (3.0 * a2)};
}

double Frente::curvatura(long k) const
{
  // O CIRCUNRAIO de `k-1, k, k+1`:  `kappa = 4 A / (a b c)`.
  //
  // Para tres pontos SOBRE UM CIRCULO o circunraio E' o raio, qualquer que
  // seja o espacamento -- logo isto e' EXATO para um circulo, e nao aproximado
  // de segunda ordem.  A diferenca finita parametrica daria so' ordem.
  const Ponto2 &p = marcador(k - 1);
  const Ponto2 &q = marcador(k);
  const Ponto2 &r = marcador(k + 1);

  const double a = std::hypot(q[0] - p[0], q[1] - p[1]);
  const double b = std::hypot(r[0] - q[0], r[1] - q[1]);
  const double c = std::hypot(r[0] - p[0], r[1] - p[1]);
  if (a < EPS || b < EPS || c < EPS) return 0.0;

  // O DOBRO DA AREA COM SINAL: positivo no anti-horario, que e' o sentido em
  // que a frente e' montada.  E' o sinal que distingue uma saliencia de uma
  // reentrancia -- tomar o modulo aqui faria a tensao superficial empurrar
  // para o mesmo lado nos dois casos.
  const double cruz = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
  return 2.0 * cruz / (a * b * c);
}

Ponto2 Frente::normal(long k) const
{
  // A TANGENTE pela diferenca central entre os vizinhos, e a normal EXTERIOR e'
  // ela girada de -90 graus -- com a frente no sentido anti-horario, `(ty, -tx)`
  // aponta para fora.
  const Ponto2 &p = marcador(k - 1);
  const Ponto2 &r = marcador(k + 1);
  const double tx = r[0] - p[0], ty = r[1] - p[1];
  const double L = std::hypot(tx, ty);
  if (L < EPS) return {0.0, 0.0};
  return {ty / L, -tx / L};
}

void Frente::avanca(const CampoDeAdvecao &u, double t, double dt, bool rk2)
{
  if (!u)
    throw std::invalid_argument("Frente::avanca: campo de adveccao nao instalado");

  std::vector<Ponto2> novo(m_.size());
  for (size_t k = 0; k < m_.size(); k++) {
    const Ponto2 &x = m_[k];
    const Ponto2 v1 = u(x, t);
    if (!rk2) {
      novo[k] = {x[0] + dt * v1[0], x[1] + dt * v1[1]};
      continue;
    }
    // RK2 DO PONTO MEDIO.  Com Euler, um campo de rotacao pura afasta o
    // marcador do circulo em `O(dt)` por passo e a gota incha -- deriva
    // radial que se acumula.  Com RK2 e' `O(dt^2)` e some no refino.
    const Ponto2 meio{x[0] + 0.5 * dt * v1[0], x[1] + 0.5 * dt * v1[1]};
    const Ponto2 v2 = u(meio, t + 0.5 * dt);
    novo[k] = {x[0] + dt * v2[0], x[1] + dt * v2[1]};
  }
  m_.swap(novo);
}

void Frente::move(const std::vector<Ponto2> &vel, double dt)
{
  if (vel.size() != m_.size())
    throw std::invalid_argument(
        "Frente::move: " + std::to_string(vel.size()) + " velocidades para " +
        std::to_string(m_.size()) + " marcadores.  A cirurgia muda a contagem, "
        "e interpolar antes dela deixa as duas listas em desacordo -- o que "
        "nao falharia se os tamanhos fossem so' parecidos");
  for (size_t k = 0; k < m_.size(); k++) {
    m_[k][0] += dt * vel[k][0];
    m_[k][1] += dt * vel[k][1];
  }
}

Frente::Cirurgia Frente::cirurgia(double l_min, double l_max)
{
  if (!(l_max > 0.0) || !(l_min > 0.0))
    throw std::invalid_argument("Frente::cirurgia: limites tem de ser positivos");
  // O LIMITE INFERIOR TEM DE SER MENOS DA METADE DO SUPERIOR.
  //
  // Com `l_min >= l_max/2`, partir uma aresta de comprimento pouco acima de
  // `l_max` produz duas de pouco mais de `l_max/2`, que ja' nascem abaixo do
  // minimo -- e a cirurgia seguinte as funde de volta.  A frente oscila entre
  // dois estados a cada passo, sem convergir, e nada acusa: os contadores
  // mostram atividade e a area mal se move.
  if (!(l_min < l_max / 2.0))
    throw std::invalid_argument(
        "Frente::cirurgia: l_min = " + std::to_string(l_min) +
        " nao e' menor que l_max/2 = " + std::to_string(l_max / 2.0) +
        " -- partir uma aresta produziria duas abaixo do minimo, e a cirurgia "
        "seguinte as fundiria de volta: a frente oscilaria entre dois estados "
        "a cada passo sem convergir");

  Cirurgia c;

  // PARTIR primeiro, e num percurso so'.  Partir e fundir na mesma varredura
  // faria o resultado depender da ordem de visita.
  {
    std::vector<Ponto2> novo;
    novo.reserve(m_.size() * 2);
    const long n = static_cast<long>(m_.size());
    for (long k = 0; k < n; k++) {
      const Ponto2 &a = m_[static_cast<size_t>(k)];
      const Ponto2 &b = marcador(k + 1);
      novo.push_back(a);
      if (std::hypot(b[0] - a[0], b[1] - a[1]) > l_max) {
        novo.push_back({0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])});
        c.partidas++;
      }
    }
    m_.swap(novo);
  }

  // FUNDIR: remove o marcador cuja aresta seguinte ficou curta, substituindo
  // os dois pelo ponto medio.  Nunca abaixo de tres -- uma frente de dois
  // marcadores nao tem curvatura.
  {
    std::vector<Ponto2> novo;
    novo.reserve(m_.size());
    const long n = static_cast<long>(m_.size());
    long k = 0;
    while (k < n) {
      const Ponto2 &a = m_[static_cast<size_t>(k)];
      const Ponto2 &b = marcador(k + 1);
      const bool curta = std::hypot(b[0] - a[0], b[1] - a[1]) < l_min;
      const bool cabe = (n - c.fundidas) > 3;
      if (curta && cabe && k + 1 < n) {
        novo.push_back({0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])});
        c.fundidas++;
        k += 2;                       // os dois viraram um
      } else {
        novo.push_back(a);
        k++;
      }
    }
    m_.swap(novo);
  }
  return c;
}

namespace {

//! O comprimento da ESTRELA do marcador: metade da soma dos dois segmentos
//! incidentes.  E' o analogo plano do volume da estrela que a versao
//! tridimensional usa como peso, e e' o que faz a altura comum NAO ser a media
//! simples quando os dois lados tem resolucao diferente.
double estrela_de(const Frente &f, long i)
{
  const Ponto2 &a = f.marcador(i - 1);
  const Ponto2 &p = f.marcador(i);
  const Ponto2 &b = f.marcador(i + 1);
  return 0.5 * (std::hypot(p[0] - a[0], p[1] - a[1]) +
                std::hypot(b[0] - p[0], b[1] - p[1]));
}

}  // namespace

long Frente::tsur(double rlx, double limiar)
{
  const long n = static_cast<long>(m_.size());
  if (n < 4) return 0;

  // OS DESLOCAMENTOS SAO ACUMULADOS E APLICADOS DE UMA VEZ.  Aplicar em
  // seguida faria o segundo segmento ler posicoes ja' movidas pelo primeiro, e
  // o resultado passaria a depender da ORDEM DE PERCURSO -- que nao e'
  // propriedade da frente.
  std::vector<Ponto2> desl(static_cast<size_t>(n), Ponto2{0.0, 0.0});
  std::vector<char> tocado(static_cast<size_t>(n), 0);

  for (long e = 0; e < n; e++) {
    const long v1 = e, v2 = (e + 1) % n;
    const Ponto2 n1 = normal(v1), n2 = normal(v2);
    Ponto2 ne{n1[0] + n2[0], n1[1] + n2[1]};
    const double mn = std::hypot(ne[0], ne[1]);
    if (mn < 1e-7) continue;          // normais opostas: nada a definir
    ne[0] /= mn;
    ne[1] /= mn;

    const Ponto2 &p1 = m_[static_cast<size_t>(v1)];
    const Ponto2 &p2 = m_[static_cast<size_t>(v2)];
    const Ponto2 meio{0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1])};
    const double h1 = (p1[0] - meio[0]) * ne[0] + (p1[1] - meio[1]) * ne[1];
    const double h2 = (p2[0] - meio[0]) * ne[0] + (p2[1] - meio[1]) * ne[1];

    // O LIMIAR compara `|h1|` com a MEIA-ARESTA, que e' o maximo que ele pode
    // valer -- de modo que `limiar` e' uma fracao adimensional e nao depende
    // da resolucao da frente.
    const double meia = 0.5 * std::hypot(p2[0] - p1[0], p2[1] - p1[1]);
    if (meia <= 0.0 || std::fabs(h1) < limiar * meia) continue;

    const double w1 = estrela_de(*this, v1), w2 = estrela_de(*this, v2);
    if (w1 + w2 <= 0.0) continue;
    const double hc = (w1 * h1 + w2 * h2) / (w1 + w2);

    const double g = 1.0 - rlx;
    desl[static_cast<size_t>(v1)][0] += g * (hc - h1) * ne[0];
    desl[static_cast<size_t>(v1)][1] += g * (hc - h1) * ne[1];
    desl[static_cast<size_t>(v2)][0] += g * (hc - h2) * ne[0];
    desl[static_cast<size_t>(v2)][1] += g * (hc - h2) * ne[1];
    tocado[static_cast<size_t>(v1)] = 1;
    tocado[static_cast<size_t>(v2)] = 1;
  }

  long movidos = 0;
  for (long i = 0; i < n; i++) {
    if (!tocado[static_cast<size_t>(i)]) continue;
    m_[static_cast<size_t>(i)][0] += desl[static_cast<size_t>(i)][0];
    m_[static_cast<size_t>(i)][1] += desl[static_cast<size_t>(i)][1];
    movidos++;
  }
  return movidos;
}

}  // namespace t8flow2d
