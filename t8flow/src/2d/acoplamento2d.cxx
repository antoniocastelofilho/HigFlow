#include "acoplamento2d.hxx"

#include "mls2d.hxx"

#include <cmath>
#include <set>
#include <unordered_map>
#include <utility>
#include <stdexcept>

namespace t8flow2d {

namespace {

constexpr double PI = 3.14159265358979323846;

//! Recorta um poligono por UM semiplano, na forma de Sutherland-Hodgman.
//! `lado(p) >= 0` e' o lado que fica.
template <typename F>
std::vector<Ponto2> recorta_semiplano(const std::vector<Ponto2> &p, F lado)
{
  std::vector<Ponto2> out;
  if (p.empty()) return out;
  out.reserve(p.size() + 4);
  for (size_t k = 0; k < p.size(); k++) {
    const Ponto2 &a = p[k];
    const Ponto2 &b = p[(k + 1) % p.size()];
    const double la = lado(a), lb = lado(b);
    if (la >= 0.0) out.push_back(a);
    if ((la >= 0.0) != (lb >= 0.0)) {
      // A INTERSECAO COM A RETA, pelo parametro.  O denominador nao se anula
      // porque os sinais diferem.
      const double t = la / (la - lb);
      out.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])});
    }
  }
  return out;
}

double cadarco(const std::vector<Ponto2> &p)
{
  if (p.size() < 3) return 0.0;
  double s = 0.0;
  for (size_t k = 0; k < p.size(); k++) {
    const Ponto2 &a = p[k];
    const Ponto2 &b = p[(k + 1) % p.size()];
    s += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * s;
}

//! O nucleo de Peskin de quatro pontos, em uma dimensao.  `r` em unidades de h.
double peskin(double r)
{
  const double a = std::fabs(r);
  if (a >= 2.0) return 0.0;
  if (a >= 1.0)
    return (5.0 - 2.0 * a - std::sqrt(-7.0 + 12.0 * a - 4.0 * a * a)) / 8.0;
  return (3.0 - 2.0 * a + std::sqrt(1.0 + 4.0 * a - 4.0 * a * a)) / 8.0;
}

}  // namespace

double area_recortada(const std::vector<Ponto2> &poligono, const Ponto2 &lo,
                      const Ponto2 &hi)
{
  std::vector<Ponto2> p = poligono;
  p = recorta_semiplano(p, [&](const Ponto2 &q) { return q[0] - lo[0]; });
  p = recorta_semiplano(p, [&](const Ponto2 &q) { return hi[0] - q[0]; });
  p = recorta_semiplano(p, [&](const Ponto2 &q) { return q[1] - lo[1]; });
  p = recorta_semiplano(p, [&](const Ponto2 &q) { return hi[1] - q[1]; });
  return std::fabs(cadarco(p));
}

void escreve_fracao(const Frente &f, const Malha &m, CampoEscalar *c)
{
  if (c == nullptr) throw std::invalid_argument("escreve_fracao: campo nulo");
  const std::vector<Ponto2> &p = f.marcadores();

  // A ENVOLTORIA DA FRENTE limita o trabalho: fora dela a fracao e' zero, e
  // recortar o poligono contra cada celula do dominio desperdicaria o grosso
  // das chamadas.  A folga de uma celula cobre o arredondamento.
  Ponto2 elo{1e300, 1e300}, ehi{-1e300, -1e300};
  for (const Ponto2 &q : p)
    for (int d = 0; d < 2; d++) {
      if (q[d] < elo[d]) elo[d] = q[d];
      if (q[d] > ehi[d]) ehi[d] = q[d];
    }

  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const Ponto2 cc = m.centroide(b, i);
      const Ponto2 hh = m.lados_da_folha(b, i);
      const Ponto2 lo{cc[0] - hh[0] / 2, cc[1] - hh[1] / 2};
      const Ponto2 hi{cc[0] + hh[0] / 2, cc[1] + hh[1] / 2};
      const long id = m.id_global(b, i);

      if (hi[0] < elo[0] || lo[0] > ehi[0] || hi[1] < elo[1] || lo[1] > ehi[1]) {
        c->em(id) = 0.0;
        continue;
      }
      c->em(id) = area_recortada(p, lo, hi) / (hh[0] * hh[1]);
    }
}

void interpola_propriedade(const Malha &m, const CampoEscalar &c,
                           double dentro, double fora, CampoEscalar *prop)
{
  if (prop == nullptr)
    throw std::invalid_argument("interpola_propriedade: campo nulo");
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const long id = m.id_global(b, i);
      prop->em(id) = fora + (dentro - fora) * c.em(id);
    }
}

void forca_tensao(const Frente &f, const Malha &m, const DominioFacetas &fd,
                  double sigma, CampoVetorial *forca)
{
  if (forca == nullptr) throw std::invalid_argument("forca_tensao: campo nulo");
  forca->preenche(0.0);

  // `(bloco, folha)` POR ID, uma vez.  Procurar o par de cada celula por
  // varredura dentro do laco dos marcadores daria quadratico -- o mesmo que
  // ja' custou caro na montagem do Poisson.
  std::unordered_map<long, std::pair<long, long>> bi;
  bi.reserve(static_cast<size_t>(m.num_folhas()) * 2);
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi[m.id_global(b, i)] = {b, i};

  const long nm = f.num_marcadores();
  for (long k = 0; k < nm; k++) {
    const Ponto2 &X = f.marcador(k);

    // O `h` DA CELULA QUE CONTEM O MARCADOR.  Numa malha adaptada nao ha' um
    // `h` do dominio, e distribuir com o errado espalha a forca por largura
    // errada -- muda a curvatura efetiva e desloca o salto de pressao sem que
    // nada acuse.
    const long cel = m.id_do_ponto(X);
    const auto itc = (cel >= 0) ? bi.find(cel) : bi.end();
    if (itc == bi.end()) continue;              // marcador fora do dominio
    const Ponto2 hh = m.lados_da_folha(itc->second.first, itc->second.second);

    // O COMPRIMENTO DE ARCO do marcador: metade de cada aresta vizinha.
    const Ponto2 &A = f.marcador(k - 1);
    const Ponto2 &B = f.marcador(k + 1);
    const double ds = 0.5 * (std::hypot(X[0] - A[0], X[1] - A[1]) +
                             std::hypot(B[0] - X[0], B[1] - X[1]));

    const double kap = f.curvatura(k);
    const Ponto2 nrm = f.normal(k);
    // A FORCA APONTA PARA DENTRO numa gota convexa: `-sigma kappa n`, com `n`
    // exterior.  O sinal trocado faria a gota explodir em vez de se conter --
    // e isso nao falha, a corrida roda.
    const double comp[2] = {-sigma * kap * nrm[0] * ds,
                            -sigma * kap * nrm[1] * ds};

    // AS FACETAS DO SUPORTE, pelas CELULAS da vizinhanca e nao por varredura
    // global.  O nucleo de Peskin tem suporte de dois `h` em cada direcao,
    // logo bastam as celulas de uma caixa de raio 2.
    std::set<std::pair<int, long>> vistas;
    for (int jj = -2; jj <= 2; jj++)
      for (int ii = -2; ii <= 2; ii++) {
        const Ponto2 q{X[0] + ii * hh[0], X[1] + jj * hh[1]};
        const long id = m.id_do_ponto(q);
        const auto it = (id >= 0) ? bi.find(id) : bi.end();
        if (it == bi.end()) continue;
        for (int d = 0; d < 2; d++)
          for (int lado = 0; lado < 2; lado++) {
            const long fid = fd.id_faceta(it->second.first, it->second.second,
                                          d, lado);
            if (fid >= 0) vistas.insert({d, fid});
          }
      }

    for (const auto &par : vistas) {
      const int d = par.first;
      const long q = par.second;
      if (comp[d] == 0.0) continue;
      const Ponto2 xf = fd.centro(d, q);
      const double rx = (xf[0] - X[0]) / hh[0];
      const double ry = (xf[1] - X[1]) / hh[1];
      if (std::fabs(rx) >= 2.0 || std::fabs(ry) >= 2.0) continue;
      const double w = peskin(rx) * peskin(ry) / (hh[0] * hh[1]);
      forca->em(d, q) += comp[d] * w;
    }
  }
}

std::vector<Ponto2> velocidade_nos_marcadores(const Frente &f, const Malha &m,
                                              const DominioFacetas &fd,
                                              const CampoVetorial &u)
{
  std::unordered_map<long, std::pair<long, long>> bi;
  bi.reserve(static_cast<size_t>(m.num_folhas()) * 2);
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi[m.id_global(b, i)] = {b, i};

  const long nm = f.num_marcadores();
  std::vector<Ponto2> U(static_cast<size_t>(nm), Ponto2{0.0, 0.0});

  for (long k = 0; k < nm; k++) {
    const Ponto2 &X = f.marcador(k);
    const long cel = m.id_do_ponto(X);
    const auto itc = (cel >= 0) ? bi.find(cel) : bi.end();
    if (itc == bi.end()) continue;
    const Ponto2 hh = m.lados_da_folha(itc->second.first, itc->second.second);

    // AS MESMAS FACETAS que a distribuicao alcanca, pelas mesmas celulas.  Se
    // os dois lacos nao percorressem o mesmo conjunto, o par deixaria de ser
    // adjunto -- e e' do par que vem a conservacao.
    std::set<std::pair<int, long>> vistas;
    for (int jj = -2; jj <= 2; jj++)
      for (int ii = -2; ii <= 2; ii++) {
        const Ponto2 q{X[0] + ii * hh[0], X[1] + jj * hh[1]};
        const long id = m.id_do_ponto(q);
        const auto it = (id >= 0) ? bi.find(id) : bi.end();
        if (it == bi.end()) continue;
        for (int d = 0; d < 2; d++)
          for (int lado = 0; lado < 2; lado++) {
            const long fid = fd.id_faceta(it->second.first, it->second.second,
                                          d, lado);
            if (fid >= 0) vistas.insert({d, fid});
          }
      }

    for (const auto &par : vistas) {
      const int d = par.first;
      const long q = par.second;
      const Ponto2 xf = fd.centro(d, q);
      const double rx = (xf[0] - X[0]) / hh[0];
      const double ry = (xf[1] - X[1]) / hh[1];
      if (std::fabs(rx) >= 2.0 || std::fabs(ry) >= 2.0) continue;
      // O MESMO NUCLEO, com `h^2` no lugar de `1/h^2`: e' o que torna a
      // interpolacao adjunta da distribuicao, e nao so' parecida com ela.
      U[static_cast<size_t>(k)][d] += u.em(d, q) * peskin(rx) * peskin(ry);
    }
  }
  return U;
}

std::vector<Ponto2> velocidade_nos_marcadores_mls(const Frente &f,
                                                  const Malha &m,
                                                  const DominioFacetas &fd,
                                                  const CampoVetorial &u,
                                                  int ordem)
{
  std::unordered_map<long, std::pair<long, long>> bi;
  bi.reserve(static_cast<size_t>(m.num_folhas()) * 2);
  for (long b = 0; b < m.num_blocos(); b++)
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++)
      bi[m.id_global(b, i)] = {b, i};

  const long nm = f.num_marcadores();
  std::vector<Ponto2> U(static_cast<size_t>(nm), Ponto2{0.0, 0.0});
  Mls mls(ordem, 64);

  for (long k = 0; k < nm; k++) {
    const Ponto2 &X = f.marcador(k);
    const long cel = m.id_do_ponto(X);
    const auto itc = (cel >= 0) ? bi.find(cel) : bi.end();
    if (itc == bi.end()) continue;
    const Ponto2 hh = m.lados_da_folha(itc->second.first, itc->second.second);

    for (int d = 0; d < 2; d++) {
      // O SUPORTE SAO AS FACETAS DA DIRECAO `d` na vizinhanca.  O raio cresce
      // ate' juntar o minimo que a ordem exige -- o mesmo criterio do
      // estencil de celulas, e pela mesma razao: numero de pontos suficiente
      // nao e' posicionamento suficiente, mas numero INsuficiente e' recusa
      // certa.
      std::set<std::pair<int, long>> vistas;
      std::vector<Ponto2> sup;
      std::vector<double> val;
      for (int raio = 2; raio <= 4 &&
           static_cast<int>(sup.size()) < Mls::min_pontos(ordem); raio++) {
        sup.clear();
        val.clear();
        vistas.clear();
        for (int jj = -raio; jj <= raio; jj++)
          for (int ii = -raio; ii <= raio; ii++) {
            const Ponto2 q{X[0] + ii * hh[0], X[1] + jj * hh[1]};
            const long id = m.id_do_ponto(q);
            const auto it = (id >= 0) ? bi.find(id) : bi.end();
            if (it == bi.end()) continue;
            for (int lado = 0; lado < 2; lado++) {
              const long fid = fd.id_faceta(it->second.first,
                                            it->second.second, d, lado);
              if (fid >= 0 && vistas.insert({d, fid}).second) {
                sup.push_back(fd.centro(d, fid));
                val.push_back(u.em(d, fid));
              }
            }
          }
      }
      if (static_cast<int>(sup.size()) < Mls::min_pontos(ordem)) continue;
      U[static_cast<size_t>(k)][d] = mls.interpola(X, sup, val);
    }
  }
  return U;
}

Ponto2 forca_total(const DominioFacetas &fd, const CampoVetorial &forca)
{
  Ponto2 s{0.0, 0.0};
  for (int d = 0; d < 2; d++)
    for (long q = 0; q < fd.num_facetas(d); q++) {
      // O VOLUME DE CONTROLE da faceta e' `medida x h_normal`; com a medida
      // transversal so', a soma nao teria dimensao de forca.
      const Ponto2 c = fd.centro(d, q);
      (void) c;
      s[d] += forca.em(d, q) * fd.medida(d, q);
    }
  return s;
}

}  // namespace t8flow2d
