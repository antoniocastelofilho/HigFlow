#include "criterio2d.hxx"

#include <algorithm>

#include <cmath>
#include <stdexcept>

namespace t8flow2d {

std::vector<double> amr_limiares(double h_base, int niveis, int cel_min,
                                 double largura_fina)
{
  if (!(h_base > 0.0))
    throw std::invalid_argument("amr_limiares: h_base tem de ser positivo");
  if (niveis < 1)
    throw std::invalid_argument("amr_limiares: niveis tem de ser ao menos 1");
  if (cel_min < 1)
    throw std::invalid_argument(
        "amr_limiares: cel_min " + std::to_string(cel_min) +
        " -- uma banda de zero celulas finas nao refina coisa nenhuma, e "
        "aceitar isso calado deixaria a malha PARECER refinada");

  double h_fino = h_base;
  for (int l = 0; l < niveis; l++) h_fino *= 0.5;

  // UM `h_base` INTEIRO DE MARGEM, e a razao e' geometrica: o criterio decide
  // por celula base, logo a borda da banda e' incerta em `h_base/2` pelo teste
  // no centro, E a celula se estende outro `h_base/2` alem dele.
  const double margem = h_base + h_fino;

  // A BANDA FINA TEM DUAS EXIGENCIAS, E SO' UMA DELAS E' DA MALHA.
  //
  //   `cel_min * h_fino`  e' da DISCRETIZACAO: o esquema precisa de tantas
  //                       celulas finas de cada lado da semente, e isso escala
  //                       com `h` porque celula e' que escala com `h`;
  //
  //   `largura_fina`      e' do PROBLEMA: a camada em volta da interface tem
  //                       espessura FISICA, que nao diminui porque a malha
  //                       refinou.
  //
  // O MAIOR DOS DOIS, e nao a soma nem so' o primeiro.
  //
  // SO' O PRIMEIRO -- que e' o que estava aqui -- FAZ O ESTUDO ANTICONVERGIR, e
  // isso foi medido na bolha de Hysing: a banda valia 0,50 R em `h_min =
  // 0,0125` e 0,25 R em 0,00625, de modo que refinar ENCOLHEU PELA METADE a
  // regiao resolvida em volta da interface.  `v_max` foi de -6,5% para -10,5%
  // ao refinar, contra +0,4% da malha uniforme na mesma `h_min` -- o AMR
  // piorava com refino, e nao por razao de discretizacao nenhuma.
  //
  // E ERA RESULTADO CONHECIDO: a secao companheira ja' havia registrado que
  // banda em numero fixo de celulas refina largura fisica diferente a cada
  // resolucao, e que o estudo nao converge por razoes alheias a'
  // discretizacao.  Estava anotado, e este codigo fazia o contrario -- ter
  // anotado nao substituiu conferir o codigo contra a anotacao.
  //
  // O PADRAO E' ZERO, e nao um valor conveniente: com zero isto reduz
  // EXATAMENTE ao que era, e os casos que nao passam largura -- a caixa
  // manufaturada, o arreio -- continuam bit a bit.  Quem tem escala fisica a
  // declarar e' quem conhece o problema, e nao este arquivo: o criterio nao
  // sabe de que forma vem a semente.
  if (largura_fina < 0.0)
    throw std::invalid_argument(
        "amr_limiares: largura_fina negativa nao quer dizer nada");
  const double banda_fina = std::max(largura_fina, cel_min * h_fino) + margem;

  std::vector<double> thr(static_cast<size_t>(niveis), 0.0);
  for (int l = niveis - 1; l >= 0; l--) {
    thr[static_cast<size_t>(l)] = banda_fina;
    double h_nivel = h_base;
    for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
    for (int k = 0; k < niveis - 1 - l; k++)
      thr[static_cast<size_t>(l)] += cel_min * h_nivel;
  }
  return thr;
}

namespace {

//! O nivel que uma distancia pede: o maior `l+1` cujo limiar ainda a cobre.
//
// MONOTONO DECRESCENTE na distancia, e e' dessa monotonia que a dilatacao tira
// a exatidao -- o maximo sobre as sementes e' o nivel da semente mais proxima.
int nivel_da_distancia(double dist, const std::vector<double> &thr)
{
  int nv = 0;
  for (size_t l = 0; l < thr.size(); l++)
    if (dist <= thr[l]) nv = static_cast<int>(l) + 1;
  return nv;
}

}  // namespace

std::vector<signed char> amr_tabela(const std::vector<Ponto2> &sem,
                                    const Ponto2 &lo, double hx, double hy,
                                    int nx, int ny,
                                    const std::vector<double> &thr)
{
  std::vector<signed char> tab(static_cast<size_t>(nx) * ny, 0);
  if (sem.empty() || thr.empty()) return tab;

  // O RAIO DO CARIMBO vem do MAIOR limiar, que e' `thr[0]`.  O `+1` cobre o
  // arredondamento do indice: sem ele a ultima coluna da banda fica de fora, e
  // o que se perde e' exatamente a borda -- onde a banda deveria terminar.
  const int Rx = static_cast<int>(thr[0] / hx) + 1;
  const int Ry = static_cast<int>(thr[0] / hy) + 1;

  for (const Ponto2 &p : sem) {
    const int ic = static_cast<int>((p[0] - lo[0]) / hx);
    const int jc = static_cast<int>((p[1] - lo[1]) / hy);
    for (int j = jc - Ry; j <= jc + Ry; j++) {
      if (j < 0 || j >= ny) continue;
      for (int i = ic - Rx; i <= ic + Rx; i++) {
        if (i < 0 || i >= nx) continue;
        const double cx = lo[0] + (i + 0.5) * hx;
        const double cy = lo[1] + (j + 0.5) * hy;
        const double dist = std::hypot(cx - p[0], cy - p[1]);
        const int nv = nivel_da_distancia(dist, thr);
        const size_t q = static_cast<size_t>(i) + static_cast<size_t>(nx) * j;
        if (nv > tab[q]) tab[q] = static_cast<signed char>(nv);
      }
    }
  }
  return tab;
}

std::vector<signed char> amr_tabela_exaustiva(const std::vector<Ponto2> &sem,
                                              const Ponto2 &lo, double hx,
                                              double hy, int nx, int ny,
                                              const std::vector<double> &thr)
{
  std::vector<signed char> tab(static_cast<size_t>(nx) * ny, 0);
  if (sem.empty() || thr.empty()) return tab;

  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++) {
      const double cx = lo[0] + (i + 0.5) * hx;
      const double cy = lo[1] + (j + 0.5) * hy;
      double dmin = 1e300;
      for (const Ponto2 &p : sem)
        dmin = std::fmin(dmin, std::hypot(cx - p[0], cy - p[1]));
      tab[static_cast<size_t>(i) + static_cast<size_t>(nx) * j] =
          static_cast<signed char>(nivel_da_distancia(dmin, thr));
    }
  return tab;
}

int amr_nivel_em(const std::vector<signed char> &tab, const Ponto2 &lo,
                 double hx, double hy, int nx, int ny, const Ponto2 &x)
{
  const int i = static_cast<int>((x[0] - lo[0]) / hx);
  const int j = static_cast<int>((x[1] - lo[1]) / hy);
  if (i < 0 || i >= nx || j < 0 || j >= ny) return 0;
  return tab[static_cast<size_t>(i) + static_cast<size_t>(nx) * j];
}

}  // namespace t8flow2d
