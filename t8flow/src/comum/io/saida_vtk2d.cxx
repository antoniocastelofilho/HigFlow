#include "saida_vtk2d.hxx"

#include <cstdio>
#include <stdexcept>

namespace t8flow {

using t8flow2d::CampoEscalar;
using t8flow2d::CampoTensor;
using t8flow2d::CampoVetorial;
using t8flow2d::DominioFacetas;
using t8flow2d::Malha;
using t8flow2d::Ponto2;

SaidaVtk2d::SaidaVtk2d(const Malha &m) : m_(&m) {}

void SaidaVtk2d::limpa()
{
  campos_.clear();
  fronteiras_.clear();
}

void SaidaVtk2d::escalar(const std::string &nome, const CampoEscalar &c)
{
  Campo f;
  f.nome = nome;
  f.componentes = 1;
  f.v.reserve(static_cast<size_t>(m_->num_folhas()));
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++)
      f.v.push_back(c.em(m_->id_global(b, i)));
  campos_.push_back(std::move(f));
}

void SaidaVtk2d::vetorial(const std::string &nome, const DominioFacetas &fd,
                          const CampoVetorial &v)
{
  // A MEDIA DAS DUAS FACETAS OPOSTAS, por direcao.  E' a interpolacao que o
  // cabecalho avisa: ela e' de segunda ordem no interior e de PRIMEIRA onde
  // uma das duas nao existe, e e' por isso que o arquivo nao serve de oraculo.
  Campo f;
  f.nome = nome;
  f.componentes = 3;                 // o VTK quer tres, mesmo em 2D
  f.v.reserve(static_cast<size_t>(m_->num_folhas()) * 3);
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      double comp[3] = {0.0, 0.0, 0.0};
      for (int d = 0; d < 2; d++) {
        const long fm = fd.id_faceta(b, i, d, 0);
        const long fp = fd.id_faceta(b, i, d, 1);
        // A faceta "mais" so' pertence a esta folha quando e' contorno; no
        // interior ela e' a "menos" do vizinho, e e' preciso ir busca-la pela
        // geometria.
        double soma = 0.0;
        int n = 0;
        if (fm >= 0) { soma += v.em(d, fm); n++; }
        if (fp >= 0) { soma += v.em(d, fp); n++; }
        if (n == 0) continue;
        comp[d] = soma / n;
      }
      f.v.push_back(comp[0]);
      f.v.push_back(comp[1]);
      f.v.push_back(comp[2]);
    }
  campos_.push_back(std::move(f));
}

void SaidaVtk2d::tensorial(const std::string &nome, const CampoTensor &t)
{
  Campo f;
  f.nome = nome;
  f.componentes = 9;                 // o VTK quer 3x3
  f.v.reserve(static_cast<size_t>(m_->num_folhas()) * 9);
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      const double a00 = t.em(id, 0, 0), a01 = t.em(id, 0, 1);
      const double a10 = t.em(id, 1, 0), a11 = t.em(id, 1, 1);
      const double z[9] = {a00, a01, 0.0, a10, a11, 0.0, 0.0, 0.0, 0.0};
      for (double q : z) f.v.push_back(q);
    }
  campos_.push_back(std::move(f));
}

void SaidaVtk2d::fronteira(const std::string &nome,
                           const std::vector<std::vector<Ponto2>> &curvas)
{
  fronteiras_.push_back({nome, curvas});
}

void SaidaVtk2d::escreve(const std::string &prefixo, int quadro) const
{
  const int rank = m_->processo();
  char nome[512];
  std::snprintf(nome, sizeof nome, "%s-%06d-%04d.vtk", prefixo.c_str(), quadro,
                rank);
  std::FILE *f = std::fopen(nome, "w");
  if (f == nullptr)
    throw std::runtime_error("SaidaVtk2d: nao consegui abrir \"" +
                             std::string(nome) + "\" para escrita");

  const long n = m_->num_folhas();
  std::fprintf(f, "# vtk DataFile Version 3.0\n");
  std::fprintf(f, "t8Flow, quadro %d, processo %d\n", quadro, rank);
  std::fprintf(f, "ASCII\nDATASET UNSTRUCTURED_GRID\n");

  // QUATRO PONTOS POR CELULA, sem compartilhar vertices entre celulas.
  //
  // Gastar quatro vezes mais pontos e' deliberado: compartilhar exigiria
  // casar vertices por coordenada, e numa malha adaptada dois vertices
  // "iguais" de niveis diferentes nao sao exatamente iguais em ponto
  // flutuante.  O custo e' disco de visualizacao; o beneficio e' que o
  // arquivo nunca costura celulas que nao se tocam.
  std::fprintf(f, "POINTS %ld double\n", 4 * n);
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      const Ponto2 h = m_->lados_da_folha(b, i);
      const double x0 = c[0] - h[0] / 2, x1 = c[0] + h[0] / 2;
      const double y0 = c[1] - h[1] / 2, y1 = c[1] + h[1] / 2;
      std::fprintf(f, "%.17g %.17g 0\n%.17g %.17g 0\n%.17g %.17g 0\n%.17g %.17g 0\n",
                   x0, y0, x1, y0, x1, y1, x0, y1);
    }

  std::fprintf(f, "CELLS %ld %ld\n", n, 5 * n);
  for (long k = 0; k < n; k++)
    std::fprintf(f, "4 %ld %ld %ld %ld\n", 4 * k, 4 * k + 1, 4 * k + 2, 4 * k + 3);
  std::fprintf(f, "CELL_TYPES %ld\n", n);
  for (long k = 0; k < n; k++) std::fprintf(f, "9\n");    // VTK_QUAD

  std::fprintf(f, "CELL_DATA %ld\n", n);
  // O ID GLOBAL VAI JUNTO, sempre.  E' o que permite casar uma celula do VTK
  // com a mesma celula do despejo de precisao cheia -- sem ele, "a celula do
  // canto" e' identificacao por coordenada, e coordenada em ASCII ja' chega
  // arredondada.
  std::fprintf(f, "SCALARS id_global int 1\nLOOKUP_TABLE default\n");
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++)
      std::fprintf(f, "%ld\n", m_->id_global(b, i));

  for (const Campo &c : campos_) {
    if (static_cast<long>(c.v.size()) != n * c.componentes)
      throw std::runtime_error(
          "SaidaVtk2d: o campo \"" + c.nome + "\" tem " +
          std::to_string(c.v.size()) + " valores para " + std::to_string(n) +
          " celulas de " + std::to_string(c.componentes) + " componentes");

    if (c.componentes == 1)
      std::fprintf(f, "SCALARS %s double 1\nLOOKUP_TABLE default\n", c.nome.c_str());
    else if (c.componentes == 3)
      std::fprintf(f, "VECTORS %s double\n", c.nome.c_str());
    else
      std::fprintf(f, "TENSORS %s double\n", c.nome.c_str());

    for (long k = 0; k < n; k++) {
      for (int q = 0; q < c.componentes; q++)
        std::fprintf(f, "%.17g%s", c.v[static_cast<size_t>(k * c.componentes + q)],
                     (q + 1 == c.componentes) ? "" : " ");
      std::fprintf(f, "\n");
    }
  }
  std::fclose(f);

  // AS FRONTEIRAS EM ARQUIVO PROPRIO: segmentos e quadrilateros sao topologias
  // diferentes, e o VTK legado nao as mistura num arquivo sem ambiguidade.
  for (const auto &fr : fronteiras_) {
    std::snprintf(nome, sizeof nome, "%s-%s-%06d-%04d.vtk", prefixo.c_str(),
                  fr.first.c_str(), quadro, rank);
    std::FILE *g = std::fopen(nome, "w");
    if (g == nullptr)
      throw std::runtime_error("SaidaVtk2d: nao consegui abrir \"" +
                               std::string(nome) + "\"");

    long np = 0, nseg = 0;
    for (const auto &curva : fr.second) {
      np += static_cast<long>(curva.size());
      if (curva.size() >= 2) nseg += static_cast<long>(curva.size()) - 1;
    }
    std::fprintf(g, "# vtk DataFile Version 3.0\nt8Flow, fronteira %s, quadro %d\n"
                    "ASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %ld double\n",
                 fr.first.c_str(), quadro, np);
    for (const auto &curva : fr.second)
      for (const Ponto2 &p : curva)
        std::fprintf(g, "%.17g %.17g 0\n", p[0], p[1]);

    std::fprintf(g, "CELLS %ld %ld\n", nseg, 3 * nseg);
    long base = 0;
    for (const auto &curva : fr.second) {
      for (size_t q = 0; q + 1 < curva.size(); q++)
        std::fprintf(g, "2 %ld %ld\n", base + static_cast<long>(q),
                     base + static_cast<long>(q) + 1);
      base += static_cast<long>(curva.size());
    }
    std::fprintf(g, "CELL_TYPES %ld\n", nseg);
    for (long q = 0; q < nseg; q++) std::fprintf(g, "3\n");   // VTK_LINE
    std::fclose(g);
  }
}

}  // namespace t8flow
