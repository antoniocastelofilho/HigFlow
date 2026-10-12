// Implementacao do arranjo deslocado.  O PORQUE de cada propriedade morar onde
// mora, e a regra de posse das facetas, estao no `campos2d.hxx`.  Aqui ficam as
// decisoes de implementacao:
//
//   1. os campos de centro sao indexados pelo id GLOBAL da folha, deslocado por
//      `primeiro_id()` -- o MESMO id do solver.  Nao ha' indice proprio de campo,
//      e isso e' deliberado: indice paralelo e' como os dois dominios do VOF
//      passaram anos trocando posicao entre si;
//   2. a numeracao das facetas e' pre-computada na construcao, em O(folhas), e
//      consultada em O(1).  Recomputar por consulta seria O(folhas) dentro de
//      laco quente;
//   3. o criterio de "face mais e' contorno" e' `num_vizinhos == 0`, que vale em
//      np=1 -- o limite esta' dito no cabecalho.
#include "campos2d.hxx"

#include <stdexcept>
#include <string>

namespace t8flow2d {

// ---------------------------------------------------------------------------
// CampoEscalar
// ---------------------------------------------------------------------------

CampoEscalar::CampoEscalar(const Malha &m)
    : m_(&m), primeiro_(m.primeiro_id()), v_(static_cast<size_t>(m.num_folhas()), 0.0)
{
}

double &CampoEscalar::em(long id)
{
  const long k = id - primeiro_;
  if (k < 0 || k >= static_cast<long>(v_.size()))
    throw std::out_of_range("CampoEscalar::em: id " + std::to_string(id) + " fora de [" +
                            std::to_string(primeiro_) + ", " +
                            std::to_string(primeiro_ + static_cast<long>(v_.size())) + ")");
  return v_[static_cast<size_t>(k)];
}

double CampoEscalar::em(long id) const
{
  return const_cast<CampoEscalar *>(this)->em(id);
}

double CampoEscalar::em_ponto(const Ponto2 &p) const
{
  const long id = m_->id_do_ponto(p);
  if (id < 0)
    throw std::out_of_range("CampoEscalar::em_ponto: ponto fora do dominio");
  return em(id);
}

void CampoEscalar::preenche(double valor)
{
  for (double &x : v_) x = valor;
}

// ---------------------------------------------------------------------------
// CampoTensor
// ---------------------------------------------------------------------------

CampoTensor::CampoTensor(const Malha &m)
{
  comp_.reserve(4);
  for (int k = 0; k < 4; k++) comp_.emplace_back(m);
}

static int indice_componente(int lin, int col)
{
  if (lin < 0 || lin > 1 || col < 0 || col > 1)
    throw std::out_of_range("CampoTensor: componente (" + std::to_string(lin) + ", " +
                            std::to_string(col) + ") -- em 2D sao 2x2");
  return 2 * lin + col;
}

double &CampoTensor::em(long id, int lin, int col)
{
  return comp_[static_cast<size_t>(indice_componente(lin, col))].em(id);
}

double CampoTensor::em(long id, int lin, int col) const
{
  return comp_[static_cast<size_t>(indice_componente(lin, col))].em(id);
}

CampoEscalar &CampoTensor::componente(int lin, int col)
{
  return comp_[static_cast<size_t>(indice_componente(lin, col))];
}

const CampoEscalar &CampoTensor::componente(int lin, int col) const
{
  return comp_[static_cast<size_t>(indice_componente(lin, col))];
}

void CampoTensor::preenche(double valor)
{
  for (CampoEscalar &c : comp_) c.preenche(valor);
}

// ---------------------------------------------------------------------------
// DominioFacetas
// ---------------------------------------------------------------------------

DominioFacetas::DominioFacetas(const Malha &m) : m_(&m)
{
  const long n = m.num_folhas();
  const long p0 = m.primeiro_id();

  for (int d = 0; d < 2; d++) {
    menos_[d].assign(static_cast<size_t>(n), -1);
    mais_[d].assign(static_cast<size_t>(n), -1);
  }

  // UMA passada pelas folhas, na ordem da curva.  Para cada direcao: a faceta
  // "menos" sempre existe e e' desta folha; a "mais" so' existe se for contorno
  // do dominio -- se tiver vizinho, pertence a ele (e' a "menos" dele).
  for (long b = 0; b < m.num_blocos(); b++) {
    for (long i = 0; i < m.num_folhas_do_bloco(b); i++) {
      const long k = m.id_global(b, i) - p0;
      const Ponto2 c = m.centroide(b, i);
      const Ponto2 h = m.lados_da_folha(b, i);

      for (int d = 0; d < 2; d++) {
        // lado MENOS: sempre desta folha
        {
          Ponto2 cf = c;
          cf[d] -= h[d] / 2.0;
          menos_[d][static_cast<size_t>(k)] = static_cast<long>(centros_[d].size());
          centros_[d].push_back(cf);
          // a medida de uma faceta de normal d e' o lado TRANSVERSAL
          medidas_[d].push_back(h[1 - d]);
        }
        // lado MAIS: so' se nao houver vizinho, isto e', se for contorno
        if (m.num_vizinhos(b, i, 2 * d + 1) == 0) {
          Ponto2 cf = c;
          cf[d] += h[d] / 2.0;
          mais_[d][static_cast<size_t>(k)] = static_cast<long>(centros_[d].size());
          centros_[d].push_back(cf);
          medidas_[d].push_back(h[1 - d]);
        }
      }
    }
  }
}

long DominioFacetas::num_facetas(int d) const
{
  if (d < 0 || d > 1)
    throw std::out_of_range("DominioFacetas: direcao " + std::to_string(d) + " -- em 2D sao 2");
  return static_cast<long>(centros_[d].size());
}

long DominioFacetas::id_faceta(long b, long i, int d, int lado) const
{
  if (d < 0 || d > 1)
    throw std::out_of_range("DominioFacetas::id_faceta: direcao " + std::to_string(d));
  if (lado < 0 || lado > 1)
    throw std::out_of_range("DominioFacetas::id_faceta: lado " + std::to_string(lado) +
                            " -- 0 e' menos, 1 e' mais");
  const long k = m_->id_global(b, i) - m_->primeiro_id();
  if (k < 0 || k >= static_cast<long>(menos_[d].size()))
    throw std::out_of_range("DominioFacetas::id_faceta: folha fora da faixa local");
  return (lado == 0) ? menos_[d][static_cast<size_t>(k)] : mais_[d][static_cast<size_t>(k)];
}

Ponto2 DominioFacetas::centro(int d, long id) const
{
  if (d < 0 || d > 1)
    throw std::out_of_range("DominioFacetas::centro: direcao " + std::to_string(d));
  if (id < 0 || id >= num_facetas(d))
    throw std::out_of_range("DominioFacetas::centro: faceta " + std::to_string(id) +
                            " fora de [0, " + std::to_string(num_facetas(d)) + ")");
  return centros_[d][static_cast<size_t>(id)];
}

double DominioFacetas::medida(int d, long id) const
{
  if (d < 0 || d > 1)
    throw std::out_of_range("DominioFacetas::medida: direcao " + std::to_string(d));
  if (id < 0 || id >= num_facetas(d))
    throw std::out_of_range("DominioFacetas::medida: faceta " + std::to_string(id) +
                            " fora de [0, " + std::to_string(num_facetas(d)) + ")");
  return medidas_[d][static_cast<size_t>(id)];
}

// ---------------------------------------------------------------------------
// CampoVetorial
// ---------------------------------------------------------------------------

CampoVetorial::CampoVetorial(const DominioFacetas &fd) : fd_(&fd)
{
  for (int d = 0; d < 2; d++)
    v_[d].assign(static_cast<size_t>(fd.num_facetas(d)), 0.0);
}

double &CampoVetorial::em(int d, long id)
{
  if (d < 0 || d > 1)
    throw std::out_of_range("CampoVetorial::em: direcao " + std::to_string(d));
  if (id < 0 || id >= static_cast<long>(v_[d].size()))
    throw std::out_of_range("CampoVetorial::em: faceta " + std::to_string(id) +
                            " fora de [0, " + std::to_string(v_[d].size()) +
                            ") na direcao " + std::to_string(d));
  return v_[d][static_cast<size_t>(id)];
}

double CampoVetorial::em(int d, long id) const
{
  return const_cast<CampoVetorial *>(this)->em(d, id);
}

void CampoVetorial::preenche(double valor)
{
  for (int d = 0; d < 2; d++)
    for (double &x : v_[d]) x = valor;
}

}  // namespace t8flow2d
