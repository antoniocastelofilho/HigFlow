// Implementacao do construtor de estencil.  O PORQUE da forma -- `rhs` separado,
// parada ao fechar o anel, suporte local -- esta' no `estencil2d.hxx`.  Aqui ficam
// as decisoes de implementacao:
//
//   1. os aneis sao construidos por VIZINHANCA DE FACE, nao por raio geometrico.
//      Raio ignora a topologia: numa transicao 2:1 ele pega celulas que nao sao
//      vizinhas e perde as que sao;
//   2. o anel e' SEMPRE fechado por inteiro antes de decidir parar.  Decidir no
//      meio faz o suporte depender da ordem de visita;
//   3. o ponto de CONTORNO entra no suporte com id -1 e com o seu valor de
//      Dirichlet guardado a parte.  O peso dele multiplica o valor e vai para o
//      `rhs`; ele NAO aparece em `ids`;
//   4. o conjunto de ja'-vistos e' por ID, nao por posicao: comparar posicao em
//      ponto flutuante para decidir "ja' tenho esta celula" e' como se duplica
//      celula perto de interface de bloco.
#include "estencil2d.hxx"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace t8flow2d {

ConstrutorEstencil::ConstrutorEstencil(const Malha &m, int ordem, int max_pontos,
                                       int max_aneis)
    : m_(&m), ordem_(ordem), max_pontos_(max_pontos), max_aneis_(max_aneis),
      mls_(ordem, max_pontos)
{
  if (max_aneis < 1)
    throw std::invalid_argument("ConstrutorEstencil: max_aneis tem de ser >= 1");
}

namespace {

//! (bloco, folha) de um id global -- a malha indexa por par, o vizinho vem por id.
//
// Varredura: O(folhas).  E' aceitavel aqui porque a selecao de suporte ja' e' um
// caminho de montagem, nao de laco quente; quando for, isto vira um indice.  Dito
// para nao ser suposto rapido.
bool localiza(const Malha &m, long id, long *b, long *i)
{
  for (long bb = 0; bb < m.num_blocos(); bb++)
    for (long ii = 0; ii < m.num_folhas_do_bloco(bb); ii++)
      if (m.id_global(bb, ii) == id) { *b = bb; *i = ii; return true; }
  return false;
}

}  // namespace

ConstrutorEstencil::Suporte ConstrutorEstencil::suporte_em(const Ponto2 &x) const
{
  const long id0 = m_->id_do_ponto(x);
  if (id0 < 0)
    throw std::out_of_range("ConstrutorEstencil: ponto fora do dominio");

  const int minimo = Mls::min_pontos(ordem_);

  Suporte s;
  std::set<long> vistos;        // por ID, nao por posicao -- ver o cabecalho
  std::vector<long> fronteira;  // os ids do anel corrente

  long b0 = 0, i0 = 0;
  if (!localiza(*m_, id0, &b0, &i0))
    throw std::runtime_error("ConstrutorEstencil: id " + std::to_string(id0) +
                             " nao localizado na malha");

  vistos.insert(id0);
  fronteira.push_back(id0);
  s.pontos.push_back(m_->centroide(b0, i0));
  s.ids.push_back(id0);
  s.valores_contorno.push_back(0.0);

  for (int anel = 1; anel <= max_aneis_; anel++) {
    // O ANEL E' FECHADO POR INTEIRO antes de qualquer decisao de parada.  Parar no
    // meio faria o suporte -- e portanto os pesos -- depender da ordem de visita.
    std::vector<long> proximo;
    for (long id : fronteira) {
      long b = 0, i = 0;
      if (!localiza(*m_, id, &b, &i)) continue;
      for (int fa = 0; fa < 4; fa++) {
        const std::vector<long> vz = m_->ids_vizinhos(b, i, fa);
        if (vz.empty()) {
          // FACE DE CONTORNO: o ponto que entra e' o centro da FACETA, e o valor
          // dele e' conhecido.  Sem condicao instalada, o contorno nao contribui --
          // e o suporte vai faltar pontos, que e' melhor que contribuir com zero.
          if (!dirichlet_) continue;
          const Ponto2 c = m_->centroide(b, i);
          const Ponto2 h = m_->lados_da_folha(b, i);
          const int d = fa / 2;
          Ponto2 cf = c;
          cf[d] += (fa % 2 == 0 ? -h[d] / 2.0 : h[d] / 2.0);
          s.pontos.push_back(cf);
          s.ids.push_back(-1);                        // -1 = ponto de contorno
          s.valores_contorno.push_back(dirichlet_(cf));
          continue;
        }
        for (long v : vz) {
          if (vistos.count(v) != 0) continue;
          long vb = 0, vi = 0;
          if (!localiza(*m_, v, &vb, &vi)) continue;
          vistos.insert(v);
          proximo.push_back(v);
          s.pontos.push_back(m_->centroide(vb, vi));
          s.ids.push_back(v);
          s.valores_contorno.push_back(0.0);
        }
      }
    }
    s.aneis = anel;
    if (proximo.empty()) break;            // a malha acabou (dominio pequeno)
    fronteira = proximo;
    if (static_cast<int>(s.pontos.size()) >= minimo) break;   // anel FECHADO, agora para
  }

  // TETO: o ajuste nao aceita mais que `max_pontos`, e cortar o excesso tem de ser
  // por DISTANCIA ao alvo -- cortar pela ordem de insercao devolveria um suporte
  // que depende do percurso.
  if (static_cast<int>(s.pontos.size()) > max_pontos_) {
    std::vector<int> ordem(s.pontos.size());
    for (size_t k = 0; k < ordem.size(); k++) ordem[k] = static_cast<int>(k);
    const auto d2 = [&](int k) {
      const double dx = s.pontos[k][0] - x[0], dy = s.pontos[k][1] - x[1];
      return dx * dx + dy * dy;
    };
    std::stable_sort(ordem.begin(), ordem.end(),
                     [&](int a, int b) { return d2(a) < d2(b); });
    Suporte t;
    t.aneis = s.aneis;
    for (int k = 0; k < max_pontos_; k++) {
      t.pontos.push_back(s.pontos[ordem[k]]);
      t.ids.push_back(s.ids[ordem[k]]);
      t.valores_contorno.push_back(s.valores_contorno[ordem[k]]);
    }
    s = t;
  }

  if (static_cast<int>(s.pontos.size()) < minimo)
    throw std::runtime_error(
        "ConstrutorEstencil: suporte de " + std::to_string(s.pontos.size()) +
        " pontos em " + std::to_string(max_aneis_) + " aneis, minimo " +
        std::to_string(minimo) + " para ordem " + std::to_string(ordem_) +
        (dirichlet_ ? "" : " (nenhuma condicao de Dirichlet instalada: o contorno"
                           " nao contribuiu)"));

  for (size_t a = 0; a < s.pontos.size(); a++)
    for (size_t b = a + 1; b < s.pontos.size(); b++) {
      const double dx = s.pontos[a][0] - s.pontos[b][0];
      const double dy = s.pontos[a][1] - s.pontos[b][1];
      s.diametro = std::fmax(s.diametro, std::sqrt(dx * dx + dy * dy));
    }
  return s;
}

Estencil ConstrutorEstencil::em(const Ponto2 &x) const
{
  const Suporte s = suporte_em(x);
  const std::vector<double> w = mls_.pesos(x, s.pontos);

  // O `rhs` RECEBE a contribuicao dos pontos de contorno, e eles NAO entram em
  // `ids`.  Quem consome o estencil monta a linha do sistema com `ids`/`pesos` e
  // soma `rhs` ao lado direito -- sem precisar saber qual id era parede.
  Estencil e;
  e.rhs = 0.0;
  for (size_t k = 0; k < w.size(); k++) {
    if (s.ids[k] < 0) e.rhs += w[k] * s.valores_contorno[k];
    else { e.ids.push_back(s.ids[k]); e.pesos.push_back(w[k]); }
  }
  return e;
}

void ConstrutorEstencil::em_lote(const std::vector<Ponto2> &alvos, LoteEstenceis &lote,
                                 std::vector<Estencil> &saida) const
{
  // A SELECAO DO SUPORTE FICA NA CPU, sempre: e' navegacao de malha -- busca por
  // ponto, vizinhanca, conjuntos -- e nao aritmetica densa.  O que vai para a
  // placa sao os PESOS, que e' onde esta' a conta.  Levar a navegacao para a GPU
  // seria portar a malha, nao acelerar o calculo.
  std::vector<Suporte> sups;
  sups.reserve(alvos.size());
  std::vector<ConsultaEstencil> consultas;
  consultas.reserve(alvos.size());

  for (const Ponto2 &x : alvos) {
    Suporte s = suporte_em(x);
    ConsultaEstencil c;
    c.alvo = x;
    c.suporte = s.pontos;
    // Os ids vao com o -1 do contorno PRESERVADO: o lote nao precisa saber o que
    // -1 significa, e a traducao para `rhs` acontece aqui, depois.
    c.ids = s.ids;
    consultas.push_back(std::move(c));
    sups.push_back(std::move(s));
  }

  std::vector<Estencil> brutos;
  lote.calcula(consultas, brutos);

  saida.assign(alvos.size(), Estencil{});
  for (size_t q = 0; q < brutos.size(); q++) {
    const Suporte &s = sups[q];
    Estencil e;
    e.rhs = 0.0;
    for (size_t k = 0; k < brutos[q].pesos.size(); k++) {
      if (s.ids[k] < 0) e.rhs += brutos[q].pesos[k] * s.valores_contorno[k];
      else { e.ids.push_back(s.ids[k]); e.pesos.push_back(brutos[q].pesos[k]); }
    }
    saida[q] = std::move(e);
  }
}

}  // namespace t8flow2d
