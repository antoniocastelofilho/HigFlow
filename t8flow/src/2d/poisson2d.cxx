// Implementacao do Poisson.  A DERIVACAO do fechamento de parede e a razao de o
// operador ser um balanco de fluxos estao no `poisson2d.hxx`.  Aqui ficam as
// decisoes de implementacao:
//
//   1. a linha de cada celula e' acumulada num mapa `id -> coeficiente` antes de
//      ir ao solver, porque o MESMO id pode receber contribuicao de mais de uma
//      faceta -- o vizinho interior usado no fechamento de uma parede e' tambem
//      vizinho por outra face.  Somar e' o correto; sobrescrever perderia um termo
//      em silencio;
//   2. o termo conhecido (o valor de Dirichlet) e' movido para o lado direito com
//      sinal trocado, uma vez, no ponto em que e' calculado.  Acumular e trocar o
//      sinal depois e' onde o sinal se perde;
//   3. a montagem RECUSA o que esta versao nao sabe: faceta com dois vizinhos
//      (transicao 2:1), celula de parede sem vizinho interior, e contorno que nao
//      seja Dirichlet.  Recusar e' o que impede um resultado plausivel de sair de
//      um caso nao suportado.
#include "poisson2d.hxx"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace t8flow2d {

Poisson::Poisson(const Malha &m, ContornoPorRegiao contorno)
    : m_(&m), contorno_(std::move(contorno)),
      primeiro_(m.primeiro_id()), n_(m.num_folhas())
{
  if (contorno_.vazio())
    throw std::invalid_argument("Poisson: nenhuma condicao de contorno instalada");
  sol_.assign(static_cast<size_t>(n_), 0.0);
}

namespace {

//! (bloco, folha) de um id global.  O(folhas) -- caminho de montagem, nao laco
//! quente; dito para nao ser suposto rapido.
bool localiza(const Malha &m, long id, long *b, long *i)
{
  for (long bb = 0; bb < m.num_blocos(); bb++)
    for (long ii = 0; ii < m.num_folhas_do_bloco(bb); ii++)
      if (m.id_global(bb, ii) == id) { *b = bb; *i = ii; return true; }
  return false;
}

}  // namespace

void Poisson::resolve(Motor motor, double tol, int max_it)
{
  if (!f_) throw std::invalid_argument("Poisson::resolve: termo de fonte nao instalado");

  // 5 por linha no interior; o fechamento de tres pontos nao acrescenta id novo
  // (o vizinho interior ja' esta' na linha), mas folga nao custa.
  Solver s(motor, static_cast<int>(primeiro_), static_cast<std::size_t>(n_), 9);
  s.tolerancia(tol);
  s.max_iteracoes(max_it);

  for (long b = 0; b < m_->num_blocos(); b++) {
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const long id = m_->id_global(b, i);
      const Ponto2 c = m_->centroide(b, i);
      const Ponto2 h = m_->lados_da_folha(b, i);
      const double area = h[0] * h[1];

      // `int_V lap(u) = sum_faces du/dn |face|`, logo a linha e' a soma dos fluxos
      // e o lado direito e' `f(c) * area`.
      std::map<long, double> linha;
      double rhs = f_(c) * area;

      for (int d = 0; d < 2; d++) {
        const double hd = h[d];
        const double larg = h[1 - d];     // comprimento da faceta normal a d
        for (int lado = 0; lado < 2; lado++) {
          const int fa = 2 * d + lado;
          const std::vector<long> vz = m_->ids_vizinhos(b, i, fa);

          if (vz.size() > 1)
            throw std::runtime_error(
                "Poisson: a faceta " + std::to_string(fa) + " da celula " +
                std::to_string(id) + " tem " + std::to_string(vz.size()) +
                " vizinhos (transicao 2:1).  Esta versao monta so' malha UNIFORME:"
                " na face grossa o fluxo e' a soma dos dois lados finos, e a linha"
                " deixa de ser esta");

          if (vz.size() == 1) {
            // FACETA INTERNA: derivada central entre os dois centros, exata para
            // grau <= 3 em malha uniforme -- e' o que faz o interior nao limitar a
            // exatidao do grau 2.
            const double coef = larg / hd;
            linha[vz[0]] += coef;
            linha[id]    -= coef;
            continue;
          }

          // PAREDE.  O ponto onde a condicao e' avaliada e' o centro da faceta.
          const double sinal = (lado == 0) ? -1.0 : 1.0;
          Ponto2 parede = c;
          parede[d] += sinal * hd / 2.0;

          const Contorno &cond = contorno_.em(parede);

          if (cond.tipo == TipoContorno::Neumann) {
            // NEUMANN E' O CASO TRIVIAL NO BALANCO DE FLUXOS, e isso vale dizer:
            // o que a condicao da' e' exatamente `du/dn`, que e' o proprio fluxo.
            // Nao ha' aproximacao nenhuma -- nem de dois pontos, nem de tres -- e
            // nenhum coeficiente novo entra na linha.  E' o lado em que o arranjo
            // deslocado paga: na forma pontual, Neumann exigiria fantasma.
            //
            // Com a normal EXTERIOR, `du/dn` e' o fluxo para fora direto.
            rhs -= cond.g(parede) * larg;
            continue;
          }

          if (cond.tipo != TipoContorno::Dirichlet)
            throw std::runtime_error(
                "Poisson: a parede em (" + std::to_string(parede[0]) + ", " +
                std::to_string(parede[1]) + ") e' Robin, e esta versao monta"
                " Dirichlet e Neumann.  Robin mistura valor e derivada na MESMA"
                " faceta: a linha ganha um coeficiente que depende de alfa, beta e"
                " h, e montar isso como se fosse Neumann daria resultado plausivel"
                " e errado");
          const double g = cond.g(parede);

          if (fechamento_ == FechamentoParede::DoisPontos) {
            // `du/dn` para fora = (g - u_c)/(h/2).  Classico, segunda ordem global,
            // NAO exato para grau 2 -- esta' aqui como contraste do fechamento.
            const double coef = larg / (hd / 2.0);
            linha[id] -= coef;
            rhs -= coef * g;                 // o conhecido vai para o lado direito
            continue;
          }

          // TRES PONTOS: parede (0, g), centro (h/2, u_c), vizinho (3h/2, u_v).
          // `du/ds(0) = -8g/(3h) + 3u_c/h - u_v/(3h)` com s para DENTRO, logo o
          // fluxo para FORA e' o negativo disso.
          const int fa_oposta = 2 * d + (1 - lado);
          const std::vector<long> vo = m_->ids_vizinhos(b, i, fa_oposta);
          if (vo.size() != 1)
            throw std::runtime_error(
                "Poisson: a celula " + std::to_string(id) + " nao tem vizinho"
                " interior na direcao " + std::to_string(d) + " para o fechamento de"
                " tres pontos (dominio de uma celula de espessura?).  Use"
                " FechamentoParede::DoisPontos, que precisa de dois");

          linha[id]    += -3.0 / hd * larg;
          linha[vo[0]] += (1.0 / (3.0 * hd)) * larg;
          rhs          -= (8.0 * g / (3.0 * hd)) * larg;
        }
      }

      std::vector<int> cols;
      std::vector<double> vals;
      cols.reserve(linha.size());
      vals.reserve(linha.size());
      for (const auto &kv : linha) {
        cols.push_back(static_cast<int>(kv.first));
        vals.push_back(kv.second);
      }
      s.linha(static_cast<int>(id), cols, vals);
      s.lado_direito(static_cast<int>(id), rhs);
    }
  }

  s.monta();
  s.resolve();

  convergiu_ = s.convergiu();
  iteracoes_ = s.iteracoes();
  residuo_ = s.residuo();
  sol_ = s.solucao();
}

double Poisson::u(long id) const
{
  const long k = id - primeiro_;
  if (k < 0 || k >= static_cast<long>(sol_.size()))
    throw std::out_of_range("Poisson::u: id " + std::to_string(id) + " fora da faixa local");
  return sol_[static_cast<size_t>(k)];
}

double Poisson::erro_inf(const std::function<double(const Ponto2 &)> &exata) const
{
  double pior = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      pior = std::fmax(pior, std::fabs(u(m_->id_global(b, i)) - exata(c)));
    }
  return pior;
}

double Poisson::erro_l2(const std::function<double(const Ponto2 &)> &exata) const
{
  double soma = 0.0, area_total = 0.0;
  for (long b = 0; b < m_->num_blocos(); b++)
    for (long i = 0; i < m_->num_folhas_do_bloco(b); i++) {
      const Ponto2 c = m_->centroide(b, i);
      const double e = u(m_->id_global(b, i)) - exata(c);
      const double a = m_->area(b, i);
      soma += e * e * a;
      area_total += a;
    }
  return (area_total > 0.0) ? std::sqrt(soma / area_total) : 0.0;
}

}  // namespace t8flow2d
