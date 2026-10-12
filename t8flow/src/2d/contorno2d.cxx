// Implementacao das condicoes de fronteira.  A DERIVACAO de (A, B) para cada tipo
// esta' no `contorno2d.hxx`, com as contas -- aqui ficam as decisoes de
// implementacao, que sao duas:
//
//   1. UMA conta para os tres tipos, parametrizada por (alfa, beta).  Tres contas
//      separadas divergiriam, e divergir aqui significa Robin certo e Neumann
//      errado -- ou pior, os dois plausiveis;
//   2. o caso DEGENERADO e' recusado, nao aproximado.  `alfa*h/2 + beta == 0`
//      anula o denominador: a condicao nao determina o fantasma naquela celula, e
//      devolver qualquer numero ali seria inventar dado.
#include "contorno2d.hxx"

#include <cmath>
#include <string>

namespace t8flow2d {

void Contorno::afim(const Ponto2 &x_parede, double h, double *A, double *B) const
{
  if (!(h > 0.0))
    throw std::invalid_argument("Contorno::afim: h = " + std::to_string(h) +
                                " -- a distancia ao fantasma tem de ser positiva");
  if (!g)
    throw std::invalid_argument("Contorno::afim: condicao sem dado `g` instalado");

  const double gv = g(x_parede);

  // Dirichlet tem caminho PROPRIO no construtor de estencil -- o ponto de suporte
  // e' a parede, com valor conhecido.  Aqui ele aparece como o limite beta=0, que
  // e' a reflexao antissimetrica (A=-1, B=2g); serve de verificacao da formula, e
  // e' o que o teste compara.
  const double den = alfa * h / 2.0 + beta;
  if (std::fabs(den) <= 1e-300)
    throw std::invalid_argument(
        "Contorno::afim: condicao degenerada -- alfa*h/2 + beta = 0 com alfa=" +
        std::to_string(alfa) + ", beta=" + std::to_string(beta) + ", h=" +
        std::to_string(h) + "; ela nao determina o valor do fantasma");

  *A = (beta - alfa * h / 2.0) / den;
  *B = (gv * h) / den;
}

// ---------------------------------------------------------------------------
// ContornoPorRegiao
// ---------------------------------------------------------------------------

void ContornoPorRegiao::registra(std::string nome, OndeFica onde, Contorno c)
{
  if (!onde)
    throw std::invalid_argument("ContornoPorRegiao: regiao \"" + nome +
                                "\" sem predicado de posicao");
  if (!c.g)
    throw std::invalid_argument("ContornoPorRegiao: regiao \"" + nome +
                                "\" sem dado `g` instalado");
  for (const Regiao &r : regioes_)
    if (r.nome == nome)
      throw std::invalid_argument("ContornoPorRegiao: regiao \"" + nome +
                                  "\" registrada duas vezes -- nome repetido torna o"
                                  " diagnostico ambiguo justamente quando ele importa");
  regioes_.push_back(Regiao{std::move(nome), std::move(onde), std::move(c)});
}

void ContornoPorRegiao::padrao(Contorno c)
{
  if (!c.g)
    throw std::invalid_argument("ContornoPorRegiao::padrao: condicao sem dado `g`");
  padrao_ = std::move(c);
  tem_padrao_ = true;
}

const ContornoPorRegiao::Regiao *ContornoPorRegiao::acha(const Ponto2 &x,
                                                         std::string *ambigua) const
{
  // VARRE TODAS, nao para na primeira.  Parar na primeira esconderia a
  // sobreposicao e faria a condicao depender da ordem de registro -- e' o ponto
  // inteiro desta classe, e custa O(regioes), que sao poucas.
  const Regiao *achada = nullptr;
  for (const Regiao &r : regioes_) {
    if (!r.onde(x)) continue;
    if (achada != nullptr) {
      if (ambigua != nullptr) *ambigua = achada->nome + "\" e \"" + r.nome;
      return nullptr;
    }
    achada = &r;
  }
  return achada;
}

const Contorno &ContornoPorRegiao::em(const Ponto2 &x_parede) const
{
  std::string ambigua;
  const Regiao *r = acha(x_parede, &ambigua);
  if (r != nullptr) return r->cond;

  const std::string onde = "(" + std::to_string(x_parede[0]) + ", " +
                           std::to_string(x_parede[1]) + ")";
  if (!ambigua.empty())
    throw std::invalid_argument("ContornoPorRegiao: o ponto " + onde +
                                " esta' em DUAS regioes, \"" + ambigua +
                                "\" -- a sobreposicao faria a condicao depender da"
                                " ordem de registro");
  if (tem_padrao_) return padrao_;
  throw std::invalid_argument(
      "ContornoPorRegiao: nenhuma regiao cobre o ponto de parede " + onde +
      " e nao ha' padrao instalado.  Contorno sem condicao faz a conta rodar com a"
      " parede ausente -- resultado plausivel e errado");
}

std::string ContornoPorRegiao::nome_em(const Ponto2 &x_parede) const
{
  std::string ambigua;
  const Regiao *r = acha(x_parede, &ambigua);
  if (r != nullptr) return r->nome;
  if (!ambigua.empty()) return "(ambigua: " + ambigua + ")";
  return tem_padrao_ ? "(padrao)" : "(nenhuma)";
}

}  // namespace t8flow2d
