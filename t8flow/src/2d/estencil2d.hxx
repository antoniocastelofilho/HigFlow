// O ESTENCIL COM A ORDEM DESEJADA, como na HiGTree.
//
// La' a interface e' `sd_set_interpolator_order` + `sd_get_stencil(sd, centro, x,
// alpha, stn)`, devolvendo um `sim_stencil` com `ids[]`, `vals[]` (os pesos) e um
// `rhs` SEPARADO para o que nao se associa a nenhuma celula -- o valor de
// Dirichlet.  Aqui a forma e' a mesma, e a separacao do `rhs` tambem: guardar o
// valor de contorno junto com os pesos obriga o chamador a saber qual id e' "a
// parede", e e' assim que o fechamento se perde.
//
// O QUE ESTE MODULO ACRESCENTA ao que ja' existia -- o `Mls` sabia calcular pesos
// de um suporte DADO, e o `LoteEstenceis` sabia faze-lo em lote -- e' a ESCOLHA DO
// SUPORTE:
//
//   1. acha a celula que contem o ponto (busca geometrica da malha);
//   2. cresce por ANEIS de vizinhanca ate' ter pontos para a ordem pedida;
//   3. PARA AO FECHAR O ANEL, nunca no meio dele.
//
// O PASSO 3 NAO E' DETALHE.  Parar no meio de um anel faz o suporte depender da
// ORDEM em que os vizinhos foram visitados, e a resposta passa a mudar com o
// percurso -- que e' exatamente como a clausula de desempate da HiGTree passava
// por sorte.  E crescer SEM PARADA e' pior: busca em largura sem limite devolve a
// malha inteira (medido na HiGTree: 64 pontos numa malha de 64), o ajuste
// reproduz polinomio mesmo assim, e o teste passa enquanto a topologia deixou de
// importar.  Daqui sai um suporte LOCAL, e a clausula F3 afirma isso.
//
// O FECHAMENTO DE CONTORNO.  Quando o anel alcanca o contorno do dominio, os
// pontos que faltam vem das FACETAS de contorno, nao de celulas: o valor ali e'
// conhecido (Dirichlet), entao o peso daquele ponto multiplica o valor e vai para
// o `rhs`.  O suporte tem dois tipos de ponto, e a diferenca nao esta' na
// geometria -- esta' em quem responde pelo valor.
//
// LIMITE DE HOJE, declarado: so' Dirichlet.  Neumann precisa de derivada na
// faceta, que muda a LINHA do ajuste e nao so' o lado direito -- e' trabalho
// proprio, nao um `else` aqui.
#ifndef T8FLOW_ESTENCIL2D_HXX
#define T8FLOW_ESTENCIL2D_HXX

#include <functional>
#include <vector>

#include "campos2d.hxx"
#include "lote2d.hxx"
#include "malha2d.hxx"
#include "mls2d.hxx"

namespace t8flow2d {

//! Condicao de contorno de Dirichlet: o valor na faceta de contorno.
//
// Uma funcao do ponto, nao uma constante: valor de contorno CONSTANTE e' o caso
// em que quase todo defeito de fechamento fica invisivel -- com parede de
// nao-deslizamento, fechar contra a parede errada da' o mesmo numero.  Foi assim
// que nove defeitos de fechamento sobreviveram cinco anos na HiGTree.
using ValorDeContorno = std::function<double(const Ponto2 &)>;

class ConstrutorEstencil {
 public:
  //! `ordem` e' o grau reproduzido exatamente.  `max_pontos` e' o teto do
  //! suporte; `max_aneis` limita o crescimento -- sem ele, um ponto numa regiao
  //! grossa cercada de contorno cresceria indefinidamente.
  ConstrutorEstencil(const Malha &m, int ordem, int max_pontos, int max_aneis = 4);

  //! Instala a condicao de Dirichlet.  Sem ela, o contorno nao contribui e o
  //! suporte perto da parede fica com menos pontos -- o que faz `em()` LANCAR em
  //! vez de devolver um ajuste de ordem menor sem avisar.
  void dirichlet(ValorDeContorno v) { dirichlet_ = std::move(v); }

  int ordem() const { return ordem_; }
  int max_pontos() const { return max_pontos_; }

  //! O estencil em `x`: ids e pesos das celulas, e `rhs` com a contribuicao do
  //! contorno.
  //
  // Lanca `std::out_of_range` se `x` estiver fora do dominio, e
  // `std::runtime_error` se o suporte nao alcancar o minimo da ordem dentro de
  // `max_aneis` -- recusar e' deliberado: um ajuste de ordem menor que a pedida
  // devolve numeros plausiveis e errados.
  Estencil em(const Ponto2 &x) const;

  //! O mesmo para muitos pontos, por um `LoteEstenceis` -- que e' onde a GPU
  //! entra.  O suporte continua sendo escolhido na CPU (e' navegacao de malha,
  //! nao aritmetica densa); o que vai para a placa sao os PESOS.
  void em_lote(const std::vector<Ponto2> &alvos, LoteEstenceis &lote,
               std::vector<Estencil> &saida) const;

  //! Diagnostico do suporte escolhido em `x`: quantos pontos de celula, quantos
  //! de contorno, quantos aneis, e o diametro.  Existe para a clausula F3 -- "o
  //! suporte cai na grade e e' local" -- poder ser afirmada de fora.
  struct Suporte {
    std::vector<Ponto2> pontos;
    std::vector<long> ids;        //!< -1 para ponto de CONTORNO
    std::vector<double> valores_contorno;  //!< so' para os de contorno; 0 nos outros
    int aneis = 0;
    double diametro = 0.0;
  };
  Suporte suporte_em(const Ponto2 &x) const;

 private:
  const Malha *m_;
  int ordem_;
  int max_pontos_;
  int max_aneis_;
  ValorDeContorno dirichlet_;
  mutable Mls mls_;
};

}  // namespace t8flow2d

#endif
