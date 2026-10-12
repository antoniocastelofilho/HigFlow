// O DISCRETIZADOR DO T8FLOW 2D: aproximacao por minimos quadrados moveis.
//
// Fachada sobre o `wls` herdado da HiGTree (decisao 1).  O que o herdado faz
// esta' certo e e' independente de malha -- `wls.h` depende so' de `coord.h` e
// `types.h`, e a API inteira e' geometrica: pontos de suporte e um ponto alvo
// entram, pesos saem.  O que esta fachada acrescenta e' o que falta nele:
//
//   * TIPOS 2D CONCRETOS.  O herdado usa `Point`, que e' `real[DIM]` -- array
//     cru, nao copiavel, que NAO entra em container da STL (`std::vector<Point>`
//     nao compila, e o erro sai dentro do <bits/stl_construct.h> sem dizer isso).
//     Aqui o ponto e' `Ponto2`, um `std::array<double,2>`: copiavel, com o MESMO
//     layout, entao a travessia para o herdado e' conversao de TIPO e nao de
//     dados;
//   * TEMPO DE VIDA.  `wls_create`/`wls_destroy` viram construtor e destrutor;
//   * A CONDICAO DE USO AFIRMADA.  O ajuste precisa de um minimo de pontos para
//     a ordem pedida (`wls_num_min_points`), e abaixo dele o resultado e'
//     silenciosamente ruim.  `pesos()` RECUSA, em vez de aproximar.
//
// A SEPARACAO 2D/3D (decisao 7) e' do codigo NOVO: este arquivo e' 2D e nao tem
// `#if`.  As unidades herdadas continuam compiladas com `-DDIM=2` para este
// lado, porque a interface delas e' por macro -- isso nao muda, e nao e' custo
// novo: a HiGTree sempre construiu duas bibliotecas.
#ifndef T8FLOW_MLS2D_HXX
#define T8FLOW_MLS2D_HXX

#include <array>
#include <cstddef>
#include <vector>

struct wls_interpolator;   // o herdado, opaco aqui

namespace t8flow2d {

using Ponto2 = std::array<double, 2>;

class Mls {
 public:
  //! `ordem` e' o grau que o ajuste reproduz exatamente; `max_pontos` dimensiona
  //! o interpolador e e' o teto do suporte aceito.
  Mls(int ordem, int max_pontos);
  ~Mls();

  Mls(const Mls &) = delete;
  Mls &operator=(const Mls &) = delete;

  //! Minimo de pontos para a ordem pedida em 2D.
  static int min_pontos(int ordem);

  int ordem() const { return ordem_; }
  //! Teto de pontos do suporte, com que o interpolador foi dimensionado.
  int max_pontos() const { return max_pontos_; }

  //! Pesos de interpolacao em `alvo` a partir de `suporte`.
  //
  // Lanca `std::invalid_argument` se o suporte for menor que o minimo da ordem
  // ou maior que `max_pontos`.  Recusar e' deliberado: com suporte insuficiente o
  // ajuste devolve numeros, e numeros ruins em silencio sao o modo de falha que
  // este projeto mais paga.
  std::vector<double> pesos(const Ponto2 &alvo, const std::vector<Ponto2> &suporte) const;

  //! PESOS DE DERIVADA: `sum w_k u_k` aproxima `du/dx_dir` no alvo, exato para
  //! polinomio de grau <= `ordem`.
  //
  // A ALGEBRA E' A MESMA DO VALOR, com um vetor de avaliacao diferente -- e e' por
  // isso que isto cabe aqui em vez de ser um modulo novo.  O ajuste resolve
  //
  //     c = (B^T W B)^-1 B^T W u
  //
  // e o valor no alvo e' `b(x)^T c`; a derivada e' `(db/dx_dir)(x)^T c`.  Com a
  // base RELATIVA ao alvo, `b(x) = (1,0,0,...)` e `db/dx(x) = (0,1,0,...)`: troca-se
  // o lado direito `e_0` por `e_1`, e mais nada.
  //
  // POR QUE ISTO NAO VEM DO `wls` HERDADO: ele nao expoe a derivada -- `wls_calc`
  // devolve pesos de VALOR e a matriz fatorada fica interna.  Logo esta e' uma
  // implementacao NOVA (equacoes normais com o mesmo peso WLSMOVING), e por isso
  // os pesos de VALOR dela sao conferidos contra o herdado no teste: se as duas
  // implementacoes do mesmo ajuste divergem, nao ha' oraculo para a derivada.
  std::vector<double> pesos_derivada(const Ponto2 &alvo,
                                     const std::vector<Ponto2> &suporte, int dir) const;

  //! Pesos de VALOR pela implementacao propria -- existe para ser comparada com
  //! `pesos()`, que vem do nucleo herdado.  Mesma conta, dois caminhos.
  std::vector<double> pesos_proprios(const Ponto2 &alvo,
                                     const std::vector<Ponto2> &suporte) const;

  //! Valor interpolado em `alvo`, dados os valores nos pontos do suporte.
  //
  // Exige `valores.size() == suporte.size()`.  Existe para que o chamador nao
  // escreva o produto interno a mao -- foi onde a troca de indices entre dois
  // dominios se esconde no HiGFlow.
  double interpola(const Ponto2 &alvo, const std::vector<Ponto2> &suporte,
                   const std::vector<double> &valores) const;

 private:
  wls_interpolator *wls_ = nullptr;
  int ordem_ = 0;
  int max_pontos_ = 0;
};

}  // namespace t8flow2d

#endif
