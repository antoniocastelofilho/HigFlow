// AJUSTE POR MINIMOS QUADRADOS MOVEIS EM TRES DIMENSOES.
//
// Dada uma nuvem de pontos e um alvo, devolve os PESOS `w_k` tais que
// `soma w_k u_k` aproxima um funcional de `u` no alvo -- o valor, uma derivada
// primeira, ou o laplaciano.  Nao devolve o polinomio: devolve os pesos, porque
// e' deles que a linha da matriz e' feita.
//
// O QUE MUDA DE DUAS PARA TRES DIMENSOES, e nao e' so' um indice a mais:
//
//   A BASE EXPLODE.  Monomios por ordem:
//
//                      2D      3D
//        ordem 1        3       4
//        ordem 2        6      10
//        ordem 3       10      20
//        ordem 4       15      35
//
//   O suporte minimo e' o numero de monomios, logo a ordem 3 em 3D pede 20
//   pontos onde em 2D pedia 10 -- e a ordem 4 pede 35.  Com o suporte crescendo
//   por aneis de vizinhos de face, isso significa mais aneis, que significa mais
//   camadas de franja, que em tres dimensoes custa o cubo do raio.
//
//   A MATRIZ NORMAL E' 35x35 na ordem 4, contra 15x15 em duas dimensoes, e o
//   custo de resolve-la e' cubico nisso.  O ajuste deixa de ser barato.
//
// O QUE NAO MUDA, e vem com os consertos que duas dimensoes custou descobrir:
//
//   A BASE E' ADIMENSIONALIZADA pelo raio do suporte.  Sem isso a matriz normal
//   da ordem 4 tem condicionamento da ordem de `h^-8` e o ajuste degrada EM
//   SILENCIO -- devolvendo pesos, nao erro;
//
//   O REBAIXAMENTO DE ORDEM E' REGISTRADO.  Quando o suporte nao sustenta a
//   ordem pedida, o ajuste cai para a que sustenta e um contador conta.  A
//   versao que devolvia um vetor de ZEROS produziu, em duas dimensoes, erro de
//   3,3e-01 num nivel cercado por 1e-10 nos vizinhos;
//
//   O PESO DO LAPLACIANO E' REGULARIZADO.  No laplaciano pontual o alvo e' o
//   centro da celula e o centro ESTA' no suporte -- sempre.  Com o peso herdado
//   `1/(d^2 + EPSMACH)` esse ponto recebe peso da ordem de 1e32 e o ajuste passa
//   a INTERPOLAR, o que e' certo para o VALOR e fatal para a segunda derivada.

#ifndef T8FLOW_3D_MLS3D_HXX
#define T8FLOW_3D_MLS3D_HXX

#include <array>
#include <vector>

namespace t8flow3d {

using Ponto3 = std::array<double, 3>;

class Mls3 {
 public:
  //! `ordem` e' o grau do polinomio reproduzido exatamente; `max_pontos` e' o
  //! teto do suporte, e tem de alcancar o minimo da ordem.
  Mls3(int ordem, int max_pontos);
  ~Mls3();

  Mls3(const Mls3 &) = delete;
  Mls3 &operator=(const Mls3 &) = delete;

  //! Quantos pontos a ordem exige em TRES dimensoes: 4, 10, 20, 35.
  static int min_pontos(int ordem);

  int ordem() const { return ordem_; }
  int max_pontos() const { return max_pontos_; }

  //! Pesos de VALOR pelo nucleo HERDADO (decomposicao QR).
  //
  // Mantido ao lado da implementacao propria de proposito: os dois resolvem o
  // mesmo problema por caminhos diferentes, e a concordancia entre eles e' um
  // oraculo que nenhum dos dois produz sozinho.  Em duas dimensoes eles
  // concordam em 1,249e-16.
  std::vector<double> pesos(const Ponto3 &alvo, const std::vector<Ponto3> &suporte) const;

  //! Pesos de VALOR pela implementacao propria.
  std::vector<double> pesos_proprios(const Ponto3 &alvo,
                                     const std::vector<Ponto3> &suporte) const;

  //! Pesos da derivada `d/dx_dir`, com `dir` em 0, 1, 2.
  std::vector<double> pesos_derivada(const Ponto3 &alvo,
                                     const std::vector<Ponto3> &suporte, int dir) const;

  //! Pesos do LAPLACIANO pontual.  Exige ordem >= 2 -- um ajuste linear nao tem
  //! segunda derivada, e devolver zero seria responder a uma pergunta mal posta.
  std::vector<double> pesos_laplaciano(const Ponto3 &alvo,
                                       const std::vector<Ponto3> &suporte) const;

  //! Interpola um campo amostrado no suporte.
  double interpola(const Ponto3 &alvo, const std::vector<Ponto3> &suporte,
                   const std::vector<double> &valores) const;

  //! Contadores: quantos ajustes, quantos REBAIXARAM a ordem, quantos foram
  //! nulos.  Zero rebaixamentos e' o que se quer poder afirmar -- um esquema que
  //! diz "ordem 3" com rebaixamento entrega ordem menor em parte do dominio.
  long ajustes() const { return ajustes_; }
  long rebaixamentos() const { return rebaixos_; }
  long ajustes_nulos() const { return nulos_; }
  void zera_contadores() const { ajustes_ = rebaixos_ = nulos_ = 0; }

 private:
  const std::vector<double> &conta(const std::vector<double> &w) const;

  int ordem_ = 0;
  int max_pontos_ = 0;
  void *wls_ = nullptr;      //!< `wls_interpolator*` do nucleo herdado

  mutable int ord_usada_ = 0;
  mutable long ajustes_ = 0;
  mutable long rebaixos_ = 0;
  mutable long nulos_ = 0;
};

}  // namespace t8flow3d

#endif  // T8FLOW_3D_MLS3D_HXX
