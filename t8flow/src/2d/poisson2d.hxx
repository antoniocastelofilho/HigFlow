// A EQUACAO DE POISSON 2D, para o metodo das SOLUCOES MANUFATURADAS.
//
//     laplaciano(u) = f   no dominio,   com condicao de contorno por regiao
//
// MONTADA COMO BALANCO DE FLUXOS, que e' o que o arranjo deslocado pede: pela
// forma integral, `int_V lap(u) = sum_faces du/dn * |face|`, entao a linha de cada
// celula e' a soma dos fluxos pelas SUAS facetas igualada a `f(c) * area`.  Nao e'
// o laplaciano pontual montado no centro -- e essa diferenca e' a razao de a
// clausula C5 existir: no centro da celula o estencil MLS e' a IDENTIDADE (peso
// 1,0 nela, rhs ~0), e montar operador ali daria uma identidade, nao um
// laplaciano.
//
// O FECHAMENTO NA PAREDE E' DE SEGUNDA ORDEM, E ISSO NAO E' REFINAMENTO -- e' o
// que o pedido exige.  Para um `u` de grau 2 sair com erro de MAQUINA, a derivada
// normal na parede tem de ser exata para grau 2, e isso precisa de TRES pontos:
//
//     a parede (distancia 0, valor g conhecido)
//     o centro da celula       (distancia h/2)
//     o centro do vizinho      (distancia 3h/2)
//
// Pela interpolacao de Lagrange nesses tres pontos, a derivada na parede, no
// sentido PARA DENTRO, e':
//
//     du/ds(0) = -8g/(3h) + 3u_c/h - u_viz/(3h)
//
// (a soma dos coeficientes e' -8/3 + 3 - 1/3 = 0, como tem de ser: a derivada de
// uma constante e' zero -- essa conferencia pega erro de sinal e de escala.)
//
// A ALTERNATIVA DE DOIS PONTOS -- `(u_c - g)/(h/2)` -- da' a derivada no ponto
// MEDIO entre a parede e o centro, nao NA parede: erro O(h) no fluxo e O(h^2) na
// solucao.  E' o fechamento classico de volumes finitos, e e' suficiente para
// segunda ordem GLOBAL, mas NAO da' erro de maquina para grau 2.  O teste mede os
// dois, e o contraste entre eles e' o que mostra que o fechamento esta' sendo
// exercitado.
//
// NEUMANN E' O CASO TRIVIAL AQUI, e vale dizer por que: no balanco de fluxos, o
// que a condicao da' e' exatamente `du/dn` -- o proprio fluxo.  Nenhuma
// aproximacao, nenhum coeficiente novo na linha.  E' o lado em que o arranjo
// deslocado paga o que a forma pontual cobraria com fantasma.  ROBIN nao: ele
// mistura valor e derivada na MESMA faceta, a linha ganha coeficiente que depende
// de alfa, beta e h, e esta versao RECUSA em vez de montar como se fosse Neumann.
//
// E ATENCAO ao problema so' com Neumann: ele e' SINGULAR -- a solucao fica a menos
// de uma constante, e o solver nao converge ou converge para qualquer uma.  Isso
// nao e' defeito da montagem, e' propriedade do problema; o exemplo avisa.
//
// A TRANSICAO 2:1 ESTA' MONTADA, pela MESMA formula da faceta uniforme: para cada
// vizinho, a face compartilhada tem o tamanho do MENOR dos dois e a distancia e' a
// soma dos meios-passos.  Uma formula e nao tres, porque tres divergem -- e
// divergir ali quebra a CONSERVACAO sem quebrar nada visivel.  A igualdade
// conferida:
//
//     grossa->fina:  larg_C/2 / (3h_C/4)  =  (2/3) larg_C/h_C
//     fina->grossa:  larg_F   / (3h_F/2)  =  (2/3) larg_C/h_C
//
// Logo o que sai de uma entra na outra com o mesmo valor, e a matriz fica
// simetrica.
//
// O QUE ISSO NAO PROMETE, e a analise da quadratura ja' diz por que: na face de
// transicao os centros NAO estao alinhados com a normal -- o centro da sub-face
// fina esta' deslocado um quarto de celula do centro da grossa --, e a distancia
// entre centros (3h/4) nao e' a distancia do centro da grossa a' face (h/2).  Os
// termos de deslocamento das DUAS sub-faces se cancelam por simetria, mas os de
// curvatura NAO.  Entao nao se espera exatidao para grau 2 na malha adaptada, e
// sim segunda ordem -- e o teste MEDE em vez de supor.
#ifndef T8FLOW_POISSON2D_HXX
#define T8FLOW_POISSON2D_HXX

#include <functional>
#include <map>
#include <vector>

#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "solver2d.hxx"

namespace t8flow2d {

//! Como a derivada normal e' fechada na parede.
//
// Existe como OPCAO, e nao como escolha fixa, porque o contraste entre as duas e'
// o oraculo do fechamento: a de tres pontos tem de dar erro de maquina para grau
// 2, e a de dois pontos NAO -- se as duas derem o mesmo, o fechamento nao esta'
// influenciando o resultado.
enum class FechamentoParede {
  TresPontos,   //!< exato para grau 2; e' o padrao
  DoisPontos,   //!< `(u_c - g)/(h/2)`: classico, segunda ordem global, nao exato
};

class Poisson {
 public:
  Poisson(const Malha &m, ContornoPorRegiao contorno);

  //! O termo de fonte `f`.  Na manufatura, `f = lap(u_exata)`, derivado a mao.
  void fonte(std::function<double(const Ponto2 &)> f) { f_ = std::move(f); }

  void fechamento(FechamentoParede fp) { fechamento_ = fp; }
  FechamentoParede fechamento() const { return fechamento_; }

  //! Monta e resolve.  `tol` vai para o solver linear; o residuo relatado por ele
  //! e' o PISO do erro que se pode exigir da solucao -- pedir erro de maquina com
  //! tolerancia frouxa seria medir o solver, nao o esquema.
  void resolve(Motor motor = Motor::Petsc, double tol = 1e-14, int max_it = 20000);

  //! A solucao na celula de id global `id`.
  double u(long id) const;

  bool convergiu() const { return convergiu_; }
  int iteracoes() const { return iteracoes_; }
  double residuo() const { return residuo_; }

  //! Erro contra a solucao exata, na norma do MAXIMO.
  double erro_inf(const std::function<double(const Ponto2 &)> &exata) const;

  //! Erro na norma 2 DISCRETA, ponderada por area: `sqrt(sum e^2 A / sum A)`.
  //
  // Ponderada por area e nao media simples: em malha nao-uniforme a media simples
  // daria peso igual a celulas de tamanhos diferentes, e o numero deixaria de ser
  // uma norma de funcao.  Em malha uniforme as duas coincidem, e e' de proposito
  // que a formula ja' esta' certa para quando a malha deixar de ser.
  double erro_l2(const std::function<double(const Ponto2 &)> &exata) const;

  //! Incognitas do sistema -- uma por celula.
  long num_incognitas() const { return n_; }

  //! `max |A*u_exata - b| / area`: o ORACULO DIRETO DO OPERADOR, sem o solver.
  //
  // Mede a CONSISTENCIA do esquema.  O erro da SOLUCAO mistura consistencia,
  // condicionamento e tolerancia do solver, e quando ele e' grande nao se sabe
  // qual dos tres -- foi o que faltou quando o esquema de transicao 2:1 deu ordem
  // negativa.  Aqui, se o esquema reproduz a exata, isto e' zero a menos de
  // arredondamento, qualquer que seja o solver.
  double residuo_da_exata(const std::function<double(const Ponto2 &)> &exata) const;

 private:
  const Malha *m_;
  ContornoPorRegiao contorno_;
  std::function<double(const Ponto2 &)> f_;
  FechamentoParede fechamento_ = FechamentoParede::TresPontos;

  //! Monta a linha de uma celula.  `indice` e' o mapa id -> (bloco, folha), opaco
  //! aqui para o cabecalho nao depender do tipo interno.
  void monta_linha(long b, long i, const void *indice,
                   std::map<long, double> *linha, double *rhs) const;

  long primeiro_ = 0;
  long n_ = 0;
  std::vector<double> sol_;
  bool convergiu_ = false;
  int iteracoes_ = 0;
  double residuo_ = 0.0;
};

}  // namespace t8flow2d

#endif
