// AS CONDICOES DE FRONTEIRA DO T8FLOW 2D: Dirichlet, Neumann e Robin.
//
// POR CELULAS FANTASMAS, e a escolha e' do usuario -- ele levantou que "Neumann
// com o t8code utilize celulas fantasmas, que e' diferente de como o HiGFlow faz".
// Esta' certo nas duas metades, e vale registrar o que cada um faz:
//
//   O HIGFLOW NAO USA FANTASMA.  Lido em `domain.c` (~2010): ele toma um ponto
//   MEIA CELULA PARA DENTRO do dominio ("Take a point half Δx into the domain"),
//   pega o estencil DAQUELE ponto, e corrige por extrapolacao linear com a
//   derivada ("Approximating the function as a line, get the line slope").  Nao
//   existe celula fora do dominio em lugar nenhum.
//
//   O GHOST DO T8CODE TAMBEM NAO SERVE para isto: ele e' de PARTICAO -- uma camada
//   de face-vizinhos de OUTRO rank --, nao de contorno.  Medido: profundidade 1 e
//   so' `T8_GHOST_FACES`.  Fantasma de contorno e' construcao nossa.
//
// POR QUE FANTASMA CASA BEM COM MLS, e esta e' a razao tecnica de aceitar a
// sugestao: o fantasma entra como um PONTO DE SUPORTE comum, e o valor dele e'
// uma funcao AFIM da incognita da celula espelhada --
//
//     u_fantasma = A * u_celula + B,   com A e B conhecidos da condicao
//
// Logo o peso que o ajuste der ao fantasma se REDISTRIBUI: `A*w` soma no id da
// celula espelhada e `B*w` vai para o `rhs`.  O ajuste nao precisa saber que
// aquele ponto e' especial, e os tres tipos de condicao viram a MESMA conta com
// (A, B) diferentes.  Nao ha' linha de derivada na matriz, nao ha' caso especial
// no laco, e o `Mls` herdado serve sem modificacao.
//
// (A, B) por tipo, com `h` a distancia entre o centro da celula e o do fantasma
// (dois meios-passos: h = lado da celula) e `g` o dado da condicao na parede:
//
//     DIRICHLET   u_parede = g
//                 o fantasma nao e' usado: o ponto de suporte e' a PROPRIA parede,
//                 com valor conhecido -- A = 0, B = g.  E' o que ja' estava
//                 implementado, e continua sendo o caminho mais direto.
//
//     NEUMANN     du/dn = g
//                 u_f = u_c + h*g        ->  A = 1,  B = h*g
//                 (n aponta para FORA; o fantasma esta' do lado de fora)
//
//     ROBIN       alfa*u + beta*du/dn = g   na parede
//                 com u_parede ~ (u_c + u_f)/2  e  du/dn ~ (u_f - u_c)/h:
//                 u_f * (alfa/2 + beta/h) = g - u_c * (alfa/2 - beta/h)
//                 ->  A = (beta/h - alfa/2) / (alfa/2 + beta/h)
//                     B =            g      / (alfa/2 + beta/h)
//                 Neumann e' o caso alfa=0 (A=1, B=h*g/beta) e Dirichlet e' o
//                 limite beta=0 (A=-1, B=2g) -- que e' a reflexao ANTISSIMETRICA,
//                 e por isso Dirichlet continua tendo caminho proprio: refletir
//                 e' mais fragil que amostrar a parede.
//
// LIMITE: a aproximacao de u_parede e da derivada e' de SEGUNDA ordem na faceta
// (media e diferenca central), o que casa com ordem 1 e 2 do ajuste.  Para ordem
// mais alta o fantasma unico deixa de bastar -- precisaria de dois --, e isso esta'
// afirmado no teste pelo CONTRASTE, nao suposto.
#ifndef T8FLOW_CONTORNO2D_HXX
#define T8FLOW_CONTORNO2D_HXX

#include <functional>
#include <stdexcept>
#include <vector>

#include "mls2d.hxx"

namespace t8flow2d {

enum class TipoContorno { Dirichlet, Neumann, Robin };

//! O dado da condicao num ponto da parede: `g` em todos os casos.
using DadoDeContorno = std::function<double(const Ponto2 &)>;

//! Uma condicao de fronteira, com o seu tipo e os seus coeficientes.
//
// `alfa` e `beta` so' valem para Robin.  Para Dirichlet e Neumann sao derivados
// (1,0) e (0,1) -- guardados assim de proposito, para a conta de (A,B) ser UMA
// so' e nao tres: tres contas divergem.
struct Contorno {
  TipoContorno tipo = TipoContorno::Dirichlet;
  DadoDeContorno g;
  double alfa = 1.0;
  double beta = 0.0;

  static Contorno dirichlet(DadoDeContorno g)
  {
    Contorno c; c.tipo = TipoContorno::Dirichlet; c.g = std::move(g);
    c.alfa = 1.0; c.beta = 0.0; return c;
  }
  static Contorno neumann(DadoDeContorno g)
  {
    Contorno c; c.tipo = TipoContorno::Neumann; c.g = std::move(g);
    c.alfa = 0.0; c.beta = 1.0; return c;
  }
  static Contorno robin(double alfa, double beta, DadoDeContorno g)
  {
    if (alfa == 0.0 && beta == 0.0)
      throw std::invalid_argument("Contorno::robin: alfa e beta ambos zero -- nao e' condicao");
    Contorno c; c.tipo = TipoContorno::Robin; c.g = std::move(g);
    c.alfa = alfa; c.beta = beta; return c;
  }

  //! Os coeficientes (A, B) da relacao afim `u_fantasma = A*u_celula + B`.
  //
  // `x_parede` e' onde a condicao e' avaliada; `h` e' a distancia entre o centro
  // da celula e o do fantasma.  Lanca se `h <= 0` ou se o denominador do Robin
  // anular -- `alfa*h/2 + beta == 0` significa condicao degenerada naquela celula,
  // e resolver com ela daria numero plausivel e errado.
  void afim(const Ponto2 &x_parede, double h, double *A, double *B) const;
};

}  // namespace t8flow2d

#endif
