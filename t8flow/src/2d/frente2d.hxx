// A FRENTE: poligonal fechada de marcadores, com conectividade explicita.
//
// E' o lado "front-tracking" do par que este projeto compara, e a escolha de
// representacao e' O UNICO ponto em que ele difere do VOF -- mesma malha,
// mesmos campos interpolados, mesmo solver de quantidade de movimento, mesma
// projecao, acoplamento no mesmo ponto do passo.  Tudo o mais tem de ser
// compartilhado, senao a diferenca medida entre os dois metodos deixa de ser
// atribuivel a' representacao.
//
// ESTE ARQUIVO E' A FASE B1: advecao CINEMATICA num campo prescrito, sem
// solver.  O acoplamento -- fracao volumetrica, densidade, viscosidade, tensao
// superficial distribuida aos marcadores -- e' a fase B2, e vem depois.  A
// separacao e' deliberada: os oraculos de B1 sao de forma fechada (circulo,
// elipse, reversao de LeVeque) e nao dependem de nada do escoamento, de modo
// que o que eles medem e' o transporte da frente e mais nada.
//
// A CURVATURA SAI DO CIRCUNRAIO de tres marcadores consecutivos,
//
//     kappa = 4 A / (a b c)
//
// com `A` a area do triangulo e `a, b, c` os lados.  A escolha nao e' de
// gosto: para tres pontos SOBRE UM CIRCULO o circunraio E' o raio do circulo,
// qualquer que seja o espacamento entre eles -- logo a curvatura de um circulo
// e' EXATA, e nao aproximada de segunda ordem.  Isso da' um oraculo de
// exatidao onde a diferenca finita parametrica so' daria ordem.
//
// E A AREA SAI DA FORMULA DO CADARCO, que e' exata para um poligono.  A area
// de um CIRCULO nao e' exata -- o poligono inscrito e' menor --, mas o erro e'
// de forma fechada, `R^2 (pi - n sin(2pi/n)/2)`, e e' isso que o arreio cobra.

#ifndef T8FLOW_2D_FRENTE2D_HXX
#define T8FLOW_2D_FRENTE2D_HXX

#include <functional>
#include <vector>

#include "malha2d.hxx"

namespace t8flow2d {

//! O campo que advecta: posicao e instante.
using CampoDeAdvecao = std::function<Ponto2(const Ponto2 &, double)>;

class Frente {
 public:
  //! Um circulo com `n` marcadores, no sentido ANTI-HORARIO.
  //
  // O SENTIDO IMPORTA e nao e' convencao arbitraria: a area do cadarco e'
  // positiva no anti-horario, e a normal "para fora" e' a tangente girada de
  // -90 graus.  Uma frente montada no sentido errado da' area negativa e
  // normal para dentro -- e a tensao superficial empurraria a gota para fora
  // em vez de a conter, o que nao falha, inverte.
  static Frente circulo(const Ponto2 &centro, double raio, int n);

  //! Uma elipse, para o oraculo de curvatura que o circulo nao discrimina.
  static Frente elipse(const Ponto2 &centro, double a, double b, int n);

  explicit Frente(std::vector<Ponto2> marcadores);

  long num_marcadores() const { return static_cast<long>(m_.size()); }
  const std::vector<Ponto2> &marcadores() const { return m_; }
  const Ponto2 &marcador(long k) const;

  //! AREA pelo cadarco, exata para o poligono.
  double area() const;

  //! PERIMETRO: a soma dos lados.
  double perimetro() const;

  //! O CENTROIDE da regiao delimitada, nao a media dos marcadores.
  //
  // A media dos marcadores e' o centroide de uma distribuicao de MASSA nos
  // vertices, e so' coincide com o da regiao quando eles sao equidistantes.
  // Depois da cirurgia eles nao sao.
  Ponto2 centroide() const;

  //! A CURVATURA no marcador `k`, pelo circunraio de `k-1, k, k+1`.
  //! Positiva quando o centro de curvatura esta' do lado de dentro.
  double curvatura(long k) const;

  //! A NORMAL EXTERIOR no marcador `k`, unitaria.
  Ponto2 normal(long k) const;

  //! AVANCA os marcadores no campo, por Euler ou Runge-Kutta 2.
  //
  // RK2 e' o padrao e nao e' luxo: com Euler, um campo de rotacao pura afasta
  // os marcadores do circulo em `O(dt)` POR PASSO, e o erro se acumula como
  // deriva radial -- a gota incha.  Com RK2 a deriva e' de `O(dt^2)` e some no
  // refino do passo.
  void avanca(const CampoDeAdvecao &u, double t, double dt, bool rk2 = true);

  //! MOVE os marcadores por velocidades JA' CALCULADAS, uma por marcador.
  //
  // Separado de `avanca` porque o campo acoplado nao e' funcao do ponto: ele
  // vem interpolado das facetas, e so' existe nos marcadores.  Integrar com
  // RK2 exigiria avaliar o campo numa posicao intermediaria, o que pediria
  // resolver o escoamento de novo -- por isso aqui e' Euler, consistente com
  // o avanco de primeira ordem do solver.
  void move(const std::vector<Ponto2> &vel, double dt);

  //! CIRURGIA LOCAL: parte a aresta que passa de `l_max`, funde a que fica
  //! abaixo de `l_min`.  Devolve quantas partiu e quantas fundiu.
  //
  // SEM ELA A FRENTE MORRE em qualquer escoamento com estiramento: as arestas
  // crescem onde o campo estica e encolhem onde comprime, e em poucas dezenas
  // de passos a resolucao da frente deixa de ter relacao com a da malha.
  //
  // O LIMITE INFERIOR TEM DE SER MENOS DA METADE DO SUPERIOR.  Com
  // `l_min >= l_max/2`, partir uma aresta produz duas que ja' nascem abaixo do
  // minimo, e a proxima cirurgia as funde de volta: a frente oscila entre dois
  // estados a cada passo, sem convergir e sem nada acusar.
  struct Cirurgia { long partidas = 0, fundidas = 0; };
  Cirurgia cirurgia(double l_min, double l_max);

  //! TSUR -- *Trapezoidal Subgrid Undulations Removal*.
  //
  // Remove a ondulacao de comprimento de onda SUBMALHA que a adveccao
  // acumula: marcadores vizinhos oscilando em torno da curva lisa, com
  // amplitude menor que `h`.  Ela nao e' fisica -- e' ruido do transporte --,
  // e alimenta a curvatura, que alimenta a tensao superficial, que alimenta a
  // ondulacao.
  //
  // COMO FUNCIONA, segmento a segmento: o ponto de base e' o MEIO do segmento,
  // de modo que as alturas dos dois extremos acima dele, medidas ao longo da
  // normal do segmento, sao iguais e opostas por construcao -- e `|h1|` mede a
  // inclinacao do segmento em relacao a' normal local, que E' a ondulacao.  Os
  // dois extremos vao para uma altura COMUM, ponderada pelas ESTRELAS (metade
  // da soma dos segmentos incidentes).
  //
  // O LIMIAR NAO E' AFINACAO, e este e' o ponto que nao pode ser perdido: sem
  // ele o filtro move TODO marcador em TODA aplicacao, e reduzir perimetro a
  // area fixa e' EXATAMENTE o que a tensao superficial faz -- o filtro passa a
  // competir com a fisica que o metodo resolve.
  //
  // E CONSERVAR AREA NAO PROTEGE DISSO.  Na circunferencia o efeito e'
  // invisivel, porque o estado liso ja' e' de perimetro minimo; e' a ELIPSE
  // que expoe, e e' nela que o arreio mede.
  //
  // `rlx` e' a relaxacao: 0 move o deslocamento inteiro, 1 nao move nada.
  // `limiar` e' a inclinacao minima que faz o segmento mover, em fracao da
  // MEIA-ARESTA -- o maximo que `|h1|` pode valer.  Devolve quantos moveu.
  long tsur(double rlx, double limiar);

 private:
  std::vector<Ponto2> m_;
};

}  // namespace t8flow2d

#endif  // T8FLOW_2D_FRENTE2D_HXX
