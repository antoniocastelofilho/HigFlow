// NAVIER-STOKES INCOMPRESSIVEL EM TRES DIMENSOES, volumes finitos, arranjo
// deslocado (MAC) sobre a floresta do t8code.
//
// ===========================================================================
// O QUE ESTE ARQUIVO E', HOJE
// ===========================================================================
//
// O SEGUNDO INCREMENTO: divergencia, gradiente e PROJECAO.  Nao ha' ainda termo
// convectivo nem viscoso -- isto nao integra Navier-Stokes.  Ha' o par
// `div`/`grad` e o Poisson que os liga, que e' o esqueleto sobre o qual os
// termos assentam.
//
// O ORACULO DA PROJECAO NAO PRECISA DE FISICA: dado um campo de velocidade
// QUALQUER, projeta-lo tem de deixar divergencia zero de maquina.  Isso
// exercita a cadeia inteira -- divergencia, montagem do Poisson, solucao,
// gradiente -- sem que exista um passo de tempo.
//
// Isto e' deliberado, e repete o caminho do lado 2D.  O risco dominante em
// tres dimensoes nao e' a fisica -- ela e' a mesma -- e' a CONTABILIDADE DE
// FACETAS com tres componentes: seis faces por celula, tres familias de
// facetas, e cada faceta pertencendo a duas celulas.  Escrever o solver
// inteiro e so' entao descobrir que a varredura erra significa procurar o
// defeito atraves de quatro operadores.
//
// O ORACULO E' DE EXATIDAO, nao de ordem.  A divergencia de volumes finitos e'
// a regra do ponto medio sobre a superficie, EXATA para campo linear.  Logo um
// campo linear solenoidal da' zero de maquina em qualquer malha, e um campo de
// divergencia conhecida da' o valor exato.  Ordem certa com constante errada
// existe; exatidao nao.
//
// ===========================================================================
// O QUE E' RECUSADO, E DITO
// ===========================================================================
//
// `np > 1` e TRANSICAO 2:1.  As duas lancam com a razao, e nao produzem numero.
//
// A transicao 2:1 em tres dimensoes NAO e' a de duas com um indice a mais: a
// face da celula grossa e' coberta por QUATRO facetas finas, nao duas, e a
// quina onde a viscosidade se interpola toca OITO celulas, nao quatro.  Tratar
// a face por restricao e a quina por media com as regras de 2D produz
// divergencia que nao fecha -- defeito que aparece como o solver de pressao
// nao convergindo, longe da causa.  Esse modo de falha ja' se manifestou no
// lado 2D e esta' registrado em `doc/relatorio-ns-fv-2d.md`.
//
// ===========================================================================
// POR QUE CODIGO SEPARADO DO 2D, E NAO `#if DIM`
// ===========================================================================
//
// E' a decisao 7 do projeto, e aqui ela se paga: `T8_ECLASS_QUAD` e
// `T8_ECLASS_HEX` sao classes distintas do t8code, e as diferencas do solver
// nao sao de contagem de indices.  A media harmonica da viscosidade passa de
// quatro celulas para oito; o estencil convectivo de segunda ordem ganha uma
// direcao; a faceta deixa de ter comprimento e passa a ter area.  Um corpo
// unico com macro esconderia essas tres como se fossem a mesma coisa.

#ifndef T8FLOW_3D_NAVIERSTOKES3D_HXX
#define T8FLOW_3D_NAVIERSTOKES3D_HXX

#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

#include "campos3d.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"

namespace t8flow3d {

class NavierStokes {
 public:
  //! Constroi sobre a malha e o dominio de facetas, com uma condicao de
  //! contorno por COMPONENTE da velocidade e uma para a pressao.
  //
  // TRES CONTORNOS DE VELOCIDADE, e nao um vetorial: cada componente tem a sua
  // condicao em cada parede, e sao diferentes -- numa parede de deslizamento
  // livre a normal e' Dirichlet e as duas tangenciais sao Neumann.  Um
  // contorno vetorial unico obrigaria quem chama a codificar essa distincao
  // noutro lugar.
  //
  // Lanca se a malha tiver transicao 2:1, ou se rodar em `np > 1`.
  NavierStokes(const Malha3 &m, const DominioFacetas &fd,
               const ContornoPorRegiao &bc_u0, const ContornoPorRegiao &bc_u1,
               const ContornoPorRegiao &bc_u2, const ContornoPorRegiao &bc_p);

  NavierStokes(const NavierStokes &) = delete;
  NavierStokes &operator=(const NavierStokes &) = delete;

  CampoVetorial &velocidade() { return u_; }
  const CampoVetorial &velocidade() const { return u_; }
  CampoEscalar &pressao() { return p_; }
  const CampoEscalar &pressao() const { return p_; }

  //! O passo de tempo, que a projecao usa como escala.
  void passo_de_tempo(double dt);
  double passo_de_tempo() const { return dt_; }

  //! PROJETA O CAMPO ATUAL: resolve `lap(dp) = div(u)/dt` e corrige
  //! `u <- u - dt grad(dp)`.
  //
  // Trata o campo que esta' la' como `u*`.  Nao ha' preditor ainda, e e' por
  // isso que isto e' um metodo proprio e nao um `passo()`: chamar de `passo`
  // o que nao avanca a equacao seria mentir no nome.
  void projeta_campo();

  //! `G q` numa faceta, para um campo de celula qualquer.
  double gradiente_em(const CampoEscalar &campo, int d, long fid) const;

  //! O VOLUME DE CONTROLE DUAL da faceta: em tres dimensoes `area x
  //! distancia entre os centros`.
  //
  // E' o peso em que a adjuncao `D = -G^T` pode valer -- ver a clausula O1, no
  // lado 2D.  Em duas dimensoes era `comprimento x distancia`; a diferenca nao
  // e' de nome.
  double volume_dual_da_faceta(int d, long fid) const;

  int iteracoes_da_pressao() const { return it_pressao_; }
  bool pressao_convergiu() const { return conv_pressao_; }

  //! `D u` numa celula: o fluxo liquido pelas SEIS faces, por unidade de
  //! volume.
  //
  // EM TRES DIMENSOES A MEDIDA DA FACETA E' UMA AREA, e e' por isso que o
  // fluxo e' `u . area` e nao `u . comprimento`.  Trocar isso da' um resultado
  // com dimensao errada que ainda assim converge em malha uniforme -- `h` e
  // `h^2` diferem por um fator constante numa malha so'.
  double divergencia_em(long b, long i) const;

  //! O maximo de `|divergencia_em|` sobre as celulas locais.
  double divergencia_maxima() const;

 private:
  //! A faceta de uma face da celula, VENHA ELA DE QUEM VIER.
  //
  // `DominioFacetas::id_faceta` devolve -1 quando a folha nao e' dona da
  // faceta -- informacao, nao falha.  O balanco de fluxos precisa das SEIS
  // faces mesmo assim, e as que nao sao dela vem do vizinho.
  long faceta_da_face(long b, long i, int d, int lado) const;

  //! `(bloco, folha)` por id global.  So' para achar o vizinho acima.
  std::unordered_map<long, std::pair<long, long>> bi_por_id_;

  //! AS DUAS CELULAS DE CADA FACETA, por direcao: `[0]` do lado menos e `[1]`
  //! do lado mais, cada uma `(bloco, folha)` ou `(-1, -1)` se for contorno.
  //
  // Construido uma vez.  Sem ele, achar as celulas de uma faceta exigiria
  // busca por posicao a cada avaliacao do gradiente, e busca por posicao numa
  // faceta e' justamente onde o ponto cai EM CIMA da fronteira entre duas
  // celulas -- o caso em que a regra semiaberta decide por desempate.
  std::vector<std::array<std::pair<long, long>, 2>> celulas_da_faceta_[3];

  CampoEscalar dp_;
  double dt_ = 0.0;
  int it_pressao_ = 0;
  bool conv_pressao_ = false;
  bool p_tem_dirichlet_ = false;

  void resolve_pressao();
  void projeta();

  const Malha3 *m_ = nullptr;
  const DominioFacetas *fd_ = nullptr;
  const ContornoPorRegiao *bcu_[3] = {nullptr, nullptr, nullptr};
  const ContornoPorRegiao *bcp_ = nullptr;

  CampoVetorial u_;
  CampoEscalar p_;
};

}  // namespace t8flow3d

#endif  // T8FLOW_3D_NAVIERSTOKES3D_HXX
