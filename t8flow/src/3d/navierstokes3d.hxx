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
// face da celula grossa e' coberta por QUATRO facetas finas, nao duas, e cada
// faceta tem DUAS direcoes transversais, nao uma.  Tratar a face por restricao
// com as regras de 2D produz divergencia que nao fecha -- defeito que aparece
// como o solver de pressao nao convergindo, longe da causa.  Esse modo de
// falha ja' se manifestou no lado 2D e esta' registrado em
// `doc/relatorio-ns-fv-2d.md`.
//
// NAO SAO OITO CELULAS NA MEDIA DA VISCOSIDADE, e este comentario ja' disse
// que eram.  O ponto onde `mu` e' avaliado para a direcao transversal fica
// sobre uma ARESTA -- dois planos de face cruzando, com a terceira coordenada
// num centro de celula --, e aresta e' tocada por QUATRO celulas, como a quina
// em duas dimensoes.  Oito seria um VERTICE, tres planos de face, que este
// termo nao usa.  A afirmacao errada sobreviveu a uma clausula de contrato e a
// dois commits antes de alguem fazer a conta.
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
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "campos3d.hxx"
#include "contorno3d.hxx"
#include "malha3d.hxx"

namespace t8flow3d {

//! Os esquemas do termo convectivo.  Os mesmos do lado 2D e do HiGFlow.
enum class Convectivo { Central, Upwind1, Cubista, Quick };

//! O nome, para configuracao por texto.  Lanca com a lista dos validos.
Convectivo convectivo_do_nome(const std::string &s);

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

  //! O ESQUEMA CONVECTIVO.
  void convectivo(Convectivo c) { conv_ = c; }
  Convectivo convectivo() const { return conv_; }

  //! `sum_e d(u_e u_d)/dx_e` na faceta, em FORMA DIVERGENTE.
  double conveccao(int d, long fid) const;

  //! Quantos termos convectivos TIVERAM DE CAIR para primeira ordem por falta
  //! de ponto a montante, desde a ultima vez que isto foi zerado.
  //
  // IMPORTA MAIS EM TRES DIMENSOES do que em duas: sao tres direcoes por
  // faceta em vez de duas, logo mais oportunidades de o ponto a montante cair
  // fora do dominio.  Um esquema que se anuncia de segunda ordem e cai para
  // primeira em parte do dominio entrega menos do que diz, e sem o contador
  // isso e' invisivel.
  long rebaixamentos_convectivos() const { return rebaixos_conv_; }
  void zera_rebaixamentos() { rebaixos_conv_ = 0; }

  //! O numero de Reynolds, quando NAO ha' campo de viscosidade.
  void reynolds(double re);

  //! O CAMPO DE VISCOSIDADE, nos centros.  Sem ele, `mu = 1/Re` uniforme.
  //
  // Ponteiro e nao copia: quem o atualiza e' o acoplamento bifasico, a cada
  // passo, e copiar obrigaria a reinstalar.
  void viscosidade(const CampoEscalar *mu) { mu_ = mu; }

  //! `div(mu grad u) + (grad u^T).grad(mu)` na faceta -- o TENSOR COMPLETO.
  //
  // E' a forma do HiGFlow, e isto custou ser aprendido no lado 2D: ele omite
  // so' `mu grad(div u)`, que some num campo solenoidal, e NAO o transposto
  // inteiro.  Descartar o transposto joga fora a tensao cisalhante atraves da
  // interface, onde o caso bifasico vive.
  double tensao_viscosa(int d, long fid) const;

  //! A parte de FLUXO, `div(mu grad u)`.
  double tensao_viscosa_fluxo(int d, long fid) const;

  //! O TRANSPOSTO QUE SOBREVIVE, `(d_e mu)(d_d u_e)`.
  double tensao_viscosa_transposto(int d, long fid) const;

  //! `mu` num ponto, pelo valor da CELULA que o contem.
  double mu_em(const Ponto3 &x) const;

  //! `mu` sobre uma ARESTA, media HARMONICA das QUATRO celulas que a tocam.
  //
  // QUATRO, e nao oito.  O ponto fica sobre dois planos de face com a terceira
  // coordenada num centro de celula -- isso e' uma aresta, e aresta e' tocada
  // por quatro.  Oito seria um VERTICE, tres planos de face, que este termo
  // nao usa.
  //
  // HARMONICA, e nao aritmetica: e' ela que preserva a continuidade da tensao
  // cisalhante atraves de uma interface plana.  Com a aritmetica a tensao
  // salta na interface, e isso nao falha -- muda a taxa de deformacao
  // exatamente onde os dois fluidos se tocam.
  double mu_na_aresta(const Ponto3 &xa, int d, int e) const;

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

  //! A faceta de direcao `d` cujo centro e' `x`, ou -1.
  long faceta_em(int d, const Ponto3 &x) const;

  //! A componente `d` da velocidade num ponto, por interpolacao linear ao
  //! longo de `d`.  Em malha uniforme o ponto cai num centro de celula.
  double u_no_ponto(int d, const Ponto3 &x, bool *dentro) const;

  //! Valor e distancia da faceta vizinha, andando em `e` no sentido `s`.
  struct Vizinha {
    double valor = 0.0;
    double dist = 0.0;
    bool dentro = false;
  };
  Vizinha vizinha_de(int d, long fid, int e, int s) const;

  //! O fechamento de contorno: o fantasma, AFIM na propria incognita.
  double contorno_vizinho(int d, long fid, int e, int s) const;

  //! A velocidade que ATRAVESSA a face do volume de controle.
  double velocidade_que_advecta(int d, long fid, int e, int s) const;

  //! O valor da componente `d` interpolado na face, pelo esquema.
  double valor_na_face(int d, long fid, int e, int s, double V) const;

  //! O comprimento do volume de controle da faceta na direcao `e`.
  double comprimento_do_volume(int d, long fid, int e) const;

  Convectivo conv_ = Convectivo::Cubista;
  mutable long rebaixos_conv_ = 0;
  const CampoEscalar *mu_ = nullptr;
  double re_ = 1.0;

  //! Centro quantizado -> id, por direcao.  So' funciona em malha uniforme,
  //! que e' a unica que este solver aceita hoje.
  std::unordered_map<long long, long> faceta_por_ponto_[3];
  double quantum_ = 1.0;
  long long chave(const Ponto3 &x) const;

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
