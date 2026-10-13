// NAVIER-STOKES INCOMPRESSIVEL EM VOLUMES FINITOS, bidimensional, newtoniano.
//
// E' o `t8SolverNS-FV` da decisao de arquitetura: o modulo de discretizacao
// escolhe a familia, e este resolve a equacao nela.  O arranjo e' DESLOCADO
// (MAC) -- escalares e tensores nos centros, vetores nas facetas --, que e' o
// que torna o acoplamento pressao-velocidade estavel sem estabilizacao.
//
// AS ESTRATEGIAS SAO AS DO HIGFLOW, por decisao: o objetivo declarado e'
// comparar os dois sistemas no mesmo caso, e um metodo diferente faz qualquer
// diferenca medida ser inconclusiva -- nao se sabe se veio da malha, que e' o
// que se quer medir, ou do metodo, que nao e'.
//
//   SEGREGADO, por PROJECAO.  Preditor sem (ou com) a pressao antiga, Poisson
//   para a correcao, projecao da velocidade.  Nao incremental e incremental;
//   a incremental e' de segunda ordem na pressao e a nao incremental de
//   primeira, e a diferenca entre as duas e' um oraculo barato;
//
//   AVANCO por EULER EXPLICITO e EULER SEMI-IMPLICITO.  Os de ordem mais alta
//   do HiGFlow -- RK2, RK3, Crank-Nicolson, BDF2 -- cabem na mesma estrutura e
//   nao estao aqui: entram quando houver caso que os exercite;
//
//   CONVECTIVO DE SEGUNDA ORDEM.  Central, upwind de primeira (para contraste,
//   nao para producao) e os de segunda do HiGFlow.  O padrao e' CUBISTA, que e'
//   o do caso de referencia.
//
// O QUE ESTE ARQUIVO NAO TRAZ, declarado aqui e nao descoberto depois: malha
// ADAPTADA.  A transicao 2:1 numa faceta de velocidade nao e' o mesmo problema
// que numa faceta de fluxo escalar -- a componente normal e a tangencial se
// restringem de formas diferentes, e errar isso produz divergencia que nao
// fecha em vez de erro visivel.  O construtor RECUSA malha com transicao, em
// vez de rodar e dar numero.
//
// O CASO QUE O VERIFICA E' O CANAL, e a escolha nao e' de conveniencia: o
// escoamento de Poiseuille tem solucao analitica `u = (3/2)(1 - y^2)`, `v = 0`,
// `dp/dx = -2 u_max / Re`, e ela e' QUADRATICA em `y`.  O laplaciano de tres
// pontos e' exato para quadraticas, o termo convectivo se anula identicamente
// (`u` nao depende de `x` e `v` e' nulo) e o gradiente de pressao e' constante.
// Logo o estacionario discreto tem de bater com o analitico em PRECISAO DE
// MAQUINA -- que e' afirmacao muito mais forte que ordem de convergencia, pela
// mesma razao de sempre: a ordem pode estar certa com a constante errada.

#ifndef T8FLOW_2D_NAVIERSTOKES2D_HXX
#define T8FLOW_2D_NAVIERSTOKES2D_HXX

#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "campos2d.hxx"
#include "contorno2d.hxx"
#include "malha2d.hxx"
#include "solverlin.hxx"

namespace t8flow2d {

//! O avanco no tempo.
enum class Tempo {
  EulerExplicito,      //!< `u* = u + dt (-conv + visc - grad p)`; estabilidade por `dt`
  EulerSemiImplicito,  //!< viscoso IMPLICITO; remove a restricao difusiva de `dt`
};

//! O acoplamento pressao-velocidade.
enum class Projecao {
  NaoIncremental,  //!< o preditor ignora a pressao; `p^{n+1} = dp`.  1a ordem em `p`
  Incremental,     //!< o preditor usa `p^n`; `p^{n+1} = p^n + dp`.  2a ordem em `p`
};

//! A discretizacao do termo convectivo.
//
// A ORDEM NAO E' O UNICO EIXO: entre os de segunda ordem, o que muda e' o
// limitador.  Central nao limita e oscila onde o gradiente e' forte; CUBISTA e
// QUICK limitam de formas diferentes.  Upwind de primeira esta' aqui para
// CONTRASTE -- um caso que da' o mesmo resultado com upwind de primeira e com
// CUBISTA nao esta' exercitando o termo convectivo.
enum class Convectivo {
  Central,   //!< 2a ordem, sem limitador
  Upwind1,   //!< 1a ordem.  Para contraste, nao para producao
  Cubista,   //!< 2a ordem limitado.  O padrao, e o do caso de referencia
  Quick,     //!< 2a ordem limitado
};

std::string nome_de(Tempo t);
std::string nome_de(Projecao p);
std::string nome_de(Convectivo c);
Tempo tempo_do_nome(const std::string &s);
Projecao projecao_do_nome(const std::string &s);
Convectivo convectivo_do_nome(const std::string &s);

class NavierStokes {
 public:
  //! A malha, o dominio de facetas e o contorno.
  //
  // O CONTORNO DA VELOCIDADE E' POR COMPONENTE, como o `velocity_0` e
  // `velocity_1` do HiGFlow, e nao um so' para o vetor: na saida do canal a
  // componente normal tem Neumann e a tangencial tambem, mas numa parede movel
  // a normal e' Dirichlet nulo e a tangencial Dirichlet nao nulo.  Um contorno
  // unico para as duas obrigaria a codificar essa diferenca dentro do dado.
  //
  // O DA PRESSAO E' SEPARADO dos dois porque e' mesmo diferente: no canal a
  // entrada tem velocidade de Dirichlet com pressao de Neumann, e a saida o
  // inverso.  Essa troca e' a condicao de contorno do metodo de projecao, nao
  // um detalhe do caso.
  NavierStokes(const Malha &m, const DominioFacetas &fd,
               const ContornoPorRegiao &bc_u0, const ContornoPorRegiao &bc_u1,
               const ContornoPorRegiao &bc_p);

  NavierStokes(const NavierStokes &) = delete;
  NavierStokes &operator=(const NavierStokes &) = delete;

  // ---- parametros ----------------------------------------------------------
  void reynolds(double re);
  void passo_de_tempo(double dt);
  void tempo(Tempo t) { tempo_ = t; }
  void projecao(Projecao p) { proj_ = p; }
  void convectivo(Convectivo c) { conv_ = c; }

  double reynolds() const { return re_; }
  double passo_de_tempo() const { return dt_; }
  Tempo tempo() const { return tempo_; }
  Projecao projecao() const { return proj_; }
  Convectivo convectivo() const { return conv_; }

  // ---- os termos fonte, que sao FUNCOES DO USUARIO -------------------------
  //
  // Como no HiGFlow, e separados em DOIS como la': o de FACETA entra na equacao
  // de quantidade de movimento (forca de corpo -- gravidade, Lorentz, fronteira
  // imersa) e o de CELULA entra na continuidade (fonte de massa).
  //
  // Sao `std::function`, e nao ponteiros de funcao com contexto: quem monta o
  // caso escreve uma lambda que captura o que precisar, e nao ha' um `void*`
  // atravessando a interface.  Nao instalado vale ZERO -- e vale zero de
  // verdade, sem ramo, porque o laco testa o `std::function` uma vez por passo
  // e nao uma vez por faceta.
  //
  // O MESMO VALE PARA AS CONDICOES DE CONTORNO: `Contorno::dirichlet` e as
  // irmas tomam `DadoDeContorno`, que tambem e' `std::function` do ponto da
  // parede, de modo que um perfil de entrada parabolico e um dependente do
  // tempo sao o mesmo mecanismo.

  //! Forca de corpo na faceta: `f(x, d, t)` devolve a componente `d` em `x`
  //! no instante `t`.
  //
  // O INSTANTE ENTRA NA ASSINATURA, e nao e' capturado de fora: o termo fonte
  // de uma solucao manufaturada depende do tempo, e deixar quem chama
  // capturar um `double` por referencia faria o valor correto depender de
  // QUANDO alguem o atualizou -- acoplamento invisivel entre o exemplo e a
  // ordem interna do passo.
  void termo_fonte_faceta(std::function<double(const Ponto2 &, int, double)> f)
  {
    fonte_faceta_ = std::move(f);
  }

  // ---- propriedades variaveis ---------------------------------------------
  //
  // POR OMISSAO NAO HA' CAMPO, e o solver usa `rho = 1` e `mu = 1/Re`.  E' o
  // que o lado monofasico sempre fez, e manter esse caminho EXATO -- e nao
  // "equivalente" -- e' o que permite os casos monofasicos continuarem
  // verificando o codigo depois desta mudanca.
  //
  // COM CAMPO, a equacao passa a ser
  //
  //     rho (du/dt + div(uu)) = -grad p + div(mu (grad u + grad u^T))
  //                             + rho g + f
  //
  // e tres coisas mudam de verdade: o preditor divide por `rho` NA FACETA, o
  // Poisson da projecao ganha coeficiente `1/rho` -- deixa de ser laplaciano e
  // vira operador de coeficiente variavel --, e o termo viscoso deixa de ser
  // `mu lap(u)`.
  //
  // O TENSOR COMPLETO NAO E' LUXO.  Com `mu` constante `div(mu grad u^T)` e'
  // `mu grad(div u)`, que some num campo solenoidal -- mas com `mu` variando
  // ele NAO some, e deixa-lo de fora erra a tensao cisalhante exatamente
  // atraves da interface, que e' onde o problema vive.
  void densidade(const CampoEscalar *rho) { rho_ = rho; }
  void viscosidade(const CampoEscalar *mu) { mu_ = mu; }

  //! A GRAVIDADE, como aceleracao.  Entra como `g`, e nao como `rho g / rho`:
  //! e' a variacao de `rho` que impede a pressao de a equilibrar em toda
  //! parte, e e' dai' que vem o empuxo.
  void gravidade(const Ponto2 &g) { g_ = g; }

  //! FORCA POR FACETA, como CAMPO e nao como funcao do ponto.
  //
  // A tensao superficial do front-tracking nasce assim: ela e' distribuida dos
  // MARCADORES para as facetas por um nucleo, e o resultado e' um campo -- nao
  // ha' `f(x)` que a descreva, porque ela depende de onde os marcadores estao
  // e nao so' de onde a faceta esta'.
  //
  // Embrulhar o campo numa lambda `f(x, d, t)` obrigaria a busca da faceta
  // pelo ponto DENTRO do laco do preditor, e isso e' justamente o quadratico
  // que ja' custou caro aqui.  O ponteiro e' observado, nao possuido: quem o
  // instala garante que o campo vive enquanto o solver o usa.
  void forca_de_faceta(const CampoVetorial *f) { forca_campo_ = f; }

  //! Fonte de massa na celula: `f(x)`.  Entra no lado direito do Poisson da
  //! projecao, e um valor nao nulo significa escoamento NAO solenoidal -- o
  //! que e' legitimo e raro, e por isso vale zero por omissao.
  void termo_fonte_celula(std::function<double(const Ponto2 &)> f)
  {
    fonte_celula_ = std::move(f);
  }

  // ---- estado --------------------------------------------------------------
  CampoVetorial &velocidade() { return u_; }
  const CampoVetorial &velocidade() const { return u_; }
  CampoEscalar &pressao() { return p_; }
  const CampoEscalar &pressao() const { return p_; }

  double instante() const { return t_; }

  //! ONDE FICA O RELOGIO QUE AS CONDICOES DE CONTORNO LEEM.
  //
  // `DadoDeContorno` e' funcao so' do ponto -- e' assim no modulo de contorno,
  // e mudar a assinatura atingiria o Poisson, os campos e tudo o mais.  Uma
  // condicao dependente do tempo captura um `double` por referencia, e o
  // PASSO o move: `t^n` antes do preditor, `t^{n+1}` antes da aplicacao final.
  //
  // QUEM DECLARA O RELOGIO E' O CASO, e nao o solver, por ordem de
  // construcao: a condicao precisa do relogio e o solver precisa da condicao.
  // Com o relogio do lado de fora, a lambda o captura antes de o solver
  // existir.
  //
  // E QUEM O MOVE E' O SOLVER, e nao o caso.  Se o caso o atualizasse antes de
  // chamar `passo`, a condicao valeria no instante errado em metade das
  // aplicacoes -- e isso nao falha: degrada a ordem no tempo, que e' o tipo de
  // erro que se atribui ao esquema.
  void relogio_do_contorno(double *t) { relogio_ = t; }
  long passos_dados() const { return n_passos_; }

  // ---- o passo -------------------------------------------------------------
  //! Um passo de tempo completo.
  void passo();

  //! A DIVERGENCIA MAXIMA do campo de velocidade, reduzida entre processos.
  //
  // E' o oraculo mais barato que existe aqui e o unico que nao precisa de
  // solucao de referencia: o metodo de projecao EXISTE para zerar isto, e um
  // valor que nao cai ao nivel da tolerancia do solver diz que a projecao nao
  // esta' fazendo o que promete -- antes de qualquer comparacao com perfil.
  double divergencia_maxima() const;

  //! A variacao da velocidade no ultimo passo, `max|u^{n+1} - u^n| / dt`.
  //! E' o criterio de estacionariedade, e e' dele que o caso do canal sabe
  //! quando parar.
  double residuo_estacionario() const { return resid_estac_; }

  //! A SOMA DOS FLUXOS no ultimo Poisson -- a condicao de compatibilidade.
  //
  // E' so' o fluxo de `u*` pelo contorno das celulas; os termos de condicao de
  // contorno tambem vao para o lado direito e NAO sao divergencia.
  //
  // Fisicamente e' "o que entra e' o que sai": numa caixa fechada tem de ser
  // nula.  Com Dirichlet na pressao em
  // algum lugar o sistema resolve mesmo incompativel, e a vizinhanca do
  // Dirichlet ABSORVE a incompatibilidade -- o resultado sai plausivel.  Por
  // isso o numero fica a' vista.
  double soma_do_lado_direito() const { return soma_rhs_; }

  int iteracoes_da_pressao() const { return it_pressao_; }
  bool pressao_convergiu() const { return conv_pressao_; }

  //! Quantos termos convectivos TIVERAM DE CAIR para primeira ordem por falta
  //! de ponto a montante, no ultimo passo.
  //
  // Nao e' defeito -- perto da parede o ponto `U` do esquema de segunda ordem
  // cai fora do dominio e nao ha' o que fazer --, mas tem de ser VISIVEL: um
  // esquema que se anuncia de segunda ordem e cai para primeira em parte do
  // dominio entrega menos do que diz, e sem o contador isso e' invisivel.
  long rebaixamentos_convectivos() const { return rebaixos_conv_; }

 private:
  void aplica_contorno_velocidade();
  void preditor();
  void resolve_pressao();
  void projeta();

  //! A faceta da direcao `d` cujo centro e' `x`, ou -1 se nao houver.
  long faceta_em(int d, const Ponto2 &x) const;

  //! UMA SUB-FACETA da face `(d, lado)` de uma celula.
  struct SubFaceta {
    long id;           //!< id da faceta na direcao `d`
    double medida;     //!< comprimento dela (o da face da vizinha MAIS FINA)
    Ponto2 centro;
  };

  //! De que CELULA uma faceta e', e de que LADO.
  //
  // A faceta pertence ao lado mais fino, logo a celula dona e' a mais fina das
  // duas -- e e' a partir dela que todo o estencil se orienta: o passo ao
  // longo da normal e' o lado DELA, nao um `h` global que a malha adaptada nao
  // tem.
  struct Dono { long b, i; int lado; };

  //! `u_d` NO PONTO `x`, por interpolacao dentro da celula que o contem.
  //
  // E' a primitiva de que todo o resto depende, e existe por um fato da
  // transicao 2:1: a face de uma celula FINA cai no MEIO da GROSSA vizinha, e
  // ali nao ha' faceta nenhuma.  Nao e' caso excepcional -- acontece em metade
  // das celulas finas de qualquer interface.
  //
  // Quando ha' faceta no ponto, devolve o valor dela.  Quando nao ha',
  // interpola LINEARMENTE entre as duas faces de direcao `d` da celula que
  // contem o ponto, cada uma tomada como a media de fluxo das suas
  // sub-facetas.  Em malha uniforme o primeiro caso e' o unico, e o segundo
  // nunca roda -- o caminho geral REDUZ ao particular, que e' o que faz os
  // casos uniformes verificarem o codigo adaptado.
  double u_no_ponto(int d, const Ponto2 &x, bool *dentro) const;

  //! A MEDIA DE FLUXO das sub-facetas de uma face: `soma(u . medida) / soma(medida)`.
  //
  // E' a RESTRICAO da componente normal, e a ponderacao nao e' detalhe: e' o
  // que faz a vazao que sai das duas finas ser a que entra na grossa.  Media
  // aritmetica daria o mesmo aqui, onde as duas sub-facetas tem a mesma
  // medida, e daria DIFERENTE se algum dia a razao deixar de ser 2:1 -- e
  // estaria errada sem avisar.
  double media_de_fluxo(long b, long i, int d, int lado) const;

  //! O gradiente TRANSVERSAL da pressao numa celula, por diferenca central
  //! entre os vizinhos em `t` (ou lateral, quando so' ha' um).
  //! OS PESOS do gradiente transversal: a UNICA fonte deles.
  //
  // Aparecem em dois lugares -- coeficientes na matriz do Poisson e avaliacao
  // na projecao -- e ter DUAS implementacoes fez o par `div`/`grad` deixar de
  // fechar na parede, com a divergencia indo de 1e-14 para 9,5e-06.
  void pesos_grad_transversal(long b, long i, int t,
                              std::map<long, double> *w) const;
  double grad_transversal(const CampoEscalar &campo, long b, long i,
                          int t) const;

  //! O vizinho de uma faceta na direcao `e`, lado `s`: VALOR e DISTANCIA.
  //
  // A distancia nao e' `h`: ao longo da normal (`e == d`) e' o lado da celula
  // que se atravessa, e transversalmente e' a meia-soma dos lados das duas.
  // Supor `h` uniforme aqui e' o que quebra o esquema em malha adaptada, e
  // quebra sem erro visivel -- produz um laplaciano plausivel e errado.
  struct Vizinha {
    double valor = 0.0;
    double dist = 0.0;
    bool dentro = false;     //!< falso quando saiu do dominio (contorno)
  };
  Vizinha vizinha_de(int d, long fid, int e, int s) const;

  //! AS SUB-FACETAS de uma face.
  //
  // Em transicao 2:1 a face da celula GROSSA e' coberta por DUAS facetas, cada
  // uma pertencente a uma vizinha fina -- e a grossa nao tem faceta propria
  // ali, porque `DominioFacetas` so' cria a faceta "mais" quando ela e'
  // contorno.
  //
  // E' POR ISSO QUE A BUSCA POR UM PONTO SO' NAO SERVE: o centro da face
  // grossa nao e' o centro de nenhuma das duas sub-facetas, e a busca devolve
  // -1.  O fluxo pela face seria DESCARTADO, e descartado em silencio -- a
  // divergencia daquela celula sairia errada sem nada acusar.
  std::vector<SubFaceta> sub_facetas(long b, long i, int d, int lado) const;

  //! O valor FANTASMA fora do dominio, pela condicao de contorno, com a parede
  //! a MEIO LADO DA CELULA DONA -- e nao a meio `h` global, que a malha
  //! adaptada nao tem.
  double contorno_vizinho(int d, long fid, int e, int s) const;
  double dist_ate_parede(int d, long fid, int e, int s) const;

  //! O lado do VOLUME DE CONTROLE da faceta na direcao `e`: meia celula de
  //! cada lado ao longo da normal, a celula dona inteira na transversal.
  double comprimento_do_volume(int d, long fid, int e) const;

  //! O valor transportado na face do volume, pelo esquema escolhido.  `V` e' a
  //! velocidade que advecta, e o SINAL dela decide quem e' `C`, `D` e `U`.
  double valor_na_face(int d, long fid, int e, int s, double V) const;

  //! `rho` e `mu` num PONTO, interpolados do campo de celulas.  Sem campo
  //! instalado devolvem `1` e `1/Re`, que e' o caminho monofasico EXATO.
  double rho_em(const Ponto2 &x) const;
  double mu_em(const Ponto2 &x) const;

  //! `mu` NA QUINA: media HARMONICA das quatro celulas que a cercam, como
  //! no HiGFlow.  Harmonica porque e' ela que preserva a continuidade da
  //! tensao cisalhante atraves de uma interface plana.
  double mu_na_quina(const Ponto2 &xq) const;

  //! `rho` na FACE entre a celula de centro `c` e a vizinha no lado `s2` da
  //! direcao `d`: media aritmetica, como no HiGFlow.  O Poisson a usa.
  double rho_na_face(const Ponto2 &c, int d, int s2, double hd) const;

  //! `rho` na FACETA.  TEM de dar o mesmo que `rho_na_face` na mesma face --
  //! o preditor e a projecao a usam, e o Poisson usa a outra.
  double rho_na_faceta(int d, long fid) const;

  //! `div(mu grad u) + (grad u^T).grad(mu)` na faceta.  A soma das duas abaixo.
  double tensao_viscosa(int d, long fid) const;

  //! A PARTE DE FLUXO, `div(mu grad u)`: a que o semi-implicito trata
  //! implicitamente.
  double tensao_viscosa_fluxo(int d, long fid) const;

  //! O TRANSPOSTO QUE SOBREVIVE, `(d_j mu)(d_i u_j)`: fica SEMPRE explicito,
  //! como no HiGFlow.  Ele depende de `grad mu` e da derivada CRUZADA; poe-lo
  //! na matriz acoplaria as duas componentes e o sistema deixaria de ser dois
  //! problemas escalares.
  double tensao_viscosa_transposto(int d, long fid) const;

  double conveccao(int d, long fid) const;
  double difusao(int d, long fid) const;
  //! O gradiente de um campo de CELULA numa faceta.  Recebe o campo porque
  //! o preditor o quer de `p` e a projecao de `dp` -- e tem de ser o MESMO
  //! operador: e' o par `div`/`grad` que fecha, nao cada um.
  double grad_campo(const CampoEscalar &campo, int d, long fid) const;
  double velocidade_que_advecta(int d, long fid, int e, int s) const;

  const Malha *m_ = nullptr;
  const DominioFacetas *fd_ = nullptr;
  const ContornoPorRegiao *bcu_[2] = {nullptr, nullptr};
  const ContornoPorRegiao *bcp_ = nullptr;

  CampoVetorial u_, ustar_;
  CampoEscalar p_, dp_;

  double re_ = 1.0;
  double dt_ = 0.0;
  double t_ = 0.0;
  double t_bc_ = 0.0;
  double *relogio_ = nullptr;
  bool p_tem_dirichlet_ = false;
  mutable double soma_rhs_ = 0.0;
  long n_passos_ = 0;
  double resid_estac_ = 0.0;
  int it_pressao_ = 0;
  bool conv_pressao_ = false;

  //! O lado da folha MAIS FINA.  A chave do indice de facetas quantiza por
  //! `h_min/2`, e nao pelo lado "da" folha, porque numa malha adaptada nao ha'
  //! um.  Como os niveis sao potencias de dois e os blocos se dividem
  //! diadicamente, TODA posicao de faceta de QUALQUER nivel cai exatamente
  //! nessa grade -- a quantizacao continua exata, e nao aproximada.
  Ponto2 h_min_{0.0, 0.0};
  Ponto2 h_{0.0, 0.0};                 //!< lado da folha quando a malha e' uniforme
  bool uniforme_ = true;
  //! Indice de FACETA POR PONTO, por direcao.  Nao existe em `DominioFacetas`
  //! porque o consumidor dele e' o esquema, nao a numeracao: quem monta a
  //! equacao pergunta "qual e' a faceta meio passo adiante", e quem numera nao
  //! tem por que saber disso.  A chave quantiza em `h/2`, que e' exato para as
  //! posicoes de faceta de uma malha uniforme.
  std::unordered_map<unsigned long long, long> fac_por_ponto_[2];

  //! Dono de cada faceta, por direcao.  Montado uma vez.
  std::vector<Dono> dono_[2];

  //! `(bloco, folha)` por id global.  `id_do_ponto` devolve id, e todo o
  //! estencil precisa de `(b, i)` para perguntar lado, nivel e sub-facetas.
  std::unordered_map<long, std::pair<long, long>> bi_por_id_;
  mutable long rebaixos_conv_ = 0;

  Tempo tempo_ = Tempo::EulerExplicito;
  Projecao proj_ = Projecao::Incremental;
  Convectivo conv_ = Convectivo::Cubista;
  std::function<double(const Ponto2 &, int, double)> fonte_faceta_;
  const CampoVetorial *forca_campo_ = nullptr;
  const CampoEscalar *rho_ = nullptr;
  const CampoEscalar *mu_ = nullptr;
  Ponto2 g_{0.0, 0.0};
  std::function<double(const Ponto2 &)> fonte_celula_;
};

}  // namespace t8flow2d

#endif  // T8FLOW_2D_NAVIERSTOKES2D_HXX
