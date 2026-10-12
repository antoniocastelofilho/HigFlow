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
#include "discretizacao2d.hxx"
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

//! O esquema de discretizacao do laplaciano.
enum class Esquema {
  //! BALANCO DE FLUXOS de segunda ordem, com cinco pontos.  E' o que funciona em
  //! malha ADAPTADA (transicao 2:1) e com Neumann.
  Fluxo2,

  //! COMPACTO DE QUARTA ORDEM (Mehrstellen), nove pontos, so' em malha UNIFORME.
  //
  //     (1/(6h^2)) [ 4(u_E+u_W+u_N+u_S) + (u_NE+u_NW+u_SE+u_SW) - 20 u_C ]
  //         =  f_C + (h^2/12) lap(f)_C
  //
  // O TERMO `h^2/12 lap(f)` E' O QUE FAZ A QUARTA ORDEM.  Sem ele, o estencil de
  // nove pontos e' apenas uma media de dois laplacianos de segunda ordem -- o de
  // cinco pontos e o diagonal -- e continua de segunda ordem.  "Compacto" quer
  // dizer que a quarta ordem sai de um estencil de largura UM, sem precisar de
  // vizinhos a dois passos: o preco e' o lado direito deixar de ser `f` e passar a
  // envolver `lap(f)`.
  //
  // O erro e' proporcional a derivadas de SEXTA ordem, logo o esquema e' EXATO
  // para polinomio de grau <= 5.
  //
  // SO' MALHA UNIFORME DE UM BLOCO, e a montagem recusa o resto: o estencil usa
  // vizinhos DIAGONAIS, e em transicao 2:1 o vizinho diagonal pode nem existir
  // como celula unica.  Nao e' limitacao de implementacao -- e' o esquema que nao
  // esta' definido ali sem uma reconstrucao propria.
  Compacto4,

  //! O MESMO estencil de nove pontos, com o lado direito do esquema COMPACTO DE
  //! LELE -- que e' implicito na segunda derivada:
  //!
  //!     (1/10) u''(i-1) + u''(i) + (1/10) u''(i+1)
  //!          =  (6/5) (u(i+1) - 2u(i) + u(i-1)) / h^2
  //!
  //! Aplicando o operador tridiagonal nas duas direcoes a `lap(u) = f`, o lado
  //! esquerdo da' os mesmos nove pontos em `u` e o lado direito fica `A_x A_y f` --
  //! `f` nos NOVE pontos vizinhos, sem `lap(f)`.
  //!
  //! E OS DOIS SAO O MESMO ESQUEMA, conferido na algebra: `A_x f / (6/5)` expande
  //! em `f + (h^2/12) f_xx + O(h^4)`, porque `(1/10)/(6/5) = 1/12`.  Logo
  //! `A_x A_y f` e' `(6/5)^2 [f + (h^2/12) lap(f)] + O(h^4)`, e os estenceis em `u`
  //! diferem pelo mesmo fator `(6/5)^2 = 1,44` -- medido contra o de Mehrstellen:
  //! `-4,8 / -3,333 = 1,44` e `0,96 / 0,667 = 1,44`.
  //!
  //! A DIFERENCA PRATICA, e e' ela que faz valer ter os dois: Mehrstellen exige
  //! `lap(f)` ANALITICO, que existe na manufatura e nao existe em problema real.
  //! Lele usa apenas `f`, avaliado onde se precisa.  Ter os dois e comparar e' o
  //! que verifica a algebra acima em vez de confiar nela.
  Compacto4Lele,

  //! RECONSTRUCAO MLS COM QUADRATURA DE GAUSS -- o unico que serve em malha
  //! ADAPTADA com ordem acima de dois.
  //
  // O fluxo pela faceta e' `integral du/dn ds`, e este esquema o calcula com as
  // duas aproximacoes levadas a ordem alta ao mesmo tempo -- porque melhorar uma
  // sozinha nao passa de segunda ordem:
  //
  //   a DERIVADA em cada ponto de quadratura vem do ajuste MLS de ordem 3, exato
  //   para polinomio de grau 3 (`Mls::pesos_derivada`);
  //
  //   a QUADRATURA e' de Gauss com DOIS pontos por faceta, exata para integrando
  //   cubico -- contra o ponto medio, exato so' para linear.  Foi a quadratura do
  //   ponto medio que tirou a exatidao do grau 2 na malha adaptada.
  //
  // CONSERVATIVO POR CONSTRUCAO, e isso sai do desenho e nao de um acerto: o
  // suporte do estencil depende apenas do PONTO de quadratura, nao da celula que
  // pergunta.  Logo as duas celulas que compartilham uma faceta obtem o mesmo
  // fluxo, sem ninguem combinar nada.
  //
  // Serve em malha uniforme e adaptada, com Dirichlet e com Neumann -- o
  // fechamento vem do `ConstrutorEstencil`, que ja' trata os tres tipos.
  MlsGauss,

  //! LAPLACIANO PONTUAL por ajuste MLS -- diferencas finitas em malha qualquer.
  //
  // E' a variante que torna a escolha da FAMILIA realmente livre.  Sem ela,
  // escolher diferencas finitas em malha ADAPTADA nao leva a lugar nenhum: as
  // duas variantes compactas recusam malha adaptada por construcao, porque o
  // estencil de nove pontos pressupoe espacamento uniforme.
  //
  // A linha e' o laplaciano avaliado NO CENTRO, sem quadratura nenhuma, e o lado
  // direito e' `f(centro)` -- sem integrar.  Comparar esta variante com a
  // `MlsGauss` e' comparar as duas familias com TUDO O MAIS IGUAL: mesmo ajuste,
  // mesma ordem, mesmo suporte, mesmo fechamento de contorno.  A diferenca
  // medida e' a diferenca entre integrar e avaliar.
  MlsPontual,
};

class Poisson {
 public:
  Poisson(const Malha &m, ContornoPorRegiao contorno);

  //! O termo de fonte `f`.  Na manufatura, `f = lap(u_exata)`, derivado a mao.
  void fonte(std::function<double(const Ponto2 &)> f) { f_ = std::move(f); }

  void fechamento(FechamentoParede fp) { fechamento_ = fp; }
  FechamentoParede fechamento() const { return fechamento_; }

  //! Troca o esquema.  `Compacto4` exige malha uniforme de um bloco e Dirichlet
  //! em todo o contorno -- `resolve` recusa, dizendo qual condicao falhou.
  void esquema(Esquema e) { esquema_ = e; }
  Esquema esquema() const { return esquema_; }

  //! Instala a escolha de DISCRETIZACAO -- familia e variante de uma vez.
  //
  // E' a interface preferida; `esquema()` continua para quem ja' a usa.  A
  // diferenca nao e' de gosto: `Discretizacao` carrega junto as CONSEQUENCIAS da
  // escolha -- o que a incognita significa, como a fonte entra, como o residuo
  // normaliza --, e e' isso que impede que as quatro decisoes se separem.
  void discretizacao(const Discretizacao &d);

  //! A escolha instalada.
  Discretizacao discretizacao() const;

  //! Como `integral f dA` e' calculada no esquema `MlsGauss`.
  //
  // EXISTE PARA DAR CONTRASTE, e o contraste aqui mediu algo que eu nao esperava:
  // com o fluxo da faceta levado a ordem alta e a fonte ainda no ponto medio, o
  // esquema TRAVA em segunda ordem (grau 4 deu 1,99 em vez de exato).  A regra do
  // ponto medio e' um erro `O(h^2)` relativo em `integral f dA`, e no esquema de
  // fluxo de segunda ordem ele se escondia atras do erro do proprio fluxo.
  //
  // Sem este interruptor, "Gauss 2x2 na fonte" seria uma afirmacao sem prova:
  // verde que nunca foi vermelho nao prova nada.
  enum class QuadFonte {
    PontoMedio,   //!< `f(c) * area` -- exata para `f` linear
    Gauss2x2,     //!< exata para `f` bicubico (o padrao)
  };
  void quadratura_fonte(QuadFonte q) { quad_fonte_ = q; }

  //! Ordem do ajuste MLS e teto do suporte, usados pelo esquema `MlsGauss`.
  //! Ordem do ajuste, teto do suporte e quantos ANEIS ele pode abrir.
  //
  // O numero de aneis nao e' so' um teto de busca: ele e' o alcance que a FRANJA
  // tem de cobrir em paralelo, e o `resolve` pede essa cobertura a' malha
  // (`aneis + 1` camadas).  Cada camada custa `(2K+1)^2` folhas publicadas por
  // folha local, logo o custo cresce com o quadrado -- e por isso o numero sai da
  // ORDEM em vez de um padrao generoso.
  //
  // `aneis <= 0` pede o automatico, que e' derivado do minimo de pontos medido:
  // a ordem 3 pede 10 monomios e dois aneis entregam de 12 a 16 pontos; a ordem 4
  // pede 15 e dois aneis entregaram 10 -- insuficiente --, logo ela leva tres.
  void mls(int ordem, int max_pontos, int aneis = 0)
  {
    mls_ordem_ = ordem;
    mls_maxpts_ = max_pontos;
    mls_aneis_ = (aneis > 0) ? aneis : ((ordem >= 4) ? 3 : 2);
  }

  //! A fonte do termo de correcao do esquema compacto: `lap(f)`.
  //
  // Separado de `fonte` de proposito.  Ele e' EXATO na manufatura (derivado a
  // mao), e aproxima-lo por diferencas introduziria no resultado um erro da mesma
  // ordem do que o esquema promete remover.  Sem ele instalado, `Compacto4`
  // RECUSA -- rodar com `lap(f) = 0` daria segunda ordem e pareceria quarta
  // falhando.
  void fonte_laplaciano(std::function<double(const Ponto2 &)> lf) { lap_f_ = std::move(lf); }

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
  //! Os dois tempos de `resolve`, SEPARADOS: montar a matriz e resolve-la.
  //
  // SEPARADOS PORQUE ELES ESCALAM DE FORMAS DIFERENTES, e somados escondem isso.
  // A montagem paraleliza SEM COMUNICACAO NENHUMA -- cada processo monta as suas
  // linhas --, e a resolucao paraleliza mal.  Se a montagem domina em np=1, o
  // tempo total pode PARECER superlinear sem que nada seja superlinear.
  //
  // E' a ultima hipotese de pe' para o ganho de 39x por iteracao medido de np=1
  // para np=8 no esquema MLS: o precondicionador foi refutado (as iteracoes
  // AUMENTARAM 63%) e a hipotese de cache nao se sustentou (a matriz tem 19,6
  // nao-zeros por linha, 3,85 MB no nivel 7, que cabe no L3 de 16 MiB).
  //
  // Em segundos; valem zero antes do primeiro `resolve`.
  double tempo_de_montagem() const { return t_monta_; }
  double tempo_de_resolucao() const { return t_resolve_; }

  //! Quantos ajustes MLS tiveram de REBAIXAR a ordem na ultima montagem.
  //
  // Zero e' o que se quer afirmar: um esquema que diz "ordem 4" com rebaixamento
  // entrega ordem menor em parte do dominio, e sem este numero isso e' invisivel.
  long rebaixamentos_mls() const { return rebaixos_mls_; }
  long ajustes_mls() const { return ajustes_mls_; }

  //! Quantos suportes a FRANJA truncou, somado sobre os processos.
  //
  // Zero em np=1 por construcao.  Diferente de zero em paralelo significa que os
  // pesos perto das fronteiras de particao nao sao os do serial -- a solucao
  // passa a depender do numero de processos, e isso tem de ser visivel.
  long truncados_pela_franja() const { return static_cast<long>(truncados_franja_); }

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
  Esquema esquema_ = Esquema::Fluxo2;
  int mls_ordem_ = 3;
  int mls_maxpts_ = 40;
  int mls_aneis_ = 2;
  QuadFonte quad_fonte_ = QuadFonte::Gauss2x2;
  mutable long rebaixos_mls_ = 0;
  mutable long ajustes_mls_ = 0;
  mutable double truncados_franja_ = 0.0;
  double t_monta_ = 0.0;
  double t_resolve_ = 0.0;
  std::function<double(const Ponto2 &)> lap_f_;

  //! Monta a linha de uma celula.  `indice` e' o mapa id -> (bloco, folha), opaco
  //! aqui para o cabecalho nao depender do tipo interno.
  void monta_linha(long b, long i, const void *indice,
                   std::map<long, double> *linha, double *rhs) const;

 public:
  //! A LINHA MONTADA de uma celula, para oraculo e diagnostico.
  //
  // Existe porque comparar a linha entre np=1 e np>1 e' o unico oraculo que
  // distingue "o esquema montou outra coisa" de "o solver convergiu diferente".
  // Sem ela, um erro de O(1) em paralelo e' indistinguivel de um problema de
  // precondicionamento, e foi exatamente onde este trabalho ficou preso.
  //
  // Devolve os coeficientes por id GLOBAL e o lado direito -- a mesma coisa que
  // vai para o solver, pelo mesmo caminho de codigo.
  std::map<long, double> linha_montada(long b, long i, double *rhs) const;

 private:
  //! Reducoes entre processos.  Identidade em np=1, e sem ramo separado para ele.
  double reduz_max(double v) const;
  double reduz_soma(double v) const;

  //! Monta a linha pelo laplaciano PONTUAL do MLS (diferencas finitas).
  void monta_linha_mls_pontual(long b, long i, const void *construtor,
                               std::map<long, double> *linha, double *rhs) const;

  //! Monta a linha pelo esquema de reconstrucao MLS com quadratura de Gauss.
  void monta_linha_mls(long b, long i, const void *construtor,
                       std::map<long, double> *linha, double *rhs) const;

  //! Monta a linha pelo esquema COMPACTO.  Exige a grade logica (so' uniforme).
  void monta_linha_compacta(long ix, long iy, const void *grade,
                            std::map<long, double> *linha, double *rhs) const;

  //! Valida e constroi a grade logica `(ix, iy) -> id` da malha uniforme.
  //! Lanca se a malha nao for uniforme de um bloco.
  std::vector<long> grade_logica(long *n) const;

  long primeiro_ = 0;
  long n_ = 0;
  std::vector<double> sol_;
  bool convergiu_ = false;
  int iteracoes_ = 0;
  double residuo_ = 0.0;
};

}  // namespace t8flow2d

#endif
