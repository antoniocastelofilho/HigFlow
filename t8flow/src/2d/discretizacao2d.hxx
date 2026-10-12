// A ESCOLHA DA DISCRETIZACAO: diferencas finitas ou volumes finitos.
//
// POR QUE ISTO E' UM MODULO, e nao um `enum` a mais no `Poisson`: a escolha entre
// as duas familias nao e' de implementacao, e' de SIGNIFICADO DA INCOGNITA, e
// essa diferenca se propaga por quatro lugares que precisam concordar entre si.
// Espalha-los por `if`s dentro do operador foi o que esta frente encontrou, e e'
// o que esta classe existe para impedir.
//
//     a equacao que se impoe    VF: INTEGRADA sobre a celula
//                               DF: PONTUAL no centro
//
//     o que `u[i]` significa    e' um EIXO SEPARADO, e colapsa-lo no primeiro e'
//                               um erro -- ver abaixo
//
// OS DOIS EIXOS SAO INDEPENDENTES, e foi uma pergunta do usuario que expos isto:
// "tem como discretizar por FV com ordem maior que 2?".  Tem, e a resposta
// mostrou que a tabela acima estava errada quando dizia que volumes finitos
// implica incognita-media.
//
// Em SEGUNDA ordem os dois eixos sao indistinguiveis: a media sobre a celula e o
// valor no centro diferem por `O(h^2)`, que e' a propria ordem do esquema, e
// nenhuma medida os separa.  Acima de segunda ordem eles se separam, e aparecem
// tres combinacoes uteis:
//
//   equacao integrada + incognita MEDIA     o volumes finitos classico de alta
//                                           ordem; exige um passo de
//                                           RECONSTRUCAO media -> pontual, que e'
//                                           onde ele gasta o esforco
//
//   equacao integrada + incognita PONTUAL   o que `vf-mlsgauss` faz: conservativo
//                                           -- o fluxo da faceta e' o mesmo visto
//                                           dos dois lados --, e sem o passo de
//                                           reconstrucao, porque o ajuste MLS ja'
//                                           trabalha com valores pontuais nos
//                                           centros
//
//   equacao pontual + incognita PONTUAL     diferencas finitas
//
// A MEDIDA QUE DECIDE QUAL E' O CASO: `vf-mlsgauss` de ordem 3 e' exato para
// polinomio de grau 3 quando comparado com o VALOR NO CENTRO (erro 1e-14).  Se a
// incognita fosse a media, essa comparacao daria `O(h^2)` -- a diferenca entre a
// media e o valor central de um cubico nao e' zero.  Logo a incognita e' pontual,
// e dizer o contrario no modulo seria descrever um esquema que nao e' o que roda.
//
//     como a fonte entra        VF: INTEGRADA sobre a celula -- e a quadratura
//                                   dessa integral importa, ver a clausula D6:
//                                   a regra do ponto medio travou um esquema
//                                   inteiro em segunda ordem
//                               DF: AVALIADA no centro, sem quadratura
//
//     como o residuo normaliza  VF: pela area, porque a linha foi multiplicada
//                                   por `|face|`
//                               DF: por nada, porque a linha ja' esta' na forma
//                                   do laplaciano pontual
//
// MISTURAR OS QUATRO E' UM ERRO QUE NAO FALHA, E' UM ERRO QUE CONVERGE: um
// operador de volumes finitos com fonte pontual, por exemplo, continua
// convergindo -- em ordem reduzida, e com a constante errada.  Nenhum teste de
// convergencia o pega; o que o pega e' a exatidao polinomial, que e' o oraculo
// que esta suite usa.
//
// A SEGUNDA RAZAO PARA O MODULO e' a proxima etapa.  Com Navier-Stokes entram
// outros operadores -- divergente, gradiente, adveccao --, e cada um tera' a
// mesma pergunta.  Responde-la uma vez, aqui, e' o que impede que o divergente
// seja de volumes finitos enquanto o laplaciano e' de diferencas finitas sem que
// ninguem perceba.

#ifndef T8FLOW_2D_DISCRETIZACAO2D_HXX
#define T8FLOW_2D_DISCRETIZACAO2D_HXX

#include <string>

namespace t8flow2d {

//! A familia da discretizacao -- o que a incognita significa.
enum class Familia {
  //! VOLUMES FINITOS: a incognita e' a MEDIA sobre a celula, e a equacao e'
  //! integrada sobre ela.  E' a familia conservativa: o fluxo que sai de uma
  //! celula entra na vizinha, exatamente, e a soma telescopa.
  //!
  //! E' tambem a unica das duas que trata malha ADAPTADA sem caso especial, e a
  //! razao e' geometrica: o balanco de fluxos nao pressupoe espacamento
  //! uniforme, so' pede a medida de cada faceta.
  VolumesFinitos,

  //! DIFERENCAS FINITAS: a incognita e' o VALOR no centro, e a equacao e'
  //! imposta pontualmente ali.  Nao e' conservativa por construcao -- nada
  //! garante que o que "sai" de um ponto "entre" no vizinho --, e em troca
  //! alcanca ordem alta com estenceis menores.
  DiferencasFinitas,
};

//! A variante dentro da familia.
//
// A lista e' curta de proposito.  Uma variante entra aqui quando esta' MEDIDA --
// grau de exatidao e ordem de convergencia, em malha uniforme e adaptada --, e
// nao quando esta' implementada: variante sem medida e' promessa, e promessa numa
// lista de escolhas e' a forma mais eficiente de alguem escolher errado.
enum class Variante {
  //! VF, segunda ordem: balanco de fluxos de cinco pontos.  Exato ate' grau 2 em
  //! malha uniforme e grau 1 em adaptada.
  FluxoBalanco,

  //! VF, ordem alta: derivada normal por ajuste MLS mais quadratura de Gauss na
  //! faceta e no volume.  Exato ate' o GRAU DA ORDEM do ajuste, em malha
  //! uniforme E adaptada -- a unica variante com essa propriedade.
  MlsGauss,

  //! DF, quarta ordem: compacto de nove pontos na forma de Mehrstellen, com a
  //! correcao `h^2 lap(f)/12` no lado direito.  SO' MALHA UNIFORME.
  CompactoMehrstellen,

  //! DF, quarta ordem: o mesmo esquema na forma de derivada segunda implicita de
  //! Lele.  Verificado equivalente ao de Mehrstellen, nao suposto.  SO' MALHA
  //! UNIFORME.
  CompactoLele,

  //! DF, ordem alta: laplaciano PONTUAL por pesos de ajuste MLS.
  //
  // E' a variante que torna a escolha da familia realmente livre: sem ela,
  // escolher diferencas finitas em malha ADAPTADA nao leva a lugar nenhum,
  // porque as duas variantes compactas recusam malha adaptada por construcao.
  //
  // EXPERIMENTAL, e `estavel_em_malha_adaptada()` diz onde.  Ao contrario do
  // balanco de fluxos, o laplaciano pontual NAO garante estrutura de matriz: o
  // balanco produz diagonal dominante com sinais certos por construcao, e o
  // ajuste MLS nao produz nada em particular.  Medido em malha adaptada com
  // ordem 3: o solver nao converge em 20 mil iteracoes e o erro vai a 1,3e+06.
  // Com ordem 4 converge, em 1,4 mil iteracoes contra 59 da ordem 3 em malha
  // uniforme.
  MlsPontual,
};

//! ONDE UMA GRANDEZA VIVE.  O arranjo deslocado (MAC) e' a decisao do usuario
//! para os DOIS solvers de Navier-Stokes, e ela so' significa algo se cada
//! operador disser de onde le' e para onde escreve.
enum class Local {
  Centro,    //!< escalares e tensores: pressao, viscosidade, tensao
  FacetaX,   //!< a componente `x` da velocidade, na faceta normal a `x`
  FacetaY,   //!< a componente `y` da velocidade, na faceta normal a `y`
};

//! Os operadores que o passo de Navier-Stokes precisa, e as localizacoes de cada.
//
// POR QUE A TABELA VIVE AQUI e nao em cada solver: em MAC a localizacao de
// entrada e a de saida sao DIFERENTES para quase todo operador, e e' justamente
// isso que o arranjo compra -- o gradiente de pressao nasce onde a velocidade
// mora, sem interpolacao, e o divergente da velocidade nasce onde a pressao
// mora.  Essa e' a razao de o MAC nao precisar de estabilizacao: nao ha' par
// pressao-velocidade colocalizado para produzir modo espurio.
//
// Escrever a tabela uma vez, aqui, e' o que impede que o t8SolverNS-FD e o
// t8SolverNS-FV discordem sobre onde uma grandeza vive -- discordancia que nao
// falharia, produziria escoamento.
enum class Operador {
  //! `grad p`: do CENTRO para a FACETA.  Sem interpolacao -- a diferenca entre
  //! as duas celulas vizinhas cai exatamente sobre a faceta entre elas.
  GradientePressao,

  //! `div u`: da FACETA para o CENTRO.  Tambem sem interpolacao, e e' o par
  //! exato do anterior: o divergente discreto e' o transposto (a menos de sinal)
  //! do gradiente discreto, e e' dessa dualidade que vem a conservacao de massa
  //! na projecao.
  DivergenteVelocidade,

  //! `lap u`: da FACETA para a MESMA faceta.  O termo viscoso.
  LaplacianoVelocidade,

  //! `u . grad u`: da FACETA para a MESMA faceta, mas lendo as DUAS componentes.
  //
  // E' O UNICO QUE EXIGE INTERPOLACAO, e por isso e' o unico lugar onde o
  // arranjo deslocado cobra: para avancar a componente `x` e' preciso da
  // componente `y` na faceta de `x`, e elas nao moram juntas.  A interpolacao
  // que se escolhe ali decide a conservacao de quantidade de movimento, e a
  // escolha NAO esta' feita -- ver `doc/projeto-t8flow.md`.
  Adveccao,
};

//! De onde o operador le'.
Local origem_de(Operador op, int direcao);

//! Para onde o operador escreve.
Local destino_de(Operador op, int direcao);

//! O operador precisa interpolar entre localizacoes diferentes da mesma grandeza?
//
// Verdadeiro so' para a adveccao, e saber disso de antemao e' o que permite a um
// solver declarar que nao a implementou ainda em vez de produzir numero.
bool exige_interpolacao(Operador op);

//! A escolha, e o que ela implica.
//
// Objeto de valor: sem estado mutavel, barato de copiar, e cada pergunta tem
// UMA resposta derivada da escolha -- nao um `if` no ponto de uso.
class Discretizacao {
 public:
  //! Constroi a partir da variante; a familia sai dela, nao e' informada.
  //
  // DE PROPOSITO: informar as duas abriria a possibilidade de informa-las
  // INCONSISTENTES, e um par (familia, variante) incoerente e' precisamente o
  // defeito que este modulo existe para tornar impossivel.
  explicit Discretizacao(Variante v, int ordem_mls = 3);

  //! A escolha padrao, e ela e' declarada em vez de implicita: volumes finitos
  //! de segunda ordem, que e' o que funciona em qualquer malha deste sistema.
  static Discretizacao padrao() { return Discretizacao(Variante::FluxoBalanco); }

  Familia familia() const { return familia_; }
  Variante variante() const { return variante_; }

  //! Ordem do ajuste MLS.  So' tem sentido nas variantes que usam MLS; nas
  //! outras vale zero, e perguntar nao e' erro -- a resposta e' "nao se aplica".
  int ordem_mls() const { return ordem_mls_; }

  //! A equacao e' INTEGRADA sobre a celula (VF) ou imposta no ponto (DF)?
  //
  // E' o eixo que DEFINE a familia, e o que decide se o esquema e' conservativo:
  // integrar sobre a celula e aplicar a divergencia transforma o laplaciano em
  // balanco de fluxos, e o fluxo que sai de uma celula entra na vizinha.
  bool equacao_integrada() const { return familia_ == Familia::VolumesFinitos; }

  //! A incognita e' a MEDIA sobre a celula, ou o VALOR no centro?
  //
  // EIXO SEPARADO DO ANTERIOR, e colapsa-los foi um erro deste modulo na
  // primeira versao.  Em segunda ordem eles coincidem -- media e valor central
  // diferem por `O(h^2)`, a propria ordem do esquema --, e so' acima de segunda
  // ordem a diferenca e' observavel.
  //
  // `vf-mlsgauss` integra a equacao e mantem a incognita PONTUAL: medido, ele e'
  // exato para grau 3 comparado com o valor no centro, o que nao aconteceria se
  // a incognita fosse a media.
  bool incognita_e_media() const;

  //! A fonte entra INTEGRADA sobre a celula?  Segue a equacao, nao a incognita.
  bool fonte_integrada() const { return equacao_integrada(); }

  //! Por quanto dividir o residuo do operador para que ele seja comparavel
  //! entre esquemas: a area da celula em VF, e 1 em DF.
  //
  // Em VF a linha foi multiplicada por `|face|`, logo o residuo bruto cresce com
  // o tamanho da celula; em DF a linha ja' esta' na forma do laplaciano pontual.
  // Dividir pela area em DF daria um numero que CRESCE com o refino sem
  // significado -- e esse foi um engano real deste projeto.
  double normalizacao_do_residuo(double area_da_celula) const
  {
    return equacao_integrada() ? area_da_celula : 1.0;
  }

  //! Esta discretizacao trata malha ADAPTADA (transicao 2:1)?
  //
  // As duas variantes compactas NAO tratam, e a resposta aqui e' o que permite
  // recusar cedo e com explicacao, em vez de produzir numeros plausiveis: o
  // estencil de nove pontos pressupoe espacamento uniforme, e aplica-lo numa
  // transicao da' um resultado que converge para a coisa errada.
  bool trata_malha_adaptada() const;

  //! Esta discretizacao produz sistema que o solver RESOLVE em malha adaptada?
  //
  // SEPARADA DE `trata_malha_adaptada` de proposito, porque sao perguntas
  // diferentes e confundi-las esconde a segunda: "trata" e' sobre o OPERADOR
  // estar definido; "estavel" e' sobre o SISTEMA ser resoluvel.
  //
  // O balanco de fluxos garante a estrutura por construcao -- diagonal dominante,
  // sinais certos --, e o laplaciano pontual por MLS nao garante nada: ele monta
  // uma linha que aproxima o laplaciano, e so'.  Medido, ordem 3 em malha
  // adaptada: 20 mil iteracoes sem convergir, erro 1,3e+06.
  //
  // Afirmar isto aqui e' o que impede alguem de escolher a variante para malha
  // adaptada e descobrir sozinho, numa corrida longa, que ela nao fecha.
  bool estavel_em_malha_adaptada() const;

  //! Grau de polinomio ate' o qual esta discretizacao e' EXATA, por tipo de
  //! malha.  Medido, nao deduzido -- ver o relatorio do Poisson.
  //
  // Serve para os testes afirmarem o grau certo sem repetir a tabela, e para o
  // exemplo nao MENTIR na leitura que imprime: dizer "grau 2" para todas as
  // discretizacoes faria o texto errar exatamente onde ele serve para algo.
  int grau_exato(bool malha_adaptada) const;

  //! O ARRANJO e' o deslocado (MAC) nas duas familias -- decisao do usuario em
  //! 12/10/2026, e ela vale para o t8SolverNS-FD e o t8SolverNS-FV.
  //
  // E' uma funcao e nao uma constante porque a pergunta tem de poder ser feita a
  // uma `Discretizacao`: quem recebe uma nao deveria precisar saber que a
  // resposta e' a mesma para todas.  Se um dia uma variante colocalizada entrar,
  // ela entra aqui e nao em cada ponto de uso.
  //
  // MAC EM DIFERENCAS FINITAS e' menos comum que colocalizado com estabilizacao,
  // e foi escolhido pela mesma razao que em volumes finitos: sem par
  // pressao-velocidade colocalizado nao ha' modo espurio de pressao, logo nao ha'
  // o que estabilizar.  O preco e' a interpolacao na adveccao.
  bool arranjo_deslocado() const { return true; }

  //! A ORDEM NOMINAL do esquema -- a ordem de convergencia que ele promete.
  //
  // Nao e' a mesma coisa que `ordem_mls()`, e distingui-las importa: a ordem do
  // AJUSTE e' um detalhe de como a derivada e' obtida, e a ordem do ESQUEMA e' o
  // que se promete a quem escolhe.  Nas variantes de MLS as duas coincidem; no
  // balanco de fluxos a primeira nao existe e a segunda e' 2.
  int ordem() const;

  //! Nome curto, para mensagens e para a linha de comando dos exemplos.
  //
  // INCLUI A ORDEM quando ela e' parte da escolha, e isso nao e' cosmetica:
  // "volumes finitos de quarta ordem" tem de ser uma opcao NOMEAVEL, nao um
  // parametro separado que alguem precise descobrir que existe.  `vf-mlsgauss4`
  // diz a familia, o metodo e a ordem; `vf-mlsgauss` mais um `--mls-ordem=4`
  // escondido em outro lugar diz menos e erra mais.
  std::string nome() const;

  //! A partir do nome; lanca com a LISTA dos validos quando nao reconhece.
  //
  // Aceita com e sem o sufixo de ordem (`vf-mlsgauss4` ou `vf-mlsgauss`); sem o
  // sufixo vale `ordem_mls`, que e' o padrao de quem chama.
  static Discretizacao do_nome(const std::string &nome, int ordem_mls = 3);

  //! Os nomes validos, separados por virgula -- para mensagens de erro e ajuda.
  static std::string nomes_validos();

 private:
  Variante variante_;
  Familia familia_;
  int ordem_mls_;
};

}  // namespace t8flow2d

#endif  // T8FLOW_2D_DISCRETIZACAO2D_HXX
