// A MALHA MULTIBLOCOS 2D DO T8FLOW.
//
// Multiblocos e' a razao pela qual o t8code casa com o desenho do HiGFlow: o
// `cmesh` e' nativamente MULTIARVORE, e uma arvore por bloco e' o mesmo conceito,
// nao uma emulacao.  Cada bloco tem a SUA caixa e o SEU nivel de refino -- e' o
// multiblocos de verdade, resolucao por regiao.
//
// A FLORESTA E' A ESTRUTURA (decisao 2), nao um produtor de malha para outra
// estrutura.  Esta classe e' dona dela.
//
// AS INTERFACES ENTRE BLOCOS SAO UNIDAS (`t8_cmesh_set_join`), e isso NAO e'
// detalhe: sem o join, cada bloco e' uma ilha e toda face na borda de um bloco
// aparece como "face sem vizinho" -- que e' o criterio de CONTORNO.  Dois blocos
// lado a lado do mesmo dominio dariam contorno fisico no meio do escoamento, em
// silencio.  E e' a decisao D3 que a HiGTree ja' tomou: interface entre blocos do
// MESMO dominio nao e' contorno, porque ha' malha dos dois lados.
//
// LIMITE DECLARADO, e ele RECUSA em vez de aproximar: so' faces EXATAMENTE
// coincidentes sao unidas -- mesmo plano, mesmo intervalo transversal.  Faces
// que se sobrepoem PARCIALMENTE (um bloco encostando no meio da face de outro)
// nao tem join simples no t8code, e aceita-las calada produziria contorno onde
// nao ha'.  O construtor lanca, nomeando o par.
//
// CODIGO SEPARADO POR DIMENSAO (decisao 7): aqui o elemento e' `T8_ECLASS_QUAD`,
// 4 vertices e 4 faces, escrito explicitamente.  A ordem dos vertices e' a do
// t8code -- z-order, x variando mais rapido: (lo,lo), (hi,lo), (lo,hi), (hi,hi)
// --, e as faces do QUAD sao 0:x-, 1:x+, 2:y-, 3:y+.
#ifndef T8FLOW_MALHA2D_HXX
#define T8FLOW_MALHA2D_HXX

#include <array>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

#include "mls2d.hxx"   // Ponto2

namespace t8flow2d {

//! Um bloco: caixa alinhada aos eixos, com o seu proprio nivel de refino.
struct Bloco {
  std::array<double, 2> lo;
  std::array<double, 2> hi;
  int nivel;   //!< refino uniforme DENTRO do bloco: 4^nivel folhas
};

//! A malha: uma floresta do t8code com uma arvore por bloco.
class Malha {
 public:
  //! Monta o cmesh multiarvore, une as interfaces coincidentes e refina cada
  //! bloco ate' o nivel pedido.
  //
  // Lanca `std::invalid_argument` se a lista estiver vazia, se algum bloco for
  // degenerado (`hi <= lo` em qualquer direcao), se algum nivel for negativo, se
  // dois blocos se SOBREPUSEREM em area, ou se duas faces se sobrepuserem
  // parcialmente sem coincidir.
  Malha(const std::vector<Bloco> &blocos, sc_MPI_Comm comm);
  ~Malha();

  Malha(const Malha &) = delete;
  Malha &operator=(const Malha &) = delete;

  //! Numero de blocos pedidos -- uma arvore por bloco, em serie.
  long num_blocos() const { return static_cast<long>(blocos_.size()); }
  const Bloco &bloco(long b) const;

  //! Interfaces unidas entre blocos.  Zero e' legitimo (blocos disjuntos) e e'
  //! exatamente o que um teste tem de poder distinguir de "o join nao rodou".
  long num_juncoes() const { return juncoes_; }

  //! Folhas locais, somadas sobre os blocos.
  long num_folhas() const;

  //! Folhas da arvore do bloco `b`.
  long num_folhas_do_bloco(long b) const;

  //! Centroide da folha `i` da arvore do bloco `b`.
  Ponto2 centroide(long b, long i) const;

  //! Centroides de TODAS as folhas locais, na ordem de percurso da floresta --
  //! que e' a da curva de preenchimento, a propriedade que da' numeracao
  //! contigua por processo.
  std::vector<Ponto2> centroides() const;

  //! Quantos vizinhos a folha `i` do bloco `b` tem pela face `fa` (0:x-, 1:x+,
  //! 2:y-, 3:y+).  Zero significa face sem vizinho, que e' o criterio de
  //! contorno do dominio.
  int num_vizinhos(long b, long i, int fa) const;

  //! Nivel de refino da folha.
  int nivel(long b, long i) const;

  //! Area da folha.  Existe para o oraculo de CONSERVACAO: a soma das areas das
  //! folhas de um bloco tem de ser a area da caixa dele, antes e depois de
  //! qualquer refino.  Contagem de folhas nao pega celula perdida ou duplicada
  //! com geometria errada; area pega.
  double area(long b, long i) const;

  //! Ids GLOBAIS dos vizinhos da folha pela face `fa`.  Vazio quando a face e'
  //! contorno do dominio; dois ou mais quando o lado de la' e' mais fino (2:1).
  //
  // O id sai do indice de FLORESTA que o t8code devolve, convertido para (arvore,
  // folha) -- e nao do indice tratado como se fosse de arvore.  Essa confusao
  // entrega um elemento EXISTENTE, so' que outro, e o erro e' silencioso.
  std::vector<long> ids_vizinhos(long b, long i, int fa) const;

  //! Niveis dos vizinhos da folha pela face `fa`.
  //
  // E' o insumo do oraculo 2:1 -- e o oraculo fica no TESTE, nao aqui: codigo que
  // verifica a si mesmo nao discrimina.
  std::vector<int> niveis_vizinhos(long b, long i, int fa) const;

  //! REFINO LOCAL: refina ate' `nivel_alvo` toda folha cujo CENTRO caia na caixa
  //! [lo, hi], e reequilibra a floresta.
  //
  // Depois desta chamada a razao entre niveis de celulas vizinhas e' no maximo
  // 2:1 -- nao porque o criterio seja suave, mas porque o `t8_forest_set_balance`
  // insere o refino de transicao.  SEM o balance o salto e' o que o criterio
  // pedir: um refino de dois niveis numa caixa da' 4:1 na borda dela, e isso
  // atravessa o estencil do discretizador sem nenhum aviso.
  //
  // A ORDEM DOS COMMITS NAO E' LIVRE: adapt num commit, balance em OUTRO.  Pedir
  // `set_adapt`, `set_ghost` e `set_partition` no mesmo commit do t8code derruba
  // com np>1 e FUNCIONA em np=1 -- armadilha medida, e das piores, porque o teste
  // serial nao a ve.
  //
  // Lanca `std::invalid_argument` se a caixa for degenerada, se nao intersectar
  // nenhum bloco, ou se `nivel_alvo` nao for maior que o nivel de algum bloco que
  // a caixa alcanca -- pedir refino que nao refina nada e' quase sempre engano de
  // quem chamou, e aceitar calado deixaria a malha parecer refinada.
  void refina_caixa(const Ponto2 &lo, const Ponto2 &hi, int nivel_alvo);

  //! REFINO ALEATORIO: refina uma fracao das folhas, escolhidas por sorteio, e
  //! reequilibra em 2:1.
  //
  // E' a malha que o teste de convergencia em quadtree precisa: refinamento sem
  // padrao, com transicoes 2:1 espalhadas, que e' o caso em que um esquema errado
  // na face de transicao aparece -- numa malha refinada por caixa, as transicoes
  // ficam todas alinhadas aos eixos e um erro que dependa da orientacao pode
  // passar.
  //
  // REPRODUTIVEL, E ISSO E' O PONTO DELICADO.  O sorteio NAO vem de um gerador
  // sequencial: vem de um HASH DA POSICAO da celula, misturado com a semente.  A
  // diferenca importa por duas razoes, e as duas sao de erro silencioso:
  //
  //   * com gerador sequencial, a decisao depende da ORDEM EM QUE AS CELULAS SAO
  //     VISITADAS -- e em paralelo, de quantas celulas cada rank tem.  A mesma
  //     semente daria malhas diferentes em np=1 e np=4, e comparar resultados
  //     entre eles deixaria de ser possivel;
  //   * com hash da posicao, a mesma celula recebe a MESMA decisao sempre, em
  //     qualquer ordem e em qualquer particao.  A malha passa a ser funcao de
  //     (blocos, semente, fracao), e so'.
  //
  // `nivel_max` limita o refino: sem teto, chamadas repetidas levam a floresta ao
  // limite do t8code, e o sintoma e' memoria, nao erro.
  //
  // Lanca se `fracao` nao estiver em (0, 1] ou se `nivel_max` nao for maior que o
  // nivel atual -- refino que nao refina nada e' quase sempre engano de quem
  // chamou, e aceito calado deixaria a malha PARECER adaptada.
  void refina_aleatorio(double fracao, unsigned long semente, int nivel_max);

  //! Quantas folhas tem vizinho de nivel DIFERENTE -- o tamanho da transicao 2:1
  //! na malha.
  //
  // Existe para o teste poder exigir que o caso EXERCITE o que afirma: uma malha
  // "aleatoria" que por acaso nao tenha transicao nenhuma passaria em qualquer
  // teste de esquema na face de transicao sem testar nada.
  long num_folhas_em_transicao() const;

  //! Maior e menor nivel presentes na malha.
  int nivel_minimo() const;
  int nivel_maximo() const;

  // -------------------------------------------------------------------------
  // A NUMERACAO GLOBAL -- a unica coisa, alem do suporte, que o nucleo herdado
  // precisa da malha.
  //
  // O solver pede `slv_create(engine, primeiro_id, n_local)`: um id global por
  // incognita, CONTIGUO por processo.  E' exatamente o que a curva de
  // preenchimento do t8code da' naturalmente -- as folhas de um rank ocupam um
  // intervalo [primeiro, primeiro + n_local) na ordem da curva --, e e' por isso
  // que esta parte do t8Flow e' leitura e nao construcao: o t8code ja' numerou.
  //
  // O QUE NAO PODE SER O ORACULO DISTO, e esta' medido no HiGFlow: comparar a
  // numeracao com o percurso que a gerou e' TAUTOLOGIA -- trocar `mp_lookup` por
  // `n_visitadas++` deixava a suite inteira verde, nos dois DIM e nos tres np.  A
  // numeracao tem de ser alcancada por outra via, e e' para isso que existe
  // `id_do_ponto`.
  // -------------------------------------------------------------------------

  //! Id global da PRIMEIRA folha local.  Coletiva (o t8code precisa das contagens
  //! de todos os ranks para montar o deslocamento).
  long primeiro_id() const;

  //! Total de folhas na floresta, somando todos os ranks.
  long num_folhas_globais() const;

  //! Id global da folha `i` do bloco `b`.
  //
  // E' `primeiro_id()` mais a posicao da folha na ordem de percurso da floresta
  // (as arvores em ordem, as folhas de cada arvore na ordem da curva).
  long id_global(long b, long i) const;

  //! Lados da folha em cada direcao.  Nao e' `sqrt(area)`: bloco nao-quadrado tem
  //! celula nao-quadrada -- o bloco [1,3]x[0,1] em nivel 1 tem celulas 1,0 x 0,5.
  Ponto2 lados_da_folha(long b, long i) const;

  //! Id global da folha que CONTEM o ponto, ou -1 se nenhum bloco o contem.
  //
  // E' a VIA INDEPENDENTE da numeracao: acha a celula pela geometria, nao por
  // contagem, e e' o que permite verificar que o id usado para escrever uma
  // incognita e' o mesmo que a geometria daquela posicao devolve.  Sem isso, a
  // numeracao so' poderia ser comparada consigo mesma.
  //
  // DESEMPATE EXPLICITO, porque a geometria sozinha nao o da': o intervalo de cada
  // celula e' SEMIABERTO, [lo, hi) em cada direcao.  Logo um ponto exatamente
  // sobre a face comum pertence a' celula de coordenada MAIOR, sempre, e a
  // resposta nao depende da ordem em que as celulas sao visitadas -- que e' como
  // a C8 da HiGTree passava por sorte.  O limite superior do dominio e' a unica
  // excecao: la' o semiaberto deixaria o ponto fora, e a ultima celula o inclui.
  //
  // Custo O(folhas): varredura. A busca em arvore (`t8_forest_search`) e' trabalho
  // futuro, e esta' dito aqui para ninguem a supor.
  long id_do_ponto(const Ponto2 &p) const;

  // =========================================================================
  // A FRANJA -- o que o PARALELO exige e o serial nunca exercita
  // =========================================================================
  //
  // Em np=1 nada aqui e' usado, e e' exatamente por isso que precisa existir
  // ANTES de rodar em paralelo: o caminho serial nao distingue um vizinho local
  // de um vizinho que vive em outro processo, e a conta errada da' um elemento
  // EXISTENTE -- so' que outro.  Erro silencioso, e e' o modo de falha que esta
  // frente herdou documentado do sistema antecessor.
  //
  // O INDICE DE FLORESTA e' a chave de tudo aqui.  O t8code numera as folhas
  // locais de 0 a `num_folhas()-1` e as da franja logo em seguida, de
  // `num_folhas()` a `num_folhas()+num_fantasmas()-1`.  As consultas de
  // vizinhanca devolvem indices nessa faixa unica, e tratar um indice de franja
  // como local -- somando o deslocamento do processo, por exemplo -- e' o erro
  // que esta API existe para tornar impossivel.

  //! Quantas folhas de FRANJA este processo recebeu.  Zero em np=1.
  long num_fantasmas() const { return n_fantasmas_; }

  //! Id GLOBAL de uma folha pelo indice de floresta, local ou de franja.
  //
  // Para as locais e' `primeiro_id() + indice`.  Para as de franja NAO e': o id
  // delas pertence a outro processo, e vem da TROCA de dados da franja, nao de
  // aritmetica.  Esta funcao e' o unico lugar onde essa distincao mora.
  long id_global_do_indice(long idx) const;

  //! Geometria de uma folha pelo indice de floresta, local ou de franja.
  //
  // Os valores da franja vem da troca, nao de recalculo: recalcular exigiria a
  // geometria da arvore do dono, que este processo pode nao ter.
  Ponto2 centroide_do_indice(long idx) const;
  Ponto2 lados_do_indice(long idx) const;
  int nivel_do_indice(long idx) const;

  //! Vizinhos de face pelo INDICE de floresta -- local ou de franja.
  //
  // Devolve vazio em duas situacoes DIFERENTES, e quem pergunta precisa das
  // duas: a face e' contorno do DOMINIO, ou a folha esta' na ULTIMA camada da
  // franja e este processo nao sabe com que ela faz fronteira.  Use
  // `expansivel_do_indice` para distinguir -- confundir as duas faria uma folha
  // da borda da franja parecer borda do dominio, e o esquema aplicaria condicao
  // de contorno no MEIO do dominio.
  std::vector<long> vizinhos_do_indice(long idx, int fa) const;

  //! Niveis dos vizinhos de face pelo INDICE -- viajam junto com os ids.
  //
  // Existe para que o nivel de um vizinho nao exija que ELE esteja na tabela com
  // geometria, o que custaria uma camada inteira de franja.  Medido: sem isso,
  // nem quatro camadas bastavam, porque a carga satura antes.
  std::vector<int> niveis_vizinhos_do_indice(long idx, int fa) const;

  //! Quantas vezes a carga da troca SATUROU na ultima montagem da franja.
  //
  // Diferente de zero significa que alguma folha publicou menos vizinhanca do
  // que conhecia, e portanto que algum processo sabe menos do que deveria.  Nao
  // e' erro por si -- o canal e' de tamanho fixo --, mas e' a explicacao de
  // qualquer ausencia inesperada, e por isso nao pode ficar invisivel.
  long saturacoes_da_franja() const { return saturacoes_; }

  //! Esta folha pode abrir mais um anel, ou esta' na borda da franja?
  bool expansivel_do_indice(long idx) const;

  //! Quantas CAMADAS de franja esta malha mantem.
  //
  // A franja do t8code tem UMA camada, e uma camada nao basta: medido, um ponto
  // de quadratura sobre uma faceta de fronteira de particao cai numa celula de
  // franja, e o suporte que cresce dali tem UM ponto -- contra os 10 que a
  // ordem 3 do MLS exige.  A alternativa de crescer o suporte a partir da
  // celula local adjacente foi rejeitada: faria o suporte depender de QUEM
  // PERGUNTA, e os dois fluxos pela mesma faceta deixariam de se cancelar,
  // perdendo a conservacao que a clausula D7 afirma em `0,0e+00`.
  //
  // As camadas adicionais vem de RODADAS DE TROCA: cada processo publica, por
  // folha sua, a vizinhanca que ja' conhece, e cada rodada estende o alcance em
  // um anel.  O custo e' `(2K-1)^2` folhas publicadas por folha local, logo
  // cresce com o quadrado de K -- e e' por isso que K e' escolha de quem monta
  // o esquema, nao um padrao generoso.
  int camadas_de_franja() const { return camadas_; }

  //! Garante ao menos `k` camadas de franja, refazendo a troca se preciso.
  //
  // Idempotente: pedir menos do que ja' existe nao desfaz nada.  Em np=1 e'
  // no-op, porque nao ha' franja -- e nao ha' ramo separado para o serial.
  void exige_camadas(int k);

  //! Indice de floresta de um id GLOBAL -- ou `-1` se este processo nao o tem.
  //
  // Cobre as locais e as de franja.  O `-1` e' informacao, nao falha: um id que
  // nao esta' nem nas locais nem na franja esta' longe demais para este processo
  // montar estencil com ele, e quem pergunta tem de decidir o que fazer.
  long indice_do_id(long id) const;

  //! Quantos processos, e qual e' este.  Para as reducoes de quem consome.
  int num_processos() const { return nprocs_; }
  int processo() const { return rank_; }
  sc_MPI_Comm comunicador() const { return comm_; }

  //! REPARTICIONA pela curva de preenchimento, refazendo a franja.
  //
  // O refino local desbalanceia: medido 2,8x nesta maquina.  Reparticionar
  // custou de 3 a 6% do tempo de passo e devolveu o balanco -- numeros da
  // secao de escalabilidade do relatorio.
  //
  // Em np=1 e' no-op por construcao (a particao de um processo e' trivial), e
  // nao cobra nada.
  void reparticiona();


  //! A floresta crua, para quem precisa da API do t8code direto.  A posse fica
  //! com a Malha: nao chamar `unref` nisto.
  t8_forest_t forest() const { return f_; }

 private:
  //! Monta a franja: ids, geometria e niveis dos fantasmas, por TROCA.
  //
  // Chamada depois de CADA commit que muda a floresta.  Esquecer uma chamada
  // deixa a franja velha apontando para elementos que mudaram de dono -- e isso
  // nao falha, apenas responde errado.
  void monta_franja();

  sc_MPI_Comm comm_ = sc_MPI_COMM_SELF;
  int rank_ = 0;
  int nprocs_ = 1;
  long n_locais_ = 0;
  long n_fantasmas_ = 0;
  int camadas_ = 1;
  long saturacoes_ = 0;

  // Indexados pelo INDICE DE FLORESTA: locais primeiro, franja depois.
  std::vector<long> id_g_;
  std::vector<double> cen_;     // 2 por folha
  std::vector<double> lad_;     // 2 por folha
  std::vector<int> niv_;
  std::unordered_map<long, long> idx_por_id_;

  //! Vizinhos de face por indice: 4 faces, ate' 2 vizinhos cada (razao 2:1).
  //
  // Preenchido pelo t8code para as locais e pelas RODADAS DE TROCA para as de
  // franja.  Vazio para as folhas da ultima camada, que nao sao expansiveis.
  std::vector<std::array<std::vector<long>, 4>> viz_;
  std::vector<std::array<std::vector<int>, 4>> viz_niv_;
  std::vector<char> expansivel_;

  //! Monta a franja com `k` camadas.  Chamada por `monta_franja` e por
  //! `exige_camadas`.
  void monta_franja_camadas(int k);

  std::vector<Bloco> blocos_;
  t8_forest_t f_ = nullptr;
  long juncoes_ = 0;
};

}  // namespace t8flow2d

#endif
