// A MALHA EM TRES DIMENSOES: floresta de hexaedros, multibloco, com franja.
//
// POR QUE UM ARQUIVO SEPARADO e nao um `#if DIM == 3` sobre o de duas dimensoes:
// decisao 7 do projeto, e a razao e' que cada laco sobre direcoes, cada contagem
// de vizinhos de face e cada formula de quadratura ganharia um caso que nao e'
// exercitado no regime em que se esta' trabalhando.  Com bases separadas, o
// codigo tridimensional e' lido e testado como codigo tridimensional.
//
// O QUE MUDA DE VERDADE EM RELACAO A DUAS DIMENSOES -- e e' mais que trocar um
// 2 por um 3:
//
//   HEX tem SEIS faces, nao quatro, e oito vertices, nao quatro.  A ordem dos
//   vertices que o `cmesh` do t8code espera e' a de bits (x mais rapido, depois
//   y, depois z), e errar essa ordem nao falha -- produz um hexaedro torcido com
//   volume positivo;
//
//   na transicao 2:1 uma face tem ate' QUATRO vizinhos, nao dois, porque a face
//   de um hexaedro e' partida em quatro quando o vizinho refina;
//
//   a franja custa mais, e foi MEDIDA antes deste arquivo existir
//   (`bench/franja3d`, e a secao 18.3 do projeto).  O resultado mudou uma
//   decisao de arquitetura: o teto teorico de 24 vizinhos por celula (seis faces
//   vezes quatro) NUNCA e' atingido -- o maximo medido e' 9, mesmo em malha
//   adaptada --, e dimensionar a carga pelo teto desperdicaria 2,7 vezes.  Esta
//   classe mede os vizinhos na construcao e dimensiona pelo que mediu.
//
// OS QUATRO MODOS DE FALHA DO PARALELO JA' VEM CONSERTADOS, e nao redescobertos:
// eles estao mapeados na secao 15 do projeto e nas clausulas X1 a X9 do contrato
// de duas dimensoes.  Em ordem de quanto custaram:
//
//   1. a franja NAO pode ser ordenada por id.  O t8code numera as folhas de
//      franja em `num_local .. num_local+num_ghosts-1`, e so' as camadas 2 e
//      acima podem ser reordenadas;
//   2. face de franja NAO e' face de contorno.  As duas dao lista de vizinhos
//      vazia, e confundi-las aplica condicao de contorno no meio do dominio;
//   3. o nivel do vizinho VIAJA JUNTO com o id.  Busca-lo depois custa uma
//      camada inteira de franja, e o valor de reserva "nivel da propria celula"
//      e' pior que a recusa;
//   4. NENHUMA chamada coletiva em caminho que divirja entre processos.  Vale
//      para dentro desta classe: `id_do_ponto` usa os niveis LOCAIS.

#ifndef T8FLOW_3D_MALHA3D_HXX
#define T8FLOW_3D_MALHA3D_HXX

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

namespace t8flow3d {

//! Um ponto em tres dimensoes.
using Ponto3 = std::array<double, 3>;

//! Um BLOCO: a caixa alinhada aos eixos e o refino uniforme dentro dela.
//
// Um bloco vira uma ARVORE do `cmesh` do t8code.  Multibloco e multiarvore sao o
// mesmo conceito, e foi a peca que encaixou sem atrito entre este sistema e a
// biblioteca.
struct Bloco {
  Ponto3 lo;   //!< canto menor
  Ponto3 hi;   //!< canto maior
  int nivel;   //!< refino uniforme DENTRO do bloco: 8^nivel folhas
};

class Malha3 {
 public:
  //! Constroi a floresta a partir dos blocos.  Uma arvore por bloco, unidas
  //! pelas faces que coincidem EXATAMENTE.
  //
  // A franja nasce com uma camada (a nativa do t8code).  Pedir mais e' com
  // `exige_camadas`, e quem monta o esquema decide quantas -- a carga cresce
  // com o cubo do raio em tres dimensoes.
  Malha3(const std::vector<Bloco> &blocos, sc_MPI_Comm comm);
  ~Malha3();

  Malha3(const Malha3 &) = delete;
  Malha3 &operator=(const Malha3 &) = delete;

  // ---- os blocos -----------------------------------------------------------
  long num_blocos() const { return static_cast<long>(blocos_.size()); }
  const Bloco &bloco(long b) const;

  //! Quantas faces foram unidas entre blocos vizinhos.
  long num_juncoes() const { return juncoes_; }

  // ---- as folhas -----------------------------------------------------------
  long num_folhas() const;
  long num_folhas_do_bloco(long b) const;
  long num_folhas_globais() const;
  long primeiro_id() const;
  long id_global(long b, long i) const;

  Ponto3 centroide(long b, long i) const;
  std::vector<Ponto3> centroides() const;

  //! Os tres lados da folha.  Saem do BLOCO e do NIVEL, e nao do diametro do
  //! elemento: a mesma formula aqui e na franja e' o que impede que a folha
  //! local e a mesma folha vista como franja tenham medidas diferentes.
  Ponto3 lados_da_folha(long b, long i) const;

  //! O VOLUME da folha -- em tres dimensoes, o produto dos tres lados.
  double volume(long b, long i) const;

  int nivel(long b, long i) const;
  int nivel_minimo() const;
  int nivel_maximo() const;

  // ---- vizinhanca ----------------------------------------------------------
  //
  // O HEX TEM SEIS FACES: 0 e 1 sao `x-` e `x+`, 2 e 3 sao `y-` e `y+`, 4 e 5
  // sao `z-` e `z+`.  A convencao e' a do t8code, e segui-la em vez de inventar
  // outra e' o que permite passar a face direto para a biblioteca.
  int num_vizinhos(long b, long i, int fa) const;
  std::vector<long> ids_vizinhos(long b, long i, int fa) const;
  std::vector<int> niveis_vizinhos(long b, long i, int fa) const;

  //! Quantas folhas tem vizinho de nivel diferente em alguma face.
  long num_folhas_em_transicao() const;

  //! O MAIOR numero de vizinhos que uma folha desta malha tem, somando as seis
  //! faces.  Medido na construcao, e nao presumido.
  //
  // O TETO TEORICO E' 24 -- seis faces vezes quatro, porque a face de um
  // hexaedro e' partida em quatro quando o vizinho refina.  Medido
  // (`bench/franja3d`): o maximo real e' 9 em malha adaptada e 6 em uniforme,
  // porque uma celula tem no maximo uma face com vizinhos refinados em cada
  // direcao, e nao todas as seis.
  //
  // E' por isso que este numero existe: a carga da troca da franja e'
  // dimensionada por ele, e dimensiona-la pelo teto custaria 2,7 vezes mais
  // memoria -- em tres dimensoes, a diferenca entre 55 KB e 24 KB por folha.
  int max_vizinhos_por_folha() const { return max_viz_; }

  // ---- refino --------------------------------------------------------------
  void refina_caixa(const Ponto3 &lo, const Ponto3 &hi, int nivel_alvo);
  void refina_aleatorio(double fracao, unsigned long semente, int nivel_max);

  // ---- busca geometrica ----------------------------------------------------
  //
  // Devolve o id GLOBAL da folha que contem o ponto, ou `-1`.  A regra de
  // desempate e' SEMIABERTA `[lo, hi)`, e ela e' global: nao depende de qual
  // processo pergunta.  E' dela que vem a conservacao do esquema de fluxo em
  // paralelo -- os dois processos que compartilham uma faceta chegam a' mesma
  // celula.
  //
  // A busca e' O(niveis), por indice espacial.  Uma varredura linear aqui torna
  // a montagem do operador QUADRATICA, e isso ja' aconteceu em duas dimensoes:
  // 8,74 s contra 0,679 s no mesmo problema, com o sintoma se disfarcando de
  // ganho de paralelismo (secao 17.1 do projeto).
  long id_do_ponto(const Ponto3 &p) const;

  //! O ponto esta' dentro de algum BLOCO, geometricamente?
  //
  // SEPARADA DE `id_do_ponto` porque as duas respostas negativas sao
  // DIFERENTES e confundi-las esconde a segunda: um ponto fora do dominio e'
  // legitimo -- o contorno trata --, e um ponto DENTRO do dominio que
  // `id_do_ponto` nao acha significa que falta FRANJA.
  //
  // Sem a distincao, o suporte por caixa descarta os dois em silencio, e a
  // medida de "quantas camadas bastam" da' zero truncamento com qualquer
  // numero de camadas -- que foi exatamente o resultado bom demais que
  // apareceu antes desta funcao existir.
  bool dentro_do_dominio(const Ponto3 &p) const;

  // ---- a franja ------------------------------------------------------------
  long num_fantasmas() const { return n_fantasmas_; }
  int camadas_de_franja() const { return camadas_; }
  void exige_camadas(int k);

  long id_global_do_indice(long idx) const;
  Ponto3 centroide_do_indice(long idx) const;
  Ponto3 lados_do_indice(long idx) const;
  int nivel_do_indice(long idx) const;
  long indice_do_id(long id) const;

  std::vector<long> vizinhos_do_indice(long idx, int fa) const;
  std::vector<int> niveis_vizinhos_do_indice(long idx, int fa) const;

  //! Esta folha pode abrir mais um anel, ou esta' na borda da franja?
  //
  // Distingue "face no contorno do DOMINIO" de "face sobre a qual nao me
  // contaram".  As duas dao lista de vizinhos vazia, e confundi-las faz o
  // esquema aplicar condicao de contorno no MEIO do dominio.
  bool expansivel_do_indice(long idx) const;

  long bytes_publicados_na_franja() const { return bytes_pub_; }
  long saturacoes_da_franja() const { return saturacoes_; }

  //! Quantas RODADAS de troca a franja precisou, e quantos pontos de caixa
  //! ficaram faltando ao parar.
  //
  // UMA RODADA SO' ALCANCA VIZINHO DE FACE.  A franja nativa do t8code e'
  // `T8_GHOST_FACES` -- a documentacao diz que so' esse tipo e' suportado --,
  // de modo que o processo diagonalmente oposto numa quina nao e' vizinho de
  // face de ninguem que fale comigo.  Em oito processos isso dava 192 suportes
  // truncados no centro do dominio, e AUMENTAR AS CAMADAS NAO RESOLVIA: a
  // informacao nao precisava ir mais longe, precisava ir por outro caminho.
  //
  // `faltando_na_franja` diferente de zero depois de `exige_camadas` significa
  // que o esquema NAO tera' o suporte que pediu em parte do dominio -- e e'
  // melhor ler isso aqui do que descobrir pelo resultado que muda com `np`.
  //! Os tempos das fases da fase B, para saber O QUE otimizar.  Em segundos.
  struct TemposFranja {
    double reindexa, reduz, zera, publica, troca, funde;
    double prologo, achata;
  };
  TemposFranja tempos_da_franja() const {
    return {t_reindexa_, t_reduz_, t_zera_, t_publica_, t_troca_,
            t_funde_, t_prologo_, t_achata_};
  }

  int rodadas_de_franja() const { return rodadas_franja_; }
  long faltando_na_franja() const { return faltando_franja_; }

  // ---- particao ------------------------------------------------------------
  int num_processos() const { return nprocs_; }
  int processo() const { return rank_; }
  sc_MPI_Comm comunicador() const { return comm_; }

  //! `folhas_do_maior / folhas_do_menor`.  Coletiva.
  double desbalanceamento() const;

  //! Reparticiona pela curva de preenchimento, refazendo a franja.
  //
  // EXPLICITA, e nao automatica, porque o t8code JA' tenta balancear e o
  // comportamento dele nao e' previsivel -- quatro medidas do mesmo
  // `set_balance(..., 0)` deram 1,00, 1,11, 2,8 e a promessa da documentacao.
  // Ver a secao 17.3 do projeto.
  void reparticiona();

 private:
  void monta_franja();
  void monta_franja_camadas(int k);
  void mede_vizinhos();

  t8_forest_t f_ = nullptr;
  std::vector<Bloco> blocos_;
  long juncoes_ = 0;

  sc_MPI_Comm comm_ = sc_MPI_COMM_SELF;
  int rank_ = 0;
  int nprocs_ = 1;
  long n_locais_ = 0;
  long n_fantasmas_ = 0;
  int camadas_ = 1;
  int max_viz_ = 6;        //!< medido em `mede_vizinhos`; 6 e' o de malha uniforme
  long saturacoes_ = 0;
  long bytes_pub_ = 0;
  int rodadas_franja_ = 0;
  double t_reindexa_ = 0, t_reduz_ = 0, t_zera_ = 0;
  double t_publica_ = 0, t_troca_ = 0, t_funde_ = 0;
  double t_prologo_ = 0, t_achata_ = 0;
  long faltando_franja_ = 0;

  // Indexados pelo INDICE DE FLORESTA: locais primeiro, franja depois.
  std::vector<long> id_g_;
  std::vector<double> cen_;     //!< 3 por folha
  std::vector<double> lad_;     //!< 3 por folha
  std::vector<int> niv_;
  std::vector<std::array<std::vector<long>, 6>> viz_;
  std::vector<std::array<std::vector<int>, 6>> viz_niv_;
  //! A MASCARA de faces conhecidas, um bit por face; `0x3F` sao as seis.
  //! `expansivel_do_indice` pergunta por `== 0x3F`, e a mascara inteira fica
  //! guardada para que a tabela possa ser reconstruida sem refazer a fase A.
  std::vector<char> expansivel_;

  std::unordered_map<long, long> idx_por_id_;

  //! Indice espacial `(bloco, nivel, i, j, k)` -> indice de floresta.
  std::unordered_map<unsigned long long, long> ponto_por_celula_;

  //! Niveis das folhas que ESTE processo conhece.  Sem reducao, de proposito:
  //! `id_do_ponto` e' chamada de caminhos que divergem entre os processos.
  int nivel_local_min_ = -1;
  int nivel_local_max_ = -1;
};

}  // namespace t8flow3d

#endif  // T8FLOW_3D_MALHA3D_HXX
