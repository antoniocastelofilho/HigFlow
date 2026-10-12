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

  //! Maior e menor nivel presentes na malha.
  int nivel_minimo() const;
  int nivel_maximo() const;

  //! A floresta crua, para quem precisa da API do t8code direto.  A posse fica
  //! com a Malha: nao chamar `unref` nisto.
  t8_forest_t forest() const { return f_; }

 private:
  std::vector<Bloco> blocos_;
  t8_forest_t f_ = nullptr;
  long juncoes_ = 0;
};

}  // namespace t8flow2d

#endif
