// O SOLVER LINEAR DO T8FLOW 2D.
//
// Fachada sobre o `solver` herdado da HiGTree (decisao 1), que e' uma vtable
// sobre o sistema linear e ja' nao sabe o que e' uma malha: fala em indices
// (`slv_create(engine, first_gid, local_size)`, `set_Ai(i, numjs, j, v)`) e tem
// ZERO referencias a `hig_cell`, `higcit`, `sim_domain` ou `psd_` em todos os
// backends.  O proprio cabecalho declara o contrato: "everything above this file
// assembles rows and asks for a solution without knowing which library solves
// it."
//
// O QUE A FACHADA ACRESCENTA, e cada item vem de uma armadilha do herdado:
//
//   * TEMPO DE VIDA.  `slv_create`/`slv_destroy` viram construtor e destrutor.
//   * A ORDEM DE MONTAGEM.  O herdado expoe `slv_assemble_matrix`,
//     `slv_assemble_rhs` e `slv_assemble_output` separadamente, e a ORDEM entre
//     elas importa: a pre-alocacao da matriz do PETSc consome as linhas
//     acumuladas.  Chamar na ordem errada custou um SIGSEGV dentro do PETSc.
//     Aqui ha' `monta()`, que faz a ordem certa (rhs -> output -> matriz, como o
//     `slv_assemble` do herdado), e nao ha' como fazer a errada.
//   * O NOME DO MOTOR, com a aritmetica em UM lugar.  As duas listas do
//     `solver.h` -- `engine_names[]` e o enum `SolverEngine` -- estao DESLOCADAS
//     DE UMA POSICAO (`engine_names[i]` corresponde a `(SolverEngine)(i+1)`), e
//     nada confere isso em tempo de compilacao: errar faz o programa ANUNCIAR um
//     motor e rodar outro.  `nome_do_motor()` concentra essa conta, e um teste a
//     cobra.
//
//     LACUNA DO HERDADO, e ela e' real: NAO ha' como perguntar qual motor foi
//     INSTALADO -- nao existe `slv_get_engine`, e as vtables dos backends sao
//     estaticas em cada unidade.  Com `Motor::Qualquer` a escolha sai de uma
//     variavel de ambiente ou de heuristica, e a unica pista e' a linha que o
//     herdado imprime em stdout ("Solver engine used: ...").  Logo
//     `nome_do_motor()` devolve o que foi PEDIDO, nao o que roda; quem precisa de
//     garantia pede o motor explicitamente em vez de `Qualquer`.
//   * TODA LINHA PREENCHIDA.  O caminho de pre-alocacao do PETSc afirma
//     `col_ids[i] && vals[i]` para TODA linha local: deixar uma linha sem
//     `linha()` nao e' sistema singular, e' abort dentro do herdado.  `monta()`
//     recusa antes, dizendo QUAL linha falta.
//
// O QUE A FACHADA NAO ESCONDE: a numeracao.  `primeiro_id` e `n_local` sao a
// unica coisa que o solver precisa da malha -- um id global por incognita,
// contiguo por processo --, e e' exatamente o que a curva de preenchimento do
// t8code da' naturalmente.  Manter isso explicito e' deliberado.
#ifndef T8FLOW_SOLVER2D_HXX
#define T8FLOW_SOLVER2D_HXX

#include <cstddef>
#include <string>
#include <vector>

struct solver;   // o herdado, opaco aqui

namespace t8flow {

//! Motor do sistema linear.  `Qualquer` deixa o herdado decidir em execucao
//! (variavel de ambiente `HIGTREE_SOLVER_ENGINE`, ou heuristica).
enum class Motor { Qualquer, Petsc, Hypre };

class Solver {
 public:
  //! `primeiro_id` e' o menor id global deste processo; `n_local` quantas
  //! incognitas ele possui.  Os ids locais vao de `primeiro_id` a
  //! `primeiro_id + n_local - 1`.
  Solver(Motor motor, int primeiro_id, std::size_t n_local, int max_nao_zeros);
  ~Solver();

  Solver(const Solver &) = delete;
  Solver &operator=(const Solver &) = delete;

  //! Nome do motor PEDIDO, pela tabela do herdado.  Ver a lacuna no cabecalho:
  //! o motor instalado nao e' consultavel.  Devolve "qualquer" para
  //! `Motor::Qualquer`, que nao tem entrada na tabela.
  std::string nome_do_motor() const;

  //! Define a linha `id_global`: colunas e valores.  Uma vez por linha.
  void linha(int id_global, const std::vector<int> &colunas,
             const std::vector<double> &valores);

  //! Lado direito da equacao `id_global`.
  void lado_direito(int id_global, double v);

  //! Monta o sistema na ordem correta.  Recusa se alguma linha local nao foi
  //! definida.
  void monta();

  void resolve();

  //! Valor da incognita `id_global` na solucao.
  double x(int id_global) const;

  //! A solucao local inteira, na ordem dos ids (`primeiro_id` em diante).
  std::vector<double> solucao() const;

  //! Se a iteracao atingiu a tolerancia.  E' o UNICO relato de que ela atingiu:
  //! um motor que para na iteracao maxima tambem devolve, com residuo e tudo.
  //! Por isso nenhum oraculo deste projeto se apoia so' nisto -- o criterio e' o
  //! valor da solucao.
  bool convergiu() const;

  //! Iteracoes gastas na ultima resolucao.
  int iteracoes() const;

  //! Residuo relatado pelo motor na ultima resolucao.  Diagnostico, nao oraculo:
  //! quem o calcula e' o proprio solver.
  double residuo() const;

  //! Tolerancia da iteracao (a mesma para relativa e absoluta, como no herdado).
  void tolerancia(double eps);

  //! Teto de iteracoes.  Atingi-lo NAO e' erro para o herdado: ele retorna, e so'
  //! `convergiu()` distingue.
  void max_iteracoes(int n);

 private:
  solver *s_ = nullptr;
  Motor motor_ = Motor::Qualquer;
  int primeiro_ = 0;
  std::size_t n_ = 0;
  bool montado_ = false;
  std::vector<bool> linha_definida_;
};

}  // namespace t8flow

#endif
