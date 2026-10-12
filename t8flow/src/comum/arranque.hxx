// O ARRANQUE DO T8FLOW, e a ordem dele NAO E' LIVRE.
//
// Descoberto com um SIGSEGV dentro do PETSc, sem uma linha que apontasse para o
// nosso codigo.  O mecanismo, lido do fonte herdado:
//
//   `_slv_create_petsc` chama `_try_initialize_petsc` (higtree/src/utils.c:34),
//   que desreferencia um `static char **_argv` -- e QUEM preenche esse estatico
//   e' o `higtree_initialize`.  Sem ele, `*_argv` e' desreferencia de nulo.  E a
//   flag `is_petsc_initialized` tambem e' `static`: inicializar o PETSc por fora
//   NAO a desarma, entao o solver tenta de novo e cai do mesmo jeito.
//
// Logo o nucleo herdado tem estado global de arranque, e a consequencia e' de
// projeto: existe UM arranque no t8Flow, com esta ordem --
//
//   1. `higtree_initialize`  MPI_Init, e guarda argc/argv para o PETSc;
//   2. `sc_init` + `t8_init` por processo, com o MPI JA' inicializado;
//   3. o PETSc entra preguicoso, pelo proprio solver, agora com argv valido.
//
// NAO chamar `sc_MPI_Init` (dois MPI_Init sao erro) nem os Finalize: o
// `higtree_initialize` registra `_higtree_finalize` no `atexit`, que fecha PETSc
// e MPI.  Daqui so' sai o `sc_finalize`.
//
// Por isso isto e' um objeto e nao uma funcao: a ordem fica no construtor, uma
// vez, e nao em cada `main` que alguem escrever.  Vive em `comum/` porque nao
// tem dimensao -- e' o mesmo arranque para o 2D e para o 3D.
#ifndef T8FLOW_ARRANQUE_HXX
#define T8FLOW_ARRANQUE_HXX

namespace t8flow {

class Arranque {
 public:
  //! Recebe argc/argv do `main` porque o PETSc herdado os consome.
  Arranque(int *argc, char ***argv);
  ~Arranque();

  Arranque(const Arranque &) = delete;
  Arranque &operator=(const Arranque &) = delete;

  //! Rank deste processo em MPI_COMM_WORLD, lido uma vez no arranque.
  int rank() const { return rank_; }

  //! Numero de processos em MPI_COMM_WORLD.
  int num_ranks() const { return num_ranks_; }

 private:
  int rank_ = 0;
  int num_ranks_ = 1;
};

}  // namespace t8flow

#endif
