#include "arranque.hxx"

#include <mpi.h>

#include <t8.h>

// O nucleo herdado e' C compilado como C++ na MESMA norma (decisao 5), entao
// nada de `extern "C"` -- os simbolos vem com encadeamento C++ dos dois lados.
#include "utils.h"

namespace t8flow {

Arranque::Arranque(int *argc, char ***argv)
{
  // 1. MPI, e o estado global que o PETSc herdado vai precisar.
  higtree_initialize(argc, argv);

  // 2. libsc e t8code: por processo, nao por chamada, e com o MPI ja' de pe'.
  sc_init(sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ERROR);
  t8_init(SC_LP_ERROR);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);

  // 3. o PETSc nao entra aqui: o solver herdado o inicializa preguicoso, e
  //    fazer isso agora so' adiantaria um custo sem remover nenhum.
}

Arranque::~Arranque()
{
  // `MPI_Finalize` e `PetscFinalize` sao do `atexit` que o `higtree_initialize`
  // registrou.  Chamar aqui seria chamar duas vezes.
  sc_finalize();
}

}  // namespace t8flow
