// Small dense linear algebra on DIM x DIM matrices, for the constitutive equations.
// Everything here is fixed-size and local to a cell -- nothing distributed.
//
// hig_flow_jacobi() assumes A is SYMMETRIC and does not check it; on a non-symmetric
// input it converges to something that is not an eigendecomposition.  That is the
// intended use (the conformation tensor is symmetric by construction), but a path
// that loses symmetry to round-off and feeds it here gets no warning.

// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Linear Algebra
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_LINEAR_ALGEBRA
#define HIG_FLOW_LINEAR_ALGEBRA

#include "types.h"

// Autovalores (d) e autovetores (V) de uma matriz simetrica A, por rotacoes
// de Jacobi.  Usada pelas variantes viscoelasticas para diagonalizar o tensor
// conformacao.
void hig_flow_jacobi(real A[DIM][DIM], real d[DIM], real V[DIM][DIM]);

// Resolve A x = b por eliminacao gaussiana com pivotamento parcial, com b
// armazenado na ultima coluna de A.
void hig_flow_solve_system_constitutive_equation(int n, real A[DIM*DIM][DIM*DIM+1], real x[DIM*DIM]);

// B = A R  e  B = A^T R, para matrizes DIMxDIM.
void hig_flow_matrix_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]);
void hig_flow_matrix_transpose_product (real A[DIM][DIM], real R[DIM][DIM], real B[DIM][DIM]);

#endif
