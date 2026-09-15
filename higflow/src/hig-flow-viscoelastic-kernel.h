// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Viscoelastic Kernel
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_VISCOELASTIC_KERNEL
#define HIG_FLOW_VISCOELASTIC_KERNEL

#include "types.h"

// Maquinaria do tensor conformacao compartilhada pelas variantes
// viscoelasticas (viscoelastic, elastoviscoplastic, variable-viscosity).
void hig_flow_kernel_rhs (real De, real K[DIM][DIM], real O[DIM][DIM], real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]);
void hig_flow_implicit_kernel_rhs (real De, real B[DIM][DIM], real M[DIM][DIM], real RHS[DIM][DIM]);
void hig_flow_kernel_system_matrix (real w[DIM*DIM][DIM*DIM+1], real Omega[DIM][DIM], real dt);

#endif
