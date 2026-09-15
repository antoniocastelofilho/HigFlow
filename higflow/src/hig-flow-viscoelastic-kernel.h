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

// Omega = R Omega_aux R^t.  Quando lambda[i] e lambda[j] coincidem a menos de
// `small` o termo de rotacao e' zerado: nessa faixa Omega multiplicaria uma
// diferenca de kernel nula, entao o valor dele nao chega ao RHS.  O termo
// simetrico que sobreviveria ao limite, (M_ij+M_ji)*lambda*jlambda, nao e'
// recuperado aqui — o hig_flow_calculate_b zera a off-diagonal.
void hig_flow_calculate_omega (real lambda[DIM], real R[DIM][DIM], real M[DIM][DIM], real Omega[DIM][DIM], real small);

#endif
