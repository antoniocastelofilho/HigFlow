#include "ns-example-2d.h"


/******************************************************************************************/
/******************************************************************************************/
/***************************** viscoelastic user functions ********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial kernel conformation tensor
real get_tensor(Point center, int i, int j, real t) {
    return get_kernel(i, 1.0, 0.0) * (i == j);
}

// interpolating quadratic - coincides up to second order derivative with the log, is invertible and differentiable at 0
real log_interp_quadratic(real lambda, real tol) { 
    return -lambda*lambda/(2.0*tol*tol) + 2.0*lambda/tol + log(tol) - 1.5;
}
real log_interp_quadratic_inverse(real lambda, real tol) { // inverse of the interpolating quadratic
    return tol*(2.0 - sqrt(1.0 + 2.0*log(tol) - 2.0*lambda)); 
}
real log_interp_quadratic_jacobian(real lambda, real tol) { // derivative of the interpolating quadratic
    return -lambda/(tol*tol) + 2.0/tol; 
}
// interpolating line - coincides up to first order derivative with the log, is invertible and differentiable at 0
real log_interp_linear(real lambda, real tol) { 
    return lambda/tol + log(tol) - 1.0;
}
real log_interp_linear_inverse(real lambda, real tol) { // inverse of the interpolating line
    return tol*(lambda - log(tol) + 1.0); 
}
real log_interp_linear_jacobian(real lambda, real tol) { // derivative of the interpolating line
    return 1.0/tol; 
}

real get_kernel(int dim, real lambda, real tol) {
    if (FLT_EQ(lambda, 0.0)) {
        if (lambda < max_neg_lambda) max_neg_lambda = lambda;
        num_neg_lambda++;
    }
    real value;
    // if (lambda < tol)
    //   value = log_interp_quadratic(lambda, tol);
    // else
    // value = log(lambda);
    //if (lambda < tol)
    //   value = sqrt(tol);
    //else
    // value = sqrt(lambda);
    value = lambda;
    return value;
}

real get_kernel_inverse(int dim, real lambda, real tol) {
    // if(lambda < 1.0e-12){
    //     printf("lambda = %lf\n", lambda);
    // }
    real value;
    // if (lambda < log(tol))
    //   value = log_interp_quadratic_inverse(lambda, tol);
    // else
    // value = exp(lambda);
    // value = lambda*lambda;
    value = lambda;
    return value;
}

real get_kernel_jacobian(int dim, real lambda, real tol) {
    real value;
    // if (lambda < tol)
    //   value = log_interp_quadratic_jacobian(lambda, tol);
    // else
    // value = 1.0/lambda;
    //if (lambda < tol)
    //   value = 0.5/sqrt(tol);
    //else
    // value = 0.5/sqrt(lambda);
    value = 1.0;
    return value;
}
// Define the user function for viscoelastic flow
void calculate_m_user(real lambda[DIM], real jlambda[DIM],  real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par) {
    // // Calculate the matrix MM and BB for Oldroyd-B model
    // for (int i = 0; i < DIM; i++) {
    //     for (int j = i+1; j < DIM; j++) {
    //         M_aux[i][j] = 0.0;
    //         M_aux[j][i] = 0.0;
    //     }
    //     M_aux[i][i]  = (1.0-lambda[i])*jlambda[i];
    // }
}



/******************************************************************************************/
/******************************************************************************************/
/*********************** multiphase viscoelastic user functions ***************************/
/******************************************************************************************/
/******************************************************************************************/

// initial multiphase kernel conformation tensor
real get_tensor_multiphase(real fracvol, Point center, int i, int j, real t) {
    return get_kernel(i, 1.0, 0.0) * (i == j);
}

void calculate_m_user_multiphase(real fracvol, real lambda[DIM], real jlambda[DIM],  real M_aux[DIM][DIM], real Re, real trS, ve_parameters *par0, ve_parameters *par1) {
    // // Calculate the matrix MM and BB for Oldroyd-B model
    // for (int i = 0; i < DIM; i++) {
    //     for (int j = i+1; j < DIM; j++) {
    //         M_aux[i][j] = 0.0;
    //         M_aux[j][i] = 0.0;
    //     }
    //     real jlambda = get_kernel_jacobian(i, lambda[i], tol);
    //     M_aux[i][i]  = (1.0-lambda[i])*jlambda;
    // }
}



/******************************************************************************************/
/******************************************************************************************/
/******************** viscoelastic integral user functions ********************************/
/******************************************************************************************/
/******************************************************************************************/


// Value of the Tensor
real get_tensor_integral(Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}