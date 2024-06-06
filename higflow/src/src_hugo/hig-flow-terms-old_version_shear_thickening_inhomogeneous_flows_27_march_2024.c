// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-terms.h"

// *******************************************************************
// Navier-Stokes Equation Components
// *******************************************************************

// Pressure term contribution for the Navier-Stokes equation
real higflow_pressure_term(higflow_solver *ns)
{
    // Set the pressure term
    real value = 0.0;
    if (ns->contr.projtype == 1)
    {
        value = ns->cc.dpdx;
        if (ns->contr.flowtype == 2)
        {
            value /= ns->cc.dens;
        }
    }
    return value;
}

// Cell source term contribution for the Navier-Stokes equation
real higflow_source_term(higflow_solver *ns)
{
    // Set the source term
    real value = ns->cc.F;
    if (ns->contr.flowtype == 2)
    {
        value /= ns->cc.dens;
    }
    return value;
}

// Cell term contribution for the interfacial tension
real higflow_interfacial_tension_term(higflow_solver *ns)
{
    real Bo = 0.18;
    real value = ns->cc.IF;
    // real value = 24.5*ns->cc.IF;
    value = value / (Bo * ns->cc.dens);
    return value;
}

// Cell term contribution for the gravity
real higflow_gravity_term(higflow_solver *ns)
{
    // real Fr = 1.0;
    //  real value = ns->cc.curv;
    // real value = 1.0/pow(Fr,2.0);
    // real value = 9.8;
    real value = 0.0;
    //Bousinessq aprooximation for non-isothermal flows
    if (ns->contr.modelflowtype == 5)
    {
        value = -1.0;
        value *= ns->cc.Bous;
    }
    return value;
}

// Tensor term contribution for the Navier-Stokes equation
real higflow_tensor_term(higflow_solver *ns)
{
    // Set the tensor term
    real value = 0.0;
    // Non Newtonian contribution
    if (ns->contr.flowtype == 3)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }
    // Non Newtonian contribution
    if (ns->contr.flowtype == 4)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }

    if (ns->contr.flowtype == 2)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            // value += ns->cc.dSdx[dim2];
        }
    }
    // Viscoelastic flow with variable viscosity
    if (ns->contr.flowtype == 6)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }
    // Viscoelastic flow with shear-banding
    if (ns->contr.flowtype == 7)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }
    // Elastoviscoplastic flow
    if (ns->contr.flowtype == 8)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }
    // Shear-thickening suspensions
    if (ns->contr.flowtype == 9)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.dSdx[dim2];
        }
    }
    return value;
}

// Convective term: Quick scheme
real higflow_convective_term_quick(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        value += (ns->cc.v[dim2] + fabs(ns->cc.v[dim2])) * 0.1875 * ns->cc.dudxl[dim2];
        value += (ns->cc.v[dim2] + fabs(ns->cc.v[dim2])) * 0.3750 * ns->cc.dudxr[dim2];
        value -= (ns->cc.v[dim2] + fabs(ns->cc.v[dim2])) * 0.0625 * ns->cc.dudxrr[dim2];
        value += (ns->cc.v[dim2] - fabs(ns->cc.v[dim2])) * 0.3750 * ns->cc.dudxl[dim2];
        value += (ns->cc.v[dim2] - fabs(ns->cc.v[dim2])) * 0.1875 * ns->cc.dudxr[dim2];
        value -= (ns->cc.v[dim2] - fabs(ns->cc.v[dim2])) * 0.0625 * ns->cc.dudxll[dim2];
    }
    return value;
}

// Convective term: Cubista scheme
real higflow_convective_term_cubista(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        value += ns->cc.vc[dim2];
    }
    return value;
}

// Convective term: Modified Coefficient Upwind scheme
real higflow_convective_term_mcupwind(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        value += 0.5 * (ns->cc.v[dim2] + fabs(ns->cc.v[dim2])) * ns->cc.dudxl[dim2];
        value += 0.5 * (ns->cc.v[dim2] - fabs(ns->cc.v[dim2])) * ns->cc.dudxr[dim2];
        if (ns->cc.v[dim2] > 0.0)
        {
            if (fabs(ns->cc.dudxl[dim2]) > 1.0e-8)
            {
                value += 0.5 * ns->cc.v[dim2] * delta[dim2] * ns->cc.du2dx2[dim2] / ns->cc.dudxl[dim2];
            }
        }
        else
        {
            if (fabs(ns->cc.dudxr[dim2]) > 1.0e-8)
            {
                value -= 0.5 * ns->cc.v[dim2] * delta[dim2] * ns->cc.du2dx2[dim2] / ns->cc.dudxr[dim2];
            }
        }
    }
    return value;
}

// Convective term: second order scheme
real higflow_convective_term_second_order(higflow_solver *ns, Point delta, int dim)
{
    // Set the second order convective term
    real value;
    switch (ns->contr.secondconvecdiscrtype)
    {
    // Modified Coefficient Upwind scheme
    case 0:
        value = higflow_convective_term_mcupwind(ns, delta, dim);
        break;
    // Cubista scheme
    case 1:
        value = higflow_convective_term_cubista(ns, delta, dim);
        break;
    // Quick scheme
    case 2:
        value = higflow_convective_term_quick(ns, delta, dim);
        break;
    }
    return value;
}

// Convective term: first order scheme
real higflow_convective_term_first_order(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        value += 0.5 * (ns->cc.v[dim2] + fabs(ns->cc.v[dim2])) * ns->cc.dudxl[dim2];
        value += 0.5 * (ns->cc.v[dim2] - fabs(ns->cc.v[dim2])) * ns->cc.dudxr[dim2];
    }
    return value;
}

// Convective term: central scheme
real higflow_convective_term_central(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        value += ns->cc.v[dim2] * ns->cc.dudxc[dim2];
    }
    return value;
}

// Convective term contribution for the Navier-Stokes equation
real higflow_convective_term(higflow_solver *ns, Point delta, int dim)
{
    // Set the convective term according choice
    real value;
    switch (ns->cc.convec_type)
    {
    // Cental differening scheme
    case 0:
        value = higflow_convective_term_central(ns, delta, dim);
        break;
    // First order upwind scheme
    case 1:
        value = higflow_convective_term_first_order(ns, delta, dim);
        break;
    // Second order upwind scheme
    case 2:
        value = higflow_convective_term_second_order(ns, delta, dim);
        break;
    }
    return value;
}

// Difusive term contribution for the Navier-Stokes equation
real higflow_difusive_term(higflow_solver *ns, Point delta)
{
    real value = 0.0;
    switch (ns->contr.flowtype)
    {
    // Newtonian
    case 0:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        break;
    // Generalized Newtonian
    case 1:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        break;
    // Multiphase
    case 2:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        value /= ns->cc.dens;
        break;
    // Viscoelastic
    case 3:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        if (ns->contr.modelflowtype == 5)
        {
        value *= 1.0*(ns->cc.FT)*(ns->ed.ve.par.beta) + 1.0*(1.0-ns->ed.ve.par.beta);
        }
        break;
        // Viscoelastic
    case 4:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        break;
    // Viscoelastic with variable viscosity
    case 6:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        if (ns->contr.modelflowtype == 5)
        {
        value *= 1.0*(ns->cc.FT)*(ns->ed.vevv.par.beta) + 1.0*(1.0-ns->ed.vevv.par.beta);
        }
        break;
    // Viscoelastic with shear-banding
    case 7:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        // value /= ns->par.Re;
        // value *= ns->ed.vesb.par.De;
        value *= (1.0 + ns->ed.vesb.par.beta) * (ns->ed.vesb.par.De) / ns->par.Re;
        break;
    // Elastoviscoplastic
    case 8:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value /= ns->par.Re;
        break;
    // Shear-thickening suspension
    case 9:
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            value += ns->cc.du2dx2[dim2];
        }
        value *= 2.0 * (ns->ed.stsp.par.eta0);
        break;
    }
    return value;
}

// Cell electroosmotic source term contribution for the Navier-Stokes equation
real higflow_electroosmotic_source_term(higflow_solver *ns)
{
    // Set the source term
    real value = ns->cc.Feo;
    return value;
}

// Cell electroosmotic diffusive ionic term contribution for the Navier-Stokes equation
real higflow_diffusive_ionic_term(higflow_solver *ns)
{
    // Set the diffusive ionic term
    real value = ns->cc.d2ndx2 / ns->ed.eo.par.Pe;
    return value;
}

// Cell electroosmotic potential ionic term contribution for the Navier-Stokes equation
real higflow_potential_ionic_term(higflow_solver *ns)
{
    // Set the potential ionic term
    real value;
    value = ns->cc.dndx * (ns->cc.dphidx + ns->cc.dpsidx) + ns->cc.ncell * (ns->cc.d2psidx2 + ns->cc.d2phidx2);
    value *= ns->ed.eo.par.alpha / ns->ed.eo.par.Pe;
    return value;
}

// Cell the diffusive term contribution for the transport equation of the specie A
real higflow_diffusive_shear_banding_nA_term(higflow_solver *ns)
{
    // Set the diffusive term
    real value;
    value = ns->cc.d2nABdx2;
    value *= 2.0;
    value /= ns->ed.vesb.par.PeA;
    return value;
}

// Cell the diffusive term contribution for the transport equation of the specie B
real higflow_diffusive_shear_banding_nB_term(higflow_solver *ns)
{
    // Set the diffusive term
    real value;
    value = ns->cc.d2nABdx2;
    value *= 2.0;
    value /= ns->ed.vesb.par.PeB;
    return value;
}

// Cell the diffusive term contribution for the transport equation of the specie A
// real higflow_diffusive_shear_banding_conformation_tensor_A_term(higflow_solver *ns) {
// Set the diffusive ionic term
// real value = ns->cc.d2ABdx2/ns->ed.vesb.par.PeA;
// return value;
//}

// Cell the diffusive term contribution for the transport equation of the specie B
// real higflow_diffusive_shear_banding_conformation_tensor_B_term(higflow_solver *ns) {
// Set the diffusive ionic term
// real value = ns->cc.d2ABdx2/ns->ed.vesb.par.PeB;
// return value;
//}

// Cell particle migration equation: calculation of the diffusive term -4*(1-phi)3*(dphidx)*(dTdx)
real higflow_vol_frac_term1(higflow_solver *ns)
{
    // Set the terms
    
    real value = 0.0;
    for (int dim2 = 0; dim2 < DIM; dim2++)
    {
        //if (dim2 == 1) {
        //    value += ns->cc.dTdx[dim2];
        //} else {
        //    value +=  0.0;
        //}
        value += ns->cc.dTdx[dim2];
        //value += dim2;
    }
    
    //value = dim;
    //value = ns->cc.dTdx[dim];
    //real value = -ns->cc.phivfcell;
    value *= (ns->cc.dphivfdx);
    // value   = ns->cc.dndx*(ns->cc.dphidx + ns->cc.dpsidx) + ns->cc.ncell*(ns->cc.d2psidx2 + ns->cc.d2phidx2);
    value *= -4.0*(pow((1.0-ns->cc.phivfcell),3.0));
    value *= -2.0*(ns->par.Re)*(pow(ns->ed.stsp.par.apsize,2.0)) / (9.0*(ns->ed.stsp.par.eta0));
    //value = 8.0;
    return value;
}

// Cell particle migration equation: calculation of the diffusive term (1-phi)4*(d2Tdx2)
real higflow_vol_frac_term2(higflow_solver *ns, real varphic)
{
    // Set the potential ionic term
    real value = 0.0;
    //value += ns->cc.d2Tdx2;
    //value = (pow((1.0-0.1),4.0));
    value = (pow((1.0-varphic),4.0))*(ns->cc.d2Tdx2);
    value *= -2.0*(ns->par.Re)*(pow(ns->ed.stsp.par.apsize,2.0)) / (9.0*(ns->ed.stsp.par.eta0));
    return value;
}

// Cell energy equation diffusion term
real higflow_difussion_energy_equation_term(higflow_solver *ns)
{
    // Set the potential ionic term
    real value;
    real PrRe = (ns->par.Re)*(ns->ed.nif.par.Pr);
    //value = (ns->cc.dKTdx) * (ns->cc.dTempdx) + (ns->cc.KT)*(ns->cc.d2Tempdx2);
    value = ns->cc.d2Tempdx2;
    value /= PrRe;
    return value;
}

// Cell energy equation diffusion term
real higflow_implicit_difussion_energy_equation_term(higflow_solver *ns)
{
    // Set the potential ionic term
    real value;
    real PrRe = (ns->par.Re)*(ns->ed.nif.par.Pr);
    value = (ns->cc.dKTdx) * (ns->cc.dTempdx);
    value /= PrRe;
    return value;
}
