// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step Electro-osmotic - version 02/02/2018
// *******************************************************************
// *******************************************************************

#include "hig-flow-step-electroosmotic.h"
// *******************************************************************
// Navier-Stokes step elements
// *******************************************************************

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_nplus(higflow_solver *ns) {
    real bcval, fracvol;
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnplus);

    int num_bc_types = 2;
    bc_type bc_t;
    for(int i = 0; i < num_bc_types; i++) {
        if(i == 0) bc_t = DIRICHLET;
        else if(i == 1) bc_t = NEUMANN;

        // Get the number of boundaries of type
        int numbcs = sd_get_num_bcs(sdp, bc_t);
        // For each boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the boundary
            sim_boundary *bc = sd_get_bc(sdp, bc_t, i);
            // Get the id defined by the user
            int userid = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bclid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to get the pressure
                    real t = ns->par.t + ns->par.dt;
                    // Get the pressure defined by the user
                    if(ns->contr.flowtype == MULTIPHASE) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nplus(fracvol, userid, bccenter, t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_nplus(userid, bccenter, t);
                    // Set the value 
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns) {
    real bcval, fracvol;
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnminus);
    int num_bc_types = 2;
    bc_type bc_t;
    for(int i = 0; i < num_bc_types; i++) {
        if(i == 0) bc_t = DIRICHLET;
        else if(i == 1) bc_t = NEUMANN;

        // Get the number of boundaries of type
        int numbcs = sd_get_num_bcs(sdp, bc_t);
        // For each boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the boundary
            sim_boundary *bc = sd_get_bc(sdp, bc_t, i);
            // Get the id defined by the user
            int userid = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bclid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to get the pressure
                    real t = ns->par.t + ns->par.dt;
                    // Get the pressure defined by the user
                    if(ns->contr.flowtype == MULTIPHASE) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_nminus(fracvol, userid, bccenter, t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_nminus(userid, bccenter, t);
                    // Set the value 
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}


void higflow_boundary_condition_for_phi(higflow_solver *ns) {
    real bcval, fracvol;
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOphi);

    int num_bc_types = 2;
    bc_type bc_t;
    for(int i = 0; i < num_bc_types; i++) {
        if(i == 0) bc_t = DIRICHLET;
        else if(i == 1) bc_t = NEUMANN;

        // Get the number of boundaries of type
        int numbcs = sd_get_num_bcs(sdp, bc_t);
        // For each boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the boundary
            sim_boundary *bc = sd_get_bc(sdp, bc_t, i);
            // Get the id defined by the user
            int userid = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bclid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to get the pressure
                    real t = ns->par.t + ns->par.dt;
                    // Get the pressure defined by the user
                    if(ns->contr.flowtype == MULTIPHASE) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_phi(fracvol, userid, bccenter, t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_phi(userid, bccenter, t);
                    // Set the value 
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}


void higflow_boundary_condition_for_psi(higflow_solver *ns) {
    real bcval, fracvol;
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);

    int num_bc_types = 2;
    bc_type bc_t;
    for(int i = 0; i < num_bc_types; i++) {
        if(i == 0) bc_t = DIRICHLET;
        else if(i == 1) bc_t = NEUMANN;

        // Get the number of boundaries of type
        int numbcs = sd_get_num_bcs(sdp, bc_t);
        // For each boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the boundary
            sim_boundary *bc = sd_get_bc(sdp, bc_t, i);
            // Get the id defined by the user
            int userid = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bclid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to get the pressure
                    real t = ns->par.t + ns->par.dt;
                    // Get the pressure defined by the user
                    if(ns->contr.flowtype == MULTIPHASE) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_psi(fracvol, userid, bccenter, t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_psi(userid, bccenter, t);
                    // Set the value 
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}


// *******************************************************************
// Ionic transport equation 
// *******************************************************************

//Explicit ionic transport equation
void higflow_explicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == PNP) {
        real alphaeo = ns->ed.eo.par.alpha;
        real Pe = ns->ed.eo.par.Pe;
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnplus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnplus); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            // real u[DIM], dnplusdx[DIM];
            // hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            real nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                // Set the computational cell in this particular direction
                higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
                // Compute the diffusive ionic term rhs
                rhs    += higflow_diffusive_ionic_term(ns, Pe);
                // convective term
                switch (ns->ed.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
                        //hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                        rhs    -= ns->cc.ucell * ns->cc.dndx;
                        rhs    += higflow_electric_convective_ionic_term_central(ns, alphaeo, Pe);
                        rhs    += higflow_electric_divergence_ionic_term(ns, alphaeo, Pe);
                        break;
                        
                    case CELL_CUBISTA: // CUBISTA scheme
                        rhs += hig_flow_convective_ionic_cell_term_cubista(ns, nplus, ccenter, cdelta, dim, alphaeo, Pe, POSITIVE);
                        break;
                }
            }
            // Compute the final value step time
            real newnplus =  nplus + ns->par.dt * rhs;
            // Set property value  
            dp_set_value(ns->ed.eo.dpnplus_temp, clid, newnplus);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed ionic property
        dp_sync(ns->ed.eo.dpnplus_temp);
    }
}

//Explicit ionic transport equation
void higflow_explicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == PNP) {
        real alphaeo = ns->ed.eo.par.alpha;
        real Pe = ns->ed.eo.par.Pe;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        // Get the map for the domain property nminus
        mp_mapper *m = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *itt;
        for (itt = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(itt); higcit_nextcell(itt)) {
            // Get the cell
            hig_cell *c = higcit_getcell(itt);
            // Get the cell identifier
            int clid    = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            // real u[DIM], dnminusdx[DIM];
            // hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            real nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                // Set the computational cell in this particular direction
                higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
                // Compute the diffusive ionic term rhs
                rhs += higflow_diffusive_ionic_term(ns, Pe);
                // convective term
                switch (ns->ed.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
                        //hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                        rhs -= ns->cc.ucell * ns->cc.dndx;
                        rhs    -= higflow_electric_convective_ionic_term_central(ns, alphaeo, Pe);
                        rhs    -= higflow_electric_divergence_ionic_term(ns, alphaeo, Pe);
                        break;
                        
                    case CELL_CUBISTA: // CUBISTA scheme
                        rhs += hig_flow_convective_ionic_cell_term_cubista(ns, nminus, ccenter, cdelta, dim, alphaeo, Pe, NEGATIVE);
                        break;
                }
            }
            // Compute the final value step time
            real newnminus =  nminus + ns->par.dt * rhs;
            // Set ionic property value  
            dp_set_value(ns->ed.eo.dpnminus_temp, clid, newnminus);
        }
        // Destroy the iterator
        higcit_destroy(itt);
        // Sync the distributed ionic property
        dp_sync(ns->ed.eo.dpnminus_temp);
    }
}

//Semi-implicit ionic transport equation
void higflow_semi_implicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == PNP) {
        real alphaeo = ns->ed.eo.par.alpha;
        real Pe = ns->ed.eo.par.Pe;
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnplus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnplus); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            // real u[DIM], dnplusdx[DIM];
            // hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            real nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                // Set the computational cell in this particular direction
                higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
                // convective term
                switch (ns->ed.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
                        //hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                        rhs -= ns->cc.ucell * ns->cc.dndx;
                        rhs    += higflow_electric_convective_ionic_term_central(ns, alphaeo, Pe);
                        rhs    += higflow_electric_divergence_ionic_term(ns, alphaeo, Pe);
                        break;  
                    case CELL_CUBISTA: // CUBISTA scheme
                        rhs += hig_flow_convective_ionic_cell_term_cubista(ns, nplus, ccenter, cdelta, dim, alphaeo, Pe, POSITIVE);
                        break;
                }
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.eo.stnnplus);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.eo.stnnplus, rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->ed.eo.par.Pe*cdelta[dim2]*cdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.eo.stnnplus);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.eo.stnnplus);                
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnplus, ccenter, ccenter, alpha, ns->ed.eo.stnnplus);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnplus, ns->ed.eo.stnnplus);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.eo.stnnplus);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.eo.stnnplus);
            int cgid = psd_get_global_id(ns->ed.eo.psdEOnplus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnplus, cgid, stn_get_rhs(ns->ed.eo.stnnplus));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.eo.slvnplus, cgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.eo.slvnplus);
        // Solve the linear system
        slv_solve(ns->ed.eo.slvnplus);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.eo.dpnplus_temp, ns->ed.eo.slvnplus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnplus_temp);
    }
}

//Semi-implicit ionic transport equation
void higflow_semi_implicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == PNP) {
        real alphaeo = ns->ed.eo.par.alpha;
        real Pe = ns->ed.eo.par.Pe;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        // Get the map for the domain property nminus
        mp_mapper *m = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(m, hig_get_cid(c));
            // Get the center of the cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity at cell center 
            // real u[DIM], dnminusdx[DIM];
            // hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            real nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                // Set the computational cell in this particular direction
                higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
                // convective term
                switch (ns->ed.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
                        //hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                        rhs -= ns->cc.ucell * ns->cc.dndx;
                        rhs    -= higflow_electric_convective_ionic_term_central(ns, alphaeo, Pe);
                        rhs    -= higflow_electric_divergence_ionic_term(ns, alphaeo, Pe);
                        break;
                    case CELL_CUBISTA: // CUBISTA scheme
                        rhs += hig_flow_convective_ionic_cell_term_cubista(ns, nminus, ccenter, cdelta, dim, alphaeo, Pe, NEGATIVE);
                        break;
                }
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.eo.stnnminus);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.eo.stnnminus, rhs);
            // Calculate the point and weight of the stencil
            real alpha2 = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w   = -ns->par.dt/(ns->ed.eo.par.Pe*cdelta[dim2]*cdelta[dim2]);
                alpha2  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.eo.stnnminus);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.eo.stnnminus);
            }
            alpha2 = 1.0 + alpha2;
            // Get the stencil
            sd_get_stencil(sdnminus, ccenter, ccenter,alpha2, ns->ed.eo.stnnminus);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnminus, ns->ed.eo.stnnminus);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.eo.stnnminus);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.eo.stnnminus);
            int cgid = psd_get_global_id(ns->ed.eo.psdEOnminus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnminus, cgid, stn_get_rhs(ns->ed.eo.stnnminus));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->ed.eo.slvnminus, cgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Assemble the solver
        slv_assemble(ns->ed.eo.slvnminus);
        // Solve the linear system
        slv_solve(ns->ed.eo.slvnminus);
        //Load property from solver
        dp_slv_load_from_solver(ns->ed.eo.dpnminus_temp, ns->ed.eo.slvnminus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnminus_temp);
    }
}



// *******************************************************************
// Calculate convective ionic cell-term CUBISTA
// *******************************************************************
real hig_flow_convective_ionic_cell_term_cubista(higflow_solver *ns, real nc, Point ccenter, Point cdelta, int dim, real alphaeo, real Pe, charge_type charge_sign) {
    
    distributed_property *dpphi = ns->ed.eo.dpphi;
    distributed_property *dppsi = ns->ed.eo.dppsi;
    sim_domain *sdphi = ns->ed.eo.sdEOphi;
    sim_domain *sdpsi = ns->ed.eo.sdEOpsi;
    sim_stencil *stnphi = ns->ed.eo.stnphi;
    sim_stencil *stnpsi = ns->ed.eo.stnpsi;

    distributed_property *dpu = ns->dpu[dim];
    sim_facet_domain *sfdu = ns->sfdu[dim];
    sim_stencil *stn = ns->stn;

    distributed_property *dpn; 
    sim_domain *sdn;
    sim_stencil *stnn;
    real transport_factor;

    switch (charge_sign) {
        case POSITIVE:
            dpn = ns->ed.eo.dpnplus;
            sdn = ns->ed.eo.sdEOnplus;
            stnn = ns->ed.eo.stnnplus;
            transport_factor = alphaeo / Pe;
        break;
        case NEGATIVE:
            dpn = ns->ed.eo.dpnminus;
            sdn = ns->ed.eo.sdEOnminus;
            stnn = ns->ed.eo.stnnminus;
            transport_factor = -alphaeo/ Pe;
        break;
    }
    
    real vbar, vr, vl, nr, nrr, nl, nll, a, b, c, d, e, fi, conv1, conv2;
    real ur, ul, dphidxr, dphidxl, dpsidxr, dpsidxl, phir, phil, phic, psir, psil, psic;
    a = 1.7500;
    b = 0.3750;
    c = 0.7500;
    d = 0.1250;
    e = 0.2500;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;

    // Get the low, high, lowlow, highhigh component kernel at center cell
    nl = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 1.0, dpn, stnn, &incell_l);
    nr = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 1.0, dpn, stnn, &incell_r);
    nll = compute_center_p_left_22(sdn, ccenter, cdelta, dim, 2.0, dpn, stnn, &incell_ll);
    nrr = compute_center_p_right_22(sdn, ccenter, cdelta, dim, 2.0, dpn, stnn, &incell_rr);

    phil = compute_center_p_left(sdphi, ccenter, cdelta, dim, 1.0, dpphi, stnphi);
    phir = compute_center_p_right(sdphi, ccenter, cdelta, dim, 1.0, dpphi, stnphi);
    phic = compute_value_at_point(sdphi, ccenter, ccenter, 1.0, dpphi, stnphi);
    psil = compute_center_p_left(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, stnpsi);
    psir = compute_center_p_right(sdpsi, ccenter, cdelta, dim, 1.0, dppsi, stnpsi);
    psic = compute_value_at_point(sdpsi, ccenter, ccenter, 1.0, dppsi, stnpsi);

    dphidxl = compute_dpdxl_at_point(cdelta, dim, 1.0, phil, phic);
    dphidxr = compute_dpdxr_at_point(cdelta, dim, 1.0, phic, phir);
    dpsidxl = compute_dpdxl_at_point(cdelta, dim, 1.0, psil, psic);
    dpsidxr = compute_dpdxr_at_point(cdelta, dim, 1.0, psic, psir);

    ul = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
    ur = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
    
    vl = -ul + (dphidxl + dpsidxl) * transport_factor;
    vr = -ur + (dphidxr + dpsidxr) * transport_factor;

    //////////////////////// Upwind ////////////////////////////////////
    /*
    real philf = compute_center_p_left(sdphi, ccenter, cdelta, dim, 0.5, dpphi, stnphi);
    real phirf = compute_center_p_right(sdphi, ccenter, cdelta, dim, 0.5, dpphi, stnphi);
    real psilf = compute_center_p_left(sdpsi, ccenter, cdelta, dim, 0.5, dppsi, stnpsi);
    real psirf = compute_center_p_right(sdpsi, ccenter, cdelta, dim, 0.5, dppsi, stnpsi);

    real dphidx = compute_dpdx_at_point(cdelta, dim, 1.0, phil, phir);
    real dpsidx = compute_dpdx_at_point(cdelta, dim, 1.0, psil, psir);
    real v = (dphidx + dpsidx) * transport_factor;

    real dkdxl = compute_dpdxl_at_point(cdelta, dim, 1.0, nl, nc);
    real dkdxr = compute_dpdxr_at_point(cdelta, dim, 1.0, nc, nr);
    //real dkdx = compute_dpdx_at_point(cdelta, dim, 1.0, nl, nr);
    
    conv1 = (v+fabs(v))*dkdxl;
    conv2 = (v-fabs(v))*dkdxr;
    return 0.5*(conv1+conv2);
    */
    ///////////////////////////////////////////////////////////////////////


    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar = vr;
    if (vbar > 0.0) {
        if (FLT_EQ(nr, nl)) {
            conv1 = vbar * nc;
        }
        else {
            fi = (nc - nl) / (nr - nl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar * nc;
            }
            else {
                if (fi < b) {
                    if (incell_l == 1)                    conv1 = vbar * (a * nc - c * nl);
                    else                                  conv1 = vbar * nc;
                }
                if ((fi >= b) && (fi <= c)) {
                    if ((incell_l == 1) && (incell_r == 1)) conv1 = vbar * (c * nc + b * nr - d * nl);
                    else                                  conv1 = vbar * nc;
                }
                if (fi > c) {
                    if (incell_r == 1)                    conv1 = vbar * (e * nc + c * nr);
                    else                                  conv1 = vbar * nc;
                }

            }
        }
        //v1bar < 0.0
    }
    else {
        if ((incell_r == 1) && (incell_rr == 1)) {
            if (FLT_EQ(nc, nrr)) {
                conv1 = vbar * nr;
            }
            else {
                fi = (nr - nrr) / (nc - nrr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar * nr;
                }
                else {
                    if (fi < b)
                        conv1 = vbar * (a * nr - c * nrr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar * (c * nr + b * nc - d * nrr);
                    if (fi > c)
                        conv1 = vbar * (c * nc + e * nr);
                }
            }
        }
        else if ((incell_r == 1) && (incell_rr == 0)) { 
            if (FLT_EQ(nc, nrr)) {
                conv1 = vbar * nr;
            }
            else {
                fi = (nr - nrr) / (nc - nrr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar * nr;
                }
                else {
                    if (fi <= c)
                        conv1 = vbar * nr;
                    if (fi > c)
                        conv1 = vbar * (c * nc + e * nr);
                }
            }
        }
        else { //Return upwind value at boundary
            vbar = vr;
            if (vbar > 0.0) conv1 = vbar * nc;
            else                 conv1 = vbar * nc;
            vbar = vl;
            if (vbar > 0.0) conv2 = vbar * nl;
            else                 conv2 = vbar * nc;
            return ((conv1 - conv2) / cdelta[dim]);
        }

    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar = vl;
    if (vbar > 0.0) {
        if ((incell_l == 1) && (incell_ll == 1)) {
            if (FLT_EQ(nc, nll)) {
                conv2 = vbar * nl;
            }
            else {
                fi = (nl - nll) / (nc - nll);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv2 = vbar * nl;
                }
                else {
                    if (fi < b)
                        conv2 = vbar * (a * nl - c * nll);
                    if ((fi >= b) && (fi <= c))
                        conv2 = vbar * (b * nc + c * nl - d * nll);
                    if (fi > c)
                        conv2 = vbar * (c * nc + e * nl);
                }
            }
        }
        else if ((incell_l == 1) && (incell_ll == 0)) {
            if (FLT_EQ(nc, nll)) {
                conv2 = vbar * nl;
            }
            else {
                fi = (nl - nll) / (nc - nll);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv2 = vbar * nl;
                }
                else {
                    if (fi <= c)
                        conv2 = vbar * nl;
                    if (fi > c)
                        conv2 = vbar * (c * nc + e * nl);
                }
            }
        }
        else { //Return upwind value at boundary
            vbar = vr;
            if (vbar > 0.0) conv1 = vbar * nc;
            else                 conv1 = vbar * nr;
            vbar = vl;
            if (vbar > 0.0) conv2 = vbar * nc;
            else                 conv2 = vbar * nc;
            return ((conv1 - conv2) / cdelta[dim]);
        }
    }
    else {
        //v2bar < 0.0 
        if (FLT_EQ(nl, nr)) {
            conv2 = vbar * nc;
        }
        else {
            fi = (nc - nr) / (nl - nr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar * nc;
            }
            else {
                if (fi < b) {
                    if (incell_r == 1)                    conv2 = vbar * (a * nc - c * nr);
                    else                                  conv2 = vbar * nc;
                }
                if ((fi >= b) && (fi <= c)) {
                    if ((incell_l == 1) && (incell_r == 1)) conv2 = vbar * (c * nc + b * nl - d * nr);
                    else                                  conv2 = vbar * nc;
                }
                if (fi > c) {
                    if (incell_l == 1)                    conv2 = vbar * (e * nc + c * nl);
                    else                                  conv2 = vbar * nc;
                }
            }
        }
    }
    return ((conv1 - conv2) / cdelta[dim]);
}


// Electroosmotic induced potential psi 
real higflow_electroosmotic_psi(higflow_solver *ns) {
    real alphaeo = ns->ed.eo.par.alpha;
    real delta   = ns->ed.eo.par.delta;
    real nminus_temp, nplus_temp, psi;
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    // Get the map for the domain property nminus
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        int clid    = mp_lookup(m, hig_get_cid(c));
        // Reset the stencil
        stn_reset(ns->ed.eo.stnpsi);
        // Initialize rhs
        real rhs = 0.0;
        // Set the right side of stencil
        if (ns->ed.eo.contr.eo_model == PNP) {
            // Poisson-Nernst-Planck model 
            // Get the ionic concentration n- at center cell
            nplus_temp    = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus_temp, ns->ed.eo.stnnplus);
            nminus_temp   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus_temp, ns->ed.eo.stnnminus);
            rhs      = -delta*(nplus_temp - nminus_temp);
        } else if (ns->ed.eo.contr.eo_model == PB) {
            // Poisson-Boltzmann model 
            psi      = dp_get_value(ns->ed.eo.dppsi, clid);
            rhs      = 2.0*delta*sinh(alphaeo*psi) - 2.0*delta*alphaeo*psi*cosh(alphaeo*psi);
        }
        // Calculate the point and weight of the stencil
        stn_set_rhs(ns->ed.eo.stnpsi, rhs);
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            real perm = ns->ed.eo.get_permittivity(p, ns->par.t);
            // Stencil point update: right point
            p[dim] = ccenter[dim] + cdelta[dim];
            real permr = ns->ed.eo.get_permittivity(p, ns->par.t);
            real permrc = 0.5*(perm + permr);
            real wr = permrc/(cdelta[dim]*cdelta[dim]);
            sd_get_stencil(sdp, ccenter, p, wr, ns->ed.eo.stnpsi);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            real perml = ns->ed.eo.get_permittivity(p, ns->par.t);
            real permlc = 0.5*(perm + perml);
            real wl = permlc/(cdelta[dim]*cdelta[dim]);
            sd_get_stencil(sdp, ccenter, p, wl, ns->ed.eo.stnpsi);

            alpha -= (wr + wl);
        }
        switch (ns->ed.eo.contr.eo_model) {
        case PB:
            // Poisson-Boltzmann model 
            alpha -= 2.0*alphaeo*delta*cosh(alphaeo*psi);
            break;
        case PBDH:
            // Debye-Hückel model (solving poisson equation for psi) 
            alpha -= 2.0*alphaeo*delta;
            break;
        case PBDH_ANALYTIC:
            // Debye-Hückel model (using the analytic solution for psi)
            printf("???????This model should not solve poisson equation for potential psi!???????\n");
            exit(1);
            break;
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.eo.stnpsi);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOpsi, ns->ed.eo.stnpsi);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.eo.stnpsi);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.eo.stnpsi);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOpsi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvpsi, cgid, stn_get_rhs(ns->ed.eo.stnpsi));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->ed.eo.slvpsi, cgid, numelems, ids, vals);
    }
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->ed.eo.slvpsi);
    // Solve the linear system
    slv_solve(ns->ed.eo.slvpsi);
    // Set the solver solution in distributed property
    //dp_slv_load_from_solver(ns->ed.eo.dppsi, ns->ed.eo.slvpsi);
    real max_psi_res = 0.0;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int clid    = mp_lookup(m, hig_get_cid(c));
        int cgid    = psd_lid_to_gid(ns->ed.eo.psdEOpsi, clid);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the value of psi
        real psi_new = slv_get_xi(ns->ed.eo.slvpsi, cgid);
        real psi_old = dp_get_value(ns->ed.eo.dppsi, clid);
        UPDATE_RESIDUAL_BUFFER_CELL(ns, psi_old, psi_new, c, ccenter)
        real res = fabs(psi_new - psi_old);
        if(res > max_psi_res) max_psi_res = res;

        // Store psi
        dp_set_value(ns->ed.eo.dppsi, clid, psi_new);   
    }
    // Destroy the iterator
    higcit_destroy(it);

    UPDATE_RESIDUALS(ns, ns->residuals->psi)

    dp_sync(ns->ed.eo.dppsi);

    // get the global residual
    real max_psi_res_global = INFINITY;
    MPI_Allreduce(&max_psi_res, &max_psi_res_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max_psi_res_global;
}



void higflow_electroosmotic_solve_pb(higflow_solver *ns) {
    real max_psi_res = INFINITY, max_psi_res_global = INFINITY;
    real pb_tol = EPSMACH;
    int iter = 0, maxiter = 50;
    
    while (max_psi_res_global > pb_tol && iter < maxiter) {
        max_psi_res = higflow_electroosmotic_psi(ns);
        MPI_Allreduce(&max_psi_res, &max_psi_res_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        print0f("Poisson-Boltzmann model residual: %e\n", max_psi_res_global);
        iter++;
    }
    if (iter < maxiter) {
        print0f("Poisson-Boltzmann model converged after %d iterations\n", iter);
    }
    else {
        print0f("Poisson-Boltzmann model did not converge after %d iterations\n", iter);
        exit(1);
    }

}

// Electroosmotic applied potential phi 
void higflow_electroosmotic_phi(higflow_solver *ns) {
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOphi);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        // Get the delta of the cell
        Point cdelta;
        hig_get_delta(c, cdelta);
        // Reset the stencil
        stn_reset(ns->ed.eo.stnphi);
        // Set the right side of stencil
        stn_set_rhs(ns->ed.eo.stnphi, 0.0);
        // Calculate the point and weight of the stencil
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil weight update
            real w = 1.0/(cdelta[dim]*cdelta[dim]);
            alpha -= 2.0 * w;
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            // Stencil point update: right point
            p[dim] = ccenter[dim] + cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.eo.stnphi);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.eo.stnphi);
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.eo.stnphi);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOphi, ns->ed.eo.stnphi);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.eo.stnphi);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.eo.stnphi);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOphi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvphi, cgid, stn_get_rhs(ns->ed.eo.stnphi));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->ed.eo.slvphi, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->ed.eo.slvphi);
    // Solve the linear system
    slv_solve(ns->ed.eo.slvphi);
    // Set the solver solution in distributed property
    dp_slv_load_from_solver(ns->ed.eo.dpphi, ns->ed.eo.slvphi);
    dp_sync(ns->ed.eo.dpphi);
}

// *******************************************************************
// Electro-osmotic source term
// *******************************************************************
void higflow_calculate_electroosmotic_source_term( higflow_solver *ns) {
    real phil, phir, psil, psir, dphidx, dpsidx, npl, npr, nml, nmr, nplus, nminus, rho, psi, Feo;
    // Get the necessary parameters
    real alphaeo = ns->ed.eo.par.alpha;
    real delta   = ns->ed.eo.par.delta;
    // Get the local sub-domain
    sim_domain *sdpsi = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi  = psd_get_local_domain(ns->ed.eo.psdEOphi);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the derivative of applied potential
            phil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            phir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, phil, phir);
            // Get the derivative of induced potential
            psil        = compute_center_p_left(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            psir        = compute_center_p_right(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            dpsidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            
            switch (ns->ed.eo.contr.eo_model) {
                case PNP:
                    // Poisson-Nernst-Planck model 
                    // Get the ionic concentration n+ at center cell
                    npl     = compute_center_p_left(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
                    npr     = compute_center_p_right(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
                    nplus    = compute_value_at_mid_point(npl, npr);
                    // Get the ionic concentration n- at center cell
                    nml     = compute_center_p_left(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
                    nmr     = compute_center_p_right(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
                    nminus   = compute_value_at_mid_point(nml, nmr);
                    // Compute the electric density
                    rho = (nplus - nminus)*delta;
                    break;
                case PB:
                    // Poisson-Boltzmann model 
                    // Get the induced potential and its derivative
                    psi    = compute_value_at_mid_point(psil, psir);
                    // Compute the charge density 
                    rho    = -2.0*delta*sinh(alphaeo*psi);
                    break;
                case PBDH:
                    // Debye-Hückel model (solving poisson equation for psi) 
                    // Get the induced potential and its derivative
                    psi    = compute_value_at_mid_point(psil, psir);
                    // Compute the charge density 
                    rho    = -2.0*alphaeo*delta*psi;
                    break;
            }
            // Compute the electro-osmotic source term
            //Feo = -rho*(dphidx + dpsidx) ;
            Feo = -rho*(dphidx) ; // large normal psi gradients make the above formulation incompatible with pressure neumann conditions
            real normE2 = dphidx*dphidx;
            for (int dim2 = 0; dim2 < DIM; dim2++) {
                if(dim2 != dim) {
                    real phil_ = compute_center_p_left_2(sdphi, fcenter, fdelta, dim, dim2, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
                    real phir_ = compute_center_p_right_2(sdphi, fcenter, fdelta, dim, dim2, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
                    real dphidx_ = compute_dpdx_at_point(fdelta, dim2, 1.0, phil_, phir_);
                    normE2 += dphidx_*dphidx_;
                }
            }
            Point ccenterl, ccenterr;
            POINT_ASSIGN(ccenterl, fcenter); POINT_ASSIGN(ccenterr, fcenter);
            ccenterl[dim] -= 0.5*fdelta[dim]; ccenterr[dim] += 0.5*fdelta[dim];
            real perml = ns->ed.eo.get_permittivity(ccenterl, ns->par.t);
            real permr = ns->ed.eo.get_permittivity(ccenterr, ns->par.t);
            real dpermdx = compute_dpdx_at_point(fdelta, dim, 0.5, perml, permr);

            // Extra term arising due to (possibly non-uniform) permittivity gradient
            // Comes from Korteweg-Helmholtz force for linear dielectric isotropic media with instantaneous polarization response
            Feo += -0.5*normE2*dpermdx;
            
            // Get the electroosmotic extra source term defined by user
            Feo   += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_analytic_pbdh( higflow_solver *ns) {
    real eps    = 7.0e-10;
    real E      = -0.912;
    real H      = 1.0e-5;
    real zeta   = 0.025;
    real kappa  = 10.0e5;
    // Get the local sub-domain
    //sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdF[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdF[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            real y   = fcenter[1];
            // Set the induced potential
            real psi = cosh(kappa*H*y)/cosh(kappa*H);
            // Set the charge density 
            real rho = -eps*zeta*kappa*kappa*psi;
            // Set the electro-osmotic source term
            real Feo;
            if(dim == 0) Feo = rho*E;
            else         Feo = 0.0;
            //DEBUG_INSPECT(Feo, %lf);
            // Get the electroosmotic source term defined by user
            Feo     += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns) {
    real bcval, fracvol;
    // Facet iterator
    higcit_celliterator *it;
    // Local sub-domain
    sim_facet_domain *sfdFeo[DIM];
    // For each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local sub-domain
        sfdFeo[dim]      = psfd_get_local_domain(ns->ed.eo.psfdEOFeo[dim]);
        sim_domain *sd = sfdFeo[dim]->cdom;
        // Get the number of Dirichlet boundary
        int numbcs = sd_get_num_bcs(sd, DIRICHLET);
        // For each Dirichlet boundary
        for (int i = 0; i < numbcs; i++) {
            // Get the Dirichlet boundary
            sim_boundary *bc = sd_get_bc(sd, DIRICHLET, i);
            // Get the id defined by the user
            int userid       = sb_get_userid(bc);
            // Get the value type of the boundary condition
            bc_valuetype valuetype = sb_get_valuetype(bc);
            // Apply the boundary condition if time dependent
            if (valuetype == timedependent) {
                // Get the mapper
                mp_mapper *bm    = sb_get_mapper(bc);
                // For each cell of the boundary
                for(it = sb_get_celliterator(bc); !higcit_isfinished(it); higcit_nextcell(it)) {
                    // Get the cell
                    hig_cell *bcell = higcit_getcell(it);
                    // Get the cell center
                    Point bccenter;
                    hig_get_center(bcell, bccenter);
                    // Get the id of the cell
                    int bclid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to apply the boundary condition
                    real t   = ns->par.t + ns->par.dt;
                    // Get the velocity defined by the user
                    if(ns->contr.flowtype == MULTIPHASE) {
                        fracvol = compute_value_at_point(ns->ed.mult.sdmult, bccenter, bccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                        bcval = ns->ed.mult.eo.get_boundary_multiphase_electroosmotic_source_term(fracvol, userid, bccenter, dim, t);
                    } else
                        bcval = ns->ed.eo.get_boundary_electroosmotic_source_term(userid, bccenter, dim, t);
                    // Set the value
                    sb_set_value(bc, bclid, bcval);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }
}

// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp    = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
        // Get the local domain for facet cell
        for (int dim2 = 0; dim2 < DIM; dim2++) {
            sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
        }
        // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns, ns->ed.eo.par.Ex);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_diffusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Compute the intermediate velocity
            real ustar = ns->cc.ufacet + ns->par.dt * rhs;
            // Update the distributed property intermediate velocity
            dp_set_value(dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(dpustar[dim]);
        // Set the velocity at outflow
   //set_outflow(ns->psfdu[dim], ns->dpustar[dim], 20.0);
    }
}


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit Euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk2  = 0.5*(u + ustar);
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk2);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for third Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Calculate the order 2 Runge-Kutta method using the euler method
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = 0.75*u + 0.25*ustar;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpuaux[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpuaux[dim]);
    }
    // Calculate the order 2 Runge-Kutta method using the euler method
    higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpuaux, ns->dpustar);
    // Loop for each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
    // Calculate the third stage velocity by the explicit euler method
        for(fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the intermediate velocity
            real u     = dp_get_value(ns->dpu[dim], flid);
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the intermediate velocity
            real urk3  = u/3.0 + 2.0*ustar/3.0;
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpustar[dim], flid, urk3);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit Euler Method
// *******************************************************************
void higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
   // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns, ns->ed.eo.par.Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
       int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Semi-Implicit Crank-Nicolson Method
// *******************************************************************
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
   // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns, ns->ed.eo.par.Ex);
            // Diffusive term term contribution
            rhs += 0.5*higflow_diffusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -0.5*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
       int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// *******************************************************************
// Navier-Stokes Step for the Implicit BDF2 Method
// *******************************************************************
void higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(higflow_solver *ns) {
    // Firt stage of Tr-BDF2 method
    // Get the facet iterator
    higfit_facetiterator *fit;
    // Get the local domain for cell
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Get the local domain for facet cell
    for (int dim2 = 0; dim2 < DIM; dim2++) {
        sfdu[dim2] = psfd_get_local_domain(ns->psfdu[dim2]);
    }
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
   // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns, ns->ed.eo.par.Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= 0.5*ns->par.dt;
            // Difusive term contribution
            rhs += 0.25*ns->par.dt*higflow_diffusive_term(ns, fdelta);
            // Velocity term contribution
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = -0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -= 2.0*w; //divide po 4 para usar regra trapezio em t(n+1/2)
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter, alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
       int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real uaux = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpuaux[dim], flid, uaux);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpuaux[dim]);
    }
    // Second Stage of Tr-BDF2
    // Looping for the velocity
    for (int dim = 0; dim < DIM; dim++) {
   // Get the map of domain
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
       int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Set the computational cell
            higflow_computational_cell_electroosmotic(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            rhs += higflow_electroosmotic_source_term(ns, ns->ed.eo.par.Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Tensor term contribution
            rhs += higflow_tensor_term(ns);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt/3.0;
            rhs += (4.0*uaux - ns->cc.ufacet)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w  = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha  -=  2.0*w ;
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, w, ns->stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sfd_get_stencil(sfdu[dim], fcenter, fcenter,alpha, ns->stn);
            // Get the index of the stencil
            int *ids   = psfd_stn_get_gids(ns->psfdu[dim], ns->stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->stn);
       int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Set the right side of solver linear system
            slv_set_bi(ns->slvu[dim], fgid, stn_get_rhs(ns->stn));
            // Set the line of matrix of the solver linear system
            slv_set_Ai(ns->slvu[dim], fgid, numelems, ids, vals);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Assemble the solver
        slv_assemble(ns->slvu[dim]);
        // Solve the linear system
        slv_solve(ns->slvu[dim]);
        // Get the solution of linear system
        //Vec *vecu = slv_get_solution_vec(ns->slvu[dim]);
        // Gets the values of the solution
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet cell identifier
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            int fgid = psfd_lid_to_gid(ns->psfdu[dim], flid);
            // Get the value of ustar
            real ustar = slv_get_xi(ns->slvu[dim], fgid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

void print_minmax_properties(higflow_solver *ns) {

    real u_min[DIM]={INFINITY}, u_max[DIM]={-INFINITY}, p_min=INFINITY, p_max=-INFINITY;
    real phi_min=INFINITY, phi_max=-INFINITY, psi_min=INFINITY, psi_max=-INFINITY;
    real nplus_min=INFINITY, nplus_max=-INFINITY, nminus_min=INFINITY, nminus_max=-INFINITY;
    real u[DIM], p, phi, psi, nplus, nminus;

    real ustar_min[DIM]={INFINITY}, ustar_max[DIM]={-INFINITY}, deltap_min=INFINITY, deltap_max=-INFINITY;
    real ustar[DIM], deltap;

    real Feo_min[DIM] = {INFINITY}, Feo_max[DIM] = {-INFINITY};
    real Feo[DIM];

    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    //mapper for domain
    mp_mapper *m = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        //get cell and its id
        hig_cell *c = higcit_getcell(it);
        int clid    = mp_lookup(m, hig_get_cid(c));
        Point ccenter;
        hig_get_center(c, ccenter);
        //get distributed properties and update min and max
        for(int dim=0; dim<DIM; dim++){
            // u[dim] = dp_get_value(ns->dpu[dim], clid);
            // if(u[dim]>u_max[dim]) u_max[dim] = u[dim];
            // if(u[dim]<u_min[dim]) u_min[dim] = u[dim];
            // ustar[dim] = dp_get_value(ns->dpustar[dim], clid);
            // if(ustar[dim]>ustar_max[dim]) ustar_max[dim] = ustar[dim];
            // if(ustar[dim]<ustar_min[dim]) ustar_min[dim] = ustar[dim];
            Feo[dim] = compute_facet_value_at_point(ns->ed.eo.sfdEOFeo[dim], ccenter, ccenter, 1.0, ns->ed.eo.dpFeo[dim], ns->ed.eo.stnpsi);
            if(Feo[dim]>Feo_max[dim]) Feo_max[dim] = Feo[dim];
            if(Feo[dim]<Feo_min[dim]) Feo_min[dim] = Feo[dim];
        }
        // p = dp_get_value(ns->dpp, clid);
        // if(p>p_max) p_max = p;
        // if(p<p_min) p_min = p;
        phi = compute_value_at_point(ns->ed.eo.sdEOphi, ccenter, ccenter, 1.0, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
        if(phi>phi_max) phi_max = phi;
        if(phi<phi_min) phi_min = phi;
        psi = dp_get_value(ns->ed.eo.dppsi, clid);
        if(psi>psi_max) psi_max = psi;
        if(psi<psi_min) psi_min = psi;
        nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
        if(nplus>nplus_max) nplus_max = nplus;
        if(nplus<nplus_min) nplus_min = nplus;
        nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
        if(nminus>nminus_max) nminus_max = nminus;
        if(nminus<nminus_min) nminus_min = nminus;
        // deltap = dp_get_value(ns->ddeltap, clid);
        // if(deltap>deltap_max) deltap_max = deltap;
        // if(deltap<deltap_min) deltap_min = deltap;
    }
    higcit_destroy(it);

    // printf("####### MINMAX properties - SUBSTEP%5.3lf ########\n", substep_id);
    // for(int dim=0; dim<DIM; dim++){
    //     printf("ustar%dmin=%15.10lf | ustar%dmax=%15.10lf | ", dim, ustar_min[dim], dim, ustar_max[dim]);
    // }
    // printf("deltapmin=%15.10lf | deltapmax=%15.10lf\n", deltap_min, deltap_max);
    // for(int dim=0; dim<DIM; dim++){
    //     printf("u%dmin=%15.10lf | u%dmax=%15.10lf | ", dim, u_min[dim], dim, u_max[dim]);

    // }
    // printf("pmin=%15.10lf | pmax=%15.10lf\n", p_min, p_max);
    // printf("phimin=%15.10lf | phimax=%15.10lf | psimin=%15.10lf | psimax=%15.10lf\n",
    //     phi_min, phi_max, psi_min, psi_max);
    // printf("nplusmin=%15.10lf | nplusmax=%15.10lf | nminusmin=%15.10lf | nminusmax=%15.10lf\n",
    //     nplus_min, nplus_max, nminus_min, nminus_max);
    // for(int dim=0; dim<DIM; dim++){
    //     printf("Feo%dmin=%15.10lf | Feo%dmax=%15.10lf |", dim, Feo_min[dim], dim, Feo_max[dim]);
    // }
    //printf("\n####################################\n");

    real phi_min_global, phi_max_global, psi_min_global, psi_max_global;
    real nplus_min_global, nplus_max_global, nminus_min_global, nminus_max_global;
    real Feo_min_global[DIM], Feo_max_global[DIM];
    MPI_Allreduce(&phi_min, &phi_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&phi_max, &phi_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&psi_min, &psi_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&psi_max, &psi_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nplus_min, &nplus_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&nplus_max, &nplus_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nminus_min, &nminus_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&nminus_max, &nminus_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&Feo_min, &Feo_min_global, DIM, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&Feo_max, &Feo_max_global, DIM, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // print psi, nplus and Feo
    print0f("===> \u03A8min = %15.10lf <===> \u03A8max = %15.10lf <===\n", psi_min_global, psi_max_global);
    if(ns->ed.eo.contr.eo_model == PNP)
        print0f("===> n+min = %15.10lf <===> n+max = %15.10lf <===\n", nplus_min_global, nplus_max_global);
    for(int dim=0; dim<DIM; dim++)
        print0f("===> %d: Feomin = %15.10lf <===> Feomax = %15.10lf <===\n", dim, Feo_min_global[dim], Feo_max_global[dim]);

}

// One step of the Navier-Stokes the projection method
void higflow_solver_step_electroosmotic(higflow_solver *ns) {
    // Boundary conditions and source terms
    higflow_boundary_condition_for_velocity(ns);
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    higflow_boundary_condition_for_pressure(ns);
    higflow_calculate_source_term(ns);
    higflow_calculate_facet_source_term(ns);

    // Calculate the electro-osmotic terms
    higflow_boundary_condition_for_phi(ns);
    higflow_boundary_condition_for_psi(ns);
    switch (ns->ed.eo.contr.eo_model) {
    case PNP: // Poisson-Nernst-Planck model
        for(int k=0; k<2; k++) {
            higflow_boundary_condition_for_electroosmotic_nplus(ns);
            higflow_boundary_condition_for_electroosmotic_nminus(ns);
            if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
            switch (ns->ed.eo.contr.tempdiscrtype) {
                case EXPLICIT_EULER:
                    higflow_explicit_euler_ionic_transport_equation_nplus(ns);
                    higflow_explicit_euler_ionic_transport_equation_nminus(ns);
                    break;
                case SEMI_IMPLICIT_EULER:
                    higflow_semi_implicit_euler_ionic_transport_equation_nplus(ns);
                    higflow_semi_implicit_euler_ionic_transport_equation_nminus(ns);
                    break;
            }
            real max_psi_res = higflow_electroosmotic_psi(ns);
            if(k>0) print0f("=+=+=+ psi residual in inner iteration %d: %15.10lf =+=+=+\n", k+1, max_psi_res);
        }
        dp_copy_values(ns->ed.eo.dpnplus, ns->ed.eo.dpnplus_temp);
        dp_copy_values(ns->ed.eo.dpnminus, ns->ed.eo.dpnminus_temp);
        higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PB: // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
        higflow_electroosmotic_solve_pb(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true
                                || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PBDH: // Debye-Hückel model (solving poisson equation for psi) 
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_electroosmotic_psi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_phibc_timedependent == true
                                || ns->ed.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_electroosmotic_source_term(ns);
        break;
    case PBDH_ANALYTIC: // Debye-Hückel model (using the analytic solution for psi)
        higflow_calculate_electroosmotic_source_term_analytic_pbdh(ns);
        break;
    } 

    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
    case EXPLICIT_EULER:
        higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
        break;
    case EXPLICIT_RK2:
        higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(ns);
        break;
    case EXPLICIT_RK3:
        higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_CN:
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_BDF2:
        higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(ns);
        break;
    }

    // Projection
    higflow_pressure(ns);
    higflow_final_velocity(ns);
    higflow_final_pressure(ns);
}


