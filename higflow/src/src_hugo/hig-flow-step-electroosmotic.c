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
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnplus);
    // Get the number of Dirichlet boundary
    int numbcs = sd_get_num_bcs(sdp, NEUMANN);
    // For each Dirichlet boundary
    for (int i = 0; i < numbcs; i++) {
        // Get the Dirichlet boundary
        sim_boundary *bc = sd_get_bc(sdp, NEUMANN, i);
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
                int bcgid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get the pressure defined by the user
                real val = ns->ed.eo.get_boundary_electroosmotic_nplus(userid, bccenter, t);
                // Set the value 
                sb_set_value(bc, bcgid, val);
             }
             // Destroy the iterator
             higcit_destroy(it);
         }
     }
}

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_nminus(higflow_solver *ns) {
    // cell iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOnminus);
    // Get the number of Dirichlet boundary
    int numbcs = sd_get_num_bcs(sdp, NEUMANN);
    // For each Dirichlet boundary
    for (int i = 0; i < numbcs; i++) {
        // Get the Dirichlet boundary
        sim_boundary *bc = sd_get_bc(sdp, NEUMANN, i);
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
                int bcgid = mp_lookup(bm, hig_get_cid(bcell));
                // Set the time to get the pressure
                real t = ns->par.t + ns->par.dt;
                // Get the pressure defined by the user
                real val = ns->ed.eo.get_boundary_electroosmotic_nminus(userid, bccenter, t);
                // Set the value 
                sb_set_value(bc, bcgid, val);
             }
             // Destroy the iterator
             higcit_destroy(it);
         }
     }
}

// *******************************************************************
// Ionic transport equation 
// *******************************************************************
void higflow_explicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        real alphaeo = ns->ed.eo.par.alpha;
        // Get the cosntants
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
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
            real u[DIM], dnplusdx[DIM], nplus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnplusdx[dim];
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnplus, ns->ed.eo.sdEOnplus, ns->ed.stn, nplus, ccenter, cdelta, dim);
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            // Compute the final value step time
            real newnplus =  nplus + ns->par.dt * rhs;
            // Set property value  
            dp_set_value(ns->ed.eo.dpnplus, clid, newnplus);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.eo.dpnplus);
    }
}

void higflow_explicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        real alphaeo = ns->ed.eo.par.alpha;
        // Get the cosntants
        real tol   = ns->ed.ve.par.kernel_tol;
        real small = 1.0e-14;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
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
            real u[DIM], dnminusdx[DIM], nminus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnminusdx[dim];
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    //Compute convective ionic term CUBISTA in rhs
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnminus, ns->ed.eo.sdEOnminus, ns->ed.stn, nminus, ccenter, cdelta, dim);
                        // Compute the diffusive ionic term rhs
                        rhs    += higflow_diffusive_ionic_term(ns);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            // Compute the final value step time
            real newnminus =  nminus + ns->par.dt * rhs;
            // Set ionic property value  
            dp_set_value(ns->ed.eo.dpnminus, clid, newnminus);
        }
        // Destroy the iterator
        higcit_destroy(itt);
        // Sync the ditributed ionic property
        dp_sync(ns->ed.eo.dpnminus);
    }
}

//Semi-implicit ionic transport equation
void higflow_implicit_euler_ionic_transport_equation_nplus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        // Get the local sub-domain for the cells
        sim_domain *sdnplus  = psd_get_local_domain(ns->ed.eo.psdEOnplus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
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
            real u[DIM], dnplusdx[DIM], nplus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nplus value at point 
            nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.stn);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nplus_at_center_cell(ns, ccenter, cdelta, nplus, dnplusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnplusdx[dim];
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnplus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnplus, ns->ed.eo.sdEOnplus, ns->ed.stn, nplus, ccenter, cdelta, dim);
                        // Compute the potential ionic term rhs
                        rhs    += higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
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
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, w, ns->ed.stn);
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnplus, ccenter, ccenter,alpha, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnplus, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	    int cgid = psd_get_global_id(ns->ed.eo.psdEOnplus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnplus, cgid, stn_get_rhs(ns->ed.stn));
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
        dp_slv_load_from_solver(ns->ed.eo.dpnplus, ns->ed.eo.slvnplus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnplus);
    }
}

//Semi-implicit ionic transport equation
void higflow_implicit_euler_ionic_transport_equation_nminus(higflow_solver *ns) {
    if (ns->ed.eo.contr.eo_model == 0) {
        // Get the local sub-domain for the cells
        sim_domain *sdnminus = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        sim_domain *sdpsi    = psd_get_local_domain(ns->ed.eo.psdEOpsi);
        sim_domain *sdphi    = psd_get_local_domain(ns->ed.eo.psdEOphi);
        // Get the local sub-domain for the facets
        sim_facet_domain *sfdu[DIM];
        for(int i = 0; i < DIM; i++) {
            sfdu[i] = psfd_get_local_domain(ns->psfdu[i]);
        }
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
            real u[DIM], dnminusdx[DIM], nminus;
            hig_flow_velocity_at_center_cell(ns, ccenter, cdelta, u);
            // compute nminus value at point 
            nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.stn);
            // Right hand side equation
            real rhs = 0.0;
            switch (ns->ed.eo.contr.convecdiscrtype) {
                // Central scheme
                case 0: 
                    // Ionic derivative at cell center
                    hig_flow_derivative_nminus_at_center_cell(ns, ccenter, cdelta, nminus, dnminusdx);
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective inonic term in rhs
                        rhs    -= u[dim]*dnminusdx[dim];
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
                // CUBISTA scheme
                case 1: 
                    //Compute convective ionic term CUBISTA in rhs
                    for (int dim = 0; dim < DIM; dim++) {
                        // Set the computational cell 
                        higflow_computational_cell_electroosmotic_ionic(ns, sdnminus, sdpsi, sdphi, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.dppsi, ns->ed.eo.dpphi);
                        //Compute convective ionic term CUBISTA in rhs
                        rhs    -= higflow_convective_ionic_term_cubista(ns, ns->dpu[dim], ns->ed.eo.dpnminus, ns->ed.eo.sdEOnminus, ns->ed.stn, nminus, ccenter, cdelta, dim);
                        // Compute the potential ionic term rhs
                        rhs    -= higflow_potential_ionic_term(ns);
                    }
                    break;
            }
            rhs *= ns->par.dt;
            rhs += ns->cc.ncell;
            // Reset the stencil
            stn_reset(ns->ed.stn);
            // Set the right side of stencil
            stn_set_rhs(ns->ed.stn,rhs);
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
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.stn);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, w, ns->ed.stn);
            }
            alpha2 = 1.0 + alpha2;
            // Get the stencil
            sd_get_stencil(sdnminus, ccenter, ccenter,alpha2, ns->ed.stn);
            // Get the index of the stencil
            int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOnminus, ns->ed.stn);
            // Get the value of the stencil
            real *vals = stn_get_vals(ns->ed.stn);
            // Get the number of elements of the stencil
            int numelems = stn_get_numelems(ns->ed.stn);
	    int cgid = psd_get_global_id(ns->ed.eo.psdEOnminus, c);
            // Set the right side of solver linear system
            slv_set_bi(ns->ed.eo.slvnminus, cgid, stn_get_rhs(ns->ed.stn));
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
        dp_slv_load_from_solver(ns->ed.eo.dpnminus, ns->ed.eo.slvnminus);
        // Syncing the distributed property
        dp_sync(ns->ed.eo.dpnminus);
    }
}

//Calculate the velocity at center cell
void hig_flow_velocity_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real u[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        // Verity if is in facet
        int infacet;
        // Get the velocity in the left facet center
        real ul = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Get the velocity in the right facet center
        real ur = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
        // Setting the velocity at cell center
        u[dim]  = 0.5*(ul + ur);
    }
}

// *******************************************************************
// Calculate convective ionic term CUBISTA
// *******************************************************************
real higflow_convective_ionic_term_cubista(higflow_solver *ns, distributed_property *dpu, distributed_property *dpn, sim_domain *sdp, sim_stencil *stn, real n, Point ccenter, Point cdelta, int dim) {
    real  vbar[DIM], dKdx[dim], kr, krr, kl, kll, kc, a, b, c, d, e, tol, fi, conv1,conv2;
    a     = 1.7500;
    b     = 0.3750;
    c     = 0.7500;
    d     = 0.1250;
    e     = 0.2500;
    tol   = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;
    // Get the kernel at center cell
    kc  = n;
    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl  = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 1.0, dpn, stn, &incell_l); 
    kr  = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 1.0, dpn, stn, &incell_r); 
    kll = compute_center_p_left_22(sdp, ccenter, cdelta, dim, 2.0, dpn, stn, &incell_ll);
    krr = compute_center_p_right_22(sdp, ccenter, cdelta, dim, 2.0, dpn, stn, &incell_rr);
    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (fabs(kr - kl) <= tol){
            conv1 = vbar[dim]*kc;
        }else {
            fi = (kc - kl)/(kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar[dim]*kc;
            }else {
                if (fi < b){ 
                    conv1 = vbar[dim]*(a*kc - c*kl);
                }
	        if ((fi >= b) && (fi <= c)){
                    conv1 = vbar[dim]*(c*kc + b*kr -d*kl);
                }
	        if (fi > c){ 
                    conv1 = vbar[dim]*(e*kc + c*kr);
                }
                    
            }    
        }
    //v1bar < 0.0
    }else {
        if (incell_r == 1){
            if (fabs(kc - krr) <= tol){
                conv1 = vbar[dim]*kr;
            }else {
                fi = (kr- krr)/(kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar[dim]*kr;
                }else {
		    if (fi < b) 
                        conv1 = vbar[dim]*(a*kr - c*krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar[dim]*(c*kr + b*kc -d*krr);
	            if (fi > c) 
                        conv1 = vbar[dim]*(c*kc + e*kr);
                }
            }
        //Return upwind value at boundary
        }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        }
        
    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
    if (vbar[dim] > 0.0){
        if (incell_l == 1){
            if (fabs(kc-kll) <= tol) {
	        conv2 = vbar[dim]*kl;
            }else {
	        fi = (kl - kll)/(kc - kll);
	        if ((fi <= 0.0) || (fi >= 1.0)) {
	            conv2 = vbar[dim]*kl;
	        }else {
	            if (fi < b)
	                conv2 = vbar[dim]*(a*kl - c*kll);
	            if ((fi >= b) && (fi <= c))
	                conv2 = vbar[dim]*(b*kc + c*kl - d*kll);
	            if (fi > c)  
	                conv2 = vbar[dim]*(c*kc + e*kl);
	        }
	    }
       }else {
                vbar[dim] = compute_facet_u_right(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv1 = vbar[dim]*kc;
                else                 conv1 = vbar[dim]*kr;
                vbar[dim] = compute_facet_u_left(ns->sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpu[dim], ns->stn, &infacet);
                if (vbar[dim] > 0.0) conv2 = vbar[dim]*kl;
                else                 conv2 = vbar[dim]*kc;
                return ((conv1 - conv2)/cdelta[dim]); 
        } 
    }else {
    //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar[dim]*kc;
        }else {
            fi = (kc - kr)/(kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar[dim]*kc;
            }else {
	        if (fi < b){
                    conv2 = vbar[dim]*(a*kc - c*kr);
                }
	        if ((fi >= b) && (fi <= c)){
                    conv2 = vbar[dim]*(c*kc + b*kl -d*kr);
                }
	        if (fi > c){ 
                    conv2 = vbar[dim]*(e*kc + c*kl);
                }
	    }
        }
    }
    return ((conv1-conv2)/cdelta[dim]);
}

// Get the derivative of nminus 
void hig_flow_derivative_nminus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.eo.sdEOnminus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnminus, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.eo.sdEOnminus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnminus, ns->ed.stn, &incell_left);
        // Compute the concentration derivative
           dndx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, nleft, nright);
    }
}

// Get the derivative of nplus 
void hig_flow_derivative_nplus_at_center_cell (higflow_solver *ns, Point ccenter, Point cdelta, real ncenter, real dndx[DIM]) {
    for (int dim = 0; dim < DIM; dim++) {
        int incell_left, incell_right;
        // Get the concentration in the left cell
        real nleft = compute_center_p_left_22(ns->ed.eo.sdEOnplus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnplus, ns->ed.stn, &incell_left);
        // Get the concentration in the right cell
        real nright = compute_center_p_right_22(ns->ed.eo.sdEOnplus, ccenter, cdelta, dim, 1.0, ns->ed.eo.dpnplus, ns->ed.stn, &incell_left);
        // Compute the concentration derivative
           dndx[dim] = compute_dpdx_at_point(cdelta, dim, 1.0, nleft, nright);
    }
}

// Electroosmotic induced potential psi 
void higflow_electroosmotic_psi(higflow_solver *ns) {
    real alphaeo = ns->ed.eo.par.alpha;
    real delta   = ns->ed.eo.par.delta;
    real nminus, nplus;
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
        stn_reset(ns->ed.stn);
        // Initialize rhs
        real rhs = 0.0;
        // Set the right side of stencil
        if (ns->ed.eo.contr.eo_model == 0) {
            // Poisson-Nernst-Planck model 
            // Get the ionic concentration n- at center cell
            nplus    = dp_get_value(ns->ed.eo.dpnplus, clid);
            nminus   = dp_get_value(ns->ed.eo.dpnminus, clid);
            rhs      = delta*(nminus - nplus);
        } 
        // Calculate the point and weight of the stencil
        stn_set_rhs(ns->ed.stn, rhs);
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
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
        }
        switch (ns->ed.eo.contr.eo_model) {
            case 1:
               // Poisson-Boltzmann model 
               printf("xxxxxx Please, you must implement this model xxxxxxx\n");
               exit(1);
               break;
            case 2:
               // Debye-Hückel model (solving poisson equation for psi) 
               alpha -= 2.0*alphaeo*delta;
               break;
            case 3:
               // Debye-Hückel model (using the analytic solution for psi)
               printf("???????This model should not solve poisson equation for potential psi!???????\n");
               exit(1);
               break;
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOpsi, ns->ed.stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOpsi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvpsi, cgid, stn_get_rhs(ns->ed.stn));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->ed.eo.slvpsi, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->ed.eo.slvpsi);
    // Solve the linear system
    slv_solve(ns->ed.eo.slvpsi);
    // Set the solver solution in distributed property
    dp_slv_load_from_solver(ns->ed.eo.dppsi, ns->ed.eo.slvpsi);
    dp_sync(ns->ed.eo.dppsi);
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
        stn_reset(ns->ed.stn);
        // Set the right side of stencil
        stn_set_rhs(ns->ed.stn, 0.0);
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
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->ed.stn);
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->ed.stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->ed.eo.psdEOphi, ns->ed.stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->ed.stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->ed.stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->ed.eo.psdEOphi, c);
        // Set the right side of solver linear system
        slv_set_bi(ns->ed.eo.slvphi, cgid, stn_get_rhs(ns->ed.stn));
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
void higflow_calculate_electroosmotic_source_term_pnp( higflow_solver *ns) {
    real psil, psir, nplus, nminus;
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
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Set the derivative of applied potential
            real phi    = compute_value_at_mid_point(psil, psir);
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Set the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Get the ionic concentration n+ in the left cell
            psil     = compute_center_p_left(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.stn);
            // Get the ionic concentration n+ in the right cell
            psir     = compute_center_p_right(ns->ed.eo.sdEOnplus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnplus, ns->ed.stn);
            // Get the ionic concentration n+ at center cell
            nplus    = compute_value_at_mid_point(psil, psir);
            // Get the ionic concentration n- in the left cell
            psil     = compute_center_p_left(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.stn);
            // Get the ionic concentration n- in the right cell
            psir     = compute_center_p_right(ns->ed.eo.sdEOnminus, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpnminus, ns->ed.stn);
            // Get the ionic concentration n- at center cell
            nminus   = compute_value_at_mid_point(psil, psir);
            // Compute the electrical density
            real rho = (nplus - nminus)*delta;   
            // Compute the electro-osmotic source term
            real Feo;
            Feo   = -rho*(dphidx + dpsidy); 
            // Get the electroosmotic extra source term defined by user
            Feo   += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_pb( higflow_solver *ns) {
    real psil, psir;
    real alphaeo= ns->ed.eo.par.alpha;
    real delta  = ns->ed.eo.par.delta;
    real eps    = 7.0e-10;
    real zeta   = 0.025;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi = psd_get_local_domain(ns->ed.eo.psdEOphi);
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
            real y      = fcenter[1];
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Set the derivative of applied potential
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
     //       DEBUG_INSPECT(dphidx,%lf);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Set the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Set the charge density 
            real rho    = -eps*zeta*2.0*alphaeo*delta*psi;
            // Set the electro-osmotic source term
            real Feo;
            if (dim == 0){
                  Feo   = -rho*dphidx;
            }else Feo   = -rho*dpsidy;
            //DEBUG_INSPECT(Feo, %lf);
            // Get the electroosmotic source term defined by user
            Feo        += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

void higflow_calculate_electroosmotic_source_term_pbdh( higflow_solver *ns) {
    real psil, psir;
    real alphaeo= ns->ed.eo.par.alpha;
    real delta  = ns->ed.eo.par.delta;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->ed.eo.psdEOpsi);
    sim_domain *sdphi = psd_get_local_domain(ns->ed.eo.psdEOphi);
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
            real y      = fcenter[1];
            // Get the applied potential in the left cell
            psil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Get the applied potential in the right cell
            psir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.stn);
            // Set the derivative of applied potential
            real dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
  //          DEBUG_INSPECT(dphidx,%lf);
            // Get the induced potential in the left cell
            psil        = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Get the induced potential in the right cell
            psir        = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.stn);
            // Set the induced potential
            real psi    = compute_value_at_mid_point(psil, psir);
            // Set the derivative of induced potential
            real dpsidy = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            // Set the charge density where  "epsez = epsilon/ez" 
            real rho   = -2.0*alphaeo*delta*psi;
            // Set the electro-osmotic source term
            real Feo;
            Feo   = -rho*(dphidx + dpsidy) ; 
            // Get the electroosmotic source term defined by user
            Feo  += ns->ed.eo.get_electroosmotic_source_term(fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the ditributed velocity property
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
        // Sync the ditributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}

// Apply the boundary condition for source term 
void higflow_boundary_condition_for_electroosmotic_source_term(higflow_solver *ns) {
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
                    int bcgid = mp_lookup(bm, hig_get_cid(bcell));
                    // Set the time to apply the boundary condition
                    real t   = ns->par.t + ns->par.dt;
                    // Get the velocity defined by the user
                    real val = ns->ed.eo.get_boundary_electroosmotic_source_term(userid, bccenter, dim, t);
                    // Set the value
                    sb_set_value(bc, bcgid, val);
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
            rhs += higflow_electroosmotic_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Compute the intermediate velocity
            real ustar = ns->cc.ucell + ns->par.dt * rhs;
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
        // Sync the ditributed velocity property
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
        // Sync the ditributed velocity property
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
        // Sync the ditributed velocity property
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
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms by delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
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
            rhs += higflow_electroosmotic_source_term(ns);
            // Diffusive term term contribution
            rhs += 0.5*higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt;
            // Velocity term contribution
            rhs += ns->cc.ucell;
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
void higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 0.5*ns->par.dt;
            // Difusive term contribution
            rhs += 0.25*ns->par.dt*higflow_difusive_term(ns, fdelta);
            // Velocity term contribution
            rhs += ns->cc.ucell;
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
            rhs += higflow_electroosmotic_source_term(ns);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= ns->par.dt/3.0;
            rhs += (4.0*uaux - ns->cc.ucell)/3.0;
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
            real ustar = slv_get_xi(ns->slvu[dim], flid);
            // Set the value of ustar
            dp_set_value(ns->dpustar[dim], flid, ustar);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Syncing the intermediate velocity
        dp_sync(ns->dpustar[dim]);
    }
}

// One step of the Navier-Stokes the projection method
void higflow_solver_step_electroosmotic(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Boundary condition for n+ 
//    higflow_boundary_condition_for_electroosmotic_nplus(ns);
    // Boundary condition for n- 
//    higflow_boundary_condition_for_electroosmotic_nminus(ns);
    // Set outflow for velocity
    higflow_outflow_u_step(ns);
    // Calculate the electro-osmotic source term
    switch (ns->ed.eo.contr.eo_model) {
        case 0:
           // Poisson-Nernst-Planck model 
 //          higflow_electroosmotic_phi(ns);  // está sendo calculado no ¨ns-example-2d.c¨
           higflow_implicit_euler_ionic_transport_equation_nplus(ns);
           higflow_implicit_euler_ionic_transport_equation_nminus(ns);
    //       higflow_explicit_euler_ionic_transport_equation_nplus(ns);
    //       higflow_explicit_euler_ionic_transport_equation_nminus(ns);
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pnp(ns); 
           break;
        case 1:
           // Poisson-Boltzmann model 
           //higflow_electroosmotic_phi(ns);
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pb(ns); 
           break;
        case 2:
           // Debye-Hückel model (solving poisson equation for psi) 
 //          higflow_electroosmotic_phi(ns);  // está sendo calculado no ¨ns-example-2d.c¨
           higflow_electroosmotic_psi(ns);
           higflow_calculate_electroosmotic_source_term_pbdh(ns); 
           break;
        case 3:
           // Debye-Hückel model (using the analytic solution for psi)
           higflow_calculate_electroosmotic_source_term_analytic_pbdh(ns); 
           break;
    }
    // Calculate the intermediated velocity
    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
        case 0:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
           break;
        case 1:
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity_electroosmotic(ns);
           break;
        case 2:
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity_electroosmotic(ns);
           break;
        case 3:
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity_electroosmotic(ns);
           break;
        case 4:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity_electroosmotic(ns);
           break;
        case 5:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity_electroosmotic(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Set outflow for ustar velocity 
    higflow_outflow_ustar_step(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Set outflow for u velocity 
    higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
}


