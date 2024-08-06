// *******************************************************************
//  HiG-Flow Solver Step Multiphase Electroosmotic - version 07/05/2024
// *******************************************************************

#include "hig-flow-step-multiphase-electroosmotic.h"


real higflow_interp_alpha_multiphase_electroosmotic(real alpha0, real alpha1, real fracvol) {
    return (1.0 - fracvol) * alpha0 + fracvol * alpha1;
}

real higflow_interp_Pe_multiphase_electroosmotic(real Pe0, real Pe1, real fracvol) {
    real D0 = 1.0/Pe0;
    real D1 = 1.0/Pe1;
    real D = (1.0 - fracvol) * D0 + fracvol * D1;
    return 1.0/D;
}

// real higflow_interp_delta_multiphase_electroosmotic(real delta0, real delta1, real perm0, real perm1, real fracvol) {
//     real perm = (1.0 - fracvol) * perm0 + fracvol * perm1;
//     return ((1 - fracvol) * perm0 * delta0 + fracvol * perm1 * delta1) / perm;
// }

real higflow_interp_delta_multiphase_electroosmotic(real delta0, real delta1, real fracvol) {
    return (1.0 - fracvol) * delta0 + fracvol * delta1;
}

real higflow_interp_Ex_multiphase_electroosmotic(real Ex0, real Ex1, real fracvol) {
    return (1.0 - fracvol) * Ex0 + fracvol * Ex1;
}


// *******************************************************************
// Ionic transport equations multiphase
// *******************************************************************


void higflow_explicit_euler_ionic_transport_equation_nplus_multiphase(higflow_solver *ns) {
    if (ns->ed.mult.eo.contr.eo_model == PNP) {
        flow_type eoflow0 = ns->ed.mult.contr.eoflow0;
        flow_type eoflow1 = ns->ed.mult.contr.eoflow1;
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

            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            // pure newtonian in this phase
            if ((FLT_EQ(fracvol, 0.0) && eoflow0 != true) ||
                (FLT_EQ(fracvol, 1.0) && eoflow1 != true)) {
                dp_set_value(ns->ed.eo.dpnplus, clid, 1.0);
                continue;
            }
            real alphaeo0 = ns->ed.mult.eo.par0.alpha; real alphaeo1 = ns->ed.mult.eo.par1.alpha;
            real Pe0 = ns->ed.mult.eo.par0.Pe; real Pe1 = ns->ed.mult.eo.par1.Pe;
            real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, 1.0);
            real Pe = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, 1.0);
            
            // compute nplus value at point 
            real nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                real fracvoll = compute_center_p_left(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real fracvolr = compute_center_p_right(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real alphaeol = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvoll);
                real alphaeor = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvolr);
                real Pel = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvoll);
                real Per = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvolr);
                // Set the computational cell in this particular direction
                higflow_computational_cell_multiphase_electroosmotic_ionic(ns, sdnplus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus, 
                                                                           alphaeo, alphaeol, alphaeor, Pe, Pel, Per);
                // Compute the diffusive ionic term rhs
                rhs    += higflow_diffusive_ionic_term(ns, Pe);
                // convective term
                switch (ns->ed.mult.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
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


void higflow_explicit_euler_ionic_transport_equation_nminus_multiphase(higflow_solver *ns) {
    if (ns->ed.mult.eo.contr.eo_model == PNP) {
        flow_type eoflow0 = ns->ed.mult.contr.eoflow0;
        flow_type eoflow1 = ns->ed.mult.contr.eoflow1;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus  = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(it); higcit_nextcell(it)) {
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

            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            // pure newtonian in this phase
            if ((FLT_EQ(fracvol, 0.0) && eoflow0 != true) ||
                (FLT_EQ(fracvol, 1.0) && eoflow1 != true)) {
                dp_set_value(ns->ed.eo.dpnplus, clid, 1.0);
                continue;
            }
            real alphaeo0 = ns->ed.mult.eo.par0.alpha; real alphaeo1 = ns->ed.mult.eo.par1.alpha;
            real Pe0 = ns->ed.mult.eo.par0.Pe; real Pe1 = ns->ed.mult.eo.par1.Pe;
            real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, 1.0);
            real Pe = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, 1.0);

            // compute nminus value at point 
            real nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                real fracvoll = compute_center_p_left(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real fracvolr = compute_center_p_right(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real alphaeol = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvoll);
                real alphaeor = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvolr);
                real Pel = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvoll);
                real Per = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvolr);
                // Set the computational cell in this particular direction
                higflow_computational_cell_multiphase_electroosmotic_ionic(ns, sdnminus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus,
                                                                           alphaeo, alphaeol, alphaeor, Pe, Pel, Per);
                // Compute the diffusive ionic term rhs
                rhs    += higflow_diffusive_ionic_term(ns, Pe);
                // convective term
                switch (ns->ed.mult.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
                        rhs    -= ns->cc.ucell * ns->cc.dndx;
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
            // Set property value  
            dp_set_value(ns->ed.eo.dpnminus_temp, clid, newnminus);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed ionic property
        dp_sync(ns->ed.eo.dpnminus_temp);
    }
}

void higflow_semi_implicit_euler_ionic_transport_equation_nplus_multiphase(higflow_solver *ns) {
    if (ns->ed.mult.eo.contr.eo_model == PNP) {
        flow_type eoflow0 = ns->ed.mult.contr.eoflow0;
        flow_type eoflow1 = ns->ed.mult.contr.eoflow1;
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

            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            // pure newtonian in this phase
            if ((FLT_EQ(fracvol, 0.0) && eoflow0 != true) ||
                (FLT_EQ(fracvol, 1.0) && eoflow1 != true)) {
                int cgid = psd_get_global_id(ns->ed.eo.psdEOnplus, c);
                slv_set_bi(ns->ed.eo.slvnplus, cgid, 1.0);
                int j = cgid;
                real val = 1.0;
                slv_set_Ai(ns->ed.eo.slvnplus, cgid, 1, &j, &val);
                continue;
            }
            real alphaeo0 = ns->ed.mult.eo.par0.alpha; real alphaeo1 = ns->ed.mult.eo.par1.alpha;
            real Pe0 = ns->ed.mult.eo.par0.Pe; real Pe1 = ns->ed.mult.eo.par1.Pe;
            real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, 1.0);
            real Pe = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, 1.0);
            
            // compute nplus value at point 
            real nplus = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            real fracvoll[DIM], fracvolr[DIM], Pel[DIM], Per[DIM];
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                fracvoll[dim] = compute_center_p_left(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                fracvolr[dim] = compute_center_p_right(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real alphaeol = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvoll[dim]);
                real alphaeor = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvolr[dim]);
                Pel[dim] = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvoll[dim]);
                Per[dim] = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvolr[dim]);
                // Set the computational cell in this particular direction
                higflow_computational_cell_multiphase_electroosmotic_ionic(ns, sdnplus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnplus, ns->ed.eo.stnnplus,
                                                                           alphaeo, alphaeol, alphaeor, Pe, Pel[dim], Per[dim]);
                // convective term
                switch (ns->ed.mult.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
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
            real alpha = 0.0, Pelc, Perc;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // consider non-uniform diffusivity around the interface
                Pelc = 0.5 * (Pel[dim2] + Pe);
                Perc = 0.5 * (Per[dim2] + Pe);
                // Stencil weight update
                real wl = -ns->par.dt/(Pelc*cdelta[dim2]*cdelta[dim2]);
                real wr = -ns->par.dt/(Perc*cdelta[dim2]*cdelta[dim2]);
                alpha  -= (wl + wr) ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, wr, ns->ed.eo.stnnplus);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnplus, ccenter, p, wl, ns->ed.eo.stnnplus);                
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

void higflow_semi_implicit_euler_ionic_transport_equation_nminus_multiphase(higflow_solver *ns) {
    if (ns->ed.mult.eo.contr.eo_model == PNP) {
        flow_type eoflow0 = ns->ed.mult.contr.eoflow0;
        flow_type eoflow1 = ns->ed.mult.contr.eoflow1;
        // Get the local sub-domain for the cells
        sim_domain *sdnminus  = psd_get_local_domain(ns->ed.eo.psdEOnminus);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdnminus);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdnminus); !higcit_isfinished(it); higcit_nextcell(it)) {
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

            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            // pure newtonian in this phase
            if ((FLT_EQ(fracvol, 0.0) && eoflow0 != true) ||
                (FLT_EQ(fracvol, 1.0) && eoflow1 != true)) {
                int cgid = psd_get_global_id(ns->ed.eo.psdEOnminus, c);
                slv_set_bi(ns->ed.eo.slvnminus, cgid, 1.0);
                int j = cgid;
                real val = 1.0;
                slv_set_Ai(ns->ed.eo.slvnminus, cgid, 1, &j, &val);
                continue;
            }
            real alphaeo0 = ns->ed.mult.eo.par0.alpha; real alphaeo1 = ns->ed.mult.eo.par1.alpha;
            real Pe0 = ns->ed.mult.eo.par0.Pe; real Pe1 = ns->ed.mult.eo.par1.Pe;
            real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, 1.0);
            real Pe = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, 1.0);
            
            // compute nminus value at point 
            real nminus = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus);
            // Solving the Transport Equation using the Euler Method
            // Right hand side equation
            real rhs = 0.0;
            real fracvoll[DIM], fracvolr[DIM], Pel[DIM], Per[DIM];
            for (int dim = 0; dim < DIM; dim++) { // sum in both directions
                fracvoll[dim] = compute_center_p_left(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                fracvolr[dim] = compute_center_p_right(ns->ed.mult.sdmult, ccenter, cdelta, dim, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
                real alphaeol = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvoll[dim]);
                real alphaeor = higflow_interp_alpha_multiphase_electroosmotic(alphaeo0, alphaeo1, fracvolr[dim]);
                Pel[dim] = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvoll[dim]);
                Per[dim] = higflow_interp_Pe_multiphase_electroosmotic(Pe0, Pe1, fracvolr[dim]);
                // Set the computational cell in this particular direction
                higflow_computational_cell_multiphase_electroosmotic_ionic(ns, sdnminus, clid, ccenter, cdelta, dim, ns->ed.eo.dpnminus, ns->ed.eo.stnnminus,
                                                                           alphaeo, alphaeol, alphaeor, Pe, Pel[dim], Per[dim]);
                // convective term
                switch (ns->ed.mult.eo.contr.convecdiscrtype) {
                    case CELL_CENTRAL: // Central scheme
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
            real alpha = 0.0, Pelc, Perc;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // consider non-uniform diffusivity around the interface
                Pelc = 0.5 * (Pel[dim2] + Pe);
                Perc = 0.5 * (Per[dim2] + Pe);
                // Stencil weight update
                real wl = -ns->par.dt/(Pelc*cdelta[dim2]*cdelta[dim2]);
                real wr = -ns->par.dt/(Perc*cdelta[dim2]*cdelta[dim2]);
                alpha  -= (wl + wr) ;
                Point p;
                POINT_ASSIGN(p, ccenter);
                // Stencil point update: right point
                p[dim2] = ccenter[dim2] + cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, wr, ns->ed.eo.stnnminus);
                // Stencil point update: left point
                p[dim2] = ccenter[dim2] - cdelta[dim2];
                sd_get_stencil(sdnminus, ccenter, p, wl, ns->ed.eo.stnnminus);                
            }
            alpha = 1.0 + alpha;
            // Get the stencil
            sd_get_stencil(sdnminus, ccenter, ccenter, alpha, ns->ed.eo.stnnminus);
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


// Electroosmotic induced potential psi multiphase
real higflow_multiphase_electroosmotic_psi(higflow_solver *ns) {
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

        real fracvol = compute_value_at_point(ns->ed.mult.sdmult, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        real deltaeo = higflow_interp_delta_multiphase_electroosmotic(ns->ed.mult.eo.par0.delta, ns->ed.mult.eo.par1.delta, fracvol);
        real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(ns->ed.mult.eo.par0.alpha, ns->ed.mult.eo.par1.alpha, 1.0);

        // Reset the stencil
        stn_reset(ns->ed.eo.stnpsi);
        // Initialize rhs
        real rhs = 0.0;
        // Set the right side of stencil
        if (ns->ed.mult.eo.contr.eo_model == PNP) {
            // Poisson-Nernst-Planck model 
            // Get the ionic concentration n- at center cell
            nplus_temp    = compute_value_at_point(ns->ed.eo.sdEOnplus, ccenter, ccenter, 1.0, ns->ed.eo.dpnplus_temp, ns->ed.eo.stnnplus);
            nminus_temp   = compute_value_at_point(ns->ed.eo.sdEOnminus, ccenter, ccenter, 1.0, ns->ed.eo.dpnminus_temp, ns->ed.eo.stnnminus);
            rhs      = -deltaeo*(nplus_temp - nminus_temp);
        } else if (ns->ed.mult.eo.contr.eo_model == PB) {
            // Poisson-Boltzmann model 
            psi      = dp_get_value(ns->ed.eo.dppsi, clid);
            rhs      = 2.0*deltaeo*sinh(alphaeo*psi) - 2.0*deltaeo*alphaeo*psi*cosh(alphaeo*psi);
        }
        // Calculate the point and weight of the stencil
        stn_set_rhs(ns->ed.eo.stnpsi, rhs);
        real alpha = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            // Stencil point update
            Point p;
            POINT_ASSIGN(p, ccenter);
            real perm = ns->ed.mult.eo.get_permittivity(fracvol, p, ns->par.t);
            // Stencil point update: right point
            p[dim] = ccenter[dim] + cdelta[dim];
            real fracvolr = compute_value_at_point(ns->ed.mult.sdmult, p, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real permr = ns->ed.mult.eo.get_permittivity(fracvolr, p, ns->par.t);
            real permrc = 0.5*(perm + permr);
            real wr = permrc/(cdelta[dim]*cdelta[dim]);
            sd_get_stencil(sdp, ccenter, p, wr, ns->ed.eo.stnpsi);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            real fracvoll = compute_value_at_point(ns->ed.mult.sdmult, p, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real perml = ns->ed.mult.eo.get_permittivity(fracvoll, p, ns->par.t);
            real permlc = 0.5*(perm + perml);
            real wl = permlc/(cdelta[dim]*cdelta[dim]);
            sd_get_stencil(sdp, ccenter, p, wl, ns->ed.eo.stnpsi);

            alpha -= (wr + wl);
        }
        switch (ns->ed.mult.eo.contr.eo_model) {
        case PB:
            // Poisson-Boltzmann model
            alpha -= 2.0*alphaeo*deltaeo*cosh(alphaeo*psi); 
            break;
        case PBDH:
            // Debye-Hückel model (solving poisson equation for psi) 
            alpha -= 2.0*alphaeo*deltaeo;
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


void higflow_multiphase_electroosmotic_solve_pb(higflow_solver *ns) {
    real max_psi_res = INFINITY, max_psi_res_global = INFINITY;
    real pb_tol = EPSMACH;
    int iter = 0, maxiter = 50;
    
    while (max_psi_res_global > pb_tol && iter < maxiter) {
        max_psi_res = higflow_multiphase_electroosmotic_psi(ns);
        MPI_Allreduce(&max_psi_res, &max_psi_res_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        print0f("Poisson-Boltzmann model residual: %e\n", max_psi_res_global);
        iter++;
    }
    if (iter < maxiter) {
        printf("Poisson-Boltzmann model converged after %d iterations\n", iter);
    }
    else {
        printf("Poisson-Boltzmann model did not converge after %d iterations\n", iter);
        exit(1);
    }
}

// *******************************************************************
// Electro-osmotic source term
// *******************************************************************
void higflow_calculate_multiphase_electroosmotic_source_term(higflow_solver *ns) {
    real phil, phir, psil, psir, dphidx, dpsidx, npl, npr, nml, nmr, nplus, nminus, rho, psi, Feo;
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

            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real deltaeo = higflow_interp_delta_multiphase_electroosmotic(ns->ed.mult.eo.par0.delta, ns->ed.mult.eo.par1.delta, fracvol);
            real alphaeo = higflow_interp_alpha_multiphase_electroosmotic(ns->ed.mult.eo.par0.alpha, ns->ed.mult.eo.par1.alpha, 1.0);

            // Get the derivative of applied potential
            phil        = compute_center_p_left(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            phir        = compute_center_p_right(sdphi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dpphi, ns->ed.eo.stnphi);
            dphidx = compute_dpdx_at_point(fdelta, dim, 0.5, phil, phir);
            // Get the derivative of induced potential
            psil        = compute_center_p_left(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            psir        = compute_center_p_right(sdpsi, fcenter, fdelta, dim, 0.5, ns->ed.eo.dppsi, ns->ed.eo.stnpsi);
            dpsidx = compute_dpdx_at_point(fdelta, dim, 0.5, psil, psir);
            
            switch (ns->ed.mult.eo.contr.eo_model) {
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
                    rho = (nplus - nminus)*deltaeo;
                    break;
                case PB:
                    // Poisson-Boltzmann model 
                    // Get the induced potential and its derivative
                    psi    = compute_value_at_mid_point(psil, psir);
                    // Compute the charge density 
                    rho    = -2.0*deltaeo*sinh(alphaeo*psi);
                    break;
                case PBDH:
                    // Debye-Hückel model (solving poisson equation for psi) 
                    // Get the induced potential and its derivative
                    psi    = compute_value_at_mid_point(psil, psir);
                    // Compute the charge density 
                    rho    = -2.0*alphaeo*deltaeo*psi;
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
            real fracvoll = compute_value_at_point(ns->ed.mult.sdmult, ccenterl, ccenterl, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real perml = ns->ed.mult.eo.get_permittivity(fracvoll, ccenterl, ns->par.t);
            real fracvolr = compute_value_at_point(ns->ed.mult.sdmult, ccenterr, ccenterr, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real permr = ns->ed.mult.eo.get_permittivity(fracvolr, ccenterr, ns->par.t);
            real dpermdx = compute_dpdx_at_point(fdelta, dim, 0.5, perml, permr);

            // Extra term arising due to (possibly non-uniform) permittivity gradient
            // Comes from Korteweg-Helmholtz force for linear dielectric isotropic media with instantaneous polarization response
            Feo += -0.5*normE2*dpermdx;
            
            // Get the electroosmotic extra source term defined by user
            Feo   += ns->ed.mult.eo.get_multiphase_electroosmotic_source_term(fracvol, fcenter, dim, ns->par.t);
            // Set the distributed source term property
            dp_set_value(ns->ed.eo.dpFeo[dim], flid, Feo);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->ed.eo.dpFeo[dim]);
    }
}



// *******************************************************************
// Navier-Stokes Step for the Explicit Euler Method
// *******************************************************************
void higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real Ex = higflow_interp_Ex_multiphase_electroosmotic(ns->ed.mult.eo.par0.Ex, ns->ed.mult.eo.par1.Ex, fracvol);
            rhs += higflow_electroosmotic_source_term(ns, Ex);
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
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
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
void higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns) {
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
void higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns) {
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
void higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns) {
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
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real Ex = higflow_interp_Ex_multiphase_electroosmotic(ns->ed.mult.eo.par0.Ex, ns->ed.mult.eo.par1.Ex, fracvol);
            rhs += higflow_electroosmotic_source_term(ns, Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
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
                real wr = - ns->par.dt*ns->cc.viscr[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - ns->par.dt*ns->cc.viscl[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns) {
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
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real Ex = higflow_interp_Ex_multiphase_electroosmotic(ns->ed.mult.eo.par0.Ex, ns->ed.mult.eo.par1.Ex, fracvol);
            rhs += higflow_electroosmotic_source_term(ns, Ex);
            // Diffusive term term contribution
            rhs += 0.5*higflow_diffusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
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
                real wr = - 0.5*ns->par.dt*ns->cc.viscr[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 0.5*ns->par.dt*ns->cc.viscl[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
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
void higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(higflow_solver *ns) {
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
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real Ex = higflow_interp_Ex_multiphase_electroosmotic(ns->ed.mult.eo.par0.Ex, ns->ed.mult.eo.par1.Ex, fracvol);
            rhs += higflow_electroosmotic_source_term(ns, Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
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
                real wr = - 0.25*ns->par.dt*ns->cc.viscr[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 0.25*ns->par.dt*ns->cc.viscl[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
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
            higflow_computational_cell_multiphase(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            // Electo-osmotic source term
            real fracvol = compute_value_at_point(ns->ed.mult.sdmult, fcenter, fcenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
            real Ex = higflow_interp_Ex_multiphase_electroosmotic(ns->ed.mult.eo.par0.Ex, ns->ed.mult.eo.par1.Ex, fracvol);
            rhs += higflow_electroosmotic_source_term(ns, Ex);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Interfacial tension contribution
            rhs += higflow_interfacial_tension_term(ns);
            // Cell term contribution for the gravity
            if (dim == 1) rhs -= higflow_gravity_term(ns);
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
                real wr = - 1.0/3.0*ns->par.dt*ns->cc.viscr[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                real wl = - 1.0/3.0*ns->par.dt*ns->cc.viscl[dim2]/(ns->par.Re*fdelta[dim2]*fdelta[dim2])/ns->cc.dens;
                alpha -= (wr + wl);
                Point p;
                POINT_ASSIGN(p, fcenter);
                // Stencil point update: right point
                p[dim2] = fcenter[dim2] + fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wr, ns->stn);
                // Stencil point update: left point
                p[dim2] = fcenter[dim2] - fdelta[dim2];
                sfd_get_stencil(sfdu[dim], fcenter, p, wl, ns->stn);
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


// One step of the Navier-Stokes the projection method
void higflow_solver_step_multiphase_electroosmotic(higflow_solver *ns) {
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
    switch (ns->ed.mult.eo.contr.eo_model) {
    case PNP: // Poisson-Nernst-Planck model
        for(int k=0; k<3; k++) {
            higflow_boundary_condition_for_electroosmotic_nplus(ns);
            higflow_boundary_condition_for_electroosmotic_nminus(ns);
            if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
            higflow_electroosmotic_phi(ns);
            switch (ns->ed.mult.eo.contr.tempdiscrtype) {
                case EXPLICIT_EULER:
                    higflow_explicit_euler_ionic_transport_equation_nplus_multiphase(ns);
                    higflow_explicit_euler_ionic_transport_equation_nminus_multiphase(ns);
                    break;
                case SEMI_IMPLICIT_EULER:
                    higflow_semi_implicit_euler_ionic_transport_equation_nplus_multiphase(ns);
                    higflow_semi_implicit_euler_ionic_transport_equation_nminus_multiphase(ns);
                    break;
            }
            real max_psi_res = higflow_multiphase_electroosmotic_psi(ns);
            if(k>0) print0f("=+=+=+ psi residual in inner iteration %d: %15.10lf =+=+=+\n", k+1, max_psi_res);
        }
        dp_copy_values(ns->ed.eo.dpnplus, ns->ed.eo.dpnplus_temp);
        dp_copy_values(ns->ed.eo.dpnminus, ns->ed.eo.dpnminus_temp);
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PB: // Poisson-Boltzmann model 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.eo.contr.is_psibc_timedependent == true)
        higflow_multiphase_electroosmotic_solve_pb(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true
                              || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
        higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PBDH: // Debye-Hückel model (solving poisson equation for psi) 
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true)
        higflow_electroosmotic_phi(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
            higflow_multiphase_electroosmotic_psi(ns);
        if( (ns->par.step == 0) || ns->ed.mult.eo.contr.is_phibc_timedependent == true
                              || ns->ed.mult.eo.contr.is_psibc_timedependent == true)
            higflow_calculate_multiphase_electroosmotic_source_term(ns);
        break;
    case PBDH_ANALYTIC:
        printf("PBDH_ANALYTIC not implemented\n");
        break;
    }

    // Interpolate the viscosity and density
    higflow_compute_viscosity_multiphase(ns);
    higflow_compute_density_multiphase(ns);
    // Calculate the curvature, interfacial force and normal
    higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
    higflow_compute_distance_multiphase_2D(ns);
    higflow_compute_plic_lines_2d(ns);

    // Calculate the intermediate velocity
    switch (ns->contr.tempdiscrtype) {
    case EXPLICIT_EULER:
        higflow_explicit_euler_intermediate_velocity_multiphase_electroosmotic(ns, ns->dpu, ns->dpustar);
        break;
    case EXPLICIT_RK2:
        higflow_explicit_runge_kutta_2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case EXPLICIT_RK3:
        higflow_explicit_runge_kutta_3_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_EULER:
        higflow_semi_implicit_euler_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_CN:
        higflow_semi_implicit_crank_nicolson_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    case SEMI_IMPLICIT_BDF2:
        higflow_semi_implicit_bdf2_intermediate_velocity_multiphase_electroosmotic(ns);
        break;
    }

    // Projection
    higflow_pressure_multiphase(ns);
    higflow_final_velocity_multiphase(ns);
    higflow_final_pressure(ns);

    // Calculate the volume fraction
    higflow_plic_advection_volume_fraction(ns);
}

