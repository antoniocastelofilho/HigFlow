// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Step - version 10/11/2016
// *******************************************************************
// *******************************************************************

#include "hig-flow-step.h"

// *******************************************************************
// Navier-Stokes step elements
// *******************************************************************

// Remove singularity for the pressure
void remove_pressure_singularity(higflow_solver *ns, solver *slvp) {
    if (ns->contr.desingpressure == true) {
        int cgid, clid;
        real value;
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        real t, dt;
        if (rank == 0) {
            sim_domain *sdp = psd_get_local_domain(ns->psdp);
            mp_mapper *mp = sd_get_domain_mapper(sdp);
            Point ccenter, ccenter2;
            // Toma a primeira celula que aparecer
            hig_cell *c = sd_get_higtree(sdp, 0);
            hig_get_center(c, ccenter);
            hig_cell *c2;
            c2 = hig_get_cell_with_point(c, ccenter);
            clid = mp_lookup(mp, hig_get_cid(c2));
            hig_get_center(c2, ccenter2);
            t = ns->par.t;
            dt = ns->par.dt;
            // Calcula valor do lado direito
            if (ns->contr.projtype == INCREMENTAL) {
                value = ns->func.get_pressure(ccenter2,t+dt) - dp_get_value(ns->dpp, clid);
            } else {
                value = ns->func.get_pressure(ccenter2,t+dt);
            }
            cgid = psd_lid_to_gid(ns->psdp, clid);
       slv_impose_value(slvp, cgid, value);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Navier-Stokes pressure using the projection method
void higflow_pressure(higflow_solver *ns) {
    real uoutflow;
    // Get the local sub-domain for the cells
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    // Get the local sub-domain for the facets
    sim_facet_domain *sfdu[DIM];
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
    }
    // aqui é pra remover
    remove_pressure_singularity(ns, ns->slvp);
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
        // Calculate the divergence of the intermediate velocity
        real sumdudx = 0.0;
        for(int dim = 0; dim < DIM; dim++) {
            int infacet;
            // Get the left facet intermediate velocity
            real dudx;
            real ustarl = compute_facet_u_left(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpustar[dim], ns->stn, &infacet);
            // Get the left facet intermediate velocity
            real ustarr = compute_facet_u_right(sfdu[dim], ccenter, cdelta, dim, 0.5, ns->dpustar[dim], ns->stn, &infacet);
            // Compute the derivative of intermediate velocity
            dudx   = compute_facet_dudxc(cdelta, dim, 0.5, ustarl, ustarl, ustarr);
            sumdudx += dudx;
        }
        // Reset the stencil
        stn_reset(ns->stn);
        // Set the right side of stencil
        stn_set_rhs(ns->stn, sumdudx / ns->par.dt);
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
            sd_get_stencil(sdp, ccenter, p, w, ns->stn);
            // Stencil point update: left point
            p[dim] = ccenter[dim] - cdelta[dim];
            sd_get_stencil(sdp, ccenter, p, w, ns->stn);
        }
        // Get the stencil
        sd_get_stencil(sdp, ccenter, ccenter, alpha, ns->stn);
        // Get the index of the stencil
        int *ids   = psd_stn_get_gids(ns->psdp, ns->stn);
        // Get the value of the stencil
        real *vals = stn_get_vals(ns->stn);
        // Get the number of elements of the stencil
        int numelems = stn_get_numelems(ns->stn);
        // Get the cell identifier of the cell
        int cgid = psd_get_global_id(ns->psdp, c);
        if (ns->contr.desingpressure == true) { // cell with pressure desingularized
            if(cgid==ns->slvp->imposed_line) continue;
        }
        // Set the right side of solver linear system
        slv_set_bi(ns->slvp, cgid, stn_get_rhs(ns->stn));
        // Set the line of matrix of the solver linear system
        slv_set_Ai(ns->slvp, cgid, numelems, ids, vals);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Assemble the solver
    slv_assemble(ns->slvp);
    // Solve the linear system
    slv_solve(ns->slvp);
    // Set the solver solution in the pressure ou pressure difference distributed property
    distributed_property *dp = (ns->contr.projtype == INCREMENTAL) ? ns->ddeltap : ns->dpp;
    dp_slv_load_from_solver(dp, ns->slvp);
    dp_sync(dp);
}

// Navier-Stokes outflow for u velocity
void higflow_outflow_u_step(higflow_solver *ns) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Set the velocity at outflow
   set_outflow(ns->psfdu[dim], ns->dpu[dim], 1.0);
    }
}

// Navier-Stokes outflow for ustar velocity
void higflow_outflow_ustar_step(higflow_solver *ns) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Set the velocity at outflow
   set_outflow(ns->psfdu[dim], ns->dpustar[dim], 1.0);
    }
}

void set_outflow(psim_facet_domain *psfdu, distributed_property *dpu, real alpha) {
   sim_facet_domain *sfd = psfd_get_local_domain(psfdu);
   int dim = sfd_get_dim(sfd);
   int numnbcs = sfd_get_num_bcs(sfd, NEUMANN);
   int numhigs = sfd_get_num_higtrees(sfd);
   for (int i = 0; i < numnbcs; i++) {
      sim_boundary *sb = sfd_get_bc(sfd, NEUMANN, i);
      hig_cell *bchig = sb_get_higtree(sb);
      Point bcl, bch;
      hig_get_lowpoint(bchig, bcl);
      hig_get_highpoint(bchig, bch);
      int bcdim = hig_get_narrowest_dim(bchig);
      POINT_SUB_SCALAR(bcl, bcl, EPSDELTA);
      POINT_ADD_SCALAR(bch, bch, EPSDELTA);
      for (int j = 0; j < numhigs; j++) {
         hig_cell *root = sfd_get_higtree(sfd, j);
         mp_mapper *m = sfd_get_domain_mapper(sfd);
         Point delta;
         Point center;
         higcit_celliterator *it =
             higcit_create_bounding_box(root, bcl, bch);
         if (higcit_isfinished(it)) {
            higcit_destroy(it);
            continue;
         }
         hig_cell *c = higcit_getcell(it);
         hig_get_delta(c, delta);
         hig_get_center(c, center);
         higcit_destroy(it);
         Point backwarddelta;
         POINT_ASSIGN_SCALAR(backwarddelta, 0.0);
         Point flowl, flowh;
         POINT_ASSIGN(flowl, bcl);
         POINT_ASSIGN(flowh, bch);
         if (center[bcdim] < bcl[bcdim]) {
            backwarddelta[bcdim] = -alpha * delta[bcdim];
            flowl[bcdim] -= alpha * delta[bcdim];
            //printf("Positivo: center[%d]=%5.10f; bcl[%d]=%5.10f; bch[%d]=%5.10f\n",bcdim,center[bcdim],bcdim,bcl[bcdim],bcdim,bch[bcdim]);
         } else {
            backwarddelta[bcdim] = alpha * delta[bcdim];
            flowh[bcdim] += alpha * delta[bcdim];
            //printf("Negativo: center[%d]=%5.10f; bcl[%d]=%5.10f; bch[%d]=%5.10f\n",bcdim,center[bcdim],bcdim,bcl[bcdim],bcdim,bch[bcdim]);
         }
         int dimofinterest[DIM];
         POINT_ASSIGN_SCALAR(dimofinterest, 0);
         dimofinterest[dim] = 1;
         higfit_facetiterator *fit;
         for (fit =
              higfit_create_bounding_box_facets(root,
                            dimofinterest,
                            flowl, flowh);
              !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(m, hig_get_fid(f));
            Point fcenter;
            hig_get_facet_center(f, fcenter);

            hig_facet fc;
            Point fcopy;
            POINT_ADD(fcopy, fcenter, backwarddelta);

            hig_get_facet_with_point(root, dim, fcopy, &fc);
            int fclid = mp_lookup(m, hig_get_fid(&fc));

                      //DEBUG_INSPECT(fcenter[0], %lf);
                      //DEBUG_INSPECT(fcenter[1], %lf);
            //DEBUG_INSPECT(flid, %d);

            if (flid >= 0 && fclid >= 0) {

               dp_set_value(dpu, flid,
                       dp_get_value(dpu, fclid));
            }
         }
              higfit_destroy(fit);
      }
   }
}

// Navier-Stokes final velocity using the projection method
void higflow_final_velocity(higflow_solver *ns) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdu[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Initialize the min and max velocity
        real velmax = -1.0e16;
        real velmin =  1.0e16;
        // Get the local partitioned domain for facets
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);

        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            real pl, pr;
            if (ns->contr.projtype == INCREMENTAL) {
                // Get the pressure in the left cell
                pl    = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->ddeltap, ns->stn);
                // Get the pressure in the right cell
                pr    = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->ddeltap, ns->stn);
            } else {
                // Get the pressure in the left cell
                pl    = compute_center_p_left(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
                // Get the pressure in the right cell
                pr    = compute_center_p_right(sdp, fcenter, fdelta, dim, 0.5, ns->dpp, ns->stn);
            }
            // Compute the pressure derivative
            real dpdx  = compute_dpdx_at_point(fdelta, dim, 0.5, pl, pr);
            // Get the intermediate velocity
            real ustar = dp_get_value(ns->dpustar[dim], flid);
            // Compute the final velocity
            real u  = ustar - ns->par.dt * dpdx;

            UPDATE_RESIDUAL_BUFFER_FACET(ns, dp_get_value(ns->dpu[dim], flid), u, f, fcenter)

            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpu[dim], flid, u);
            // Update the min and max velocity
            if (u > velmax) velmax = u;
            if (u < velmin) velmin = u;
        }
        // Destroy the iterator
        higfit_destroy(fit);

        UPDATE_RESIDUALS(ns, ns->residuals->u[dim])

        // Sync the distributed velocity property
        dp_sync(ns->dpu[dim]);
        // Printing the min and max velocity
        real velmin_global, velmax_global;
        MPI_Allreduce(&velmin, &velmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&velmax, &velmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        print0f("===> %d: Vmin = %15.10lf <===> Vmax = %15.10lf <===\n",dim,velmin_global,velmax_global);
    }
}

// Navier-Stokes final pressure using the projection method
void higflow_final_pressure(higflow_solver *ns) {
    // Projection method
    if (ns->contr.projtype == INCREMENTAL) {
        // Incremental projection method
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the map of the distributd properties in the cells
        mp_mapper  *mp  = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the pressure in the distributed pressure property
            real p      = dp_get_value(ns->dpp, clid);
            // Get the pressure difference in the distributed difference pressure property
            real deltap = dp_get_value(ns->ddeltap, clid);
            // Calculate the final pressure
            real newp   = p + deltap;
            // Set the final pressure in the distributed pressure property
            dp_set_value(ns->dpp, clid, newp);
        }
        // Destroy the iterator
        higcit_destroy(it);
        // Sync the distributed pressure property
        dp_sync(ns->dpp);
    }
}

// Navier-Stokes calculate the source term
void higflow_calculate_source_term(higflow_solver *ns) {
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdF);
    // Get the map of the distributd properties in the cells
    mp_mapper  *mp  = sd_get_domain_mapper(sdp);
    // Loop for each cell
    higcit_celliterator *it;
    for(it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        // Get the cell identifier
        int clid    = mp_lookup(mp, hig_get_cid(c));
        // Get the cell center
        Point ccenter;
        hig_get_center(c, ccenter);
        // Set the source term
        real F      = ns->func.get_source_term(ccenter, ns->par.t);
        // Set the final pressure in the distributed pressure property
        dp_set_value(ns->dpF, clid, F);
    }
    // Destroy the iterator
    higcit_destroy(it);
    // Sync the distributed pressure property
    dp_sync(ns->dpF);
}

// Navier-Stokes calculate the facet source term
void higflow_calculate_facet_source_term(higflow_solver *ns) {
    // Get the local sub-domain
    //sim_domain *sdp = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfdF[DIM];
    // Loop for each dimension
    higfit_facetiterator *fit;
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local partitioned domain for facets
        sfdF[dim] = psfd_get_local_domain(ns->psfdF[dim]);
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
            // Get the pressure in the left cell
            real F = ns->func.get_facet_source_term(fcenter, dim, ns->par.t);
            // Set the final velocity in the distributed velocity property
            dp_set_value(ns->dpFU[dim], flid, F);
        }
        // Destroy the iterator
        higfit_destroy(fit);
        // Sync the distributed velocity property
        dp_sync(ns->dpFU[dim]);
    }
}

// Apply the boundary condition for the velocity
void higflow_boundary_condition_for_velocity(higflow_solver *ns) {
    // Facet iterator
    higcit_celliterator *it;
    // Local sub-domain
    sim_facet_domain *sfdu[DIM];
    // For each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local sub-domain
        sfdu[dim]      = psfd_get_local_domain(ns->psfdu[dim]);
        sim_domain *sd = sfdu[dim]->cdom;

        int num_bc_types = 2;
        bc_type bc_t;
        for(int i = 0; i < num_bc_types; i++) {
            if(i == 0) bc_t = DIRICHLET;
            else if(i == 1) bc_t = NEUMANN;

            // Get the number of boundaries of type
            int numbcs = sd_get_num_bcs(sd, bc_t);
            // For each boundary
            for (int i = 0; i < numbcs; i++) {
                // Get the boundary
                sim_boundary *bc = sd_get_bc(sd, bc_t, i);
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
                        real val = ns->func.get_boundary_velocity(userid, bccenter, dim, t);
                        // Set the value
                        sb_set_value(bc, bclid, val);
                    }
                    // Destroy the iterator
                    higcit_destroy(it);
                }
            }
        }
    }
}

// Apply the boundary condition for the pressure
void higflow_boundary_condition_for_pressure(higflow_solver *ns) {
    // Facet iterator
    higcit_celliterator *it;
    // Get the local sub-domain
    sim_domain *sdp = psd_get_local_domain(ns->psdp);

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
                    real val;
                    if (ns->contr.projtype == 0) {
                        // Non incremental projection method
                        val = ns->func.get_boundary_pressure(userid, bccenter, t);
                    } else {
                        // Incremental projection method
                        val = ns->func.get_boundary_pressure(userid, bccenter, t) -
                            ns->func.get_boundary_pressure(userid, bccenter, ns->par.t);
                    }
                    // Set the value
                    sb_set_value(bc, bclid, val);
                }
                // Destroy the iterator
                higcit_destroy(it);
            }
        }
    }

}


// Apply the boundary condition for the cell source term
void higflow_boundary_condition_for_cell_source_term(higflow_solver *ns) {
    // Facet iterator
    higcit_celliterator *it;
    // Local sub-domain
    sim_domain *sd = psd_get_local_domain(ns->psdF);
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
                real val = ns->func.get_boundary_source_term(userid, bccenter, t);
                // Set the value
                sb_set_value(bc, bclid, val);
            }
            // Destroy the iterator
            higcit_destroy(it);
        }
    }
}


// Apply the boundary condition for the facet source term
void higflow_boundary_condition_for_facet_source_term(higflow_solver *ns) {
    // Facet iterator
    higcit_celliterator *it;
    // Local sub-domain
    sim_facet_domain *sfdF[DIM];
    // For each dimension
    for(int dim = 0; dim < DIM; dim++) {
        // Get the local sub-domain
        sfdF[dim]      = psfd_get_local_domain(ns->psfdF[dim]);
        sim_domain *sd = sfdF[dim]->cdom;
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
                    real val = ns->func.get_boundary_facet_source_term(userid, bccenter, dim, t);
                    // Set the value
                    sb_set_value(bc, bclid, val);
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
void higflow_explicit_euler_intermediate_velocity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Difusive term contribution
            rhs += higflow_difusive_term(ns, fdelta);
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
    }
}


// *******************************************************************
// Navier-Stokes Step for the Second Order Explicit Runge-Kutta Method
// *******************************************************************
void higflow_explicit_runge_kutta_2_intermediate_velocity(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit Euler method
    higflow_explicit_euler_intermediate_velocity(ns, ns->dpu, ns->dpuaux);
    // Calculate the star velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity(ns, ns->dpuaux, ns->dpustar);
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
void higflow_explicit_runge_kutta_3_intermediate_velocity(higflow_solver *ns) {
    // Calculate the auxiliar velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity(ns, ns->dpu, ns->dpuaux);
    // Calculate the second stage velocity by the explicit euler method
    higflow_explicit_euler_intermediate_velocity(ns, ns->dpuaux, ns->dpustar);
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
    higflow_explicit_euler_intermediate_velocity(ns, ns->dpuaux, ns->dpustar);
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
void higflow_semi_implicit_euler_intermediate_velocity(higflow_solver *ns) {
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
            higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
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
            real w = - ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
            alpha -= 2.0 * w ;
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
void higflow_semi_implicit_crank_nicolson_intermediate_velocity(higflow_solver *ns) {
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
            higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
            // Diffusive term term contribution
            rhs += 0.5 * higflow_difusive_term(ns, fdelta);
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
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
            real w = - 0.5 * ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
            alpha -= 2.0 * w ;
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
void higflow_semi_implicit_bdf2_intermediate_velocity(higflow_solver *ns, distributed_property *dpu[DIM], distributed_property *dpustar[DIM]) {
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
            higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            // Right hand side equation
            real rhs = 0.0;
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
            rhs += ns->cc.ufacet;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 0.25*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -= 2.0 * w ; //divide po 4 para usar regra trapezio em t(n+1/2)
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
            higflow_computational_cell(ns, sdp, sfdu, flid, fcenter, fdelta, dim, ns->dpu);
            //Get the uaux
            real uaux = dp_get_value(ns->dpuaux[dim], flid);
            // Right hand side equation
            real rhs = 0.0;
            // Source term contribution
            rhs += higflow_source_term(ns);
            // Pressure term contribution
            rhs -= higflow_pressure_term(ns);
            // Convective term contribution
            rhs -= higflow_convective_term(ns, fdelta, dim);
            // Total contribuition terms times delta t
            rhs *= 1.0/3.0*ns->par.dt;
            rhs += (4.0*uaux - ns->cc.ufacet)/3.0;
            // Reset the stencil
            stn_reset(ns->stn);
            // Set the right side of stencil
            stn_set_rhs(ns->stn,rhs);
            // Calculate the point and weight of the stencil
            real alpha = 0.0;
            for(int dim2 = 0; dim2 < DIM; dim2++) {
                // Stencil weight update
                real w = - 1.0/3.0*ns->par.dt/(ns->par.Re*fdelta[dim2]*fdelta[dim2]);
                alpha -=  2.0 * w ;
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

// One step of the Navier-Stokes the projection method
void higflow_solver_step(higflow_solver *ns) {
    // Boundary condition for velocity
    higflow_boundary_condition_for_velocity(ns);
    // Boundary conditions for source term
    higflow_boundary_condition_for_cell_source_term(ns);
    higflow_boundary_condition_for_facet_source_term(ns);
    // Boundary condition for pressure
    higflow_boundary_condition_for_pressure(ns);
    // Calculate the source term
    higflow_calculate_source_term(ns);
    // Calculate the facet source term
    higflow_calculate_facet_source_term(ns);
    
    // Calculate the intermediated velocity
    switch (ns->contr.tempdiscrtype) {
        case EXPLICIT_EULER:
           // Explicit Euler method
           higflow_explicit_euler_intermediate_velocity(ns, ns->dpu, ns->dpustar);
           break;
        case EXPLICIT_RK2:
           // Explicit RK2 method
           higflow_explicit_runge_kutta_2_intermediate_velocity(ns);
           break;
        case EXPLICIT_RK3:
           // Explicit RK3 method
           higflow_explicit_runge_kutta_3_intermediate_velocity(ns);
           break;
        case SEMI_IMPLICIT_EULER:
           // Semi-Implicit Euler Method
           higflow_semi_implicit_euler_intermediate_velocity(ns);
           break;
        case SEMI_IMPLICIT_CN:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_crank_nicolson_intermediate_velocity(ns);
           break;
        case SEMI_IMPLICIT_BDF2:
           // Semi-Implicit Crank-Nicolson Method
           higflow_semi_implicit_bdf2_intermediate_velocity(ns, ns->dpu, ns->dpustar);
           break;
    }
    // Set outflow for ustar velocity 
    //higflow_outflow_ustar_step(ns);
    // Calculate the pressure
    higflow_pressure(ns);
    // Calculate the final velocity
    higflow_final_velocity(ns);
    // Set outflow for u velocity 
    //higflow_outflow_u_step(ns);
    // Calculate the final pressure
    higflow_final_pressure(ns);
}

// *******************************************************************
// Calculate convective cell term CUBISTA
// *******************************************************************
real hig_flow_convective_cell_term_cubista(distributed_property* dpu, sim_facet_domain* sfdu, sim_stencil* stn, distributed_property* dpK, sim_domain* sdED, sim_stencil* stnED, real kc, Point ccenter, Point cdelta, int dim) {
    real  vbar, vl, vr, kr, krr, kl, kll, a, b, c, d, e, tol, fi, conv1, conv2;
    a = 1.7500;
    b = 0.3750;
    c = 0.7500;
    d = 0.1250;
    e = 0.2500;
    tol = 1.0e-14;
    conv1 = 0.0;
    conv2 = 0.0;
    int   incell_r, incell_l, incell_ll, incell_rr, infacet;

    // Get the low, high, lowlow, highhigh component kernel at center cell
    kl = compute_center_p_left_22(sdED, ccenter, cdelta, dim, 1.0, dpK, stnED, &incell_l);
    kr = compute_center_p_right_22(sdED, ccenter, cdelta, dim, 1.0, dpK, stnED, &incell_r);
    kll = compute_center_p_left_22(sdED, ccenter, cdelta, dim, 2.0, dpK, stnED, &incell_ll);
    krr = compute_center_p_right_22(sdED, ccenter, cdelta, dim, 2.0, dpK, stnED, &incell_rr);

    vl = compute_facet_u_left(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);
    vr = compute_facet_u_right(sfdu, ccenter, cdelta, dim, 0.5, dpu, stn, &infacet);

    // Get the velocity  v1bar(i+1/2,j) in the facet center
    vbar = vr;
    if (vbar > 0.0) {
        if (fabs(kr - kl) <= tol) {
            conv1 = vbar * kc;
        }
        else {
            fi = (kc - kl) / (kr - kl);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv1 = vbar * kc;
            }
            else {
                if (fi < b) {
                    if (incell_l == 1)                    conv1 = vbar * (a * kc - c * kl);
                    else                                  conv1 = vbar * kc;
                }
                if ((fi >= b) && (fi <= c)) {
                    if ((incell_l == 1) && (incell_r == 1)) conv1 = vbar * (c * kc + b * kr - d * kl);
                    else                                  conv1 = vbar * kc;
                }
                if (fi > c) {
                    if (incell_r == 1)                    conv1 = vbar * (e * kc + c * kr);
                    else                                  conv1 = vbar * kc;
                }

            }
        }
        //v1bar < 0.0
    }
    else {
        if ((incell_r == 1) && (incell_rr == 1)) {
            if (fabs(kc - krr) <= tol) {
                conv1 = vbar * kr;
            }
            else {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar * kr;
                }
                else {
                    if (fi < b)
                        conv1 = vbar * (a * kr - c * krr);
                    if ((fi >= b) && (fi <= c))
                        conv1 = vbar * (c * kr + b * kc - d * krr);
                    if (fi > c)
                        conv1 = vbar * (c * kc + e * kr);
                }
            }
        }
        else if ((incell_r == 1) && (incell_rr == 0)) { 
            if (fabs(kc - krr) <= tol) {
                conv1 = vbar * kr;
            }
            else {
                fi = (kr - krr) / (kc - krr);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv1 = vbar * kr;
                }
                else {
                    if (fi <= c)
                        conv1 = vbar * kr;
                    if (fi > c)
                        conv1 = vbar * (c * kc + e * kr);
                }
            }
        }
        else { //Return upwind value at boundary
            vbar = vr;
            if (vbar > 0.0) conv1 = vbar * kc;
            else                 conv1 = vbar * kc;
            vbar = vl;
            if (vbar > 0.0) conv2 = vbar * kl;
            else                 conv2 = vbar * kc;
            return ((conv1 - conv2) / cdelta[dim]);
        }

    }
    // Get the velocity  v2bar(i-1/2,j) in the facet center
    vbar = vl;
    if (vbar > 0.0) {
        if ((incell_l == 1) && (incell_ll == 1)) {
            if (fabs(kc - kll) <= tol) {
                conv2 = vbar * kl;
            }
            else {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv2 = vbar * kl;
                }
                else {
                    if (fi < b)
                        conv2 = vbar * (a * kl - c * kll);
                    if ((fi >= b) && (fi <= c))
                        conv2 = vbar * (b * kc + c * kl - d * kll);
                    if (fi > c)
                        conv2 = vbar * (c * kc + e * kl);
                }
            }
        }
        else if ((incell_l == 1) && (incell_ll == 0)) {
            if (fabs(kc - kll) <= tol) {
                conv2 = vbar * kl;
            }
            else {
                fi = (kl - kll) / (kc - kll);
                if ((fi <= 0.0) || (fi >= 1.0)) {
                    conv2 = vbar * kl;
                }
                else {
                    if (fi <= c)
                        conv2 = vbar * kl;
                    if (fi > c)
                        conv2 = vbar * (c * kc + e * kl);
                }
            }
        }
        else { //Return upwind value at boundary
            vbar = vr;
            if (vbar > 0.0) conv1 = vbar * kc;
            else                 conv1 = vbar * kr;
            vbar = vl;
            if (vbar > 0.0) conv2 = vbar * kc;
            else                 conv2 = vbar * kc;
            return ((conv1 - conv2) / cdelta[dim]);
        }
    }
    else {
        //v2bar < 0.0 
        if (fabs(kl - kr) <= tol) {
            conv2 = vbar * kc;
        }
        else {
            fi = (kc - kr) / (kl - kr);
            if ((fi <= 0.0) || (fi >= 1.0)) {
                conv2 = vbar * kc;
            }
            else {
                if (fi < b) {
                    if (incell_r == 1)                    conv2 = vbar * (a * kc - c * kr);
                    else                                  conv2 = vbar * kc;
                }
                if ((fi >= b) && (fi <= c)) {
                    if ((incell_l == 1) && (incell_r == 1)) conv2 = vbar * (c * kc + b * kl - d * kr);
                    else                                  conv2 = vbar * kc;
                }
                if (fi > c) {
                    if (incell_l == 1)                    conv2 = vbar * (e * kc + c * kl);
                    else                                  conv2 = vbar * kc;
                }
            }
        }
    }
    return ((conv1 - conv2) / cdelta[dim]);
}


