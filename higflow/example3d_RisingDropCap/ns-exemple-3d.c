#include "ns-exemple-3d.h"
#include <string.h>

real func(Point p) {
    real cx = 0.5, cy = 0.5, cz = 0.5, r = 0.25;
    real dx = p[0] - cx, dy = p[1] - cy, dz = p[2] - cz;
    return (dx * dx + dy * dy + dz * dz) / (r * r) - 1.0;
}

int square_case(real f0, real f1, real f10, real f3) { return 0; }

void Intersec(Point p0, Point p1, real f0, real f1, Point p) {
    real lambda = -f0/(f1-f0);
    p[0] = p0[0] + lambda*(p1[0]-p0[0]);
    p[1] = p0[1] + lambda*(p1[1]-p0[1]);
    p[2] = p0[2] + lambda*(p1[2]-p0[2]);
}

real get_fracvol(Point center, Point delta, real t) {
    Point p; int nsub = 8; real vol = 0.0;
    for (int i = 0; i < nsub; i++) {
        p[0] = center[0] - 0.5*delta[0] + (i+0.5)*delta[0]/nsub;
        for (int j = 0; j < nsub; j++) {
            p[1] = center[1] - 0.5*delta[1] + (j+0.5)*delta[1]/nsub;
            for (int k = 0; k < nsub; k++) {
                p[2] = center[2] - 0.5*delta[2] + (k+0.5)*delta[2]/nsub;
                if (func(p) <= 0.0) vol += 1.0;
            }
        }
    }
    return vol/(nsub*nsub*nsub);
}

real get_pressure(Point center, real t) { return 0.0; }
real get_velocity(Point center, int dim, real t) { return 0.0; }
real get_source_term(Point center, real t) { return 0.0; }
real get_facet_source_term(Point center, int dim, real t) { return 0.0; }
real get_boundary_pressure(int id, Point center, real t) { return 0.0; }
real get_boundary_velocity(int id, Point center, int dim, real t) { return 0.0; }
real get_boundary_source_term(int id, Point center, real t) { return 0.0; }
real get_boundary_facet_source_term(int id, Point center, int dim, real t) { return 0.0; }
/* Density/viscosity follow the solver mixing rule
 *   value = (1 - fracvol)*x0 + fracvol*x1
 * so phase 1 (x1) is the region where fracvol == 1, i.e. the drop (the
 * sphere defined by func() <= 0).  The drop is the LIGHT phase that rises;
 * the ambient (fracvol == 0, phase 0) is the heavy phase.  Values are
 * non-dimensionalised by the ambient (OpenFOAM risingDrop_cap: drop rho=1,
 * mu=0.1; ambient rho=1000, mu=10 -> ratios 1/1000 and 1/100). */
real get_density0(Point center, real t) { return 1.0; }    /* ambient (heavy) */
real get_density1(Point center, real t) { return 0.001; }  /* drop (light)    */
real get_viscosity0(Point center, real t) { return 1.0; }  /* ambient         */
real get_viscosity1(Point center, real t) { return 0.01; } /* drop            */

#if ADAPT_ENABLED
real REFINE_THRESHOLDS[] = {0.05, 0.03, -1.0};
#include "mesh_adapt_function.c"

void higflow_interpolate_all_cells_3d(higflow_solver *ns, higflow_solver *ns2) {
    sim_domain *sdp  = psd_get_local_domain(ns->psdp);
    sim_domain *sdm  = psd_get_local_domain(ns->ed.mult.psdmult);
    sim_domain *sdp2 = psd_get_local_domain(ns2->psdp);
    sim_domain *sdm2 = psd_get_local_domain(ns2->ed.mult.psdmult);
    mp_mapper *mpp2 = sd_get_domain_mapper(sdp2);
    mp_mapper *mpm2 = sd_get_domain_mapper(sdm2);
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm2); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        Point ccenter;
        hig_get_center(c, ccenter);
        int clid_m = mp_lookup(mpm2, hig_get_cid(c));
        int clid_p = mp_lookup(mpp2, hig_get_cid(c));
        real fracvol = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpfracvol, ns->ed.mult.stn);
        real visc = compute_value_at_point(sdm, ccenter, ccenter, 1.0, ns->ed.mult.dpvisc, ns->ed.mult.stn);
        real dens0 = ns->ed.mult.get_density0(ccenter, ns->par.t);
        real dens1 = ns->ed.mult.get_density1(ccenter, ns->par.t);
        real dens = (1.0 - fracvol) * dens0 + fracvol * dens1;
        real p = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->dpp, ns->stn);
        dp_set_value(ns2->ed.mult.dpfracvol, clid_m, fracvol);
        dp_set_value(ns2->ed.mult.dpvisc,   clid_m, visc);
        dp_set_value(ns2->ed.mult.dpdens,   clid_m, dens);
        dp_set_value(ns2->dpp,              clid_p, p);
    }
    higcit_destroy(it);
    dp_sync(ns2->ed.mult.dpfracvol);
    dp_sync(ns2->ed.mult.dpvisc);
    dp_sync(ns2->ed.mult.dpdens);
    dp_sync(ns2->dpp);
}

void higflow_interpolate_velocity_3d(higflow_solver *ns, higflow_solver *ns2) {
    for (int dim = 0; dim < DIM; dim++) {
        sim_facet_domain *sfdu  = psfd_get_local_domain(ns->psfdu[dim]);
        sim_facet_domain *sfdu2 = psfd_get_local_domain(ns2->psfdu[dim]);
        mp_mapper *mu2 = sfd_get_domain_mapper(sfdu2);
        higfit_facetiterator *fit;
        for (fit = sfd_get_domain_facetiterator(sfdu2); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu2, hig_get_fid(f));
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            stn_reset(ns->stn);
            sfd_get_stencil(sfdu, fcenter, fcenter, 1, ns->stn);
            real u = dp_interpolate_from_stencil(ns->dpu[dim], ns->stn);
            dp_set_value(ns2->dpu[dim], flid, u);
        }
        higfit_destroy(fit);
        dp_sync(ns2->dpu[dim]);
    }
}
#endif

int main(int argc, char *argv[]) {

    START_CLOCK(total);
    int ntasks, myrank;
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    higflow_solver *ns = higflow_create();
    memset(ns, 0, sizeof(higflow_solver));
    higflow_load_data_file_names(argc, argv, ns);
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
    higflow_set_external_functions(ns, get_pressure, get_velocity,
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term);
    int order_center = 2, order_facet = 2, cache = 1;
    higflow_create_domain(ns, cache, order_center);
    higflow_create_domain_multiphase(ns, cache, order_center,
        get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet);
    higflow_initialize_boundaries_yaml(ns);


#if ADAPT_ENABLED
    if (ns->par.step == 0) {
        int cache2 = 1;
        higflow_solver *ns2 = higflow_create();
        higflow_load_data_file_names(argc, argv, ns2);
        higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
        higflow_set_external_functions(ns2, get_pressure, get_velocity,
            get_source_term, get_facet_source_term,
            get_boundary_pressure, get_boundary_velocity,
            get_boundary_source_term, get_boundary_facet_source_term);
        higflow_create_domain(ns2, cache2, order_center);
        higflow_create_domain_multiphase(ns2, cache2, order_center,
            get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
        hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);
        partition_graph *pg = pg_create(MPI_COMM_WORLD);
        pg_set_fringe_size(pg, 5);
        sd_add_higtree(ns2->sdp, root);
        sd_add_higtree(ns2->sdF, root);
        sd_add_higtree(ns2->ed.mult.sdmult, root);
        higflow_create_partitioned_domain(ns2, pg, order_center);
        higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);
        higflow_create_stencil(ns2);
        higflow_create_stencil_multiphase(ns2);
        ns2->par = ns->par; ns2->contr = ns->contr;
        ns = ns2;
        printf("===> Initial analytic interface AMR applied\n");
    }
#endif

    higflow_create_distributed_properties(ns);

    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);

    higflow_create_solver(ns);

    if (ns->par.step > 0) higflow_load_properties(ns, myrank, ntasks);

    higflow_save_domain_yaml(ns, myrank, ntasks);
    higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
    if (ns->par.step == 0) {
        higflow_print_vtk(ns, myrank);
        ns->par.tp += ns->par.dtp; ns->par.frame++;
        higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }

    DECL_CLOCK(firstiter);
    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
        if (ns->par.step == step0) START_CLOCK(firstiter);
        higflow_solver_step_multiphase(ns);
        ns->par.t += ns->par.dt;
        if (ns->par.step == step0) STOP_CLOCK(firstiter);

#if ADAPT_ENABLED
        if (ns->par.step % ADAPT_FREQ == 0) {
            int cache2 = 1;
            higflow_solver *ns2 = higflow_create();
            ns2->par = ns->par; ns2->contr = ns->contr;
            ns2->sdp = ns->sdp; ns2->sdF = ns->sdF;
            higflow_set_external_functions(ns2, get_pressure, get_velocity,
                get_source_term, get_facet_source_term,
                get_boundary_pressure, get_boundary_velocity,
                get_boundary_source_term, get_boundary_facet_source_term);
            higflow_create_domain(ns2, cache2, order_center);
            higflow_create_domain_multiphase(ns2, cache2, order_center,
                get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
            hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);
            partition_graph *pg = pg_create(MPI_COMM_WORLD);
            pg_set_fringe_size(pg, 5);
            sd_add_higtree(ns2->sdp, root);
            sd_add_higtree(ns2->sdF, root);
            sd_add_higtree(ns2->ed.mult.sdmult, root);
            higflow_create_partitioned_domain(ns2, pg, order_center);
            higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);
            higflow_create_stencil(ns2);
            higflow_create_stencil_multiphase(ns2);
            higflow_create_distributed_properties(ns2);
            higflow_initialize_boundaries_yaml(ns2);
            higflow_interpolate_velocity_3d(ns, ns2);
            higflow_interpolate_all_cells_3d(ns, ns2);
            higflow_create_solver(ns2);
            ns2->par = ns->par; ns2->contr = ns->contr;
            higflow_destroy(ns);
            ns = ns2;
            psd_synced_mapper(ns->psdp);
            for (int dim = 0; dim < DIM; dim++)
                psfd_synced_mapper(ns->psfdu[dim]);
        }
#endif

        if (ns->par.t >= ns->par.tp) {
            higflow_print_vtk(ns, myrank);
            ns->par.tp += ns->par.dtp; ns->par.frame++;
        }
        if (ns->par.t >= ns->par.ts) {
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            ns->par.ts += ns->par.dts;
        }
    }
    higflow_destroy(ns);
    STOP_CLOCK(total);
    return 0;
}
