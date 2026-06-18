// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/2023
// *******************************************************************
// *******************************************************************

#include "ns-example-2d.h"
#include <string.h>

// -----------------------------------------------------------------------
// BC higtree refinement for AMR cases.
// Registered via higflow_set_bc_refine_hook() so each boundary higtree
// is refined in place to match the adjacent internal mesh before the
// sim_boundary is created.  _bc_domain_root is set just before calling
// higflow_initialize_boundaries_yaml() and cleared immediately after.
// -----------------------------------------------------------------------
static hig_cell *_bc_domain_root = NULL;

static int _bc_level(hig_cell *c) {
    int l = 0;
    while (c) { c = hig_get_parent(c); l++; }
    return l;
}

static void _refine_bc_tree(hig_cell *bc_root, int bc_id) {
    if (!_bc_domain_root) return;
    const real eps = 1e-7;
    bool changed = true;
    while (changed) {
        changed = false;
        higcit_celliterator *it = higcit_create_all_leaves(bc_root);
        for (; !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *bc_leaf = higcit_getcell(it);
            Point center; hig_get_center(bc_leaf, center);
            Point q; q[0] = center[0]; q[1] = center[1];
            if      (bc_id == 0) q[0] += eps;
            else if (bc_id == 1) q[1] -= eps;
            else if (bc_id == 2) q[0] -= eps;
            else                 q[1] += eps;
            hig_cell *dom = hig_get_cell_with_point(_bc_domain_root, q);
            if (!dom) continue;
            int bc_lev  = _bc_level(bc_leaf);
            int dom_lev = _bc_level(dom);
            if (dom_lev > bc_lev) {
                int nc[DIM];
                if (bc_id == 0 || bc_id == 2) { nc[0]=1; nc[1]=2; }
                else                            { nc[0]=2; nc[1]=1; }
                hig_refine_uniform(bc_leaf, nc);
                changed = true;
                higcit_destroy(it);
                break;
            }
        }
        if (!changed) higcit_destroy(it);
    }
}

/************************************ user functions **************************************/

#include "ns-user-functions-vof.c"

#include "ns-user-functions-electroosmotic.c"

#include "ns-user-functions-viscoelastic.c"
#include "ns-user-functions-newtonian-gn.c"
#include <stdlib.h>

/******************************************************************************************/
/******************************************************************************************/
/********************************* main user functions ************************************/
/******************************************************************************************/
/******************************************************************************************/


void create_initialize_all_domains(higflow_solver* ns, int myrank, int ntasks) {
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 1;
    int order_facet = 1;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;

    // Set Initial conditions and Boundary Conditions from user defined functions
    higflow_create_domain(ns, cache, order_center);
    higflow_set_external_functions(ns, get_pressure, get_velocity,
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term);

    switch (ns->contr.flowtype) {
        case GENERALIZED_NEWTONIAN:
            higflow_create_domain_generalized_newtonian(ns, cache, order_center, 
                                                        get_viscosity_gn);
            break;
        case MULTIPHASE:
            higflow_create_domain_multiphase(ns, cache, order_center, get_viscosity0, get_viscosity1, 
                                             get_density0, get_density1, get_fracvol);
            if(ns->ed.mult.contr.viscoelastic_either == true) {
                higflow_create_domain_multiphase_viscoelastic(ns, get_tensor_multiphase, get_kernel,
                                                             get_kernel_inverse, get_kernel_jacobian);
                higflow_define_user_function_multiphase_viscoelastic(ns, calculate_m_user_multiphase);
            }
            if(ns->ed.mult.contr.eoflow_either == true) 
                higflow_create_domain_multiphase_electroosmotic(ns, cache, order_center, get_multiphase_electroosmotic_source_term,
                                                               get_multiphase_electroosmotic_phi, get_multiphase_electroosmotic_psi, 
                                                               get_multiphase_electroosmotic_nplus, get_multiphase_electroosmotic_nminus, 
                                                               get_boundary_multiphase_electroosmotic_source_term, 
                                                               get_boundary_multiphase_electroosmotic_phi, get_boundary_multiphase_electroosmotic_psi, 
                                                               get_boundary_multiphase_electroosmotic_nplus, get_boundary_multiphase_electroosmotic_nminus,
                                                               get_multiphase_electroosmotic_permittivity);
            break;
        case VISCOELASTIC:
            higflow_create_domain_viscoelastic(ns, cache, order_center, get_tensor, get_kernel,
                                               get_kernel_inverse, get_kernel_jacobian);
            higflow_define_user_function_viscoelastic(ns, calculate_m_user);
            break;
        case VISCOELASTIC_INTEGRAL:
            higflow_create_domain_viscoelastic_integral(ns, cache, order_center, get_tensor_integral);
            break;
    }
        
    if (ns->contr.eoflow == true)
        higflow_create_domain_electroosmotic(ns, cache, order_center, get_electroosmotic_source_term,
                                             get_electroosmotic_phi, get_electroosmotic_psi, 
                                             get_electroosmotic_nplus, get_electroosmotic_nminus, 
                                             get_boundary_electroosmotic_source_term, 
                                             get_boundary_electroosmotic_phi, get_boundary_electroosmotic_psi, 
                                             get_boundary_electroosmotic_nplus, get_boundary_electroosmotic_nminus,
                                             get_electroosmotic_permittivity);
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet);
}

void solver_step(higflow_solver* ns) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        switch (ns->contr.flowtype) {
            case NEWTONIAN:
                higflow_solver_step_electroosmotic(ns);
                break;
            case MULTIPHASE:
                if(ns->ed.mult.contr.viscoelastic_either == true)
                    higflow_solver_step_multiphase_electroosmotic_viscoelastic(ns);
                else
                    higflow_solver_step_multiphase_electroosmotic(ns);
                break;
            case VISCOELASTIC:
                higflow_solver_step_electroosmotic_viscoelastic(ns);
                break;
        }
    }
    else {
        switch (ns->contr.flowtype) {
            case NEWTONIAN:
                higflow_solver_step(ns);
                break;
            case MULTIPHASE:
                if(ns->ed.mult.contr.viscoelastic_either == true) 
                    higflow_solver_step_multiphase_viscoelastic(ns);
                else
                    higflow_solver_step_multiphase(ns);
                ns->par.stepaux=ns->par.stepaux+1;
                break;
            case GENERALIZED_NEWTONIAN:
                higflow_solver_step_gen_newt(ns);
                break;
            case VISCOELASTIC:
                higflow_solver_step_viscoelastic(ns);
                break;
            case VISCOELASTIC_INTEGRAL:
                higflow_solver_step_viscoelastic_integral(ns);
                break;
        }
    }
}

int errors(higflow_solver* ns, sim_residuals* sim_res, int myrank) {
    int errcode = 0;

    if (sim_res != NULL) {
        real dudt_norm = 0.0;
        write_residuals(sim_res, ns);
        if (myrank == 0) dudt_norm = sim_res->u[0]->midrange->res_max->avg[0] / ns->par.dt;

        MPI_Bcast(&dudt_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        print0f("|    dudt_norm = %15.10lf", dudt_norm);
        if(flowtype != MULTIPHASE) {
            real tol = 1.0e-5;
            if (ns->contr.eoflow == true) tol = 5.0e-5;
            if (dudt_norm < max(1.0e-10 / ns->par.dt, tol)) {
                print0f("\nsteady state reached\n");
                errcode = 1;
            }
            if (dudt_norm > 1.0e8) {
                print0f("\nsimulation 'diverged'\n");
                errcode = -1;
            }
        }
    }

    return errcode;
}

void save(higflow_solver* ns, int myrank, int ntasks) {
    if (myrank == 0) printf("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
    //higflow_save_properties(ns, myrank, ntasks);
    ns->par.ts += ns->par.dts;
}

void print(higflow_solver* ns, int myrank, int ntasks){
    if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
            printf("===> Printing frame: %4d <====> tp = %15.10lf <===\n", ns->par.frame, ns->par.tp);
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
    }
    
    higflow_print_vtk(ns, myrank);
    // if(ns->contr.flowtype == MULTIPHASE) higflow_print_vtk2D_multiphase(ns, myrank);
    // higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
    // write_xdmf(ns);
    if(ns->contr.flowtype == MULTIPHASE) {
        //higflow_print_vtk2D_multiphase_parallel_single(ns, myrank, ntasks);
        if(ns->par.step==0) {
            higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns);
            higflow_compute_distance_multiphase_2D(ns);
            higflow_compute_plic_lines_2d(ns);
        }
        save_deformation_parameter(ns, myrank);
        higflow_print_vtk2d_multiphase_plic_lines_serial_single(ns, myrank, ntasks);
    }
    ns->par.frame++;
    ns->par.tp += ns->par.dtp;
}

void init_global_var(higflow_solver* ns) {
    flowtype = ns->contr.flowtype;
    sdp_ptr = &ns->sdp;
    stn_ptr = &ns->stn;
    dpp_ptr = &ns->dpp;
    sfdv_ptr = &(ns->sfdu[1]);
    dpvstar_ptr = &ns->dpustar[1];
    if(flowtype == VISCOELASTIC) visc_model = ns->ed.ve.contr.model;
    if(flowtype == MULTIPHASE) {
        dpfracvol_ptr = &ns->ed.mult.dpfracvol;
        sdmult_ptr = &ns->ed.mult.sdmult;
        stnmult_ptr = &ns->ed.mult.stn;
        flowtype0 = ns->ed.mult.contr.flowtype0;
        flowtype1 = ns->ed.mult.contr.flowtype1;
        viscoelastic_either = ns->ed.mult.contr.viscoelastic_either;
        eoflow0 = ns->ed.mult.contr.eoflow0;
        eoflow1 = ns->ed.mult.contr.eoflow1;
        eoflow_either = ns->ed.mult.contr.eoflow_either;
    }
    eoflow = ns->contr.eoflow;
    if(eoflow == true) {
        sfdFeoy_ptr = &ns->ed.eo.sfdEOFeo[1];
        dpFeoy_ptr = &ns->ed.eo.dpFeo[1];
        stnFeoy_ptr = &ns->ed.eo.stnpsi;
    }
}

void get_inlet_types(higflow_solver* ns) {
    char namefile[1024];
    sprintf(namefile,"%s.bc.yaml",ns->par.nameload);
    
    FILE *fbc = fopen(namefile, "r");
    struct fy_document *fyd = NULL;
    fyd = fy_document_build_from_file(NULL, namefile);
     
    if (fyd == NULL) {
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }

    char aux[1024];
    int ifd = fy_document_scanf(fyd,"/bc/bc0/velocity_0/type %s",aux);
    if (strcmp(aux,"dirichlet") == 0) u_inlet = DIRICHLET;
    else if (strcmp(aux,"neumann") == 0) u_inlet = NEUMANN;
    else {
        printf("=+=+=+= Error loading boundary condition type for the inlet velocity\n");
        exit(1);
    }

    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        ifd = fy_document_scanf(fyd,"/bc_electroosmotic/bc0/psi/type %s",aux);
        if (strcmp(aux,"dirichlet") == 0) psi_inlet = DIRICHLET;
        else if (strcmp(aux,"neumann") == 0) psi_inlet = NEUMANN;
        else {
            printf("=+=+=+= Error loading boundary condition type for the inlet velocity\n");
            exit(1);
        }
    }
    
    fy_document_destroy(fyd);
    fclose(fbc);
}

/******************************************************************************************/
/******************************************************************************************/
/************************************** main program **************************************/
/******************************************************************************************/
/******************************************************************************************/

#if ADAPT_ENABLED
// Refinement thresholds (sentinel -1.0 marks end of array)
real REFINE_THRESHOLDS[] = {0.05, 0.03, -1.0};

#include "mesh_adapt_function.c"
#endif

real compute_total_fracvol(higflow_solver *ns) {
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    real total = 0.0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        Point delta;
        hig_get_delta(c, delta);
        real volcell = delta[0] * delta[1];
        real fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
        total += fracvol * volcell;
    }
    higcit_destroy(it);
    real global_total;
    MPI_Allreduce(&total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_total;
}

// Interpolate all cell-centered properties in a single iterator pass
void higflow_interpolate_all_cells(higflow_solver *ns, higflow_solver *ns2) {
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

        real fracvol = compute_value_at_point(sdm, ccenter, ccenter, 1.0,
                                               ns->ed.mult.dpfracvol, ns->ed.mult.stn);

        if (ns->par.step == 5 || ns->par.step == 10) {
            real val = 0.4;
            if (fracvol >= (1.0 - val))      fracvol = 1.0;
            else if (fracvol <= val)         fracvol = 0.0;
        }

        real visc = compute_value_at_point(sdm, ccenter, ccenter, 1.0,
                                            ns->ed.mult.dpvisc, ns->ed.mult.stn);

        real dens0 = ns->ed.mult.get_density0(ccenter, ns->par.t);
        real dens1 = ns->ed.mult.get_density1(ccenter, ns->par.t);
        real dens  = (1.0 - fracvol) * dens0 + fracvol * dens1;

        real p = compute_value_at_point(sdp, ccenter, ccenter, 1.0,
                                         ns->dpp, ns->stn);

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

// Interpola velocidade (facets) da malha velha para a nova
void higflow_interpolate_velocity(higflow_solver *ns, higflow_solver *ns2) {
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

void higflow_interpolate_all_bcs(higflow_solver *ns, higflow_solver *ns2) {
    sim_domain *sdp  = psd_get_local_domain(ns->psdp);
    sim_domain *sdp2 = psd_get_local_domain(ns2->psdp);
    int num_bc_types = 2;
    bc_type bc_t;

    for (int i = 0; i < num_bc_types; i++) {
        if (i == 0) bc_t = DIRICHLET; else bc_t = NEUMANN;
        int numbcs = sd_get_num_bcs(sdp2, bc_t);
        for (int h = 0; h < numbcs; h++) {
            sim_boundary *bc2 = sd_get_bc(sdp2, bc_t, h);
            if (bc_t == NEUMANN || sb_get_valuetype(bc2) == fixedValue) continue;
            mp_mapper *bm2 = sb_get_mapper(bc2);
            higcit_celliterator *it;
            for (it = sb_get_celliterator(bc2); !higcit_isfinished(it); higcit_nextcell(it)) {
                hig_cell *bcell = higcit_getcell(it);
                Point bc; hig_get_center(bcell, bc);
                int bclid = mp_lookup(bm2, hig_get_cid(bcell));
                stn_reset(ns->stn);
                sd_get_stencil(sdp, bc, bc, 1, ns->stn);
                sb_set_value(bc2, bclid, compute_value_at_point(sdp, bc, bc, 1.0, ns->dpp, ns->stn));
            }
            higcit_destroy(it);
        }
    }

    for (int dim = 0; dim < DIM; dim++) {
        sim_facet_domain *sfdu2 = psfd_get_local_domain(ns2->psfdu[dim]);
        sim_domain *sd2 = sfdu2->cdom;
        for (int i = 0; i < num_bc_types; i++) {
            if (i == 0) bc_t = DIRICHLET; else bc_t = NEUMANN;
            int numbcs2 = sd_get_num_bcs(sd2, bc_t);
            for (int h = 0; h < numbcs2; h++) {
                sim_boundary *bc2 = sd_get_bc(sd2, bc_t, h);
                int uid = sb_get_userid(bc2);
                mp_mapper *bm2 = sb_get_mapper(bc2);
                higcit_celliterator *it;
                for (it = sb_get_celliterator(bc2); !higcit_isfinished(it); higcit_nextcell(it)) {
                    hig_cell *bcell = higcit_getcell(it);
                    Point bc; hig_get_center(bcell, bc);
                    int bclid = mp_lookup(bm2, hig_get_cid(bcell));
                    sb_set_value(bc2, bclid, ns2->func.get_boundary_velocity(uid, bc, dim, ns2->par.t + ns2->par.dt));
                }
                higcit_destroy(it);
            }
        }
    }
}


int main(int argc, char* argv[]) {
    int errcode = 0;
    START_CLOCK(total);
    int ntasks; // Number of tasks
    int myrank; // Identifier of the process
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 1;
    int order_facet = 1;
    // Set the cache: Reuse interpolation, 0 on, 1 off
    int cache = 1;
    printf("=+=+ Initializing Navier-Stokes Solver =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    higflow_solver* ns = higflow_create();
    memset(ns, 0, sizeof(higflow_solver));

    // Set data file type names
    higflow_load_data_file_names(argc, argv, ns);

    print0f("=+=+ Loading Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);

    p_par = create_initialize_physical_parameters(ns, myrank);
    init_global_var(ns);

    print0f("=+=+ Loading, Creating, Partitioning and Initializing Domains =+=+=+=+=+=\n");
    create_initialize_all_domains(ns, myrank, ntasks);

    print0f("=+=+ Creating and Initializing Distributed Properties =+=+=+=+=+=+=+=+=+=+=+=+=\n");

    // =====================================================
    // INITIAL ADAPTATION BASED ON ANALYTICAL INTERFACE
    // =====================================================
#if ADAPT_ENABLED
    if (ns->par.step == 0) {
        int order_center = 1;
        int cache = 1;

        higflow_solver *ns2 = higflow_create();
        higflow_load_data_file_names(argc, argv, ns2);
        higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
        higflow_set_external_functions(ns2,
            get_pressure, get_velocity,
            get_source_term, get_facet_source_term,
            get_boundary_pressure, get_boundary_velocity,
            get_boundary_source_term, get_boundary_facet_source_term);

        higflow_create_domain(ns2, cache, order_center);
        higflow_create_domain_multiphase(ns2, cache, order_center,
            get_viscosity0, get_viscosity1,
            get_density0, get_density1, get_fracvol);

        hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);

        sd_add_higtree(ns2->sdp, root);
        sd_add_higtree(ns2->sdF, root);
        sd_add_higtree(ns2->ed.mult.sdmult, root);

        partition_graph *pg = pg_create(MPI_COMM_WORLD);
        pg_set_fringe_size(pg, 5);
        load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
        lb_destroy(lb);

        higflow_create_partitioned_domain(ns2, pg, order_center);
        higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);

        higflow_create_stencil(ns2);
        higflow_create_stencil_multiphase(ns2);

        ns2->par = ns->par;
        ns2->contr = ns->contr;
        ns = ns2;

        print0f("===> Initial analytic interface AMR applied\n");
    }
#endif

    higflow_create_distributed_properties(ns);
    higflow_print_vtk(ns, myrank);
    if (ns->par.step == 0) {
        higflow_initialize_distributed_properties(ns);
        real v0 = compute_total_fracvol(ns);
        print0f("=+= Total volume at step 0 = %lf =+=\n", v0);
    }


    // get inlet boundary types to set boundary conditions correctly
    get_inlet_types(ns);

    print0f("=+=+ Initializing Boundaries =+=+\n");
#if ADAPT_ENABLED
    _bc_domain_root = sd_get_higtree(psd_get_local_domain(ns->psdp), 0);
    higflow_set_bc_refine_hook(_refine_bc_tree);
#endif
    higflow_initialize_boundaries_yaml(ns);
#if ADAPT_ENABLED
    higflow_set_bc_refine_hook(NULL);
    _bc_domain_root = NULL;
#endif

    if(ns->par.step > 0){
        if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
            printf("===> Reloading properties from previous simulation <====> step = %d <====> t = %15.10lf <===\n", ns->par.step, ns->par.t);
            printf("*********************************************************************************\n");
            printf("*********************************************************************************\n");
        }
        //higflow_load_properties(ns, myrank, ntasks);
    }

    print0f("=+=+ Creating Linear System Solvers =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_create_solver(ns);

    MPI_Barrier(MPI_COMM_WORLD);
    print0f("=+=+ Saving Domain and Boundary Properties =+=+\n");
    higflow_save_domain_yaml(ns, myrank, ntasks);
    higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank); //copying necessary yamls

    write_init(ns);
    
    if (ns->par.step == 0) {
        print(ns, myrank, ntasks);
        save(ns, myrank, ntasks);
    }
    // =====================================================
    // INITIAL ADAPTATION BASED ON ANALYTICAL INTERFACE
    // =====================================================

    sim_residuals *sim_res = create_initialize_sim_residuals(ns);
    real u_center;

    // ********************************************************
    // Begin Loop for the Navier-Stokes equations integration
    // ********************************************************
    for (; ns->par.step <= ns->par.finalstep; ns->par.step++) {
        //////////////////////////////////////////// IO ///////////////////////////////////////////////
        GET_NSEC_CLOCK(iter_total) = 0.0; START_CLOCK(iter_total);

        if (FLT_GE(ns->par.t, ns->par.tp))
            if(ns->par.step - ns->par.initstep > 0) print(ns, myrank, ntasks);
            // higflow_print_vtk(ns, myrank);

        if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf(" ===> Step:   %7d <====> t   = %15.10lf <===\n", ns->par.step, ns->par.t);
            printf("*********************************************************************************\n");
        }
        GET_NSEC_CLOCK(currentiter) = 0.0; START_CLOCK(currentiter);

        ///////////////////////////////////////////////////////
        solver_step(ns);

#if ADAPT_ENABLED
        if (ns->par.step % ADAPT_FREQ == 0) {
              // real vol_before = compute_total_fracvol(ns);
              // print0f("=+= Volume before interpolation = %16.10lf =+=\n", vol_before);

              // Initializing Navier-Stokes solver
              // Create Navier-Stokes solver
              higflow_solver *ns2 = higflow_create();

              // Inherit parameters and controllers from the old solver (avoids re-reading YAML)
              ns2->par = ns->par;
              ns2->contr = ns->contr;
              ns2->sdp = ns->sdp;
              ns2->sdF = ns->sdF;

              // set the external functions
              higflow_set_external_functions(ns2, get_pressure, get_velocity, 
                                            get_source_term, get_facet_source_term,
                                            get_boundary_pressure, get_boundary_velocity,
                                            get_boundary_source_term, get_boundary_facet_source_term); 
              // Reset simulation domain
              higflow_create_domain(ns2, cache, order_center); 

              // // case MULTIPHASE:
              higflow_create_domain_multiphase(ns2, cache, order_center, get_viscosity0, get_viscosity1, 
                                              get_density0, get_density1, get_fracvol);

              // // Initialize the domain
              // print0f("=+=+=+= Load Domain (ns2) =+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
              //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
              hig_cell *root = higflow_make_adapted_tree_params(ns, REFINE_THRESHOLDS);

              partition_graph *pg = pg_create(MPI_COMM_WORLD);
              pg_set_fringe_size(pg, 5);
              load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
              lb_destroy(lb);

              sd_add_higtree(ns2->sdp, root);
              sd_add_higtree(ns2->sdF, root);
              sd_add_higtree(ns2->ed.mult.sdmult, root);

              // // Creating the partitioned sub-domain to simulation
              higflow_create_partitioned_domain(ns2, pg, order_center);
              higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);

              // Creating the stencil for properties interpolation
              higflow_create_stencil(ns2);
              higflow_create_stencil_multiphase(ns2);

              // Creating distributed property  
              print0f("=+=+=+= Creating distributed property (ns2) +=+=+=+=+=\n");
              higflow_create_distributed_properties(ns2);

              // Refine BC higtrees to match the adapted internal mesh before building boundaries
              _bc_domain_root = root;
              higflow_set_bc_refine_hook(_refine_bc_tree);
              higflow_initialize_boundaries_yaml(ns2);
              higflow_set_bc_refine_hook(NULL);
              _bc_domain_root = NULL;

              // Interpolar (3 chamadas fundidas em vez de 6+)
              print0f("=+=+=+= Interpolation (ns2) +=+=+=+=+=\n");
              higflow_interpolate_velocity(ns, ns2);
              higflow_interpolate_all_cells(ns, ns2);

              higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns2);
              higflow_compute_distance_multiphase_2D(ns2);
              higflow_compute_plic_lines_2d(ns2);

              higflow_interpolate_all_bcs(ns, ns2);

              //{
              //    real vol_after = compute_total_fracvol(ns2);
              //    print0f("=+= Volume after interpolation  = %16.10lf =+=\n", vol_after);
              //    print0f("=+= Volume change               = %16.10lf =+=\n", vol_after - vol_before);
              //}

              higflow_create_solver(ns2); 

              // Copy essential parameters from the old solver, then destroy it
              ns2->par = ns->par;
              ns2->contr = ns->contr;
              higflow_destroy(ns);
              ns = (higflow_solver *) ns2;

              // Sync mappers so PETSc global IDs are consistent
              psd_synced_mapper(ns->psdp); 
              for(int dim = 0; dim < DIM; dim++) {
                psfd_synced_mapper(ns->psfdu[dim]); 
              }

              // higflow_print_vtk(ns, myrank);
        }
#endif
        /////////////////////////////////////////////////

        write_mem_usage(ns, "after step");

        ////////////////// errors //////////////////////
        u_center = get_fdp_value_at_point(ns, ns->dpu[0], ns->psfdu[0], (Point) { p_par->center[0].val, p_par->center[1].val });
        print0f("u_center = %15.10lf    ", u_center);
        errcode = errors(ns, sim_res, myrank);
        print0f("\n");

        /////// check if A is SPD ////////////////////////////////////////////
        real max_neg_lambda_global; int num_neg_lambda_global;
        MPI_Allreduce(&max_neg_lambda, &max_neg_lambda_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&num_neg_lambda, &num_neg_lambda_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (num_neg_lambda_global > 0) print0f("Warning! %d neg eigenvalues <====> max neg = %15.10lf \n", num_neg_lambda_global, max_neg_lambda_global);

        /////////////////////////////////////////////////////////
        // Time update 
        ns->par.t += ns->par.dt;
        STOP_CLOCK(currentiter);
        print0f("------------------------ time of iteration = %lf s -----------------------\n", GET_NSEC_CLOCK(currentiter) / 1.0e9);

        if (FLT_GE(ns->par.t, ns->par.ts))
            save(ns, myrank, ntasks);

        STOP_CLOCK(iter_total);
        if (FLT_GE(ns->par.t, ns->par.ts) || FLT_GE(ns->par.t, ns->par.tp))
            print0f("------------------------ time of IO in iteration = %lf s -----------------------\n", (GET_NSEC_CLOCK(iter_total) - GET_NSEC_CLOCK(currentiter)) / 1.0e9);

        if (errcode != 0) break;
    }
    if(errcode > 0) {
        print(ns, myrank, ntasks);
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    free_sim_residuals(sim_res);
    STOP_CLOCK(total);
    print0f("------------------------ total time = %lf s -----------------------\n", GET_NSEC_CLOCK(total) / 1.0e9);
}

