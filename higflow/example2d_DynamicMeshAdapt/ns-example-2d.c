// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/2023
// *******************************************************************
// *******************************************************************

#include "ns-example-2d.h"

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

// Função para criar um snapshot da malha refinada sem afetar a simulação
// Certifique-se de que estes includes estão no topo do arquivo ns-example-2d.c

#include "mesh_adapt_function.c"

// Navier-Stokes final pressure using the projection method
void higflow_interpolate_pressure(higflow_solver *ns, higflow_solver *ns2) {
  // Incremental projection method
  // Get the local sub-domain
  sim_domain *sdp = psd_get_local_domain(ns->psdp);
  sim_domain *sdp2 = psd_get_local_domain(ns2->psdp);
  // Get the map of the distributd properties in the cells
  mp_mapper  *mp2  = sd_get_domain_mapper(sdp2);
  // Loop for each cell
  higcit_celliterator *it;
  for(it = sd_get_domain_celliterator(sdp2); !higcit_isfinished(it); higcit_nextcell(it)) {
    // Get the cell
    hig_cell *c = higcit_getcell(it);
    // Get the cell identifier
    int clid    = mp_lookup(mp2, hig_get_cid(c));
    // Get the cell center
    Point ccenter;
    hig_get_center(c, ccenter);
    // Get the pressure in the distributed pressure property
    real interpolated_p;
    // compute value of p in a point
    interpolated_p = compute_value_at_point(sdp, ccenter,
                                            ccenter, 1.0,
                                            ns->dpp, ns->stn);

    // Coloca o valor interplolado na malha 2
    dp_set_value(ns2->dpp, clid, interpolated_p);
  }
  // Destroy the iterator
  higcit_destroy(it);

  // Sync the distributed pressure property
  dp_sync(ns2->dpp);
}

// Navier-Stokes final pressure using the projection method
void higflow_interpolate_velocity(higflow_solver *ns, higflow_solver *ns2) {
  // Get the local sub-domain
  sim_facet_domain *sfdu[DIM];
  sim_facet_domain *sfdu2[DIM];
  // Loop for each dimension
  higfit_facetiterator *fit;
  for(int dim = 0; dim < DIM; dim++) {
    // Get the local partitioned domain for facets
    sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
    sfdu2[dim] = psfd_get_local_domain(ns2->psfdu[dim]);
    // Get the map of the distributd properties in the facets
    mp_mapper *mu2 = sfd_get_domain_mapper(sfdu2[dim]);

    // Loop for each facet
    for (fit = sfd_get_domain_facetiterator(sfdu2[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
      // Get the facet
      hig_facet *f = higfit_getfacet(fit);
      int flid = mp_lookup(mu2, hig_get_fid(f));
      // Get the center of the facet
      Point fcenter;
      hig_get_facet_center(f, fcenter);
      real interpolated_u;
      stn_reset(ns->stn);
      // Get the stencil parameters
      sfd_get_stencil(sfdu[dim], fcenter, fcenter, 1, ns->stn);
      interpolated_u = dp_interpolate_from_stencil(ns->dpu[dim], ns->stn);
      dp_set_value(ns2->dpu[dim], flid, interpolated_u);
    }
    // Destroy the iterator
    higfit_destroy(fit);
    // Sync the distributed velocity property
    dp_sync(ns2->dpu[dim]);
  }
}

int signn(double x) {
    return (x >= 0) - (x <= 0);
}

void higflow_interpolate_viscosity(higflow_solver *ns, higflow_solver *ns2) {
    // Get the local sub-domain for the cells (SOURCE - Old Mesh)
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    
    // Get the local sub-domain for the cells (TARGET - New Mesh)
    sim_domain *sdm2 = psd_get_local_domain(ns2->ed.mult.psdmult);
    
    // Get the map for the domain properties (Target)
    mp_mapper *mp2 = sd_get_domain_mapper(sdm2);
    
    // Loop for each cell in the NEW domain (ns2)
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm2); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *c = higcit_getcell(it);
        
        // Get the cell identifier in the new mapper
        int clid = mp_lookup(mp2, hig_get_cid(c));
        
        // Get the center of the cell
        Point ccenter;
        hig_get_center(c, ccenter);
        
        // Interpolate fracvol from the OLD solver (ns)
        // Usamos o domínio 'sdm' (velho) e a propriedade 'dpfracvol' do 'ns'
        // (velho)
        real fracvol = compute_value_at_point(sdm, ccenter,
                                              ccenter, 1.0,
                                              ns->ed.mult.dpfracvol,
                                              ns->ed.mult.stn);

        // Valor definido empiricamente (altamente testado)
        if (ns->par.step == 5 || ns->par.step == 10 || ns->par.step == 15){
          printf("Step %d\n", ns->par.step);
          real val = 0.4;
          // Operador ternário para calcular a fração de volume sharp
          // fracvol = (fracvol > 1 - val) ? 1 : (fracvol < val ? 0 : fracvol);
          // if (fracvol > (1 - 0.4)) {
          //     fracvol = 1;
          // } 
          // else if (fracvol <= 0.4) {
          //     fracvol = 0;
          // }
        }

        real visc = compute_value_at_point(sdm, ccenter,
                                              ccenter, 1.0,
                                              ns->ed.mult.dpvisc,
                                              ns->ed.mult.stn);
        
        // Set the viscosity in the distributed viscosity property of the NEW solver (ns2)
        dp_set_value(ns2->ed.mult.dpvisc, clid, visc);
        
        // É importante atualizar também o fracvol no novo solver para manter consistência
        dp_set_value(ns2->ed.mult.dpfracvol, clid, fracvol);
    }
    
    // Destroy the iterator
    higcit_destroy(it);
    
    // Sync the distributed properties in the new solver
    dp_sync(ns2->ed.mult.dpvisc);
    dp_sync(ns2->ed.mult.dpfracvol);
}

void higflow_interpolate_density(higflow_solver *ns, higflow_solver *ns2) {
    // Obter o subdomínio local para as células da malha ANTIGA (Fonte)
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    
    // Obter o subdomínio local para as células da NOVA malha (Destino)
    sim_domain *sdm2 = psd_get_local_domain(ns2->ed.mult.psdmult);
    
    // Obter o mapa de propriedades para a nova malha
    mp_mapper *mp2 = sd_get_domain_mapper(sdm2);
    
    // Iterar sobre cada célula da NOVA malha
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm2); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Obter a célula atual
        hig_cell *c = higcit_getcell(it);
        
        // Obter o identificador da célula no novo mapa
        int clid = mp_lookup(mp2, hig_get_cid(c));
        
        // Obter o centro da célula
        Point ccenter;
        hig_get_center(c, ccenter);
        
        // INTERPOLAÇÃO: Buscar a fração de volume na malha ANTIGA (sdm) usando
        // o stencil antigo (ns->stn)
        real fracvol = compute_value_at_point(sdm2, ccenter,
                                              ccenter, 1.0,
                                              ns2->ed.mult.dpfracvol,
                                              ns2->ed.mult.stn);
        
        // Recalcular a densidade com base na fração interpolada e nas propriedades dos fluidos
        // Usa o tempo atual (ns->par.t) para propriedades que variam no tempo
        real dens0 = ns->ed.mult.get_density0(ccenter, ns->par.t);
        real dens1 = ns->ed.mult.get_density1(ccenter, ns->par.t);
        
        // Mistura linear baseada na fração de volume (regra da mistura)
        real dens = (1.0 - fracvol) * dens0 + fracvol * dens1;
        
        // Armazenar a densidade calculada na propriedade distribuída do NOVO
        // solver (ns2)
        dp_set_value(ns2->ed.mult.dpdens, clid, dens);
    }
    
    // Destruir o iterador
    higcit_destroy(it);
    
    // Sincronizar a propriedade de densidade distribuída na nova malha
    dp_sync(ns2->ed.mult.dpdens);
}

void higflow_interpolate_bc_for_pressure(higflow_solver *ns, higflow_solver *ns2) {
  // Facet iterator
  higcit_celliterator *it;
  // Get the local sub-domain
  sim_domain *sdp = psd_get_local_domain(ns->psdp);
  sim_domain *sdp2 = psd_get_local_domain(ns2->psdp);

  int num_bc_types = 2;
  bc_type bc_t;
  for(int i = 0; i < num_bc_types; i++) {
    if(i == 0) bc_t = DIRICHLET;
    else if(i == 1) bc_t = NEUMANN;

    // Get the number of boundaries of type
    int numbcs = sd_get_num_bcs(sdp2, bc_t);
    // For each boundary
    for (int h = 0; h < numbcs; h++) {
      // Get the boundary
      // sim_boundary *bc = sd_get_bc(sdp, bc_t, h); // Cuidado: pode não corresponder ao h do sdp2
      sim_boundary *bc2 = sd_get_bc(sdp2, bc_t, h);
      
      // ===> CORREÇÃO 1: Filtrar tipos de valor <===
      // Se for Neumann, não podemos interpolar P escalar.
      // Se for Valor Fixo (fixedValue), mantemos o valor analítico (ex: 0.0) já setado.
      bc_valuetype valuetype = sb_get_valuetype(bc2);
      
      if (bc_t == NEUMANN || valuetype == fixedValue) {
          continue; // Pula para a próxima fronteira
      }

      // Get the mapper
      // mp_mapper *bm = sb_get_mapper(bc); // Não use o mapper antigo se a topologia mudou
      mp_mapper *bm2 = sb_get_mapper(bc2);

      // For each cell of the boundary
      for(it = sb_get_celliterator(bc2); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Get the cell
        hig_cell *bcell = higcit_getcell(it);
        // Get the cell center
        Point bccenter;
        hig_get_center(bcell, bccenter);
        // Get the id of the cell
        int bclid = mp_lookup(bm2, hig_get_cid(bcell));
        
        // Get the pressure in the distributed pressure property
        real interpolated_p;
        stn_reset(ns->stn);
        
        // Usa o domínio antigo para pegar o valor
        sd_get_stencil(sdp, bccenter, bccenter, 1, ns->stn);
        interpolated_p = compute_value_at_point(sdp, bccenter,
                                                bccenter, 1.0,
                                                ns->dpp, ns->stn);
        
        // Set the value
        sb_set_value(bc2, bclid, interpolated_p);
      }
      // Destroy the iterator
      higcit_destroy(it);
    }
  }
}

void higflow_interpolate_fracvolaux(higflow_solver *ns, higflow_solver *ns2) {
    // Obter o subdomínio local para as células da malha ANTIGA (Fonte)
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    
    // Obter o subdomínio local para as células da NOVA malha (Destino)
    sim_domain *sdm2 = psd_get_local_domain(ns2->ed.mult.psdmult);
    
    // Obter o mapa de propriedades para a nova malha
    mp_mapper *mp2 = sd_get_domain_mapper(sdm2);
    
    // Iterar sobre cada célula da NOVA malha
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm2); !higcit_isfinished(it); higcit_nextcell(it)) {
        // Obter a célula atual
        hig_cell *c = higcit_getcell(it);
        
        // Obter o identificador da célula no novo mapa
        int clid = mp_lookup(mp2, hig_get_cid(c));
        
        // Obter o centro da célula
        Point ccenter;
        hig_get_center(c, ccenter);
        
        // INTERPOLAÇÃO: Buscar a fração de volume na malha ANTIGA (sdm) usando o stencil antigo (ns->stn)
        real fracvolaux = compute_value_at_point(sdm, ccenter, ccenter, 1.0, 
                                              ns->ed.mult.dpfracvolaux, ns->ed.mult.stn);
        
        // Se desejar inicializar também a variável auxiliar com o mesmo valor interpolado:
        dp_set_value(ns2->ed.mult.dpfracvolaux, clid, fracvolaux);
    }
    
    // Destruir o iterador
    higcit_destroy(it);
    
    // Sincronizar as propriedades distribuídas na nova malha
    dp_sync(ns2->ed.mult.dpfracvolaux);
}

void higflow_interpolate_bc_for_velocity(higflow_solver *ns, higflow_solver *ns2) {
  // Facet iterator
  higcit_celliterator *it;
  // Local sub-domain
  sim_facet_domain *sfdu[DIM];
  sim_facet_domain *sfdu2[DIM];
  // For each dimension
  for(int dim = 0; dim < DIM; dim++) {
    // Get the local sub-domain
    sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
    sfdu2[dim] = psfd_get_local_domain(ns2->psfdu[dim]);
    sim_domain *sd = sfdu[dim]->cdom;
    sim_domain *sd2 = sfdu2[dim]->cdom;

    int num_bc_types = 2;
    bc_type bc_t;
    for(int i = 0; i < num_bc_types; i++) {
      if(i == 0) bc_t = DIRICHLET;
      else if(i == 1) bc_t = NEUMANN;

      // Get the number of boundaries of type
      // int numbcs = sd_get_num_bcs(sd, bc_t);
      int numbcs2 = sd_get_num_bcs(sd2, bc_t);
      // For each boundary
      for (int h = 0; h < numbcs2; h++) {
        // Get the boundary
        // sim_boundary *bc = sd_get_bc(sd, bc_t, i);
        sim_boundary *bc2 = sd_get_bc(sd2, bc_t, h);
        // Get the id defined by the user
        int userid       = sb_get_userid(bc2);
        // Get the value type of the boundary condition
        bc_valuetype valuetype = sb_get_valuetype(bc2);
        // Get the mapper
        mp_mapper *bm2    = sb_get_mapper(bc2);
        // For each cell of the boundary
        for(it = sb_get_celliterator(bc2); !higcit_isfinished(it); higcit_nextcell(it)) {
          // Get the cell
          hig_cell *bcell = higcit_getcell(it);
          // Get the cell center
          Point bccenter;
          hig_get_center(bcell, bccenter);
          // Get the id of the cell
          int bclid = mp_lookup(bm2, hig_get_cid(bcell));
          stn_reset(ns->stn);
          sfd_get_stencil(sfdu[dim], bccenter, bccenter, 1, ns->stn);
          // real interpolated_u = sfd_dp_interpolate(sfdu[dim], ns->dpu[dim],bccenter, bccenter, ns->stn);
          real interpolated_u = dp_interpolate_from_stencil(ns->dpu[dim], ns->stn);
          real t   = ns2->par.t + ns2->par.dt;
          // Get the velocity defined by the user
          real val = ns2->func.get_boundary_velocity(userid, bccenter, dim, t);
          sb_set_value(bc2, bclid, val);
        }
        // Destroy the iterator
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

    // Set data file type names
    higflow_load_data_file_names(argc, argv, ns);

    print0f("=+=+ Loading Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);

    p_par = create_initialize_physical_parameters(ns, myrank);
    init_global_var(ns);

    print0f("=+=+ Loading, Creating, Partitioning and Initializing Domains =+=+=+=+=+=\n");
    create_initialize_all_domains(ns, myrank, ntasks);

    print0f("=+=+ Creating and Initializing Distributed Properties =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_create_distributed_properties(ns);
    if(ns->par.step == 0) higflow_initialize_distributed_properties(ns);

    // get inlet boundary types to set boundary conditions correctly
    get_inlet_types(ns);

    print0f("=+=+ Initializing Boundaries =+=+\n");
    higflow_initialize_boundaries_yaml(ns);
    
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

        if (ns->par.step % 5 == 0) {
              // Initializing Navier-Stokes solver
              // Create Navier-Stokes solver
              higflow_solver *ns2 = higflow_create();

              // Load the data files
              higflow_load_data_file_names(argc, argv, ns2); 
              print0f("=+=+=+= Load Controllers and Parameters (ns2) =+=+=+=+=+\n");
              higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
              // // set the external functions
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
              partition_graph *pg = pg_create(MPI_COMM_WORLD);
              // Initializing partition table
              // Setting the fringe size of the sub-domain
              // The fringe is a buffer around the cells of a given node
              pg_set_fringe_size(pg, 5);
              /* Partitioning the grid from AMR information */
              load_balancer *lb = lb_create(MPI_COMM_WORLD, 1);
              /* Creating the distributed HigTree data structure */

              real thresholds[] = {0.08, 0.04, 0.02}; 
              // 2 para passo de tempo == 5 e 3 otherwise 
                
              int num_levels;
              if(ns->par.step <= 5) {
                num_levels = 1;
              } else if( ns->par.step > 5 && ns->par.step <= 10 ) {
                num_levels = 2;
              } else if(ns->par.step > 10 && ns->par.step <= 15){
                num_levels = 3;
              } else {
                num_levels = 3;
              }

              hig_cell *root = higflow_make_adapted_tree_params(ns, num_levels, thresholds);
              // hig_cell *root = higflow_save_refined_mesh_preview(ns, ns->par.step);
              // hig_cell *root = sd_get_higtree(ns->sdp, 0);
              // Add higtree for SDs
              sd_add_higtree(ns2->sdp, root);
              sd_add_higtree(ns2->sdF, root);
              sd_add_higtree(ns2->ed.mult.sdmult, root);
              if (ns2->contr.flowtype == MULTIPHASE) {
                  if(ns2->ed.mult.contr.viscoelastic_either == true) {
                      sd_add_higtree(ns2->ed.sdED, root);
                  }
              }
              lb_destroy(lb);

              // // Creating the partitioned sub-domain to simulation
              higflow_create_partitioned_domain(ns2, pg, order_center);
              higflow_create_partitioned_domain_multiphase(ns2, pg, order_center);

              // Creating the stencil for properties interpolation
              higflow_create_stencil(ns2);
              higflow_create_stencil_multiphase(ns2);

              // Creating distributed property  
              print0f("=+=+=+= Creating distributed property (ns2) +=+=+=+=+=\n");
              higflow_create_distributed_properties(ns2);

              // Treatment for boundary conditions
              higflow_initialize_boundaries_yaml(ns2);

              // Interpolar
              print0f("=+=+=+= Interpolation (ns2) +=+=+=+=+=\n");
              // dpu
              higflow_interpolate_velocity(ns, ns2);
              // dpp
              higflow_interpolate_viscosity(ns, ns2);
              higflow_interpolate_density(ns, ns2);

              higflow_compute_curvature_interfacial_force_normal_multiphase_2D_hf_shirani(ns2);
              higflow_compute_distance_multiphase_2D(ns2);
              higflow_compute_plic_lines_2d(ns2);

              higflow_interpolate_pressure(ns, ns2);
              higflow_interpolate_bc_for_velocity(ns, ns2);
              higflow_interpolate_bc_for_pressure(ns, ns2);
              higflow_create_solver(ns2); 

              // Destroy the Navier-Stokes object
   
              // 2. Copia parâmetros essenciais do solver antigo
              ns2->par = ns->par;
              ns2->contr = ns->contr;
              higflow_destroy(ns);
              ns = (higflow_solver *) ns2;

              // 1. Atualizar Mappers (Garante que os IDs globais do PETSc estejam certos)
              // É boa prática sincronizar pressão também
              psd_synced_mapper(ns->psdp); 
              for(int dim = 0; dim < DIM; dim++) {
                psfd_synced_mapper(ns->psfdu[dim]); 
              }

              // ===> FIM DA INSERÇÃO <===
              // higflow_print_vtk(ns, myrank);
        }
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

