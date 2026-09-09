// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************

#include "ns-example-2d.h"
#include "utils.h"
#include <unistd.h> // Required for sleep()

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************


real dpdx = 3.0;
real L = 8.0;
//higflow_solver *nsaux;

// Value of the pressure
real get_pressure(Point center, real t) {
  real value = 0.0;
  return value; 
}

// Value of the velocity
real get_velocity(Point center, int dim, real t) {
  real value = 0.0;
  return value; 
}

// Value of the cell source term
real get_source_term(Point center, real t) {
  real value = 0.0;
  return value; 
}

// Value of the facet source term
real get_facet_source_term(Point center, int dim, real t) {
  real value = 0.0;
  return value; 
}

// Value of the viscosity
real get_viscosity(Point center, real q, real t) {
  real value = 1.0;
  return value; 
}

// Value of the pressure at boundary
real get_boundary_pressure(int id, Point center, real t) {
  real value;
  switch (id) {
    case 0:
      value = 0.0;      
      break;
    case 1:
      value = 0.0;        
      break;
    case 2:
      value = 0.0;        
      break;
    case 3:
      value = 0.0;        
      break;
  }
  return value; 
}

// Value of the velocity at boundary
real get_boundary_velocity(int id, Point center, int dim, real t) {
  real value;
  switch (id) {
    case 0:
      switch (dim) {
        case 0: ;
          value = 1.5*(1.0 - center[1]*center[1]);
          // value = 1.5;
          break;
        case 1:
          value = 0.0;
          break;
      }
      break;
    case 1:
      switch (dim) {
        case 0:
          value = 0.0;
          //value = 1.0;
          break;
        case 1:
          value = 0.0;
          break;
      }
      break;
    case 2:
      switch (dim) {
        case 0:
          value = 0.0;
          break;
        case 1:
          value = 0.0;
          break;
      }
      break;
    case 3:
      switch (dim) {
        case 0:
          value = 0.0;
          break;
        case 1:
          value = 0.0;
          break;
      }
      break;
  }
  return value; 
}

// Value of the cell source term at boundary
real get_boundary_source_term(int id, Point center, real t) {
  real value = 0.0;
  return value; 
}

// Value of the facet source term at boundary
real get_boundary_facet_source_term(int id, Point center, int dim, real t) {
  real value = 0.0;
  return value; 
}

// Value of the boundary viscosity
real get_boundary_viscosity(int id, Point center, real q, real t) {
  real value = 1.0;
  return value; 
}

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
    interpolated_p = compute_value_at_point(sdp, ccenter, ccenter, 1.0, ns->dpp, ns->stn);

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

// *******************************************************************
// Navier-Stokes main program
// *******************************************************************

// Main program for the Navier-Stokes simulation 
int main (int argc, char *argv[]) {
  // Initialize the total time counting
  START_CLOCK(total);
  // Number of tasks
  int ntasks;
  // Identifier of the process
  int myrank;
  // Initializing Navier-Stokes solver
  higflow_initialize(&argc, &argv, &myrank, &ntasks);
  // Create Navier-Stokes solver
  higflow_solver *ns = higflow_create();
  // Load the data files
  higflow_load_data_file_names(argc, argv, ns); 
  print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+n");
  higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
  // set the external functions
  higflow_set_external_functions(ns, get_pressure, get_velocity, 
                                 get_source_term, get_facet_source_term,
                                 get_boundary_pressure, get_boundary_velocity,
                                 get_boundary_source_term, get_boundary_facet_source_term); 
  // Set the order of the interpolation to be used in the SD. 
  int order_center = 2;
  int order_facet = 2;
  // Set the cache: Reuse interpolation, 0 on, 1 off
  int cache = 1;

  // Create the simulation domain
  higflow_create_domain(ns, cache, order_center); 

  // Initialize the domain
  print0f("=+=+=+= Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
  //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
  higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 

  // Initialize the boundaries
  print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
  //higflow_initialize_boundaries(ns);
  higflow_initialize_boundaries_yaml(ns);

  // Creating distributed property  
  higflow_create_distributed_properties(ns);
  // Initialize distributed properties
  if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);
  // Create the linear system solvers
  higflow_create_solver(ns);

  // Load the properties form 
  if (ns->par.step > 0) {
    // Loading the velocities 
    if (myrank == 0) {
      printf("*********************************************************************************\n");
      printf("*********************************************************************************\n");
      printf("===> Reloading properties from previous simulation <====> step = %d <====> t = %15.10lf <===\n", ns->par.step, ns->par.t);
      printf("*********************************************************************************\n");
      printf("*********************************************************************************\n");
    }
    higflow_load_properties(ns, myrank, ntasks);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  print0f("=+=+ Saving Domain and Boundary Properties =+=+\n");
  higflow_save_domain_yaml(ns, myrank, ntasks);
  higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
  higflow_save_all_controllers_and_parameters_yaml(ns, myrank); //copying necessary yamls

  // Printing the properties to visualize: first step
  if (ns->par.step == 0) {
    print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
    higflow_print_vtk(ns, myrank);
    //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
    ns->par.tp += ns->par.dtp;
    ns->par.frame++;
    print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
    higflow_save_properties(ns, myrank, ntasks);
    ns->par.ts += ns->par.dts;
  }

  // ********************************************************
  // Begin Loop for the Navier-Stokes equations integration
  // ********************************************************

  for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
    // Print the step
    print0f("===> Step:        %7d <====> t  = %15.10lf <===\n", ns->par.step, ns->par.t);
    // Start the first step time
    if (ns->par.step == step0)  START_CLOCK(firstiter); 
    // Update velocities and pressure using the projection method 
    higflow_solver_step(ns);
    // Time update 
    ns->par.t += ns->par.dt;
    // Stop the first step time
    if (ns->par.step == step0) STOP_CLOCK(firstiter); 
    // Printing
    if (ns->par.t >= ns->par.tp) {
      print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
      higflow_print_vtk(ns, myrank);
      //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
      ns->par.tp += ns->par.dtp;
      ns->par.frame++;
    }
    // Saving the properties
    if (ns->par.t >= ns->par.ts) {
      print0f("===> Saving <===> ts = %15.10lf <===\n", ns->par.ts);
      higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
      higflow_save_properties(ns, myrank, ntasks);
      ns->par.ts += ns->par.dts;
    }

    // Condition for testing, do this step only for the half time
    if (ns->par.step == 50) {
      // Initializing Navier-Stokes solver
      // Create Navier-Stokes solver
      higflow_solver *ns2 = higflow_create();
      // Load the data files
      // argv_new 
      char *argv_new[argc];
      argv_new[1] = "input/fine_mesh";
      argv_new[2] = "output/fine_mesh.save";
      argv_new[3] = "VTKS/fine_mesh.print";
      higflow_load_data_file_names(argc, argv_new, ns2); 
      print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+\n");
      higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
      // set the external functions
      higflow_set_external_functions(ns2, get_pressure, get_velocity, 
                                     get_source_term, get_facet_source_term,
                                     get_boundary_pressure, get_boundary_velocity,
                                     get_boundary_source_term, get_boundary_facet_source_term); 

      // Create the simulation domain
      higflow_create_domain(ns2, cache, order_center); 

      // Initialize the domain
      print0f("=+=+=+= Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
      //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
      higflow_initialize_domain_yaml(ns2, ntasks, myrank, order_facet); 

      // Initialize the boundaries
      print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+\n");
      //higflow_initialize_boundaries(ns);
      higflow_initialize_boundaries_yaml(ns2);

      // Creating distributed property  
      print0f("=+=+=+= Creating distributed property  +=+=+=+=+=\n");
      higflow_create_distributed_properties(ns2);

      // Interpolar
      print0f("=+=+=+= Interpolation  +=+=+=+=+=\n");
      higflow_interpolate_velocity(ns, ns2);
      higflow_interpolate_pressure(ns, ns2);
      higflow_interpolate_bc_for_velocity(ns, ns2);
      higflow_interpolate_bc_for_pressure(ns, ns2);

      // Keep parameters like time  and steps
      ns2->par = ns->par;
      // Destroy the Navier-Stokes object

      higflow_destroy(ns);
      ns = (higflow_solver *) ns2;

      // 1. Atualizar Mappers (Garante que os IDs globais do PETSc estejam certos)
      // É boa prática sincronizar pressão também
      psd_synced_mapper(ns->psdp); 
      for(int dim = 0; dim < DIM; dim++) {
        psfd_synced_mapper(ns->psfdu[dim]); 
      }

      // 2. Criar Stencils (ns2 é novo, não tem stencils calculados ainda)
      higflow_create_stencil(ns);

      // 3. CRIAR OS SOLVERS 
      // Atenção: Use higflow_create_solver em vez de realloc, 
      // pois ns2 é um objeto novo que nunca teve solver.
      higflow_create_solver(ns); 

      // ===> FIM DA INSERÇÃO <===
      higflow_print_vtk(ns, myrank);
    }
  }
  // ********************************************************
  // End Loop for the Navier-Stokes equations integration
  // ********************************************************
  // adapt_mesh()

  // Stop the total time
  STOP_CLOCK(total);
  // Getting the execution time 
  if(myrank == 0) {
    DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
    DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
  }
}
