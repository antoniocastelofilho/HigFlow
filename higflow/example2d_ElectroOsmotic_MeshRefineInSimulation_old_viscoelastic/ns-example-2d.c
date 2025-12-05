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

/******************************************************************************************/
/******************************************************************************************/
/********************************* main user functions ************************************/
/******************************************************************************************/
/******************************************************************************************/


void create_initialize_all_domains(higflow_solver* ns, int myrank, int ntasks) {
    // Set the order of the interpolation to be used in the SD. 
    int order_center = 2;
    int order_facet = 2;
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
    
    //higflow_print_vtk(ns, myrank);
    // if(ns->contr.flowtype == MULTIPHASE) higflow_print_vtk2D_multiphase(ns, myrank);
    //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
    write_xdmf(ns);
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


// Defina um nível máximo para evitar refinamento infinito
#define MAX_REF_LEVEL 3

// Função auxiliar para pegar o nível da célula
int get_cell_level(hig_cell *c) {
    int level = 0;
    hig_cell *p = hig_get_parent(c);
    while (p != NULL) {
        level++;
        p = hig_get_parent(p);
    }
    return level;
}

// Função para criar um snapshot da malha refinada sem afetar a simulação
// Certifique-se de que estes includes estão no topo do arquivo ns-example-2d.c
#include "higtree.h"
#include "higtree-io.h" 
#include "higtree-iterator.h"
#include "pdomain.h"

void higflow_save_refined_mesh_preview(higflow_solver *ns, int frame_id) {
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 1. Obtém a raiz da árvore local
    hig_cell *root_original = sd_get_higtree(sdm, 0);
    if (root_original == NULL) return;

    // 2. Clone a estrutura (Topologia apenas)
    hig_cell *root_copy = hig_clone(root_original);
    if (root_copy == NULL) return;

    if(myrank == 0) printf("Preview: Gerando malha virtual refinada (Frame %d)...\n", frame_id);

    // 3. Preparação para refino
    int capacity = 1000;
    int count = 0;
    hig_cell **cells_to_refine = (hig_cell **)malloc(capacity * sizeof(hig_cell*));

    higcit_celliterator *it = higcit_create_all_higtree(root_copy);
    while (!higcit_isfinished(it)) {
        hig_cell *c_copy = higcit_getcell(it);
        
        Point ccenter;
        hig_get_center(c_copy, ccenter);
        
        // CORREÇÃO AQUI: Usando _with_point em vez de _at_point
        hig_cell *c_orig = sd_get_cell_with_point(sdm, ccenter);
        
        if (c_orig != NULL) {
            int clid = mp_lookup(mp, hig_get_cid(c_orig));
            if (clid >= 0) {
                real fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
                
                if (fracvol > 0.001 && fracvol < 0.999) {
                    if (get_cell_level(c_copy) < MAX_REF_LEVEL) {
                        if (count >= capacity) {
                            capacity *= 2;
                            cells_to_refine = (hig_cell **)realloc(cells_to_refine, capacity * sizeof(hig_cell*));
                        }
                        cells_to_refine[count++] = c_copy;
                    }
                }
            }
        }
        higcit_nextcell(it);
    }
    higcit_destroy(it);

    // 4. Aplica Refino na Cópia
    for (int i = 0; i < count; i++) {
        int numcells[DIM];
        hig_get_cells_per_dim(cells_to_refine[i], numcells);
        POINT_ASSIGN_INTS(numcells, 2, 2); 
        hig_refine_uniform(cells_to_refine[i], numcells);
    }
    free(cells_to_refine);

    // 5. Salva VTK da Cópia
    char filename[256];
    sprintf(filename, "preview_mesh_r%d_f%d.vtk", myrank, frame_id);
    
    FILE *fd = fopen(filename, "w");
    if (fd) {
        // Tenta hig_print_vtk. Se der erro de linker, tente higio_print_vtk
        higio_print_in_vtk2d(fd, root_copy);
        fclose(fd);
    }

    // 6. Limpeza
    hig_destroy(root_copy);
}

void adapt_mesh_and_update_solver(higflow_solver *ns) {
    if (ns->contr.flowtype != MULTIPHASE) return;

    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    
    higcit_celliterator *it;
    int changes = 0;

    // 1. Refinamento
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        
        // Pega o centro e a fração de volume
        Point ccenter;
        hig_get_center(c, ccenter);
        real fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
        
        // Critério de refinamento: Interface presente (0 < fracvol < 1)
        // Usamos uma tolerância para evitar erros numéricos
        if (fracvol > 0.001 && fracvol < 0.999) {
            int level = get_cell_level(c);
            if (level < MAX_REF_LEVEL) {
                // Refina a célula uniformemente
                // hig_refine_uniform cria filhos com o mesmo padrão de divisão do pai
                int numcells[DIM];
                hig_get_cells_per_dim(c, numcells);
                POINT_ASSIGN_INTS(numcells, 4, 4);
                // printf("c ---- %d\n", *c->numcells);
                hig_refine_uniform(c, numcells);
                // printf("c ---- %d\n", NC);
                changes++;
            }
        }
    }
    higcit_destroy(it);   

    // higflow_print_vtk(ns, 0);
    
    // int nc[DIM];
    // POINT_ASSIGN_INTS(nc, 100, 100);
    // hig_refine_uniform(&ns->sdp->higtrees[0][0], nc);

    // Se não houve mudanças, retorna e economiza tempo
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) printf("===> Adapt Mesh: %d cells refined. Rebuilding solver structures...\n", changes);

    // 2. ATUALIZAÇÃO CRÍTICA DAS ESTRUTURAS (O que faltava no seu código)
    
    // A. Sincroniza os Mapeadores (IDs Globais vs Locais)
    psd_synced_mapper(ns->ed.mult.psdmult);
    psd_synced_mapper(ns->psdp); // Pressão
    for(int dim = 0; dim < DIM; dim++) {
         psfd_synced_mapper(ns->psfdu[dim]); // Velocidade
    }

    // B. Atualiza Stencils (Interpolações vizinhas)
    higflow_create_stencil(ns);

    // C. Realoca os Solvers Lineares e Propriedades Distribuídas
    // ATENÇÃO: Isso pode zerar a solução se não houver interpolação de dados (ver nota abaixo)
    higflow_realloc_solver(ns); 
    
    // Se sua versão do HigFlow tiver suporte a interpolação automática na realocação, ótimo.
    // Caso contrário, você precisaria interpolar manualmente os campos U e P da malha velha para a nova antes deste passo.
}

void copy_domain(higflow_solver *ns, higflow_solver *ns2, int argc, char*
                 argv[], int myrank, int cache, int order_center, int
                 order_facet, int ntasks)
{
  // Create the new Navier-Stokes solver
  ns2 = higflow_create();

  // Load the data files
  // argv_new 
  char *argv_new[argc];
  argv_new[1] = "input_fine_mesh/load";
  argv_new[2] = argv[2];
  argv_new[3] = argv[3];
  higflow_load_data_file_names(argc, argv_new, ns2); 
  print0f("=+=+=+= Load Controllers and Parameters (ns2) =+=+=+=+=+\n");
  higflow_load_all_controllers_and_parameters_yaml(ns2, myrank);
  // set the external functions
  higflow_set_external_functions(ns2, get_pressure, get_velocity, 
                                 get_source_term, get_facet_source_term,
                                 get_boundary_pressure, get_boundary_velocity,
                                 get_boundary_source_term, get_boundary_facet_source_term); 

  // Create the simulation domain
  higflow_create_domain(ns2, cache, order_center); 

  // Initialize the domain
  print0f("=+=+=+= Load Domain (ns2) =+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
  //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
  higflow_initialize_domain_yaml(ns2, ntasks, myrank, order_facet); 

  // Initialize the boundaries
  print0f("=+=+=+= Load Bondary Condtions (ns2) =+=+=+=+=+=+=+=+=+\n");
  higflow_initialize_boundaries_yaml(ns2);

  // Creating distributed property  
  print0f("=+=+=+= Creating distributed property (ns2)  +=+=+=+=+=\n");
  higflow_create_distributed_properties(ns2);

  // Interpolar
  print0f("=+=+=+= Interpolation (ns2)  +=+=+=+=+=\n");
  higflow_interpolate_velocity(ns, ns2);
  higflow_interpolate_pressure(ns, ns2);
  higflow_interpolate_viscosity(ns, ns2);
  higflow_interpolate_density(ns, ns2);
  higflow_interpolate_bc_for_velocity(ns, ns2);
  higflow_interpolate_bc_for_pressure(ns, ns2);

  // Keep parameters like time and steps
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

  // Para debug
  // higflow_print_vtk(ns, myrank);
}

int main(int argc, char* argv[]) {
    int errcode = 0;
    START_CLOCK(total);
    int ntasks; // Number of tasks
    int myrank; // Identifier of the process
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

        if (myrank == 0) {
            printf("*********************************************************************************\n");
            printf(" ===> Step:   %7d <====> t   = %15.10lf <===\n", ns->par.step, ns->par.t);
            printf("*********************************************************************************\n");
        }
        GET_NSEC_CLOCK(currentiter) = 0.0; START_CLOCK(currentiter);

        solver_step(ns);

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

        if (ns->par.step == 5) {
          higflow_solver *ns2;
          int order_center = 2;
          int order_facet = 2;
          int cache = 1;
          copy_domain(ns, ns2, argc, argv, myrank, cache, order_center, order_facet, ntasks);
         }

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

