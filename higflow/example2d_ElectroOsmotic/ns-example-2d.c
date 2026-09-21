// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 03/2023
// *******************************************************************
// *******************************************************************
//
// Electroosmotic flow: an applied potential drives the fluid through the charge in
// the double layer.  The ionic distribution model is chosen by `model:` under the
// electroosmotic block (pnp here).
//
// The suite runs this directory as TWO cases -- example2d_ElectroOsmotic and
// example2d_ElectroOsmotic_eo -- which differ only in configuration.

#include "ns-example-2d.h"

/************************************ user functions **************************************/

#include "ns-user-functions-vof.c"

#include "ns-user-functions-electroosmotic.c"

#include "ns-user-functions-viscoelastic.c"

#include "ns-user-functions-newtonian-gn.c"

// Infraestrutura identica nos dois exemplos que usam este conjunto de funcoes
// de usuario; o que distingue cada um fica no main e nas funcoes de modelo.
#include "../examples-common/eo-droplet-2d.c"

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
    // Registro por objeto: a interface substitui os oito ponteiros.
    higflow_set_problem(ns, &problema);

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
                higflow_create_domain_multiphase_electroosmotic(ns, cache, order_center, &problema_eo);
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
        higflow_create_domain_electroosmotic(ns, cache, order_center, &problema_eo);
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet);
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
    
    higflow_print_vtk(ns, myrank);  // alem do XDMF: o verificador da suite le VTK
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

        ///////////////////////////////////////////////////////
        solver_step(ns);
        // adapt_mesh_and_update_solver(ns);
        higflow_save_refined_mesh_preview(ns, ns->par.step);
        //higflow_adjust_timestep(ns);
        // printf("Num cells HT[0][0] %d\n", *ns->sdp->higtrees[0][0].numcells);
        // int nc[DIM];
        // POINT_ASSIGN_INTS(nc, 100, 100);
        // hig_refine_uniform(&ns->sdp->higtrees[0][0], nc);
        // printf("Num cells HT[0][0] %d\n", *ns->sdp->higtrees[0][0].numcells);
        ///////////////////////////////////////////////////////

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

