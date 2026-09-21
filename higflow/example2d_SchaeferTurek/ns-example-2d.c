// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
//
// Newtonian channel flow, 2D.  The simplest case in the suite and the one to read
// first: it shows the minimum an example has to provide -- boundary values, initial
// state, and the main loop.
//
// It is also the case wired to the t8code mesh backends (HIGFLOW_MALHA selects the
// source); the default path reads the .amr files as before.
//
// Config: singlephase / newtonian / semi_implicit_euler, in
// input/example-2d.load.par.contr.yaml.

#include "ns-example-2d.h"

#ifdef HIGFLOW_COM_T8CODE
// Definida em malha-t8.cxx, compilada so' quando T8CODE esta' no ambiente.
// Definida em ../examples-common/malha-t8.cxx, compilada so' com T8CODE.
extern "C" void malha_t8_instala(higflow_solver *ns, int myrank);
#endif

// A fronteira imersa NAO depende do t8code -- fica fora do ifdef acima.
#include "../src/hig-flow-fronteira-imersa.h"
extern "C" void fronteira_imersa_instala(higflow_solver *ns, fi_corpo *corpo);

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

real dpdx = 3.0;
real L = 8.0;
//higflow_solver *nsaux;

// ---------------------------------------------------------------------------
// O problema deste exemplo, agora como um tipo em vez de oito funcoes soltas.
// As oito eram registradas de uma vez por higflow_set_external_functions; o
// registro passa a ser higflow_set_problem com a instancia abaixo.  Os corpos
// sao os mesmos, so' mudaram de lugar e perderam o prefixo get_.
// ---------------------------------------------------------------------------
class NewtProblem : public HigFlowProblem {
public:
    // Value of the pressure
    real pressure(Point center, real t) {
        real value = 0.0;
        return value; 
    }
    // Value of the velocity
    real velocity(Point center, int dim, real t) {
        real value = 0.0;
        return value; 
    }
    // Value of the cell source term
    real source_term(Point center, real t) {
        real value = 0.0;
        return value; 
    }
    // Value of the facet source term
    real facet_source_term(Point center, int dim, real t) {
        real value = 0.0;
        return value; 
    }
    // Value of the pressure at boundary
    real boundary_pressure(int id, Point center, real t) {
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
    real boundary_velocity(int id, Point center, int dim, real t) {
        real value;
        switch (id) {
            case 0:
                switch (dim) {
                    case 0: ;
                        //set max velocity = 8.0e-4 
                        //value = 3.2e-3*(-center[1]*center[1] + center[1]);
                        //value = 4.0*(-center[1]*center[1] + 0.25);
                        //set max velocity = 1.5 
                        //value = -4.0*center[1]*(center[1] - 1.0);
                        // SCHAEFER-TUREK 2D-1, adimensionalizado por D e pela
                        // velocidade MEDIA.  No original: canal H = 0,41 e
                        // u(y) = 4 Um y (H-y) / H^2 com Um = 0,3, media 2Um/3 =
                        // 0,2, e Re = media * D / nu = 20.
                        //
                        // Dividindo todo comprimento por D = 0,1 e a velocidade
                        // pela media, fica H = 4,1, media = 1, Um = 1,5 -- e o
                        // Re do arquivo passa a ser o Re do benchmark, que e'
                        // exatamente o que o HiGFlow espera (ele aplica 1/Re no
                        // termo viscoso).  Cd e Cl sao adimensionais e nao mudam
                        // com essa escala.
                        value = 6.0*center[1]*(4.1 - center[1])/(4.1*4.1);
                        //value = 1.25*(1.0 - center[1]*center[1]*center[1]*center[1]);
                        //value = 1.0*(1.0 - fabs(center[1]));
                        //value = 2.0*(1.0 - sqrt(fabs(center[1])));
                        //value = 1.0;   
                        //value = 0.0;
                                            // if(t == 0.0) value = 0.0;
                        // else{ // set stagnation pressure
                        //     hig_cell *c = sd_get_cell_with_point(nsaux->sdp, center);
                        //     Point ccenter;
                        //     hig_get_center(c, ccenter);
                        //     sim_stencil *stn = stn_create();
                        //     real p = compute_value_at_point(nsaux->sdp, ccenter, center, 1.0, nsaux->dpp, stn);
                        //     real p_0 = dpdx * L;
                        //     value = sqrt(2.0*fabs(p_0-p));
                        //     printf("p = %f, p_0 = %f, value = %f\n", p, p_0, value);
                        //     stn_destroy(stn);
                        // }
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
    real boundary_source_term(int id, Point center, real t) {
        real value = 0.0;
        return value; 
    }
    // Value of the facet source term at boundary
    real boundary_facet_source_term(int id, Point center, int dim, real t) {
        real value = 0.0;
        return value; 
    }
};

static NewtProblem problema;

// Value of the viscosity
real get_viscosity(Point center, real q, real t) {
    real value = 1.0;
    return value; 
}

// Value of the boundary viscosity
real get_boundary_viscosity(int id, Point center, real q, real t) {
    real value = 1.0;
    return value; 
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
	print0f("=+=+=+= Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
        // set the external functions
    // Registro por objeto: a interface substitui os oito ponteiros.
    higflow_set_problem(ns, &problema); 
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
    // O EXEMPLO ESCOLHE A FONTE DA MALHA: uma linha, e so' com T8CODE no build.
    // Sem HIGFLOW_MALHA no ambiente nada muda, e a suite padrao afirma isso.
#ifdef HIGFLOW_COM_T8CODE
    malha_t8_instala(ns, myrank);
#endif
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 

    // O OBSTACULO.  Circulo de raio 0,25 em (2,0), no canal [0,8]x[-1,1] com
    // h = 0,05: diametro de 10 celulas, obstrucao de 25% -- a mesma ordem do
    // benchmark de Schaefer-Turek, e nao 50%, que e' o que raio 0,5 daria.
    //
    // O poligono tem 64 lados; com perimetro 2*pi*0,25 = 1,571 isso da' lado de
    // 0,0245, pouco abaixo de h, que e' o espacamento que o nucleo regularizado
    // quer.  A curva entra por SEGMENTOS -- `fi_cria_circulo` so' constroi o
    // poligono e chama `fi_cria_curva`.
    fi_corpo *obstaculo = NULL;
    {
        // O cilindro do benchmark: centro (0,2 ; 0,2) e D = 0,1 no original,
        // que na escala de D vira centro (2,2) e raio 0,5.
        //
        // O centro em y = 2 num canal [0 ; 4,1] esta' LIGEIRAMENTE FORA DO EIXO
        // -- 2 do fundo, 2,1 do topo.  Isso e' de proposito no benchmark: a
        // assimetria e' o que fixa o Cl publicado.  Centralizar "para ficar
        // bonito" invalida a comparacao.
        Point centro; centro[0] = 2.0; centro[1] = 2.0;
        for (int d = 2; d < DIM; d++) centro[d] = 0.0;
        // 128 lados: perimetro pi, lado 0,0245, abaixo de h = 0,05.
        obstaculo = fi_cria_circulo(ns->sfdu[0], centro, 0.5, 128, 0.05);
        fronteira_imersa_instala(ns, obstaculo);
        print0f("=+=+=+= Schaefer-Turek 2D-1: cilindro D=1 em (2,2), 20 celulas "
                "no diametro =+=+=+=\n");
    }

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
            print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            ns->par.ts += ns->par.dts;
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    // Stop the total time
    STOP_CLOCK(total);
    // Getting the execution time 
    // O PRIMEIRO ORACULO.  Com forca defasada o residuo nao vai a zero, vai a
    // O(dt): o teste e' refinar dt e ver a TAXA, nao olhar o valor.
    if (obstaculo != NULL)
        print0f("=+=+=+= RESIDUO_NAO_ESCORREGAMENTO %.6e  (dt = %.6e) =+=+=+=\n",
                (double) fi_residuo_max(obstaculo), (double) ns->par.dt);

    if(myrank == 0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
