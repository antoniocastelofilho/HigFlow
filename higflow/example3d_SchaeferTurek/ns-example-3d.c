// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************
//
// Lid-driven cavity in 3D: the classic verification case, a box with one moving
// wall.
//
// Mesh is 10x10x10, too coarse to divide, so the suite pins it at np=1 (max_np).

#include "ns-example-3d.h"

#include "../src/hig-flow-fronteira-imersa.h"
extern "C" void fronteira_imersa_instala(higflow_solver *ns, fi_corpo *corpo);

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// ---------------------------------------------------------------------------
// O problema deste exemplo, como um tipo em vez de oito funcoes soltas.
// Os corpos sao os mesmos; so' mudaram de lugar e perderam o prefixo get_.
// ---------------------------------------------------------------------------
class LidDrivenProblem : public HigFlowProblem {
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
        real value = 0.0;
        return value; 
    }
    // Value of the velocity at boundary
    real boundary_velocity(int id, Point center, int dim, real t) {
        real value;
        switch (id) {
            case 0:
                // ENTRADA (x = 0), Schaefer-Turek 3D-1Z, adimensionalizado por
                // D e pela velocidade media -- mesma escala do caso 2D, e pelo
                // mesmo motivo: o HiGFlow aplica 1/Re no termo viscoso.
                //
                // No original H = 0,41 e u = 16 Um y z (H-y)(H-z) / H^4 com
                // Um = 0,45; a media e' 4Um/9 = 0,2 e Re = media*D/nu = 20.
                // Na escala de D: H = 4,1, media = 1, Um = 2,25.
                switch (dim) {
                    case 0:
                        value = 36.0 * center[1] * center[2]
                                * (4.1 - center[1]) * (4.1 - center[2])
                                / (4.1*4.1*4.1*4.1);
                        break;
                    case 1:
                    case 2:
                        value = 0.0;
                        break;
                }
                break;
            case 1:
                switch (dim) {
                    case 0:
                        value = 0.0;
                        break;
                    case 1:
                        value = 0.0;
                        break;
                    case 2:
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
                    case 2:
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
                        value = 1.0;
                        break;
                    case 2:
                        value = 0.0;
                        break;
                }
                break;
            case 4:
                switch (dim) {
                    case 0:
                        value = 0.0;
                        break;
                    case 1:
                        value = 0.0;
                        break;
                    case 2:
                        value = 0.0;
                        break;
                }
                break;
            case 5:
                switch (dim) {
                    case 0:
                        value = 0.0;
                        break;
                    case 1:
                        value = 0.0;
                        break;
                    case 2:
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

static LidDrivenProblem problema;

// Value of the Tensor
real get_tensor(Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Tensor
real get_boundary_tensor(int id, Point center, int i, int j, real t) {
    real value = 0.0;
    return value; 
}

// Value of the Kernel
real get_kernel(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = log(tol);
    //else
    //   value = log(lambda);
    //if (lambda < tol)
    //   value = sqrt(tol);
    //else
    //   value = sqrt(lambda);
    value = lambda;
    return value; 
}

// Value of the Kernel inverse
real get_kernel_inverse(int dim, real lambda, real tol) {
    real value;
    //real value = exp(lambda);
    //real value = lambda*lambda;
    value = lambda;
    return value; 
}

// Value of the Kernel Jacobian
real get_kernel_jacobian(int dim, real lambda, real tol) {
    real value;
    //if (lambda < tol)
    //   value = 1.0/tol;
    //else
    //   value = 1.0/lambda;
    //if (lambda < tol)
    //   value = 0.5/sqrt(tol);
    //else
    //   value = 0.5/sqrt(lambda);
    value = 1.0;
    return value; 
}

// Impressao de perfis: identica nos dois exemplos 3-D (ver o arquivo).
#include "../examples-common/print-3d.c"

// Print the velocity

// Print the Polymeric Tensor at point

// Print the velocity at point

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
    print0f("=+=+=+= Irrrrraaaaaaa... Load Controllers and Parameters =+=+=+=+=+=+=+=+=+=+=+=+=\n");
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
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 

    // O CILINDRO DO BENCHMARK 3D-1Z.  No original: eixo em z, centro
    // (0,2 ; 0,2), D = 0,1, atravessando todo o canal.  Na escala de D: centro
    // (5 ; 2), raio 0,5, de z = 0 a z = 4,1.
    //
    // MAS a extrusao vai de 0,10 a 4,00, UMA celula afastada de cada parede.
    // Uma e nao meia: com meia, os centros dos retalhos caem em 0,10; 0,20; ...
    // que sao FACES de celula, e a localizacao degenera ali.  Com uma celula os
    // centros ficam em 0,15; 0,25; ... -- centros de celula, como devem ser.
    // Marcador encostado na parede tem suporte caindo fora da regiao mapeada, e
    // a contribuicao se perde em silencio -- medido: 46.080 pontos de suporte
    // com peso nao nulo descartados.  As paredes ja' impoem nao escorregamento
    // por conta propria, entao afastar meia celula nao tira fisica nenhuma.
    //
    // Em 3D o corpo e' SUPERFICIE, nao curva: `fi_cria_cilindro` extruda o
    // poligono em z e o peso de cada marcador e' a AREA do retalho.  Nao precisa
    // de topologia -- corpo rigido tem peso fixo na criacao; topologia so' faz
    // falta quando for preciso curvatura, que e' o caso de interface.
    //
    // HIGFLOW_SEM_CORPO=1 roda a MESMA geometria sem o obstaculo.  Nao e'
    // andaime: e' a bifurcacao que separa "a geometria nao carrega" de "a
    // fronteira imersa quebra", e ela se paga toda vez que algo aborta.
    fi_corpo *cilindro = NULL;
    if (getenv("HIGFLOW_SEM_CORPO") != NULL) {
        print0f("=+=+=+= HIGFLOW_SEM_CORPO: canal vazio, sem o cilindro =+=+=+=\n");
    } else {
        cilindro = fi_cria_cilindro(ns->sfdu[0], 5.0, 2.0, 0.5, 64, 0.1, 4.0, 0.1);
        fronteira_imersa_instala(ns, cilindro);
        print0f("=+=+=+= Schaefer-Turek 3D-1Z: cilindro D=1 em (5;2), eixo z, "
                "10 celulas no diametro =+=+=+=\n");
    }
    print0f("=+=+=+= Irrrrrraaaaaaa... Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    // Initialize the boundaries
    print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_boundaries(ns);
    higflow_initialize_boundaries_yaml(ns);
    print0f("=+=+=+= Irrrrrraaaaaaa... Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
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
    //print0f("=+=+ Saving Domain and Boundary Properties =+=+\n");
    //higflow_save_domain_yaml(ns, myrank, ntasks);
    //higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    //higflow_save_all_controllers_and_parameters_yaml(ns, myrank); //copying necessary yamls

    // Printing the properties to visualize: first step
    if (ns->par.step == 0) {
        print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        higflow_print_vtk(ns, myrank);
        //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
        ns->par.tp += ns->par.dtp;
        ns->par.frame++;
        print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
        //higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        //higflow_save_properties(ns, myrank, ntasks);
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
            //print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
            //higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            //higflow_save_properties(ns, myrank, ntasks);
            //ns->par.ts += ns->par.dts;
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // ANTES do higflow_destroy e FORA do `if (myrank == 0)`, e nenhum dos dois
    // e' questao de estilo.
    //
    // `fi_residuo_max` faz um MPI_Allreduce, que e' COLETIVO: dentro de um ramo
    // so' do rank 0, os outros sete nunca entram na chamada e o rank 0 gira nela
    // para sempre.  Custou duas horas de uma corrida que ja' tinha terminado, e
    // o backtrace do rank ocupado foi o que resolveu -- a inferencia a partir do
    // log apontava para o lugar errado.
    //
    // E depois do destroy, `ns->par.dt` seria leitura de memoria liberada.
    if (cilindro != NULL) {
        print0f("=+=+=+= RESIDUO_NAO_ESCORREGAMENTO %.6e  (dt = %.6e) =+=+=+=\n",
                (double) fi_residuo_max(cilindro), (double) ns->par.dt);
        // DEVE SER ZERO.  Diferente de zero e' nucleo somando menos que 1 em
        // algum marcador -- corpo levemente poroso, que se le' como malha ruim.
        printf("=+=+=+= SUPORTE_PERDIDO rank %d: total %ld = espelho %ld + mapa %ld"
               "  (sem faceta: %ld) =+=+=+=\n",
               myrank, fi_suporte_perdidos(), fi_perdidos_espelho(),
               fi_perdidos_mapa(), fi_perdidos_sem_faceta());
        fflush(stdout);
        // A malha lagrangeana em VTK, ao lado da euleriana.  Sem isto o corpo
        // nao aparece em lugar nenhum.
        fi_escreve_vtk(cilindro, argv[3], 0);
    }

    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    // Stop the total time
    STOP_CLOCK(total);
    // Getting the execution time 
    if(myrank == 0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
