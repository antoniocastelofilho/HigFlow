// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************
//
// FASE B2 3D: GOTA ESTATICA, lei de Laplace.  Caixa fechada em repouso, uma
// esfera de front-tracking parada no centro, e a tensao superficial impondo o
// salto de pressao.  O oraculo e' Dp = 2*sigma/R -- DOIS sobre R, porque em 3D
// a curvatura media da esfera e' 2/R.  O B2 bidimensional deste repositorio
// fechou a 0,02% contra sigma/R; copiar aquele numero aqui validaria o errado.
//
// MONOFASICO de proposito: rho e mu uniformes.  Isola o termo de tensao do
// salto de propriedade, que e' o que esta fase existe para verificar.  A fracao
// de volume da superficie triangulada (necessaria para o bifasico) e' trabalho
// da fase seguinte.
//
// A superficie fica PARADA: FT3_ADVECTA=0.  Advectar sem oraculo de movimento
// seria exatamente o que os degraus anteriores existem para evitar.
//
// Mesh is 10x10x10, too coarse to divide, so the suite pins it at np=1 (max_np).

#include "ns-example-3d.h"
#include "hig-flow-front-tracking-3d.h"
#include <stdlib.h>

extern "C" void front_tracking_3d_instala(higflow_solver *ns, ft3_superficie *sup,
                                          real sigma);

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
                        value = 0.0;   // caixa em REPOUSO (nao ha' tampa movel)
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

    // ---- FASE B2 3D: a esfera e a tensao superficial -------------------
    // Parametros por ambiente, como nos exemplos 2D, para a serie de refino
    // nao exigir recompilar.
    real R_gota = 0.25, cx = 0.5, cy = 0.5, cz = 0.5, sigma = 1.0;
    int  nsub = 3;
    { const char *e;
      if ((e = getenv("FT3_R"))     != NULL) R_gota = atof(e);
      if ((e = getenv("FT3_CX"))    != NULL) cx = atof(e);
      if ((e = getenv("FT3_CY"))    != NULL) cy = atof(e);
      if ((e = getenv("FT3_CZ"))    != NULL) cz = atof(e);
      if ((e = getenv("FT3_SIGMA")) != NULL) sigma = atof(e);
      if ((e = getenv("FT3_NSUB"))  != NULL) nsub = atoi(e); }

    Point centro_gota; centro_gota[0] = cx; centro_gota[1] = cy; centro_gota[2] = cz;
    ft3_superficie *gota = ft3_cria_esfera(centro_gota, R_gota, nsub);
    if (gota == NULL) { print0f("FT3: falha ao criar a esfera\n"); exit(1); }

    // Conferencia de que a superficie chegou sa' ao solver: area, volume e
    // topologia contra os valores fechados, ANTES de qualquer passo.
    {
        real A = ft3_area(gota), V = ft3_volume(gota);
        real Aex = 4.0 * M_PI * R_gota * R_gota;
        real Vex = 4.0 / 3.0 * M_PI * R_gota * R_gota * R_gota;
        print0f("=+=+=+= FT3 esfera: nsub=%d nv=%d nt=%d  chi=%d\n"
                "        area=%.8f (exata %.8f, erro %.3e)\n"
                "        volume=%.8f (exato %.8f, erro %.3e)\n"
                "        R=%.4f sigma=%.4f  ==>  Laplace esperado Dp = 2*sigma/R = %.6f =+=+=+=\n",
                nsub, ft3_num_vertices(gota), ft3_num_triangulos(gota), ft3_euler(gota),
                (double) A, (double) Aex, (double) (fabs(A-Aex)/Aex),
                (double) V, (double) Vex, (double) (fabs(V-Vex)/Vex),
                (double) R_gota, (double) sigma, (double) (2.0*sigma/R_gota));
    }
    front_tracking_3d_instala(ns, gota, sigma);
    // --------------------------------------------------------------------
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

    // ---- O ORACULO DO B2 3D: o salto de Laplace ------------------------
    // Amostra a pressao no centro da gota (dentro) e num ponto bem fora, e
    // compara Dp com 2*sigma/R.  Tambem mede a CORRENTE ESPURIA: num equilibrio
    // exato a velocidade seria zero, e o que sobra e' o erro do acoplamento.
    {
        Point p_in, p_out;
        p_in[0] = cx; p_in[1] = cy; p_in[2] = cz;
        p_out[0] = cx; p_out[1] = cy + 2.4 * R_gota; p_out[2] = cz;
        { const char *e = getenv("FT3_POUT_Y"); if (e != NULL) p_out[1] = atof(e); }

        sim_stencil *stn = stn_create();
        real pin = 0.0, pout = 0.0;
        hig_cell *ci = sd_get_cell_with_point(ns->sdp, p_in);
        hig_cell *co = sd_get_cell_with_point(ns->sdp, p_out);
        if (ci != NULL) { Point cc; hig_get_center(ci, cc);
            pin  = compute_value_at_point(ns->sdp, cc, p_in,  1.0, ns->dpp, stn); }
        if (co != NULL) { Point cc; hig_get_center(co, cc);
            pout = compute_value_at_point(ns->sdp, cc, p_out, 1.0, ns->dpp, stn); }
        stn_destroy(stn);

        // Corrente espuria: maximo de |u| sobre as facetas proprias.
        real umax = 0.0;
        for (int dim = 0; dim < DIM; dim++) {
            higfit_facetiterator *fit;
            sim_facet_domain *sfdu = psfd_get_local_domain(ns->psfdu[dim]);
            for (fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit);
                 higfit_nextfacet(fit)) {
                hig_facet *f = higfit_getfacet(fit);
                int flid = mp_lookup(sfd_get_domain_mapper(sfdu), hig_get_fid(f));
                if (flid < 0) continue;
                real v = fabs(dp_get_value(ns->dpu[dim], flid));
                if (v > umax) umax = v;
            }
            higfit_destroy(fit);
        }
        real umax_g = umax;
        MPI_Allreduce(&umax, &umax_g, 1, MPI_HIGREAL, MPI_MAX, MPI_COMM_WORLD);

        const real dp_exato = 2.0 * sigma / R_gota;
        print0f("=+=+=+= LAPLACE 3D  p_in=%.6f  p_out=%.6f  Dp=%.6f  "
                "2*sigma/R=%.6f  erro_rel=%.4f  |u|max=%.3e =+=+=+=\n",
                (double) pin, (double) pout, (double) (pin - pout),
                (double) dp_exato,
                (double) (fabs((pin - pout) - dp_exato) / dp_exato),
                (double) umax_g);
    }
    ft3_destroi(gota);

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
