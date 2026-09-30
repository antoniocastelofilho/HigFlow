// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************
//
// CASO MINIMO VISCOELASTICO 3D -- a pergunta e' se a maquinaria RODA em tres
// dimensoes, nao se algum escoamento esta' certo.
//
// MONTAGEM.  Caixa fechada em repouso, esfera de fase 1 no centro, fase 0
// (ambiente) VISCOELASTICA com De e beta nao nulos.  Sem gravidade, sem entrada,
// sem tensao superficial: o estado inicial e' um EQUILIBRIO EXATO.
//
// O ORACULO E' ESSE EQUILIBRIO.  Partindo do repouso com o tensor de conformacao
// na identidade -- que e' o equilibrio do Oldroyd-B --, tudo tem de FICAR ali.
// Qualquer velocidade ou desvio da identidade e' espurio, e mede o erro da
// maquinaria, nao a fisica.
//
// E' teste fraco de fisica e forte de encanamento: pega NaN, laco [DIM][DIM] que
// supoe duas dimensoes, componente z nao inicializada, e indice trocado -- que e'
// exatamente o que pode ter sobrado num caminho nunca compilado em 3D.
//
// POR QUE ESTE CASO EXISTE.  Nem esta arvore nem a versao de junho jamais
// rodaram 3D exercitando viscoelasticidade: a de junho tinha a viscosidade
// polimerica FIXADA EM ZERO no codigo.  Aqui De e beta vem do YAML e sao nao
// nulos.

#include "ns-example-3d.h"
#include <stdlib.h>

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// ---------------------------------------------------------------------------
// O problema deste exemplo, como um tipo em vez de oito funcoes soltas.
// Os corpos sao os mesmos; so' mudaram de lugar e perderam o prefixo get_.
// ---------------------------------------------------------------------------
class VE3DProblem : public HigFlowProblem, public HigFlowMultiphaseProblem,
                    public HigFlowMultiphaseViscoelasticProblem {
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
    //
    // TAMPA MOVEL, por ambiente (FT3_TAMPA, padrao 0).  Com zero a caixa fica em
    // repouso e a gota estatica nao muda em bit nenhum -- o caso do B2 continua
    // sendo o padrao.  Com valor nao nulo, a face y=L (id 3) desliza em x e gera
    // um vortice que CARREGA a gota: e' o que poe a adveccao acoplada a prova,
    // porque com a gota parada a velocidade e' ~1e-3 e ela mal se move.
    real boundary_velocity(int id, Point center, int dim, real t) {
        if (id == 3 && dim == 0) {
            const char *e = getenv("FT3_TAMPA");
            return (e != NULL) ? atof(e) : 0.0;
        }
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
    // --- interface multifasica -------------------------------------------
    // Propriedades por ambiente.  Os padroes sao IGUAIS nas duas fases de
    // proposito: e' o estagio que move o encanamento para o caminho multifasico
    // SEM mudar a fisica, de modo que qualquer diferenca no Laplace venha do
    // encanamento e nao do salto de propriedade.  Foi assim que o 2D fez.
    //
    // CONVENCAO DO SOLVER, conferida no codigo e NAO suposta:
    //   dens = (1 - fracvol)*dens0 + fracvol*dens1
    // logo fracvol=1 e' a FASE 1, e fase 1 e' DENTRO da gota.
    static real _amb(const char *nome, real padrao) {
        const char *s = getenv(nome);
        return (s != NULL) ? atof(s) : padrao;
    }
    real viscosity0(Point center, real t) { return _amb("VE_MU0", 1.0); }
    real viscosity1(Point center, real t) { return _amb("VE_MU1", 1.0); }
    real density0(Point center, real t)   { return _amb("VE_RHO0", 1.0); }
    real density1(Point center, real t)   { return _amb("VE_RHO1", 1.0); }

    // Fracao INICIAL.  Do segundo passo em diante quem manda e' a superficie,
    // pelo gancho; esta funcao so' serve a' condicao inicial, e usa a MESMA
    // esfera para as duas nascerem coerentes.
    real fracvol(Point center, Point delta, real t) {
        real cx = _amb("VE_CX", 0.5), cy = _amb("VE_CY", 0.5);
        real cz = _amb("VE_CZ", 0.5), R  = _amb("VE_R", 0.25);
        real d[3] = { center[0]-cx, center[1]-cy, center[2]-cz };
        real r = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        // degrau suavizado na espessura de uma celula: a condicao inicial nao
        // precisa ser exata (o gancho recalcula), mas um degrau abrupto deixaria
        // o primeiro passo com um gradiente artificial
        real e = 0.5 * delta[0];
        if (r < R - e) return 1.0;
        if (r > R + e) return 0.0;
        return 0.5 * (1.0 - (r - R) / e);
    }

    // --- multifasico viscoelastico ---
    // Conformacao INICIAL: identidade, que e' o equilibrio do Oldroyd-B em
    // repouso.  E' o estado cujo desvio o oraculo mede.
    real tensor_multiphase(real fracvol, Point center, int i, int j, real t) {
        return kernel(i, 1.0, 0.0) * (i == j);
    }
    real kernel(int dim, real lambda, real tol)          { return lambda; }
    real kernel_inverse(int dim, real lambda, real tol)  { return lambda; }
    real kernel_jacobian(int dim, real lambda, real tol) { return 1.0; }

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

static VE3DProblem problema;

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
    // O multifasico e' ADICAO ao dominio base, nao substituto: sem
    // higflow_create_domain_multiphase os dominios sdmult e sdED nascem NULOS e
    // higflow_initialize_domain_yaml segmenta em sd_add_higtree.  Foi o que
    // aconteceu aqui, e o comentario do exemplo 2D ja' avisava.
    higflow_create_domain(ns, cache, order_center);
    if (ns->contr.flowtype == MULTIPHASE)
        higflow_create_domain_multiphase(ns, cache, order_center, &problema);
    if (ns->contr.flowtype == MULTIPHASE)
        higflow_create_domain_multiphase_viscoelastic(ns, &problema);
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

    // Parametros do caso, por ambiente.
    real R_gota = 0.25, cx = 0.5, cy = 0.5, cz = 0.5;
    { const char *e;
      if ((e = getenv("VE_R"))  != NULL) R_gota = atof(e);
      if ((e = getenv("VE_CX")) != NULL) cx = atof(e);
      if ((e = getenv("VE_CY")) != NULL) cy = atof(e);
      if ((e = getenv("VE_CZ")) != NULL) cz = atof(e); }
    print0f("=+=+=+= VE3D: esfera R=%.3f em (%.2f,%.2f,%.2f); ambiente viscoelastico =+=+=+=\n",
            (double) R_gota, (double) cx, (double) cy, (double) cz);

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
        // DESPACHO POR FASE.  higflow_solver_step e' o passo MONOFASICO; com
        // flowphase: multiphase e' higflow_solver_step_multiphase que tem de
        // rodar.  Chamar o errado nao quebra nem avisa -- roda 200 passos e
        // entrega velocidade e pressao IDENTICAMENTE NULAS, que foi o que
        // aconteceu aqui.
        // VISCOELASTICO multifasico: e' este o caminho que nunca rodou em 3D.
        if (ns->contr.flowtype == MULTIPHASE)
            higflow_solver_step_multiphase_viscoelastic(ns);
        else
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

    // ---- O ORACULO: o equilibrio trivial se manteve? -------------------
    // Partindo do repouso, sem forcante, a velocidade tem de continuar nula e o
    // tensor de conformacao na identidade.  Desvio aqui e' erro da maquinaria.
    {
        real umax = 0.0;
        int  nan_u = 0;
        for (int dim = 0; dim < DIM; dim++) {
            higfit_facetiterator *fit;
            sim_facet_domain *sfdu = psfd_get_local_domain(ns->psfdu[dim]);
            for (fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit);
                 higfit_nextfacet(fit)) {
                hig_facet *f = higfit_getfacet(fit);
                int flid = mp_lookup(sfd_get_domain_mapper(sfdu), hig_get_fid(f));
                if (flid < 0) continue;
                real v = dp_get_value(ns->dpu[dim], flid);
                if (v != v) { nan_u++; continue; }
                if (fabs(v) > umax) umax = fabs(v);
            }
            higfit_destroy(fit);
        }

        // Desvio do tensor de conformacao em relacao a' identidade, e NaN nele.
        real pior_S = 0.0;
        int  nan_S = 0, ncel = 0;
        {
            sim_domain *sdm = psd_get_local_domain(ns->ed.psdED);
            mp_mapper  *mp  = sd_get_domain_mapper(sdm);
            higcit_celliterator *it;
            for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
                 higcit_nextcell(it)) {
                hig_cell *c = higcit_getcell(it);
                int clid = mp_lookup(mp, hig_get_cid(c));
                if (clid < 0) continue;
                ncel++;
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++) {
                        real a = dp_get_value(ns->ed.ve.dpKernel[i][j], clid);
                        if (a != a) { nan_S++; continue; }
                        real alvo = (i == j) ? 1.0 : 0.0;
                        if (fabs(a - alvo) > pior_S) pior_S = fabs(a - alvo);
                    }
            }
            higcit_destroy(it);
        }

        print0f("=+=+=+= VE3D EQUILIBRIO  |u|max=%.3e  pior |A-I|=%.3e  "
                "NaN(u)=%d  NaN(A)=%d  celulas=%d =+=+=+=\n",
                (double) umax, (double) pior_S, nan_u, nan_S, ncel);
        print0f("=+=+=+= VE3D %s =+=+=+=\n",
                (nan_u == 0 && nan_S == 0 && umax < 1e-6 && pior_S < 1e-6)
                  ? "PASSOU: o equilibrio se manteve"
                  : "FALHOU: o equilibrio nao se manteve");
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
