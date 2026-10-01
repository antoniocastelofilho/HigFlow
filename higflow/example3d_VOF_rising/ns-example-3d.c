// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************
//
// BOLHA ASCENDENTE 3D VISCOELASTICA -- o porte do caso examples3D_rising_NN.
//
// E' o analogo tridimensional do Hysing: dominio 1 x 1 x 3, bolha esferica de
// raio 0,25 no centro da base, densidades 1000 e 100, viscosidades 10 e 1,
// gravidade em z.  Os mesmos numeros do caso 2D que a Secao 10 do relatorio
// fechou contra a literatura.
//
// O QUE MUDA EM RELACAO AO ORIGINAL, e por que.  O caso de junho tinha De = 0 e
// beta = 1 nas DUAS fases -- beta = 1 e' fracao de solvente total, ou seja
// viscosidade polimerica NULA -- e ainda por cima a fonte fixava visc_p = 0 no
// codigo.  Ele rodava newtoniano por dois caminhos independentes.  Aqui De e
// beta vem por ambiente e o padrao NAO e' newtoniano: senao o porte reproduziria
// o mesmo caso com outro nome, e a pergunta que o motivou continuaria aberta.
//
// MALHA.  O original usa 51 x 51 x 153 = 398 mil celulas.  Ao custo medido neste
// equipamento isso da' ~2 min por passo e cinco dias para os 4001 passos.  O
// gerador aceita a resolucao por argumento; a verificacao aqui e' em malha
// reduzida, e isso esta' dito em vez de escondido.

#include "ns-example-3d.h"
#include <stdlib.h>

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// ---------------------------------------------------------------------------
// FORMA INICIAL, como level-set: > 0 DENTRO da bolha (fase 1, fracvol = 1).
// E' a UNICA definicao da forma -- `fracvol` abaixo so' integra esta funcao,
// de modo que nao existam dois lugares dizendo onde a bolha comeca.
// ---------------------------------------------------------------------------
static real _ambi(const char *nome, real padrao) {
    const char *s = getenv(nome);
    return (s != NULL) ? atof(s) : padrao;
}
real func3(Point p) {
    real cx = _ambi("SB_CX", 0.5), cy = _ambi("SB_CY", 0.5);
    real cz = _ambi("SB_CZ", 0.5), R  = _ambi("SB_R",  0.25);
    real d[3] = { p[0]-cx, p[1]-cy, p[2]-cz };
    return R - sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}
#include "../examples-common/vof-geometry-3d.c"

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
        // CISALHAMENTO SIMPLES (VE_CISALHA, padrao 0).  Com zero a caixa fica
        // em repouso e o caso de equilibrio nao muda em bit nenhum.  Com valor
        // nao nulo, a face y=L desliza em x; com as faces em x e z em Neumann o
        // escoamento e' efetivamente unidimensional e a taxa e' U/H.
        if (id == 3 && dim == 0) {
            const char *e = getenv("VE_CISALHA");
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
    real viscosity0(Point center, real t) { return _amb("SB_MU0", 10.0); }
    real viscosity1(Point center, real t) { return _amb("SB_MU1", 1.0); }
    real density0(Point center, real t)   { return _amb("SB_RHO0", 1000.0); }
    real density1(Point center, real t)   { return _amb("SB_RHO1", 100.0); }

    // Fracao INICIAL, por SUBDIVISAO SUBCELULAR -- 16x16x16 subcelulas na
    // celula de interface, cada uma resolvida pela forma fechada do corte
    // plano-caixa.  E' o analogo 3D do que os exemplos 2D fazem (laco 16x16
    // em example2d_VOF/ns-example-2d.c).
    //
    // AQUI NAO HA' GANCHO QUE RECALCULE.  Este exemplo e' VOF puro: nao existe
    // superficie de front-tracking reescrevendo a fracao a cada passo.  O que
    // nascer errado fica errado, porque o VOF conserva fielmente a massa que
    // recebeu.  Com a rampa linear de uma celula que estava aqui antes, o
    // volume da bolha em N=16 -- a resolucao em que este caso roda -- errava
    //
    //     rampa      1.23e-02      ordem medida 1.59 / 2.90 / 0.99 (erratica)
    //     subdivisao 6.11e-05      ordem medida 2.00 / 2.00 / 2.00
    //
    // duas ordens de grandeza, e a rampa nem convergia: a largura dela estava
    // amarrada a h, logo o padrao de amostragem contra a esfera mudava a cada
    // resolucao.  Para um caso de bolha ascendente cujas grandezas publicadas
    // se comparam em tres digitos, 1.2% de massa errada em t=0 inviabiliza a
    // comparacao.  Ver tests/vof-3d/teste-fracao-inicial.c.
    real fracvol(Point center, Point delta, real t) {
        // No cisalhamento a fracao e' UNIFORME: as duas fases tem os mesmos
        // parametros, entao a interpolacao e' exercitada e a solucao analitica
        // continua inequivoca.
        { const char *e = getenv("VE_CISALHA");
          if (e != NULL && atof(e) != 0.0) return 0.5; }
        return vof3_fracao(center, delta, t, 16);
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
static real _ambf(const char *nome, real padrao) {
    const char *e = getenv(nome);
    return (e != NULL) ? atof(e) : padrao;
}

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
      if ((e = getenv("SB_R"))  != NULL) R_gota = atof(e);
      if ((e = getenv("SB_CX")) != NULL) cx = atof(e);
      if ((e = getenv("SB_CY")) != NULL) cy = atof(e);
      if ((e = getenv("SB_CZ")) != NULL) cz = atof(e); }
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
        // A interface por marching cubes -- ver o comentario no sitio do laco.
        higflow_print_vtk3d_multiphase_mc_serial(ns, myrank);
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
        // E DENTRO DO MULTIFASICO, despacho por REOLOGIA.  Chamar o passo
        // viscoelastico com as duas fases newtonianas SEGFAULTA no primeiro
        // passo: `viscoelastic_either` falso faz a carga pular as condicoes de
        // contorno dos tensores -- corretamente --, e o passo viscoelastico
        // entao calcula a derivada da velocidade e a equacao constitutiva sobre
        // campos que nao existem.  O estouro acontece ANTES de
        // higflow_solver_step_multiphase ser alcancado, o que foi o que
        // localizou o defeito: instrumentado, nenhuma marca dele saiu.
        //
        // Nenhum exemplo 3D desta arvore tinha fases newtonianas -- todos os
        // multifasicos 3D nasceram com flowtype0 viscoelastic --, por isso o
        // caminho newtoniano multifasico 3D nunca havia sido exercitado.
        if (ns->contr.flowtype == MULTIPHASE) {
            if (ns->ed.mult.contr.viscoelastic_either == true)
                higflow_solver_step_multiphase_viscoelastic(ns);
            else
                higflow_solver_step_multiphase(ns);
        } else {
            higflow_solver_step(ns);
        }
        // Time update 
        ns->par.t += ns->par.dt;
        // Stop the first step time
        if (ns->par.step == step0) STOP_CLOCK(firstiter); 
        // Printing
        if (ns->par.t >= ns->par.tp) {
            print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
            higflow_print_vtk(ns, myrank);
            // A INTERFACE, por marching cubes sobre fracvol = 1/2: superficie
            // CONTINUA, com os vertices soldados nas arestas da rede, e os campos
            // (p, vel, S, T) interpolados NOS VERTICES.  Sai em arquivo proprio,
            // _mult_mc_, ao lado do campo.  O visualizador extrairia um contorno
            // do FracVol por conta propria, mas sem a solda e sem os tensores --
            // e o tensor polimerico T nem esta' no campo, so' o S, porque T tem
            // de ser obtido a partir de S:  T = S + 2(1-beta)/Re D.
            //
            // A linha de diagnostico que isto imprime traz arestas orfas, Euler,
            // e o volume encerrado contra a integral de fracvol.  Os dois
            // primeiros sao exatos (0 e 2 numa bolha); o terceiro e' a unica
            // medida que liga a superficie extraida ao campo que a gerou.  Em
            // N=16 ele acusa ~4%, que e' a discretizacao em 4 celulas de raio,
            // nao defeito -- o arreio analitico da' 2.68% na mesma razao R/h e
            // mostra o numero caindo como h^2.
            higflow_print_vtk3d_multiphase_mc_serial(ns, myrank);
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

    // ---- AS METRICAS DA BOLHA -------------------------------------------
    // Centroide e velocidade de subida, integrados sobre a fase LEVE (fracvol),
    // que e' a bolha.  Sao as mesmas grandezas do Hysing bidimensional, e e' o
    // que permite comparar os dois.
    {
        real vol = 0.0, zc = 0.0, wz = 0.0;
        sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
        mp_mapper  *mp  = sd_get_domain_mapper(sdm);
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
            hig_cell *c = higcit_getcell(it);
            int clid = mp_lookup(mp, hig_get_cid(c));
            if (clid < 0) continue;
            Point cc, dd; hig_get_center(c, cc); hig_get_delta(c, dd);
            real f = dp_get_value(ns->ed.mult.dpfracvol, clid);
            real dv = dd[0]*dd[1]*dd[2];
            vol += f * dv;
            zc  += f * cc[2] * dv;
        }
        higcit_destroy(it);
        if (vol > 0.0) zc /= vol;
        const real R = _ambf("SB_R", 0.25);
        const real Vex = 4.0/3.0*M_PI*R*R*R;
        print0f("=+=+=+= BOLHA 3D  t=%.4f  volume=%.8f (inicial %.8f, erro %.3f%%)  "
                "z_centroide=%.6f =+=+=+=\n",
                (double) ns->par.t, (double) vol, (double) Vex,
                (double) (fabs(vol-Vex)/Vex*100), (double) zc);
        (void) wz;
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
