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
#include "../examples-common/malha-adaptativa.h"   // B4: medida de banda
#include "hig-flow-fronteira-imersa.h"              // fi_suporte_nivel_trocado
#include "hig-flow-front-tracking-3d.h"
#include <stdlib.h>

extern "C" void front_tracking_3d_instala(higflow_solver *ns, ft3_superficie *sup,
                                          real sigma);
extern "C" void ft3_escreve_amr_criterio(higflow_solver *ns, ft3_superficie *sup,
                                         real lx, real ly, real lz,
                                         int nx, int ny, int nz,
                                         int niveis, int cel_min,
                                         const char *caminho);

// Converte uma banda FISICA em numero de celulas finas.  Existe porque comparar
// resolucoes com `celmin` FIXO compara regioes refinadas de extensao fisica
// DIFERENTE em cada malha -- a mais fina refina metade do espaco --, e a esteira,
// onde o centroide e a velocidade se decidem, fica de fora na malha fina.  O
// piso de 5 do VOF vence sempre: abaixo dele a reconstrucao nao fecha.
static int _celmin_da_banda_3d(const char *var, int celmin,
                               real lx, int nx, int niveis)
{
    const char *b = getenv(var);
    if (b == NULL) return celmin;
    real banda = atof(b);
    if (!(banda > 0.0)) return celmin;
    real h_fino = lx / (real) nx;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    int novo = (int) (banda / h_fino + 0.5);
    if (novo < 5) {
        print0f("===> %s=%g pede %d celulas finas em h=%.5f; o piso de 5 vence. "
                "A banda efetiva sera' %.5f, nao %.5f\n",
                var, (double) banda, novo, (double) h_fino,
                (double) (5 * h_fino), (double) banda);
        novo = 5;
    }
    print0f("===> banda fisica %.5f em h_fino=%.5f  ==>  celmin = %d (era %d)\n",
            (double) banda, (double) h_fino, novo, celmin);
    return novo;
}

extern "C" void front_tracking_3d_fracvol_inicial(higflow_solver *ns,
                                                  ft3_superficie *sup);

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

// ---------------------------------------------------------------------------
// O problema deste exemplo, como um tipo em vez de oito funcoes soltas.
// Os corpos sao os mesmos; so' mudaram de lugar e perderam o prefixo get_.
// ---------------------------------------------------------------------------
class GotaProblem : public HigFlowProblem, public HigFlowMultiphaseProblem,
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
    real viscosity0(Point center, real t) { return _amb("FT3_MU0", 1.0); }
    real viscosity1(Point center, real t) { return _amb("FT3_MU1", 1.0); }
    real density0(Point center, real t)   { return _amb("FT3_RHO0", 1.0); }
    real density1(Point center, real t)   { return _amb("FT3_RHO1", 1.0); }

    // Fracao INICIAL.  Do segundo passo em diante quem manda e' a superficie,
    // pelo gancho; esta funcao so' serve a' condicao inicial, e usa a MESMA
    // esfera para as duas nascerem coerentes.
    real fracvol(Point center, Point delta, real t) {
        real cx = _amb("FT3_CX", 0.5), cy = _amb("FT3_CY", 0.5);
        real cz = _amb("FT3_CZ", 0.5), R  = _amb("FT3_R", 0.25);
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
    // Conformacao inicial: identidade, o equilibrio do Oldroyd-B em repouso.
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

static GotaProblem problema;

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

// MARCA DE PONTO DE SALVAMENTO COMPLETO.
//
// Os `._dp*` sao gravados em VARIOS arquivos, um por campo, SOBRESCREVENDO o
// conjunto anterior.  Se a maquina morrer no meio da gravacao o conjunto fica
// metade novo e metade velho -- estado remendado de dois instantes, que carrega
// adiante sem queixa nenhuma.  E' a unica falha deste mecanismo que CORROMPE em
// vez de so' perder, e por isso ela precisa ser DETECTAVEL.
//
// A marca e' gravada DEPOIS de todos eles, num temporario que e' depois
// renomeado -- o rename e' atomico no mesmo sistema de arquivos, de modo que a
// marca existe com o passo certo ou nao existe.  Quem retoma compara o passo da
// marca com o do init.yaml: se diferirem, a gravacao foi interrompida e aquele
// conjunto nao serve.  O `fsync` antes do rename e' o que impede a marca de
// chegar ao disco antes dos dados que ela atesta.
static void marca_salvamento(const char *base, int passo, double t) {
    char tmp[1200], fim[1200];
    snprintf(tmp, sizeof tmp, "%s.ok.tmp", base);
    snprintf(fim, sizeof fim, "%s.ok",     base);
    FILE *f = fopen(tmp, "w");
    if (f == NULL) return;
    fprintf(f, "passo %d\nt %.10f\n", passo + 1, t);   // +1: a convencao do init.yaml
    fflush(f);
    fsync(fileno(f));
    fclose(f);
    if (rename(tmp, fim) != 0) remove(tmp);
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

    // -------------------------------------------------------------------
    // B4 3D: REFINO ADAPTATIVO.  Parametros por ambiente, e a malha INICIAL ja'
    // adaptada -- ela tem de existir ANTES de higflow_create_domain, que a le'.
    //
    // FT3_AMR_NIVEIS  quantos niveis de refino (0 desliga tudo)
    // FT3_AMR_CELMIN  piso de celulas finas de cada lado da interface
    // FT3_AMR_BANDA   banda FISICA; converte-se em celmin, e vence o piso de 5
    // FT3_AMR_CADA    remalha a cada quantos passos (0 desliga o dinamico)
    // FT3_AMR_NX/NY/NZ e FT3_AMR_LX/LY/LZ   a malha BASE do criterio
    // FT3_AMR_CAMINHO o .amr que o criterio escreve e o solver rele'
    // -------------------------------------------------------------------
    // Remocao de ondulacoes submalha (BALANCE + TSUR).  O padrao LIGADO a
    // cada 20 passos e' deliberado -- foi o pedido --, e por isso o exemplo
    // ANUNCIA o intervalo no arranque: nenhuma corrida deve ter a superficie
    // suavizada sem que o log diga.
    int  tsur_cada = 20;
    real tsur_rlx  = 0.0;        // 0 reproduz o Freeflow (relaxacao total)
    // LIMIAR DE ONDULACAO, e ele NAO e' afinacao.  Sem limiar o filtro move
    // todo vertice em toda aplicacao e introduz tensao superficial artificial:
    // medido, 75 aplicacoes num elipsoide a volume fixo mudam a esfericidade em
    // 0,025 -- nove vezes a diferenca entre os dois metodos que a comparacao
    // mede.  Com limiar >= 0,02 o elipsoide liso fica INTOCADO e a ondulacao
    // ainda cai 82%.  0,05 fica com margem da borda, a 80%.
    real tsur_lim  = 0.05;
    { const char *e;
      if ((e = getenv("FT3_TSUR_CADA")) != NULL) tsur_cada = atoi(e);
      if ((e = getenv("FT3_TSUR_RLX"))  != NULL) tsur_rlx  = atof(e);
      if ((e = getenv("FT3_TSUR_LIM"))  != NULL) tsur_lim  = atof(e); }
    print0f("=+=+=+= Ondulacoes submalha: BALANCE+TSUR a cada %d passo(s), "
            "rlx=%g limiar=%g %s=+=+=+=\n", tsur_cada, (double) tsur_rlx,
            (double) tsur_lim,
            tsur_cada > 0 ? "" : "(DESLIGADO) ");

    int  amr_niveis = 0, amr_celmin = 5, amr_cada = 0;
    int  amr_nx = 25, amr_ny = 25, amr_nz = 50;
    real amr_lx = 1.0, amr_ly = 1.0, amr_lz = 2.0;
    const char *amr_caminho = "amrs/criterio.amr";
    {
        const char *e;
        if ((e = getenv("FT3_AMR_NIVEIS")) != NULL) amr_niveis = atoi(e);
        if ((e = getenv("FT3_AMR_CELMIN")) != NULL) amr_celmin = atoi(e);
        if ((e = getenv("FT3_AMR_CADA"))   != NULL) amr_cada   = atoi(e);
        if ((e = getenv("FT3_AMR_NX"))     != NULL) amr_nx = atoi(e);
        if ((e = getenv("FT3_AMR_NY"))     != NULL) amr_ny = atoi(e);
        if ((e = getenv("FT3_AMR_NZ"))     != NULL) amr_nz = atoi(e);
        if ((e = getenv("FT3_AMR_LX"))     != NULL) amr_lx = atof(e);
        if ((e = getenv("FT3_AMR_LY"))     != NULL) amr_ly = atof(e);
        if ((e = getenv("FT3_AMR_LZ"))     != NULL) amr_lz = atof(e);
        if ((e = getenv("FT3_AMR_CAMINHO"))!= NULL) amr_caminho = e;
        // DEPOIS de amr_lx e amr_nx, nao antes: a conversao da banda fisica
        // depende dos dois.  No 2D esta ordem foi trocada em UM dos dois
        // sitios de chamada, e o resultado foi celmin=10 em toda resolucao --
        // o que invalidou duas series inteiras sem dar erro nenhum.
        amr_celmin = _celmin_da_banda_3d("FT3_AMR_BANDA", amr_celmin,
                                         amr_lx, amr_nx, amr_niveis);
        if (amr_niveis > 0) {
            // RETOMADA: NAO reescrever o .amr quando a corrida esta' sendo
            // continuada.  O arquivo em disco JA' e' a malha correta -- ele e'
            // reescrito a cada remalhamento, e a malha so' muda ali.  Reescreve-lo
            // do estado INICIAL poria os campos salvos (gravados com a posicao
            // junto) numa malha diferente da que os gerou.
            if (ns->par.step == 0) {
                // Uma esfera TEMPORARIA, com os mesmos parametros da definitiva:
                // a de verdade so' nasce adiante, e a malha tem de existir antes.
                real cx0 = 0.5, cy0 = 0.5, cz0 = 0.5, r0 = 0.25;
                int  ns0 = 3;
                const char *q;
                if ((q = getenv("FT3_CX"))   != NULL) cx0 = atof(q);
                if ((q = getenv("FT3_CY"))   != NULL) cy0 = atof(q);
                if ((q = getenv("FT3_CZ"))   != NULL) cz0 = atof(q);
                if ((q = getenv("FT3_R"))    != NULL) r0  = atof(q);
                if ((q = getenv("FT3_NSUB")) != NULL) ns0 = atoi(q);
                Point c0; c0[0] = cx0; c0[1] = cy0; c0[2] = cz0;
                ft3_superficie *s0 = ft3_cria_esfera(c0, r0, ns0);
                if (s0 != NULL) {
                    ft3_escreve_amr_criterio(NULL, s0, amr_lx, amr_ly, amr_lz,
                                             amr_nx, amr_ny, amr_nz,
                                             amr_niveis, amr_celmin, amr_caminho);
                    ft3_destroi(s0);
                    print0f("===> malha inicial JA' adaptada, do criterio: %s\n",
                            amr_caminho);
                } else {
                    print0f("FT3: falha ao criar a esfera do criterio inicial\n");
                }
            } else {
                print0f("===> retomada: mantendo a malha adaptada de %s\n",
                        amr_caminho);
            }
        }
    }

    // Ponto de salvamento: DESLIGADO por omissao (ver o bloco no laco).
    int salva = 0;
    { const char *e = getenv("FT3_SALVA"); if (e != NULL) salva = atoi(e); }
    print0f("=+=+=+= Pontos de salvamento: %s =+=+=+=\n",
            salva ? "LIGADOS (corrida retomavel)" : "desligados");

    // Create the simulation domain
    // O multifasico e' ADICAO ao dominio base, nao substituto: sem
    // higflow_create_domain_multiphase os dominios sdmult e sdED nascem NULOS e
    // higflow_initialize_domain_yaml segmenta em sd_add_higtree.  Foi o que
    // aconteceu aqui, e o comentario do exemplo 2D ja' avisava.
    higflow_create_domain(ns, cache, order_center);
    if (ns->contr.flowtype == MULTIPHASE)
        higflow_create_domain_multiphase(ns, cache, order_center, &problema);
    if (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.viscoelastic_either == true)
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
    // RETOMADA: a superficie vem do DISCO, nao de uma esfera nova.  Criar a
    // esfera aqui numa corrida retomada poria a interface em t=0 e os campos no
    // tempo salvo -- dois instantes diferentes, e nada avisaria.
    ft3_superficie *gota = NULL;
    if (ns->par.step > 0) {
        char cam[1024];
        snprintf(cam, sizeof cam, "%s.frente", ns->par.namesave);
        gota = ft3_le(cam);
        if (gota != NULL)
            print0f("===> retomada: frente lida de %s, nv=%d nt=%d chi=%d\n",
                    cam, ft3_num_vertices(gota), ft3_num_triangulos(gota),
                    ft3_euler(gota));
        else
            print0f("FT3: nao li a frente de %s -- a retomada NAO pode "
                    "continuar com uma esfera nova\n", cam);
    }
    if (gota == NULL && ns->par.step > 0) exit(1);
    if (gota == NULL) gota = ft3_cria_esfera(centro_gota, R_gota, nsub);
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
    // A FRACAO INICIAL VEM DA SUPERFICIE, nao do degrau suavizado que a funcao
    // `fracvol` acima devolve.  Sem isto o quadro 0 -- o que fixa a massa de
    // referencia -- nasce com o erro do degrau (1.2e-02 em volume), e a
    // comparacao com o VOF, que inicializa por subdivisao subcelular (6.1e-05),
    // mediria a diferenca das condicoes iniciais em vez da dos metodos.
    front_tracking_3d_fracvol_inicial(ns, gota);
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

    // CADENCIA DO CAMPO, separada da da interface.  O VTK de campo e' um
    // UNSTRUCTURED_GRID com oito cantos por celula, sem partilhar vertices:
    // 250 mil celulas dao 2 milhoes de pontos.  Medido, em N=25 o campo sai com
    // 28 MB por quadro contra 124 KB da interface -- fator 226 --, e em N=50 ele
    // chega a ~300 MB.  Medido com o processo SOZINHO: 130 MB em menos de 10
    // minutos incluindo a montagem, da ordem de 20 MINUTOS por quadro de campo.
    // (A primeira estimativa, 75-100 min, saiu de uma taxa colhida com dois
    // processos competindo pela banda de memoria, e a medicao sozinho a
    // refutou.)  Ainda assim: 120 quadros dariam 40 horas so' de impressao,
    // contra ~16 de conta.
    //
    // A interface NAO precisa do campo: ela ja' sai com p, vel, S e T
    // interpolados nos vertices, que e' o que a animacao da bolha com tensoes
    // pede.  O campo serve a corte e a visualizacao volumetrica, para o que
    // tres amostras bastam (t=0, 1.5, 3.0 com dtp=0.025 e t final 3).  As
    // variaveis SB_CAMPO_CADA / FT3_CAMPO_CADA controlam: 1
    // imprime todo quadro (o comportamento antigo), 0 nunca imprime.
    int campo_cada = 60;
    { const char *e = getenv("FT3_CAMPO_CADA"); if (e != NULL) campo_cada = atoi(e); }
    print0f("=+=+=+= VTK de campo a cada %d quadro(s); interface em todos =+=+=+=\n",
            campo_cada);

    // Printing the properties to visualize: first step
    if (ns->par.step == 0) {
        print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
        if (campo_cada > 0 && ns->par.frame % campo_cada == 0) higflow_print_vtk(ns, myrank);
        // A SUPERFICIE, para a animacao.  O FT a carrega explicitamente, logo
        // nao ha' nada a extrair: grava-se a triangulacao como ela esta'.
        // (O analogo no VOF e' a isosuperficie por marching cubes, que precisa
        // ser extraida do campo porque la' a interface nao existe como objeto.)
        if (myrank == 0) {
            char cam[1024];
            snprintf(cam, sizeof cam, "%s_frente_%d-%d.vtk",
                     ns->par.nameprint, myrank, ns->par.frame);
            if (ft3_escreve_vtk(gota, cam) != 0)
                print0f("FT3: falha ao gravar %s\n", cam);
        }
        //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
        ns->par.tp += ns->par.dtp;
        ns->par.frame++;
        print0f("===> Saving               <====> ts = %15.10lf <===\n", ns->par.ts);
        //higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        //higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }
    
    // PRIMEIRO PONTO DE SALVAMENTO, e ele NAO e' opcional para quem quer
    // retomar.  `higflow_save_all_controllers_and_parameters_yaml` so' copia o
    // par.contr.yaml quando `step == initstep`; em qualquer passo posterior ela
    // verifica um flag global, descobre que a copia nunca foi feita, IMPRIME UM
    // AVISO E RETORNA SEM GRAVAR o init.yaml.
    //
    // Com as duas chamadas do bloco do passo 0 comentadas, era exatamente isso
    // que acontecia: os `._dp*` eram gravados a cada `dts` e o init.yaml --
    // que carrega step, t, frame, ts e tp -- NUNCA.  Ponto de salvamento sem
    // instante nenhum: nao se retoma de nada, e o unico sinal era uma linha de
    // aviso no meio do log.
    //
    // E tem de ser AQUI, fora do `if (step == 0)`, porque numa corrida RETOMADA
    // o bloco do passo 0 nao executa e initstep e' o passo carregado.  Se a
    // copia dependesse dele, a segunda retomada gravaria propriedades novas
    // sobre um init.yaml velho -- estado remendado de dois instantes, que e'
    // pior que nao ter ponto de salvamento.
    if (salva) {
        print0f("===> Ponto de salvamento inicial <====> step = %d  t = %15.10lf <===\n",
                ns->par.step, ns->par.t);
        higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        higflow_save_properties(ns, myrank, ntasks);
        // A FRENTE junto, senao este ponto fica incompleto: a retomada LE a
        // superficie do disco e ABORTA se nao achar (ver o ft3_le acima).
        if (myrank == 0) {
            char cam[1024];
            snprintf(cam, sizeof cam, "%s.frente", ns->par.namesave);
            if (ft3_grava(gota, cam) != 0)
                print0f("FT3: falha ao gravar a frente em %s\n", cam);
            marca_salvamento(ns->par.namesave, ns->par.step, ns->par.t);
        }
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
        // VISCOELASTICO: higflow_solver_step_multiphase_viscoelastic resolve a
        // equacao constitutiva, monta o tensor polimerico e ENTAO chama
        // higflow_solver_step_multiphase -- que e' onde o gancho da frente vive.
        // Logo a frente continua mandando na fracao, e a maquinaria polimerica
        // ja' a usa.  Front-tracking viscoelastico e' esta linha.
        //
        // UMA DEFASAGEM, dita para nao ser descoberta depois: a equacao
        // constitutiva roda ANTES do gancho, entao os parametros De e beta sao
        // interpolados com o fracvol do passo ANTERIOR.  E' O(dt), da mesma
        // ordem da defasagem que a forca de tensao ja' tem.
        // ONDULACOES SUBMALHA: BALANCE e TSUR, portadas do Freeflow.
        //
        // A adveccao lagrangiana acumula zigue-zague na escala da malha, e ele
        // aparece sobretudo na parte INFERIOR da bolha, onde o escoamento
        // comprime a superficie contra a esteira.  As duas rotinas o removem
        // PRESERVANDO VOLUME -- medido no arreio: a ondulacao cai 84% em tres
        // passadas com variacao de volume de 0,000%, contra 9,74% de perda da
        // media laplaciana que suaviza o mesmo tanto.
        //
        // ANTES do passo, de proposito: o gancho do front-tracking regenera a
        // fracao volumetrica da GEOMETRIA da superficie no passo seguinte, e
        // suavizar depois dele deixaria o campo um passo atrasado em relacao a'
        // superficie que o gerou.
        //
        // FT3_TSUR_CADA controla o intervalo; 0 desliga.  O padrao e' 20.
        if (tsur_cada > 0 && ns->par.step > 0 &&
            ns->par.step % tsur_cada == 0) {
            const real vol_antes = ft3_volume(gota);
            const real area_antes = ft3_area(gota);
            const int nb = ft3_balance(gota, tsur_lim);
            const int nt2 = ft3_tsur(gota, tsur_rlx, tsur_lim);
            const real vol_depois = ft3_volume(gota);
            const real area_depois = ft3_area(gota);
            print0f("===> TSUR passo %d: BALANCE moveu %d, TSUR %d, "
                    "chi=%d, volume %+.3e, AREA %+.3e\n",
                    ns->par.step, nb, nt2, ft3_euler(gota),
                    (double) ((vol_depois - vol_antes) / vol_antes),
                    (double) ((area_depois - area_antes) / area_antes));
        }

        // REMALHAMENTO ADAPTATIVO, pelo caminho verificado da fronteira
        // imersa.  A superficie escreve o .amr multinivel pelo criterio, e
        // higflow_reconstroi_dominio o rele' pelo caminho normal de arranque.
        //
        // A SUPERFICIE NAO E' TRANSFERIDA: ela e' replicada e independente de
        // malha, sobrevive intocada ao remalhamento e regenera o fracvol no
        // passo seguinte, pelo gancho.  So' u e p viajam.
        if (amr_niveis > 0 && amr_cada > 0 && ns->par.step > 0 &&
            ns->par.step % amr_cada == 0) {
            const bool crono = (getenv("HIGFLOW_CRONO_REMALHA") != NULL);
            double t0 = 0.0, t1 = 0.0, t2 = 0.0;
            if (crono) { MPI_Barrier(MPI_COMM_WORLD); t0 = MPI_Wtime(); }
            ft3_escreve_amr_criterio(ns, gota, amr_lx, amr_ly, amr_lz,
                                     amr_nx, amr_ny, amr_nz,
                                     amr_niveis, amr_celmin, amr_caminho);
            if (crono) { MPI_Barrier(MPI_COMM_WORLD); t1 = MPI_Wtime(); }
            long faltam = higflow_reconstroi_dominio(ns, ntasks, myrank,
                                                     cache, order_center,
                                                     order_facet);
            if (crono) {
                MPI_Barrier(MPI_COMM_WORLD); t2 = MPI_Wtime();
                int nt = 1; MPI_Comm_size(MPI_COMM_WORLD, &nt);
                print0f("CRONO np=%d passo %d: criterio %.4f s  "
                        "reconstrucao %.4f s  total %.4f s\n",
                        nt, ns->par.step, t1 - t0, t2 - t1, t2 - t0);
            }
            if (faltam != 0)
                print0f("===> remalha passo %d: %ld posicoes sem valor\n",
                        ns->par.step, faltam);
        }

        if (ns->contr.flowtype == MULTIPHASE) {
            if (ns->ed.mult.contr.viscoelastic_either == true)
                higflow_solver_step_multiphase_viscoelastic(ns);
            else
                higflow_solver_step_multiphase(ns);
        } else {
            higflow_solver_step(ns);
        }
        // ORACULOS DA MALHA ADAPTADA, DEPOIS do passo -- antes dele o fracvol
        // esta' zerado pelo remalhamento e nao ha' semente que colher.
        //   banda   -- a regra das celulas minimas foi CUMPRIDA de fato?  Tem de
        //              dar >= FT3_AMR_CELMIN.  O criterio garante no PAPEL; isto
        //              mede na malha que o solver de fato montou.
        //   nivel   -- pontos de suporte do nucleo caidos em celula de TAMANHO
        //              DIFERENTE do h do marcador.  DEVE SER ZERO: o nucleo so'
        //              e' normalizado para um h, e atravessar nivel quebra a
        //              particao da unidade EM SILENCIO, porque a soma dos pesos
        //              pode continuar perto de 1.
        if (getenv("FT3_DIAG_MALHA") != NULL && amr_niveis > 0 &&
            ns->par.step % (amr_cada > 0 ? amr_cada : 25) == 0) {
            real banda = malha_adapt_mede_banda(ns, amr_niveis,
                                                amr_lx / amr_nx);
            print0f("===> MALHA passo %d: banda=%.2f celulas finas "
                    "(exigido >=%d)  suporte_nivel_trocado=%ld  perdidos=%ld\n",
                    ns->par.step, (double) banda, amr_celmin,
                    fi_suporte_nivel_trocado(), fi_suporte_perdidos());
        }

        // Time update 
        ns->par.t += ns->par.dt;
        // Stop the first step time
        if (ns->par.step == step0) STOP_CLOCK(firstiter); 
        // Printing
        if (ns->par.t >= ns->par.tp) {
            print0f("===> Printing frame: %4d <====> tp = %15.10lf <===\n",ns->par.frame, ns->par.tp);
            if (campo_cada > 0 && ns->par.frame % campo_cada == 0) higflow_print_vtk(ns, myrank);
            if (myrank == 0) {
                char cam[1024];
                snprintf(cam, sizeof cam, "%s_frente_%d-%d.vtk",
                         ns->par.nameprint, myrank, ns->par.frame);
                if (ft3_escreve_vtk(gota, cam) != 0)
                    print0f("FT3: falha ao gravar %s\n", cam);
            }
            //higflow_print_vtk2D_parallel_single(ns, myrank, ntasks);
            ns->par.tp += ns->par.dtp;
            ns->par.frame++;
        }
        // PONTO DE SALVAMENTO, para a corrida poder ser RETOMADA.
        //
        // Estava comentado, e por isso 37 horas de uma corrida de referencia se
        // perderam quando ela teve de ser interrompida: o estado existia so' em
        // memoria.  O mecanismo de retomada JA' estava pronto do outro lado --
        // higflow_load_properties e' chamado quando init_par.step > 0 -- e o
        // yaml de retomada sai do save_all_controllers.  Faltava so' gravar.
        //
        // Os TRES juntos, e a ordem importa: o yaml leva step, t, frame, ts e
        // tp DESTE instante, e as propriedades vao em seguida.  Retomar com um
        // yaml de um instante e campos de outro poria os valores num tempo que
        // nao os gerou.
        //
        // Controlado por ambiente porque gravar u, p e os tensores de 250 mil
        // celulas a cada `dts` custa disco e tempo, e ligar por omissao mudaria
        // o custo de toda corrida existente sem aviso.
        if (salva && ns->par.t >= ns->par.ts) {
            print0f("===> Salvando (retomavel)  <====> ts = %15.10lf <===\n",
                    ns->par.ts);
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            // A SUPERFICIE TAMBEM, e sem ela a retomada nao existe: ela nao e'
            // parte do estado do solver -- e' replicada, independente de malha,
            // e e' quem regenera a fracao volumetrica a cada passo.  Retomar com
            // u e p salvos e uma esfera nova poria o campo e a interface em
            // instantes diferentes.
            if (myrank == 0) {
                char cam[1024];
                snprintf(cam, sizeof cam, "%s.frente", ns->par.namesave);
                if (ft3_grava(gota, cam) != 0)
                    print0f("FT3: falha ao gravar a frente em %s\n", cam);
            }
            if (myrank == 0)
                marca_salvamento(ns->par.namesave, ns->par.step, ns->par.t);
            ns->par.ts += ns->par.dts;
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
