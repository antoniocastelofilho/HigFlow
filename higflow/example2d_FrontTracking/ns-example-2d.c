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

// FRONT-TRACKING (fase B2, gota estatica / lei de Laplace).  NAO e' o corpo
// rigido: a frente ORDENADA carrega tensao superficial sigma*kappa*n, espalhada
// pelo adaptador examples-common/front-tracking.c.
#include "../src/hig-flow-front-tracking.h"
#include "../examples-common/malha-adaptativa.h"   // B4: refino adaptativo
#include "../examples-common/hysing-metricas.h"      // grandezas do benchmark
#include "../src/hig-flow-fronteira-imersa.h"       // contadores de suporte
extern "C" void front_tracking_passo_inicial(int passo);
extern "C" void front_tracking_instala(higflow_solver *ns, ft_frente *frente,
                                        real sigma);
extern "C" void ft_escreve_amr_criterio(higflow_solver *ns, ft_frente *frente,
                                        real lx, real ly,
                                        int nx, int ny, int niveis, int cel_min,
                                        const char *caminho);

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
class NewtProblem : public HigFlowProblem, public HigFlowMultiphaseProblem {
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
                        // GOTA ESTATICA (Laplace): sem entrada.  A unica
                        // dinamica e' a tensao superficial; a velocidade na
                        // parada tem de ficar ~0 (correntes parasitas).
                        value = 0.0;
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

    // --- interface multifasica (fase B4) ---
    // Propriedades por ambiente: FT_RHO0/FT_RHO1 e FT_MU0/FT_MU1.  Os padroes
    // sao IGUAIS nas duas fases, que e' o estagio B4.0: mover o encanamento para
    // o caminho multifasico sem mudar a fisica, de modo que qualquer diferenca
    // no resultado seja do encanamento e nao do salto de propriedade.
    //
    // CONVENCAO DO SOLVER, conferida no codigo e NAO suposta:
    //   dens = (1 - fracvol)*dens0 + fracvol*dens1   (hig-flow-step-multiphase.c)
    // ou seja fracvol=1 -> fase 1.  FASE 1 E' DENTRO da gota; FASE 0 e' FORA.
    // Eu havia suposto o contrario, e o resultado foi uma bolha PESADA num meio
    // leve: ela afundou ate' encostar na base, e a quebra apareceu primeiro na
    // particao da unidade (6,7e-15 -> 0,53), porque os marcadores chegaram a'
    // fronteira do dominio.  Quem localizou foi o despejo das FORMAS da frente;
    // os numeros sozinhos so' diziam que o deslocamento parava de crescer.
    static real _amb(const char *nome, real padrao) {
        const char *s = getenv(nome);
        return (s != NULL) ? atof(s) : padrao;
    }
    real viscosity0(Point center, real t) { return _amb("FT_MU0", 1.0); }
    real viscosity1(Point center, real t) { return _amb("FT_MU1", 1.0); }
    real density0(Point center, real t)   { return _amb("FT_RHO0", 1.0); }
    real density1(Point center, real t)   { return _amb("FT_RHO1", 1.0); }

    // A fracao INICIAL.  Depois do primeiro passo quem manda e' a frente, pelo
    // gancho -- esta funcao so' serve a' condicao inicial.  Usa a mesma elipse
    // (ou circulo) que a frente, para que as duas nascam coerentes.
    real fracvol(Point center, Point delta, real t) {
        real cx = _amb("FT_CX", 2.0), cy = _amb("FT_CY", 0.0);
        real a  = _amb("FT_A", 0.0),  b  = _amb("FT_B", 0.0);
        if (a <= 0.0 || b <= 0.0) { a = _amb("FT_R", 0.25); b = a; }
        // Subamostragem 16x16, como o exemplo VOF faz: a celula cortada recebe a
        // fracao integrada, nao o teste do centro.
        const int N = 16;
        int dentro = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                real x = center[0] - 0.5*delta[0] + (i + 0.5)*delta[0]/N;
                real y = center[1] - 0.5*delta[1] + (j + 0.5)*delta[1]/N;
                real dx = (x - cx)/a, dy = (y - cy)/b;
                if (dx*dx + dy*dy <= 1.0) dentro++;
            }
        return (real) dentro / (real) (N*N);
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
    // MALHA INICIAL JA' ADAPTADA.  A frente e' GEOMETRIA PURA -- nao depende de
    // dominio nenhum --, entao pode nascer antes e ditar a malha do arranque.
    //
    // Sem isto o primeiro remalhamento saltava da malha uniforme para DOIS
    // niveis de uma vez, e a transferencia por posicao nao preenchia: medido,
    // 8820 posicoes sem valor, com a area da gota caindo 1,4% em 100 passos.
    // A funcao so' promete preenchimento completo quando a malha muda de UM
    // nivel -- que e' o que uma adaptacao por passo produz.
    {
        const char *e; int niv = 0;
        if ((e = getenv("FT_AMR_NIVEIS")) != NULL) niv = atoi(e);
        if (niv > 0) {
            real cx0 = 2.0, cy0 = 0.0, r0 = 0.25;
            int nm0 = 128;
            if ((e = getenv("FT_CX"))    != NULL) cx0 = atof(e);
            if ((e = getenv("FT_CY"))    != NULL) cy0 = atof(e);
            if ((e = getenv("FT_R"))     != NULL) r0  = atof(e);
            if ((e = getenv("FT_NMARC")) != NULL) nm0 = atoi(e);
            int cel0 = 5;   if ((e = getenv("FT_AMR_CELMIN")) != NULL) cel0 = atoi(e);
            int nx0 = 40;   if ((e = getenv("FT_AMR_NX")) != NULL) nx0 = atoi(e);
            int ny0 = 80;   if ((e = getenv("FT_AMR_NY")) != NULL) ny0 = atoi(e);
            real lx0 = 1.0; if ((e = getenv("FT_AMR_LX")) != NULL) lx0 = atof(e);
            real ly0 = 2.0; if ((e = getenv("FT_AMR_LY")) != NULL) ly0 = atof(e);
            const char *cam0 = getenv("FT_AMR_CAMINHO");
            if (cam0 == NULL) cam0 = "amrs-hysing/criterio/dominio.amr";
            // RETOMADA: NAO reescrever o .amr quando a corrida esta' sendo
            // continuada.  O arquivo em disco JA' e' a malha correta: ele e'
            // reescrito a cada remalhamento, e a malha so' muda ali -- entao em
            // qualquer ponto de salvamento o .amr em disco E' a malha em
            // memoria.  Reescreve-lo do estado INICIAL poria os campos salvos
            // (que sao gravados com a posicao junto) numa malha diferente da
            // que os gerou: celula que sumiu perde o valor, celula que apareceu
            // nao recebe nenhum.
            if (ns->par.step == 0) {
                Point c0; c0[0] = cx0; c0[1] = cy0;
                for (int d = 2; d < DIM; d++) c0[d] = 0.0;
                ft_frente *f0 = ft_cria_circulo(c0, r0, nm0);
                ft_escreve_amr_criterio(NULL, f0, lx0, ly0, nx0, ny0, niv, cel0, cam0);
                ft_destroi(f0);
                print0f("===> malha inicial JA' adaptada, do criterio\n");
            } else {
                print0f("===> retomada: mantendo a malha adaptada de %s\n", cam0);
            }
        }
    }

    // B4: o multifasico e' uma ADICAO ao dominio base, nao um substituto -- e'
    // assim que o example2d_VOF faz, e trocar um pelo outro custou um SEGV.
    higflow_create_domain(ns, cache, order_center);
    if (ns->contr.flowtype == MULTIPHASE)
        higflow_create_domain_multiphase(ns, cache, order_center, &problema);
    
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
    // A GOTA.  Frente circular de raio 0,25 em (2,0), no mesmo canal e h=0,05 do
    // exemplo do corpo rigido.  A tensao superficial sigma impoe o salto de
    // Laplace Delta p = sigma/R atraves da interface.  A frente e' REPLICADA e
    // ORDENADA (o front-tracking precisa da ordem para a curvatura).
    //
    // Para o teste de Laplace a gota fica em REPOUSO: a forca sigma*kappa*n tem
    // de equilibrar o gradiente de pressao, com u ~ 0 (correntes parasitas
    // pequenas).  sigma via codigo (nao ha' campo bifasico aqui -- e' um fluido
    // so' com forca singular na interface).
    // Geometria da gota por ambiente, para poder rodar o MESMO binario no setup
    // do example2d_VOF e comparar os dois metodos no mesmo problema:
    //   FT_R, FT_CX, FT_CY, FT_SIGMA, FT_NMARC
    // Os padroes reproduzem o caso original (R=0,25 em (2,0), sigma=1).
    ft_frente *gota = NULL;
    real R_gota = 0.25, cx_gota = 2.0, cy_gota = 0.0;
    {
        const char *s;
        if ((s = getenv("FT_R"))     != NULL) R_gota  = atof(s);
        if ((s = getenv("FT_CX"))    != NULL) cx_gota = atof(s);
        if ((s = getenv("FT_CY"))    != NULL) cy_gota = atof(s);
        real sigma = 1.0;
        if ((s = getenv("FT_SIGMA")) != NULL) sigma = atof(s);
        real ea = 0.0, eb = 0.0;
        if ((s = getenv("FT_A")) != NULL) ea = atof(s);
        if ((s = getenv("FT_B")) != NULL) eb = atof(s);
        // NUMERO DE MARCADORES DERIVADO DA MALHA, nao escolhido a mao.
        //
        // O front-tracking exige Delta s ~ h.  Fixar o numero de marcadores e
        // depois refinar a malha quebra essa relacao em silencio: MEDIDO, a
        // mesma gota com 64 marcadores perde 2,1e-2 de area em h=1/160
        // (Delta s/h = 3,93) e 2,3e-5 com 256 (Delta s/h = 0,98) -- 920 vezes.
        // Nao era o AMR, nem o remalhamento, nem o passo de tempo: era a frente
        // ficar relativamente grossa ao refinar, e a massa vazar entre
        // marcadores.
        //
        // FT_NMARC continua existindo e VENCE quando dado, para reproduzir
        // corridas antigas; sem ele o numero sai de FT_DS_SOBRE_H (padrao 1,0)
        // vezes o h da malha mais fina que a frente vai encontrar.
        int nmarc = 0;
        if ((s = getenv("FT_NMARC")) != NULL) nmarc = atoi(s);
        if (nmarc <= 0) {
            real h_base = 0.025, ds_h = 1.0;
            int niv = 0;
            if ((s = getenv("FT_AMR_LX")) != NULL && getenv("FT_AMR_NX") != NULL)
                h_base = atof(s) / atoi(getenv("FT_AMR_NX"));
            if ((s = getenv("FT_AMR_NIVEIS"))   != NULL) niv  = atoi(s);
            if ((s = getenv("FT_DS_SOBRE_H"))   != NULL) ds_h = atof(s);
            real h_fino = h_base;
            for (int l = 0; l < niv; l++) h_fino *= 0.5;
            real raio = (ea > 0.0 && eb > 0.0) ? sqrt(ea * eb) : R_gota;
            nmarc = (int) ceil(2.0 * M_PI * raio / (ds_h * h_fino));
            if (nmarc < 16) nmarc = 16;
            print0f("===> frente: %d marcadores, de Delta s/h=%.2f com "
                    "h_fino=%.5f\n", nmarc, (double) ds_h, (double) h_fino);
        }

        Point centro; centro[0] = cx_gota; centro[1] = cy_gota;
        for (int d = 2; d < DIM; d++) centro[d] = 0.0;

        // ELIPSE (FT_A e FT_B): a gota deixa de estar em equilibrio, e a tensao
        // superficial a relaxa para o circulo de MESMA AREA.  E' o teste que
        // mede o que a gota circular esconde: com kappa VARIAVEL ao longo da
        // frente, sigma*kappa*grad(H) deixa de ser exatamente grad(sigma*kappa*H)
        // e sobra um residuo de balanco.  Tambem e' o primeiro caso em que a
        // frente deforma, e portanto o primeiro que exercita a CIRURGIA com o
        // solver no circuito.
        if (ea > 0.0 && eb > 0.0) {
            Point *v = (Point *) malloc((size_t) nmarc * sizeof(Point));
            for (int k = 0; k < nmarc; k++) {
                real th = 2.0 * M_PI * (real) k / (real) nmarc;
                v[k][0] = centro[0] + ea * cos(th);
                v[k][1] = centro[1] + eb * sin(th);
                for (int d = 2; d < DIM; d++) v[k][d] = 0.0;
            }
            // ds_alvo = o espacamento medio dos vertices que acabei de gerar.
            // Assim `ft_cria_curva` NAO reamostra (nsub=1 em todo segmento) e os
            // marcadores ficam EXATAMENTE sobre a elipse -- tres pontos sobre
            // uma corda dariam curvatura zero --, mas o alvo fica registrado
            // para a cirurgia usar quando a frente deformar.
            real per = 0.0;
            for (int k = 0; k < nmarc; k++) {
                const real *p = v[k], *q = v[(k + 1) % nmarc];
                per += sqrt((q[0]-p[0])*(q[0]-p[0]) + (q[1]-p[1])*(q[1]-p[1]));
            }
            gota = ft_cria_curva(v, nmarc, per / nmarc * 1.0001);
            free(v);
            R_gota = sqrt(ea * eb);       // raio de equilibrio: mesma area
            print0f("=+=+=+= Front-tracking: ELIPSE a=%.5f b=%.5f, area=%.6f, "
                    "R_eq=%.5f, Laplace final esperado Dp=%.4f =+=+=+=\n",
                    (double) ea, (double) eb, (double) ft_area(gota),
                    (double) R_gota, (double)(sigma / R_gota));
        } else {
            gota = ft_cria_circulo(centro, R_gota, nmarc);
        }
        // RETOMADA.  Se a corrida esta' sendo continuada de um `h.save`
        // (init_par.step > 0), a frente vem do ARQUIVO, e nao da geometria
        // inicial que acabou de ser construida.
        //
        // Sem isto a retomada seria INCOERENTE e silenciosa: os campos voltam no
        // instante salvo e os marcadores no circulo de t=0.  O solver nao tem
        // como notar -- a frente regenera o fracvol da propria geometria, o
        // arranque fica plausivel, e a corrida segue medindo outra coisa.  Por
        // isso a falta do arquivo ABORTA em vez de cair na geometria inicial:
        // continuar aqui e' pior do que nao continuar.
        if (ns->par.step > 0) {
            char nomef[1024];
            snprintf(nomef, sizeof nomef, "%s.frente", ns->par.namesave);
            ft_frente *lida = ft_le(nomef);
            if (lida == NULL) {
                print0f("=+=+=+= ERRO: retomada no passo %d sem a frente "
                        "gravada em %s.  A frente nao e' campo e nao sai do "
                        "h.save; continuar poria os marcadores em t=0 com os "
                        "campos em t=%g. =+=+=+=\n",
                        ns->par.step, nomef, (double) ns->par.t);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            ft_destroi(gota);
            gota = lida;
            print0f("===> frente RETOMADA de %s: %d marcadores, "
                    "area=%.8f, ds_alvo=%.6f\n", nomef, ft_num(gota),
                    (double) ft_area(gota), (double) ft_ds_alvo(gota));
        }
        // O passo em que a corrida comeca, para o despejo nomear os arquivos
        // pelo passo real e nao reiniciar a contagem numa retomada.
        front_tracking_passo_inicial(ns->par.step);
        front_tracking_instala(ns, gota, sigma);
        // A FRACAO TEM FONTE EXTERNA: a frente a regenera da geometria a cada
        // passo, e e' replicada -- sobrevive intocada ao remalhamento.  E' a
        // condicao para higflow_reconstroi_dominio aceitar MULTIPHASE.
        if (ns->contr.flowtype == MULTIPHASE)
            higflow_set_fracvol_externo(ns, true);
        print0f("=+=+=+= Front-tracking: gota R=%.5f em (%.3f,%.3f), sigma=%.3f, "
                "%d marcadores, Laplace esperado Dp=%.4f =+=+=+=\n",
                (double) R_gota, (double) cx_gota, (double) cy_gota,
                (double) sigma, nmarc, (double)(sigma / R_gota));
    }

    // Initialize the boundaries
    print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    //higflow_initialize_boundaries(ns);
    higflow_initialize_boundaries_yaml(ns);

    // Parametros do refino adaptativo, por ambiente.
    int  amr_niveis = 0, amr_celmin = 5, amr_cada = 0, amr_nx = 40, amr_ny = 80;
    real amr_lx = 1.0, amr_ly = 2.0;
    const char *amr_caminho = "amrs-hysing/criterio/dominio.amr";
    {
        const char *e;
        if ((e = getenv("FT_AMR_NIVEIS")) != NULL) amr_niveis = atoi(e);
        if ((e = getenv("FT_AMR_CELMIN")) != NULL) amr_celmin = atoi(e);
        if ((e = getenv("FT_AMR_CADA"))   != NULL) amr_cada   = atoi(e);
        if ((e = getenv("FT_AMR_NX"))     != NULL) amr_nx     = atoi(e);
        if ((e = getenv("FT_AMR_NY"))     != NULL) amr_ny     = atoi(e);
        if ((e = getenv("FT_AMR_LX"))     != NULL) amr_lx     = atof(e);
        if ((e = getenv("FT_AMR_LY"))     != NULL) amr_ly     = atof(e);
        if ((e = getenv("FT_AMR_CAMINHO"))!= NULL) amr_caminho = e;
    }

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
        if (gota != NULL) {   // a frente vai junto: o h.save nao a cobre
            char nomef[1024];
            snprintf(nomef, sizeof nomef, "%s.frente", ns->par.namesave);
            ft_grava(gota, nomef);
        }
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
        // B4: o passo multifasico e' outra rotina, e e' nela que o gancho do
        // front-tracking foi instalado -- chamar o monofasico com config
        // multifasica roda sem forca nenhuma, em silencio.
        // REMALHAMENTO ADAPTATIVO, pelo caminho verificado da fronteira imersa.
        // A frente escreve o .amr multinivel pelo criterio (o MESMO do VOF:
        // celulas de interface como sementes, distancia define o nivel), com os
        // limiares DERIVADOS da regra das celulas minimas; e entao o
        // higflow_reconstroi_dominio o rele' pelo caminho normal de arranque.
        //
        // A frente NAO e' transferida: ela e' replicada e independente de malha,
        // e regenera o fracvol no passo seguinte.  So' u e p viajam.
        if (amr_niveis > 0 && amr_cada > 0 && ns->par.step > 0 &&
            ns->par.step % amr_cada == 0) {
            // FASE 0 do escalonamento (HIGFLOW_CRONO_REMALHA): cronometra o
            // criterio e a reconstrucao SEPARADOS do resto.  A pergunta e' se o
            // remalhamento domina, e ela so' se responde medindo os dois lados.
            // Barreira antes de cada marca: sem ela mede-se desbalanceamento
            // acumulado, nao o custo da etapa.
            const bool crono = (getenv("HIGFLOW_CRONO_REMALHA") != NULL);
            double t0 = 0.0, t1 = 0.0, t2 = 0.0;
            if (crono) { MPI_Barrier(MPI_COMM_WORLD); t0 = MPI_Wtime(); }
            ft_escreve_amr_criterio(ns, gota, amr_lx, amr_ly, amr_nx, amr_ny,
                                    amr_niveis, amr_celmin, amr_caminho);
                        if (crono) { MPI_Barrier(MPI_COMM_WORLD); t1 = MPI_Wtime(); }
            long faltam = higflow_reconstroi_dominio(ns, ntasks, myrank,
                                                     cache, order_center,
                                                     order_facet);
            if (crono) {
                MPI_Barrier(MPI_COMM_WORLD); t2 = MPI_Wtime();
                int nt = 1; MPI_Comm_size(MPI_COMM_WORLD, &nt);
                print0f("CRONO np=%d passo %d: criterio %.4f s  reconstrucao %.4f s  total %.4f s\n",
                        nt, ns->par.step, t1 - t0, t2 - t1, t2 - t0);
            }
            if (faltam != 0)
                print0f("===> remalha passo %d: %ld posicoes sem valor\n",
                        ns->par.step, faltam);
        }

        if (ns->contr.flowtype == MULTIPHASE)
            higflow_solver_step_multiphase(ns);
        else
            higflow_solver_step(ns);
        // ORACULOS DA MALHA ADAPTADA, DEPOIS do passo -- antes dele o fracvol esta'
        // zerado pelo remalhamento e nao ha' semente que colher (FT_DIAG_MALHA):
        //   banda   -- a regra das celulas minimas foi CUMPRIDA de fato?  Tem de
        //              dar >= FT_AMR_CELMIN.  Nunca foi medida ate' agora.
        //   nivel   -- pontos de suporte do nucleo caidos em celula de TAMANHO
        //              DIFERENTE do h do marcador.  DEVE SER ZERO: o nucleo so'
        //              e' normalizado para um h, e atravessar nivel quebra a
        //              particao da unidade -- em silencio, porque a soma dos
        //              pesos pode continuar perto de 1.
        if (getenv("FT_DIAG_MALHA") != NULL && amr_niveis > 0 &&
            ns->par.step % (amr_cada > 0 ? amr_cada : 25) == 0) {
            real banda = malha_adapt_mede_banda(ns, amr_niveis,
                                                amr_lx / amr_nx);
            print0f("===> MALHA passo %d: banda=%.2f celulas finas "
                    "(exigido >=%d)  suporte_nivel_trocado=%ld  perdidos=%ld\n",
                    ns->par.step, (double) banda, amr_celmin,
                    fi_suporte_nivel_trocado(), fi_suporte_perdidos());
        }

        // GRANDEZAS DO BENCHMARK DE HYSING (FT_HYSING=1), nas definicoes do
        // artigo: yc e Vc integrados no campo de fracao -- a MESMA formula que o
        // VOF usa --, e a circularidade Pa/Pb com o perimetro exato do poligono.
        if (getenv("FT_HYSING") != NULL) {
            static int cada_h = 0;
            if (cada_h == 0) { const char *c = getenv("FT_HYSING_CADA");
                               cada_h = (c != NULL) ? atoi(c) : 50;
                               if (cada_h < 1) cada_h = 50; }
            if (ns->par.step % cada_h == 0) {
                real Ah, ych, vch;
                hysing_medidas(ns, &Ah, &ych, &vch);
                real Ap = ft_area(gota), Pp = ft_perimetro(gota);
                print0f("HYSING %10.5f %12.6f %12.6f %12.6f %12.8f\n",
                        (double) ns->par.t, (double) ych, (double) vch,
                        (double) hysing_circularidade(Ap, Pp), (double) Ah);
            }
        }

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
            if (gota != NULL) {   // a frente vai junto: o h.save nao a cobre
                char nomef[1024];
                snprintf(nomef, sizeof nomef, "%s.frente", ns->par.namesave);
                ft_grava(gota, nomef);
            }
            ns->par.ts += ns->par.dts;
        }
    }
    // ********************************************************
    // End Loop for the Navier-Stokes equations integration
    // ********************************************************

    // O ORACULO DO B2: o salto de Laplace.  ANTES de destruir o solver: amostra
    // a pressao no centro da gota (dentro) e num ponto longe (fora), e compara
    // Dp = p_in - p_out com sigma/R.
    {
        // Dentro = centro da gota; fora = deslocado 2,4R em y (bem fora da gota,
        // e dentro do dominio nos dois setups).  FT_POUT_Y sobrepoe se preciso.
        Point p_in, p_out;
        p_in[0]  = cx_gota;  p_in[1]  = cy_gota;
        p_out[0] = cx_gota;  p_out[1] = cy_gota + 2.4 * R_gota;
        { const char *s = getenv("FT_POUT_Y"); if (s != NULL) p_out[1] = atof(s); }
        for (int d = 2; d < DIM; d++) { p_in[d] = 0.0; p_out[d] = 0.0; }
        sim_stencil *stn = stn_create();
        real pin = 0.0, pout = 0.0;
        hig_cell *ci = sd_get_cell_with_point(ns->sdp, p_in);
        hig_cell *co = sd_get_cell_with_point(ns->sdp, p_out);
        if (ci != NULL) {
            Point cc; hig_get_center(ci, cc);
            pin = compute_value_at_point(ns->sdp, cc, p_in, 1.0, ns->dpp, stn);
        }
        if (co != NULL) {
            Point cc; hig_get_center(co, cc);
            pout = compute_value_at_point(ns->sdp, cc, p_out, 1.0, ns->dpp, stn);
        }
        stn_destroy(stn);
        real sigma = 1.0;
        { const char *s = getenv("FT_SIGMA"); if (s != NULL) sigma = atof(s); }
        const real dp_exato = sigma / R_gota;
        print0f("=+=+=+= LAPLACE  p_in=%.6f  p_out=%.6f  Dp=%.6f  "
                "sigma/R=%.6f  erro_rel=%.4f =+=+=+=\n",
                (double) pin, (double) pout, (double)(pin - pout),
                (double) dp_exato,
                (double) fabs((pin - pout) - dp_exato) / dp_exato);
    }
    (void) gota;

    // Destroy the Navier-Stokes object
    higflow_destroy(ns);
    // Stop the total time
    STOP_CLOCK(total);

    if(myrank == 0) {
        DEBUG_INSPECT(GET_NSEC_CLOCK(total)/1.0e9, %lf);
        DEBUG_INSPECT(GET_NSEC_CLOCK(firstiter)/1.0e9, %lf);
    }
}
