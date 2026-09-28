// *******************************************************************
// *******************************************************************
//  Example for HiG-Flow Solver - version 10/11/2016
// *******************************************************************
// *******************************************************************
//
// Two Newtonian fluids with an interface tracked by volume fraction.
//
// READING THE CONFIG: what makes this multiphase is `flowphase: multiphase`.  Under
// it the `flowtype:` field of the singlephase block is INERT -- the ones that count
// are `flowtype0` and `flowtype1`, one per phase.  Reading `flowtype:` alone
// misleads; here both phases are newtonian.

#include "ns-example-2d.h"
#include "../examples-common/malha-adaptativa.h"   // oraculo da banda
#include "../examples-common/hysing-metricas.h"      // grandezas do benchmark

#ifdef HIGFLOW_COM_T8CODE
// Definida em ../examples-common/malha-t8.cxx, compilada so' com T8CODE.
extern "C" void malha_t8_instala(higflow_solver *ns, int myrank);
#endif
#define pi 3.1415926535897932384626434

// *******************************************************************
// Extern functions for the Navier-Stokes program
// *******************************************************************

real func (Point p) {
	real value;
	Point c, pc, r, cav;
	// c[0]     = 0.5;
	// c[1]     = 0.75;
	// r[0]     = 0.15;
	// r[1]     = 0.15;
	// cav[0]   = 0.05;
	// cav[1]   = 0.2479;
	// pc[0]    = p[0] - c[0];
	// pc[1]    = p[1] - c[1];
	// // negative when in the circle
	// real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
	// real incav1  = (fabs(pc[0]) - 0.5*cav[0])/r[0]; ;
	// real incav2  = (pc[1] - (cav[1] - r[1]))/r[1];
	// // negative when in the cavity - both incav1 and incav2 must be negative
	// real incav   = max(incav1, incav2);
	// value = min(outcirc, incav);

	c[0] = 0.5; c[1] = 0.5;
	// Raios por ambiente (VOF_RA, VOF_RB): a formula ja' e' de ELIPSE -- os dois
	// raios so' estavam iguais.  Com raios distintos a gota nao esta' em
	// equilibrio e a tensao superficial a relaxa para o circulo de mesma area,
	// que e' o teste de curvatura VARIAVEL.  Sem as variaveis, o caso original.
	r[0] = 1.0/6.0; r[1] = 1.0/6.0;
	{ const char *e;
	  if ((e = getenv("VOF_RA")) != NULL) r[0] = atof(e);
	  if ((e = getenv("VOF_RB")) != NULL) r[1] = atof(e); }
	pc[0] = p[0] - c[0]; pc[1] = p[1] - c[1];
	real outcirc = 1.0 - (pc[0]*pc[0]/(r[0]*r[0]) + pc[1]*pc[1]/(r[1]*r[1]));
	value = outcirc;
	return value;
}

// Geometria 2-D comum aos exemplos VOF: as seis funcoes de recorte de poligono
// eram identicas em cinco exemplos e passaram a viver num lugar so'.
#include "../examples-common/vof-geometry-2d.c"

// MODO GOTA ESTATICA (VOF_ESTATICO=1), para a comparacao com o front-tracking.
//
// O contorno id=1 deste exemplo IMPOE velocidade que cresce no tempo
// (8*(1+tanh(8t-4))*x^2*(1-x)^2): e' escoamento forcado, nao repouso.  Com ele
// ligado a velocidade medida NAO e' corrente parasita -- e' parasita mais
// escoamento imposto, e comparar isso com o front-tracking estatico mediria
// coisas diferentes.  Este modo zera o forcamento para que o caso seja a gota
// estatica de Laplace de verdade.
//
// ATRAS DE UM FLAG de proposito: o padrao preserva o comportamento do exemplo, e
// portanto a referencia da suite.
static int _vof_estatico(void) {
	static int lido = -1;
	if (lido < 0) {
		const char *s = getenv("VOF_ESTATICO");
		lido = (s != NULL) ? atoi(s) : 0;
	}
	return lido;
}

// Deformacao da gota pelos MOMENTOS de area do campo de fracao volumetrica:
//   A = int f dV,  Mxx = int f (x-xc)^2 dV,  Myy = int f (y-yc)^2 dV
//   a = 2 sqrt(Mxx/A),  b = 2 sqrt(Myy/A),  D = (a-b)/(a+b)
// E' EXATAMENTE a grandeza que ft_semieixos mede do lado do front-tracking --
// mesma definicao, nao duas parecidas --, o que e' a condicao para comparar a
// oscilacao dos dois metodos.  Imprime com VOF_DIAG_D a cada VOF_DIAG_CADA passos.
static void _vof_deformacao(higflow_solver *ns, int passo) {
	sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
	mp_mapper *mp = sd_get_domain_mapper(sdm);
	real A = 0.0, sx = 0.0, sy = 0.0;
	higcit_celliterator *it;
	for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int clid = mp_lookup(mp, hig_get_cid(c));
		if (clid < 0) continue;
		Point delta, cc;
		hig_get_delta(c, delta); hig_get_center(c, cc);
		real w = dp_get_value(ns->ed.mult.dpfracvol, clid) * delta[0] * delta[1];
		A += w; sx += w * cc[0]; sy += w * cc[1];
	}
	higcit_destroy(it);
	if (A <= 0.0) return;
	real xc = sx / A, yc = sy / A, Mxx = 0.0, Myy = 0.0;
	for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int clid = mp_lookup(mp, hig_get_cid(c));
		if (clid < 0) continue;
		Point delta, cc;
		hig_get_delta(c, delta); hig_get_center(c, cc);
		real w = dp_get_value(ns->ed.mult.dpfracvol, clid) * delta[0] * delta[1];
		// Momento proprio da celula (+delta^2/12) para nao subestimar a largura
		// nas celulas cortadas -- a celula tem extensao, nao e' um ponto.
		Mxx += w * ((cc[0]-xc)*(cc[0]-xc) + delta[0]*delta[0]/12.0);
		Myy += w * ((cc[1]-yc)*(cc[1]-yc) + delta[1]*delta[1]/12.0);
	}
	higcit_destroy(it);
	real a = 2.0*sqrt(fabs(Mxx/A)), b = 2.0*sqrt(fabs(Myy/A));
	print0f("VOF passo %5d: A=%.8f  a=%.6f  b=%.6f  D=%+.6e\n",
	        passo, (double) A, (double) a, (double) b,
	        (double)((a-b)/(a+b)));
}


// ---------------------------------------------------------------------------
// REFINO ADAPTATIVO (VOF), pelo MESMO criterio e pela MESMA maquina do
// front-tracking: celulas com 0<fracvol<1 sao SEMENTES, a distancia define o
// nivel, os limiares saem da REGRA DAS CELULAS MINIMAS, e o
// higflow_reconstroi_dominio rele' o .amr pelo caminho normal de arranque.
//
// A diferenca para o front-tracking e' de onde vem a fracao: la' da geometria da
// frente, aqui do proprio campo advectado -- que e' justamente o que a
// comparacao entre os dois isola.
//
// AQUI A FRACAO E' ESTADO e viaja no remalhamento (o reconstroi a transporta).
// ---------------------------------------------------------------------------
static void vof_escreve_amr_criterio(higflow_solver *ns, real lx, real ly,
                                     int nx, int ny, int niveis, int cel_min,
                                     const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny;
    const real h_base = (hx < hy) ? hx : hy;
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    const real margem = h_base + h_fino;          // ver malha-adaptativa.c
    const real banda_fina = cel_min * h_fino + margem;
    real thr[8];
    for (int l = niveis - 1; l >= 0; l--) {
        thr[l] = banda_fina;
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        for (int k = 0; k < niveis - 1 - l; k++) thr[l] += cel_min * h_nivel;
    }

    // Sementes: celulas de interface do CAMPO de fracao.
    sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
    mp_mapper *mp = sd_get_domain_mapper(sdm);
    const long ncel = (long) nx * ny;
    Point *sem = (Point *) malloc((size_t) ncel * 4 * sizeof(Point));
    long nsem = 0;
    higcit_celliterator *it;
    for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it);
         higcit_nextcell(it)) {
        hig_cell *c = higcit_getcell(it);
        int clid = mp_lookup(mp, hig_get_cid(c));
        if (clid < 0) continue;
        real f;
        if (ns->par.step == 0) {
            Point xc, d; hig_get_center(c, xc); hig_get_delta(c, d);
            f = ns->ed.mult.problem->fracvol(xc, d, ns->par.t);
        } else {
            f = dp_get_value(ns->ed.mult.dpfracvol, clid);
        }
        if (f > 0.001 && f < 0.999 && nsem < ncel * 4) {
            hig_get_center(c, sem[nsem]); nsem++;
        }
    }
    higcit_destroy(it);

    signed char *tab = (signed char *) calloc((size_t) ncel, 1);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            real cx = (i + 0.5) * hx, cy = (j + 0.5) * hy, d2 = 1e300;
            for (long k = 0; k < nsem; k++) {
                real a = cx - sem[k][0], b = cy - sem[k][1];
                real q = a*a + b*b; if (q < d2) d2 = q;
            }
            int nivel = 0;
            for (int l = niveis - 1; l >= 0; l--)
                if (d2 <= thr[l]*thr[l]) { nivel = l + 1; break; }
            tab[i + (long) j * nx] = (signed char) nivel;
        }

    FILE *f = fopen(caminho, "w");
    if (f == NULL) { perror(caminho); free(tab); free(sem); return; }
    fprintf(f, "0.0 %.10g 0.0 %.10g\n", (double) lx, (double) ly);
    long cont[8] = {0};
    for (long q = 0; q < ncel; q++)
        for (int l = 1; l <= niveis; l++) if (tab[q] >= l) cont[l]++;
    int niv_esc = 1;
    for (int l = 1; l <= niveis; l++) if (cont[l] > 0) niv_esc = l + 1;
    fprintf(f, "%d\n", niv_esc);
    // indice base 1, como os .amr do repositorio
    fprintf(f, "%.10g %.10g\n1\n1 1 %d %d\n", (double) hx, (double) hy, nx, ny);
    for (int l = 1; l < niv_esc; l++) {
        int fator = 1 << l;
        fprintf(f, "%.10g %.10g\n%ld\n", (double)(hx/fator), (double)(hy/fator), cont[l]);
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                if (tab[i + (long) j * nx] >= l)
                    fprintf(f, "%d %d %d %d\n", fator*i+1, fator*j+1, fator, fator);
    }
    fclose(f);
    fprintf(stderr, "VOF amr: %ld sementes, niveis=%d, limiares=", nsem, niv_esc-1);
    for (int l = 0; l < niveis; l++) fprintf(stderr, " %.5f", (double) thr[l]);
    fprintf(stderr, "  (>=%d celulas finas por lado)\n", cel_min);
    free(tab); free(sem);
}


// Criterio de ARRANQUE: sem dominio ainda, a fracao sai da forma analitica
// (func(), a mesma que inicializa o fracvol).  Mesma regra e mesmos limiares.
static void vof_escreve_amr_inicial(real lx, real ly, int nx, int ny,
                                    int niveis, int cel_min, const char *caminho)
{
    const real hx = lx / nx, hy = ly / ny;
    const real h_base = (hx < hy) ? hx : hy;
    real h_fino = h_base;
    for (int l = 0; l < niveis; l++) h_fino *= 0.5;
    const real margem = h_base + h_fino;
    const real banda_fina = cel_min * h_fino + margem;
    real thr[8];
    for (int l = niveis - 1; l >= 0; l--) {
        thr[l] = banda_fina;
        real h_nivel = h_base;
        for (int k = 0; k < l + 1; k++) h_nivel *= 0.5;
        for (int k = 0; k < niveis - 1 - l; k++) thr[l] += cel_min * h_nivel;
    }
    const long ncel = (long) nx * ny;
    Point *sem = (Point *) malloc((size_t) ncel * sizeof(Point));
    long nsem = 0;
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            // celula cortada = os quatro cantos nao concordam em sinal de func()
            int dentro = 0;
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++) {
                    Point q; q[0] = (i + a) * hx; q[1] = (j + b) * hy;
                    if (func(q) > 0.0) dentro++;
                }
            if (dentro > 0 && dentro < 4) {
                sem[nsem][0] = (i + 0.5) * hx; sem[nsem][1] = (j + 0.5) * hy;
                nsem++;
            }
        }
    signed char *tab = (signed char *) calloc((size_t) ncel, 1);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            real cx = (i + 0.5) * hx, cy = (j + 0.5) * hy, d2 = 1e300;
            for (long k = 0; k < nsem; k++) {
                real a = cx - sem[k][0], b = cy - sem[k][1];
                real q = a*a + b*b; if (q < d2) d2 = q;
            }
            int nivel = 0;
            for (int l = niveis - 1; l >= 0; l--)
                if (d2 <= thr[l]*thr[l]) { nivel = l + 1; break; }
            tab[i + (long) j * nx] = (signed char) nivel;
        }
    FILE *f = fopen(caminho, "w");
    if (f == NULL) { perror(caminho); free(tab); free(sem); return; }
    fprintf(f, "0.0 %.10g 0.0 %.10g\n", (double) lx, (double) ly);
    long cont[8] = {0};
    for (long q = 0; q < ncel; q++)
        for (int l = 1; l <= niveis; l++) if (tab[q] >= l) cont[l]++;
    int niv_esc = 1;
    for (int l = 1; l <= niveis; l++) if (cont[l] > 0) niv_esc = l + 1;
    fprintf(f, "%d\n", niv_esc);
    fprintf(f, "%.10g %.10g\n1\n1 1 %d %d\n", (double) hx, (double) hy, nx, ny);
    for (int l = 1; l < niv_esc; l++) {
        int fator = 1 << l;
        fprintf(f, "%.10g %.10g\n%ld\n", (double)(hx/fator), (double)(hy/fator), cont[l]);
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                if (tab[i + (long) j * nx] >= l)
                    fprintf(f, "%d %d %d %d\n", fator*i+1, fator*j+1, fator, fator);
    }
    fclose(f);
    fprintf(stderr, "VOF amr inicial: %ld sementes, niveis=%d\n", nsem, niv_esc-1);
    free(tab); free(sem);
}

// Volume fraction

// ---------------------------------------------------------------------------
// O problema deste exemplo, como um tipo em vez de oito funcoes soltas.
// Os corpos sao os mesmos; so' mudaram de lugar e perderam o prefixo get_.
// ---------------------------------------------------------------------------
class VofProblem : public HigFlowProblem, public HigFlowMultiphaseProblem {
public:
    // Value of the pressure
    real pressure(Point center, real t) {
        real value = 0.0;
    	return value;
    }
    // Value of the velocity
    real velocity(Point center, int dim, real t) {
    	real value;
    	real x = center[0];
    	real y = center[1];

    	// switch(dim){
    	// 	case 0: ;
    	// 		value = 0.5 - y;
    	// 		// real sx = sin(pi*x);
    	// 		// real cy = cos(pi*y);
    	// 		// value = sx*cy;
    	// 	break;
    	// 	case 1: ;
    	// 		value = x - 0.5;
    	// 		// real sy = sin(pi*y);
    	// 		// real cx = cos(pi*x);
    	// 		// value = -sy*cx;
    	// 	break;
    	// }

    	value = 0.0;
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
    	real x = center[0];
    	real y = center[1];

    	switch (dim) {
    		case 0: ;
    			// value = 0.5 - y;
    			// // real sx = sin(pi*x);
    			// // real cy = cos(pi*y);
    			// // if(t<pi) value = sx*cy;
    			// // else value = -sx*cy;
    			switch (id) {
    				case 0: ;
    					value = 0.0;
    					break;
    				case 1: ;
    					// Forcamento da tampa; zerado no modo gota estatica.
    					value = _vof_estatico() ? 0.0
    					      : 8.0*(1.0 + tanh(8.0*t - 4.0))*x*x*(1.0 - x)*(1.0 - x);
    					break;
    				case 2: ;
    					value = 0.0;
    					break;
    				case 3: ;
    					value = 0.0;
    					break;
    			}
    			break;
    		case 1: ;
    			// value = x - 0.5;
    			// // real sy = sin(pi*y);
    			// // real cx = cos(pi*x);
    			// // if(t<pi) value = -sy*cx;
    			// // else value = sy*cx;
    			value = 0.0;
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

    // --- modelo multifasico ---
    // Value of the viscosity
    real viscosity0(Point center, real t) {
    	// Propriedade por ambiente (VOF_MU0), para o VOF rodar o MESMO caso de
    	// Hysing que o front-tracking.  Padrao 1.0 -- o comportamento original.
    	const char *e = getenv("VOF_MU0");
    	return (e != NULL) ? atof(e) : 1.0;
    }
    // Value of the viscosity
    real viscosity1(Point center, real t) {
    	// Propriedade por ambiente (VOF_MU1), para o VOF rodar o MESMO caso de
    	// Hysing que o front-tracking.  Padrao 1.0 -- o comportamento original.
    	const char *e = getenv("VOF_MU1");
    	return (e != NULL) ? atof(e) : 1.0;
    }
    // Value of the density
    real density0(Point center, real t) {
    	// Propriedade por ambiente (VOF_RHO0), para o VOF rodar o MESMO caso de
    	// Hysing que o front-tracking.  Padrao 1.0 -- o comportamento original.
    	const char *e = getenv("VOF_RHO0");
    	return (e != NULL) ? atof(e) : 1.0;
    }
    // Value of the density
    real density1(Point center, real t) {
    	// Propriedade por ambiente (VOF_RHO1), para o VOF rodar o MESMO caso de
    	// Hysing que o front-tracking.  Padrao 1.0 -- o comportamento original.
    	const char *e = getenv("VOF_RHO1");
    	return (e != NULL) ? atof(e) : 1.0;
    }
    real fracvol(Point center, Point delta, real t) {
    	Point p0, p1, p2, p3;
    	real  f0, f1, f2, f3;
    	real  value;
    	int var = 0;

    	// Canto inferior esquerdo
    	p0[0] = center[0] - 0.5*delta[0];
    	p0[1] = center[1] - 0.5*delta[1];
    	f0    = func(p0);
    	if (f0 > 0.0) var += 1;

    	// Canto inferior direito
    	p1[0] = center[0] + 0.5*delta[0];
    	p1[1] = center[1] - 0.5*delta[1];
    	f1    = func(p1);
    	if (f1 > 0.0) var += 1;

    	// Canto superior esquerdo
    	p2[0] = center[0] - 0.5*delta[0];
    	p2[1] = center[1] + 0.5*delta[1];
    	f2    = func(p2);
    	if (f2 > 0.0) var += 1;

    	// Canto superior direito
    	p3[0] = center[0] + 0.5*delta[0];
    	p3[1] = center[1] + 0.5*delta[1];
    	f3    = func(p3);
    	if (f3 > 0.0) var += 1;

    	if (var == 0){
    		value = 0.0;
    	} else if (var == 4){
    		value = delta[0]*delta[1];
    	}
    	else{
    		value =0.0;
    		int N=16;
    		real delta_new[2];
    		delta_new[0]=delta[0]/N;delta_new[1]=delta[1]/N;
    		real center_new[2];
    		for(int i=0;i<N;i++)
    		{	
    			for (int j=0;j<N;j++)
    			{
    				center_new[0]=p0[0]+(i+0.5)*delta_new[0];
    				center_new[1]=p0[1]+(j+0.5)*delta_new[1];
    				value = value + get_fracvolN(center_new,delta_new, t);
    			}
    		}
    	}	

    	value = value/delta[0]/delta[1];
    	return value;
    }
};

static VofProblem problema;

// Massa da gota: integral da fracao volumetrica sobre o dominio.  Mesma conta do
// example2d_DynamicMeshAdapt (compute_total_fracvol), trazida para ca' porque nao
// esta' na biblioteca.  Serve a comparacao com o front-tracking, onde o analogo e'
// a area fechada pela frente (ft_area, formula do laco).
real compute_total_fracvol_local(higflow_solver *ns) {
	sim_domain *sdm = psd_get_local_domain(ns->ed.mult.psdmult);
	mp_mapper *mp = sd_get_domain_mapper(sdm);
	real total = 0.0;
	higcit_celliterator *it;
	for (it = sd_get_domain_celliterator(sdm); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int clid = mp_lookup(mp, hig_get_cid(c));
		if (clid < 0) continue;
		Point delta;
		hig_get_delta(c, delta);
		real volcell = delta[0] * delta[1];
		real fracvol = dp_get_value(ns->ed.mult.dpfracvol, clid);
		total += fracvol * volcell;
	}
	higcit_destroy(it);
	real global_total;
	MPI_Allreduce(&total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return global_total;
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
	// MALHA INICIAL JA' ADAPTADA, e aqui ela e' obrigatoria, nao otimizacao.
	//
	// O sdmult do HiGFlow e' sempre UNIFORME no nivel mais fino.  Comecando na
	// malha base e remalhando para 2 niveis, ele salta de h para h/4 DE UMA VEZ
	// -- e a transferencia por posicao so' promete salto de UM nivel.  MEDIDO:
	// a massa caiu a 1/4 exato (0,19635 -> 0,04909) no primeiro remalhamento.
	// Nascendo ja' no nivel final, o sdmult nao muda mais e a fracao viaja
	// sempre em salto de um nivel -- ou em nenhum.
	{
		const char *e; int niv = 0;
		if ((e = getenv("VOF_AMR_NIVEIS")) != NULL) niv = atoi(e);
		if (niv > 0) {
			int nx0 = 40, ny0 = 80, cel0 = 5;
			real lx0 = 1.0, ly0 = 2.0;
			if ((e = getenv("VOF_AMR_NX")) != NULL) nx0 = atoi(e);
			if ((e = getenv("VOF_AMR_NY")) != NULL) ny0 = atoi(e);
			if ((e = getenv("VOF_AMR_CELMIN")) != NULL) cel0 = atoi(e);
			if ((e = getenv("VOF_AMR_LX")) != NULL) lx0 = atof(e);
			if ((e = getenv("VOF_AMR_LY")) != NULL) ly0 = atof(e);
			const char *cam0 = getenv("VOF_AMR_CAMINHO");
			if (cam0 == NULL) cam0 = "amrs-hysing/criterio/dominio.amr";
			vof_escreve_amr_inicial(lx0, ly0, nx0, ny0, niv, cel0, cam0);
			print0f("===> malha inicial JA' adaptada, do criterio analitico\n");
		}
	}

	higflow_create_domain(ns, cache, order_center);
	// Create the simulation domain for non newtonian simulation
	// higflow_create_domain_generalized_newtonian(ns, cache, order_center, get_viscosity);
	higflow_create_domain_multiphase(ns, cache, order_center, &problema);
	
	// Initialize the domain
    print0f("=+=+=+= Load Domain =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
   //higflow_initialize_domain(ns, ntasks, myrank, order_facet); 
#ifdef HIGFLOW_COM_T8CODE
    malha_t8_instala(ns, myrank);
#endif
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet); 
    
	// Initialize the boundaries
    print0f("=+=+=+= Load Bondary Condtions =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
	//higflow_initialize_boundaries(ns);
	higflow_initialize_boundaries_yaml(ns);
   
	// Creating distributed property
	higflow_create_distributed_properties(ns);
	// Initialize distributed properties
    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);

	// Volume INICIAL, para separar DERIVA de erro de INICIALIZACAO.  Sem isto a
	// comparacao com o front-tracking mediria coisas diferentes: o desvio contra
	// a area analitica mistura a discretizacao inicial do circulo com a perda de
	// massa ao longo do tempo, e e' a segunda que caracteriza o metodo.
	real vol_inicial = -1.0;
	if (getenv("VOF_DIAG_LAPLACE") != NULL) {
		vol_inicial = compute_total_fracvol_local(ns);
		print0f("=+=+=+= VOF MASSA INICIAL  vol0=%.8f =+=+=+=\n",
		        (double) vol_inicial);
	}
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

    // Parametros do refino adaptativo, por ambiente.
    int  vof_amr_niveis = 0, vof_amr_celmin = 5, vof_amr_cada = 0;
    int  vof_amr_nx = 40, vof_amr_ny = 80;
    real vof_amr_lx = 1.0, vof_amr_ly = 2.0;
    const char *vof_amr_caminho = "amrs-hysing/criterio/dominio.amr";
    {
        const char *e;
        if ((e = getenv("VOF_AMR_NIVEIS")) != NULL) vof_amr_niveis = atoi(e);
        if ((e = getenv("VOF_AMR_CELMIN")) != NULL) vof_amr_celmin = atoi(e);
        if ((e = getenv("VOF_AMR_CADA"))   != NULL) vof_amr_cada   = atoi(e);
        if ((e = getenv("VOF_AMR_NX"))     != NULL) vof_amr_nx     = atoi(e);
        if ((e = getenv("VOF_AMR_NY"))     != NULL) vof_amr_ny     = atoi(e);
        if ((e = getenv("VOF_AMR_LX"))     != NULL) vof_amr_lx     = atof(e);
        if ((e = getenv("VOF_AMR_LY"))     != NULL) vof_amr_ly     = atof(e);
        if ((e = getenv("VOF_AMR_CAMINHO"))!= NULL) vof_amr_caminho = e;
    }

    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
		// Print the step
        print0f("===> Step:        %7d <====> t  = %15.10lf <===\n", ns->par.step, ns->par.t);
		// Start the first step time
        if (ns->par.step == step0)  START_CLOCK(firstiter); 
		// Update velocities and pressure using the projection method
		//printf("===> calling step multiphase at step: %d <=== \n",ns->par.stepaux);
		// real max_IF0 = max_dp(ns, ns->ed.mult.dpIF[0]);
		// real max_IF1 = max_dp(ns, ns->ed.mult.dpIF[1]);
		// printf("===> max_IF0 = %f <=== \n",max_IF0);
		// printf("===> max_IF1 = %f <=== \n",max_IF1);
		// REMALHAMENTO ADAPTATIVO (VOF_AMR_NIVEIS>0), mesma maquina do
		// front-tracking.  Aqui a fracao E' estado: o reconstroi a transporta
		// junto com u e p, com interpolacao conservativa para salto de um nivel.
		if (vof_amr_niveis > 0 && vof_amr_cada > 0 && ns->par.step > 0 &&
		    ns->par.step % vof_amr_cada == 0) {
			vof_escreve_amr_criterio(ns, vof_amr_lx, vof_amr_ly,
			                         vof_amr_nx, vof_amr_ny,
			                         vof_amr_niveis, vof_amr_celmin,
			                         vof_amr_caminho);
			long faltam = higflow_reconstroi_dominio(ns, ntasks, myrank,
			                                         cache, order_center,
			                                         order_facet);
			if (faltam != 0)
				print0f("===> remalha passo %d: %ld posicoes sem valor\n",
				        ns->par.step, faltam);
			if (getenv("VOF_DIAG_MALHA") != NULL)
				print0f("===> MALHA passo %d: banda=%.2f celulas finas "
				        "(exigido >=%d)\n", ns->par.step,
				        (double) malha_adapt_mede_banda(ns, vof_amr_niveis,
				                                        vof_amr_lx/vof_amr_nx),
				        vof_amr_celmin);
		}

		higflow_solver_step_multiphase(ns);

		// GRANDEZAS DO BENCHMARK DE HYSING (VOF_HYSING=1).  yc e Vc pela MESMA
		// formula do front-tracking.  O PERIMETRO nao sai de graca do campo: a
		// reconstrucao PLIC do solver nao esta' exposta, entao a circularidade
		// do lado VOF NAO e' reportada -- e' melhor faltar do que sair por uma
		// aproximacao diferente da do outro lado, que tornaria a comparacao
		// entre os dois uma comparacao de definicoes.
		if (getenv("VOF_HYSING") != NULL) {
			static int cada_h = 0;
			if (cada_h == 0) { const char *c = getenv("VOF_HYSING_CADA");
			                   cada_h = (c != NULL) ? atoi(c) : 50;
			                   if (cada_h < 1) cada_h = 50; }
			if (ns->par.step % cada_h == 0) {
				real Ah, ych, vch;
				hysing_medidas(ns, &Ah, &ych, &vch);
				print0f("HYSING %10.5f %12.6f %12.6f %12.8f\n",
				        (double) ns->par.t, (double) ych, (double) vch,
				        (double) Ah);
			}
		}
		ns->par.stepaux=ns->par.stepaux+1;
		
		// Time update
		ns->par.t += ns->par.dt;
		// Stop the first step time
        if (ns->par.step == step0) STOP_CLOCK(firstiter); 
        { static int cada = 0;
          if (cada == 0) { const char *c = getenv("VOF_DIAG_CADA");
                           cada = (c != NULL) ? atoi(c) : 50; if (cada < 1) cada = 50; }
          if (getenv("VOF_DIAG_D") != NULL && ns->par.step % cada == 0)
              _vof_deformacao(ns, ns->par.step); }

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

	// ------------------------------------------------------------------
	// SONDAS DE COMPARACAO COM O FRONT-TRACKING (mesmas metricas, mesma
	// forma de medir).  Ligadas por VOF_DIAG_LAPLACE.
	//
	// Este caso E' a gota estatica de Laplace: gota circular R=1/6 em
	// (0,5;0,5), newtoniana bifasica com PROPRIEDADES IGUAIS (rho=mu=1 nas
	// duas fases), u=0 inicial, tensao superficial ligada.  A forca entra
	// como IF/(We*rho) com We=Ca*Re=1, entao sigma_efetivo = 1 e o salto
	// esperado e' sigma/R = 6,0 -- exatamente a mesma fisica que o
	// front-tracking resolve em example2d_FrontTracking, so' mudando a
	// representacao da interface.
	//
	// As correntes parasitas saem do proprio print do solver (Vmin/Vmax),
	// que e' o MESMO codigo nos dois exemplos -- comparaveis sem ajuste.
	// ------------------------------------------------------------------
	if (getenv("VOF_DIAG_LAPLACE") != NULL) {
		// Massa: integral da fracao volumetrica.  Conservacao por construcao
		// e' a vantagem classica do VOF -- aqui ela e' MEDIDA, nao suposta.
		real vol = compute_total_fracvol_local(ns);
		// Raio de EQUILIBRIO: a gota relaxa para o circulo de mesma area, entao
		// R_eq = sqrt(ra*rb) -- que para ra=rb devolve o proprio raio.
		real ra = 1.0/6.0, rb = 1.0/6.0;
		{ const char *e;
		  if ((e = getenv("VOF_RA")) != NULL) ra = atof(e);
		  if ((e = getenv("VOF_RB")) != NULL) rb = atof(e); }
		const real R = sqrt(ra * rb);
		const real vol_exato = M_PI * ra * rb;

		// Salto de pressao: dentro (centro da gota) e fora (longe dela).
		Point p_in, p_out;
		p_in[0]  = 0.5; p_in[1]  = 0.5;
		p_out[0] = 0.5; p_out[1] = 0.9;
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

		const real We = ns->ed.mult.spar.Ca * ns->par.Re;   // = 1
		const real sigma = 1.0 / We;
		const real dp_exato = sigma / R;
		print0f("=+=+=+= VOF LAPLACE  p_in=%.6f  p_out=%.6f  Dp=%.6f  "
		        "sigma/R=%.6f  erro_rel=%.4f =+=+=+=\n",
		        (double) pin, (double) pout, (double)(pin - pout),
		        (double) dp_exato,
		        (double) fabs((pin - pout) - dp_exato) / dp_exato);
		// DUAS medidas distintas, e a distincao importa:
		//   deriva      |vol(t) - vol(0)|/vol(0) -- perda de massa NO TEMPO, que
		//               e' o que caracteriza o metodo.  Comparavel com o dA/A do
		//               front-tracking.
		//   inicializacao |vol(0) - pi R^2|/(pi R^2) -- quao bem a discretizacao
		//               representa o circulo no instante zero.  No front-tracking
		//               o analogo e' o deficit de area do poligono de N lados.
		print0f("=+=+=+= VOF MASSA  vol=%.8f  vol0=%.8f  deriva_rel=%.3e  "
		        "exato=%.8f  inicializacao_rel=%.3e =+=+=+=\n",
		        (double) vol, (double) vol_inicial,
		        (vol_inicial > 0.0)
		            ? (double) (fabs(vol - vol_inicial) / vol_inicial) : -1.0,
		        (double) vol_exato,
		        (vol_inicial > 0.0)
		            ? (double) (fabs(vol_inicial - vol_exato) / vol_exato) : -1.0);
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
