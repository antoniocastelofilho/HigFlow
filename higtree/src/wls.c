#include "utils.h"
#include "coord.h"
#include "rng.h"
#include "string.h"
#include "Debug-c.h"
#include "wls.h"

#define USE_MLS

/****************************************************************************/
/* Esta rotina aplica a matriz de rotacao de Givens nas matrizes A e Q para */
/* a decomposicao A = QR.                                                   */
/****************************************************************************/

void wls_givens ( int m, int n, int i, int k, mreal A, mreal Q, real s1, real s2 ) {
	int   j;   /* variavel auxiliar para controle de lacos  */
	real  s,t; /* variaveis auxiliares para armazenar dados */

	if (fabs(s2)+fabs(s1) > 0.0) {
		if (fabs(s2) >= fabs(s1))
			s  = sqrt(1.0+(s1/s2)*(s1/s2))*fabs(s2);
		else
			s  = sqrt(1.0+(s2/s1)*(s2/s1))*fabs(s1);
		s1 = s1/s; s2 = s2/s;

		/* aplicando a matriz de rotacao em A */
		for (j = 0; j < n; j++) {
			s       =  s1*A[i][j]+s2*A[k][j];
			t       = -s2*A[i][j]+s1*A[k][j];
			A[i][j] =  s;
			A[k][j] =  t;
		}

		/* aplicando a matriz de rotacao em Q */
		for (j = 0; j < m; j++) {
			s       =  s1*Q[i][j]+s2*Q[k][j];
			t       = -s2*Q[i][j]+s1*Q[k][j];
			Q[i][j] =  s;
			Q[k][j] =  t;
		}
	}

	return;
}



/**************************************************************************/
/* Esta rotina retorna a decomposicao QR de uma matriz A (n+1)x(n) usando */
/* o metodo de Givens. A matriz R volta em A e em Q volta a matriz Q^t da */
/* decomposicao A = QR. Esta rotina tambem calcula uma estimativa para a  */
/* condicao da matriz A.                                                  */
/**************************************************************************/

void wls_qr ( int m, int n, mreal A, mreal Q) {
	int   i,j; /* variaveis auxiliares para controle de lacos       */

	/* inicializando Q com a identidade */
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) Q[i][j] = 0.0;
		Q[i][i] = 1.0;
	}

	/* decomposicao QR por Givens */
	for (i = 0; i < n; i++) {
		for (j = i+1; j < m; j++) {
			/* escolha da matriz de rotacao */
			wls_givens(m,n,i,j,A,Q,A[i][i],A[j][i]);
		}
	}

	return;
}


/**************************************************************************/
/* Esta rotina retorna a decomposicao QR de uma matriz A (n+1)x(n) usando */
/* o metodo de Givens. A matriz R volta em A e em Q volta a matriz Q^t da */
/* decomposicao A = QR. Esta rotina tambem calcula uma estimativa para a  */
/* condicao da matriz A.                                                  */
/**************************************************************************/

void wls_qr_update ( int m, int k, int n, mreal A, mreal Q) {
	int   i,j; /* variaveis auxiliares para controle de lacos       */

	/* decomposicao QR por Givens */
	for (i = 0; i < n; i++) {
		for (j = m; j < m+k; j++) {
			/* escolha da matriz de rotacao */
			wls_givens(m+k,n,i,j,A,Q,A[i][i],A[j][i]);
		}
	}

	return;
}



/**************************************************************************/
/* Esta obtem a solucao de Ax = b, onde A (m)x(n) usando o metodo dos     */
/* minimos quadrados, sabendo a decomposicao QR de A armazenados em Q e A */
/**************************************************************************/

void wls_solve_system ( int m, int n, mreal A, mreal Q, vreal x, vreal y, real eps ) {
	int   i,j;

	for (i = 0; i < n; i++) {
                if (fabs(A[i][i]) > eps) {
			for (j = 0, x[i] = y[i]; j < i; j++) x[i] -= A[j][i]*x[j];
			x[i] = x[i]/A[i][i];
		} else {
			x[i] = 0.0;
		}
	}

	for (i = 0; i < n; i++) y[i] = x[i];
	for (i = n; i < m; i++) y[i] = 0.0;

	for (i = 0; i < m; i++)
		for (j = 0, x[i] = 0.0; j < m; j++)
			x[i] += Q[j][i]*y[j];

	return;
}


/********************************************/
/* Funcao para alocar um vetor de reais     */
/********************************************/
/* PARAMETROS DE ENTRADA: dimensao do vetor */
/* VALOR RETORNADO: ponteiro para o vetor   */
/* FUNCOES ATIVADAS: malloc()               */
/********************************************/

vreal wls_allocvreal ( int n ) {
	int  i;   /* variavel auxiliar     */

	/* alocando o vetor */
	DECL_AND_ALLOC(real, v, n+1);
	for (i = 0; i <= n; i++) v[i] = 0.0;

	/* retornando o ponteiro */
	return (v);
}


/*********************************************************/
/* Funcao para liberar um vetor de reais                 */
/*********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para vetor e dimensao */
/* VALOR RETORNADO: ponteiro para o vetor                */
/* FUNCOES ATIVADAS: free()                             */
/*********************************************************/

vreal wls_freevreal ( int n, vreal v ) {
	/* liberando o vetor */
	free(v);

	/* retornando o ponteiro */
	return (NULL);
}


/**********************************************/
/* Funcao para alocar uma matriz de reais     */
/**********************************************/
/* PARAMETROS DE ENTRADA: dimensoes da matriz */
/* VALOR RETORNADO: ponteiro para a matriz    */
/* FUNCOES ATIVADAS: free()                  */
/**********************************************/

mreal wls_allocmreal ( int m, int n ) {
	int   i,j; /* variaveis auxiliares   */

	/* alocando as linhas da matriz */
	DECL_AND_ALLOC(vreal, v, m+1);

	/* alocando as colunas da matriz */
	for (i = 0; i <= m; i++) {
  		ALLOC(real, v[i], n+1);
  		for (j = 0; j <= n; j++) v[i][j] = 0.0;
	}

	/* retornando o ponteiro */
	return (v);
}


/***********************************************************/
/* Funcao para liberar uma matriz de reais                 */
/***********************************************************/
/* PARAMETROS DE ENTRADA: ponteiro para matriz e dimensoes */
/* VALOR RETORNADO: ponteiro para matriz                   */
/* FUNCOES ATIVADAS: free()                               */
/***********************************************************/

mreal wls_freemreal ( int m, int n, mreal v ) {
	int      i;  /* variavel auxiliar */

	/* liberando as linhas da matriz */
	for (i = 0; i <= m; i++) free(v[i]);

	/* libreando a matriz */
	free(v);

	/* retornando o ponteiro */
	return (NULL);
}


// Combination
static int combinations_with_repetition ( int n, int k ) {
	int i, c;

	for (i = n, c = 1; i > n-k; i--) c *= i;
	for (i = 1; i <= k; i++)  c /= i;

	return c;
}


// Minimal number of points for dimension n and order ord
int wls_num_min_points (unsigned dim, int ord )
{
	int i, npts_min;

	npts_min = 1;
	for (i = 1; i <= ord; i++) npts_min +=
		combinations_with_repetition(dim+i-1,i);

	return npts_min;
}


// Assembly matrix line for a dimension dim point x and order ord
static void _assembly_pts_matrix (unsigned dim, int ord, int m, vreal a, const Point x ) {
	int  k, i1, i2, i3, i4;

	k    = 0;

	if (ord >= 0) {
		a[k] = 1;
		k++;
		if (k > m) return;
	}

	if (ord >= 1) {
		for (i1 = 0; i1 < dim; i1++) {
			a[k] = x[i1];
			k++;
			if (k > m) return;
		}
	}

	if (ord >= 2) {
		for (i1 = 0; i1 < dim; i1++) {
			for (i2 = i1; i2 < dim; i2++) {
				a[k] = x[i1]*x[i2];
				k++;
				if (k > m) return;
			}
		}
	}

	if (ord >= 3) {
		for (i1 = 0; i1 < dim; i1++) {
			for (i2 = i1; i2 < dim; i2++) {
				for (i3 = i2; i3 < dim; i3++) {
					a[k] = x[i1]*x[i2]*x[i3];
					k++;
					if (k > m) return;
				}
			}
		}
	}

	if (ord >= 4) {
		for (i1 = 0; i1 < dim; i1++) {
			for (i2 = i1; i2 < dim; i2++) {
				for (i3 = i2; i3 < dim; i3++) {
					for (i4 = i3; i4 < dim; i4++) {
						a[k] = x[i1]*x[i2]*x[i3]*x[i4];
						k++;
						if (k > m) return;
					}
				}
			}
		}
	}

	return;
}

static bool inc_vector(unsigned dim, int degree, int i[degree])
{
	for(int j = 0; j < degree; ++j) {
		if(++i[j] != dim) {
			int val = i[j];
			while(--j >= 0) {
				i[j] = val;
			}
			return true;
		}
	}
	return false;
}

static real eval_poly(CPPoint x, unsigned dim, int degree, real* coefs)
{
	assert(degree >= 0);

	real ret = 0.0;

	if(degree > 0) {
		int i[degree];
		memset(i, 0, sizeof i);
		do {
			real accum = *coefs++;
			for(unsigned j = 0; j < degree; ++j) {
				accum *= x[i[j]];
			}
			ret += accum;
		} while(inc_vector(dim, degree, i));

		ret += eval_poly(x, dim, degree - 1, coefs);
	} else {
		ret += *coefs;
	}

	return ret;
}

wls_interpolator * wls_create(unsigned dim, int ord, int maxnumpts)
{
	DECL_AND_ALLOC(wls_interpolator, wls, 1);
	wls->type = WLSMOVING;
	//wls->type = WLSSTATIC;
	wls->numpoly = wls_num_min_points(dim, ord);
	wls->maxnumpts = maxnumpts;
	wls->ord = ord;
	wls->dim = dim;
	int n = wls->numpoly;
	int m = wls->maxnumpts;
	wls->A = (mreal) wls_allocmreal(m,n);
	wls->Q = (mreal) wls_allocmreal(m,m);
	wls->D = (vreal) wls_allocvreal(m);
	wls->b = (vreal) wls_allocvreal(m);
	wls->w = (vreal) wls_allocvreal(m);

	// Random number seed for poly coefs,
	// Merkle root of Bitcoin Block 0:
	const char seed[] = {
		0x4a, 0x5e, 0x1e, 0x4b, 0xaa, 0xb8, 0x9f,
		0x3a, 0x32, 0x51, 0x8a, 0x88, 0xc3, 0x1b,
		0xc8, 0x7f, 0x61, 0x8f, 0x76, 0x67, 0x3e,
		0x2c, 0xc7, 0x7a, 0xb2, 0x12, 0x7b, 0x7a,
		0xfd, 0xed, 0xa3, 0x3b
	};
	RNG *rng = rng_create(seed, sizeof seed);
	wls->poly_coefs = wls_allocvreal(n);
	for(int i = 0; i < n; ++i) {
		wls->poly_coefs[i] = rng_uniform(rng, -1.0, 1.0);
	}
	rng_destroy(rng);

	return wls;
}

void wls_set_type(wls_interpolator *wls, wls_type t) {
	wls->type = t;
}

real wls_set_points(wls_interpolator *wls, int numpts, Point pts[], Point x) {
        for (int i = 0; i < numpts; i++) {
		_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->A[i],pts[i]);
	}
	if (wls->type == WLSMOVING) {
        	for (int i = 0; i < numpts; i++) {
			real mlsw = 0.0;
			for (int j = 0; j < wls->dim; j++) {
				mlsw += (x[j]-pts[i][j])*(x[j]-pts[i][j]);
			}
			wls->D[i] = sqrt(1.0/(mlsw+1.0e-10));
			for (int j = 0; j < wls->numpoly; j++) {
				wls->A[i][j] *= wls->D[i];
			}
		}
	}
	wls_qr(numpts,wls->numpoly,wls->A,wls->Q);
        real mindiag = fabs(wls->A[1][1]);
	for (int i = 2; i < wls->numpoly; i++) {
		real z = fabs(wls->A[i][i]);
		if (z < mindiag) {
			mindiag = z;
		}
	}
	return mindiag;
}

real wls_add_points(wls_interpolator *wls, int numpts, int addnumpts, Point pts[], Point x) {
        for (int i = 0; i < addnumpts; i++) {
		_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->A[numpts+i],pts[i]);
	}
	if (wls->type == WLSMOVING) {
        	for (int i = numpts; i < numpts+addnumpts; i++) {
			real mlsw = 0.0;
			for (int j = 0; j < wls->dim; j++) {
				mlsw += (x[j]-pts[i][j])*(x[j]-pts[i][j]);
			}
			wls->D[i] = sqrt(1.0/(mlsw+1.0e-10));
			for (int j = 0; j < wls->numpoly; j++) {
				wls->A[i][j] *= wls->D[i];
			}
		}
	}
        for (int i = 0; i < numpts; i++) {
		for (int j = numpts; j < numpts+addnumpts; j++) {
			wls->Q[i][j] = 0.0;
		}
	}
        for (int i = numpts; i < numpts+addnumpts; i++) {
		for (int j = 0; j < numpts; j++) {
			wls->Q[i][j] = 0.0;
		}
	}
        for (int i = numpts; i < numpts+addnumpts; i++) {
		for (int j = numpts; j < numpts+addnumpts; j++) {
			wls->Q[i][j] = 0.0;
		}
		wls->Q[i][i] = 1.0;
	}
	wls_qr_update(numpts,addnumpts,wls->numpoly,wls->A,wls->Q);
        real mindiag = fabs(wls->A[1][1]);
	for (int i = 2; i < wls->numpoly; i++) {
		real z = fabs(wls->A[i][i]);
		if (z < mindiag) {
			mindiag = z;
		}
	}
	return mindiag;
}

void wls_calc(wls_interpolator *wls, Point x, int numpts, real w[]) {
	_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->b,x);
	wls_solve_system(numpts,wls->numpoly,wls->A,wls->Q,wls->w,wls->b,EPSMACH);
	if (wls->type == WLSMOVING) {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->D[i]*wls->w[i];
		}
	} else {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->w[i];
		}
	}
}

real wls_set_points_and_calc(wls_interpolator *wls, int numpts, Point pts[], Point x, real w[]) {
        for (int i = 0; i < numpts; i++) {
		_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->A[i],pts[i]);
	}
	if (wls->type == WLSMOVING) {
        	for (int i = 0; i < numpts; i++) {
			real mlsw = 0.0;
			for (int j = 0; j < wls->dim; j++) {
				mlsw += (x[j]-pts[i][j])*(x[j]-pts[i][j]);
			}
			wls->D[i] = sqrt(1.0/(mlsw+1.0e-10));
			for (int j = 0; j < wls->numpoly; j++) {
				wls->A[i][j] *= wls->D[i];
			}
		}
	}
	wls_qr(numpts,wls->numpoly,wls->A,wls->Q);
	_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->b,x);
        real mindiag = fabs(wls->A[1][1]);
	for (int i = 2; i < wls->numpoly; i++) {
		real z = fabs(wls->A[i][i]);
		if (z < mindiag) {
			mindiag = z;
		}
	}
	wls_solve_system(numpts,wls->numpoly,wls->A,wls->Q,wls->w,wls->b,EPSMACH);
	if (wls->type == WLSMOVING) {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->D[i]*wls->w[i];
		}
	} else {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->w[i];
		}
	}
	return mindiag;
}

real wls_set_samples_and_calc(wls_interpolator *wls, int numpts, wls_item items[], const Point x, real w[])
{
	for (int i = 0; i < numpts; i++) {
		_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->A[i], items[i].x);
	}
	if (wls->type == WLSMOVING) {
        	for (int i = 0; i < numpts; i++) {
			real mlsw = 0.0;
			for (int j = 0; j < wls->dim; j++) {
				mlsw += (x[j]-items[i].x[j])*(x[j]-items[i].x[j]);
			}
			wls->D[i] = sqrt(1.0/(mlsw+1.0e-10));
			for (int j = 0; j < wls->numpoly; j++) {
				wls->A[i][j] *= wls->D[i];
			}
		}
	}
	wls_qr(numpts,wls->numpoly,wls->A,wls->Q);
	_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->b,x);
        real mindiag = fabs(wls->A[1][1]);
	for (int i = 2; i < wls->numpoly; i++) {
		real z = fabs(wls->A[i][i]);
		if (z < mindiag) {
			mindiag = z;
		}
	}
	wls_solve_system(numpts,wls->numpoly,wls->A,wls->Q,wls->w,wls->b,EPSMACH);
	if (wls->type == WLSMOVING) {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->D[i]*wls->w[i];
		}
	} else {
		for (int i = 0; i < numpts; i++) {
			w[i] = wls->w[i];
		}
	}
	return mindiag;
}

real wls_add_points_and_calc(wls_interpolator *wls, int numpts, int addnumpts, Point pts[], Point x, real w[]) {
        for (int i = 0; i < addnumpts; i++) {
		_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->A[numpts+i],pts[i]);
	}
	if (wls->type == WLSMOVING) {
        	for (int i = numpts; i < numpts+addnumpts; i++) {
			real mlsw = 0.0;
			for (int j = 0; j < wls->dim; j++) {
				mlsw += (x[j]-pts[i][j])*(x[j]-pts[i][j]);
			}
			wls->D[i] = sqrt(1.0/(mlsw+1.0e-10));
			for (int j = 0; j < wls->numpoly; j++) {
				wls->A[i][j] *= wls->D[i];
			}
		}
	}
        for (int i = 0; i < numpts; i++) {
		for (int j = numpts; j < numpts+addnumpts; j++) {
			wls->Q[i][j] = 0.0;
		}
	}
        for (int i = numpts; i < numpts+addnumpts; i++) {
		for (int j = 0; j < numpts; j++) {
			wls->Q[i][j] = 0.0;
		}
	}
        for (int i = numpts; i < numpts+addnumpts; i++) {
		for (int j = numpts; j < numpts+addnumpts; j++) {
			wls->Q[i][j] = 0.0;
		}
		wls->Q[i][i] = 1.0;
	}
	wls_qr_update(numpts,addnumpts,wls->numpoly,wls->A,wls->Q);
	_assembly_pts_matrix(wls->dim,wls->ord,numpts,wls->b,x);
        real mindiag = fabs(wls->A[1][1]);
	for (int i = 2; i < wls->numpoly; i++) {
		real z = fabs(wls->A[i][i]);
		if (z < mindiag) {
			mindiag = z;
		}
	}
	wls_solve_system(numpts+addnumpts,wls->numpoly,wls->A,wls->Q,wls->w,wls->b,EPSMACH);
	if (wls->type == WLSMOVING) {
		for (int i = 0; i < numpts+addnumpts; i++) {
			w[i] = wls->D[i]*wls->w[i];
		}
	} else {
		for (int i = 0; i < numpts+addnumpts; i++) {
			w[i] = wls->w[i];
		}
	}
	return mindiag;
}

bool
wls_is_good_enough(wls_interpolator *wls, int numusedpts, wls_item items[], const Point x, real w[])
{
	real interpolated_value = -eval_poly(x, wls->dim, wls->ord+1, wls->poly_coefs);
	for (int i = 0; i < numusedpts; i++) {
		interpolated_value += w[i]
			* eval_poly(items[i].x, wls->dim, wls->ord+1, wls->poly_coefs);
	}
	return FLT_EQ(interpolated_value, 0);
}

void wls_destroy(wls_interpolator *wls) {
	int n = wls->numpoly;
	int m = wls->maxnumpts;

        wls_freemreal(m,n,wls->A);
        wls_freemreal(m,m,wls->Q);
        wls_freevreal(m,wls->D);
        wls_freevreal(m,wls->b);
        wls_freevreal(m,wls->w);
	wls_freevreal(n,wls->poly_coefs);
	free(wls);
}

