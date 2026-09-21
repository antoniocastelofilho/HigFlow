// PLIC reconstruction in 3D: the plane that cuts the cell with the right volume.
// See the family header in hig-flow-vof-plic.c.
#if DIM == 3
#include "hig-flow-vof-plic-3D.h"
#include "hig-mesh-snapshot.h"
#include "hig-flow-vof-plic.h"
#include "hig-flow-vof-finite-difference-normal-curvature.h"

// Removidas daqui seis funcoes que ficaram sem chamador em 0a90511, quando o
// caso geral passou a inverter por bisseccao sobre a formula fechada de volume
// em vez de escolher entre dezesseis ramos:
//
//   newton_raphson                 -- orfa desde db52df8
//   raiz_cubica_bisseccao          -- introduzida e orfanada no mesmo dia
//   solver_equation_second_order   -- a quadratica de um corte
//   solver_equation_third_order    -- a cubica de dois cortes
//   solver_equation_third_order_2  -- a cubica de tres cortes
//   solver_equation_third_order_3  -- a forma fechada do regime de laje
//
// A aritmetica delas estava CORRETA; o que estava errado era a escolha de ramo
// que decidia qual chamar.  Ficam no historico, nao no arquivo.
//======================================================================

//if only one argument of the normal vector is non-zero ================
real solver_equation_n_nonzero(real volume, real n1, real dy, real dz){

	real Value=n1*volume/(dy*dz);
	
	return Value;
}
//======================================================================

//calculating the volume when all arguments of the normal vector is non-zero
real volume_3D(real n_x, real n_y, real n_z, real dx, real dy, real dz, real d, real aux_d){
	
	real aux_x, aux_y, aux_z, aux_x_y, aux_x_z, aux_y_z;
	n_x = fabs(n_x);
	n_y = fabs(n_y);
	n_z = fabs(n_z);
	real n1 = n_x+1e-14;
	real n2 = n_y+1e-14;
	real n3 = n_z+1e-14;	
	d = fabs(d);
	real volume;
	volume = 0.0;
	real tol_n = 1e-8;
	
	aux_x = d - n1*dx;
	aux_y = d - n2*dy;
	aux_z = d - n3*dz;
	aux_x_y = d - n1*dx - n2*dy;
	aux_x_z = d - n1*dx - n3*dz;
	aux_y_z = d - n2*dy - n3*dz;
			
	volume = (1/(6.0*n1*n2*n3))*(pow(d,3) - H(aux_x)*pow(d-n1*dx,3) - H(aux_y)*pow(d-n2*dy,3) 
	- H(aux_z)*pow(d-n3*dz,3) + H(aux_x_y)*pow(d-n1*dx-n2*dy,3)
	+ H(aux_x_z)*pow(d-n1*dx-n3*dz,3) + H(aux_y_z)*pow(d-n2*dy-n3*dz,3)); 
	
	if(aux_d<0){
		volume=dx*dy*dz-volume;
	}
	
	if(n_x<=tol_n && (n_y>tol_n && n_z>tol_n) || n_y<=tol_n  && 
	(n_x>tol_n && n_z>tol_n) || n_z<=tol_n  && (n_y>tol_n && n_x>tol_n)){
		volume = (dz/(2.0*n1*n2))*(pow(d,2) - H(aux_x)*pow(d-n1*dx,2) 
		- H(aux_y)*pow(d-n2*dy,2));
		if(aux_d<0){
			volume = dx*dy*dz-volume;
		}
	}
	return volume;
}
//======================================================================
//checking if there is an extra triangle.
real H(real x){
	
	real value = (x <= 0.0) ? 0.0 : 1.0;
	
	return value;
}
//======================================================================
real trans_center_to_p0(Point Delta,Point Normal,real d_from_center){
	
		real d_from_p0 =  -Normal[0] * 0.5 * Delta[0] - Normal[1] * 0.5 * Delta[1]
						- Normal[2] * 0.5 * Delta[2] - d_from_center;
		
	return d_from_p0;
}
//======================================================================
real trans_p0_to_center(Point Delta, Point Normal, real d_from_p0) {
	
	real nx = Normal[0];
	real ny = Normal[1];
	real nz = Normal[2];
	real n1 = fabs(nx);
	real n2 = fabs(ny);
	real n3 = fabs(nz);
	real tol_n = 1e-8;
	real n;
	
	real d_from_center =  - Normal[0] * 0.5 * Delta[0] - Normal[1] * 0.5 * Delta[1]
						- Normal[2] * 0.5 * Delta[2] - d_from_p0;
	
	return d_from_center;
}
//======================================================================
// Volume da parte da celula com n.x <= d_do_centro, exato, para |n_i| todos > 0.
// Inclusao-exclusao nos cantos: cada termo desconta a piramide que ja' saiu por
// uma face, cada par devolve o que foi descontado duas vezes, e assim por diante.
// Vale em TODOS os regimes -- e' isso que torna desnecessaria a escolha de ramo.
static real P3_(real x){ return (x > 0.0) ? x*x*x : 0.0; }

// G(u) = [P3(u) - P3(u-b)] / b, avaliado SEM subtrair termos quase iguais.
// E' o fator que permite tirar o menor b da inclusao-exclusao analiticamente.
static real G_(real u, real b){
	if (u <= 0.0)  return 0.0;
	if (u <= b)    return (b > 0.0) ? (u*u*u)/b : 3.0*u*u;
	return 3.0*u*u - 3.0*u*b + b*b;                 // b -> 0 da' 3u^2, o caso 2-D
}

static real volume_abaixo_do_plano(const real m[3], const real D[3], real d_centro){
	real b[3] = { m[0]*D[0], m[1]*D[1], m[2]*D[2] };
	real a = d_centro + 0.5*(b[0] + b[1] + b[2]);   // centro -> canto
	real vmax = D[0]*D[1]*D[2];

	// Ordena para que b[i0] seja o MENOR.  A forma direta da inclusao-exclusao
	// subtrai termos de ordem a^3 para obter um resultado de ordem 6.m0.m1.m2.V,
	// entao perde toda a precisao quando alguma componente e' pequena: medido,
	// com duas componentes em 1e-7 a resposta virava constante em relacao a
	// fracao.  Fatorando b[i0] o m[i0] cancela contra o denominador e some da
	// conta, e o caso b -> 0 cai continuamente na formula 2-D.
	int i0 = 0;
	if (b[1] < b[i0]) i0 = 1;
	if (b[2] < b[i0]) i0 = 2;
	int i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;

	real soma = G_(a, b[i0])
	          - G_(a - b[i1], b[i0]) - G_(a - b[i2], b[i0])
	          + G_(a - b[i1] - b[i2], b[i0]);
	real v = D[i0]*soma/(6.0*m[i1]*m[i2]);
	if (v < 0.0)  v = 0.0;
	if (v > vmax) v = vmax;
	return v;
}

//======================================================================
// O vetor normal deve apontar para fora da interface.
real distance_from_center_3D(Point Normal,Point Delta,real VOLUME){
	
	real tol_volume = 1e-8, tol_n = 1e-8;
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real n_x = fabs(nx), n_y = fabs(ny), n_z = fabs(nz);
	real dx = Delta[0], dy = Delta[1], dz = Delta[2];
	Point Value;
	real volume, d;
		
	real n;
	n = maximo(n_x, n_y, n_z);
	n = (n==n_x) ? nx : ((n==n_y) ? ny : nz);

	// Celula exatamente meio cheia: a caixa e' centralmente simetrica, entao um
	// plano pelo centro a corta em duas metades congruentes qualquer que seja a
	// normal -- a distancia ao centro e' 0, exata.  Sem esta saida a cadeia
	// estrita abaixo nao casa com VOLUME == 0.5*dx*dy*dz e a funcao cai fora do
	// fim sem retornar.  A versao 2D resolve o mesmo pelos FLT_EQ do topo
	// (hig-flow-vof-plic.c:89 e :92).
	if (FLT_EQ(VOLUME, 0.5*(dx*dy*dz))) return 0.0;
	//==================================================================
	//first case========================================================
	//==================================================================
	// Guardas com <=, alinhadas a' porta do chamador volume_left_line_origin_center
	// (que testa <= tol_n) e a' de parallel_case_volume.  Com `<`, uma componente
	// EXATAMENTE igual a tol_n escapava dos ramos degenerados e caia no caso geral,
	// onde a inclusao-exclusao sofre cancelamento catastrofico: com duas componentes
	// em 1e-8 o numerador e' da ordem de 1e-19 obtido subtraindo termos de 1e-4, e
	// a precisao relativa cai para ~15%.  Os ramos degenerados nao tem esse
	// problema, porque tratam a direcao pequena como livre em vez de dividir por ela.
	if ((n_x <= tol_n && n_y <= tol_n) || (n_x <= tol_n && 
	n_z <= tol_n) || (n_y <= tol_n && n_z <= tol_n)){
		
		real n1, n_1;	
		n1 = (n_x <= tol_n && n_y <= tol_n) ? nz : ((n_x <= tol_n && n_z <= tol_n) ? ny : nx);
		n_1 = (n_x <= tol_n && n_y <= tol_n) ? -n_z : ((n_x <= tol_n && n_z <= tol_n) ? -n_y : -n_x);
		
		// kd: eixo dominante.  A altura do corte e' VOLUME dividido pelas duas
		// extensoes PERPENDICULARES a ele, e o deslocamento e' meia extensao DELE.
		// Antes eram dy*dz e Delta[0] fixos -- certo so' quando o dominante era x,
		// e invisivel em celula cubica, onde os tres Delta coincidem.
		int kd = (n_x <= tol_n && n_y <= tol_n) ? 2 : ((n_x <= tol_n && n_z <= tol_n) ? 1 : 0);
		real perp1 = Delta[(kd+1)%3], perp2 = Delta[(kd+2)%3];
		d  = solver_equation_n_nonzero(VOLUME, n1, perp1, perp2);
		d = 0.5*Delta[kd] - sign(n)*d;
		return d;
		
	}else if(n_x <= tol_n || n_y <= tol_n || n_z <= tol_n){
	//==================================================================
	//Second case=======================================================
	//==================================================================
		// Exatamente UMA componente e' nula -- duas ou mais ja' foram tratadas no
		// ramo acima.  Entao o problema e' 2-D no plano dos dois eixos vivos, com
		// espessura L0 ao longo do eixo nulo, e a distancia ao centro da celula 3-D
		// e' a mesma distancia ao centro do retangulo, porque a normal nao tem
		// componente na direcao nula.
		//
		// distance_from_center (hig-flow-vof-plic.c:68) ja' resolve exatamente essa
		// geometria, com os tres regimes e os testes de igualdade no topo.  O que
		// havia aqui era so' `d = sqrt(2*volume*n_1*n_2/L0)`, a formula do TRIANGULO
		// DE CANTO, valida enquanto o triangulo cabe no retangulo -- sem verificacao
		// de regime alguma.  Media em celula cubica com eixo nulo z: acertava de
		// frac 0,05 a 0,35 e errava em 0,45.
		//
		// A 2-D foi conferida contra area analitica antes de ser adotada aqui: 0
		// discordancias em 1.248 casos por geometria, em quatro formas de celula,
		// com a MESMA convencao de lado que a versao 3-D.  Ela nao escreve em
		// Normal, e nem ela nem os seus auxiliares leem o indice [2], entao receber
		// um Point de tres posicoes e' seguro.
		int  k0 = (n_x <= tol_n) ? 0 : ((n_y <= tol_n) ? 1 : 2);
		int  ka = (k0 + 1) % 3, kb = (k0 + 2) % 3;
		real L0 = Delta[k0];
		Point N2, D2;
		N2[0] = Normal[ka];  N2[1] = Normal[kb];  N2[2] = 0.0;
		D2[0] = Delta[ka];   D2[1] = Delta[kb];   D2[2] = L0;
		return distance_from_center(N2, D2, VOLUME / L0);
	}else{
	//==================================================================
	//Third case =======================================================
	//==================================================================
		// Nenhuma componente nula.  Antes havia aqui dois blocos espelhados com
		// oito ramos cada, escolhidos comparando c = 6.n1n2n3.V com (n_i.D_i)^3.
		// Essa comparacao equivale a comparar d^3 com a_i^3, o que so' vale no
		// PRIMEIRO regime, onde a piramide cabe inteira na celula.  Fora dele a
		// formula da piramide superestima o volume para um dado d, logo cbrt(c) <= d,
		// e a assimetria machuca: `c > a_i^3` garante `d > a_i`, mas `c <= a_i^3` NAO
		// garante `d <= a_i`.  A classificacao subcontava cortes sistematicamente.
		// Medido: onde ela acertava o ramo, 2% de erro; onde errava, 83%.
		//
		// A funcao volume e' monotona em d e tem forma fechada valida em TODOS os
		// regimes, entao a escolha de ramo e' desnecessaria: basta inverter por
		// bisseccao no intervalo fisico [-soma(a_i)/2, +soma(a_i)/2].  Custa ~100
		// avaliacoes de um polinomio por celula de interface, fora de laco quente.
		//
		// Convencao, a mesma dos outros ramos: VOLUME e' o lado n.x >= d.
		real m[3]  = { n_x, n_y, n_z };          // ja' sao os modulos
		real amax  = 0.5*(m[0]*dx + m[1]*dy + m[2]*dz);
		real celv  = dx*dy*dz;
		real lo = -amax, hi = amax;
		for (int it = 0; it < 200; it++) {
			real mid = 0.5*(lo + hi);
			if (celv - volume_abaixo_do_plano(m, Delta, mid) > VOLUME) lo = mid;
			else                                                       hi = mid;
			if (hi - lo < 1e-15*amax) break;
		}
		return 0.5*(lo + hi);
	}
}
//======================================================================

real parallel_case_volume(Point Normal,Point Delta,real d_from_center,real tol_n){
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real dx = Delta[0], dy = Delta[1], dz = Delta[2];
	int  n_x, n_y, n_z;
	// n_1..n_3 guardam |Normal[i]|, que e' real em [0,1]: declarados como int,
	// a componente dominante truncava para 0 sempre que nao fosse exatamente
	// 1.0, e Normal virava o vetor nulo.  trans_center_to_p0 devolvia entao
	// -d em vez de |n|.Delta/2 - d.  A inversa trans_p0_to_center (linha 155)
	// ja' usa Normal[i] direto, sem truncar.
	real n_1, n_2, n_3;
	real volume, d;
	
	n_1 = fabs(Normal[0]);
	n_2 = fabs(Normal[1]);
	n_3 = fabs(Normal[2]);
	Normal[0] = -n_1;
	Normal[1] = -n_2;
	Normal[2] = -n_3;
	
	d = fabs(d_from_center);
	d = trans_center_to_p0(Delta,Normal,d);
	
	if(fabs(nx)<=tol_n && fabs(ny)<=tol_n){
		n_x = 0;
		n_y = 0;
		n_z = 1;
	}else if(fabs(nx)<=tol_n && fabs(nz)<=tol_n){
		n_x = 0;
		n_y = 1;
		n_z = 0;
	}else if(fabs(ny)<=tol_n && fabs(nz)<=tol_n){
		n_x = 1;
		n_y = 0;
		n_z = 0;
	}else{
		printf("=+=+=+= Error in parallel_case_volume: normal (%g %g %g) nao e'"
		       " paralela a eixo algum com tol_n=%g =+=+=+=\n", nx, ny, nz, tol_n);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
		
	
	volume = (n_x==1) ? d*(dy*dz) : ((n_y==1) ? d*(dx*dz) : d*(dy*dx));

	return volume;
}
//======================================================================

real volume_left_line_origin_center(Point Normal,Point Delta,real d_from_center){
	
	real aux_d = d_from_center;
	d_from_center = fabs(d_from_center);
	real d, volume;
	real dx  = Delta[0], dy  = Delta[1], dz  = Delta[2];
	real nx = Normal[0], ny = Normal[1], nz = Normal[2];
	real n_x = fabs(nx), n_y = fabs(ny), n_z = fabs(nz);
	real tol_n = 1e-8;
	
	real dmax   = n_x * 0.5 * Delta[0] + n_y * 0.5 * Delta[1] + n_z * 0.5 * Delta[2];
	real volmax = dx*dy*dz;
	
	if((d_from_center >= dmax) && aux_d<0) {
		 volume = dx*dy*dz;
	} else if((d_from_center >= dmax) && aux_d>0) {
		 volume = 0.0;
	}else{
		if((n_x<=tol_n && n_y<=tol_n ) || (n_x<=tol_n && n_z<=tol_n ) || (n_y<=tol_n && n_z<=tol_n )) {
			volume = parallel_case_volume(Normal, Delta, d_from_center, tol_n);
			if(aux_d<0){
				volume=dx*dy*dz-volume;
			}
		}else if(n_x<=tol_n && (n_y>tol_n && n_z>tol_n) || n_y<=tol_n  && 
		(n_x>tol_n && n_z>tol_n) || n_z<=tol_n  && (n_y>tol_n && n_x>tol_n)){
			
			real n1, n2, n3, n_1, n_2, n_3;
			n1 = (n_x < tol_n) ? ny : nx;
			n2 = (n_x < tol_n || n_y < tol_n) ? nz : ny;
			n3 = 0.0;
			n_1 = (n_x < tol_n) ? -n_y : -n_x;
			n_2 = (n_x < tol_n || n_y < tol_n) ? -n_z : -n_y;
			n_3 = 0.0;
			Normal[0] = n_1;
			Normal[1] = n_2;
			Normal[2] = n_3;
			
			d = trans_center_to_p0(Delta, Normal, d_from_center);
			volume = volume_3D(n_1, n_2, n_3, dx, dy, dz, d, aux_d);
		}else if(n_x>tol_n && n_y>tol_n && n_z>tol_n){
			
			real n1, n2, n3, n_1, n_2, n_3;
			n_1 = -n_x;
			n_2 = -n_y;
			n_3 = -n_z;
			Normal[0] = n_1;
			Normal[1] = n_2;
			Normal[2] = n_3;
			
			d = trans_center_to_p0(Delta, Normal, d_from_center);
			volume = volume_3D(nx, ny, nz, dx, dy, dz, d, aux_d);
		}
	}
	
	volume = minimo(volume,volmax);
	
	return volume;
}
//======================================================================
real maximo(real nx, real ny, real nz){
		real max;		
		max = (nx > ny && nx > nz) ? nx : ((ny > nx && ny > nz) ? ny : nz);
		
		return max;
}
//======================================================================
real minimo(real vol1, real vol2){
		real min;
		min = (vol1 < vol2) ? vol1 : vol2;
		
		return min;
}
//======================================================================
void higflow_compute_distance_multiphase_3D(higflow_solver *ns) {
	if (ns->contr.flowtype == 2) {
		real IF[DIM];
		// Get the local sub-domain for the cells
		sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);

		// Get the map for the domain properties
		mp_mapper *mp = sd_get_domain_mapper(sdp);

		// Loop for each cell
		higcit_celliterator *it;

		{
		const hig_mesh_snapshot *hms = sd_get_snapshot(sdp);
		for(int clid = 0; clid < hms->n; clid++) {
			// Get the cell

			// Get the cell identifier

			// Get the center of the cell
			Point center;
			hms_center(hms, clid, center);

			// Get the delta of the cell
			Point delta;
			hms_delta(hms, clid, delta);

			Point p;
			p[0] = center[0];
			p[1] = center[1];
			p[2] = center[2];
			real frac = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			
			Point Normal;
			Normal[0] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[0], ns->ed.stn);
			Normal[1] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[1], ns->ed.stn);
			Normal[2] = compute_value_at_point(sdp, center, p, 1.0, ns->ed.mult.dpnormal[2], ns->ed.stn);
						
			real tol_n = 1e-8;
			if (fabs(Normal[0]) < tol_n && fabs(Normal[1]) < tol_n && fabs(Normal[2]) < tol_n){
				continue;
			}
			
			//real frac = compute_value_at_point(sdp, center, center, 1.0, ns->ed.mult.dpfracvol, ns->ed.stn);
			real volume = frac*delta[0]*delta[1]*delta[2];
			
			real distance3D = distance_from_center_3D(Normal,delta,volume);
			
			dp_set_value(ns->ed.mult.dpdistance, clid, distance3D);

			//arquivoDist(NULL,center[0],center[1],distance);
		}
		}
		// Destroy the iterator
		// Sync the distributed pressure property
		dp_sync(ns->ed.mult.dpdistance);
	}
}
#endif
