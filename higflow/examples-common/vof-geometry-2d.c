// ---------------------------------------------------------------------------
// Geometria 2-D para a fracao volumetrica inicial dos exemplos VOF.
//
// Estas seis funcoes eram IDENTICAS -- byte a byte -- em cinco exemplos:
// VOF, VOF_Gptt, VOF_Oldroyd, ElectroOsmotic e DynamicMeshAdapt.  Sao recorte de
// poligono e calculo de area: nao dependem de nada do exemplo, entao repeti-las
// era custo sem contrapartida.
//
// O que NAO veio para ca': get_fracvolN.  Ele tem a mesma forma nos cinco, mas
// avalia func(p) nos cantos da celula -- e func e' a definicao da FORMA, que cada
// exemplo escreve a seu modo (disco de Zalesak, circulo, ...).  Traze-lo criaria
// dependencia implicita de ordem de #include em troca de 168 linhas; a ponte fica
// onde esta', com o exemplo que a define.
//
// Incluido por #include, nao compilado a' parte: os Makefiles dos exemplos so'
// alcancam os .c do proprio diretorio.
// ---------------------------------------------------------------------------
#ifndef HIGFLOW_EXAMPLES_VOF_GEOMETRY_2D
#define HIGFLOW_EXAMPLES_VOF_GEOMETRY_2D

void Intersec(Point p0, Point p1, real f0, real f1, Point p) {
	real lambda = -f0/(f1-f0);
	p[0] = p0[0] + lambda*(p1[0]-p0[0]);
	p[1] = p0[1] + lambda*(p1[1]-p0[1]);
	return;
}
real volume_convex_facet(int n, real P[][DIM]) {
	real area   = 0.0;
	real altura = 0.0;
	for (int i = 0; i< n-1; i++) {
		area   += P[i+1][0] - P[i][0];
		altura += P[i][1];
	}
	altura += P[n-1][1];

	return altura*area/n;
}
real volume_convex_13(Point p0, Point p01, Point p02){
	real facet0[2][DIM], facet1[2][DIM], facet2[2][DIM];
	real volume = 0.0;

	facet0[0][0] = p0[0];
	facet0[0][1] = p0[1];
	facet0[1][0] = p01[0];
	facet0[1][1] = p01[1];

	volume += volume_convex_facet(2,facet0);

	facet1[0][0] = p01[0];
	facet1[0][1] = p01[1];
	facet1[1][0] = p02[0];
	facet1[1][1] = p02[1];

	volume += volume_convex_facet(2,facet1);

	facet2[0][0] = p02[0];
	facet2[0][1] = p02[1];
	facet2[1][0] = p0[0];
	facet2[1][1] = p0[1];

	volume += volume_convex_facet(2,facet2);

	return volume;
}
real volume_convex_22(Point p0, Point p1, Point p13, Point p02){
	real facet0[2][DIM], facet1[2][DIM], facet2[2][DIM], facet3[2][DIM];
	real volume = 0.0;

	facet0[0][0] = p0[0];
	facet0[0][1] = p0[1];
	facet0[1][0] = p1[0];
	facet0[1][1] = p1[1];

	volume += volume_convex_facet(2,facet0);

	facet1[0][0] = p1[0];
	facet1[0][1] = p1[1];
	facet1[1][0] = p13[0];
	facet1[1][1] = p13[1];

	volume += volume_convex_facet(2,facet1);

	facet2[0][0] = p13[0];
	facet2[0][1] = p13[1];
	facet2[1][0] = p02[0];
	facet2[1][1] = p02[1];

	volume += volume_convex_facet(2,facet2);

	facet3[0][0] = p02[0];
	facet3[0][1] = p02[1];
	facet3[1][0] = p0[0];
	facet3[1][1] = p0[1];

	volume += volume_convex_facet(2,facet3);

	return volume;
}
real square_case_13(Point center, Point delta, Point p0, Point p1, Point p2, Point p3, real f0, real f1, real f2, real f3, int var){
	Point p01, p02, p13, p23;
	real value = 0.0;
	int caso;
	if (var == 1) {
		if (f0 > 0.0) {
			//caso = 2;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p0, p2, f0, f2, p02);
			value = fabs(volume_convex_13(p0, p01, p02));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f1 > 0.0) {
			//caso = 3;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p1, p3, f1, f3, p13);
			value = fabs(volume_convex_13(p1, p13, p01));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f2 > 0.0) {
			//caso = 4;
			Intersec(p0, p2, f0, f2, p02);
			Intersec(p2, p3, f2, f3, p23);
			value = fabs(volume_convex_13(p2, p02, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f3 > 0.0) {
			//caso = 5;
			Intersec(p1, p3, f1, f3, p13);
			Intersec(p2, p3, f2, f3, p23);
			value = fabs(volume_convex_13(p3, p23, p13));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		}
	} else {
		if (f0 <= 0.0) {
			//caso = 6;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p0, p2, f0, f2, p02);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p0, p01, p02));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f1 <= 0.0) {
			//caso = 7;
			Intersec(p0, p1, f0, f1, p01);
			Intersec(p1, p3, f1, f3, p13);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p1, p13, p01));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f2 <= 0.0) {
			//caso = 8;
			Intersec(p0, p2, f0, f2, p02);
			Intersec(p2, p3, f2, f3, p23);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p2, p02, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, value);
			return value;
		} else if (f3 <= 0.0) {
			//caso = 9;
			Intersec(p1, p3, f1, f3, p13);
			Intersec(p2, p3, f2, f3, p23);
			value = delta[0] * delta[1] - fabs(volume_convex_13(p3, p13, p23));
			//printf("Caso %d, volume = %.18lf\n", caso, volume);
			return value;
		}
	}
	return value;
}
real square_case_22(Point center, Point delta, Point p0, Point p1, Point p2, Point p3, real f0, real f1, real f2, real f3){
	Point p01, p02, p13, p23;
	real value  = 0.0;
	int caso;
	if ((f0 > 0.0) && (f1 > 0.0)){
		//caso = 10;
		Intersec(p0,p2,f0,f2,p02);
		Intersec(p1,p3,f1,f3,p13);
		value = fabs(volume_convex_22(p0, p1, p13, p02));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f0 > 0.0) && (f2 > 0.0)){
		//caso = 11;
		Intersec(p0,p1,f0,f1,p01);
		Intersec(p2,p3,f2,f3,p23);
		value = fabs(volume_convex_22(p2, p0, p01, p23));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f1 > 0.0) && (f3 > 0.0)){
		//caso = 12;
		Intersec(p0,p1,f0,f1,p01);
		Intersec(p2,p3,f2,f3,p23);
		value = fabs(volume_convex_22(p1, p3, p23, p01));
		//printf("Caso %d, volume = %.18lf\n", caso, value);
		return value;
	} else if ((f2 > 0.0) && (f3 > 0.0)){
		//caso = 13;
		Intersec(p0,p2,f0,f2,p02);
		Intersec(p1,p3,f1,f3,p13);
		value = fabs(volume_convex_22(p3, p2, p02, p13));
		//printf("Caso %d, volume = %.18lf", caso, value);
		return value;
	} else {
		// caso = 14 e caso = 15
		printf("######## Malha pouco refinada para este problema ########\n");
		exit(1);
	}
}

#endif

// ---------------------------------------------------------------------------
// Fracao volumetrica de uma celula 2-D cortada pela interface.
//
// Os cinco exemplos VOF tinham este corpo identico byte a byte.  O que varia
// entre eles NAO esta' aqui: e' o `func` abaixo -- a funcao level-set que
// define a forma inicial --, que cada exemplo continua definindo por conta
// propria.  Esta funcao so' amostra `func` nos quatro cantos e despacha para
// os casos analiticos acima.
// ---------------------------------------------------------------------------
real func(Point p);   // definida por cada exemplo: a forma inicial da interface

real get_fracvolN(Point center, Point delta, real t) {
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
	} else if ((var == 1)||(var == 3)){
		value = square_case_13(center, delta, p0, p1, p2, p3, f0, f1, f2, f3, var);
	} else {
		value = square_case_22(center, delta, p0, p1, p2, p3, f0, f1, f2, f3);
	}

	return value;
}
