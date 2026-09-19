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
