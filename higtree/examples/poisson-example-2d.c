/**
 * @file poisson-example-2d.c
 * @brief Solucionador da equação de Poisson 2D usando métodos de diferenças finitas
 * 
 * Este programa resolve a equação de Poisson 2D: ∇²u = g(x,y) com várias condições de contorno
 * e funções de teste diferentes. Utiliza uma estrutura de dados de árvore hierárquica (higtree)
 * para suportar malhas uniformes e não-uniformes com refinamento adaptativo.
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "utils.h"
#include "higtree.h"
#include "higtree-io.h"
#include "mapper.h"
#include "wls.h"
#include "glib.h"
#include "domain.h"
#include "solver.h"

#define DEBUG
#include "Debug-c.h"

/* Variável global para seleção da função teste */
int func = 0;

/**
 * @brief Função analítica u(x,y) - solução exata
 * @param x Coordenada x
 * @param y Coordenada y
 * @return Valor da função no ponto (x,y)
 */
real f(real x, real y) {
	switch(func) {
		case 0:	return sin(x)*cos(y);                    // Solução suave
		case 1: return 4*x*x*y+2*x*x-4*y*x-5*x;        // Polinômio
		case 2: return cos(x*y)+sin(x+y);              // Combinação trigonométrica
		case 3: return cos(x-y)*sin(x+y);              // Produto trigonométrico
		case 4: return cos(x+sin(y));                  // Composição trigonométrica
		case 5: return x*x+y*y-x*y-1.1234;            // Forma quadrática
	}
	printf("undefined function! Should be in 0..5\n");
	exit(0);
}

/**
 * @brief Termo fonte g(x,y) = ∇²u (lado direito da equação de Poisson)
 * @param x Coordenada x
 * @param y Coordenada y
 * @return Valor do Laplaciano no ponto (x,y)
 */
real g(real x, real y) {
	switch(func) {
		case 0: return -2.0*sin(x)*cos(y);            // Laplaciano de sin(x)cos(y)
		case 1: return 4 + 8*y;                       // Laplaciano do polinômio
		case 2: return (-x*x-y*y)*cos(x*y)-2*sin(x+y); // Laplaciano complexo
		case 3: return -4*cos(x-y)*sin(x+y);          // Laplaciano do produto
		case 4: return -(1+cos(y)*cos(y))*cos(x+sin(y))+sin(y)*(sin(x+sin(y))); // Laplaciano composto
		case 5: return 4.0;                           // Laplaciano constante
	}
}

/**
 * @brief Derivada parcial ∂u/∂y
 * @param x Coordenada x
 * @param y Coordenada y
 * @return Valor da derivada em y
 */
real dfdy(real x, real y) {
	switch(func) {
		case 0: return -sin(x)*sin(y);
		case 1: return -4*x+4*x*x;
		case 2: return cos(x+y) - x*sin(x*y);
		case 3: return cos(x-y)*cos(x+y) + sin(x-y)*sin(x+y);
		case 4: return -cos(y)*sin(x+sin(y));
		case 5: return 2.0*y-x;
	}
}

/**
 * @brief Derivada parcial ∂u/∂x
 * @param x Coordenada x
 * @param y Coordenada y
 * @return Valor da derivada em x
 */
real dfdx(real x, real y) {
	switch(func) {
		case 0: return cos(x)*cos(y);
		case 1: return -5 + 4*x - 4*y + 8*x*y;
		case 2: return cos(x+y) - y*sin(x*y);
		case 3: return cos(x-y)*cos(x+y) - sin(x-y)*sin(x+y);
		case 4: return -sin(x+sin(y));
		case 5: return 2.0*x-y;
	}
}

/**
 * @brief Retorna string representando a função atual
 * @return String com expressão matemática da função
 */
const char * f_str() {
	switch(func) {
		case 0:	return "sin(x)*cos(y)";
		case 1: return "4*x*x*y+2*x*x-4*y*x-5*x";
		case 2: return "cos(x*y)+sin(x+y)";
		case 3: return "cos(x-y)*sin(x+y)";
		case 4: return "cos(x+sin(y))";
		case 5: return "x*x+y*y-x*y-1.1234";
	}
	return "undefined function!";
}

/**
 * @brief Refina regiões específicas da malha para criar não-uniformidade
 * @param root Célula raiz da árvore hierárquica
 * @param nc Fator de refinamento em cada direção
 * 
 * Aplica refinamento em duas regiões retangulares aninhadas para criar
 * uma malha com resolução variável
 */
void __refine_non_uniform(hig_cell *root, int nc[]) {
	Point bl, bh;
	/* Primeira região de refinamento: [-2.59,-0.59]×[0.59,0.59] */
	POINT_ASSIGN_REALS(bl, -2.59, -0.59);
	POINT_ASSIGN_REALS(bh, 0.59, 0.59);
	higcit_celliterator *it;
	for(it = higcit_create_bounding_box(root, bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		hig_refine_uniform(c, nc);
	}
	higcit_destroy(it);
	
	/* Segunda região de refinamento (mais interna): [-2.19,-0.19]×[0.19,0.19] */
	POINT_ASSIGN_REALS(bl, -2.19, -0.19);
	POINT_ASSIGN_REALS(bh, 0.19, 0.19);
	for(it = higcit_create_bounding_box(root, bl, bh); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		hig_refine_uniform(c, nc);
	}
	higcit_destroy(it);
}

/**
 * @brief Exibe instruções de uso do programa
 */
void __usage() {
	printf("./poisson-example \n\t[--side s (default 10)] \n\t[--refinement-level r (default 0)]\n\t[--uniform u (default: 1)]\n\t[--bctypes (default: DDDD)]\n\t[--writevtk vtkfilename]\n\n");
}

/**
 * @brief Processa argumentos de linha de comando
 * @param argc Contador de argumentos
 * @param argv Vetor de argumentos
 * @param side Número de células por lado (saída)
 * @param refinementLevel Níveis de refinamento (saída)
 * @param isUniform Flag para malha uniforme (saída)
 * @param bctypes Tipos de condições de contorno (saída)
 * @param refineFactor Fator de refinamento (saída)
 * @param writevtk Flag para escrita VTK (saída)
 * @param vtkfilename Nome do arquivo VTK (saída)
 * @param maxlen Tamanho máximo do nome do arquivo
 */
void __get_opts(int argc, char *argv[], int *side, int *refinementLevel, int *isUniform, int *bctypes, int *refineFactor, int *writevtk, char vtkfilename[], int maxlen) {
	while(argc > 2) {
		if(strcmp(argv[1], "--side") == 0) {
			argv++; argc--;
			*side = atoi(argv[1]);
			argv++; argc--;
		} else if(strcmp(argv[1], "--refinement-level") == 0) {
			argv++; argc--;
			*refinementLevel = atoi(argv[1]);
			argv++; argc--;
		} else if(strcmp(argv[1], "--uniform") == 0) {
			argv++; argc--;
			*isUniform = atoi(argv[1]);
			argv++; argc--;
		} else if(strcmp(argv[1], "--refinement-factor") == 0) {
			argv++; argc--;
			*refineFactor = atoi(argv[1]);
			argv++; argc--;
		} else if(strcmp(argv[1], "--bctypes") == 0) {
			argv++; argc--;
			*bctypes = 0;
			/* Processa string de 4 caracteres (N/D para cada borda) */
			for(int i = 0; i < 4; i++) {
				if (argv[1][i] == 'N' || argv[1][i] == 'n') {
					*bctypes |= (1 << i);
				} else if (argv[1][i] != 'D' && argv[1][i] != 'd') {
					__usage();
				}
			}
			argv++; argc--;
		} else if(strcmp(argv[1], "--writevtk") == 0) {
			argv++; argc--;
			*writevtk = 1;
			strncpy(vtkfilename, argv[1], maxlen);
			argv++; argc--;
		} else if(strcmp(argv[1], "--func") == 0) {
			argv++; argc--;
			func = atoi(argv[1]);
			argv++; argc--;
		} else {
			__usage();
			break;
		}
	}
}

/**
 * @brief Exibe configurações atuais do problema
 * @param side Número de células por lado
 * @param size Número total de células
 * @param refinementLevel Níveis de refinamento
 * @param isUniform Flag de malha uniforme
 * @param bctypes Tipos de condições de contorno
 * @param refineFactor Fator de refinamento
 * @param writevtk Flag de escrita VTK
 * @param vtkfilename Nome do arquivo VTK
 */
void __print_opts(int side, int size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, char vtkfilename[]) {
	printf("func = %s\n", f_str());
	printf("side = %d\n", side);
	printf("number of cells = %d\n", size);
	printf("refinement levels = %d\n", refinementLevel);
	printf("%s\n", (isUniform?"uniform":"non-uniform"));
	printf("refinement factor = %d\n", refineFactor);
	printf("bc types:\n\t  Left = %s\n\t Rigth = %s\n\tBottom = %s\n\t   Top = %s\n",
		((bctypes & (1<<0))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<1))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<2))?"NEUMANN":"DIRICHLET"),
		((bctypes & (1<<3))?"NEUMANN":"DIRICHLET"));
	if(writevtk) {
		printf("vtkfilename = %s\n", vtkfilename);
	} else {
		printf("no vtk\n");
	}
}

/**
 * @brief Constrói o domínio computacional com malha e condições de contorno
 * @param side Número inicial de células por lado
 * @param size Ponteiro para armazenar número total de células (saída)
 * @param refinementLevel Níveis de refinamento
 * @param isUniform Flag para malha uniforme
 * @param bctypes Tipos de condições de contorno
 * @param refineFactor Fator de refinamento
 * @param writevtk Flag para escrita VTK
 * @param vtkfilename Nome do arquivo VTK
 * @return Ponteiro para o domínio simulação construído
 */
sim_domain * __build_domain(int side, int *size, int refinementLevel, int isUniform, int bctypes, int refineFactor, int writevtk, char vtkfilename[]) {
	Point l, h;
	/* Domínio principal: [-1,1]×[-1,1] */
	POINT_ASSIGN_REALS(l, -1.0, -1.0);
	POINT_ASSIGN_REALS(h, 1.0, 1.0);
	
	/* Cria árvore hierárquica raiz */
	hig_cell *root = hig_create_root(l, h);
	int nc[DIM];
	POINT_ASSIGN_INTS(nc, side, side);
	hig_refine_uniform(root, nc);
	
	/* Aplica refinamento não-uniforme se solicitado */
	POINT_ASSIGN_INTS(nc, refineFactor, refineFactor);
	higcit_celliterator *it;
	if (!isUniform) {
		__refine_non_uniform(root, nc);
	}
	
	/* Refinamento uniforme adicional */
	POINT_ASSIGN_INTS(nc, 2, 2);
	for (int i = 0; i < refinementLevel; i++) {
		for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			hig_refine_uniform(c, nc);
		}
		higcit_destroy(it);
		side *= 2;  // Atualiza contagem de células
	}
	
	/* Exporta malha para VTK se solicitado */
	if (writevtk) {
		FILE *fd = fopen(vtkfilename, "w");
		higio_print_in_vtk2d(fd, root);
		fclose(fd);
	}

	/* Cria mapeamento do domínio */
	mp_mapper *m = mp_create();
	it = higcit_create_all_leaves(root);
	*size = mp_assign_from_celliterator(m, it, 0);
	higcit_destroy(it);

	/* Cria domínio de simulação */
	sim_domain *sd = sd_create(m);
	sd_add_higtree(sd, root);

	/* Configura condições de contorno para as 4 bordas */
	Point ls[4], hs[4];
	int ncs[4][DIM];
	hig_cell *bcs[4];
	mp_mapper *ms[4];
	sim_boundary *sbs[4];
	bc_type types[4];

	/* Borda esquerda */
	POINT_ASSIGN_REALS(ls[0], l[0] - EPSDELTA, l[1]);
	POINT_ASSIGN_REALS(hs[0], l[0] + EPSDELTA, h[1]);
	POINT_ASSIGN_INTS(ncs[0], 1, side);
	types[0] = ((bctypes&(1<<0))?NEUMANN:DIRICHLET);

	/* Borda direita */
	POINT_ASSIGN_REALS(ls[1], h[0] - EPSDELTA, l[1]);
	POINT_ASSIGN_REALS(hs[1], h[0] + EPSDELTA, h[1]);
	POINT_ASSIGN_INTS(ncs[1], 1, side);
	types[1] = ((bctypes&(1<<1))?NEUMANN:DIRICHLET);

	/* Borda inferior */
	POINT_ASSIGN_REALS(ls[2], l[0], l[1] - EPSDELTA);
	POINT_ASSIGN_REALS(hs[2], h[0], l[1] + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[2], side, 1);
	types[2] = ((bctypes&(1<<2))?NEUMANN:DIRICHLET);

	/* Borda superior */
	POINT_ASSIGN_REALS(ls[3], l[0], h[1] - EPSDELTA);
	POINT_ASSIGN_REALS(hs[3], h[0], h[1] + EPSDELTA);
	POINT_ASSIGN_INTS(ncs[3], side, 1);
	types[3] = ((bctypes&(1<<3))?NEUMANN:DIRICHLET);

	/* Cria e configura cada borda */
	for(int i = 0; i < 4; i++) {
		bcs[i] = hig_create_root(ls[i], hs[i]);
		hig_refine_uniform(bcs[i], ncs[i]);
		
		/* Aplica refinamento não-uniforme na borda esquerda se necessário */
		POINT_ASSIGN_INTS(nc, 1, refineFactor);
		if (!isUniform && i == 0) {
			__refine_non_uniform(bcs[i], nc);
			/* Debug: exporta borda para VTK */
			FILE *fd = fopen("/tmp/bc.vtk", "w");
			higio_print_in_vtk2d(fd, bcs[i]);
			fclose(fd);
		}
		
		ms[i] = mp_create();
		it = higcit_create_all_leaves(bcs[i]);
		mp_assign_from_celliterator(ms[i], it, 0);
		higcit_destroy(it);
		sbs[i] = sb_create(bcs[i], types[i], ms[i]);
		sd_add_boundary(sd, sbs[i]);
	}

	/* Define valores das condições de contorno */
	for(int i = 0; i < 4; i++) {
		for(it = higcit_create_all_leaves(bcs[i]); !higcit_isfinished(it); higcit_nextcell(it)) {
			hig_cell *c = higcit_getcell(it);
			Point ccenter;
			hig_get_center(c, ccenter);
			int cgid = mp_lookup(ms[i], hig_get_cid(c));
			
			if(types[i] == DIRICHLET) {
				/* Condição Dirichlet: u = f(x,y) */
				sb_set_value(sbs[i], cgid, f(ccenter[0], ccenter[1]));
			} else {
				/* Condição Neumann: ∂u/∂n = valor */
				switch (i) {
					case 0: sb_set_value(sbs[i], cgid, -dfdx(ccenter[0], ccenter[1]));  // Borda esquerda: -∂u/∂x
						break;
					case 1: sb_set_value(sbs[i], cgid, dfdx(ccenter[0], ccenter[1]));   // Borda direita: ∂u/∂x
						break;
					case 2: sb_set_value(sbs[i], cgid, -dfdy(ccenter[0], ccenter[1]));  // Borda inferior: -∂u/∂y
						break;
					case 3: sb_set_value(sbs[i], cgid, dfdy(ccenter[0], ccenter[1]));   // Borda superior: ∂u/∂y
						break;
					default : exit(0);
				}
			}
		}
		higcit_destroy(it);
	}
	return sd;
}

/**
 * @brief Função para cálculo de estêncil (atualmente vazia)
 */
void __calc_poisson_stencil(sim_domain *sd, real cx, real cy, real dx, real dy, sim_stencil *stn) {
}

/**
 * @brief Função principal
 * @param argc Contador de argumentos
 * @param argv Vetor de argumentos
 * @return 0 em caso de sucesso
 */
int main(int argc, char *argv[]) {
	int nc[DIM];
	int side = 10;
	int refinementLevel = 0;
	int isUniform = 1;
	int bctypes = 0;
	int writevtk = 0;
	int refineFactor = 2;
	int size;
	char vtkfilename[100];

	/* Inicialização do sistema higtree */
	higtree_initialize(&argc, &argv);
	DEBUG_DIFF_TIME;
	
	/* Processamento de argumentos */
	__get_opts(argc, argv, &side, &refinementLevel, &isUniform, &bctypes, &refineFactor, &writevtk, vtkfilename, 100);

	/* Construção do domínio */
	sim_domain *sd = __build_domain(side, &size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, vtkfilename);

	__print_opts(side, size, refinementLevel, isUniform, bctypes, refineFactor, writevtk, vtkfilename);

	/* Preparação para solução do sistema linear */
	sim_stencil *stn = stn_create();

	printf("generating matrix...\n");
	DEBUG_DIFF_TIME;
	solver *slv = slv_create(SOLVER_ANY, 0, size);
	slv_set_maxnonzeros(slv, 300);

	/* Montagem da matriz do sistema usando diferenças finitas */
	higcit_celliterator *it;
	hig_cell *root = sd_get_higtree(sd, 0);
	mp_mapper * m = sd_get_domain_mapper(sd);
	sd_set_interpolator_order(sd, 3);
	
	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		Point ccenter;
		Point cdelta;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];
		hig_get_delta(c, cdelta);
		real dx = cdelta[0];
		real dy = cdelta[1];
		dy *= 1.0;

		/* Define pontos para estêncil de 5 pontos */
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		/* Monta estêncil do Laplaciano discreto */
		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));  // Lado direito
		
		/* Coeficientes do esquema de diferenças finitas:
		   ∇²u ≈ (u_{i-1,j} + u_{i+1,j} - 2u_{i,j})/dx² + (u_{i,j-1} + u_{i,j+1} - 2u_{i,j})/dy² */
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);  // Ponto central
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);  // Vizinho esquerdo
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);  // Vizinho direito
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);  // Vizinho superior
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);  // Vizinho inferior

		/* Adiciona linha na matriz do sistema */
		int *ids = stn_get_ids(stn);
		real *vals = stn_get_vals(stn);
		int numelems = stn_get_numelems(stn);

		int cgid = mp_lookup(m, hig_get_cid(c));
		slv_set_Ai(slv, cgid, numelems, ids, vals);
		slv_set_bi(slv, cgid, stn_get_rhs(stn));
	}
	higcit_destroy(it);
	DEBUG_DIFF_TIME;
	
	/* Montagem e solução do sistema linear */
	printf("assembling...\n");
	slv_assemble(slv);
	DEBUG_DIFF_TIME;

	printf("solving...\n");
	slv_solve(slv);
	DEBUG_DIFF_TIME;

	/* Cálculo de erros */
	real norm_inf = -1.0;
	real norm_2 = 0.0;

	for(it = higcit_create_all_leaves(root); !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		int cgid = mp_lookup(m, hig_get_cid(c));
		Point ccenter;
		hig_get_center(c, ccenter);
		real cx = ccenter[0];
		real cy = ccenter[1];

		/* Erro absoluto no ponto */
		real error = fabs(f(cx, cy) - slv_get_xi(slv, cgid));
		norm_2 += error*error;
		norm_inf = ((norm_inf<error) ? error : norm_inf);
	}
	norm_2 /= (real) size;
	norm_2 = sqrt(norm_2);

	/* Exibe resultados */
	printf("norm_inf = %f\n", norm_inf);
	printf("norm_2 = %f\n", norm_2);
	printf("{%e, %e, %e}\n", 2.0/(side*pow(2.0, refinementLevel)), norm_inf, norm_2);
	higcit_destroy(it);
	return 0;
}

/* COMENTÁRIOS SOBRE ESQUEMAS ALTERNATIVOS DE DIFERENÇAS FINITAS:

O código inclui várias alternativas comentadas para esquemas de diferenças finitas:

1. Esquema de alta ordem com 3 pontos em cada direção
2. Esquema com passo diferente em y (assimétrico)
3. Esquema padrão de 5 pontos (utilizado)
4. Esquema com passo duplo em y
5. Outra variante com passo duplo

Estes podem ser ativados substituindo o bloco atual do estêncil.
*/
/*
		Point center = {cx         , cy     };
		Point px1    = {cx - 1.0*dx, cy};
		Point px2    = {cx + 0.5*dx, cy};
		Point px3    = {cx + 2.0*dx, cy};
		Point py1    = {cx,          cy - 1.0*dy};
		Point py2    = {cx,          cy + 0.5*dy};
		Point py3    = {cx,          cy + 2.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		real wc = -27.0/9.0, w1 = 10.0/9.0, w2 = 16.0/9.0, w3 = 1.0/9.0;

		sd_get_stencil(sd, center, center,  (wc/(dx*dx)), stn);
		sd_get_stencil(sd, center, px1,     (w1/(dx*dx)),  stn);
		sd_get_stencil(sd, center, px2,     (w2/(dx*dx)),  stn);
		sd_get_stencil(sd, center, px3,     (w3/(dx*dx)), stn);
		sd_get_stencil(sd, center, center,  (wc/(dy*dy)), stn);
		sd_get_stencil(sd, center, py1,     (w1/(dy*dy)),  stn);
		sd_get_stencil(sd, center, py2,     (w2/(dy*dy)),  stn);
		sd_get_stencil(sd, center, py3,     (w3/(dy*dy)), stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point p1    = {cx,       cy - 1.0*dy};
		Point p2    = {cx,       cy + 2.0*dy};
		Point p3    = {cx,       cy + 3.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  ((-8.0/6.0)/(dy*dy)), stn);
		sd_get_stencil(sd, center, p1,  ((5.0/6.0)/(dy*dy)),  stn);
		sd_get_stencil(sd, center, p2,  ((4.0/6.0)/(dy*dy)),  stn);
		sd_get_stencil(sd, center, p3,  ((-1.0/6.0)/(dy*dy)),  stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center, (-2.0/(dx*dx)-2.0/(dy*dy)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),              stn);
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),              stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),              stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + dy};
		Point bottom = {cx,      cy - dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  (-2.0/(dy*dy)), stn);
		sd_get_stencil(sd, center, top,     (1.0/(dy*dy)),  stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(dy*dy)),  stn);
*/

/*
		Point center = {cx     , cy     };
		Point left   = {cx - dx, cy     };
		Point right  = {cx + dx, cy     };
		Point top    = {cx,      cy + 2.0*dy};
		Point bottom = {cx,      cy - 2.0*dy};

		stn_reset(stn);
		stn_set_rhs(stn, g(cx, cy));
		sd_get_stencil(sd, center, center,  (-2.0/(dx*dx)), stn);
		sd_get_stencil(sd, center, left,    (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, right,   (1.0/(dx*dx)),  stn);
		sd_get_stencil(sd, center, center,  (-2.0/(4.0*dy*dy)), stn);
		sd_get_stencil(sd, center, top,     (1.0/(4.0*dy*dy)),  stn);
		sd_get_stencil(sd, center, bottom,  (1.0/(4.0*dy*dy)),  stn);
*/
