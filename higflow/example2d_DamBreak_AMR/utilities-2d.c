// *******************************************************************
// *******************************************************************
//  Utility functions for directory - version 03/2023
// *******************************************************************
// *******************************************************************

#include "utilities-2d.h"

// *******************************************************************

// *******************************************************************

dimensional_unit sum_units(dimensional_unit a, dimensional_unit b) {
    dimensional_unit c;
    c.kg = a.kg + b.kg;
    c.m = a.m + b.m;
    c.s = a.s + b.s;
    c.K = a.K + b.K;
    c.mol = a.mol + b.mol;
    c.A = a.A + b.A;
    c.cd = a.cd + b.cd;
    return c;
}

dimensional_unit sub_units(dimensional_unit a, dimensional_unit b) {
    dimensional_unit c;
    c.kg = a.kg - b.kg;
    c.m = a.m - b.m;
    c.s = a.s - b.s;
    c.K = a.K - b.K;
    c.mol = a.mol - b.mol;
    c.A = a.A - b.A;
    c.cd = a.cd - b.cd;
    return c;
}

bool are_units_equal(dimensional_unit a, dimensional_unit b) {
    return (a.kg == b.kg) && (a.m == b.m) && (a.s == b.s) && (a.K == b.K) && (a.mol == b.mol) && (a.A == b.A) && (a.cd == b.cd);
}

physical_quantity mult_pq(physical_quantity a, physical_quantity b) {
    physical_quantity c;
    c.val = a.val * b.val;
    c.unit = sum_units(a.unit, b.unit);
    return c;
}

physical_quantity mult_scalar_pq(real scalar, physical_quantity a) {
    physical_quantity c;
    c.val = scalar * a.val;
    c.unit = a.unit;
    return c;
}

physical_quantity div_pq(physical_quantity a, physical_quantity b) {
    if(b.val != 0.0) {
        physical_quantity c;
        c.val = a.val / b.val;
        c.unit = sub_units(a.unit, b.unit);
        return c;
    }
    printf("Error: division by zero in div_pq");
    exit(1);
}

physical_quantity sum_pq(physical_quantity a, physical_quantity b) {
    if(are_units_equal(a.unit, b.unit)){ 
        physical_quantity c;
        c.unit = a.unit;
        c.val = a.val + b.val;
        return c;
    }
    printf("Error: units are not equal in sum_pq");
    exit(1);
}

// subtract physical quantities
physical_quantity sub_pq(physical_quantity a, physical_quantity b) {
    if(are_units_equal(a.unit, b.unit)){ 
        physical_quantity c;
        c.unit = a.unit;
        c.val = a.val - b.val;
        return c;
    }
    printf("Error: units are not equal in sub_pq");
    exit(1);
}

physical_quantity sqrt_pq(physical_quantity a) {
    
    if(a.unit.kg % 2 == 0 && a.unit.m % 2 == 0 && a.unit.s % 2 == 0 && a.unit.K % 2 == 0 && a.unit.mol % 2 == 0 && a.unit.A % 2 == 0 && a.unit.cd % 2 == 0) {
        physical_quantity c;
        c.val = sqrt(a.val);
        c.unit = (dimensional_unit) {a.unit.kg/2, a.unit.m/2, a.unit.s/2, a.unit.K/2, a.unit.mol/2, a.unit.A/2, a.unit.cd/2};
        return c;
    }
    printf("Error: sqrt_pq: units are not even");
    exit(1);
}

void print_pq(physical_quantity a) {
    if(a.val >= 1.0 && a.val < 1.0e4) {
        printf("%7.2f", a.val);
    }
    else if(a.val < 1.0 && a.val >= 1.0e-3) {
        printf("%7.5f", a.val);
    }
    else printf("%.4e", a.val);

    if(a.unit.kg != 0) {
        if(a.unit.kg != 1) printf(" kg^%d", a.unit.kg);
        else printf(" kg");
    }
    if(a.unit.m != 0) {
        if(a.unit.m != 1) printf(" m^%d", a.unit.m);
        else printf(" m");
    }
    if(a.unit.s != 0) {
        if(a.unit.s != 1) printf(" s^%d", a.unit.s);
        else printf(" s");
    }
    if(a.unit.K != 0) {
        if(a.unit.K != 1) printf(" K^%d", a.unit.K);
        else printf(" K");
    }
    if(a.unit.mol != 0) {
        if(a.unit.mol != 1) printf(" mol^%d", a.unit.mol);
        else printf(" mol");
    }
    if(a.unit.A != 0) {
        if(a.unit.A != 1) printf(" A^%d", a.unit.A);
        else printf(" A");
    }
    if(a.unit.cd != 0) {
        if(a.unit.cd != 1) printf(" cd^%d", a.unit.cd);
        else printf(" cd");
    }
}


// this function returns the dimensions of the domain using the mesh file
bool get_lengths(physical_parameters *p_par, char *nameload) {
    bool is_domain_cuboid = true;
    char namefile[1024];
    sprintf(namefile,"%s.domain", nameload);
    FILE *fdomain = fopen(namefile, "r");
    if (fdomain == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
    // Number of HigTrees
    int numhigs;
    int ifd = fscanf(fdomain,"%d\n",&numhigs);
    if(numhigs!=1) {
        printf("Number of HigTrees: %d ==============> EXPECTED TO BE 1 - domain may not be cuboid\n", numhigs);
        is_domain_cuboid = false;
    }
    char amrfilename[1024];
    __higflow_readstring(amrfilename,1024,fdomain);
    // Open the AMR format file
    FILE *fd = fopen(amrfilename, "r");
    if (fd == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",amrfilename);
        exit(1);
    }
    // low and high points of the domain
    real l[DIM], h[DIM];
    // Reading the higtree information from the file 
    for(int dim = 0; dim < DIM; dim++) {
		int status = fscanf(fd, "%lf", &l[dim]);
		status = fscanf(fd, "%lf", &h[dim]);
        p_par->center[dim].val = (h[dim] + l[dim])/2.0;
        p_par->center[dim].unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->L[dim].val = h[dim] - l[dim];
        p_par->L[dim].unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->L_physical[dim] = mult_pq(p_par->L[dim], p_par->H);
	}
    // Close the AMR format file
    fclose(fd);

    fclose(fdomain);
    return is_domain_cuboid;
}

bool get_lengths_yaml(physical_parameters *p_par, char *nameload) {
    bool is_domain_cuboid = true;
    char namefile[1024];
    sprintf(namefile,"%s.domain.yaml", nameload);

    struct fy_document *fyd = fy_document_build_from_file(NULL, namefile);
    if (fyd == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }
    // Number of HigTrees
    int numhigs;
    int ifd = fy_document_scanf(fyd,"domain/number_domains %d",&numhigs);
    if(numhigs!=1) {
        printf("Number of HigTrees: %d ==============> EXPECTED TO BE 1 - domain may not be cuboid\n", numhigs);
        is_domain_cuboid = false;
    }
    char amrfilename[1024];
    ifd = fy_document_scanf(fyd,"domain/domain0/path %s",amrfilename);
    fy_document_destroy(fyd);
    // Open the AMR format file
    FILE *fd = fopen(amrfilename, "r");
    if (fd == NULL) {
        // Error in open the file
        printf("=+=+=+= Error loading file %s =+=+=+=\n",amrfilename);
        exit(1);
    }
    // low and high points of the domain
    real l[DIM], h[DIM];
    // Reading the higtree information from the file 
    for(int dim = 0; dim < DIM; dim++) {
		int status = fscanf(fd, "%lf", &l[dim]);
		status = fscanf(fd, "%lf", &h[dim]);
        p_par->center[dim].val = (h[dim] + l[dim])/2.0;
        p_par->center[dim].unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->L[dim].val = h[dim] - l[dim];
        p_par->L[dim].unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->L_physical[dim] = mult_pq(p_par->L[dim], p_par->H);
	}
    // Close the AMR format file
    fclose(fd);
    return is_domain_cuboid;
}


void initialize_physical_parameters(physical_parameters *p_par, higflow_controllers hig_contr, higflow_parameters hig_par, 
ve_controllers ve_contr, ve_parameters ve_par, eo_controllers eo_contr, 
eo_parameters eo_par, mult_controllers mult_contr, mult_surf_parameters mult_spar, 
mult_parameters mult_par0, mult_parameters mult_par1, 
mult_ve_controllers mult_ve_contr, ve_parameters mult_ve_par0, 
ve_parameters mult_ve_par1, eo_controllers mult_eo_contr,
eo_parameters mult_eo_par0, eo_parameters mult_eo_par1, 
int myrank) {

    print0f("\n=+=+=+=+=+=+=+=+=+=+=+=+ Physical parameters =+=+=+=+=+=+=+=+=+=+=+=+\n");

    // density of water
    p_par->rho_ref.val = 1000.0;
    p_par->rho_ref.unit = (dimensional_unit) {1, -3, 0, 0, 0, 0, 0};
    // reference length of 1 meter - dam break scale
    p_par->H.val = 1.0;
    p_par->H.unit = (dimensional_unit) {0, 1, 0, 0, 0, 0, 0};
    // approximate viscosity of water
    p_par->mu_ref.val = 1.0e-3;
    p_par->mu_ref.unit = (dimensional_unit) {1, -1, -1, 0, 0, 0, 0};
    // dynamic viscosity
    p_par->nu = div_pq(p_par->mu_ref, p_par->rho_ref);
    // obtain reynolds number from higflow parameters
    p_par->Re.val = hig_par.Re;
    p_par->Re.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
    // compute characteristic velocity from reynolds number
    p_par->U = div_pq (p_par->Re , div_pq ( mult_pq(p_par->rho_ref, p_par->H), p_par->mu_ref ) );
    // get time step
    p_par->dt.val = hig_par.dt;
    p_par->dt.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
    // compute characteristic time
    p_par->tau_time.val = p_par->H.val / p_par->U.val;
    p_par->tau_time.unit = (dimensional_unit) {0, 0, 1, 0, 0, 0, 0};
    // real time step
    p_par->dt_physical = mult_pq(p_par->dt, p_par->tau_time);
    // gravitational acceleration
    p_par->g.val = 9.80665;
    p_par->g.unit = (dimensional_unit) {0, 1, -2, 0, 0, 0, 0};
    if(hig_contr.add_gravity ==  true) {
        // compute Froude number
        p_par->Fr.val = hig_par.Fr;
        p_par->Fr.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
    }

    bool is_domain_cuboid = get_lengths_yaml(p_par, hig_par.nameload);
    if(myrank == 0) {
        if(is_domain_cuboid == true){
            printf("+=+=+= Cuboid Domain Dimensions: +=+=+=\n");
        // print domain dimensions
            for(int dim = 0; dim < DIM; dim++) {
                printf("L_%d = ", dim);
                print_pq(p_par->L[dim]);
                printf("\n");
                printf("L_physical_%d = ", dim);
                print_pq(p_par->L_physical[dim]);
                printf("\n");
                printf("center_%d = ", dim);
                print_pq(p_par->center[dim]);
                printf("\n");
            }
        }
    }

    if(myrank == 0) {
        //print all physical parameters yet defined
        printf("rho_ref = ");
        print_pq(p_par->rho_ref);
        printf("\n");
        printf("mu_ref = ");
        print_pq(p_par->mu_ref);
        printf("\n");
        printf("nu = ");
        print_pq(p_par->nu);
        printf("\n");
        printf("Re = ");
        print_pq(p_par->Re);
        printf("\n");
        printf("U = ");
        print_pq(p_par->U);
        printf("\n");
        printf("dt = ");
        print_pq(p_par->dt);
        printf("\n");
        printf("tau_time = ");
        print_pq(p_par->tau_time);
        printf("\n");
        printf("dt_physical = ");
        print_pq(p_par->dt_physical);
        printf("\n");
        printf("g = ");
        print_pq(p_par->g);
        printf("\n");
        if(hig_contr.add_gravity ==  true) {
            printf("Fr = ");
            print_pq(p_par->Fr);
            printf("\n");
        }
    }

    if(hig_contr.eoflow == true) {
        print0f("\n+=+=+= Electroosmotic parameters +=+=+=\n");
        // elementary charge of the electron = 1.602176634e-19 C = 9.64853321233100184 C/mol
        p_par->e.val = 96485.3321233100184;
        p_par->e.unit = (dimensional_unit) {0, 0, 1, 0, -1, 1, 0};
        // assuming monovalent ionic solute
        p_par->Z.val = 1.0;
        p_par->Z.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // permittivity of free space
        p_par->epsilon_0.val = 8.8541878128e-12;
        p_par->epsilon_0.unit = (dimensional_unit) {-1, -3, 4, 0, 0, 2, 0};
        // relative permittivity of water around normal temperatures
        p_par->epsilon_r.val = 80.1;
        p_par->epsilon_r.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // permittivity of water around normal temperatures
        p_par->epsilon_e = mult_pq(p_par->epsilon_0, p_par->epsilon_r);
        if(eo_contr.eo_model == PNP) {
            // obtain Péclet number
            p_par->Pe.val = eo_par.Pe;
            p_par->Pe.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // compute diffusivity using Péclet number
            p_par->D = div_pq (mult_pq (p_par->Re, p_par->nu), p_par->Pe);
        }
        // boltzmann constant = 1.380649e-23 J/K  =  8.31446261815324 J/(mol K)
        p_par->k_B.val = 8.31446261815324;
        p_par->k_B.unit = (dimensional_unit) {1, 2, -2, -1, -1, 0, 0};
        // get alpha
        p_par->alpha_eo.val = eo_par.alpha;
        p_par->alpha_eo.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // temperature
        p_par->T.val = 298;
        p_par->T.unit = (dimensional_unit) {0, 0, 0, 1, 0, 0, 0};
        // reference potential calculated from alpha = zeta_ref e Z / (k_B T)
        p_par->zeta_ref = div_pq(mult_pq(p_par->alpha_eo, mult_pq(p_par->k_B, p_par->T)), mult_pq(p_par->e, p_par->Z));
        // get delta
        p_par->delta_eo.val = eo_par.delta;
        p_par->delta_eo.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // determine reference concentration from delta = n_0 H^2 e z/ (ep_e zeta_ref)
        p_par->n_ref = div_pq(mult_pq(p_par->delta_eo, mult_pq(p_par->epsilon_e, p_par->zeta_ref)), mult_pq(mult_pq(p_par->H, p_par->H), mult_pq(p_par->e, p_par->Z)) );
        // get Debye number kappa = sqrt(2 alpha delta)
        p_par->kappa_eo = sqrt_pq( mult_scalar_pq(2.0, mult_pq(p_par->alpha_eo, p_par->delta_eo)) );
        // Debye length
        p_par->lambda_D = div_pq(p_par->H, p_par->kappa_eo);
        // get potential difference for every unit of reference length
        p_par->Ex.val = eo_par.Ex;
        p_par->Ex.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // compute real potential difference per unit of length
        p_par->Ex_physical = div_pq(mult_pq(p_par->Ex, p_par->zeta_ref), p_par->H);
        // get physical Helmholtz-Smoluchowski velocity (ep_e zeta_ref / mu) * E
        p_par->u_hs_physical = mult_pq (div_pq(mult_pq(p_par->zeta_ref, p_par->epsilon_e), p_par->mu_ref), p_par->Ex_physical);
        // compute conversion factor between eletro-osmotic and momentum systems
        physical_quantity one = (physical_quantity) {1.0, (dimensional_unit) {0, 0, 0, 0, 0, 0, 0}};
        p_par->G_x = div_pq(one, mult_pq(p_par->Ex, p_par->Re));
        
        if(myrank == 0) {
            // print all electroosmotic physical parameters
            printf("e = "); 
            print_pq(p_par->e);
            printf("\n");
            printf("Z = ");
            print_pq(p_par->Z);
            printf("\n");
            printf("epsilon_0 = ");
            print_pq(p_par->epsilon_0);
            printf("\n");
            printf("epsilon_r = ");
            print_pq(p_par->epsilon_r);
            printf("\n");
            printf("epsilon_e = ");
            print_pq(p_par->epsilon_e);
            printf("\n");
            printf("zeta_ref = ");
            print_pq(p_par->zeta_ref);
            printf("\n");
            if(eo_contr.eo_model == PNP) {
                printf("Pe = ");
                print_pq(p_par->Pe);
                printf("\n");
                printf("D = ");
                print_pq(p_par->D);
                printf("\n");
            }
            printf("k_B = ");
            print_pq(p_par->k_B);
            printf("\n");
            printf("alpha_eo = ");
            print_pq(p_par->alpha_eo);
            printf("\n");
            printf("T = ");
            print_pq(p_par->T);
            printf("\n");
            printf("delta_eo = ");
            print_pq(p_par->delta_eo);
            printf("\n");
            printf("n_0 = ");
            print_pq(p_par->n_ref);
            printf("\n");
            printf("kappa_eo = ");
            print_pq(p_par->kappa_eo);
            printf("\n");
            printf("lambda_D = ");
            print_pq(p_par->lambda_D);
            printf("\n");
            printf("Ex = ");
            print_pq(p_par->Ex);
            printf("\n");
            printf("Ex_physical = ");
            print_pq(p_par->Ex_physical);
            printf("\n");
            printf("u_hs_physical = ");
            print_pq(p_par->u_hs_physical);
            printf("\n");
            printf("G_x = ");
            print_pq(p_par->G_x);
            printf("\n");
            
        }
    }

    if(hig_contr.flowtype == VISCOELASTIC) { //viscoelastic flow
        print0f("\n+=+=+= Viscoelastic parameters +=+=+=\n");
        // get Deborah number
        p_par->De.val = ve_par.De;
        p_par->De.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // find relaxation time using Deborah number (De = lambda U / H)
        p_par->lambda_time = mult_pq(p_par->De, div_pq(p_par->H, p_par->U));
        // get beta
        p_par->beta.val = ve_par.beta;
        p_par->beta.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // find viscosity of solvent using beta
        p_par->mu_s = mult_pq(p_par->mu_ref, p_par->beta);
        // find polymeric viscosity (mu = mu_p + mu_s)
        p_par->mu_p = sub_pq(p_par->mu_ref, p_par->mu_s);
        // find spring constant (lambda = K_s/mu_p)
        p_par->K_s = div_pq(p_par->mu_p, p_par->lambda_time);
        // get Giesekus alpha
        p_par->alpha_giesekus.val = ve_par.alpha;
        p_par->alpha_giesekus.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get PTT epsilon
        p_par->epsilon_ptt.val = ve_par.epsilon;
        p_par->epsilon_ptt.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get PTT xi
        p_par->xi_ptt.val = ve_par.xi;
        p_par->xi_ptt.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get gptt alpha
        p_par->alpha_gptt.val = ve_par.alpha;
        p_par->alpha_gptt.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get gptt beta
        p_par->beta_gptt.val = ve_par.beta;
        p_par->beta_gptt.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get fene L^2
        p_par->L2_fene.val = ve_par.L2_fene;
        p_par->L2_fene.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get e-fene lambda
        p_par->lambda_fene.val = ve_par.lambda_fene;
        p_par->lambda_fene.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // get e-fene E
        p_par->E_fene.val = ve_par.E_fene;
        p_par->E_fene.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        // compute steady state fully developed newtonian velocity normalized by flow rate
        p_par->Un_U.val = solve_un_u(ve_contr.model, p_par);
        p_par->Un_U.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};

        if(myrank == 0) {
            // print all physical parameters
            printf("De = ");
            print_pq(p_par->De);
            printf("\n");
            printf("lambda_time = ");
            print_pq(p_par->lambda_time);
            printf("\n");
            printf("beta = ");
            print_pq(p_par->beta);
            printf("\n");
            printf("mu_s = ");
            print_pq(p_par->mu_s);
            printf("\n");
            printf("mu_p = ");
            print_pq(p_par->mu_p);
            printf("\n");
            printf("K_s = ");
            print_pq(p_par->K_s);
            printf("\n");
            printf("alpha_giesekus = ");
            print_pq(p_par->alpha_giesekus);
            printf("\n");
            printf("epsilon_ptt = ");
            print_pq(p_par->epsilon_ptt);
            printf("\n");
            printf("xi_ptt = ");
            print_pq(p_par->xi_ptt);
            printf("\n");
            printf("alpha_gptt = ");
            print_pq(p_par->alpha_gptt);
            printf("\n");
            printf("beta_gptt = ");
            print_pq(p_par->beta_gptt);
            printf("\n");
            printf("L2_fene = ");
            print_pq(p_par->L2_fene);
            printf("\n");
            printf("lambda_fene = ");
            print_pq(p_par->lambda_fene);
            printf("\n");
            printf("E_fene = ");
            print_pq(p_par->E_fene);
            printf("\n");
            printf("Un_U = ");
            print_pq(p_par->Un_U);
            printf("\n");
        }
    }

    if(hig_contr.flowtype == MULTIPHASE) { // multiphase
        if(mult_contr.add_surface_tension == true) {
            // Capillary number
            p_par->Ca.val = mult_spar.Ca;
            p_par->Ca.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // interfacial tension coefficient
            p_par->sigma = div_pq(mult_pq(p_par->mu_ref,p_par->U), p_par->Ca);
            // Bond number
            p_par->Bo = div_pq(mult_pq(mult_pq(p_par->rho_ref,p_par->g), mult_pq(p_par->H, p_par->H)), p_par->sigma);
            // Weber number
            p_par->We = mult_pq(p_par->Ca, p_par->Re);

            if(myrank == 0) {
                printf("\n+=+=+= Multiphase parameters +=+=+=\n");
                // print all physical parameters
                printf("sigma = ");
                print_pq(p_par->sigma);
                printf("\n");
                printf("Bo = ");
                print_pq(p_par->Bo);
                printf("\n");
                printf("Ca = ");
                print_pq(p_par->Ca);
                printf("\n");
                printf("We = ");
                print_pq(p_par->We);
                printf("\n");
            }

        }
        // reference density given by user in phase 0
        p_par->rho0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->rho0.val = mult_par0.rho;
        // reference density given by user in phase 1
        p_par->rho1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->rho1.val = mult_par1.rho;
        // reference viscosity given by user in phase 0
        p_par->mu0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->mu0.val = mult_par0.mu;
        // reference viscosity given by user in phase 1
        p_par->mu1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
        p_par->mu1.val = mult_par1.mu;

        if(myrank == 0) {
            printf("\n+=+=+= Multiphase parameters - Phase 0 +=+=+=\n");
            printf("rho0 = ");
            print_pq(p_par->rho0);
            printf("\n");
            printf("mu0 = ");
            print_pq(p_par->mu0);
            printf("\n");
            printf("\n+=+=+= Multiphase parameters - Phase 1 +=+=+=\n");
            printf("rho1 = ");
            print_pq(p_par->rho1);
            printf("\n");
            printf("mu1 = ");
            print_pq(p_par->mu1);
            printf("\n");
        }

        if(mult_contr.viscoelastic_either == true)  {
            // get Deborah number
            p_par->De0.val = mult_ve_par0.De;
            p_par->De0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->De1.val = mult_ve_par1.De;
            p_par->De1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // find relaxation time using Deborah number (De = lambda U / H)
            p_par->lambda_time0 = mult_pq(p_par->De0, div_pq(p_par->H, p_par->U));
            p_par->lambda_time1 = mult_pq(p_par->De1, div_pq(p_par->H, p_par->U));
            // get beta
            p_par->beta0.val = mult_ve_par0.beta;
            p_par->beta0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->beta1.val = mult_ve_par1.beta;
            p_par->beta1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // find viscosity of solvent using beta
            p_par->mu_s0 = mult_pq(p_par->mu_ref, p_par->beta0);
            p_par->mu_s1 = mult_pq(p_par->mu_ref, p_par->beta1);
            // find polymeric viscosity (mu = mu_p + mu_s)
            p_par->mu_p0 = sub_pq(p_par->mu_ref, p_par->mu_s0);
            p_par->mu_p1 = sub_pq(p_par->mu_ref, p_par->mu_s1);
            // find spring constant (lambda = K_s/mu_p)
            if(FLT_EQ(p_par->lambda_time0.val,0.0)){
                p_par->K_s0.unit = (dimensional_unit) {1, -1, -2, 0, 0, 0, 0};
                p_par->K_s0.val = 0.0;
            }
            else p_par->K_s0 = div_pq(p_par->mu_p0, p_par->lambda_time0);
            if(FLT_EQ(p_par->lambda_time1.val,0.0)){
                p_par->K_s1.unit = (dimensional_unit) {1, -1, -2, 0, 0, 0, 0};
                p_par->K_s1.val = 0.0;
            }
            else p_par->K_s1 = div_pq(p_par->mu_p1, p_par->lambda_time1);
            // get Giesekus alpha
            p_par->alpha_giesekus0.val = mult_ve_par0.alpha;
            p_par->alpha_giesekus0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->alpha_giesekus1.val = mult_ve_par1.alpha;
            p_par->alpha_giesekus1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get PTT epsilon
            p_par->epsilon_ptt0.val = mult_ve_par0.epsilon;
            p_par->epsilon_ptt0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->epsilon_ptt1.val = mult_ve_par1.epsilon;
            p_par->epsilon_ptt1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get PTT xi
            p_par->xi_ptt0.val = mult_ve_par0.xi;
            p_par->xi_ptt0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->xi_ptt1.val = mult_ve_par1.xi;
            p_par->xi_ptt1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get gptt alpha
            p_par->alpha_gptt0.val = mult_ve_par0.alpha;
            p_par->alpha_gptt0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->alpha_gptt1.val = mult_ve_par1.alpha;
            p_par->alpha_gptt1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get gptt beta
            p_par->beta_gptt0.val = mult_ve_par0.beta;
            p_par->beta_gptt0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->beta_gptt1.val = mult_ve_par1.beta;
            p_par->beta_gptt1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get fene L^2
            p_par->L2_fene0.val = mult_ve_par0.L2_fene;
            p_par->L2_fene0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->L2_fene1.val = mult_ve_par1.L2_fene;
            p_par->L2_fene1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get e-fene lambda
            p_par->lambda_fene0.val = mult_ve_par0.lambda_fene;
            p_par->lambda_fene0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->lambda_fene1.val = mult_ve_par1.lambda_fene;
            p_par->lambda_fene1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // get e-fene E
            p_par->E_fene0.val = mult_ve_par0.E_fene;
            p_par->E_fene0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->E_fene1.val = mult_ve_par1.E_fene;
            p_par->E_fene1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            
            if(myrank == 0) {

                // print all physical parameters
                printf("\n+=+= Multiphase viscoelastic parameters - Phase 0 +=+=\n");
                if(mult_contr.flowtype0 != VISCOELASTIC)
                    printf("Warning: Since phase 0 is not viscoelastic, these parameters are only a means of interpolating in the simulation\n");
                printf("De = ");
                print_pq(p_par->De0);
                printf("\n");
                printf("lambda_time = ");
                print_pq(p_par->lambda_time0);
                printf("\n");
                printf("beta = ");
                print_pq(p_par->beta0);
                printf("\n");
                printf("mu_s = ");
                print_pq(p_par->mu_s0);
                printf("\n");
                printf("mu_p = ");
                print_pq(p_par->mu_p0);
                printf("\n");
                printf("K_s = ");
                print_pq(p_par->K_s0);
                printf("\n");
                printf("alpha_giesekus = ");
                print_pq(p_par->alpha_giesekus0);
                printf("\n");
                printf("epsilon_ptt = ");
                print_pq(p_par->epsilon_ptt0);
                printf("\n");
                printf("xi_ptt = ");
                print_pq(p_par->xi_ptt0);
                printf("\n");
                printf("alpha_gptt = ");
                print_pq(p_par->alpha_gptt0);
                printf("\n");
                printf("beta_gptt = ");
                print_pq(p_par->beta_gptt0);
                printf("\n");
                printf("L2_fene = ");
                print_pq(p_par->L2_fene0);
                printf("\n");
                printf("lambda_fene = ");
                print_pq(p_par->lambda_fene0);
                printf("\n");
                printf("E_fene = ");
                print_pq(p_par->E_fene0);
                printf("\n");

                printf("\n+=+= Multiphase viscoelastic parameters - Phase 1 +=+=\n");
                if(mult_contr.flowtype1 != VISCOELASTIC)
                    printf("Warning: Since phase 1 is not viscoelastic, these parameters are only a means of interpolating in the simulation\n");
                printf("Phase 1\n");
                printf("De = ");
                print_pq(p_par->De1);
                printf("\n");
                printf("lambda_time = ");
                print_pq(p_par->lambda_time1);
                printf("\n");
                printf("beta = ");
                print_pq(p_par->beta1);
                printf("\n");
                printf("mu_s = ");
                print_pq(p_par->mu_s1);
                printf("\n");
                printf("mu_p = ");
                print_pq(p_par->mu_p1);
                printf("\n");
                printf("K_s = ");
                print_pq(p_par->K_s1);
                printf("\n");
                printf("alpha_giesekus = ");
                print_pq(p_par->alpha_giesekus1);
                printf("\n");
                printf("epsilon_ptt = ");
                print_pq(p_par->epsilon_ptt1);
                printf("\n");
                printf("xi_ptt = ");
                print_pq(p_par->xi_ptt1);
                printf("\n");
                printf("alpha_gptt = ");
                print_pq(p_par->alpha_gptt1);
                printf("\n");
                printf("beta_gptt = ");
                print_pq(p_par->beta_gptt1);
                printf("\n");
                printf("L2_fene = ");
                print_pq(p_par->L2_fene1);
                printf("\n");
                printf("lambda_fene = ");
                print_pq(p_par->lambda_fene1);
                printf("\n");
                printf("E_fene = ");
                print_pq(p_par->E_fene1);
                printf("\n");
            }
        }


        if(mult_contr.eoflow_either ==  true) {
            print0f("\n+=+=+= Multiphase Electroosmotic parameters +=+=+=\n");
            // elementary charge of the electron = 1.602176634e-19 C = 9.64853321233100184 C/mol
            p_par->e.val = 96485.3321233100184;
            p_par->e.unit = (dimensional_unit) {0, 0, 1, 0, -1, 1, 0};
            // assuming monovalent ionic solute
            p_par->Z.val = 1.0;
            p_par->Z.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // permittivity of free space
            p_par->epsilon_0.val = 8.8541878128e-12;
            p_par->epsilon_0.unit = (dimensional_unit) {-1, -3, 4, 0, 0, 2, 0};
            // relative permittivity of water around normal temperatures
            p_par->epsilon_r.val = 80.1;
            p_par->epsilon_r.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // reference relative adimensional permittivity given by the user in phase 0
            p_par->perm0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->perm0.val = mult_eo_par0.perm;
            p_par->epsilon_r0 = mult_pq(p_par->epsilon_r, p_par->perm0);
            // reference relative adimensional permittivity given by the user in phase 1
            p_par->perm1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->perm1.val = mult_eo_par1.perm;
            p_par->epsilon_r1 = mult_pq(p_par->epsilon_r, p_par->perm1);
            // permittivity of water around normal temperatures
            p_par->epsilon_e = mult_pq(p_par->epsilon_0, p_par->epsilon_r);
            if(mult_eo_contr.eo_model == PNP) {
                // obtain Péclet number
                p_par->Pe0.val = mult_eo_par0.Pe;
                p_par->Pe0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
                p_par->Pe1.val = mult_eo_par1.Pe;
                p_par->Pe1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
                // compute diffusivity using Péclet number
                if(FLT_EQ(p_par->Pe0.val,0.0)){
                    p_par->D0.unit = (dimensional_unit) {0, 2, 1, 0, 0, 0, 0};
                    p_par->D0.val = 0.0;
                }
                else p_par->D0 = div_pq (mult_pq (p_par->Re, p_par->nu), p_par->Pe0);
                if(FLT_EQ(p_par->Pe1.val,0.0)){
                p_par->D1.unit = (dimensional_unit) {0, 2, 1, 0, 0, 0, 0};
                p_par->D1.val = 0.0;
                }
                else p_par->D1 = div_pq (mult_pq (p_par->Re, p_par->nu), p_par->Pe1);
            }
            // boltzmann constant = 1.380649e-23 J/K  =  8.31446261815324 J/(mol K)
            p_par->k_B.val = 8.31446261815324;
            p_par->k_B.unit = (dimensional_unit) {1, 2, -2, -1, -1, 0, 0};
            // get alpha
            p_par->alpha_eo0.val = mult_eo_par0.alpha;
            p_par->alpha_eo0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->alpha_eo1.val = mult_eo_par1.alpha;
            p_par->alpha_eo1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // temperature
            p_par->T.val = 298;
            p_par->T.unit = (dimensional_unit) {0, 0, 0, 1, 0, 0, 0};
            // reference potential calculated from alpha = zeta_ref e Z / (k_B T)
            p_par->zeta_ref0 = div_pq(mult_pq(p_par->alpha_eo0, mult_pq(p_par->k_B, p_par->T)), mult_pq(p_par->e, p_par->Z));
            p_par->zeta_ref1 = div_pq(mult_pq(p_par->alpha_eo1, mult_pq(p_par->k_B, p_par->T)), mult_pq(p_par->e, p_par->Z));
            // get delta
            p_par->delta_eo0.val = mult_eo_par0.delta;
            p_par->delta_eo0.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->delta_eo1.val = mult_eo_par1.delta;
            p_par->delta_eo1.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // determine reference concentration from delta = n_0 H^2 e z/ (ep_e zeta_ref)
            p_par->n_ref0 = div_pq(mult_pq(p_par->delta_eo0, mult_pq(p_par->epsilon_e, p_par->zeta_ref0)), mult_pq(mult_pq(p_par->H, p_par->H), mult_pq(p_par->e, p_par->Z)) );
            p_par->n_ref1 = div_pq(mult_pq(p_par->delta_eo1, mult_pq(p_par->epsilon_e, p_par->zeta_ref1)), mult_pq(mult_pq(p_par->H, p_par->H), mult_pq(p_par->e, p_par->Z)) );
            // get Debye number kappa = sqrt(2 alpha delta)
            p_par->kappa_eo0 = sqrt_pq( mult_scalar_pq(2.0, mult_pq(p_par->alpha_eo0, p_par->delta_eo0)) );
            p_par->kappa_eo1 = sqrt_pq( mult_scalar_pq(2.0, mult_pq(p_par->alpha_eo1, p_par->delta_eo1)) );
            // Debye length
            if(FLT_EQ(p_par->kappa_eo0.val,0.0)){
                p_par->lambda_D0.unit = (dimensional_unit) {0, 1, 0, 0, 0, 0, 0};
                p_par->lambda_D0.val = 0.0;
            }
            else p_par->lambda_D0 = div_pq(p_par->H, p_par->kappa_eo0);
            if(FLT_EQ(p_par->kappa_eo1.val,0.0)){
                p_par->lambda_D1.unit = (dimensional_unit) {0, 1, 0, 0, 0, 0, 0};
                p_par->lambda_D1.val = 0.0;
            }
            else p_par->lambda_D1 = div_pq(p_par->H, p_par->kappa_eo1);
            // get potential difference for every unit of reference length
            p_par->Ex0.val = mult_eo_par0.Ex;
            p_par->Ex.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            p_par->Ex1.val = mult_eo_par1.Ex;
            p_par->Ex.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
            // compute real potential difference per unit of length
            p_par->Ex_physical0 = div_pq(mult_pq(p_par->Ex0, p_par->zeta_ref0), p_par->H);
            p_par->Ex_physical1 = div_pq(mult_pq(p_par->Ex1, p_par->zeta_ref1), p_par->H);
            // get physical Helmholtz-Smoluchowski velocity (ep_e zeta_ref / mu) * E
            p_par->u_hs_physical0 = mult_pq (div_pq(mult_pq(p_par->zeta_ref0, p_par->epsilon_e), p_par->mu_ref), p_par->Ex_physical0);
            p_par->u_hs_physical1 = mult_pq (div_pq(mult_pq(p_par->zeta_ref1, p_par->epsilon_e), p_par->mu_ref), p_par->Ex_physical1);
            // compute conversion factor between eletro-osmotic and momentum systems
            physical_quantity one = (physical_quantity) {1.0, (dimensional_unit) {0, 0, 0, 0, 0, 0, 0}};
            if(FLT_EQ(p_par->Ex0.val,0.0)){
                p_par->G_x.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
                p_par->G_x.val = 0.0;
            }
            else p_par->G_x = div_pq(one, mult_pq(p_par->Ex0, p_par->Re));
            if(FLT_EQ(p_par->Ex1.val,0.0)){
                p_par->G_x.unit = (dimensional_unit) {0, 0, 0, 0, 0, 0, 0};
                p_par->G_x.val = 0.0;
            }
            else p_par->G_x = div_pq(one, mult_pq(p_par->Ex1, p_par->Re));
            
            if(myrank == 0) {
                // print all electroosmotic physical parameters
                printf("\n+=+=+= Multiphase Electroosmotic parameters - both phases +=+=+=\n");
                printf("e = "); 
                print_pq(p_par->e);
                printf("\n");
                printf("Z = ");
                print_pq(p_par->Z);
                printf("\n");
                printf("epsilon_0 = ");
                print_pq(p_par->epsilon_0);
                printf("\n");
                printf("epsilon_r = ");
                print_pq(p_par->epsilon_r);
                printf("\n");
                printf("epsilon_e = ");
                print_pq(p_par->epsilon_e);
                printf("\n");
                printf("k_B = ");
                print_pq(p_par->k_B);
                printf("\n");
                printf("T = ");
                print_pq(p_par->T);
                printf("\n");

                printf("\n+=+=+= Multiphase Electroosmotic parameters - Phase 0 +=+=+=\n");
                printf("zeta_ref = ");
                print_pq(p_par->zeta_ref0);
                printf("\n");
                if(mult_eo_contr.eo_model == PNP) {
                    printf("Pe = ");
                    print_pq(p_par->Pe0);
                    printf("\n");
                    printf("D = ");
                    print_pq(p_par->D0);
                    printf("\n");
                }
                printf("alpha_eo = ");
                print_pq(p_par->alpha_eo0);
                printf("\n");
                printf("delta_eo = ");
                print_pq(p_par->delta_eo0);
                printf("\n");
                printf("n_0 = ");
                print_pq(p_par->n_ref0);
                printf("\n");
                printf("kappa_eo = ");
                print_pq(p_par->kappa_eo0);
                printf("\n");
                printf("lambda_D = ");
                print_pq(p_par->lambda_D0);
                printf("\n");
                printf("Ex = ");
                print_pq(p_par->Ex0);
                printf("\n");
                printf("Ex_physical = ");
                print_pq(p_par->Ex_physical0);
                printf("\n");
                printf("u_hs_physical = ");
                print_pq(p_par->u_hs_physical0);
                printf("\n");
                printf("G_x = ");
                print_pq(p_par->G_x0);
                printf("\n");
                printf("perm = ");
                print_pq(p_par->perm0);
                printf("\n");
                printf("epsilon_r = ");
                print_pq(p_par->epsilon_r0);
                printf("\n");

                printf("\n+=+=+= Multiphase Electroosmotic parameters - Phase 1 +=+=+=\n");
                printf("zeta_ref = ");
                print_pq(p_par->zeta_ref1);
                printf("\n");
                if(mult_eo_contr.eo_model == PNP) {
                    printf("Pe = ");
                    print_pq(p_par->Pe1);
                    printf("\n");
                    printf("D = ");
                    print_pq(p_par->D1);
                    printf("\n");
                }
                printf("alpha_eo = ");
                print_pq(p_par->alpha_eo1);
                printf("\n");
                printf("delta_eo = ");
                print_pq(p_par->delta_eo1);
                printf("\n");
                printf("n_0 = ");
                print_pq(p_par->n_ref1);
                printf("\n");
                printf("kappa_eo = ");
                print_pq(p_par->kappa_eo1);
                printf("\n");
                printf("lambda_D = ");
                print_pq(p_par->lambda_D1);
                printf("\n");
                printf("Ex = ");
                print_pq(p_par->Ex1);
                printf("\n");
                printf("Ex_physical = ");
                print_pq(p_par->Ex_physical1);
                printf("\n");
                printf("u_hs_physical = ");
                print_pq(p_par->u_hs_physical1);
                printf("\n");
                printf("G_x = ");
                print_pq(p_par->G_x1);
                printf("\n");
                printf("perm = ");
                print_pq(p_par->perm1);
                printf("\n");
                printf("epsilon_r = ");
                print_pq(p_par->epsilon_r1);
                printf("\n");
            }
        }
    }

    // end physical parameters
    print0f("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+\n\n");
}

physical_parameters *create_initialize_physical_parameters(higflow_solver *ns, int myrank) {
    higflow_controllers hig_contr = ns->contr;
    higflow_parameters hig_par = ns->par;
    ve_controllers ve_contr = ns->ed.ve.contr;
    ve_parameters ve_par = ns->ed.ve.par;
    eo_controllers eo_contr = ns->ed.eo.contr;
    eo_parameters eo_par = ns->ed.eo.par;
    mult_controllers mult_contr = ns->ed.mult.contr;
    mult_surf_parameters mult_spar = ns->ed.mult.spar;
    mult_parameters mult_par0 = ns->ed.mult.par0;
    mult_parameters mult_par1 = ns->ed.mult.par1;
    mult_ve_controllers mult_ve_contr = ns->ed.mult.ve.contr;
    ve_parameters mult_ve_par0 = ns->ed.mult.ve.par0;
    ve_parameters mult_ve_par1 = ns->ed.mult.ve.par1;
    eo_controllers mult_eo_contr = ns->ed.mult.eo.contr;
    eo_parameters mult_eo_par0 = ns->ed.mult.eo.par0;
    eo_parameters mult_eo_par1 = ns->ed.mult.eo.par1;
    physical_parameters *p_par = (physical_parameters *) malloc(sizeof(physical_parameters));
    initialize_physical_parameters(p_par, hig_contr, hig_par, ve_contr, ve_par, eo_contr, eo_par, mult_contr, mult_spar, mult_par0, mult_par1,
    mult_ve_contr, mult_ve_par0, mult_ve_par1, mult_eo_contr, mult_eo_par0, mult_eo_par1, myrank);
    return p_par;
}

void free_physical_parameters(physical_parameters *p_par) {
    free(p_par);
    p_par = NULL;
}


distributed_property  **create_velocity_copy(higflow_solver *ns){
    distributed_property **u_copy;
    for(int dim=0; dim<DIM; dim++){
        //psfd_compute_sfbi(ns->psfdu[dim]);
	    //psfd_synced_mapper(ns->psfdu[dim]);
        u_copy[dim] = psfd_create_property(ns->psfdu[dim]);
    }
    return u_copy;
}

void free_velocity_copy(higflow_solver *ns, distributed_property *u_copy[DIM]){
    for(int dim=0; dim<DIM; dim++){
        dp_destroy(u_copy[dim]);
    }
}

void copy_velocity_facet(higflow_solver *ns, distributed_property *u_copy[DIM], distributed_property *dpu[DIM]){
    hig_facet *f;
    int flid;
    real dp_value;
    hig_cell *cell_with_facet;
    for(int dim=0; dim<DIM; dim++){
        higfit_facetiterator *fit;
        sim_facet_domain *sdu = psfd_get_local_domain(ns->psfdu[dim]);
        mp_mapper *m = sfd_get_domain_mapper(sdu);
        for(fit = sfd_get_domain_facetiterator(sdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            f = higfit_getfacet(fit);
            flid = mp_lookup(m, hig_get_fid(f));
            dp_value = dp_get_value(dpu[dim], flid);
            dp_set_value(u_copy[dim], flid, dp_value);
        }
        higfit_destroy(fit);
    }

}

real get_fdp_value_at_point(higflow_solver *ns, distributed_property *dp, psim_facet_domain *psfd, Point p){
    hig_cell *c;
    real dp_value_local, dp_value;
    Point center;
    sim_facet_domain *sfd = psfd_get_local_domain(psfd);
    c = sfd_get_cell_with_point(sfd, p);
    if(c==NULL) dp_value_local =  -INFINITY;
    else{
        hig_get_center(c, center);
        dp_value_local = compute_facet_value_at_point(sfd, center, p, 1.0, dp, ns->stn);
    }
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    if (ntasks > 1)
        MPI_Allreduce(&dp_value_local, &dp_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    else
        dp_value = dp_value_local;
    return dp_value;
}

real max_dif_fdp(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd){
    real max_dif = 0.0;
    real dif;
    hig_facet *f;
    int flid;
    real fdp1_value;
    real fdp2_value;
    
    higfit_facetiterator *fit;
    sim_facet_domain *sdu = psfd_get_local_domain(psfd);
    mp_mapper *m = sfd_get_domain_mapper(sdu);
    for(fit = sfd_get_domain_facetiterator(sdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        f = higfit_getfacet(fit);
        flid = mp_lookup(m, hig_get_fid(f));
        fdp1_value = dp_get_value(fdp1, flid);
        fdp2_value = dp_get_value(fdp2, flid);
        dif = fabs(fdp1_value - fdp2_value);
        if(dif > max_dif){
            max_dif = dif;
        }
    }
    higfit_destroy(fit);
    return max_dif;
}

real max_dif_fdp_midline(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd, real midlinex){
    real max_dif = 0.0;
    real dif;
    hig_facet *f;
    int flid;
    real fdp1_value;
    real fdp2_value;
    hig_cell *cell_with_facet;
    Point lowpoint, highpoint;
    
    higfit_facetiterator *fit;
    sim_facet_domain *sdu = psfd_get_local_domain(psfd);
    mp_mapper *m = sfd_get_domain_mapper(sdu);
    for(fit = sfd_get_domain_facetiterator(sdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        f = higfit_getfacet(fit);
        cell_with_facet = hig_get_facet_cell(f);
        hig_get_lowpoint(cell_with_facet, lowpoint);
        hig_get_highpoint(cell_with_facet, highpoint);
        if(POS_GE(midlinex, lowpoint[0]) && POS_LE(midlinex, highpoint[0])) {
            flid = mp_lookup(m, hig_get_fid(f));
            fdp1_value = dp_get_value(fdp1, flid);
            fdp2_value = dp_get_value(fdp2, flid);
            dif = fabs(fdp1_value - fdp2_value);
            if(dif > max_dif){
                max_dif = dif;
            }
        }
    }
    higfit_destroy(fit);

    return max_dif;
}

real max_dif_fdp_midrange(higflow_solver *ns, distributed_property *fdp1, distributed_property *fdp2, psim_facet_domain *psfd, real midlinex, real Lx){
    real max_dif = 0.0;
    real dif;
    hig_facet *f;
    int flid;
    real fdp1_value;
    real fdp2_value;
    Point fcenter;
    real leftx = midlinex - 0.25*Lx;
    real rightx = midlinex + 0.25*Lx;
    
    higfit_facetiterator *fit;
    sim_facet_domain *sdu = psfd_get_local_domain(psfd);
    mp_mapper *m = sfd_get_domain_mapper(sdu);
    for(fit = sfd_get_domain_facetiterator(sdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        f = higfit_getfacet(fit);
        hig_get_facet_center(f, fcenter);
        if(POS_GE(fcenter[0], leftx) && POS_LE(fcenter[0], rightx)) {
            flid = mp_lookup(m, hig_get_fid(f));
            fdp1_value = dp_get_value(fdp1, flid);
            fdp2_value = dp_get_value(fdp2, flid);
            dif = fabs(fdp1_value - fdp2_value);
            if(dif > max_dif){
                max_dif = dif;
            }
        }
    }
    higfit_destroy(fit);
    return max_dif;
}

real max_dp(higflow_solver *ns, distributed_property *dp){
    real max = 0.0;
    real dp_value;
    hig_cell *c;
    int clid;
    higcit_celliterator *cit;
    sim_domain *sd = psd_get_local_domain(ns->psdp);
    mp_mapper *m = sd_get_domain_mapper(sd);
    for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        clid = mp_lookup(m, hig_get_cid(c));
        dp_value = dp_get_value(dp, clid);
        if(dp_value > max) max = dp_value;
    }
    higcit_destroy(cit);
    return max;
}

real max_fdp(higflow_solver *ns, distributed_property *fdp){
    real max = 0.0;
    real fdp_value;
    hig_facet *f;
    int flid;
    higfit_facetiterator *fit;
    sim_facet_domain *sfdu = psfd_get_local_domain(ns->psfdu[0]);
    mp_mapper *m = sfd_get_domain_mapper(sfdu);
    for(fit = sfd_get_domain_facetiterator(sfdu); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
        f = higfit_getfacet(fit);
        flid = mp_lookup(m, hig_get_fid(f));
        fdp_value = dp_get_value(fdp, flid);
        if(fdp_value > max) max = fdp_value;
    }
    higfit_destroy(fit);
    return max;
}

dp_copies *create_dp_copies(higflow_solver *ns){
    dp_copies* copies = (dp_copies *)malloc(sizeof(dp_copies)); 
    for(int dim=0; dim<DIM; dim++){
        copies->u[dim] = psfd_create_property(ns->psfdu[dim]);
    }
    if((ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.viscoelastic_either == true) || ns->contr.flowtype == VISCOELASTIC) {
        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                copies->Kernel[dim][dim2] = psd_create_property(ns->ed.psdED);
            }
        }
    }
    if(ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        copies->psi = psd_create_property(ns->ed.eo.psdEOpsi);
    }
    if(ns->contr.flowtype == MULTIPHASE) {
        copies->fracvol = psd_create_property(ns->ed.psdED);
    }
    return copies;
}

void free_dp_copies(dp_copies *copies, higflow_solver *ns){
    for(int dim=0; dim<DIM; dim++){
        dp_destroy(copies->u[dim]);
    }
    if((ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.viscoelastic_either == true) || ns->contr.flowtype == VISCOELASTIC) {
        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                dp_destroy(copies->Kernel[dim][dim2]);
            }
        }
    }
    if(ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        dp_destroy(copies->psi);
    }
    if(ns->contr.flowtype == MULTIPHASE) {
        dp_destroy(copies->fracvol);
    }
}

void copy_dps(higflow_solver *ns, dp_copies *copies){
    real dp_value;
    for(int dim=0; dim<DIM; dim++){
        higfit_facetiterator *fit;
        hig_facet *f;
        int flid;
        sim_facet_domain *sfd = psfd_get_local_domain(ns->psfdu[dim]);
        mp_mapper *m = sfd_get_domain_mapper(sfd);
        for(fit = sfd_get_domain_facetiterator(sfd); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            f = higfit_getfacet(fit);
            flid = mp_lookup(m, hig_get_fid(f));
            dp_value = dp_get_value(ns->dpu[dim], flid);
            dp_set_value(copies->u[dim], flid, dp_value);
        }
        higfit_destroy(fit);
    }
    
    if(ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.viscoelastic_either == true) {
        hig_cell *c;
        int clid;
        higcit_celliterator *cit;
        sim_domain *sd = psd_get_local_domain(ns->ed.psdED);
        mp_mapper *m = sd_get_domain_mapper(sd);
        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
                    c = higcit_getcell(cit);
                    clid = mp_lookup(m, hig_get_cid(c));
                    dp_value = dp_get_value(ns->ed.ve.dpKernel[dim][dim2], clid);
                    dp_set_value(copies->Kernel[dim][dim2], clid, dp_value);
                }
                higcit_destroy(cit);
            }
        }
    }

    if(ns->contr.flowtype == VISCOELASTIC) {
        hig_cell *c;
        int clid;
        higcit_celliterator *cit;
        sim_domain *sd = psd_get_local_domain(ns->ed.psdED);
        mp_mapper *m = sd_get_domain_mapper(sd);
        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
                    c = higcit_getcell(cit);
                    clid = mp_lookup(m, hig_get_cid(c));
                    dp_value = dp_get_value(ns->ed.ve.dpKernel[dim][dim2], clid);
                    dp_set_value(copies->Kernel[dim][dim2], clid, dp_value);
                }
                higcit_destroy(cit);
            }
        }
    }

    if(ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true) {
        hig_cell *c;
        int clid;
        higcit_celliterator *cit;
        sim_domain *sd = psd_get_local_domain(ns->ed.psdED);
        mp_mapper *m = sd_get_domain_mapper(sd);
        for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
            c = higcit_getcell(cit);
            clid = mp_lookup(m, hig_get_cid(c));
            dp_value = dp_get_value(ns->ed.eo.dppsi, clid);
            dp_set_value(copies->psi, clid, dp_value);
        }
        higcit_destroy(cit);
    }

    if(ns->contr.eoflow == true) {
        hig_cell *c;
        int clid;
        higcit_celliterator *cit;
        sim_domain *sd = psd_get_local_domain(ns->ed.psdED);
        mp_mapper *m = sd_get_domain_mapper(sd);
        for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
            c = higcit_getcell(cit);
            clid = mp_lookup(m, hig_get_cid(c));
            dp_value = dp_get_value(ns->ed.eo.dppsi, clid);
            dp_set_value(copies->psi, clid, dp_value);
        }
        higcit_destroy(cit);
    }

    if(ns->contr.flowtype == MULTIPHASE) {
        hig_cell *c;
        int clid;
        higcit_celliterator *cit;
        sim_domain *sd = psd_get_local_domain(ns->ed.mult.psdmult);
        mp_mapper *m = sd_get_domain_mapper(sd);
        for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
            c = higcit_getcell(cit);
            clid = mp_lookup(m, hig_get_cid(c));
            dp_value = dp_get_value(ns->ed.mult.dpfracvol, clid);
            dp_set_value(copies->fracvol, clid, dp_value);
        }
        higcit_destroy(cit);
    }

}

// not parallel yet
real compute_mid_end_err_x(higflow_solver *ns, real midlinex, real Lx){
    real near_endx = midlinex + 0.75*Lx/2.0;
    real err, max_err = 0.0;
    
    Point lp, hp, ccenter, cdelta, fcenter;
    real lx, hx, dl, dh;
    Point midline_pt, near_end_pt;
    midline_pt[0] = midlinex; near_end_pt[0] = near_endx;
    real midline_val, near_end_val;
    real succ;

    hig_cell *c;
    //hig_facet *f;
    higcit_celliterator *cit;
    sim_domain *sd = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfd = psfd_get_local_domain(ns->psfdu[0]);
    distributed_property *dpu = ns->dpu[0];
    mp_mapper *m = sd_get_domain_mapper(sd);

    for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        hig_get_highpoint(c,hp); hig_get_lowpoint(c,lp);
        hx = hp[0]; lx = lp[0];
        if(POS_LE(midlinex, hx) && POS_GE(midlinex, lx)){ // centerline of channel in cell
            ////////// get midline value
            hig_get_center(c,ccenter);
            fcenter[1] = ccenter[1];
            midline_pt[1] = ccenter[1]; 
            near_end_pt[1] = ccenter[1];
            hig_get_delta(c, cdelta);
            dl = midlinex - lx; dh = hx - midlinex;
            if(dl < dh) fcenter[0] = ccenter[0] - 0.5*cdelta[0];
            else fcenter[0] = ccenter[0] + 0.5*cdelta[0];
            // succ = sfd_get_facet_with_point(sfd, fcenter, f); // get facet closest to point
            midline_val = compute_facet_value_at_point(sfd, fcenter, midline_pt, 1.0, dpu, ns->stn);
            /////////// get near-end value
            c = sd_get_cell_with_point(sd, near_end_pt);
            hig_get_highpoint(c,hp); hig_get_lowpoint(c,lp);
            hx = hp[0]; lx = lp[0];
            hig_get_center(c,ccenter);
            fcenter[1] = ccenter[1]; 
            near_end_pt[1] = ccenter[1];
            hig_get_delta(c, cdelta);
            dl = midlinex - lx; dh = hx - midlinex;
            if(dl < dh) fcenter[0] = ccenter[0] - 0.5*cdelta[0];
            else fcenter[0] = ccenter[0] + 0.5*cdelta[0];
            // succ = sfd_get_facet_with_point(sfd, fcenter, f); // get facet closest to point
            near_end_val = compute_facet_value_at_point(sfd, fcenter, near_end_pt, 1.0, dpu, ns->stn);
            /////////// compute error
            err = fabs(midline_val - near_end_val);
            if(err > max_err) max_err = err;
        }
    }
    higcit_destroy(cit);

    return max_err;
}

// not parallel yet
real compute_inlet_mid_err_x(higflow_solver *ns, real midlinex, real Lx){
    real inletx = midlinex - Lx/2.0;
    real err, max_err = 0.0;
    
    Point lp, hp, ccenter, cdelta, fcenter;
    real lx, hx, dl, dh;
    Point midline_pt, inlet_pt;
    midline_pt[0] = midlinex; inlet_pt[0] = inletx;
    real midline_val, inlet_val;
    real succ;

    hig_cell *c;
    //hig_facet *f;
    higcit_celliterator *cit;
    sim_domain *sd = psd_get_local_domain(ns->psdp);
    sim_facet_domain *sfd = psfd_get_local_domain(ns->psfdu[0]);
    distributed_property *dpu = ns->dpu[0];
    mp_mapper *m = sd_get_domain_mapper(sd);

    for(cit = sd_get_domain_celliterator(sd); !higcit_isfinished(cit); higcit_nextcell(cit)) {
        c = higcit_getcell(cit);
        hig_get_highpoint(c,hp); hig_get_lowpoint(c,lp);
        hx = hp[0]; lx = lp[0];
        if(POS_LE(midlinex, hx) && POS_GE(midlinex, lx)){ // centerline of channel in cell
            ////////// get midline value
            hig_get_center(c,ccenter);
            fcenter[1] = ccenter[1];
            midline_pt[1] = ccenter[1]; 
            inlet_pt[1] = ccenter[1];
            hig_get_delta(c, cdelta);
            dl = midlinex - lx; dh = hx - midlinex;
            if(dl < dh) fcenter[0] = ccenter[0] - 0.5*cdelta[0];
            else fcenter[0] = ccenter[0] + 0.5*cdelta[0];
            // succ = sfd_get_facet_with_point(sfd, fcenter, f); // get facet closest to point
            midline_val = compute_facet_value_at_point(sfd, fcenter, midline_pt, 1.0, dpu, ns->stn);
            /////////// get inlet value
            inlet_val = ns->func.get_boundary_velocity(0, inlet_pt, 0, ns->par.t);
            /////////// compute error
            err = fabs(midline_val - inlet_val);
            if(err > max_err) max_err = err;
        }
    }
    higcit_destroy(cit);

    return max_err;
}

long int get_current_mem_usage() {
    struct rusage ru_mem;
    getrusage(RUSAGE_SELF, &ru_mem);
    long int local_usage = ru_mem.ru_maxrss;
    long int total_usage;
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    if (ntasks > 1)
        MPI_Allreduce(&local_usage, &total_usage, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    else
        total_usage = local_usage;
    return total_usage;
}

void fill_column_string_util(const char *literal, char *str, size_t size) {
    // Copy the literal to the beginning of the string
    size_t literalLength = strlen(literal);
    size_t copyLength = (literalLength < size) ? literalLength : size;
    
    strncpy(str, literal, copyLength);

    // Fill the remaining space with spaces
    size_t remainingSpace = size - copyLength;
    if (remainingSpace > 0) {
        memset(str + copyLength, ' ', remainingSpace);
    }
    memset(str + size, '\0', 1);
}

void get_namemem(char **source, char *destination) {
    // Find the second to last '/' character
    char *last_slash = NULL;
    char *current = *source;
    
    while (*current != '\0') {
        if (*current == '/') last_slash = current;
        current++;
    }

    // Copy the portion of the string up to the second to last '/' character
    if (last_slash != NULL) {
        strncpy(destination, *source, last_slash - *source);
        destination[last_slash - *source] = '\0';
    } else {
        strcpy(destination, *source); // No '/' found, copy the entire string
    }

    // Append "/res/" to the destination string
    strcat(destination, "/mem.txt");
}

void write_mem_usage(higflow_solver *ns, char *brief_desc) {
    
    long int usage = get_current_mem_usage();
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank == 0)  {
        int step = ns->par.step;
        FILE *f;
        char filename[1024];
        strcpy(filename, ns->par.namesave);
        get_namemem(&ns->par.namesave, filename);

        int step_len = 9;
        int desc_len = 20;
        int num_len = 10;
        char step_str[step_len + 1];
        char desc_str[desc_len + 1];
        char num_str[num_len + 1];
        
        if(step == 0){
            f = fopen(filename, "w");

            fill_column_string_util("step", step_str, step_len);
            fprintf(f, "%s", step_str);
            fill_column_string_util("desc", desc_str, desc_len);
            fprintf(f, "%s", desc_str);
            fill_column_string_util("mem", num_str, num_len);
            fprintf(f, "%s", num_str);
            fprintf(f, "%s", "\n");

            fclose(f);
        }

        f = fopen(filename, "a");

        sprintf(step_str, "%d", step);
        fill_column_string_util(step_str, step_str, step_len);
        fprintf(f, "%s", step_str);
        fill_column_string_util(brief_desc, desc_str, desc_len);
        fprintf(f, "%s", desc_str);
        sprintf(num_str, "%ld", usage);
        fill_column_string_util(num_str, num_str, num_len);
        fprintf(f, "%s", num_str);
        fprintf(f, "%s", "\n");
        
        fclose(f);
    }

}

real solve_un_u(visc_model_type visc_model, physical_parameters *p_par) {
    real k = 1.5;
    real De = p_par->De.val;
    real beta = p_par->beta.val;
    real adim_group;
    if(visc_model == LPTT){ //ptt
        real eps = p_par->epsilon_ptt.val;
        adim_group = sqrt(eps)*De;
    } else if(visc_model == FENE_P) { //FENE-P
        real L2 = p_par->L2_fene.val;
        real a = L2/(L2-3);
        adim_group = sqrt(1.0/L2)*(De/a);
    } else if(visc_model == OLDROYD_B) { // Oldroyd-B
        return beta;
    } else{
        print0f("will not calculate un/U - no suitable viscoelastic model\n");
        return -1.0;
    }

    real un_u;
    real tol = 1.0e-9;
    real ag2 = adim_group*adim_group;

    if(beta > 1.0 - tol) return 1.0; // newtonian case

    if(beta < tol) { // purely polymeric fluid
        real b = 24.0/5.0*k*k*ag2;
        real d = sqrt(27.0*b + 4) + sqrt(27.0*b);
        un_u = pow(432.0,1.0/6.0)*(pow(d,2.0/3.0)-pow(2.0,2.0/3.0))
                                      /(6.0*sqrt(b)*cbrt(d));
        return un_u;
    }

    real x0 = 0.5;
    real delta = 1.0;
    int maxiter = 20;

    real A, A3, C, d_C, C2, C3, S, d_S, Fp, d_Fp, Fm, d_Fm, Gp, d_Gp, Gm, d_Gm,
         Hp, d_Hp, Hm, d_Hm, Ic, d_Ic, Ip, d_Ip, Im, d_Im, Ix, d_Ix;
    
    real x = x0;
    for (int i = 0; i < maxiter; i++) { // Newton iterations
        A = 1.0/(6.0*ag2*beta); A3 = A*A*A;
        C = -3.0*k*A*x; C2 = C*C; C3 = C2*C; d_C = -3.0*k*A;
        S = sqrt(A3 + C2); d_S = C*d_C/S;
        Fp = cbrt(C + S); d_Fp = (d_C + d_S)/(3.0*Fp*Fp);
        Fm = cbrt(C - S); d_Fm = (d_C - d_S)/(3.0*Fm*Fm);
        Gp = 3.0*C + S; d_Gp = 3.0*d_C + d_S;
        Gm = 3.0*C - S; d_Gm = 3.0*d_C - d_S;
        Hp = 8.0*A3 + C*(-19.0*C + 9.0*S); d_Hp = -38.0*C*d_C + 9.0*(d_C*S + C*d_S);
        Hm = 8.0*A3 + C*(-19.0*C - 9.0*S); d_Hm = -38.0*C*d_C - 9.0*(d_C*S + C*d_S);
        Ic = Fp*Gm + Fm*Gp; d_Ic = d_Fp*Gm + Fp*d_Gm + d_Fm*Gp + Fm*d_Gp;
        Ip = Fp*Hp - 8.0*pow(A, 3.5); d_Ip = d_Fp*Hp + Fp*d_Hp;
        Im = Fm*Hm + 8.0*pow(A, 3.5); d_Im = d_Fm*Hm + Fm*d_Hm;
        Ix = x/beta + 3.0/8.0*(1.0-beta)/beta * (Ic/C + 3.0/35.0*(Ip + Im)/C2) - 1.0;
        d_Ix = 1.0/beta + 3.0/8.0*(1.0-beta)/beta * 
               (-Ic*d_C/C2 + d_Ic/C - 6.0/35.0*d_C*(Ip+Im)/C3 + 3.0/35.0*(d_Ip+d_Im)/C2);

        delta = Ix/d_Ix;
        x = x - delta;
        if(fabs(delta) < tol) break;
    }

    un_u = x;
    return un_u;
}

real calc_u_ptt_fene(visc_model_type visc_model, physical_parameters *p_par, real y){
    real un_u = p_par->Un_U.val;
    
    real k = 1.5;
    real De = p_par->De.val;
    real beta = p_par->beta.val;
    real adim_group;
    if(visc_model == LPTT){ //ptt
        real eps = p_par->epsilon_ptt.val;
        adim_group = sqrt(eps)*De;
    } else if(visc_model == FENE_P) { //FENE-P
        real L2 = p_par->L2_fene.val;
        real a = L2/(L2-3);
        adim_group = sqrt(1.0/L2)*(De/a);
    } else if(visc_model == OLDROYD_B) { // Oldroyd-B
        return 1.5 * (1.0 - y*y);
    } else{
        print0f("will not calculate u/U - no suitable viscoelastic model\n");
        return -1.0;
    }

    real uy_u;
    real y2 = y*y;
    real tol = 1.0e-7;
    real ag2 = adim_group*adim_group;

    if(beta > 1.0 - tol) { // newtonian case
        real uy_u = k*un_u*(1.0 - y2);
        return uy_u;
    }
    if(beta < tol) {
        real uy_u = k*un_u*(1.0 - y2)*(1.0 + 4.0*ag2*k*k*un_u*un_u*(1.0 + y2));
        return uy_u;
    }

    real A, A3, C, C2, Cy, Cy2, S, Sy, Fp, Fm, Gp, Gm, Fyp, Fym, Gyp, Gym;

    A = 1.0/(6.0*ag2*beta); A3 = A*A*A;
    C = -3.0*k*A*un_u; C2 = C*C; Cy = C*y; Cy2 = Cy*Cy;
    S = sqrt(A3 + C2); Sy = sqrt(A3 + Cy2);
    Fp = cbrt(C + S); Fm = cbrt(C - S);
    Gp = 3.0*C + S; Gm = 3.0*C - S;
    Fyp = cbrt(Cy + Sy); Fym = cbrt(Cy - Sy);
    Gyp = 3.0*Cy + Sy; Gym = 3.0*Cy - Sy;

    uy_u = k*un_u/beta * (1.0 - y2) + 3.0/8.0 * (1.0 - beta)/beta * 1.0/C * 
           (Fp*Gm - Fyp*Gym + Fm*Gp - Fym*Gyp);
    
    return uy_u;
}

typedef struct nonlinear_pb_info {
    real alphaeo;
    real deltaeo;
    real dy;
} nonlinear_pb_info;

PetscErrorCode nonlinear_pb(SNES snes, Vec psi, Vec F, void *ctx) {

    PetscErrorCode ierr;
    const PetscScalar  *lpsi;
    PetscScalar   *f;
    PetscInt      n, ng, start, end;
    PetscScalar lpsi_minusone, lpsi_plusone;

    nonlinear_pb_info *info = (nonlinear_pb_info *) ctx;
    PetscScalar alphaeo = info->alphaeo;
    PetscScalar deltaeo = info->deltaeo;
    PetscScalar dy = info->dy;
    PetscScalar idy2 = 1.0/(dy*dy);

    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
    
    VecGetArrayRead(psi, &lpsi);
    VecGetArray(F, &f);
    VecGetLocalSize(psi, &n);
    VecGetSize(psi, &ng);

    VecGetOwnershipRange(psi, &start, &end);

    if(start == 0) f[0] = lpsi[0] + 1.0;
    else{
        ierr = MPI_Send(&lpsi[0], 1, MPI_HIGREAL, myrank-1, 0, MPI_COMM_WORLD);
        ierr = MPI_Recv(&lpsi_minusone, 1, MPI_HIGREAL, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        f[0] = (lpsi[1] - 2.0*lpsi[0] + lpsi_minusone)*idy2 - 2.0*deltaeo * sinh(alphaeo*lpsi[0]);
    }

    for (PetscInt i = 1; i < n-1; i++) {
        f[i] = (lpsi[i+1] - 2.0*lpsi[i] + lpsi[i-1])*idy2 - 2.0*deltaeo * sinh(alphaeo*lpsi[i]);
    }

    if(end==ng) f[n-1] = lpsi[n-1] + 1.0;
    else{
        ierr = MPI_Send(&lpsi[n-1], 1, MPI_HIGREAL, myrank+1, 0, MPI_COMM_WORLD);
        ierr = MPI_Recv(&lpsi_plusone, 1, MPI_HIGREAL, myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        f[n-1] = (lpsi_plusone - 2.0*lpsi[n-1] + lpsi[n-2])*idy2 - 2.0*deltaeo * sinh(alphaeo*lpsi[n-1]);
    }

    VecRestoreArrayRead(psi, &lpsi);
    ierr = VecRestoreArray(F, &f);

    return ierr;
}

PetscErrorCode jac_nonlinear_pb(SNES snes,Vec psi, Mat jac,Mat B,void *ctx)
{   
    PetscErrorCode    ierr;
    const PetscScalar  *lpsi;
    PetscInt      n, ng, start, end, j[3], i;
    PetscScalar       A[3];

    nonlinear_pb_info *info = (nonlinear_pb_info *) ctx;
    PetscScalar alphaeo = info->alphaeo;
    PetscScalar deltaeo = info->deltaeo;
    PetscScalar dy = info->dy;
    PetscScalar idy2 = 1.0/(dy*dy);

    VecGetArrayRead(psi, &lpsi);
    VecGetLocalSize(psi, &n);
    VecGetSize(psi, &ng);

    VecGetOwnershipRange(psi, &start, &end);

    i = start;
    if(start == 0){ 
        A[0] = 1.0; //be careful
        MatSetValues(jac,1,&i,1,&i,&(A[0]),INSERT_VALUES);
    }
    else{
        A[0] = idy2; A[2] = idy2;
        A[1] = -2.0*idy2 - 2.0*deltaeo*alphaeo*cosh(alphaeo*lpsi[0]);
        j[0] = start - 1; j[1] = start; j[2] = start + 1;
        MatSetValues(jac,1,&i,3,j,A,INSERT_VALUES);
    }

    for (i = start+1; i < end-1; i++) {
        A[0] = idy2; A[2] = idy2;
        A[1] = -2.0*idy2 - 2.0*deltaeo*alphaeo*cosh(alphaeo*lpsi[i-start]);
        j[0] = i - 1; j[1] = i; j[2] = i + 1;
        MatSetValues(jac,1,&i,3,j,A,INSERT_VALUES);
    }

    i = end - 1;
    if(end==ng) { 
        A[0] = 1.0; //be careful
        MatSetValues(jac,1,&i,1,&i,&(A[0]),INSERT_VALUES);
    }
    else{
        A[0] = idy2; A[2] = idy2;
        A[1] = -2.0*idy2 - 2.0*deltaeo*alphaeo*cosh(alphaeo*lpsi[n-1]);
        j[0] = end - 2; j[1] = end - 1; j[2] = end;
        MatSetValues(jac,1,&i,3,j,A,INSERT_VALUES);
    }

    VecRestoreArrayRead(psi, &lpsi);

    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
    return ierr;
}

void solve_psi_in(real *psi_in, physical_parameters *p_par, int npoints) {

    int myrank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    void _try_initialize_petsc(void);
    _try_initialize_petsc();

    if(myrank == 0) printf("=+=+ Solving psi_in =+=+=+=+=+=+=+=+=+=+=+=+=+=+=\n");
    real alphaeo = p_par->alpha_eo.val;
    real deltaeo = p_par->delta_eo.val;
    real kappaeo = p_par->kappa_eo.val;
    real dy = 2.0/(npoints-1);
    real y;

    Vec            psi_sol;
    Mat            J;             // Jacobian matrix
    SNES           snes;          // Nonlinear solver context
    PetscInt      start, end;

    VecCreate(PETSC_COMM_WORLD, &psi_sol);
    VecSetSizes(psi_sol, PETSC_DECIDE, npoints);
    VecSetFromOptions(psi_sol);

    PetscScalar   *lpsi;
    VecGetArray(psi_sol, &lpsi);
    VecGetOwnershipRange(psi_sol, &start, &end);

    for(int i=start; i<end; i++) {
        y = -1.0 + i*dy;
        lpsi[i-start] = -cosh(kappaeo*y)/cosh(kappaeo); //initial guess
    }

    SNESCreate(PETSC_COMM_WORLD, &snes);
    nonlinear_pb_info info={alphaeo, deltaeo, dy};
    SNESSetFunction(snes, NULL, &nonlinear_pb, &info);

    SNESSetType(snes, SNESNEWTONLS);
    SNESSetFromOptions(snes);

    MatCreate(PETSC_COMM_WORLD,&J);
    MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,npoints,npoints);
    MatSetFromOptions(J);
    MatSeqAIJSetPreallocation(J,3,NULL);
    MatMPIAIJSetPreallocation(J,3,NULL,3,NULL);
    SNESSetJacobian(snes,J,J,jac_nonlinear_pb,&info);

    SNESSetTolerances(snes, 1.0e-10, 1.0e-14, 1.0e-11, PETSC_DEFAULT, PETSC_DEFAULT);
    if(myrank == 0){
        PetscInt       maxit,maxf;
        PetscReal      abstol,rtol,stol;
        SNESGetTolerances(snes,&abstol,&rtol,&stol,&maxit,&maxf);
        PetscPrintf(PETSC_COMM_WORLD,"atol=%g, rtol=%g, stol=%g, maxit=%D, maxf=%D\n",(double)abstol,(double)rtol,(double)stol,maxit,maxf);
    }
    SNESSolve(snes, NULL, psi_sol);

    int vals_per_proc[nprocs];
    int n = end-start;
    MPI_Allgather(&n, 1, MPI_INT, vals_per_proc, 1, MPI_INT, MPI_COMM_WORLD);
    int displs[nprocs];
    displs[0] = 0;
    for(int i=1; i<nprocs; i++) displs[i] = displs[i-1] + vals_per_proc[i-1];

    MPI_Allgatherv(lpsi, end-start, MPI_HIGREAL, psi_in, vals_per_proc, displs, MPI_HIGREAL, MPI_COMM_WORLD);

    //prints psi to file
    if(myrank == 0) {

        PetscInt       its;
        SNESGetIterationNumber(snes,&its);
        PetscPrintf(PETSC_COMM_WORLD,"Number of SNES iterations = %D\n",its);

        SNESConvergedReason reason;
        SNESGetConvergedReason(snes,&reason);
        PetscPrintf(PETSC_COMM_WORLD,"Converged reason = %D\n",reason);

        /*FILE *fpsi = fopen("psi.dat", "w");
        for(int i=0; i<npoints; i++) {
            y = -1.0 + i*dy;
            fprintf(fpsi, "%20.16lf %20.16lf\n", y, psi_in[i]);
        }
        fclose(fpsi);
        FILE *fpsi_err = fopen("psi_lin_dif.dat", "w");
        for(int i=0; i<npoints; i++) {
            y = -1.0 + i*dy;
            fprintf(fpsi_err, "%20.16lf %20.16lf\n", y, psi_in[i] - (-cosh(kappaeo*y)/cosh(kappaeo)) );
        }
        fclose(fpsi_err);
        FILE *fpsi_nl = fopen("visualization/1dpsi_nonlinear_solution", "r");
        FILE *fpsi_dif = fopen("psi_nl_dif.dat", "w");
        for(int i=0; i<npoints; i++) {
            real pos,nl;
            y = -1.0 + i*dy;
            int err = fscanf(fpsi_nl, "%lf %lf\n", &pos, &nl);
            fprintf(fpsi_dif, "%20.16lf %20.16lf\n", y, psi_in[i]-nl);
        }
        fclose(fpsi_nl);
        fclose(fpsi_dif);*/
    }

    VecDestroy(&psi_sol);
    MatDestroy(&J);
    SNESDestroy(&snes);
}

real get_psi_val_from_array(real *psi_in, real y, int npoints){
    real dy = 2.0/(npoints-1);
    int loc = (int) ((y+1.0)/dy);
    return psi_in[loc];
}

real get_psi_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y){
    real psi_pb, ek1my, ek1py, tanhma4;

    ek1my = exp(kappa*(1.0-y));
    ek1py = exp(kappa*(1.0+y));
    tanhma4 = tanh(-alpha/4.0);
    real infsol_m1 = 2.0/alpha * log((ek1my+tanhma4)/(ek1my-tanhma4));
    real infsol_p1 = 2.0/alpha * log((ek1py+tanhma4)/(ek1py-tanhma4));
    psi_pb = (-psi_up) * infsol_m1 + (-psi_down) * infsol_p1;
    
    return psi_pb;
}

real get_nplus_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y){
    real ek1my, ek1py, tanhma4;

    ek1my = exp(kappa*(1.0-y));
    ek1py = exp(kappa*(1.0+y));
    tanhma4 = tanh(-alpha/4.0);
    real infsol_m1 = (ek1my-tanhma4)/(ek1my+tanhma4);
    real infsol_p1 = (ek1py-tanhma4)/(ek1py+tanhma4);
    real sqrt_np_pb = pow(infsol_m1, -psi_up) * pow(infsol_p1, -psi_down);
    real np_pb = sqrt_np_pb*sqrt_np_pb;
    
    return np_pb;
}

real get_nminus_channel_sol(real alpha, real kappa, real psi_up, real psi_down, real y){
    real ek1my, ek1py, tanhma4;

    ek1my = exp(kappa*(1.0-y));
    ek1py = exp(kappa*(1.0+y));
    tanhma4 = tanh(-alpha/4.0);
    real infsol_m1 = (ek1my+tanhma4)/(ek1my-tanhma4);
    real infsol_p1 = (ek1py+tanhma4)/(ek1py-tanhma4);
    real sqrt_nm_pb = pow(infsol_m1, -psi_up) * pow(infsol_p1, -psi_down);
    real nm_pb = sqrt_nm_pb*sqrt_nm_pb;
    
    return nm_pb;
}





xdmf_output *output;
xdmf_output *output_mult;
bool reuse_grid;
bool reuse_grid_mult;


void write_init(higflow_solver *ns) {
    higio_set_chunk_size(1024);

    output = xdmf_init(ns->par.nameprint, ns->sdp);
    xdmf_register_cell_property(output, "p", 1, ns->dpp);
    xdmf_register_facet_property_array(output, "u", DIM, ns->sfdu, ns->dpu);

    if((ns->contr.flowtype != MULTIPHASE && ns->contr.flowtype == VISCOELASTIC) ||
        (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.viscoelastic_either == true)) {
        xdmf_register_cell_property(output, "Tp", 3, 
        ns->ed.ve.dpTaup[0][0], ns->ed.ve.dpTaup[0][1], ns->ed.ve.dpTaup[1][1]);
    }

    if((ns->contr.flowtype != MULTIPHASE && ns->contr.eoflow == true) ||
       (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        xdmf_register_cell_property(output, "phi", 1, ns->ed.eo.dpphi);
        xdmf_register_cell_property(output, "psi", 1, ns->ed.eo.dppsi);
        xdmf_register_facet_property_array(output, "Feo", DIM, ns->ed.eo.sfdEOFeo, ns->ed.eo.dpFeo);
        if((ns->contr.flowtype != MULTIPHASE && ns->ed.eo.contr.eo_model == PNP) ||
           (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.eo.contr.eo_model == PNP)) {
            xdmf_register_cell_property(output, "n+", 1, ns->ed.eo.dpnplus);
            xdmf_register_cell_property(output, "n-", 1, ns->ed.eo.dpnminus);
        }
    }

    if(ns->contr.flowtype == MULTIPHASE) {
        char mul_nameprint[1024];
        strcpy(mul_nameprint, ns->par.nameprint);
        strcat(mul_nameprint, "_mult");

        output_mult = xdmf_init(mul_nameprint, ns->ed.mult.sdmult);
        xdmf_register_cell_property(output_mult, "fvol", 1, ns->ed.mult.dpfracvol);
        xdmf_register_cell_property(output_mult, "kappa", 1, ns->ed.mult.dpcurvature);
        reuse_grid_mult = false;
    }

    reuse_grid = false;

}

void write_xdmf(higflow_solver *ns) {

    xdmf_write_timestep(output, ns->par.step, reuse_grid);
    if(reuse_grid == false) reuse_grid = true;

    if(ns->contr.flowtype == MULTIPHASE) {
        xdmf_write_timestep(output_mult, ns->par.step, reuse_grid_mult);
        if(reuse_grid_mult == false) reuse_grid_mult = true;
    }
    
}
