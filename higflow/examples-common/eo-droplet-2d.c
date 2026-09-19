// ---------------------------------------------------------------------------
// Infraestrutura comum ao example2d_DynamicMeshAdapt e ao
// example2d_ElectroOsmotic.
//
// Os dois sao exemplos distintos e os `main` deles continuam separados -- e' o
// `main` que carrega a identidade de cada demonstracao.  O que esta' aqui e' o
// maquinario que os dois tinham identico byte a byte e que nao diz nada sobre
// qual problema cada um resolve:
//
//   get_fracvol                    fracao volumetrica por subamostragem 16x16,
//                                  delegando o caso analitico a get_fracvolN
//   compute_deformation_parameter  diagnostico da gota a partir das retas PLIC
//   save_deformation_parameter     gravacao desse diagnostico
//   init_global_var                cacheia ponteiros de `ns` nos globais
//   get_inlet_types                le os tipos de entrada do .bc.yaml
//   solver_step                    despacho do passo pelo flowtype
//   errors                         norma do residuo e criterio de parada
//
// O que NAO veio para ca', mesmo estando duplicado: get_kernel, get_viscosity,
// calculate_m_user e companhia.  Sao curtas, mas e' nelas que cada exemplo
// declara o modelo que usa; compartilha-las esconderia essa escolha.
//
// Este arquivo e' INCLUIDO, nao compilado a parte: get_fracvol chama func() e
// get_fracvolN, que cada exemplo define por conta propria.
// ---------------------------------------------------------------------------

real get_fracvol(Point center, Point delta, real t) {
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

real compute_deformation_parameter(higflow_solver *ns) {
    int num_plic = ns->ed.mult.num_plic_lines;
    int num_plic_global;
    Point center, center_global, midseg;
    
    center[0] = 0.0; center[1] = 0.0;
    for(int i=0; i<num_plic; i++) {
        midseg[0] = 0.5*(ns->ed.mult.plic_lines[i][0][0] + ns->ed.mult.plic_lines[i][1][0]);
        midseg[1] = 0.5*(ns->ed.mult.plic_lines[i][0][1] + ns->ed.mult.plic_lines[i][1][1]);
        center[0] += midseg[0];
        center[1] += midseg[1];
    }
    MPI_Allreduce(center, center_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&num_plic, &num_plic_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    center_global[0] /= num_plic_global;
    center_global[1] /= num_plic_global;
    if(psi_down == 0.0 && ns->ed.mult.contr.eoflow_either == true) {
        center_global[1] = 0.0;
    }

    real Rmin=INFINITY, Rmax=0.0, dist2;
    real Rmin_global, Rmax_global;
    for(int i=0; i<num_plic; i++) {
        midseg[0] = 0.5*(ns->ed.mult.plic_lines[i][0][0] + ns->ed.mult.plic_lines[i][1][0]);
        midseg[1] = 0.5*(ns->ed.mult.plic_lines[i][0][1] + ns->ed.mult.plic_lines[i][1][1]);
        dist2 = (midseg[0] - center_global[0])*(midseg[0] - center_global[0]) 
              + (midseg[1] - center_global[1])*(midseg[1] - center_global[1]);
        if(dist2 < Rmin) Rmin = dist2;
        if(dist2 > Rmax) Rmax = dist2;
    }
    Rmin = sqrt(Rmin); Rmax = sqrt(Rmax);
    MPI_Allreduce(&Rmin, &Rmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&Rmax, &Rmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return (Rmax_global-Rmin_global)/(Rmax_global+Rmin_global);
}

void save_deformation_parameter(higflow_solver *ns, int myrank) {
    real defpar = compute_deformation_parameter(ns);
    if(myrank==0) {
        char filename[1024];
        sprintf(filename, "%s.deformation_parameter.txt", ns->par.namesave);
        FILE *f;
        if(ns->par.step==0) f = fopen(filename, "w");
        else f = fopen(filename, "a");
        fprintf(f, "%d %.10lf %.10lf\n", ns->par.step, ns->par.t, defpar);
        fclose(f);
    }
}

void init_global_var(higflow_solver* ns) {
    flowtype = ns->contr.flowtype;
    sdp_ptr = &ns->sdp;
    stn_ptr = &ns->stn;
    dpp_ptr = &ns->dpp;
    sfdv_ptr = &(ns->sfdu[1]);
    dpvstar_ptr = &ns->dpustar[1];
    if(flowtype == VISCOELASTIC) visc_model = ns->ed.ve.contr.model;
    if(flowtype == MULTIPHASE) {
        dpfracvol_ptr = &ns->ed.mult.dpfracvol;
        sdmult_ptr = &ns->ed.mult.sdmult;
        stnmult_ptr = &ns->ed.mult.stn;
        flowtype0 = ns->ed.mult.contr.flowtype0;
        flowtype1 = ns->ed.mult.contr.flowtype1;
        viscoelastic_either = ns->ed.mult.contr.viscoelastic_either;
        eoflow0 = ns->ed.mult.contr.eoflow0;
        eoflow1 = ns->ed.mult.contr.eoflow1;
        eoflow_either = ns->ed.mult.contr.eoflow_either;
    }
    eoflow = ns->contr.eoflow;
    if(eoflow == true) {
        sfdFeoy_ptr = &ns->ed.eo.sfdEOFeo[1];
        dpFeoy_ptr = &ns->ed.eo.dpFeo[1];
        stnFeoy_ptr = &ns->ed.eo.stnpsi;
    }
}

void get_inlet_types(higflow_solver* ns) {
    char namefile[1024];
    sprintf(namefile,"%s.bc.yaml",ns->par.nameload);
    
    FILE *fbc = fopen(namefile, "r");
    struct fy_document *fyd = NULL;
    fyd = fy_document_build_from_file(NULL, namefile);
     
    if (fyd == NULL) {
        printf("=+=+=+= Error loading file %s =+=+=+=\n",namefile);
        exit(1);
    }

    char aux[1024];
    int ifd = fy_document_scanf(fyd,"/bc/bc0/velocity_0/type %s",aux);
    if (strcmp(aux,"dirichlet") == 0) u_inlet = DIRICHLET;
    else if (strcmp(aux,"neumann") == 0) u_inlet = NEUMANN;
    else {
        printf("=+=+=+= Error loading boundary condition type for the inlet velocity\n");
        exit(1);
    }

    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        ifd = fy_document_scanf(fyd,"/bc_electroosmotic/bc0/psi/type %s",aux);
        if (strcmp(aux,"dirichlet") == 0) psi_inlet = DIRICHLET;
        else if (strcmp(aux,"neumann") == 0) psi_inlet = NEUMANN;
        else {
            printf("=+=+=+= Error loading boundary condition type for the inlet velocity\n");
            exit(1);
        }
    }
    
    fy_document_destroy(fyd);
    fclose(fbc);
}

void solver_step(higflow_solver* ns) {
    if (ns->contr.eoflow == true || (ns->contr.flowtype == MULTIPHASE && ns->ed.mult.contr.eoflow_either == true)) {
        switch (ns->contr.flowtype) {
            case NEWTONIAN:
                higflow_solver_step_electroosmotic(ns);
                break;
            case MULTIPHASE:
                if(ns->ed.mult.contr.viscoelastic_either == true)
                    higflow_solver_step_multiphase_electroosmotic_viscoelastic(ns);
                else
                    higflow_solver_step_multiphase_electroosmotic(ns);
                break;
            case VISCOELASTIC:
                higflow_solver_step_electroosmotic_viscoelastic(ns);
                break;
        }
    }
    else {
        switch (ns->contr.flowtype) {
            case NEWTONIAN:
                higflow_solver_step(ns);
                break;
            case MULTIPHASE:
                if(ns->ed.mult.contr.viscoelastic_either == true) 
                    higflow_solver_step_multiphase_viscoelastic(ns);
                else
                    higflow_solver_step_multiphase(ns);
                ns->par.stepaux=ns->par.stepaux+1;
                break;
            case GENERALIZED_NEWTONIAN:
                higflow_solver_step_gen_newt(ns);
                break;
            case VISCOELASTIC:
                higflow_solver_step_viscoelastic(ns);
                break;
            case VISCOELASTIC_INTEGRAL:
                higflow_solver_step_viscoelastic_integral(ns);
                break;
        }
    }
}

int errors(higflow_solver* ns, sim_residuals* sim_res, int myrank) {
    int errcode = 0;

    if (sim_res != NULL) {
        real dudt_norm = 0.0;
        write_residuals(sim_res, ns);
        if (myrank == 0) dudt_norm = sim_res->u[0]->midrange->res_max->avg[0] / ns->par.dt;

        MPI_Bcast(&dudt_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        print0f("|    dudt_norm = %15.10lf", dudt_norm);
        if(flowtype != MULTIPHASE) {
            real tol = 1.0e-5;
            if (ns->contr.eoflow == true) tol = 5.0e-5;
            if (dudt_norm < max(1.0e-10 / ns->par.dt, tol)) {
                print0f("\nsteady state reached\n");
                errcode = 1;
            }
            if (dudt_norm > 1.0e8) {
                print0f("\nsimulation 'diverged'\n");
                errcode = -1;
            }
        }
    }

    return errcode;
}
