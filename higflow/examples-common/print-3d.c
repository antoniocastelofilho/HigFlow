// ---------------------------------------------------------------------------
// Impressao de perfis 3-D em cortes e em pontos.
//
// As quatro funcoes eram identicas byte a byte no example3d_complex e no
// example3d_lid_driven.  Nao ha' nada de cada exemplo nelas: sao percursos de
// dominio da biblioteca gravando .dat.  O que cada exemplo escolhe -- qual
// corte, qual componente, com que frequencia -- e' argumento na chamada, e
// continua no main de cada um.
//
// Incluido, nao compilado a parte: assim nenhum dos dois Makefiles precisa de
// regra nova, e o objeto continua sendo um so' por exemplo.
// ---------------------------------------------------------------------------

void print_tensor(higflow_solver *ns, int myrank, int i, int j, int dimprint1, real pprint1, int dimprint2, real pprint2) {
    char filename[1024];
    snprintf(filename,sizeof filename,"tensor-%d-%d-%d.dat",i,j,myrank);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpDu[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpDu[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpTaup[i][j], ns->ed.stn);
            real T = S + 2.0*(1-beta)*D/Re; 
            //Print polymeric stress data file
            if ((fabs(ccenter[dimprint1] - pprint1) < 0.5*cdelta[dimprint1]) &&
                (fabs(ccenter[dimprint2] - pprint2) < 0.5*cdelta[dimprint2])) { 
                fprintf(fd,"%10.6f  %10.6f  %10.6f %15.10f\n", ccenter[0], ccenter[1], ccenter[2], T); 
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
}

void print_velocity (higflow_solver *ns, int myrank, int dim, int dimprint1, real pprint1, int dimprint2, real pprint2) {
    char filename[1024];
    snprintf(filename,sizeof filename,"velocity-%d-%d.dat",dim,myrank);
    FILE *fd = fopen(filename, "w");

    if (fd != NULL) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        higfit_facetiterator *fit;
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the velocity
            real u = dp_get_value(ns->dpu[dim], flid);
            if ((fabs(fcenter[dimprint1] - pprint1) < 0.5*fdelta[dimprint1]) &&
                (fabs(fcenter[dimprint2] - pprint2) < 0.5*fdelta[dimprint2])) { 
                fprintf(fd,"%10.6f  %10.6f %10.6f %15.10f\n",fcenter[0],fcenter[1],fcenter[2],u);
            }
        }
        // Destroy the iterator
        higfit_destroy(fit);
        fclose(fd);
    } else {
        printf("Arquivo %s nao aberto\n",filename);
        exit(1);
    }
}

void print_tensor_at_point(higflow_solver *ns, FILE *fd, int i, int j, real time, real pprint1, real pprint2, real pprint3) {
        // Get the cosntants
        real Re   = ns->par.Re;
        real De   = ns->ed.ve.par.De;
        real beta = ns->ed.ve.par.beta;
        real tol  = ns->ed.ve.par.kernel_tol;
        // Get the local sub-domain for the cells
        sim_domain *sdp = psd_get_local_domain(ns->ed.psdED);
        // Get the map for the domain properties
        mp_mapper *mp = sd_get_domain_mapper(sdp);
        // Loop for each cell
        higcit_celliterator *it;
        for (it = sd_get_domain_celliterator(sdp); !higcit_isfinished(it); higcit_nextcell(it)) {
            // Get the cell
            hig_cell *c = higcit_getcell(it);
            // Get the cell identifier
            int clid    = mp_lookup(mp, hig_get_cid(c));
            // Get the inside/outside inflow point cell
            Point ccenter;
            hig_get_center(c, ccenter);
            // Get the delta of the cell
            Point cdelta;
            hig_get_delta(c, cdelta);
            // Get the velocity derivative tensor Du and the Kernel tensor
            real Du[DIM][DIM];
            // Get Du
            Du[i][j] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpDu[i][j], ns->ed.stn);
            Du[j][i] = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpDu[j][i], ns->ed.stn);
            // Get T tensor
            real D  = 0.5*(Du[i][j]+Du[j][i]);
            real S = compute_value_at_point(ns->ed.sdED, ccenter, ccenter, 1.0, ns->ed.ve.dpTaup[i][j], ns->ed.stn);
            real T = S + 2.0*(1-beta)*D/Re; 
            //Print polymeric stress data file
            if ((fabs(ccenter[0] - pprint1) < 0.5*cdelta[0]) &&
                (fabs(ccenter[1] - pprint2) < 0.5*cdelta[1]) &&
                (fabs(ccenter[2] - pprint3) < 0.5*cdelta[2])) { 
                fprintf(fd,"%10.6f %10.6f %10.6f  %10.6f  %15.10f\n", time, ccenter[0], ccenter[1], ccenter[2], T); 
            }
        }
        // Destroy the iterator
        higcit_destroy(it);
}

void print_velocity_at_point (higflow_solver *ns, FILE *fd, int dim, real time, real pprint1, real pprint2, real pprint3) {
        // Get the local sub-domain
        sim_domain *sdp = psd_get_local_domain(ns->psdp);
        // Get the local partitioned domain for facets
        sim_facet_domain *sfdu[DIM];
        sfdu[dim] = psfd_get_local_domain(ns->psfdu[dim]);
        // Get the map of the distributd properties in the facets
        mp_mapper *mu = sfd_get_domain_mapper(sfdu[dim]);
        higfit_facetiterator *fit;
        // Loop for each facet
        for (fit = sfd_get_domain_facetiterator(sfdu[dim]); !higfit_isfinished(fit); higfit_nextfacet(fit)) {
            // Get the facet
            hig_facet *f = higfit_getfacet(fit);
            int flid = mp_lookup(mu, hig_get_fid(f));
            // Get the center of the facet
            Point fcenter;
            hig_get_facet_center(f, fcenter);
            // Get the delta of the facet
            Point fdelta;
            hig_get_facet_delta(f, fdelta);
            // Get the velocity
            real u = dp_get_value(ns->dpu[dim], flid);
            if ((fabs(fcenter[0] - pprint1) < 0.5*fdelta[0]) &&
                (fabs(fcenter[1] - pprint2) < 0.5*fdelta[1]) &&
                (fabs(fcenter[2] - pprint3) < 0.5*fdelta[2])) { 
                fprintf(fd,"%10.6f %10.6f  %10.6f %10.6f %15.10f\n",time, fcenter[0],fcenter[1],fcenter[2],u);
            }
        }
        // Destroy the iterator
        higfit_destroy(fit);
}
