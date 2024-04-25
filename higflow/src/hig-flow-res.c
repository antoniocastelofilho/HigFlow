#include "hig-flow-res.h"

// sim_residuals *residuals_global; // has to be initialized in the main program using create_sim_residuals
//                                  // and freed using free_sim_residuals


residual_timeinfo *create_residual_timeinfo() {
    residual_timeinfo *res = (residual_timeinfo *)calloc(1, sizeof(residual_timeinfo));
    return res;
}

void free_residual_timeinfo(residual_timeinfo *res) {
    free(res);
    res = NULL;
}

residual_normgroup *create_residual_normgroup() {
    residual_normgroup *resn = (residual_normgroup *)malloc(sizeof(residual_normgroup));
    resn->res_max = create_residual_timeinfo();
    resn->res_2 = create_residual_timeinfo();
    resn->res_1 = create_residual_timeinfo();
    return resn;
}

void free_residual_normgroup(residual_normgroup *resn) {
    free_residual_timeinfo(resn->res_max);
    free_residual_timeinfo(resn->res_2);
    free_residual_timeinfo(resn->res_1);
    free(resn);
    resn = NULL;
}

dp_residuals *create_dp_residuals() {
    dp_residuals *dp_res = (dp_residuals *)malloc(sizeof(dp_residuals));
    dp_res->all = create_residual_normgroup();
    #ifdef COMPUTE_MIDRANGE
        dp_res->midrange = create_residual_normgroup();
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        dp_res->midline = create_residual_normgroup();
    #endif // COMPUTE_MIDLINE
    return dp_res;
}

void free_dp_residuals(dp_residuals *dp_res) {
    free_residual_normgroup(dp_res->all);
    #ifdef COMPUTE_MIDRANGE
        free_residual_normgroup(dp_res->midrange);
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        free_residual_normgroup(dp_res->midline);
    #endif // COMPUTE_MIDLINE
    free(dp_res);
    dp_res = NULL;
}

sim_residuals *create_sim_residuals(higflow_controllers hig_contr, mult_controllers mult_contr, int myrank) {
    sim_residuals *sim_res = (sim_residuals *)malloc(sizeof(sim_residuals));
    if(myrank==0) {
        ////////////////// Make dps NULL by default //////////////////
        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                sim_res->Kernel[dim][dim2] = NULL;
            }
        }

        sim_res->psi = NULL;

        for (int dim=0; dim<DIM; dim++){
            sim_res->fracvol_adv[dim] = NULL;
        }
        ///////////////////////////////////////////////////////////////

        for(int dim=0; dim<DIM; dim++){
            sim_res->u[dim] = create_dp_residuals();
        }
        
        switch(hig_contr.flowtype) {
            case MULTIPHASE:
                for(int dim=0; dim<DIM; dim++){
                    sim_res->fracvol_adv[dim] = create_dp_residuals();
                }

                if(mult_contr.flowtype_either == VISCOELASTIC) {
                    for(int dim=0; dim<DIM; dim++){
                        for(int dim2=dim; dim2<DIM; dim2++){
                            sim_res->Kernel[dim][dim2] = create_dp_residuals();
                        }
                    }
                }
                break;
            case VISCOELASTIC:
                for(int dim=0; dim<DIM; dim++){
                    for(int dim2=dim; dim2<DIM; dim2++){
                        sim_res->Kernel[dim][dim2] = create_dp_residuals();
                    }
                }
                break;
        }

        if(hig_contr.eoflow == true) {
            sim_res->psi = create_dp_residuals();
        }

    }
    sim_res->res_buffer = (residual_buffer *)calloc(1, sizeof(residual_buffer));
    return sim_res;
}

void free_sim_residuals(sim_residuals *sim_res) {
    if(sim_res == NULL) return;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank==0) {
        for(int dim=0; dim<DIM; dim++){
            free_dp_residuals(sim_res->u[dim]);
        }

        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                if(sim_res->Kernel[dim][dim2] != NULL) {
                    free_dp_residuals(sim_res->Kernel[dim][dim2]);
                }
            }
        }

        if(sim_res->psi != NULL) {
            free_dp_residuals(sim_res->psi);
        }

        for(int dim=0; dim<DIM; dim++){
            if(sim_res->fracvol_adv[dim] != NULL) {
                free_dp_residuals(sim_res->fracvol_adv[dim]);
            }
        }
    }
    free(sim_res->res_buffer);
    free(sim_res);
    sim_res = NULL;
}

void get_channel_lengths(Point center, Point L, char *nameload) {
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
    if(numhigs!=1) printf("Number of HigTrees: %d ==============> EXPECTED TO BE 1\n", numhigs);
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
        center[dim] = (h[dim] + l[dim])/2.0;
        L[dim] = h[dim] - l[dim];
	}
    // Close the AMR format file
    fclose(fd);
    fclose(fdomain);
}


sim_residuals *create_initialize_sim_residuals(higflow_solver *ns) {
    sim_residuals *sim_res = NULL;
    #ifdef COMPUTE_RESIDUALS
        higflow_controllers hig_contr = ns->contr;
        mult_controllers mult_contr;
        if(ns->contr.flowtype == MULTIPHASE) {
            mult_contr = ns->ed.mult.contr;
        }
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        sim_res = create_sim_residuals(hig_contr, mult_contr, myrank);

        /////////////////// get channel dimensions if needed ///////////////
        #if defined(COMPUTE_MIDRANGE) || defined(COMPUTE_MIDLINE)
            Point center; Point L;
            get_channel_lengths(center, L, ns->par.nameload);
            #ifdef COMPUTE_MIDRANGE
                real leftx = center[0] - 0.25*L[0];
                real rightx = center[0] + 0.25*L[0];
                sim_res->res_buffer->leftx = leftx;
                sim_res->res_buffer->rightx = rightx;
            #endif // COMPUTE_MIDRANGE
            #ifdef COMPUTE_MIDLINE
                real midlinex = center[0];
                sim_res->res_buffer->midlinex = midlinex;
            #endif // COMPUTE_MIDLINE
        #endif // COMPUTE_MIDRANGE
    #endif // COMPUTE_RESIDUALS

    ns->residuals = sim_res;

    return sim_res;
}


void update_residual_timeinfo(residual_timeinfo *res, real r, int step, int initstep) {

    int rstep = step - initstep;

    if(rstep==0){
        for(int i=0; i<RES_STORE_NUM_POW; i++){
            res->avg[i] = r;
            res->geoavg[i] = r;
            //res->med[i] = r;
            res->max[i] = r;
            res->min[i] = r;
        }
        res->stored[RES_STORE_NUM-1] = r;
        res->last = r;
        return;
    }
    
    int num_last = RES_STORE_NUM_BASE;
    real popped;
    for(int i=0; i<RES_STORE_NUM_POW; i++){
        if(rstep < num_last) { //does not yet have base, base^2... steps necessary to calculate avg
            res->avg[i] = (res->avg[i]*rstep + r)/(rstep+1);
            res->geoavg[i] = pow(res->geoavg[i], rstep/(rstep+1)) * pow(r, 1.0/(rstep+1));
            if(r > res->max[i]) res->max[i] = r;
            if(r < res->min[i]) res->min[i] = r;
        }
        else{
            popped = res->stored[RES_STORE_NUM-num_last];
            res->avg[i] -= popped/num_last;
            res->avg[i] += r/num_last;
            if(popped!=0.0) res->geoavg[i] /= pow(popped, 1.0/num_last);
            res->geoavg[i] *= pow(r, 1.0/num_last);
            if(popped < res->max[i]) {
                if(r > res->max[i]) res->max[i] = r;
            }
            else { // current maximum is being popped - finds new maximum
                res->max[i] = r;
                for (int j=RES_STORE_NUM-num_last+1; j<RES_STORE_NUM; j++)
                    if(res->stored[j] > res->max[i]) res->max[i] = res->stored[j];
            }
            if(popped > res->min[i]) {
                if(r < res->min[i]) res->min[i] = r;
            }
            else { // current minimum is being popped - finds new minimum
                res->min[i] = r;
                for (int j=RES_STORE_NUM-num_last+1; j<RES_STORE_NUM; j++)
                    if(res->stored[j] < res->min[i]) res->min[i] = res->stored[j];
            }
        }
        num_last*=RES_STORE_NUM_BASE;
    }

    for (int j=RES_STORE_NUM-min(RES_STORE_NUM-1,rstep); j<RES_STORE_NUM; j++)
        res->stored[j-1] = res->stored[j];
    res->stored[RES_STORE_NUM-1] = r;
    res->last = r;
    return;
}

void reduce_residual_normgroup_buffer(residual_normgroup_buffer *resn_buff) {
    real rmax, r2, r1, area;
    real rmax_global, r2_global, r1_global, area_global;

    rmax = resn_buff->res_max;
    r2 = resn_buff->res_2;
    r1 = resn_buff->res_1;
    area = resn_buff->area;
    MPI_Allreduce(&rmax, &rmax_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&area, &area_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&r2, &r2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r2_global = sqrt(r2_global) / area_global; // vector 2-norm
    MPI_Allreduce(&r1, &r1_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r1_global = r1_global / area_global; // vector 1-norm

    resn_buff->res_max = rmax_global;
    resn_buff->res_2 = r2_global;
    resn_buff->res_1 = r1_global;
    resn_buff->area = area_global;
}

void reduce_residual_buffer(residual_buffer *res_buff) {
    
    reduce_residual_normgroup_buffer(&res_buff->all);

    #ifdef COMPUTE_MIDRANGE
        reduce_residual_normgroup_buffer(&res_buff->midrange);
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        reduce_residual_normgroup_buffer(&res_buff->midline);
    #endif // COMPUTE_MIDLINE
}

void update_residual_normgroup(residual_normgroup *resn, residual_normgroup_buffer *resn_buff, int step, int initstep) {
    real rmax = resn_buff->res_max;
    real r2 = resn_buff->res_2;
    real r1 = resn_buff->res_1;

    update_residual_timeinfo(resn->res_max, rmax, step, initstep);
    update_residual_timeinfo(resn->res_2, r2, step, initstep);
    update_residual_timeinfo(resn->res_1, r1, step, initstep);
}

void update_dp_residuals(dp_residuals *dp_res, residual_buffer *res_buff, int step, int initstep) {
    update_residual_normgroup(dp_res->all, &res_buff->all, step, initstep);
    #ifdef COMPUTE_MIDRANGE
        update_residual_normgroup(dp_res->midrange, &res_buff->midrange, step, initstep);
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        update_residual_normgroup(dp_res->midline, &res_buff->midline, step, initstep);
    #endif // COMPUTE_MIDLINE
}

void clear_residual_normgroup_buffer(residual_normgroup_buffer *resn_buff) {
    resn_buff->res_max = 0.0;
    resn_buff->res_2 = 0.0;
    resn_buff->res_1 = 0.0;
    resn_buff->area = 0.0;
}

void clear_residual_buffer(residual_buffer *res_buff) {
    clear_residual_normgroup_buffer(&res_buff->all);
    #ifdef COMPUTE_MIDRANGE
        clear_residual_normgroup_buffer(&res_buff->midrange);
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        clear_residual_normgroup_buffer(&res_buff->midline);
    #endif // COMPUTE_MIDLINE
}

void update_residual_buffer_facet(residual_buffer *res_buff, real val1, real val2, hig_facet *f, Point fcenter){
    Point cdelta;
    hig_cell *cell_with_facet = hig_get_facet_cell(f);
    hig_get_delta(cell_with_facet, cdelta);
    real cell_area = cdelta[0]*cdelta[1];

    real dif = fabs(val1 - val2);
    real difc = dif*cell_area;
    real dif2c = dif*dif*cell_area;

    res_buff->all.area += cell_area;
    if(dif > res_buff->all.res_max) res_buff->all.res_max = dif;
    res_buff->all.res_1 += difc;
    res_buff->all.res_2 += dif2c;


    #ifdef COMPUTE_MIDRANGE
    if(POS_GE(fcenter[0], res_buff->leftx) && POS_LE(fcenter[0], res_buff->rightx)) { //midrange
        res_buff->midrange.area += cell_area;
        if(dif > res_buff->midrange.res_max) res_buff->midrange.res_max = dif;
        res_buff->midrange.res_1 += difc;
        res_buff->midrange.res_2 += dif2c;
    }
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        Point lowpoint, highpoint;
        hig_get_lowpoint(cell_with_facet, lowpoint);
        hig_get_highpoint(cell_with_facet, highpoint);
        if(POS_GE(res_buff->midlinex, lowpoint[0]) && POS_LE(res_buff->midlinex, highpoint[0])){ // midline
            res_buff->midline.area += cell_area;
            if(dif > res_buff->midline.res_max) res_buff->midline.res_max = dif;
            res_buff->midline.res_1 += difc;
            res_buff->midline.res_2 += dif2c;
        }
    #endif // COMPUTE_MIDLINE
}

void update_residual_buffer_cell(residual_buffer *res_buff, real val1, real val2, hig_cell *c, Point ccenter){
    Point cdelta;
    hig_get_delta(c, cdelta);
    real cell_area = cdelta[0]*cdelta[1];

    real dif = fabs(val1 - val2);
    real difc = dif*cell_area;
    real dif2c = dif*dif*cell_area;

    res_buff->all.area += cell_area;
    if(dif > res_buff->all.res_max) res_buff->all.res_max = dif;
    res_buff->all.res_1 += difc;
    res_buff->all.res_2 += dif2c;

    #ifdef COMPUTE_MIDRANGE
    if(POS_GE(ccenter[0], res_buff->leftx) && POS_LE(ccenter[0], res_buff->rightx)) { //midrange
        res_buff->midrange.area += cell_area;
        if(dif > res_buff->midrange.res_max) res_buff->midrange.res_max = dif;
        res_buff->midrange.res_1 += difc;
        res_buff->midrange.res_2 += dif2c;
    }
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        Point lowpoint, highpoint;
        hig_get_lowpoint(c, lowpoint);
        hig_get_highpoint(c, highpoint);
        if(POS_GE(res_buff->midlinex, lowpoint[0]) && POS_LE(res_buff->midlinex, highpoint[0])){ // midline
            res_buff->midline.area += cell_area;
            if(dif > res_buff->midline.res_max) res_buff->midline.res_max = dif;
            res_buff->midline.res_1 += difc;
            res_buff->midline.res_2 += dif2c;
        }
    #endif // COMPUTE_MIDLINE
}

void fill_column_string(const char *literal, char *str, size_t size) {
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

/*void get_nameres(char **source, char *destination) {
    // Find the second to last '/' character
    char *second_last_slash = NULL;
    char *last_slash = NULL;
    char *current = *source;
    
    while (*current != '\0') {
        if (*current == '/') {
            second_last_slash = last_slash;
            last_slash = current;
        }
        current++;
    }

    // Copy the portion of the string up to the second to last '/' character
    if (second_last_slash != NULL) {
        strncpy(destination, *source, second_last_slash - *source);
        destination[second_last_slash - *source] = '\0';
    } else {
        strcpy(destination, *source); // No '/' found, copy the entire string
    }

    // Append "/res/" to the destination string
    strcat(destination, "/res/");
}*/


void write_residual_normgroup(residual_normgroup *resn, char *nameres_dp, char *resn_name, real t, int step) {
    char nameres_dp_resn[1024];
    strcpy(nameres_dp_resn, nameres_dp);
    strcat(nameres_dp_resn, resn_name);
    strcat(nameres_dp_resn,".txt");

    int num_size = 16;
    int num_dec_size = 12;
    double num_max = pow(10, num_size - num_dec_size - 1) - pow(10, -(num_dec_size - 1));
    int num_space = 4;
    int str_col_size = num_size + num_space;
    char curr_col[str_col_size+1];

    // #define NUM_TIMEINFOS_WRITE NUM_TIMEINFOS;
    // #define NUM_TIMESTATS_WRITE NUM_TIMESTATS;
    
    /* opted to reduce the information written */
    #define NUM_TIMEINFOS_WRITE 1
    #define NUM_TIMESTATS_WRITE 2
    

    FILE *f;
    
    if(step == 0) {
        char curr_name[str_col_size+1];
        f = fopen(nameres_dp_resn, "w");

        fill_column_string("step", curr_col, str_col_size);
        fprintf(f, "%s", curr_col);
        fill_column_string("t", curr_col, str_col_size);
        fprintf(f, "%s", curr_col);

        // char normgroup_names[NUM_TIMEINFOS_WRITE][8] = {"maxnorm", "2norm", "1norm"};
        // char timeinfo_names[NUM_TIMESTATS_WRITE][5] = {"avg", "gavg", "max", "min"};

        /* opted to reduce the information written */
        char normgroup_names[NUM_TIMEINFOS_WRITE][8] = {"maxnorm"};
        char timeinfo_names[NUM_TIMESTATS_WRITE][5] = {"avg", "max"};

        // header
        for(int i=0; i<NUM_TIMEINFOS_WRITE; i++){
            sprintf(curr_name, "%s", normgroup_names[i]);
            fill_column_string(curr_name, curr_col, str_col_size);
            fprintf(f, "%s", curr_col);
            for(int j=0; j<NUM_TIMESTATS_WRITE; j++){
                int num_last = RES_STORE_NUM_BASE;
                for(int k=0; k<RES_STORE_NUM_POW; k++){
                    sprintf(curr_name, "%s_%s%d", normgroup_names[i], timeinfo_names[j], num_last);
                    fill_column_string(curr_name, curr_col, str_col_size);
                    fprintf(f, "%s", curr_col);
                    num_last*=RES_STORE_NUM_BASE;
                }
            }
        }
        fprintf(f, "\n");

        fclose(f);
    }

    f = fopen(nameres_dp_resn, "a");
    
    char format[10];
    sprintf(format, "%%%d.%dlf", num_size, num_dec_size);
    char curr_num[str_col_size+1];
    real num;

    sprintf(curr_num, "%d", step);
    fill_column_string(curr_num, curr_col, str_col_size);
    fprintf(f, "%s", curr_col);
    sprintf(curr_num, "%lf", t);
    fill_column_string(curr_num, curr_col, str_col_size);
    fprintf(f, "%s", curr_col);

    // residual_timeinfo *timeinfos[NUM_TIMEINFOS_WRITE] = {resn->res_max, resn->res_2, resn->res_1};

    /* opted to reduce the information written */
    residual_timeinfo *timeinfos[NUM_TIMEINFOS_WRITE] = {resn->res_max};
    for(int i=0; i<NUM_TIMEINFOS_WRITE; i++){
        residual_timeinfo *timeinfo = timeinfos[i];
        num = min(timeinfo->last, num_max);
        sprintf(curr_num, format, num);
        fill_column_string(curr_num, curr_col, str_col_size);
        fprintf(f, "%s", curr_col);

        // real *timestats[NUM_TIMESTATS_WRITE] = {timeinfo->avg, timeinfo->geoavg, timeinfo->max, timeinfo->min};

        /* opted to reduce the information written */
        real *timestats[NUM_TIMESTATS_WRITE] = {timeinfo->avg, timeinfo->max};
        for(int j=0; j<NUM_TIMESTATS_WRITE; j++){
            real *timestat_over_time  = timestats[j];
            int num_last = RES_STORE_NUM_BASE;
            for(int k=0; k<RES_STORE_NUM_POW; k++){
                num = min(timestat_over_time[k], num_max);
                sprintf(curr_num, format, num);
                fill_column_string(curr_num, curr_col, str_col_size);
                fprintf(f, "%s", curr_col);
                num_last*=RES_STORE_NUM_BASE;
            }
        }
    }
    fprintf(f, "\n");

    fclose(f);

}


void write_dp_residuals(dp_residuals *dp_res, char *nameres, char *dp_name, real t, int step) {
    char nameres_dp[1024];
    strcpy(nameres_dp, nameres);
    strcat(nameres_dp, dp_name);
    write_residual_normgroup(dp_res->all, nameres_dp, "", t, step);
    #ifdef COMPUTE_MIDRANGE
        write_residual_normgroup(dp_res->midrange, nameres_dp, "_midrange", t, step);
    #endif // COMPUTE_MIDRANGE
    #ifdef COMPUTE_MIDLINE
        write_residual_normgroup(dp_res->midline, nameres_dp, "_midline", t, step);
    #endif // COMPUTE_MIDLINE
}


void write_residuals(sim_residuals *sim_res, higflow_solver *ns) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int step = ns->par.step;
    double t = ns->par.t;
    char *nameres = ns->par.nameres;
    if(myrank == 0){
        char dp_name[64];
        for(int dim=0; dim<DIM; dim++){
            sprintf(dp_name, "u%d", dim);
            write_dp_residuals(sim_res->u[dim], nameres, dp_name, t, step);
        }

        for(int dim=0; dim<DIM; dim++){
            for(int dim2=dim; dim2<DIM; dim2++){
                if(sim_res->Kernel[dim][dim2] != NULL){
                    sprintf(dp_name, "Kernel%d%d", dim, dim2);
                    write_dp_residuals(sim_res->Kernel[dim][dim2], nameres, dp_name, t, step);
                }
            }
        }

        if(sim_res->psi != NULL){
            sprintf(dp_name, "psi");
            write_dp_residuals(sim_res->psi, nameres, dp_name, t, step);
        }
        
        for(int dim=0; dim<DIM; dim++){
            if(sim_res->fracvol_adv[dim] != NULL){
                sprintf(dp_name, "fracvoladv%d", dim);
                write_dp_residuals(sim_res->fracvol_adv[dim], nameres, dp_name, t, step);
            }
        }
        
    }
}