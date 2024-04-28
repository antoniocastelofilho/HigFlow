// *******************************************************************
// *******************************************************************
//  HiG-Flow Solver Residuals - version 13/05/2024
// *******************************************************************
// *******************************************************************

#ifndef HIG_FLOW_RES
#define HIG_FLOW_RES

#include "hig-flow-kernel.h"

#define RES_STORE_NUM_BASE 5
#define RES_STORE_NUM_POW 3
#define RES_STORE_NUM 125 //access to the last num residuals

#define NUM_TIMEINFOS 3
#define NUM_TIMESTATS 4

#define COMPUTE_MIDRANGE // channel only
#define COMPUTE_MIDLINE // channel only

#ifdef COMPUTE_RESIDUALS
    #define UPDATE_RESIDUAL_BUFFER_FACET(ns, val1, val2, f, fcenter) \
        do { \
            update_residual_buffer_facet((ns)->residuals->res_buffer, (val1), (val2), (f), (fcenter)); \
        } while (0);
#else
    #define UPDATE_RESIDUAL_BUFFER_FACET(ns, val1, val2, f, fcenter) do {} while (0);
#endif // COMPUTE_RESIDUALS

#ifdef COMPUTE_RESIDUALS
    #define UPDATE_RESIDUAL_BUFFER_CELL(ns, val1, val2, c, ccenter) \
        do { \
            update_residual_buffer_cell((ns)->residuals->res_buffer, (val1), (val2), (c), (ccenter)); \
        } while (0);
#else
    #define UPDATE_RESIDUAL_BUFFER_CELL(ns, val1, val2, c, ccenter) do {} while (0);
#endif // COMPUTE_RESIDUALS

#ifdef COMPUTE_RESIDUALS
    #define UPDATE_RESIDUALS(ns, property) \
    do { \
        int myrank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank); \
        reduce_residual_buffer((ns)->residuals->res_buffer); \
        if (myrank == 0) \
            update_dp_residuals(property, (ns)->residuals->res_buffer, (ns)->par.step, (ns)->par.initstep); \
        clear_residual_buffer((ns)->residuals->res_buffer); \
    } while (0);
#else
    #define UPDATE_RESIDUALS(ns, property) do {} while (0);
#endif // COMPUTE_RESIDUALS


typedef struct residual_normgroup_buffer{
    real res_max;
    real res_2;
    real res_1;
    real area;
} residual_normgroup_buffer;

typedef struct residual_buffer{
    #ifdef COMPUTE_MIDLINE
        residual_normgroup_buffer midline; // for channels only
        real midlinex;
    #endif // COMPUTE_MIDLINE
    #ifdef COMPUTE_MIDRANGE
        residual_normgroup_buffer midrange; // for channels only
        real leftx;
        real rightx;
    #endif // COMPUTE_MIDRANGE
    residual_normgroup_buffer all;
}residual_buffer;

typedef struct residual_timeinfo{
    real stored[RES_STORE_NUM]; // last residual, last base residuals, last base^2 residuals and so on;
    real avg[RES_STORE_NUM_POW]; // average of last base residuals, average of last base^2 residuals and so on;
    real geoavg[RES_STORE_NUM_POW]; // geometric average of last base residuals, geometric average of last 100 residuals and so on;
    //real med[RES_STORE_NUM_POW]; // median of last base residuals, median of last base^2 residuals and so on;
    real max[RES_STORE_NUM_POW]; // max of last base residuals, max of last base^2 residuals and so on;
    real min[RES_STORE_NUM_POW]; // min of last base residuals, min of last base^2 residuals and so on;
    real last; // last residual
} residual_timeinfo;

typedef struct residual_normgroup{
    residual_timeinfo *res_max;
    residual_timeinfo *res_2;
    residual_timeinfo *res_1;
} residual_normgroup;

typedef struct dp_residuals{
    #ifdef COMPUTE_MIDLINE
        residual_normgroup *midline; // for channels only
    #endif // COMPUTE_MIDLINE
    #ifdef COMPUTE_MIDRANGE
        residual_normgroup *midrange; // for channels only
    #endif // COMPUTE_MIDRANGE
    residual_normgroup *all;
}dp_residuals;

struct sim_residuals{
    // dp_residuals *u_star[DIM];
    dp_residuals *u[DIM];
    // dp_residuals *p;

    // dp_residuals *D[DIM][DIM];
    // dp_residuals *S[DIM][DIM];
    dp_residuals *Kernel[DIM][DIM];
    // dp_residuals *B0_integral[DIM][DIM];

    // dp_residuals *phi;
    dp_residuals *psi;
    // dp_residuals *nplus;
    // dp_residuals *nminus;
    // dp_residuals *Feo[DIM];

    dp_residuals *fracvol_adv[DIM];
    // dp_residuals *curvature;
    // dp_residuals *normal[DIM];
    // dp_residuals *IF[DIM];

    residual_buffer *res_buffer;
};

// extern sim_residuals *(RES_VAR_NAME); // has to be initialized in the main program using create_sim_residuals
//                                         // and freed using free_sim_residuals

sim_residuals *create_initialize_sim_residuals(higflow_solver *ns);
void free_sim_residuals(sim_residuals *sim_res);

// void get_nameres(char **source, char *destination);
// write the residuals files
void write_residuals(sim_residuals *sim_res, higflow_solver *ns);

// the routines below are called to update the residual buffer according to the type of iterator
void update_residual_buffer_facet(residual_buffer *res_buff, real val1, real val2, hig_facet *f, Point fcenter);
void update_residual_buffer_cell(residual_buffer *res_buff, real val1, real val2, hig_cell *c, Point ccenter);

/*
  the routines below are called to update the residuals after the iterator finishes
  update_dp_residuals is ONLY called from process 0
*/
void reduce_residual_buffer(residual_buffer *res_buff); // called first to gather info from all processes
void update_dp_residuals(dp_residuals *dp_res, residual_buffer *res_buff, int step, int initstep); // only called to update the residuals structure in process 0
void clear_residual_buffer(residual_buffer *res_buff); // called after update_dp_residuals

#endif // HIG_FLOW_RES