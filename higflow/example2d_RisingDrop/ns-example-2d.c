#include "ns-example-2d.h"

// 2D rising drop (cap regime) — analogue of example3d_RisingDropCap.
// A light circular drop rises through a denser ambient fluid in a [0,1]x[0,2]
// box (Y is the tall/gravity direction; the solver applies gravity on dim==1).
// Initial drop: circle of radius 0.25 centred at (0.5, 0.5).  All walls no-slip.
// Non-dimensionalised like OpenFOAM risingDrop_cap: Re=35, Ca=3.57, Fr=1.

// Signed indicator of the drop: < 0 inside the circle, > 0 outside.
real func(Point p) {
    real cx = 0.5, cy = 0.5, r = 0.25;
    real dx = p[0] - cx, dy = p[1] - cy;
    return (dx * dx + dy * dy) / (r * r) - 1.0;
}

int square_case(real f0, real f1, real f10, real f3) { return 0; }

void Intersec(Point p0, Point p1, real f0, real f1, Point p) {
    real lambda = -f0 / (f1 - f0);
    p[0] = p0[0] + lambda * (p1[0] - p0[0]);
    p[1] = p0[1] + lambda * (p1[1] - p0[1]);
}

// Volume fraction by subcell sampling (fracvol == 1 inside the drop).
real get_fracvol(Point center, Point delta, real t) {
    Point p;
    int nsub = 16;
    real vol = 0.0;
    for (int i = 0; i < nsub; i++) {
        p[0] = center[0] - 0.5 * delta[0] + (i + 0.5) * delta[0] / nsub;
        for (int j = 0; j < nsub; j++) {
            p[1] = center[1] - 0.5 * delta[1] + (j + 0.5) * delta[1] / nsub;
            if (func(p) <= 0.0) vol += 1.0;
        }
    }
    return vol / (nsub * nsub);
}

real get_pressure(Point center, real t) { return 0.0; }
real get_velocity(Point center, int dim, real t) { return 0.0; }
real get_source_term(Point center, real t) { return 0.0; }
real get_facet_source_term(Point center, int dim, real t) { return 0.0; }
real get_boundary_pressure(int id, Point center, real t) { return 0.0; }
real get_boundary_velocity(int id, Point center, int dim, real t) { return 0.0; }
real get_boundary_source_term(int id, Point center, real t) { return 0.0; }
real get_boundary_facet_source_term(int id, Point center, int dim, real t) { return 0.0; }

// Mixing rule is value = (1 - fracvol)*x0 + fracvol*x1, so phase 1 is the
// region where fracvol == 1 (the drop).  Drop = light, ambient = heavy.
real get_density0(Point center, real t) { return 1.0; }    /* ambient (heavy) */
real get_density1(Point center, real t) { return 0.001; }  /* drop (light)    */
real get_viscosity0(Point center, real t) { return 1.0; }  /* ambient         */
real get_viscosity1(Point center, real t) { return 0.01; } /* drop            */

int main(int argc, char *argv[]) {
    START_CLOCK(total);
    int ntasks, myrank;
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    higflow_solver *ns = higflow_create();
    higflow_load_data_file_names(argc, argv, ns);
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
    higflow_set_external_functions(ns, get_pressure, get_velocity,
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term);
    int order_center = 2;
    int order_facet = 2;
    int cache = 1;
    higflow_create_domain(ns, cache, order_center);
    higflow_create_domain_multiphase(ns, cache, order_center,
        get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet);
    higflow_initialize_boundaries_yaml(ns);
    higflow_create_distributed_properties(ns);
    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);
    higflow_create_solver(ns);
    if (ns->par.step > 0) {
        higflow_load_properties(ns, myrank, ntasks);
    }
    higflow_save_domain_yaml(ns, myrank, ntasks);
    higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
    if (ns->par.step == 0) {
        higflow_print_vtk(ns, myrank);
        ns->par.tp += ns->par.dtp;
        ns->par.frame++;
        higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }
    DECL_CLOCK(firstiter);
    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
        if (ns->par.step == step0) START_CLOCK(firstiter);
        higflow_solver_step_multiphase(ns);
        ns->par.t += ns->par.dt;
        if (ns->par.step == step0) STOP_CLOCK(firstiter);
        if (ns->par.t >= ns->par.tp) {
            higflow_print_vtk(ns, myrank);
            ns->par.tp += ns->par.dtp;
            ns->par.frame++;
        }
        if (ns->par.t >= ns->par.ts) {
            higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
            higflow_save_properties(ns, myrank, ntasks);
            ns->par.ts += ns->par.dts;
        }
    }
    higflow_destroy(ns);
    STOP_CLOCK(total);
    return 0;
}
