#include "ns-exemple-3d.h"
#include <string.h>

real func(Point p) {
    real cx = 0.5, cy = 0.5, cz = 0.5, r = 0.25;
    real dx = p[0] - cx, dy = p[1] - cy, dz = p[2] - cz;
    return (dx * dx + dy * dy + dz * dz) / (r * r) - 1.0;
}

int square_case(real f0, real f1, real f10, real f3) { return 0; }

void Intersec(Point p0, Point p1, real f0, real f1, Point p) {
    real lambda = -f0/(f1-f0);
    p[0] = p0[0] + lambda*(p1[0]-p0[0]);
    p[1] = p0[1] + lambda*(p1[1]-p0[1]);
    p[2] = p0[2] + lambda*(p1[2]-p0[2]);
}

real get_fracvol(Point center, Point delta, real t) {
    Point p; int nsub = 8; real vol = 0.0;
    for (int i = 0; i < nsub; i++) {
        p[0] = center[0] - 0.5*delta[0] + (i+0.5)*delta[0]/nsub;
        for (int j = 0; j < nsub; j++) {
            p[1] = center[1] - 0.5*delta[1] + (j+0.5)*delta[1]/nsub;
            for (int k = 0; k < nsub; k++) {
                p[2] = center[2] - 0.5*delta[2] + (k+0.5)*delta[2]/nsub;
                if (func(p) <= 0.0) vol += 1.0;
            }
        }
    }
    return vol/(nsub*nsub*nsub);
}

real get_pressure(Point center, real t) { return 0.0; }
real get_velocity(Point center, int dim, real t) { return 0.0; }
real get_source_term(Point center, real t) { return 0.0; }
real get_facet_source_term(Point center, int dim, real t) { return 0.0; }
real get_boundary_pressure(int id, Point center, real t) { return 0.0; }
real get_boundary_velocity(int id, Point center, int dim, real t) { return 0.0; }
real get_boundary_source_term(int id, Point center, real t) { return 0.0; }
real get_boundary_facet_source_term(int id, Point center, int dim, real t) { return 0.0; }
/* Densidade/viscosidade seguem regra de mistura:
 *   valor = (1 - fracvol)*x0 + fracvol*x1
 * Phase 0 = ambiente (pesado, fracvol==0)
 * Phase 1 = bolha (leve, fracvol==1)
 * Colleague's values: ambient rho=1000 mu=10; bubble rho=100 mu=1
 * Non-dimensionalizados pelo ambiente: razao 10:1 */
real get_density0(Point center, real t) { return 1.0; }   /* ambiente (pesado) */
real get_density1(Point center, real t) { return 0.1; }   /* bolha (leve)      */
real get_viscosity0(Point center, real t) { return 1.0; } /* ambiente          */
real get_viscosity1(Point center, real t) { return 0.1; } /* bolha             */

/* Viscoelastic kernel functions (Oldroyd-B: identity transformations) */
real get_tensor_multiphase(real fracvol, Point center, int i, int j, real t) { return 0.0; }
real get_kernel(int dim, real lambda, real tol) { return lambda; }
real get_kernel_inverse(int dim, real lambda, real tol) { return lambda; }
real get_kernel_jacobian(int dim, real lambda, real tol) { return 1.0; }

int main(int argc, char *argv[]) {

    START_CLOCK(total);
    int ntasks, myrank;
    higflow_initialize(&argc, &argv, &myrank, &ntasks);
    higflow_solver *ns = higflow_create();
    memset(ns, 0, sizeof(higflow_solver));
    higflow_load_data_file_names(argc, argv, ns);
    higflow_load_all_controllers_and_parameters_yaml(ns, myrank);
    higflow_set_external_functions(ns, get_pressure, get_velocity,
        get_source_term, get_facet_source_term,
        get_boundary_pressure, get_boundary_velocity,
        get_boundary_source_term, get_boundary_facet_source_term);
    int order_center = 2, order_facet = 2, cache = 1;
    higflow_create_domain(ns, cache, order_center);
    higflow_create_domain_multiphase(ns, cache, order_center,
        get_viscosity0, get_viscosity1, get_density0, get_density1, get_fracvol);
    higflow_create_domain_multiphase_viscoelastic(ns,
        get_tensor_multiphase, get_kernel, get_kernel_inverse, get_kernel_jacobian);
    higflow_initialize_domain_yaml(ns, ntasks, myrank, order_facet);
    higflow_initialize_boundaries_yaml(ns);

    higflow_create_distributed_properties(ns);

    if (ns->par.step == 0) higflow_initialize_distributed_properties(ns);

    higflow_create_solver(ns);

    if (ns->par.step > 0) higflow_load_properties(ns, myrank, ntasks);

    higflow_save_domain_yaml(ns, myrank, ntasks);
    higflow_save_all_boundaries_yaml(ns, myrank, ntasks);
    higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
    if (ns->par.step == 0) {
        if (ntasks > 1) higflow_print_vtk3D_parallel_single(ns, myrank, ntasks); else higflow_print_vtk(ns, myrank);
        ns->par.tp += ns->par.dtp; ns->par.frame++;
        higflow_save_all_controllers_and_parameters_yaml(ns, myrank);
        higflow_save_properties(ns, myrank, ntasks);
        ns->par.ts += ns->par.dts;
    }

    DECL_CLOCK(firstiter);
    for (int step0 = ns->par.initstep; ns->par.step <= ns->par.finalstep; ns->par.step++) {
        if (ns->par.step == step0) START_CLOCK(firstiter);
        higflow_solver_step_multiphase_viscoelastic(ns);
        ns->par.t += ns->par.dt;
        if (ns->par.step == step0) STOP_CLOCK(firstiter);

        if (ns->par.t >= ns->par.tp) {
            if (ntasks > 1) higflow_print_vtk3D_parallel_single(ns, myrank, ntasks); else higflow_print_vtk(ns, myrank);
            ns->par.tp += ns->par.dtp; ns->par.frame++;
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
