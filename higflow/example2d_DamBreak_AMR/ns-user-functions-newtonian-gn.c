#include "ns-example-2d.h"


/******************************************************************************************/
/******************************************************************************************/
/***************************** newtonian user functions ***********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial pressure
real get_pressure(Point center, real t) {
    return 0.0;
}
/*!
 * \brief Initial velocity — dam break starts from rest.
 *
 * Returns zero everywhere (the fluid is initially stationary).
 * Gravity and the hydrostatic pressure gradient drive the flow.
 */
real
get_velocity(Point center, int dim, real t)
{
    return 0.0;
}
// initial source term
real get_source_term(Point center, real t) {
    return 0.0;
}
// initial facet source term
real get_facet_source_term(Point center, int dim, real t) {
    return 0.0;
}

real get_boundary_pressure(int id, Point center, real t) {
    real value;
    switch (id) {
    case 0:
        value = 0.0;
        break;
    case 1:
        value = 0.0;
        break;
    case 2:
        value = 0.0;
        break;
    case 3:
        value = 0.0;
        break;
    }
    return value;
}

/*!
 * \brief Boundary velocity — all walls are no-slip.
 *
 * bc0 (left), bc2 (right), bc3 (bottom): u = v = 0 (Dirichlet).
 * bc1 (top): Neumann in the YAML, so this function is never called.
 */
real
get_boundary_velocity(int id, Point center, int dim, real t)
{
    (void)id;
    (void)center;
    (void)dim;
    (void)t;
    return 0.0;
}

real get_boundary_source_term(int id, Point center, real t) {
    return 0.0;
}

real get_boundary_facet_source_term(int id, Point center, int dim, real t) {
    return 0.0;
}


/******************************************************************************************/
/******************************************************************************************/
/******************** generalized newtonian user functions ********************************/
/******************************************************************************************/
/******************************************************************************************/

// initial viscosity
real get_viscosity_gn(Point center, real q, real t) {
    return 1.0;
}
