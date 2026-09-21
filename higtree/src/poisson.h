// DEAD HEADER -- kept only because nothing has decided to delete it yet.
//
// There is no poisson.c in this tree, nothing includes this file, and it could not
// compile if anything did: it has no include guard, and it names two types that do
// not exist in the library (`poisson`, and `higtree` -- the cell type is `hig_cell`).
// Every one of its doc comments is the literal word TODO.
//
// It is not a stub for planned work: the Poisson solve the system actually performs
// goes through the generic solver interface in solver.h, assembled by higflow's
// pressure step.  Do not extend this file; remove it when convenient.


//! TODO
poisson * pss_create(int order, higtree *root, mp_mapper *m);

//! TODO
void pss_set_rho(poisson *p, real *rho);

//! TODO
void pss_set_rhs(poisson *p, real *rhs);

//! TODO

void pss_set_bc(poisson *p,
