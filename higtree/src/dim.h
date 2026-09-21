// DIM is fixed at COMPILE TIME, and the library is built once per value: the tree
// layout, Point, the facet count per cell and the solver's stride are all baked in.
//
// That is why there is a libhig2d.a and a libhig3d.a, and why an object compiled at
// one DIM linked against the library of the other produces garbage rather than a
// link error -- the symbols match, the structure sizes do not.  Anything that builds
// against HiGTree must pass the same -DDIM as the library it links.

#ifndef __DIM_H
#define __DIM_H

#ifndef DIM
//! Defines the dimension of the model. It can be defined by commandline as well (e.g. using DIM=3 the make invocation)
#define DIM 2
#endif


#endif
