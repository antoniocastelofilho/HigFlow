// The scalar type of the whole system, and its MPI counterpart.
//
// `real` and `MPI_HIGREAL` MUST be changed together.  Narrowing `real` to float and
// leaving MPI_HIGREAL as MPI_DOUBLE compiles cleanly and then corrupts every halo
// exchange in the run, because each side reads twice the bytes it was sent.  There
// is no build-time check for this; the pairing is a convention held by hand.

#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>

//! Defines the size of a real value.
typedef double real;

//! Defines a boolean.
typedef bool boolean;

//! Defines a vector of integers.
typedef int     *vint;

//! Defines a vector of reals.
typedef real    *vreal;

//! Defines a matrix of reals.
typedef real   **mreal;

#define MPI_HIGREAL MPI_DOUBLE


#endif
