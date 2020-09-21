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
