#ifndef RNG_H
#define RNG_H

#include "coord.h"

typedef struct RNG RNG;

//! Creates a random number generator from the given seed.
// If seed_size is 0, the seed is random.
RNG *rng_create(const char* seed, size_t seed_size);

//! Generates a value uniformly distributed in the range [a, b).
real rng_uniform(RNG* rng, real a, real b);

//! Destroys the random number generator.
void rng_destroy(RNG* rng);

#endif
