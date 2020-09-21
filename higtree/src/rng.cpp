#include <random>

extern "C" {
#include "rng.h"
}

struct RNG
{
	//! Mersenne Twistter engine.
	std::mt19937_64 mt;
};

RNG *rng_create(const char* seed, size_t seed_size)
{
	RNG* ret = new RNG;

	if(seed_size) {
		std::seed_seq ss(seed, seed + seed_size);
		ret->mt.seed(ss);
	} else {
		std::random_device seeder;
		std::vector<uint64_t> seed(ret->mt.state_size);
		for(uint64_t &s: seed) {
			s = seeder();
		}

		std::seed_seq ss(seed.begin(), seed.end());
		ret->mt.seed(ss);
	}

	return ret;
}

real rng_uniform(RNG* rng, real a, real b)
{
	std::uniform_real_distribution<real> dist(a, b);
	return dist(rng->mt);
}

//! Destroys the random number generator.
void rng_destroy(RNG* rng)
{
	delete rng;
}
