#include "point-mapper.h"

static void discretize_point(CPPoint x, int64_t *dx)
{
	for(unsigned dim = 0; dim < DIM; ++dim) {
		dx[dim] = x[dim] / EPSMACH + 0.5;
	}
}

static guint hash_func(gconstpointer key)
{
	static const unsigned iters = DIM * 8 / sizeof(guint);
	assert(8 % sizeof(guint) == 0);

	const guint *p = key;

	guint hash = p[0];
	for(unsigned dim = 1; dim < iters; ++dim)
	{
    		hash ^= p[dim] + 0x9e3779b9 + (hash<<6) + (hash>>2);
	}
	return hash;
}

static gboolean key_equal_func(gconstpointer a, gconstpointer b)
{
	const int64_t *va = a;
	const int64_t *vb = b;

	for(unsigned dim = 0; dim < DIM; ++dim) {
		if(va[dim] != vb[dim]) {
			return false;
		}
	}
	return true;
}

point_mapper *ptm_create(void (*value_free)(void *))
{
	point_mapper* ptm = g_hash_table_new_full(hash_func, key_equal_func,
		free, value_free
	);

	return ptm;
}

bool ptm_lookup(point_mapper *ptm, CPPoint x, void **value)
{
	int64_t key[DIM];
	discretize_point(x, key);

	return g_hash_table_lookup_extended(ptm, key,
		NULL, (gpointer *)value);
}

bool ptm_remove(point_mapper *ptm, CPPoint x)
{
	int64_t key[DIM];
	discretize_point(x, key);

	return g_hash_table_remove(ptm, key);
}

void ptm_insert(point_mapper *ptm, CPPoint x, void *value)
{
	DECL_AND_ALLOC(int64_t, key, DIM);
	discretize_point(x, key);

	g_hash_table_insert(ptm, key, (gpointer)value);
}

void ptm_destroy(point_mapper *ptm)
{
	g_hash_table_destroy(ptm);
}
