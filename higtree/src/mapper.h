#ifndef MAPPER_H
#define MAPPER_H

#include "uniqueid.h"
#include "higtree-iterator.h"
#include <stdint.h>
#include <glib.h>

typedef uint32_t mp_value_t;

typedef void * mpp_value_t;

#define MP_UNDEF ((mp_value_t)(GPOINTER_TO_INT(NULL)-1))
#define MP_ALL_BITS -1
#define MP_1_BIT 1
#define MP_2_BITS 2
#define MP_4_BITS 4
#define MP_8_BITS 8
#define MP_16_BITS 16
#define MP_32_BITS 32

//! \brief Structure which represents a mapper from keys to values. A value is an integer with
//! numbits. It internally uses the glib hashtable.
typedef struct mp_mapper {
	GHashTable *map;
	int largest_key;
	int numbits;
	int numvals_per_pos;
	int mask;
} mp_mapper;

//! Creates a mapper from integer keys to integer values.
mp_mapper * mp_create();

//! Sets how many bits each value has.
void mp_set_numbits(mp_mapper *m, int numbits);

//! Resets the mapper, i.e., erases all information.
void mp_reset(mp_mapper *m);

//! Destroys the mapper.
void mp_destroy(mp_mapper *m);

//! Assigns each consecutive value, starting at firstid, to the id of each cell returned by the iterator it.
mp_value_t mp_assign_from_celliterator(mp_mapper *m, higcit_celliterator *it, mp_value_t firstid);

//! Assigns each consecutive value, starting at firstid, to the id of each facet returned by the iterator fit.
mp_value_t mp_assign_from_facetiterator(mp_mapper *m, higfit_facetiterator *fit, mp_value_t firstid);

//! Assigns each consecutive value, starting at firstid, to the id of each leaf of the higtree rooted root.
mp_value_t mp_assign_standard_for_hig(mp_mapper *m, hig_cell *root, mp_value_t firstid);

//! Assigns val to id.
void mp_assign(mp_mapper *m, uniqueid id, mp_value_t val);

//! Gets the largest key assigned in the mapper.
int mp_get_largest_key(mp_mapper *m);

//! Gets the value assigned to id.
mp_value_t mp_lookup(mp_mapper *m, uniqueid id);

//! Defines a mapper from integer keys to pointers.
typedef mp_mapper mpp_mapper;

//! Creates a mapper to pointers.
mpp_mapper * mpp_create();

//! Destroys a mapper.
void mpp_destroy(mpp_mapper *m);

//! Assigns the pointer v to id.
void mpp_assign(mpp_mapper *m, uniqueid id, mpp_value_t v);

//! Gets the pointer assigned to id.
mpp_value_t mpp_lookup(mpp_mapper *m, uniqueid id);

//! Assigns to the mapper mt the value assigned to mapper mf, in the cell ids corresponding to the respective hig-trees.
void mp_copy(hig_cell *from, hig_cell *to, mp_mapper *mf, mp_mapper *mt);


#endif
