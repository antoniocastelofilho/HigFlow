
#include <stdlib.h>
#include <string.h>
#include "Debug-c.h"
#include "mapper.h"

mp_mapper * mp_create() {
	mp_mapper *m = g_new(mp_mapper, 1);
	m->numbits = MP_ALL_BITS;
	m->numvals_per_pos = 1;
	m->map = g_hash_table_new(NULL, NULL);
	m->largest_key = 0;
	return m;
}

void mp_set_numbits(mp_mapper *m, int numbits) {
	assert(numbits <= MP_32_BITS);
	assert(
		numbits == MP_1_BIT ||
		numbits == MP_2_BITS ||
		numbits == MP_4_BITS ||
		numbits == MP_8_BITS ||
		numbits == MP_16_BITS ||
		numbits == MP_32_BITS
	);
	m->numbits = numbits;
	m->numvals_per_pos = sizeof(GPOINTER_TO_INT(NULL))*8/numbits;
	m->mask = 0;
	for(int i = 0, b = 1; i < m->numbits; i++, b <<= 1) {
		m->mask |= b;
	}
}

void mp_assign(mp_mapper *m, uniqueid id, mp_value_t val) {
	if (m->numbits == MP_ALL_BITS) {
		g_hash_table_insert(m->map, GINT_TO_POINTER(id), GINT_TO_POINTER(val+1));
	} else {
		DEBUG_PASS;
		assert(val < (1 << m->numbits));
		uniqueid new_id = id / m->numvals_per_pos;
		mp_value_t v = GPOINTER_TO_INT(g_hash_table_lookup(m->map, GINT_TO_POINTER(new_id)));
		int pos_id = id % m->numvals_per_pos;
		int shift = m->numbits * pos_id;
		mp_value_t new_v = (v & ~(m->mask << shift)) | (val << shift);
		g_hash_table_insert(m->map, GINT_TO_POINTER(new_id), GINT_TO_POINTER(new_v));
	}
	if (id > m->largest_key) {
		m->largest_key = id;
	}
}

mp_value_t mp_lookup(mp_mapper *m, uniqueid id) {
	if (m->numbits == MP_ALL_BITS) {
		gpointer *col = g_hash_table_lookup(m->map, GINT_TO_POINTER(id));
		return (mp_value_t) (GPOINTER_TO_INT(col) - 1);
	} else {
		uniqueid new_id = id / m->numvals_per_pos;
		mp_value_t v = GPOINTER_TO_INT(g_hash_table_lookup(m->map, GINT_TO_POINTER(new_id)));
		int pos_id = id % m->numvals_per_pos;
		int shift = m->numbits * pos_id;
		mp_value_t new_v = (v >> shift) % (1 << m->numbits);
		return new_v;
	}
}

int mp_get_largest_key(mp_mapper *m) {
	return m->largest_key;
}

mp_value_t mp_assign_from_celliterator(mp_mapper *m, higcit_celliterator *it, mp_value_t firstid) {
	mp_value_t id = firstid;
	const int redef = 0;
	for(; !higcit_isfinished(it); higcit_nextcell(it)) {
		hig_cell *c = higcit_getcell(it);
		uniqueid cid = hig_get_cid(c);
		if (redef || m->numbits != MP_ALL_BITS || mp_lookup(m, cid) == MP_UNDEF) {
			mp_assign(m, cid, id);
			id++;
		}
	}
	return id;
}

mp_value_t mp_assign_from_facetiterator(mp_mapper *m,
	higfit_facetiterator *fit, mp_value_t firstid)
{
	mp_value_t id = firstid;
	const int redef = 0;
	for(; !higfit_isfinished(fit); higfit_nextfacet(fit)) {
		hig_facet *f = higfit_getfacet(fit);
		uniqueid fid = hig_get_fid(f);
		if (redef || m->numbits != MP_ALL_BITS || mp_lookup(m, fid) == MP_UNDEF) {
			mp_assign(m, fid, id);
			id++;
		}
	}
	return id;
}

mp_value_t mp_assign_standard_for_hig(mp_mapper *m, hig_cell *root, mp_value_t firstid) {
	higcit_celliterator *it = higcit_create_all_leaves(root);
	mp_value_t id = mp_assign_from_celliterator(m, it, firstid);
	higcit_destroy(it);
	return id;
}

void mp_destroy(mp_mapper *m) {
	g_hash_table_destroy(m->map);
	free(m);
}

void mp_reset(mp_mapper *m) {
	g_hash_table_remove_all(m->map);
}

mpp_mapper * mpp_create() {
	return mp_create();
}

void mpp_destroy(mpp_mapper *m) {
	mp_destroy(m);
}

void mpp_assign(mpp_mapper *m, uniqueid id, mpp_value_t v) {
	g_hash_table_insert(m->map, GINT_TO_POINTER(id), v);
}

mpp_value_t mpp_lookup(mpp_mapper *m, uniqueid id) {
	return g_hash_table_lookup(m->map, GINT_TO_POINTER(id));
}

void mp_copy(hig_cell *from, hig_cell *to, mp_mapper *mf, mp_mapper *mt) {
	higcit_celliterator *itf, *itt;
	for(itf = higcit_create_all_leaves(from), itt = higcit_create_all_leaves(to);
	    !higcit_isfinished(itf) && !higcit_isfinished(itt);
	    higcit_nextcell(itf), higcit_nextcell(itt)) {
		hig_cell *cf = higcit_getcell(itf);
		hig_cell *ct = higcit_getcell(itt);
		int cfid = hig_get_cid(cf);
		int ctid = hig_get_cid(ct);
		int gfid = mp_lookup(mf, cfid);
		mp_assign(mt, ctid, gfid);
	}
	higcit_destroy(itf);
	higcit_destroy(itt);
}

