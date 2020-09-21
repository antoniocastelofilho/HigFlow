#include <string.h>
#include <stdint.h>

#include "higtree-iterator.h"
#include "higtree-serialize.h"

void hs_displaced_serialization_start(hig_serialization_ctx *hss,
	hig_cell *tree, CPPoint disp)
{
	hss->it = higcit_create_all_higtree(tree);
	POINT_ASSIGN(hss->displacement, disp);
}

void hs_serialization_start(hig_serialization_ctx *hss, hig_cell *tree)
{
	Point disp = {0};
	hs_displaced_serialization_start(hss, tree, disp);
}

size_t hs_serialize_n(hig_serialization_ctx *hss,
	hig_serial_tree *buf, size_t count)
{
	size_t i = 0;
	for(; i < count; ++i)
	{
		if(higcit_isfinished(hss->it))
			break;

		hig_cell *c = higcit_getcell(hss->it);

		POINT_ADD(buf[i].lowpoint, c->lowpoint, hss->displacement);
		POINT_ADD(buf[i].highpoint, c->highpoint, hss->displacement);
		POINT_ASSIGN(buf[i].numcells, c->numcells);

		higcit_nextcell(hss->it);
	}

	return i;
}

void hs_serialization_finish(hig_serialization_ctx *hss)
{
	higcit_destroy(hss->it);
	hss->it = NULL;
}

void hs_deserialization_start(hig_deserialization_ctx *hsd)
{
	hsd->tree = NULL;
}

size_t hs_deserialize_n(hig_deserialization_ctx *hsd,
	hig_serial_tree *buf, size_t count)
{
	if(count == 0)
		return 0;

	if(!hsd->tree) {
		hsd->tree = hig_create_root(buf[0].lowpoint, buf[0].highpoint);
		hsd->it = higcit_create_all_higtree(hsd->tree);
	}

	size_t i = 0;
	for(; i < count; ++i)
	{
		if(higcit_isfinished(hsd->it))
			break;

		hig_cell *cell = higcit_getcell(hsd->it);
		hig_fill_empty(cell, buf[i].lowpoint, buf[i].highpoint);
		hig_refine_empty(cell, buf[i].numcells);

		unsigned numcells = hig_get_number_of_children(cell);
		for(unsigned j = 0; j < numcells; ++j) {
			hig_create_empty_leaf(cell, j);
		}

		higcit_nextcell(hsd->it);
	}

	return i;
}

hig_cell* hs_deserialization_finish(hig_deserialization_ctx *hsd)
{
	higcit_destroy(hsd->it);
	return hsd->tree;
}

size_t hs_serialize(hig_serial_tree *buf, hig_cell *tree)
{
	Point disp = {0};
	return hs_serialize_displaced(buf, tree, disp);
}

size_t hs_serialize_displaced(hig_serial_tree *buf, hig_cell *tree, CPPoint disp)
{
	hig_serialization_ctx hss;
	hs_displaced_serialization_start(&hss, tree, disp);
	size_t ret = hs_serialize_n(&hss, buf, SIZE_MAX);
	hs_serialization_finish(&hss);

	return ret;
}

size_t hs_deserialize(hig_cell **tree, hig_serial_tree *buf)
{
	hig_deserialization_ctx hsd;
	hs_deserialization_start(&hsd);
	size_t ret = hs_deserialize_n(&hsd, buf, SIZE_MAX);
	*tree = hs_deserialization_finish(&hsd);

	return ret;
}
