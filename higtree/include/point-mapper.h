#ifndef POINT_MAPPER
#define POINT_MAPPER

#include <stdbool.h>
#include <stdint.h>
#include <glib.h>

#include "coord.h"

typedef GHashTable point_mapper;

/*! Creates a best effort map from a Point to a pointer or int.
 *
 * The position is discretized by a delta of value EPSMACH. This means if two
 * points are farther apart by more than EPSMACH^(1/DIM) (box diagonal), they
 * will never match as equals. But the closer they are within the magnitude of
 * EPSMACH, greater the chance they are considered as equals.
 *
 * If value_free is not NULL, it will be used to free the value when an entry
 * is removed from the map. If NULL, no dealocation will be performed.
 */
point_mapper *ptm_create(void (*value_free)(void *));

/*! Find a value with key x.
 *
 * @return True if found, in which case, value is filled accordingly.
 */
bool ptm_lookup(point_mapper *ptm, CPPoint x, void **value);

/*! Removes the entry with key x from the map.
 *
 * If value_free was set on ptm_create(), value will be freed.
 *
 * @return True if found and removed, false otherwise.
 */
bool ptm_remove(point_mapper *ptm, CPPoint x);

/*! Inserts a value with key x.
 *
 * If key is already present, value will be replaced. Old value will be freed if
 * value_free was set on ptm_create().
 */
void ptm_insert(point_mapper *ptm, CPPoint x, void *value);

/*! Destroys the map.
 *
 * If value_free was set on ptm_create(), values will be freed.
 */
void ptm_destroy(point_mapper *ptm);

#endif
