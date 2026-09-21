// The source of cell and facet identity.
//
// Ids come from a single process-local counter, so they are unique WITHIN a process
// and say nothing across ranks: two processes independently number their own cells
// starting from the same place.  The global numbering that the solver needs is built
// on top, by the mappers in pdomain.c, from the local ids plus the partition.

#ifndef UNIQUEID_H
#define UNIQUEID_H

#include<stdint.h>

//! Defines the type of a unique id.
typedef int32_t uniqueid;

//! Gets a unique id which time it is called.
uniqueid uid_getuniqueid();

//! \brief. Gets the current id. It is is the same as the id returned by the previous call to uid_getuniqueid.
uniqueid uid_getcurrentid();


#endif
