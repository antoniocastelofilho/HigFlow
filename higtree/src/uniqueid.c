// The counter behind `uniqueid`, the identity a cell keeps regardless of which rank
// holds it or how the domain is divided.
//
// Contrast with the local index from mapper.c, which is per-process and per-partition.


#include <stdlib.h>
#include "higtree.h"
#include "uniqueid.h"

#include "Debug-c.h"

static uniqueid curid = 1;

uniqueid uid_getuniqueid() {
	uniqueid id = 0;
	//unsigned short int d = 4;
	//id = (curid << d) + (rand() % (1<<d));
	id = curid;
	curid++;
	return id;
}

uniqueid uid_getcurrentid() {
	return curid;
}
