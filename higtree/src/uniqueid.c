
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
