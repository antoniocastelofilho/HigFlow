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
