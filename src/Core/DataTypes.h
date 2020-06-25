#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>

#include "VirtualFluidsDefinitions.h"

#ifdef VF_DOUBLE_ACCURACY
typedef double real;
#else
typedef float  real;
#endif

typedef unsigned int uint;
#define INVALID_INDEX 4294967295 //max uint

#endif
