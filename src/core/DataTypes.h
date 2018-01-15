/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DATATYPES_H
#define DATATYPES_H

#include "VirtualFluidsDefinitions.h"

#ifdef VF_DOUBLE_ACCURACY
typedef double real;
#else
typedef float  real;
#endif

typedef unsigned int uint;

#endif
