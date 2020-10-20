#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>

#ifdef VF_DOUBLE_ACCURACY
typedef double real;
#else
using real = float;
#endif

using uint = unsigned int;
#define INVALID_INDEX 4294967295 // max uint

#endif
