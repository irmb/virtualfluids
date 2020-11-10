#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>

#ifdef VF_DOUBLE_ACCURACY
using real = double;
#else
using real = float;
#endif

using uint = unsigned int;
#define INVALID_INDEX 4294967295 // max uint

#endif
