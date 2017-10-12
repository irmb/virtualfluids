#ifndef global_h
#define global_h

typedef float doubflo;

#define DEBUG 1
#define GLOB_NODE 5

#define DIMENSION 3

#ifdef __unix__
#define PATH_TO_DATA "/home/soeren/gridgen/DATA/"
#else
#define PATH_TO_DATA "C:/Users/Soeren/Documents/Development/LBM_GRIG/DATA/"
#endif

#define TESTSUITE "TESTSUITE/"
#define VTK_OUTPUT "VTK_OUTPUT/"
#define STL_OUTPUT "STL_OUTPUT/"
#define STL "STL/"

#define MASTERRANK 0

#include <sstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


#include <GridGenerator/utilities/cuda/cudaDefines.h>

#endif // !global_h
