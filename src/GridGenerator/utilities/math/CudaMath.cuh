#ifndef CudaMath_H
#define CudaMath_H

#include "GridGenerator/global.h"


#include "cuda.h"
#include <cuda_runtime.h>
#include "float.h"

#define EPSILON FLT_EPSILON

class VF_PUBLIC CudaMath
{
public:
	HOSTDEVICE static bool equal(const real& val1, const real& val2, real maxRelDiff = EPSILON);
	HOSTDEVICE static bool lessEqual(const real& val1, const real& val2, real maxRelDiff = EPSILON);
	HOSTDEVICE static bool greaterEqual(const real& val1, const real& val2, real maxRelDiff = EPSILON);

	HOSTDEVICE static real sqrtReal(const real& val);
	HOSTDEVICE static real acosReal(const real& val);
    
    HOSTDEVICE static real getDecimalPart(real number) {
        return number - real(int(number));
    }


};

#endif
