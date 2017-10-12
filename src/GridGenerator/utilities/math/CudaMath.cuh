#ifndef CudaMath_H
#define CudaMath_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include "cuda.h"
#include <cuda_runtime.h>
#include "float.h"

#define EPSILON FLT_EPSILON

struct GridGenerator_EXPORT CudaMath
{
	HOSTDEVICE static bool equal(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff = EPSILON);
	HOSTDEVICE static bool lessEqual(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff = EPSILON);
	HOSTDEVICE static bool greaterEqual(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff = EPSILON);

	HOSTDEVICE static doubflo sqrt(const doubflo& val);
	HOSTDEVICE static doubflo acos(const doubflo& val);

};

#endif
