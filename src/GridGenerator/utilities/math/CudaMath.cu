#include "CudaMath.cuh"

#include <cmath>

HOSTDEVICE bool CudaMath::equal(const real& val1, const real& val2, real maxRelDiff)
{
	real diff = std::fabs(val1 - val2);
	real val1_abs = std::fabs(val1);
	real val2_abs = std::fabs(val2);

	real largest = (val2_abs > val1_abs) ? val2_abs : val1_abs;
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}

HOSTDEVICE bool CudaMath::lessEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 < val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE bool CudaMath::greaterEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 > val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE real CudaMath::sqrt(const real& val)
{
    return sqrtf(val);
}

HOSTDEVICE real CudaMath::acos(const real& val)
{
    return acosf(val);
}

