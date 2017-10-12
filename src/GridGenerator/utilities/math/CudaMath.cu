#include "CudaMath.cuh"

#include <cmath>

HOSTDEVICE bool CudaMath::equal(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff)
{
	doubflo diff = std::fabs(val1 - val2);
	doubflo val1_abs = std::fabs(val1);
	doubflo val2_abs = std::fabs(val2);

	doubflo largest = (val2_abs > val1_abs) ? val2_abs : val1_abs;
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}

HOSTDEVICE bool CudaMath::lessEqual(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff)
{
	if (val1 < val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE bool CudaMath::greaterEqual(const doubflo& val1, const doubflo& val2, doubflo maxRelDiff)
{
	if (val1 > val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE doubflo CudaMath::sqrt(const doubflo& val)
{
    return sqrtf(val);
}

HOSTDEVICE doubflo CudaMath::acos(const doubflo& val)
{
    return acosf(val);
}

