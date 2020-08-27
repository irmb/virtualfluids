#include "Math.h"

#include <cmath>

#include "VirtualFluidsDefinitions.h"

HOSTDEVICE bool vf::Math::equal(const real& val1, const real& val2, real maxRelDiff)
{
	const real diff = std::fabs(val1 - val2);
	const real val1_abs = std::fabs(val1);
	const real val2_abs = std::fabs(val2);

	const real largest = (val2_abs > val1_abs) ? val2_abs : val1_abs;
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}

HOSTDEVICE bool vf::Math::lessEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 < val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE bool vf::Math::greaterEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 > val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

HOSTDEVICE real vf::Math::sqrtReal(const real& val)
{
#ifdef VF_DOUBLE_ACCURACY
    return sqrt(val);
#else
    return sqrtf(val);
#endif
}

HOSTDEVICE real vf::Math::acosReal(const real& val)
{
#ifdef VF_DOUBLE_ACCURACY
    return acos(val);
#else
    return acosf(val);
#endif
}


