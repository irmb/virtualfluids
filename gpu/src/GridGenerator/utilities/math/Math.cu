//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Math.cu
//! \ingroup utilities
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "Math.h"

#include <cmath>

bool vf::Math::equal(const real& val1, const real& val2, real maxRelDiff)
{
	const real diff = std::fabs(val1 - val2);
	const real val1_abs = std::fabs(val1);
	const real val2_abs = std::fabs(val2);

	const real largest = (val2_abs > val1_abs) ? val2_abs : val1_abs;
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}

bool vf::Math::lessEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 < val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

bool vf::Math::greaterEqual(const real& val1, const real& val2, real maxRelDiff)
{
	if (val1 > val2 || equal(val1, val2, maxRelDiff))
		return true;
	return false;
}

real vf::Math::sqrtReal(const real& val)
{
#ifdef VF_DOUBLE_ACCURACY
    return sqrt(val);
#else
    return sqrtf(val);
#endif
}

real vf::Math::acosReal(const real& val)
{
#ifdef VF_DOUBLE_ACCURACY
    return acos(val);
#else
    return acosf(val);
#endif
}


