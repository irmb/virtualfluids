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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
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
//! \file RheologyPowellEyringModelLBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RheologyPowellEyringModelLBMKernel_H
#define RheologyPowellEyringModelLBMKernel_H

#include "RheologyK17LBMKernel.h"
#include "Rheology.h"

//! \brief    Cumulant LBM kernel + Herschel-Bulkley plastic model 
//! \author K. Kutscher, M. Geier
class RheologyPowellEyringModelLBMKernel : public RheologyK17LBMKernel
{
public:
	RheologyPowellEyringModelLBMKernel() {};
	~RheologyPowellEyringModelLBMKernel() {};
	SPtr<LBMKernel> clone() override
	{
		SPtr<LBMKernel> kernel(new RheologyPowellEyringModelLBMKernel());
		kernel->setNX(nx);
		kernel->setCollisionFactor(collFactor);
		dynamicPointerCast<RheologyPowellEyringModelLBMKernel>(kernel)->initDataSet();
		kernel->setBCSet(bcSet->clone(kernel));
		kernel->setWithForcing(withForcing);
		kernel->setForcingX1(muForcingX1);
		kernel->setForcingX2(muForcingX2);
		kernel->setForcingX3(muForcingX3);
		kernel->setIndex(ix1, ix2, ix3);
		kernel->setDeltaT(deltaT);

		return kernel;
	}
protected:
	real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const override
	{
		return Rheology::getPowellEyringCollFactor(omegaInf, shearRate, drho);
	}
};


#endif