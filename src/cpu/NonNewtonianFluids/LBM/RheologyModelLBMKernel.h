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
//! \file RheologyModelLBMKernel.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RheologyModelLBMKernel_H
#define RheologyModelLBMKernel_H

#include "LBMKernel.h"
#include "BCSet.h"
#include "D3Q27System.h"
#include "basics/utilities/UbTiming.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"

class RheologyModelLBMKernel;

//! \brief Base class for model of thixotropy based on K16. Use Template Method design pattern for Implementation of different models. 
//! \author K. Kutscher, M. Geier
class RheologyModelLBMKernel : public LBMKernel
{
public:
	RheologyModelLBMKernel();
	virtual ~RheologyModelLBMKernel();
	void calculate(int step);
	virtual SPtr<LBMKernel> clone() { UB_THROW(UbException("SPtr<LBMKernel> clone() - belongs in the derived class")); };
	real getCalculationTime();

	void swapDistributions();

protected:
	void initDataSet();

	virtual real getRheologyCollFactor(real omegaInf, real shearRate, real drho) const { UB_THROW(UbException("real getRheologyCollFactor() - belongs in the derived class")); }

	real f[D3Q27System::ENDF + 1];

	UbTimer timer;

	real OxyyMxzz;
	
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsF;
	CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsF;
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsF;

	mu::value_type muX1, muX2, muX3;
	mu::value_type muDeltaT;
	mu::value_type muNu;
	real forcingX1;
	real forcingX2;
	real forcingX3;

	bool test;
};

#endif
