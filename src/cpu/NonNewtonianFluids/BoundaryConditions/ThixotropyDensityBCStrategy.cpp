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
//! \file ThixotropyDensityBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#include "ThixotropyDensityBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

ThixotropyDensityBCStrategy::ThixotropyDensityBCStrategy()
{
	BCStrategy::type = BCStrategy::ThixotropyDensityBCStrategy;
	BCStrategy::preCollision = false;
	BCStrategy::thixotropy = true;
	lambdaBC = 0.0;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyDensityBCStrategy::~ThixotropyDensityBCStrategy()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> ThixotropyDensityBCStrategy::clone()
{
	SPtr<BCStrategy> bc(new ThixotropyDensityBCStrategy());
	dynamicPointerCast<ThixotropyDensityBCStrategy>(bc)->setLambdaBC(lambdaBC);
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void ThixotropyDensityBCStrategy::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCStrategy::applyBC()
{
	using namespace vf::lbm::dir;
    using namespace D3Q27System;

	real f[D3Q27System::ENDF + 1];
	real feq[D3Q27System::ENDF + 1];
	real h[D3Q27System::ENDF + 1];
	real heq[D3Q27System::ENDF + 1];
	distributions->getPostCollisionDistribution(f, x1, x2, x3);
	distributionsH->getPostCollisionDistribution(h, x1, x2, x3);
	
	real rho, vx1, vx2, vx3;
	
	calcMacrosFct(f, rho, vx1, vx2, vx3);
	calcFeqFct(feq, rho, vx1, vx2, vx3);

	real lambda = D3Q27System::getDensity(h);
	D3Q27System::calcCompFeq(heq, lambda, vx1, vx2, vx3);


	int nx1 = x1;
	int nx2 = x2;
	int nx3 = x3;

	//flag points in direction of fluid
	if (bcPtr->hasDensityBoundaryFlag(DIR_P00)) { nx1 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(DIR_M00)) { nx1 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(DIR_0P0)) { nx2 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(DIR_0M0)) { nx2 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(DIR_00P)) { nx3 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(DIR_00M)) { nx3 += 1; }
	else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

	real rhoBC = bcPtr->getBoundaryDensity();

	for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++)
	{
		if (bcPtr->hasDensityBoundaryFlag(fdir))
		{
			real ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
			ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
			distributions->setPostCollisionDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);

			real htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
			htemp = D3Q27System::getCompFeqForDirection(fdir,lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
			distributionsH->setPostCollisionDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
		}
	}
}
