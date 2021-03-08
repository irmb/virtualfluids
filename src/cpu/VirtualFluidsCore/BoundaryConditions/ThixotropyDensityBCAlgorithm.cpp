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
//! \file ThixotropyDensityBCAlgorithm.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#include "ThixotropyDensityBCAlgorithm.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

ThixotropyDensityBCAlgorithm::ThixotropyDensityBCAlgorithm()
{
	BCAlgorithm::type = BCAlgorithm::ThixotropyDensityBCAlgorithm;
	BCAlgorithm::preCollision = false;
	BCAlgorithm::thixotropy = true;
	lambdaBC = 0.0;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyDensityBCAlgorithm::~ThixotropyDensityBCAlgorithm()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCAlgorithm> ThixotropyDensityBCAlgorithm::clone()
{
	SPtr<BCAlgorithm> bc(new ThixotropyDensityBCAlgorithm());
	dynamicPointerCast<ThixotropyDensityBCAlgorithm>(bc)->setLambdaBC(lambdaBC);
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCAlgorithm::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void ThixotropyDensityBCAlgorithm::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCAlgorithm::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyDensityBCAlgorithm::applyBC()
{
   using namespace D3Q27System;

	LBMReal f[D3Q27System::ENDF + 1];
	LBMReal feq[D3Q27System::ENDF + 1];
	LBMReal h[D3Q27System::ENDF + 1];
	LBMReal heq[D3Q27System::ENDF + 1];
	distributions->getDistributionInv(f, x1, x2, x3);
	distributionsH->getDistributionInv(h, x1, x2, x3);
	
	LBMReal rho, vx1, vx2, vx3;
	
	calcMacrosFct(f, rho, vx1, vx2, vx3);
	calcFeqFct(feq, rho, vx1, vx2, vx3);

	LBMReal lambda = D3Q27System::getDensity(h);
	D3Q27System::calcCompFeq(heq, lambda, vx1, vx2, vx3);


	int nx1 = x1;
	int nx2 = x2;
	int nx3 = x3;

	//flag points in direction of fluid
	if (bcPtr->hasDensityBoundaryFlag(D3Q27System::E)) { nx1 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::W)) { nx1 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::N)) { nx2 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::S)) { nx2 += 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::T)) { nx3 -= 1; }
	else if (bcPtr->hasDensityBoundaryFlag(D3Q27System::B)) { nx3 += 1; }
	else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

	LBMReal rhoBC = bcPtr->getBoundaryDensity();

	for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++)
	{
		if (bcPtr->hasDensityBoundaryFlag(fdir))
		{
			LBMReal ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
			ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
			distributions->setDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);

			LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
			htemp = D3Q27System::getCompFeqForDirection(fdir,lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
			distributionsH->setDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
		}
	}
}
