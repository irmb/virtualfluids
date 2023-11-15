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
//! \file ThixotropyVelocityBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "ThixotropyVelocityBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

ThixotropyVelocityBCStrategy::ThixotropyVelocityBCStrategy()
{
	BCStrategy::type = BCStrategy::ThixotropyVelocityBCStrategy;
	BCStrategy::preCollision = false;
	BCStrategy::thixotropy = true;
	lambdaBC = vf::basics::constant::c0o1;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyVelocityBCStrategy::~ThixotropyVelocityBCStrategy()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> ThixotropyVelocityBCStrategy::clone()
{
	SPtr<BCStrategy> bc(new ThixotropyVelocityBCStrategy());
	dynamicPointerCast<ThixotropyVelocityBCStrategy>(bc)->setLambdaBC(lambdaBC);
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void ThixotropyVelocityBCStrategy::addDistributionsF(DistributionArray3DPtr distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyVelocityBCStrategy::applyBC()
{
	using namespace vf::lbm::dir;

	real f[D3Q27System::ENDF + 1];
	real feq[D3Q27System::ENDF + 1];
	real h[D3Q27System::ENDF + 1];

	distributions->getPostCollisionDistribution(f, x1, x2, x3);
	distributionsH->getPostCollisionDistribution(h, x1, x2, x3);
	
	real rho, vx1, vx2, vx3, drho;
	calcMacrosFct(f, drho, vx1, vx2, vx3);
	calcFeqFct(feq, drho, vx1, vx2, vx3);

	rho = vf::basics::constant::c1o1 + drho * compressibleFactor;

	//calcDiffusionMacrosFctPost(h, concentration, fl1, fl2, fl3, m100, collFactor);
	real lambda = D3Q27System::getDensity(h);

	int nx1 = x1;
	int nx2 = x2;
	int nx3 = x3;

	//flag points in direction of fluid
	if (bcPtr->hasVelocityBoundaryFlag(dP00)) { nx1 -= 1; }
	else if (bcPtr->hasVelocityBoundaryFlag(dM00)) { nx1 += 1; }
	else if (bcPtr->hasVelocityBoundaryFlag(d0P0)) { nx2 -= 1; }
	else if (bcPtr->hasVelocityBoundaryFlag(d0M0)) { nx2 += 1; }
	else if (bcPtr->hasVelocityBoundaryFlag(d00P)) { nx3 -= 1; }
	else if (bcPtr->hasVelocityBoundaryFlag(d00M)) { nx3 += 1; }
	else	 UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));

	//lambdaBC = bcPtr->getBoundaryThixotropy();

	//LBMReal rhoBC = bcPtr->getBoundaryDensity();

	//for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++)
	//{
	//	if (bcPtr->hasDensityBoundaryFlag(fdir))
	//	{
	//		LBMReal ftemp = calcFeqsForDirFct(fdir, rho, vx1, vx2, vx3);
	//		ftemp = calcFeqsForDirFct(fdir, rhoBC, vx1, vx2, vx3) + f[fdir] - ftemp;
	//		distributions->setPostCollisionDistributionForDirection(ftemp, nx1, nx2, nx3, fdir);

	//		LBMReal htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
	//		htemp = D3Q27System::getCompFeqForDirection(fdir,lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
	//		distributionsH->setPostCollisionDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
	//	}
	//}

	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
	{
		if (bcPtr->hasVelocityBoundaryFlag(fdir))
		{
			const int invDir = D3Q27System::INVDIR[fdir];
			real q = bcPtr->getQ(invDir);// m+m q=0 stabiler
			real velocity = bcPtr->getBoundaryVelocity(invDir);
			real fReturn = ((vf::basics::constant::c1o1 - q) / (vf::basics::constant::c1o1 + q)) * ((f[invDir] - feq[invDir]) / (vf::basics::constant::c1o1 - collFactor) + feq[invDir]) + ((q * (f[invDir] + f[fdir]) - velocity * rho) / (vf::basics::constant::c1o1 + q));
			distributions->setPostCollisionDistributionForDirection(fReturn, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

			real htemp = D3Q27System::getCompFeqForDirection(fdir, lambda, vx1, vx2, vx3);
			htemp = D3Q27System::getCompFeqForDirection(fdir, lambdaBC, vx1, vx2, vx3) + h[fdir] - htemp;
			distributionsH->setPostCollisionDistributionForDirection(htemp, nx1, nx2, nx3, fdir);
		}
	}
}
