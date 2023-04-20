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
//! \file ThixotropyNoSlipBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#include "ThixotropyNoSlipBCStrategy.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

ThixotropyNoSlipBCStrategy::ThixotropyNoSlipBCStrategy()
{
	BCStrategy::type = BCStrategy::ThixotropyNoSlipBCStrategy;
	BCStrategy::preCollision = false;
	BCStrategy::thixotropy = true;
	
}
//////////////////////////////////////////////////////////////////////////
ThixotropyNoSlipBCStrategy::~ThixotropyNoSlipBCStrategy()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> ThixotropyNoSlipBCStrategy::clone()
{
	SPtr<BCStrategy> bc(new ThixotropyNoSlipBCStrategy());
	return bc;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNoSlipBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
	this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
//void ThixotropyNoSlipBCStrategy::addDistributionsF(SPtr<DistributionArray3D> distributions)
//{
//	this->distributionsf = distributions;
//}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNoSlipBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributions)
{
	this->distributionsH = distributions;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyNoSlipBCStrategy::applyBC()
{
	real f[D3Q27System::ENDF + 1];
	real feq[D3Q27System::ENDF + 1];
	real h[D3Q27System::ENDF + 1];
	real heq[D3Q27System::ENDF + 1];
	distributions->getDistributionInv(f, x1, x2, x3);
	distributionsH->getDistributionInv(h, x1, x2, x3);
	real rho, vx1, vx2, vx3;//, concentration, fl1, fl2, fl3, m100;
	calcMacrosFct(f, rho, vx1, vx2, vx3);
	calcFeqFct(feq, rho, vx1, vx2, vx3);

	//calcDiffusionMacrosFctPost(h, concentration, fl1, fl2, fl3, m100, collFactor);
	real lambda = D3Q27System::getDensity(h);
	D3Q27System::calcCompFeq(heq, lambda, 0., 0., 0.);

	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++)
	{
		if (bcPtr->hasNoSlipBoundaryFlag(fdir))
		{
			//quadratic bounce back
			const int invDir = D3Q27System::INVDIR[fdir];
			real q = bcPtr->getQ(invDir);
			real fReturnf = ((vf::basics::constant::c1o1 - q) / (vf::basics::constant::c1o1 + q))*((f[invDir] - feq[invDir]) / (vf::basics::constant::c1o1 - collFactor) + feq[invDir]) + ((q / (vf::basics::constant::c1o1 + q))*(f[invDir] + f[fdir]));
			real fReturnh = ((vf::basics::constant::c1o1 - q) / (vf::basics::constant::c1o1 + q))*((h[invDir] - heq[invDir]) / (vf::basics::constant::c1o1 - collFactor) + heq[invDir]) + ((q / (vf::basics::constant::c1o1 + q))*(h[invDir] + h[fdir]));

			distributions->setDistributionForDirection(fReturnf, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);
			distributionsH->setDistributionForDirection(fReturnh, x1 + D3Q27System::DX1[invDir], x2 + D3Q27System::DX2[invDir], x3 + D3Q27System::DX3[invDir], fdir);

		}
	}
}
