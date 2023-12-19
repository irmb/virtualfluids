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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#include "SlipBounceBack.h"
#include "DistributionArray3D.h"
#include "BoundaryConditions.h"

SlipBounceBack::SlipBounceBack()
{
   BCStrategy::preCollision = false;
}
//////////////////////////////////////////////////////////////////////////
SlipBounceBack::~SlipBounceBack()
{

}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> SlipBounceBack::clone()
{
   SPtr<BCStrategy> bc(new SlipBounceBack());
   return bc;
}
//////////////////////////////////////////////////////////////////////////
void SlipBounceBack::addDistributions(SPtr<DistributionArray3D> distributions)
{
   this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void SlipBounceBack::applyBC()
{
    using namespace vf::lbm::dir;

   real f[D3Q27System::ENDF+1];
   real feq[D3Q27System::ENDF+1];
   distributions->getPostCollisionDistribution(f, x1, x2, x3);
   real vx1, vx2, vx3, drho, rho;
   calcMacrosFct(f, drho, vx1, vx2, vx3);
   calcFeqFct(feq, drho, vx1, vx2, vx3);

   rho = vf::basics::constant::c1o1 + drho * compressibleFactor;

   UbTupleFloat3 normale = bcPtr->getNormalVector();
   real amp = vx1*val<1>(normale)+vx2*val<2>(normale)+vx3*val<3>(normale);

   vx1 = vx1 - amp * val<1>(normale); //normale zeigt von struktur weg!
   vx2 = vx2 - amp * val<2>(normale); //normale zeigt von struktur weg!
   vx3 = vx3 - amp * val<3>(normale); //normale zeigt von struktur weg!

   for (int fdir = D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
   {
      if (bcPtr->hasSlipBoundaryFlag(fdir))
      {
         //quadratic bounce back
         const int invDir = D3Q27System::INVDIR[fdir];
         real velocity = vf::basics::constant::c0o1;
         switch (invDir)
         {
         case dP00: velocity = (vf::basics::constant::c4o9*(+vx1)); break;      //(2/cs^2)(=6)*rho_0(=1 bei imkompr)*wi*u*ei mit cs=1/sqrt(3)
         case dM00: velocity = (vf::basics::constant::c4o9*(-vx1)); break;      //z.B. aus paper manfred MRT LB models in three dimensions (2002)   
         case d0P0: velocity = (vf::basics::constant::c4o9*(+vx2)); break;
         case d0M0: velocity = (vf::basics::constant::c4o9*(-vx2)); break;
         case d00P: velocity = (vf::basics::constant::c4o9*(+vx3)); break;
         case d00M: velocity = (vf::basics::constant::c4o9*(-vx3)); break;
         case dPP0: velocity = (vf::basics::constant::c1o9*(+vx1+vx2)); break;
         case dMM0: velocity = (vf::basics::constant::c1o9*(-vx1-vx2)); break;
         case dPM0: velocity = (vf::basics::constant::c1o9*(+vx1-vx2)); break;
         case dMP0: velocity = (vf::basics::constant::c1o9*(-vx1+vx2)); break;
         case dP0P: velocity = (vf::basics::constant::c1o9*(+vx1+vx3)); break;
         case dM0M: velocity = (vf::basics::constant::c1o9*(-vx1-vx3)); break;
         case dP0M: velocity = (vf::basics::constant::c1o9*(+vx1-vx3)); break;
         case dM0P: velocity = (vf::basics::constant::c1o9*(-vx1+vx3)); break;
         case d0PP: velocity = (vf::basics::constant::c1o9*(+vx2+vx3)); break;
         case d0MM: velocity = (vf::basics::constant::c1o9*(-vx2-vx3)); break;
         case d0PM: velocity = (vf::basics::constant::c1o9*(+vx2-vx3)); break;
         case d0MP: velocity = (vf::basics::constant::c1o9*(-vx2+vx3)); break;
         case dPPP: velocity = (vf::basics::constant::c1o36*(+vx1+vx2+vx3)); break;
         case dMMM: velocity = (vf::basics::constant::c1o36*(-vx1-vx2-vx3)); break;
         case dPPM: velocity = (vf::basics::constant::c1o36*(+vx1+vx2-vx3)); break;
         case dMMP: velocity = (vf::basics::constant::c1o36*(-vx1-vx2+vx3)); break;
         case dPMP: velocity = (vf::basics::constant::c1o36*(+vx1-vx2+vx3)); break;
         case dMPM: velocity = (vf::basics::constant::c1o36*(-vx1+vx2-vx3)); break;
         case dPMM: velocity = (vf::basics::constant::c1o36*(+vx1-vx2-vx3)); break;
         case dMPP: velocity = (vf::basics::constant::c1o36*(-vx1+vx2+vx3)); break;
         default: throw UbException(UB_EXARGS, "unknown error");
         }
         real fReturn = f[invDir] - velocity * rho;
         distributions->setPostCollisionDistributionForDirection(fReturn, x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);
      }
   }
}
//! \}
