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
//! \file VelocityNonReflecting.cpp
//! \ingroup BoundarConditions
//! \author Hussein Alihussein
//=======================================================================================
#include "VelocityNonReflecting.h"

#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"

VelocityNonReflecting::VelocityNonReflecting()
{
    BCStrategy::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
VelocityNonReflecting::VelocityNonReflecting(real relaxationRate)
{
    BCStrategy::preCollision = true;
    this->BCVeloWeight = relaxationRate;
}
//////////////////////////////////////////////////////////////////////////
VelocityNonReflecting::~VelocityNonReflecting() = default;
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> VelocityNonReflecting::clone()
{
    SPtr<BCStrategy> bc(new VelocityNonReflecting(this->BCVeloWeight));
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void VelocityNonReflecting::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void VelocityNonReflecting::applyBC()
{
    using namespace vf::lbm::dir;
    using namespace D3Q27System;
 //   using namespace UbMath;
    using namespace vf::basics::constant;

    real f[ENDF + 1];
    real ftemp[ENDF + 1];

    int nx1       = x1;
    int nx2       = x2;
    int nx3       = x3;
    int direction = -1;

    // flag points in direction of fluid
    if (bcPtr->hasVelocityBoundaryFlag(dP00)) {
        nx1 += 1;
        direction = dP00;
        this->velocity = bcPtr->getBoundaryVelocityX1();
    } else if (bcPtr->hasVelocityBoundaryFlag(dM00)) {
        nx1 -= 1;
        direction = dM00;
        this->velocity = bcPtr->getBoundaryVelocityX1();
    } else if (bcPtr->hasVelocityBoundaryFlag(d0P0)) {
        nx2 += 1;
        direction = d0P0;
        this->velocity = bcPtr->getBoundaryVelocityX2();
    } else if (bcPtr->hasVelocityBoundaryFlag(d0M0)) {
        nx2 -= 1;
        direction = d0M0;
        this->velocity = bcPtr->getBoundaryVelocityX2();
    } else if (bcPtr->hasVelocityBoundaryFlag(d00P)) {
        nx3 += 1;
        direction = d00P;
        this->velocity = bcPtr->getBoundaryVelocityX3();
    } else if (bcPtr->hasVelocityBoundaryFlag(d00M)) {
        nx3 -= 1;
        direction = d00M;
        this->velocity = bcPtr->getBoundaryVelocityX3();
    } else
        UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on velocity boundary..."));

    distributions->getPreCollisionDistribution(f, x1, x2, x3);
    distributions->getPreCollisionDistribution(ftemp, nx1, nx2, nx3);

    real rho, vx1, vx2, vx3;
    calcMacrosFct(f, rho, vx1, vx2, vx3);
    //vx1                  = 0.;
    real BCVeloWeight = c1o2;
    // real velocity     = 0.004814077025232405; 
     // real velocity     = 0.00057735;
    //real velocity = 0.04; 
       real velocity = 0.01; 
     // real velocity = 1./112.; 
    // real velocity = 1./126.; 
     //real velocity = c1o100/2;
     // real velocity = 0.005; 
    //real delf         =(-velocity+vx1)*0.5 ;
    real delf; 

    switch (direction) {
        case dP00:
            delf = (-velocity + vx1) * BCVeloWeight; 
            // delf = (-velocity ) * BCVeloWeight;
            f[dP00]   = ftemp[dP00] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP00] - delf* WEIGTH[dP00];
            f[dPP0]  = ftemp[dPP0] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPP0]- delf* WEIGTH[dPP0];
            f[dPM0]  = ftemp[dPM0] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPM0]- delf* WEIGTH[dPM0];
            f[dP0P]  = ftemp[dP0P] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP0P]- delf* WEIGTH[dP0P];
            f[dP0M]  = ftemp[dP0M] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dP0M]- delf* WEIGTH[dP0M];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPPP]- delf* WEIGTH[dPPP];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPMP]- delf* WEIGTH[dPMP];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPPM]- delf* WEIGTH[dPPM];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 + vx1) + (c1o1 - c1oSqrt3 - vx1) * f[dPMM]- delf* WEIGTH[dPMM];
            //f[dP00] = (ftemp[dP00] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dP00]) *
            //           (1 - BCVeloWeight) +
            //       (ftemp[dM00] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dM00] +
            //       velocity*(6)*WEIGTH[dP00]/* bcPtr->getBoundaryVelocity(INVDIR[dM00])*/) *
            //           (BCVeloWeight)  ;
            //f[dPP0] = (ftemp[dPP0] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPP0]) *
            //            (1 - BCVeloWeight) +
            //        (ftemp[dMM0] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMM0] +
            //         velocity * (6) * WEIGTH[dPP0] /*bcPtr->getBoundaryVelocity(INVDIR[dMM0])*/) *
            //            (BCVeloWeight); 
            //f[dPM0] = (ftemp[dPM0] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPM0]) *
            //            (1 - BCVeloWeight) +
            //        (ftemp[dMP0] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMP0] +
            //        velocity*(6)*WEIGTH[dPP0]/* bcPtr->getBoundaryVelocity(INVDIR[dMP0])*/) *
            //            (BCVeloWeight); 
            //f[dP0P] = (ftemp[dP0P] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dP0P]) *
            //            (1 - BCVeloWeight) +
            //        (ftemp[dM0M] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dM0M] +
            //        velocity*(6)*WEIGTH[dP0P]/* bcPtr->getBoundaryVelocity(INVDIR[dM0M])*/) *
            //            (BCVeloWeight); 
            //f[dP0M] = (ftemp[dP0M] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dP0M])*
            //            (1 - BCVeloWeight) +
            //        (ftemp[dM0P] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dM0P] +
            //        velocity*(6)*WEIGTH[dP0M]/* bcPtr->getBoundaryVelocity(INVDIR[dM0P])*/) *
            //            (BCVeloWeight); 
            //f[dPPP] = (ftemp[dPPP] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPPP])*
            //            (1 - BCVeloWeight) +
            //        (ftemp[dMMM] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMMM] +
            //     velocity * (6) * WEIGTH[dPPP] /* bcPtr->getBoundaryVelocity(INVDIR[dMMM])*/) *
            //            (BCVeloWeight); 
            //f[dPMP] = (ftemp[dPMP] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPMP]) *
            //             (1 - BCVeloWeight) +
            //         (ftemp[dMPM] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMPM] +
            //     velocity * (6) * WEIGTH[dPPP] /*bcPtr->getBoundaryVelocity(INVDIR[dMPM])*/) *
            //             (BCVeloWeight); 
            //f[dPPM] = (ftemp[dPPM] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPPM]) *
            //             (1 - BCVeloWeight) +
            //         (ftemp[dMMP] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMMP] +
            //     velocity * (6) * WEIGTH[dPPP] /* bcPtr->getBoundaryVelocity(INVDIR[dMMP])*/) *
            //             (BCVeloWeight); 
            //f[dPMM] = (ftemp[dPMM] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dPMM]) *
            //             (1 - BCVeloWeight) +
            //         (ftemp[dMPP] * (c1oSqrt3 + vx1) + (1.0 - c1oSqrt3 - vx1) * f[dMPP] +
            //     velocity * (6) * WEIGTH[dPPP] /* bcPtr->getBoundaryVelocity(INVDIR[dMPP])*/) *
            //             (BCVeloWeight); 

            distributions->setPreCollisionDistributionForDirection(f[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            break;
        case dM00:
            delf = (-velocity - vx1) * BCVeloWeight;
            f[dM00] = ftemp[dM00] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM00] -
                   delf * WEIGTH[dM00];
            f[dMP0] = ftemp[dMP0] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMP0] -
                    delf * WEIGTH[dMP0];
            f[dMM0] = ftemp[dMM0] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMM0] -
                    delf * WEIGTH[dMM0];
            f[dM0P] = ftemp[dM0P] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM0P] -
                    delf * WEIGTH[dM0P];
            f[dM0M] = ftemp[dM0M] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dM0M] -
                    delf * WEIGTH[dM0M];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMPP] -
                     delf * WEIGTH[dMPP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMMP] -
                     delf * WEIGTH[dMMP];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMPM] -
                     delf * WEIGTH[dMPM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx1) + (c1o1 - c1oSqrt3 + vx1) * f[dMMM] -
                     delf * WEIGTH[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        case d0P0:
            delf = (-velocity + vx2) * BCVeloWeight;
            f[d0P0] = ftemp[d0P0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0P0] -
                   delf * WEIGTH[d0P0];
            f[dPP0] = ftemp[dPP0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPP0] -
                    delf * WEIGTH[dPP0];
            f[dMP0] = ftemp[dMP0] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMP0] -
                    delf * WEIGTH[dMP0];
            f[d0PP] = ftemp[d0PP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0PP] -
                    delf * WEIGTH[d0PP];
            f[d0PM] = ftemp[d0PM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[d0PM] -
                    delf * WEIGTH[d0PM];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPPP] -
                     delf * WEIGTH[dPPP];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMPP] -
                     delf * WEIGTH[dMPP];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dPPM] -
                     delf * WEIGTH[dPPM];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 + vx2) + (c1o1 - c1oSqrt3 - vx2) * f[dMPM] -
                     delf * WEIGTH[dMPM];

            distributions->setPreCollisionDistributionForDirection(f[d0P0], x1 + DX1[d0M0], x2 + DX2[d0M0], x3 + DX3[d0M0], d0M0);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            break;
        case d0M0:
            delf = (-velocity - vx2) * BCVeloWeight;
            f[d0M0] = ftemp[d0M0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0M0] -
                   delf * WEIGTH[d0M0];
            f[dPM0] = ftemp[dPM0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPM0] -
                    delf * WEIGTH[dPM0];
            f[dMM0] = ftemp[dMM0] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMM0] -
                    delf * WEIGTH[dMM0];
            f[d0MP] = ftemp[d0MP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0MP] -
                    delf * WEIGTH[d0MP];
            f[d0MM] = ftemp[d0MM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[d0MM] -
                    delf * WEIGTH[d0MM];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPMP] -
                     delf * WEIGTH[dPMP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMMP] -
                     delf * WEIGTH[dMMP];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dPMM] -
                     delf * WEIGTH[dPMM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx2) + (c1o1 - c1oSqrt3 + vx2) * f[dMMM] -
                     delf * WEIGTH[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d0M0], x1 + DX1[d0P0], x2 + DX2[d0P0], x3 + DX3[d0P0], d0P0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        case d00P:
            delf = (-velocity + vx3) * BCVeloWeight;
            f[d00P] = ftemp[d00P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d00P] -
                   delf * WEIGTH[d00P];
            f[dP0P] = ftemp[dP0P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dP0P] -
                    delf * WEIGTH[dP0P];
            f[dM0P] = ftemp[dM0P] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dM0P] -
                    delf * WEIGTH[dM0P];
            f[d0PP] = ftemp[d0PP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d0PP] -
                    delf * WEIGTH[d0PP];
            f[d0MP] = ftemp[d0MP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[d0MP] -
                    delf * WEIGTH[d0MP];
            f[dPPP] = ftemp[dPPP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dPPP] -
                     delf * WEIGTH[dPPP];
            f[dMPP] = ftemp[dMPP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dMPP] -
                     delf * WEIGTH[dMPP];
            f[dPMP] = ftemp[dPMP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dPMP] -
                     delf * WEIGTH[dPMP];
            f[dMMP] = ftemp[dMMP] * (c1oSqrt3 + vx3) + (c1o1 - c1oSqrt3 - vx3) * f[dMMP] -
                     delf * WEIGTH[dMMP];

            distributions->setPreCollisionDistributionForDirection(f[d00P], x1 + DX1[d00M], x2 + DX2[d00M], x3 + DX3[d00M], d00M);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            break;
        case d00M:
            delf = (-velocity - vx3) * BCVeloWeight;
            f[d00M] = ftemp[d00M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d00M] -
                   delf * WEIGTH[d00M];
            f[dP0M] = ftemp[dP0M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dP0M] -
                    delf * WEIGTH[dP0M];
            f[dM0M] = ftemp[dM0M] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dM0M] -
                    delf * WEIGTH[dM0M];
            f[d0PM] = ftemp[d0PM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d0PM] -
                    delf * WEIGTH[d0PM];
            f[d0MM] = ftemp[d0MM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[d0MM] -
                    delf * WEIGTH[d0MM];
            f[dPPM] = ftemp[dPPM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dPPM] -
                     delf * WEIGTH[dPPM];
            f[dMPM] = ftemp[dMPM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dMPM] -
                     delf * WEIGTH[dMPM];
            f[dPMM] = ftemp[dPMM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dPMM] -
                     delf * WEIGTH[dPMM];
            f[dMMM] = ftemp[dMMM] * (c1oSqrt3 - vx3) + (c1o1 - c1oSqrt3 + vx3) * f[dMMM] -
                     delf * WEIGTH[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d00M], x1 + DX1[d00P], x2 + DX2[d00P], x3 + DX3[d00P], d00P);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        default:
            UB_THROW(
                UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
    }
}
