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
//! \file MultiphaseNonReflectingOutflowBCStrategy.cpp
//! \ingroup BoundarConditions
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseNonReflectingOutflowBCStrategy.h"
#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "DistributionArray3D.h"

MultiphaseNonReflectingOutflowBCStrategy::MultiphaseNonReflectingOutflowBCStrategy()
{
    BCStrategy::type = BCStrategy::NonReflectingOutflowBCStrategy;
    BCStrategy::preCollision = true;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseNonReflectingOutflowBCStrategy::~MultiphaseNonReflectingOutflowBCStrategy()
{
}
//////////////////////////////////////////////////////////////////////////
SPtr<BCStrategy> MultiphaseNonReflectingOutflowBCStrategy::clone()
{
    SPtr<BCStrategy> bc(new MultiphaseNonReflectingOutflowBCStrategy());
    return bc;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCStrategy::addDistributions(SPtr<DistributionArray3D> distributions)
{
    this->distributions = distributions;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCStrategy::addDistributionsH(SPtr<DistributionArray3D> distributionsH)
{
    this->distributionsH = distributionsH;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCStrategy::addDistributionsH2(SPtr<DistributionArray3D> distributionsH2)
{
    this->distributionsH2 = distributionsH2;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseNonReflectingOutflowBCStrategy::applyBC()
{
    using namespace D3Q27System;
//    using namespace UbMath;
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    real f[ENDF + 1];
    real ftemp[ENDF + 1];
    real h[D3Q27System::ENDF + 1];
    real htemp[ENDF + 1];
    real h2[D3Q27System::ENDF + 1];
    real h2temp[ENDF + 1];

    int nx1 = x1;
    int nx2 = x2;
    int nx3 = x3;
    int direction = -1;

    // flag points in direction of fluid
    if (bcPtr->hasDensityBoundaryFlag(dP00)) {
        nx1 += 1;
        direction = dP00;
    } else if (bcPtr->hasDensityBoundaryFlag(dM00)) {
        nx1 -= 1;
        direction = dM00;
    } else if (bcPtr->hasDensityBoundaryFlag(d0P0)) {
        nx2 += 1;
        direction = d0P0;
    } else if (bcPtr->hasDensityBoundaryFlag(d0M0)) {
        nx2 -= 1;
        direction = d0M0;
    } else if (bcPtr->hasDensityBoundaryFlag(d00P)) {
        nx3 += 1;
        direction = d00P;
    } else if (bcPtr->hasDensityBoundaryFlag(d00M)) {
        nx3 -= 1;
        direction = d00M;
    } else
        UB_THROW(UbException(UB_EXARGS, "Danger...no orthogonal BC-Flag on density boundary..."));

    distributions->getPreCollisionDistribution(f, x1, x2, x3);
    distributions->getPreCollisionDistribution(ftemp, nx1, nx2, nx3);
    distributionsH->getPreCollisionDistribution(h, x1, x2, x3);
    distributionsH->getPreCollisionDistribution(htemp, nx1, nx2, nx3);
    distributionsH2->getPreCollisionDistribution(h2, x1, x2, x3);
    distributionsH2->getPreCollisionDistribution(h2temp, nx1, nx2, nx3);

    real /* phi,*/ p1, vx1, vx2, vx3;

    // D3Q27System::calcDensity(h, phi);

    calcMacrosFct(f, p1, vx1, vx2, vx3);

    switch (direction) {
        case dP00:
            f[dP00] = ftemp[dP00] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dP00];
            f[dPP0] = ftemp[dPP0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPP0];
            f[dPM0] = ftemp[dPM0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPM0];
            f[dP0P] = ftemp[dP0P] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dP0P];
            f[dP0M] = ftemp[dP0M] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dP0M];
            f[dPPP] = ftemp[dPPP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPPP];
            f[dPMP] = ftemp[dPMP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPMP];
            f[dPPM] = ftemp[dPPM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPPM];
            f[dPMM] = ftemp[dPMM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * f[dPMM];

            distributions->setPreCollisionDistributionForDirection(f[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);

            h[dP00] = htemp[dP00] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dP00];
            h[dPP0] = htemp[dPP0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPP0];
            h[dPM0] = htemp[dPM0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPM0];
            h[dP0P] = htemp[dP0P] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dP0P];
            h[dP0M] = htemp[dP0M] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dP0M];
            h[dPPP] = htemp[dPPP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPPP];
            h[dPMP] = htemp[dPMP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPMP];
            h[dPPM] = htemp[dPPM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPPM];
            h[dPMM] = htemp[dPMM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h[dPMM];

            distributionsH->setPreCollisionDistributionForDirection(h[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributionsH->setPreCollisionDistributionForDirection(h[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributionsH->setPreCollisionDistributionForDirection(h[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributionsH->setPreCollisionDistributionForDirection(h[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributionsH->setPreCollisionDistributionForDirection(h[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);

            h2[dP00] = c1o2 * (h2temp[dP00] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dP00]);
            h2[dPP0] = c1o2 * (h2temp[dPP0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPP0]);
            h2[dPM0] = c1o2 * (h2temp[dPM0] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPM0]);
            h2[dP0P] = c1o2 * (h2temp[dP0P] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dP0P]);
            h2[dP0M] = c1o2 * (h2temp[dP0M] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dP0M]);
            h2[dPPP] = c1o2 * (h2temp[dPPP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPPP]);
            h2[dPMP] = c1o2 * (h2temp[dPMP] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPMP]);
            h2[dPPM] = c1o2 * (h2temp[dPPM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPPM]);
            h2[dPMM] = c1o2 * (h2temp[dPMM] * (one_over_sqrt3 + vx1) + (c1o1 - one_over_sqrt3 - vx1) * h2[dPMM]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[dP00], x1 + DX1[dM00], x2 + DX2[dM00], x3 + DX3[dM00], dM00);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);

            break;
        case dM00:
            f[dM00] = ftemp[dM00] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dM00];
            f[dMP0] = ftemp[dMP0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMP0];
            f[dMM0] = ftemp[dMM0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMM0];
            f[dM0P] = ftemp[dM0P] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dM0P];
            f[dM0M] = ftemp[dM0M] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dM0M];
            f[dMPP] = ftemp[dMPP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMPP];
            f[dMMP] = ftemp[dMMP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMMP];
            f[dMPM] = ftemp[dMPM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMPM];
            f[dMMM] = ftemp[dMMM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h[dM00] = htemp[dM00] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dM00];
            h[dMP0] = htemp[dMP0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMP0];
            h[dMM0] = htemp[dMM0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMM0];
            h[dM0P] = htemp[dM0P] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dM0P];
            h[dM0M] = htemp[dM0M] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dM0M];
            h[dMPP] = htemp[dMPP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMPP];
            h[dMMP] = htemp[dMMP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMMP];
            h[dMPM] = htemp[dMPM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMPM];
            h[dMMM] = htemp[dMMM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h[dMMM];

            distributionsH->setPreCollisionDistributionForDirection(h[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributionsH->setPreCollisionDistributionForDirection(h[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributionsH->setPreCollisionDistributionForDirection(h[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributionsH->setPreCollisionDistributionForDirection(h[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributionsH->setPreCollisionDistributionForDirection(h[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h2[dM00] = c1o2 * (htemp[dM00] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dM00]);
            h2[dMP0] = c1o2 * (htemp[dMP0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMP0]);
            h2[dMM0] = c1o2 * (htemp[dMM0] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMM0]);
            h2[dM0P] = c1o2 * (htemp[dM0P] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dM0P]);
            h2[dM0M] = c1o2 * (htemp[dM0M] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dM0M]);
            h2[dMPP] = c1o2 * (htemp[dMPP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMPP]);
            h2[dMMP] = c1o2 * (htemp[dMMP] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMMP]);
            h2[dMPM] = c1o2 * (htemp[dMPM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMPM]);
            h2[dMMM] = c1o2 * (htemp[dMMM] * (one_over_sqrt3 - vx1) + (c1o1 - one_over_sqrt3 + vx1) * h2[dMMM]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[dM00], x1 + DX1[dP00], x2 + DX2[dP00], x3 + DX3[dP00], dP00);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);
            break;
        case d0P0:
            f[d0P0] = ftemp[d0P0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[d0P0];
            f[dPP0] = ftemp[dPP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dPP0];
            f[dMP0] = ftemp[dMP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dMP0];
            f[d0PP] = ftemp[d0PP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[d0PP];
            f[d0PM] = ftemp[d0PM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[d0PM];
            f[dPPP] = ftemp[dPPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dPPP];
            f[dMPP] = ftemp[dMPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dMPP];
            f[dPPM] = ftemp[dPPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dPPM];
            f[dMPM] = ftemp[dMPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * f[dMPM];

            distributions->setPreCollisionDistributionForDirection(f[d0P0], x1 + DX1[d0M0], x2 + DX2[d0M0], x3 + DX3[d0M0], d0M0);
            distributions->setPreCollisionDistributionForDirection(f[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributions->setPreCollisionDistributionForDirection(f[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);

            h[d0P0] = htemp[d0P0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[d0P0];
            h[dPP0] = htemp[dPP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dPP0];
            h[dMP0] = htemp[dMP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dMP0];
            h[d0PP] = htemp[d0PP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[d0PP];
            h[d0PM] = htemp[d0PM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[d0PM];
            h[dPPP] = htemp[dPPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dPPP];
            h[dMPP] = htemp[dMPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dMPP];
            h[dPPM] = htemp[dPPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dPPM];
            h[dMPM] = htemp[dMPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h[dMPM];

            distributionsH->setPreCollisionDistributionForDirection(h[d0P0], x1 + DX1[d0M0], x2 + DX2[d0M0], x3 + DX3[d0M0], d0M0);
            distributionsH->setPreCollisionDistributionForDirection(h[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributionsH->setPreCollisionDistributionForDirection(h[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributionsH->setPreCollisionDistributionForDirection(h[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributionsH->setPreCollisionDistributionForDirection(h[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);

            h2[d0P0] = c1o2 * (htemp[d0P0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[d0P0]);
            h2[dPP0] = c1o2 * (htemp[dPP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dPP0]);
            h2[dMP0] = c1o2 * (htemp[dMP0] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dMP0]);
            h2[d0PP] = c1o2 * (htemp[d0PP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[d0PP]);
            h2[d0PM] = c1o2 * (htemp[d0PM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[d0PM]);
            h2[dPPP] = c1o2 * (htemp[dPPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dPPP]);
            h2[dMPP] = c1o2 * (htemp[dMPP] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dMPP]);
            h2[dPPM] = c1o2 * (htemp[dPPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dPPM]);
            h2[dMPM] = c1o2 * (htemp[dMPM] * (one_over_sqrt3 + vx2) + (c1o1 - one_over_sqrt3 - vx2) * h2[dMPM]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[d0P0], x1 + DX1[d0M0], x2 + DX2[d0M0], x3 + DX3[d0M0], d0M0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPP0], x1 + DX1[dMM0], x2 + DX2[dMM0], x3 + DX3[dMM0], dMM0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMP0], x1 + DX1[dPM0], x2 + DX2[dPM0], x3 + DX3[dPM0], dPM0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);

            break;
        case d0M0:
            f[d0M0] = ftemp[d0M0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[d0M0];
            f[dPM0] = ftemp[dPM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dPM0];
            f[dMM0] = ftemp[dMM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dMM0];
            f[d0MP] = ftemp[d0MP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[d0MP];
            f[d0MM] = ftemp[d0MM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[d0MM];
            f[dPMP] = ftemp[dPMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dPMP];
            f[dMMP] = ftemp[dMMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dMMP];
            f[dPMM] = ftemp[dPMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dPMM];
            f[dMMM] = ftemp[dMMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d0M0], x1 + DX1[d0P0], x2 + DX2[d0P0], x3 + DX3[d0P0], d0P0);
            distributions->setPreCollisionDistributionForDirection(f[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributions->setPreCollisionDistributionForDirection(f[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h[d0M0] = htemp[d0M0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[d0M0];
            h[dPM0] = htemp[dPM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dPM0];
            h[dMM0] = htemp[dMM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dMM0];
            h[d0MP] = htemp[d0MP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[d0MP];
            h[d0MM] = htemp[d0MM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[d0MM];
            h[dPMP] = htemp[dPMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dPMP];
            h[dMMP] = htemp[dMMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dMMP];
            h[dPMM] = htemp[dPMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dPMM];
            h[dMMM] = htemp[dMMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h[dMMM];

            distributionsH->setPreCollisionDistributionForDirection(h[d0M0], x1 + DX1[d0P0], x2 + DX2[d0P0], x3 + DX3[d0P0], d0P0);
            distributionsH->setPreCollisionDistributionForDirection(h[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributionsH->setPreCollisionDistributionForDirection(h[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributionsH->setPreCollisionDistributionForDirection(h[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributionsH->setPreCollisionDistributionForDirection(h[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h2[d0M0] = c1o2 * (htemp[d0M0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[d0M0]);
            h2[dPM0] = c1o2 * (htemp[dPM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dPM0]);
            h2[dMM0] = c1o2 * (htemp[dMM0] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dMM0]);
            h2[d0MP] = c1o2 * (htemp[d0MP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[d0MP]);
            h2[d0MM] = c1o2 * (htemp[d0MM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[d0MM]);
            h2[dPMP] = c1o2 * (htemp[dPMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dPMP]);
            h2[dMMP] = c1o2 * (htemp[dMMP] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dMMP]);
            h2[dPMM] = c1o2 * (htemp[dPMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dPMM]);
            h2[dMMM] = c1o2 * (htemp[dMMM] * (one_over_sqrt3 - vx2) + (c1o1 - one_over_sqrt3 + vx2) * h2[dMMM]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[d0M0], x1 + DX1[d0P0], x2 + DX2[d0P0], x3 + DX3[d0P0], d0P0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPM0], x1 + DX1[dMP0], x2 + DX2[dMP0], x3 + DX3[dMP0], dMP0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMM0], x1 + DX1[dPP0], x2 + DX2[dPP0], x3 + DX3[dPP0], dPP0);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            break;
        case d00P:
            f[d00P] = ftemp[d00P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[d00P];
            f[dP0P] = ftemp[dP0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dP0P];
            f[dM0P] = ftemp[dM0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dM0P];
            f[d0PP] = ftemp[d0PP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[d0PP];
            f[d0MP] = ftemp[d0MP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[d0MP];
            f[dPPP] = ftemp[dPPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dPPP];
            f[dMPP] = ftemp[dMPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dMPP];
            f[dPMP] = ftemp[dPMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dPMP];
            f[dMMP] = ftemp[dMMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * f[dMMP];

            distributions->setPreCollisionDistributionForDirection(f[d00P], x1 + DX1[d00M], x2 + DX2[d00M], x3 + DX3[d00M], d00M);
            distributions->setPreCollisionDistributionForDirection(f[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributions->setPreCollisionDistributionForDirection(f[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributions->setPreCollisionDistributionForDirection(f[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributions->setPreCollisionDistributionForDirection(f[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributions->setPreCollisionDistributionForDirection(f[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributions->setPreCollisionDistributionForDirection(f[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributions->setPreCollisionDistributionForDirection(f[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributions->setPreCollisionDistributionForDirection(f[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);

            h[d00P] = htemp[d00P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[d00P];
            h[dP0P] = htemp[dP0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dP0P];
            h[dM0P] = htemp[dM0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dM0P];
            h[d0PP] = htemp[d0PP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[d0PP];
            h[d0MP] = htemp[d0MP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[d0MP];
            h[dPPP] = htemp[dPPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dPPP];
            h[dMPP] = htemp[dMPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dMPP];
            h[dPMP] = htemp[dPMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dPMP];
            h[dMMP] = htemp[dMMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h[dMMP];

            distributionsH->setPreCollisionDistributionForDirection(h[d00P], x1 + DX1[d00M], x2 + DX2[d00M], x3 + DX3[d00M], d00M);
            distributionsH->setPreCollisionDistributionForDirection(h[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributionsH->setPreCollisionDistributionForDirection(h[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributionsH->setPreCollisionDistributionForDirection(h[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributionsH->setPreCollisionDistributionForDirection(h[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);

            h2[d00P] = c1o2 * (htemp[d00P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[d00P]);
            h2[dP0P] = c1o2 * (htemp[dP0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dP0P]);
            h2[dM0P] = c1o2 * (htemp[dM0P] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dM0P]);
            h2[d0PP] = c1o2 * (htemp[d0PP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[d0PP]);
            h2[d0MP] = c1o2 * (htemp[d0MP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[d0MP]);
            h2[dPPP] = c1o2 * (htemp[dPPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dPPP]);
            h2[dMPP] = c1o2 * (htemp[dMPP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dMPP]);
            h2[dPMP] = c1o2 * (htemp[dPMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dPMP]);
            h2[dMMP] = c1o2 * (htemp[dMMP] * (one_over_sqrt3 + vx3) + (c1o1 - one_over_sqrt3 - vx3) * h2[dMMP]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[d00P], x1 + DX1[d00M], x2 + DX2[d00M], x3 + DX3[d00M], d00M);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dP0P], x1 + DX1[dM0M], x2 + DX2[dM0M], x3 + DX3[dM0M], dM0M);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dM0P], x1 + DX1[dP0M], x2 + DX2[dP0M], x3 + DX3[dP0M], dP0M);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0PP], x1 + DX1[d0MM], x2 + DX2[d0MM], x3 + DX3[d0MM], d0MM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0MP], x1 + DX1[d0PM], x2 + DX2[d0PM], x3 + DX3[d0PM], d0PM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPP], x1 + DX1[dMMM], x2 + DX2[dMMM], x3 + DX3[dMMM], dMMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPP], x1 + DX1[dPMM], x2 + DX2[dPMM], x3 + DX3[dPMM], dPMM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMP], x1 + DX1[dMPM], x2 + DX2[dMPM], x3 + DX3[dMPM], dMPM);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMP], x1 + DX1[dPPM], x2 + DX2[dPPM], x3 + DX3[dPPM], dPPM);

            break;
        case d00M:
            f[d00M] = ftemp[d00M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[d00M];
            f[dP0M] = ftemp[dP0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dP0M];
            f[dM0M] = ftemp[dM0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dM0M];
            f[d0PM] = ftemp[d0PM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[d0PM];
            f[d0MM] = ftemp[d0MM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[d0MM];
            f[dPPM] = ftemp[dPPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dPPM];
            f[dMPM] = ftemp[dMPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dMPM];
            f[dPMM] = ftemp[dPMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dPMM];
            f[dMMM] = ftemp[dMMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * f[dMMM];

            distributions->setPreCollisionDistributionForDirection(f[d00M], x1 + DX1[d00P], x2 + DX2[d00P], x3 + DX3[d00P], d00P);
            distributions->setPreCollisionDistributionForDirection(f[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributions->setPreCollisionDistributionForDirection(f[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributions->setPreCollisionDistributionForDirection(f[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributions->setPreCollisionDistributionForDirection(f[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributions->setPreCollisionDistributionForDirection(f[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributions->setPreCollisionDistributionForDirection(f[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributions->setPreCollisionDistributionForDirection(f[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributions->setPreCollisionDistributionForDirection(f[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h[d00M] = htemp[d00M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[d00M];
            h[dP0M] = htemp[dP0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dP0M];
            h[dM0M] = htemp[dM0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dM0M];
            h[d0PM] = htemp[d0PM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[d0PM];
            h[d0MM] = htemp[d0MM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[d0MM];
            h[dPPM] = htemp[dPPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dPPM];
            h[dMPM] = htemp[dMPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dMPM];
            h[dPMM] = htemp[dPMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dPMM];
            h[dMMM] = htemp[dMMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h[dMMM];

            distributionsH->setPreCollisionDistributionForDirection(h[d00M], x1 + DX1[d00P], x2 + DX2[d00P], x3 + DX3[d00P], d00P);
            distributionsH->setPreCollisionDistributionForDirection(h[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributionsH->setPreCollisionDistributionForDirection(h[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributionsH->setPreCollisionDistributionForDirection(h[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributionsH->setPreCollisionDistributionForDirection(h[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributionsH->setPreCollisionDistributionForDirection(h[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH->setPreCollisionDistributionForDirection(h[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributionsH->setPreCollisionDistributionForDirection(h[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributionsH->setPreCollisionDistributionForDirection(h[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            h2[d00M] = c1o2 * (htemp[d00M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[d00M]);
            h2[dP0M] = c1o2 * (htemp[dP0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dP0M]);
            h2[dM0M] = c1o2 * (htemp[dM0M] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dM0M]);
            h2[d0PM] = c1o2 * (htemp[d0PM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[d0PM]);
            h2[d0MM] = c1o2 * (htemp[d0MM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[d0MM]);
            h2[dPPM] = c1o2 * (htemp[dPPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dPPM]);
            h2[dMPM] = c1o2 * (htemp[dMPM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dMPM]);
            h2[dPMM] = c1o2 * (htemp[dPMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dPMM]);
            h2[dMMM] = c1o2 * (htemp[dMMM] * (one_over_sqrt3 - vx3) + (c1o1 - one_over_sqrt3 + vx3) * h2[dMMM]);

            distributionsH2->setPreCollisionDistributionForDirection(h2[d00M], x1 + DX1[d00P], x2 + DX2[d00P], x3 + DX3[d00P], d00P);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dP0M], x1 + DX1[dM0P], x2 + DX2[dM0P], x3 + DX3[dM0P], dM0P);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dM0M], x1 + DX1[dP0P], x2 + DX2[dP0P], x3 + DX3[dP0P], dP0P);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0PM], x1 + DX1[d0MP], x2 + DX2[d0MP], x3 + DX3[d0MP], d0MP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[d0MM], x1 + DX1[d0PP], x2 + DX2[d0PP], x3 + DX3[d0PP], d0PP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPPM], x1 + DX1[dMMP], x2 + DX2[dMMP], x3 + DX3[dMMP], dMMP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMPM], x1 + DX1[dPMP], x2 + DX2[dPMP], x3 + DX3[dPMP], dPMP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dPMM], x1 + DX1[dMPP], x2 + DX2[dMPP], x3 + DX3[dMPP], dMPP);
            distributionsH2->setPreCollisionDistributionForDirection(h2[dMMM], x1 + DX1[dPPP], x2 + DX2[dPPP], x3 + DX3[dPPP], dPPP);

            break;
        default:
            UB_THROW(UbException(UB_EXARGS, "It isn't implemented non reflecting density boundary for this direction!"));
    }
}
