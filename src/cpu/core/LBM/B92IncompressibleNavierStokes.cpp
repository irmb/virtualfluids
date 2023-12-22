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
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "B92IncompressibleNavierStokes.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "EsoSplit.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "Block3D.h"
#include "basics/constants/NumericConstants.h"

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
B92IncompressibleNavierStokes::B92IncompressibleNavierStokes() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
B92IncompressibleNavierStokes::~B92IncompressibleNavierStokes(void) = default;
//////////////////////////////////////////////////////////////////////////
void B92IncompressibleNavierStokes::initDataSet()
{
    SPtr<DistributionArray3D> d(new EsoSplit(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> B92IncompressibleNavierStokes::clone()
{
    SPtr<LBMKernel> kernel(new B92IncompressibleNavierStokes());
    kernel->setNX(nx);
    std::dynamic_pointer_cast<B92IncompressibleNavierStokes>(kernel)->initDataSet();
    kernel->setCollisionFactor(this->collFactor);
    kernel->setBCSet(bcSet->clone(kernel));
    kernel->setWithForcing(withForcing);
    kernel->setForcingX1(muForcingX1);
    kernel->setForcingX2(muForcingX2);
    kernel->setForcingX3(muForcingX3);
    kernel->setIndex(ix1, ix2, ix3);
    kernel->setDeltaT(deltaT);
    kernel->setBlock(block.lock());
    return kernel;
}
//////////////////////////////////////////////////////////////////////////
void B92IncompressibleNavierStokes::calculate(int step)
{
    using namespace D3Q27System;
 //   using namespace UbMath;
   using namespace vf::basics::constant;
   using namespace vf::lbm::dir;

    // initializing of forcing stuff
    if (withForcing) {
        muForcingX1.DefineVar("x1", &muX1);
        muForcingX1.DefineVar("x2", &muX2);
        muForcingX1.DefineVar("x3", &muX3);
        muForcingX2.DefineVar("x1", &muX1);
        muForcingX2.DefineVar("x2", &muX2);
        muForcingX2.DefineVar("x3", &muX3);
        muForcingX3.DefineVar("x1", &muX1);
        muForcingX3.DefineVar("x2", &muX2);
        muForcingX3.DefineVar("x3", &muX3);
        forcingX1 = 0;
        forcingX2 = 0;
        forcingX3 = 0;
    }
    /////////////////////////////////////

    localDistributions =
        std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())->getLocalDistributions();
    nonLocalDistributions = std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())
                                ->getNonLocalDistributions();
    zeroDistributions =
        std::dynamic_pointer_cast<EsoSplit>(dataSet->getFdistributions())->getZeroDistributions();

    SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();
    real f[D3Q27System::ENDF + 1];
    real feq[D3Q27System::ENDF + 1];
    real drho, vx1, vx2, vx3;
    const int bcArrayMaxX1 = (int)bcArray->getNX1();
    const int bcArrayMaxX2 = (int)bcArray->getNX2();
    const int bcArrayMaxX3 = (int)bcArray->getNX3();

    int minX1 = ghostLayerWidth;
    int minX2 = ghostLayerWidth;
    int minX3 = ghostLayerWidth;
    int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
    int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
    int maxX3 = bcArrayMaxX3 - ghostLayerWidth;

    for (int x3 = minX3; x3 < maxX3; x3++) {
        for (int x2 = minX2; x2 < maxX2; x2++) {
            for (int x1 = minX1; x1 < maxX1; x1++) {
                if (bcArray->isUnvalidForCollision(x1, x2, x3)) {
                    continue;
                }

                int x1p = x1 + 1;
                int x2p = x2 + 1;
                int x3p = x3 + 1;
                //////////////////////////////////////////////////////////////////////////
                // read distribution
                ////////////////////////////////////////////////////////////////////////////
                f[d000] = (*this->zeroDistributions)(x1, x2, x3);

                f[dP00]   = (*this->localDistributions)(eP00, x1, x2, x3);
                f[d0P0]   = (*this->localDistributions)(e0P0, x1, x2, x3);
                f[d00P]   = (*this->localDistributions)(e00P, x1, x2, x3);
                f[dPP0]  = (*this->localDistributions)(ePP0, x1, x2, x3);
                f[dMP0]  = (*this->localDistributions)(eMP0, x1p, x2, x3);
                f[dP0P]  = (*this->localDistributions)(eP0P, x1, x2, x3);
                f[dM0P]  = (*this->localDistributions)(eM0P, x1p, x2, x3);
                f[d0PP]  = (*this->localDistributions)(e0PP, x1, x2, x3);
                f[d0MP]  = (*this->localDistributions)(e0MP, x1, x2p, x3);
                f[dPPP] = (*this->localDistributions)(ePPP, x1, x2, x3);
                f[dMPP] = (*this->localDistributions)(eMPP, x1p, x2, x3);
                f[dPMP] = (*this->localDistributions)(ePMP, x1, x2p, x3);
                f[dMMP] = (*this->localDistributions)(eMMP, x1p, x2p, x3);

                f[dM00]   = (*this->nonLocalDistributions)(eM00, x1p, x2, x3);
                f[d0M0]   = (*this->nonLocalDistributions)(e0M0, x1, x2p, x3);
                f[d00M]   = (*this->nonLocalDistributions)(e00M, x1, x2, x3p);
                f[dMM0]  = (*this->nonLocalDistributions)(eMM0, x1p, x2p, x3);
                f[dPM0]  = (*this->nonLocalDistributions)(ePM0, x1, x2p, x3);
                f[dM0M]  = (*this->nonLocalDistributions)(eM0M, x1p, x2, x3p);
                f[dP0M]  = (*this->nonLocalDistributions)(eP0M, x1, x2, x3p);
                f[d0MM]  = (*this->nonLocalDistributions)(e0MM, x1, x2p, x3p);
                f[d0PM]  = (*this->nonLocalDistributions)(e0PM, x1, x2, x3p);
                f[dMMM] = (*this->nonLocalDistributions)(eMMM, x1p, x2p, x3p);
                f[dPMM] = (*this->nonLocalDistributions)(ePMM, x1, x2p, x3p);
                f[dMPM] = (*this->nonLocalDistributions)(eMPM, x1p, x2, x3p);
                f[dPPM] = (*this->nonLocalDistributions)(ePPM, x1, x2, x3p);
                //////////////////////////////////////////////////////////////////////////

                drho = f[d000] + f[dP00] + f[dM00] + f[d0P0] + f[d0M0] + f[d00P] + f[d00M] + f[dPP0] + f[dMM0] + f[dPM0] + f[dMP0] + f[dP0P] +
                        f[dM0M] + f[dP0M] + f[dM0P] + f[d0PP] + f[d0MM] + f[d0PM] + f[d0MP] + f[dPPP] + f[dMMP] + f[dPMP] + f[dMPP] +
                        f[dPPM] + f[dMMM] + f[dPMM] + f[dMPM];

                vx1 = f[dP00] - f[dM00] + f[dPP0] - f[dMM0] + f[dPM0] - f[dMP0] + f[dP0P] - f[dM0M] + f[dP0M] - f[dM0P] + f[dPPP] -
                        f[dMMP] + f[dPMP] - f[dMPP] + f[dPPM] - f[dMMM] + f[dPMM] - f[dMPM];

                vx2 = f[d0P0] - f[d0M0] + f[dPP0] - f[dMM0] - f[dPM0] + f[dMP0] + f[d0PP] - f[d0MM] + f[d0PM] - f[d0MP] + f[dPPP] -
                        f[dMMP] - f[dPMP] + f[dMPP] + f[dPPM] - f[dMMM] - f[dPMM] + f[dMPM];

                vx3 = f[d00P] - f[d00M] + f[dP0P] - f[dM0M] - f[dP0M] + f[dM0P] + f[d0PP] - f[d0MM] - f[d0PM] + f[d0MP] + f[dPPP] +
                        f[dMMP] + f[dPMP] + f[dMPP] - f[dPPM] - f[dMMM] - f[dPMM] - f[dMPM];

                real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

                feq[d000] = c8o27 * (drho - cu_sq);
                feq[dP00]    = c2o27 * (drho + c3o1 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq);
                feq[dM00]    = c2o27 * (drho + c3o1 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq);
                feq[d0P0]    = c2o27 * (drho + c3o1 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq);
                feq[d0M0]    = c2o27 * (drho + c3o1 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq);
                feq[d00P]    = c2o27 * (drho + c3o1 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq);
                feq[d00M]    = c2o27 * (drho + c3o1 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq);
                feq[dPP0]   = c1o54 * (drho + c3o1 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
                feq[dMM0]   = c1o54 * (drho + c3o1 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
                feq[dPM0]   = c1o54 * (drho + c3o1 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
                feq[dMP0]   = c1o54 * (drho + c3o1 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
                feq[dP0P]   = c1o54 * (drho + c3o1 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
                feq[dM0M]   = c1o54 * (drho + c3o1 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
                feq[dP0M]   = c1o54 * (drho + c3o1 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
                feq[dM0P]   = c1o54 * (drho + c3o1 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
                feq[d0PP]   = c1o54 * (drho + c3o1 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
                feq[d0MM]   = c1o54 * (drho + c3o1 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
                feq[d0PM]   = c1o54 * (drho + c3o1 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
                feq[d0MP]   = c1o54 * (drho + c3o1 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
                feq[dPPP]  = c1o216 *
                            (drho + c3o1 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
                feq[dMMM] = c1o216 * (drho + c3o1 * (-vx1 - vx2 - vx3) +
                                        c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
                feq[dPPM] = c1o216 *
                            (drho + c3o1 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
                feq[dMMP] = c1o216 * (drho + c3o1 * (-vx1 - vx2 + vx3) +
                                        c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
                feq[dPMP] = c1o216 *
                            (drho + c3o1 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
                feq[dMPM] = c1o216 * (drho + c3o1 * (-vx1 + vx2 - vx3) +
                                        c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
                feq[dPMM] = c1o216 *
                            (drho + c3o1 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
                feq[dMPP] = c1o216 * (drho + c3o1 * (-vx1 + vx2 + vx3) +
                                        c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);

                // Relaxation
                f[d000] += (feq[d000] - f[d000]) * collFactor;
                f[dP00] += (feq[dP00] - f[dP00]) * collFactor;
                f[dM00] += (feq[dM00] - f[dM00]) * collFactor;
                f[d0P0] += (feq[d0P0] - f[d0P0]) * collFactor;
                f[d0M0] += (feq[d0M0] - f[d0M0]) * collFactor;
                f[d00P] += (feq[d00P] - f[d00P]) * collFactor;
                f[d00M] += (feq[d00M] - f[d00M]) * collFactor;
                f[dPP0] += (feq[dPP0] - f[dPP0]) * collFactor;
                f[dMM0] += (feq[dMM0] - f[dMM0]) * collFactor;
                f[dPM0] += (feq[dPM0] - f[dPM0]) * collFactor;
                f[dMP0] += (feq[dMP0] - f[dMP0]) * collFactor;
                f[dP0P] += (feq[dP0P] - f[dP0P]) * collFactor;
                f[dM0M] += (feq[dM0M] - f[dM0M]) * collFactor;
                f[dP0M] += (feq[dP0M] - f[dP0M]) * collFactor;
                f[dM0P] += (feq[dM0P] - f[dM0P]) * collFactor;
                f[d0PP] += (feq[d0PP] - f[d0PP]) * collFactor;
                f[d0MM] += (feq[d0MM] - f[d0MM]) * collFactor;
                f[d0PM] += (feq[d0PM] - f[d0PM]) * collFactor;
                f[d0MP] += (feq[d0MP] - f[d0MP]) * collFactor;

                f[dPPP] += (feq[dPPP] - f[dPPP]) * collFactor;
                f[dMMM] += (feq[dMMM] - f[dMMM]) * collFactor;
                f[dPPM] += (feq[dPPM] - f[dPPM]) * collFactor;
                f[dMMP] += (feq[dMMP] - f[dMMP]) * collFactor;
                f[dPMP] += (feq[dPMP] - f[dPMP]) * collFactor;
                f[dMPM] += (feq[dMPM] - f[dMPM]) * collFactor;
                f[dPMM] += (feq[dPMM] - f[dPMM]) * collFactor;
                f[dMPP] += (feq[dMPP] - f[dMPP]) * collFactor;

                //////////////////////////////////////////////////////////////////////////
                // forcing
                if (withForcing) {
                    muX1 = x1 + ix1 * bcArrayMaxX1;
                    muX2 = x2 + ix2 * bcArrayMaxX2;
                    muX3 = x3 + ix3 * bcArrayMaxX3;

                    forcingX1 = muForcingX1.Eval();
                    forcingX2 = muForcingX2.Eval();
                    forcingX3 = muForcingX3.Eval();

                    f[d000] += c0o1;
                    f[dP00] += c3o1 * c2o27 * (forcingX1);
                    f[dM00] += c3o1 * c2o27 * (-forcingX1);
                    f[d0P0] += c3o1 * c2o27 * (forcingX2);
                    f[d0M0] += c3o1 * c2o27 * (-forcingX2);
                    f[d00P] += c3o1 * c2o27 * (forcingX3);
                    f[d00M] += c3o1 * c2o27 * (-forcingX3);
                    f[dPP0] += c3o1 * c1o54 * (forcingX1 + forcingX2);
                    f[dMM0] += c3o1 * c1o54 * (-forcingX1 - forcingX2);
                    f[dPM0] += c3o1 * c1o54 * (forcingX1 - forcingX2);
                    f[dMP0] += c3o1 * c1o54 * (-forcingX1 + forcingX2);
                    f[dP0P] += c3o1 * c1o54 * (forcingX1 + forcingX3);
                    f[dM0M] += c3o1 * c1o54 * (-forcingX1 - forcingX3);
                    f[dP0M] += c3o1 * c1o54 * (forcingX1 - forcingX3);
                    f[dM0P] += c3o1 * c1o54 * (-forcingX1 + forcingX3);
                    f[d0PP] += c3o1 * c1o54 * (forcingX2 + forcingX3);
                    f[d0MM] += c3o1 * c1o54 * (-forcingX2 - forcingX3);
                    f[d0PM] += c3o1 * c1o54 * (forcingX2 - forcingX3);
                    f[d0MP] += c3o1 * c1o54 * (-forcingX2 + forcingX3);
                    f[dPPP] += c3o1 * c1o216 * (forcingX1 + forcingX2 + forcingX3);
                    f[dMMM] += c3o1 * c1o216 * (-forcingX1 - forcingX2 - forcingX3);
                    f[dPPM] += c3o1 * c1o216 * (forcingX1 + forcingX2 - forcingX3);
                    f[dMMP] += c3o1 * c1o216 * (-forcingX1 - forcingX2 + forcingX3);
                    f[dPMP] += c3o1 * c1o216 * (forcingX1 - forcingX2 + forcingX3);
                    f[dMPM] += c3o1 * c1o216 * (-forcingX1 + forcingX2 - forcingX3);
                    f[dPMM] += c3o1 * c1o216 * (forcingX1 - forcingX2 - forcingX3);
                    f[dMPP] += c3o1 * c1o216 * (-forcingX1 + forcingX2 + forcingX3);
                }
                //////////////////////////////////////////////////////////////////////////
#ifdef PROOF_CORRECTNESS
                real rho_post = f[d000] + f[dP00] + f[dM00] + f[d0P0] + f[d0M0] + f[d00P] + f[d00M] + f[dPP0] + f[dMM0] + f[dPM0] +
                                    f[dMP0] + f[dP0P] + f[dM0M] + f[dP0M] + f[dM0P] + f[d0PP] + f[d0MM] + f[d0PM] + f[d0MP] + f[dPPP] +
                                    f[dMMP] + f[dPMP] + f[dMPP] + f[dPPM] + f[dMMM] + f[dPMM] + f[dMPM];
                real dif = drho - rho_post;
#ifdef SINGLEPRECISION
                if (dif > 10.0E-7 || dif < -10.0E-7)
#else
                if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
                {
                    UB_THROW(UbException(UB_EXARGS, "rho="+UbSystem::toString(drho)+", rho_post="+UbSystem::toString(rho_post)
                        +" dif="+UbSystem::toString(dif)
                        +" rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)
                        +" in " + block.lock()->toString()+" step = "+UbSystem::toString(step)));
                }
#endif
                //////////////////////////////////////////////////////////////////////////
                // write distribution
                //////////////////////////////////////////////////////////////////////////
                (*this->localDistributions)(eP00, x1, x2, x3)     = f[iP00];
                (*this->localDistributions)(e0P0, x1, x2, x3)     = f[i0P0];
                (*this->localDistributions)(e00P, x1, x2, x3)     = f[i00P];
                (*this->localDistributions)(ePP0, x1, x2, x3)    = f[iPP0];
                (*this->localDistributions)(eMP0, x1p, x2, x3)   = f[iMP0];
                (*this->localDistributions)(eP0P, x1, x2, x3)    = f[iP0P];
                (*this->localDistributions)(eM0P, x1p, x2, x3)   = f[iM0P];
                (*this->localDistributions)(e0PP, x1, x2, x3)    = f[i0PP];
                (*this->localDistributions)(e0MP, x1, x2p, x3)   = f[i0MP];
                (*this->localDistributions)(ePPP, x1, x2, x3)   = f[iPPP];
                (*this->localDistributions)(eMPP, x1p, x2, x3)  = f[iMPP];
                (*this->localDistributions)(ePMP, x1, x2p, x3)  = f[iPMP];
                (*this->localDistributions)(eMMP, x1p, x2p, x3) = f[iMMP];

                (*this->nonLocalDistributions)(eM00, x1p, x2, x3)     = f[iM00];
                (*this->nonLocalDistributions)(e0M0, x1, x2p, x3)     = f[i0M0];
                (*this->nonLocalDistributions)(e00M, x1, x2, x3p)     = f[i00M];
                (*this->nonLocalDistributions)(eMM0, x1p, x2p, x3)   = f[iMM0];
                (*this->nonLocalDistributions)(ePM0, x1, x2p, x3)    = f[iPM0];
                (*this->nonLocalDistributions)(eM0M, x1p, x2, x3p)   = f[iM0M];
                (*this->nonLocalDistributions)(eP0M, x1, x2, x3p)    = f[iP0M];
                (*this->nonLocalDistributions)(e0MM, x1, x2p, x3p)   = f[i0MM];
                (*this->nonLocalDistributions)(e0PM, x1, x2, x3p)    = f[i0PM];
                (*this->nonLocalDistributions)(eMMM, x1p, x2p, x3p) = f[iMMM];
                (*this->nonLocalDistributions)(ePMM, x1, x2p, x3p)  = f[iPMM];
                (*this->nonLocalDistributions)(eMPM, x1p, x2, x3p)  = f[iMPM];
                (*this->nonLocalDistributions)(ePPM, x1, x2, x3p)   = f[iPPM];

                (*this->zeroDistributions)(x1, x2, x3) = f[d000];
                //////////////////////////////////////////////////////////////////////////
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
real B92IncompressibleNavierStokes::getCalculationTime() { return vf::basics::constant::c0o1; }

//! \}
