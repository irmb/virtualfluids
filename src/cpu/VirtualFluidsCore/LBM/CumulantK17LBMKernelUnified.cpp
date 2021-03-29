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
//! \file CumulantK17LBMKernel.cpp
//! \ingroup LBM
//! \author Konstantin Kutscher, Martin Geier
//=======================================================================================
#include "CumulantK17LBMKernelUnified.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <cmath>
#include "DataSet3D.h"
#include "LBMKernel.h"
#include "Block3D.h"
#include "BCArray3D.h"

//#include <lbm/CumulantChimera.h>
#include <lbm/CumulantChimeraPreCompiled.h>

//#define PROOF_CORRECTNESS

using namespace UbMath;

//////////////////////////////////////////////////////////////////////////
CumulantK17LBMKernelUnified::CumulantK17LBMKernelUnified()
{
    this->compressible = true;
}
//////////////////////////////////////////////////////////////////////////
void CumulantK17LBMKernelUnified::initDataSet()
{
    SPtr<DistributionArray3D> d(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> CumulantK17LBMKernelUnified::clone()
{
    SPtr<LBMKernel> kernel(new CumulantK17LBMKernelUnified());
    kernel->setNX(nx);
    std::dynamic_pointer_cast<CumulantK17LBMKernelUnified>(kernel)->initDataSet();
    kernel->setCollisionFactor(this->collFactor);
    kernel->setBCProcessor(bcProcessor->clone(kernel));
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
void CumulantK17LBMKernelUnified::calculate(int step)
{
    //////////////////////////////////////////////////////////////////////////
    //! Cumulant K17 Kernel is based on
    //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
    //! and
    //! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
    //!
    //! The cumulant kernel is executed in the following steps
    //!
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from thredIdx, blockIdx, blockDim and gridDim.
    //!

    using namespace std;

    //initializing of forcing stuff
    if (withForcing)
    {
        muForcingX1.DefineVar("x1", &muX1); muForcingX1.DefineVar("x2", &muX2); muForcingX1.DefineVar("x3", &muX3);
        muForcingX2.DefineVar("x1", &muX1); muForcingX2.DefineVar("x2", &muX2); muForcingX2.DefineVar("x3", &muX3);
        muForcingX3.DefineVar("x1", &muX1); muForcingX3.DefineVar("x2", &muX2); muForcingX3.DefineVar("x3", &muX3);

        muDeltaT = deltaT;

        muForcingX1.DefineVar("dt", &muDeltaT);
        muForcingX2.DefineVar("dt", &muDeltaT);
        muForcingX3.DefineVar("dt", &muDeltaT);

        muNu = (1.0 / 3.0) * (1.0 / collFactor - 1.0 / 2.0);

        muForcingX1.DefineVar("nu", &muNu);
        muForcingX2.DefineVar("nu", &muNu);
        muForcingX3.DefineVar("nu", &muNu);
    }
    /////////////////////////////////////

    localDistributions = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
    nonLocalDistributions = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
    restDistributions = dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

    SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();

    const int bcArrayMaxX1 = (int)bcArray->getNX1();
    const int bcArrayMaxX2 = (int)bcArray->getNX2();
    const int bcArrayMaxX3 = (int)bcArray->getNX3();

    int minX1 = ghostLayerWidth;
    int minX2 = ghostLayerWidth;
    int minX3 = ghostLayerWidth;
    int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
    int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
    int maxX3 = bcArrayMaxX3 - ghostLayerWidth;

    LBMReal omega = collFactor;

    for (int x3 = minX3; x3 < maxX3; x3++)
    {
        for (int x2 = minX2; x2 < maxX2; x2++)
        {
            for (int x1 = minX1; x1 < maxX1; x1++)
            {
                if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3))
                {
                    int x1p = x1 + 1;
                    int x2p = x2 + 1;
                    int x3p = x3 + 1;
                    //////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////
                    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm
                    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
                    //!
                    ////////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////

                    //E   N  T
                    //c   c  c
                    //////////
                    //W   S  B
                    //a   a  a

                    //Rest is b

                    //mfxyz
                    //a - negative
                    //b - null
                    //c - positive

                    // a b c
                    //-1 0 1

                    LBMReal mfcbb = (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
                    LBMReal mfbcb = (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
                    LBMReal mfbbc = (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
                    LBMReal mfccb = (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
                    LBMReal mfacb = (*this->localDistributions)(D3Q27System::ET_NW, x1p, x2, x3);
                    LBMReal mfcbc = (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
                    LBMReal mfabc = (*this->localDistributions)(D3Q27System::ET_TW, x1p, x2, x3);
                    LBMReal mfbcc = (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
                    LBMReal mfbac = (*this->localDistributions)(D3Q27System::ET_TS, x1, x2p, x3);
                    LBMReal mfccc = (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
                    LBMReal mfacc = (*this->localDistributions)(D3Q27System::ET_TNW, x1p, x2, x3);
                    LBMReal mfcac = (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2p, x3);
                    LBMReal mfaac = (*this->localDistributions)(D3Q27System::ET_TSW, x1p, x2p, x3);

                    LBMReal mfabb = (*this->nonLocalDistributions)(D3Q27System::ET_W, x1p, x2, x3);
                    LBMReal mfbab = (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2p, x3);
                    LBMReal mfbba = (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3p);
                    LBMReal mfaab = (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1p, x2p, x3);
                    LBMReal mfcab = (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2p, x3);
                    LBMReal mfaba = (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1p, x2, x3p);
                    LBMReal mfcba = (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3p);
                    LBMReal mfbaa = (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2p, x3p);
                    LBMReal mfbca = (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3p);
                    LBMReal mfaaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                    LBMReal mfcaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2p, x3p);
                    LBMReal mfaca = (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1p, x2, x3p);
                    LBMReal mfcca = (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3p);

                    LBMReal mfbbb = (*this->restDistributions)(x1, x2, x3);

                    
                    LBMReal forces[3] = {0., 0., 0.};
                    if (withForcing)
                    {
                        muX1 = static_cast<double>(x1 - 1 + ix1 * maxX1);
                        muX2 = static_cast<double>(x2 - 1 + ix2 * maxX2);
                        muX3 = static_cast<double>(x3 - 1 + ix3 * maxX3);

                        forcingX1 = muForcingX1.Eval();
                        forcingX2 = muForcingX2.Eval();
                        forcingX3 = muForcingX3.Eval();

                        forces[0] = forcingX1 * deltaT;
                        forces[1] = forcingX2 * deltaT;
                        forces[2] = forcingX3 * deltaT;
                    }

                    VF::LBM::Distribution27 distribution;

                    distribution.f[VF::LBM::DIR::PZZ] = mfcbb;
                    distribution.f[VF::LBM::DIR::MZZ] = mfabb;
                    distribution.f[VF::LBM::DIR::ZPZ] = mfbcb;
                    distribution.f[VF::LBM::DIR::ZMZ] = mfbab;
                    distribution.f[VF::LBM::DIR::ZZP] = mfbbc;
                    distribution.f[VF::LBM::DIR::ZZM] = mfbba;
                    distribution.f[VF::LBM::DIR::PPZ] = mfccb;
                    distribution.f[VF::LBM::DIR::MMZ] = mfaab;
                    distribution.f[VF::LBM::DIR::PMZ] = mfcab;
                    distribution.f[VF::LBM::DIR::MPZ] = mfacb;
                    distribution.f[VF::LBM::DIR::PZP] = mfcbc;
                    distribution.f[VF::LBM::DIR::MZM] = mfaba;
                    distribution.f[VF::LBM::DIR::PZM] = mfcba;
                    distribution.f[VF::LBM::DIR::MZP] = mfabc;
                    distribution.f[VF::LBM::DIR::ZPP] = mfbcc;
                    distribution.f[VF::LBM::DIR::ZMM] = mfbaa;
                    distribution.f[VF::LBM::DIR::ZPM] = mfbca;
                    distribution.f[VF::LBM::DIR::ZMP] = mfbac;
                    distribution.f[VF::LBM::DIR::PPP] = mfccc;
                    distribution.f[VF::LBM::DIR::MPP] = mfacc;
                    distribution.f[VF::LBM::DIR::PMP] = mfcac;
                    distribution.f[VF::LBM::DIR::MMP] = mfaac;
                    distribution.f[VF::LBM::DIR::PPM] = mfcca;
                    distribution.f[VF::LBM::DIR::MPM] = mfaca;
                    distribution.f[VF::LBM::DIR::PMM] = mfcaa;
                    distribution.f[VF::LBM::DIR::MMM] = mfaaa;
                    distribution.f[VF::LBM::DIR::ZZZ] = mfbbb;

                    VF::LBM::cumulantChimera(distribution, omega, forces);

                    mfcbb = distribution.f[VF::LBM::DIR::PZZ];
                    mfabb = distribution.f[VF::LBM::DIR::MZZ];
                    mfbcb = distribution.f[VF::LBM::DIR::ZPZ];
                    mfbab = distribution.f[VF::LBM::DIR::ZMZ];
                    mfbbc = distribution.f[VF::LBM::DIR::ZZP];
                    mfbba = distribution.f[VF::LBM::DIR::ZZM];
                    mfccb = distribution.f[VF::LBM::DIR::PPZ];
                    mfaab = distribution.f[VF::LBM::DIR::MMZ];
                    mfcab = distribution.f[VF::LBM::DIR::PMZ];
                    mfacb = distribution.f[VF::LBM::DIR::MPZ];
                    mfcbc = distribution.f[VF::LBM::DIR::PZP];
                    mfaba = distribution.f[VF::LBM::DIR::MZM];
                    mfcba = distribution.f[VF::LBM::DIR::PZM];
                    mfabc = distribution.f[VF::LBM::DIR::MZP];
                    mfbcc = distribution.f[VF::LBM::DIR::ZPP];
                    mfbaa = distribution.f[VF::LBM::DIR::ZMM];
                    mfbca = distribution.f[VF::LBM::DIR::ZPM];
                    mfbac = distribution.f[VF::LBM::DIR::ZMP];
                    mfccc = distribution.f[VF::LBM::DIR::PPP];
                    mfacc = distribution.f[VF::LBM::DIR::MPP];
                    mfcac = distribution.f[VF::LBM::DIR::PMP];
                    mfaac = distribution.f[VF::LBM::DIR::MMP];
                    mfcca = distribution.f[VF::LBM::DIR::PPM];
                    mfaca = distribution.f[VF::LBM::DIR::MPM];
                    mfcaa = distribution.f[VF::LBM::DIR::PMM];
                    mfaaa = distribution.f[VF::LBM::DIR::MMM];
                    mfbbb = distribution.f[VF::LBM::DIR::ZZZ];

                    //////////////////////////////////////////////////////////////////////////
                    //proof correctness
                    //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
                    LBMReal drho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
                                        + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
                                        + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
                    LBMReal dif = distribution.getDensity_() - drho_post;
#ifdef SINGLEPRECISION
                    if (dif > 10.0E-7 || dif < -10.0E-7)
#else
                    if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
                    {
                        UB_THROW(UbException(UB_EXARGS, "rho=" + UbSystem::toString(distribution.getDensity_()) + ", rho_post=" + UbSystem::toString(drho_post)
                                                        + " dif=" + UbSystem::toString(dif)
                                                        + " rho is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)
                                                        + " in " + block.lock()->toString() + " step = " + UbSystem::toString(step)));
                    }
#endif
                    
                    (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3) = mfcbb;
                    (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3) = mfbcb;
                    (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3) = mfbbc;
                    (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3) = mfccb;
                    (*this->localDistributions)(D3Q27System::ET_NW, x1p, x2, x3) = mfacb;
                    (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3) = mfcbc;
                    (*this->localDistributions)(D3Q27System::ET_TW, x1p, x2, x3) = mfabc;
                    (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3) = mfbcc;
                    (*this->localDistributions)(D3Q27System::ET_TS, x1, x2p, x3) = mfbac;
                    (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3) = mfccc;
                    (*this->localDistributions)(D3Q27System::ET_TNW, x1p, x2, x3) = mfacc;
                    (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2p, x3) = mfcac;
                    (*this->localDistributions)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfaac;

                    (*this->nonLocalDistributions)(D3Q27System::ET_W, x1p, x2, x3) =  mfabb;
                    (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2p, x3) =  mfbab;
                    (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3p) =  mfbba;
                    (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1p, x2p, x3) = mfaab;
                    (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2p, x3) = mfcab;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1p, x2, x3p) = mfaba;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3p) = mfcba;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbaa;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3p) = mfbca;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfaaa;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfcaa;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfaca;
                    (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3p) = mfcca;
                    (*this->restDistributions)(x1, x2, x3) = mfbbb;
                    //////////////////////////////////////////////////////////////////////////

                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////

