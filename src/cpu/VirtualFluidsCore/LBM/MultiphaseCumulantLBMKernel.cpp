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
//! \file MultiphaseCumulantLBMKernel.cpp
//! \ingroup LBMKernel
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseCumulantLBMKernel.h"
#include "BCArray3D.h"
#include "Block3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
#include <cmath>

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
MultiphaseCumulantLBMKernel::MultiphaseCumulantLBMKernel() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
void MultiphaseCumulantLBMKernel::initDataSet()
{
    SPtr<DistributionArray3D> f(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    SPtr<DistributionArray3D> h(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9)); // For phase-field
    SPtr<PhaseFieldArray3D> divU(new PhaseFieldArray3D(nx[0] + 2, nx[1] + 2, nx[2] + 2, 0.0));
    dataSet->setFdistributions(f);
    dataSet->setHdistributions(h); // For phase-field
    dataSet->setPhaseField(divU);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> MultiphaseCumulantLBMKernel::clone()
{
    SPtr<LBMKernel> kernel(new MultiphaseCumulantLBMKernel());
    kernel->setNX(nx);
    dynamicPointerCast<MultiphaseCumulantLBMKernel>(kernel)->initDataSet();
    kernel->setCollisionFactorMultiphase(this->collFactorL, this->collFactorG);
    kernel->setDensityRatio(this->densityRatio);
    kernel->setMultiphaseModelParameters(this->beta, this->kappa);
    kernel->setContactAngle(this->contactAngle);
    kernel->setPhiL(this->phiL);
    kernel->setPhiH(this->phiH);
    kernel->setPhaseFieldRelaxation(this->tauH);
    kernel->setMobility(this->mob);

    kernel->setBCProcessor(bcProcessor->clone(kernel));
    kernel->setWithForcing(withForcing);
    kernel->setForcingX1(muForcingX1);
    kernel->setForcingX2(muForcingX2);
    kernel->setForcingX3(muForcingX3);
    kernel->setIndex(ix1, ix2, ix3);
    kernel->setDeltaT(deltaT);

    return kernel;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseCumulantLBMKernel::calculate(int step)
{
    using namespace D3Q27System;
    using namespace UbMath;

    forcingX1 = 0.0;
    forcingX2 = 0.0;
    forcingX3 = 0.0;
    /////////////////////////////////////

    localDistributionsF    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
    nonLocalDistributionsF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
    zeroDistributionsF     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

    localDistributionsH    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getLocalDistributions();
    nonLocalDistributionsH = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getNonLocalDistributions();
    zeroDistributionsH     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getZeroDistributions();

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

        CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr phaseField(
            new CbArray3D<LBMReal, IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, -999.0));
        CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr divU(
            new CbArray3D<LBMReal, IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, 0.0));


        for (int x3 = 0; x3 <= maxX3; x3++) {
            for (int x2 = 0; x2 <= maxX2; x2++) {
                for (int x1 = 0; x1 <= maxX1; x1++) {
                    if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                        int x1p = x1 + 1;
                        int x2p = x2 + 1;
                        int x3p = x3 + 1;

                        LBMReal mfcbb = (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3);
                        LBMReal mfbcb = (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3);
                        LBMReal mfbbc = (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3);
                        LBMReal mfccb = (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3);
                        LBMReal mfacb = (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3);
                        LBMReal mfcbc = (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3);
                        LBMReal mfabc = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3);
                        LBMReal mfbcc = (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3);
                        LBMReal mfbac = (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3);
                        LBMReal mfccc = (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3);
                        LBMReal mfacc = (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3);
                        LBMReal mfcac = (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3);
                        LBMReal mfaac = (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3);
                        LBMReal mfabb = (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3);
                        LBMReal mfbab = (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3);
                        LBMReal mfbba = (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p);
                        LBMReal mfaab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3);
                        LBMReal mfcab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3);
                        LBMReal mfaba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p);
                        LBMReal mfcba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p);
                        LBMReal mfbaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p);
                        LBMReal mfbca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p);
                        LBMReal mfaaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                        LBMReal mfcaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p);
                        LBMReal mfaca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p);
                        LBMReal mfcca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p);

                        LBMReal mfbbb = (*this->zeroDistributionsH)(x1, x2, x3);
                        (*phaseField)(x1, x2, x3) = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
                                                    (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
                                                    (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) +
                                                    (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
                    }
                }
            }
        }

        LBMReal collFactorM;
        LBMReal forcingTerm[D3Q27System::ENDF + 1];

        for (int x3 = minX3; x3 < maxX3; x3++) {
            for (int x2 = minX2; x2 < maxX2; x2++) {
                for (int x1 = minX1; x1 < maxX1; x1++) {
                    if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
                        int x1p = x1 + 1;
                        int x2p = x2 + 1;
                        int x3p = x3 + 1;

                        //////////////////////////////////////////////////////////////////////////
                        // Read distributions and phase field
                        ////////////////////////////////////////////////////////////////////////////
                        //////////////////////////////////////////////////////////////////////////

                        // E   N  T
                        // c   c  c
                        //////////
                        // W   S  B
                        // a   a  a

                        // Rest ist b

                        // mfxyz
                        // a - negative
                        // b - null
                        // c - positive

                        // a b c
                        //-1 0 1

                        findNeighbors(phaseField, x1, x2, x3);

                        LBMReal mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
                        LBMReal mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
                        LBMReal mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
                        LBMReal mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
                        LBMReal mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
                        LBMReal mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
                        LBMReal mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
                        LBMReal mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
                        LBMReal mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
                        LBMReal mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
                        LBMReal mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
                        LBMReal mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
                        LBMReal mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
                        LBMReal mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
                        LBMReal mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
                        LBMReal mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
                        LBMReal mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
                        LBMReal mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
                        LBMReal mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
                        LBMReal mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
                        LBMReal mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
                        LBMReal mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
                        LBMReal mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                        LBMReal mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
                        LBMReal mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
                        LBMReal mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

                        LBMReal mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);

                        LBMReal rhoH = 1.0;
                        LBMReal rhoL = 1.0 / densityRatio;

                        LBMReal rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

                        LBMReal dX1_phi = gradX1_phi();
                        LBMReal dX2_phi = gradX2_phi();
                        LBMReal dX3_phi = gradX3_phi();

                        LBMReal denom = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi) + 1e-9;
                        collFactorM = collFactorL + (collFactorL - collFactorG) * (phi[DIR_000] - phiH) / (phiH - phiL);


                        LBMReal mu = 2 * beta * phi[DIR_000] * (phi[DIR_000] - 1) * (2 * phi[DIR_000] - 1) - kappa * nabla2_phi();

                        //----------- Calculating Macroscopic Values -------------
                        LBMReal rho = rhoH + rhoToPhi * (phi[DIR_000] - phiH);

                        if (withForcing) {
                            // muX1 = static_cast<double>(x1-1+ix1*maxX1);
                            // muX2 = static_cast<double>(x2-1+ix2*maxX2);
                            // muX3 = static_cast<double>(x3-1+ix3*maxX3);

                            muForcingX1.DefineVar("rho",&muRho); 
				            muForcingX2.DefineVar("rho",&muRho); 
				            muForcingX3.DefineVar("rho",&muRho); 

				            muRho = rho;

                            forcingX1 = muForcingX1.Eval();
                            forcingX2 = muForcingX2.Eval();
                            forcingX3 = muForcingX3.Eval();

                            LBMReal rho_m = 1.0 / densityRatio;
                            forcingX1     = forcingX1 * (rho - rho_m);
                            forcingX2     = forcingX2 * (rho - rho_m);
                            forcingX3     = forcingX3 * (rho - rho_m);

                            // ux += forcingX1*deltaT*0.5; // X
                            // uy += forcingX2*deltaT*0.5; // Y
                            // uz += forcingX3*deltaT*0.5; // Z
                        }

                        LBMReal ux = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                                      (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                                      (mfcbb - mfabb)) /
                                         (rho * c1o3) +
                                     (mu * dX1_phi + forcingX1) / (2 * rho);

                        LBMReal uy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                                      (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                                      (mfbcb - mfbab)) /
                                         (rho * c1o3) +
                                     (mu * dX2_phi + forcingX2) / (2 * rho);

                        LBMReal uz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                                      (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                                      (mfbbc - mfbba)) /
                                         (rho * c1o3) +
                                     (mu * dX3_phi + forcingX3) / (2 * rho);

                        //--------------------------------------------------------

                        LBMReal ux2 = ux * ux;
                        LBMReal uy2 = uy * uy;
                        LBMReal uz2 = uz * uz;

                        //----------- Calculating Forcing Terms * -------------
                        for (int dir = STARTF; dir <= (FENDDIR); dir++) {
                            LBMReal velProd = DX1[dir] * ux + DX2[dir] * uy + DX3[dir] * uz;
                            LBMReal velSq1  = velProd * velProd;
                            LBMReal gamma = WEIGTH[dir] * (1.0 + 3 * velProd + 4.5 * velSq1 - 1.5 * (ux2 + uy2 + uz2));

                            LBMReal fac1 = (gamma - WEIGTH[dir]) * c1o3 * rhoToPhi;

                            forcingTerm[dir] = ((-ux) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
                                                (-uy) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
                                                (-uz) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3))) +
                                               (DX1[dir]) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
                                               (DX2[dir]) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
                                               (DX3[dir]) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3));
                        }

                        LBMReal gamma = WEIGTH[DIR_000] * (1.0 - 1.5 * (ux2 + uy2 + uz2));
                        LBMReal fac1      = (gamma - WEIGTH[DIR_000]) * c1o3 * rhoToPhi;
                        forcingTerm[DIR_000] = (-ux) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
                                            (-uy) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
                                            (-uz) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3));

                        //--------------------------------------------------------

                        mfcbb = 3.0 * (mfcbb + 0.5 * forcingTerm[DIR_P00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[E  ];
                        mfbcb = 3.0 * (mfbcb + 0.5 * forcingTerm[DIR_0P0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[N  ];
                        mfbbc = 3.0 * (mfbbc + 0.5 * forcingTerm[DIR_00P]) / rho;    //-(3.0*p1 - rho)*WEIGTH[T  ];
                        mfccb = 3.0 * (mfccb + 0.5 * forcingTerm[DIR_PP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NE ];
                        mfacb = 3.0 * (mfacb + 0.5 * forcingTerm[DIR_MP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NW ];
                        mfcbc = 3.0 * (mfcbc + 0.5 * forcingTerm[DIR_P0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TE ];
                        mfabc = 3.0 * (mfabc + 0.5 * forcingTerm[DIR_M0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TW ];
                        mfbcc = 3.0 * (mfbcc + 0.5 * forcingTerm[DIR_0PP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TN ];
                        mfbac = 3.0 * (mfbac + 0.5 * forcingTerm[DIR_0MP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TS ];
                        mfccc = 3.0 * (mfccc + 0.5 * forcingTerm[DIR_PPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNE];
                        mfacc = 3.0 * (mfacc + 0.5 * forcingTerm[DIR_MPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNW];
                        mfcac = 3.0 * (mfcac + 0.5 * forcingTerm[DIR_PMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSE];
                        mfaac = 3.0 * (mfaac + 0.5 * forcingTerm[DIR_MMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSW];
                        mfabb = 3.0 * (mfabb + 0.5 * forcingTerm[DIR_M00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[W  ];
                        mfbab = 3.0 * (mfbab + 0.5 * forcingTerm[DIR_0M0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[S  ];
                        mfbba = 3.0 * (mfbba + 0.5 * forcingTerm[DIR_00M]) / rho;    //-(3.0*p1 - rho)*WEIGTH[B  ];
                        mfaab = 3.0 * (mfaab + 0.5 * forcingTerm[DIR_MM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SW ];
                        mfcab = 3.0 * (mfcab + 0.5 * forcingTerm[DIR_PM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SE ];
                        mfaba = 3.0 * (mfaba + 0.5 * forcingTerm[DIR_M0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BW ];
                        mfcba = 3.0 * (mfcba + 0.5 * forcingTerm[DIR_P0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BE ];
                        mfbaa = 3.0 * (mfbaa + 0.5 * forcingTerm[DIR_0MM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BS ];
                        mfbca = 3.0 * (mfbca + 0.5 * forcingTerm[DIR_0PM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BN ];
                        mfaaa = 3.0 * (mfaaa + 0.5 * forcingTerm[DIR_MMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSW];
                        mfcaa = 3.0 * (mfcaa + 0.5 * forcingTerm[DIR_PMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSE];
                        mfaca = 3.0 * (mfaca + 0.5 * forcingTerm[DIR_MPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNW];
                        mfcca = 3.0 * (mfcca + 0.5 * forcingTerm[DIR_PPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNE];
                        mfbbb = 3.0 * (mfbbb + 0.5 * forcingTerm[DIR_000]) / rho; //- (3.0*p1 - rho)*WEIGTH[REST];

                        LBMReal rho1 = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
                                       (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
                                       (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) +
                                       (mfbba + mfbbc) + mfbbb;


                        LBMReal oMdrho, m0, m1, m2;

                        oMdrho = mfccc + mfaaa;
                        m0     = mfaca + mfcac;
                        m1     = mfacc + mfcaa;
                        m2     = mfaac + mfcca;
                        oMdrho += m0;
                        m1 += m2;
                        oMdrho += m1;
                        m0 = mfbac + mfbca;
                        m1 = mfbaa + mfbcc;
                        m0 += m1;
                        m1 = mfabc + mfcba;
                        m2 = mfaba + mfcbc;
                        m1 += m2;
                        m0 += m1;
                        m1 = mfacb + mfcab;
                        m2 = mfaab + mfccb;
                        m1 += m2;
                        m0 += m1;
                        oMdrho += m0;
                        m0 = mfabb + mfcbb;
                        m1 = mfbab + mfbcb;
                        m2 = mfbba + mfbbc;
                        m0 += m1 + m2;
                        m0 += mfbbb; // hat gefehlt
                        oMdrho = 1. - (oMdrho + m0);
                        // oMdrho = rho - (oMdrho + m0);

                        ////////////////////////////////////////////////////////////////////////////////////
                        LBMReal wadjust;
                        LBMReal qudricLimit = 0.01;
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Hin
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Z - Dir
                        m2    = mfaaa + mfaac;
                        m1    = mfaac - mfaaa;
                        m0    = m2 + mfaab;
                        mfaaa = m0;
                        m0 += c1o36 * oMdrho;
                        mfaab = m1 - m0 * uz;
                        mfaac = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaba + mfabc;
                        m1    = mfabc - mfaba;
                        m0    = m2 + mfabb;
                        mfaba = m0;
                        m0 += c1o9 * oMdrho;
                        mfabb = m1 - m0 * uz;
                        mfabc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaca + mfacc;
                        m1    = mfacc - mfaca;
                        m0    = m2 + mfacb;
                        mfaca = m0;
                        m0 += c1o36 * oMdrho;
                        mfacb = m1 - m0 * uz;
                        mfacc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbaa + mfbac;
                        m1    = mfbac - mfbaa;
                        m0    = m2 + mfbab;
                        mfbaa = m0;
                        m0 += c1o9 * oMdrho;
                        mfbab = m1 - m0 * uz;
                        mfbac = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbba + mfbbc;
                        m1    = mfbbc - mfbba;
                        m0    = m2 + mfbbb;
                        mfbba = m0;
                        m0 += c4o9 * oMdrho;
                        mfbbb = m1 - m0 * uz;
                        mfbbc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbca + mfbcc;
                        m1    = mfbcc - mfbca;
                        m0    = m2 + mfbcb;
                        mfbca = m0;
                        m0 += c1o9 * oMdrho;
                        mfbcb = m1 - m0 * uz;
                        mfbcc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcaa + mfcac;
                        m1    = mfcac - mfcaa;
                        m0    = m2 + mfcab;
                        mfcaa = m0;
                        m0 += c1o36 * oMdrho;
                        mfcab = m1 - m0 * uz;
                        mfcac = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcba + mfcbc;
                        m1    = mfcbc - mfcba;
                        m0    = m2 + mfcbb;
                        mfcba = m0;
                        m0 += c1o9 * oMdrho;
                        mfcbb = m1 - m0 * uz;
                        mfcbc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcca + mfccc;
                        m1    = mfccc - mfcca;
                        m0    = m2 + mfccb;
                        mfcca = m0;
                        m0 += c1o36 * oMdrho;
                        mfccb = m1 - m0 * uz;
                        mfccc = m2 - 2. * m1 * uz + uz2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Y - Dir
                        m2    = mfaaa + mfaca;
                        m1    = mfaca - mfaaa;
                        m0    = m2 + mfaba;
                        mfaaa = m0;
                        m0 += c1o6 * oMdrho;
                        mfaba = m1 - m0 * uy;
                        mfaca = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaab + mfacb;
                        m1    = mfacb - mfaab;
                        m0    = m2 + mfabb;
                        mfaab = m0;
                        mfabb = m1 - m0 * uy;
                        mfacb = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaac + mfacc;
                        m1    = mfacc - mfaac;
                        m0    = m2 + mfabc;
                        mfaac = m0;
                        m0 += c1o18 * oMdrho;
                        mfabc = m1 - m0 * uy;
                        mfacc = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbaa + mfbca;
                        m1    = mfbca - mfbaa;
                        m0    = m2 + mfbba;
                        mfbaa = m0;
                        m0 += c2o3 * oMdrho;
                        mfbba = m1 - m0 * uy;
                        mfbca = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbab + mfbcb;
                        m1    = mfbcb - mfbab;
                        m0    = m2 + mfbbb;
                        mfbab = m0;
                        mfbbb = m1 - m0 * uy;
                        mfbcb = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfbac + mfbcc;
                        m1    = mfbcc - mfbac;
                        m0    = m2 + mfbbc;
                        mfbac = m0;
                        m0 += c2o9 * oMdrho;
                        mfbbc = m1 - m0 * uy;
                        mfbcc = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcaa + mfcca;
                        m1    = mfcca - mfcaa;
                        m0    = m2 + mfcba;
                        mfcaa = m0;
                        m0 += c1o6 * oMdrho;
                        mfcba = m1 - m0 * uy;
                        mfcca = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcab + mfccb;
                        m1    = mfccb - mfcab;
                        m0    = m2 + mfcbb;
                        mfcab = m0;
                        mfcbb = m1 - m0 * uy;
                        mfccb = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfcac + mfccc;
                        m1    = mfccc - mfcac;
                        m0    = m2 + mfcbc;
                        mfcac = m0;
                        m0 += c1o18 * oMdrho;
                        mfcbc = m1 - m0 * uy;
                        mfccc = m2 - 2. * m1 * uy + uy2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // X - Dir
                        m2    = mfaaa + mfcaa;
                        m1    = mfcaa - mfaaa;
                        m0    = m2 + mfbaa;
                        mfaaa = m0;
                        m0 += 1. * oMdrho;
                        mfbaa = m1 - m0 * ux;
                        mfcaa = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaba + mfcba;
                        m1    = mfcba - mfaba;
                        m0    = m2 + mfbba;
                        mfaba = m0;
                        mfbba = m1 - m0 * ux;
                        mfcba = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaca + mfcca;
                        m1    = mfcca - mfaca;
                        m0    = m2 + mfbca;
                        mfaca = m0;
                        m0 += c1o3 * oMdrho;
                        mfbca = m1 - m0 * ux;
                        mfcca = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaab + mfcab;
                        m1    = mfcab - mfaab;
                        m0    = m2 + mfbab;
                        mfaab = m0;
                        mfbab = m1 - m0 * ux;
                        mfcab = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfabb + mfcbb;
                        m1    = mfcbb - mfabb;
                        m0    = m2 + mfbbb;
                        mfabb = m0;
                        mfbbb = m1 - m0 * ux;
                        mfcbb = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfacb + mfccb;
                        m1    = mfccb - mfacb;
                        m0    = m2 + mfbcb;
                        mfacb = m0;
                        mfbcb = m1 - m0 * ux;
                        mfccb = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfaac + mfcac;
                        m1    = mfcac - mfaac;
                        m0    = m2 + mfbac;
                        mfaac = m0;
                        m0 += c1o3 * oMdrho;
                        mfbac = m1 - m0 * ux;
                        mfcac = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfabc + mfcbc;
                        m1    = mfcbc - mfabc;
                        m0    = m2 + mfbbc;
                        mfabc = m0;
                        mfbbc = m1 - m0 * ux;
                        mfcbc = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m2    = mfacc + mfccc;
                        m1    = mfccc - mfacc;
                        m0    = m2 + mfbcc;
                        mfacc = m0;
                        m0 += c1o9 * oMdrho;
                        mfbcc = m1 - m0 * ux;
                        mfccc = m2 - 2. * m1 * ux + ux2 * m0;
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Cumulants
                        ////////////////////////////////////////////////////////////////////////////////////
                        LBMReal OxxPyyPzz = 1.; // omega2 or bulk viscosity
                        LBMReal OxyyPxzz  = 1.; //-s9;//2+s9;//
                        LBMReal OxyyMxzz  = 1.; // 2+s9;//
                        LBMReal O4        = 1.;
                        LBMReal O5        = 1.;
                        LBMReal O6        = 1.;

                        // Cum 4.
                        LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
                        LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
                        LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

                        LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) +
                                                  c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);
                        LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) +
                                                  c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);
                        LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) +
                                                  c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);

                        // Cum 5.
                        LBMReal CUMbcc = mfbcc -
                                         (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb +
                                          2. * (mfbab * mfacb + mfbba * mfabc)) -
                                         c1o3 * (mfbca + mfbac) * oMdrho;
                        LBMReal CUMcbc = mfcbc -
                                         (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb +
                                          2. * (mfabb * mfcab + mfbba * mfbac)) -
                                         c1o3 * (mfcba + mfabc) * oMdrho;
                        LBMReal CUMccb = mfccb -
                                         (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb +
                                          2. * (mfbab * mfbca + mfabb * mfcba)) -
                                         c1o3 * (mfacb + mfcab) * oMdrho;

                        // Cum 6.
                        LBMReal CUMccc =
                            mfccc +
                            ((-4. * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                              4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
                              2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) +
                             (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                              2. * (mfcaa * mfaca * mfaac) + 16. * mfbba * mfbab * mfabb) -
                             c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho -
                             c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) -
                             c1o27 * oMdrho * oMdrho * (-2. * oMdrho) +
                             (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                              (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) *
                                 c2o3 * oMdrho) +
                            c1o27 * oMdrho;

                        // 2.
                        // linear combinations
                        LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
                        LBMReal mxxMyy    = mfcaa - mfaca;
                        LBMReal mxxMzz    = mfcaa - mfaac;

                        LBMReal dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
                        LBMReal dyuy = dxux + collFactorM * c3o2 * mxxMyy;
                        LBMReal dzuz = dxux + collFactorM * c3o2 * mxxMzz;

                        (*divU)(x1, x2, x3) = dxux + dyuy + dzuz;

                        // relax
                        mxxPyyPzz += OxxPyyPzz * (mfaaa - mxxPyyPzz) -
                                     3. * (1. - c1o2 * OxxPyyPzz) * (ux2 * dxux + uy2 * dyuy + uz2 * dzuz);
                        mxxMyy += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uy2 * dyuy);
                        mxxMzz += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uz2 * dzuz);

                        mfabb += collFactorM * (-mfabb);
                        mfbab += collFactorM * (-mfbab);
                        mfbba += collFactorM * (-mfbba);

                        // linear combinations back
                        mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
                        mfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
                        mfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);

                        // 3.
                        // linear combinations
                        LBMReal mxxyPyzz = mfcba + mfabc;
                        LBMReal mxxyMyzz = mfcba - mfabc;

                        LBMReal mxxzPyyz = mfcab + mfacb;
                        LBMReal mxxzMyyz = mfcab - mfacb;

                        LBMReal mxyyPxzz = mfbca + mfbac;
                        LBMReal mxyyMxzz = mfbca - mfbac;

                        // relax
                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mfbbb) / (fabs(mfbbb) + qudricLimit);
                        mfbbb += wadjust * (-mfbbb);
                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxyPyzz) / (fabs(mxxyPyzz) + qudricLimit);
                        mxxyPyzz += wadjust * (-mxxyPyzz);
                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxyMyzz) / (fabs(mxxyMyzz) + qudricLimit);
                        mxxyMyzz += wadjust * (-mxxyMyzz);
                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxzPyyz) / (fabs(mxxzPyyz) + qudricLimit);
                        mxxzPyyz += wadjust * (-mxxzPyyz);
                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxzMyyz) / (fabs(mxxzMyyz) + qudricLimit);
                        mxxzMyyz += wadjust * (-mxxzMyyz);
                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxyyPxzz) / (fabs(mxyyPxzz) + qudricLimit);
                        mxyyPxzz += wadjust * (-mxyyPxzz);
                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxyyMxzz) / (fabs(mxyyMxzz) + qudricLimit);
                        mxyyMxzz += wadjust * (-mxyyMxzz);

                        // linear combinations back
                        mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
                        mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
                        mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
                        mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
                        mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
                        mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

                        // 4.
                        CUMacc += O4 * (-CUMacc);
                        CUMcac += O4 * (-CUMcac);
                        CUMcca += O4 * (-CUMcca);

                        CUMbbc += O4 * (-CUMbbc);
                        CUMbcb += O4 * (-CUMbcb);
                        CUMcbb += O4 * (-CUMcbb);

                        // 5.
                        CUMbcc += O5 * (-CUMbcc);
                        CUMcbc += O5 * (-CUMcbc);
                        CUMccb += O5 * (-CUMccb);

                        // 6.
                        CUMccc += O6 * (-CUMccc);

                        // back cumulants to central moments
                        // 4.
                        mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
                        mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
                        mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

                        mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho +
                                c1o9 * (oMdrho - 1) * oMdrho;
                        mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho +
                                c1o9 * (oMdrho - 1) * oMdrho;
                        mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho +
                                c1o9 * (oMdrho - 1) * oMdrho;

                        // 5.
                        mfbcc = CUMbcc +
                                (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb +
                                 2. * (mfbab * mfacb + mfbba * mfabc)) +
                                c1o3 * (mfbca + mfbac) * oMdrho;
                        mfcbc = CUMcbc +
                                (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb +
                                 2. * (mfabb * mfcab + mfbba * mfbac)) +
                                c1o3 * (mfcba + mfabc) * oMdrho;
                        mfccb = CUMccb +
                                (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb +
                                 2. * (mfbab * mfbca + mfabb * mfcba)) +
                                c1o3 * (mfacb + mfcab) * oMdrho;

                        // 6.
                        mfccc = CUMccc -
                                ((-4. * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
                                  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) -
                                  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) +
                                 (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
                                  2. * (mfcaa * mfaca * mfaac) + 16. * mfbba * mfbab * mfabb) -
                                 c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho -
                                 c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) -
                                 c1o27 * oMdrho * oMdrho * (-2. * oMdrho) +
                                 (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
                                  (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) *
                                     c2o3 * oMdrho) -
                                c1o27 * oMdrho;

                        ////////////////////////////////////////////////////////////////////////////////////
                        // forcing
                        mfbaa = -mfbaa;
                        mfaba = -mfaba;
                        mfaab = -mfaab;
                        //////////////////////////////////////////////////////////////////////////////////////

                        ////////////////////////////////////////////////////////////////////////////////////
                        // back
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Z - Dir
                        m0    = mfaac * c1o2 + mfaab * (uz - c1o2) + (mfaaa + 1. * oMdrho) * (uz2 - uz) * c1o2;
                        m1    = -mfaac - 2. * mfaab * uz + mfaaa * (1. - uz2) - 1. * oMdrho * uz2;
                        m2    = mfaac * c1o2 + mfaab * (uz + c1o2) + (mfaaa + 1. * oMdrho) * (uz2 + uz) * c1o2;
                        mfaaa = m0;
                        mfaab = m1;
                        mfaac = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfabc * c1o2 + mfabb * (uz - c1o2) + mfaba * (uz2 - uz) * c1o2;
                        m1    = -mfabc - 2. * mfabb * uz + mfaba * (1. - uz2);
                        m2    = mfabc * c1o2 + mfabb * (uz + c1o2) + mfaba * (uz2 + uz) * c1o2;
                        mfaba = m0;
                        mfabb = m1;
                        mfabc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfacc * c1o2 + mfacb * (uz - c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
                        m1    = -mfacc - 2. * mfacb * uz + mfaca * (1. - uz2) - c1o3 * oMdrho * uz2;
                        m2    = mfacc * c1o2 + mfacb * (uz + c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
                        mfaca = m0;
                        mfacb = m1;
                        mfacc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfbac * c1o2 + mfbab * (uz - c1o2) + mfbaa * (uz2 - uz) * c1o2;
                        m1    = -mfbac - 2. * mfbab * uz + mfbaa * (1. - uz2);
                        m2    = mfbac * c1o2 + mfbab * (uz + c1o2) + mfbaa * (uz2 + uz) * c1o2;
                        mfbaa = m0;
                        mfbab = m1;
                        mfbac = m2;
                        /////////b//////////////////////////////////////////////////////////////////////////
                        m0    = mfbbc * c1o2 + mfbbb * (uz - c1o2) + mfbba * (uz2 - uz) * c1o2;
                        m1    = -mfbbc - 2. * mfbbb * uz + mfbba * (1. - uz2);
                        m2    = mfbbc * c1o2 + mfbbb * (uz + c1o2) + mfbba * (uz2 + uz) * c1o2;
                        mfbba = m0;
                        mfbbb = m1;
                        mfbbc = m2;
                        /////////b//////////////////////////////////////////////////////////////////////////
                        m0    = mfbcc * c1o2 + mfbcb * (uz - c1o2) + mfbca * (uz2 - uz) * c1o2;
                        m1    = -mfbcc - 2. * mfbcb * uz + mfbca * (1. - uz2);
                        m2    = mfbcc * c1o2 + mfbcb * (uz + c1o2) + mfbca * (uz2 + uz) * c1o2;
                        mfbca = m0;
                        mfbcb = m1;
                        mfbcc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcac * c1o2 + mfcab * (uz - c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
                        m1    = -mfcac - 2. * mfcab * uz + mfcaa * (1. - uz2) - c1o3 * oMdrho * uz2;
                        m2    = mfcac * c1o2 + mfcab * (uz + c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
                        mfcaa = m0;
                        mfcab = m1;
                        mfcac = m2;
                        /////////c//////////////////////////////////////////////////////////////////////////
                        m0    = mfcbc * c1o2 + mfcbb * (uz - c1o2) + mfcba * (uz2 - uz) * c1o2;
                        m1    = -mfcbc - 2. * mfcbb * uz + mfcba * (1. - uz2);
                        m2    = mfcbc * c1o2 + mfcbb * (uz + c1o2) + mfcba * (uz2 + uz) * c1o2;
                        mfcba = m0;
                        mfcbb = m1;
                        mfcbc = m2;
                        /////////c//////////////////////////////////////////////////////////////////////////
                        m0    = mfccc * c1o2 + mfccb * (uz - c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 - uz) * c1o2;
                        m1    = -mfccc - 2. * mfccb * uz + mfcca * (1. - uz2) - c1o9 * oMdrho * uz2;
                        m2    = mfccc * c1o2 + mfccb * (uz + c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 + uz) * c1o2;
                        mfcca = m0;
                        mfccb = m1;
                        mfccc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // Y - Dir
                        m0    = mfaca * c1o2 + mfaba * (uy - c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfaca - 2. * mfaba * uy + mfaaa * (1. - uy2) - c1o6 * oMdrho * uy2;
                        m2    = mfaca * c1o2 + mfaba * (uy + c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
                        mfaaa = m0;
                        mfaba = m1;
                        mfaca = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfacb * c1o2 + mfabb * (uy - c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfacb - 2. * mfabb * uy + mfaab * (1. - uy2) - c2o3 * oMdrho * uy2;
                        m2    = mfacb * c1o2 + mfabb * (uy + c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 + uy) * c1o2;
                        mfaab = m0;
                        mfabb = m1;
                        mfacb = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfacc * c1o2 + mfabc * (uy - c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfacc - 2. * mfabc * uy + mfaac * (1. - uy2) - c1o6 * oMdrho * uy2;
                        m2    = mfacc * c1o2 + mfabc * (uy + c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
                        mfaac = m0;
                        mfabc = m1;
                        mfacc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfbca * c1o2 + mfbba * (uy - c1o2) + mfbaa * (uy2 - uy) * c1o2;
                        m1    = -mfbca - 2. * mfbba * uy + mfbaa * (1. - uy2);
                        m2    = mfbca * c1o2 + mfbba * (uy + c1o2) + mfbaa * (uy2 + uy) * c1o2;
                        mfbaa = m0;
                        mfbba = m1;
                        mfbca = m2;
                        /////////b//////////////////////////////////////////////////////////////////////////
                        m0    = mfbcb * c1o2 + mfbbb * (uy - c1o2) + mfbab * (uy2 - uy) * c1o2;
                        m1    = -mfbcb - 2. * mfbbb * uy + mfbab * (1. - uy2);
                        m2    = mfbcb * c1o2 + mfbbb * (uy + c1o2) + mfbab * (uy2 + uy) * c1o2;
                        mfbab = m0;
                        mfbbb = m1;
                        mfbcb = m2;
                        /////////b//////////////////////////////////////////////////////////////////////////
                        m0    = mfbcc * c1o2 + mfbbc * (uy - c1o2) + mfbac * (uy2 - uy) * c1o2;
                        m1    = -mfbcc - 2. * mfbbc * uy + mfbac * (1. - uy2);
                        m2    = mfbcc * c1o2 + mfbbc * (uy + c1o2) + mfbac * (uy2 + uy) * c1o2;
                        mfbac = m0;
                        mfbbc = m1;
                        mfbcc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcca * c1o2 + mfcba * (uy - c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfcca - 2. * mfcba * uy + mfcaa * (1. - uy2) - c1o18 * oMdrho * uy2;
                        m2    = mfcca * c1o2 + mfcba * (uy + c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
                        mfcaa = m0;
                        mfcba = m1;
                        mfcca = m2;
                        /////////c//////////////////////////////////////////////////////////////////////////
                        m0    = mfccb * c1o2 + mfcbb * (uy - c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfccb - 2. * mfcbb * uy + mfcab * (1. - uy2) - c2o9 * oMdrho * uy2;
                        m2    = mfccb * c1o2 + mfcbb * (uy + c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 + uy) * c1o2;
                        mfcab = m0;
                        mfcbb = m1;
                        mfccb = m2;
                        /////////c//////////////////////////////////////////////////////////////////////////
                        m0    = mfccc * c1o2 + mfcbc * (uy - c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
                        m1    = -mfccc - 2. * mfcbc * uy + mfcac * (1. - uy2) - c1o18 * oMdrho * uy2;
                        m2    = mfccc * c1o2 + mfcbc * (uy + c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
                        mfcac = m0;
                        mfcbc = m1;
                        mfccc = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
                        ////////////////////////////////////////////////////////////////////////////////////
                        // X - Dir
                        m0    = mfcaa * c1o2 + mfbaa * (ux - c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcaa - 2. * mfbaa * ux + mfaaa * (1. - ux2) - c1o36 * oMdrho * ux2;
                        m2    = mfcaa * c1o2 + mfbaa * (ux + c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
                        mfaaa = m0;
                        mfbaa = m1;
                        mfcaa = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcba * c1o2 + mfbba * (ux - c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcba - 2. * mfbba * ux + mfaba * (1. - ux2) - c1o9 * oMdrho * ux2;
                        m2    = mfcba * c1o2 + mfbba * (ux + c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
                        mfaba = m0;
                        mfbba = m1;
                        mfcba = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcca * c1o2 + mfbca * (ux - c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcca - 2. * mfbca * ux + mfaca * (1. - ux2) - c1o36 * oMdrho * ux2;
                        m2    = mfcca * c1o2 + mfbca * (ux + c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
                        mfaca = m0;
                        mfbca = m1;
                        mfcca = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcab * c1o2 + mfbab * (ux - c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcab - 2. * mfbab * ux + mfaab * (1. - ux2) - c1o9 * oMdrho * ux2;
                        m2    = mfcab * c1o2 + mfbab * (ux + c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
                        mfaab = m0;
                        mfbab = m1;
                        mfcab = m2;
                        ///////////b////////////////////////////////////////////////////////////////////////
                        m0    = mfcbb * c1o2 + mfbbb * (ux - c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcbb - 2. * mfbbb * ux + mfabb * (1. - ux2) - c4o9 * oMdrho * ux2;
                        m2    = mfcbb * c1o2 + mfbbb * (ux + c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 + ux) * c1o2;
                        mfabb = m0;
                        mfbbb = m1;
                        mfcbb = m2;
                        ///////////b////////////////////////////////////////////////////////////////////////
                        m0    = mfccb * c1o2 + mfbcb * (ux - c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfccb - 2. * mfbcb * ux + mfacb * (1. - ux2) - c1o9 * oMdrho * ux2;
                        m2    = mfccb * c1o2 + mfbcb * (ux + c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
                        mfacb = m0;
                        mfbcb = m1;
                        mfccb = m2;
                        ////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////////
                        m0    = mfcac * c1o2 + mfbac * (ux - c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcac - 2. * mfbac * ux + mfaac * (1. - ux2) - c1o36 * oMdrho * ux2;
                        m2    = mfcac * c1o2 + mfbac * (ux + c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
                        mfaac = m0;
                        mfbac = m1;
                        mfcac = m2;
                        ///////////c////////////////////////////////////////////////////////////////////////
                        m0    = mfcbc * c1o2 + mfbbc * (ux - c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfcbc - 2. * mfbbc * ux + mfabc * (1. - ux2) - c1o9 * oMdrho * ux2;
                        m2    = mfcbc * c1o2 + mfbbc * (ux + c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
                        mfabc = m0;
                        mfbbc = m1;
                        mfcbc = m2;
                        ///////////c////////////////////////////////////////////////////////////////////////
                        m0    = mfccc * c1o2 + mfbcc * (ux - c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
                        m1    = -mfccc - 2. * mfbcc * ux + mfacc * (1. - ux2) - c1o36 * oMdrho * ux2;
                        m2    = mfccc * c1o2 + mfbcc * (ux + c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
                        mfacc = m0;
                        mfbcc = m1;
                        mfccc = m2;

                        ///////////////////////////////////////////////////////////////////////////

                        //////////////////////////////////////////////////////////////////////////
                        // proof correctness
                        //////////////////////////////////////////////////////////////////////////
#ifdef PROOF_CORRECTNESS
                        LBMReal rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
                                           (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
                                           (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) +
                                           (mfbba + mfbbc) + mfbbb;

                        LBMReal dif = rho1 - rho_post;
#ifdef SINGLEPRECISION
                        if (dif > 10.0E-7 || dif < -10.0E-7)
#else
                        if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
                        {
                            UB_THROW(UbException(UB_EXARGS,
                                                 "rho=" + UbSystem::toString(rho) + ", rho_post=" +
                                                     UbSystem::toString(rho_post) + " dif=" + UbSystem::toString(dif) +
                                                     " rho is not correct for node " + UbSystem::toString(x1) + "," +
                                                     UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
                        }
#endif

                        mfcbb = rho * c1o3 * (mfcbb) + 0.5 * forcingTerm[DIR_P00];
                        mfbcb = rho * c1o3 * (mfbcb) + 0.5 * forcingTerm[DIR_0P0];
                        mfbbc = rho * c1o3 * (mfbbc) + 0.5 * forcingTerm[DIR_00P];
                        mfccb = rho * c1o3 * (mfccb) + 0.5 * forcingTerm[DIR_PP0];
                        mfacb = rho * c1o3 * (mfacb) + 0.5 * forcingTerm[DIR_MP0];
                        mfcbc = rho * c1o3 * (mfcbc) + 0.5 * forcingTerm[DIR_P0P];
                        mfabc = rho * c1o3 * (mfabc) + 0.5 * forcingTerm[DIR_M0P];
                        mfbcc = rho * c1o3 * (mfbcc) + 0.5 * forcingTerm[DIR_0PP];
                        mfbac = rho * c1o3 * (mfbac) + 0.5 * forcingTerm[DIR_0MP];
                        mfccc = rho * c1o3 * (mfccc) + 0.5 * forcingTerm[DIR_PPP];
                        mfacc = rho * c1o3 * (mfacc) + 0.5 * forcingTerm[DIR_MPP];
                        mfcac = rho * c1o3 * (mfcac) + 0.5 * forcingTerm[DIR_PMP];
                        mfaac = rho * c1o3 * (mfaac) + 0.5 * forcingTerm[DIR_MMP];
                        mfabb = rho * c1o3 * (mfabb) + 0.5 * forcingTerm[DIR_M00];
                        mfbab = rho * c1o3 * (mfbab) + 0.5 * forcingTerm[DIR_0M0];
                        mfbba = rho * c1o3 * (mfbba) + 0.5 * forcingTerm[DIR_00M];
                        mfaab = rho * c1o3 * (mfaab) + 0.5 * forcingTerm[DIR_MM0];
                        mfcab = rho * c1o3 * (mfcab) + 0.5 * forcingTerm[DIR_PM0];
                        mfaba = rho * c1o3 * (mfaba) + 0.5 * forcingTerm[DIR_M0M];
                        mfcba = rho * c1o3 * (mfcba) + 0.5 * forcingTerm[DIR_P0M];
                        mfbaa = rho * c1o3 * (mfbaa) + 0.5 * forcingTerm[DIR_0MM];
                        mfbca = rho * c1o3 * (mfbca) + 0.5 * forcingTerm[DIR_0PM];
                        mfaaa = rho * c1o3 * (mfaaa) + 0.5 * forcingTerm[DIR_MMM];
                        mfcaa = rho * c1o3 * (mfcaa) + 0.5 * forcingTerm[DIR_PMM];
                        mfaca = rho * c1o3 * (mfaca) + 0.5 * forcingTerm[DIR_MPM];
                        mfcca = rho * c1o3 * (mfcca) + 0.5 * forcingTerm[DIR_PPM];
                        mfbbb = rho * c1o3 * (mfbbb) + 0.5 * forcingTerm[DIR_000];

                        //////////////////////////////////////////////////////////////////////////
                        // write distribution for F
                        //////////////////////////////////////////////////////////////////////////

                        (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3)     = mfabb;
                        (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3)     = mfbab;
                        (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)     = mfbba;
                        (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3)    = mfaab;
                        (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3)   = mfcab;
                        (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3)    = mfaba;
                        (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3)   = mfcba;
                        (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3)    = mfbaa;
                        (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3)   = mfbca;
                        (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3)   = mfaaa;
                        (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3)  = mfcaa;
                        (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3)  = mfaca;
                        (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

                        (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)     = mfcbb;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3)     = mfbcb;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p)     = mfbbc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfccb;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3)    = mfacb;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfcbc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p)    = mfabc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfbcc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p)    = mfbac;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfacc;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfcac;
                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfaac;

                        (*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;

                        /////////////////////  P H A S E - F I E L D   S O L V E R
                        ////////////////////////////////////////////

                        /////////////////////   PHASE-FIELD BGK SOLVER ///////////////////////////////

                        h[DIR_P00]   = (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3);
                        h[DIR_0P0]   = (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3);
                        h[DIR_00P]   = (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3);
                        h[DIR_PP0]  = (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3);
                        h[DIR_MP0]  = (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3);
                        h[DIR_P0P]  = (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3);
                        h[DIR_M0P]  = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3);
                        h[DIR_0PP]  = (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3);
                        h[DIR_0MP]  = (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3);
                        h[DIR_PPP] = (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3);
                        h[DIR_MPP] = (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3);
                        h[DIR_PMP] = (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3);
                        h[DIR_MMP] = (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3);

                        h[DIR_M00]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3);
                        h[DIR_0M0]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3);
                        h[DIR_00M]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p);
                        h[DIR_MM0]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3);
                        h[DIR_PM0]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3);
                        h[DIR_M0M]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p);
                        h[DIR_P0M]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p);
                        h[DIR_0MM]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p);
                        h[DIR_0PM]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p);
                        h[DIR_MMM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                        h[DIR_PMM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p);
                        h[DIR_MPM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p);
                        h[DIR_PPM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p);

                        h[DIR_000] = (*this->zeroDistributionsH)(x1, x2, x3);

                        for (int dir = STARTF; dir < (ENDF + 1); dir++) {
                            LBMReal velProd = DX1[dir] * ux + DX2[dir] * uy + DX3[dir] * uz;
                            LBMReal velSq1  = velProd * velProd;
                            LBMReal hEq; //, gEq;

                            if (dir != DIR_000) {
                                LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]]) / 2.0;
                                LBMReal hSource     = (tauH - 0.5) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * (dirGrad_phi) / denom; 
                                hEq = phi[DIR_000] * WEIGTH[dir] * (1.0 + 3.0 * velProd + 4.5 * velSq1 - 1.5 * (ux2 + uy2 + uz2)) +                                 hSource * WEIGTH[dir];

                                // This corresponds with the collision factor of 1.0 which equals (tauH + 0.5).
                                h[dir] = h[dir] - (h[dir] - hEq) / (tauH); 

                            } else {
                                hEq = phi[DIR_000] * WEIGTH[DIR_000] * (1.0 - 1.5 * (ux2 + uy2 + uz2));
                                h[DIR_000] = h[DIR_000] - (h[DIR_000] - hEq) / (tauH); 
                            }
                        }

                        (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3)     = h[D3Q27System::INV_P00];
                        (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3)     = h[D3Q27System::INV_0P0];
                        (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3)     = h[D3Q27System::INV_00P];
                        (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3)    = h[D3Q27System::INV_PP0];
                        (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3)   = h[D3Q27System::INV_MP0];
                        (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3)    = h[D3Q27System::INV_P0P];
                        (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3)   = h[D3Q27System::INV_M0P];
                        (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3)    = h[D3Q27System::INV_0PP];
                        (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3)   = h[D3Q27System::INV_0MP];
                        (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3)   = h[D3Q27System::INV_PPP];
                        (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3)  = h[D3Q27System::INV_MPP];
                        (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3)  = h[D3Q27System::INV_PMP];
                        (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = h[D3Q27System::INV_MMP];

                        (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3)     = h[D3Q27System::INV_M00];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3)     = h[D3Q27System::INV_0M0];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p)     = h[D3Q27System::INV_00M];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3)   = h[D3Q27System::INV_MM0];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3)    = h[D3Q27System::INV_PM0];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p)   = h[D3Q27System::INV_M0M];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p)    = h[D3Q27System::INV_P0M];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p)   = h[D3Q27System::INV_0MM];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p)    = h[D3Q27System::INV_0PM];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = h[D3Q27System::INV_MMM];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p)  = h[D3Q27System::INV_PMM];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p)  = h[D3Q27System::INV_MPM];
                        (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p)   = h[D3Q27System::INV_PPM];

                        (*this->zeroDistributionsH)(x1, x2, x3) = h[D3Q27System::DIR_000];

                        /////////////////////   END OF OLD BGK SOLVER ///////////////////////////////
                    }
                }
            }
        }
        dataSet->setPhaseField(divU);
    }

//////////////////////////////////////////////////////////////////////////

LBMReal MultiphaseCumulantLBMKernel::gradX1_phi()
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX1[k] * phi[k];
    }
    return 3.0 * sum;
}

LBMReal MultiphaseCumulantLBMKernel::gradX2_phi()
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX2[k] * phi[k];
    }
    return 3.0 * sum;
}

LBMReal MultiphaseCumulantLBMKernel::gradX3_phi()
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * DX3[k] * phi[k];
    }
    return 3.0 * sum;
}

LBMReal MultiphaseCumulantLBMKernel::nabla2_phi()
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
    for (int k = FSTARTDIR; k <= FENDDIR; k++) {
        sum += WEIGTH[k] * (phi[k] - phi[DIR_000]);
    }
    return 6.0 * sum;
}

void MultiphaseCumulantLBMKernel::computePhasefield()
{
    using namespace D3Q27System;
    SPtr<DistributionArray3D> distributionsH = dataSet->getHdistributions();

    int minX1 = ghostLayerWidth;
    int minX2 = ghostLayerWidth;
    int minX3 = ghostLayerWidth;
    int maxX1 = (int)distributionsH->getNX1() - ghostLayerWidth;
    int maxX2 = (int)distributionsH->getNX2() - ghostLayerWidth;
    int maxX3 = (int)distributionsH->getNX3() - ghostLayerWidth;

    //------------- Computing the phase-field ------------------
    for (int x3 = minX3; x3 < maxX3; x3++) {
        for (int x2 = minX2; x2 < maxX2; x2++) {
            for (int x1 = minX1; x1 < maxX1; x1++) {
                // if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
                {
                    int x1p = x1 + 1;
                    int x2p = x2 + 1;
                    int x3p = x3 + 1;

                    h[DIR_P00]   = (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3);
                    h[DIR_0P0]   = (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3);
                    h[DIR_00P]   = (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3);
                    h[DIR_PP0]  = (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3);
                    h[DIR_MP0]  = (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3);
                    h[DIR_P0P]  = (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3);
                    h[DIR_M0P]  = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3);
                    h[DIR_0PP]  = (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3);
                    h[DIR_0MP]  = (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3);
                    h[DIR_PPP] = (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3);
                    h[DIR_MPP] = (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3);
                    h[DIR_PMP] = (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3);
                    h[DIR_MMP] = (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3);

                    h[DIR_M00]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3);
                    h[DIR_0M0]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3);
                    h[DIR_00M]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p);
                    h[DIR_MM0]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3);
                    h[DIR_PM0]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3);
                    h[DIR_M0M]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p);
                    h[DIR_P0M]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p);
                    h[DIR_0MM]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p);
                    h[DIR_0PM]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p);
                    h[DIR_MMM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p);
                    h[DIR_PMM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p);
                    h[DIR_MPM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p);
                    h[DIR_PPM] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p);

                    h[DIR_000] = (*this->zeroDistributionsH)(x1, x2, x3);
                }
            }
        }
    }
}

void MultiphaseCumulantLBMKernel::findNeighbors(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
                                                int x3)
{
    using namespace D3Q27System;

    SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();

    phi[DIR_000] = (*ph)(x1, x2, x3);

    for (int k = FSTARTDIR; k <= FENDDIR; k++) {

        if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
            phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
        } else {

            phi[k] = 0.;//16.03.2021 quick fix for uninitialized variables, might influence contact angle!
         }
    }
}

void MultiphaseCumulantLBMKernel::swapDistributions()
{
    LBMKernel::swapDistributions();
    dataSet->getHdistributions()->swap();
}