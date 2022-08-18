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
//! \file MultiphaseScratchCumulantLBMKernel.cpp
//! \ingroup LBMKernel
//! \author Hesameddin Safari
//=======================================================================================

#include "MultiphaseScratchCumulantLBMKernel.h"
#include "BCArray3D.h"
#include "Block3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
#include <cmath>
#include <iostream>

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
MultiphaseScratchCumulantLBMKernel::MultiphaseScratchCumulantLBMKernel() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
void MultiphaseScratchCumulantLBMKernel::initDataSet()
{
    SPtr<DistributionArray3D> f(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    SPtr<DistributionArray3D> h(new D3Q27EsoTwist3DSplittedVector(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9)); // For phase-field
    SPtr<PhaseFieldArray3D> divU(new PhaseFieldArray3D(nx[0] + 2, nx[1] + 2, nx[2] + 2, 0.0));
    dataSet->setFdistributions(f);
    dataSet->setHdistributions(h); // For phase-field
    dataSet->setPhaseField(divU);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> MultiphaseScratchCumulantLBMKernel::clone()
{
    SPtr<LBMKernel> kernel(new MultiphaseScratchCumulantLBMKernel());
    kernel->setNX(nx);
    dynamicPointerCast<MultiphaseScratchCumulantLBMKernel>(kernel)->initDataSet();
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
 void  MultiphaseScratchCumulantLBMKernel::forwardInverseChimeraWithKincompressible(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2, LBMReal Kinverse, LBMReal K, LBMReal oneMinusRho) {
	using namespace UbMath;
    LBMReal m2 = mfa + mfc;
	LBMReal m1 = mfc - mfa;
	LBMReal m0 = m2 + mfb;
	mfa = m0;
	m0 *= Kinverse;
	m0 += oneMinusRho;
	mfb = (m1 * Kinverse - m0 * vv) * K;
	mfc = ((m2 - c2 * m1 * vv) * Kinverse + v2 * m0) * K;
}

////////////////////////////////////////////////////////////////////////////////
 void  MultiphaseScratchCumulantLBMKernel::backwardInverseChimeraWithKincompressible(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2, LBMReal Kinverse, LBMReal K, LBMReal oneMinusRho) {
	using namespace UbMath;
    LBMReal m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 - vv) * c1o2) * K;
	LBMReal m1 = (((mfa - mfc) - c2 * mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (-v2)) * K;
	mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 + vv) * c1o2) * K;
	mfa = m0;
	mfb = m1;
}


////////////////////////////////////////////////////////////////////////////////
 void  MultiphaseScratchCumulantLBMKernel::forwardChimera(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2) {
	using namespace UbMath;
    LBMReal m1 = (mfa + mfc) + mfb;
	LBMReal m2 = mfc - mfa;
	mfc = (mfc + mfa) + (v2 * m1 - c2 * vv * m2);
	mfb = m2 - vv * m1;
	mfa = m1;
}


 void  MultiphaseScratchCumulantLBMKernel::backwardChimera(LBMReal& mfa, LBMReal& mfb, LBMReal& mfc, LBMReal vv, LBMReal v2) {
	using namespace UbMath;
    LBMReal ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
	LBMReal mb = ((mfa - mfc) - mfa * v2) - c2 * mfb * vv;
	mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
	mfb = mb;
	mfa = ma;
}


void MultiphaseScratchCumulantLBMKernel::calculate(int step)
{
    using namespace D3Q27System;
    using namespace UbMath;

    forcingX1 = 0.0;
    forcingX2 = 0.0;
    forcingX3 = 0.0;

	LBMReal oneOverInterfaceScale = 1.0;// 1.0 / 3.0;
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


		/////For velocity filter

		//CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr velocityX(
		//	new CbArray3D<LBMReal, IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, 0.0));
		//CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr velocityY(
		//	new CbArray3D<LBMReal, IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, 0.0));
		//CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr velocityZ(
		//	new CbArray3D<LBMReal, IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, 0.0));


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
                        (*phaseField)(x1, x2, x3) = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca)  + (mfcaa + mfacc))  ) +
                                                    (((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
                                                    ((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
                                                    (mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;
						//(*phaseField)(x1, x2, x3) = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
						//	(mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
						//	(mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) +
						//	(mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

						///Velocity filter


						LBMReal rhoH = 1.0;
						LBMReal rhoL = 1.0 / densityRatio;

						LBMReal rhoToPhi = (rhoH - rhoL) / (phiH - phiL);


						LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH);

						mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) / rho * c3;
						mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) / rho * c3;
						mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) / rho * c3;
						mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) / rho * c3;
						mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) / rho * c3;
						mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) / rho * c3;
						mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) / rho * c3;
						mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) / rho * c3;
						mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) / rho * c3;
						mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) / rho * c3;
						mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) / rho * c3;
						mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) / rho * c3;
						mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) / rho * c3;

						mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) / rho * c3;
						mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) / rho * c3;
						mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) / rho * c3;
						mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) / rho * c3;
						mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) / rho * c3;
						mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) / rho * c3;
						mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) / rho * c3;
						mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) / rho * c3;
						mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) / rho * c3;
						mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) / rho * c3;
						mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) / rho * c3;
						mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) / rho * c3;
						mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) / rho * c3;

						mfbbb = (*this->zeroDistributionsF)(x1, x2, x3) / rho * c3;

						//(*velocityX)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
						//	(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
						//	(mfcbb - mfabb)) ;
						//(*velocityY)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
						//	(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
						//	(mfbcb - mfbab)) ;
						//(*velocityZ)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
						//	(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
						//	(mfbbc - mfbba)) ;




                    }
					else { (*phaseField)(x1, x2, x3) = 0; }
                }
            }
        }

        LBMReal collFactorM;
        //LBMReal forcingTerm[D3Q27System::ENDF + 1];

        for (int x3 = minX3; x3 < maxX3; x3++) {
            for (int x2 = minX2; x2 < maxX2; x2++) {
                for (int x1 = minX1; x1 < maxX1; x1++) {

					//for (int x3 = minX3+1; x3 < maxX3-1; x3++) {
					//	for (int x2 = minX2+1; x2 < maxX2-1; x2++) {
					//		for (int x1 = minX1+1; x1 < maxX1-1; x1++) {
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
						//// reading distributions here appears to be unnecessary!
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

						//LBMReal dX1_phi = 3.0*((
						//	WEIGTH[TNE]*((((*phaseField)(x1 + 1, x2+1, x3+1)- (*phaseField)(x1 - 1, x2 - 1, x3 - 1))+ ((*phaseField)(x1 + 1, x2 - 1, x3 + 1) - (*phaseField)(x1 - 1, x2 + 1, x3 - 1)))
						//	+ (((*phaseField)(x1 + 1, x2 - 1, x3 - 1) - (*phaseField)(x1 - 1, x2 + 1, x3 + 1)) + ((*phaseField)(x1 + 1, x2 + 1, x3 - 1) - (*phaseField)(x1 - 1, x2 - 1, x3 + 1))))
						//	+WEIGTH[NE]* ((((*phaseField)(x1 + 1, x2 + 1, x3) - (*phaseField)(x1 - 1, x2 - 1, x3)) + ((*phaseField)(x1 + 1, x2 - 1, x3) - (*phaseField)(x1 - 1, x2 + 1, x3 )))
						//	+ (((*phaseField)(x1 + 1, x2, x3 - 1) - (*phaseField)(x1 - 1, x2, x3 + 1)) + ((*phaseField)(x1 + 1, x2, x3 + 1) - (*phaseField)(x1 - 1, x2, x3 - 1)))))
						//	+WEIGTH[N]*((*phaseField)(x1 + 1, x2, x3 ) - (*phaseField)(x1 - 1, x2, x3))
						//	); 
						////if (dX1_phi != NdX1_phi) {std::cout<<dX1_phi<<" "<< NdX1_phi<<std::endl;}

						//LBMReal dX2_phi = 3.0 * ((
						//	WEIGTH[TNE] * ((((*phaseField)(x1 + 1, x2 + 1, x3 + 1) - (*phaseField)(x1 - 1, x2 - 1, x3 - 1)) + ((*phaseField)(x1 -1, x2 + 1, x3 + 1) - (*phaseField)(x1 + 1, x2 - 1, x3 - 1)))
						//	+ (((*phaseField)(x1 - 1, x2 + 1, x3 - 1) - (*phaseField)(x1 + 1, x2 - 1, x3 + 1)) + ((*phaseField)(x1 + 1, x2 + 1, x3 - 1) - (*phaseField)(x1 - 1, x2 - 1, x3 + 1))))
						//	+ WEIGTH[NE] * ((((*phaseField)(x1 + 1, x2 + 1, x3) - (*phaseField)(x1 - 1, x2 - 1, x3)) + ((*phaseField)(x1 - 1, x2 + 1, x3) - (*phaseField)(x1 + 1, x2 - 1, x3)))
						//		+ (((*phaseField)(x1, x2+1, x3 - 1) - (*phaseField)(x1 , x2-1, x3 + 1)) + ((*phaseField)(x1 , x2+1, x3 + 1) - (*phaseField)(x1 , x2-1, x3 - 1)))))
						//	+ WEIGTH[N] * ((*phaseField)(x1 , x2+1, x3) - (*phaseField)(x1 , x2-1, x3))
						//	);

						//LBMReal dX3_phi = 3.0 * ((
						//	WEIGTH[TNE] * ((((*phaseField)(x1 + 1, x2 + 1, x3 + 1) - (*phaseField)(x1 - 1, x2 - 1, x3 - 1)) + ((*phaseField)(x1 - 1, x2 + 1, x3 + 1) - (*phaseField)(x1 + 1, x2 - 1, x3 - 1)))
						//	+ (((*phaseField)(x1 - 1, x2 - 1, x3 + 1) - (*phaseField)(x1 + 1, x2 + 1, x3 - 1)) + ((*phaseField)(x1 + 1, x2 - 1, x3 + 1) - (*phaseField)(x1 - 1, x2 + 1, x3 - 1))))
						//	+ WEIGTH[NE] * ((((*phaseField)(x1 + 1, x2, x3+1) - (*phaseField)(x1 - 1, x2, x3-1)) + ((*phaseField)(x1 - 1, x2, x3+1) - (*phaseField)(x1 + 1, x2, x3-1)))
						//		+ (((*phaseField)(x1, x2 - 1, x3 + 1) - (*phaseField)(x1, x2 + 1, x3 - 1)) + ((*phaseField)(x1, x2 + 1, x3 + 1) - (*phaseField)(x1, x2 - 1, x3 - 1)))))
						//	+ WEIGTH[N] * ((*phaseField)(x1, x2, x3+1) - (*phaseField)(x1, x2, x3-1))
						//	);

						///////////////////////////////////////

						//LBMReal dX1_phi2 = 1.5 * ((
						//	WEIGTH[TNE] * ((((*phaseField)(x1 + 2, x2 + 2, x3 + 2) - (*phaseField)(x1 - 2, x2 - 2, x3 - 2)) + ((*phaseField)(x1 + 2, x2 - 2, x3 + 2) - (*phaseField)(x1 - 2, x2 + 2, x3 - 2)))
						//		+ (((*phaseField)(x1 + 2, x2 - 2, x3 - 2) - (*phaseField)(x1 - 2, x2 + 2, x3 + 2)) + ((*phaseField)(x1 + 2, x2 + 2, x3 - 2) - (*phaseField)(x1 - 2, x2 - 2, x3 + 2))))
						//	+ WEIGTH[NE] * ((((*phaseField)(x1 + 2, x2 + 2, x3) - (*phaseField)(x1 - 2, x2 - 2, x3)) + ((*phaseField)(x1 + 2, x2 - 2, x3) - (*phaseField)(x1 - 2, x2 + 2, x3)))
						//		+ (((*phaseField)(x1 + 2, x2, x3 - 2) - (*phaseField)(x1 - 2, x2, x3 + 2)) + ((*phaseField)(x1 + 2, x2, x3 + 2) - (*phaseField)(x1 - 2, x2, x3 - 2)))))
						//	+ WEIGTH[N] * ((*phaseField)(x1 + 2, x2, x3) - (*phaseField)(x1 - 2, x2, x3))
						//	);
						////if (dX1_phi != NdX1_phi) {std::cout<<dX1_phi<<" "<< NdX1_phi<<std::endl;}

						//LBMReal dX2_phi2 = 1.5 * ((
						//	WEIGTH[TNE] * ((((*phaseField)(x1 + 2, x2 + 2, x3 + 2) - (*phaseField)(x1 - 2, x2 - 2, x3 - 2)) + ((*phaseField)(x1 - 2, x2 + 2, x3 + 2) - (*phaseField)(x1 + 2, x2 - 2, x3 - 2)))
						//		+ (((*phaseField)(x1 - 2, x2 + 2, x3 - 2) - (*phaseField)(x1 + 2, x2 - 2, x3 + 2)) + ((*phaseField)(x1 + 2, x2 + 2, x3 - 2) - (*phaseField)(x1 - 2, x2 - 2, x3 + 2))))
						//	+ WEIGTH[NE] * ((((*phaseField)(x1 + 2, x2 + 2, x3) - (*phaseField)(x1 - 2, x2 - 2, x3)) + ((*phaseField)(x1 - 2, x2 + 2, x3) - (*phaseField)(x1 + 2, x2 - 2, x3)))
						//		+ (((*phaseField)(x1, x2 + 2, x3 - 2) - (*phaseField)(x1, x2 - 2, x3 + 2)) + ((*phaseField)(x1, x2 + 2, x3 + 2) - (*phaseField)(x1, x2 - 2, x3 - 2)))))
						//	+ WEIGTH[N] * ((*phaseField)(x1, x2 + 2, x3) - (*phaseField)(x1, x2 - 2, x3))
						//	);

						//LBMReal dX3_phi2 = 1.5 * ((
						//	WEIGTH[TNE] * ((((*phaseField)(x1 + 2, x2 + 2, x3 + 2) - (*phaseField)(x1 - 2, x2 - 2, x3 - 2)) + ((*phaseField)(x1 - 2, x2 + 2, x3 + 2) - (*phaseField)(x1 + 2, x2 - 2, x3 - 2)))
						//		+ (((*phaseField)(x1 - 2, x2 - 2, x3 + 2) - (*phaseField)(x1 + 2, x2 + 2, x3 - 2)) + ((*phaseField)(x1 + 2, x2 - 2, x3 + 2) - (*phaseField)(x1 - 2, x2 + 2, x3 - 2))))
						//	+ WEIGTH[NE] * ((((*phaseField)(x1 + 2, x2, x3 + 2) - (*phaseField)(x1 - 2, x2, x3 - 2)) + ((*phaseField)(x1 - 2, x2, x3 + 2) - (*phaseField)(x1 + 2, x2, x3 - 2)))
						//		+ (((*phaseField)(x1, x2 - 2, x3 + 2) - (*phaseField)(x1, x2 + 2, x3 - 2)) + ((*phaseField)(x1, x2 + 2, x3 + 2) - (*phaseField)(x1, x2 - 2, x3 - 2)))))
						//	+ WEIGTH[N] * ((*phaseField)(x1, x2, x3 + 2) - (*phaseField)(x1, x2, x3 - 2))
						//	);

						//dX1_phi = (2*dX1_phi -1*dX1_phi2);// 2 * dX1_phi - dX1_phi2;
						//dX2_phi = (2*dX2_phi -1*dX2_phi2);// 2 * dX2_phi - dX2_phi2;
						//dX3_phi = (2*dX3_phi -1*dX3_phi2);// 2 * dX3_phi - dX3_phi2;


                        LBMReal denom = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi) + 1e-9;
                        LBMReal normX1 = dX1_phi/denom;
						LBMReal normX2 = dX2_phi/denom;
						LBMReal normX3 = dX3_phi/denom; 


						///test for magnitude of gradient from phase indicator directly
						//if (fabs((1.0 - phi[REST]) * (phi[REST]) */* c4*/ - (denom- 1e-9)) / denom > 1e-3 &&phi[REST]>0.4 &&phi[REST]<0.6) {
						//	std::cout << (1.0 - phi[REST]) * (phi[REST])  // *c4 
						//		<< " " << denom <<" "<< ((1.0 - phi[REST]) * (phi[REST]) * c4 ) / denom << std::endl;
						//}
						//dX1_phi = (1.0 - phi[REST]) * (phi[REST]) /* c4 */* normX1;
						//dX2_phi = (1.0 - phi[REST]) * (phi[REST]) /* c4 */* normX2;
						//dX3_phi = (1.0 - phi[REST]) * (phi[REST]) /* c4 */* normX3;

						//denom = 1.0;

						///!test

						collFactorM = collFactorL + (collFactorL - collFactorG) * (phi[DIR_000] - phiH) / (phiH - phiL);
						//collFactorM = phi[REST] - phiL < (phiH - phiL) * 0.05 ? collFactorG : collFactorL;

                        LBMReal mu = 2 * beta * phi[DIR_000] * (phi[DIR_000] - 1) * (2 * phi[DIR_000] - 1) - kappa * nabla2_phi();

                        //----------- Calculating Macroscopic Values -------------
                        LBMReal rho = rhoH + rhoToPhi * (phi[DIR_000] - phiH);

						if (withForcing) {
							// muX1 = static_cast<double>(x1-1+ix1*maxX1);
							// muX2 = static_cast<double>(x2-1+ix2*maxX2);
							// muX3 = static_cast<double>(x3-1+ix3*maxX3);

							forcingX1 = muForcingX1.Eval();
							forcingX2 = muForcingX2.Eval();
							forcingX3 = muForcingX3.Eval();

							LBMReal rho_m = 1.0 / densityRatio;
							forcingX1 = forcingX1 * (rho - rho_m);
							forcingX2 = forcingX2 * (rho - rho_m);
							forcingX3 = forcingX3 * (rho - rho_m);
						}
                            			   ////Incompressible Kernal

			    mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)/rho*c3;
			    mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) / rho * c3;
			    mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) / rho * c3;
			    mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) / rho * c3;
			    mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) / rho * c3;
			    mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) / rho * c3;
			    mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) / rho * c3;
			    mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) / rho * c3;
			    mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) / rho * c3;
			    mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) / rho * c3;
			    mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) / rho * c3;
			    mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) / rho * c3;
			    mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) / rho * c3;

			    mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) / rho * c3;
			    mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) / rho * c3;
			    mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) / rho * c3;
			    mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) / rho * c3;
			    mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) / rho * c3;
			    mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) / rho * c3;
			    mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) / rho * c3;
			    mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) / rho * c3;
			    mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) / rho * c3;
			    mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) / rho * c3;
			    mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) / rho * c3;
			    mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) / rho * c3;
			    mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) / rho * c3;

			    mfbbb = (*this->zeroDistributionsF)(x1, x2, x3) / rho * c3;





			   LBMReal m0, m1, m2;
			   LBMReal rhoRef=c1;

			  //LBMReal
			  // FIXME: warning: unused variable 'drho'
//			   LBMReal drho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
//				   + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
//				   + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

			   LBMReal vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				   (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
				   (mfcbb - mfabb))/rhoRef;
			   LBMReal vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				   (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
				   (mfbcb - mfbab))/rhoRef;
			   LBMReal vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				   (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
				   (mfbbc - mfbba))/rhoRef;

			   ///surface tension force
			   vvx += mu * dX1_phi*c1o2;
			   vvy += mu * dX2_phi * c1o2;
			   vvz += mu * dX3_phi * c1o2;
			  


			   ////Velocity filter 14.04.2021
			  // LBMReal lap_vx, lap_vy,lap_vz;
			  // {
				 //  LBMReal sum = 0.0;
				 //  sum += WEIGTH[TNE] * (((((*velocityX)(x1+1, x2+1, x3+1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2 - 1, x3 - 1) - (*velocityX)(x1, x2, x3))) + (((*velocityX)(x1 + 1, x2 + 1, x3 - 1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 + 1, x2 - 1, x3 + 1) - (*velocityX)(x1, x2, x3))))
					//   + ((((*velocityX)(x1 + 1, x2 - 1, x3 + 1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2 + 1, x3 - 1) - (*velocityX)(x1, x2, x3))) + (((*velocityX)(x1 - 1, x2 + 1, x3 + 1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 + 1, x2 - 1, x3 - 1) - (*velocityX)(x1, x2, x3)))));
				 //  sum += WEIGTH[TN] * (
					//   ((((*velocityX)(x1 + 1, x2 + 1, x3 ) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2 - 1, x3) - (*velocityX)(x1, x2, x3))) + (((*velocityX)(x1 + 1, x2 - 1, x3) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2 + 1, x3) - (*velocityX)(x1, x2, x3))))
					//   + ((((*velocityX)(x1 + 1, x2 , x3+1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2 , x3-1) - (*velocityX)(x1, x2, x3))) + (((*velocityX)(x1 +1 , x2 , x3-1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 - 1, x2, x3 + 1) - (*velocityX)(x1, x2, x3))))
					//   + ((((*velocityX)(x1 , x2+1, x3 + 1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1, x2 - 1, x3 - 1) - (*velocityX)(x1, x2, x3))) + (((*velocityX)(x1, x2 + 1, x3 - 1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1, x2 - 1, x3 + 1) - (*velocityX)(x1, x2, x3))))
					//   );
				 //  sum += WEIGTH[T] * (
					//   (((*velocityX)(x1-1, x2 , x3 ) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1 + 1, x2, x3) - (*velocityX)(x1, x2, x3)))
					//   + (((*velocityX)(x1 , x2-1, x3) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1, x2 + 1, x3) - (*velocityX)(x1, x2, x3)))
					//   + (((*velocityX)(x1, x2, x3-1) - (*velocityX)(x1, x2, x3)) + ((*velocityX)(x1, x2, x3+1) - (*velocityX)(x1, x2, x3)))
					//   );
				 //  //for (int k = FSTARTDIR; k <= FENDDIR; k++) {
				 //  //    sum += WEIGTH[k] * (phi[k] - phi[REST]);
				 //  //}
				 //   lap_vx=6.0 * sum;

					//sum = 0.0;
					//sum += WEIGTH[TNE] * (((((*velocityY)(x1 + 1, x2 + 1, x3 + 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2 - 1, x3 - 1) - (*velocityY)(x1, x2, x3))) + (((*velocityY)(x1 + 1, x2 + 1, x3 - 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 + 1, x2 - 1, x3 + 1) - (*velocityY)(x1, x2, x3))))
					//	+ ((((*velocityY)(x1 + 1, x2 - 1, x3 + 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2 + 1, x3 - 1) - (*velocityY)(x1, x2, x3))) + (((*velocityY)(x1 - 1, x2 + 1, x3 + 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 + 1, x2 - 1, x3 - 1) - (*velocityY)(x1, x2, x3)))));
					//sum += WEIGTH[TN] * (
					//	((((*velocityY)(x1 + 1, x2 + 1, x3) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2 - 1, x3) - (*velocityY)(x1, x2, x3))) + (((*velocityY)(x1 + 1, x2 - 1, x3) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2 + 1, x3) - (*velocityY)(x1, x2, x3))))
					//	+ ((((*velocityY)(x1 + 1, x2, x3 + 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2, x3 - 1) - (*velocityY)(x1, x2, x3))) + (((*velocityY)(x1 + 1, x2, x3 - 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 - 1, x2, x3 + 1) - (*velocityY)(x1, x2, x3))))
					//	+ ((((*velocityY)(x1, x2 + 1, x3 + 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1, x2 - 1, x3 - 1) - (*velocityY)(x1, x2, x3))) + (((*velocityY)(x1, x2 + 1, x3 - 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1, x2 - 1, x3 + 1) - (*velocityY)(x1, x2, x3))))
					//	);
					//sum += WEIGTH[T] * (
					//	(((*velocityY)(x1 - 1, x2, x3) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1 + 1, x2, x3) - (*velocityY)(x1, x2, x3)))
					//	+ (((*velocityY)(x1, x2 - 1, x3) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1, x2 + 1, x3) - (*velocityY)(x1, x2, x3)))
					//	+ (((*velocityY)(x1, x2, x3 - 1) - (*velocityY)(x1, x2, x3)) + ((*velocityY)(x1, x2, x3 + 1) - (*velocityY)(x1, x2, x3)))
					//	);

					//lap_vy = 6.0 * sum;

					//sum = 0.0;
					//sum += WEIGTH[TNE] * (((((*velocityZ)(x1 + 1, x2 + 1, x3 + 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2 - 1, x3 - 1) - (*velocityZ)(x1, x2, x3))) + (((*velocityZ)(x1 + 1, x2 + 1, x3 - 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 + 1, x2 - 1, x3 + 1) - (*velocityZ)(x1, x2, x3))))
					//	+ ((((*velocityZ)(x1 + 1, x2 - 1, x3 + 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2 + 1, x3 - 1) - (*velocityZ)(x1, x2, x3))) + (((*velocityZ)(x1 - 1, x2 + 1, x3 + 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 + 1, x2 - 1, x3 - 1) - (*velocityZ)(x1, x2, x3)))));
					//sum += WEIGTH[TN] * (
					//	((((*velocityZ)(x1 + 1, x2 + 1, x3) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2 - 1, x3) - (*velocityZ)(x1, x2, x3))) + (((*velocityZ)(x1 + 1, x2 - 1, x3) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2 + 1, x3) - (*velocityZ)(x1, x2, x3))))
					//	+ ((((*velocityZ)(x1 + 1, x2, x3 + 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2, x3 - 1) - (*velocityZ)(x1, x2, x3))) + (((*velocityZ)(x1 + 1, x2, x3 - 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 - 1, x2, x3 + 1) - (*velocityZ)(x1, x2, x3))))
					//	+ ((((*velocityZ)(x1, x2 + 1, x3 + 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1, x2 - 1, x3 - 1) - (*velocityZ)(x1, x2, x3))) + (((*velocityZ)(x1, x2 + 1, x3 - 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1, x2 - 1, x3 + 1) - (*velocityZ)(x1, x2, x3))))
					//	);
					//sum += WEIGTH[T] * (
					//	(((*velocityZ)(x1 - 1, x2, x3) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1 + 1, x2, x3) - (*velocityZ)(x1, x2, x3)))
					//	+ (((*velocityZ)(x1, x2 - 1, x3) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1, x2 + 1, x3) - (*velocityZ)(x1, x2, x3)))
					//	+ (((*velocityZ)(x1, x2, x3 - 1) - (*velocityZ)(x1, x2, x3)) + ((*velocityZ)(x1, x2, x3 + 1) - (*velocityZ)(x1, x2, x3)))
					//	);

					//lap_vz = 6.0 * sum;

			  // }

			  // if (lap_vx != 0.0) {
				 //  lap_vx = lap_vx;
			  // }

			   ///----Classic source term 8.4.2021

			   LBMReal vvxF, vvyF, vvzF;
			   vvxF = vvx;//-2*c1o24 * lap_vx;// 
			   vvyF = vvy;//-2*c1o24 * lap_vy;// 
			   vvzF = vvz;//-2*c1o24 * lap_vz;// 

//			   vvxF = 1.2* vvx- 0.2*0.5 * ((*velocityX)(x1 - 1, x2, x3) + (*velocityX)(x1 + 1, x2, x3));
//			   vvyF = 1.2 *vvy- 0.2*0.5* ((*velocityY)(x1 , x2-1, x3) + (*velocityY)(x1 , x2+1, x3));
//			   vvzF = 1.2 *vvz-0.2*0.5* ((*velocityZ)(x1 , x2, x3-1) + (*velocityZ)(x1 , x2, x3+1));
			   //if (vvxF != vvx) {
				  // vvxF = vvxF;
			   //}
			   LBMReal weightGrad =  1.0-denom*denom/(denom*denom+0.0001*0.001);
			   LBMReal dX1_phiF = dX1_phi * weightGrad + (1.0 - weightGrad) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * normX1;
			   LBMReal dX2_phiF = dX2_phi * weightGrad + (1.0 - weightGrad) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * normX2;
			   LBMReal dX3_phiF = dX3_phi * weightGrad + (1.0 - weightGrad) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * normX3;

			   //dX1_phiF *= 1.2;
			   //dX2_phiF *= 1.2;
			   //dX3_phiF *= 1.2;

			   //LBMReal gradFD = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);
			   //LBMReal gradPhi = (1.0 - phi[REST]) * (phi[REST]);
			   //gradPhi = (gradPhi > gradFD) ? gradPhi : gradFD;
			   //dX1_phiF = gradPhi * normX1;
				  // dX2_phiF = gradPhi * normX2;
				  // dX3_phiF = gradPhi * normX3;

			   LBMReal ux2;
			   LBMReal uy2;
			   LBMReal uz2;
			   ux2 = vvxF * vvxF;
			   uy2 = vvyF * vvyF;
			   uz2 = vvzF * vvzF;
			   LBMReal forcingTerm[D3Q27System::ENDF + 1];
			   for (int dir = STARTF; dir <= (FENDDIR); dir++) {
				   LBMReal velProd = DX1[dir] * vvxF + DX2[dir] * vvyF + DX3[dir] * vvzF;
				   LBMReal velSq1 = velProd * velProd;
				   LBMReal gamma = WEIGTH[dir] * (1.0 + 3 * velProd + (4.5 * velSq1 - 1.5 * (ux2 + uy2 + uz2)));

				   LBMReal fac1 = (gamma - WEIGTH[dir]) * c1o3 * rhoToPhi;

				   forcingTerm[dir] = 
					   (-vvxF) * (fac1 * dX1_phiF ) +
					   (-vvyF) * (fac1 * dX2_phiF ) +
					   (-vvzF) * (fac1 * dX3_phiF ) +
					   (DX1[dir]) * (fac1 * dX1_phiF ) +
					   (DX2[dir]) * (fac1 * dX2_phiF ) +
					   (DX3[dir]) * (fac1 * dX3_phiF );

				   //LBMReal biDif= (-((*phaseField)(x1 + 2 * DX1[dir], x2 + 2 * DX2[dir], x3 + 2 * DX3[dir])) + 4 * ((*phaseField)(x1 + DX1[dir], x2 + DX2[dir], x3 + DX3[dir]))
					  // - 3*((*phaseField)(x1 , x2 , x3 )) )*0.5;
				   //LBMReal ceDif = (((*phaseField)(x1 + DX1[dir], x2 + DX2[dir], x3 + DX3[dir])) - ((*phaseField)(x1 - DX1[dir], x2 - DX2[dir], x3 - DX3[dir]))) * 0.5;

				   ////ceDif = ((((*phaseField)(x1 + 2*DX1[dir], x2 + 2*DX2[dir], x3 + 2*DX3[dir])) - ((*phaseField)(x1 , x2 , x3 ))) * biDif < 0) ?
					  //// (!bcArray->isSolid(x1+2*DX1[dir], x2+2*DX2[dir], x3+2*DX3[dir]) && !bcArray->isUndefined(x1 + 2 * DX1[dir], x2 + 2 * DX2[dir], x3 + 2 * DX3[dir]) && !bcArray->isSolid(x1 + DX1[dir], x2 +  DX2[dir], x3 +  DX3[dir]) && !bcArray->isUndefined(x1 +  DX1[dir], x2 + DX2[dir], x3 + DX3[dir]) && !bcArray->isSolid(x1 - DX1[dir], x2 - DX2[dir], x3 - DX3[dir]) && !bcArray->isUndefined(x1 - DX1[dir], x2 - DX2[dir], x3 - DX3[dir])) ?
					  //// (biDif+ceDif)*0.5 : ceDif: ceDif;

				   //ceDif = ((((*phaseField)(x1 + 2 * DX1[dir], x2 + 2 * DX2[dir], x3 + 2 * DX3[dir])) - ((*phaseField)(x1, x2, x3))) * biDif < 0) ? biDif : ceDif;

				   //forcingTerm[dir] =
					  // (-vvxF) * (fac1 * dX1_phiF) +
					  // (-vvyF) * (fac1 * dX2_phiF) +
					  // (-vvzF) * (fac1 * dX3_phiF) +
					  // fac1 * ceDif;//(((*phaseField)(x1 + DX1[dir], x2 + DX2[dir], x3 + DX3[dir])) -  ((*phaseField)(x1 - DX1[dir], x2 - DX2[dir], x3 - DX3[dir]))) * 0.5;
					  // //( -((*phaseField)(x1 +2* DX1[dir], x2 + 2 * DX2[dir], x3 + 2 * DX3[dir])) + 5*((*phaseField)(x1 + DX1[dir], x2 +  DX2[dir], x3 +  DX3[dir])) 
						 //  //- 3*((*phaseField)(x1 , x2 , x3 )) - ((*phaseField)(x1 - DX1[dir], x2 - DX2[dir], x3 - DX3[dir])) )*0.25;


			   }

			   LBMReal gamma = WEIGTH[DIR_000] * (1.0 - 1.5 * (ux2 + uy2 + uz2));
			   LBMReal fac1 = (gamma - WEIGTH[DIR_000]) * c1o3 * rhoToPhi;
			   forcingTerm[DIR_000] = (-vvxF) * (fac1 * dX1_phiF ) +
				   (-vvyF) * (fac1 * dX2_phiF ) +
				   (-vvzF) * (fac1 * dX3_phiF );

			   ////////
			  // LBMReal divAfterSource=
			  //( mfcbb + 3.0 * (0.5 * forcingTerm[DIR_P00]) / rho	) *((vvxF-1)*(vvxF-1)+(vvyF)  *(vvyF)  +(vvzF)  *(vvzF)-1)+
			  //( mfbcb + 3.0 * (0.5 * forcingTerm[N]) / rho	) *((vvxF)  *(vvxF)  +(vvyF-1)*(vvyF-1)+(vvzF)  *(vvzF)-1)+
			  //( mfbbc + 3.0 * (0.5 * forcingTerm[T]) / rho	) *((vvxF)  *(vvxF)  +(vvyF)  *(vvyF)  +(vvzF-1)*(vvzF-1)-1)+
			  //( mfccb + 3.0 * (0.5 * forcingTerm[NE]) / rho	) *((vvxF-1)*(vvxF-1)+(vvyF-1)*(vvyF-1)+(vvzF)  *(vvzF)-1)+
			  //( mfacb + 3.0 * (0.5 * forcingTerm[NW]) / rho	) *((vvxF+1)*(vvxF+1)+(vvyF-1)*(vvyF-1)+(vvzF)  *(vvzF)-1)+
			  //( mfcbc + 3.0 * (0.5 * forcingTerm[TE]) / rho	) *((vvxF-1)*(vvxF-1)+(vvyF)  *(vvyF)  +(vvzF-1)*(vvzF-1)-1)+
			  //( mfabc + 3.0 * (0.5 * forcingTerm[TW]) / rho	) *((vvxF+1)*(vvxF+1)+(vvyF)  *(vvyF)  +(vvzF-1)*(vvzF-1)-1)+
			  //( mfbcc + 3.0 * (0.5 * forcingTerm[TN]) / rho	) *((vvxF)  *(vvxF)  +(vvyF-1)*(vvyF-1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfbac + 3.0 * (0.5 * forcingTerm[TS]) / rho	) *((vvxF)  *(vvxF)  +(vvyF+1)*(vvyF+1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfccc + 3.0 * (0.5 * forcingTerm[TNE]) / rho) *((vvxF-1)*(vvxF-1)+(vvyF-1)*(vvyF-1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfacc + 3.0 * (0.5 * forcingTerm[TNW]) / rho) *((vvxF+1)*(vvxF+1)+(vvyF-1)*(vvyF-1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfcac + 3.0 * (0.5 * forcingTerm[TSE]) / rho) *((vvxF-1)*(vvxF-1)+(vvyF+1)*(vvyF+1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfaac + 3.0 * (0.5 * forcingTerm[TSW]) / rho) *((vvxF+1)*(vvxF+1)+(vvyF+1)*(vvyF+1)+(vvzF-1)*(vvzF-1)-1)+
			  //( mfabb + 3.0 * (0.5 * forcingTerm[W]) / rho	) *((vvxF+1)*(vvxF+1)+(vvyF)  *(vvyF)  +(vvzF)  *(vvzF)-1)+
			  //( mfbab + 3.0 * (0.5 * forcingTerm[S]) / rho	) *((vvxF)  *(vvxF)  +(vvyF+1)*(vvyF+1)+(vvzF)  *(vvzF)-1)+
			  //( mfbba + 3.0 * (0.5 * forcingTerm[B]) / rho	) *((vvxF)  *(vvxF)  +(vvyF)  *(vvyF)  +(vvzF+1)*(vvzF+1)-1)+
			  //( mfaab + 3.0 * (0.5 * forcingTerm[SW]) / rho	) *((vvxF+1)*(vvxF+1)+(vvyF+1)*(vvyF+1)+(vvzF)  *(vvzF)-1)+
			  //( mfcab + 3.0 * (0.5 * forcingTerm[SE]) / rho	) *((vvxF-1)*(vvxF-1)+(vvyF+1)*(vvyF+1)+(vvzF)  *(vvzF)-1)+
			  //( mfaba + 3.0 * (0.5 * forcingTerm[BW]) / rho	) *((vvxF+1)*(vvxF+1)+(vvyF)  *(vvyF)  +(vvzF+1)*(vvzF+1)-1)+
			  //( mfcba + 3.0 * (0.5 * forcingTerm[BE]) / rho	) *((vvxF-1)*(vvxF-1)+(vvyF)  *(vvyF)  +(vvzF+1)*(vvzF+1)-1)+
			  //( mfbaa + 3.0 * (0.5 * forcingTerm[BS]) / rho	) *((vvxF)  *(vvxF)  +(vvyF+1)*(vvyF+1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfbca + 3.0 * (0.5 * forcingTerm[BN]) / rho	) *((vvxF)  *(vvxF)  +(vvyF-1)*(vvyF-1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfaaa + 3.0 * (0.5 * forcingTerm[BSW]) / rho) *((vvxF+1)*(vvxF+1)+(vvyF+1)*(vvyF+1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfcaa + 3.0 * (0.5 * forcingTerm[BSE]) / rho) *((vvxF-1)*(vvxF-1)+(vvyF+1)*(vvyF+1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfaca + 3.0 * (0.5 * forcingTerm[BNW]) / rho) *((vvxF+1)*(vvxF+1)+(vvyF-1)*(vvyF-1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfcca + 3.0 * (0.5 * forcingTerm[BNE]) / rho) *((vvxF-1)*(vvxF-1)+(vvyF-1)*(vvyF-1)+(vvzF+1)*(vvzF+1)-1)+
			  //( mfbbb + 3.0 * (0.5 * forcingTerm[REST]) / rho)*((vvxF)*(vvxF)+(vvyF)*(vvyF)+(vvzF)*(vvzF)-1);

			  // LBMReal divBeforeSource =
				 //  (mfcbb)    * ((vvxF - 1) * (vvxF - 1) + (vvyF) * (vvyF)+(vvzF) * (vvzF)-1) +
				 //  (mfbcb)    * ((vvxF) * (vvxF)+(vvyF - 1) * (vvyF - 1) + (vvzF) * (vvzF)-1) +
				 //  (mfbbc)    * ((vvxF) * (vvxF)+(vvyF) * (vvyF)+(vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfccb)   * ((vvxF - 1) * (vvxF - 1) + (vvyF - 1) * (vvyF - 1) + (vvzF) * (vvzF)-1) +
				 //  (mfacb)   * ((vvxF + 1) * (vvxF + 1) + (vvyF - 1) * (vvyF - 1) + (vvzF) * (vvzF)-1) +
				 //  (mfcbc)   * ((vvxF - 1) * (vvxF - 1) + (vvyF) * (vvyF)+(vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfabc)   * ((vvxF + 1) * (vvxF + 1) + (vvyF) * (vvyF)+(vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfbcc)   * ((vvxF) * (vvxF)+(vvyF - 1) * (vvyF - 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfbac)   * ((vvxF) * (vvxF)+(vvyF + 1) * (vvyF + 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfccc)  * ((vvxF - 1) * (vvxF - 1) + (vvyF - 1) * (vvyF - 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfacc)  * ((vvxF + 1) * (vvxF + 1) + (vvyF - 1) * (vvyF - 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfcac)  * ((vvxF - 1) * (vvxF - 1) + (vvyF + 1) * (vvyF + 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfaac)  * ((vvxF + 1) * (vvxF + 1) + (vvyF + 1) * (vvyF + 1) + (vvzF - 1) * (vvzF - 1)-1) +
				 //  (mfabb)    * ((vvxF + 1) * (vvxF + 1) + (vvyF) * (vvyF)+(vvzF) * (vvzF)-1) +
				 //  (mfbab)    * ((vvxF) * (vvxF)+(vvyF + 1) * (vvyF + 1) + (vvzF) * (vvzF)-1) +
				 //  (mfbba)    * ((vvxF) * (vvxF)+(vvyF) * (vvyF)+(vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfaab)   * ((vvxF + 1) * (vvxF + 1) + (vvyF + 1) * (vvyF + 1) + (vvzF) * (vvzF)-1) +
				 //  (mfcab)   * ((vvxF - 1) * (vvxF - 1) + (vvyF + 1) * (vvyF + 1) + (vvzF) * (vvzF)-1) +
				 //  (mfaba)   * ((vvxF + 1) * (vvxF + 1) + (vvyF) * (vvyF)+(vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfcba)   * ((vvxF - 1) * (vvxF - 1) + (vvyF) * (vvyF)+(vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfbaa)   * ((vvxF) * (vvxF)+(vvyF + 1) * (vvyF + 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfbca)   * ((vvxF) * (vvxF)+(vvyF - 1) * (vvyF - 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfaaa)  * ((vvxF + 1) * (vvxF + 1) + (vvyF + 1) * (vvyF + 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfcaa)  * ((vvxF - 1) * (vvxF - 1) + (vvyF + 1) * (vvyF + 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfaca)  * ((vvxF + 1) * (vvxF + 1) + (vvyF - 1) * (vvyF - 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfcca)  * ((vvxF - 1) * (vvxF - 1) + (vvyF - 1) * (vvyF - 1) + (vvzF + 1) * (vvzF + 1)-1) +
				 //  (mfbbb) * ((vvxF) * (vvxF)+(vvyF) * (vvyF)+(vvzF) * (vvzF)-1);
			   //if (divAfterSource - divBeforeSource != 0 && phi[REST]>0.0001 && phi[REST]<0.999) {
				  // std::cout << phi[REST]<<" "<< divAfterSource << " " << divBeforeSource <<" "<< divAfterSource/ divBeforeSource << std::endl;
			   //}

			   //if (fabs(divAfterSource - divBeforeSource)/(fabs(divAfterSource) + fabs(divBeforeSource)+1e-10) > 1e-5) {
				  // LBMReal scaleDiv =0.95+(1-0.95)* (divBeforeSource) / (divBeforeSource - divAfterSource);

				  // forcingTerm[DIR_P00]	 *=scaleDiv;
				  // forcingTerm[N]	 *=scaleDiv;
				  // forcingTerm[T]	 *=scaleDiv;
				  // forcingTerm[NE]	 *=scaleDiv;
				  // forcingTerm[NW]	 *=scaleDiv;
				  // forcingTerm[TE]	 *=scaleDiv;
				  // forcingTerm[TW]	 *=scaleDiv;
				  // forcingTerm[TN]	 *=scaleDiv;
				  // forcingTerm[TS]	 *=scaleDiv;
				  // forcingTerm[TNE]	 *=scaleDiv;
				  // forcingTerm[TNW]	 *=scaleDiv;
				  // forcingTerm[TSE]	 *=scaleDiv;
				  // forcingTerm[TSW]	 *=scaleDiv;
				  // forcingTerm[W]	 *=scaleDiv;
				  // forcingTerm[S]	 *=scaleDiv;
				  // forcingTerm[B]	 *=scaleDiv;
				  // forcingTerm[SW]	 *=scaleDiv;
				  // forcingTerm[SE]	 *=scaleDiv;
				  // forcingTerm[BW]	 *=scaleDiv;
				  // forcingTerm[BE]	 *=scaleDiv;
				  // forcingTerm[BS]	 *=scaleDiv;
				  // forcingTerm[BN]	 *=scaleDiv;
				  // forcingTerm[BSW]	 *=scaleDiv;
				  // forcingTerm[BSE]	 *=scaleDiv;
				  // forcingTerm[BNW]	 *=scaleDiv;
				  // forcingTerm[BNE]	 *=scaleDiv;
				  // forcingTerm[REST] *=scaleDiv;
			   //}
			   ////////


			   mfcbb +=3.0 * ( 0.5 * forcingTerm[DIR_P00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[E  ];
			   mfbcb +=3.0 * ( 0.5 * forcingTerm[DIR_0P0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[N  ];
			   mfbbc +=3.0 * ( 0.5 * forcingTerm[DIR_00P]) / rho;    //-(3.0*p1 - rho)*WEIGTH[T  ];
			   mfccb +=3.0 * ( 0.5 * forcingTerm[DIR_PP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NE ];
			   mfacb +=3.0 * ( 0.5 * forcingTerm[DIR_MP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NW ];
			   mfcbc +=3.0 * ( 0.5 * forcingTerm[DIR_P0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TE ];
			   mfabc +=3.0 * ( 0.5 * forcingTerm[DIR_M0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TW ];
			   mfbcc +=3.0 * ( 0.5 * forcingTerm[DIR_0PP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TN ];
			   mfbac +=3.0 * ( 0.5 * forcingTerm[DIR_0MP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TS ];
			   mfccc +=3.0 * ( 0.5 * forcingTerm[DIR_PPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNE];
			   mfacc +=3.0 * ( 0.5 * forcingTerm[DIR_MPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNW];
			   mfcac +=3.0 * ( 0.5 * forcingTerm[DIR_PMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSE];
			   mfaac +=3.0 * ( 0.5 * forcingTerm[DIR_MMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSW];
			   mfabb +=3.0 * ( 0.5 * forcingTerm[DIR_M00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[W  ];
			   mfbab +=3.0 * ( 0.5 * forcingTerm[DIR_0M0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[S  ];
			   mfbba +=3.0 * ( 0.5 * forcingTerm[DIR_00M]) / rho;    //-(3.0*p1 - rho)*WEIGTH[B  ];
			   mfaab +=3.0 * ( 0.5 * forcingTerm[DIR_MM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SW ];
			   mfcab +=3.0 * ( 0.5 * forcingTerm[DIR_PM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SE ];
			   mfaba +=3.0 * ( 0.5 * forcingTerm[DIR_M0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BW ];
			   mfcba +=3.0 * ( 0.5 * forcingTerm[DIR_P0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BE ];
			   mfbaa +=3.0 * ( 0.5 * forcingTerm[DIR_0MM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BS ];
			   mfbca +=3.0 * ( 0.5 * forcingTerm[DIR_0PM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BN ];
			   mfaaa +=3.0 * ( 0.5 * forcingTerm[DIR_MMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSW];
			   mfcaa +=3.0 * ( 0.5 * forcingTerm[DIR_PMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSE];
			   mfaca +=3.0 * ( 0.5 * forcingTerm[DIR_MPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNW];
			   mfcca +=3.0 * ( 0.5 * forcingTerm[DIR_PPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNE];
			   mfbbb +=3.0 * ( 0.5 * forcingTerm[DIR_000]) / rho; //- (3.0*p1 - rho)*WEIGTH[REST]

			   //--------------------------------------------------------


			   //////////End classic source term
			   //forcing 
			   ///////////////////////////////////////////////////////////////////////////////////////////
			   if (withForcing)
			   {
				   muX1 = static_cast<double>(x1 - 1 + ix1 * maxX1);
				   muX2 = static_cast<double>(x2 - 1 + ix2 * maxX2);
				   muX3 = static_cast<double>(x3 - 1 + ix3 * maxX3);

				   forcingX1 = muForcingX1.Eval();
				   forcingX2 = muForcingX2.Eval();
				   forcingX3 = muForcingX3.Eval();

				   vvx += forcingX1 * deltaT * 0.5; // X
				   vvy += forcingX2 * deltaT * 0.5; // Y
				   vvz += forcingX3 * deltaT * 0.5; // Z
			   }

			   LBMReal vx2;
			   LBMReal vy2;
			   LBMReal vz2;
			   vx2 = vvx * vvx;
			   vy2 = vvy * vvy;
			   vz2 = vvz * vvz;

			   ///////



			   ///////////////////////////////////////////////////////////////////////////////////////////               
			   LBMReal oMdrho;


			   oMdrho = mfccc + mfaaa;
			   m0 = mfaca + mfcac;
			   m1 = mfacc + mfcaa;
			   m2 = mfaac + mfcca;
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
			   m0 += mfbbb; //hat gefehlt
			   oMdrho = (rhoRef - (oMdrho + m0))/rhoRef;// 12.03.21 check derivation!!!!


			   ////////////////////////////////////////////////////////////////////////////////////
			   LBMReal wadjust;
			   LBMReal qudricLimit = 0.01;
			   ////////////////////////////////////////////////////////////////////////////////////
			   //Hin
			   ////////////////////////////////////////////////////////////////////////////////////
			   // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   m2 = mfaaa + mfaac;
			   m1 = mfaac - mfaaa;
			   m0 = m2 + mfaab;
			   mfaaa = m0;
			   m0 += c1o36 * oMdrho;
			   mfaab = m1 - m0 * vvz;
			   mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaba + mfabc;
			   m1 = mfabc - mfaba;
			   m0 = m2 + mfabb;
			   mfaba = m0;
			   m0 += c1o9 * oMdrho;
			   mfabb = m1 - m0 * vvz;
			   mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaca + mfacc;
			   m1 = mfacc - mfaca;
			   m0 = m2 + mfacb;
			   mfaca = m0;
			   m0 += c1o36 * oMdrho;
			   mfacb = m1 - m0 * vvz;
			   mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbaa + mfbac;
			   m1 = mfbac - mfbaa;
			   m0 = m2 + mfbab;
			   mfbaa = m0;
			   m0 += c1o9 * oMdrho;
			   mfbab = m1 - m0 * vvz;
			   mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbba + mfbbc;
			   m1 = mfbbc - mfbba;
			   m0 = m2 + mfbbb;
			   mfbba = m0;
			   m0 += c4o9 * oMdrho;
			   mfbbb = m1 - m0 * vvz;
			   mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbca + mfbcc;
			   m1 = mfbcc - mfbca;
			   m0 = m2 + mfbcb;
			   mfbca = m0;
			   m0 += c1o9 * oMdrho;
			   mfbcb = m1 - m0 * vvz;
			   mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcaa + mfcac;
			   m1 = mfcac - mfcaa;
			   m0 = m2 + mfcab;
			   mfcaa = m0;
			   m0 += c1o36 * oMdrho;
			   mfcab = m1 - m0 * vvz;
			   mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcba + mfcbc;
			   m1 = mfcbc - mfcba;
			   m0 = m2 + mfcbb;
			   mfcba = m0;
			   m0 += c1o9 * oMdrho;
			   mfcbb = m1 - m0 * vvz;
			   mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcca + mfccc;
			   m1 = mfccc - mfcca;
			   m0 = m2 + mfccb;
			   mfcca = m0;
			   m0 += c1o36 * oMdrho;
			   mfccb = m1 - m0 * vvz;
			   mfccc = m2 - 2. * m1 * vvz + vz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Y - Dir
			   m2 = mfaaa + mfaca;
			   m1 = mfaca - mfaaa;
			   m0 = m2 + mfaba;
			   mfaaa = m0;
			   m0 += c1o6 * oMdrho;
			   mfaba = m1 - m0 * vvy;
			   mfaca = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaab + mfacb;
			   m1 = mfacb - mfaab;
			   m0 = m2 + mfabb;
			   mfaab = m0;
			   mfabb = m1 - m0 * vvy;
			   mfacb = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaac + mfacc;
			   m1 = mfacc - mfaac;
			   m0 = m2 + mfabc;
			   mfaac = m0;
			   m0 += c1o18 * oMdrho;
			   mfabc = m1 - m0 * vvy;
			   mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbaa + mfbca;
			   m1 = mfbca - mfbaa;
			   m0 = m2 + mfbba;
			   mfbaa = m0;
			   m0 += c2o3 * oMdrho;
			   mfbba = m1 - m0 * vvy;
			   mfbca = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbab + mfbcb;
			   m1 = mfbcb - mfbab;
			   m0 = m2 + mfbbb;
			   mfbab = m0;
			   mfbbb = m1 - m0 * vvy;
			   mfbcb = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbac + mfbcc;
			   m1 = mfbcc - mfbac;
			   m0 = m2 + mfbbc;
			   mfbac = m0;
			   m0 += c2o9 * oMdrho;
			   mfbbc = m1 - m0 * vvy;
			   mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcaa + mfcca;
			   m1 = mfcca - mfcaa;
			   m0 = m2 + mfcba;
			   mfcaa = m0;
			   m0 += c1o6 * oMdrho;
			   mfcba = m1 - m0 * vvy;
			   mfcca = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcab + mfccb;
			   m1 = mfccb - mfcab;
			   m0 = m2 + mfcbb;
			   mfcab = m0;
			   mfcbb = m1 - m0 * vvy;
			   mfccb = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcac + mfccc;
			   m1 = mfccc - mfcac;
			   m0 = m2 + mfcbc;
			   mfcac = m0;
			   m0 += c1o18 * oMdrho;
			   mfcbc = m1 - m0 * vvy;
			   mfccc = m2 - 2. * m1 * vvy + vy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   m2 = mfaaa + mfcaa;
			   m1 = mfcaa - mfaaa;
			   m0 = m2 + mfbaa;
			   mfaaa = m0;
			   m0 += 1. * oMdrho;
			   mfbaa = m1 - m0 * vvx;
			   mfcaa = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaba + mfcba;
			   m1 = mfcba - mfaba;
			   m0 = m2 + mfbba;
			   mfaba = m0;
			   mfbba = m1 - m0 * vvx;
			   mfcba = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaca + mfcca;
			   m1 = mfcca - mfaca;
			   m0 = m2 + mfbca;
			   mfaca = m0;
			   m0 += c1o3 * oMdrho;
			   mfbca = m1 - m0 * vvx;
			   mfcca = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaab + mfcab;
			   m1 = mfcab - mfaab;
			   m0 = m2 + mfbab;
			   mfaab = m0;
			   mfbab = m1 - m0 * vvx;
			   mfcab = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfabb + mfcbb;
			   m1 = mfcbb - mfabb;
			   m0 = m2 + mfbbb;
			   mfabb = m0;
			   mfbbb = m1 - m0 * vvx;
			   mfcbb = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfacb + mfccb;
			   m1 = mfccb - mfacb;
			   m0 = m2 + mfbcb;
			   mfacb = m0;
			   mfbcb = m1 - m0 * vvx;
			   mfccb = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaac + mfcac;
			   m1 = mfcac - mfaac;
			   m0 = m2 + mfbac;
			   mfaac = m0;
			   m0 += c1o3 * oMdrho;
			   mfbac = m1 - m0 * vvx;
			   mfcac = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfabc + mfcbc;
			   m1 = mfcbc - mfabc;
			   m0 = m2 + mfbbc;
			   mfabc = m0;
			   mfbbc = m1 - m0 * vvx;
			   mfcbc = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfacc + mfccc;
			   m1 = mfccc - mfacc;
			   m0 = m2 + mfbcc;
			   mfacc = m0;
			   m0 += c1o9 * oMdrho;
			   mfbcc = m1 - m0 * vvx;
			   mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Cumulants
			   ////////////////////////////////////////////////////////////////////////////////////
			   LBMReal OxxPyyPzz = 1.; //omega2 or bulk viscosity
			   //LBMReal OxyyPxzz = 2.0 - collFactorM;// 1.;//-s9;//2+s9;//
			   //LBMReal OxyyMxzz  = 2.0 - collFactorM;// 1.;//2+s9;//
			   LBMReal O4 = 1.0;//collFactorM;// 1.;
			   LBMReal O5 = 1.;
			   LBMReal O6 = 1.;


			   /////fourth order parameters; here only for test. Move out of loop!

			   LBMReal OxyyPxzz =  8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
			   LBMReal OxyyMxzz =  8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
			   LBMReal Oxyz =  24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 * collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0));
			   LBMReal A =  (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
			   //FIXME:  warning C4459: declaration of 'B' hides global declaration (message : see declaration of 'D3Q27System::DIR_00M' )
			   LBMReal BB =   (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));


			   //Cum 4.
			   //LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
			   //LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
			   //LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

			   LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
			   LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
			   LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

			   LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho);
			   LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho);
			   LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho);

			   //Cum 5.
			   LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			   LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			   LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			   //Cum 6.
			   LBMReal CUMccc = mfccc + ((-4. * mfbbb * mfbbb
				   - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				   - 4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
				   - 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				   + (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					   + 2. * (mfcaa * mfaca * mfaac)
					   + 16. * mfbba * mfbab * mfabb)
				   - c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho
				   - c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) - c1o27 * oMdrho * oMdrho * (-2. * oMdrho)
				   + (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					   + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 * oMdrho) + c1o27 * oMdrho;

			   //2.
			   // linear combinations
			   LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
				mxxPyyPzz-=mfaaa;//12.03.21 shifted by mfaaa
			   LBMReal mxxMyy = mfcaa - mfaca;
			   LBMReal mxxMzz = mfcaa - mfaac;

			   //applying phase field gradients first part:
			  // mxxPyyPzz += c2o3 * rhoToPhi * (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz);

			   //17.03.2021 attempt for statililization by assymptotically vanishing bias
			   //LBMReal correctionScaling = rhoToPhi /rho;// +0.5;// (vx2 + vy2 + vz2) * 100;// +0.5;//(vx2 + vy2 + vz2)*1000;
			   //mxxPyyPzz += (1.0/3.0) * (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz)* correctionScaling; // As in Hesam's code
			   //mxxMyy += c1o3 * (dX1_phi * vvx - dX2_phi * vvy)* correctionScaling;
			   //mxxMzz += c1o3 * (dX1_phi * vvx - dX3_phi * vvz) * correctionScaling;
			   //mfabb += c1o6 * (dX2_phi * vvz + dX3_phi * vvy) * correctionScaling;
			   //mfbab += c1o6 * (dX1_phi * vvz + dX3_phi * vvx) * correctionScaling;
			   //mfbba += c1o6 * (dX1_phi * vvy + dX2_phi * vvx) * correctionScaling;


			   //14.04.2021 filtered velocity

			   //LBMReal correctionScaling =  rhoToPhi / rho;// +0.5;// (vx2 + vy2 + vz2) * 100;// +0.5;//(vx2 + vy2 + vz2)*1000;
			   //mxxPyyPzz += (1.0 / 3.0) * (dX1_phi * vvxF + dX2_phi * vvyF + dX3_phi * vvzF) * correctionScaling; // As in Hesam's code
			   //mxxMyy += c1o3 * (dX1_phi * vvxF - dX2_phi * vvyF) * correctionScaling;
			   //mxxMzz += c1o3 * (dX1_phi * vvxF - dX3_phi * vvzF) * correctionScaling;
			   //mfabb += c1o6 * (dX2_phi * vvzF + dX3_phi * vvyF) * correctionScaling;
			   //mfbab += c1o6 * (dX1_phi * vvzF + dX3_phi * vvxF) * correctionScaling;
			   //mfbba += c1o6 * (dX1_phi * vvyF + dX2_phi * vvxF) * correctionScaling;


			   LBMReal dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz);
			   LBMReal dyuy =  dxux + collFactorM * c3o2 * mxxMyy;
			   LBMReal dzuz =  dxux + collFactorM * c3o2 * mxxMzz;

			   LBMReal Dxy = -three * collFactorM * mfbba;
			   LBMReal Dxz = -three * collFactorM * mfbab;
			   LBMReal Dyz = -three * collFactorM * mfabb;

			   ////relax unfiltered
			   //! divergenceFilter 10.05.2021
			   LBMReal divMag= (1.0 - phi[DIR_000]) * (phi[DIR_000])*10*5*sqrt(fabs((OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))));
			  // LBMReal divMag = 500 *500* 50*(fabs((OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))))* (fabs((OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))));
			   //LBMReal divMag = (dX1_phi * dxux) > 0 ? (dX1_phi * dxux) : 0;
			   //divMag += (dX2_phi * dyuy) > 0 ? (dX2_phi * dyuy) : 0;
			   //divMag += (dX3_phi * dzuz) > 0 ? (dX3_phi * dzuz) : 0;
			   //divMag *= 5000;
			   //divMag+= denom * 10 * 5 * sqrt(fabs((OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz))));
			   //LBMReal divMag = 5000 * (fabs(dX1_phi * dxux)+fabs(dX2_phi * dyuy)+fabs(dX3_phi * dzuz));
			   collFactorM = collFactorM / (1.0 + 3.0 * divMag);

			   collFactorM = (collFactorM > 1.0) ? collFactorM : 1.0;


			   mxxPyyPzz += OxxPyyPzz * (/*mfaaa*/ - mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
			   mxxMyy += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vy2 * dyuy);
			   mxxMzz += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vz2 * dzuz);

			   mfabb += collFactorM * (-mfabb);
			   mfbab += collFactorM * (-mfbab);
			   mfbba += collFactorM * (-mfbba);


			   //relax filtered
			   //LBMReal interfaceFilter=0.001;
			   //LBMReal interfaceFactor = c1;// (dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);

			   //mxxPyyPzz += OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
			   //
			   //wadjust = collFactorM + (1. - collFactorM) * fabs(mxxMyy) / (fabs(mxxMyy) * interfaceFactor + interfaceFilter)* interfaceFactor;
			   //mxxMyy += wadjust * (-mxxMyy);// -3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vy2 * dyuy);
			   //wadjust = collFactorM + (1. - collFactorM) * fabs(mxxMzz) / (fabs(mxxMzz) * interfaceFactor + interfaceFilter) * interfaceFactor;
			   //mxxMzz += wadjust * (-mxxMzz);// -3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vz2 * dzuz);

			   //wadjust = collFactorM + (1. - collFactorM) * fabs(mfabb) / (fabs(mfabb) * interfaceFactor + interfaceFilter) * interfaceFactor;
			   //mfabb += wadjust * (-mfabb);
			   //wadjust = collFactorM + (1. - collFactorM) * fabs(mfbab) / (fabs(mfbab) * interfaceFactor + interfaceFilter) * interfaceFactor;
			   //mfbab += wadjust * (-mfbab);
			   //wadjust = collFactorM + (1. - collFactorM) * fabs(mfbba) / (fabs(mfbba) * interfaceFactor + interfaceFilter) * interfaceFactor;
			   //mfbba += wadjust * (-mfbba);

			   //applying phase field gradients second part:
			   //mxxPyyPzz += c2o3 * rhoToPhi * (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz);
			   //mxxPyyPzz += (1.0 / 3.0) * (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling; // As in Hesam's code
			   //mxxMyy += c1o3 * (dX1_phi * vvx - dX2_phi * vvy) * correctionScaling;
			   //mxxMzz += c1o3 * (dX1_phi * vvx - dX3_phi * vvz) * correctionScaling;
			   //mfabb += c1o6 * (dX2_phi * vvz + dX3_phi * vvy) * correctionScaling;
			   //mfbab += c1o6 * (dX1_phi * vvz + dX3_phi * vvx) * correctionScaling;
			   //mfbba += c1o6 * (dX1_phi * vvy + dX2_phi * vvx) * correctionScaling;


			   //////updated pressure
			   //mfaaa += (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling;


			   //mxxPyyPzz += (1.0 / 3.0) * (dX1_phi * vvxF + dX2_phi * vvyF + dX3_phi * vvzF) * correctionScaling; // As in Hesam's code
			   //mxxMyy += c1o3 * (dX1_phi * vvxF - dX2_phi * vvyF) * correctionScaling;
			   //mxxMzz += c1o3 * (dX1_phi * vvxF - dX3_phi * vvzF) * correctionScaling;
			   //mfabb += c1o6 * (dX2_phi * vvzF + dX3_phi * vvyF) * correctionScaling;
			   //mfbab += c1o6 * (dX1_phi * vvzF + dX3_phi * vvxF) * correctionScaling;
			   //mfbba += c1o6 * (dX1_phi * vvyF + dX2_phi * vvxF) * correctionScaling;


			   //////updated pressure
			   //mfaaa += (dX1_phi * vvxF + dX2_phi * vvyF + dX3_phi * vvzF) * correctionScaling;


			   mxxPyyPzz += mfaaa;//12.03.21 shifted by mfaaa
			 //  mxxPyyPzz = mfaaa; //12.03.21 reguarized pressure !?
			   // linear combinations back
			   mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
			   mfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
			   mfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);

			   //3.
			   // linear combinations
			   LBMReal mxxyPyzz = mfcba + mfabc;
			   LBMReal mxxyMyzz = mfcba - mfabc;

			   LBMReal mxxzPyyz = mfcab + mfacb;
			   LBMReal mxxzMyyz = mfcab - mfacb;

			   LBMReal mxyyPxzz = mfbca + mfbac;
			   LBMReal mxyyMxzz = mfbca - mfbac;

			   //relax
			   wadjust = Oxyz + (1. - Oxyz) * fabs(mfbbb) / (fabs(mfbbb) + qudricLimit);
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

			   //4.
			   //CUMacc += O4 * (-CUMacc);
			   //CUMcac += O4 * (-CUMcac);
			   //CUMcca += O4 * (-CUMcca);

			   //CUMbbc += O4 * (-CUMbbc);
			   //CUMbcb += O4 * (-CUMbcb);
			   //CUMcbb += O4 * (-CUMcbb);


			   CUMacc = -O4 * (one / collFactorM - c1o2) * (dyuy + dzuz) * c2o3 * A + (one - O4) * (CUMacc);
			   CUMcac = -O4 * (one / collFactorM - c1o2) * (dxux + dzuz) * c2o3 * A + (one - O4) * (CUMcac);
			   CUMcca = -O4 * (one / collFactorM - c1o2) * (dyuy + dxux) * c2o3 * A + (one - O4) * (CUMcca);
			   CUMbbc = -O4 * (one / collFactorM - c1o2) * Dxy * c1o3 * BB + (one - O4) * (CUMbbc);
			   CUMbcb = -O4 * (one / collFactorM - c1o2) * Dxz * c1o3 * BB + (one - O4) * (CUMbcb);
			   CUMcbb = -O4 * (one / collFactorM - c1o2) * Dyz * c1o3 * BB + (one - O4) * (CUMcbb);




			   //CUMacc -= (one / collFactorM - c1o2) * (dyuy + dzuz) * c2o3 * A ;
			   //CUMcac -= (one / collFactorM - c1o2) * (dxux + dzuz) * c2o3 * A ;
			   //CUMcca -= (one / collFactorM - c1o2) * (dyuy + dxux) * c2o3 * A ;
			   //CUMbbc -= (one / collFactorM - c1o2) * Dxy * c1o3 * B ;
			   //CUMbcb -= (one / collFactorM - c1o2) * Dxz * c1o3 * B ;
			   //CUMcbb -= (one / collFactorM - c1o2) * Dyz * c1o3 * B ;

			   //wadjust = O4 + (1. - O4) * fabs(CUMacc) / (fabs(CUMacc) + qudricLimit);
			   //CUMacc += wadjust * (-CUMacc);
			   //wadjust = O4 + (1. - O4) * fabs(CUMcac) / (fabs(CUMcac) + qudricLimit);
			   //CUMcac += wadjust * (-CUMcac);
			   //wadjust = O4 + (1. - O4) * fabs(CUMcca) / (fabs(CUMcca) + qudricLimit);
			   //CUMcca += wadjust * (-CUMcca);
			   //wadjust = O4 + (1. - O4) * fabs(CUMbbc) / (fabs(CUMbbc) + qudricLimit);
			   //CUMbbc += wadjust * (-CUMbbc);
			   //wadjust = O4 + (1. - O4) * fabs(CUMbcb) / (fabs(CUMbcb) + qudricLimit);
			   //CUMbcb += wadjust * (-CUMbcb);
			   //wadjust = O4 + (1. - O4) * fabs(CUMcbb) / (fabs(CUMcbb) + qudricLimit);
			   //CUMcbb += wadjust * (-CUMcbb);






			   //CUMacc += (one / collFactorM - c1o2) * (dyuy + dzuz) * c2o3 * A;
			   //CUMcac += (one / collFactorM - c1o2) * (dxux + dzuz) * c2o3 * A;
			   //CUMcca += (one / collFactorM - c1o2) * (dyuy + dxux) * c2o3 * A;
			   //CUMbbc += (one / collFactorM - c1o2) * Dxy * c1o3 * B;
			   //CUMbcb += (one / collFactorM - c1o2) * Dxz * c1o3 * B;
			   //CUMcbb += (one / collFactorM - c1o2) * Dyz * c1o3 * B;

			   //5.
			   CUMbcc += O5 * (-CUMbcc);
			   CUMcbc += O5 * (-CUMcbc);
			   CUMccb += O5 * (-CUMccb);


			   //wadjust = O5 + (1. - O5) * fabs(CUMbcc) / (fabs(CUMbcc) + qudricLimit);
			   //CUMbcc += wadjust * (-CUMbcc);
			   //wadjust = O5 + (1. - O5) * fabs(CUMcbc) / (fabs(CUMcbc) + qudricLimit);
			   //CUMbcc += wadjust * (-CUMcbc);
			   //wadjust = O5 + (1. - O5) * fabs(CUMccb) / (fabs(CUMccb) + qudricLimit);
			   //CUMbcc += wadjust * (-CUMccb);


			   //6.
			   CUMccc += O6 * (-CUMccc);

			   //back cumulants to central moments
			   //4.
			   //mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
			   //mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
			   //mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

			   mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
			   mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
			   mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

			   mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho;
			   mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho;
			   mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - c1) * oMdrho;

			   //5.
			   mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
			   mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
			   mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

			   //6.
			   mfccc = CUMccc - ((-4. * mfbbb * mfbbb
				   - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
				   - 4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
				   - 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
				   + (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
					   + 2. * (mfcaa * mfaca * mfaac)
					   + 16. * mfbba * mfbab * mfabb)
				   - c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho
				   - c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) - c1o27 * oMdrho * oMdrho * (-2. * oMdrho)
				   + (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
					   + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 * oMdrho) - c1o27 * oMdrho;


			   ////////


			   ////////////////////////////////////////////////////////////////////////////////////
			   //forcing
			   mfbaa = -mfbaa;
			   mfaba = -mfaba;
			   mfaab = -mfaab;
			   //////////////////////////////////////////////////////////////////////////////////////

			   ////////////////////////////////////////////////////////////////////////////////////
			   //back
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + 1. * oMdrho) * (vz2 - vvz) * c1o2;
			   m1 = -mfaac - 2. * mfaab * vvz + mfaaa * (1. - vz2) - 1. * oMdrho * vz2;
			   m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + 1. * oMdrho) * (vz2 + vvz) * c1o2;
			   mfaaa = m0;
			   mfaab = m1;
			   mfaac = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
			   m1 = -mfabc - 2. * mfabb * vvz + mfaba * (1. - vz2);
			   m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
			   mfaba = m0;
			   mfabb = m1;
			   mfabc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			   m1 = -mfacc - 2. * mfacb * vvz + mfaca * (1. - vz2) - c1o3 * oMdrho * vz2;
			   m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			   mfaca = m0;
			   mfacb = m1;
			   mfacc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
			   m1 = -mfbac - 2. * mfbab * vvz + mfbaa * (1. - vz2);
			   m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
			   mfbaa = m0;
			   mfbab = m1;
			   mfbac = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
			   m1 = -mfbbc - 2. * mfbbb * vvz + mfbba * (1. - vz2);
			   m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
			   mfbba = m0;
			   mfbbb = m1;
			   mfbbc = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
			   m1 = -mfbcc - 2. * mfbcb * vvz + mfbca * (1. - vz2);
			   m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
			   mfbca = m0;
			   mfbcb = m1;
			   mfbcc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
			   m1 = -mfcac - 2. * mfcab * vvz + mfcaa * (1. - vz2) - c1o3 * oMdrho * vz2;
			   m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
			   mfcaa = m0;
			   mfcab = m1;
			   mfcac = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
			   m1 = -mfcbc - 2. * mfcbb * vvz + mfcba * (1. - vz2);
			   m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
			   mfcba = m0;
			   mfcbb = m1;
			   mfcbc = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
			   m1 = -mfccc - 2. * mfccb * vvz + mfcca * (1. - vz2) - c1o9 * oMdrho * vz2;
			   m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
			   mfcca = m0;
			   mfccb = m1;
			   mfccc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Y - Dir
			   m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfaca - 2. * mfaba * vvy + mfaaa * (1. - vy2) - c1o6 * oMdrho * vy2;
			   m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfaaa = m0;
			   mfaba = m1;
			   mfaca = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfacb - 2. * mfabb * vvy + mfaab * (1. - vy2) - c2o3 * oMdrho * vy2;
			   m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfaab = m0;
			   mfabb = m1;
			   mfacb = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfacc - 2. * mfabc * vvy + mfaac * (1. - vy2) - c1o6 * oMdrho * vy2;
			   m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfaac = m0;
			   mfabc = m1;
			   mfacc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
			   m1 = -mfbca - 2. * mfbba * vvy + mfbaa * (1. - vy2);
			   m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
			   mfbaa = m0;
			   mfbba = m1;
			   mfbca = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
			   m1 = -mfbcb - 2. * mfbbb * vvy + mfbab * (1. - vy2);
			   m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
			   mfbab = m0;
			   mfbbb = m1;
			   mfbcb = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
			   m1 = -mfbcc - 2. * mfbbc * vvy + mfbac * (1. - vy2);
			   m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
			   mfbac = m0;
			   mfbbc = m1;
			   mfbcc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfcca - 2. * mfcba * vvy + mfcaa * (1. - vy2) - c1o18 * oMdrho * vy2;
			   m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfcaa = m0;
			   mfcba = m1;
			   mfcca = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfccb - 2. * mfcbb * vvy + mfcab * (1. - vy2) - c2o9 * oMdrho * vy2;
			   m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfcab = m0;
			   mfcbb = m1;
			   mfccb = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
			   m1 = -mfccc - 2. * mfcbc * vvy + mfcac * (1. - vy2) - c1o18 * oMdrho * vy2;
			   m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
			   mfcac = m0;
			   mfcbc = m1;
			   mfccc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcaa - 2. * mfbaa * vvx + mfaaa * (1. - vx2) - c1o36 * oMdrho * vx2;
			   m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfaaa = m0;
			   mfbaa = m1;
			   mfcaa = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcba - 2. * mfbba * vvx + mfaba * (1. - vx2) - c1o9 * oMdrho * vx2;
			   m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfaba = m0;
			   mfbba = m1;
			   mfcba = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcca - 2. * mfbca * vvx + mfaca * (1. - vx2) - c1o36 * oMdrho * vx2;
			   m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfaca = m0;
			   mfbca = m1;
			   mfcca = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcab - 2. * mfbab * vvx + mfaab * (1. - vx2) - c1o9 * oMdrho * vx2;
			   m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfaab = m0;
			   mfbab = m1;
			   mfcab = m2;
			   ///////////b////////////////////////////////////////////////////////////////////////
			   m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcbb - 2. * mfbbb * vvx + mfabb * (1. - vx2) - c4o9 * oMdrho * vx2;
			   m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfabb = m0;
			   mfbbb = m1;
			   mfcbb = m2;
			   ///////////b////////////////////////////////////////////////////////////////////////
			   m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfccb - 2. * mfbcb * vvx + mfacb * (1. - vx2) - c1o9 * oMdrho * vx2;
			   m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfacb = m0;
			   mfbcb = m1;
			   mfccb = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcac - 2. * mfbac * vvx + mfaac * (1. - vx2) - c1o36 * oMdrho * vx2;
			   m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfaac = m0;
			   mfbac = m1;
			   mfcac = m2;
			   ///////////c////////////////////////////////////////////////////////////////////////
			   m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfcbc - 2. * mfbbc * vvx + mfabc * (1. - vx2) - c1o9 * oMdrho * vx2;
			   m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfabc = m0;
			   mfbbc = m1;
			   mfcbc = m2;
			   ///////////c////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
			   m1 = -mfccc - 2. * mfbcc * vvx + mfacc * (1. - vx2) - c1o36 * oMdrho * vx2;
			   m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
			   mfacc = m0;
			   mfbcc = m1;
			   mfccc = m2;

			   //////////////////////////////////////////////////////////////////////////
			   //proof correctness
			   //////////////////////////////////////////////////////////////////////////
//#ifdef  PROOF_CORRECTNESS
//			   LBMReal rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
//				   + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
//				   + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
//			   //LBMReal dif = fabs(drho - rho_post);
//			   LBMReal dif = drho+  (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) *correctionScaling - rho_post;
//#ifdef SINGLEPRECISION
//			   if (dif > 10.0E-7 || dif < -10.0E-7)
//#else
//			   if (dif > 10.0E-15 || dif < -10.0E-15)
//#endif
//			   {
//				   UB_THROW(UbException(UB_EXARGS, "drho=" + UbSystem::toString(drho) + ", rho_post=" + UbSystem::toString(rho_post)
//					   + " dif=" + UbSystem::toString(dif)
//					   + " drho is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
//				   //UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): drho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
//				   //exit(EXIT_FAILURE);
//			   }
//#endif
			   //////////////////////////////////////////////////////////////////////////
			   //write distribution
			   //////////////////////////////////////////////////////////////////////////

			   /////classical source term 8.4.2021

			   mfcbb += 3.0 * (0.5 * forcingTerm[DIR_P00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[E  ];
			   mfbcb += 3.0 * (0.5 * forcingTerm[DIR_0P0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[N  ];
			   mfbbc += 3.0 * (0.5 * forcingTerm[DIR_00P]) / rho;    //-(3.0*p1 - rho)*WEIGTH[T  ];
			   mfccb += 3.0 * (0.5 * forcingTerm[DIR_PP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NE ];
			   mfacb += 3.0 * (0.5 * forcingTerm[DIR_MP0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NW ];
			   mfcbc += 3.0 * (0.5 * forcingTerm[DIR_P0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TE ];
			   mfabc += 3.0 * (0.5 * forcingTerm[DIR_M0P]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TW ];
			   mfbcc += 3.0 * (0.5 * forcingTerm[DIR_0PP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TN ];
			   mfbac += 3.0 * (0.5 * forcingTerm[DIR_0MP]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TS ];
			   mfccc += 3.0 * (0.5 * forcingTerm[DIR_PPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNE];
			   mfacc += 3.0 * (0.5 * forcingTerm[DIR_MPP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNW];
			   mfcac += 3.0 * (0.5 * forcingTerm[DIR_PMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSE];
			   mfaac += 3.0 * (0.5 * forcingTerm[DIR_MMP]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSW];
			   mfabb += 3.0 * (0.5 * forcingTerm[DIR_M00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[W  ];
			   mfbab += 3.0 * (0.5 * forcingTerm[DIR_0M0]) / rho;    //-(3.0*p1 - rho)*WEIGTH[S  ];
			   mfbba += 3.0 * (0.5 * forcingTerm[DIR_00M]) / rho;    //-(3.0*p1 - rho)*WEIGTH[B  ];
			   mfaab += 3.0 * (0.5 * forcingTerm[DIR_MM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SW ];
			   mfcab += 3.0 * (0.5 * forcingTerm[DIR_PM0]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SE ];
			   mfaba += 3.0 * (0.5 * forcingTerm[DIR_M0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BW ];
			   mfcba += 3.0 * (0.5 * forcingTerm[DIR_P0M]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BE ];
			   mfbaa += 3.0 * (0.5 * forcingTerm[DIR_0MM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BS ];
			   mfbca += 3.0 * (0.5 * forcingTerm[DIR_0PM]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BN ];
			   mfaaa += 3.0 * (0.5 * forcingTerm[DIR_MMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSW];
			   mfcaa += 3.0 * (0.5 * forcingTerm[DIR_PMM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSE];
			   mfaca += 3.0 * (0.5 * forcingTerm[DIR_MPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNW];
			   mfcca += 3.0 * (0.5 * forcingTerm[DIR_PPM]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNE];
			   mfbbb += 3.0 * (0.5 * forcingTerm[DIR_000]) / rho; //- (3.0*p1 - rho)*WEIGTH[REST]



			   ////////////////////


			   (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = mfabb * rho*c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = mfbab * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = mfbba * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = mfaab * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = mfcab * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = mfaba * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = mfcba * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = mfbaa * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = mfbca * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = mfaaa * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = mfcaa * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = mfaca * rho * c1o3;
			   (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca * rho * c1o3;

			   (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac * rho * c1o3;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac * rho * c1o3;

			   (*this->zeroDistributionsF)(x1, x2, x3) = mfbbb * rho * c1o3;
			   //////////////////////////////////////////////////////////////////////////

			   ////!Incompressible Kernal

                            
//                            ///////Old Kernel \|/
//                            // ux += forcingX1*deltaT*0.5; // X
//                            // uy += forcingX2*deltaT*0.5; // Y
//                            // uz += forcingX3*deltaT*0.5; // Z
//                        }
//
//                        LBMReal ux = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
//                                      (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
//                                      (mfcbb - mfabb)) /
//                                         (rho * c1o3) +
//                                     (mu * dX1_phi + forcingX1) / (2 * rho);
//
//                        LBMReal uy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
//                                      (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
//                                      (mfbcb - mfbab)) /
//                                         (rho * c1o3) +
//                                     (mu * dX2_phi + forcingX2) / (2 * rho);
//
//                        LBMReal uz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
//                                      (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
//                                      (mfbbc - mfbba)) /
//                                         (rho * c1o3) +
//                                     (mu * dX3_phi + forcingX3) / (2 * rho);
//
//                        //--------------------------------------------------------
//
//                        LBMReal ux2 = ux * ux;
//                        LBMReal uy2 = uy * uy;
//                        LBMReal uz2 = uz * uz;
//
//                        //----------- Calculating Forcing Terms * -------------
//                        for (int dir = STARTF; dir <= (FENDDIR); dir++) {
//                            LBMReal velProd = DX1[dir] * ux + DX2[dir] * uy + DX3[dir] * uz;
//                            LBMReal velSq1  = velProd * velProd;
//                            LBMReal gamma = WEIGTH[dir] * (1.0 + 3 * velProd + 4.5 * velSq1 - 1.5 * (ux2 + uy2 + uz2));
//
//                            LBMReal fac1 = (gamma - WEIGTH[dir]) * c1o3 * rhoToPhi;
//
//                            forcingTerm[dir] = ((-ux) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
//                                                (-uy) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
//                                                (-uz) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3))) +
//                                               (DX1[dir]) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
//                                               (DX2[dir]) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
//                                               (DX3[dir]) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3));
//                        }
//
//                        LBMReal gamma = WEIGTH[REST] * (1.0 - 1.5 * (ux2 + uy2 + uz2));
//                        LBMReal fac1      = (gamma - WEIGTH[REST]) * c1o3 * rhoToPhi;
//                        forcingTerm[REST] = (-ux) * (fac1 * dX1_phi + gamma * (mu * dX1_phi + forcingX1)) +
//                                            (-uy) * (fac1 * dX2_phi + gamma * (mu * dX2_phi + forcingX2)) +
//                                            (-uz) * (fac1 * dX3_phi + gamma * (mu * dX3_phi + forcingX3));
//
//                        //--------------------------------------------------------
//
//                        mfcbb = 3.0 * (mfcbb + 0.5 * forcingTerm[DIR_P00]) / rho;    //-(3.0*p1 - rho)*WEIGTH[E  ];
//                        mfbcb = 3.0 * (mfbcb + 0.5 * forcingTerm[N]) / rho;    //-(3.0*p1 - rho)*WEIGTH[N  ];
//                        mfbbc = 3.0 * (mfbbc + 0.5 * forcingTerm[T]) / rho;    //-(3.0*p1 - rho)*WEIGTH[T  ];
//                        mfccb = 3.0 * (mfccb + 0.5 * forcingTerm[NE]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NE ];
//                        mfacb = 3.0 * (mfacb + 0.5 * forcingTerm[NW]) / rho;   //-(3.0*p1 - rho)*WEIGTH[NW ];
//                        mfcbc = 3.0 * (mfcbc + 0.5 * forcingTerm[TE]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TE ];
//                        mfabc = 3.0 * (mfabc + 0.5 * forcingTerm[TW]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TW ];
//                        mfbcc = 3.0 * (mfbcc + 0.5 * forcingTerm[TN]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TN ];
//                        mfbac = 3.0 * (mfbac + 0.5 * forcingTerm[TS]) / rho;   //-(3.0*p1 - rho)*WEIGTH[TS ];
//                        mfccc = 3.0 * (mfccc + 0.5 * forcingTerm[TNE]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNE];
//                        mfacc = 3.0 * (mfacc + 0.5 * forcingTerm[TNW]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TNW];
//                        mfcac = 3.0 * (mfcac + 0.5 * forcingTerm[TSE]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSE];
//                        mfaac = 3.0 * (mfaac + 0.5 * forcingTerm[TSW]) / rho;  //-(3.0*p1 - rho)*WEIGTH[TSW];
//                        mfabb = 3.0 * (mfabb + 0.5 * forcingTerm[W]) / rho;    //-(3.0*p1 - rho)*WEIGTH[W  ];
//                        mfbab = 3.0 * (mfbab + 0.5 * forcingTerm[S]) / rho;    //-(3.0*p1 - rho)*WEIGTH[S  ];
//                        mfbba = 3.0 * (mfbba + 0.5 * forcingTerm[B]) / rho;    //-(3.0*p1 - rho)*WEIGTH[B  ];
//                        mfaab = 3.0 * (mfaab + 0.5 * forcingTerm[SW]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SW ];
//                        mfcab = 3.0 * (mfcab + 0.5 * forcingTerm[SE]) / rho;   //-(3.0*p1 - rho)*WEIGTH[SE ];
//                        mfaba = 3.0 * (mfaba + 0.5 * forcingTerm[BW]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BW ];
//                        mfcba = 3.0 * (mfcba + 0.5 * forcingTerm[BE]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BE ];
//                        mfbaa = 3.0 * (mfbaa + 0.5 * forcingTerm[BS]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BS ];
//                        mfbca = 3.0 * (mfbca + 0.5 * forcingTerm[BN]) / rho;   //-(3.0*p1 - rho)*WEIGTH[BN ];
//                        mfaaa = 3.0 * (mfaaa + 0.5 * forcingTerm[BSW]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSW];
//                        mfcaa = 3.0 * (mfcaa + 0.5 * forcingTerm[BSE]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BSE];
//                        mfaca = 3.0 * (mfaca + 0.5 * forcingTerm[BNW]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNW];
//                        mfcca = 3.0 * (mfcca + 0.5 * forcingTerm[BNE]) / rho;  //-(3.0*p1 - rho)*WEIGTH[BNE];
//                        mfbbb = 3.0 * (mfbbb + 0.5 * forcingTerm[REST]) / rho; //- (3.0*p1 - rho)*WEIGTH[REST];
//
//                        LBMReal rho1 = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
//                                       (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
//                                       (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) +
//                                       (mfbba + mfbbc) + mfbbb;
//
//
//                        LBMReal oMdrho, m0, m1, m2;
//
//                        oMdrho = mfccc + mfaaa;
//                        m0     = mfaca + mfcac;
//                        m1     = mfacc + mfcaa;
//                        m2     = mfaac + mfcca;
//                        oMdrho += m0;
//                        m1 += m2;
//                        oMdrho += m1;
//                        m0 = mfbac + mfbca;
//                        m1 = mfbaa + mfbcc;
//                        m0 += m1;
//                        m1 = mfabc + mfcba;
//                        m2 = mfaba + mfcbc;
//                        m1 += m2;
//                        m0 += m1;
//                        m1 = mfacb + mfcab;
//                        m2 = mfaab + mfccb;
//                        m1 += m2;
//                        m0 += m1;
//                        oMdrho += m0;
//                        m0 = mfabb + mfcbb;
//                        m1 = mfbab + mfbcb;
//                        m2 = mfbba + mfbbc;
//                        m0 += m1 + m2;
//                        m0 += mfbbb; // hat gefehlt
//                        oMdrho = 1. - (oMdrho + m0);
//                        // oMdrho = rho - (oMdrho + m0);
//
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        LBMReal wadjust;
//                        LBMReal qudricLimit = 0.01;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Hin
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Z - Dir
//                        m2    = mfaaa + mfaac;
//                        m1    = mfaac - mfaaa;
//                        m0    = m2 + mfaab;
//                        mfaaa = m0;
//                        m0 += c1o36 * oMdrho;
//                        mfaab = m1 - m0 * uz;
//                        mfaac = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaba + mfabc;
//                        m1    = mfabc - mfaba;
//                        m0    = m2 + mfabb;
//                        mfaba = m0;
//                        m0 += c1o9 * oMdrho;
//                        mfabb = m1 - m0 * uz;
//                        mfabc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaca + mfacc;
//                        m1    = mfacc - mfaca;
//                        m0    = m2 + mfacb;
//                        mfaca = m0;
//                        m0 += c1o36 * oMdrho;
//                        mfacb = m1 - m0 * uz;
//                        mfacc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbaa + mfbac;
//                        m1    = mfbac - mfbaa;
//                        m0    = m2 + mfbab;
//                        mfbaa = m0;
//                        m0 += c1o9 * oMdrho;
//                        mfbab = m1 - m0 * uz;
//                        mfbac = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbba + mfbbc;
//                        m1    = mfbbc - mfbba;
//                        m0    = m2 + mfbbb;
//                        mfbba = m0;
//                        m0 += c4o9 * oMdrho;
//                        mfbbb = m1 - m0 * uz;
//                        mfbbc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbca + mfbcc;
//                        m1    = mfbcc - mfbca;
//                        m0    = m2 + mfbcb;
//                        mfbca = m0;
//                        m0 += c1o9 * oMdrho;
//                        mfbcb = m1 - m0 * uz;
//                        mfbcc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcaa + mfcac;
//                        m1    = mfcac - mfcaa;
//                        m0    = m2 + mfcab;
//                        mfcaa = m0;
//                        m0 += c1o36 * oMdrho;
//                        mfcab = m1 - m0 * uz;
//                        mfcac = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcba + mfcbc;
//                        m1    = mfcbc - mfcba;
//                        m0    = m2 + mfcbb;
//                        mfcba = m0;
//                        m0 += c1o9 * oMdrho;
//                        mfcbb = m1 - m0 * uz;
//                        mfcbc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcca + mfccc;
//                        m1    = mfccc - mfcca;
//                        m0    = m2 + mfccb;
//                        mfcca = m0;
//                        m0 += c1o36 * oMdrho;
//                        mfccb = m1 - m0 * uz;
//                        mfccc = m2 - 2. * m1 * uz + uz2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Y - Dir
//                        m2    = mfaaa + mfaca;
//                        m1    = mfaca - mfaaa;
//                        m0    = m2 + mfaba;
//                        mfaaa = m0;
//                        m0 += c1o6 * oMdrho;
//                        mfaba = m1 - m0 * uy;
//                        mfaca = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaab + mfacb;
//                        m1    = mfacb - mfaab;
//                        m0    = m2 + mfabb;
//                        mfaab = m0;
//                        mfabb = m1 - m0 * uy;
//                        mfacb = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaac + mfacc;
//                        m1    = mfacc - mfaac;
//                        m0    = m2 + mfabc;
//                        mfaac = m0;
//                        m0 += c1o18 * oMdrho;
//                        mfabc = m1 - m0 * uy;
//                        mfacc = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbaa + mfbca;
//                        m1    = mfbca - mfbaa;
//                        m0    = m2 + mfbba;
//                        mfbaa = m0;
//                        m0 += c2o3 * oMdrho;
//                        mfbba = m1 - m0 * uy;
//                        mfbca = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbab + mfbcb;
//                        m1    = mfbcb - mfbab;
//                        m0    = m2 + mfbbb;
//                        mfbab = m0;
//                        mfbbb = m1 - m0 * uy;
//                        mfbcb = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfbac + mfbcc;
//                        m1    = mfbcc - mfbac;
//                        m0    = m2 + mfbbc;
//                        mfbac = m0;
//                        m0 += c2o9 * oMdrho;
//                        mfbbc = m1 - m0 * uy;
//                        mfbcc = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcaa + mfcca;
//                        m1    = mfcca - mfcaa;
//                        m0    = m2 + mfcba;
//                        mfcaa = m0;
//                        m0 += c1o6 * oMdrho;
//                        mfcba = m1 - m0 * uy;
//                        mfcca = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcab + mfccb;
//                        m1    = mfccb - mfcab;
//                        m0    = m2 + mfcbb;
//                        mfcab = m0;
//                        mfcbb = m1 - m0 * uy;
//                        mfccb = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfcac + mfccc;
//                        m1    = mfccc - mfcac;
//                        m0    = m2 + mfcbc;
//                        mfcac = m0;
//                        m0 += c1o18 * oMdrho;
//                        mfcbc = m1 - m0 * uy;
//                        mfccc = m2 - 2. * m1 * uy + uy2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // X - Dir
//                        m2    = mfaaa + mfcaa;
//                        m1    = mfcaa - mfaaa;
//                        m0    = m2 + mfbaa;
//                        mfaaa = m0;
//                        m0 += 1. * oMdrho;
//                        mfbaa = m1 - m0 * ux;
//                        mfcaa = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaba + mfcba;
//                        m1    = mfcba - mfaba;
//                        m0    = m2 + mfbba;
//                        mfaba = m0;
//                        mfbba = m1 - m0 * ux;
//                        mfcba = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaca + mfcca;
//                        m1    = mfcca - mfaca;
//                        m0    = m2 + mfbca;
//                        mfaca = m0;
//                        m0 += c1o3 * oMdrho;
//                        mfbca = m1 - m0 * ux;
//                        mfcca = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaab + mfcab;
//                        m1    = mfcab - mfaab;
//                        m0    = m2 + mfbab;
//                        mfaab = m0;
//                        mfbab = m1 - m0 * ux;
//                        mfcab = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfabb + mfcbb;
//                        m1    = mfcbb - mfabb;
//                        m0    = m2 + mfbbb;
//                        mfabb = m0;
//                        mfbbb = m1 - m0 * ux;
//                        mfcbb = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfacb + mfccb;
//                        m1    = mfccb - mfacb;
//                        m0    = m2 + mfbcb;
//                        mfacb = m0;
//                        mfbcb = m1 - m0 * ux;
//                        mfccb = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfaac + mfcac;
//                        m1    = mfcac - mfaac;
//                        m0    = m2 + mfbac;
//                        mfaac = m0;
//                        m0 += c1o3 * oMdrho;
//                        mfbac = m1 - m0 * ux;
//                        mfcac = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfabc + mfcbc;
//                        m1    = mfcbc - mfabc;
//                        m0    = m2 + mfbbc;
//                        mfabc = m0;
//                        mfbbc = m1 - m0 * ux;
//                        mfcbc = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m2    = mfacc + mfccc;
//                        m1    = mfccc - mfacc;
//                        m0    = m2 + mfbcc;
//                        mfacc = m0;
//                        m0 += c1o9 * oMdrho;
//                        mfbcc = m1 - m0 * ux;
//                        mfccc = m2 - 2. * m1 * ux + ux2 * m0;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Cumulants
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        LBMReal OxxPyyPzz = 1.; // omega2 or bulk viscosity
//                        LBMReal OxyyPxzz  = 1.; //-s9;//2+s9;//
//                        LBMReal OxyyMxzz  = 1.; // 2+s9;//
//                        LBMReal O4        = 1.;
//                        LBMReal O5        = 1.;
//                        LBMReal O6        = 1.;
//
//                        // Cum 4.
//                        LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
//                        LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
//                        LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
//
//                        LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) +
//                                                  c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);
//                        LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) +
//                                                  c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);
//                        LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) +
//                                                  c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho);
//
//                        // Cum 5.
//                        LBMReal CUMbcc = mfbcc -
//                                         (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb +
//                                          2. * (mfbab * mfacb + mfbba * mfabc)) -
//                                         c1o3 * (mfbca + mfbac) * oMdrho;
//                        LBMReal CUMcbc = mfcbc -
//                                         (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb +
//                                          2. * (mfabb * mfcab + mfbba * mfbac)) -
//                                         c1o3 * (mfcba + mfabc) * oMdrho;
//                        LBMReal CUMccb = mfccb -
//                                         (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb +
//                                          2. * (mfbab * mfbca + mfabb * mfcba)) -
//                                         c1o3 * (mfacb + mfcab) * oMdrho;
//
//                        // Cum 6.
//                        LBMReal CUMccc =
//                            mfccc +
//                            ((-4. * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
//                              4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc) -
//                              2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) +
//                             (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
//                              2. * (mfcaa * mfaca * mfaac) + 16. * mfbba * mfbab * mfabb) -
//                             c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho -
//                             c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) -
//                             c1o27 * oMdrho * oMdrho * (-2. * oMdrho) +
//                             (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
//                              (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) *
//                                 c2o3 * oMdrho) +
//                            c1o27 * oMdrho;
//
//                        // 2.
//                        // linear combinations
//                        LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
//                        LBMReal mxxMyy    = mfcaa - mfaca;
//                        LBMReal mxxMzz    = mfcaa - mfaac;
//
//                        LBMReal dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
//                        LBMReal dyuy = dxux + collFactorM * c3o2 * mxxMyy;
//                        LBMReal dzuz = dxux + collFactorM * c3o2 * mxxMzz;
//
//                        (*divU)(x1, x2, x3) = dxux + dyuy + dzuz;
//
//                        // relax
//                        mxxPyyPzz += OxxPyyPzz * (mfaaa - mxxPyyPzz) -
//                                     3. * (1. - c1o2 * OxxPyyPzz) * (ux2 * dxux + uy2 * dyuy + uz2 * dzuz);
//                        mxxMyy += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uy2 * dyuy);
//                        mxxMzz += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uz2 * dzuz);
//
//                        mfabb += collFactorM * (-mfabb);
//                        mfbab += collFactorM * (-mfbab);
//                        mfbba += collFactorM * (-mfbba);
//
//                        // linear combinations back
//                        mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
//                        mfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
//                        mfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);
//
//                        // 3.
//                        // linear combinations
//                        LBMReal mxxyPyzz = mfcba + mfabc;
//                        LBMReal mxxyMyzz = mfcba - mfabc;
//
//                        LBMReal mxxzPyyz = mfcab + mfacb;
//                        LBMReal mxxzMyyz = mfcab - mfacb;
//
//                        LBMReal mxyyPxzz = mfbca + mfbac;
//                        LBMReal mxyyMxzz = mfbca - mfbac;
//
//                        // relax
//                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mfbbb) / (fabs(mfbbb) + qudricLimit);
//                        mfbbb += wadjust * (-mfbbb);
//                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxyPyzz) / (fabs(mxxyPyzz) + qudricLimit);
//                        mxxyPyzz += wadjust * (-mxxyPyzz);
//                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxyMyzz) / (fabs(mxxyMyzz) + qudricLimit);
//                        mxxyMyzz += wadjust * (-mxxyMyzz);
//                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxzPyyz) / (fabs(mxxzPyyz) + qudricLimit);
//                        mxxzPyyz += wadjust * (-mxxzPyyz);
//                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxzMyyz) / (fabs(mxxzMyyz) + qudricLimit);
//                        mxxzMyyz += wadjust * (-mxxzMyyz);
//                        wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxyyPxzz) / (fabs(mxyyPxzz) + qudricLimit);
//                        mxyyPxzz += wadjust * (-mxyyPxzz);
//                        wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxyyMxzz) / (fabs(mxyyMxzz) + qudricLimit);
//                        mxyyMxzz += wadjust * (-mxyyMxzz);
//
//                        // linear combinations back
//                        mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
//                        mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
//                        mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
//                        mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
//                        mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
//                        mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
//
//                        // 4.
//                        CUMacc += O4 * (-CUMacc);
//                        CUMcac += O4 * (-CUMcac);
//                        CUMcca += O4 * (-CUMcca);
//
//                        CUMbbc += O4 * (-CUMbbc);
//                        CUMbcb += O4 * (-CUMbcb);
//                        CUMcbb += O4 * (-CUMcbb);
//
//                        // 5.
//                        CUMbcc += O5 * (-CUMbcc);
//                        CUMcbc += O5 * (-CUMcbc);
//                        CUMccb += O5 * (-CUMccb);
//
//                        // 6.
//                        CUMccc += O6 * (-CUMccc);
//
//                        // back cumulants to central moments
//                        // 4.
//                        mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
//                        mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
//                        mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
//
//                        mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho +
//                                c1o9 * (oMdrho - 1) * oMdrho;
//                        mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho +
//                                c1o9 * (oMdrho - 1) * oMdrho;
//                        mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho +
//                                c1o9 * (oMdrho - 1) * oMdrho;
//
//                        // 5.
//                        mfbcc = CUMbcc +
//                                (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb +
//                                 2. * (mfbab * mfacb + mfbba * mfabc)) +
//                                c1o3 * (mfbca + mfbac) * oMdrho;
//                        mfcbc = CUMcbc +
//                                (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb +
//                                 2. * (mfabb * mfcab + mfbba * mfbac)) +
//                                c1o3 * (mfcba + mfabc) * oMdrho;
//                        mfccb = CUMccb +
//                                (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb +
//                                 2. * (mfbab * mfbca + mfabb * mfcba)) +
//                                c1o3 * (mfacb + mfcab) * oMdrho;
//
//                        // 6.
//                        mfccc = CUMccc -
//                                ((-4. * mfbbb * mfbbb - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) -
//                                  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) -
//                                  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) +
//                                 (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) +
//                                  2. * (mfcaa * mfaca * mfaac) + 16. * mfbba * mfbab * mfabb) -
//                                 c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho -
//                                 c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) -
//                                 c1o27 * oMdrho * oMdrho * (-2. * oMdrho) +
//                                 (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) +
//                                  (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) *
//                                     c2o3 * oMdrho) -
//                                c1o27 * oMdrho;
//
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // forcing
//                        mfbaa = -mfbaa;
//                        mfaba = -mfaba;
//                        mfaab = -mfaab;
//                        //////////////////////////////////////////////////////////////////////////////////////
//
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // back
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Z - Dir
//                        m0    = mfaac * c1o2 + mfaab * (uz - c1o2) + (mfaaa + 1. * oMdrho) * (uz2 - uz) * c1o2;
//                        m1    = -mfaac - 2. * mfaab * uz + mfaaa * (1. - uz2) - 1. * oMdrho * uz2;
//                        m2    = mfaac * c1o2 + mfaab * (uz + c1o2) + (mfaaa + 1. * oMdrho) * (uz2 + uz) * c1o2;
//                        mfaaa = m0;
//                        mfaab = m1;
//                        mfaac = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfabc * c1o2 + mfabb * (uz - c1o2) + mfaba * (uz2 - uz) * c1o2;
//                        m1    = -mfabc - 2. * mfabb * uz + mfaba * (1. - uz2);
//                        m2    = mfabc * c1o2 + mfabb * (uz + c1o2) + mfaba * (uz2 + uz) * c1o2;
//                        mfaba = m0;
//                        mfabb = m1;
//                        mfabc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfacc * c1o2 + mfacb * (uz - c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
//                        m1    = -mfacc - 2. * mfacb * uz + mfaca * (1. - uz2) - c1o3 * oMdrho * uz2;
//                        m2    = mfacc * c1o2 + mfacb * (uz + c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
//                        mfaca = m0;
//                        mfacb = m1;
//                        mfacc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfbac * c1o2 + mfbab * (uz - c1o2) + mfbaa * (uz2 - uz) * c1o2;
//                        m1    = -mfbac - 2. * mfbab * uz + mfbaa * (1. - uz2);
//                        m2    = mfbac * c1o2 + mfbab * (uz + c1o2) + mfbaa * (uz2 + uz) * c1o2;
//                        mfbaa = m0;
//                        mfbab = m1;
//                        mfbac = m2;
//                        /////////b//////////////////////////////////////////////////////////////////////////
//                        m0    = mfbbc * c1o2 + mfbbb * (uz - c1o2) + mfbba * (uz2 - uz) * c1o2;
//                        m1    = -mfbbc - 2. * mfbbb * uz + mfbba * (1. - uz2);
//                        m2    = mfbbc * c1o2 + mfbbb * (uz + c1o2) + mfbba * (uz2 + uz) * c1o2;
//                        mfbba = m0;
//                        mfbbb = m1;
//                        mfbbc = m2;
//                        /////////b//////////////////////////////////////////////////////////////////////////
//                        m0    = mfbcc * c1o2 + mfbcb * (uz - c1o2) + mfbca * (uz2 - uz) * c1o2;
//                        m1    = -mfbcc - 2. * mfbcb * uz + mfbca * (1. - uz2);
//                        m2    = mfbcc * c1o2 + mfbcb * (uz + c1o2) + mfbca * (uz2 + uz) * c1o2;
//                        mfbca = m0;
//                        mfbcb = m1;
//                        mfbcc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcac * c1o2 + mfcab * (uz - c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
//                        m1    = -mfcac - 2. * mfcab * uz + mfcaa * (1. - uz2) - c1o3 * oMdrho * uz2;
//                        m2    = mfcac * c1o2 + mfcab * (uz + c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
//                        mfcaa = m0;
//                        mfcab = m1;
//                        mfcac = m2;
//                        /////////c//////////////////////////////////////////////////////////////////////////
//                        m0    = mfcbc * c1o2 + mfcbb * (uz - c1o2) + mfcba * (uz2 - uz) * c1o2;
//                        m1    = -mfcbc - 2. * mfcbb * uz + mfcba * (1. - uz2);
//                        m2    = mfcbc * c1o2 + mfcbb * (uz + c1o2) + mfcba * (uz2 + uz) * c1o2;
//                        mfcba = m0;
//                        mfcbb = m1;
//                        mfcbc = m2;
//                        /////////c//////////////////////////////////////////////////////////////////////////
//                        m0    = mfccc * c1o2 + mfccb * (uz - c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 - uz) * c1o2;
//                        m1    = -mfccc - 2. * mfccb * uz + mfcca * (1. - uz2) - c1o9 * oMdrho * uz2;
//                        m2    = mfccc * c1o2 + mfccb * (uz + c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 + uz) * c1o2;
//                        mfcca = m0;
//                        mfccb = m1;
//                        mfccc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // Y - Dir
//                        m0    = mfaca * c1o2 + mfaba * (uy - c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfaca - 2. * mfaba * uy + mfaaa * (1. - uy2) - c1o6 * oMdrho * uy2;
//                        m2    = mfaca * c1o2 + mfaba * (uy + c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfaaa = m0;
//                        mfaba = m1;
//                        mfaca = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfacb * c1o2 + mfabb * (uy - c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfacb - 2. * mfabb * uy + mfaab * (1. - uy2) - c2o3 * oMdrho * uy2;
//                        m2    = mfacb * c1o2 + mfabb * (uy + c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfaab = m0;
//                        mfabb = m1;
//                        mfacb = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfacc * c1o2 + mfabc * (uy - c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfacc - 2. * mfabc * uy + mfaac * (1. - uy2) - c1o6 * oMdrho * uy2;
//                        m2    = mfacc * c1o2 + mfabc * (uy + c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfaac = m0;
//                        mfabc = m1;
//                        mfacc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfbca * c1o2 + mfbba * (uy - c1o2) + mfbaa * (uy2 - uy) * c1o2;
//                        m1    = -mfbca - 2. * mfbba * uy + mfbaa * (1. - uy2);
//                        m2    = mfbca * c1o2 + mfbba * (uy + c1o2) + mfbaa * (uy2 + uy) * c1o2;
//                        mfbaa = m0;
//                        mfbba = m1;
//                        mfbca = m2;
//                        /////////b//////////////////////////////////////////////////////////////////////////
//                        m0    = mfbcb * c1o2 + mfbbb * (uy - c1o2) + mfbab * (uy2 - uy) * c1o2;
//                        m1    = -mfbcb - 2. * mfbbb * uy + mfbab * (1. - uy2);
//                        m2    = mfbcb * c1o2 + mfbbb * (uy + c1o2) + mfbab * (uy2 + uy) * c1o2;
//                        mfbab = m0;
//                        mfbbb = m1;
//                        mfbcb = m2;
//                        /////////b//////////////////////////////////////////////////////////////////////////
//                        m0    = mfbcc * c1o2 + mfbbc * (uy - c1o2) + mfbac * (uy2 - uy) * c1o2;
//                        m1    = -mfbcc - 2. * mfbbc * uy + mfbac * (1. - uy2);
//                        m2    = mfbcc * c1o2 + mfbbc * (uy + c1o2) + mfbac * (uy2 + uy) * c1o2;
//                        mfbac = m0;
//                        mfbbc = m1;
//                        mfbcc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcca * c1o2 + mfcba * (uy - c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfcca - 2. * mfcba * uy + mfcaa * (1. - uy2) - c1o18 * oMdrho * uy2;
//                        m2    = mfcca * c1o2 + mfcba * (uy + c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfcaa = m0;
//                        mfcba = m1;
//                        mfcca = m2;
//                        /////////c//////////////////////////////////////////////////////////////////////////
//                        m0    = mfccb * c1o2 + mfcbb * (uy - c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfccb - 2. * mfcbb * uy + mfcab * (1. - uy2) - c2o9 * oMdrho * uy2;
//                        m2    = mfccb * c1o2 + mfcbb * (uy + c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfcab = m0;
//                        mfcbb = m1;
//                        mfccb = m2;
//                        /////////c//////////////////////////////////////////////////////////////////////////
//                        m0    = mfccc * c1o2 + mfcbc * (uy - c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
//                        m1    = -mfccc - 2. * mfcbc * uy + mfcac * (1. - uy2) - c1o18 * oMdrho * uy2;
//                        m2    = mfccc * c1o2 + mfcbc * (uy + c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
//                        mfcac = m0;
//                        mfcbc = m1;
//                        mfccc = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        // X - Dir
//                        m0    = mfcaa * c1o2 + mfbaa * (ux - c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcaa - 2. * mfbaa * ux + mfaaa * (1. - ux2) - c1o36 * oMdrho * ux2;
//                        m2    = mfcaa * c1o2 + mfbaa * (ux + c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfaaa = m0;
//                        mfbaa = m1;
//                        mfcaa = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcba * c1o2 + mfbba * (ux - c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcba - 2. * mfbba * ux + mfaba * (1. - ux2) - c1o9 * oMdrho * ux2;
//                        m2    = mfcba * c1o2 + mfbba * (ux + c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfaba = m0;
//                        mfbba = m1;
//                        mfcba = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcca * c1o2 + mfbca * (ux - c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcca - 2. * mfbca * ux + mfaca * (1. - ux2) - c1o36 * oMdrho * ux2;
//                        m2    = mfcca * c1o2 + mfbca * (ux + c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfaca = m0;
//                        mfbca = m1;
//                        mfcca = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcab * c1o2 + mfbab * (ux - c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcab - 2. * mfbab * ux + mfaab * (1. - ux2) - c1o9 * oMdrho * ux2;
//                        m2    = mfcab * c1o2 + mfbab * (ux + c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfaab = m0;
//                        mfbab = m1;
//                        mfcab = m2;
//                        ///////////b////////////////////////////////////////////////////////////////////////
//                        m0    = mfcbb * c1o2 + mfbbb * (ux - c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcbb - 2. * mfbbb * ux + mfabb * (1. - ux2) - c4o9 * oMdrho * ux2;
//                        m2    = mfcbb * c1o2 + mfbbb * (ux + c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfabb = m0;
//                        mfbbb = m1;
//                        mfcbb = m2;
//                        ///////////b////////////////////////////////////////////////////////////////////////
//                        m0    = mfccb * c1o2 + mfbcb * (ux - c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfccb - 2. * mfbcb * ux + mfacb * (1. - ux2) - c1o9 * oMdrho * ux2;
//                        m2    = mfccb * c1o2 + mfbcb * (ux + c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfacb = m0;
//                        mfbcb = m1;
//                        mfccb = m2;
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        ////////////////////////////////////////////////////////////////////////////////////
//                        m0    = mfcac * c1o2 + mfbac * (ux - c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcac - 2. * mfbac * ux + mfaac * (1. - ux2) - c1o36 * oMdrho * ux2;
//                        m2    = mfcac * c1o2 + mfbac * (ux + c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfaac = m0;
//                        mfbac = m1;
//                        mfcac = m2;
//                        ///////////c////////////////////////////////////////////////////////////////////////
//                        m0    = mfcbc * c1o2 + mfbbc * (ux - c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfcbc - 2. * mfbbc * ux + mfabc * (1. - ux2) - c1o9 * oMdrho * ux2;
//                        m2    = mfcbc * c1o2 + mfbbc * (ux + c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfabc = m0;
//                        mfbbc = m1;
//                        mfcbc = m2;
//                        ///////////c////////////////////////////////////////////////////////////////////////
//                        m0    = mfccc * c1o2 + mfbcc * (ux - c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
//                        m1    = -mfccc - 2. * mfbcc * ux + mfacc * (1. - ux2) - c1o36 * oMdrho * ux2;
//                        m2    = mfccc * c1o2 + mfbcc * (ux + c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
//                        mfacc = m0;
//                        mfbcc = m1;
//                        mfccc = m2;
//
//                        ///////////////////////////////////////////////////////////////////////////
//
//                        //////////////////////////////////////////////////////////////////////////
//                        // proof correctness
//                        //////////////////////////////////////////////////////////////////////////
//#ifdef PROOF_CORRECTNESS
//                        LBMReal rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) +
//                                           (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) +
//                                           (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) +
//                                           (mfbba + mfbbc) + mfbbb;
//
//                        LBMReal dif = rho1 - rho_post;
//#ifdef SINGLEPRECISION
//                        if (dif > 10.0E-7 || dif < -10.0E-7)
//#else
//                        if (dif > 10.0E-15 || dif < -10.0E-15)
//#endif
//                        {
//                            UB_THROW(UbException(UB_EXARGS,
//                                                 "rho=" + UbSystem::toString(rho) + ", rho_post=" +
//                                                     UbSystem::toString(rho_post) + " dif=" + UbSystem::toString(dif) +
//                                                     " rho is not correct for node " + UbSystem::toString(x1) + "," +
//                                                     UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
//                        }
//#endif
//
//                        mfcbb = rho * c1o3 * (mfcbb) + 0.5 * forcingTerm[DIR_P00];
//                        mfbcb = rho * c1o3 * (mfbcb) + 0.5 * forcingTerm[N];
//                        mfbbc = rho * c1o3 * (mfbbc) + 0.5 * forcingTerm[T];
//                        mfccb = rho * c1o3 * (mfccb) + 0.5 * forcingTerm[NE];
//                        mfacb = rho * c1o3 * (mfacb) + 0.5 * forcingTerm[NW];
//                        mfcbc = rho * c1o3 * (mfcbc) + 0.5 * forcingTerm[TE];
//                        mfabc = rho * c1o3 * (mfabc) + 0.5 * forcingTerm[TW];
//                        mfbcc = rho * c1o3 * (mfbcc) + 0.5 * forcingTerm[TN];
//                        mfbac = rho * c1o3 * (mfbac) + 0.5 * forcingTerm[TS];
//                        mfccc = rho * c1o3 * (mfccc) + 0.5 * forcingTerm[TNE];
//                        mfacc = rho * c1o3 * (mfacc) + 0.5 * forcingTerm[TNW];
//                        mfcac = rho * c1o3 * (mfcac) + 0.5 * forcingTerm[TSE];
//                        mfaac = rho * c1o3 * (mfaac) + 0.5 * forcingTerm[TSW];
//                        mfabb = rho * c1o3 * (mfabb) + 0.5 * forcingTerm[W];
//                        mfbab = rho * c1o3 * (mfbab) + 0.5 * forcingTerm[S];
//                        mfbba = rho * c1o3 * (mfbba) + 0.5 * forcingTerm[B];
//                        mfaab = rho * c1o3 * (mfaab) + 0.5 * forcingTerm[SW];
//                        mfcab = rho * c1o3 * (mfcab) + 0.5 * forcingTerm[SE];
//                        mfaba = rho * c1o3 * (mfaba) + 0.5 * forcingTerm[BW];
//                        mfcba = rho * c1o3 * (mfcba) + 0.5 * forcingTerm[BE];
//                        mfbaa = rho * c1o3 * (mfbaa) + 0.5 * forcingTerm[BS];
//                        mfbca = rho * c1o3 * (mfbca) + 0.5 * forcingTerm[BN];
//                        mfaaa = rho * c1o3 * (mfaaa) + 0.5 * forcingTerm[BSW];
//                        mfcaa = rho * c1o3 * (mfcaa) + 0.5 * forcingTerm[BSE];
//                        mfaca = rho * c1o3 * (mfaca) + 0.5 * forcingTerm[BNW];
//                        mfcca = rho * c1o3 * (mfcca) + 0.5 * forcingTerm[BNE];
//                        mfbbb = rho * c1o3 * (mfbbb) + 0.5 * forcingTerm[REST];
//
//                        //////////////////////////////////////////////////////////////////////////
//                        // write distribution for F
//                        //////////////////////////////////////////////////////////////////////////
//
//                        (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3)     = mfabb;
//                        (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3)     = mfbab;
//                        (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)     = mfbba;
//                        (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3)    = mfaab;
//                        (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3)   = mfcab;
//                        (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3)    = mfaba;
//                        (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3)   = mfcba;
//                        (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3)    = mfbaa;
//                        (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3)   = mfbca;
//                        (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3)   = mfaaa;
//                        (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3)  = mfcaa;
//                        (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3)  = mfaca;
//                        (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;
//
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)     = mfcbb;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3)     = mfbcb;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p)     = mfbbc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfccb;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3)    = mfacb;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfcbc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p)    = mfabc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfbcc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p)    = mfbac;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfacc;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfcac;
//                        (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfaac;
//
//                        (*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;
// !Old Kernel
                        /////////////////////  P H A S E - F I E L D   S O L V E R
                        ////////////////////////////////////////////
		/////CUMULANT PHASE-FIELD
				LBMReal omegaD =1.0/( 3.0 * mob + 0.5);

			   mfcbb = (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3);
			   mfbcb = (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3);
			   mfbbc = (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3);
			   mfccb = (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3);
			   mfacb = (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3);
			   mfcbc = (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3);
			   mfabc = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3);
			   mfbcc = (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3);
			   mfbac = (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3);
			   mfccc = (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3);
			   mfacc = (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3);
			   mfcac = (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3);
			   mfaac = (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3);
			   mfabb = (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3);
			   mfbab = (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3);
			   mfbba = (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p);
			   mfaab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3);
			   mfcab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3);
			   mfaba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p);
			   mfcba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p);
			   mfbaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p);
			   mfbca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p);
			   mfaaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p);
			   mfcaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p);
			   mfaca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p);
			   mfcca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p);
			   mfbbb = (*this->zeroDistributionsH)(x1, x2, x3);


					////////////////////////////////////////////////////////////////////////////////////
		//! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
		//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
		//!
		////////////////////////////////////////////////////////////////////////////////////
		// fluid component
			   //LBMReal drhoFluid =
				  // ((((fccc + faaa) + (faca + fcac)) + ((facc + fcaa) + (faac + fcca))) +
				  // (((fbac + fbca) + (fbaa + fbcc)) + ((fabc + fcba) + (faba + fcbc)) + ((facb + fcab) + (faab + fccb))) +
					 //  ((fabb + fcbb) + (fbab + fbcb) + (fbba + fbbc))) + fbbb;

			   //LBMReal rhoFluid = c1 + drhoFluid;
			   //LBMReal OOrhoFluid = c1 / rhoFluid;


			   //LBMReal vvx =
				  // ((((fccc - faaa) + (fcac - faca)) + ((fcaa - facc) + (fcca - faac))) +
				  // (((fcba - fabc) + (fcbc - faba)) + ((fcab - facb) + (fccb - faab))) +
					 //  (fcbb - fabb)) * OOrhoFluid;
			   //LBMReal vvy =
				  // ((((fccc - faaa) + (faca - fcac)) + ((facc - fcaa) + (fcca - faac))) +
				  // (((fbca - fbac) + (fbcc - fbaa)) + ((facb - fcab) + (fccb - faab))) +
					 //  (fbcb - fbab)) * OOrhoFluid;
			   //LBMReal vvz =
				  // ((((fccc - faaa) + (fcac - faca)) + ((facc - fcaa) + (faac - fcca))) +
				  // (((fbac - fbca) + (fbcc - fbaa)) + ((fabc - fcba) + (fcbc - faba))) +
					 //  (fbbc - fbba)) * OOrhoFluid;

			 //  LBMReal vvx = ux;
			 //  LBMReal vvy = uy;
			 //  LBMReal vvz = uz;
			   ////////////////////////////////////////////////////////////////////////////////////
			   // second component
			   LBMReal concentration =
				   ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				   (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
					   ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
			   ////////////////////////////////////////////////////////////////////////////////////
			   //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
			   //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
			   //!
			  // LBMReal fx = forces[0];
			  // LBMReal fy = forces[1];
			  // LBMReal fz = -concentration * forces[2];
			  // vvx += fx * c1o2;
			  // vvy += fy * c1o2;
			  // vvz += fz * c1o2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   LBMReal oneMinusRho = c1- concentration;

			   LBMReal cx =
				   ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
				   (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					   (mfcbb - mfabb));
			   LBMReal cy =
				   ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
				   (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					   (mfbcb - mfbab));
			   LBMReal cz =
				   ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
				   (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					   (mfbbc - mfbba));

			   ////////////////////////////////////////////////////////////////////////////////////
			   // calculate the square of velocities for this lattice node
			   LBMReal cx2 = cx * cx;
			   LBMReal cy2 = cy * cy;
			   LBMReal cz2 = cz * cz;
			   ////////////////////////////////////////////////////////////////////////////////////
			   //! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
			   //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
			   //! see also Eq. (6)-(14) in \ref
			   //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
			   //!
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   forwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, cz, cz2, c36, c1o36, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, cz, cz2, c9, c1o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, cz, cz2, c36, c1o36, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, cz, cz2, c9, c1o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, cz, cz2, c9o4, c4o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, cz, cz2, c9, c1o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, cz, cz2, c36, c1o36, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, cz, cz2, c9, c1o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, cz, cz2, c36, c1o36, oneMinusRho);

			   ////////////////////////////////////////////////////////////////////////////////////
			   // Y - Dir
			   forwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, cy, cy2, c6, c1o6, oneMinusRho);
			   forwardChimera(mfaab, mfabb, mfacb, cy, cy2);
			   forwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, cy, cy2, c18, c1o18, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, cy, cy2, c3o2, c2o3, oneMinusRho);
			   forwardChimera(mfbab, mfbbb, mfbcb, cy, cy2);
			   forwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, cy, cy2, c9o2, c2o9, oneMinusRho);
			   forwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, cy, cy2, c6, c1o6, oneMinusRho);
			   forwardChimera(mfcab, mfcbb, mfccb, cy, cy2);
			   forwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, cy, cy2, c18, c1o18, oneMinusRho);

			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   forwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, cx, cx2, c1, c1, oneMinusRho);
			   forwardChimera(mfaba, mfbba, mfcba, cx, cx2);
			   forwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, cx, cx2, c3, c1o3, oneMinusRho);
			   forwardChimera(mfaab, mfbab, mfcab, cx, cx2);
			   forwardChimera(mfabb, mfbbb, mfcbb, cx, cx2);
			   forwardChimera(mfacb, mfbcb, mfccb, cx, cx2);
			   forwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, cx, cx2, c3, c1o3, oneMinusRho);
			   forwardChimera(mfabc, mfbbc, mfcbc, cx, cx2);
			   forwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, cx, cx2, c3, c1o9, oneMinusRho);

			   ////////////////////////////////////////////////////////////////////////////////////
			   //! - experimental Cumulant ... to be published ... hopefully
			   //!

			   // linearized orthogonalization of 3rd order central moments
			   LBMReal Mabc = mfabc - mfaba * c1o3;
			   LBMReal Mbca = mfbca - mfbaa * c1o3;
			   LBMReal Macb = mfacb - mfaab * c1o3;
			   LBMReal Mcba = mfcba - mfaba * c1o3;
			   LBMReal Mcab = mfcab - mfaab * c1o3;
			   LBMReal Mbac = mfbac - mfbaa * c1o3;
			   // linearized orthogonalization of 5th order central moments
			   LBMReal Mcbc = mfcbc - mfaba * c1o9;
			   LBMReal Mbcc = mfbcc - mfbaa * c1o9;
			   LBMReal Mccb = mfccb - mfaab * c1o9;

			   // collision of 1st order moments
			   cx = cx * (c1 - omegaD) + omegaD * vvx * concentration + normX1 * (c1 - 0.5 * omegaD) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * c1o3 * oneOverInterfaceScale;
			   cy = cy * (c1 - omegaD) + omegaD * vvy * concentration + normX2 * (c1 - 0.5 * omegaD) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * c1o3 * oneOverInterfaceScale;
			   cz = cz * (c1 - omegaD) + omegaD * vvz * concentration + normX3 * (c1 - 0.5 * omegaD) * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * c1o3 * oneOverInterfaceScale;

			   //mhx = (ux * phi[REST] + normX1 * (tauH - 0.5) * (1.0 - phi[REST]) * (phi[REST])) / tauH + (1.0 - 1.0 / tauH) * mhx;
			   //mhy = (uy * phi[REST] + normX2 * (tauH - 0.5) * (1.0 - phi[REST]) * (phi[REST])) / tauH + (1.0 - 1.0 / tauH) * mhy;
			//mhz = (uz * phi[REST] + normX3 * (tauH - 0.5) * (1.0 - phi[REST]) * (phi[REST])) / tauH + (1.0 - 1.0 / tauH) * mhz;


			   cx2 = cx * cx;
			   cy2 = cy * cy;
			   cz2 = cz * cz;

			   // equilibration of 2nd order moments
			   mfbba = zeroReal;
			   mfbab = zeroReal;
			   mfabb = zeroReal;

			   mfcaa = c1o3 * concentration;
			   mfaca = c1o3 * concentration;
			   mfaac = c1o3 * concentration;


			   //LBMReal omega2 = 1.0f;// omegaD;
			   //mfbba *= (c1 - omega2);
			   //mfbab *= (c1 - omega2);
			   //mfabb *= (c1 - omega2);

			   //mfcaa = mfcaa*(c1 - omega2) + omega2*c1o3 * concentration;
			   //mfaca = mfaca*(c1 - omega2) + omega2*c1o3 * concentration;
			   //mfaac = mfaac*(c1 - omega2) + omega2*c1o3 * concentration;

			   // equilibration of 3rd order moments
			   Mabc = zeroReal;
			   Mbca = zeroReal;
			   Macb = zeroReal;
			   Mcba = zeroReal;
			   Mcab = zeroReal;
			   Mbac = zeroReal;
			   mfbbb = zeroReal;

			   // from linearized orthogonalization 3rd order central moments to central moments
			   mfabc = Mabc + mfaba * c1o3;
			   mfbca = Mbca + mfbaa * c1o3;
			   mfacb = Macb + mfaab * c1o3;
			   mfcba = Mcba + mfaba * c1o3;
			   mfcab = Mcab + mfaab * c1o3;
			   mfbac = Mbac + mfbaa * c1o3;

			   // equilibration of 4th order moments
			   mfacc = c1o9 * concentration;
			   mfcac = c1o9 * concentration;
			   mfcca = c1o9 * concentration;

			   mfcbb = zeroReal;
			   mfbcb = zeroReal;
			   mfbbc = zeroReal;

			   // equilibration of 5th order moments
			   Mcbc = zeroReal;
			   Mbcc = zeroReal;
			   Mccb = zeroReal;

			   // from linearized orthogonalization 5th order central moments to central moments
			   mfcbc = Mcbc + mfaba * c1o9;
			   mfbcc = Mbcc + mfbaa * c1o9;
			   mfccb = Mccb + mfaab * c1o9;

			   // equilibration of 6th order moment
			   mfccc = c1o27 * concentration;

			   ////////////////////////////////////////////////////////////////////////////////////
			   //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
			   //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
			   //! see also Eq. (88)-(96) in
			   //! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
			   //!
			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   backwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, cx, cx2, c1, c1, oneMinusRho);
			   backwardChimera(mfaba, mfbba, mfcba, cx, cx2);
			   backwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, cx, cx2, c3, c1o3, oneMinusRho);
			   backwardChimera(mfaab, mfbab, mfcab, cx, cx2);
			   backwardChimera(mfabb, mfbbb, mfcbb, cx, cx2);
			   backwardChimera(mfacb, mfbcb, mfccb, cx, cx2);
			   backwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, cx, cx2, c3, c1o3, oneMinusRho);
			   backwardChimera(mfabc, mfbbc, mfcbc, cx, cx2);
			   backwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, cx, cx2, c9, c1o9, oneMinusRho);

			   ////////////////////////////////////////////////////////////////////////////////////
			   // Y - Dir
			   backwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, cy, cy2, c6, c1o6, oneMinusRho);
			   backwardChimera(mfaab, mfabb, mfacb, cy, cy2);
			   backwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, cy, cy2, c18, c1o18, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, cy, cy2, c3o2, c2o3, oneMinusRho);
			   backwardChimera(mfbab, mfbbb, mfbcb, cy, cy2);
			   backwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, cy, cy2, c9o2, c2o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, cy, cy2, c6, c1o6, oneMinusRho);
			   backwardChimera(mfcab, mfcbb, mfccb, cy, cy2);
			   backwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, cy, cy2, c18, c1o18, oneMinusRho);

			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   backwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, cz, cz2, c36, c1o36, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, cz, cz2, c9, c1o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, cz, cz2, c36, c1o36, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, cz, cz2, c9, c1o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, cz, cz2, c9o4, c4o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, cz, cz2, c9, c1o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, cz, cz2, c36, c1o36, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, cz, cz2, c9, c1o9, oneMinusRho);
			   backwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, cz, cz2, c36, c1o36, oneMinusRho);



			   (*this->localDistributionsH)(D3Q27System::ET_E,   x1,  x2,  x3) = mfabb;
   (*this->localDistributionsH)(D3Q27System::ET_N,   x1,  x2,  x3) = mfbab;
   (*this->localDistributionsH)(D3Q27System::ET_T,   x1,  x2,  x3) = mfbba;
   (*this->localDistributionsH)(D3Q27System::ET_NE,  x1,  x2,  x3) = mfaab;
   (*this->localDistributionsH)(D3Q27System::ET_NW,  x1p, x2,  x3) = mfcab;
   (*this->localDistributionsH)(D3Q27System::ET_TE,  x1,  x2,  x3) = mfaba;
   (*this->localDistributionsH)(D3Q27System::ET_TW,  x1p, x2,  x3) = mfcba;
   (*this->localDistributionsH)(D3Q27System::ET_TN,  x1,  x2,  x3) = mfbaa;
   (*this->localDistributionsH)(D3Q27System::ET_TS,  x1,  x2p, x3) = mfbca;
   (*this->localDistributionsH)(D3Q27System::ET_TNE, x1,  x2,  x3) = mfaaa;
   (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2,  x3) = mfcaa;
   (*this->localDistributionsH)(D3Q27System::ET_TSE, x1,  x2p, x3) = mfaca;
   (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

   (*this->nonLocalDistributionsH)(D3Q27System::ET_W,   x1p, x2,  x3 ) = mfcbb;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_S,   x1,  x2p, x3 ) = mfbcb;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_B,   x1,  x2,  x3p) = mfbbc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = mfccb;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = mfacb;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,  x1p, x2,  x3p) = mfcbc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,  x1,  x2,  x3p) = mfabc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,  x1,  x2p, x3p) = mfbcc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,  x1,  x2,  x3p) = mfbac;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1,  x2p, x3p) = mfacc;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2,  x3p) = mfcac;
   (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1,  x2,  x3p) = mfaac;

   (*this->zeroDistributionsH)(x1,x2,x3) = mfbbb;

		/////!CUMULANT PHASE-FIELD



                        ///////////////////   PHASE-FIELD BGK SOLVER ///////////////////////////////
//using namespace D3Q27System;

      //                  h[DIR_P00]   = (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3);
      //                  h[N]   = (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3);
      //                  h[T]   = (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3);
      //                  h[NE]  = (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3);
      //                  h[NW]  = (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3);
      //                  h[TE]  = (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3);
      //                  h[TW]  = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3);
      //                  h[TN]  = (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3);
      //                  h[TS]  = (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3);
      //                  h[TNE] = (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3);
      //                  h[TNW] = (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3);
      //                  h[TSE] = (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3);
      //                  h[TSW] = (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3);

      //                  h[W]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3);
      //                  h[S]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3);
      //                  h[B]   = (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p);
      //                  h[SW]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3);
      //                  h[SE]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3);
      //                  h[BW]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p);
      //                  h[BE]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p);
      //                  h[BS]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p);
      //                  h[BN]  = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p);
      //                  h[BSW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p);
      //                  h[BSE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p);
      //                  h[BNW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p);
      //                  h[BNE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p);

      //                  h[REST] = (*this->zeroDistributionsH)(x1, x2, x3);
						////vvx *= 3;
						////vvy *= 3;
						////vvz *= 3;
						////vx2 = vvx * vvx;
						////vy2 = vvy * vvy;
						////vz2 = vvz * vvz;

      //                  for (int dir = STARTF; dir < (ENDF + 1); dir++) {
      //                      LBMReal velProd = DX1[dir] * vvx + DX2[dir] * vvy + DX3[dir] * vvz;
      //                      LBMReal velSq1  = velProd * velProd;
      //                      LBMReal hEq; //, gEq;

      //                      if (dir != REST) {
      //                          LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]]) / 2.0;
      //                          LBMReal hSource     = (tauH - 0.5) * (1.0 - phi[REST]) * (phi[REST]) * (dirGrad_phi) / denom; 
      //                          hEq = phi[REST] * WEIGTH[dir] * (1.0 + 3.0 * velProd + 4.5 * velSq1 - 1.5 * (vx2 + vy2 + vz2)) +                                 hSource * WEIGTH[dir];

      //                          // This corresponds with the collision factor of 1.0 which equals (tauH + 0.5).
      //                          h[dir] = h[dir] - (h[dir] - hEq) / (tauH); 

      //                      } else {
      //                          hEq = phi[REST] * WEIGTH[REST] * (1.0 - 1.5 * (vx2 + vy2 + vz2));
      //                          h[REST] = h[REST] - (h[REST] - hEq) / (tauH); 
      //                      }
      //                  }

      //                  (*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3)     = h[D3Q27System::INV_E];
      //                  (*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3)     = h[D3Q27System::INV_N];
      //                  (*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3)     = h[D3Q27System::INV_T];
      //                  (*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3)    = h[D3Q27System::INV_NE];
      //                  (*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3)   = h[D3Q27System::INV_NW];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3)    = h[D3Q27System::INV_TE];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3)   = h[D3Q27System::INV_TW];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3)    = h[D3Q27System::INV_TN];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3)   = h[D3Q27System::INV_TS];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3)   = h[D3Q27System::INV_TNE];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3)  = h[D3Q27System::INV_TNW];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3)  = h[D3Q27System::INV_TSE];
      //                  (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = h[D3Q27System::INV_TSW];

      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3)     = h[D3Q27System::INV_W];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3)     = h[D3Q27System::INV_S];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p)     = h[D3Q27System::INV_B];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3)   = h[D3Q27System::INV_SW];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3)    = h[D3Q27System::INV_SE];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p)   = h[D3Q27System::INV_BW];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p)    = h[D3Q27System::INV_BE];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p)   = h[D3Q27System::INV_BS];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p)    = h[D3Q27System::INV_BN];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = h[D3Q27System::INV_BSW];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p)  = h[D3Q27System::INV_BSE];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p)  = h[D3Q27System::INV_BNW];
      //                  (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p)   = h[D3Q27System::INV_BNE];

      //                  (*this->zeroDistributionsH)(x1, x2, x3) = h[D3Q27System::REST];

                        ///////////////////   END OF OLD BGK SOLVER ///////////////////////////////
                    }
                }
            }
        
        dataSet->setPhaseField(divU);
		}
}
//////////////////////////////////////////////////////////////////////////

LBMReal MultiphaseScratchCumulantLBMKernel::gradX1_phi()
{
    using namespace D3Q27System;
	return 3.0* ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) + (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) + (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) + (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_PM0] - phi[DIR_MP0]) + (phi[DIR_PP0] - phi[DIR_MM0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_P00] - phi[DIR_M00]));
    //LBMReal sum = 0.0;
    //for (int k = FSTARTDIR; k <= FENDDIR; k++) {
    //    sum += WEIGTH[k] * DX1[k] * phi[k];
    //}
    //return 3.0 * sum;
}

LBMReal MultiphaseScratchCumulantLBMKernel::gradX2_phi()
{
    using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PPM] - phi[DIR_MMP])- (phi[DIR_PMP] - phi[DIR_MPM])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_0PP] - phi[DIR_0MM]) + (phi[DIR_0PM] - phi[DIR_0MP])) + ((phi[DIR_PP0] - phi[DIR_MM0])- (phi[DIR_PM0] - phi[DIR_MP0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_0P0] - phi[DIR_0M0]));
    //LBMReal sum = 0.0;
    //for (int k = FSTARTDIR; k <= FENDDIR; k++) {
    //    sum += WEIGTH[k] * DX2[k] * phi[k];
    //}
    //return 3.0 * sum;
}

LBMReal MultiphaseScratchCumulantLBMKernel::gradX3_phi()
{
    using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) - (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) - (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_0MP] - phi[DIR_0PM]) + (phi[DIR_0PP] - phi[DIR_0MM])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_00P] - phi[DIR_00M]));
    //LBMReal sum = 0.0;
    //for (int k = FSTARTDIR; k <= FENDDIR; k++) {
    //    sum += WEIGTH[k] * DX3[k] * phi[k];
    //}
    //return 3.0 * sum;
}

LBMReal MultiphaseScratchCumulantLBMKernel::nabla2_phi()
{
    using namespace D3Q27System;
    LBMReal sum = 0.0;
	sum += WEIGTH[DIR_PPP] * ((((phi[DIR_PPP] - phi[DIR_000]) + (phi[DIR_MMM] - phi[DIR_000])) + ((phi[DIR_MMP] - phi[DIR_000]) + (phi[DIR_PPM] - phi[DIR_000])))
		+ (((phi[DIR_MPP] - phi[DIR_000]) + (phi[DIR_PMM] - phi[DIR_000])) + ((phi[DIR_PMP] - phi[DIR_000]) + (phi[DIR_MPM] - phi[DIR_000]))));
	sum += WEIGTH[DIR_0PP] * (
			(((phi[DIR_0PP] - phi[DIR_000]) + (phi[DIR_0MM] - phi[DIR_000])) + ((phi[DIR_0MP] - phi[DIR_000]) + (phi[DIR_0PM] - phi[DIR_000])))
		+	(((phi[DIR_P0P] - phi[DIR_000]) + (phi[DIR_M0M] - phi[DIR_000])) + ((phi[DIR_M0P] - phi[DIR_000]) + (phi[DIR_P0M] - phi[DIR_000])))
		+	(((phi[DIR_PP0] - phi[DIR_000]) + (phi[DIR_MM0] - phi[DIR_000])) + ((phi[DIR_MP0] - phi[DIR_000]) + (phi[DIR_PM0] - phi[DIR_000])))
		);
	sum += WEIGTH[DIR_00P] * (
			((phi[DIR_00P] - phi[DIR_000]) + (phi[DIR_00M] - phi[DIR_000]))
		+	((phi[DIR_0P0] - phi[DIR_000]) + (phi[DIR_0M0] - phi[DIR_000]))
		+	((phi[DIR_P00] - phi[DIR_000]) + (phi[DIR_M00] - phi[DIR_000]))
		);
    //for (int k = FSTARTDIR; k <= FENDDIR; k++) {
    //    sum += WEIGTH[k] * (phi[k] - phi[REST]);
    //}
    return 6.0 * sum;
}

void MultiphaseScratchCumulantLBMKernel::computePhasefield()
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

void MultiphaseScratchCumulantLBMKernel::findNeighbors(CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
                                                int x3)
{
    using namespace D3Q27System;

    SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();

    phi[DIR_000] = (*ph)(x1, x2, x3);

    for (int k = FSTARTDIR; k <= FENDDIR; k++) {

        if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
            phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
        } else {
			phi[k] = 0.0;
         }
    }
}

void MultiphaseScratchCumulantLBMKernel::swapDistributions()
{
    LBMKernel::swapDistributions();
    dataSet->getHdistributions()->swap();
}