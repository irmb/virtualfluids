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
//! \file MultiphaseScaleDistributionLBMKernel.cpp
//! \ingroup LBMKernel
//! \author M. Geier, K. Kutscher, Hesameddin Safari
//=======================================================================================

#include "MultiphaseScaleDistributionLBMKernel.h"
#include "BCArray3D.h"
#include "Block3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
#include <cmath>
#include <iostream>
#include <string>
#include "NonNewtonianFluids/LBM/Rheology.h"

using namespace vf::lbm::dir;
using namespace vf::basics::constant;

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
MultiphaseScaleDistributionLBMKernel::MultiphaseScaleDistributionLBMKernel() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
void MultiphaseScaleDistributionLBMKernel::initDataSet()
{
	SPtr<DistributionArray3D> f(new D3Q27EsoTwist3DSplittedVector( nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9));
	SPtr<DistributionArray3D> h(new D3Q27EsoTwist3DSplittedVector( nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9)); // For phase-field
	SPtr<DistributionArray3D> h2(new D3Q27EsoTwist3DSplittedVector(nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9));
	SPtr<PhaseFieldArray3D> divU1(new PhaseFieldArray3D(            nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressure(new  CbArray3D<real, IndexerX3X2X1>(    nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	pressureOld = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	p1Old = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));

	rhoNode = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	vxNode = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	vyNode = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	vzNode = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	dataSet->setFdistributions(f);
	dataSet->setHdistributions(h); // For phase-field
	dataSet->setH2distributions(h2);
	dataSet->setPhaseField(divU1);
	dataSet->setPressureField(pressure);

	phaseField = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.0));
	phaseFieldOld = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 999.0));

	divU = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> MultiphaseScaleDistributionLBMKernel::clone()
{
	SPtr<LBMKernel> kernel(new MultiphaseScaleDistributionLBMKernel());
	kernel->setNX(nx);
	dynamicPointerCast<MultiphaseScaleDistributionLBMKernel>(kernel)->initDataSet();
	kernel->setCollisionFactorMultiphase(this->collFactorL, this->collFactorG);
	kernel->setDensityRatio(this->densityRatio);
	kernel->setMultiphaseModelParameters(this->beta, this->kappa);
	kernel->setContactAngle(this->contactAngle);
	kernel->setPhiL(this->phiL);
	kernel->setPhiH(this->phiH);
	kernel->setPhaseFieldRelaxation(this->tauH);
	kernel->setMobility(this->mob);
	kernel->setInterfaceWidth(this->interfaceWidth);
    kernel->setSigma(this->sigma);

	kernel->setBCSet(bcSet->clone(kernel));
	kernel->setWithForcing(withForcing);
	kernel->setForcingX1(muForcingX1);
	kernel->setForcingX2(muForcingX2);
	kernel->setForcingX3(muForcingX3);
	kernel->setIndex(ix1, ix2, ix3);
	kernel->setDeltaT(deltaT);
	kernel->setGhostLayerWidth(2);
	dynamicPointerCast<MultiphaseScaleDistributionLBMKernel>(kernel)->initForcing();

	return kernel;
}
//////////////////////////////////////////////////////////////////////////
void  MultiphaseScaleDistributionLBMKernel::forwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho) {
	//using namespace UbMath;
	real m2 = mfa + mfc;
	real m1 = mfc - mfa;
	real m0 = m2 + mfb;
	mfa = m0;
	m0 *= Kinverse;
	m0 += oneMinusRho;
	mfb = (m1 * Kinverse - m0 * vv) * K;
	mfc = ((m2 - c2o1 * m1 * vv) * Kinverse + v2 * m0) * K;
}

////////////////////////////////////////////////////////////////////////////////
void  MultiphaseScaleDistributionLBMKernel::backwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho) {
	//using namespace UbMath;
	real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 - vv) * c1o2) * K;
	real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (-v2)) * K;
	mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 + vv) * c1o2) * K;
	mfa = m0;
	mfb = m1;
}


////////////////////////////////////////////////////////////////////////////////
void  MultiphaseScaleDistributionLBMKernel::forwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2) {
	//using namespace UbMath;
	real m1 = (mfa + mfc) + mfb;
	real m2 = mfc - mfa;
	mfc = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
	mfb = m2 - vv * m1;
	mfa = m1;
}


void  MultiphaseScaleDistributionLBMKernel::backwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2) {
	//using namespace UbMath;
	real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
	real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;
	mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
	mfb = mb;
	mfa = ma;
}


void MultiphaseScaleDistributionLBMKernel::calculate(int step)
{
	using namespace D3Q27System;
	//using namespace UbMath;

	forcingX1 = 0.0;
	forcingX2 = 0.0;
	forcingX3 = 0.0;
    real phiLim = 0.5;

	real oneOverInterfaceScale = c4o1 / interfaceWidth; //1.0;//1.5;
														 /////////////////////////////////////

	localDistributionsF    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
	nonLocalDistributionsF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
	zeroDistributionsF     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

	localDistributionsH1    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getLocalDistributions();
	nonLocalDistributionsH1 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getNonLocalDistributions();
	zeroDistributionsH1     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getZeroDistributions();

	localDistributionsH2    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getH2distributions())->getLocalDistributions();
	nonLocalDistributionsH2 = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getH2distributions())->getNonLocalDistributions();
	zeroDistributionsH2     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getH2distributions())->getZeroDistributions();


	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressure = dataSet->getPressureField();

	SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

	const int bcArrayMaxX1 = (int)bcArray->getNX1();
	const int bcArrayMaxX2 = (int)bcArray->getNX2();
	const int bcArrayMaxX3 = (int)bcArray->getNX3();

	int minX1 = ghostLayerWidth;
	int minX2 = ghostLayerWidth;
	int minX3 = ghostLayerWidth;
	int maxX1 = bcArrayMaxX1 - ghostLayerWidth;
	int maxX2 = bcArrayMaxX2 - ghostLayerWidth;
	int maxX3 = bcArrayMaxX3 - ghostLayerWidth;
	real omegaDRho = 1.0;// 1.25;// 1.3;
	for (int x3 = minX3 - ghostLayerWidth; x3 < maxX3 + ghostLayerWidth; x3++) {
		for (int x2 = minX2 - ghostLayerWidth; x2 < maxX2 + ghostLayerWidth; x2++) {
			for (int x1 = minX1 - ghostLayerWidth; x1 < maxX1 + ghostLayerWidth; x1++) {
				if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
					int x1p = x1 + 1;
					int x2p = x2 + 1;
					int x3p = x3 + 1;



					real mfcbb = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);
					real mfbcb = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);
					real mfbbc = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);
					real mfccb = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);
					real mfacb = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);
					real mfcbc = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);
					real mfabc = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);
					real mfbcc = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);
					real mfbac = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);
					real mfccc = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);
					real mfacc = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);
					real mfcac = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);
					real mfaac = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);
					real mfabb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);
					real mfbab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);
					real mfbba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);
					real mfaab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);
					real mfcab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);
					real mfaba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);
					real mfcba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);
					real mfbaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);
					real mfbca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);
					real mfaaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					real mfcaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);
					real mfaca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);
					real mfcca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);
					real mfbbb = (*this->zeroDistributionsH1)(x1, x2, x3);

					omegaDRho = 2.0;// 1.5;
					real phiOld = (*phaseField)(x1, x2, x3);

					(*phaseField)(x1, x2, x3) = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca) + (mfcaa + mfacc))) +
						(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
							((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
								(mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;


					if ((*phaseField)(x1, x2, x3) > 1) {
						(*phaseField)(x1, x2, x3) = c1o1;
					}

					if ((*phaseField)(x1, x2, x3) < 0) {
						(*phaseField)(x1, x2, x3) = 0;

					
				}
			}
		}
	}
	}

	this->swapDistributions();
	for (int x3 = minX3 - ghostLayerWidth+1; x3 < maxX3 + ghostLayerWidth-1; x3++) {
		for (int x2 = minX2 - ghostLayerWidth+1; x2 < maxX2 + ghostLayerWidth-1; x2++) {
			for (int x1 = minX1 - ghostLayerWidth+1; x1 < maxX1 + ghostLayerWidth-1; x1++) {
				if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
					int x1p = x1 + 1;
					int x2p = x2 + 1;
					int x3p = x3 + 1;

					//real mfabb = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);//* rho * c1o3;
     //               real mfbab = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);//* rho * c1o3;
     //               real mfbba = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);//* rho * c1o3;
     //               real mfaab = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);//* rho * c1o3;
     //               real mfcab = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);//* rho * c1o3;
     //               real mfaba = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);//* rho * c1o3;
     //               real mfcba = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);//* rho * c1o3;
     //               real mfbaa = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);//* rho * c1o3;
     //               real mfbca = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);//* rho * c1o3;
     //               real mfaaa = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);//* rho * c1o3;
     //               real mfcaa = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);//* rho * c1o3;
     //               real mfaca = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);//* rho * c1o3;
     //               real mfcca = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);//* rho * c1o3;
     //               real mfcbb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);//* rho * c1o3;
     //               real mfbcb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);//* rho * c1o3;
     //               real mfbbc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);//* rho * c1o3;
     //               real mfccb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);//* rho * c1o3;
     //               real mfacb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);//* rho * c1o3;
     //               real mfcbc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);//* rho * c1o3;
     //               real mfabc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);//* rho * c1o3;
     //               real mfbcc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);//* rho * c1o3;
     //               real mfbac = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);//* rho * c1o3;
     //               real mfccc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);//* rho * c1o3;
     //               real mfacc = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);//* rho * c1o3;
     //               real mfcac = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);//* rho * c1o3;
     //               real mfaac = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);//* rho * c1o3;
     //               real mfbbb = (*this->zeroDistributionsH1)(x1, x2, x3);
					
					SPtr<DistributionArray3D> distributionH = this->getDataSet()->getHdistributions();
					real hh[27];
					distributionH->getDistributionInv(hh, x1, x2, x3);
					real phiD, vxP, vyP, vzP;

					D3Q27System::calcIncompMacroscopicValues(hh, phiD, vxP, vyP, vzP);
					(*phaseFieldOld)(x1, x2, x3) = phiD;
					
					//real mfcbb = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);
					//real mfbcb = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);
					//real mfbbc = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);
					//real mfccb = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);
					//real mfacb = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);
					//real mfcbc = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);
					//real mfabc = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);
					//real mfbcc = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);
					//real mfbac = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);
					//real mfccc = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);
					//real mfacc = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);
					//real mfcac = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);
					//real mfaac = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//real mfabb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);
					//real mfbab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);
					//real mfbba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);
					//real mfaab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);
					//real mfcab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);
					//real mfaba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);
					//real mfcba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);
					//real mfbaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);
					//real mfbca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);
					//real mfaaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//real mfcaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//real mfaca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//real mfcca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);

					//real mfbbb = (*this->zeroDistributionsH1)(x1, x2, x3);
					//(*phaseField)(x1, x2, x3) = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca) + (mfcaa + mfacc))) +
					//	(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
					//		((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
					//			(mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;
					//if ((*phaseField)(x1, x2, x3) > 1) {
					//	(*phaseField)(x1, x2, x3) = c1o1;
					//}

					//if ((*phaseField)(x1, x2, x3) < 0) {
					//	(*phaseField)(x1, x2, x3) = 0;
					//}
					////// read F-distributions for velocity formalism
						 //mfabb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);//* rho * c1o3;
						 //mfbab = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);//* rho * c1o3;
						 //mfbba = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);//* rho * c1o3;
						 //mfaab = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);//* rho * c1o3;
						 //mfcab = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);//* rho * c1o3;
						 //mfaba = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);//* rho * c1o3;
						 //mfcba = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);//* rho * c1o3;
						 //mfbaa = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);//* rho * c1o3;
						 //mfbca = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);//* rho * c1o3;
						 //mfaaa = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);//* rho * c1o3;
						 //mfcaa = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);//* rho * c1o3;
						 //mfaca = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);//* rho * c1o3;
						 //mfcca = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);//* rho * c1o3;
						 //mfcbb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);//* rho * c1o3;
						 //mfbcb = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);//* rho * c1o3;
						 //mfbbc = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);//* rho * c1o3;
						 //mfccb = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);//* rho * c1o3;
						 //mfacb = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);//* rho * c1o3;
						 //mfcbc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);//* rho * c1o3;
						 //mfabc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);//* rho * c1o3;
						 //mfbcc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);//* rho * c1o3;
						 //mfbac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);//* rho * c1o3;
						 //mfccc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);//* rho * c1o3;
						 //mfacc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);//* rho * c1o3;
						 //mfcac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);//* rho * c1o3;
						 //mfaac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);//* rho * c1o3;
						 //mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);

					//mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
					//mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
					//mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
					//mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
					//mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
					//mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
					//mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
					//mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
					//mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
					//mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
					//mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
					//mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
					//mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
					//mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
					//mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
					//mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
					//mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
					//mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
					//mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
					//mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
					//mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
					//mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

					//mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);


					//real drho = (((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

					//(*rhoNode)(x1, x2, x3) = drho;
					//(*vxNode)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
					//	(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					//	(mfcbb - mfabb));
					//(*vyNode)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
					//	(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					//	(mfbcb - mfbab));
					//(*vzNode)(x1, x2, x3) = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
					//	(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					//	(mfbbc - mfbba));

					SPtr<DistributionArray3D> distribution = this->getDataSet()->getFdistributions();
					real ff[27];
					distribution->getDistributionInv(ff, x1, x2, x3);
					real rhoG,vx,vy,vz;
					//real rhoGG = (((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//				+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//				+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;


					//vx= ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
					//	(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					//	(mfcbb - mfabb));
					//vy	 = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
					//		(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					//		(mfbcb - mfbab));
					//vz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
					//		(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					//		(mfbbc - mfbba));
					D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
                    //if (withForcing) {

                    //    real forcingX1 = muForcingX1.Eval();
                    //    real forcingX2 = muForcingX2.Eval();
                    //    real forcingX3 = muForcingX3.Eval();

                    //    vx += (forcingX1)*deltaT * c1o2;
                    //    vy += (forcingX2)*deltaT * c1o2;
                    //    vz += (forcingX3)*deltaT * c1o2;
                    //}

					//if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {  }
					//else { rhoG = 0.0; vx = 0.0; vy = 0.0; vz = 0.0; }
					//// very bad save the world procedure!!!!
					//vx = (vx > 1 || vx < -1) ? 0 : vx;
					//vy = (vy > 1 || vy < -1) ? 0 : vy;
					//vz = (vz > 1 || vz < -1) ? 0 : vz;
					//rhoG = (rhoG > 10 || rhoG < -10) ? 0 : rhoG;
					(*rhoNode)(x1, x2, x3) = rhoG;// *((*phaseField)(x1, x2, x3) > c1o2 ? densityRatio : c1o1);
					(*vxNode)(x1, x2, x3) = vx;
					(*vyNode)(x1, x2, x3) = vy;
					(*vzNode)(x1, x2, x3) = vz;
					//if (fabsf(vx) > 0 && fabsf(vx) < 0.01) {
					//	int test = 0;
					//}



					//if ((*vzNode)(x1, x2, x3) != 0) {
					//	real vvvv = (*vzNode)(x1, x2, x3);
					//	real pppp = vvvv / (*phaseField)(x1, x2, x3);
					//	int ii = 1;
					//}

				}
			}
		}
	}

	SPtr<DistributionArray3D> distribution = this->getDataSet()->getFdistributions();
	real ff[27];
	for (int x3 = minX3 - 1; x3 < maxX3 + 1; x3++) {
		for (int x2 = minX2 - 1; x2 < maxX2 + 1; x2++) {
			for (int x1 = minX1 - 1; x1 < maxX1 + 1; x1++) {
				if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
					int x1p = x1 + 1;
					int x2p = x2 + 1;
					int x3p = x3 + 1;
					findNeighbors(phaseFieldOld, x1, x2, x3);
                    findNeighbors2(phaseField, x1, x2, x3);
					////////////////////////////////Momentum conservation experiment 06.03.2023
					//surfacetension
 

					if ((phi2[DIR_000] <= phiLim) || (phi[DIR_000] <= phiLim) &&
                        (
						(phi[DIR_P00] > phiLim) ||
						(phi[DIR_M00] > phiLim) ||
						(phi[DIR_00P] > phiLim) ||
						(phi[DIR_00M] > phiLim) ||
						(phi[DIR_0M0] > phiLim) ||
						(phi[DIR_0P0] > phiLim) ||
						(phi[DIR_PP0] > phiLim) ||
						(phi[DIR_PM0] > phiLim) ||
						(phi[DIR_P0P] > phiLim) ||
						(phi[DIR_P0M] > phiLim) ||
						(phi[DIR_MP0] > phiLim) ||
						(phi[DIR_MM0] > phiLim) ||
						(phi[DIR_M0P] > phiLim) ||
						(phi[DIR_M0M] > phiLim) ||
						(phi[DIR_0PM] > phiLim) ||
						(phi[DIR_0MM] > phiLim) ||
						(phi[DIR_0PP] > phiLim) ||
						(phi[DIR_0MP] > phiLim) ||
						(phi[DIR_PPP] > phiLim) ||
						(phi[DIR_PMP] > phiLim) ||
						(phi[DIR_MPP] > phiLim) ||
						(phi[DIR_MMP] > phiLim) ||
						(phi[DIR_PPM] > phiLim) ||
						(phi[DIR_PMM] > phiLim) ||
						(phi[DIR_MPM] > phiLim) ||
						(phi[DIR_MMM] > phiLim)
						)) {
						real vx = (*vxNode)(x1, x2, x3);
						real vy = (*vyNode)(x1, x2, x3);
						real vz = (*vzNode)(x1, x2, x3);
                        real rho = (*rhoNode)(x1, x2, x3);
						findNeighbors(phaseField, x1, x2, x3);
                        real dX1_phi = gradX1_phi();
                        real dX2_phi = gradX2_phi();
                        real dX3_phi = gradX3_phi();
						//real curv = computeCurvature_phi();
                        real laplacePressure = c12o1 * sigma * computeCurvature_phi();
						findNeighbors(phaseFieldOld, x1, x2, x3);


						//real sigma = c3o1*c2o1*1e-3;

                        real flowDirection = vx * dX1_phi + vy * dX2_phi + vy * dX3_phi;


//16.03.23 c: BB gas side with updated boundary velocity

						distribution->getDistributionInv(ff, x1, x2, x3);
						real rhoG;
                        if (phi[DIR_000] > phiLim) { // initialization necessary
							real sumRho = 0;
							real sumWeight = 1.e-100;
							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                if ((phi[fdir] <= phiLim)) {
									sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
									sumWeight += WEIGTH[fdir];
								}

							}

							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
                            }
							rhoG = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neighbor?
							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                if ((phi[fdir] > phiLim)) {
									real vxBC = ((*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real vyBC = ((*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real vzBC = ((*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                    //vx = vxBC;
                                    //vy = vyBC;
                                    //vz = vzBC;
									real fPEQNeighbor = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir],c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real fPEQNeighborInv = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                    //real vBC = (fPEQNeighborInv - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
									//real fPEQHere = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , vx,vy,vz);
                                    //real vBC = (fPEQHere - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
									real vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
									real vDir = (D3Q27System::DX1[fdir] * vx + D3Q27System::DX2[fdir] * vy + D3Q27System::DX3[fdir] * vz);
									//vBC = (vBC + vDir) / (c2o1 + vBC - vDir);
                                   // real dvDir = vBC - vDir;
                                    // 27.04.23
                                    real vxI = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                    real vyI = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                    real vzI = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                    real vIDir = (D3Q27System::DX1[fdir] * vxI + D3Q27System::DX2[fdir] * vyI + D3Q27System::DX3[fdir] * vzI);
                                   // real dvDir = (vBC - vIDir) * c1o2;
                                    real dvDir = (vBC - vDir) ;

									//// 3.7.23
                                   // vIDir = (vIDir + vDir) * c1o2;
                                   // real qq = (c1o2 - (*phaseField)(x1, x2, x3)) / ((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) - (*phaseField)(x1, x2, x3));
                                    //vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);

         //                           // vBC = (qq > c1o2) ? vIDir + (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq)) * c3o2
         //                           //                                                : vBC + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * -c1o2 * (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq));

                                    //dvDir = (vBC - vIDir) /(c1o2+qq);
                                    //real fGEQInv = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vxI, vyI, vzI); 

         //                           ///!03.07.2023


									real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
									
									//if ((phi[D3Q27System::INVDIR[fdir]] > phiLim))
									{
										///here we need reconstruction from scrach
									real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									//real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

										//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
										//real fGEQOld = D3Q27System::getIncompFeqForDirection(fdir, (*rhoNode)(x1, x2, x3), vx, vy, vz);
										//real fGEQNew = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy, vz);
                                    real fGEQ = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy, vz); 
									//real fBC = ( fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
									//3.7.23
                                    //real fBC = ((fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - qq)) - c6o1 * WEIGTH[fdir] * (vBC))*c1o2 / (qq + c1o2) + ((fGEQInv - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o2)))*qq/(qq+c1o2);
                                    //qq = 0;   
									//real fBC = ((distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - collFactorG * fGEQ) / (c1o1 - collFactorG) * (1-qq)+qq*distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - c6o1 * WEIGTH[fdir]*vBC) / (qq + c1o1) +
                                    //           (distribution->getDistributionInvForDirection(x1 , x2 , x3, D3Q27System::INVDIR[fdir])) * qq / (qq + c1o1);
								//real fBC = (fGEQ - c3o1*WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
								//13.07.2023
                                real fBC = feqNew + (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1);
                                    //real fBC = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG )) - c6o1 * WEIGTH[fdir] * (vDir);
                                    //real vNG = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
									//real fBC = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG)) - c6o1 * WEIGTH[fdir] * (vDir)/(c1o1-vDir+vNG);
                                    // 15.5.23
                                    //real fBC = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG)) - c6o1 * WEIGTH[fdir] * (vDir);
                                    //real fBC = (c2o1 * (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
                                    //real fBC = (distribution->getDistributionInvForDirection(x1, x2, x3, D3Q27System::INVDIR[fdir]) - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) ;
                                    //real fBC = (fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (collFactorG )/(c3o1-collFactorG)) - c6o1 * WEIGTH[fdir] * (vBC);
										//real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;// fL -feqOLD + feqNew;
										//real fBC = fGG - c6o1 * WEIGTH[fdir] * (vBC);
									distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
                                    ff[D3Q27System::INVDIR[fdir]] = fBC;
									///// other possibility is tor replace the node itself instead of the neighbor (only c1o1 of them is allowed!)
									//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
									//real feqOLD = D3Q27System::getIncompFeqForDirection(fdir, (*rhoNode)(x1 , x2 , x3 ), (*vxNode)(x1 , x2 , x3 ), (*vyNode)(x1 , x2 , x3 ), (*vzNode)(x1 , x2 , x3 ));
									//real feqNew = D3Q27System::getIncompFeqForDirection(fdir, rhoG, (*vxNode)(x1 , x2 , x3 ), (*vyNode)(x1 , x2 , x3 ), (*vzNode)(x1, x2, x3 ));
									//real fBC = fG - feqOLD + feqNew;
									//distribution->setDistributionForDirection(fBC, x1, x2, x3, fdir);


									}
								}
							}
							//distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1, x2, x3, DIR_000);
							{
								real fL = distribution->getDistributionInvForDirection(x1, x2, x3, DIR_000);
								real feqOLD = D3Q27System::getIncompFeqForDirection(DIR_000, (*rhoNode)(x1, x2, x3), vx,vy,vz);
								real feqNew = D3Q27System::getIncompFeqForDirection(DIR_000, rhoG,vx,vy,vz);
								distribution->setDistributionForDirection(fL-feqOLD+feqNew, x1, x2, x3, DIR_000);
							}
                            D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
                            ff[DIR_000] = vx * vx + vy * vy + vz * vz +
                                          (((ff[DIR_MM0] + ff[DIR_PP0]) + (ff[DIR_MP0] + ff[DIR_PM0])) + ((ff[DIR_0MM] + ff[DIR_0PP]) + (ff[DIR_0MP] + ff[DIR_0PM])) + ((ff[DIR_M0M] + ff[DIR_P0P]) + (ff[DIR_M0P] + ff[DIR_P0M])) +
                                           c2o1 * ((((ff[DIR_MMM] + ff[DIR_PPP]) + (ff[DIR_MMP] + ff[DIR_PPM]))) + (((ff[DIR_MPM] + ff[DIR_PMP]) + (ff[DIR_MPP] + ff[DIR_PMM])))));
                            distribution->setDistributionForDirection(ff[DIR_000], x1, x2, x3, DIR_000);

						}
						else {//no refill of gas required
							rhoG = (*rhoNode)(x1, x2, x3);
                            if (phi2[DIR_000] <= phiLim) { // no refill liquid
								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                    if ((phi[fdir] > phiLim)) {
										real vxBC = ((*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real vyBC = ((*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real vzBC = ((*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                        //vx = vxBC;
                                        //vy = vyBC;
                                        //vz = vzBC;
										real vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
									real fPEQNeighbor = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir],c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real fPEQNeighborInv = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                    //real vBC = (fPEQNeighborInv - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
                                    //									real fPEQHere = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , vx,vy,vz);
//                                  real vBC = (fPEQHere - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
										real vDir = (D3Q27System::DX1[fdir] * vx + D3Q27System::DX2[fdir] * vy + D3Q27System::DX3[fdir] * vz);
										//real dvDir = vBC - vDir;
                                        // real dvDir = vBC - vDir;
                                        // 27.04.23
                                        real vxI = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vyI = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vzI = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vIDir = (D3Q27System::DX1[fdir] * vxI + D3Q27System::DX2[fdir] * vyI + D3Q27System::DX3[fdir] * vzI);
                                       // real dvDir = (vBC - vIDir) * c1o2;
                                        real dvDir = (vBC - vDir) ;

										//vBC = (vBC + vDir) / (c2o1 + vBC - vDir);
										real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
										real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
										//real fBC = fG - c6o1 * WEIGTH[fdir] * (vBC);
										//real fGInv = distribution->getDistributionInvForDirection(x1, x2, x3, D3Q27System::INVDIR[fdir]);
										//real fGInvEQ = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz);
										real fGEQ = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy, vz);
										//real fBC = (-fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) )- c6o1 * WEIGTH[fdir] * (vBC);
										//// 3.7.23
                                      //  vIDir = (vIDir + vDir) * c1o2;
                                      //  real qq = (c1o2 - (*phaseField)(x1, x2, x3)) / ((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) - (*phaseField)(x1, x2, x3));
                                      //  vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);

										//
										////vBC = (qq > c1o2) ? vIDir + (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq)) * c3o2
          ////                                                : vBC + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * -c1o2 * (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq));

										//dvDir = (vBC - vIDir)*c2o3;
                                       // dvDir = (vBC - vIDir) / (c1o2 + qq);
                                       // real fGEQInv = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vxI, vyI, vzI); 

										///!03.07.2023
                                        // 3.7.23
                                        //real fBC = ((fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - qq)) - c6o1 * WEIGTH[fdir] * (vBC))*c1o2 / (qq + c1o2) + ((fGEQInv - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o2))) * qq / (qq + c1o2);
                                        //qq = 0;
                                        //real fBC = ((distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - collFactorG * fGEQ) / (c1o1 - collFactorG) * (1 - qq) + qq * distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - c6o1 * WEIGTH[fdir] * vBC) / (qq + c1o1) +
                                        //       (distribution->getDistributionInvForDirection(x1 , x2 , x3, D3Q27System::INVDIR[fdir])) * qq / (qq + c1o1);

										//real fBC = ( fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        // 13.07.2023
                                        real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                        real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

                                        real fBC = feqNew + (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1);

										//real fBC = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        
										//real qq = c1o1 - ((c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) /
                                        //                  (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])));
                                        //real vNG = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
                                        //real fBC = fGEQ + (-WEIGTH[fdir] * dvDir * ((c1o1 + qq) / collFactorG - c2o1 * qq) - c6o1 * WEIGTH[fdir] * (vDir + qq * vNG)) / (c1o1 + qq);
                                        //real fBC = fGEQ + (-WEIGTH[fdir] * dvDir * ((c1o1 + qq) / collFactorG - c2o1 * qq) - c6o1 * WEIGTH[fdir] * (vDir + qq * (vDir + qq * (vNG - vDir))) / (c1o1 - vDir + vNG)) / (c1o1 + qq);
										
										//real fBC = (c2o1*(fGEQ -   WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
                                        //real fBC = (c2o1*(fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
										//real fBC = (fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (collFactorG) / (c3o1 - collFactorG)) - c6o1 * WEIGTH[fdir] * (vBC);
										//26.04.23 flux BC:
                                        //real fBC = (c2o1*(fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
                                        //if (flowDirection > 0) fBC = (fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
										//if (fabsf(-fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) - fGEQ) >1000* (fabsf(fG - fGEQ))) fBC = fG - c6o1 * WEIGTH[fdir] * (vBC);
										//if (fGEQ > 1.0e-8&& step>30&& vyBC!=0) {
										//	std::cout << D3Q27System::DX1[fdir] <<","<< D3Q27System::DX2[fdir] << "," << D3Q27System::DX3[fdir] <<" " << -fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) - fGEQ << " fg:" << fG - fGEQ << " ratio=" << (-fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) - fGEQ) / (fG - fGEQ) << " feq" << fGEQ << " vy =" << vy << "vyBC=" << vyBC << "\n";
										//}

										//real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
										//real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

										//if ((*phaseField)(x1, x2, x3) <= c1o2) 
										distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
                                        //if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > phiLim)
                                        if (phi2[fdir] > phiLim)
										{
											//real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
											real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
											real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
											//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

											//distribution->setDistributionForDirection((fBC + fG) / densityRatio*0 - fL  - (feqG - feqL) * (c1o1 / densityRatio*0 - c1o1) * vBC, x1, x2, x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
											//distribution->setDistributionForDirection((fBC + fG) / densityRatio * 0 - fL - (feqG - feqL-2*fL+2*feqL) * (c1o1 / densityRatio - c1o1) * vBC, x1, x2, x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
											//real flW = (fBC + fG) / densityRatio * 0 - fL - (feqG - feqL) * (c1o1 / densityRatio*0 - c1o1) * vBC;
											//real flWW = (fBC + fG) / densityRatio * 0 - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio*0 - c1o1) * vBC;
											//real fLi = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], fdir);
											//real number = 666;
											//distribution->setDistributionForDirection((fBC + fG) / densityRatio * 0 - fL - (feqG - feqL) * (c1o1 / densityRatio * 0 - c1o1) * vBC, x1, x2, x3, fdir);
											real eqBC= D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, vx, vy, vz);
											real eqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
											real eqBCN = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real eqGN = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

										//real flNew = (fBC + fG-eqBC-eqG) / densityRatio +eqBC+eqG - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio  - c1o1) * vBC;
										


											
                                            // real flNew = (fBC + fG-eqBC-eqG) / densityRatio +eqBC+eqG - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio  - c1o1) * vBC;
                                            real laplacePressureBC;
                                            if ((x1 + D3Q27System::DX1[fdir] > 0) && (x1 + D3Q27System::DX1[fdir] < maxX1 + 1) && (x2 + D3Q27System::DX2[fdir] > 0) && (x2 + D3Q27System::DX2[fdir] < maxX2 + 1) && (x3 + D3Q27System::DX3[fdir] > 0) && (x3 + D3Q27System::DX3[fdir] < maxX3 + 1) &&
                                                phi2[DIR_000] != phi2[fdir]) {
                                                findNeighbors(phaseField, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
                                                laplacePressureBC = c6o1 * c2o1 * computeCurvature_phi() * sigma;
                                                findNeighbors(phaseFieldOld, x1, x2, x3);
                                            } else
                                                laplacePressureBC = laplacePressure; 
												//if (UbMath::isNaN(laplacePressureBC) || UbMath::isInfinity(laplacePressureBC)) {
            //                                    laplacePressureBC = laplacePressure;
            //                                }
											// curv; // reset to the above
                                            if (phi2[DIR_000] != phi2[fdir])
                                                {

                                                    laplacePressureBC = laplacePressure * (c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) /
                                                                            (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) +
                                                                        laplacePressureBC * (-c1o1 + c2o1 * (*phaseField)(x1, x2, x3)) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                                }
                                            else
                                                laplacePressureBC = laplacePressure;
                                            // laplacePressureBC *= sigma;
                                            // eqBCN = eqBC;
                                            // distribution->setDistributionForDirection(LaplacePressure* WEIGTH[fdir] + (fBC + fG - eqBC - eqG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio*0) - fL - 0*(feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio - c1o1) * vBC, x1, x2,
                                            // x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
                                            distribution->setDistributionForDirection(laplacePressureBC * WEIGTH[fdir] + (fBC + fG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio) - fL, x1, x2, x3, fdir);




										//real laplacePressureBC;
          //                                  if ((x1 + D3Q27System::DX1[fdir] > 0) && (x1 + D3Q27System::DX1[fdir] < maxX1 + 1) && (x2 + D3Q27System::DX2[fdir] > 0) && (x2 + D3Q27System::DX2[fdir] < maxX2 + 1) && (x3 + D3Q27System::DX3[fdir] > 0) && (x3 + D3Q27System::DX3[fdir] < maxX3 + 1)) {
          //                                      findNeighbors(phaseField, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
          //                                      laplacePressureBC = c6o1 * c2o1 * computeCurvature_phi() * sigma;
          //                                      findNeighbors(phaseFieldOld, x1, x2, x3);
          //                                  } else
          //                                      laplacePressureBC = laplacePressure; // curv; // reset to the above
          //                                  laplacePressureBC = laplacePressure * (c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) /
          //                                                          (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) +
          //                                                      laplacePressureBC * (-c1o1 + c2o1 * (*phaseField)(x1, x2, x3)) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
          //                                  // laplacePressureBC *= sigma;
          //                                 
          //                              }

										//	
										//	//real curvBC;
										//	//if ((x1 + D3Q27System::DX1[fdir] > 0) && (x1 + D3Q27System::DX1[fdir] < maxX1 + 1) && (x2 + D3Q27System::DX2[fdir] > 0) && (x2 + D3Q27System::DX2[fdir] < maxX2 + 1) && (x3 + D3Q27System::DX3[fdir] > 0) && (x3 + D3Q27System::DX3[fdir] < maxX3 + 1)) {
										//	//	findNeighbors(phaseField, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
										//	//		 curvBC = computeCurvature_phi();
										//	//		findNeighbors(phaseFieldOld, x1, x2, x3);
										//	//}
										//	//else curvBC = curv;//reset to the above
										//	//real LaplacePressure = curv *(c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) + curvBC * (-c1o1 + c2o1 * (*phaseField)(x1, x2, x3)) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//	////16.04.23
										//	//real eqLL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//	////fL = fL*0.99 +0.01*(eqLL - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1));
										//	//LaplacePressure *= sigma;
										//	//eqBCN = eqBC;
										//	//distribution->setDistributionForDirection(LaplacePressure * WEIGTH[fdir] +(fBC + fG - eqBC - eqG) / densityRatio + (eqBCN + eqGN) * (c1o1-c1o1 / densityRatio*0 ) - fL -0* (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio  - c1o1) * vBC, x1, x2, x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
										//	distribution->setDistributionForDirection(laplacePressureBC* WEIGTH[fdir] + (fBC + fG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio ) - fL, x1, x2, x3, fdir);
											//if (vxBC != 0) {
											//	int set = 0;
											//}

										}

									}


								}
							}
							else {//refill liquid

								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                    if ((phi[fdir] > phiLim)) {
										real vxBC = ((*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real vyBC = ((*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real vzBC = ((*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
									real fPEQNeighbor = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir],c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real fPEQNeighborInv = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                    //real vBC = (fPEQNeighborInv - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
									//real fPEQHere = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , vx,vy,vz);
                                    //real vBC = (fPEQHere - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
										real vDir = (D3Q27System::DX1[fdir] * vx + D3Q27System::DX2[fdir] * vy + D3Q27System::DX3[fdir] * vz);
										//real dvDir = vBC - vDir;
										//27.04.23
                                        real vxI = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vyI = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vzI = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                        real vIDir = (D3Q27System::DX1[fdir] * vxI + D3Q27System::DX2[fdir] * vyI + D3Q27System::DX3[fdir] * vzI);
                                      //  real dvDir = (vBC - vIDir)*c1o2;
                                        real dvDir = (vBC - vDir) ;

										//// 3.7.23
                                       // vIDir = (vIDir + vDir) * c1o2;
                                       // real qq = (c1o2 - (*phaseField)(x1, x2, x3)) / ((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) - (*phaseField)(x1, x2, x3));
                                       // vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);

          //                              //vBC = (qq > c1o2) ? vIDir + (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq)) * c3o2
          //                              //                  : vBC + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * -c1o2 * (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq));

          //                          
										//dvDir = (vBC - vIDir) * c2o3;
                                        //dvDir = (vBC - vIDir) / (c1o2 + qq);
                                        //real fGEQInv = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vxI, vyI, vzI); 
          //                              ///!03.07.2023

										//vBC = (vBC + vDir) / (c2o1 + vBC - vDir);
										real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
										real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
										//real fBC = fG - c6o1 * WEIGTH[fdir] * (vBC);
										//alternative way to bounce back by recovering fG from the opiste direction
										//real fGInv= distribution->getDistributionInvForDirection(x1, x2, x3, D3Q27System::INVDIR[fdir]);
										//real fGInvEQ = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz);
										real fGEQ = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy, vz);
										//real fBC = (-fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        // 3.7.23
                                        //real fBC = ((fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - qq)) - c6o1 * WEIGTH[fdir] * (vBC))*c1o2 / (qq + c1o2) + ((fGEQInv - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o2))) * qq / (qq + c1o2);
                                        //qq = 0;
                                        //real fBC = ((distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - collFactorG * fGEQ) / (c1o1 - collFactorG) * (1 - qq) + qq * distribution->getDistributionInvForDirection(x1, x2, x3, fdir) - c6o1 * WEIGTH[fdir] * vBC) / (qq + c1o1) +
                                        //       (distribution->getDistributionInvForDirection(x1 , x2 , x3, D3Q27System::INVDIR[fdir])) * qq / (qq + c1o1);

										//real fBC = (fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        //  13.07.2023
                                        real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                        real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                            (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

                                        real fBC = feqNew + (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1);

										//real fBC = (fGEQ -  WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);

										//real qq = c1o1-((c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) /
                                         //         (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])));
                                        //real vNG = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
                                        //real fBC = fGEQ + (-WEIGTH[fdir] * dvDir * ((c1o1 + qq) / collFactorG - c2o1 * qq) - c6o1 * WEIGTH[fdir] * (vDir + qq * vNG))/(c1o1+qq);
                                        //real fBC = fGEQ + (-WEIGTH[fdir] * dvDir * ((c1o1 + qq) / collFactorG - c2o1 * qq) - c6o1 * WEIGTH[fdir] * (vDir + qq *( vDir+qq*( vNG-vDir)))/(c1o1-vDir+vNG)) / (c1o1 + qq);
										//real fBC = (c2o1 * (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
										//real fBC = (distribution->getDistributionInvForDirection(x1, x2, x3, D3Q27System::INVDIR[fdir]) - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) ;
                                        //real fBC = (fGEQ - c6o1 * WEIGTH[fdir] * dvDir * ( collFactorG )/(c3o1-collFactorG)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        // 26.04.23 flux BC:
                                        //real fBC = (c2o1 * (fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1))) - c6o1 * WEIGTH[fdir] * (vBC)-fG;
                                        //if (flowDirection > 0) fBC = (fGEQ - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
										//if (fabsf(-fGInv + fGInvEQ + fGEQ - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) - fGEQ) > 1000*(fabsf(fG - fGEQ))) fBC = fG - c6o1 * WEIGTH[fdir] * (vBC);
										//real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
										//real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

										ff[D3Q27System::INVDIR[fdir]] = fBC;
                                        if (phi2[fdir] > phiLim) {
											//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
											real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
											real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
											//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

											//distribution->setDistributionForDirection((fBC + fG) / densityRatio*0 - fL- (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vBC), x1, x2, x3, fdir);
											//distribution->setDistributionForDirection((fBC + fG) / densityRatio * 0 - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio*0 - c1o1) * vBC, x1, x2, x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
											//distribution->setDistributionForDirection(0, x1, x2, x3, fdir);
											//real flW = (fBC + fG) / densityRatio * 0 - fL - (feqG - feqL) * (c1o1 / densityRatio * 0 - c1o1) * vBC;
											//real flWW = (fBC + fG) / densityRatio * 0 - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio * 0 - c1o1) * vBC;
											//real fLi = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], fdir);
											real eqBCN = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real eqGN = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real eqBC = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, vx, vy, vz);
											real eqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
											////real flNew = (fBC + fG - eqBC - eqG) / densityRatio + eqBC + eqG - fL - (feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio - c1o1) * vBC;
											//real curvBC;
											//if ((x1 + D3Q27System::DX1[fdir] > 0) && (x1 + D3Q27System::DX1[fdir] < maxX1 + 1) && (x2 + D3Q27System::DX2[fdir] > 0) && (x2 + D3Q27System::DX2[fdir] < maxX2 + 1) && (x3 + D3Q27System::DX3[fdir] > 0) && (x3 + D3Q27System::DX3[fdir] < maxX3 + 1)) {
											//	findNeighbors(phaseField, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
											//	curvBC = computeCurvature_phi();
											//	findNeighbors(phaseFieldOld, x1, x2, x3);
											//}
											//else curvBC = curv;//reset to the above
											////16.04.23
											//real eqLL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											////fL = fL * 0.99 + 0.01 * (eqLL -  c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1));
											//real LaplacePressure = curv *(c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) + curvBC * (-c1o1 + c2o1 * (*phaseField)(x1, x2, x3)) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											//LaplacePressure *= sigma;
											
											
											real laplacePressureBC;
                                            if ((x1 + D3Q27System::DX1[fdir] > 0) && (x1 + D3Q27System::DX1[fdir] < maxX1 + 1) && (x2 + D3Q27System::DX2[fdir] > 0) && (x2 + D3Q27System::DX2[fdir] < maxX2 + 1) && (x3 + D3Q27System::DX3[fdir] > 0) && (x3 + D3Q27System::DX3[fdir] < maxX3 + 1) ) {
                                                findNeighbors(phaseField, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
                                                laplacePressureBC = c6o1 * c2o1 * computeCurvature_phi() * sigma;
                                                findNeighbors(phaseFieldOld, x1, x2, x3);
                                            } else
                                                laplacePressureBC = laplacePressure; // curv; // reset to the above

											//if (UbMath::isNaN(laplacePressureBC) || UbMath::isInfinity(laplacePressureBC)) {
           //                                     laplacePressureBC = laplacePressure;
           //                                 }
                                            if (phi2[DIR_000] != phi2[fdir])
                                            {

                                                laplacePressureBC = laplacePressure * (c1o1 - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) /
                                                                        (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) +
                                                                    laplacePressureBC * (-c1o1 + c2o1 * (*phaseField)(x1, x2, x3)) / (c2o1 * (*phaseField)(x1, x2, x3) - c2o1 * (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                            }
                                            else laplacePressureBC = laplacePressure;

											                                 //               real pp1 = (*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
                                            //real pp2 = (*phaseField)(x1, x2, x3);

											//eqBCN = eqBC;
											//distribution->setDistributionForDirection(LaplacePressure* WEIGTH[fdir] + (fBC + fG - eqBC - eqG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio*0) - fL - 0*(feqG - feqL - 2 * fL + 2 * feqL) * (c1o1 / densityRatio - c1o1) * vBC, x1, x2, x3, fdir);// (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
                                           // fBC = (fG) / (densityRatio - c1o1) +
                                           //       ((densityRatio) / (densityRatio - c1o1)) * ((eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio) - c2o1 * fL + (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1)) + laplacePressureBC * WEIGTH[fdir]);
                                            // 13.07.2023
                                            real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                                (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                                (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                            real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
                                                                                                (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

                                            real fBC = feqNew + (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1);

											
											distribution->setDistributionForDirection(laplacePressureBC* WEIGTH[fdir] + (fBC + fG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio)  - fL , x1, x2, x3, fdir);
										//	real number = 666;



										}

									}
									else {
										ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);;
									}


								}

								real sum2 = 1e-100;
								real sumRho = 0;
								real sumVx = 0;
								real sumVy = 0;
								real sumVz = 0;
								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                    if ((phi[fdir] > phiLim)) {

										sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);// * tempRho;
										sumVx += WEIGTH[fdir] * (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
										sumVy += WEIGTH[fdir] * (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
										sumVz += WEIGTH[fdir] * (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
										sum2 += WEIGTH[fdir];
									}
								}
								real rhoL;
								D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
								rhoL = sumRho / sum2;
								//vx = sumVx / sum2;
								//vy = sumVy / sum2;
								//vz = sumVz / sum2;
								//rhoL = (*rhoNode)(x1, x2, x3)/densityRatio;

								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
									ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
								}

								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
									//if (!((phi[fdir] > c1o2) && (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2))) {
                                    if (!((phi[fdir] > phiLim))) {
											real vxBC = ((*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real vyBC = ((*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real vzBC = ((*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
											real vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
									real fPEQNeighbor = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir],c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
									real fPEQNeighborInv = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                    //real vBC = (fPEQNeighborInv - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
									//real fPEQHere = D3Q27System::getIncompFeqForDirection(fdir,c0o1 , vx,vy,vz);
                                    //real vBC = (fPEQHere - fPEQNeighbor) / WEIGTH[fdir] * c1o6;
											real vDir = (D3Q27System::DX1[fdir] * vx + D3Q27System::DX2[fdir] * vy + D3Q27System::DX3[fdir] * vz);
											//real dvDir = vBC - vDir;
                                            // 27.04.23
                                            real vxI = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                            real vyI = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                            real vzI = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
                                            real vIDir = (D3Q27System::DX1[fdir] * vxI + D3Q27System::DX2[fdir] * vyI + D3Q27System::DX3[fdir] * vzI);
                                            //real dvDir = (vBC - vIDir) * c1o2;
                                            real dvDir = (vBC - vDir) ;


											vBC = (vBC + vDir) / (c2o1 + vBC - vDir);

											//// 3.7.23
           //                                 real qq = (c1o2 - (*phaseField)(x1, x2, x3)) / ((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) - (*phaseField)(x1, x2, x3));
           //                                 vBC = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);

           //                                 //vBC = (qq > c1o2) ? vIDir + (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq)) * c3o2
           //                                 //                  : vBC + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * -c1o2 * (vBC - vIDir) / (c1o1 + qq + (c1o1 / collFactorG - c1o2) / (c1o1 / collFactorL - c1o2) / densityRatio * (c1o1 - qq));

                                            //dvDir = (vBC - vIDir) * c2o3;
                                            ///!03.07.2023


										//real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz);
										//real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoL, vx, vy, vz);
										//ff[D3Q27System::INVDIR[fdir]]=(feqNew - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1));
                                        //ff[D3Q27System::INVDIR[fdir]] = (feqNew - WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1));
                                        real fGEQ = D3Q27System::getIncompFeqForDirection(fdir, rhoL, vx, vy, vz);
                                        //ff[D3Q27System::INVDIR[fdir]] = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1)) - c6o1 * WEIGTH[fdir] * (vBC);
                                        //real vNG = (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX3[fdir] * vzBC);
                                        //ff[D3Q27System::INVDIR[fdir]] = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorL)) - c6o1 * WEIGTH[fdir] * (vDir) / (c1o1 - vDir + vNG);
          //                              real fG, fBCPseudo;
          //                              if (phi2[fdir] <= phiLim)
          //                                  {

          //                                  fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
          //                                  fBCPseudo = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
          //                              }
										//else {
          //                                  // 13.07.2023
          // //                                 real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
										//	//real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                     (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                     (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
          // //                                 real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                     (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

          // //                                 fBCPseudo = feqNew + (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1);


          //                                  //unfortunately we have to reconstruct our own populations because they are overwritten by our neighbor. This is realy ugly but doing it the right way would require a completly other code structure.
										//	//fBCPseudo = -c1o3*WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) + D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                                                                         (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                                                                         (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                                                                                         (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
          // //                                 fG = -c1o3*WEIGTH[fdir] * dvDir * (c1o1 / collFactorG - c1o1) +
          // //                                      D3Q27System::getIncompFeqForDirection(fdir, (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          // //                                                                            (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
          //                              // 16.07.23 attempt to make it at least momentum conserving
          //                                  fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
          //                                  fBCPseudo = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 , x2 , x3 ),
          //                                                                                    (*vxNode)(x1 , x2 , x3 ), (*vyNode)(x1 , x2 , x3 ),
          //                                                                                    (*vzNode)(x1 , x2 , x3 )) + fG - D3Q27System::getIncompFeqForDirection(fdir, (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          //                                                                                    (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										//}

										//real eqBCN = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          //                                                                                 (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
          //                              real eqGN = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]),
          //                                                                                (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

										//ff[D3Q27System::INVDIR[fdir]] = (laplacePressure * WEIGTH[fdir] + (fBCPseudo + fG) / densityRatio + (eqBCN + eqGN) * (c1o1 - c1o1 / densityRatio))*c1o2 +(fG-fBCPseudo)*c1o2;
										//16.07.2023: alternative: In attempt to balance the momentum loss in the generation of new gas nodes the inverse operation is performed here with the liquid phase:
                                        real fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
										real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
										real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoL, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
                                        real fBC = feqNew + (fG - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1);
										//not adding laplace pressure as it should be included in rhoL;
                                        ff[D3Q27System::INVDIR[fdir]] = fBC;
										//15.5.23
										//ff[D3Q27System::INVDIR[fdir]] = (fGEQ - WEIGTH[fdir] * dvDir * (c1o1 / collFactorL)) - c6o1 * WEIGTH[fdir] * (vDir);
										//real fBC = (distribution->getDistributionInvForDirection(x1, x2, x3, D3Q27System::INVDIR[fdir]) - c6o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1)) ;
                                        //ff[D3Q27System::INVDIR[fdir]] = (feqNew - c6o1 * WEIGTH[fdir] * dvDir * (collFactorL)/(c3o1-collFactorL));
										//ff[D3Q27System::INVDIR[fdir]] = (ff[D3Q27System::INVDIR[fdir]] - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew;
										distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
									}
								}
                                real eqRest = D3Q27System::getIncompFeqForDirection(DIR_000, 0, (*vxNode)(x1, x2 , x3 ),
                                                                                   (*vyNode)(x1, x2 , x3 ), (*vzNode)(x1 , x2 , x3 ));
                                real fRest = distribution->getDistributionInvForDirection(x1 , x2 , x3 , DIR_000);
                                distribution->setDistributionForDirection((laplacePressure * WEIGTH[DIR_000] + c2o1*(fRest) / densityRatio + (eqRest) * (c1o1 - c1o1 / densityRatio))  , x1, x2, x3, DIR_000);

                                //03.04.2023 alternative initialization of liquid nodes based on FD
								//for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
								//	//if (!((phi[fdir] > c1o2) && (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2))) {
								//	if (!((phi[fdir] > c1o2))) {
								//		real vxBC = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real vyBC = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real vzBC = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real vBC = (-D3Q27System::DX1[fdir] * vxBC - D3Q27System::DX2[fdir] * vyBC - D3Q27System::DX3[fdir] * vzBC);
								//		real vDir = (-D3Q27System::DX1[fdir] * vx - D3Q27System::DX2[fdir] * vy - D3Q27System::DX3[fdir] * vz);
								//		real dvDir = vBC - vDir;

								//		//real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz);
								//		real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoL, vx, vy, vz);
								//		ff[D3Q27System::INVDIR[fdir]] =  feqNew - c3o1 * WEIGTH[fdir] * dvDir * (c1o1 / collFactorL - c1o1);
								//		distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
								//	}
								//}

								//for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
								//	if ((phi[D3Q27System::INVDIR[fdir]] <= c1o2) && (phi[fdir] > c1o2)) {
								//		//real vxBC = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		//real vyBC = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		//real vzBC = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		//real vBC = -(D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX2[fdir] * vzBC);
								//		real vDir = -(D3Q27System::DX1[fdir] * vx + D3Q27System::DX2[fdir] * vy + D3Q27System::DX2[fdir] * vz);
								//		//vBC = (vBC + vDir) / (c2o1 -( vBC - vDir));
								//		//real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]) - c6o1 * WEIGTH[fdir] * vDir;
								//		//real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]) + c6o1 * WEIGTH[fdir] * (vx * D3Q27System::DX1[fdir] + vy * D3Q27System::DX2[fdir] + vz * D3Q27System::DX3[fdir]);
								//		real fL= D3Q27System::getIncompFeqForDirection(fdir, rhoL, vx, vy, vz);
								//		distribution->setDistributionForDirection(fL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
								//		ff[fdir] = fL;
								//	}
								//	if (!(phi[fdir] > c1o2)) {
								//		//std::cout << "Eq at dir=" << fdir << "\n";
								//		real vxBC = ((*vxNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real vyBC = ((*vyNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real vzBC = ((*vzNode)(x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir]));
								//		real feqL = D3Q27System::getIncompFeqForDirection(fdir, rhoL, vx, vy, vz);
								//		distribution->setDistributionForDirection(feqL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
								//		ff[fdir] = feqL;
								//	}
								//}
						//real sumRho2= 0;
						//for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//	sumRho2 += ff[fdir];// -D3Q27System::getIncompFeqForDirection(fdir, 0, sumVx, sumVy, sumVz);
						//}
						//ff[DIR_000] = rhoL - sumRho2;
						//rhoL = 27.0 / 18.0 * sumRho2;
						//std::cout << "rhoL=" << rhoL <<" sumRho="<< 27.0 / 18.0 * sumRho2 << " vx=" << vx << " vy=" << vy << "\n";
						//D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
						//std::cout << "RecalCrhoL=" << rhoL << " sumRho=" << 27.0 / 18.0 * sumRho2 << " vx=" << vx << " vy=" << vy << "ffRest="<<ff[DIR_000]<<"\n";
//						distribution->setDistributionForDirection(ff[DIR_000], x1, x2, x3, DIR_000);
						//{
						//	real fG = distribution->getDistributionInvForDirection(x1, x2, x3, DIR_000);
						//	real feqOLD = D3Q27System::getIncompFeqForDirection(DIR_000, (*rhoNode)(x1, x2, x3), vx, vy, vz);
						//	real feqNew = D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz);
						//	distribution->setDistributionForDirection(fG - feqOLD + feqNew, x1, x2, x3, DIR_000);
						//}
                        //ff[DIR_000] = vx * vx + vy * vy + vz * vz +
                        //              (((ff[DIR_MM0] + ff[DIR_PP0]) + (ff[DIR_MP0] + ff[DIR_PM0])) + ((ff[DIR_0MM] + ff[DIR_0PP]) + (ff[DIR_0MP] + ff[DIR_0PM])) + ((ff[DIR_M0M] + ff[DIR_P0P]) + (ff[DIR_M0P] + ff[DIR_P0M])) +
                        //               c2o1 * ((((ff[DIR_MMM] + ff[DIR_PPP]) + (ff[DIR_MMP] + ff[DIR_PPM]))) + (((ff[DIR_MPM] + ff[DIR_PMP]) + (ff[DIR_MPP] + ff[DIR_PMM])))));
                        //distribution->setDistributionForDirection(ff[DIR_000], x1, x2, x3, DIR_000);

                        //for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//	ff[D3Q27System::INVDIR[fdir]]=distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//}
						//D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
						//std::cout << "AfterRead rhoL=" << rhoL << " rhoGToL=" << rhoG/densityRatio << " vx=" << vx << " vy=" << vy << "ffRest=" << ff[DIR_000] <<" x="<<x1<<" y="<<x2<<" z="<<x3<< "\n";

								//real feqL = D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz);
								//distribution->setDistributionForDirection(feqL, x1, x2, x3, DIR_000);



							}



						}


						//if ((*phaseField)(x1, x2, x3) <= c1o2) {
						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//		
						//		real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, 0, 0.0001, 0);
						//		ff[D3Q27System::INVDIR[fdir]] = feqNew;
						//		distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//	}
						//}
						//16.03.23 B: Bounce Back gas side
						//distribution->getDistributionInv(ff, x1, x2, x3);
						//real rhoG;
						//if (phi[DIR_000] > c1o2) { //initialization necessary
						//	real sumRho = 0;
						//	real sumWeight = 0;
						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//		if ((phi[fdir] <= c1o2)) {
						//			sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
						//			sumWeight += WEIGTH[fdir];
						//		}

						//	}
						//	rhoG = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neigbor?
						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//		if ((phi[fdir] > c1o2)) {
						//			real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//			real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//			real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//			real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//			real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
						//			real fBC = fG - c6o1 * WEIGTH[fdir] * (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX2[fdir] * vzBC);

						//			distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);

						//		}
						//	}
						//	distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1, x2, x3, DIR_000);



						//}
						//else {//no refill of gas required
						//	rhoG = (*rhoNode)(x1, x2, x3);
						//	if ((*phaseField)(x1, x2, x3) <= c1o2) {//no refill liquid
						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {
						//				real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
						//				real fBC = fG - c6o1 * WEIGTH[fdir] * (D3Q27System::DX1[fdir] * vxBC + D3Q27System::DX2[fdir] * vyBC + D3Q27System::DX2[fdir] * vzBC);

						//				//if ((*phaseField)(x1, x2, x3) <= c1o2) 
						//				distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
						//					//real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
						//					real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

						//					distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
						//				}

						//			}


						//		}
						//	}
						//	else {//refill liquid

						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {
						//				real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
						//				real fBC = fG-c6o1*WEIGTH[fdir]*(D3Q27System::DX1[fdir]*vxBC+ D3Q27System::DX2[fdir] * vyBC+ D3Q27System::DX2[fdir] * vzBC);

						//				ff[D3Q27System::INVDIR[fdir]] = fBC;
						//				if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
						//					//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
						//					real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

						//					distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
						//				}

						//			}
						//			else {
						//				ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);;
						//			}


						//		}

						//		real sum2 = 1e-100;
						//		real sumRho = 0;
						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {

						//				sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);// * tempRho;
						//				sum2 += WEIGTH[fdir];
						//			}
						//		}
						//		real rhoL;
						//		D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
						//		rhoL = sumRho / sum2;


						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[D3Q27System::INVDIR[fdir]] <= c1o2) && (phi[fdir] > c1o2)) {
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]) + c6o1 * WEIGTH[fdir] * (vx * D3Q27System::DX1[fdir] + vy * D3Q27System::DX2[fdir] + vz * D3Q27System::DX3[fdir]);;
						//				distribution->setDistributionForDirection(fL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
						//			}
						//			if ((phi[D3Q27System::INVDIR[fdir]] <= c1o2) && (phi[fdir] <= c1o2)) {
						//				real feqL = D3Q27System::getIncompFeqForDirection(fdir, rhoL, vx, vy, vz);
						//				distribution->setDistributionForDirection(feqL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
						//			}
						//		}

						//		real feqL = D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz);
						//		distribution->setDistributionForDirection(feqL, x1, x2, x3, DIR_000);



						//	}



						//}




						//16.03.23 A: scaled pressure
						//distribution->getDistributionInv(ff, x1, x2, x3);
						//real rhoG;
						//if (phi[DIR_000] > c1o2) { //initialization necessary
						//	real sumRho = 0;
						//	real sumWeight = 0;
						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//		if ((phi[fdir] <= c1o2)) {
						//			sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
						//			sumWeight += WEIGTH[fdir];
						//		}

						//	}
						//	rhoG = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neigbor?
						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//		if ((phi[fdir] > c1o2)) {
						//			real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//			real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//			real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//			real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

						//			distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);

						//		}
						//	}
						//	distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1, x2, x3, DIR_000);



						//}
						//else {//no refill of gas required
						//	rhoG = (*rhoNode)(x1, x2, x3);
						//	if ((*phaseField)(x1, x2, x3) <= c1o2) {//no refill liquid
						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
						//				real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

						//				//if ((*phaseField)(x1, x2, x3) <= c1o2) 
						//					distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
						//					real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
						//					real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

						//					distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
						//				}

						//			}


						//		}
						//	}
						//	else {//refill liquid

						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
						//				real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//				real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
						//				real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

						//				ff[D3Q27System::INVDIR[fdir]] = fBC;
						//				if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
						//					real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
						//					real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
						//					//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

						//					distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
						//				}

						//			}
						//			else {
						//				ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);;
						//			}


						//		}

						//		real sum2 = 1e-100;
						//		real sumRho = 0;
						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[fdir] > c1o2)) {
						//				
						//					sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);// * tempRho;
						//					sum2 += WEIGTH[fdir];									
						//			}
						//		}
						//		real rhoL;
						//		D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
						//		rhoL=sumRho/sum2;


						//		for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
						//			if ((phi[D3Q27System::INVDIR[fdir]] <= c1o2) && (phi[fdir] > c1o2)) {
						//				real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]) + c6o1 * WEIGTH[fdir] * (vx * D3Q27System::DX1[fdir] + vy * D3Q27System::DX2[fdir] + vz * D3Q27System::DX3[fdir]);;
						//				distribution->setDistributionForDirection(fL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
						//			}
						//			if ((phi[D3Q27System::INVDIR[fdir]] <= c1o2) && (phi[fdir] <= c1o2)) {
						//				real feqL = D3Q27System::getIncompFeqForDirection(fdir, rhoL, vx, vy, vz);
						//				distribution->setDistributionForDirection(feqL, x1 - D3Q27System::DX1[fdir], x2 - D3Q27System::DX2[fdir], x3 - D3Q27System::DX3[fdir], fdir);
						//			}
						//		}

						//		real feqL = D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz);
						//		distribution->setDistributionForDirection(feqL, x1 , x2, x3 , DIR_000);



						//	}



						//}










					
}//end Loop
					




	//for (int x3 = minX3-1; x3 < maxX3+1; x3++) {
	//	for (int x2 = minX2-1; x2 < maxX2+1; x2++) {
	//		for (int x1 = minX1-1; x1 < maxX1+1; x1++) {
	//			if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
	//				int x1p = x1 + 1;
	//				int x2p = x2 + 1;
	//				int x3p = x3 + 1;
	//				findNeighbors(phaseFieldOld, x1, x2, x3);
	//				////////////////////////////////Momentum conservation experiment 06.03.2023
	//				//surfacetension
	//				real  kapkap = 0*1.0e-5;
	//				//real scalRefill = 0.0;
	//				real slowerFactor = 1.0e6;
	//				if (((*phaseField)(x1, x2, x3) <= c1o2) && (
	//					(phi[DIR_P00] > c1o2) ||
	//					(phi[DIR_M00] > c1o2) ||
	//					(phi[DIR_00P] > c1o2) ||
	//					(phi[DIR_00M] > c1o2) ||
	//					(phi[DIR_0M0] > c1o2) ||
	//					(phi[DIR_0P0] > c1o2) ||
	//					(phi[DIR_PP0] > c1o2) ||
	//					(phi[DIR_PM0] > c1o2) ||
	//					(phi[DIR_P0P] > c1o2) ||
	//					(phi[DIR_P0M] > c1o2) ||
	//					(phi[DIR_MP0] > c1o2) ||
	//					(phi[DIR_MM0] > c1o2) ||
	//					(phi[DIR_M0P] > c1o2) ||
	//					(phi[DIR_M0M] > c1o2) ||
	//					(phi[DIR_0PM] > c1o2) ||
	//					(phi[DIR_0MM] > c1o2) ||
	//					(phi[DIR_0PP] > c1o2) ||
	//					(phi[DIR_0MP] > c1o2) ||
	//					(phi[DIR_PPP] > c1o2) ||
	//					(phi[DIR_PMP] > c1o2) ||
	//					(phi[DIR_MPP] > c1o2) ||
	//					(phi[DIR_MMP] > c1o2) ||
	//					(phi[DIR_PPM] > c1o2) ||
	//					(phi[DIR_PMM] > c1o2) ||
	//					(phi[DIR_MPM] > c1o2) ||
	//					(phi[DIR_MMM] > c1o2)
	//					)) {
	//						real vx = (*vxNode)(x1, x2, x3);
	//						real vy =  (*vyNode)(x1, x2, x3);
	//						real vz = (*vzNode)(x1, x2, x3);


	//						distribution->getDistributionInv(ff, x1, x2, x3);
	//						real rhoG;
	//						if (phi[DIR_000] > c1o2) { //initialization necessary
	//							real sumRho = 0;
	//							real sumWeight = 0;
	//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
	//								if ((phi[fdir] <= c1o2)) {
	//									sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
	//									sumWeight += WEIGTH[fdir];
	//								}

	//							}
	//							rhoG = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neigbor?
	//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
	//								if ((phi[fdir] > c1o2)) {
	//									real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//									real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//									real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//									real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

	//									distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);

	//								}
	//							}
	//							distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1, x2, x3, DIR_000);



	//						}
	//						else {//no refill required

	//							rhoG = (*rhoNode)(x1, x2, x3);
	//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
	//								if ((phi[fdir] > c1o2)) {
	//									real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//									real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//									real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//									real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//									real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;

	//									distribution->setDistributionForDirection(fBC, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//									if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
	//										real vxBC =c1o2*(vx+ (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										real vyBC =c1o2*(vy+ (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										real vzBC =c1o2*(vz+ (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
	//										real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])*(D3Q27System::DX1[fdir])* (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//										real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//										//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

	//										distribution->setDistributionForDirection((fBC+fG) / densityRatio-fL+kapkap* WEIGTH[fdir]* computeCurvature_phi() -(feqG-feqL)*(c1o1/densityRatio-c1o1)*(vxBC* D3Q27System::DX1[fdir]+vyBC* D3Q27System::DX2[fdir]+vzBC* D3Q27System::DX3[fdir]), x1 , x2 , x3 , fdir);
	//									}

	//									}
	//								else {
	//									if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {

	//										real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//										real feqOLD = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy,vz);
	//										real slower = c1o1/(c1o1+slowerFactor * (vx * vx + vy * vy + vz * vz));
	//										real feqNew = D3Q27System::getIncompFeqForDirection(fdir, rhoG / densityRatio, slower * vx, slower * vy, slower * vz);
	//										real fBC = (fG - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew;

	//										distribution->setDistributionForDirection(fBC, x1, x2 , x3 , fdir);

	//										////inverse refill from here
	//										//int xn1 = x1 + D3Q27System::DX1[fdir];
	//										//int xn2 = x2 + D3Q27System::DX2[fdir];
	//										//int xn3 = x3 + D3Q27System::DX3[fdir];
	//										//real sumRho = 0;
	//										//real sumWeight = 0;
	//										//for (int nfdir = D3Q27System::STARTF; nfdir < D3Q27System::ENDF; nfdir++) {
	//										//	if ((phi[nfdir] > c1o2)) {
	//										//		sumRho += WEIGTH[nfdir] * (*rhoNode)(xn1 + D3Q27System::DX1[nfdir], xn2 + D3Q27System::DX2[nfdir], xn3 + D3Q27System::DX3[nfdir]);
	//										//		sumWeight += WEIGTH[nfdir];
	//										//	}
	//										//}
	//										////real rhoL = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neigbor?
	//										//real rhoL = c1o2*(sumRho / sumWeight * scalRefill + (c1o1 - scalRefill) * rhoG / densityRatio);//

	//										//// what comes next is the inversion of BC for the gas phase which is only used to derive the liquid BC
	//										//real fBC = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//										////real feqOld = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoL, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										////Dirty
	//										//real feqOld = D3Q27System::getIncompFeqForDirection(fdir, rhoL, -(*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), -(*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), -(*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

	//										//real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										//real fL = (fBC - feqNew) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqOld;


	//										//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//										//real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										//real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										//real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										////real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//										////real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
	//										//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//										//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//										////real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//										//distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL + kapkap * WEIGTH[fdir] * computeCurvature_phi() - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
	//										////distribution->setDistributionForDirection(( fG) / densityRatio - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
	//									
	//									
	//									}
	//								}
	//								
	//							}

	//						}






	//						

	//				
	//				
	//				}
	//				if (((*phaseField)(x1, x2, x3) > c1o2) && ((*phaseFieldOld)(x1, x2, x3) <= c1o2)) {
	//					real vx = (*vxNode)(x1, x2, x3);
	//					real vy = (*vyNode)(x1, x2, x3);
	//					real vz = (*vzNode)(x1, x2, x3);
	//					for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
	//						if (((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) > c1o2) {
	//							real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//							real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//							real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//							real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//							//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//							//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
	//							real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//							real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

	//							if (((*phaseFieldOld)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) <= c1o2) {
	//								real rhoG = (*rhoNode)(x1, x2, x3);
	//								real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//								real feqOLD = D3Q27System::getIncompFeqForDirection(fdir, rhoG, vx, vy, vz);
	//								real slower = c1o1 / (c1o1 + slowerFactor * (vx * vx + vy * vy + vz * vz));
	//								real feqNew = D3Q27System::getIncompFeqForDirection(fdir, rhoG / densityRatio, slower*vx, slower*vy, slower*vz);
	//								real fBC = (fG - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew;

	//								distribution->setDistributionForDirection(fBC, x1, x2, x3, fdir);

	//								/////reverse liquid
	//								//int xn1 = x1 + D3Q27System::DX1[fdir];
	//								//int xn2 = x2 + D3Q27System::DX2[fdir];
	//								//int xn3 = x3 + D3Q27System::DX3[fdir];
	//								//real sumRho = 0;
	//								//real sumWeight = 0;
	//								//for (int nfdir = D3Q27System::STARTF; nfdir < D3Q27System::ENDF; nfdir++) {
	//								//	if ((phi[nfdir] > c1o2)) {
	//								//		sumRho += WEIGTH[nfdir] * (*rhoNode)(xn1 + D3Q27System::DX1[nfdir], xn2 + D3Q27System::DX2[nfdir], xn3 + D3Q27System::DX3[nfdir]);
	//								//		sumWeight += WEIGTH[nfdir];
	//								//	}
	//								//}
	//								////real rhoL = sumRho / sumWeight;// uncheck excpetion: what if there is no adequate neigbor?
	//								//real rhoL = (sumRho / sumWeight*scalRefill+(c1o1-scalRefill)*(*rhoNode)(x1, x2, x3) / densityRatio);//
	//								//// what comes next is the inversion of BC for the gas phase which is only used to derive the liquid BC
	//								//real fBC = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//								////real feqOld = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoL, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								////Dirty
	//								//real feqOld = D3Q27System::getIncompFeqForDirection(fdir, rhoL, -(*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), -(*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), -(*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));

	//								//real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1, x2, x3), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								//real fL = (fBC - feqNew) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqOld;


	//								//real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//								//real vxBC = c1o2 * (vx + (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								//real vyBC = c1o2 * (vy + (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								//real vzBC = c1o2 * (vz + (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								////real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								////real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
	//								//real feqL = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//								//real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));
	//								////real feqG = D3Q27System::getIncompFeqForDirection(fdir, 0, vx * (D3Q27System::DX1[fdir]) * (D3Q27System::DX1[fdir]), vy * (D3Q27System::DX2[fdir]) * (D3Q27System::DX2[fdir]), vz * (D3Q27System::DX3[fdir]) * (D3Q27System::DX3[fdir]));

	//								//distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL + kapkap * WEIGTH[fdir] * computeCurvature_phi() - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);

	//								/////!reverse liquid
	//								//
	//								////distribution->setDistributionForDirection((fG) / densityRatio - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
	//							}
	//							else {
	//								real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
	//								real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1, x2, x3), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
	//								real fG = distribution->getDistributionInvForDirection(x1, x2, x3, fdir);
	//								real fBC = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;
	//								distribution->setDistributionForDirection((fBC + fG) / densityRatio - fL + kapkap * WEIGTH[fdir] * computeCurvature_phi() - (feqG - feqL) * (c1o1 / densityRatio - c1o1) * (vxBC * D3Q27System::DX1[fdir] + vyBC * D3Q27System::DX2[fdir] + vzBC * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
	//							}
	//						}
	//					
	//					}



	//				}





					//////////////////////////////////////

//					//if ((phi[DIR_000] > c1o2) && (
//					//	(phi[DIR_P00] <= c1o2) ||
//					//	(phi[DIR_M00] <= c1o2) ||
//					//	(phi[DIR_00P] <= c1o2) ||
//					//	(phi[DIR_00M] <= c1o2) ||
//					//	(phi[DIR_0M0] <= c1o2) ||
//					//	(phi[DIR_0P0] <= c1o2) ||
//					//	(phi[DIR_PP0] <= c1o2) ||
//					//	(phi[DIR_PM0] <= c1o2) ||
//					//	(phi[DIR_P0P] <= c1o2) ||
//					//	(phi[DIR_P0M] <= c1o2) ||
//					//	(phi[DIR_MP0] <= c1o2) ||
//					//	(phi[DIR_MM0] <= c1o2) ||
//					//	(phi[DIR_M0P] <= c1o2) ||
//					//	(phi[DIR_M0M] <= c1o2) ||
//					//	(phi[DIR_0PM] <= c1o2) ||
//					//	(phi[DIR_0MM] <= c1o2) ||
//					//	(phi[DIR_0PP] <= c1o2) ||
//					//	(phi[DIR_0MP] <= c1o2) ||
//					//	(phi[DIR_PPP] <= c1o2) ||
//					//	(phi[DIR_PMP] <= c1o2) ||
//					//	(phi[DIR_MPP] <= c1o2) ||
//					//	(phi[DIR_MMP] <= c1o2) ||
//					//	(phi[DIR_PPM] <= c1o2) ||
//					//	(phi[DIR_PMM] <= c1o2) ||
//					//	(phi[DIR_MPM] <= c1o2) ||
//					//	(phi[DIR_MMM] <= c1o2)
//					//	)) {
//
//					//	real vx = (*vxNode)(x1, x2, x3);
//					//	real vy =  (*vyNode)(x1, x2, x3);
//					//	real vz = (*vzNode)(x1, x2, x3);
//
//
//					//	distribution->getDistributionInv(ff, x1, x2, x3);
//					//	real rhoG;
//					//	//D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//					//	real sumRhoG = 0.0;
//					//	int countRhoG = 0;
//					//	for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {
//					//		if ((phi[fdir] <= c1o2)) {
//					//			//BC version
//					//			// rhoG =  (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//					//			//real ftemp = D3Q27System::getCompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz) + D3Q27System::getCompFeqForDirection(fdir, rhoG, vx, vy, vz);
//					//			//
//					//			//real fBB;
//					//			//fBB = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//					//			//distribution->setDistributionForDirection((ftemp - ff[fdir]), x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//					//			//distribution->setDistributionForDirection(fBB - c6o1 * D3Q27System::WEIGTH[fdir] * (-vx * D3Q27System::DX1[fdir] - vy * D3Q27System::DX2[fdir] - vz * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
//					//		//scaled Version
//
//					//			real fG;
//					//			fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//
//					//			//Liquid
//					//			real ssrho = 1;
//					//			real rhoLL = (*rhoNode)(x1, x2, x3);
//
//					//			//real rhoFilter = (*rhoNode)(x1, x2, x3)*c8o27
//					//			//	+ c2o27*(((*rhoNode)(x1 + 1, x2, x3) + (*rhoNode)(x1 - 1, x2, x3)) + ((*rhoNode)(x1, x2 + 1, x3) + (*rhoNode)(x1, x2 - 1, x3)) + ((*rhoNode)(x1, x2, x3 + 1) + (*rhoNode)(x1, x2, x3 - 1)))
//					//			//	+ c1o54*((((*rhoNode)(x1 + 1, x2 + 1, x3) + (*rhoNode)(x1 - 1, x2 - 1, x3)) + ((*rhoNode)(x1 - 1, x2 + 1, x3) + (*rhoNode)(x1 + 1, x2 - 1, x3)))
//					//			//		+ (((*rhoNode)(x1 + 1, x2, x3 + 1) + (*rhoNode)(x1 - 1, x2, x3 - 1)) + ((*rhoNode)(x1 - 1, x2, x3 + 1) + (*rhoNode)(x1 + 1, x2, x3 - 1)))
//					//			//		+ (((*rhoNode)(x1, x2 + 1, x3 + 1) + (*rhoNode)(x1, x2 - 1, x3 - 1)) + ((*rhoNode)(x1, x2 - 1, x3 + 1) + (*rhoNode)(x1, x2 + 1, x3 - 1)))
//					//			//		)
//					//			//	+ c1o216*(
//					//			//		(((*rhoNode)(x1 + 1, x2 + 1, x3 + 1) + (*rhoNode)(x1 - 1, x2 - 1, x3 - 1)) + ((*rhoNode)(x1 + 1, x2 - 1, x3 + 1) + (*rhoNode)(x1 - 1, x2 + 1, x3 - 1)))
//					//			//		+ (((*rhoNode)(x1 + 1, x2 + 1, x3 - 1) + (*rhoNode)(x1 - 1, x2 - 1, x3 + 1)) + ((*rhoNode)(x1 + 1, x2 - 1, x3 - 1) + (*rhoNode)(x1 - 1, x2 + 1, x3 + 1)))
//					//			//		);
//					//			real rhoGG = (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//					//			real feqOLD = (D3Q27System::getIncompFeqForDirection(fdir, rhoLL/densityRatio, vx, vy, vz));
//					//			real feqNew = (D3Q27System::getIncompFeqForDirection(fdir, rhoLL*(c1o1-ssrho)+ssrho*rhoGG, vx, vy, vz));
//					//			//real feqNew = (D3Q27System::getIncompFeqForDirection(fdir, rhoFilter, vx, vy, vz));
//					//			distribution->setDistributionForDirection( (ff[fdir] - feqOLD)*(c1o1/collFactorG-c1o1)/(c1o1/collFactorL-c1o1) + feqNew, x1, x2, x3, fdir);
//
//					//			feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//					//			feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) / densityRatio, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//					//			distribution->setDistributionForDirection((fG - feqOLD)* (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//					//			sumRhoG += (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//					//			countRhoG++;
//
//
//					//		}
//					//	}
//					//	(*rhoNode)(x1, x2, x3) = sumRhoG / countRhoG;
//
//					if ((phi[DIR_000] > c1o2) && (
//						(phi[DIR_P00] <= c1o2) ||
//						(phi[DIR_M00] <= c1o2) ||
//						(phi[DIR_00P] <= c1o2) ||
//						(phi[DIR_00M] <= c1o2) ||
//						(phi[DIR_0M0] <= c1o2) ||
//						(phi[DIR_0P0] <= c1o2) ||
//						(phi[DIR_PP0] <= c1o2) ||
//						(phi[DIR_PM0] <= c1o2) ||
//						(phi[DIR_P0P] <= c1o2) ||
//						(phi[DIR_P0M] <= c1o2) ||
//						(phi[DIR_MP0] <= c1o2) ||
//						(phi[DIR_MM0] <= c1o2) ||
//						(phi[DIR_M0P] <= c1o2) ||
//						(phi[DIR_M0M] <= c1o2) ||
//						(phi[DIR_0PM] <= c1o2) ||
//						(phi[DIR_0MM] <= c1o2) ||
//						(phi[DIR_0PP] <= c1o2) ||
//						(phi[DIR_0MP] <= c1o2) ||
//						(phi[DIR_PPP] <= c1o2) ||
//						(phi[DIR_PMP] <= c1o2) ||
//						(phi[DIR_MPP] <= c1o2) ||
//						(phi[DIR_MMP] <= c1o2) ||
//						(phi[DIR_PPM] <= c1o2) ||
//						(phi[DIR_PMM] <= c1o2) ||
//						(phi[DIR_MPM] <= c1o2) ||
//						(phi[DIR_MMM] <= c1o2)
//						)) {
//							real vx = (*vxNode)(x1, x2, x3);
//							real vy =  (*vyNode)(x1, x2, x3);
//							real vz = (*vzNode)(x1, x2, x3);
//
//
//						//distribution->getDistributionInv(ff, x1, x2, x3);
//
//						if ((*phaseField)(x1, x2, x3) > c1o2) {
//						
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {// populations without DIR_000
//								if ((phi[fdir] <= c1o2)) {
//										real fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//										real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//										real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) / densityRatio, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//										distribution->setDistributionForDirection((fG - feqOLD)* (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//
//								}
//							}
//
//
//						}
//						else {
//						//refill necessary
//							real sumRho = 0;
//							real sumWeight = 0;
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								if ((phi[fdir] > c1o2)) {
//									sumRho += WEIGTH[fdir]*(*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//									sumWeight += WEIGTH[fdir];
//									}
//							}
//							sumRho /= sumWeight;
//
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								if ((phi[fdir] > c1o2)) {
//									real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//									real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//									real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], sumRho, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//									ff[D3Q27System::INVDIR[fdir]] = (fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew;
//								}
//								else { ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//								}
//							}
//
//							real rhoG;
//							D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//							sumRho = 0;
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								sumRho = ff[fdir] - D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
//							}
//							rhoG = 27.0 / 19.0 * sumRho;
//							distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1 , x2 , x3 , DIR_000);
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								//if ((phi[fdir] > c1o2)) {
//									distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//								//}
//							}
//
//
//						}
//
//
//					}
//					else if ((phi[DIR_000] <= c1o2) && (
//						(phi[DIR_P00] > c1o2) ||
//						(phi[DIR_M00] > c1o2) ||
//						(phi[DIR_00P] > c1o2) ||
//						(phi[DIR_00M] > c1o2) ||
//						(phi[DIR_0M0] > c1o2) ||
//						(phi[DIR_0P0] > c1o2) ||
//						(phi[DIR_PP0] > c1o2) ||
//						(phi[DIR_PM0] > c1o2) ||
//						(phi[DIR_P0P] > c1o2) ||
//						(phi[DIR_P0M] > c1o2) ||
//						(phi[DIR_MP0] > c1o2) ||
//						(phi[DIR_MM0] > c1o2) ||
//						(phi[DIR_M0P] > c1o2) ||
//						(phi[DIR_M0M] > c1o2) ||
//						(phi[DIR_0PM] > c1o2) ||
//						(phi[DIR_0MM] > c1o2) ||
//						(phi[DIR_0PP] > c1o2) ||
//						(phi[DIR_0MP] > c1o2) ||
//						(phi[DIR_PPP] > c1o2) ||
//						(phi[DIR_PMP] > c1o2) ||
//						(phi[DIR_MPP] > c1o2) ||
//						(phi[DIR_MMP] > c1o2) ||
//						(phi[DIR_PPM] > c1o2) ||
//						(phi[DIR_PMM] > c1o2) ||
//						(phi[DIR_MPM] > c1o2) ||
//						(phi[DIR_MMM] > c1o2)
//						)) {
//						real vx = (*vxNode)(x1, x2, x3);
//						real vy = (*vyNode)(x1, x2, x3);
//						real vz = (*vzNode)(x1, x2, x3);
//
//
//						//distribution->getDistributionInv(ff, x1, x2, x3);
//						if ((*phaseField)(x1, x2, x3) <= c1o2) {
//						////explicit way:
//						////real ppph = (*phaseField)(x1, x2, x3);
//						//	for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//						//		if ((phi[fdir] > c1o2)) {
//						//			//vx = (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//						//			//vy = (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//						//			//vz = (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//						//			//real rhorho = (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//						//			//int xx1 = x1 + D3Q27System::DX1[fdir];
//						//			//int xx2 = x2 + D3Q27System::DX2[fdir];
//						//			//int xx3 = x3 + D3Q27System::DX3[fdir];
//
//						//			real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//						//			real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//						//			real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 , x2 , x3 ), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//						//			distribution->setDistributionForDirection((fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//						//		}
//						//	}
///////iterative way:
//							real rhoG = (*rhoNode)(x1, x2, x3);
//							//real sumWeight=0;
//							//for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//							//	if ((phi[fdir] > c1o2)) {
//							//		rhoG += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//							//		sumWeight += WEIGTH[fdir];
//							//	}
//							//}
//							//rhoG = rhoG/sumWeight*densityRatio;
//
//							for (int itter = 0; itter < 5; itter++) {
//								for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//									if ((phi[fdir] > c1o2)) {
//										real fL = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//										real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//										real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//										ff[D3Q27System::INVDIR[fdir]] = ((fL - feqOLD) * (c1o1 / collFactorG - c1o1) / (c1o1 / collFactorL - c1o1) + feqNew);
//									}
//									else { ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]); }
//								}
//								ff[DIR_000]= distribution->getDistributionInvForDirection(x1, x2, x3, DIR_000);
//								D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//								//real sumRho = 0;
//								//for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								//	sumRho = ff[fdir] - D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
//								//}
//								//rhoG = 27.0 / 19.0 * sumRho;
//							}
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								
//								distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//								
//							}
//							//distribution->setDistributionForDirection(ff[DIR_000], x1, x2, x3, DIR_000);
//
//
//
//
//						}
//						else {
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								if ((phi[fdir] <= c1o2)) {
//									real fG = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//									real feqOLD = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//									real feqNew = D3Q27System::getIncompFeqForDirection(D3Q27System::INVDIR[fdir], (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) / densityRatio, (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]), (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]));
//									ff[D3Q27System::INVDIR[fdir]] = (fG - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew;
//									distribution->setDistributionForDirection(ff[D3Q27System::INVDIR[fdir]], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);
//								}
//								else { ff[D3Q27System::INVDIR[fdir]] = distribution->getDistributionInvForDirection(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], D3Q27System::INVDIR[fdir]);}
//							}
//							real rhoG;
//							D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//							real sumRho = 0;
//							for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//								sumRho = ff[fdir] - D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
//							}
//							rhoG = 27.0 / 19.0 * sumRho;
//							distribution->setDistributionForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoG, vx, vy, vz), x1, x2, x3, DIR_000);
//
//
//
//						}
//
//					}
//
//
					}
				}
			}
		}
	


	this->swapDistributions();

	for (int x3 = minX3; x3 < maxX3; x3++) {
	for (int x2 = minX2; x2 < maxX2; x2++) {
		for (int x1 = minX1; x1 < maxX1; x1++) {
			 
				int x1p = x1 + 1;
				int x2p = x2 + 1;
				int x3p = x3 + 1;
				findNeighbors(phaseFieldOld, x1, x2, x3);

				//if (((*phaseField)(x1, x2, x3) > c1o2) && (((*phaseFieldOld)(x1, x2, x3) <= c1o2)))
				{//Refill liquid
					real vx;
					real vy;
					real vz;


					distribution->getDistribution(ff, x1, x2, x3);
					real rhoL;
					D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
					//if (vz != 0) {

					//	std::cout << "precol: rhoL=" << rhoL << " vx=" << vx << " vy=" << vy << " vz=" << vz << "ffRest=" << ff[DIR_000] << " x=" << x1 << " y=" << x2 << " z=" << x3 << "\n";
					//}
				}
			
		}
	}
}



	////////momentum balance 06.03.2023
	//for (int x3 = minX3 - ghostLayerWidth; x3 < maxX3 + ghostLayerWidth; x3++) {
	//	for (int x2 = minX2 - ghostLayerWidth; x2 < maxX2 + ghostLayerWidth; x2++) {
	//		for (int x1 = minX1 - ghostLayerWidth; x1 < maxX1 + ghostLayerWidth; x1++) {
	//			if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
	//				int x1p = x1 + 1;
	//				int x2p = x2 + 1;
	//				int x3p = x3 + 1;
	//				if (((*phaseField)(x1, x2, x3) > c1o2) && (((*phaseFieldOld)(x1, x2, x3) <= c1o2)))
	//				{//Refill liquid
	//					real vx;
	//					real vy;
	//					real vz;


	//					distribution->getDistribution(ff, x1, x2, x3);
	//					real rhoL;
	//					D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
	//					std::cout << "precol: rhoL=" << rhoL << " vx=" << vx << " vy=" << vy << "ffRest=" << ff[DIR_000] << " x=" << x1 << " y=" << x2 << " z=" << x3 << "\n";
	//				}
	//			}
	//		}
	//	}
	//}


	//						real sumRho = 0;
	//						for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
	//							sumRho = ff[fdir] - D3Q27System::getIncompFeqForDirection(fdir, 0, vx, vy, vz);
	//						}
	//						rhoL = 27.0 / 19.0 * sumRho;
	//						distribution->setDistributionInvForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz), x1, x2, x3, DIR_000);


	//				}
	//			}
	//		}
	//	}
	//}


	//////rescaling new liquid nodes 10.03.2023
//	for (int x3 = minX3 ; x3 < maxX3 ; x3++) {
//		for (int x2 = minX2 ; x2 < maxX2; x2++) {
//			for (int x1 = minX1 ; x1 < maxX1 ; x1++) {
//				if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
//					int x1p = x1 + 1;
//					int x2p = x2 + 1;
//					int x3p = x3 + 1;
//					if (((*phaseField)(x1, x2, x3) > c1o2) && (((*phaseFieldOld)(x1, x2, x3) <= c1o2)))
//					{//Refill liquid
//						real vx;
//						real vy;
//						real vz;
//
//						findNeighbors(phaseFieldOld, x1, x2, x3);
//
//						//real rhoG;
//						//D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//
//						//vx = (*vxNode)(x1, x2, x3);
//						//vy = (*vyNode)(x1, x2, x3);
//						//vz = (*vzNode)(x1, x2, x3);
//
//
//						//for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {// loop includes DIR_000 position, different from all the others
//						//	real feqOLD = D3Q27System::getIncompFeqForDirection(fdir, rhoG,vx,vy,vz);
//						//	real feqNew = D3Q27System::getIncompFeqForDirection(fdir, rhoG/densityRatio,vx,vy,vz);
//						//	real fBC = (ff[fdir] - feqOLD) * (c1o1 / collFactorL - c1o1) / (c1o1 / collFactorG - c1o1) + feqNew;
//
//						//	distribution->setDistributionInvForDirection(fBC, x1 , x2 , x3 ,fdir);
//
//
//
//						//}
////15.03.2023
//						real sumVx=0, sumVy=0, sumVz=0;
//						real tempRho, tempVx, tempVy, tempVz;
//						real sumRho = 0;
//						real sum = 1e-100;
//						real sum2 = 1e-100;
//						for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//							if (!(((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) > c1o2) && (((*phaseFieldOld)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) <= c1o2)))&& !bcArray->isSolid(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) && !bcArray->isUndefined(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir])) {
//								//distribution->getDistribution(ff, x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
//								//D3Q27System::calcIncompMacroscopicValues(ff, tempRho, tempVx, tempVy, tempVz);
//								sum += WEIGTH[fdir];
//								sumVx += WEIGTH[fdir] * (*vxNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);//*tempVx;
//								sumVy += WEIGTH[fdir] * (*vyNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);//*tempVy;
//								sumVz += WEIGTH[fdir] * (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);//*tempVz;
//								if ((*phaseField)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) > c1o2) {
//									sumRho += WEIGTH[fdir] * (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);// * tempRho;
//									sum2 += WEIGTH[fdir];
//								}
//								if (tempVz != 0) {
//									std::cout << "vz=" << tempVz << " " << "x=" << x1 << " " << "y=" << x2 << " " << "z=" << x3 << " fdir=" << fdir << " " << "xn=" << x1 + D3Q27System::DX1[fdir] << " " << "yn=" << x2 + D3Q27System::DX2[fdir] << " " << "zn=" << x3 + D3Q27System::DX3[fdir]<<"vzold="<< (*vzNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]) << "\n";
//								}
//							}
//						}
//						sumRho/=sum2 ;
//						sumVx /= sum;
//						sumVy /= sum;
//						sumVz /= sum;
//						distribution->setDistributionInvForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, sumRho, sumVx, sumVy, sumVz), x1, x2, x3, DIR_000);
//
//						std::cout << "x=" << x1 << " " << "y=" << x2 << " " << "z=" << x3 <<" sumVx="<<sumVx<< " sumVy=" << sumVy << " sumVz=" << sumVz << " sumRho=" << sumRho << "\n";
//
////14.03.2023 
//						distribution->getDistribution(ff, x1, x2, x3);
//						real rhoG= (*rhoNode)(x1, x2, x3);
//
//						vx = (*vxNode)(x1, x2, x3);
//						vy = (*vyNode)(x1, x2, x3);
//						vz = (*vzNode)(x1, x2, x3);
//						std::cout << " Vx=" << vx << " Vy=" << vy << " Vz=" << vz << " rhoL=" << (*rhoNode)(x1, x2, x3) / densityRatio << "\n";
//
//						for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//							if (phi[D3Q27System::INVDIR[fdir]] <= c1o2 && phi[fdir] > c1o2) {
//								//ff[fdir] = ff[D3Q27System::INVDIR[fdir]] + c6o1 * WEIGTH[fdir] * (vx * D3Q27System::DX1[fdir] + vy * D3Q27System::DX2[fdir] + vz * D3Q27System::DX3[fdir]);
//								ff[fdir] = ff[D3Q27System::INVDIR[fdir]] + c6o1 * WEIGTH[fdir] * (sumVx * D3Q27System::DX1[fdir] + sumVy * D3Q27System::DX2[fdir] + sumVz * D3Q27System::DX3[fdir]);
//								//ff[fdir] = D3Q27System::getIncompFeqForDirection(fdir, sumRho, sumVx, sumVy, sumVz);
//								distribution->setDistributionInvForDirection(ff[fdir], x1, x2, x3, fdir);
//							}
//							if (phi[fdir] <= c1o2 && phi[D3Q27System::INVDIR[fdir]] <= c1o2) {
//								//ff[fdir] = D3Q27System::getIncompFeqForDirection(fdir, rhoG / densityRatio, vx, vy, vz);
//								ff[fdir]= D3Q27System::getIncompFeqForDirection(fdir, sumRho, sumVx, sumVy, sumVz);
//								distribution->setDistributionInvForDirection(ff[fdir], x1, x2, x3, fdir);
//							}
//						}
//					real rhoL;
//					//D3Q27System::calcIncompMacroscopicValues(ff, rhoL, vx, vy, vz);
//
//
//					//real sumRho;
//					real sumRho2= 0;
//						for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
//							sumRho2 += ff[fdir];// -D3Q27System::getIncompFeqForDirection(fdir, 0, sumVx, sumVy, sumVz);
//						}
//						rhoL = 27.0 / 19.0 * sumRho;
//						//distribution->setDistributionInvForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, vx, vy, vz), x1, x2, x3, DIR_000);
//						//distribution->setDistributionInvForDirection(D3Q27System::getIncompFeqForDirection(DIR_000, rhoL, sumVx, sumVy, sumVz), x1, x2, x3, DIR_000);
//						ff[DIR_000] = sumRho - sumRho2;
//						distribution->setDistributionInvForDirection(ff[DIR_000], x1, x2, x3, DIR_000);
//						D3Q27System::calcIncompMacroscopicValues(ff, rhoG, vx, vy, vz);
//						std::cout << " calcVx=" << vx << " calcVy=" << vy << " calcVz=" << vz << " rhoG=" << rhoG << "\n";
//
//
//
//					}
//				}
//			}
//		}
//	}


	//for (int x3 = minX3-ghostLayerWidth; x3 < maxX3+ghostLayerWidth; x3++) {
	//	for (int x2 = minX2-ghostLayerWidth; x2 < maxX2+ghostLayerWidth; x2++) {
	//		for (int x1 = minX1-ghostLayerWidth; x1 < maxX1+ghostLayerWidth; x1++) {
	//			if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
	//				int x1p = x1 + 1;
	//				int x2p = x2 + 1;
	//				int x3p = x3 + 1;



	//				real mfcbb = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);
	//				real mfbcb = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);
	//				real mfbbc = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);
	//				real mfccb = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);
	//				real mfacb = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);
	//				real mfcbc = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);
	//				real mfabc = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);
	//				real mfbcc = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);
	//				real mfbac = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);
	//				real mfccc = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);
	//				real mfacc = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);
	//				real mfcac = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);
	//				real mfaac = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);
	//				real mfabb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);
	//				real mfbab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);
	//				real mfbba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);
	//				real mfaab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);
	//				real mfcab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);
	//				real mfaba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);
	//				real mfcba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);
	//				real mfbaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);
	//				real mfbca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);
	//				real mfaaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);
	//				real mfcaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);
	//				real mfaca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);
	//				real mfcca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);
	//				real mfbbb = (*this->zeroDistributionsH1)(x1, x2, x3);

	//				omegaDRho = 2.0;// 1.5;
	//				real phiOld = (*phaseField)(x1, x2, x3);

	//				(*phaseField)(x1, x2, x3) = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca)  + (mfcaa + mfacc))  ) +
	//					(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
	//						((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
	//							(mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;
	//				//if (phiOld > 0.49 && phiOld < 0.501) {
	//				//	real ppppppppp = (*phaseField)(x1, x2, x3);
	//				//	int ist = 1;
	//				//}
	//				//if (phiOld > 0.5 && (*phaseField)(x1, x2, x3) <= 0.5) {
	//				//	real ppppppppp = (*phaseField)(x1, x2, x3);
	//				//	int ist = 1;
	//				//}


	//				if ((*phaseField)(x1, x2, x3) > 1 ) {
	//					(*phaseField)(x1, x2, x3) = c1o1;
	//				}

	//				if ((*phaseField)(x1, x2, x3) < 0) {
	//					(*phaseField)(x1, x2, x3) = 0;
	//				}
	//				////// read F-distributions for velocity formalism
	//				if (((phiOld <= 0.5) && ((*phaseField)(x1, x2, x3) <= 0.5)) || ((phiOld > 0.5) && ((*phaseField)(x1, x2, x3) > 0.5))) {}
	//				else {
	//					real scaleDistribution = densityRatio;// gas turn liquid
	//					real scaleStress = (c1o1/collFactorG)/(c1o1/collFactorL);
	//					if ((phiOld > 0.5) && ((*phaseField)(x1, x2, x3) <= 0.5)) {
	//						scaleDistribution = 1.0 / densityRatio;
	//						scaleStress = (c1o1 / collFactorL) / (c1o1 / collFactorG);
	//						//liquid turned gas
	//					}

	//					mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
	//					mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
	//					mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
	//					mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
	//					mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
	//					mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
	//					mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
	//					mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
	//					mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
	//					mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
	//					mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
	//					mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
	//					mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
	//					mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
	//					mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
	//					mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
	//					mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
	//					mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
	//					mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
	//					mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
	//					mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
	//					mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
	//					mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
	//					mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
	//					mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
	//					mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

	//					mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);

	//					distribution->getDistribution(ff, x1, x2, x3);
	//					real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
	//						(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
	//						(mfcbb - mfabb));
	//					real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
	//						(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
	//						(mfbcb - mfbab)) ;
	//					real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
	//						(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
	//						(mfbbc - mfbba)) ;
	//					real drho = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca) + (mfcaa + mfacc))) +
	//						(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
	//							((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
	//								(mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;
	//					//real mp= c3o1*(((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca) + (mfcaa + mfacc))) +
	//					//	c2o1*(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
	//					//		((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
	//					//			(mfbab + mfbcb) + (mfbba + mfbbc));
	//					//mp -= vvx * vvx - vvy * vvy - vvz * vvz;
	//					real drhoScaled = drho / scaleDistribution;
	//					if (((*phaseField)(x1, x2, x3) <= 0.5)) { drhoScaled = (*rhoNode)(x1, x2, x3); }

	//					//mp = 2 * drho - mp;
	//					for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {
	//						real feqOLD = (D3Q27System::getIncompFeqForDirection(fdir, drho, vvx, vvy, vvz));
	//						real feqNew = (D3Q27System::getIncompFeqForDirection(fdir, drhoScaled, vvx, vvy, vvz));
	//						distribution->setDistributionInvForDirection((ff[fdir]-feqOLD)* scaleStress +feqNew,x1,x2,x3,fdir);
	//					}





	//				}
	//			}
	//		}
	//	}
	//}

	real collFactorM;






	for (int x3 = minX3; x3 < maxX3; x3++) {
		for (int x2 = minX2; x2 < maxX2; x2++) {
			for (int x1 = minX1; x1 < maxX1; x1++) {
				if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
					int x1p = x1 + 1;
					int x2p = x2 + 1;
					int x3p = x3 + 1;



//					real mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
//					real mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
//					real mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
//					real mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
//					real mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
//					real mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
//					real mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
//					real mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
//					real mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
//					real mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
//					real mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
//					real mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
//					real mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
//
//					real mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
//					real mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
//					real mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
//					real mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
//					real mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
//					real mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
//					real mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
//					real mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
//					real mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
//					real mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
//					real mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
//					real mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
//					real mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);
//
//					real mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);
//
//					real m0, m1, m2;
//
//					real rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
//						+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
//						+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
//
//					real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
//						(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
//						(mfcbb - mfabb));
//					real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
//						(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
//						(mfbcb - mfbab));
//					real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
//						(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
//						(mfbbc - mfbba));
//
//					//forcing 
//					///////////////////////////////////////////////////////////////////////////////////////////
//					if (withForcing)
//					{
//						muX1 = static_cast<double>(x1 - 1 + ix1 * maxX1);
//						muX2 = static_cast<double>(x2 - 1 + ix2 * maxX2);
//						muX3 = static_cast<double>(x3 - 1 + ix3 * maxX3);
//
//						forcingX1 = muForcingX1.Eval();
//						forcingX2 = muForcingX2.Eval();
//						forcingX3 = muForcingX3.Eval();
//
//						vvx += forcingX1 * deltaT * 0.5; // X
//						vvy += forcingX2 * deltaT * 0.5; // Y
//						vvz += forcingX3 * deltaT * 0.5; // Z
//					}
//					///////////////////////////////////////////////////////////////////////////////////////////               
//					real oMdrho;
//
//					oMdrho = mfccc + mfaaa;
//					m0 = mfaca + mfcac;
//					m1 = mfacc + mfcaa;
//					m2 = mfaac + mfcca;
//					oMdrho += m0;
//					m1 += m2;
//					oMdrho += m1;
//					m0 = mfbac + mfbca;
//					m1 = mfbaa + mfbcc;
//					m0 += m1;
//					m1 = mfabc + mfcba;
//					m2 = mfaba + mfcbc;
//					m1 += m2;
//					m0 += m1;
//					m1 = mfacb + mfcab;
//					m2 = mfaab + mfccb;
//					m1 += m2;
//					m0 += m1;
//					oMdrho += m0;
//					m0 = mfabb + mfcbb;
//					m1 = mfbab + mfbcb;
//					m2 = mfbba + mfbbc;
//					m0 += m1 + m2;
//					m0 += mfbbb; //hat gefehlt
//					oMdrho = 1. - (oMdrho + m0);
//
//					real vx2;
//					real vy2;
//					real vz2;
//					vx2 = vvx * vvx;
//					vy2 = vvy * vvy;
//					vz2 = vvz * vvz;
//					////////////////////////////////////////////////////////////////////////////////////
//					real wadjust;
//					real qudricLimit = 0.01;
//					////////////////////////////////////////////////////////////////////////////////////
//					//Hin
//					////////////////////////////////////////////////////////////////////////////////////
//					// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// Z - Dir
//					m2 = mfaaa + mfaac;
//					m1 = mfaac - mfaaa;
//					m0 = m2 + mfaab;
//					mfaaa = m0;
//					m0 += c1o36 * oMdrho;
//					mfaab = m1 - m0 * vvz;
//					mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaba + mfabc;
//					m1 = mfabc - mfaba;
//					m0 = m2 + mfabb;
//					mfaba = m0;
//					m0 += c1o9 * oMdrho;
//					mfabb = m1 - m0 * vvz;
//					mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaca + mfacc;
//					m1 = mfacc - mfaca;
//					m0 = m2 + mfacb;
//					mfaca = m0;
//					m0 += c1o36 * oMdrho;
//					mfacb = m1 - m0 * vvz;
//					mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbaa + mfbac;
//					m1 = mfbac - mfbaa;
//					m0 = m2 + mfbab;
//					mfbaa = m0;
//					m0 += c1o9 * oMdrho;
//					mfbab = m1 - m0 * vvz;
//					mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbba + mfbbc;
//					m1 = mfbbc - mfbba;
//					m0 = m2 + mfbbb;
//					mfbba = m0;
//					m0 += c4o9 * oMdrho;
//					mfbbb = m1 - m0 * vvz;
//					mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbca + mfbcc;
//					m1 = mfbcc - mfbca;
//					m0 = m2 + mfbcb;
//					mfbca = m0;
//					m0 += c1o9 * oMdrho;
//					mfbcb = m1 - m0 * vvz;
//					mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcaa + mfcac;
//					m1 = mfcac - mfcaa;
//					m0 = m2 + mfcab;
//					mfcaa = m0;
//					m0 += c1o36 * oMdrho;
//					mfcab = m1 - m0 * vvz;
//					mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcba + mfcbc;
//					m1 = mfcbc - mfcba;
//					m0 = m2 + mfcbb;
//					mfcba = m0;
//					m0 += c1o9 * oMdrho;
//					mfcbb = m1 - m0 * vvz;
//					mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcca + mfccc;
//					m1 = mfccc - mfcca;
//					m0 = m2 + mfccb;
//					mfcca = m0;
//					m0 += c1o36 * oMdrho;
//					mfccb = m1 - m0 * vvz;
//					mfccc = m2 - 2. * m1 * vvz + vz2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// Y - Dir
//					m2 = mfaaa + mfaca;
//					m1 = mfaca - mfaaa;
//					m0 = m2 + mfaba;
//					mfaaa = m0;
//					m0 += c1o6 * oMdrho;
//					mfaba = m1 - m0 * vvy;
//					mfaca = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaab + mfacb;
//					m1 = mfacb - mfaab;
//					m0 = m2 + mfabb;
//					mfaab = m0;
//					mfabb = m1 - m0 * vvy;
//					mfacb = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaac + mfacc;
//					m1 = mfacc - mfaac;
//					m0 = m2 + mfabc;
//					mfaac = m0;
//					m0 += c1o18 * oMdrho;
//					mfabc = m1 - m0 * vvy;
//					mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbaa + mfbca;
//					m1 = mfbca - mfbaa;
//					m0 = m2 + mfbba;
//					mfbaa = m0;
//					m0 += c2o3 * oMdrho;
//					mfbba = m1 - m0 * vvy;
//					mfbca = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbab + mfbcb;
//					m1 = mfbcb - mfbab;
//					m0 = m2 + mfbbb;
//					mfbab = m0;
//					mfbbb = m1 - m0 * vvy;
//					mfbcb = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfbac + mfbcc;
//					m1 = mfbcc - mfbac;
//					m0 = m2 + mfbbc;
//					mfbac = m0;
//					m0 += c2o9 * oMdrho;
//					mfbbc = m1 - m0 * vvy;
//					mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcaa + mfcca;
//					m1 = mfcca - mfcaa;
//					m0 = m2 + mfcba;
//					mfcaa = m0;
//					m0 += c1o6 * oMdrho;
//					mfcba = m1 - m0 * vvy;
//					mfcca = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcab + mfccb;
//					m1 = mfccb - mfcab;
//					m0 = m2 + mfcbb;
//					mfcab = m0;
//					mfcbb = m1 - m0 * vvy;
//					mfccb = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfcac + mfccc;
//					m1 = mfccc - mfcac;
//					m0 = m2 + mfcbc;
//					mfcac = m0;
//					m0 += c1o18 * oMdrho;
//					mfcbc = m1 - m0 * vvy;
//					mfccc = m2 - 2. * m1 * vvy + vy2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// X - Dir
//					m2 = mfaaa + mfcaa;
//					m1 = mfcaa - mfaaa;
//					m0 = m2 + mfbaa;
//					mfaaa = m0;
//					m0 += 1. * oMdrho;
//					mfbaa = m1 - m0 * vvx;
//					mfcaa = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaba + mfcba;
//					m1 = mfcba - mfaba;
//					m0 = m2 + mfbba;
//					mfaba = m0;
//					mfbba = m1 - m0 * vvx;
//					mfcba = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaca + mfcca;
//					m1 = mfcca - mfaca;
//					m0 = m2 + mfbca;
//					mfaca = m0;
//					m0 += c1o3 * oMdrho;
//					mfbca = m1 - m0 * vvx;
//					mfcca = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaab + mfcab;
//					m1 = mfcab - mfaab;
//					m0 = m2 + mfbab;
//					mfaab = m0;
//					mfbab = m1 - m0 * vvx;
//					mfcab = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfabb + mfcbb;
//					m1 = mfcbb - mfabb;
//					m0 = m2 + mfbbb;
//					mfabb = m0;
//					mfbbb = m1 - m0 * vvx;
//					mfcbb = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfacb + mfccb;
//					m1 = mfccb - mfacb;
//					m0 = m2 + mfbcb;
//					mfacb = m0;
//					mfbcb = m1 - m0 * vvx;
//					mfccb = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfaac + mfcac;
//					m1 = mfcac - mfaac;
//					m0 = m2 + mfbac;
//					mfaac = m0;
//					m0 += c1o3 * oMdrho;
//					mfbac = m1 - m0 * vvx;
//					mfcac = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfabc + mfcbc;
//					m1 = mfcbc - mfabc;
//					m0 = m2 + mfbbc;
//					mfabc = m0;
//					mfbbc = m1 - m0 * vvx;
//					mfcbc = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					m2 = mfacc + mfccc;
//					m1 = mfccc - mfacc;
//					m0 = m2 + mfbcc;
//					mfacc = m0;
//					m0 += c1o9 * oMdrho;
//					mfbcc = m1 - m0 * vvx;
//					mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
//					////////////////////////////////////////////////////////////////////////////////////
//					// Cumulants
//					////////////////////////////////////////////////////////////////////////////////////
//					real OxxPyyPzz = 1.; //omega2 or bulk viscosity
//					real OxyyPxzz = 1.;//-s9;//2+s9;//
//					//real OxyyMxzz  = 1.;//2+s9;//
//					real O4 = 1.;
//					real O5 = 1.;
//					real O6 = 1.;
//					real OxyyMxzz = 1.;
//					//real OxyyPxzz = 1.;
//
//					//Cum 4.
//					//real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
//					//real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
//					//real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015
//
//					real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
//					real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
//					real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
//
//					real CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - 1.) * oMdrho);
//					real CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - 1.) * oMdrho);
//					real CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - 1.) * oMdrho);
//
//					//Cum 5.
//					real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
//					real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
//					real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;
//
//					//Cum 6.
//					real CUMccc = mfccc + ((-4. * mfbbb * mfbbb
//						- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
//						- 4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
//						- 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
//						+ (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
//							+ 2. * (mfcaa * mfaca * mfaac)
//							+ 16. * mfbba * mfbab * mfabb)
//						- c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho
//						- c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) - c1o27 * oMdrho * oMdrho * (-2. * oMdrho)
//						+ (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
//							+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 * oMdrho) + c1o27 * oMdrho;
//
//					//2.
//					// linear combinations
//					real mxxPyyPzz = mfcaa + mfaca + mfaac;
//					real mxxMyy = mfcaa - mfaca;
//					real mxxMzz = mfcaa - mfaac;
//
//					real dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
//					real dyuy = dxux + collFactorM * c3o2 * mxxMyy;
//					real dzuz = dxux + collFactorM * c3o2 * mxxMzz;
//
//					//relax
//					mxxPyyPzz += OxxPyyPzz * (mfaaa - mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
//					mxxMyy += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactor) * (vx2 * dxux - vy2 * dyuy);
//					mxxMzz += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactor) * (vx2 * dxux - vz2 * dzuz);
//
//					mfabb += collFactorM * (-mfabb);
//					mfbab += collFactorM * (-mfbab);
//					mfbba += collFactorM * (-mfbba);
//
//					// linear combinations back
//					mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
//					mfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
//					mfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);
//
//					//3.
//					// linear combinations
//					real mxxyPyzz = mfcba + mfabc;
//					real mxxyMyzz = mfcba - mfabc;
//
//					real mxxzPyyz = mfcab + mfacb;
//					real mxxzMyyz = mfcab - mfacb;
//
//					real mxyyPxzz = mfbca + mfbac;
//					real mxyyMxzz = mfbca - mfbac;
//
//					//relax
//					wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mfbbb) / (fabs(mfbbb) + qudricLimit);
//					mfbbb += wadjust * (-mfbbb);
//					wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxyPyzz) / (fabs(mxxyPyzz) + qudricLimit);
//					mxxyPyzz += wadjust * (-mxxyPyzz);
//					wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxyMyzz) / (fabs(mxxyMyzz) + qudricLimit);
//					mxxyMyzz += wadjust * (-mxxyMyzz);
//					wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxxzPyyz) / (fabs(mxxzPyyz) + qudricLimit);
//					mxxzPyyz += wadjust * (-mxxzPyyz);
//					wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxxzMyyz) / (fabs(mxxzMyyz) + qudricLimit);
//					mxxzMyyz += wadjust * (-mxxzMyyz);
//					wadjust = OxyyPxzz + (1. - OxyyPxzz) * fabs(mxyyPxzz) / (fabs(mxyyPxzz) + qudricLimit);
//					mxyyPxzz += wadjust * (-mxyyPxzz);
//					wadjust = OxyyMxzz + (1. - OxyyMxzz) * fabs(mxyyMxzz) / (fabs(mxyyMxzz) + qudricLimit);
//					mxyyMxzz += wadjust * (-mxyyMxzz);
//
//					// linear combinations back
//					mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
//					mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
//					mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
//					mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
//					mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
//					mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
//
//					//4.
//					CUMacc += O4 * (-CUMacc);
//					CUMcac += O4 * (-CUMcac);
//					CUMcca += O4 * (-CUMcca);
//
//					CUMbbc += O4 * (-CUMbbc);
//					CUMbcb += O4 * (-CUMbcb);
//					CUMcbb += O4 * (-CUMcbb);
//
//					//5.
//					CUMbcc += O5 * (-CUMbcc);
//					CUMcbc += O5 * (-CUMcbc);
//					CUMccb += O5 * (-CUMccb);
//
//					//6.
//					CUMccc += O6 * (-CUMccc);
//
//					//back cumulants to central moments
//					//4.
//					//mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
//					//mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
//					//mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015
//
//					mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
//					mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
//					mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
//
//					mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho;
//					mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho;
//					mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - 1) * oMdrho;
//
//					//5.
//					mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
//					mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
//					mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;
//
//					//6.
//					mfccc = CUMccc - ((-4. * mfbbb * mfbbb
//						- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
//						- 4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
//						- 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
//						+ (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
//							+ 2. * (mfcaa * mfaca * mfaac)
//							+ 16. * mfbba * mfbab * mfabb)
//						- c1o3 * (mfacc + mfcac + mfcca) * oMdrho - c1o9 * oMdrho * oMdrho
//						- c1o9 * (mfcaa + mfaca + mfaac) * oMdrho * (1. - 2. * oMdrho) - c1o27 * oMdrho * oMdrho * (-2. * oMdrho)
//						+ (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
//							+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 * oMdrho) - c1o27 * oMdrho;
//
//					////////////////////////////////////////////////////////////////////////////////////
//					//forcing
//					mfbaa = -mfbaa;
//					mfaba = -mfaba;
//					mfaab = -mfaab;
//					//////////////////////////////////////////////////////////////////////////////////////
//
//					////////////////////////////////////////////////////////////////////////////////////
//					//back
//					////////////////////////////////////////////////////////////////////////////////////
//					//mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// Z - Dir
//					m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + 1. * oMdrho) * (vz2 - vvz) * c1o2;
//					m1 = -mfaac - 2. * mfaab * vvz + mfaaa * (1. - vz2) - 1. * oMdrho * vz2;
//					m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + 1. * oMdrho) * (vz2 + vvz) * c1o2;
//					mfaaa = m0;
//					mfaab = m1;
//					mfaac = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
//					m1 = -mfabc - 2. * mfabb * vvz + mfaba * (1. - vz2);
//					m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
//					mfaba = m0;
//					mfabb = m1;
//					mfabc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
//					m1 = -mfacc - 2. * mfacb * vvz + mfaca * (1. - vz2) - c1o3 * oMdrho * vz2;
//					m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
//					mfaca = m0;
//					mfacb = m1;
//					mfacc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
//					m1 = -mfbac - 2. * mfbab * vvz + mfbaa * (1. - vz2);
//					m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
//					mfbaa = m0;
//					mfbab = m1;
//					mfbac = m2;
//					/////////b//////////////////////////////////////////////////////////////////////////
//					m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
//					m1 = -mfbbc - 2. * mfbbb * vvz + mfbba * (1. - vz2);
//					m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
//					mfbba = m0;
//					mfbbb = m1;
//					mfbbc = m2;
//					/////////b//////////////////////////////////////////////////////////////////////////
//					m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
//					m1 = -mfbcc - 2. * mfbcb * vvz + mfbca * (1. - vz2);
//					m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
//					mfbca = m0;
//					mfbcb = m1;
//					mfbcc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
//					m1 = -mfcac - 2. * mfcab * vvz + mfcaa * (1. - vz2) - c1o3 * oMdrho * vz2;
//					m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
//					mfcaa = m0;
//					mfcab = m1;
//					mfcac = m2;
//					/////////c//////////////////////////////////////////////////////////////////////////
//					m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
//					m1 = -mfcbc - 2. * mfcbb * vvz + mfcba * (1. - vz2);
//					m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
//					mfcba = m0;
//					mfcbb = m1;
//					mfcbc = m2;
//					/////////c//////////////////////////////////////////////////////////////////////////
//					m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
//					m1 = -mfccc - 2. * mfccb * vvz + mfcca * (1. - vz2) - c1o9 * oMdrho * vz2;
//					m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
//					mfcca = m0;
//					mfccb = m1;
//					mfccc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					//mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// Y - Dir
//					m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfaca - 2. * mfaba * vvy + mfaaa * (1. - vy2) - c1o6 * oMdrho * vy2;
//					m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfaaa = m0;
//					mfaba = m1;
//					mfaca = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfacb - 2. * mfabb * vvy + mfaab * (1. - vy2) - c2o3 * oMdrho * vy2;
//					m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfaab = m0;
//					mfabb = m1;
//					mfacb = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfacc - 2. * mfabc * vvy + mfaac * (1. - vy2) - c1o6 * oMdrho * vy2;
//					m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfaac = m0;
//					mfabc = m1;
//					mfacc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
//					m1 = -mfbca - 2. * mfbba * vvy + mfbaa * (1. - vy2);
//					m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
//					mfbaa = m0;
//					mfbba = m1;
//					mfbca = m2;
//					/////////b//////////////////////////////////////////////////////////////////////////
//					m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
//					m1 = -mfbcb - 2. * mfbbb * vvy + mfbab * (1. - vy2);
//					m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
//					mfbab = m0;
//					mfbbb = m1;
//					mfbcb = m2;
//					/////////b//////////////////////////////////////////////////////////////////////////
//					m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
//					m1 = -mfbcc - 2. * mfbbc * vvy + mfbac * (1. - vy2);
//					m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
//					mfbac = m0;
//					mfbbc = m1;
//					mfbcc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfcca - 2. * mfcba * vvy + mfcaa * (1. - vy2) - c1o18 * oMdrho * vy2;
//					m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfcaa = m0;
//					mfcba = m1;
//					mfcca = m2;
//					/////////c//////////////////////////////////////////////////////////////////////////
//					m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfccb - 2. * mfcbb * vvy + mfcab * (1. - vy2) - c2o9 * oMdrho * vy2;
//					m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfcab = m0;
//					mfcbb = m1;
//					mfccb = m2;
//					/////////c//////////////////////////////////////////////////////////////////////////
//					m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
//					m1 = -mfccc - 2. * mfcbc * vvy + mfcac * (1. - vy2) - c1o18 * oMdrho * vy2;
//					m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
//					mfcac = m0;
//					mfcbc = m1;
//					mfccc = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					//mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
//					////////////////////////////////////////////////////////////////////////////////////
//					// X - Dir
//					m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcaa - 2. * mfbaa * vvx + mfaaa * (1. - vx2) - c1o36 * oMdrho * vx2;
//					m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfaaa = m0;
//					mfbaa = m1;
//					mfcaa = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcba - 2. * mfbba * vvx + mfaba * (1. - vx2) - c1o9 * oMdrho * vx2;
//					m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfaba = m0;
//					mfbba = m1;
//					mfcba = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcca - 2. * mfbca * vvx + mfaca * (1. - vx2) - c1o36 * oMdrho * vx2;
//					m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfaca = m0;
//					mfbca = m1;
//					mfcca = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcab - 2. * mfbab * vvx + mfaab * (1. - vx2) - c1o9 * oMdrho * vx2;
//					m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfaab = m0;
//					mfbab = m1;
//					mfcab = m2;
//					///////////b////////////////////////////////////////////////////////////////////////
//					m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcbb - 2. * mfbbb * vvx + mfabb * (1. - vx2) - c4o9 * oMdrho * vx2;
//					m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfabb = m0;
//					mfbbb = m1;
//					mfcbb = m2;
//					///////////b////////////////////////////////////////////////////////////////////////
//					m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfccb - 2. * mfbcb * vvx + mfacb * (1. - vx2) - c1o9 * oMdrho * vx2;
//					m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfacb = m0;
//					mfbcb = m1;
//					mfccb = m2;
//					////////////////////////////////////////////////////////////////////////////////////
//					////////////////////////////////////////////////////////////////////////////////////
//					m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcac - 2. * mfbac * vvx + mfaac * (1. - vx2) - c1o36 * oMdrho * vx2;
//					m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfaac = m0;
//					mfbac = m1;
//					mfcac = m2;
//					///////////c////////////////////////////////////////////////////////////////////////
//					m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfcbc - 2. * mfbbc * vvx + mfabc * (1. - vx2) - c1o9 * oMdrho * vx2;
//					m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfabc = m0;
//					mfbbc = m1;
//					mfcbc = m2;
//					///////////c////////////////////////////////////////////////////////////////////////
//					m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
//					m1 = -mfccc - 2. * mfbcc * vvx + mfacc * (1. - vx2) - c1o36 * oMdrho * vx2;
//					m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
//					mfacc = m0;
//					mfbcc = m1;
//					mfccc = m2;
//
//					//////////////////////////////////////////////////////////////////////////
//					//proof correctness
//					//////////////////////////////////////////////////////////////////////////
//#ifdef  PROOF_CORRECTNESS
//					real rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
//						+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
//						+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
//					//real dif = fabs(rho - rho_post);
//					real dif = rho - rho_post;
//#ifdef SINGLEPRECISION
//					if (dif > 10.0E-7 || dif < -10.0E-7)
//#else
//					if (dif > 10.0E-15 || dif < -10.0E-15)
//#endif
//					{
//						UB_THROW(UbException(UB_EXARGS, "rho=" + UbSystem::toString(rho) + ", rho_post=" + UbSystem::toString(rho_post)
//							+ " dif=" + UbSystem::toString(dif)
//							+ " rho is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)
//							+ " in " + block.lock()->toString() + " step = " + UbSystem::toString(step)));
//					}
//#endif
//					//////////////////////////////////////////////////////////////////////////
//					//write distribution
//					//////////////////////////////////////////////////////////////////////////
//					(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = mfabb;
//					(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = mfbab;
//					(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = mfbba;
//					(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = mfaab;
//					(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = mfcab;
//					(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = mfaba;
//					(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = mfcba;
//					(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = mfbaa;
//					(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = mfbca;
//					(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = mfaaa;
//					(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = mfcaa;
//					(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = mfaca;
//					(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;
//
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac;
//					(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac;
//
//					(*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;
//
//					findNeighbors(phaseField, x1, x2, x3);
//					real dX1_phi = gradX1_phi();
//					real dX2_phi = gradX2_phi();
//					real dX3_phi = gradX3_phi();
//
//					real denom = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi)+ 1.0e-20;//+ 1e-9+1e-3;
//
//					real normX1 = dX1_phi / denom;
//					real normX2 = dX2_phi / denom;
//					real normX3 = dX3_phi / denom;
					//////////////////////////////////////////////////////////////////////////

//					////////////////old kernel
//					//////////////////////////////////////////////////////////////////////////
//					// Read distributions and phase field
//					////////////////////////////////////////////////////////////////////////////
//					//////////////////////////////////////////////////////////////////////////
//
//					// E   N  T
//					// c   c  c
//					//////////
//					// W   S  B
//					// a   a  a
//
//					// Rest ist b
//
//					// mfxyz
//					// a - negative
//					// b - null
//					// c - positive
//
//					// a b c
//					//-1 0 1

					findNeighbors(phaseField, x1, x2, x3);

					real mfcbb ;//= (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
					real mfbcb ;//= (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
					real mfbbc ;//= (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
					real mfccb ;//= (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
					real mfacb ;//= (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
					real mfcbc ;//= (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
					real mfabc ;//= (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
					real mfbcc ;//= (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
					real mfbac ;//= (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
					real mfccc ;//= (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
					real mfacc ;//= (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
					real mfcac ;//= (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
					real mfaac ;//= (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
					real mfabb ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
					real mfbab ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
					real mfbba ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
					real mfaab ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
					real mfcab ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
					real mfaba ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
					real mfcba ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
					real mfbaa ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
					real mfbca ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
					real mfaaa ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					real mfcaa ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
					real mfaca ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
					real mfcca ;//= (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);
					real mfbbb ;//= (*this->zeroDistributionsF)(x1, x2, x3);


					//real mfhcbb ;//= (*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3);
					//real mfhbcb ;//= (*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3);
					//real mfhbbc ;//= (*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3);
					//real mfhccb ;//= (*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3);
					//real mfhacb ;//= (*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3);
					//real mfhcbc ;//= (*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3);
					//real mfhabc ;//= (*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3);
					//real mfhbcc ;//= (*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3);
					//real mfhbac ;//= (*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3);
					//real mfhccc ;//= (*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3);
					//real mfhacc ;//= (*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3);
					//real mfhcac ;//= (*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3);
					//real mfhaac ;//= (*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//real mfhabb ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3);
					//real mfhbab ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3);
					//real mfhbba ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p);
					//real mfhaab ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3);
					//real mfhcab ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3);
					//real mfhaba ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p);
					//real mfhcba ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p);
					//real mfhbaa ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p);
					//real mfhbca ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p);
					//real mfhaaa ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//real mfhcaa ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//real mfhaca ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//real mfhcca ;//= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p);
					//real mfhbbb ;//= (*this->zeroDistributionsH2)(x1, x2, x3);


					//if (phi[DIR_000] < c1o2)
					//{
						 mfcbb= (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
						 mfbcb= (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
						 mfbbc= (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
						 mfccb= (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
						 mfacb= (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
						 mfcbc= (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
						 mfabc= (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
						 mfbcc= (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
						 mfbac= (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
						 mfccc= (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
						 mfacc= (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
						 mfcac= (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
						 mfaac= (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
						 mfabb= (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
						 mfbab= (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
						 mfbba= (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
						 mfaab= (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
						 mfcab= (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
						 mfaba= (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
						 mfcba= (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
						 mfbaa= (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
						 mfbca= (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
						 mfaaa= (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
						 mfcaa= (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
						 mfaca= (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
						 mfcca= (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);
						 mfbbb= (*this->zeroDistributionsF)(x1, x2, x3);


					//	 mfhcbb= (*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3);
					//	 mfhbcb= (*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3);
					//	 mfhbbc= (*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3);
					//	 mfhccb= (*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3);
					//	 mfhacb= (*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3);
					//	 mfhcbc= (*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3);
					//	 mfhabc= (*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3);
					//	 mfhbcc= (*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3);
					//	 mfhbac= (*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3);
					//	 mfhccc= (*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3);
					//	 mfhacc= (*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3);
					//	 mfhcac= (*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3);
					//	 mfhaac= (*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//	 mfhabb= (*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3);
					//	 mfhbab= (*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3);
					//	 mfhbba= (*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p);
					//	 mfhaab= (*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3);
					//	 mfhcab= (*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3);
					//	 mfhaba= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p);
					//	 mfhcba= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p);
					//	 mfhbaa= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p);
					//	 mfhbca= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p);
					//	 mfhaaa= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//	 mfhcaa= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//	 mfhaca= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//	 mfhcca= (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p);
					//	 mfhbbb= (*this->zeroDistributionsH2)(x1, x2, x3);

					//}
					//else
					//{
					//	mfhcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
					//	mfhbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
					//	mfhbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
					//	mfhccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
					//	mfhacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
					//	mfhcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
					//	mfhabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
					//	mfhbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
					//	mfhbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
					//	mfhccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
					//	mfhacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
					//	mfhcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
					//	mfhaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//	mfhabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
					//	mfhbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
					//	mfhbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
					//	mfhaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
					//	mfhcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
					//	mfhaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
					//	mfhcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
					//	mfhbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
					//	mfhbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
					//	mfhaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//	mfhcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//	mfhaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//	mfhcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);
					//	mfhbbb = (*this->zeroDistributionsF)(x1, x2, x3);


					//	mfcbb = (*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3);
					//	mfbcb = (*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3);
					//	mfbbc = (*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3);
					//	mfccb = (*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3);
					//	mfacb = (*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3);
					//	mfcbc = (*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3);
					//	mfabc = (*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3);
					//	mfbcc = (*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3);
					//	mfbac = (*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3);
					//	mfccc = (*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3);
					//	mfacc = (*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3);
					//	mfcac = (*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3);
					//	mfaac = (*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3);
					//	mfabb = (*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3);
					//	mfbab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3);
					//	mfbba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p);
					//	mfaab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3);
					//	mfcab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3);
					//	mfaba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p);
					//	mfcba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p);
					//	mfbaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p);
					//	mfbca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p);
					//	mfaaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					//	mfcaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p);
					//	mfaca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p);
					//	mfcca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p);
					//	mfbbb = (*this->zeroDistributionsH2)(x1, x2, x3);

					//}

					//real rhoH = 1.0;
					//real rhoL = 1.0 / densityRatio;

					real rhoH = 1.0;
					real rhoL = 1.0/ densityRatio;

					real rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

					real dX1_phi = gradX1_phi();
					real dX2_phi = gradX2_phi();
					real dX3_phi = gradX3_phi();

					real denom = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi)+ 1.0e-20;//+ 1e-9+1e-3;
					// 01.09.2022: unclear what value we have to add to the normal: lager values better cut of in gas phase?
					real normX1 = dX1_phi / denom;
					real normX2 = dX2_phi / denom;
					real normX3 = dX3_phi / denom;


					//real pushInterface = 2.0;
					//collFactorM = collFactorL + (collFactorL - collFactorG) * (phi[DIR_000] - phiH) / (phiH - phiL);
					//collFactorM = collFactorL + (collFactorL - collFactorG) * (tanh(pushInterface * (c2o1 * phi[DIR_000] - c1o1)) / tanh(pushInterface) * c1o2 + c1o2 - phiH) / (phiH - phiL);
					collFactorM = phi[DIR_000] > phiLim ? collFactorL : collFactorG;
					//collFactorM=(((*phaseField)(x1, x2, x3) > c1o2) && ((*phaseFieldOld)(x1, x2, x3) <= c1o2)) ? 1.8 : collFactorM;
					real collFactorMInv = phi[DIR_000] > phiLim ? collFactorG : collFactorL;

					real mu = 2 * beta * phi[DIR_000] * (phi[DIR_000] - 1) * (2 * phi[DIR_000] - 1) - kappa * nabla2_phi();

					//----------- Calculating Macroscopic Values -------------
					real rho = phi[DIR_000] > phiLim ? rhoH : rhoL;//rhoH + rhoToPhi * (phi[DIR_000] - phiH); //Incompressible

					//real rho = rhoH + rhoToPhi * (tanh(pushInterface*(c2o1*phi[DIR_000]-c1o1))/tanh(pushInterface)*c1o2 +c1o2 - phiH); //Incompressible
																		///scaled phase field
					//real rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) * (*phaseField)(x1, x2, x3) / ((*phaseField)(x1, x2, x3) * (*phaseField)(x1, x2, x3) + (c1o1 - (*phaseField)(x1, x2, x3)) * (c1o1 - (*phaseField)(x1, x2, x3))) - phiH);
					///!scaled phase field
					
					//real rho = rhoH + rhoToPhi * (phi[DIR_000] - phiH)+(c1o1-phi[DIR_000])* (*pressure)(x1, x2, x3)*c3o1; //compressible
					//real rho = rhoL + (rhoH - rhoL) * phi[DIR_000] + (c1o1 - phi[DIR_000]) * (*pressure)(x1, x2, x3) * c3o1; //compressible

					real m0, m1, m2;
					real rhoRef=c1o1;

					real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
						(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
						(mfcbb - mfabb))/rhoRef;
					real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
						(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
						(mfbcb - mfbab))/rhoRef;
					real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
						(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
						(mfbbc - mfbba))/rhoRef;
					////Filter&Gradient merged
					//real pressureHere = (*pressureOld)(x1, x2, x3);
					//real pressureHere = (*pressure)(x1, x2, x3);

					//real arrayP[3][3][3] = { {{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}},
					//							{{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}},
					//							{ {pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}} };
					//real LaplaceP = 0.0;
					//LaplaceP += WEIGTH[DIR_PPP] * (((((*pressureOld)(x1+1,x2+1,x3+1) - pressureHere) + ((*pressureOld)(x1 - 1, x2 - 1, x3 - 1) - pressureHere)) + (((*pressureOld)(x1 + 1, x2 + 1, x3 - 1) - pressureHere) + ((*pressureOld)(x1 - 1, x2 - 1, x3 + 1) - pressureHere)))
					//	+ ((((*pressureOld)(x1 + 1, x2 - 1, x3 + 1) - pressureHere) + ((*pressureOld)(x1 - 1, x2 + 1, x3 - 1) - pressureHere)) + (((*pressureOld)(x1 + 1, x2 - 1, x3 - 1) - pressureHere) + ((*pressureOld)(x1 - 1, x2 + 1, x3 + 1) - pressureHere))));
					//LaplaceP += WEIGTH[DIR_0PP] * (
					//	((((*pressureOld)(x1 + 1, x2 + 1, x3) - pressureHere) + ((*pressureOld)(x1 - 1, x2 - 1, x3) - pressureHere)) + (((*pressureOld)(x1 + 1, x2 - 1, x3) - pressureHere) + ((*pressureOld)(x1 - 1, x2 + 1, x3) - pressureHere)))
					//	+ ((((*pressureOld)(x1 + 1, x2, x3 + 1) - pressureHere) + ((*pressureOld)(x1 - 1, x2, x3 -1) - pressureHere)) + (((*pressureOld)(x1 + 1, x2, x3 - 1) - pressureHere) + ((*pressureOld)(x1 - 1, x2, x3 + 1) - pressureHere)))
					//	+ ((((*pressureOld)(x1, x2 + 1, x3 + 1) - pressureHere) + ((*pressureOld)(x1, x2 - 1, x3 - 1) - pressureHere)) + (((*pressureOld)(x1, x2 + 1, x3 - 1) - pressureHere) + ((*pressureOld)(x1, x2 - 1, x3 + 1) - pressureHere)))
					//	);
					//LaplaceP += WEIGTH[DIR_00P] * (
					//	(((*pressureOld)(x1 + 1, x2, x3) - pressureHere) + ((*pressureOld)(x1, x2-1, x3) - pressureHere))
					//	+ (((*pressureOld)(x1, x2 + 1, x3) - pressureHere) + ((*pressureOld)(x1, x2 - 1, x3) - pressureHere))
					//	+ (((*pressureOld)(x1, x2, x3 + 1) - pressureHere) + ((*pressureOld)(x1, x2, x3 - 1) - pressureHere))
					//	);

					//LaplaceP= 6.0 * LaplaceP;
					
					//real sum = 0.0;

//					for (int dir1 = -1; dir1 <= 1; dir1++) {
//						for (int dir2 = -1; dir2 <= 1; dir2++) {
//							for (int dir3 = -1; dir3 <= 1; dir3++){
//								int xxx = x1 + dir1;
//								int yyy = x2 + dir2;
//								int zzz = x3 + dir3;
//								if (!bcArray->isSolid(xxx, yyy, zzz) && !bcArray->isUndefined(xxx, yyy, zzz)) arrayP[dir1 + 1][dir2 + 1][dir3 + 1] = (*pressureOld)(xxx, yyy, zzz);
//								//if (!bcArray->isSolid(xxx, yyy, zzz) && !bcArray->isUndefined(xxx, yyy, zzz)) arrayP[dir1 + 1][dir2 + 1][dir3 + 1] = (*pressure)(xxx, yyy, zzz);
//							//	sum += 64.0 / (216.0 * (c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)) * (c1o1 + c3o1 * abs(dir3))) * arrayP[dir1 + 1][dir2 + 1][dir3 + 1];
//							}
//						}
//					}
////					(*pressure)(x1, x2, x3) = sum;// *0.1 + (1.0 - 0.1) * (*pressureOld)(x1, x2, x3);
//
//
//					(*pressure)(x1, x2, x3) = (((((arrayP[0][0][0] + arrayP[2][2][2]) + (arrayP[0][2][0] + arrayP[2][0][2])) + ((arrayP[2][0][0] + arrayP[0][2][2]) + (arrayP[2][2][0] + arrayP[0][0][2]))) * c1o216
//						+ (((arrayP[0][0][1] + arrayP[2][2][1]) + (arrayP[0][1][0] + arrayP[2][1][2])) + ((arrayP[1][0][0] + arrayP[1][2][2]) + (arrayP[0][1][2] + arrayP[2][1][0])) + ((arrayP[1][0][2] + arrayP[1][2][0]) + (arrayP[0][2][1] + arrayP[2][0][1]))) * c1o54)
//						+ ((arrayP[0][1][1] + arrayP[2][1][1]) + (arrayP[1][0][1] + arrayP[1][2][1]) + (arrayP[1][1][0] + arrayP[1][1][2])) * c2o27)
//						+ arrayP[1][1][1] * c8o27;
					//real gradPx = 0.0;
					//real gradPy = 0.0;
					//real gradPz = 0.0;
					//for (int dir1 = -1; dir1 <= 1; dir1++) {
					//	for (int dir2 = -1; dir2 <= 1; dir2++) {
					//		gradPx -= arrayP[0][dir1+1][dir2+1] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		gradPx += arrayP[2][dir1+1][dir2+1] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));



					//		gradPy -= arrayP[dir1+1][0][dir2+1] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		gradPy += arrayP[dir1+1][2][dir2+1] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		

					//		gradPz -= arrayP[dir1+1][dir2+1][0] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		gradPz += arrayP[dir1+1][dir2+1][2] * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//	}
					//}

					//real gradPx = ((((arrayP[2][0][0] - arrayP[0][2][2]) + (arrayP[2][2][0] - arrayP[0][0][2])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[2][0][2] - arrayP[0][2][0]))) * c1o72
					//	+ (((arrayP[2][1][0] - arrayP[0][1][2]) + (arrayP[2][2][1] - arrayP[0][0][1])) + ((arrayP[2][0][1] - arrayP[0][2][1]) + (arrayP[2][1][2] - arrayP[0][1][0]))) * c1o18)
					//	+ (arrayP[2][1][1] - arrayP[0][1][1]) * c2o9;
					//real gradPy = ((((arrayP[0][2][0] - arrayP[2][0][2]) + (arrayP[2][2][0] - arrayP[0][0][2])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[0][2][2] - arrayP[2][0][0]))) * c1o72
					//	+ (((arrayP[1][2][0] - arrayP[1][0][2]) + (arrayP[2][2][1] - arrayP[0][0][1])) + ((arrayP[0][2][1] - arrayP[2][0][1]) + (arrayP[1][2][2] - arrayP[1][0][0]))) * c1o18)
					//	+ (arrayP[1][2][1] - arrayP[1][0][1]) * c2o9;
					//real gradPz = ((((arrayP[0][0][2] - arrayP[2][2][0]) + (arrayP[0][2][2] - arrayP[2][0][0])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[2][0][2] - arrayP[0][2][0]))) * c1o72
					//	+ (((arrayP[0][1][2] - arrayP[2][1][0]) + (arrayP[1][2][2] - arrayP[1][0][0])) + ((arrayP[1][0][2] - arrayP[1][2][0]) + (arrayP[2][1][2] - arrayP[0][1][0]))) * c1o18)
					//	+ (arrayP[1][1][2] - arrayP[1][1][0]) * c2o9;

					//gradPx *=c1o1 - (*pressure)(x1, x2, x3)+pressureHere;
					//gradPy *=c1o1 - (*pressure)(x1, x2, x3) + pressureHere;
					//gradPz *=c1o1 - (*pressure)(x1, x2, x3) + pressureHere;

					////!Filter&Gradient merged
					//real gradPx = 0.0;
					//real gradPy = 0.0;
					//real gradPz = 0.0;
					//for (int dir1 = -1; dir1 <= 1; dir1++) {
					//	for (int dir2 = -1; dir2 <= 1; dir2++) {
					//		int yyy = x2 + dir1;
					//		int zzz = x3 + dir2;
					//		if (!bcArray->isSolid(x1-1, yyy, zzz) && !bcArray->isUndefined(x1-1, yyy, zzz)) {
					//			gradPx -= (*pressure)(x1 - 1, yyy, zzz) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPx -= (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(x1 + 1, yyy, zzz) && !bcArray->isUndefined(x1 + 1, yyy, zzz)) {
					//			gradPx += (*pressure)(x1 + 1, yyy, zzz) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPx += (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}

					//		int xxx = x1 + dir1;
					//		if (!bcArray->isSolid(xxx, x2-1, zzz) && !bcArray->isUndefined(xxx, x2-1, zzz)) {
					//			gradPy -= (*pressure)(xxx, x2-1, zzz) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPy -= (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(xxx, x2+1, zzz) && !bcArray->isUndefined(xxx, x2+1, zzz)) {
					//			gradPy += (*pressure)(xxx, x2+1, zzz) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPy += (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}

					//		yyy = x2 + dir2;
					//		if (!bcArray->isSolid(xxx, yyy, x3-1) && !bcArray->isUndefined(xxx, yyy, x3-1)) {
					//			gradPz -= (*pressure)(xxx, yyy, x3-1) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPz -= (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(xxx, yyy, x3+1) && !bcArray->isUndefined(xxx, yyy, x3+1)) {
					//			gradPz += (*pressure)(xxx, yyy, x3+1) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}
					//		else {
					//			gradPz += (*pressure)(x1, x2, x3) * c2o9 / ((c1o1 + c3o1 * abs(dir1)) * (c1o1 + c3o1 * abs(dir2)));
					//		}

					//	}
					//}

					//Viscosity increase by phase field residuum
					//real errPhi = (((1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale)- denom);
					//real limVis = 0.01;// 0.0000001 * 10;//0.01;
					// collFactorM =collFactorM/(c1o1+limVis*(errPhi*errPhi)*collFactorM);
					// collFactorM = (collFactorM < 1.8) ? 1.8 : collFactorM;
					//errPhi = errPhi * errPhi* errPhi * errPhi * errPhi * errPhi;
					//collFactorM = collFactorM + (1.8 - collFactorM) * errPhi / (errPhi + limVis);

					//3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) - (phi2[DIR_PPM] - phi2[DIR_MMP])))
					//+WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) - (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_0MP] - phi2[DIR_0PM]) + (phi2[DIR_0PP] - phi2[DIR_0MM])))) +
					//+WEIGTH[DIR_0P0] * (phi2[DIR_00P] - phi2[DIR_00M]));



					////external pressure
					//forcingX1 =/* muForcingX1.Eval()/rho */- gradPx/rho;
					//forcingX2 =/* muForcingX2.Eval()/rho */- gradPy/rho;
					//forcingX3 =/* muForcingX3.Eval()/rho */- gradPz/rho;

					///////////////////////////////////////////////

					//real pBefore = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb) * c1o3;
					//pBefore = -c1o3 * (-1.0e-10)/((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) );
					////if (vvx * vvx + vvy * vvy + vvz * vvz > 1.0e-100) {
					//	mfabb -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_P00] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbab -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0P0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbba -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_00P] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfaab -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PP0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcab -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MP0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfaba -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_P0P] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcba -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_M0P] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbaa -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0PP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbca -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0MP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfaaa -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PPP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcaa -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MPP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfaca -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PMP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcca -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MMP] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcbb -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_M00] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbcb -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0M0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbbc -= pBefore * c2o9 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_00M] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfccb -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MM0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfacb -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PM0] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcbc -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_M0M] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfabc -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_P0M] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbcc -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0MM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbac -= pBefore * c1o18 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_0PM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfccc -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MMM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfacc -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PMM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfcac -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_MPM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfaac -= pBefore * c1o72 * ((rhoL + phi[DIR_000] * (rhoH - rhoL) / (phiH - phiL)) / (rhoL + phi[DIR_PPM] * (rhoH - rhoL) / (phiH - phiL)));
					//	mfbbb -= pBefore * 8.0 / 9.0;
					//}

					///////////////////////////////////////////////

					//real pStarStart = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb) * c1o3;

					//rho = rhoH + rhoToPhi * ((phi[DIR_000] - phiH)+fabs(pStarStart)*0); //Incompressible

					muRho = rho;


					/////////////////////

					 vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
						(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
						(mfcbb - mfabb)) / rhoRef;
					 vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
						(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
						(mfbcb - mfbab)) / rhoRef;
					 vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
						(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
						(mfbbc - mfbba)) / rhoRef;


					 //real dRhoInvX = -(((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhcaa - mfhacc) + (mfhcca - mfhaac))) +
						// (((mfhcba - mfhabc) + (mfhcbc - mfhaba)) + ((mfhcab - mfhacb) + (mfhccb - mfhaab))) +
						// (mfhcbb - mfhabb)));
					 //real dRhoInvY = -(((((mfhccc - mfhaaa) + (mfhaca - mfhcac)) + ((mfhacc - mfhcaa) + (mfhcca - mfhaac))) +
						// (((mfhbca - mfhbac) + (mfhbcc - mfhbaa)) + ((mfhacb - mfhcab) + (mfhccb - mfhaab))) +
						// (mfhbcb - mfhbab)));
					 //real dRhoInvZ = -(((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhacc - mfhcaa) + (mfhaac - mfhcca))) +
						// (((mfhbac - mfhbca) + (mfhbcc - mfhbaa)) + ((mfhabc - mfhcba) + (mfhcbc - mfhaba))) +
						// (mfhbbc - mfhbba)));


					 forcingX1 = 0.0;
					 forcingX2 = 0.0;
					 forcingX3 = 0.0;
					//!Abbas
					//real dX1_rhoInv = gradX1_rhoInv(rhoL, rhoH - rhoL);
					//real dX2_rhoInv = gradX2_rhoInv(rhoL, rhoH - rhoL);
					//real dX3_rhoInv = gradX3_rhoInv(rhoL, rhoH - rhoL);
					//forcingX1 =/* muForcingX1.Eval() / rho*/ +pStar * dX1_rhoInv * rho;
					//forcingX2 =/* muForcingX2.Eval() / rho*/ +pStar * dX2_rhoInv * rho;
					//forcingX3 =/* muForcingX3.Eval() / rho*/ +pStar * dX3_rhoInv * rho;

					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho + pStar * dX1_rhoInv * rho) *c1o2;
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho + pStar * dX2_rhoInv * rho) *c1o2;
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho + pStar * dX3_rhoInv * rho) *c1o2;
					 //real FdX1_phi = normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					 //real FdX2_phi = normX2 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					 //real FdX3_phi = normX3 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;


					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho ) ;
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho ) ;
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho ) ;

					//forcingX1 = (pStar * dRhoInvX* rho *c3o1) ;
					//forcingX2 = (pStar * dRhoInvY* rho *c3o1) ;
					//forcingX3 = (pStar * dRhoInvZ* rho *c3o1) ;
					//if (phi[DIR_000] > 0.1 && phi[DIR_000] < 0.9) std::cout << phi[DIR_000] << " " << dX1_phi * rhoToPhi / rho << " " << dRhoInvX * rho *3<< std::endl;
					//real forcingX1ALTERNAT = ( pStar * dX1_rhoInv * rho) ;
					//real forcingX2ALTERNAT = ( pStar * dX2_rhoInv * rho) ;
					//real forcingX3ALTERNAT = ( pStar * dX3_rhoInv * rho) ;

					//forcingX1 = (fabs(vvx + c1o2 * forcingX1) < fabs(vvx + c1o2 * forcingX1ALTERNAT)) ? forcingX1 : forcingX1ALTERNAT;
					//forcingX2 = (fabs(vvy + c1o2 * forcingX2) < fabs(vvy + c1o2 * forcingX2ALTERNAT)) ? forcingX2 : forcingX2ALTERNAT;
					//forcingX3 = (fabs(vvz + c1o2 * forcingX3) < fabs(vvz + c1o2 * forcingX3ALTERNAT)) ? forcingX3 : forcingX3ALTERNAT;

					//	 forcingX1 = -pStar * rhoToPhi / rho * normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					//	 forcingX2 = -pStar * rhoToPhi / rho * normX2 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					//	 forcingX3 = -pStar * rhoToPhi / rho * normX3 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;

					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho *(c1o1- phi[DIR_000]) + pStar * dX1_rhoInv * rho*(phi[DIR_000]));
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho *(c1o1- phi[DIR_000]) + pStar * dX2_rhoInv * rho*(phi[DIR_000]));
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho *(c1o1- phi[DIR_000]) + pStar * dX3_rhoInv * rho*(phi[DIR_000]));
						 //if (phi[DIR_000] > 0.3 && phi[DIR_000] < 0.7)
						 //{
							// int test = 1;
							// std::cout << phi[DIR_000] <<" "<< dX1_phi <<" "<< normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale<<" "<< normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale/ dX1_phi<< std::endl;
						 //}



					 //real scaleGrad = c2o1 * phi[DIR_000] * (1.0 - phi[DIR_000]) / ((phi[DIR_000] * phi[DIR_000] + (1.0 - phi[DIR_000]) * (1.0 - phi[DIR_000])) * (phi[DIR_000] * phi[DIR_000] + (1.0 - phi[DIR_000]) * (1.0 - phi[DIR_000])));
					 //dX1_phi *= scaleGrad;
					 //dX2_phi *= scaleGrad;
					 //dX3_phi *= scaleGrad;

					 ///Experimental interface sharpening force 20.06.2022

					 real scaleSharpener = 1.0;
					 //forcingX1 += scaleSharpener * (FdX1_phi - dX1_phi) * fabsf(FdX1_phi - dX1_phi)  / rho;
					 //forcingX2 += scaleSharpener * (FdX2_phi - dX2_phi) * fabsf(FdX2_phi - dX2_phi)  / rho;
					 //forcingX3 += scaleSharpener * (FdX3_phi - dX3_phi) * fabsf(FdX3_phi - dX3_phi)  / rho;
					///surface tension force
					//forcingX1 += mu * dX1_phi/rho;
					//forcingX2 += mu * dX2_phi/rho;
					//forcingX3 += mu * dX3_phi/rho;

					//real forcingBIAS = 0.5;
				//	forcingX1 += muForcingX1.Eval() / rho;//*phi[DIR_000];
                //     forcingX2 += -5.0e-7;//  *phi[DIR_000];                         // muForcingX2.Eval() / rho - 5.0e-7 * phi[DIR_000] * 0;// * phi[DIR_000];
				//	forcingX3 += muForcingX3.Eval() / rho;// * phi[DIR_000];

				//	//19.08.2022
					//vvx += vvxh / rho * c1o2;
					//vvy += vvyh / rho * c1o2;
					//vvz += vvzh / rho * c1o2;
				//	//


				//	vvx += (forcingX1) * deltaT * c1o2;
				//	vvy += (forcingX2) * deltaT * c1o2;
				//	vvz += (forcingX3) * deltaT * c1o2;

					    if (withForcing) {

                        forcingX1 += muForcingX1.Eval();
                        forcingX2 += muForcingX2.Eval();
                        forcingX3 += muForcingX3.Eval();

                        vvx += (forcingX1)*deltaT * c1o2;
                        vvy += (forcingX2)*deltaT * c1o2;
                        vvz += (forcingX3)*deltaT * c1o2;
                    }

					//vvx += (forcingX1 + muForcingX1.Eval() / rho) * deltaT *  c1o2; // X
					//vvy += (forcingX2 + muForcingX2.Eval() / rho) * deltaT *  c1o2; // Y
					//vvz += (forcingX3 + muForcingX3.Eval() / rho) * deltaT *  c1o2; // Z



				//	vvx += (forcingX1 + muForcingX1.Eval() / rho) * deltaT * forcingBIAS; // X
				//	vvy += (forcingX2 + muForcingX2.Eval() / rho) * deltaT * forcingBIAS; // Y
				//	vvz += (forcingX3 + muForcingX3.Eval() / rho) * deltaT * forcingBIAS; // Z



					real vx2;
					real vy2;
					real vz2;
					vx2 = vvx * vvx;
					vy2 = vvy * vvy;
					vz2 = vvz * vvz;
					//pStar =ppStar- (vx2 + vy2 + vz2)*pStar;
				//	pStar = (pStar + ppStar)*c1o2;
					///////////////////////////////////////////////////////////////////////////////////////////               
					real oMdrho;
					///////////////

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
					real wadjust;
					//real qudricLimit = 0.01 / (c1o1 + 1.0e5 * phi[DIR_000] * (c1o1 - phi[DIR_000])); //real qudricLimit = 0.01;
					//real qudricLimit = (((*phaseField)(x1, x2, x3) > c1o2) && (normX1 * vvx + normX2 * vvy + normX3 * vvz < 0)) ? 0.01 / (c1o1 + 1.0e5 * phi[DIR_000] * (c1o1 - phi[DIR_000])) : 0.01;
					//real qudricLimit = (((*phaseField)(x1, x2, x3) > c1o2)&& (normX1*vvx+normX2*vvy+normX3*vvz<0) ) ? 0.01 / (c1o1 + 1.0e5 * phi[DIR_000] * (c1o1 - phi[DIR_000])) : 0.01 ;
                    real qudricLimit = 0.0001;
                    /// (c1o1 + 1.0e3 * (dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi));
					//real qudricLimit = 0.01 / (c1o1 + 1.0e5 * phi[DIR_000] * (c1o1 - phi[DIR_000])) ;
					//real qudricLimit = (((*phaseField)(x1, x2, x3) > c1o2) ) ? 0.01 / (c1o1 + 1.0e5 * phi[DIR_000] * (c1o1 - phi[DIR_000])) : 0.01 ;
					//qudricLimit = (((*phaseField)(x1, x2, x3)-c1o2 ) * (normX1 * vvx + normX2 * vvy + normX3 * vvz) < 0) ? 0.01 / (c1o1 + 1.0e8 * phi[DIR_000] * (c1o1 - phi[DIR_000])) : 0.01;
				//	if (phi[DIR_000] > c1o2 && (*phaseFieldOld)(x1, x2, x3) <= c1o2) collFactorM = 1.8;
					
																													////////////////////////////////////////////////////////////////////////////////////
						//! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
						//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
						//! see also Eq. (6)-(14) in \ref
						//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
						//!
						////////////////////////////////////////////////////////////////////////////////////
						// Z - Dir
					forwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36, oMdrho);
					forwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36, oMdrho);
					forwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36, oMdrho);
					forwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36, oMdrho);

					////////////////////////////////////////////////////////////////////////////////////
					// Y - Dir
					forwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6, oMdrho);
					forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
					forwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18, oMdrho);
					forwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3, oMdrho);
					forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
					forwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9, oMdrho);
					forwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6, oMdrho);
					forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
					forwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18, oMdrho);

					////////////////////////////////////////////////////////////////////////////////////
					// X - Dir
					forwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1, oMdrho);
					forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
					forwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3, oMdrho);
					forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
					forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
					forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
					forwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3, oMdrho);
					forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
					forwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, vvx, vx2, c3o1, c1o9, oMdrho);

																							  
																							  ////////////////////////////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////////////////////////////
					// Cumulants
					////////////////////////////////////////////////////////////////////////////////////

					// mfaaa = 0.0;
					real OxxPyyPzz = 1.0; //omega2 or bulk viscosity
											//  real OxyyPxzz = 1.;//-s9;//2+s9;//
											//  real OxyyMxzz  = 1.;//2+s9;//
					real O4 = 1.;
					real O5 = 1.;
					real O6 = 1.;

					//collFactorM+= (1.7 - collFactorM) * fabs(mfaaa) / (fabs(mfaaa) + 0.001f);
					//////
					//M110 -= vvx * vvy - mfbba;
					//M101 -= vvx * vvz - mfbab;
					//M011 -= vvy * vvz - mfabb;



					//M200 -= vx2;
					//M020 -= vy2;
					//M002 -= vz2;
					//real Mp = (M200 + M020 + M002);



					//M200 -= c1o3 * Mp;
					//M020 -= c1o3 * Mp;
					//M002 -= c1o3 * Mp;

					//M200 -= -mfcaa + (mfcaa + mfaca + mfaac) * c1o3;
					//M020 -= -mfaca + (mfcaa + mfaca + mfaac) * c1o3;
					//M002 -= -mfaac + (mfcaa + mfaca + mfaac) * c1o3;

					/////

					/////fourth order parameters; here only for test. Move out of loop!

					real OxyyPxzz =  8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
					real OxyyMxzz = 8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
				    real Oxyz = 24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 * collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0));
					real A = (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//FIXME:  warning C4459: declaration of 'B' hides global declaration (message : see declaration of 'D3Q27System::B' )
					real BB =  (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//real stress = 1.0;// stress / (stress + 1.0e-10);
					//stress = 1.0;
					//OxyyPxzz += stress*(1.0-OxyyPxzz);
					//OxyyPxzz = c3o1 * (collFactorM - c2o1) / (collFactorM - c3o1);
					//OxyyMxzz += stress*(1.0-OxyyMxzz);
					//Oxyz +=  stress*(1.0-Oxyz);
					//A *= 1.0-stress;
					//BB *= 1.0-stress;

					//Cum 4.
					//real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
					//real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
					//real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

					real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
					real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
					real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

					real CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho);
					real CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho);
					real CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho);

					//Cum 5.
					real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
					real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
					real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

					//Cum 6.
					real CUMccc = mfccc + ((-4. * mfbbb * mfbbb
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
					real mxxPyyPzz = mfcaa + mfaca + mfaac;
					//pStar = (mxxPyyPzz+vx2+vy2+vz2) * c1o3;//does not work
					//pStar = (mxxPyyPzz) * c1o3;
					//pStar = pStar + 1.5 * (mxxPyyPzz * c1o3 - pStar);
					//mfaaa = mxxPyyPzz;
					//  real mfaaaS = (mfaaa * (-4 - 3 * OxxPyyPzz * (-1 + rho)) + 6 * mxxPyyPzz * OxxPyyPzz * (-1 + rho)) / (-4 + 3 * OxxPyyPzz * (-1 + rho));
					mxxPyyPzz -= mfaaa ;//12.03.21 shifted by mfaaa
										//mxxPyyPzz-=(mfaaa+mfaaaS)*c1o2;//12.03.21 shifted by mfaaa
					//dirty (04.04.2023)
					//if (phi[DIR_000] > c1o2  &&(dX1_phi*vvx+dX2_phi*vvy+dX3_phi*vvz)<0) {
					//	//collFactorM = c1o1 / (c1o1 / collFactorM +1e10* fabsf(mxxPyyPzz)*(c1o1-fabsf(phi[DIR_000])));
					//	collFactorM = c1o1 / (c1o1 / collFactorM - 1e15*(dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz)* fabsf(mxxPyyPzz) );
					//	collFactorM = (collFactorM < 1.8) ? 1.8 : collFactorM;
					//	 OxyyPxzz = 8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
					//	 OxyyMxzz = 8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
					//	 Oxyz = 24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 * collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0));
					//	 A = (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//	 BB = (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));

					//}
					real mxxMyy = mfcaa - mfaca;
					real mxxMzz = mfcaa - mfaac;

					///
					real mmfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
					real mmfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
					real mmfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);
					real mmfabb = mfabb;
					real mmfbab = mfbab;
					real mmfbba = mfbba;
					///

					real dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz);// *0;
					//real dxux = -c1o2 * (mxxMyy + mxxMzz) * collFactorM - mfaaa * c1o3* omegaDRho;
					real dyuy =  dxux + collFactorM * c3o2 * mxxMyy;
					real dzuz =  dxux + collFactorM * c3o2 * mxxMzz;
					real Dxy = -c3o1 * collFactorM * mfbba;
					real Dxz = -c3o1 * collFactorM * mfbab;
					real Dyz = -c3o1 * collFactorM * mfabb;
                    // if (phi[DIR_000] > phiLim) 
						if ((phi[DIR_000] > phiLim) && ((phi[DIR_P00] <= phiLim) || (phi[DIR_M00] <= phiLim) || (phi[DIR_00P] <= phiLim) || (phi[DIR_00M] <= phiLim) || (phi[DIR_0M0] <= phiLim) || (phi[DIR_0P0] <= phiLim) || (phi[DIR_PP0] <= phiLim) || (phi[DIR_PM0] <= phiLim) || (phi[DIR_P0P] <= phiLim) ||
                                                  (phi[DIR_P0M] <= phiLim) || (phi[DIR_MP0] <= phiLim) || (phi[DIR_MM0] <= phiLim) || (phi[DIR_M0P] <= phiLim) || (phi[DIR_M0M] <= phiLim) || (phi[DIR_0PM] <= phiLim) || (phi[DIR_0MM] <= phiLim) || (phi[DIR_0PP] <= phiLim) || (phi[DIR_0MP] <= phiLim) ||
                                                  (phi[DIR_PPP] <= phiLim) || (phi[DIR_PMP] <= phiLim) || (phi[DIR_MPP] <= phiLim) || (phi[DIR_MMP] <= phiLim) ||
                         (phi[DIR_PPM] <= phiLim) || (phi[DIR_PMM] <= phiLim) || (phi[DIR_MPM] <= phiLim) || (phi[DIR_MMM] <= phiLim))) {

					// {
                        /// QR eddyviscosity:
                        real eddyR = -(Dxy * Dxy + Dxz * Dxz + c1o3 * dxux * dxux) * (dxux) - (Dxy * Dxy + Dyz * Dyz + c1o3 * dyuy * dyuy) * dyuy - (Dxz * Dxz + Dyz * Dyz + c1o3 * dzuz * dzuz) * dzuz - c2o1 * Dxy * Dxz * Dyz;
                        real eddyQ = Dxy * Dxz + Dxy * Dyz + Dxz * Dyz + c1o2 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz);
                        //real nuEddy = 5.0e5 * (eddyR / (eddyQ + 1e-100)) * (dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);
                        real nuEddy = 5.0e3 * (eddyR / (eddyQ + 1e-100));
                        nuEddy = (nuEddy < c1o1 / collFactorM) ? c1o1 / collFactorM : nuEddy;
                        collFactorM = c1o1 / nuEddy;
                        // collFactorM = c1o1 / (c1o1 / collFactorM +1.e2*nuEddy*(dX1_phi*dX1_phi+dX2_phi*dX2_phi+dX3_phi*dX3_phi));
                        collFactorM = (collFactorM < 1.8) ? 1.8 : collFactorM;
                        OxyyPxzz = 8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
                        OxyyMxzz = 8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
                        Oxyz = 24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 *
                        collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0)); A = (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM
                        - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM)); BB = (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM
                        - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
                    }

					//if ((phi[DIR_000] > c1o2)&& (normX1 * vvx + normX2 * vvy + normX3 * vvz < 0)){//&& ((*phaseFieldOld)(x1, x2, x3) <= c1o2)) {
     //               if ((phi[DIR_000] > 0.01) && (phi[DIR_000]<0.99)){
					//	//std::cout << "new node\n";
					//	///QR eddyviscosity:
					//	real eddyR = -(Dxy * Dxy + Dxz * Dxz + c1o3 * dxux * dxux) * (dxux)-(Dxy * Dxy + Dyz * Dyz + c1o3 * dyuy * dyuy) * dyuy - (Dxz * Dxz + Dyz * Dyz + c1o3 * dzuz * dzuz) * dzuz - c2o1 * Dxy * Dxz * Dyz;
					//	real eddyQ = Dxy * Dxz + Dxy * Dyz + Dxz * Dyz + c1o2 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz);
					//	real nuEddy = 10.0e4*(eddyR / (eddyQ + 1e-100)) * (dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);
     //                   nuEddy = 1000*(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);
					//	//nuEddy=10.0e4*fabsf(dxux+dyuy+dzuz) * (dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi);
					//	//if (nuEddy > c1o1 / collFactorM) std::cout << nuEddy <<" "<< fabsf(dxux + dyuy + dzuz)<< "\n";
					//	nuEddy = (nuEddy < c1o1 / collFactorM) ? c1o1 / collFactorM : nuEddy;
					//	collFactorM = c1o1 / nuEddy;
					//	//collFactorM = 1.8;
					//	//collFactorM = c1o1 / (c1o1 / collFactorM +1.e2*nuEddy*(dX1_phi*dX1_phi+dX2_phi*dX2_phi+dX3_phi*dX3_phi));
					//	collFactorM = (collFactorM < 1.8) ? 1.8 : collFactorM;
					//	OxyyPxzz = 8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
					//	OxyyMxzz = 8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
					//	Oxyz = 24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 * collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0));
					//	A = (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//	BB = (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));

					//}
					///////

                    // non Newtonian fluid collision factor
                    if (phi[DIR_000] > phiLim) {
                        real shearRate = sqrt(c2o1 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz);
                        collFactorM = Rheology::getBinghamCollFactor(collFactorM, shearRate, c1o1);
                        collFactorM = (collFactorM < c1o1) ? c1o1 : collFactorM;
                    }

					real mxxMyyh = -c2o1 * (dxux - dyuy) / collFactorMInv * c1o3;
					real mxxMzzh = -c2o1 * (dxux - dzuz) / collFactorMInv * c1o3;
//					mfhbba = -Dxy / collFactorMInv*c1o3;
//					mfhbab = -Dxz / collFactorMInv * c1o3;
//					mfhabb = -Dyz / collFactorMInv * c1o3;

//					// attempt to improve implicit  stress computation by fixed iteration
//					real dX2_rho = (rhoToPhi)*dX2_phi;
//					real dX1_rho = (rhoToPhi)*dX1_phi;
//					real dX3_rho = (rhoToPhi)*dX3_phi;
//
//						real dfx= c1o3 * (c1o1 / collFactorM - c1o2) *(2 * dxux * dX1_rho + Dxy * dX2_rho + Dxz * dX3_rho) / (rho);
//						real dfy = c1o3 * (c1o1 / collFactorM - c1o2) *(Dxy * dX1_rho + 2 * dyuy * dX2_rho + Dyz * dX3_rho) / (rho);
//						real dfz = c1o3 * (c1o1 / collFactorM - c1o2) *(Dxz * dX1_rho + Dyz * dX2_rho + 2 * dyuy * dX3_rho) / (rho);
//
//						for (int iteration = 0; iteration < 5; iteration++) {
//							mxxMyy = (mfcaa - dfx * dfx * c1o2) - (mfaca - dfy * dfy * c1o2);
//							mxxMzz = (mfcaa - dfx * dfx * c1o2) - (mfaac - dfz * dfz * c1o2);
//						}
/////end fixed iteration
//


					//relax
					mxxPyyPzz += OxxPyyPzz * (/*mfaaa*/ - mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
					mxxMyy += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vy2 * dyuy);
					mxxMzz += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactorM) * (vx2 * dxux - vz2 * dzuz);

					mfabb += collFactorM * (-mfabb);
					mfbab += collFactorM * (-mfbab);
					mfbba += collFactorM * (-mfbba);

//					mfhaaa = (phi[DIR_000] < c1o2) ? mfaaa * rhoL / rhoH : mfaaa * rhoL / rhoH;
					mxxMyyh += collFactorMInv * (-mxxMyyh) - 3. * (1. - c1o2 * collFactorMInv) * (vx2 * dxux - vy2 * dyuy);
					mxxMzzh += collFactorMInv * (-mxxMzzh) - 3. * (1. - c1o2 * collFactorMInv) * (vx2 * dxux - vz2 * dzuz);

					//mfhabb += collFactorMInv * (-mfhabb);
					//mfhbab += collFactorMInv * (-mfhbab);
					//mfhbba += collFactorMInv * (-mfhbba);

					//mfhcaa = c1o3 * (mxxMyyh + mxxMzzh + mfhaaa);
					//mfhaca = c1o3 * (-2. * mxxMyyh + mxxMzzh + mfhaaa);
					//mfhaac = c1o3 * (mxxMyyh - 2. * mxxMzzh + mfhaaa);


					//if (fabsf(mfaaa + (dxux + dyuy + dzuz) > 1e-9)){
					//	std::cout << mfaaa <<" "<< (dxux + dyuy + dzuz)<< std::endl;
					//}


					////updated pressure
					//mfaaa += (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling;
					//mfaaa *= (c1o1-omegaDRho);// (mfaaa + (dxux + dyuy + dzuz)) * .5; // Pressure elimination as in standard velocity model
								 //  mfaaa += (rho - c1o1) * (dxux + dyuy + dzuz);
				
					mxxPyyPzz += mfaaa; // 12.03.21 shifted by mfaaa

										// mxxPyyPzz += (mfaaa + mfaaaS) * c1o2;
										//mfaaa = mfaaaS;
										// linear combinations back
					mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
					mfaca = c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz);
					mfaac = c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz);

					//3.
					// linear combinations
					real mxxyPyzz = mfcba + mfabc;
					real mxxyMyzz = mfcba - mfabc;

					real mxxzPyyz = mfcab + mfacb;
					real mxxzMyyz = mfcab - mfacb;

					real mxyyPxzz = mfbca + mfbac;
					real mxyyMxzz = mfbca - mfbac;

					 mmfcaa += c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz- mfaaa);
					 mmfaca += c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz- mfaaa);
					 mmfaac += c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz- mfaaa);
					 mmfabb += mfabb;
					 mmfbab += mfbab;
					 mmfbba += mfbba;

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
					CUMacc = -O4 * (c1o1 / collFactorM - c1o2) * (dyuy + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMacc);
					CUMcac = -O4 * (c1o1 / collFactorM - c1o2) * (dxux + dzuz) * c2o3 * A + (c1o1 - O4) * (CUMcac);
					CUMcca = -O4 * (c1o1 / collFactorM - c1o2) * (dyuy + dxux) * c2o3 * A + (c1o1 - O4) * (CUMcca);
					CUMbbc = -O4 * (c1o1 / collFactorM - c1o2) * Dxy * c1o3 * BB + (c1o1 - O4) * (CUMbbc);
					CUMbcb = -O4 * (c1o1 / collFactorM - c1o2) * Dxz * c1o3 * BB + (c1o1 - O4) * (CUMbcb);
					CUMcbb = -O4 * (c1o1 / collFactorM - c1o2) * Dyz * c1o3 * BB + (c1o1 - O4) * (CUMcbb);

					//5.
					CUMbcc += O5 * (-CUMbcc);
					CUMcbc += O5 * (-CUMcbc);
					CUMccb += O5 * (-CUMccb);

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

					mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho;
					mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho;
					mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9 * (oMdrho - c1o1) * oMdrho;

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


					////////save central moments for the phase field
					real MMxx = mfcaa - c1o3 * mfaaa;
					real MMyy = mfaca - c1o3 * mfaaa;
					real MMzz = mfaac - c1o3 * mfaaa;
					real MMxy = mfbba;
					real MMxz = mfbab;
					real MMyz = mfabb;


					////////////////////////////////////////////////////////////////////////////////////
					//forcing
					mfbaa = -mfbaa;// *(c1o1 - forcingBIAS) / forcingBIAS;
					mfaba = -mfaba;// *(c1o1 - forcingBIAS) / forcingBIAS;
					mfaab = -mfaab;// *(c1o1 - forcingBIAS) / forcingBIAS;


					//mfhbaa = mfbaa;
					//mfhaba = mfaba;
					//mfhaab = mfaab;

					//mfhcba = 0.;
					//mfhabc = 0.;
					//mfhcab = 0.;
					//mfhacb = 0.;
					//mfhbca = 0.;
					//mfhbac = 0.;
					//mfhbbb = 0.;

					//real oMdrhoInv = (rhoRef - (mfhaaa )) / rhoRef;
					//mfhcbb =/* CUMcbb + */((mfhcaa + c1o3) * mfhabb + 2. * mfhbba * mfhbab);
					//mfhbcb =/* CUMbcb + */((mfhaca + c1o3) * mfhbab + 2. * mfhbba * mfhabb);
					//mfhbbc =/* CUMbbc + */((mfhaac + c1o3) * mfhbba + 2. * mfhbab * mfhabb);

					//mfhcca = /*CUMcca + */(mfhcaa * mfhaca + 2. * mfhbba * mfhbba) + c1o3 * (mfhcaa + mfhaca) * oMdrhoInv + c1o9 * (oMdrhoInv - c1o1) * oMdrhoInv;
					//mfhcac = /*CUMcac + */(mfhcaa * mfhaac + 2. * mfhbab * mfhbab) + c1o3 * (mfhcaa + mfhaac) * oMdrhoInv + c1o9 * (oMdrhoInv - c1o1) * oMdrhoInv;
					//mfhacc = /*CUMacc + */(mfhaac * mfhaca + 2. * mfhabb * mfhabb) + c1o3 * (mfhaac + mfhaca) * oMdrhoInv + c1o9 * (oMdrhoInv - c1o1) * oMdrhoInv;

					////5.
					//mfhbcc = /*CUMbcc +*/ (mfhaac * mfhbca + mfhaca * mfhbac + 4. * mfhabb * mfhbbb + 2. * (mfhbab * mfhacb + mfhbba * mfhabc)) + c1o3 * (mfhbca + mfhbac) * oMdrhoInv;
					//mfhcbc = /*CUMcbc +*/ (mfhaac * mfhcba + mfhcaa * mfhabc + 4. * mfhbab * mfhbbb + 2. * (mfhabb * mfhcab + mfhbba * mfhbac)) + c1o3 * (mfhcba + mfhabc) * oMdrhoInv;
					//mfhccb = /*CUMccb +*/ (mfhcaa * mfhacb + mfhaca * mfhcab + 4. * mfhbba * mfhbbb + 2. * (mfhbab * mfhbca + mfhabb * mfhcba)) + c1o3 * (mfhacb + mfhcab) * oMdrhoInv;

					////6.
					//mfhccc = /*CUMccc */- ((-4. * mfhbbb * mfhbbb
					//	- (mfhcaa * mfhacc + mfhaca * mfhcac + mfhaac * mfhcca)
					//	- 4. * (mfhabb * mfhcbb + mfhbac * mfhbca + mfhbba * mfhbbc)
					//	- 2. * (mfhbca * mfhbac + mfhcba * mfhabc + mfhcab * mfhacb))
					//	+ (4. * (mfhbab * mfhbab * mfhaca + mfhabb * mfhabb * mfhcaa + mfhbba * mfhbba * mfhaac)
					//		+ 2. * (mfhcaa * mfhaca * mfhaac)
					//		+ 16. * mfhbba * mfhbab * mfhabb)
					//	- c1o3 * (mfhacc + mfhcac + mfhcca) * oMdrhoInv - c1o9 * oMdrhoInv * oMdrhoInv
					//	- c1o9 * (mfhcaa + mfhaca + mfhaac) * oMdrhoInv * (1. - 2. * oMdrhoInv) - c1o27 * oMdrhoInv * oMdrhoInv * (-2. * oMdrhoInv)
					//	+ (2. * (mfhbab * mfhbab + mfhabb * mfhabb + mfhbba * mfhbba)
					//		+ (mfhaac * mfhaca + mfhaac * mfhcaa + mfhaca * mfhcaa)) * c2o3 * oMdrhoInv) - c1o27 * oMdrhoInv;





					backwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1, oMdrho);
					backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
					backwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, vvx, vx2, c3o1, c1o3, oMdrho);
					backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
					backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
					backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
					backwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, vvx, vx2, c3o1, c1o3, oMdrho);
					backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
					backwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, vvx, vx2, c9o1, c1o9, oMdrho);

					////////////////////////////////////////////////////////////////////////////////////
					// Y - Dir
					backwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, vvy, vy2, c6o1, c1o6, oMdrho);
					backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
					backwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, vvy, vy2, c18o1, c1o18, oMdrho);
					backwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, vvy, vy2, c3o2, c2o3, oMdrho);
					backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
					backwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, vvy, vy2, c9o2, c2o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, vvy, vy2, c6o1, c1o6, oMdrho);
					backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
					backwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, vvy, vy2, c18o1, c1o18, oMdrho);

					////////////////////////////////////////////////////////////////////////////////////
					// Z - Dir
					backwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, vvz, vz2, c36o1, c1o36, oMdrho);
					backwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, vvz, vz2, c9o1, c1o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, vvz, vz2, c36o1, c1o36, oMdrho);
					backwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, vvz, vz2, c9o1, c1o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, vvz, vz2, c9o4, c4o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, vvz, vz2, c9o1, c1o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, vvz, vz2, c36o1, c1o36, oMdrho);
					backwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, vvz, vz2, c9o1, c1o9, oMdrho);
					backwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, vvz, vz2, c36o1, c1o36, oMdrho);





					//backwardInverseChimeraWithKincompressible(mfhaaa, mfhbaa, mfhcaa, vvx, vx2, c1o1, c1o1, oMdrhoInv);
					//backwardChimera(mfhaba, mfhbba, mfhcba, vvx, vx2);
					//backwardInverseChimeraWithKincompressible(mfhaca, mfhbca, mfhcca, vvx, vx2, c3o1, c1o3, oMdrhoInv);
					//backwardChimera(mfhaab, mfhbab, mfhcab, vvx, vx2);
					//backwardChimera(mfhabb, mfhbbb, mfhcbb, vvx, vx2);
					//backwardChimera(mfhacb, mfhbcb, mfhccb, vvx, vx2);
					//backwardInverseChimeraWithKincompressible(mfhaac, mfhbac, mfhcac, vvx, vx2, c3o1, c1o3, oMdrhoInv);
					//backwardChimera(mfhabc, mfhbbc, mfhcbc, vvx, vx2);
					//backwardInverseChimeraWithKincompressible(mfhacc, mfhbcc, mfhccc, vvx, vx2, c9o1, c1o9, oMdrhoInv);

					//////////////////////////////////////////////////////////////////////////////////////
					//// Y - Dir
					//backwardInverseChimeraWithKincompressible(mfhaaa, mfhaba, mfhaca, vvy, vy2, c6o1, c1o6, oMdrhoInv);
					//backwardChimera(mfhaab, mfhabb, mfhacb, vvy, vy2);
					//backwardInverseChimeraWithKincompressible(mfhaac, mfhabc, mfhacc, vvy, vy2, c18o1, c1o18, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhbaa, mfhbba, mfhbca, vvy, vy2, c3o2, c2o3, oMdrhoInv);
					//backwardChimera(mfhbab, mfhbbb, mfhbcb, vvy, vy2);
					//backwardInverseChimeraWithKincompressible(mfhbac, mfhbbc, mfhbcc, vvy, vy2, c9o2, c2o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhcaa, mfhcba, mfhcca, vvy, vy2, c6o1, c1o6, oMdrhoInv);
					//backwardChimera(mfhcab, mfhcbb, mfhccb, vvy, vy2);
					//backwardInverseChimeraWithKincompressible(mfhcac, mfhcbc, mfhccc, vvy, vy2, c18o1, c1o18, oMdrhoInv);

					//////////////////////////////////////////////////////////////////////////////////////
					//// Z - Dir
					//backwardInverseChimeraWithKincompressible(mfhaaa, mfhaab, mfhaac, vvz, vz2, c36o1, c1o36, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhaba, mfhabb, mfhabc, vvz, vz2, c9o1, c1o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhaca, mfhacb, mfhacc, vvz, vz2, c36o1, c1o36, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhbaa, mfhbab, mfhbac, vvz, vz2, c9o1, c1o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhbba, mfhbbb, mfhbbc, vvz, vz2, c9o4, c4o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhbca, mfhbcb, mfhbcc, vvz, vz2, c9o1, c1o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhcaa, mfhcab, mfhcac, vvz, vz2, c36o1, c1o36, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhcba, mfhcbb, mfhcbc, vvz, vz2, c9o1, c1o9, oMdrhoInv);
					//backwardInverseChimeraWithKincompressible(mfhcca, mfhccb, mfhccc, vvz, vz2, c36o1, c1o36, oMdrhoInv);

					/////////////////////

					//////////////////////////////////////////////////////////////////////////
					//proof correctness
					//////////////////////////////////////////////////////////////////////////
					//#ifdef  PROOF_CORRECTNESS
					real rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
						+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
						+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
					//			   //real dif = fabs(drho - rho_post);
					//               real dif = drho + (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling - rho_post;
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

					if (UbMath::isNaN(rho_post) || UbMath::isInfinity(rho_post))
						UB_THROW(UbException(UB_EXARGS, "rho_post is not a number (nan or -1.#IND) or infinity number -1.#INF, node=" + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3) + ",phi=" + UbSystem::toString(phi[DIR_000])));

					//////////////////////////////////////////////////////////////////////////
					//write distribution
					//////////////////////////////////////////////////////////////////////////
				//	if (phi[DIR_000] < c1o2) {
						(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = mfabb;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = mfbab;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = mfbba;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = mfaab;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = mfcab;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = mfaba;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = mfcba;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = mfbaa;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = mfbca;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = mfaaa;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = mfcaa;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = mfaca;//* rho * c1o3;
						(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac;//* rho * c1o3;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac;//* rho * c1o3;

						(*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;// *rho* c1o3;


						//(*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3) = mfhabb;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3) = mfhbab;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3) = mfhbba;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3) = mfhaab;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3) = mfhcab;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3) = mfhaba;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3) = mfhcba;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3) = mfhbaa;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3) = mfhbca;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3) = mfhaaa;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3) = mfhcaa;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3) = mfhaca;//* rho * c1o3;
						//(*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfhcca;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3) = mfhcbb;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3) = mfhbcb;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p) = mfhbbc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3) = mfhccb;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3) = mfhacb;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p) = mfhcbc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p) = mfhabc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p) = mfhbcc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p) = mfhbac;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfhccc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfhacc;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfhcac;//* rho * c1o3;
						//(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p) = mfhaac;//* rho * c1o3;

						//(*this->zeroDistributionsH2)(x1, x2, x3) = mfhbbb;// *rho* c1o3;
				//	}


				//	else {
					//	(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3)         = mfhabb;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3)         = mfhbab;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)         = mfhbba;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3)        = mfhaab;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3)       = mfhcab;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3)        = mfhaba;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3)       = mfhcba;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3)        = mfhbaa;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3)       = mfhbca;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3)       = mfhaaa;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3)      = mfhcaa;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3)      = mfhaca;//* rho * c1o3;
					//	(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3)     = mfhcca;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)     = mfhcbb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3)     = mfhbcb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p)     = mfhbbc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfhccb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3)    = mfhacb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfhcbc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p)    = mfhabc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfhbcc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p)    = mfhbac;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfhccc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfhacc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfhcac;//* rho * c1o3;
					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfhaac;//* rho * c1o3;

					//	(*this->zeroDistributionsF)(x1, x2, x3) = mfhbbb;// *rho* c1o3;


					//	(*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3)         = mfabb;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3)         = mfbab;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3)         = mfbba;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3)        = mfaab;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3)       = mfcab;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3)        = mfaba;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3)       = mfcba;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3)        = mfbaa;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3)       = mfbca;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3)       = mfaaa;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3)      = mfcaa;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3)      = mfaca;//* rho * c1o3;
					//	(*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3)     = mfcca;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3)     = mfcbb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3)     = mfbcb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p)     = mfbbc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfccb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3)    = mfacb;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfcbc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p)    = mfabc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfbcc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p)    = mfbac;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfacc;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfcac;//* rho * c1o3;
					//	(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfaac;//* rho * c1o3;

					//	(*this->zeroDistributionsH2)(x1, x2, x3) = mfbbb;// *rho* c1o3;
					//}
																	// !Old Kernel
/////////////////////  P H A S E - F I E L D   S O L V E R
////////////////////////////////////////////
/////CUMULANT PHASE-FIELD
					real omegaD =1.0/( 3.0 * mob + 0.5);
					{
						mfcbb = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);
						mfbcb = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);
						mfbbc = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);
						mfccb = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);
						mfacb = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);
						mfcbc = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);
						mfabc = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);
						mfbcc = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);
						mfbac = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);
						mfccc = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);
						mfacc = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);
						mfcac = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);
						mfaac = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);
						mfabb = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);
						mfbab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);
						mfbba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);
						mfaab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);
						mfcab = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);
						mfaba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);
						mfcba = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);
						mfbaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);
						mfbca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);
						mfaaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);
						mfcaa = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);
						mfaca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);
						mfcca = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);
						mfbbb = (*this->zeroDistributionsH1)(x1, x2, x3);


						////////////////////////////////////////////////////////////////////////////////////
						//! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
						//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
						//!
						////////////////////////////////////////////////////////////////////////////////////
						// second component
						real concentration =
							((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
								(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
								((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
						////////////////////////////////////////////////////////////////////////////////////
						real oneMinusRho = c1o1- concentration;

						real cx =
							((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
								(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
								(mfcbb - mfabb));
						real cy =
							((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
								(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
								(mfbcb - mfbab));
						real cz =
							((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
								(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
								(mfbbc - mfbba));

						////////////////////////////////////////////////////////////////////////////////////
						// calculate the square of velocities for this lattice node
						real cx2 = cx * cx;
						real cy2 = cy * cy;
						real cz2 = cz * cz;
						////////////////////////////////////////////////////////////////////////////////////
						//! - Chimera transform from well conditioned distributions to central moments as defined in Appendix J in \ref
						//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
						//! see also Eq. (6)-(14) in \ref
						//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
						//!
						////////////////////////////////////////////////////////////////////////////////////
						// Z - Dir
						forwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, cz, cz2, c36o1, c1o36, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, cz, cz2, c9o1, c1o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, cz, cz2, c36o1, c1o36, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, cz, cz2, c9o1, c1o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, cz, cz2, c9o4, c4o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, cz, cz2, c9o1, c1o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, cz, cz2, c36o1, c1o36, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, cz, cz2, c9o1, c1o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, cz, cz2, c36o1, c1o36, oneMinusRho);

						////////////////////////////////////////////////////////////////////////////////////
						// Y - Dir
						forwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, cy, cy2, c6o1, c1o6, oneMinusRho);
						forwardChimera(mfaab, mfabb, mfacb, cy, cy2);
						forwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, cy, cy2, c18o1, c1o18, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, cy, cy2, c3o2, c2o3, oneMinusRho);
						forwardChimera(mfbab, mfbbb, mfbcb, cy, cy2);
						forwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, cy, cy2, c9o2, c2o9, oneMinusRho);
						forwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, cy, cy2, c6o1, c1o6, oneMinusRho);
						forwardChimera(mfcab, mfcbb, mfccb, cy, cy2);
						forwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, cy, cy2, c18o1, c1o18, oneMinusRho);

						////////////////////////////////////////////////////////////////////////////////////
						// X - Dir
						forwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, cx, cx2, c1o1, c1o1, oneMinusRho);
						forwardChimera(mfaba, mfbba, mfcba, cx, cx2);
						forwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, cx, cx2, c3o1, c1o3, oneMinusRho);
						forwardChimera(mfaab, mfbab, mfcab, cx, cx2);
						forwardChimera(mfabb, mfbbb, mfcbb, cx, cx2);
						forwardChimera(mfacb, mfbcb, mfccb, cx, cx2);
						forwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, cx, cx2, c3o1, c1o3, oneMinusRho);
						forwardChimera(mfabc, mfbbc, mfcbc, cx, cx2);
						forwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, cx, cx2, c3o1, c1o9, oneMinusRho);

						////////////////////////////////////////////////////////////////////////////////////
						//! - experimental Cumulant ... to be published ... hopefully
						//!

						// linearized orthogonalization of 3rd order central moments
						real Mabc = mfabc - mfaba * c1o3;
						real Mbca = mfbca - mfbaa * c1o3;
						real Macb = mfacb - mfaab * c1o3;
						real Mcba = mfcba - mfaba * c1o3;
						real Mcab = mfcab - mfaab * c1o3;
						real Mbac = mfbac - mfbaa * c1o3;
						// linearized orthogonalization of 5th order central moments
						real Mcbc = mfcbc - mfaba * c1o9;
						real Mbcc = mfbcc - mfbaa * c1o9;
						real Mccb = mfccb - mfaab * c1o9;

						//25.03.2023 mixed normals
						real MomX1 = vvx * concentration - cx;
						real MomX2 = vvy * concentration - cy;
						real MomX3 = vvz * concentration - cz;
						real mixNormal = 0.5;

						real MomXDenom = sqrt(MomX1 * MomX1 + MomX2 * MomX2 + MomX3 * MomX3)+1.0e-100;
						real scaleNorm = (normX1 * MomX1 + normX2 * MomX2 + normX3 * MomX3) / MomXDenom;
						scaleNorm = scaleNorm * scaleNorm;// *scaleNorm* scaleNorm;

						normX1 = (normX1 * (c1o1 - mixNormal) + mixNormal * MomX1 / MomXDenom )* scaleNorm;
						normX2 = (normX2 * (c1o1 - mixNormal) + mixNormal * MomX2 / MomXDenom )* scaleNorm;
						normX3 = (normX3 * (c1o1 - mixNormal) + mixNormal * MomX3 / MomXDenom )* scaleNorm;

						//31.05.2022 addaptive mobility
						//omegaD = c1o1 + (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration))) / (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration)) + fabs((1.0 - concentration) * (concentration)) * c1o6 * oneOverInterfaceScale+1.0e-200);
						//omegaD = c2o1 * (concentration * (concentration - c1o1)) / (-c6o1 * (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration))) + (concentration * (concentration - c1o1))+1.0e-200);
						// collision of 1st order moments
						cx = cx * (c1o1 - omegaD) + omegaD * vvx * concentration +
							normX1 * (c1o1 - 0.5 * omegaD) * (1.0 - concentration) * (concentration) * c1o3 * oneOverInterfaceScale;
						cy = cy * (c1o1 - omegaD) + omegaD * vvy * concentration +
							normX2 * (c1o1 - 0.5 * omegaD) * (1.0 - concentration) * (concentration) * c1o3 * oneOverInterfaceScale;
						cz = cz * (c1o1 - omegaD) + omegaD * vvz * concentration +
							normX3 * (c1o1 - 0.5 * omegaD) * (1.0 - concentration) * (concentration) * c1o3 * oneOverInterfaceScale;

						cx2 = cx * cx;
						cy2 = cy * cy;
						cz2 = cz * cz;

						//// equilibration of 2nd order moments
						mfbba = c0o1;
						mfbab = c0o1;
						mfabb = c0o1;

						mfcaa = c1o3 * concentration;
						mfaca = c1o3 * concentration;
						mfaac = c1o3 * concentration;

						//take second moment from fluid
						//mfbba = concentration*MMxy;
						//mfbab = concentration*MMxz;
						//mfabb = concentration*MMyz;

						//mfcaa = (c1o3+MMxx) * concentration;
						//mfaca = (c1o3+MMyy) * concentration;
						//mfaac = (c1o3+MMzz) * concentration;



						// equilibration of 3rd order moments
						Mabc = c0o1;
						Mbca = c0o1;
						Macb = c0o1;
						Mcba = c0o1;
						Mcab = c0o1;
						Mbac = c0o1;
						mfbbb = c0o1;

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

						mfcbb = c0o1;
						mfbcb = c0o1;
						mfbbc = c0o1;

						// equilibration of 5th order moments
						Mcbc = c0o1;
						Mbcc = c0o1;
						Mccb = c0o1;

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
						backwardInverseChimeraWithKincompressible(mfaaa, mfbaa, mfcaa, cx, cx2, c1o1, c1o1, oneMinusRho);
						backwardChimera(mfaba, mfbba, mfcba, cx, cx2);
						backwardInverseChimeraWithKincompressible(mfaca, mfbca, mfcca, cx, cx2, c3o1, c1o3, oneMinusRho);
						backwardChimera(mfaab, mfbab, mfcab, cx, cx2);
						backwardChimera(mfabb, mfbbb, mfcbb, cx, cx2);
						backwardChimera(mfacb, mfbcb, mfccb, cx, cx2);
						backwardInverseChimeraWithKincompressible(mfaac, mfbac, mfcac, cx, cx2, c3o1, c1o3, oneMinusRho);
						backwardChimera(mfabc, mfbbc, mfcbc, cx, cx2);
						backwardInverseChimeraWithKincompressible(mfacc, mfbcc, mfccc, cx, cx2, c9o1, c1o9, oneMinusRho);

						////////////////////////////////////////////////////////////////////////////////////
						// Y - Dir
						backwardInverseChimeraWithKincompressible(mfaaa, mfaba, mfaca, cy, cy2, c6o1, c1o6, oneMinusRho);
						backwardChimera(mfaab, mfabb, mfacb, cy, cy2);
						backwardInverseChimeraWithKincompressible(mfaac, mfabc, mfacc, cy, cy2, c18o1, c1o18, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfbaa, mfbba, mfbca, cy, cy2, c3o2, c2o3, oneMinusRho);
						backwardChimera(mfbab, mfbbb, mfbcb, cy, cy2);
						backwardInverseChimeraWithKincompressible(mfbac, mfbbc, mfbcc, cy, cy2, c9o2, c2o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfcaa, mfcba, mfcca, cy, cy2, c6o1, c1o6, oneMinusRho);
						backwardChimera(mfcab, mfcbb, mfccb, cy, cy2);
						backwardInverseChimeraWithKincompressible(mfcac, mfcbc, mfccc, cy, cy2, c18o1, c1o18, oneMinusRho);

						////////////////////////////////////////////////////////////////////////////////////
						// Z - Dir
						backwardInverseChimeraWithKincompressible(mfaaa, mfaab, mfaac, cz, cz2, c36o1, c1o36, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfaba, mfabb, mfabc, cz, cz2, c9o1, c1o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfaca, mfacb, mfacc, cz, cz2, c36o1, c1o36, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfbaa, mfbab, mfbac, cz, cz2, c9o1, c1o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfbba, mfbbb, mfbbc, cz, cz2, c9o4, c4o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfbca, mfbcb, mfbcc, cz, cz2, c9o1, c1o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfcaa, mfcab, mfcac, cz, cz2, c36o1, c1o36, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfcba, mfcbb, mfcbc, cz, cz2, c9o1, c1o9, oneMinusRho);
						backwardInverseChimeraWithKincompressible(mfcca, mfccb, mfccc, cz, cz2, c36o1, c1o36, oneMinusRho);



						(*this->localDistributionsH1)(D3Q27System::ET_E,   x1,  x2,  x3) = mfabb;
						(*this->localDistributionsH1)(D3Q27System::ET_N,   x1,  x2,  x3) = mfbab;
						(*this->localDistributionsH1)(D3Q27System::ET_T,   x1,  x2,  x3) = mfbba;
						(*this->localDistributionsH1)(D3Q27System::ET_NE,  x1,  x2,  x3) = mfaab;
						(*this->localDistributionsH1)(D3Q27System::ET_NW,  x1p, x2,  x3) = mfcab;
						(*this->localDistributionsH1)(D3Q27System::ET_TE,  x1,  x2,  x3) = mfaba;
						(*this->localDistributionsH1)(D3Q27System::ET_TW,  x1p, x2,  x3) = mfcba;
						(*this->localDistributionsH1)(D3Q27System::ET_TN,  x1,  x2,  x3) = mfbaa;
						(*this->localDistributionsH1)(D3Q27System::ET_TS,  x1,  x2p, x3) = mfbca;
						(*this->localDistributionsH1)(D3Q27System::ET_TNE, x1,  x2,  x3) = mfaaa;
						(*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2,  x3) = mfcaa;
						(*this->localDistributionsH1)(D3Q27System::ET_TSE, x1,  x2p, x3) = mfaca;
						(*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

						(*this->nonLocalDistributionsH1)(D3Q27System::ET_W,   x1p, x2,  x3 ) = mfcbb;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_S,   x1,  x2p, x3 ) = mfbcb;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_B,   x1,  x2,  x3p) = mfbbc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = mfccb;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = mfacb;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BW,  x1p, x2,  x3p) = mfcbc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BE,  x1,  x2,  x3p) = mfabc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BS,  x1,  x2p, x3p) = mfbcc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BN,  x1,  x2,  x3p) = mfbac;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1,  x2p, x3p) = mfacc;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2,  x3p) = mfcac;
						(*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1,  x2,  x3p) = mfaac;

						(*this->zeroDistributionsH1)(x1,x2,x3) = mfbbb;




					}
				}
			}
		}
	}
	//Set multiphase BCs
	//for (int x3 = minX3; x3 < maxX3; x3++) {
	//	for (int x2 = minX2; x2 < maxX2; x2++) {
	//		for (int x1 = minX1; x1 < maxX1; x1++) {
	//			if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {
	//				int x1p = x1 + 1;
	//				int x2p = x2 + 1;
	//				int x3p = x3 + 1;
	//				findNeighbors(phaseField, x1, x2, x3);
	//				if ((phi[DIR_000] > c1o2) && (
	//					(phi[DIR_P00] <= c1o2) ||
	//					(phi[DIR_M00] <= c1o2) ||
	//					(phi[DIR_00P] <= c1o2) ||
	//					(phi[DIR_00M] <= c1o2) ||
	//					(phi[DIR_0M0] <= c1o2) ||
	//					(phi[DIR_0P0] <= c1o2) ||
	//					(phi[DIR_PP0] <= c1o2) ||
	//					(phi[DIR_PM0] <= c1o2) ||
	//					(phi[DIR_P0P] <= c1o2) ||
	//					(phi[DIR_P0M] <= c1o2) ||
	//					(phi[DIR_MP0] <= c1o2) ||
	//					(phi[DIR_MM0] <= c1o2) ||
	//					(phi[DIR_M0P] <= c1o2) ||
	//					(phi[DIR_M0M] <= c1o2) ||
	//					(phi[DIR_0PM] <= c1o2) ||
	//					(phi[DIR_0MM] <= c1o2) ||
	//					(phi[DIR_0PP] <= c1o2) ||
	//					(phi[DIR_0MP] <= c1o2) ||
	//					(phi[DIR_PPP] <= c1o2) ||
	//					(phi[DIR_PMP] <= c1o2) ||
	//					(phi[DIR_MPP] <= c1o2) ||
	//					(phi[DIR_MMP] <= c1o2) ||
	//					(phi[DIR_PPM] <= c1o2) ||
	//					(phi[DIR_PMM] <= c1o2) ||
	//					(phi[DIR_MPM] <= c1o2) ||
	//					(phi[DIR_MMM] <= c1o2)
	//					)) {
	//					//real mfabb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);//* rho * c1o3;
	//					//real mfbab = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);//* rho * c1o3;
	//					//real mfbba = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);//* rho * c1o3;
	//					//real mfaab = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);//* rho * c1o3;
	//					//real mfcab = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);//* rho * c1o3;
	//					//real mfaba = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);//* rho * c1o3;
	//					//real mfcba = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);//* rho * c1o3;
	//					//real mfbaa = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);//* rho * c1o3;
	//					//real mfbca = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);//* rho * c1o3;
	//					//real mfaaa = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);//* rho * c1o3;
	//					//real mfcaa = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);//* rho * c1o3;
	//					//real mfaca = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);//* rho * c1o3;
	//					//real mfcca = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);//* rho * c1o3;
	//					//real mfcbb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);//* rho * c1o3;
	//					//real mfbcb = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);//* rho * c1o3;
	//					//real mfbbc = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);//* rho * c1o3;
	//					//real mfccb = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);//* rho * c1o3;
	//					//real mfacb = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);//* rho * c1o3;
	//					//real mfcbc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);//* rho * c1o3;
	//					//real mfabc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);//* rho * c1o3;
	//					//real mfbcc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);//* rho * c1o3;
	//					//real mfbac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);//* rho * c1o3;
	//					//real mfccc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);//* rho * c1o3;
	//					//real mfacc = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);//* rho * c1o3;
	//					//real mfcac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);//* rho * c1o3;
	//					//real mfaac = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);//* rho * c1o3;
	//					//real mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);
	//					
	//					real vx = (*vxNode)(x1, x2, x3);
	//					real vy = (*vyNode)(x1, x2, x3);
	//					real vz = (*vzNode)(x1, x2, x3);
	//					SPtr<DistributionArray3D> distribution = this->getDataSet()->getFdistributions();
	//					real ff[27];
	//					distribution->getDistributionInv(ff, x1, x2 , x3 );

	//					for (int fdir = D3Q27System::STARTF; fdir <= D3Q27System::ENDF; fdir++) {
	//						if ((phi[fdir] <= c1o2)) {
	//							real rhoG = (*rhoNode)(x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir]);
	//							real ftemp= D3Q27System::getCompFeqForDirection(D3Q27System::INVDIR[fdir], rhoG, vx, vy, vz) + D3Q27System::getCompFeqForDirection(fdir, rhoG, vx, vy, vz);
	//							real fBB;
	//							fBB=distribution->getDistributionInvForDirection( x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir],fdir);
	//							distribution->setDistributionInvForDirection(ftemp - ff[fdir], x1 + D3Q27System::DX1[fdir], x2 + D3Q27System::DX2[fdir], x3 + D3Q27System::DX3[fdir], fdir);
	//							distribution->setDistributionForDirection(fBB-c6o1*D3Q27System::WEIGTH[fdir] * (-vx * D3Q27System::DX1[fdir] - vy * D3Q27System::DX2[fdir] - vz * D3Q27System::DX3[fdir]), x1, x2, x3, fdir);
	//						}
	//					}

	//					//if ((phi[DIR_P00] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(W, rhoG, vx,vy,vz )+ D3Q27System::getCompFeqForDirection(E, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_E, x1 + 1, x2, x3);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_E, x1+1, x2, x3)=ftemp-mfcbb;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)=fBB-c6o1*c2o27*(-vx);
	//					//}
	//					//if ((phi[DIR_M00] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(E, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(W, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p-1, x2, x3);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p - 1, x2, x3) = ftemp - mfabb;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = fBB - c6o1 * c2o27 * ( vx);
	//					//}
	//					//if ((phi[DIR_0P0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2+1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(S, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(N, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2+1, x3);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2+1, x3) = ftemp - mfbcb;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = fBB - c6o1 * c2o27 * (-vy);
	//					//}
	//					//if ((phi[DIR_0M0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2 - 1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(N, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(S, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p-1, x3);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p - 1, x3) = ftemp - mfbab;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = fBB - c6o1 * c2o27 * ( vy);
	//					//}
	//					//if ((phi[DIR_00P] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2 , x3+1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(B, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(T, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3+1) = ftemp - mfbbc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = fBB - c6o1 * c2o27 * (-vz);
	//					//}
	//					//if ((phi[DIR_00M] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(T, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(B, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p-1) = ftemp - mfbba;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = fBB - c6o1 * c2o27 * ( vz);
	//					//}
	//					//
	//					//if ((phi[DIR_PP0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2+1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(SW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(NE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_NE, x1+1, x2+1, x3);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_NE, x1 + 1, x2 + 1, x3) = ftemp - mfccb;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = fBB - c6o1 * c1o54 * (-vx-vy);
	//					//}
	//					//if ((phi[DIR_MM0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 - 1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(NE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(SW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p-1, x2p-1, x3);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p - 1, x2p - 1, x3) = ftemp - mfaab;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = fBB - c6o1 * c1o54 * ( vx + vy);
	//					//}
	//					//if ((phi[DIR_MP0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 + 1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(SE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(NW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p-1, x2+1, x3);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_NW, x1p - 1, x2 + 1, x3) = ftemp - mfacb;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = fBB - c6o1 * c1o54 * ( vx - vy);
	//					//}
	//					//if ((phi[DIR_PM0] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2 - 1, x3);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(NW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(SE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1+1, x2p-1, x3);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1 + 1, x2p - 1, x3) = ftemp - mfcab;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = fBB - c6o1 * c1o54 * (-vx + vy);
	//					//}
	//					//if ((phi[DIR_P0P] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2 , x3+1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TE, x1+1, x2, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TE, x1 + 1, x2, x3 + 1) = ftemp - mfcbc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = fBB - c6o1 * c1o54 * (-vx - vz);
	//					//}
	//					//if ((phi[DIR_M0P] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p-1, x2, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TW, x1p - 1, x2, x3 + 1) = ftemp - mfabc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = fBB - c6o1 * c1o54 * ( vx - vz);
	//					//}
	//					//if ((phi[DIR_P0M] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1+1, x2, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1 + 1, x2, x3p - 1) = ftemp - mfcba;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = fBB - c6o1 * c1o54 * (-vx + vz);
	//					//}
	//					//if ((phi[DIR_M0M] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p-1, x2, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p - 1, x2, x3p - 1) = ftemp - mfaba;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = fBB - c6o1 * c1o54 * ( vx + vz);
	//					//}
	//					//if ((phi[DIR_0PP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2+1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BS, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TN, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2+1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2 + 1, x3 + 1) = ftemp - mfbcc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = fBB - c6o1 * c1o54 * (-vy - vz);
	//					//}
	//					//if ((phi[DIR_0MP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2 - 1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BN, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TS, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p-1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p - 1, x3 + 1) = ftemp - mfbac;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = fBB - c6o1 * c1o54 * ( vy - vz);
	//					//}
	//					//if ((phi[DIR_0PM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2 + 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TS, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BN, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2+1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2 + 1, x3p - 1) = ftemp - mfbca;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = fBB - c6o1 * c1o54 * (-vy + vz);
	//					//}
	//					//if ((phi[DIR_0MM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1, x2 - 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TN, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BS, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p-1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p - 1, x3p - 1) = ftemp - mfbaa;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = fBB - c6o1 * c1o54 * (-vy - vz);
	//					//}

	//					//if ((phi[DIR_PPP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1+1, x2 + 1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BSW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TNE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1+1, x2+1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TNE, x1 + 1, x2 + 1, x3 + 1) = ftemp - mfccc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = fBB - c6o1 * c1o216 * (-vx -vy - vz);
	//					//}
	//					//if ((phi[DIR_MPP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 + 1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BSE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TNW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p-1, x2+1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p - 1, x2 + 1, x3 + 1) = ftemp - mfacc;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = fBB - c6o1 * c1o216 * ( vx - vy - vz);
	//					//}
	//					//if ((phi[DIR_MMP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 - 1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BNE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TSW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p-1, x2p-1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p - 1, x2p - 1, x3 + 1) = ftemp - mfaac;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = fBB - c6o1 * c1o216 * (vx + vy - vz);
	//					//}
	//					//if ((phi[DIR_PMP] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2 - 1, x3 + 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(BNW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(TSE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1+1, x2p-1, x3+1);
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TSE, x1 + 1, x2p - 1, x3 + 1) = ftemp - mfcac;
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = fBB - c6o1 * c1o216 * (-vx + vy - vz);
	//					//}
	//					//if ((phi[DIR_PMM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2 - 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TNW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BSE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1+1, x2p-1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1+1, x2p-1, x3p-1) = ftemp - mfcaa;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = fBB - c6o1 * c1o216 * (-vx + vy + vz);
	//					//}
	//					//if ((phi[DIR_MMM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 - 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TNE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BSW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p-1, x2p-1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p - 1, x2p - 1, x3p - 1) = ftemp - mfaaa;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = fBB - c6o1 * c1o216 * ( vx + vy + vz);
	//					//}
	//					//if ((phi[DIR_MPM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 - 1, x2 + 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TSE, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BNW, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p-1, x2+1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p - 1, x2 + 1, x3p - 1) = ftemp - mfaca;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = fBB - c6o1 * c1o216 * (vx - vy + vz);
	//					//}
	//					//if ((phi[DIR_PPM] <= c1o2)) {
	//					//	real rhoG = (*rhoNode)(x1 + 1, x2 + 1, x3 - 1);
	//					//	real ftemp = D3Q27System::getCompFeqForDirection(TSW, rhoG, vx, vy, vz)+ D3Q27System::getCompFeqForDirection(BNE, rhoG, vx, vy, vz);
	//					//	real fBB = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1+1, x2+1, x3p-1);
	//					//	(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1 + 1, x2 + 1, x3p - 1) = ftemp - mfcca;
	//					//	(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = fBB - c6o1 * c1o216 * (-vx - vy + vz);
	//					//}




	//				}
	//			}
	//		}
	//	}
	//}
}







//////////////////////////////////////////////////////////////////////////

real MultiphaseScaleDistributionLBMKernel::gradX1_phi()
{
	using namespace D3Q27System;
	return 3.0* ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) + (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) + (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) + (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_PM0] - phi[DIR_MP0]) + (phi[DIR_PP0] - phi[DIR_MM0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_P00] - phi[DIR_M00]));
}

real MultiphaseScaleDistributionLBMKernel::gradX2_phi()
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PPM] - phi[DIR_MMP])- (phi[DIR_PMP] - phi[DIR_MPM])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_0PP] - phi[DIR_0MM]) + (phi[DIR_0PM] - phi[DIR_0MP])) + ((phi[DIR_PP0] - phi[DIR_MM0])- (phi[DIR_PM0] - phi[DIR_MP0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_0P0] - phi[DIR_0M0]));
}

real MultiphaseScaleDistributionLBMKernel::gradX3_phi()
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) - (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) - (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_0MP] - phi[DIR_0PM]) + (phi[DIR_0PP] - phi[DIR_0MM])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_00P] - phi[DIR_00M]));
}

real MultiphaseScaleDistributionLBMKernel::gradX1_rhoInv(real rhoL,real rhoDIV)
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0/(rhoL+rhoDIV*phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_P0P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0M])) + (1.0 / (rhoL + rhoDIV * phi[DIR_P0M]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0P]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PM0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MP0])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PP0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MM0]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_P00]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M00])));
}

real MultiphaseScaleDistributionLBMKernel::gradX2_rhoInv(real rhoL,real rhoDIV)
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0 / (rhoL + rhoDIV * phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_0PP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_0PM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PP0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MM0])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PM0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MP0]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_0P0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0M0])));
}

real MultiphaseScaleDistributionLBMKernel::gradX3_rhoInv(real rhoL, real rhoDIV)
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0 / (rhoL + rhoDIV * phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_P0P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0M])) - (1.0 / (rhoL + rhoDIV * phi[DIR_P0M]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0P]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_0MP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0PM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_0PP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MM]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_00P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_00M])));
}

real MultiphaseScaleDistributionLBMKernel::gradX1_phi2()
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) + (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) + (phi2[DIR_PPM] - phi2[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) + (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_PM0] - phi2[DIR_MP0]) + (phi2[DIR_PP0] - phi2[DIR_MM0])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_P00] - phi2[DIR_M00]));
}

real MultiphaseScaleDistributionLBMKernel::gradX2_phi2()
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PPM] - phi2[DIR_MMP]) - (phi2[DIR_PMP] - phi2[DIR_MPM])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_0PP] - phi2[DIR_0MM]) + (phi2[DIR_0PM] - phi2[DIR_0MP])) + ((phi2[DIR_PP0] - phi2[DIR_MM0]) - (phi2[DIR_PM0] - phi2[DIR_MP0])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_0P0] - phi2[DIR_0M0]));
}

real MultiphaseScaleDistributionLBMKernel::gradX3_phi2()
{
	using namespace D3Q27System;
	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) - (phi2[DIR_PPM] - phi2[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) - (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_0MP] - phi2[DIR_0PM]) + (phi2[DIR_0PP] - phi2[DIR_0MM])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_00P] - phi2[DIR_00M]));
}

real MultiphaseScaleDistributionLBMKernel::nabla2_phi()
{
	using namespace D3Q27System;
	real sum = 0.0;
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

	return 6.0 * sum;
}

real MultiphaseScaleDistributionLBMKernel::computeCurvature_phi()
{
    using namespace D3Q27System;
    using namespace UbMath;

    real phiX = gradX1_phi();
    real phiY = gradX2_phi();
    real phiZ = gradX3_phi();
    real phiXX =
        c4o9 * (phi[DIR_P00] - c2o1 * phi[DIR_000] + phi[DIR_M00]) + (c1o9 * (((phi[DIR_PP0] - c2o1 * phi[DIR_0P0] + phi[DIR_MP0]) + (phi[DIR_PM0] - c2o1 * phi[DIR_0M0] + phi[DIR_MM0])) + ((phi[DIR_P0P] - c2o1 * phi[DIR_00P] + phi[DIR_M0P]) + (phi[DIR_P0M] - c2o1 * phi[DIR_00M] + phi[DIR_M0M]))) +
                                                                      c1o36 * (((phi[DIR_PPP] - c2o1 * phi[DIR_0PP] + phi[DIR_MPP]) + (phi[DIR_PMP] - c2o1 * phi[DIR_0MP] + phi[DIR_MMP])) + ((phi[DIR_PPM] - c2o1 * phi[DIR_0PM] + phi[DIR_MPM]) + (phi[DIR_PMM] - c2o1 * phi[DIR_0MM] + phi[DIR_MMM]))));
    real phiYY =
        c4o9 * (phi[DIR_0P0] - c2o1 * phi[DIR_000] + phi[DIR_0M0]) + (c1o9 * (((phi[DIR_PP0] - c2o1 * phi[DIR_P00] + phi[DIR_PM0]) + (phi[DIR_MP0] - c2o1 * phi[DIR_M00] + phi[DIR_MM0])) + ((phi[DIR_0PP] - c2o1 * phi[DIR_00P] + phi[DIR_0MP]) + (phi[DIR_0PM] - c2o1 * phi[DIR_00M] + phi[DIR_0MM]))) +
                                                                      c1o36 * (((phi[DIR_PPP] - c2o1 * phi[DIR_P0P] + phi[DIR_PMP]) + (phi[DIR_MPM] - c2o1 * phi[DIR_M0M] + phi[DIR_MMM])) + ((phi[DIR_MPP] - c2o1 * phi[DIR_M0P] + phi[DIR_MMP]) + (phi[DIR_PPM] - c2o1 * phi[DIR_P0M] + phi[DIR_PMM]))));
    real phiZZ =
        c4o9 * (phi[DIR_00P] - c2o1 * phi[DIR_000] + phi[DIR_00M]) + (c1o9 * (((phi[DIR_M0P] - c2o1 * phi[DIR_M00] + phi[DIR_M0M]) + (phi[DIR_P0P] - c2o1 * phi[DIR_P00] + phi[DIR_P0M])) + ((phi[DIR_0MP] - c2o1 * phi[DIR_0M0] + phi[DIR_0MM]) + (phi[DIR_0PP] - c2o1 * phi[DIR_0P0] + phi[DIR_0PM]))) +
                                                                      c1o36 * (((phi[DIR_MPP] - c2o1 * phi[DIR_MP0] + phi[DIR_MPM]) + (phi[DIR_PMP] - c2o1 * phi[DIR_PM0] + phi[DIR_PMM])) + ((phi[DIR_MMP] - c2o1 * phi[DIR_MM0] + phi[DIR_MMM]) + (phi[DIR_PPP] - c2o1 * phi[DIR_PP0] + phi[DIR_PPM]))));
    real phiXY = c1o4 * (c2o3 * (phi[DIR_MM0] - phi[DIR_PM0] + phi[DIR_PP0] - phi[DIR_MP0]) + c1o6 * ((phi[DIR_MMP] - phi[DIR_PMP] + phi[DIR_PPP] - phi[DIR_MPP]) + (phi[DIR_MMM] - phi[DIR_PMM] + phi[DIR_PPM] - phi[DIR_MPM])));
    real phiXZ = c1o4 * (c2o3 * (phi[DIR_M0M] - phi[DIR_P0M] + phi[DIR_P0P] - phi[DIR_M0P]) + c1o6 * ((phi[DIR_MPM] - phi[DIR_PPM] + phi[DIR_PPP] - phi[DIR_MPP]) + (phi[DIR_MMM] - phi[DIR_PMM] + phi[DIR_PMP] - phi[DIR_MMP])));
    real phiYZ = c1o4 * (c2o3 * (phi[DIR_0MM] - phi[DIR_0MP] + phi[DIR_0PP] - phi[DIR_0PM]) + c1o6 * ((phi[DIR_MMM] - phi[DIR_MMP] + phi[DIR_MPP] - phi[DIR_MPM]) + (phi[DIR_PMM] - phi[DIR_PMP] + phi[DIR_PPP] - phi[DIR_PPM])));

    // non isotropic FD (to be improved):
    // real phiX = (phi[DIR_P00] - phi[DIR_M00]) * c1o2; //gradX1_phi();
    // real phiY = (phi[DIR_0P0] - phi[DIR_0M0]) * c1o2; //gradX2_phi();
    // real phiZ = (phi[DIR_00P] - phi[DIR_00M]) * c1o2; //gradX3_phi();

    // real phiXX = phi[DIR_P00] - c2o1 * phi[DIR_000] + phi[DIR_M00];
    // real phiYY = phi[DIR_0P0] - c2o1 * phi[DIR_000] + phi[DIR_0M0];
    // real phiZZ =( phi[DIR_00P] - c2o1 * phi[DIR_000] + phi[DIR_00M]);
    // real phiXY = c1o4 * (phi[DIR_MM0] - phi[DIR_PM0] + phi[DIR_PP0] - phi[DIR_MP0]);
    // real phiXZ = c1o4 * (phi[DIR_M0M] - phi[DIR_P0M] + phi[DIR_P0P] - phi[DIR_M0P]);
    // real phiYZ = c1o4 * (phi[DIR_0MM] - phi[DIR_0MP] + phi[DIR_0PP] - phi[DIR_0PM]);
    // real back= (c2o1 * (phiX * phiY * phiXY + phiX * phiZ * phiXZ + phiY * phiZ * phiYZ) - phiXX * (phiY * phiY + phiZ * phiZ) - phiYY * (phiX * phiX + phiZ * phiZ) - phiZZ * (phiX * phiX + phiY * phiY)) / (c2o1 * pow(phiX * phiX + phiY * phiY + phiZ * phiZ, c3o2));

	return (c2o1 * (phiX * phiY * phiXY + phiX * phiZ * phiXZ + phiY * phiZ * phiYZ) - phiXX * (phiY * phiY + phiZ * phiZ) - phiYY * (phiX * phiX + phiZ * phiZ) - phiZZ * (phiX * phiX + phiY * phiY)) / (c2o1 * pow(phiX * phiX + phiY * phiY + phiZ * phiZ, c3o2)+1e-200);
}
void MultiphaseScaleDistributionLBMKernel::computePhasefield()
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

					h[DIR_P00]   = (*this->localDistributionsH1)(D3Q27System::ET_E, x1, x2, x3);
					h[DIR_0P0]   = (*this->localDistributionsH1)(D3Q27System::ET_N, x1, x2, x3);
					h[DIR_00P]   = (*this->localDistributionsH1)(D3Q27System::ET_T, x1, x2, x3);
					h[DIR_PP0]  = (*this->localDistributionsH1)(D3Q27System::ET_NE, x1, x2, x3);
					h[DIR_MP0]  = (*this->localDistributionsH1)(D3Q27System::ET_NW, x1p, x2, x3);
					h[DIR_P0P]  = (*this->localDistributionsH1)(D3Q27System::ET_TE, x1, x2, x3);
					h[DIR_M0P]  = (*this->localDistributionsH1)(D3Q27System::ET_TW, x1p, x2, x3);
					h[DIR_0PP]  = (*this->localDistributionsH1)(D3Q27System::ET_TN, x1, x2, x3);
					h[DIR_0MP]  = (*this->localDistributionsH1)(D3Q27System::ET_TS, x1, x2p, x3);
					h[DIR_PPP] = (*this->localDistributionsH1)(D3Q27System::ET_TNE, x1, x2, x3);
					h[DIR_MPP] = (*this->localDistributionsH1)(D3Q27System::ET_TNW, x1p, x2, x3);
					h[DIR_PMP] = (*this->localDistributionsH1)(D3Q27System::ET_TSE, x1, x2p, x3);
					h[DIR_MMP] = (*this->localDistributionsH1)(D3Q27System::ET_TSW, x1p, x2p, x3);

					h[DIR_M00]   = (*this->nonLocalDistributionsH1)(D3Q27System::ET_W, x1p, x2, x3);
					h[DIR_0M0]   = (*this->nonLocalDistributionsH1)(D3Q27System::ET_S, x1, x2p, x3);
					h[DIR_00M]   = (*this->nonLocalDistributionsH1)(D3Q27System::ET_B, x1, x2, x3p);
					h[DIR_MM0]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SW, x1p, x2p, x3);
					h[DIR_PM0]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_SE, x1, x2p, x3);
					h[DIR_M0M]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BW, x1p, x2, x3p);
					h[DIR_P0M]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BE, x1, x2, x3p);
					h[DIR_0MM]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BS, x1, x2p, x3p);
					h[DIR_0PM]  = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BN, x1, x2, x3p);
					h[DIR_MMM] = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					h[DIR_PMM] = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BSE, x1, x2p, x3p);
					h[DIR_MPM] = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNW, x1p, x2, x3p);
					h[DIR_PPM] = (*this->nonLocalDistributionsH1)(D3Q27System::ET_BNE, x1, x2, x3p);

					h[DIR_000] = (*this->zeroDistributionsH1)(x1, x2, x3);
				}
			}
		}
	}
}

void MultiphaseScaleDistributionLBMKernel::findNeighbors(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
	int x3)
{
	using namespace D3Q27System;

	SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

	phi[DIR_000] = (*ph)(x1, x2, x3);


	for (int k = FSTARTDIR; k <= FENDDIR; k++) {

		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
			phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
		} else {
			//phi[k] = (*ph)(x1 , x2, x3 );// neutral wetting
			phi[k] = 0.0;//unwetting
		}
	}
}

void MultiphaseScaleDistributionLBMKernel::findNeighbors2(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
	int x3)
{
	using namespace D3Q27System;

	SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

	phi2[DIR_000] = (*ph)(x1, x2, x3);


	for (int k = FSTARTDIR; k <= FENDDIR; k++) {

		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
			phi2[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
		}
		else {
			phi2[k] = 0.0;
            // phi2[k] = (*ph)(x1 , x2, x3 );// neutral wetting
		}
	}
}

void MultiphaseScaleDistributionLBMKernel::swapDistributions()
{
	LBMKernel::swapDistributions();
	dataSet->getHdistributions()->swap();
	dataSet->getH2distributions()->swap();
}

void MultiphaseScaleDistributionLBMKernel::initForcing()
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

	muForcingX1.DefineVar("rho",&muRho); 
	muForcingX2.DefineVar("rho",&muRho); 
	muForcingX3.DefineVar("rho",&muRho); 

}
