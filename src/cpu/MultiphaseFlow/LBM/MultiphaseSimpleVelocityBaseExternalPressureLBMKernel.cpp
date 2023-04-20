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
//! \file MultiphaseSimpleVelocityBaseExternalPressureLBMKernel.cpp
//! \ingroup LBMKernel
//! \author M. Geier, K. Kutscher, Hesameddin Safari
//=======================================================================================

#include "MultiphaseSimpleVelocityBaseExternalPressureLBMKernel.h"
#include "BCArray3D.h"
#include "Block3D.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
#include <cmath>
#include <iostream>
#include <string>
#include "basics/constants/NumericConstants.h"
//#include <basics/utilities/UbMath.h>

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::MultiphaseSimpleVelocityBaseExternalPressureLBMKernel() { this->compressible = false; }
//////////////////////////////////////////////////////////////////////////
void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::initDataSet()
{
	SPtr<DistributionArray3D> f(new D3Q27EsoTwist3DSplittedVector( nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9));
	SPtr<DistributionArray3D> h(new D3Q27EsoTwist3DSplittedVector( nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9)); // For phase-field
	SPtr<DistributionArray3D> h2(new D3Q27EsoTwist3DSplittedVector(nx[0] + 4, nx[1] + 4, nx[2] + 4, -999.9));
	SPtr<PhaseFieldArray3D> divU1(new PhaseFieldArray3D(            nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr pressure(new  CbArray3D<real, IndexerX3X2X1>(    nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	pressureOld = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
	p1Old = CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr(new  CbArray3D<real, IndexerX3X2X1>(nx[0] + 4, nx[1] + 4, nx[2] + 4, 0.0));
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
SPtr<LBMKernel> MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::clone()
{
	SPtr<LBMKernel> kernel(new MultiphaseSimpleVelocityBaseExternalPressureLBMKernel());
	kernel->setNX(nx);
	dynamicPointerCast<MultiphaseSimpleVelocityBaseExternalPressureLBMKernel>(kernel)->initDataSet();
	kernel->setCollisionFactorMultiphase(this->collFactorL, this->collFactorG);
	kernel->setDensityRatio(this->densityRatio);
	kernel->setMultiphaseModelParameters(this->beta, this->kappa);
	kernel->setContactAngle(this->contactAngle);
	kernel->setPhiL(this->phiL);
	kernel->setPhiH(this->phiH);
	kernel->setPhaseFieldRelaxation(this->tauH);
	kernel->setMobility(this->mob);
	kernel->setInterfaceWidth(this->interfaceWidth);

	kernel->setBCSet(bcSet->clone(kernel));
	kernel->setWithForcing(withForcing);
	kernel->setForcingX1(muForcingX1);
	kernel->setForcingX2(muForcingX2);
	kernel->setForcingX3(muForcingX3);
	kernel->setIndex(ix1, ix2, ix3);
	kernel->setDeltaT(deltaT);
	kernel->setGhostLayerWidth(2);
	dynamicPointerCast<MultiphaseSimpleVelocityBaseExternalPressureLBMKernel>(kernel)->initForcing();

	return kernel;
}
//////////////////////////////////////////////////////////////////////////
void  MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::forwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho) {
	using namespace vf::basics::constant;
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
void  MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::backwardInverseChimeraWithKincompressible(real& mfa, real& mfb, real& mfc, real vv, real v2, real Kinverse, real K, real oneMinusRho) {
	using namespace vf::basics::constant;
	real m0 = (((mfc - mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 - vv) * c1o2) * K;
	real m1 = (((mfa - mfc) - c2o1 * mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (-v2)) * K;
	mfc = (((mfc + mfb) * c1o2 + mfb * vv) * Kinverse + (mfa * Kinverse + oneMinusRho) * (v2 + vv) * c1o2) * K;
	mfa = m0;
	mfb = m1;
}


////////////////////////////////////////////////////////////////////////////////
void  MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::forwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2) {
	using namespace vf::basics::constant;
	real m1 = (mfa + mfc) + mfb;
	real m2 = mfc - mfa;
	mfc = (mfc + mfa) + (v2 * m1 - c2o1 * vv * m2);
	mfb = m2 - vv * m1;
	mfa = m1;
}


void  MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::backwardChimera(real& mfa, real& mfb, real& mfc, real vv, real v2) {
	using namespace vf::basics::constant;
	real ma = (mfc + mfa * (v2 - vv)) * c1o2 + mfb * (vv - c1o2);
	real mb = ((mfa - mfc) - mfa * v2) - c2o1 * mfb * vv;
	mfc = (mfc + mfa * (v2 + vv)) * c1o2 + mfb * (vv + c1o2);
	mfb = mb;
	mfa = ma;
}


void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::calculate(int step)
{
	using namespace D3Q27System;
	using namespace vf::basics::constant;
	using namespace vf::lbm::dir;

	forcingX1 = 0.0;
	forcingX2 = 0.0;
	forcingX3 = 0.0;

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

	for (int x3 = minX3-ghostLayerWidth; x3 < maxX3+ghostLayerWidth; x3++) {
		for (int x2 = minX2-ghostLayerWidth; x2 < maxX2+ghostLayerWidth; x2++) {
			for (int x1 = minX1-ghostLayerWidth; x1 < maxX1+ghostLayerWidth; x1++) {
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
					(*phaseField)(x1, x2, x3) = (((mfaaa + mfccc) + (mfaca + mfcac)) + ((mfaac + mfcca)  + (mfcaa + mfacc))  ) +
						(((mfaab + mfacb) + (mfcab + mfccb)) + ((mfaba + mfabc) + (mfcba + mfcbc)) +
							((mfbaa + mfbac) + (mfbca + mfbcc))) + ((mfabb + mfcbb) +
								(mfbab + mfbcb) + (mfbba + mfbbc)) + mfbbb;
					if ((*phaseField)(x1, x2, x3) > 1 ) {
						(*phaseField)(x1, x2, x3) = c1o1;
					}

					if ((*phaseField)(x1, x2, x3) < 0) {
						(*phaseField)(x1, x2, x3) = 0;
					}
					////// read F-distributions for velocity formalism

					mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
					mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
					mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
					mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
					mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
					mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
					mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
					mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
					mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
					mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
					mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
					mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
					mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
					mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
					mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
					mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
					mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
					mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
					mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
					mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
					mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
					mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
					mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
					mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
					mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

					mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);

					//LBMReal rhoH = 1.0;
					//LBMReal rhoL = 1.0 / densityRatio;

					real rhoH = 1.0*densityRatio;
					real rhoL = 1.0;

					real rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

					real drho = (((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc))   )
						+ (((mfaab + mfccb) + (mfacb + mfcab) ) + ((mfaba + mfcbc) + (mfabc + mfcba) ) + ((mfbaa + mfbcc) + (mfbac + mfbca) )))
						+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
				
					omegaDRho = 2.0;// 1.5;
					drho *= omegaDRho;
					real keepDrho = drho;
					drho = ((*p1Old)(x1, x2, x3) + drho) * c1o2;
				//	drho = ((*p1Old)(x1, x2, x3)*c2o3 + drho*c1o3) ;
					(*p1Old)(x1, x2, x3) = keepDrho;
					
					//LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH); //Incompressible
///Density correction
					//LBMReal dX1_phi = gradX1_phi();
					//LBMReal dX2_phi = gradX2_phi();
					//LBMReal dX3_phi = gradX3_phi();
					//LBMReal vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
					//	(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
					//	(mfcbb - mfabb)) ;
					//LBMReal vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
					//	(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
					//	(mfbcb - mfbab)) ;
					//LBMReal vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
					//	(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
					//	(mfbbc - mfbba)) ;
					//LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH) + (one - (*phaseField)(x1, x2, x3)) * three * (*pressure)(x1, x2, x3); //explicit Compressible
					//(*pressureOld)(x1, x2, x3) = (((*pressure)(x1, x2, x3) + rho * c1o3 * drho-(rhoH-rhoL)*(vvx*dX1_phi+vvy*dX2_phi+vvz*dX3_phi)*c1o3)) / (one - (one - (*phaseField)(x1, x2, x3)) * drho);
					
					//(*pressureOld)(x1, x2, x3) = ((*pressure)(x1, x2, x3) - c1o3 * drho * ((*phaseField)(x1, x2, x3) * (rhoH - rhoL) + rhoL)) / (c1 - ((*phaseField)(x1, x2, x3) - c1) * drho);
					//LBMReal rho=rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH) + (one - (*phaseField)(x1, x2, x3)) * three * (*pressureOld)(x1, x2, x3);
					//LBMReal tempDrho = drho;
					//drho = (drho*0.9 + (*pressureOld)(x1, x2, x3)*0.1) ;
					//(*pressureOld)(x1, x2, x3) = tempDrho;

					//Mathematica

					//LBMReal rho = ((*pressure)(x1, x2, x3) - (*phaseField)(x1, x2, x3) * (*pressure)(x1, x2, x3) + c1o3 * (rhoH + ((*phaseField)(x1, x2, x3) - phiH) * rhoToPhi)) / (c1o3 + c1o3 * drho * (-1 + (*phaseField)(x1, x2, x3)));
					(*pressureOld)(x1, x2, x3) = ((*pressure)(x1, x2, x3) + c1o3 * drho * (rhoH + ((*phaseField)(x1, x2, x3) - phiH) * rhoToPhi)) / (1 + drho * (-1 + (*phaseField)(x1, x2, x3)));
/////Full Filter
					//LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH)+(one- (*phaseField)(x1, x2, x3))*three* (*pressure)(x1, x2, x3); //explicit Compressible
					//(*pressureOld)(x1, x2, x3) = (((*pressure)(x1, x2, x3) + rho * c1o3 * drho)) / (one - (one - (*phaseField)(x1, x2, x3)) * drho);
//// reduced Filter
					//LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) - phiH) + (one - (*phaseField)(x1, x2, x3)) * three * (*pressureOld)(x1, x2, x3); //explicit Compressible
					//(*pressure)(x1, x2, x3) = (((*pressureOld)(x1, x2, x3) + rho * c1o3 * drho)) / (one - (one - (*phaseField)(x1, x2, x3)) * drho);

					//rho = (rho)/(one- (one - (*phaseField)(x1, x2, x3)) * drho); // now implicit Compressible
					
					//(*pressure)(x1, x2, x3) = (((*phaseField)(x1, x2, x3)) + ((*phaseField2)(x1, x2, x3)) - c1) * c1o3;
					////!!!!!! relplace by pointer swap!
					//(*pressureOld)(x1, x2, x3) = (*pressure)(x1, x2, x3);
				}
			}
		}
	}

	real collFactorM;

	////Periodic Filter
	//for (int x3 = minX3-1; x3 <= maxX3; x3++) {
	//	for (int x2 = minX2-1; x2 <= maxX2; x2++) {
	//		for (int x1 = minX1-1; x1 <= maxX1; x1++) {
	//			if (!bcArray->isSolid(x1, x2, x3) && !bcArray->isUndefined(x1, x2, x3)) {

	//				LBMReal sum = 0.;

	//				///Version for boundaries
	//				for (int xx = -1; xx <= 1; xx++) {
	//					//int xxx = (xx+x1 <= maxX1) ? ((xx + x1 > 0) ? xx + x1 : maxX1) : 0;
	//					int xxx = xx + x1;

	//					for (int yy = -1; yy <= 1; yy++) {
	//						//int yyy = (yy+x2 <= maxX2) ?( (yy + x2 > 0) ? yy + x2 : maxX2) : 0;
	//						int yyy = yy + x2;

	//						for (int zz = -1; zz <= 1; zz++) {
	//							//int zzz = (zz+x3 <= maxX3) ? zzz = ((zz + x3 > 0) ? zz + x3 : maxX3 ): 0;
	//							int zzz = zz + x3;

	//							if (!bcArray->isSolid(xxx, yyy, zzz) && !bcArray->isUndefined(xxx, yyy, zzz)) {
	//								sum+= 64.0/(216.0*(c1+c3*abs(xx))* (c1 + c3 * abs(yy))* (c1 + c3 * abs(zz)))*(*pressureOld)(xxx, yyy, zzz);
	//							}
	//							else{ sum+= 64.0 / (216.0 * (c1 + c3 * abs(xx)) * (c1 + c3 * abs(yy)) * (c1 + c3 * abs(zz))) * (*pressureOld)(x1, x2, x3);
	//							}


	//						}
	//					}
	//				}
	//				(*pressure)(x1, x2, x3) = sum;
	//			}
	//		}
	//	}
	//}

	////!filter

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

					real mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3);
					real mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3);
					real mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3);
					real mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3);
					real mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3);
					real mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3);
					real mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3);
					real mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3);
					real mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3);
					real mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3);
					real mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3);
					real mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3);
					real mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3);
					real mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3);
					real mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3);
					real mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p);
					real mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3);
					real mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3);
					real mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p);
					real mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p);
					real mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p);
					real mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p);
					real mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					real mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p);
					real mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p);
					real mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p);

					real mfbbb = (*this->zeroDistributionsF)(x1, x2, x3);


					real mfhcbb = (*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3);
					real mfhbcb = (*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3);
					real mfhbbc = (*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3);
					real mfhccb = (*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3);
					real mfhacb = (*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3);
					real mfhcbc = (*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3);
					real mfhabc = (*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3);
					real mfhbcc = (*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3);
					real mfhbac = (*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3);
					real mfhccc = (*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3);
					real mfhacc = (*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3);
					real mfhcac = (*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3);
					real mfhaac = (*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3);
					real mfhabb = (*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3);
					real mfhbab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3);
					real mfhbba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p);
					real mfhaab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3);
					real mfhcab = (*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3);
					real mfhaba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p);
					real mfhcba = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p);
					real mfhbaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p);
					real mfhbca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p);
					real mfhaaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p);
					real mfhcaa = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p);
					real mfhaca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p);
					real mfhcca = (*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p);

					real mfhbbb = (*this->zeroDistributionsH2)(x1, x2, x3);

					//LBMReal rhoH = 1.0;
					//LBMReal rhoL = 1.0 / densityRatio;

					real rhoH = 1.0;
					real rhoL = 1.0/ densityRatio;

					real rhoToPhi = (rhoH - rhoL) / (phiH - phiL);

					real dX1_phi = gradX1_phi();
					real dX2_phi = gradX2_phi();
					real dX3_phi = gradX3_phi();

					real denom = sqrt(dX1_phi * dX1_phi + dX2_phi * dX2_phi + dX3_phi * dX3_phi) + 1e-9+1e-3;
					// 01.09.2022: unclear what value we have to add to the normal: lager values better cut of in gas phase?
					real normX1 = dX1_phi / denom;
					real normX2 = dX2_phi / denom;
					real normX3 = dX3_phi / denom;



					collFactorM = collFactorL + (collFactorL - collFactorG) * (phi[DIR_000] - phiH) / (phiH - phiL);


					real mu = 2 * beta * phi[DIR_000] * (phi[DIR_000] - 1) * (2 * phi[DIR_000] - 1) - kappa * nabla2_phi();

					//----------- Calculating Macroscopic Values -------------
					real rho = rhoH + rhoToPhi * (phi[DIR_000] - phiH); //Incompressible

																		///scaled phase field
					//LBMReal rho = rhoH + rhoToPhi * ((*phaseField)(x1, x2, x3) * (*phaseField)(x1, x2, x3) / ((*phaseField)(x1, x2, x3) * (*phaseField)(x1, x2, x3) + (c1 - (*phaseField)(x1, x2, x3)) * (c1 - (*phaseField)(x1, x2, x3))) - phiH);
					///!scaled phase field
					
					//LBMReal rho = rhoH + rhoToPhi * (phi[DIR_000] - phiH)+(one-phi[DIR_000])* (*pressure)(x1, x2, x3)*three; //compressible
					//LBMReal rho = rhoL + (rhoH - rhoL) * phi[DIR_000] + (one - phi[DIR_000]) * (*pressure)(x1, x2, x3) * three; //compressible

					real m0, m1, m2;
					real rhoRef= c1o1;

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
					real pressureHere = (*pressureOld)(x1, x2, x3);
					//LBMReal pressureHere = (*pressure)(x1, x2, x3);

					real arrayP[3][3][3] = { {{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}},
												{{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}},
												{ {pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere},{pressureHere,pressureHere,pressureHere}} };
					//LBMReal LaplaceP = 0.0;
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
					
					//LBMReal sum = 0.0;

					for (int dir1 = -1; dir1 <= 1; dir1++) {
						for (int dir2 = -1; dir2 <= 1; dir2++) {
							for (int dir3 = -1; dir3 <= 1; dir3++){
								int xxx = x1 + dir1;
								int yyy = x2 + dir2;
								int zzz = x3 + dir3;
								if (!bcArray->isSolid(xxx, yyy, zzz) && !bcArray->isUndefined(xxx, yyy, zzz)) arrayP[dir1 + 1][dir2 + 1][dir3 + 1] = (*pressureOld)(xxx, yyy, zzz);
								//if (!bcArray->isSolid(xxx, yyy, zzz) && !bcArray->isUndefined(xxx, yyy, zzz)) arrayP[dir1 + 1][dir2 + 1][dir3 + 1] = (*pressure)(xxx, yyy, zzz);
							//	sum += 64.0 / (216.0 * (c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)) * (c1 + c3 * abs(dir3))) * arrayP[dir1 + 1][dir2 + 1][dir3 + 1];
							}
						}
					}
//					(*pressure)(x1, x2, x3) = sum;// *0.1 + (1.0 - 0.1) * (*pressureOld)(x1, x2, x3);


					(*pressure)(x1, x2, x3) = (((((arrayP[0][0][0] + arrayP[2][2][2]) + (arrayP[0][2][0] + arrayP[2][0][2])) + ((arrayP[2][0][0] + arrayP[0][2][2]) + (arrayP[2][2][0] + arrayP[0][0][2]))) * c1o216
						+ (((arrayP[0][0][1] + arrayP[2][2][1]) + (arrayP[0][1][0] + arrayP[2][1][2])) + ((arrayP[1][0][0] + arrayP[1][2][2]) + (arrayP[0][1][2] + arrayP[2][1][0])) + ((arrayP[1][0][2] + arrayP[1][2][0]) + (arrayP[0][2][1] + arrayP[2][0][1]))) * c1o54)
						+ ((arrayP[0][1][1] + arrayP[2][1][1]) + (arrayP[1][0][1] + arrayP[1][2][1]) + (arrayP[1][1][0] + arrayP[1][1][2])) * c2o27)
						+ arrayP[1][1][1] * c8o27;
					//LBMReal gradPx = 0.0;
					//LBMReal gradPy = 0.0;
					//LBMReal gradPz = 0.0;
					//for (int dir1 = -1; dir1 <= 1; dir1++) {
					//	for (int dir2 = -1; dir2 <= 1; dir2++) {
					//		gradPx -= arrayP[0][dir1+1][dir2+1] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		gradPx += arrayP[2][dir1+1][dir2+1] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));



					//		gradPy -= arrayP[dir1+1][0][dir2+1] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		gradPy += arrayP[dir1+1][2][dir2+1] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		

					//		gradPz -= arrayP[dir1+1][dir2+1][0] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		gradPz += arrayP[dir1+1][dir2+1][2] * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//	}
					//}

					//LBMReal gradPx = ((((arrayP[2][0][0] - arrayP[0][2][2]) + (arrayP[2][2][0] - arrayP[0][0][2])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[2][0][2] - arrayP[0][2][0]))) * c1o72
					//	+ (((arrayP[2][1][0] - arrayP[0][1][2]) + (arrayP[2][2][1] - arrayP[0][0][1])) + ((arrayP[2][0][1] - arrayP[0][2][1]) + (arrayP[2][1][2] - arrayP[0][1][0]))) * c1o18)
					//	+ (arrayP[2][1][1] - arrayP[0][1][1]) * c2o9;
					//LBMReal gradPy = ((((arrayP[0][2][0] - arrayP[2][0][2]) + (arrayP[2][2][0] - arrayP[0][0][2])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[0][2][2] - arrayP[2][0][0]))) * c1o72
					//	+ (((arrayP[1][2][0] - arrayP[1][0][2]) + (arrayP[2][2][1] - arrayP[0][0][1])) + ((arrayP[0][2][1] - arrayP[2][0][1]) + (arrayP[1][2][2] - arrayP[1][0][0]))) * c1o18)
					//	+ (arrayP[1][2][1] - arrayP[1][0][1]) * c2o9;
					//LBMReal gradPz = ((((arrayP[0][0][2] - arrayP[2][2][0]) + (arrayP[0][2][2] - arrayP[2][0][0])) + ((arrayP[2][2][2] - arrayP[0][0][0]) + (arrayP[2][0][2] - arrayP[0][2][0]))) * c1o72
					//	+ (((arrayP[0][1][2] - arrayP[2][1][0]) + (arrayP[1][2][2] - arrayP[1][0][0])) + ((arrayP[1][0][2] - arrayP[1][2][0]) + (arrayP[2][1][2] - arrayP[0][1][0]))) * c1o18)
					//	+ (arrayP[1][1][2] - arrayP[1][1][0]) * c2o9;

					//gradPx *=c1 - (*pressure)(x1, x2, x3)+pressureHere;
					//gradPy *=c1 - (*pressure)(x1, x2, x3) + pressureHere;
					//gradPz *=c1 - (*pressure)(x1, x2, x3) + pressureHere;

					////!Filter&Gradient merged
					//LBMReal gradPx = 0.0;
					//LBMReal gradPy = 0.0;
					//LBMReal gradPz = 0.0;
					//for (int dir1 = -1; dir1 <= 1; dir1++) {
					//	for (int dir2 = -1; dir2 <= 1; dir2++) {
					//		int yyy = x2 + dir1;
					//		int zzz = x3 + dir2;
					//		if (!bcArray->isSolid(x1-1, yyy, zzz) && !bcArray->isUndefined(x1-1, yyy, zzz)) {
					//			gradPx -= (*pressure)(x1 - 1, yyy, zzz) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPx -= (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(x1 + 1, yyy, zzz) && !bcArray->isUndefined(x1 + 1, yyy, zzz)) {
					//			gradPx += (*pressure)(x1 + 1, yyy, zzz) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPx += (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}

					//		int xxx = x1 + dir1;
					//		if (!bcArray->isSolid(xxx, x2-1, zzz) && !bcArray->isUndefined(xxx, x2-1, zzz)) {
					//			gradPy -= (*pressure)(xxx, x2-1, zzz) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPy -= (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(xxx, x2+1, zzz) && !bcArray->isUndefined(xxx, x2+1, zzz)) {
					//			gradPy += (*pressure)(xxx, x2+1, zzz) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPy += (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}

					//		yyy = x2 + dir2;
					//		if (!bcArray->isSolid(xxx, yyy, x3-1) && !bcArray->isUndefined(xxx, yyy, x3-1)) {
					//			gradPz -= (*pressure)(xxx, yyy, x3-1) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPz -= (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		if (!bcArray->isSolid(xxx, yyy, x3+1) && !bcArray->isUndefined(xxx, yyy, x3+1)) {
					//			gradPz += (*pressure)(xxx, yyy, x3+1) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}
					//		else {
					//			gradPz += (*pressure)(x1, x2, x3) * c2o9 / ((c1 + c3 * abs(dir1)) * (c1 + c3 * abs(dir2)));
					//		}

					//	}
					//}

					//Viscosity increase by phase field residuum
					//LBMReal errPhi = (((1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale)- denom);
					//LBMReal limVis = 0.01;// 0.0000001 * 10;//0.01;
					// collFactorM =collFactorM/(c1+limVis*(errPhi*errPhi)*collFactorM);
					// collFactorM = (collFactorM < 1.8) ? 1.8 : collFactorM;
					//errPhi = errPhi * errPhi* errPhi * errPhi * errPhi * errPhi;
					//collFactorM = collFactorM + (1.8 - collFactorM) * errPhi / (errPhi + limVis);

					//3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) - (phi2[DIR_PPM] - phi2[DIR_MMP])))
					//+WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) - (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_0MP] - phi2[DIR_0PM]) + (phi2[DIR_0PP] - phi2[DIR_0MM])))) +
					//+WEIGTH[DIR_0P0] * (phi2[DIR_00P] - phi2[DIR_00M]));

					muRho = rho;

					////external pressure
					//forcingX1 =/* muForcingX1.Eval()/rho */- gradPx/rho;
					//forcingX2 =/* muForcingX2.Eval()/rho */- gradPy/rho;
					//forcingX3 =/* muForcingX3.Eval()/rho */- gradPz/rho;

					///////////////////////////////////////////////

					//LBMReal pBefore = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
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

					real pStarStart = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
						+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
						+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb) * c1o3;

					/////////////////////
					//LBMReal vvxh = ((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhcaa - mfhacc) + (mfhcca - mfhaac))) +
					//	(((mfhcba - mfhabc) + (mfhcbc - mfhaba)) + ((mfhcab - mfhacb) + (mfhccb - mfhaab))) +
					//	(mfhcbb - mfhabb)) / rhoRef;
					//LBMReal vvyh = ((((mfhccc - mfhaaa) + (mfhaca - mfhcac)) + ((mfhacc - mfhcaa) + (mfhcca - mfhaac))) +
					//	(((mfhbca - mfhbac) + (mfhbcc - mfhbaa)) + ((mfhacb - mfhcab) + (mfhccb - mfhaab))) +
					//	(mfhbcb - mfhbab)) / rhoRef;
					//LBMReal vvzh = ((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhacc - mfhcaa) + (mfhaac - mfhcca))) +
					//	(((mfhbac - mfhbca) + (mfhbcc - mfhbaa)) + ((mfhabc - mfhcba) + (mfhcbc - mfhaba))) +
					//	(mfhbbc - mfhbba)) / rhoRef;

					//LBMReal deltaPP = 0*(vvxh * dX1_phi + vvyh * dX2_phi + vvzh * dX3_phi) * rhoToPhi / (rho);
					//mfhbcb += c1o6* c2o9  * deltaPP;
					//mfhbbc += c1o6* c2o9  * deltaPP;
					//mfhcbb += c1o6* c2o9  * deltaPP;
					//mfhccb += c1o6* c1o18 * deltaPP;
					//mfhacb += c1o6* c1o18 * deltaPP;
					//mfhcbc += c1o6* c1o18 * deltaPP;
					//mfhabc += c1o6* c1o18 * deltaPP;
					//mfhbcc += c1o6* c1o18 * deltaPP;
					//mfhbac += c1o6* c1o18 * deltaPP;
					//mfhccc += c1o6* c1o72 * deltaPP;
					//mfhacc += c1o6* c1o72 * deltaPP;
					//mfhcac += c1o6* c1o72 * deltaPP;
					//mfhaac += c1o6* c1o72 * deltaPP;
					//mfhabb += c1o6* c2o9  * deltaPP;
					//mfhbab += c1o6* c2o9  * deltaPP;
					//mfhbba += c1o6* c2o9  * deltaPP;
					//mfhaab += c1o6* c1o18 * deltaPP;
					//mfhcab += c1o6* c1o18 * deltaPP;
					//mfhaba += c1o6* c1o18 * deltaPP;
					//mfhcba += c1o6* c1o18 * deltaPP;
					//mfhbaa += c1o6* c1o18 * deltaPP;
					//mfhbca += c1o6* c1o18 * deltaPP;
					//mfhaaa += c1o6* c1o72 * deltaPP;
					//mfhcaa += c1o6* c1o72 * deltaPP;
					//mfhaca += c1o6* c1o72 * deltaPP;
					//mfhcca += c1o6* c1o72 * deltaPP;
					//mfhbbb += c1o6* c4 * c2o9 * deltaPP;

					//////////////////////

					/////Recovering the origin distributions
					//LBMReal mfStartcbb = mfcbb ;
					//LBMReal mfStartbcb = mfbcb ;
					//LBMReal mfStartbbc = mfbbc ;
					//LBMReal mfStartccb = mfccb ;
					//LBMReal mfStartacb = mfacb ;
					//LBMReal mfStartcbc = mfcbc ;
					//LBMReal mfStartabc = mfabc ;
					//LBMReal mfStartbcc = mfbcc ;
					//LBMReal mfStartbac = mfbac ;
					//LBMReal mfStartccc = mfccc ;
					//LBMReal mfStartacc = mfacc ;
					//LBMReal mfStartcac = mfcac ;
					//LBMReal mfStartaac = mfaac ;
					//LBMReal mfStartabb = mfabb ;
					//LBMReal mfStartbab = mfbab ;
					//LBMReal mfStartbba = mfbba ;
					//LBMReal mfStartaab = mfaab ;
					//LBMReal mfStartcab = mfcab ;
					//LBMReal mfStartaba = mfaba ;
					//LBMReal mfStartcba = mfcba ;
					//LBMReal mfStartbaa = mfbaa ;
					//LBMReal mfStartbca = mfbca ;
					//LBMReal mfStartaaa = mfaaa ;
					//LBMReal mfStartcaa = mfcaa ;
					//LBMReal mfStartaca = mfaca ;
					//LBMReal mfStartcca = mfcca ;
					//LBMReal mfStartbbb = mfbbb ;


						mfcbb += mfhcbb /rho;
						mfbcb += mfhbcb /rho;
						mfbbc += mfhbbc /rho;
						mfccb += mfhccb /rho;
						mfacb += mfhacb /rho;
						mfcbc += mfhcbc /rho;
						mfabc += mfhabc /rho;
						mfbcc += mfhbcc /rho;
						mfbac += mfhbac /rho;
						mfccc += mfhccc /rho;
						mfacc += mfhacc /rho;
						mfcac += mfhcac /rho;
						mfaac += mfhaac /rho;
						mfabb += mfhabb /rho;
						mfbab += mfhbab /rho;
						mfbba += mfhbba /rho;
						mfaab += mfhaab /rho;
						mfcab += mfhcab /rho;
						mfaba += mfhaba /rho;
						mfcba += mfhcba /rho;
						mfbaa += mfhbaa /rho;
						mfbca += mfhbca /rho;
						mfaaa += mfhaaa /rho;
						mfcaa += mfhcaa /rho;
						mfaca += mfhaca /rho;
						mfcca += mfhcca /rho;
						mfbbb += mfhbbb /rho;



					//Abbas
					real pStar = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
						+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
						+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb) * c1o3;
					//22.09.22 not yet in balance, repaire here
					//LBMReal ppStar = ((((((mfhaaa + mfhccc) + (mfhaac + mfhcca)) + ((mfhcac + mfhaca) + (mfhcaa + mfhacc)))*c3
					//	+ (((mfhaab + mfhccb) + (mfhacb + mfhcab)) + ((mfhaba + mfhcbc) + (mfhabc + mfhcba)) + ((mfhbaa + mfhbcc) + (mfhbac + mfhbca))))*c2
					//	+ ((mfhabb + mfhcbb) + (mfhbab + mfhbcb) + (mfhbba + mfhbbc))) ) * c1o3/rho;
	
					//ppStar = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc))) * c3
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca)))) * c2
					//	+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc)))) * c1o3 ;

					//(*phaseFieldOld)(x1, x2, x3) = ((*phaseFieldOld)(x1, x2, x3) > 99.0) ? phi[DIR_000] : (*phaseFieldOld)(x1, x2, x3);
					//LBMReal dtPhi = phi[DIR_000] - (*phaseFieldOld)(x1, x2, x3);
					//LBMReal deltaP = -pStar * (c1 - rho / (rho + c1o2 * rhoToPhi * dtPhi));// -pStar * pStar * pStar * 1.0e-4 * rho * rho * rho;
					//LBMReal deltaP = pStar * (c1 - mfhbbb*rho) * c1o2;//Explicit
					//LBMReal deltaP = pStar * (c1 - mfhbbb * rho) / (c1 + mfhbbb * rho);//Semi-Implicit
					//(*phaseFieldOld)(x1, x2, x3) = phi[DIR_000];

					//mfabb += c2o9 *deltaP;
					//mfbab += c2o9 *deltaP;
					//mfbba += c2o9 *deltaP;
					//mfaab += c1o18*deltaP;
					//mfcab += c1o18*deltaP;
					//mfaba += c1o18*deltaP;
					//mfcba += c1o18*deltaP;
					//mfbaa += c1o18*deltaP;
					//mfbca += c1o18*deltaP;
					//mfaaa += c1o72*deltaP;
					//mfcaa += c1o72*deltaP;
					//mfaca += c1o72*deltaP;
					//mfcca += c1o72*deltaP;
					//mfcbb += c2o9 *deltaP;
					//mfbcb += c2o9 *deltaP;
					//mfbbc += c2o9 *deltaP;
					//mfccb += c1o18*deltaP;
					//mfacb += c1o18*deltaP;
					//mfcbc += c1o18*deltaP;
					//mfabc += c1o18*deltaP;
					//mfbcc += c1o18*deltaP;
					//mfbac += c1o18*deltaP;
					//mfccc += c1o72*deltaP;
					//mfacc += c1o72*deltaP;
					//mfcac += c1o72*deltaP;
					//mfaac += c1o72*deltaP;

					//pStar = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb) * c1o3;




					//mfabb -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_P00])*rhoToPhi/rho;
					//mfbab -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_0P0])*rhoToPhi/rho;
					//mfbba -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_00P])*rhoToPhi/rho;
					//mfaab -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_PP0])*rhoToPhi/rho;
					//mfcab -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_MP0])*rhoToPhi/rho;
					//mfaba -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_P0P])*rhoToPhi/rho;
					//mfcba -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_M0P])*rhoToPhi/rho;
					//mfbaa -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_0PP])*rhoToPhi/rho;
					//mfbca -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_0MP])*rhoToPhi/rho;
					//mfaaa -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_PPP])*rhoToPhi/rho;
					//mfcaa -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_MPP])*rhoToPhi/rho;
					//mfaca -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_PMP])*rhoToPhi/rho;
					//mfcca -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_MMP])*rhoToPhi/rho;
					//mfcbb -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_M00])*rhoToPhi/rho;
					//mfbcb -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_0M0])*rhoToPhi/rho;
					//mfbbc -= c1o2 * c2o9 *pStar*(phi[DIR_000]-phi[DIR_00M])*rhoToPhi/rho;
					//mfccb -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_MM0])*rhoToPhi/rho;
					//mfacb -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_PM0])*rhoToPhi/rho;
					//mfcbc -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_M0M])*rhoToPhi/rho;
					//mfabc -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_P0M])*rhoToPhi/rho;
					//mfbcc -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_0MM])*rhoToPhi/rho;
					//mfbac -= c1o2 * c1o18*pStar*(phi[DIR_000]-phi[DIR_0PM])*rhoToPhi/rho;
					//mfccc -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_MMM])*rhoToPhi/rho;
					//mfacc -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_PMM])*rhoToPhi/rho;
					//mfcac -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_MPM])*rhoToPhi/rho;
					//mfaac -= c1o2 * c1o72*pStar*(phi[DIR_000]-phi[DIR_PPM])*rhoToPhi/rho;


					//forcingX1 =/* muForcingX1.Eval() / rho*/ - pStar * dX1_phi * rhoToPhi / rho;
					//forcingX2 =/* muForcingX2.Eval() / rho*/ - pStar * dX2_phi * rhoToPhi / rho;
					//forcingX3 =/* muForcingX3.Eval() / rho*/ - pStar * dX3_phi * rhoToPhi / rho;


					//mfabb += (-forcingX1) * c2o9;
					//mfbab += (-forcingX2) * c2o9;
					//mfbba += (-forcingX3) * c2o9;
					//mfaab += (-forcingX1 - forcingX2) * c1o16;
					//mfcab += (forcingX1 - forcingX2) * c1o16;
					//mfaba += (-forcingX1 - forcingX3) * c1o16;
					//mfcba += (forcingX1 - forcingX3) * c1o16;
					//mfbaa += (-forcingX2 - forcingX3) * c1o16;
					//mfbca += (forcingX2 - forcingX3) * c1o16;
					//mfaaa += (-forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfcaa += (forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfaca += (-forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcca += (forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcbb += (forcingX1)*c2o9;
					//mfbcb += (forcingX2)*c2o9;
					//mfbbc += (forcingX3)*c2o9;
					//mfccb += (forcingX1 + forcingX2) * c1o16;
					//mfacb += (-forcingX1 + forcingX2) * c1o16;
					//mfcbc += (forcingX1 + forcingX3) * c1o16;
					//mfabc += (-forcingX1 + forcingX3) * c1o16;
					//mfbcc += (forcingX2 + forcingX3) * c1o16;
					//mfbac += (-forcingX2 + forcingX3) * c1o16;
					//mfccc += (forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfacc += (-forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfcac += (forcingX1 - forcingX2 + forcingX3) * c1o72;
					//mfaac += (-forcingX1 - forcingX2 + forcingX3) * c1o72;

					//LBMReal saveForceX1 = forcingX1;
					//LBMReal saveForceX2 = forcingX2;
					//LBMReal saveForceX3 = forcingX3;

					 vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
						(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
						(mfcbb - mfabb)) / rhoRef;
					 vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
						(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
						(mfbcb - mfbab)) / rhoRef;
					 vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
						(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
						(mfbbc - mfbba)) / rhoRef;


					 //LBMReal dRhoInvX = -(((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhcaa - mfhacc) + (mfhcca - mfhaac))) +
						// (((mfhcba - mfhabc) + (mfhcbc - mfhaba)) + ((mfhcab - mfhacb) + (mfhccb - mfhaab))) +
						// (mfhcbb - mfhabb)));
					 //LBMReal dRhoInvY = -(((((mfhccc - mfhaaa) + (mfhaca - mfhcac)) + ((mfhacc - mfhcaa) + (mfhcca - mfhaac))) +
						// (((mfhbca - mfhbac) + (mfhbcc - mfhbaa)) + ((mfhacb - mfhcab) + (mfhccb - mfhaab))) +
						// (mfhbcb - mfhbab)));
					 //LBMReal dRhoInvZ = -(((((mfhccc - mfhaaa) + (mfhcac - mfhaca)) + ((mfhacc - mfhcaa) + (mfhaac - mfhcca))) +
						// (((mfhbac - mfhbca) + (mfhbcc - mfhbaa)) + ((mfhabc - mfhcba) + (mfhcbc - mfhaba))) +
						// (mfhbbc - mfhbba)));


					 forcingX1 = 0.0;
					 forcingX2 = 0.0;
					 forcingX3 = 0.0;
					//!Abbas
					//LBMReal dX1_rhoInv = gradX1_rhoInv(rhoL, rhoH - rhoL);
					//LBMReal dX2_rhoInv = gradX2_rhoInv(rhoL, rhoH - rhoL);
					//LBMReal dX3_rhoInv = gradX3_rhoInv(rhoL, rhoH - rhoL);
					//forcingX1 =/* muForcingX1.Eval() / rho*/ +pStar * dX1_rhoInv * rho;
					//forcingX2 =/* muForcingX2.Eval() / rho*/ +pStar * dX2_rhoInv * rho;
					//forcingX3 =/* muForcingX3.Eval() / rho*/ +pStar * dX3_rhoInv * rho;

					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho + pStar * dX1_rhoInv * rho) *c1o2;
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho + pStar * dX2_rhoInv * rho) *c1o2;
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho + pStar * dX3_rhoInv * rho) *c1o2;
					 //LBMReal FdX1_phi = normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					 //LBMReal FdX2_phi = normX2 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					 //LBMReal FdX3_phi = normX3 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;


					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho ) ;
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho ) ;
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho ) ;

					//forcingX1 = (pStar * dRhoInvX* rho *c3) ;
					//forcingX2 = (pStar * dRhoInvY* rho *c3) ;
					//forcingX3 = (pStar * dRhoInvZ* rho *c3) ;
					//if (phi[DIR_000] > 0.1 && phi[DIR_000] < 0.9) std::cout << phi[DIR_000] << " " << dX1_phi * rhoToPhi / rho << " " << dRhoInvX * rho *3<< std::endl;
					//LBMReal forcingX1ALTERNAT = ( pStar * dX1_rhoInv * rho) ;
					//LBMReal forcingX2ALTERNAT = ( pStar * dX2_rhoInv * rho) ;
					//LBMReal forcingX3ALTERNAT = ( pStar * dX3_rhoInv * rho) ;

					//forcingX1 = (fabs(vvx + c1o2 * forcingX1) < fabs(vvx + c1o2 * forcingX1ALTERNAT)) ? forcingX1 : forcingX1ALTERNAT;
					//forcingX2 = (fabs(vvy + c1o2 * forcingX2) < fabs(vvy + c1o2 * forcingX2ALTERNAT)) ? forcingX2 : forcingX2ALTERNAT;
					//forcingX3 = (fabs(vvz + c1o2 * forcingX3) < fabs(vvz + c1o2 * forcingX3ALTERNAT)) ? forcingX3 : forcingX3ALTERNAT;

					//	 forcingX1 = -pStar * rhoToPhi / rho * normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					//	 forcingX2 = -pStar * rhoToPhi / rho * normX2 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;
					//	 forcingX3 = -pStar * rhoToPhi / rho * normX3 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale;

					//forcingX1 = (-pStar * dX1_phi * rhoToPhi / rho *(c1- phi[DIR_000]) + pStar * dX1_rhoInv * rho*(phi[DIR_000]));
					//forcingX2 = (-pStar * dX2_phi * rhoToPhi / rho *(c1- phi[DIR_000]) + pStar * dX2_rhoInv * rho*(phi[DIR_000]));
					//forcingX3 = (-pStar * dX3_phi * rhoToPhi / rho *(c1- phi[DIR_000]) + pStar * dX3_rhoInv * rho*(phi[DIR_000]));
						 //if (phi[DIR_000] > 0.3 && phi[DIR_000] < 0.7)
						 //{
							// int test = 1;
							// std::cout << phi[DIR_000] <<" "<< dX1_phi <<" "<< normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale<<" "<< normX1 * (1.0 - phi[DIR_000]) * (phi[DIR_000]) * oneOverInterfaceScale/ dX1_phi<< std::endl;
						 //}



					 //LBMReal scaleGrad = c2 * phi[DIR_000] * (1.0 - phi[DIR_000]) / ((phi[DIR_000] * phi[DIR_000] + (1.0 - phi[DIR_000]) * (1.0 - phi[DIR_000])) * (phi[DIR_000] * phi[DIR_000] + (1.0 - phi[DIR_000]) * (1.0 - phi[DIR_000])));
					 //dX1_phi *= scaleGrad;
					 //dX2_phi *= scaleGrad;
					 //dX3_phi *= scaleGrad;

					 ///Experimental interface sharpening force 20.06.2022

					 //LBMReal scaleSharpener = 1.0;
					 //forcingX1 += scaleSharpener * (FdX1_phi - dX1_phi) * fabsf(FdX1_phi - dX1_phi)  / rho;
					 //forcingX2 += scaleSharpener * (FdX2_phi - dX2_phi) * fabsf(FdX2_phi - dX2_phi)  / rho;
					 //forcingX3 += scaleSharpener * (FdX3_phi - dX3_phi) * fabsf(FdX3_phi - dX3_phi)  / rho;
					///surface tension force
					forcingX1 += mu * dX1_phi/rho;
					forcingX2 += mu * dX2_phi/rho;
					forcingX3 += mu * dX3_phi/rho;

					//LBMReal forcingBIAS = 0.5;
					forcingX1 += muForcingX1.Eval() / rho;//*phi[DIR_000];
					forcingX2 += muForcingX2.Eval() / rho;// * phi[DIR_000];
					forcingX3 += muForcingX3.Eval() / rho;// * phi[DIR_000];

				//	//19.08.2022
					//vvx += vvxh / rho * c1o2;
					//vvy += vvyh / rho * c1o2;
					//vvz += vvzh / rho * c1o2;
				//	//


					vvx += (forcingX1) * deltaT * c1o2;
					vvy += (forcingX2) * deltaT * c1o2;
					vvz += (forcingX3) * deltaT * c1o2;

					//vvx += (forcingX1 + muForcingX1.Eval() / rho) * deltaT *  c1o2; // X
					//vvy += (forcingX2 + muForcingX2.Eval() / rho) * deltaT *  c1o2; // Y
					//vvz += (forcingX3 + muForcingX3.Eval() / rho) * deltaT *  c1o2; // Z



				//	vvx += (forcingX1 + muForcingX1.Eval() / rho) * deltaT * forcingBIAS; // X
				//	vvy += (forcingX2 + muForcingX2.Eval() / rho) * deltaT * forcingBIAS; // Y
				//	vvz += (forcingX3 + muForcingX3.Eval() / rho) * deltaT * forcingBIAS; // Z



					//Abbas
					//LBMReal M200 = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab)) + ((mfaba + mfcbc) + (mfabc + mfcba)) ))
					//	+ ((mfabb + mfcbb))) );
					//LBMReal M020 = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ (((mfaab + mfccb) + (mfacb + mfcab))  + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ( (mfbab + mfbcb) )) );
					//LBMReal M002 = ((((((mfaaa + mfccc) + (mfaac + mfcca)) + ((mfcac + mfaca) + (mfcaa + mfacc)))
					//	+ ( + ((mfaba + mfcbc) + (mfabc + mfcba)) + ((mfbaa + mfbcc) + (mfbac + mfbca))))
					//	+ ( (mfbba + mfbbc))));

					//LBMReal M110 = ((((((mfaaa + mfccc) + (-mfcac - mfaca)) + ((mfaac + mfcca) + (-mfcaa -mfacc)))
					//	+ (((mfaab + mfccb) + (-mfacb - mfcab))   ))
					//	) );
					//LBMReal M101 = ((((((mfaaa + mfccc) - (mfaac + mfcca)) + ((mfcac + mfaca) - (mfcaa + mfacc)))
					//	+ (((mfaba + mfcbc) + (-mfabc - mfcba))))
					//	));
					//LBMReal M011 = ((((((mfaaa + mfccc) - (mfaac + mfcca)) + ( (mfcaa + mfacc)- (mfcac + mfaca)))
					//	+ (((mfbaa + mfbcc) + (-mfbac - mfbca))))
					//	));
					real vvxI = vvx;
					real vvyI = vvy;
					real vvzI = vvz;

					//LBMReal collFactorStore=collFactorM;
					//LBMReal stress;
					//for(int iter=0;iter<5;iter++)
				 //{
					//	LBMReal OxxPyyPzz = 1.0;
					//	LBMReal mxxPyyPzz = (M200-vvxI*vvxI) + (M020-vvyI*vvyI) + (M002-vvzI*vvzI);
					//	//pStar = mxxPyyPzz * c1o3;
					//mxxPyyPzz -= c3 *pStar;

					//LBMReal mxxMyy = (M200-vvxI*vvxI) - (M020-vvyI*vvyI);
					//LBMReal mxxMzz = (M200-vvxI*vvxI) - (M002-vvzI*vvzI);
					//LBMReal mxy = M110 - vvxI * vvyI;
					//LBMReal mxz = M101 - vvxI * vvzI;
					//LBMReal myz = M011 - vvyI * vvzI;


					//mxxMyy *= c1 - collFactorM * c1o2;
					//mxxMzz *= c1 - collFactorM * c1o2;
					//mxy *= c1 - collFactorM * c1o2;
					//mxz *= c1 - collFactorM * c1o2;
					//myz *= c1 - collFactorM * c1o2;
					//mxxPyyPzz *= c1 - OxxPyyPzz * c1o2;
					////mxxPyyPzz = mxxPyyPzz*fabs(mxxPyyPzz)/(1.0e-6+fabs(mxxPyyPzz));
					////mxxPyyPzz += c3 * pStar;
					//LBMReal mxx = (mxxMyy + mxxMzz + mxxPyyPzz)*c1o3;
					//LBMReal myy = (-c2*mxxMyy + mxxMzz + mxxPyyPzz)*c1o3;
					//LBMReal mzz = (mxxMyy -c2* mxxMzz + mxxPyyPzz) * c1o3;
					//vvxI = vvx - (mxx * dX1_phi + mxy * dX2_phi + mxz * dX3_phi) * rhoToPhi / (rho);
					//vvyI = vvy - (mxy * dX1_phi + myy * dX2_phi + myz * dX3_phi) * rhoToPhi / (rho);
					//vvzI = vvz - (mxz * dX1_phi + myz * dX2_phi + mzz * dX3_phi) * rhoToPhi / (rho);


				////	vvzI = vvz + (mxz * dRhoInvX + myz * dRhoInvY + mzz * dRhoInvZ) *  (rho)*c3;
				////	vvxI = vvx + (mxx * dRhoInvX + mxy * dRhoInvY + mxz * dRhoInvZ) *  (rho)*c3;
				////	vvyI = vvy + (mxy * dRhoInvX + myy * dRhoInvY + myz * dRhoInvZ) *  (rho)*c3;


				//	//LBMReal dxux = -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz);
				//	//
				//	//LBMReal dyuy = dxux + collFactorM * c3o2 * mxxMyy;
				//	//LBMReal dzuz = dxux + collFactorM * c3o2 * mxxMzz;
				//	//LBMReal Dxy = -three * collFactorM * mxy;
				//	//LBMReal Dxz = -three * collFactorM * mxz;
				//	//LBMReal Dyz = -three * collFactorM * myz;
				//	////LBMReal stress = sqrt(sqrt((dyuy+dxux+dzuz)* (dyuy + dxux + dzuz))) * sqrt(forcingX1 * forcingX1 + forcingX2 * forcingX2 + forcingX3 * forcingX3);
				//	// stress = sqrt(dyuy * dyuy + dxux * dxux + dzuz*dzuz + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz)*sqrt(forcingX1*forcingX1+forcingX2*forcingX2+forcingX3*forcingX3);
				//	////collFactorM = collFactorStore + (1.75 - collFactorStore) * stress / (stress + 1.0e-8);
				//	//
				//	//LBMReal dX2_rho = (rhoToPhi)*dX2_phi;
				//	//LBMReal dX1_rho = (rhoToPhi)*dX1_phi;
				//	//LBMReal dX3_rho = (rhoToPhi)*dX3_phi;
				//	////vvxI= vvx+ c1o6 * (c1 / collFactorM - c1o2) * (2 * dxux * dX1_rho + Dxy * dX2_rho + Dxz * dX3_rho) / (rho);
				//	////vvyI= vvy+ c1o6 * (c1 / collFactorM - c1o2) * (Dxy * dX1_rho + 2 * dyuy * dX2_rho + Dyz * dX3_rho) / (rho);
				//	////vvzI= vvz+ c1o6 * (c1 / collFactorM - c1o2) * (Dxz * dX1_rho + Dyz * dX2_rho + 2 * dyuy * dX3_rho) / (rho);

				//	//vvxI = vvx + c1o3*forcingBIAS * (c1 / collFactorM - c1o2) * (2 * dxux * dX1_rho + Dxy * dX2_rho + Dxz * dX3_rho) / (rho);
				//	//vvyI = vvy + c1o3*forcingBIAS * (c1 / collFactorM - c1o2) * (Dxy * dX1_rho + 2 * dyuy * dX2_rho + Dyz * dX3_rho) / (rho);
				//	//vvzI = vvz + c1o3*forcingBIAS * (c1 / collFactorM - c1o2) * (Dxz * dX1_rho + Dyz * dX2_rho + 2 * dyuy * dX3_rho) / (rho);

				//	////vvxI = vvx - c1o3 * forcingBIAS * (c1 / collFactorM - c1o2) * (2 * dxux * dX1_rhoInv + Dxy * dX2_rhoInv + Dxz * dX3_rhoInv);
				//	////vvyI = vvy - c1o3 * forcingBIAS * (c1 / collFactorM - c1o2) * (Dxy * dX1_rhoInv + 2 * dyuy * dX2_rhoInv + Dyz * dX3_rhoInv);
				//	////vvzI = vvz - c1o3 * forcingBIAS * (c1 / collFactorM - c1o2) * (Dxz * dX1_rhoInv + Dyz * dX2_rhoInv + 2 * dyuy * dX3_rhoInv);


					//}
				//	//forcingX1+=(vvxI-vvx)/(deltaT* forcingBIAS) + muForcingX1.Eval() / rho;
				//	//forcingX2 += (vvyI - vvy) / (deltaT * forcingBIAS) + muForcingX2.Eval() / rho;
				//	//forcingX3 += (vvzI - vvz) / (deltaT * forcingBIAS) + muForcingX3.Eval() / rho;


				////	forcingX1 += c2 * (vvxI - vvx);
				////	forcingX2 += c2 * (vvyI - vvy);
				////	forcingX3 += c2 * (vvzI - vvz);


					//mfabb += c1o2*(-forcingX1) * c2o9;
					//mfbab += c1o2*(-forcingX2) * c2o9;
					//mfbba += c1o2*(-forcingX3) * c2o9;
					//mfaab += c1o2*(-forcingX1 - forcingX2) * c1o18;
					//mfcab += c1o2*( forcingX1 - forcingX2) * c1o18;
					//mfaba += c1o2*(-forcingX1 - forcingX3) * c1o18;
					//mfcba += c1o2*( forcingX1 - forcingX3) * c1o18;
					//mfbaa += c1o2*(-forcingX2 - forcingX3) * c1o18;
					//mfbca += c1o2*( forcingX2 - forcingX3) * c1o18;
					//mfaaa += c1o2*(-forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfcaa += c1o2*(forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfaca += c1o2*(-forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcca += c1o2*(forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcbb += c1o2*(forcingX1)*c2o9;
					//mfbcb += c1o2*(forcingX2)*c2o9;
					//mfbbc += c1o2*(forcingX3)*c2o9;
					//mfccb += c1o2*( forcingX1 + forcingX2) * c1o18;
					//mfacb += c1o2*(-forcingX1 + forcingX2) * c1o18;
					//mfcbc += c1o2*( forcingX1 + forcingX3) * c1o18;
					//mfabc += c1o2*(-forcingX1 + forcingX3) * c1o18;
					//mfbcc += c1o2*( forcingX2 + forcingX3) * c1o18;
					//mfbac += c1o2*(-forcingX2 + forcingX3) * c1o18;
					//mfccc += c1o2*(forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfacc += c1o2*(-forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfcac += c1o2*(forcingX1 - forcingX2 + forcingX3) * c1o72;
					//mfaac += c1o2*(-forcingX1 - forcingX2 + forcingX3) * c1o72;


					//forcingX1 = saveForceX1;
					//forcingX2 = saveForceX2;
					//forcingX3 = saveForceX3;
					vvx = vvxI;
					vvy = vvyI;
					vvz = vvzI;



					//!Abbas

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
						//mfabb -= pStar * c2o9;
						//mfbab -= pStar * c2o9;
						//mfbba -= pStar * c2o9;
						//mfaab -= pStar * c1o16;
						//mfcab -= pStar * c1o16;
						//mfaba -= pStar * c1o16;
						//mfcba -= pStar * c1o16;
						//mfbaa -= pStar * c1o16;
						//mfbca -= pStar * c1o16;
						//mfaaa -= pStar * c1o72;
						//mfcaa -= pStar * c1o72;
						//mfaca -= pStar * c1o72;
						//mfcca -= pStar * c1o72;
						//mfcbb -= pStar * c2o9;
						//mfbcb -= pStar * c2o9;
						//mfbbc -= pStar * c2o9;
						//mfccb -= pStar * c1o16;
						//mfacb -= pStar * c1o16;
						//mfcbc -= pStar * c1o16;
						//mfabc -= pStar * c1o16;
						//mfbcc -= pStar * c1o16;
						//mfbac -= pStar * c1o16;
						//mfccc -= pStar * c1o72;
						//mfacc -= pStar * c1o72;
						//mfcac -= pStar * c1o72;
						//mfaac -= pStar * c1o72;
						//mfbbb -= pStar * 8.0/9.0;
					///////////////////

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
					real qudricLimit = 0.01 / (c1o1 + 1.0e4 * phi[DIR_000] * (c1o1 - phi[DIR_000])); //real qudricLimit = 0.01;
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

					// mfaaa = 0.0;
					real OxxPyyPzz = 1.0; //omega2 or bulk viscosity
											//  LBMReal OxyyPxzz = 1.;//-s9;//2+s9;//
											//  LBMReal OxyyMxzz  = 1.;//2+s9;//
					real O4 = 1.;
					real O5 = 1.;
					real O6 = 1.;

					//collFactorM+= (1.7 - collFactorM) * fabs(mfaaa) / (fabs(mfaaa) + 0.001f);


					/////fourth order parameters; here only for test. Move out of loop!

					real OxyyPxzz = 8.0 * (collFactorM - 2.0) * (OxxPyyPzz * (3.0 * collFactorM - 1.0) - 5.0 * collFactorM) / (8.0 * (5.0 - 2.0 * collFactorM) * collFactorM + OxxPyyPzz * (8.0 + collFactorM * (9.0 * collFactorM - 26.0)));
					real OxyyMxzz = 8.0 * (collFactorM - 2.0) * (collFactorM + OxxPyyPzz * (3.0 * collFactorM - 7.0)) / (OxxPyyPzz * (56.0 - 42.0 * collFactorM + 9.0 * collFactorM * collFactorM) - 8.0 * collFactorM);
				    real Oxyz = 24.0 * (collFactorM - 2.0) * (4.0 * collFactorM * collFactorM + collFactorM * OxxPyyPzz * (18.0 - 13.0 * collFactorM) + OxxPyyPzz * OxxPyyPzz * (2.0 + collFactorM * (6.0 * collFactorM - 11.0))) / (16.0 * collFactorM * collFactorM * (collFactorM - 6.0) - 2.0 * collFactorM * OxxPyyPzz * (216.0 + 5.0 * collFactorM * (9.0 * collFactorM - 46.0)) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (3.0 * collFactorM - 10.0) * (15.0 * collFactorM - 28.0) - 48.0));
					real A = (4.0 * collFactorM * collFactorM + 2.0 * collFactorM * OxxPyyPzz * (collFactorM - 6.0) + OxxPyyPzz * OxxPyyPzz * (collFactorM * (10.0 - 3.0 * collFactorM) - 4.0)) / ((collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//FIXME:  warning C4459: declaration of 'B' hides global declaration (message : see declaration of 'D3Q27System::B' )
					real BB = (4.0 * collFactorM * OxxPyyPzz * (9.0 * collFactorM - 16.0) - 4.0 * collFactorM * collFactorM - 2.0 * OxxPyyPzz * OxxPyyPzz * (2.0 + 9.0 * collFactorM * (collFactorM - 2.0))) / (3.0 * (collFactorM - OxxPyyPzz) * (OxxPyyPzz * (2.0 + 3.0 * collFactorM) - 8.0 * collFactorM));
					//LBMReal stress = 1.0;// stress / (stress + 1.0e-10);
					//stress = 1.0;
					//OxyyPxzz += stress*(1.0-OxyyPxzz);
					//OxyyPxzz = c3 * (collFactorM - c2) / (collFactorM - c3);
					//OxyyMxzz += stress*(1.0-OxyyMxzz);
					//Oxyz +=  stress*(1.0-Oxyz);
					//A *= 1.0-stress;
					//BB *= 1.0-stress;

					//Cum 4.
					//LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
					//LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
					//LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

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
					//  LBMReal mfaaaS = (mfaaa * (-4 - 3 * OxxPyyPzz * (-1 + rho)) + 6 * mxxPyyPzz * OxxPyyPzz * (-1 + rho)) / (-4 + 3 * OxxPyyPzz * (-1 + rho));
					mxxPyyPzz -= mfaaa ;//12.03.21 shifted by mfaaa
										//mxxPyyPzz-=(mfaaa+mfaaaS)*c1o2;//12.03.21 shifted by mfaaa
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

					real dxux =  -c1o2 * collFactorM * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (/*mfaaa*/ -mxxPyyPzz)*0;
					//LBMReal dxux = -c1o2 * (mxxMyy + mxxMzz) * collFactorM - mfaaa * c1o3* omegaDRho;
					real dyuy =  dxux + collFactorM * c3o2 * mxxMyy;
					real dzuz =  dxux + collFactorM * c3o2 * mxxMzz;
					real Dxy = -c3o1 * collFactorM * mfbba;
					real Dxz = -c3o1 * collFactorM * mfbab;
					real Dyz = -c3o1 * collFactorM * mfabb;
//					// attempt to improve implicit  stress computation by fixed iteration
//					LBMReal dX2_rho = (rhoToPhi)*dX2_phi;
//					LBMReal dX1_rho = (rhoToPhi)*dX1_phi;
//					LBMReal dX3_rho = (rhoToPhi)*dX3_phi;
//
//						LBMReal dfx= c1o3 * (c1 / collFactorM - c1o2) *(2 * dxux * dX1_rho + Dxy * dX2_rho + Dxz * dX3_rho) / (rho);
//						LBMReal dfy = c1o3 * (c1 / collFactorM - c1o2) *(Dxy * dX1_rho + 2 * dyuy * dX2_rho + Dyz * dX3_rho) / (rho);
//						LBMReal dfz = c1o3 * (c1 / collFactorM - c1o2) *(Dxz * dX1_rho + Dyz * dX2_rho + 2 * dyuy * dX3_rho) / (rho);
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

					//if (fabsf(mfaaa + (dxux + dyuy + dzuz) > 1e-9)){
					//	std::cout << mfaaa <<" "<< (dxux + dyuy + dzuz)<< std::endl;
					//}


					////updated pressure
					//mfaaa += (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling;
					//mfaaa *= (one-omegaDRho);// (mfaaa + (dxux + dyuy + dzuz)) * .5; // Pressure elimination as in standard velocity model
								 //  mfaaa += (rho - c1) * (dxux + dyuy + dzuz);
				
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


					////////


					////////////////////////////////////////////////////////////////////////////////////
					//forcing
					mfbaa = -mfbaa;// *(c1 - forcingBIAS) / forcingBIAS;
					mfaba = -mfaba;// *(c1 - forcingBIAS) / forcingBIAS;
					mfaab = -mfaab;// *(c1 - forcingBIAS) / forcingBIAS;

					//mfbaa += c1o3 * (c1 / collFactorM - c1o2) * rhoToPhi * (2 * dxux * dX1_phi + Dxy * dX2_phi + Dxz * dX3_phi) / (rho);
					//mfaba += c1o3 * (c1 / collFactorM - c1o2) * rhoToPhi * (Dxy * dX1_phi + 2 * dyuy * dX2_phi + Dyz * dX3_phi) / (rho);
					//mfaab += c1o3 * (c1 / collFactorM - c1o2) * rhoToPhi * (Dxz * dX1_phi + Dyz * dX2_phi + 2 * dyuy * dX3_phi) / (rho);

					mfbaa -= c1o2 * rhoToPhi * (mmfcaa* dX1_phi + mmfbba * dX2_phi + mmfbab * dX3_phi) / (rho);
					mfaba -= c1o2 * rhoToPhi * (mmfbba* dX1_phi + mmfaca * dX2_phi + mmfabb * dX3_phi) / (rho);
					mfaab -= c1o2 * rhoToPhi * (mmfbab* dX1_phi + mmfabb * dX2_phi + mmfaac * dX3_phi) / (rho);
					
					vvx -= c1o4 * rhoToPhi * (mmfcaa * dX1_phi + mmfbba * dX2_phi + mmfbab * dX3_phi) / (rho);
					vvy -= c1o4 * rhoToPhi * (mmfbba * dX1_phi + mmfaca * dX2_phi + mmfabb * dX3_phi) / (rho);
					vvz -= c1o4 * rhoToPhi * (mmfbab * dX1_phi + mmfabb * dX2_phi + mmfaac * dX3_phi) / (rho);

					vx2 = vvx * vvx;
					vy2 = vvy * vvy;
					vz2 = vvz * vvz;

					//mmfcaa =0;// c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz - mfaaa);
					//mmfaca =0;// c1o3 * (-2. * mxxMyy + mxxMzz + mxxPyyPzz - mfaaa);
					//mmfaac =0;// c1o3 * (mxxMyy - 2. * mxxMzz + mxxPyyPzz - mfaaa);
					//mmfabb =0;// mfabb;
					//mmfbab =0;// mfbab;
					//mmfbba =0;// mfbba;


					//////////////////////////////////////////////////////////////////////////////////////
					//grad Rho
					//LBMReal dX1_rho = (rhoToPhi - three * (*pressure)(x1, x2, x3)) * dX1_phi - phi[DIR_000] * three * gradPx;
					//LBMReal dX2_rho = (rhoToPhi - three * (*pressure)(x1, x2, x3)) * dX2_phi - phi[DIR_000] * three * gradPy;
					//LBMReal dX3_rho = (rhoToPhi - three * (*pressure)(x1, x2, x3)) * dX3_phi - phi[DIR_000] * three * gradPz;

					//LBMReal dX2_rho = (rhoToPhi ) * dX2_phi ;
					//LBMReal dX1_rho = (rhoToPhi ) * dX1_phi ;
					//LBMReal dX3_rho = (rhoToPhi ) * dX3_phi ;
					///////////////////////////////////////////////////////////////////////////////////////
					//mfbaa += c1o3 * (c1 / collFactorM - c1o2) *(2 * dxux * dX1_rho + Dxy * dX2_rho + Dxz * dX3_rho) / (rho);
					//mfaba += c1o3 * (c1 / collFactorM - c1o2) *(Dxy * dX1_rho + 2 * dyuy * dX2_rho + Dyz * dX3_rho) / (rho);
					//mfaab += c1o3 * (c1 / collFactorM - c1o2) *(Dxz * dX1_rho + Dyz * dX2_rho + 2 * dyuy * dX3_rho) / (rho);
					
					///////Fakhari pressure correction
					//mfbaa -= mfaaa / rho * dX1_rho*c1o3;
					//mfaba -= mfaaa / rho * dX2_rho*c1o3;
					//mfaab -= mfaaa / rho * dX3_rho*c1o3;
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
					/////SimpleForce

					//mfabb += c2o9 * deltaP;
					//mfbab += c2o9 * deltaP;
					//mfbba += c2o9 * deltaP;
					//mfaab += c1o18 * deltaP;
					//mfcab += c1o18 * deltaP;
					//mfaba += c1o18 * deltaP;
					//mfcba += c1o18 * deltaP;
					//mfbaa += c1o18 * deltaP;
					//mfbca += c1o18 * deltaP;
					//mfaaa += c1o72 * deltaP;
					//mfcaa += c1o72 * deltaP;
					//mfaca += c1o72 * deltaP;
					//mfcca += c1o72 * deltaP;
					//mfcbb += c2o9 * deltaP;
					//mfbcb += c2o9 * deltaP;
					//mfbbc += c2o9 * deltaP;
					//mfccb += c1o18 * deltaP;
					//mfacb += c1o18 * deltaP;
					//mfcbc += c1o18 * deltaP;
					//mfabc += c1o18 * deltaP;
					//mfbcc += c1o18 * deltaP;
					//mfbac += c1o18 * deltaP;
					//mfccc += c1o72 * deltaP;
					//mfacc += c1o72 * deltaP;
					//mfcac += c1o72 * deltaP;
					//mfaac += c1o72 * deltaP;

					//mfabb += c1o2*(-forcingX1                        ) * c2o9;
					//mfbab += c1o2*(           - forcingX2            ) * c2o9;
					//mfbba += c1o2*(                       - forcingX3) * c2o9;
					//mfaab += c1o2*(-forcingX1 - forcingX2            ) * c1o18;
					//mfcab += c1o2*( forcingX1 - forcingX2            ) * c1o18;
					//mfaba += c1o2*(-forcingX1             - forcingX3) * c1o18;
					//mfcba += c1o2*( forcingX1             - forcingX3) * c1o18;
					//mfbaa += c1o2*(           - forcingX2 - forcingX3) * c1o18;
					//mfbca += c1o2*(             forcingX2 - forcingX3) * c1o18;
					//mfaaa += c1o2*(-forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfcaa += c1o2*( forcingX1 - forcingX2 - forcingX3) * c1o72;
					//mfaca += c1o2*(-forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcca += c1o2*( forcingX1 + forcingX2 - forcingX3) * c1o72;
					//mfcbb += c1o2*( forcingX1                        ) * c2o9;
					//mfbcb += c1o2*(             forcingX2            ) * c2o9;
					//mfbbc += c1o2*(                         forcingX3) * c2o9;
					//mfccb += c1o2*( forcingX1 + forcingX2            ) * c1o18;
					//mfacb += c1o2*(-forcingX1 + forcingX2            ) * c1o18;
					//mfcbc += c1o2*( forcingX1             + forcingX3) * c1o18;
					//mfabc += c1o2*(-forcingX1             + forcingX3) * c1o18;
					//mfbcc += c1o2*(             forcingX2 + forcingX3) * c1o18;
					//mfbac += c1o2*(           - forcingX2 + forcingX3) * c1o18;
					//mfccc += c1o2*( forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfacc += c1o2*(-forcingX1 + forcingX2 + forcingX3) * c1o72;
					//mfcac += c1o2*( forcingX1 - forcingX2 + forcingX3) * c1o72;
					//mfaac += c1o2*(-forcingX1 - forcingX2 + forcingX3) * c1o72;
					//pStarStart -= (vx2 + vy2 + vz2) * c1o3;

					///Take the diffusion part with out

					//mfStartcbb -= D3Q27System::getIncompFeqForDirection(D3Q27System::E  , zeroReal, vvx, vvy, vvz);
					//mfStartbcb -= D3Q27System::getIncompFeqForDirection(D3Q27System::N  , zeroReal, vvx, vvy, vvz);
					//mfStartbbc -= D3Q27System::getIncompFeqForDirection(D3Q27System::T  , zeroReal, vvx, vvy, vvz);
					//mfStartccb -= D3Q27System::getIncompFeqForDirection(D3Q27System::NE , zeroReal, vvx, vvy, vvz);
					//mfStartacb -= D3Q27System::getIncompFeqForDirection(D3Q27System::NW , zeroReal, vvx, vvy, vvz);
					//mfStartcbc -= D3Q27System::getIncompFeqForDirection(D3Q27System::TE , zeroReal, vvx, vvy, vvz);
					//mfStartabc -= D3Q27System::getIncompFeqForDirection(D3Q27System::TW , zeroReal, vvx, vvy, vvz);
					//mfStartbcc -= D3Q27System::getIncompFeqForDirection(D3Q27System::TN , zeroReal, vvx, vvy, vvz);
					//mfStartbac -= D3Q27System::getIncompFeqForDirection(D3Q27System::TS , zeroReal, vvx, vvy, vvz);
					//mfStartccc -= D3Q27System::getIncompFeqForDirection(D3Q27System::TNE, zeroReal, vvx, vvy, vvz);
					//mfStartacc -= D3Q27System::getIncompFeqForDirection(D3Q27System::TNW, zeroReal, vvx, vvy, vvz);
					//mfStartcac -= D3Q27System::getIncompFeqForDirection(D3Q27System::TSE, zeroReal, vvx, vvy, vvz);
					//mfStartaac -= D3Q27System::getIncompFeqForDirection(D3Q27System::TSW, zeroReal, vvx, vvy, vvz);
					//mfStartabb -= D3Q27System::getIncompFeqForDirection(D3Q27System::W  , zeroReal, vvx, vvy, vvz);
					//mfStartbab -= D3Q27System::getIncompFeqForDirection(D3Q27System::S  , zeroReal, vvx, vvy, vvz);
					//mfStartbba -= D3Q27System::getIncompFeqForDirection(D3Q27System::B  , zeroReal, vvx, vvy, vvz);
					//mfStartaab -= D3Q27System::getIncompFeqForDirection(D3Q27System::SW , zeroReal, vvx, vvy, vvz);
					//mfStartcab -= D3Q27System::getIncompFeqForDirection(D3Q27System::SE , zeroReal, vvx, vvy, vvz);
					//mfStartaba -= D3Q27System::getIncompFeqForDirection(D3Q27System::BW , zeroReal, vvx, vvy, vvz);
					//mfStartcba -= D3Q27System::getIncompFeqForDirection(D3Q27System::BE , zeroReal, vvx, vvy, vvz);
					//mfStartbaa -= D3Q27System::getIncompFeqForDirection(D3Q27System::BS , zeroReal, vvx, vvy, vvz);
					//mfStartbca -= D3Q27System::getIncompFeqForDirection(D3Q27System::BN , zeroReal, vvx, vvy, vvz);
					//mfStartaaa -= D3Q27System::getIncompFeqForDirection(D3Q27System::BSW, zeroReal, vvx, vvy, vvz);
					//mfStartcaa -= D3Q27System::getIncompFeqForDirection(D3Q27System::BSE, zeroReal, vvx, vvy, vvz);
					//mfStartaca -= D3Q27System::getIncompFeqForDirection(D3Q27System::BNW, zeroReal, vvx, vvy, vvz);
					//mfStartcca -= D3Q27System::getIncompFeqForDirection(D3Q27System::BNE, zeroReal, vvx, vvy, vvz);
					//mfStartbbb -= D3Q27System::getIncompFeqForDirection(D3Q27System::REST, zeroReal, vvx, vvy, vvz);
					//
					//pStar += pStarStart*(omegaDRho-c1);

					//mfStartcbb = c2o9 * pStar;
					//	mfStartbcb= c2o9 * pStar;
					//	mfStartbbc= c2o9 * pStar;
					//	mfStartccb= c1o18 * pStar;
					//	mfStartacb= c1o18 * pStar;
					//	mfStartcbc= c1o18 * pStar;
					//	mfStartabc= c1o18 * pStar;
					//	mfStartbcc= c1o18 * pStar;
					//	mfStartbac= c1o18 * pStar;
					//	mfStartccc= c1o72 * pStar;
					//	mfStartacc= c1o72 * pStar;
					//	mfStartcac= c1o72 * pStar;
					//	mfStartaac= c1o72 * pStar;
					//	mfStartabb= c2o9 * pStar;
					//	mfStartbab= c2o9 * pStar;
					//	mfStartbba= c2o9 * pStar;
					//	mfStartaab= c1o18 * pStar;
					//	mfStartcab= c1o18 * pStar;
					//	mfStartaba= c1o18 * pStar;
					//	mfStartcba= c1o18 * pStar;
					//	mfStartbaa= c1o18 * pStar;
					//	mfStartbca= c1o18 * pStar;
					//	mfStartaaa= c1o72 * pStar;
					//	mfStartcaa= c1o72 * pStar;
					//	mfStartaca= c1o72 * pStar;
					//	mfStartcca= c1o72 * pStar;
					//	mfStartbbb= c4 * c2o9 * pStar;

					//mfaaa -= c1o2 * (mfStartaaa + mfStartccc)+ c1o72 * (mmfaac + c3 * mmfabb + mmfaca + c3 * mmfbab + c3 * mmfbba + mmfcaa);
					//mfaab -= c1o2 * (mfStartaab + mfStartccb)+c1o36 * (-mmfaac + c2 * (mmfaca + c3 * mmfbba + mmfcaa));
					//mfaac -= c1o2 * (mfStartaac + mfStartcca)+c1o72 * (mmfaac - c3 * mmfabb + mmfaca - c3 * mmfbab + c3 * mmfbba + mmfcaa);
					//mfaba -= c1o2 * (mfStartaba + mfStartcbc)+c1o36 * (c2 * mmfaac - mmfaca + c6 * mmfbab + c2 * mmfcaa);
					//mfabb -= c1o2 * (mfStartabb + mfStartcbb)+c1o9 * (-mmfaac - mmfaca + c2 * mmfcaa);
					//mfabc -= c1o2 * (mfStartabc + mfStartcba)+c1o36 * (c2 * mmfaac - mmfaca - 6 * mmfbab + c2 * mmfcaa);
					//mfaca -= c1o2 * (mfStartaca + mfStartcac)+c1o72 * (mmfaac - c3 * mmfabb + mmfaca + c3 * mmfbab - c3 * mmfbba + mmfcaa);
					//mfacb -= c1o2 * (mfStartacb + mfStartcab)+c1o36 * (-mmfaac + c2 * (mmfaca - c3 * mmfbba + mmfcaa));
					//mfacc -= c1o2 * (mfStartacc + mfStartcaa)+c1o72 * (mmfaac + c3 * mmfabb + mmfaca - c3 * mmfbab - c3 * mmfbba + mmfcaa);
					//mfbaa -= c1o2 * (mfStartbaa + mfStartbcc)+c1o36 * (c2 * mmfaac + c6 * mmfabb + c2 * mmfaca - mmfcaa);
					//mfbab -= c1o2 * (mfStartbab + mfStartbcb)+c1o9 * (-mmfaac + c2 * mmfaca - mmfcaa);
					//mfbac -= c1o2 * (mfStartbac + mfStartbca)+c1o36 * (c2 * mmfaac - 6 * mmfabb + c2 * mmfaca - mmfcaa);
					//mfbba -= c1o2 * (mfStartbba + mfStartbbc)+c1o9 * (c2 * mmfaac - mmfaca - mmfcaa);
					//mfbbb -=  (mfStartbbb)-(c4o9 * (mmfaac + mmfaca + mmfcaa));
					//mfbbc -= c1o2 * (mfStartbbc + mfStartbba)+c1o9 * (c2 * mmfaac - mmfaca - mmfcaa);
					//mfbca -= c1o2 * (mfStartbca + mfStartbac)+c1o36 * (c2 * mmfaac - 6 * mmfabb + c2 * mmfaca - mmfcaa);
					//mfbcb -= c1o2 * (mfStartbcb + mfStartbab)+c1o9 * (-mmfaac + c2 * mmfaca - mmfcaa);
					//mfbcc -= c1o2 * (mfStartbcc + mfStartbaa)+c1o36 * (c2 * mmfaac + c6 * mmfabb + c2 * mmfaca - mmfcaa);
					//mfcaa -= c1o2 * (mfStartcaa + mfStartacc)+c1o72 * (mmfaac + c3 * mmfabb + mmfaca - c3 * mmfbab - c3 * mmfbba + mmfcaa);
					//mfcab -= c1o2 * (mfStartcab + mfStartacb)+c1o36 * (-mmfaac + c2 * (mmfaca - c3 * mmfbba + mmfcaa));
					//mfcac -= c1o2 * (mfStartcac + mfStartaca)+c1o72 * (mmfaac - c3 * mmfabb + mmfaca + c3 * mmfbab - c3 * mmfbba + mmfcaa);
					//mfcba -= c1o2 * (mfStartcba + mfStartabc)+c1o36 * (c2 * mmfaac - mmfaca - 6 * mmfbab + c2 * mmfcaa);
					//mfcbb -= c1o2 * (mfStartcbb + mfStartabb)+c1o9 * (-mmfaac - mmfaca + c2 * mmfcaa);
					//mfcbc -= c1o2 * (mfStartcbc + mfStartaba)+c1o36 * (c2 * mmfaac - mmfaca + c6 * mmfbab + c2 * mmfcaa);
					//mfcca -= c1o2 * (mfStartcca + mfStartaac)+c1o72 * (mmfaac - c3 * mmfabb + mmfaca - c3 * mmfbab + c3 * mmfbba + mmfcaa);
					//mfccb -= c1o2 * (mfStartccb + mfStartaab)+c1o36 * (-mmfaac + c2 * (mmfaca + c3 * mmfbba + mmfcaa));
					//mfccc -= c1o2 * (mfStartccc + mfStartaaa)+c1o72 * (mmfaac + c3 * mmfabb + mmfaca + c3 * mmfbab + c3 * mmfbba + mmfcaa);

					//mfhaaa =rho*( c1o2 * (mfStartaaa + mfStartccc) + c1o72 * (mmfaac + c3 * mmfabb + mmfaca + c3 * mmfbab + c3 * mmfbba + mmfcaa));
					//mfhaab =rho*( c1o2 * (mfStartaab + mfStartccb) + c1o36 * (-mmfaac + c2 * (mmfaca + c3 * mmfbba + mmfcaa)));
					//mfhaac =rho*( c1o2 * (mfStartaac + mfStartcca) + c1o72 * (mmfaac - c3 * mmfabb + mmfaca - c3 * mmfbab + c3 * mmfbba + mmfcaa));
					//mfhaba =rho*( c1o2 * (mfStartaba + mfStartcbc) + c1o36 * (c2 * mmfaac - mmfaca + c6 * mmfbab + c2 * mmfcaa));
					//mfhabb =rho*( c1o2 * (mfStartabb + mfStartcbb) + c1o9 * (-mmfaac - mmfaca + c2 * mmfcaa));
					//mfhabc =rho*( c1o2 * (mfStartabc + mfStartcba) + c1o36 * (c2 * mmfaac - mmfaca - 6 * mmfbab + c2 * mmfcaa));
					//mfhaca =rho*( c1o2 * (mfStartaca + mfStartcac) + c1o72 * (mmfaac - c3 * mmfabb + mmfaca + c3 * mmfbab - c3 * mmfbba + mmfcaa));
					//mfhacb =rho*( c1o2 * (mfStartacb + mfStartcab) + c1o36 * (-mmfaac + c2 * (mmfaca - c3 * mmfbba + mmfcaa)));
					//mfhacc =rho*( c1o2 * (mfStartacc + mfStartcaa) + c1o72 * (mmfaac + c3 * mmfabb + mmfaca - c3 * mmfbab - c3 * mmfbba + mmfcaa));
					//mfhbaa =rho*( c1o2 * (mfStartbaa + mfStartbcc) + c1o36 * (c2 * mmfaac + c6 * mmfabb + c2 * mmfaca - mmfcaa));
					//mfhbab =rho*( c1o2 * (mfStartbab + mfStartbcb) + c1o9 * (-mmfaac + c2 * mmfaca - mmfcaa));
					//mfhbac =rho*( c1o2 * (mfStartbac + mfStartbca) + c1o36 * (c2 * mmfaac - 6 * mmfabb + c2 * mmfaca - mmfcaa));
					//mfhbba =rho*( c1o2 * (mfStartbba + mfStartbbc) + c1o9 * (c2 * mmfaac - mmfaca - mmfcaa));
					//mfhbbb =rho*( (mfStartbbb)-(c4o9 * (mmfaac + mmfaca + mmfcaa)));
					//mfhbbc =rho*( c1o2 * (mfStartbbc + mfStartbba) + c1o9 * (c2 * mmfaac - mmfaca - mmfcaa));
					//mfhbca =rho*( c1o2 * (mfStartbca + mfStartbac) + c1o36 * (c2 * mmfaac - 6 * mmfabb + c2 * mmfaca - mmfcaa));
					//mfhbcb =rho*( c1o2 * (mfStartbcb + mfStartbab) + c1o9 * (-mmfaac + c2 * mmfaca - mmfcaa));
					//mfhbcc =rho*( c1o2 * (mfStartbcc + mfStartbaa) + c1o36 * (c2 * mmfaac + c6 * mmfabb + c2 * mmfaca - mmfcaa));
					//mfhcaa =rho*( c1o2 * (mfStartcaa + mfStartacc) + c1o72 * (mmfaac + c3 * mmfabb + mmfaca - c3 * mmfbab - c3 * mmfbba + mmfcaa));
					//mfhcab =rho*( c1o2 * (mfStartcab + mfStartacb) + c1o36 * (-mmfaac + c2 * (mmfaca - c3 * mmfbba + mmfcaa)));
					//mfhcac =rho*( c1o2 * (mfStartcac + mfStartaca) + c1o72 * (mmfaac - c3 * mmfabb + mmfaca + c3 * mmfbab - c3 * mmfbba + mmfcaa));
					//mfhcba =rho*( c1o2 * (mfStartcba + mfStartabc) + c1o36 * (c2 * mmfaac - mmfaca - 6 * mmfbab + c2 * mmfcaa));
					//mfhcbb =rho*( c1o2 * (mfStartcbb + mfStartabb) + c1o9 * (-mmfaac - mmfaca + c2 * mmfcaa));
					//mfhcbc =rho*( c1o2 * (mfStartcbc + mfStartaba) + c1o36 * (c2 * mmfaac - mmfaca + c6 * mmfbab + c2 * mmfcaa));
					//mfhcca =rho*( c1o2 * (mfStartcca + mfStartaac) + c1o72 * (mmfaac - c3 * mmfabb + mmfaca - c3 * mmfbab + c3 * mmfbba + mmfcaa));
					//mfhccb =rho*( c1o2 * (mfStartccb + mfStartaab) + c1o36 * (-mmfaac + c2 * (mmfaca + c3 * mmfbba + mmfcaa)));
					//mfhccc =rho*( c1o2 * (mfStartccc + mfStartaaa) + c1o72 * (mmfaac + c3 * mmfabb + mmfaca + c3 * mmfbab + c3 * mmfbba + mmfcaa));




					pStar += pStarStart*(omegaDRho- c1o1);

					mfcbb -= c2o9*pStar;
					mfbcb -= c2o9*pStar;
					mfbbc -= c2o9*pStar;
					mfccb -= c1o18*pStar;
					mfacb -= c1o18*pStar;
					mfcbc -= c1o18*pStar;
					mfabc -= c1o18*pStar;
					mfbcc -= c1o18*pStar;
					mfbac -= c1o18*pStar;
					mfccc -= c1o72*pStar;
					mfacc -= c1o72*pStar;
					mfcac -= c1o72*pStar;
					mfaac -= c1o72*pStar;
					mfabb -= c2o9*pStar;
					mfbab -= c2o9*pStar;
					mfbba -= c2o9*pStar;
					mfaab -= c1o18*pStar;
					mfcab -= c1o18*pStar;
					mfaba -= c1o18*pStar;
					mfcba -= c1o18*pStar;
					mfbaa -= c1o18*pStar;
					mfbca -= c1o18*pStar;
					mfaaa -= c1o72*pStar;
					mfcaa -= c1o72*pStar;
					mfaca -= c1o72*pStar;
					mfcca -= c1o72*pStar;
					mfbbb -= c4o1*c2o9*pStar;

					mfhbcb = rho*c2o9 * pStar;
					mfhbbc = rho*c2o9 * pStar;
					mfhcbb = rho*c2o9 * pStar;
					mfhccb = rho*c1o18 * pStar;
					mfhacb = rho*c1o18 * pStar;
					mfhcbc = rho*c1o18 * pStar;
					mfhabc = rho*c1o18 * pStar;
					mfhbcc = rho*c1o18 * pStar;
					mfhbac = rho*c1o18 * pStar;
					mfhccc = rho*c1o72 * pStar;
					mfhacc = rho*c1o72 * pStar;
					mfhcac = rho*c1o72 * pStar;
					mfhaac = rho*c1o72 * pStar;
					mfhabb = rho*c2o9 * pStar;
					mfhbab = rho*c2o9 * pStar;
					mfhbba = rho*c2o9 * pStar;
					mfhaab = rho*c1o18 * pStar;
					mfhcab = rho*c1o18 * pStar;
					mfhaba = rho*c1o18 * pStar;
					mfhcba = rho*c1o18 * pStar;
					mfhbaa = rho*c1o18 * pStar;
					mfhbca = rho*c1o18 * pStar;
					mfhaaa = rho*c1o72 * pStar;
					mfhcaa = rho*c1o72 * pStar;
					mfhaca = rho*c1o72 * pStar;
					mfhcca = rho*c1o72 * pStar;
					mfhbbb = rho* c4o1 * c2o9 * pStar;

					//mfStartbcb =  c2o9  * pStarStart;
					//mfStartbbc =  c2o9  * pStarStart;
					//mfStartcbb =  c2o9  * pStarStart;
					//mfStartccb =  c1o18 * pStarStart;
					//mfStartacb =  c1o18 * pStarStart;
					//mfStartcbc =  c1o18 * pStarStart;
					//mfStartabc =  c1o18 * pStarStart;
					//mfStartbcc =  c1o18 * pStarStart;
					//mfStartbac =  c1o18 * pStarStart;
					//mfStartccc =  c1o72 * pStarStart;
					//mfStartacc =  c1o72 * pStarStart;
					//mfStartcac =  c1o72 * pStarStart;
					//mfStartaac =  c1o72 * pStarStart;
					//mfStartabb =  c2o9  * pStarStart;
					//mfStartbab =  c2o9  * pStarStart;
					//mfStartbba =  c2o9  * pStarStart;
					//mfStartaab =  c1o18 * pStarStart;
					//mfStartcab =  c1o18 * pStarStart;
					//mfStartaba =  c1o18 * pStarStart;
					//mfStartcba =  c1o18 * pStarStart;
					//mfStartbaa =  c1o18 * pStarStart;
					//mfStartbca =  c1o18 * pStarStart;
					//mfStartaaa =  c1o72 * pStarStart;
					//mfStartcaa =  c1o72 * pStarStart;
					//mfStartaca =  c1o72 * pStarStart;
					//mfStartcca =  c1o72 * pStarStart;
					//mfStartbbb =  c4 * c2o9 * pStarStart;

					//LBMReal scaleSplit = 0.5;
					//mfStartbcb = mfStartbcb*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartbbc = mfStartbbc*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartcbb = mfStartcbb*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartccb = mfStartccb*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartacb = mfStartacb*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartcbc = mfStartcbc*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartabc = mfStartabc*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartbcc = mfStartbcc*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartbac = mfStartbac*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartccc = mfStartccc*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartacc = mfStartacc*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartcac = mfStartcac*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartaac = mfStartaac*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartabb = mfStartabb*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartbab = mfStartbab*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartbba = mfStartbba*scaleSplit+(c1-scaleSplit)* c2o9 * pStarStart;
					//mfStartaab = mfStartaab*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartcab = mfStartcab*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartaba = mfStartaba*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartcba = mfStartcba*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartbaa = mfStartbaa*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartbca = mfStartbca*scaleSplit+(c1-scaleSplit)* c1o18 * pStarStart;
					//mfStartaaa = mfStartaaa*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartcaa = mfStartcaa*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartaca = mfStartaca*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartcca = mfStartcca*scaleSplit+(c1-scaleSplit)* c1o72 * pStarStart;
					//mfStartbbb = mfStartbbb*scaleSplit+(c1-scaleSplit)* c4 * c2o9 * pStarStart;


					//mfaaa -= c1o2 * (mfStartaaa + mfStartccc);
     //               mfaab -= c1o2 * (mfStartaab + mfStartccb);
     //               mfaac -= c1o2 * (mfStartaac + mfStartcca);
     //               mfaba -= c1o2 * (mfStartaba + mfStartcbc);
     //               mfabb -= c1o2 * (mfStartabb + mfStartcbb);
     //               mfabc -= c1o2 * (mfStartabc + mfStartcba);
     //               mfaca -= c1o2 * (mfStartaca + mfStartcac);
     //               mfacb -= c1o2 * (mfStartacb + mfStartcab);
     //               mfacc -= c1o2 * (mfStartacc + mfStartcaa);
     //               mfbaa -= c1o2 * (mfStartbaa + mfStartbcc);
     //               mfbab -= c1o2 * (mfStartbab + mfStartbcb);
     //               mfbac -= c1o2 * (mfStartbac + mfStartbca);
     //               mfbba -= c1o2 * (mfStartbba + mfStartbbc);
					//mfbbb -= (mfStartbbb);
     //               mfbbc -= c1o2 * (mfStartbbc + mfStartbba);
     //               mfbca -= c1o2 * (mfStartbca + mfStartbac);
     //               mfbcb -= c1o2 * (mfStartbcb + mfStartbab);
     //               mfbcc -= c1o2 * (mfStartbcc + mfStartbaa);
     //               mfcaa -= c1o2 * (mfStartcaa + mfStartacc);
     //               mfcab -= c1o2 * (mfStartcab + mfStartacb);
     //               mfcac -= c1o2 * (mfStartcac + mfStartaca);
     //               mfcba -= c1o2 * (mfStartcba + mfStartabc);
     //               mfcbb -= c1o2 * (mfStartcbb + mfStartabb);
     //               mfcbc -= c1o2 * (mfStartcbc + mfStartaba);
     //               mfcca -= c1o2 * (mfStartcca + mfStartaac);
     //               mfccb -= c1o2 * (mfStartccb + mfStartaab);
     //               mfccc -= c1o2 * (mfStartccc + mfStartaaa);
					//												
					//mfhaaa += rho*c1o2 * (mfStartaaa + mfStartccc);
					//mfhaab += rho*c1o2 * (mfStartaab + mfStartccb);
					//mfhaac += rho*c1o2 * (mfStartaac + mfStartcca);
					//mfhaba += rho*c1o2 * (mfStartaba + mfStartcbc);
					//mfhabb += rho*c1o2 * (mfStartabb + mfStartcbb);
					//mfhabc += rho*c1o2 * (mfStartabc + mfStartcba);
					//mfhaca += rho*c1o2 * (mfStartaca + mfStartcac);
					//mfhacb += rho*c1o2 * (mfStartacb + mfStartcab);
					//mfhacc += rho*c1o2 * (mfStartacc + mfStartcaa);
					//mfhbaa += rho*c1o2 * (mfStartbaa + mfStartbcc);
					//mfhbab += rho*c1o2 * (mfStartbab + mfStartbcb);
					//mfhbac += rho*c1o2 * (mfStartbac + mfStartbca);
					//mfhbba += rho*c1o2 * (mfStartbba + mfStartbbc);
					//mfhbbb += rho*(mfStartbbb);
					//mfhbbc += rho*c1o2 * (mfStartbbc + mfStartbba);
					//mfhbca += rho*c1o2 * (mfStartbca + mfStartbac);
					//mfhbcb += rho*c1o2 * (mfStartbcb + mfStartbab);
					//mfhbcc += rho*c1o2 * (mfStartbcc + mfStartbaa);
					//mfhcaa += rho*c1o2 * (mfStartcaa + mfStartacc);
					//mfhcab += rho*c1o2 * (mfStartcab + mfStartacb);
					//mfhcac += rho*c1o2 * (mfStartcac + mfStartaca);
					//mfhcba += rho*c1o2 * (mfStartcba + mfStartabc);
					//mfhcbb += rho*c1o2 * (mfStartcbb + mfStartabb);
					//mfhcbc += rho*c1o2 * (mfStartcbc + mfStartaba);
					//mfhcca += rho*c1o2 * (mfStartcca + mfStartaac);
					//mfhccb += rho*c1o2 * (mfStartccb + mfStartaab);
					//mfhccc += rho*c1o2 * (mfStartccc + mfStartaaa);
					//mfhbcb += c1o6 * c2o9 * deltaPP;
					//mfhbbc += c1o6 * c2o9 * deltaPP;
					//mfhcbb += c1o6 * c2o9 * deltaPP;
					//mfhccb += c1o6 * c1o18 * deltaPP;
					//mfhacb += c1o6 * c1o18 * deltaPP;
					//mfhcbc += c1o6 * c1o18 * deltaPP;
					//mfhabc += c1o6 * c1o18 * deltaPP;
					//mfhbcc += c1o6 * c1o18 * deltaPP;
					//mfhbac += c1o6 * c1o18 * deltaPP;
					//mfhccc += c1o6 * c1o72 * deltaPP;
					//mfhacc += c1o6 * c1o72 * deltaPP;
					//mfhcac += c1o6 * c1o72 * deltaPP;
					//mfhaac += c1o6 * c1o72 * deltaPP;
					//mfhabb += c1o6 * c2o9 * deltaPP;
					//mfhbab += c1o6 * c2o9 * deltaPP;
					//mfhbba += c1o6 * c2o9 * deltaPP;
					//mfhaab += c1o6 * c1o18 * deltaPP;
					//mfhcab += c1o6 * c1o18 * deltaPP;
					//mfhaba += c1o6 * c1o18 * deltaPP;
					//mfhcba += c1o6 * c1o18 * deltaPP;
					//mfhbaa += c1o6 * c1o18 * deltaPP;
					//mfhbca += c1o6 * c1o18 * deltaPP;
					//mfhaaa += c1o6 * c1o72 * deltaPP;
					//mfhcaa += c1o6 * c1o72 * deltaPP;
					//mfhaca += c1o6 * c1o72 * deltaPP;
					//mfhcca += c1o6 * c1o72 * deltaPP;
					//mfhbbb += c1o6 * c4 * c2o9 * deltaPP;


					//mfhbcb = c1o3/rho * c2o9 ;
					//mfhbbc = c1o3/rho * c2o9 ;
					//mfhcbb = c1o3/rho * c2o9 ;
					//mfhccb = c1o3/rho * c1o18 ;
					//mfhacb = c1o3/rho * c1o18 ;
					//mfhcbc = c1o3/rho * c1o18 ;
					//mfhabc = c1o3/rho * c1o18 ;
					//mfhbcc = c1o3/rho * c1o18 ;
					//mfhbac = c1o3/rho * c1o18 ;
					//mfhccc = c1o3/rho * c1o72 ;
					//mfhacc = c1o3/rho * c1o72 ;
					//mfhcac = c1o3/rho * c1o72 ;
					//mfhaac = c1o3/rho * c1o72 ;
					//mfhabb = c1o3/rho * c2o9 ;
					//mfhbab = c1o3/rho * c2o9 ;
					//mfhbba = c1o3/rho * c2o9 ;
					//mfhaab = c1o3/rho * c1o18 ;
					//mfhcab = c1o3/rho * c1o18 ;
					//mfhaba = c1o3/rho * c1o18 ;
					//mfhcba = c1o3/rho * c1o18 ;
					//mfhbaa = c1o3/rho * c1o18 ;
					//mfhbca = c1o3/rho * c1o18 ;
					//mfhaaa = c1o3/rho * c1o72 ;
					//mfhcaa = c1o3/rho * c1o72 ;
					//mfhaca = c1o3/rho * c1o72 ;
					//mfhcca = c1o3/rho * c1o72 ;
					//mfhbbb = c1/rho;//c1o3/rho * c4 * c2o9 ;


					
					//mfabb += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_M00]) * rhoToPhi / rho;
					//mfbab += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_0M0]) * rhoToPhi / rho;
					//mfbba += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_00M]) * rhoToPhi / rho;
					//mfaab += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_MM0]) * rhoToPhi / rho;
					//mfcab += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_PM0]) * rhoToPhi / rho;
					//mfaba += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_M0M]) * rhoToPhi / rho;
					//mfcba += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_P0M]) * rhoToPhi / rho;
					//mfbaa += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_0MM]) * rhoToPhi / rho;
					//mfbca += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_0PM]) * rhoToPhi / rho;
					//mfaaa += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_MMM]) * rhoToPhi / rho;
					//mfcaa += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_PMM]) * rhoToPhi / rho;
					//mfaca += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_MPM]) * rhoToPhi / rho;
					//mfcca += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_PPM]) * rhoToPhi / rho;
					//mfcbb += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_P00]) * rhoToPhi / rho;
					//mfbcb += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_0P0]) * rhoToPhi / rho;
					//mfbbc += c1o2 * c2o9 * pStar * (phi[DIR_000] - phi[DIR_00P]) * rhoToPhi / rho;
					//mfccb += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_PP0]) * rhoToPhi / rho;
					//mfacb += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_MP0]) * rhoToPhi / rho;
					//mfcbc += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_P0P]) * rhoToPhi / rho;
					//mfabc += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_M0P]) * rhoToPhi / rho;
					//mfbcc += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_0PP]) * rhoToPhi / rho;
					//mfbac += c1o2 * c1o18 * pStar * (phi[DIR_000] - phi[DIR_0MP]) * rhoToPhi / rho;
					//mfccc += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_PPP]) * rhoToPhi / rho;
					//mfacc += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_MPP]) * rhoToPhi / rho;
					//mfcac += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_PMP]) * rhoToPhi / rho;
					//mfaac += c1o2 * c1o72 * pStar * (phi[DIR_000] - phi[DIR_MMP]) * rhoToPhi / rho;
					
					///////////////
					//mfabb += (pBefore-pStar) * c2o9  ;
					//mfbab += (pBefore-pStar) * c2o9  ;
					//mfbba += (pBefore-pStar) * c2o9  ;
					//mfaab += (pBefore-pStar) * c1o18 ;
					//mfcab += (pBefore-pStar) * c1o18 ;
					//mfaba += (pBefore-pStar) * c1o18 ;
					//mfcba += (pBefore-pStar) * c1o18 ;
					//mfbaa += (pBefore-pStar) * c1o18 ;
					//mfbca += (pBefore-pStar) * c1o18 ;
					//mfaaa += (pBefore-pStar) * c1o72 ;
					//mfcaa += (pBefore-pStar) * c1o72 ;
					//mfaca += (pBefore-pStar) * c1o72 ;
					//mfcca += (pBefore-pStar) * c1o72 ;
					//mfcbb += (pBefore-pStar) * c2o9  ;
					//mfbcb += (pBefore-pStar) * c2o9  ;
					//mfbbc += (pBefore-pStar) * c2o9  ;
					//mfccb += (pBefore-pStar) * c1o18 ;
					//mfacb += (pBefore-pStar) * c1o18 ;
					//mfcbc += (pBefore-pStar) * c1o18 ;
					//mfabc += (pBefore-pStar) * c1o18 ;
					//mfbcc += (pBefore-pStar) * c1o18 ;
					//mfbac += (pBefore-pStar) * c1o18 ;
					//mfccc += (pBefore-pStar) * c1o72 ;
					//mfacc += (pBefore-pStar) * c1o72 ;
					//mfcac += (pBefore-pStar) * c1o72 ;
					//mfaac += (pBefore-pStar) * c1o72 ;
					//mfbbb += (pBefore-pStar) * 8.0 / 9.0;

					//mfabb = (pBefore ) * c2o9;
					//mfbab = (pBefore ) * c2o9;
					//mfbba = (pBefore ) * c2o9;
					//mfaab = (pBefore ) * c1o16;
					//mfcab = (pBefore ) * c1o16;
					//mfaba = (pBefore ) * c1o16;
					//mfcba = (pBefore ) * c1o16;
					//mfbaa = (pBefore ) * c1o16;
					//mfbca = (pBefore ) * c1o16;
					//mfaaa = (pBefore ) * c1o72;
					//mfcaa = (pBefore ) * c1o72;
					//mfaca = (pBefore ) * c1o72;
					//mfcca = (pBefore ) * c1o72;
					//mfcbb = (pBefore ) * c2o9;
					//mfbcb = (pBefore ) * c2o9;
					//mfbbc = (pBefore ) * c2o9;
					//mfccb = (pBefore ) * c1o16;
					//mfacb = (pBefore ) * c1o16;
					//mfcbc = (pBefore ) * c1o16;
					//mfabc = (pBefore ) * c1o16;
					//mfbcc = (pBefore ) * c1o16;
					//mfbac = (pBefore ) * c1o16;
					//mfccc = (pBefore ) * c1o72;
					//mfacc = (pBefore ) * c1o72;
					//mfcac = (pBefore ) * c1o72;
					//mfaac = (pBefore ) * c1o72;
					//mfbbb = (pBefore ) * 8.0 / 9.0;
					///////////////////

					//////////////////////////////////////////////////////////////////////////
					//proof correctness
					//////////////////////////////////////////////////////////////////////////
					//#ifdef  PROOF_CORRECTNESS
					real rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
						+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
						+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
					//			   //LBMReal dif = fabs(drho - rho_post);
					//               LBMReal dif = drho + (dX1_phi * vvx + dX2_phi * vvy + dX3_phi * vvz) * correctionScaling - rho_post;
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
						UB_THROW(UbException(
							UB_EXARGS, "rho_post is not a number (nan or -1.#IND) or infinity number -1.#INF, node=" + UbSystem::toString(x1) + "," +
							UbSystem::toString(x2) + "," + UbSystem::toString(x3)));

					//////////////////////////////////////////////////////////////////////////
					//write distribution
					//////////////////////////////////////////////////////////////////////////
					(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3)         = mfabb         ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3)         = mfbab         ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3)         = mfbba         ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3)        = mfaab        ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3)       = mfcab       ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3)        = mfaba        ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3)       = mfcba       ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3)        = mfbaa        ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3)       = mfbca       ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3)       = mfaaa       ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3)      = mfcaa      ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3)      = mfaca      ;//* rho * c1o3;
					(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3)     = mfcca     ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3)     = mfcbb     ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3)     = mfbcb     ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p)     = mfbbc     ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfccb   ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3)    = mfacb    ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfcbc   ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p)    = mfabc    ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfbcc   ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p)    = mfbac    ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfacc  ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfcac  ;//* rho * c1o3;
					(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfaac   ;//* rho * c1o3;

					(*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;// *rho* c1o3;

			
					(*this->localDistributionsH2)(D3Q27System::ET_E, x1, x2, x3)         = mfhabb;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_N, x1, x2, x3)         = mfhbab;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_T, x1, x2, x3)         = mfhbba;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_NE, x1, x2, x3)        = mfhaab;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_NW, x1p, x2, x3)       = mfhcab;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TE, x1, x2, x3)        = mfhaba;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TW, x1p, x2, x3)       = mfhcba;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TN, x1, x2, x3)        = mfhbaa;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TS, x1, x2p, x3)       = mfhbca;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TNE, x1, x2, x3)       = mfhaaa;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TNW, x1p, x2, x3)      = mfhcaa;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TSE, x1, x2p, x3)      = mfhaca;//* rho * c1o3;
					(*this->localDistributionsH2)(D3Q27System::ET_TSW, x1p, x2p, x3)     = mfhcca;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_W, x1p, x2, x3)     = mfhcbb;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_S, x1, x2p, x3)     = mfhbcb;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_B, x1, x2, x3p)     = mfhbbc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_SW, x1p, x2p, x3)   = mfhccb;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_SE, x1, x2p, x3)    = mfhacb;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BW, x1p, x2, x3p)   = mfhcbc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BE, x1, x2, x3p)    = mfhabc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BS, x1, x2p, x3p)   = mfhbcc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BN, x1, x2, x3p)    = mfhbac;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfhccc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BSE, x1, x2p, x3p)  = mfhacc;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNW, x1p, x2, x3p)  = mfhcac;//* rho * c1o3;
					(*this->nonLocalDistributionsH2)(D3Q27System::ET_BNE, x1, x2, x3p)   = mfhaac;//* rho * c1o3;

					(*this->zeroDistributionsH2)(x1, x2, x3) = mfhbbb;// *rho* c1o3;

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
						real oneMinusRho = c1o1 - concentration;

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

						//31.05.2022 addaptive mobility
						//omegaD = c1 + (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration))) / (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration)) + fabs((1.0 - concentration) * (concentration)) * c1o6 * oneOverInterfaceScale+1.0e-200);
						//omegaD = c2 * (concentration * (concentration - c1)) / (-c6 * (sqrt((cx - vvx * concentration) * (cx - vvx * concentration) + (cy - vvy * concentration) * (cy - vvy * concentration) + (cz - vvz * concentration) * (cz - vvz * concentration))) + (concentration * (concentration - c1))+1.0e-200);
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

						// equilibration of 2nd order moments
						mfbba = c0o1;
						mfbab = c0o1;
						mfabb = c0o1;

						mfcaa = c1o3 * concentration;
						mfaca = c1o3 * concentration;
						mfaac = c1o3 * concentration;

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
}
//////////////////////////////////////////////////////////////////////////

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX1_phi()
{
	using namespace vf::lbm::dir;
	using namespace D3Q27System;

	return 3.0* ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) + (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) + (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) + (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_PM0] - phi[DIR_MP0]) + (phi[DIR_PP0] - phi[DIR_MM0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_P00] - phi[DIR_M00]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX2_phi()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PPM] - phi[DIR_MMP])- (phi[DIR_PMP] - phi[DIR_MPM])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_0PP] - phi[DIR_0MM]) + (phi[DIR_0PM] - phi[DIR_0MP])) + ((phi[DIR_PP0] - phi[DIR_MM0])- (phi[DIR_PM0] - phi[DIR_MP0])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_0P0] - phi[DIR_0M0]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX3_phi()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi[DIR_PPP] - phi[DIR_MMM]) - (phi[DIR_PMM] - phi[DIR_MPP])) + ((phi[DIR_PMP] - phi[DIR_MPM]) - (phi[DIR_PPM] - phi[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi[DIR_P0P] - phi[DIR_M0M]) - (phi[DIR_P0M] - phi[DIR_M0P])) + ((phi[DIR_0MP] - phi[DIR_0PM]) + (phi[DIR_0PP] - phi[DIR_0MM])))) +
		+WEIGTH[DIR_0P0] * (phi[DIR_00P] - phi[DIR_00M]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX1_rhoInv(real rhoL,real rhoDIV)
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0/(rhoL+rhoDIV*phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_P0P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0M])) + (1.0 / (rhoL + rhoDIV * phi[DIR_P0M]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0P]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PM0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MP0])) + (1.0 / (rhoL + rhoDIV * phi[DIR_PP0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MM0]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_P00]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M00])));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX2_rhoInv(real rhoL,real rhoDIV)
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0 / (rhoL + rhoDIV * phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_0PP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_0PM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PP0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MM0])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PM0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MP0]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_0P0]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0M0])));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX3_rhoInv(real rhoL, real rhoDIV)
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((1.0 / (rhoL + rhoDIV * phi[DIR_PPP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PMM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPP]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_PMP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MPM])) - (1.0 / (rhoL + rhoDIV * phi[DIR_PPM]) - 1.0 / (rhoL + rhoDIV * phi[DIR_MMP]))))
		+ WEIGTH[DIR_PP0] * (((1.0 / (rhoL + rhoDIV * phi[DIR_P0P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0M])) - (1.0 / (rhoL + rhoDIV * phi[DIR_P0M]) - 1.0 / (rhoL + rhoDIV * phi[DIR_M0P]))) + ((1.0 / (rhoL + rhoDIV * phi[DIR_0MP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0PM])) + (1.0 / (rhoL + rhoDIV * phi[DIR_0PP]) - 1.0 / (rhoL + rhoDIV * phi[DIR_0MM]))))) +
		+WEIGTH[DIR_0P0] * (1.0 / (rhoL + rhoDIV * phi[DIR_00P]) - 1.0 / (rhoL + rhoDIV * phi[DIR_00M])));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX1_phi2()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) + (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) + (phi2[DIR_PPM] - phi2[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) + (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_PM0] - phi2[DIR_MP0]) + (phi2[DIR_PP0] - phi2[DIR_MM0])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_P00] - phi2[DIR_M00]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX2_phi2()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PPM] - phi2[DIR_MMP]) - (phi2[DIR_PMP] - phi2[DIR_MPM])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_0PP] - phi2[DIR_0MM]) + (phi2[DIR_0PM] - phi2[DIR_0MP])) + ((phi2[DIR_PP0] - phi2[DIR_MM0]) - (phi2[DIR_PM0] - phi2[DIR_MP0])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_0P0] - phi2[DIR_0M0]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::gradX3_phi2()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	return 3.0 * ((WEIGTH[DIR_PPP] * (((phi2[DIR_PPP] - phi2[DIR_MMM]) - (phi2[DIR_PMM] - phi2[DIR_MPP])) + ((phi2[DIR_PMP] - phi2[DIR_MPM]) - (phi2[DIR_PPM] - phi2[DIR_MMP])))
		+ WEIGTH[DIR_PP0] * (((phi2[DIR_P0P] - phi2[DIR_M0M]) - (phi2[DIR_P0M] - phi2[DIR_M0P])) + ((phi2[DIR_0MP] - phi2[DIR_0PM]) + (phi2[DIR_0PP] - phi2[DIR_0MM])))) +
		+WEIGTH[DIR_0P0] * (phi2[DIR_00P] - phi2[DIR_00M]));
}

real MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::nabla2_phi()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

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

void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::computePhasefield()
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

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

void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::findNeighbors(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
	int x3)
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;
	using namespace vf::basics::constant;

	SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

	phi[DIR_000] = (*ph)(x1, x2, x3);
    if (phi[DIR_000] < 0) {
        phi[DIR_000] = c0o1;
    }


	for (int k = FSTARTDIR; k <= FENDDIR; k++) {

		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
			phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
		} else {
			//phi[k] = (*ph)(x1 , x2, x3 );// neutral wetting
			phi[k] = 0.0;//unwetting
		}
	}
}

void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::findNeighbors2(CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2,
	int x3)
{
	using namespace D3Q27System;
	using namespace vf::lbm::dir;

	SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

	phi2[DIR_000] = (*ph)(x1, x2, x3);


	for (int k = FSTARTDIR; k <= FENDDIR; k++) {

		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])) {
			phi2[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
		}
		else {
			phi2[k] = 0.05;
		}
	}
}

void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::swapDistributions()
{
	LBMKernel::swapDistributions();
	dataSet->getHdistributions()->swap();
	dataSet->getH2distributions()->swap();
}

void MultiphaseSimpleVelocityBaseExternalPressureLBMKernel::initForcing()
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
