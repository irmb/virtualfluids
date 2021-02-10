#include "ThixotropyExpLBMKernel.h"
#include "D3Q27System.h"
#include "InterpolationProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <math.h>
#include "DataSet3D.h"
#include "LBMKernel.h"

#define PROOF_CORRECTNESS

using namespace UbMath;

//////////////////////////////////////////////////////////////////////////
ThixotropyExpLBMKernel::ThixotropyExpLBMKernel()
{
   this->parameter = ThixotropyExpLBMKernel::NORMAL;
	this->compressible = false;
	//this->TwoDistributions = true;
}
//////////////////////////////////////////////////////////////////////////
ThixotropyExpLBMKernel::~ThixotropyExpLBMKernel(void)
{

}
//////////////////////////////////////////////////////////////////////////
void ThixotropyExpLBMKernel::initDataSet()
{
	//DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+ghostLayerWitdh*2, nx2+ghostLayerWitdh*2, nx3+ghostLayerWitdh*2, -999.0));
	SPtr<DistributionArray3D> df(new D3Q27EsoTwist3DSplittedVector(nx[0]+2, nx[1]+2, nx[2]+2, -999.0));
	SPtr<DistributionArray3D> dh(new D3Q27EsoTwist3DSplittedVector(nx[0]+2, nx[1]+2, nx[2]+2, -999.0));
	dataSet->setFdistributions(df);
	dataSet->setHdistributions(dh);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> ThixotropyExpLBMKernel::clone()
{
	SPtr<LBMKernel> kernel(new ThixotropyExpLBMKernel());
	kernel->setNX(nx);
	kernel->setCollisionFactor(collFactor);
	collFactorF = collFactor;
	collFactorH = collFactor;
	dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->initDataSet();
	dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->setCollisionFactorF(this->collFactorF);
	dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->setCollisionFactorH(this->collFactorH);
	dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->setAlpha(this->alpha);
	dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->setTheta(this->theta);
	kernel->setBCProcessor(bcProcessor->clone(kernel));
	kernel->setWithForcing(withForcing);
	kernel->setForcingX1(muForcingX1);
	kernel->setForcingX2(muForcingX2);
	kernel->setForcingX3(muForcingX3);
	kernel->setIndex(ix1, ix2, ix3);
	kernel->setDeltaT(deltaT);

	switch (parameter)
	{
	case NORMAL:
		dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->OxyyMxzz = 1.0;
		break;
	case MAGIC:
		dynamicPointerCast<ThixotropyExpLBMKernel>(kernel)->OxyyMxzz = 2.0 + (-collFactorF);
		break;
	}
	return kernel;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyExpLBMKernel::calculate(int step)
{
	using namespace D3Q27System;

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

		muNu = (1.0 / 3.0)*(1.0 / collFactor - 1.0 / 2.0);

		muForcingX1.DefineVar("nu", &muNu);
		muForcingX2.DefineVar("nu", &muNu);
		muForcingX3.DefineVar("nu", &muNu);
	}
	/////////////////////////////////////

	localDistributionsF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
	nonLocalDistributionsF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
	zeroDistributionsF = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

	localDistributionsH = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getLocalDistributions();
	nonLocalDistributionsH = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getNonLocalDistributions();
	zeroDistributionsH = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getZeroDistributions();

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


	//#pragma omp parallel num_threads(8)
	{
		//   int i = omp_get_thread_num();
		//   printf_s("Hello from thread %d\n", i);
		//}
		//#pragma omp for 
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
						//read distribution
						// Cumulant (NSE part) 
						////////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////

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

						LBMReal lambda = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
							+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
							+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

						//E   N  T
						//c   c  c
						//////////
						//W   S  B
						//a   a  a

						//Rest ist b

						//mfxyz
						//a - negative
						//b - null
						//c - positive

						// a b c
						//-1 0 1

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

						LBMReal m0, m1, m2;

						LBMReal rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
							+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
							+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

						LBMReal vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
							(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
							(mfcbb - mfabb));
						LBMReal vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
							(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
							(mfbcb - mfbab));
						LBMReal vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
							(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
							(mfbbc - mfbba));
						

//						LBMReal eta0 = (1/collFactor-c1o2)*c1o3;
//						LBMReal eta = (1 + lambda)* eta0;
						//collFactorF = one/(3*eta/(rho+one)+c1o2);
						collFactorF = collFactor;

						//forcing 
						///////////////////////////////////////////////////////////////////////////////////////////
						if (withForcing)
						{
							muX1 = static_cast<double>(x1 - 1 + ix1*maxX1);
							muX2 = static_cast<double>(x2 - 1 + ix2*maxX2);
							muX3 = static_cast<double>(x3 - 1 + ix3*maxX3);

							forcingX1 = muForcingX1.Eval();
							forcingX2 = muForcingX2.Eval();
							forcingX3 = muForcingX3.Eval();

							vvx += forcingX1*deltaT*0.5; // X
							vvy += forcingX2*deltaT*0.5; // Y
							vvz += forcingX3*deltaT*0.5; // Z
						}
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
						oMdrho = 1. - (oMdrho + m0);

						LBMReal vx2;
						LBMReal vy2;
						LBMReal vz2;
						vx2 = vvx*vvx;
						vy2 = vvy*vvy;
						vz2 = vvz*vvz;
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
						mfaac = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaba + mfabc;
						m1 = mfabc - mfaba;
						m0 = m2 + mfabb;
						mfaba = m0;
						m0 += c1o9 * oMdrho;
						mfabb = m1 - m0 * vvz;
						mfabc = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaca + mfacc;
						m1 = mfacc - mfaca;
						m0 = m2 + mfacb;
						mfaca = m0;
						m0 += c1o36 * oMdrho;
						mfacb = m1 - m0 * vvz;
						mfacc = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbaa + mfbac;
						m1 = mfbac - mfbaa;
						m0 = m2 + mfbab;
						mfbaa = m0;
						m0 += c1o9 * oMdrho;
						mfbab = m1 - m0 * vvz;
						mfbac = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbba + mfbbc;
						m1 = mfbbc - mfbba;
						m0 = m2 + mfbbb;
						mfbba = m0;
						m0 += c4o9 * oMdrho;
						mfbbb = m1 - m0 * vvz;
						mfbbc = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbca + mfbcc;
						m1 = mfbcc - mfbca;
						m0 = m2 + mfbcb;
						mfbca = m0;
						m0 += c1o9 * oMdrho;
						mfbcb = m1 - m0 * vvz;
						mfbcc = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcaa + mfcac;
						m1 = mfcac - mfcaa;
						m0 = m2 + mfcab;
						mfcaa = m0;
						m0 += c1o36 * oMdrho;
						mfcab = m1 - m0 * vvz;
						mfcac = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcba + mfcbc;
						m1 = mfcbc - mfcba;
						m0 = m2 + mfcbb;
						mfcba = m0;
						m0 += c1o9 * oMdrho;
						mfcbb = m1 - m0 * vvz;
						mfcbc = m2 - 2. *   m1 * vvz + vz2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcca + mfccc;
						m1 = mfccc - mfcca;
						m0 = m2 + mfccb;
						mfcca = m0;
						m0 += c1o36 * oMdrho;
						mfccb = m1 - m0 * vvz;
						mfccc = m2 - 2. *   m1 * vvz + vz2 * m0;
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
						mfaca = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaab + mfacb;
						m1 = mfacb - mfaab;
						m0 = m2 + mfabb;
						mfaab = m0;
						mfabb = m1 - m0 * vvy;
						mfacb = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaac + mfacc;
						m1 = mfacc - mfaac;
						m0 = m2 + mfabc;
						mfaac = m0;
						m0 += c1o18 * oMdrho;
						mfabc = m1 - m0 * vvy;
						mfacc = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbaa + mfbca;
						m1 = mfbca - mfbaa;
						m0 = m2 + mfbba;
						mfbaa = m0;
						m0 += c2o3 * oMdrho;
						mfbba = m1 - m0 * vvy;
						mfbca = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbab + mfbcb;
						m1 = mfbcb - mfbab;
						m0 = m2 + mfbbb;
						mfbab = m0;
						mfbbb = m1 - m0 * vvy;
						mfbcb = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfbac + mfbcc;
						m1 = mfbcc - mfbac;
						m0 = m2 + mfbbc;
						mfbac = m0;
						m0 += c2o9 * oMdrho;
						mfbbc = m1 - m0 * vvy;
						mfbcc = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcaa + mfcca;
						m1 = mfcca - mfcaa;
						m0 = m2 + mfcba;
						mfcaa = m0;
						m0 += c1o6 * oMdrho;
						mfcba = m1 - m0 * vvy;
						mfcca = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcab + mfccb;
						m1 = mfccb - mfcab;
						m0 = m2 + mfcbb;
						mfcab = m0;
						mfcbb = m1 - m0 * vvy;
						mfccb = m2 - 2. *   m1 * vvy + vy2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfcac + mfccc;
						m1 = mfccc - mfcac;
						m0 = m2 + mfcbc;
						mfcac = m0;
						m0 += c1o18 * oMdrho;
						mfcbc = m1 - m0 * vvy;
						mfccc = m2 - 2. *   m1 * vvy + vy2 * m0;
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
						mfcaa = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaba + mfcba;
						m1 = mfcba - mfaba;
						m0 = m2 + mfbba;
						mfaba = m0;
						mfbba = m1 - m0 * vvx;
						mfcba = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaca + mfcca;
						m1 = mfcca - mfaca;
						m0 = m2 + mfbca;
						mfaca = m0;
						m0 += c1o3 * oMdrho;
						mfbca = m1 - m0 * vvx;
						mfcca = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaab + mfcab;
						m1 = mfcab - mfaab;
						m0 = m2 + mfbab;
						mfaab = m0;
						mfbab = m1 - m0 * vvx;
						mfcab = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfabb + mfcbb;
						m1 = mfcbb - mfabb;
						m0 = m2 + mfbbb;
						mfabb = m0;
						mfbbb = m1 - m0 * vvx;
						mfcbb = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfacb + mfccb;
						m1 = mfccb - mfacb;
						m0 = m2 + mfbcb;
						mfacb = m0;
						mfbcb = m1 - m0 * vvx;
						mfccb = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfaac + mfcac;
						m1 = mfcac - mfaac;
						m0 = m2 + mfbac;
						mfaac = m0;
						m0 += c1o3 * oMdrho;
						mfbac = m1 - m0 * vvx;
						mfcac = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfabc + mfcbc;
						m1 = mfcbc - mfabc;
						m0 = m2 + mfbbc;
						mfabc = m0;
						mfbbc = m1 - m0 * vvx;
						mfcbc = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						m2 = mfacc + mfccc;
						m1 = mfccc - mfacc;
						m0 = m2 + mfbcc;
						mfacc = m0;
						m0 += c1o9 * oMdrho;
						mfbcc = m1 - m0 * vvx;
						mfccc = m2 - 2. *   m1 * vvx + vx2 * m0;
						////////////////////////////////////////////////////////////////////////////////////
						// Cumulants
						////////////////////////////////////////////////////////////////////////////////////
						LBMReal OxxPyyPzz = 1.; //omega2 or bulk viscosity
						LBMReal OxyyPxzz = 1.;//-s9;//2+s9;//
											  //LBMReal OxyyMxzz  = 1.;//2+s9;//
						LBMReal O4 = 1.;
						LBMReal O5 = 1.;
						LBMReal O6 = 1.;

						//Cum 4.
						//LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
						//LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
						//LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

						LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
						LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
						LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);

						LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho - 1)*oMdrho);
						LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho - 1)*oMdrho);
						LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho - 1)*oMdrho);

						//Cum 5.
						LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
						LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
						LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

						//Cum 6.
						LBMReal CUMccc = mfccc + ((-4. *  mfbbb * mfbbb
							- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							- 4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							- 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+ (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
								+ 2. * (mfcaa * mfaca * mfaac)
								+ 16. *  mfbba * mfbab * mfabb)
							- c1o3* (mfacc + mfcac + mfcca) * oMdrho - c1o9*oMdrho*oMdrho
							- c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(1. - 2.* oMdrho) - c1o27* oMdrho * oMdrho*(-2.* oMdrho)
							+ (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
								+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) + c1o27*oMdrho;

						//2.
						// linear combinations
						LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
						LBMReal mxxMyy = mfcaa - mfaca;
						LBMReal mxxMzz = mfcaa - mfaac;

						LBMReal dxux = -c1o2 * collFactorF *(mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz*(mfaaa - mxxPyyPzz);
						LBMReal dyuy = dxux + collFactorF * c3o2 * mxxMyy;
						LBMReal dzuz = dxux + collFactorF * c3o2 * mxxMzz;

						LBMReal Dxy =-three*collFactorF*mfbba;
                  LBMReal Dxz =-three*collFactorF*mfbab;
                  LBMReal Dyz =-three*collFactorF*mfabb;

						LBMReal gammaDot = sqrt(dxux * dxux + dyuy * dyuy + dzuz * dzuz + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz) / (rho + one);
						//collFactorF = BinghamModel::getBinghamCollFactor(collFactorF, gammaDot, rho);

						//relax
						mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
						mxxMyy += collFactorF * (-mxxMyy) - 3. * (1. - c1o2 * collFactorF) * (vx2 * dxux - vy2 * dyuy);
						mxxMzz += collFactorF * (-mxxMzz) - 3. * (1. - c1o2 * collFactorF) * (vx2 * dxux - vz2 * dzuz);

						mfabb += collFactorF * (-mfabb);
						mfbab += collFactorF * (-mfbab);
						mfbba += collFactorF * (-mfbba);

						// linear combinations back
						mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
						mfaca = c1o3 * (-2. *  mxxMyy + mxxMzz + mxxPyyPzz);
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
						wadjust = OxyyMxzz + (1. - OxyyMxzz)*fabs(mfbbb) / (fabs(mfbbb) + qudricLimit);
						mfbbb += wadjust * (-mfbbb);
						wadjust = OxyyPxzz + (1. - OxyyPxzz)*fabs(mxxyPyzz) / (fabs(mxxyPyzz) + qudricLimit);
						mxxyPyzz += wadjust * (-mxxyPyzz);
						wadjust = OxyyMxzz + (1. - OxyyMxzz)*fabs(mxxyMyzz) / (fabs(mxxyMyzz) + qudricLimit);
						mxxyMyzz += wadjust * (-mxxyMyzz);
						wadjust = OxyyPxzz + (1. - OxyyPxzz)*fabs(mxxzPyyz) / (fabs(mxxzPyyz) + qudricLimit);
						mxxzPyyz += wadjust * (-mxxzPyyz);
						wadjust = OxyyMxzz + (1. - OxyyMxzz)*fabs(mxxzMyyz) / (fabs(mxxzMyyz) + qudricLimit);
						mxxzMyyz += wadjust * (-mxxzMyyz);
						wadjust = OxyyPxzz + (1. - OxyyPxzz)*fabs(mxyyPxzz) / (fabs(mxyyPxzz) + qudricLimit);
						mxyyPxzz += wadjust * (-mxyyPxzz);
						wadjust = OxyyMxzz + (1. - OxyyMxzz)*fabs(mxyyMxzz) / (fabs(mxyyMxzz) + qudricLimit);
						mxyyMxzz += wadjust * (-mxyyMxzz);

						// linear combinations back
						mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
						mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
						mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
						mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
						mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
						mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

						//4.
						CUMacc += O4 * (-CUMacc);
						CUMcac += O4 * (-CUMcac);
						CUMcca += O4 * (-CUMcca);

						CUMbbc += O4 * (-CUMbbc);
						CUMbcb += O4 * (-CUMbcb);
						CUMcbb += O4 * (-CUMcbb);

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

						mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho - 1)*oMdrho;
						mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho - 1)*oMdrho;
						mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho - 1)*oMdrho;

						//5.
						mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
						mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
						mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

						//6.
						mfccc = CUMccc - ((-4. *  mfbbb * mfbbb
							- (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							- 4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
							- 2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+ (4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
								+ 2. * (mfcaa * mfaca * mfaac)
								+ 16. *  mfbba * mfbab * mfabb)
							- c1o3* (mfacc + mfcac + mfcca) * oMdrho - c1o9*oMdrho*oMdrho
							- c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(1. - 2.* oMdrho) - c1o27* oMdrho * oMdrho*(-2.* oMdrho)
							+ (2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
								+ (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) - c1o27*oMdrho;

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
						m1 = -mfaac - 2. * mfaab *  vvz + mfaaa                * (1. - vz2) - 1. * oMdrho * vz2;
						m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + 1. * oMdrho) * (vz2 + vvz) * c1o2;
						mfaaa = m0;
						mfaab = m1;
						mfaac = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
						m1 = -mfabc - 2. * mfabb *  vvz + mfaba * (1. - vz2);
						m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
						mfaba = m0;
						mfabb = m1;
						mfabc = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
						m1 = -mfacc - 2. * mfacb *  vvz + mfaca                  * (1. - vz2) - c1o3 * oMdrho * vz2;
						m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
						mfaca = m0;
						mfacb = m1;
						mfacc = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
						m1 = -mfbac - 2. * mfbab *  vvz + mfbaa * (1. - vz2);
						m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
						mfbaa = m0;
						mfbab = m1;
						mfbac = m2;
						/////////b//////////////////////////////////////////////////////////////////////////
						m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
						m1 = -mfbbc - 2. * mfbbb *  vvz + mfbba * (1. - vz2);
						m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
						mfbba = m0;
						mfbbb = m1;
						mfbbc = m2;
						/////////b//////////////////////////////////////////////////////////////////////////
						m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
						m1 = -mfbcc - 2. * mfbcb *  vvz + mfbca * (1. - vz2);
						m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
						mfbca = m0;
						mfbcb = m1;
						mfbcc = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
						m1 = -mfcac - 2. * mfcab *  vvz + mfcaa                  * (1. - vz2) - c1o3 * oMdrho * vz2;
						m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
						mfcaa = m0;
						mfcab = m1;
						mfcac = m2;
						/////////c//////////////////////////////////////////////////////////////////////////
						m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
						m1 = -mfcbc - 2. * mfcbb *  vvz + mfcba * (1. - vz2);
						m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
						mfcba = m0;
						mfcbb = m1;
						mfcbc = m2;
						/////////c//////////////////////////////////////////////////////////////////////////
						m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
						m1 = -mfccc - 2. * mfccb *  vvz + mfcca                  * (1. - vz2) - c1o9 * oMdrho * vz2;
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
						m1 = -mfaca - 2. * mfaba *  vvy + mfaaa                  * (1. - vy2) - c1o6 * oMdrho * vy2;
						m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
						mfaaa = m0;
						mfaba = m1;
						mfaca = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
						m1 = -mfacb - 2. * mfabb *  vvy + mfaab                  * (1. - vy2) - c2o3 * oMdrho * vy2;
						m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
						mfaab = m0;
						mfabb = m1;
						mfacb = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
						m1 = -mfacc - 2. * mfabc *  vvy + mfaac                  * (1. - vy2) - c1o6 * oMdrho * vy2;
						m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
						mfaac = m0;
						mfabc = m1;
						mfacc = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
						m1 = -mfbca - 2. * mfbba *  vvy + mfbaa * (1. - vy2);
						m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
						mfbaa = m0;
						mfbba = m1;
						mfbca = m2;
						/////////b//////////////////////////////////////////////////////////////////////////
						m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
						m1 = -mfbcb - 2. * mfbbb *  vvy + mfbab * (1. - vy2);
						m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
						mfbab = m0;
						mfbbb = m1;
						mfbcb = m2;
						/////////b//////////////////////////////////////////////////////////////////////////
						m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
						m1 = -mfbcc - 2. * mfbbc *  vvy + mfbac * (1. - vy2);
						m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
						mfbac = m0;
						mfbbc = m1;
						mfbcc = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
						m1 = -mfcca - 2. * mfcba *  vvy + mfcaa                   * (1. - vy2) - c1o18 * oMdrho * vy2;
						m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
						mfcaa = m0;
						mfcba = m1;
						mfcca = m2;
						/////////c//////////////////////////////////////////////////////////////////////////
						m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
						m1 = -mfccb - 2. * mfcbb *  vvy + mfcab                  * (1. - vy2) - c2o9 * oMdrho * vy2;
						m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
						mfcab = m0;
						mfcbb = m1;
						mfccb = m2;
						/////////c//////////////////////////////////////////////////////////////////////////
						m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
						m1 = -mfccc - 2. * mfcbc *  vvy + mfcac                   * (1. - vy2) - c1o18 * oMdrho * vy2;
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
						m1 = -mfcaa - 2. * mfbaa *  vvx + mfaaa                   * (1. - vx2) - c1o36 * oMdrho * vx2;
						m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
						mfaaa = m0;
						mfbaa = m1;
						mfcaa = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcba - 2. * mfbba *  vvx + mfaba                  * (1. - vx2) - c1o9 * oMdrho * vx2;
						m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
						mfaba = m0;
						mfbba = m1;
						mfcba = m2;
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcca - 2. * mfbca *  vvx + mfaca                   * (1. - vx2) - c1o36 * oMdrho * vx2;
						m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
						mfaca = m0;
						mfbca = m1;
						mfcca = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcab - 2. * mfbab *  vvx + mfaab                  * (1. - vx2) - c1o9 * oMdrho * vx2;
						m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
						mfaab = m0;
						mfbab = m1;
						mfcab = m2;
						///////////b////////////////////////////////////////////////////////////////////////
						m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcbb - 2. * mfbbb *  vvx + mfabb                  * (1. - vx2) - c4o9 * oMdrho * vx2;
						m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
						mfabb = m0;
						mfbbb = m1;
						mfcbb = m2;
						///////////b////////////////////////////////////////////////////////////////////////
						m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfccb - 2. * mfbcb *  vvx + mfacb                  * (1. - vx2) - c1o9 * oMdrho * vx2;
						m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
						mfacb = m0;
						mfbcb = m1;
						mfccb = m2;
						////////////////////////////////////////////////////////////////////////////////////
						////////////////////////////////////////////////////////////////////////////////////
						m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcac - 2. * mfbac *  vvx + mfaac                   * (1. - vx2) - c1o36 * oMdrho * vx2;
						m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
						mfaac = m0;
						mfbac = m1;
						mfcac = m2;
						///////////c////////////////////////////////////////////////////////////////////////
						m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfcbc - 2. * mfbbc *  vvx + mfabc                  * (1. - vx2) - c1o9 * oMdrho * vx2;
						m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
						mfabc = m0;
						mfbbc = m1;
						mfcbc = m2;
						///////////c////////////////////////////////////////////////////////////////////////
						m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
						m1 = -mfccc - 2. * mfbcc *  vvx + mfacc                   * (1. - vx2) - c1o36 * oMdrho * vx2;
						m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
						mfacc = m0;
						mfbcc = m1;
						mfccc = m2;

						//////////////////////////////////////////////////////////////////////////
						//proof correctness
						//////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
						LBMReal rho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
							+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
							+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
						//LBMReal dif = fabs(rho - rho_post);
						LBMReal dif = rho - rho_post;
#ifdef SINGLEPRECISION
						if (dif > 10.0E-7 || dif < -10.0E-7)
#else
						if (dif > 10.0E-15 || dif < -10.0E-15)
#endif					
						{
							UB_THROW(UbException(UB_EXARGS, "rho=" + UbSystem::toString(rho) + ", rho_post=" + UbSystem::toString(rho_post)
								+ " dif=" + UbSystem::toString(dif)
								+ " rho is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
							//UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
							//exit(EXIT_FAILURE);
						}
#endif
						//////////////////////////////////////////////////////////////////////////
						//write distribution
						//////////////////////////////////////////////////////////////////////////
						(*this->localDistributionsF)(D3Q27System::ET_E, x1, x2, x3) = mfabb;
						(*this->localDistributionsF)(D3Q27System::ET_N, x1, x2, x3) = mfbab;
						(*this->localDistributionsF)(D3Q27System::ET_T, x1, x2, x3) = mfbba;
						(*this->localDistributionsF)(D3Q27System::ET_NE, x1, x2, x3) = mfaab;
						(*this->localDistributionsF)(D3Q27System::ET_NW, x1p, x2, x3) = mfcab;
						(*this->localDistributionsF)(D3Q27System::ET_TE, x1, x2, x3) = mfaba;
						(*this->localDistributionsF)(D3Q27System::ET_TW, x1p, x2, x3) = mfcba;
						(*this->localDistributionsF)(D3Q27System::ET_TN, x1, x2, x3) = mfbaa;
						(*this->localDistributionsF)(D3Q27System::ET_TS, x1, x2p, x3) = mfbca;
						(*this->localDistributionsF)(D3Q27System::ET_TNE, x1, x2, x3) = mfaaa;
						(*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2, x3) = mfcaa;
						(*this->localDistributionsF)(D3Q27System::ET_TSE, x1, x2p, x3) = mfaca;
						(*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

						(*this->nonLocalDistributionsF)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac;
						(*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac;

						(*this->zeroDistributionsF)(x1, x2, x3) = mfbbb;
















						LBMReal ux, uy, uz;

						ux = vvx;						
						uy = vvy;
						uz = vvz;
							
							

				
						////////////////////////////////////////////////////////////////////////////
						// Central Factorized Moment (advection diffusion part) 
						//////////////////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////////////////
						//read distribution

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
						LBMReal drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
							(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
							((mfabb + mfcbb) + (mfbab + mfbcb)) + (mfbba + mfbbc)) + mfbbb;
					
						//flocculation
						//LBMReal lambda = drho;

						//LBMReal theta = 60 * 1.28172e+06;
						//LBMReal alpha = 0.005;// *10.0;

						//LBMReal gammaDot = sqrt(dxux * dxux + dyuy * dyuy + dzuz * dzuz + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz) / (rho + one);

						//collFactorF = collFactorF - 1e-4 / (gammaDot + one * 1e-9);

						//collFactorF = (collFactorF < 0.5) ? 0.5 : collFactorF;

						LBMReal dlambda = one / theta - alpha * lambda * gammaDot;


						//////////////////////////////////////////////////////////////////////////
						//collision Factorized Central moment Kernel Geier 2015
						//////////////////////////////////////////////////////////////////////////               
						LBMReal Mom000 = mfaaa + mfaab + mfaac + mfaba + mfabb + mfabc + mfaca +
							mfacb + mfacc + mfbaa + mfbab + mfbac + mfbba + mfbbb + mfbbc + mfbca +
							mfbcb + mfbcc + mfcaa + mfcab + mfcac + mfcba + mfcbb + mfcbc + mfcca + mfccb + mfccc;
						
						Mom000 += dlambda*c1o2;  //1

																												   //(100)//
						LBMReal Mom100 = mfaaa*(-1 - ux) + mfaab*(-1 - ux) + mfaac*(-1 - ux) + mfaba*(-1 - ux) + mfabb*(-1 - ux) + mfabc*(-1 - ux) + mfaca*(-1 - ux) +
							mfacb*(-1 - ux) + mfacc*(-1 - ux) + mfcaa*(1 - ux) + mfcab*(1 - ux) + mfcac*(1 - ux) + mfcba*(1 - ux) + mfcbb*(1 - ux) +
							mfcbc*(1 - ux) + mfcca*(1 - ux) + mfccb*(1 - ux) + mfccc*(1 - ux) - mfbaa*ux - mfbab*ux - mfbac*ux - mfbba*ux - mfbbb*ux -
							mfbbc*ux - mfbca*ux - mfbcb*ux - mfbcc*ux;

						LBMReal Mom010 = mfaaa*(-1 - uy) + mfaab*(-1 - uy) + mfaac*(-1 - uy) + mfbaa*(-1 - uy) + mfbab*(-1 - uy) + mfbac*(-1 - uy) + mfcaa*(-1 - uy) +
							mfcab*(-1 - uy) + mfcac*(-1 - uy) + mfaca*(1 - uy) + mfacb*(1 - uy) + mfacc*(1 - uy) + mfbca*(1 - uy) + mfbcb*(1 - uy) +
							mfbcc*(1 - uy) + mfcca*(1 - uy) + mfccb*(1 - uy) + mfccc*(1 - uy) - mfaba*uy - mfabb*uy - mfabc*uy - mfbba*uy - mfbbb*uy -
							mfbbc*uy - mfcba*uy - mfcbb*uy - mfcbc*uy;

						LBMReal Mom001 = mfaaa*(-1 - uz) + mfaba*(-1 - uz) + mfaca*(-1 - uz) + mfbaa*(-1 - uz) + mfbba*(-1 - uz) + mfbca*(-1 - uz) + mfcaa*(-1 - uz) +
							mfcba*(-1 - uz) + mfcca*(-1 - uz) + mfaac*(1 - uz) + mfabc*(1 - uz) + mfacc*(1 - uz) + mfbac*(1 - uz) + mfbbc*(1 - uz) +
							mfbcc*(1 - uz) + mfcac*(1 - uz) + mfcbc*(1 - uz) + mfccc*(1 - uz) - mfaab*uz - mfabb*uz - mfacb*uz - mfbab*uz - mfbbb*uz -
							mfbcb*uz - mfcab*uz - mfcbb*uz - mfccb*uz;
						////

						//(110)//
						LBMReal Mom110 = mfaaa*(-1 - ux)*(-1 - uy) + mfaab*(-1 - ux)*(-1 - uy) + mfaac*(-1 - ux)*(-1 - uy) + mfcaa*(1 - ux)*(-1 - uy) +
							mfcab*(1 - ux)*(-1 - uy) + mfcac*(1 - ux)*(-1 - uy) - mfbaa*ux*(-1 - uy) - mfbab*ux*(-1 - uy) - mfbac*ux*(-1 - uy) +
							mfaca*(-1 - ux)*(1 - uy) + mfacb*(-1 - ux)*(1 - uy) + mfacc*(-1 - ux)*(1 - uy) + mfcca*(1 - ux)*(1 - uy) + mfccb*(1 - ux)*(1 - uy) +
							mfccc*(1 - ux)*(1 - uy) - mfbca*ux*(1 - uy) - mfbcb*ux*(1 - uy) - mfbcc*ux*(1 - uy) - mfaba*(-1 - ux)*uy - mfabb*(-1 - ux)*uy -
							mfabc*(-1 - ux)*uy - mfcba*(1 - ux)*uy - mfcbb*(1 - ux)*uy - mfcbc*(1 - ux)*uy + mfbba*ux*uy + mfbbb*ux*uy + mfbbc*ux*uy;

						LBMReal Mom101 = mfaaa*(-1 - ux)*(-1 - uz) + mfaba*(-1 - ux)*(-1 - uz) + mfaca*(-1 - ux)*(-1 - uz) + mfcaa*(1 - ux)*(-1 - uz) +
							mfcba*(1 - ux)*(-1 - uz) + mfcca*(1 - ux)*(-1 - uz) - mfbaa*ux*(-1 - uz) - mfbba*ux*(-1 - uz) - mfbca*ux*(-1 - uz) +
							mfaac*(-1 - ux)*(1 - uz) + mfabc*(-1 - ux)*(1 - uz) + mfacc*(-1 - ux)*(1 - uz) + mfcac*(1 - ux)*(1 - uz) + mfcbc*(1 - ux)*(1 - uz) +
							mfccc*(1 - ux)*(1 - uz) - mfbac*ux*(1 - uz) - mfbbc*ux*(1 - uz) - mfbcc*ux*(1 - uz) - mfaab*(-1 - ux)*uz - mfabb*(-1 - ux)*uz -
							mfacb*(-1 - ux)*uz - mfcab*(1 - ux)*uz - mfcbb*(1 - ux)*uz - mfccb*(1 - ux)*uz + mfbab*ux*uz + mfbbb*ux*uz + mfbcb*ux*uz;

						LBMReal Mom011 = mfaaa*(-1 - uy)*(-1 - uz) + mfbaa*(-1 - uy)*(-1 - uz) + mfcaa*(-1 - uy)*(-1 - uz) + mfaca*(1 - uy)*(-1 - uz) +
							mfbca*(1 - uy)*(-1 - uz) + mfcca*(1 - uy)*(-1 - uz) - mfaba*uy*(-1 - uz) - mfbba*uy*(-1 - uz) - mfcba*uy*(-1 - uz) +
							mfaac*(-1 - uy)*(1 - uz) + mfbac*(-1 - uy)*(1 - uz) + mfcac*(-1 - uy)*(1 - uz) + mfacc*(1 - uy)*(1 - uz) + mfbcc*(1 - uy)*(1 - uz) +
							mfccc*(1 - uy)*(1 - uz) - mfabc*uy*(1 - uz) - mfbbc*uy*(1 - uz) - mfcbc*uy*(1 - uz) - mfaab*(-1 - uy)*uz - mfbab*(-1 - uy)*uz -
							mfcab*(-1 - uy)*uz - mfacb*(1 - uy)*uz - mfbcb*(1 - uy)*uz - mfccb*(1 - uy)*uz + mfabb*uy*uz + mfbbb*uy*uz + mfcbb*uy*uz;
						////

						LBMReal Mom111 = mfaaa*(-1 - ux)*(-1 - uy)*(-1 - uz) + mfcaa*(1 - ux)*(-1 - uy)*(-1 - uz) - mfbaa*ux*(-1 - uy)*(-1 - uz) +
							mfaca*(-1 - ux)*(1 - uy)*(-1 - uz) + mfcca*(1 - ux)*(1 - uy)*(-1 - uz) - mfbca*ux*(1 - uy)*(-1 - uz) -
							mfaba*(-1 - ux)*uy*(-1 - uz) - mfcba*(1 - ux)*uy*(-1 - uz) + mfbba*ux*uy*(-1 - uz) + mfaac*(-1 - ux)*(-1 - uy)*(1 - uz) +
							mfcac*(1 - ux)*(-1 - uy)*(1 - uz) - mfbac*ux*(-1 - uy)*(1 - uz) + mfacc*(-1 - ux)*(1 - uy)*(1 - uz) +
							mfccc*(1 - ux)*(1 - uy)*(1 - uz) - mfbcc*ux*(1 - uy)*(1 - uz) - mfabc*(-1 - ux)*uy*(1 - uz) - mfcbc*(1 - ux)*uy*(1 - uz) +
							mfbbc*ux*uy*(1 - uz) - mfaab*(-1 - ux)*(-1 - uy)*uz - mfcab*(1 - ux)*(-1 - uy)*uz + mfbab*ux*(-1 - uy)*uz -
							mfacb*(-1 - ux)*(1 - uy)*uz - mfccb*(1 - ux)*(1 - uy)*uz + mfbcb*ux*(1 - uy)*uz + mfabb*(-1 - ux)*uy*uz + mfcbb*(1 - ux)*uy*uz -
							mfbbb*ux*uy*uz;

						//(200)//
						LBMReal Mom200 = ((mfcaa + mfcab + mfcac + mfcba + mfcbb + mfcbc + mfcca + mfccb +
							mfccc)*pow(-1 + ux, 2) +
							(mfbaa + mfbab + mfbac + mfbba + mfbbb + mfbbc + mfbca +
								mfbcb + mfbcc)*pow(ux, 2) +
								(mfaaa + mfaab + mfaac + mfaba + mfabb + mfabc + mfaca +
									mfacb + mfacc)*pow(1 + ux, 2)) - Mom000 / 3;

						LBMReal Mom020 = ((mfaca + mfacb + mfacc)*pow(-1 + uy, 2) +
							(mfbca + mfbcb + mfbcc)*pow(-1 + uy, 2) +
							(mfcca + mfccb + mfccc)*pow(-1 + uy, 2) +
							(mfaba + mfabb + mfabc)*pow(uy, 2) +
							(mfbba + mfbbb + mfbbc)*pow(uy, 2) +
							(mfcba + mfcbb + mfcbc)*pow(uy, 2) +
							(mfaaa + mfaab + mfaac)*pow(1 + uy, 2) +
							(mfbaa + mfbab + mfbac)*pow(1 + uy, 2) +
							(mfcaa + mfcab + mfcac)*pow(1 + uy, 2)) - Mom000 / 3;

						LBMReal Mom002 = (mfaba + mfabc + mfaca + mfacc + mfbba + mfbbc + mfbca + mfbcc +
							mfcba + mfcbc + mfcca + mfccc + mfaac*pow(-1 + uz, 2) +
							mfbac*pow(-1 + uz, 2) + mfcac*pow(-1 + uz, 2) +
							2 * mfaba*uz - 2 * mfabc*uz + 2 * mfaca*uz - 2 * mfacc*uz +
							2 * mfbba*uz - 2 * mfbbc*uz + 2 * mfbca*uz - 2 * mfbcc*uz +
							2 * mfcba*uz - 2 * mfcbc*uz + 2 * mfcca*uz - 2 * mfccc*uz +
							mfaab*pow(uz, 2) + mfaba*pow(uz, 2) + mfabb*pow(uz, 2) +
							mfabc*pow(uz, 2) + mfaca*pow(uz, 2) + mfacb*pow(uz, 2) +
							mfacc*pow(uz, 2) + mfbab*pow(uz, 2) + mfbba*pow(uz, 2) +
							mfbbb*pow(uz, 2) + mfbbc*pow(uz, 2) + mfbca*pow(uz, 2) +
							mfbcb*pow(uz, 2) + mfbcc*pow(uz, 2) + mfcab*pow(uz, 2) +
							mfcba*pow(uz, 2) + mfcbb*pow(uz, 2) + mfcbc*pow(uz, 2) +
							mfcca*pow(uz, 2) + mfccb*pow(uz, 2) + mfccc*pow(uz, 2) +
							mfaaa*pow(1 + uz, 2) + mfbaa*pow(1 + uz, 2) +
							mfcaa*pow(1 + uz, 2)) - Mom000 / 3;
						////

						//(210)//
						LBMReal Mom210 = (pow(1 + ux, 2)*(-((mfaca + mfacb + mfacc)*(-1 + uy)) -
							(mfaba + mfabb + mfabc)*uy -
							(mfaaa + mfaab + mfaac)*(1 + uy)) +
							pow(ux, 2)*(-((mfbca + mfbcb + mfbcc)*(-1 + uy)) -
							(mfbba + mfbbb + mfbbc)*uy -
								(mfbaa + mfbab + mfbac)*(1 + uy)) +
							pow(-1 + ux, 2)*(-((mfcca + mfccb + mfccc)*(-1 + uy)) -
							(mfcba + mfcbb + mfcbc)*uy -
								(mfcaa + mfcab + mfcac)*(1 + uy))) - Mom010 / 3;

						LBMReal Mom201 = (-(pow(1 + ux, 2)*(mfaba - mfabc + mfaca - mfacc +
							mfaac*(-1 + uz) + mfaab*uz + mfaba*uz + mfabb*uz +
							mfabc*uz + mfaca*uz + mfacb*uz + mfacc*uz +
							mfaaa*(1 + uz))) -
							pow(ux, 2)*(mfbba - mfbbc + mfbca - mfbcc +
								mfbac*(-1 + uz) + mfbab*uz + mfbba*uz + mfbbb*uz +
								mfbbc*uz + mfbca*uz + mfbcb*uz + mfbcc*uz + mfbaa*(1 + uz))
							- pow(-1 + ux, 2)*(mfcba - mfcbc + mfcca - mfccc +
								mfcac*(-1 + uz) + mfcab*uz + mfcba*uz + mfcbb*uz +
								mfcbc*uz + mfcca*uz + mfccb*uz + mfccc*uz + mfcaa*(1 + uz))) - Mom001 / 3;

						LBMReal Mom120 = ((-1 - ux)*((mfaca + mfacb + mfacc)*pow(-1 + uy, 2) +
							(mfaba + mfabb + mfabc)*pow(uy, 2) +
							(mfaaa + mfaab + mfaac)*pow(1 + uy, 2)) -
							ux*((mfbca + mfbcb + mfbcc)*pow(-1 + uy, 2) +
							(mfbba + mfbbb + mfbbc)*pow(uy, 2) +
								(mfbaa + mfbab + mfbac)*pow(1 + uy, 2)) +
								(1 - ux)*((mfcca + mfccb + mfccc)*pow(-1 + uy, 2) +
							(mfcba + mfcbb + mfcbc)*pow(uy, 2) +
									(mfcaa + mfcab + mfcac)*pow(1 + uy, 2))) - Mom100 / 3;


						LBMReal Mom102 = (-((1 + ux)*(mfaba + mfabc + mfaca + mfacc +
							mfaac*pow(-1 + uz, 2) + 2 * mfaba*uz - 2 * mfabc*uz +
							2 * mfaca*uz - 2 * mfacc*uz + mfaab*pow(uz, 2) +
							mfaba*pow(uz, 2) + mfabb*pow(uz, 2) +
							mfabc*pow(uz, 2) + mfaca*pow(uz, 2) +
							mfacb*pow(uz, 2) + mfacc*pow(uz, 2) +
							mfaaa*pow(1 + uz, 2))) -
							ux*(mfbba + mfbbc + mfbca + mfbcc + mfbac*pow(-1 + uz, 2) +
								2 * mfbba*uz - 2 * mfbbc*uz + 2 * mfbca*uz - 2 * mfbcc*uz +
								mfbab*pow(uz, 2) + mfbba*pow(uz, 2) +
								mfbbb*pow(uz, 2) + mfbbc*pow(uz, 2) +
								mfbca*pow(uz, 2) + mfbcb*pow(uz, 2) +
								mfbcc*pow(uz, 2) + mfbaa*pow(1 + uz, 2)) -
								(-1 + ux)*(mfcba + mfcbc + mfcca + mfccc +
									mfcac*pow(-1 + uz, 2) + 2 * mfcba*uz - 2 * mfcbc*uz +
									2 * mfcca*uz - 2 * mfccc*uz + mfcab*pow(uz, 2) +
									mfcba*pow(uz, 2) + mfcbb*pow(uz, 2) +
									mfcbc*pow(uz, 2) + mfcca*pow(uz, 2) +
									mfccb*pow(uz, 2) + mfccc*pow(uz, 2) +
									mfcaa*pow(1 + uz, 2))) - Mom100 / 3;

						LBMReal Mom021 = (-(pow(1 + uy, 2)*(mfaac*(-1 + uz) + mfaab*uz +
							mfaaa*(1 + uz))) -
							pow(uy, 2)*(mfabc*(-1 + uz) + mfabb*uz + mfaba*(1 + uz)) -
							pow(-1 + uy, 2)*(mfacc*(-1 + uz) + mfacb*uz +
								mfaca*(1 + uz)) - pow(1 + uy, 2)*
								(mfbac*(-1 + uz) + mfbab*uz + mfbaa*(1 + uz)) -
							pow(uy, 2)*(mfbbc*(-1 + uz) + mfbbb*uz + mfbba*(1 + uz)) -
							pow(-1 + uy, 2)*(mfbcc*(-1 + uz) + mfbcb*uz +
								mfbca*(1 + uz)) - pow(1 + uy, 2)*
								(mfcac*(-1 + uz) + mfcab*uz + mfcaa*(1 + uz)) -
							pow(uy, 2)*(mfcbc*(-1 + uz) + mfcbb*uz + mfcba*(1 + uz)) -
							pow(-1 + uy, 2)*(mfccc*(-1 + uz) + mfccb*uz + mfcca*(1 + uz))) - Mom001 / 3;

						LBMReal Mom012 = (-((1 + uy)*(mfaac*pow(-1 + uz, 2) + mfaab*pow(uz, 2) +
							mfaaa*pow(1 + uz, 2))) -
							uy*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) -
								(-1 + uy)*(mfacc*pow(-1 + uz, 2) + mfacb*pow(uz, 2) +
									mfaca*pow(1 + uz, 2)) -
									(1 + uy)*(mfbac*pow(-1 + uz, 2) + mfbab*pow(uz, 2) +
										mfbaa*pow(1 + uz, 2)) -
							uy*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
								mfbba*pow(1 + uz, 2)) -
								(-1 + uy)*(mfbcc*pow(-1 + uz, 2) + mfbcb*pow(uz, 2) +
									mfbca*pow(1 + uz, 2)) -
									(1 + uy)*(mfcac*pow(-1 + uz, 2) + mfcab*pow(uz, 2) +
										mfcaa*pow(1 + uz, 2)) -
							uy*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
								mfcba*pow(1 + uz, 2)) -
								(-1 + uy)*(mfccc*pow(-1 + uz, 2) + mfccb*pow(uz, 2) +
									mfcca*pow(1 + uz, 2))) - Mom010 / 3;
						////


						//(220)//
						LBMReal Mom220 = (pow(1 + ux, 2)*((mfaca + mfacb + mfacc)*pow(-1 + uy, 2) +
							(mfaba + mfabb + mfabc)*pow(uy, 2) +
							(mfaaa + mfaab + mfaac)*pow(1 + uy, 2)) +
							pow(ux, 2)*((mfbca + mfbcb + mfbcc)*pow(-1 + uy, 2) +
							(mfbba + mfbbb + mfbbc)*pow(uy, 2) +
								(mfbaa + mfbab + mfbac)*pow(1 + uy, 2)) +
							pow(-1 + ux, 2)*((mfcca + mfccb + mfccc)*pow(-1 + uy, 2) +
							(mfcba + mfcbb + mfcbc)*pow(uy, 2) +
								(mfcaa + mfcab + mfcac)*pow(1 + uy, 2))) - Mom000 / 9;

						LBMReal Mom202 = (pow(1 + ux, 2)*(mfaba + mfabc + mfaca + mfacc +
							mfaac*pow(-1 + uz, 2) + 2 * mfaba*uz - 2 * mfabc*uz +
							2 * mfaca*uz - 2 * mfacc*uz + mfaab*pow(uz, 2) +
							mfaba*pow(uz, 2) + mfabb*pow(uz, 2) +
							mfabc*pow(uz, 2) + mfaca*pow(uz, 2) +
							mfacb*pow(uz, 2) + mfacc*pow(uz, 2) +
							mfaaa*pow(1 + uz, 2)) +
							pow(ux, 2)*(mfbba + mfbbc + mfbca + mfbcc +
								mfbac*pow(-1 + uz, 2) + 2 * mfbba*uz - 2 * mfbbc*uz +
								2 * mfbca*uz - 2 * mfbcc*uz + mfbab*pow(uz, 2) +
								mfbba*pow(uz, 2) + mfbbb*pow(uz, 2) +
								mfbbc*pow(uz, 2) + mfbca*pow(uz, 2) +
								mfbcb*pow(uz, 2) + mfbcc*pow(uz, 2) +
								mfbaa*pow(1 + uz, 2)) +
							pow(-1 + ux, 2)*(mfcba + mfcbc + mfcca + mfccc +
								mfcac*pow(-1 + uz, 2) + 2 * mfcba*uz - 2 * mfcbc*uz +
								2 * mfcca*uz - 2 * mfccc*uz + mfcab*pow(uz, 2) +
								mfcba*pow(uz, 2) + mfcbb*pow(uz, 2) +
								mfcbc*pow(uz, 2) + mfcca*pow(uz, 2) +
								mfccb*pow(uz, 2) + mfccc*pow(uz, 2) +
								mfcaa*pow(1 + uz, 2))) - Mom000 / 9;

						LBMReal Mom022 = (pow(1 + uy, 2)*(mfaac*pow(-1 + uz, 2) + mfaab*pow(uz, 2) +
							mfaaa*pow(1 + uz, 2)) +
							pow(uy, 2)*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) +
							pow(-1 + uy, 2)*(mfacc*pow(-1 + uz, 2) +
								mfacb*pow(uz, 2) + mfaca*pow(1 + uz, 2)) +
							pow(1 + uy, 2)*(mfbac*pow(-1 + uz, 2) + mfbab*pow(uz, 2) +
								mfbaa*pow(1 + uz, 2)) +
							pow(uy, 2)*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
								mfbba*pow(1 + uz, 2)) +
							pow(-1 + uy, 2)*(mfbcc*pow(-1 + uz, 2) +
								mfbcb*pow(uz, 2) + mfbca*pow(1 + uz, 2)) +
							pow(1 + uy, 2)*(mfcac*pow(-1 + uz, 2) + mfcab*pow(uz, 2) +
								mfcaa*pow(1 + uz, 2)) +
							pow(uy, 2)*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
								mfcba*pow(1 + uz, 2)) +
							pow(-1 + uy, 2)*(mfccc*pow(-1 + uz, 2) +
								mfccb*pow(uz, 2) + mfcca*pow(1 + uz, 2))) - Mom000 / 9;
						////

						//(221)//
						LBMReal Mom221 = (pow(1 + ux, 2)*(-(pow(1 + uy, 2)*
							(mfaac*(-1 + uz) + mfaab*uz + mfaaa*(1 + uz))) -
							pow(uy, 2)*(mfabc*(-1 + uz) + mfabb*uz +
								mfaba*(1 + uz)) -
							pow(-1 + uy, 2)*(mfacc*(-1 + uz) + mfacb*uz +
								mfaca*(1 + uz))) +
							pow(ux, 2)*(-(pow(1 + uy, 2)*
							(mfbac*(-1 + uz) + mfbab*uz + mfbaa*(1 + uz))) -
								pow(uy, 2)*(mfbbc*(-1 + uz) + mfbbb*uz +
									mfbba*(1 + uz)) -
								pow(-1 + uy, 2)*(mfbcc*(-1 + uz) + mfbcb*uz +
									mfbca*(1 + uz))) +
							pow(-1 + ux, 2)*(-(pow(1 + uy, 2)*
							(mfcac*(-1 + uz) + mfcab*uz + mfcaa*(1 + uz))) -
								pow(uy, 2)*(mfcbc*(-1 + uz) + mfcbb*uz +
									mfcba*(1 + uz)) -
								pow(-1 + uy, 2)*(mfccc*(-1 + uz) + mfccb*uz +
									mfcca*(1 + uz)))) - Mom001 / 9;

						LBMReal Mom212 = (pow(1 + ux, 2)*(-((1 + uy)*
							(mfaac*pow(-1 + uz, 2) + mfaab*pow(uz, 2) +
								mfaaa*pow(1 + uz, 2))) -
							uy*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) -
								(-1 + uy)*(mfacc*pow(-1 + uz, 2) + mfacb*pow(uz, 2) +
									mfaca*pow(1 + uz, 2))) +
							pow(ux, 2)*(-((1 + uy)*
							(mfbac*pow(-1 + uz, 2) + mfbab*pow(uz, 2) +
								mfbaa*pow(1 + uz, 2))) -
								uy*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
									mfbba*pow(1 + uz, 2)) -
									(-1 + uy)*(mfbcc*pow(-1 + uz, 2) + mfbcb*pow(uz, 2) +
										mfbca*pow(1 + uz, 2))) +
							pow(-1 + ux, 2)*(-((1 + uy)*
							(mfcac*pow(-1 + uz, 2) + mfcab*pow(uz, 2) +
								mfcaa*pow(1 + uz, 2))) -
								uy*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
									mfcba*pow(1 + uz, 2)) -
									(-1 + uy)*(mfccc*pow(-1 + uz, 2) + mfccb*pow(uz, 2) +
										mfcca*pow(1 + uz, 2)))) - Mom010 / 9;

						LBMReal Mom122 = ((-1 - ux)*(pow(1 + uy, 2)*
							(mfaac*pow(-1 + uz, 2) + mfaab*pow(uz, 2) +
								mfaaa*pow(1 + uz, 2)) +
							pow(uy, 2)*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) +
							pow(-1 + uy, 2)*(mfacc*pow(-1 + uz, 2) +
								mfacb*pow(uz, 2) + mfaca*pow(1 + uz, 2))) -
							ux*(pow(1 + uy, 2)*(mfbac*pow(-1 + uz, 2) +
								mfbab*pow(uz, 2) + mfbaa*pow(1 + uz, 2)) +
								pow(uy, 2)*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
									mfbba*pow(1 + uz, 2)) +
								pow(-1 + uy, 2)*(mfbcc*pow(-1 + uz, 2) +
									mfbcb*pow(uz, 2) + mfbca*pow(1 + uz, 2))) +
									(1 - ux)*(pow(1 + uy, 2)*
							(mfcac*pow(-1 + uz, 2) + mfcab*pow(uz, 2) +
								mfcaa*pow(1 + uz, 2)) +
										pow(uy, 2)*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
											mfcba*pow(1 + uz, 2)) +
										pow(-1 + uy, 2)*(mfccc*pow(-1 + uz, 2) +
											mfccb*pow(uz, 2) + mfcca*pow(1 + uz, 2)))) - Mom100 / 9;
						////

						//(211)//
						LBMReal Mom211 = (pow(1 + ux, 2)*((1 + uy)*(mfaac*(-1 + uz) + mfaab*uz +
							mfaaa*(1 + uz)) +
							uy*(mfabc*(-1 + uz) + mfabb*uz + mfaba*(1 + uz)) +
							(-1 + uy)*(mfacc*(-1 + uz) + mfacb*uz + mfaca*(1 + uz))) +
							pow(ux, 2)*((1 + uy)*(mfbac*(-1 + uz) + mfbab*uz +
								mfbaa*(1 + uz)) +
								uy*(mfbbc*(-1 + uz) + mfbbb*uz + mfbba*(1 + uz)) +
								(-1 + uy)*(mfbcc*(-1 + uz) + mfbcb*uz + mfbca*(1 + uz))) +
							pow(-1 + ux, 2)*((1 + uy)*
							(mfcac*(-1 + uz) + mfcab*uz + mfcaa*(1 + uz)) +
								uy*(mfcbc*(-1 + uz) + mfcbb*uz + mfcba*(1 + uz)) +
								(-1 + uy)*(mfccc*(-1 + uz) + mfccb*uz + mfcca*(1 + uz)))) - Mom011 / 3;

						LBMReal Mom121 = ((-1 - ux)*(-(pow(1 + uy, 2)*
							(mfaac*(-1 + uz) + mfaab*uz + mfaaa*(1 + uz))) -
							pow(uy, 2)*(mfabc*(-1 + uz) + mfabb*uz +
								mfaba*(1 + uz)) -
							pow(-1 + uy, 2)*(mfacc*(-1 + uz) + mfacb*uz +
								mfaca*(1 + uz))) -
							ux*(-(pow(1 + uy, 2)*(mfbac*(-1 + uz) + mfbab*uz +
								mfbaa*(1 + uz))) -
								pow(uy, 2)*(mfbbc*(-1 + uz) + mfbbb*uz +
									mfbba*(1 + uz)) -
								pow(-1 + uy, 2)*(mfbcc*(-1 + uz) + mfbcb*uz +
									mfbca*(1 + uz))) +
									(1 - ux)*(-(pow(1 + uy, 2)*
							(mfcac*(-1 + uz) + mfcab*uz + mfcaa*(1 + uz))) -
										pow(uy, 2)*(mfcbc*(-1 + uz) + mfcbb*uz +
											mfcba*(1 + uz)) -
										pow(-1 + uy, 2)*(mfccc*(-1 + uz) + mfccb*uz +
											mfcca*(1 + uz)))) - Mom101 / 3;

						LBMReal Mom112 = ((-1 - ux)*(-((1 + uy)*(mfaac*pow(-1 + uz, 2) +
							mfaab*pow(uz, 2) + mfaaa*pow(1 + uz, 2))) -
							uy*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) -
								(-1 + uy)*(mfacc*pow(-1 + uz, 2) + mfacb*pow(uz, 2) +
									mfaca*pow(1 + uz, 2))) -
							ux*(-((1 + uy)*(mfbac*pow(-1 + uz, 2) + mfbab*pow(uz, 2) +
								mfbaa*pow(1 + uz, 2))) -
								uy*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
									mfbba*pow(1 + uz, 2)) -
									(-1 + uy)*(mfbcc*pow(-1 + uz, 2) + mfbcb*pow(uz, 2) +
										mfbca*pow(1 + uz, 2))) +
										(1 - ux)*(-((1 + uy)*(mfcac*pow(-1 + uz, 2) +
											mfcab*pow(uz, 2) + mfcaa*pow(1 + uz, 2))) -
											uy*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
												mfcba*pow(1 + uz, 2)) -
												(-1 + uy)*(mfccc*pow(-1 + uz, 2) + mfccb*pow(uz, 2) +
													mfcca*pow(1 + uz, 2)))) - Mom110 / 3;
						////

						//(222)//
						LBMReal Mom222 = (pow(1 + ux, 2)*(pow(1 + uy, 2)*
							(mfaac*pow(-1 + uz, 2) + mfaab*pow(uz, 2) +
								mfaaa*pow(1 + uz, 2)) +
							pow(uy, 2)*(mfabc*pow(-1 + uz, 2) + mfabb*pow(uz, 2) +
								mfaba*pow(1 + uz, 2)) +
							pow(-1 + uy, 2)*(mfacc*pow(-1 + uz, 2) +
								mfacb*pow(uz, 2) + mfaca*pow(1 + uz, 2))) +
							pow(ux, 2)*(pow(1 + uy, 2)*
							(mfbac*pow(-1 + uz, 2) + mfbab*pow(uz, 2) +
								mfbaa*pow(1 + uz, 2)) +
								pow(uy, 2)*(mfbbc*pow(-1 + uz, 2) + mfbbb*pow(uz, 2) +
									mfbba*pow(1 + uz, 2)) +
								pow(-1 + uy, 2)*(mfbcc*pow(-1 + uz, 2) +
									mfbcb*pow(uz, 2) + mfbca*pow(1 + uz, 2))) +
							pow(-1 + ux, 2)*(pow(1 + uy, 2)*
							(mfcac*pow(-1 + uz, 2) + mfcab*pow(uz, 2) +
								mfcaa*pow(1 + uz, 2)) +
								pow(uy, 2)*(mfcbc*pow(-1 + uz, 2) + mfcbb*pow(uz, 2) +
									mfcba*pow(1 + uz, 2)) +
								pow(-1 + uy, 2)*(mfccc*pow(-1 + uz, 2) +
									mfccb*pow(uz, 2) + mfcca*pow(1 + uz, 2)))) - Mom000 / 27;
						////





						LBMReal Meq000 = drho+dlambda*c1o2;


						// relaxation Central Moment MRT

						Mom000 = Meq000;

						Mom000 += dlambda*c1o2;

						Mom100 = (1 - collFactorH) * Mom100;
						Mom010 = (1 - collFactorH) * Mom010;
						Mom001 = (1 - collFactorH) * Mom001;

						Mom110 = 0;
						Mom101 = 0;
						Mom011 = 0;

						Mom111 = 0;

						//(200)//
						Mom200 = Mom000 / 3;
						Mom020 = Mom000 / 3;
						Mom002 = Mom000 / 3;
						////

						//(210)//
						Mom210 = Mom010 / 3;
						Mom201 = Mom001 / 3;
						Mom120 = Mom100 / 3;


						Mom102 = Mom100 / 3;
						Mom021 = Mom001 / 3;
						Mom012 = Mom010 / 3;
						////


						//(220)//
						Mom220 = Mom000 / 9;
						Mom202 = Mom000 / 9;
						Mom022 = Mom000 / 9;
						////

						//(221)//
						Mom221 = Mom001 / 9;
						Mom212 = Mom010 / 9;
						Mom122 = Mom100 / 9;
						////

						//(211)//
						Mom211 = Mom011 / 3;
						Mom121 = Mom101 / 3;
						Mom112 = Mom110 / 3;
						////

						//(222)//
						Mom222 = Mom000 / 27;
						////
						


						//Back transformation to distributions

						mfcbb = (Mom122 + Mom222 + 2 * Mom122*ux + Mom022*ux*(1 + ux) +
							2 * (Mom112 + Mom212 + 2 * Mom112*ux + Mom012*ux*(1 + ux))*uy +
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*
							(-1 + pow(uy, 2)) -
							2 * (-Mom221 - Mom021*ux*(1 + ux) - Mom121*(1 + 2 * ux) -
								2 * (Mom111 + Mom211 + 2 * Mom111*ux + Mom011*ux*(1 + ux))*
								uy - (Mom101 + Mom201 + 2 * Mom101*ux +
									Mom001*ux*(1 + ux))*(-1 + pow(uy, 2)))*uz +
									(-Mom220 - Mom020*ux*(1 + ux) - Mom120*(1 + 2 * ux) -
										2 * (Mom110 + Mom210 + 2 * Mom110*ux + Mom010*ux*(1 + ux))*
										uy - (Mom100 + Mom200 + 2 * Mom100*ux +
											Mom000*ux*(1 + ux))*(-1 + pow(uy, 2)))*
											(1 - pow(uz, 2))) / 2.;
						mfbcb = (Mom222 + 2 * Mom122*ux + Mom022*(-1 + pow(ux, 2)) +
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*uy*
							(1 + uy) + (Mom212 + 2 * Mom112*ux +
								Mom012*(-1 + pow(ux, 2)))*(1 + 2 * uy) -
							2 * (Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*uy*
								(1 + uy) - (Mom211 + 2 * Mom111*ux +
									Mom011*(-1 + pow(ux, 2)))*(1 + 2 * uy))*uz +
									(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*uy*
										(1 + uy) - (Mom210 + 2 * Mom110*ux +
											Mom010*(-1 + pow(ux, 2)))*(1 + 2 * uy))*
											(1 - pow(uz, 2))) / 2.;
						mfbbc = (Mom222 + 2 * Mom122*ux + Mom022*(-1 + pow(ux, 2)) +
							2 * (Mom212 + 2 * Mom112*ux + Mom012*(-1 + pow(ux, 2)))*uy +
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*
							(-1 + pow(uy, 2)) +
							(Mom220 + 2 * Mom120*ux + Mom020*(-1 + pow(ux, 2)) +
								2 * (Mom210 + 2 * Mom110*ux + Mom010*(-1 + pow(ux, 2)))*
								uy + (Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*
								(-1 + pow(uy, 2)))*uz*(1 + uz) +
								(Mom221 + 2 * Mom121*ux + Mom021*(-1 + pow(ux, 2)) +
									2 * (Mom211 + 2 * Mom111*ux + Mom011*(-1 + pow(ux, 2)))*
									uy + (Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*
									(-1 + pow(uy, 2)))*(1 + 2 * uz)) / 2.;
						mfccb = (-Mom222 - Mom022*ux*(1 + ux) - Mom122*(1 + 2 * ux) -
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*uy*
							(1 + uy) - (Mom112 + Mom212 + 2 * Mom112*ux +
								Mom012*ux*(1 + ux))*(1 + 2 * uy) -
							2 * (Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*uy*
								(1 + uy) + (Mom111 + Mom211 + 2 * Mom111*ux +
									Mom011*ux*(1 + ux))*(1 + 2 * uy))*uz +
									(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*uy*
										(1 + uy) + (Mom110 + Mom210 + 2 * Mom110*ux +
											Mom010*ux*(1 + ux))*(1 + 2 * uy))*(1 - pow(uz, 2))) / 4.;
						mfacb = (Mom122 - Mom222 - 2 * Mom122*ux - Mom022*(-1 + ux)*ux -
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*uy*
							(1 + uy) - (Mom212 + Mom012*(-1 + ux)*ux +
								Mom112*(-1 + 2 * ux))*(1 + 2 * uy) -
							2 * (Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*uy*
								(1 + uy) + (Mom211 + Mom011*(-1 + ux)*ux +
									Mom111*(-1 + 2 * ux))*(1 + 2 * uy))*uz +
									(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*uy*
										(1 + uy) + (Mom210 + Mom010*(-1 + ux)*ux +
											Mom110*(-1 + 2 * ux))*(1 + 2 * uy))*(1 - pow(uz, 2))) / 4.;
						mfcbc = (-Mom222 - Mom022*ux*(1 + ux) - Mom122*(1 + 2 * ux) -
							2 * (Mom112 + Mom212 + 2 * Mom112*ux + Mom012*ux*(1 + ux))*uy -
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*
							(-1 + pow(uy, 2)) +
							(-Mom220 - Mom020*ux*(1 + ux) - Mom120*(1 + 2 * ux) -
								2 * (Mom110 + Mom210 + 2 * Mom110*ux + Mom010*ux*(1 + ux))*
								uy - (Mom100 + Mom200 + 2 * Mom100*ux +
									Mom000*ux*(1 + ux))*(-1 + pow(uy, 2)))*uz*(1 + uz) +
									(-Mom221 - Mom021*ux*(1 + ux) - Mom121*(1 + 2 * ux) -
										2 * (Mom111 + Mom211 + 2 * Mom111*ux + Mom011*ux*(1 + ux))*
										uy - (Mom101 + Mom201 + 2 * Mom101*ux +
											Mom001*ux*(1 + ux))*(-1 + pow(uy, 2)))*(1 + 2 * uz)) / 4.;
						mfabc = (Mom122 - Mom222 - 2 * Mom122*ux - Mom022*(-1 + ux)*ux -
							2 * (Mom212 + Mom012*(-1 + ux)*ux + Mom112*(-1 + 2 * ux))*uy -
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*
							(-1 + pow(uy, 2)) +
							(Mom120 - Mom220 - 2 * Mom120*ux - Mom020*(-1 + ux)*ux -
								2 * (Mom210 + Mom010*(-1 + ux)*ux + Mom110*(-1 + 2 * ux))*
								uy - (Mom200 + Mom000*(-1 + ux)*ux +
									Mom100*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*uz*(1 + uz) +
									(Mom121 - Mom221 - 2 * Mom121*ux - Mom021*(-1 + ux)*ux -
										2 * (Mom211 + Mom011*(-1 + ux)*ux + Mom111*(-1 + 2 * ux))*
										uy - (Mom201 + Mom001*(-1 + ux)*ux +
											Mom101*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*(1 + 2 * uz)) / 4.;
						mfbcc = (Mom022 - Mom222 - 2 * Mom122*ux - Mom022*pow(ux, 2) -
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*uy*
							(1 + uy) - (Mom212 + 2 * Mom112*ux +
								Mom012*(-1 + pow(ux, 2)))*(1 + 2 * uy) +
								(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*uy*
									(1 + uy) - (Mom210 + 2 * Mom110*ux +
										Mom010*(-1 + pow(ux, 2)))*(1 + 2 * uy))*uz*(1 + uz) +
										(Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*uy*
											(1 + uy) - (Mom211 + 2 * Mom111*ux +
												Mom011*(-1 + pow(ux, 2)))*(1 + 2 * uy))*(1 + 2 * uz)) / 4.;
						mfbac = (Mom022 - Mom222 - 2 * Mom122*ux - Mom022*pow(ux, 2) -
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*
							(-1 + uy)*uy - (Mom212 + 2 * Mom112*ux +
								Mom012*(-1 + pow(ux, 2)))*(-1 + 2 * uy) +
								(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*
									(-1 + uy)*uy - (Mom210 + 2 * Mom110*ux +
										Mom010*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*uz*(1 + uz) +
										(Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*
											(-1 + uy)*uy - (Mom211 + 2 * Mom111*ux +
												Mom011*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*(1 + 2 * uz)) / 4.;
						mfccc = (Mom122 + Mom222 + 2 * Mom122*ux + Mom022*ux*(1 + ux) +
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*uy*
							(1 + uy) + (Mom112 + Mom212 + 2 * Mom112*ux +
								Mom012*ux*(1 + ux))*(1 + 2 * uy) +
								(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*uy*
									(1 + uy) + (Mom110 + Mom210 + 2 * Mom110*ux +
										Mom010*ux*(1 + ux))*(1 + 2 * uy))*uz*(1 + uz) +
										(Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*uy*
											(1 + uy) + (Mom111 + Mom211 + 2 * Mom111*ux +
												Mom011*ux*(1 + ux))*(1 + 2 * uy))*(1 + 2 * uz)) / 8.;
						mfacc = (Mom222 + Mom022*(-1 + ux)*ux + Mom122*(-1 + 2 * ux) +
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*uy*
							(1 + uy) + (Mom212 + Mom012*(-1 + ux)*ux +
								Mom112*(-1 + 2 * ux))*(1 + 2 * uy) +
								(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*uy*
									(1 + uy) + (Mom210 + Mom010*(-1 + ux)*ux +
										Mom110*(-1 + 2 * ux))*(1 + 2 * uy))*uz*(1 + uz) +
										(Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*uy*
											(1 + uy) + (Mom211 + Mom011*(-1 + ux)*ux +
												Mom111*(-1 + 2 * ux))*(1 + 2 * uy))*(1 + 2 * uz)) / 8.;
						mfcac = (Mom122 + Mom222 + 2 * Mom122*ux + Mom022*ux*(1 + ux) +
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*
							(-1 + uy)*uy + (Mom112 + Mom212 + 2 * Mom112*ux +
								Mom012*ux*(1 + ux))*(-1 + 2 * uy) +
								(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*
									(-1 + uy)*uy + (Mom110 + Mom210 + 2 * Mom110*ux +
										Mom010*ux*(1 + ux))*(-1 + 2 * uy))*uz*(1 + uz) +
										(Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*
											(-1 + uy)*uy + (Mom111 + Mom211 + 2 * Mom111*ux +
												Mom011*ux*(1 + ux))*(-1 + 2 * uy))*(1 + 2 * uz)) / 8.;
						mfaac = (Mom222 + Mom022*(-1 + ux)*ux + Mom122*(-1 + 2 * ux) +
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*
							(-1 + uy)*uy + (Mom212 + Mom012*(-1 + ux)*ux +
								Mom112*(-1 + 2 * ux))*(-1 + 2 * uy) +
								(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*
									(-1 + uy)*uy + (Mom210 + Mom010*(-1 + ux)*ux +
										Mom110*(-1 + 2 * ux))*(-1 + 2 * uy))*uz*(1 + uz) +
										(Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*
											(-1 + uy)*uy + (Mom211 + Mom011*(-1 + ux)*ux +
												Mom111*(-1 + 2 * ux))*(-1 + 2 * uy))*(1 + 2 * uz)) / 8.;

						mfabb = (Mom222 + Mom022*(-1 + ux)*ux + Mom122*(-1 + 2 * ux) +
							2 * (Mom212 + Mom012*(-1 + ux)*ux + Mom112*(-1 + 2 * ux))*uy +
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*
							(-1 + pow(uy, 2)) -
							2 * (Mom121 - Mom221 - 2 * Mom121*ux - Mom021*(-1 + ux)*ux -
								2 * (Mom211 + Mom011*(-1 + ux)*ux + Mom111*(-1 + 2 * ux))*
								uy - (Mom201 + Mom001*(-1 + ux)*ux +
									Mom101*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*uz +
									(Mom120 - Mom220 - 2 * Mom120*ux - Mom020*(-1 + ux)*ux -
										2 * (Mom210 + Mom010*(-1 + ux)*ux + Mom110*(-1 + 2 * ux))*
										uy - (Mom200 + Mom000*(-1 + ux)*ux +
											Mom100*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*
											(1 - pow(uz, 2))) / 2.;
						mfbab = (Mom222 + 2 * Mom122*ux + Mom022*(-1 + pow(ux, 2)) +
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*
							(-1 + uy)*uy + (Mom212 + 2 * Mom112*ux +
								Mom012*(-1 + pow(ux, 2)))*(-1 + 2 * uy) -
							2 * (Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*
								(-1 + uy)*uy - (Mom211 + 2 * Mom111*ux +
									Mom011*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*uz +
									(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*
										(-1 + uy)*uy - (Mom210 + 2 * Mom110*ux +
											Mom010*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*
											(1 - pow(uz, 2))) / 2.;
						mfbba = (Mom222 + 2 * Mom122*ux + Mom022*(-1 + pow(ux, 2)) +
							2 * (Mom212 + 2 * Mom112*ux + Mom012*(-1 + pow(ux, 2)))*uy +
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*
							(-1 + pow(uy, 2)) +
							(Mom220 + 2 * Mom120*ux + Mom020*(-1 + pow(ux, 2)) +
								2 * (Mom210 + 2 * Mom110*ux + Mom010*(-1 + pow(ux, 2)))*
								uy + (Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*
								(-1 + pow(uy, 2)))*(-1 + uz)*uz +
								(Mom221 + 2 * Mom121*ux + Mom021*(-1 + pow(ux, 2)) +
									2 * (Mom211 + 2 * Mom111*ux + Mom011*(-1 + pow(ux, 2)))*
									uy + (Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*
									(-1 + pow(uy, 2)))*(-1 + 2 * uz)) / 2.;
						mfaab = (Mom122 - Mom222 - 2 * Mom122*ux - Mom022*(-1 + ux)*ux -
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*
							(-1 + uy)*uy - (Mom212 + Mom012*(-1 + ux)*ux +
								Mom112*(-1 + 2 * ux))*(-1 + 2 * uy) -
							2 * (Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*
								(-1 + uy)*uy + (Mom211 + Mom011*(-1 + ux)*ux +
									Mom111*(-1 + 2 * ux))*(-1 + 2 * uy))*uz +
									(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*
										(-1 + uy)*uy + (Mom210 + Mom010*(-1 + ux)*ux +
											Mom110*(-1 + 2 * ux))*(-1 + 2 * uy))*(1 - pow(uz, 2))) / 4.;
						mfcab = (-Mom222 - Mom022*ux*(1 + ux) - Mom122*(1 + 2 * ux) -
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*
							(-1 + uy)*uy - (Mom112 + Mom212 + 2 * Mom112*ux +
								Mom012*ux*(1 + ux))*(-1 + 2 * uy) -
							2 * (Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*
								(-1 + uy)*uy + (Mom111 + Mom211 + 2 * Mom111*ux +
									Mom011*ux*(1 + ux))*(-1 + 2 * uy))*uz +
									(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*
										(-1 + uy)*uy + (Mom110 + Mom210 + 2 * Mom110*ux +
											Mom010*ux*(1 + ux))*(-1 + 2 * uy))*(1 - pow(uz, 2))) / 4.;
						mfaba = (Mom122 - Mom222 - 2 * Mom122*ux - Mom022*(-1 + ux)*ux - 2 * (Mom212 + Mom012*(-1 + ux)*ux + Mom112*(-1 + 2 * ux))*uy -
							(Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*(-1 + pow(uy, 2)) +
							(Mom120 - Mom220 - 2 * Mom120*ux - Mom020*(-1 + ux)*ux - 2 * (Mom210 + Mom010*(-1 + ux)*ux + Mom110*(-1 + 2 * ux))*uy -
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*(-1 + uz)*uz +
								(Mom121 - Mom221 - 2 * Mom121*ux - Mom021*(-1 + ux)*ux - 2 * (Mom211 + Mom011*(-1 + ux)*ux + Mom111*(-1 + 2 * ux))*uy -
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*(-1 + pow(uy, 2)))*(-1 + 2 * uz)) / 4.;
						mfcba = (-Mom222 - Mom022*ux*(1 + ux) - Mom122*(1 + 2 * ux) - 2 * (Mom112 + Mom212 + 2 * Mom112*ux + Mom012*ux*(1 + ux))*uy -
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*(-1 + pow(uy, 2)) +
							(-Mom220 - Mom020*ux*(1 + ux) - Mom120*(1 + 2 * ux) - 2 * (Mom110 + Mom210 + 2 * Mom110*ux + Mom010*ux*(1 + ux))*uy -
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*(-1 + pow(uy, 2)))*(-1 + uz)*uz +
								(-Mom221 - Mom021*ux*(1 + ux) - Mom121*(1 + 2 * ux) - 2 * (Mom111 + Mom211 + 2 * Mom111*ux + Mom011*ux*(1 + ux))*uy -
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*(-1 + pow(uy, 2)))*(-1 + 2 * uz)) / 4.;
						mfbaa = (Mom022 - Mom222 - 2 * Mom122*ux - Mom022*pow(ux, 2) - (Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*(-1 + uy)*uy -
							(Mom212 + 2 * Mom112*ux + Mom012*(-1 + pow(ux, 2)))*(-1 + 2 * uy) +
							(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*(-1 + uy)*uy -
								(Mom210 + 2 * Mom110*ux + Mom010*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*(-1 + uz)*uz +
								(Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*(-1 + uy)*uy -
									(Mom211 + 2 * Mom111*ux + Mom011*(-1 + pow(ux, 2)))*(-1 + 2 * uy))*(-1 + 2 * uz)) / 4.;
						mfbca = (Mom022 - Mom222 - 2 * Mom122*ux - Mom022*pow(ux, 2) - (Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*uy*(1 + uy) -
							(Mom212 + 2 * Mom112*ux + Mom012*(-1 + pow(ux, 2)))*(1 + 2 * uy) +
							(Mom020 - Mom220 - 2 * Mom120*ux - Mom020*pow(ux, 2) -
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*uy*(1 + uy) -
								(Mom210 + 2 * Mom110*ux + Mom010*(-1 + pow(ux, 2)))*(1 + 2 * uy))*(-1 + uz)*uz +
								(Mom021 - Mom221 - 2 * Mom121*ux - Mom021*pow(ux, 2) -
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*uy*(1 + uy) -
									(Mom211 + 2 * Mom111*ux + Mom011*(-1 + pow(ux, 2)))*(1 + 2 * uy))*(-1 + 2 * uz)) / 4.;
						mfaaa = (Mom222 + Mom022*(-1 + ux)*ux + Mom122*(-1 + 2 * ux) + (Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*(-1 + uy)*uy +
							(Mom212 + Mom012*(-1 + ux)*ux + Mom112*(-1 + 2 * ux))*(-1 + 2 * uy) +
							(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*(-1 + uy)*uy +
								(Mom210 + Mom010*(-1 + ux)*ux + Mom110*(-1 + 2 * ux))*(-1 + 2 * uy))*(-1 + uz)*uz +
								(Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*(-1 + uy)*uy +
									(Mom211 + Mom011*(-1 + ux)*ux + Mom111*(-1 + 2 * ux))*(-1 + 2 * uy))*(-1 + 2 * uz)) / 8.;
						mfcaa = (Mom122 + Mom222 + 2 * Mom122*ux + Mom022*ux*(1 + ux) +
							(Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*(-1 + uy)*uy +
							(Mom112 + Mom212 + 2 * Mom112*ux + Mom012*ux*(1 + ux))*(-1 + 2 * uy) +
							(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*(-1 + uy)*uy +
								(Mom110 + Mom210 + 2 * Mom110*ux + Mom010*ux*(1 + ux))*(-1 + 2 * uy))*(-1 + uz)*uz +
								(Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*(-1 + uy)*uy +
									(Mom111 + Mom211 + 2 * Mom111*ux + Mom011*ux*(1 + ux))*(-1 + 2 * uy))*(-1 + 2 * uz)) / 8.;
						mfaca = (Mom222 + Mom022*(-1 + ux)*ux + Mom122*(-1 + 2 * ux) + (Mom202 + Mom002*(-1 + ux)*ux + Mom102*(-1 + 2 * ux))*uy*(1 + uy) +
							(Mom212 + Mom012*(-1 + ux)*ux + Mom112*(-1 + 2 * ux))*(1 + 2 * uy) +
							(Mom220 + Mom020*(-1 + ux)*ux + Mom120*(-1 + 2 * ux) +
							(Mom200 + Mom000*(-1 + ux)*ux + Mom100*(-1 + 2 * ux))*uy*(1 + uy) +
								(Mom210 + Mom010*(-1 + ux)*ux + Mom110*(-1 + 2 * ux))*(1 + 2 * uy))*(-1 + uz)*uz +
								(Mom221 + Mom021*(-1 + ux)*ux + Mom121*(-1 + 2 * ux) +
							(Mom201 + Mom001*(-1 + ux)*ux + Mom101*(-1 + 2 * ux))*uy*(1 + uy) +
									(Mom211 + Mom011*(-1 + ux)*ux + Mom111*(-1 + 2 * ux))*(1 + 2 * uy))*(-1 + 2 * uz)) / 8.;
						mfcca = (Mom122 + Mom222 + 2 * Mom122*ux + Mom022*ux*(1 + ux) + (Mom102 + Mom202 + 2 * Mom102*ux + Mom002*ux*(1 + ux))*uy*(1 + uy) +
							(Mom112 + Mom212 + 2 * Mom112*ux + Mom012*ux*(1 + ux))*(1 + 2 * uy) +
							(Mom120 + Mom220 + 2 * Mom120*ux + Mom020*ux*(1 + ux) +
							(Mom100 + Mom200 + 2 * Mom100*ux + Mom000*ux*(1 + ux))*uy*(1 + uy) +
								(Mom110 + Mom210 + 2 * Mom110*ux + Mom010*ux*(1 + ux))*(1 + 2 * uy))*(-1 + uz)*uz +
								(Mom121 + Mom221 + 2 * Mom121*ux + Mom021*ux*(1 + ux) +
							(Mom101 + Mom201 + 2 * Mom101*ux + Mom001*ux*(1 + ux))*uy*(1 + uy) +
									(Mom111 + Mom211 + 2 * Mom111*ux + Mom011*ux*(1 + ux))*(1 + 2 * uy))*(-1 + 2 * uz)) / 8.;

						mfbbb = Mom022 - Mom222 - 2 * Mom122*ux - Mom022*pow(ux, 2) - 2 * (Mom212 + 2 * Mom112*ux + Mom012*(-1 + pow(ux, 2)))*uy -
							(Mom202 + 2 * Mom102*ux + Mom002*(-1 + pow(ux, 2)))*(-1 + pow(uy, 2)) -
							2 * (Mom221 + 2 * Mom121*ux + Mom021*(-1 + pow(ux, 2)) + 2 * (Mom211 + 2 * Mom111*ux + Mom011*(-1 + pow(ux, 2)))*uy +
							(Mom201 + 2 * Mom101*ux + Mom001*(-1 + pow(ux, 2)))*(-1 + pow(uy, 2)))*uz +
								(Mom220 + 2 * Mom120*ux + Mom020*(-1 + pow(ux, 2)) + 2 * (Mom210 + 2 * Mom110*ux + Mom010*(-1 + pow(ux, 2)))*uy +
							(Mom200 + 2 * Mom100*ux + Mom000*(-1 + pow(ux, 2)))*(-1 + pow(uy, 2)))*(1 - pow(uz, 2));

						////////////////////////////////////////////////////////////////////////////////////

						//////////////////////////////////////////////////////////////////////////
						//proof correctness
						//////////////////////////////////////////////////////////////////////////
//#ifdef  PROOF_CORRECTNESS
//						LBMReal drho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
//							+ (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
//							+ (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
						//UBLOG(logINFO, "lambda ="<<drho_post);
//						//LBMReal dif = fabs(rho - rho_post);
//						dif = drho - drho_post;
//#ifdef SINGLEPRECISION
//						if (dif > 10.0E-7 || dif < -10.0E-7)
//#else
//						if (dif > 10.0E-15 || dif < -10.0E-15)
//#endif
//						{
//							UB_THROW(UbException(UB_EXARGS, "Flocculation=" + UbSystem::toString(drho) + ", flocculation post=" + UbSystem::toString(drho_post)
//								+ " dif=" + UbSystem::toString(dif)
//								+ " Flocculation is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
//							//UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
//							//exit(EXIT_FAILURE);
//						}
//#endif
						//////////////////////////////////////////////////////////////////////////
						//write distribution
						//////////////////////////////////////////////////////////////////////////
						(*this->localDistributionsH)(D3Q27System::ET_E, x1, x2, x3) = mfabb;
						(*this->localDistributionsH)(D3Q27System::ET_N, x1, x2, x3) = mfbab;
						(*this->localDistributionsH)(D3Q27System::ET_T, x1, x2, x3) = mfbba;
						(*this->localDistributionsH)(D3Q27System::ET_NE, x1, x2, x3) = mfaab;
						(*this->localDistributionsH)(D3Q27System::ET_NW, x1p, x2, x3) = mfcab;
						(*this->localDistributionsH)(D3Q27System::ET_TE, x1, x2, x3) = mfaba;
						(*this->localDistributionsH)(D3Q27System::ET_TW, x1p, x2, x3) = mfcba;
						(*this->localDistributionsH)(D3Q27System::ET_TN, x1, x2, x3) = mfbaa;
						(*this->localDistributionsH)(D3Q27System::ET_TS, x1, x2p, x3) = mfbca;
						(*this->localDistributionsH)(D3Q27System::ET_TNE, x1, x2, x3) = mfaaa;
						(*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2, x3) = mfcaa;
						(*this->localDistributionsH)(D3Q27System::ET_TSE, x1, x2p, x3) = mfaca;
						(*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

						(*this->nonLocalDistributionsH)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac;
						(*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac;

						(*this->zeroDistributionsH)(x1, x2, x3) = mfbbb;
						//////////////////////////////////////////////////////////////////////////


					}
				}
			}
		}

	}
}
//////////////////////////////////////////////////////////////////////////
double ThixotropyExpLBMKernel::getCalculationTime()
{
	//return timer.getDuration();
	return timer.getTotalTime();
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyExpLBMKernel::setCollisionFactorF(double collFactor)
{
	setCollisionFactor(collFactor);
	this->collFactorF = collFactor;

}
//////////////////////////////////////////////////////////////////////////
void ThixotropyExpLBMKernel::setCollisionFactorH(double collFactor)
{
	this->collFactorH = collFactor;
}
//////////////////////////////////////////////////////////////////////////
double ThixotropyExpLBMKernel::getCollisionFactorF() const
{
	return this->collFactorF;
}
//////////////////////////////////////////////////////////////////////////
double ThixotropyExpLBMKernel::getCollisionFactorH() const
{
	return this->collFactorH;
}
void ThixotropyExpLBMKernel::setAlpha(double alpha)
{
	this->alpha = alpha;
}
double ThixotropyExpLBMKernel::getAlpha() const
{
	return this->alpha;
}
void ThixotropyExpLBMKernel::setTheta(double theta)
{
	this->theta = theta;
}
double ThixotropyExpLBMKernel::getTheta() const
{
	return this->theta;
}
void ThixotropyExpLBMKernel::swapDistributions()
{
	LBMKernel::swapDistributions();
	dataSet->getHdistributions()->swap();
}
//////////////////////////////////////////////////////////////////////////
