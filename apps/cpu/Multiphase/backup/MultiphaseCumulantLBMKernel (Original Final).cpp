#include "MultiphaseCumulantLBMKernel.h"
#include "D3Q27System.h"
#include "Interpolator.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <math.h>
#include <omp.h>

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
MultiphaseCumulantLBMKernel::MultiphaseCumulantLBMKernel()
{
   this->nx1 = 0;
   this->nx2 = 0;
   this->nx3 = 0;
   this->parameter = NORMAL;
   this->OxyyMxzz = 1.0;
   this->compressible = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseCumulantLBMKernel::MultiphaseCumulantLBMKernel(int nx1, int nx2, int nx3, Parameter p) 
{
   this->nx1 = nx1;
   this->nx2 = nx2;
   this->nx3 = nx3;
   parameter = p;
   this->compressible = false;
}
//////////////////////////////////////////////////////////////////////////
MultiphaseCumulantLBMKernel::~MultiphaseCumulantLBMKernel(void)
{

}
//////////////////////////////////////////////////////////////////////////
void MultiphaseCumulantLBMKernel::init()
{
   //DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+ghostLayerWitdh*2, nx2+ghostLayerWitdh*2, nx3+ghostLayerWitdh*2, -999.0));
   DistributionArray3DPtr f(new D3Q27EsoTwist3DSplittedVector(nx1+2, nx2+2, nx3+2, -999.0));
   DistributionArray3DPtr h(new D3Q27EsoTwist3DSplittedVector(nx1+2, nx2+2, nx3+2, -999.0)); // For phase-field
   PhaseFieldArray3DPtr divU(new CbArray3D<LBMReal,IndexerX3X2X1>(nx1+2, nx2+2, nx3+2, 0.0));
   dataSet->setFdistributions(f);
   dataSet->setHdistributions(h); // For phase-field
   dataSet->setPhaseField(divU);
}
//////////////////////////////////////////////////////////////////////////
LBMKernelPtr MultiphaseCumulantLBMKernel::clone()
{
   LBMKernelPtr kernel(new MultiphaseCumulantLBMKernel(nx1, nx2, nx3, parameter));
   boost::dynamic_pointer_cast<MultiphaseCumulantLBMKernel>(kernel)->init();
   
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
   switch (parameter)
   {
   case NORMAL:
      boost::dynamic_pointer_cast<MultiphaseCumulantLBMKernel>(kernel)->OxyyMxzz = 1.0;
   	break;
   case MAGIC:
      boost::dynamic_pointer_cast<MultiphaseCumulantLBMKernel>(kernel)->OxyyMxzz = 2.0 +(-collFactor);
      break;
   }
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseCumulantLBMKernel::calculate()
{
   timer.resetAndStart();
   collideAll();
   timer.stop();
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseCumulantLBMKernel::collideAll()
{
   using namespace D3Q27System;

   //initializing of forcing stuff 
   /*if (withForcing)
   {
   muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
   muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
   muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);

   muDeltaT = deltaT;

   muForcingX1.DefineVar("dt",&muDeltaT);
   muForcingX2.DefineVar("dt",&muDeltaT);
   muForcingX3.DefineVar("dt",&muDeltaT);

   muNu = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

   muForcingX1.DefineVar("nu",&muNu);
   muForcingX2.DefineVar("nu",&muNu);
   muForcingX3.DefineVar("nu",&muNu);

   LBMReal forcingX1 = 0;
   LBMReal forcingX2 = 0;
   LBMReal forcingX3 = 0;
   }*/
   forcingX1 = 0.0;
   forcingX2 = 0.0;
   forcingX3 = 0.0;
   /////////////////////////////////////

   localDistributionsF = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributionsF = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributionsF = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   localDistributionsH = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getLocalDistributions();
   nonLocalDistributionsH = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getNonLocalDistributions();
   zeroDistributionsH = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getHdistributions())->getZeroDistributions();

   

   //phaseField = dataSet->getPhaseField();

   BCArray3DPtr bcArray = this->getBCProcessor()->getBCArray();

   

   const int bcArrayMaxX1 = (int)bcArray->getNX1();
   const int bcArrayMaxX2 = (int)bcArray->getNX2();
   const int bcArrayMaxX3 = (int)bcArray->getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;


//#pragma omp parallel num_threads(8)
   {
   //   int i = omp_get_thread_num();
   //   printf_s("Hello from thread %d\n", i);
   //}
//#pragma omp for 

   
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr phaseField(new CbArray3D<LBMReal,IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, -999.0));
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr divU(new CbArray3D<LBMReal,IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, 0.0));
   
   //CbArray3D<LBMReal> phaseField(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3,-999);
   

   for(int x3 = 0; x3 <= maxX3; x3++)
   {
      for(int x2 = 0; x2 <= maxX2; x2++)
      {
         for(int x1 = 0; x1 <= maxX1; x1++)
         {
            if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
            {
				int x1p = x1 + 1;
				int x2p = x2 + 1;
				int x3p = x3 + 1;

				LBMReal mfcbb = (*this->localDistributionsH)(D3Q27System::ET_E, x1,x2,x3);
				LBMReal mfbcb = (*this->localDistributionsH)(D3Q27System::ET_N,x1,x2,x3); 
				LBMReal mfbbc = (*this->localDistributionsH)(D3Q27System::ET_T,x1,x2,x3);
				LBMReal mfccb = (*this->localDistributionsH)(D3Q27System::ET_NE,x1,x2,x3);
				LBMReal mfacb = (*this->localDistributionsH)(D3Q27System::ET_NW,x1p,x2,x3);
				LBMReal mfcbc = (*this->localDistributionsH)(D3Q27System::ET_TE,x1,x2,x3);
				LBMReal mfabc = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p,x2,x3);
				LBMReal mfbcc = (*this->localDistributionsH)(D3Q27System::ET_TN,x1,x2,x3);
				LBMReal mfbac = (*this->localDistributionsH)(D3Q27System::ET_TS,x1,x2p,x3);
				LBMReal mfccc = (*this->localDistributionsH)(D3Q27System::ET_TNE,x1,x2,x3);
				LBMReal mfacc = (*this->localDistributionsH)(D3Q27System::ET_TNW,x1p,x2,x3);
				LBMReal mfcac = (*this->localDistributionsH)(D3Q27System::ET_TSE,x1,x2p,x3);
				LBMReal mfaac = (*this->localDistributionsH)(D3Q27System::ET_TSW,x1p,x2p,x3);
				LBMReal mfabb = (*this->nonLocalDistributionsH)(D3Q27System::ET_W,x1p,x2,x3  );
				LBMReal mfbab = (*this->nonLocalDistributionsH)(D3Q27System::ET_S,x1,x2p,x3  );
				LBMReal mfbba = (*this->nonLocalDistributionsH)(D3Q27System::ET_B,x1,x2,x3p  );
				LBMReal mfaab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,x1p,x2p,x3 );
				LBMReal mfcab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,x1,x2p,x3 );
				LBMReal mfaba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,x1p,x2,x3p );
				LBMReal mfcba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,x1,x2,x3p );
				LBMReal mfbaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,x1,x2p,x3p );
				LBMReal mfbca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,x1,x2,x3p );
				LBMReal mfaaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW,x1p,x2p,x3p);
				LBMReal mfcaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE,x1,x2p,x3p);
				LBMReal mfaca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW,x1p,x2,x3p);
				LBMReal mfcca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE,x1,x2,x3p);

				LBMReal mfbbb = (*this->zeroDistributionsH)(x1,x2,x3);
				//LBMReal phase = h[ZERO] + h[E] + h[W] + h[N] + h[S] + h[T] + h[B] + h[NE] + h[SW] + h[SE] + h[NW] + h[TE] + h[BW] + 
				//	h[BE] + h[TW] + h[TN] + h[BS] + h[BN] + h[TS] + h[TNE] + h[TNW] + h[TSE] + h[TSW] + h[BNE] + h[BNW] + h[BSE] + h[BSW];
				//if (phase > 1.0) phase = 1.0e0;
				//(*phaseField)(x1,x2,x3) = phase;
				(*phaseField)(x1,x2,x3) = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
					+(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
					+(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;
			}
		 }
	  }
   }
   
   
   LBMReal collFactorM;
   LBMReal forcingTerm[D3Q27System::ENDF+1];
   LBMReal m000, m100, m010, m001, m110, m101, m011, m200, m020, m002, m120, m102, m210, m012, m201, m021, m111, m220, m202, m022, m211, m121, m112, m221, m212, m122, m222;
   LBMReal k000, k100, k010, k001, k110, k101, k011, k200, k020, k002, k120, k102, k210, k012, k201, k021, k111, k220, k202, k022, k211, k121, k112, k221, k212, k122, k222;
   LBMReal c000, c100, c010, c001, c110, c101, c011, c200, c020, c002, c120, c102, c210, c012, c201, c021, c111, c220, c202, c022, c211, c121, c112, c221, c212, c122, c222;

   LBMReal k200_pl_k020_pl_k002, k200_mi_k020, k200_mi_k002, k210_pl_k012, k210_mi_k012, k201_pl_k021, k201_mi_k021, k120_pl_k102, k120_mi_k102, k220_pl_k202_pl_k022, 
	   k220_mi2_k202_pl_k022, k220_pl_k202_mi2_k022;

   LBMReal c200_pl_c020_pl_c002, c200_mi_c020, c200_mi_c002, c210_pl_c012, c210_mi_c012, c201_pl_c021, c201_mi_c021, c120_pl_c102, c120_mi_c102, c220_pl_c202_pl_c022, 
	   c220_mi2_c202_pl_c022, c220_pl_c202_mi2_c022;

   LBMReal w1, w2, w3, w4, w5, w6, w7, w8, w9, w10;
   
   w2  = 1.0;
   w3  = 1.0;
   w4  = 1.0;
   w5  = 1.0;
   w6  = 1.0;
   w7  = 1.0;
   w8  = 1.0;
   w9  = 1.0;
   w10 = 1.0;

   for(int x3 = minX3; x3 < maxX3; x3++)
   {
      for(int x2 = minX2; x2 < maxX2; x2++)
      {
         for(int x1 = minX1; x1 < maxX1; x1++)
         {
            if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
            {
               int x1p = x1 + 1;
               int x2p = x2 + 1;
               int x3p = x3 + 1;


               //////////////////////////////////////////////////////////////////////////
               //Read distributions and phase field
               ////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////

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
			   
			   /*
			   phi[ZERO] = (phaseField)(x1,x2,x3);
			   phi[E  ] = (phaseField)(x1 + DX1[E  ], x2 + DX2[E  ], x3 + DX3[E  ]);
			   phi[N  ] = (phaseField)(x1 + DX1[N  ], x2 + DX2[N  ], x3 + DX3[N  ]);
			   phi[T  ] = (phaseField)(x1 + DX1[T  ], x2 + DX2[T  ], x3 + DX3[T  ]);
			   phi[W  ] = (phaseField)(x1 + DX1[W  ], x2 + DX2[W  ], x3 + DX3[W  ]);
			   phi[S  ] = (phaseField)(x1 + DX1[S  ], x2 + DX2[S  ], x3 + DX3[S  ]);
			   phi[B  ] = (phaseField)(x1 + DX1[B  ], x2 + DX2[B  ], x3 + DX3[B  ]);
			   phi[NE ] = (phaseField)(x1 + DX1[NE ], x2 + DX2[NE ], x3 + DX3[NE ]);
			   phi[NW ] = (phaseField)(x1 + DX1[NW ], x2 + DX2[NW ], x3 + DX3[NW ]);
			   phi[TE ] = (phaseField)(x1 + DX1[TE ], x2 + DX2[TE ], x3 + DX3[TE ]);
			   phi[TW ] = (phaseField)(x1 + DX1[TW ], x2 + DX2[TW ], x3 + DX3[TW ]);
			   phi[TN ] = (phaseField)(x1 + DX1[TN ], x2 + DX2[TN ], x3 + DX3[TN ]);
			   phi[TS ] = (phaseField)(x1 + DX1[TS ], x2 + DX2[TS ], x3 + DX3[TS ]);
			   phi[SW ] = (phaseField)(x1 + DX1[SW ], x2 + DX2[SW ], x3 + DX3[SW ]);
			   phi[SE ] = (phaseField)(x1 + DX1[SE ], x2 + DX2[SE ], x3 + DX3[SE ]);
			   phi[BW ] = (phaseField)(x1 + DX1[BW ], x2 + DX2[BW ], x3 + DX3[BW ]);
			   phi[BE ] = (phaseField)(x1 + DX1[BE ], x2 + DX2[BE ], x3 + DX3[BE ]);
			   phi[BS ] = (phaseField)(x1 + DX1[BS ], x2 + DX2[BS ], x3 + DX3[BS ]);
			   phi[BN ] = (phaseField)(x1 + DX1[BN ], x2 + DX2[BN ], x3 + DX3[BN ]);
			   phi[BSW] = (phaseField)(x1 + DX1[BSW], x2 + DX2[BSW], x3 + DX3[BSW]);
			   phi[BSE] = (phaseField)(x1 + DX1[BSE], x2 + DX2[BSE], x3 + DX3[BSE]);
			   phi[BNW] = (phaseField)(x1 + DX1[BNW], x2 + DX2[BNW], x3 + DX3[BNW]);
			   phi[BNE] = (phaseField)(x1 + DX1[BNE], x2 + DX2[BNE], x3 + DX3[BNE]);
			   phi[TNE] = (phaseField)(x1 + DX1[TNE], x2 + DX2[TNE], x3 + DX3[TNE]);
			   phi[TNW] = (phaseField)(x1 + DX1[TNW], x2 + DX2[TNW], x3 + DX3[TNW]);
			   phi[TSE] = (phaseField)(x1 + DX1[TSE], x2 + DX2[TSE], x3 + DX3[TSE]);
			   phi[TSW] = (phaseField)(x1 + DX1[TSW], x2 + DX2[TSW], x3 + DX3[TSW]);
			   */
			   findNeighbors(phaseField, x1, x2, x3);

			   LBMReal mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1,x2,x3);
			   LBMReal mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N,x1,x2,x3); 
			   LBMReal mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T,x1,x2,x3);
			   LBMReal mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE,x1,x2,x3);
			   LBMReal mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW,x1p,x2,x3);
			   LBMReal mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE,x1,x2,x3);
			   LBMReal mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p,x2,x3);
			   LBMReal mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN,x1,x2,x3);
			   LBMReal mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS,x1,x2p,x3);
			   LBMReal mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE,x1,x2,x3);
			   LBMReal mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW,x1p,x2,x3);
			   LBMReal mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE,x1,x2p,x3);
			   LBMReal mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW,x1p,x2p,x3);
			   LBMReal mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W,x1p,x2,x3  );
			   LBMReal mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S,x1,x2p,x3  );
			   LBMReal mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B,x1,x2,x3p  );
			   LBMReal mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW,x1p,x2p,x3 );
			   LBMReal mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE,x1,x2p,x3 );
			   LBMReal mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW,x1p,x2,x3p );
			   LBMReal mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE,x1,x2,x3p );
			   LBMReal mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS,x1,x2p,x3p );
			   LBMReal mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN,x1,x2,x3p );
			   LBMReal mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW,x1p,x2p,x3p);
			   LBMReal mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE,x1,x2p,x3p);
			   LBMReal mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW,x1p,x2,x3p);
			   LBMReal mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE,x1,x2,x3p);

			   LBMReal mfbbb = (*this->zeroDistributionsF)(x1,x2,x3);

			   

			   LBMReal rhoH = 1.0;
			   LBMReal rhoL = 1.0/densityRatio;

			   //LBMReal rhoToPhi = (1.0 - 1.0/densityRatio);
			   LBMReal rhoToPhi = (rhoH - rhoL)/(phiH - phiL);

			   //collFactorM = phi[ZERO]*collFactorL + (1-phi[ZERO])*collFactorG;
			   //collFactorM = phi[ZERO]*collFactorG + (1-phi[ZERO])*collFactorL;
			   
			   //LBMReal tauH = 1.0;
			   LBMReal di = sqrt(8*kappa/beta);
			   
			   LBMReal dX1_phi = gradX1_phi();
			   LBMReal dX2_phi = gradX2_phi();
			   LBMReal dX3_phi = gradX3_phi();
			   
			   LBMReal denom = sqrt(dX1_phi*dX1_phi + dX2_phi*dX2_phi + dX3_phi*dX3_phi) + 1e-9;
			   LBMReal normX1 = dX1_phi/denom;
			   LBMReal normX2 = dX2_phi/denom;
			   LBMReal normX3 = dX3_phi/denom;

			   collFactorM = collFactorL + (collFactorL - collFactorG)*(phi[ZERO] - phiH)/(phiH - phiL);
			   
			   /*if ((phi[ZERO] > 0.1)||(phi[ZERO] < 0.9))
			   {
				   collFactorM*=(1.0-denom);
			   }*/

			   w1 = collFactorM;
			   
			   /*dX1_phi = -normX1*((phi[ZERO]>phiH || phi[ZERO]<phiL) ? 0.0 : 4*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)/di);
			   dX2_phi = -normX2*((phi[ZERO]>phiH || phi[ZERO]<phiL) ? 0.0 : 4*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)/di);
			   dX3_phi = -normX3*((phi[ZERO]>phiH || phi[ZERO]<phiL) ? 0.0 : 4*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)/di);*/

			   //UbTupleDouble3 coords = grid->getNodeCoordinates(block, x1, x2, x3);
			   /*Block3D bl = this->block();
			   
			   int wX1 = bl->getX1()  + x1;
			   int wX2 = bl->getX2()  + x2;
			   int wX3 = bl->getX3()  + x3;*/
			   
			   /*if (wX3 >= 30.0)
			   {
			   dX1_phi = 0.0;
			   dX2_phi = 0.0;
			   dX3_phi = 0.0;
			   }*/


			   LBMReal mu = 2*beta*phi[ZERO]*(phi[ZERO]-1)*(2*phi[ZERO]-1) - kappa*nabla2_phi();
			   
			   //LBMReal rhoToPhi = (1.0/densityRatio - 1.0);
			   
			   			   

			   //----------- Calculating Macroscopic Values -------------

			   //LBMReal rho = phi[ZERO] + (1.0 - phi[ZERO])*1.0/densityRatio;
			   LBMReal rho = rhoH + rhoToPhi*(phi[ZERO] - phiH);
			   //LBMReal rho = phi[ZERO]*1.0/densityRatio + (1.0 - phi[ZERO]);

			   if (withForcing)
			   {
				   //muX1 = static_cast<double>(x1-1+ix1*maxX1);
				   //muX2 = static_cast<double>(x2-1+ix2*maxX2);
				   //muX3 = static_cast<double>(x3-1+ix3*maxX3);

				   forcingX1 = muForcingX1.Eval();
				   forcingX2 = muForcingX2.Eval();
				   forcingX3 = muForcingX3.Eval();

				   LBMReal rho_m = 1.0/densityRatio;
				   forcingX1 = forcingX1*(rho-rho_m);
				   forcingX2 = forcingX2*(rho-rho_m);
				   forcingX3 = forcingX3*(rho-rho_m);

				   //ux += forcingX1*deltaT*0.5; // X
				   //uy += forcingX2*deltaT*0.5; // Y
				   //uz += forcingX3*deltaT*0.5; // Z
			   }

			   LBMReal ux = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
				   (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
				   (mfcbb-mfabb)) / (rho*c1o3) + (mu*dX1_phi + forcingX1)/(2*rho);

			   LBMReal uy = ((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
				   (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				   (mfbcb-mfbab)) / (rho*c1o3) + (mu*dX2_phi + forcingX2)/(2*rho);

			   LBMReal uz = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
				   (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				   (mfbbc-mfbba)) / (rho*c1o3) + (mu*dX3_phi + forcingX3)/(2*rho);

			   LBMReal p1 = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
				   +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
				   +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb +
			   (ux*rhoToPhi*dX1_phi*c1o3 + uy*rhoToPhi*dX2_phi*c1o3 + uz*rhoToPhi*dX3_phi*c1o3)/2.0;
			   
			   //vvx = 0.0; vvy = 0.0; vvz = 0.0;
			   //--------------------------------------------------------
			   



			   LBMReal ux2 = ux*ux;
			   LBMReal uy2 = uy*uy;
			   LBMReal uz2 = uz*uz;
			   LBMReal ux_uy = ux*uy;
			   LBMReal ux_uz = ux*uz;
			   LBMReal uy_uz = uy*uz;
			   LBMReal ux_uy_uz = ux*uy*uz;


/*
			   //----------- Calculating Forcing Terms -------------
			   LBMReal forcingTerm1 = (ux*mu*dX1_phi + uy*mu*dX2_phi + uz*mu*dX3_phi);
			   for (int dir = STARTF; dir < (ENDF+1); dir++)
			   {
				   if (dir != ZERO)
				   {
					   LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]])/2.0;
					   forcingTerm[dir] = (c1o3*rhoToPhi*dirGrad_phi + mu*dirGrad_phi)*(DX1[dir]*ux + DX2[dir]*uy + DX3[dir]*uz)*WEIGTH[dir]/c1o3 + mu*dirGrad_phi*WEIGTH[dir] - 
						   (forcingTerm1)*WEIGTH[dir];
				   } 
				   else
				   {
					   forcingTerm[ZERO] =  -(forcingTerm1)*WEIGTH[ZERO];
				   }
			   }
			  //--------------------------------------------------------
*/

			   //----------- Calculating Forcing Terms * -------------
			   //LBMReal forcingTerm1 = (ux*mu*dX1_phi + uy*mu*dX2_phi + uz*mu*dX3_phi);
			   for (int dir = STARTF; dir <= (FENDDIR); dir++)
			   {
				   LBMReal velProd = DX1[dir]*ux + DX2[dir]*uy + DX3[dir]*uz;
				   LBMReal velSq1 = velProd*velProd;
				   LBMReal gamma = WEIGTH[dir]*(1.0 + 3*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2));
				   
				   //forcingTerm[dir] = (DX1[dir] - ux)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX1_phi + gamma*mu*dX1_phi) + 
					//   (DX2[dir] - uy)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX2_phi + gamma*mu*dX2_phi) + 
					//   (DX3[dir] - uz)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX3_phi + gamma*mu*dX3_phi);
				   
				   LBMReal fac1 = (gamma - WEIGTH[dir])*c1o3*rhoToPhi;
				   LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]])/2.0;
				   //LBMReal dirGrad_phi = DX1[dir]*dX1_phi + DX2[dir]*dX2_phi + DX3[dir]*dX3_phi;
				   
				   /*forcingTerm[dir] =  (- (ux)*(fac1*dX1_phi + gamma*mu*dX1_phi) - 
				   (uy)*(fac1*dX2_phi + gamma*mu*dX2_phi) - 
				   (uz)*(fac1*dX3_phi + gamma*mu*dX3_phi)) + (fac1*dirGrad_phi + gamma*mu*dirGrad_phi + DX1[dir]*forcingX1 + DX2[dir]*forcingX2 + DX3[dir]*forcingX3);*/
				   
				   
				   forcingTerm[dir] = ((-ux)*(fac1*dX1_phi + gamma*(mu*dX1_phi + forcingX1)) +
					   				  (-uy)*(fac1*dX2_phi + gamma*(mu*dX2_phi + forcingX2)) +
					   				  (-uz)*(fac1*dX3_phi + gamma*(mu*dX3_phi + forcingX3))) +
									  (DX1[dir])*(fac1*dX1_phi + gamma*(mu*dX1_phi + forcingX1)) +
									  (DX2[dir])*(fac1*dX2_phi + gamma*(mu*dX2_phi + forcingX2)) +
									  (DX3[dir])*(fac1*dX3_phi + gamma*(mu*dX3_phi + forcingX3));

			   }

			   LBMReal gamma = WEIGTH[ZERO]*(1.0 - 1.5*(ux2+uy2+uz2));
			   /*forcingTerm[ZERO] = -(ux)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX1_phi + gamma*mu*dX1_phi) - 
			   (uy)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX2_phi + gamma*mu*dX2_phi) - 
			   (uz)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX3_phi + gamma*mu*dX3_phi);*/
			   LBMReal fac1 = (gamma - WEIGTH[ZERO])*c1o3*rhoToPhi;
			   forcingTerm[ZERO] = (-ux)*(fac1*dX1_phi + gamma*(mu*dX1_phi + forcingX1)) +
				   (-uy)*(fac1*dX2_phi + gamma*(mu*dX2_phi + forcingX2)) +
				   (-uz)*(fac1*dX3_phi + gamma*(mu*dX3_phi + forcingX3));

			   //--------------------------------------------------------

/*			   
			   f1[E  ] = (g[E  ] + 0.5*forcingTerm[E  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[E  ]/c1o3;
			   f1[N  ] = (g[N  ] + 0.5*forcingTerm[N  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[N  ]/c1o3;
			   f1[T  ] = (g[T  ] + 0.5*forcingTerm[T  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[T  ]/c1o3;
			   f1[NE ] = (g[NE ] + 0.5*forcingTerm[NE ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[NE ]/c1o3;
			   f1[NW ] = (g[NW ] + 0.5*forcingTerm[NW ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[NW ]/c1o3;
			   f1[TE ] = (g[TE ] + 0.5*forcingTerm[TE ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TE ]/c1o3;
			   f1[TW ] = (g[TW ] + 0.5*forcingTerm[TW ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TW ]/c1o3;
			   f1[TN ] = (g[TN ] + 0.5*forcingTerm[TN ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TN ]/c1o3;
			   f1[TS ] = (g[TS ] + 0.5*forcingTerm[TS ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TS ]/c1o3;
			   f1[TNE] = (g[TNE] + 0.5*forcingTerm[TNE])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TNE]/c1o3;
			   f1[TNW] = (g[TNW] + 0.5*forcingTerm[TNW])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TNW]/c1o3;
			   f1[TSE] = (g[TSE] + 0.5*forcingTerm[TSE])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TSE]/c1o3;
			   f1[TSW] = (g[TSW] + 0.5*forcingTerm[TSW])/c1o3 - (p1 - rho*c1o3)*WEIGTH[TSW]/c1o3;
			   f1[W  ] = (g[W  ] + 0.5*forcingTerm[W  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[W  ]/c1o3;
			   f1[S  ] = (g[S  ] + 0.5*forcingTerm[S  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[S  ]/c1o3;
			   f1[B  ] = (g[B  ] + 0.5*forcingTerm[B  ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[B  ]/c1o3;
			   f1[SW ] = (g[SW ] + 0.5*forcingTerm[SW ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[SW ]/c1o3;
			   f1[SE ] = (g[SE ] + 0.5*forcingTerm[SE ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[SE ]/c1o3;
			   f1[BW ] = (g[BW ] + 0.5*forcingTerm[BW ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BW ]/c1o3;
			   f1[BE ] = (g[BE ] + 0.5*forcingTerm[BE ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BE ]/c1o3;
			   f1[BS ] = (g[BS ] + 0.5*forcingTerm[BS ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BS ]/c1o3;
			   f1[BN ] = (g[BN ] + 0.5*forcingTerm[BN ])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BN ]/c1o3;
			   f1[BSW] = (g[BSW] + 0.5*forcingTerm[BSW])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BSW]/c1o3;
			   f1[BSE] = (g[BSE] + 0.5*forcingTerm[BSE])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BSE]/c1o3;
			   f1[BNW] = (g[BNW] + 0.5*forcingTerm[BNW])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BNW]/c1o3;
			   f1[BNE] = (g[BNE] + 0.5*forcingTerm[BNE])/c1o3 - (p1 - rho*c1o3)*WEIGTH[BNE]/c1o3;
			   f1[ZERO] = (g[ZERO] + 0.5*forcingTerm[ZERO])/c1o3 - (p1 - rho*c1o3)*WEIGTH[ZERO]/c1o3;
*/
			   
			   mfcbb = 3.0*(mfcbb + 0.5*forcingTerm[E  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[E  ];
			   mfbcb = 3.0*(mfbcb + 0.5*forcingTerm[N  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[N  ];
			   mfbbc = 3.0*(mfbbc + 0.5*forcingTerm[T  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[T  ];
			   mfccb = 3.0*(mfccb + 0.5*forcingTerm[NE ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[NE ];
			   mfacb = 3.0*(mfacb + 0.5*forcingTerm[NW ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[NW ];
			   mfcbc = 3.0*(mfcbc + 0.5*forcingTerm[TE ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TE ];
			   mfabc = 3.0*(mfabc + 0.5*forcingTerm[TW ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TW ];
			   mfbcc = 3.0*(mfbcc + 0.5*forcingTerm[TN ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TN ];
			   mfbac = 3.0*(mfbac + 0.5*forcingTerm[TS ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TS ];
			   mfccc = 3.0*(mfccc + 0.5*forcingTerm[TNE] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TNE];
			   mfacc = 3.0*(mfacc + 0.5*forcingTerm[TNW] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TNW];
			   mfcac = 3.0*(mfcac + 0.5*forcingTerm[TSE] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TSE];
			   mfaac = 3.0*(mfaac + 0.5*forcingTerm[TSW] )/rho ;//-(3.0*p1 - rho)*WEIGTH[TSW];
			   mfabb = 3.0*(mfabb + 0.5*forcingTerm[W  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[W  ];
			   mfbab = 3.0*(mfbab + 0.5*forcingTerm[S  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[S  ];
			   mfbba = 3.0*(mfbba + 0.5*forcingTerm[B  ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[B  ];
			   mfaab = 3.0*(mfaab + 0.5*forcingTerm[SW ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[SW ];
			   mfcab = 3.0*(mfcab + 0.5*forcingTerm[SE ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[SE ];
			   mfaba = 3.0*(mfaba + 0.5*forcingTerm[BW ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BW ];
			   mfcba = 3.0*(mfcba + 0.5*forcingTerm[BE ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BE ];
			   mfbaa = 3.0*(mfbaa + 0.5*forcingTerm[BS ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BS ];
			   mfbca = 3.0*(mfbca + 0.5*forcingTerm[BN ] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BN ];
			   mfaaa = 3.0*(mfaaa + 0.5*forcingTerm[BSW] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BSW];
			   mfcaa = 3.0*(mfcaa + 0.5*forcingTerm[BSE] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BSE];
			   mfaca = 3.0*(mfaca + 0.5*forcingTerm[BNW] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BNW];
			   mfcca = 3.0*(mfcca + 0.5*forcingTerm[BNE] )/rho ;//-(3.0*p1 - rho)*WEIGTH[BNE];
			   mfbbb = 3.0*(mfbbb + 0.5*forcingTerm[ZERO])/rho ;//- (3.0*p1 - rho)*WEIGTH[ZERO];
			   
			   LBMReal rho1 = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
				   +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
				   +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;


			   /*
               //forcing 
               ///////////////////////////////////////////////////////////////////////////////////////////
               if (withForcing)
               {
                  muX1 = static_cast<double>(x1-1+ix1*maxX1);
                  muX2 = static_cast<double>(x2-1+ix2*maxX2);
                  muX3 = static_cast<double>(x3-1+ix3*maxX3);

                  forcingX1 = muForcingX1.Eval();
                  forcingX2 = muForcingX2.Eval();
                  forcingX3 = muForcingX3.Eval();

                  vvx += forcingX1*deltaT*0.5; // X
                  vvy += forcingX2*deltaT*0.5; // Y
                  vvz += forcingX3*deltaT*0.5; // Z
               }
               /////////////////////////////////////////////////////////////////////////////////////////// 
			   */
              
			   LBMReal oMdrho, m0, m1, m2;
               
			   oMdrho=mfccc+mfaaa;
               m0=mfaca+mfcac;
               m1=mfacc+mfcaa;
               m2=mfaac+mfcca;
               oMdrho+=m0;
               m1+=m2;
               oMdrho+=m1;
               m0=mfbac+mfbca;
               m1=mfbaa+mfbcc;
               m0+=m1;
               m1=mfabc+mfcba;
               m2=mfaba+mfcbc;
               m1+=m2;
               m0+=m1;
               m1=mfacb+mfcab;
               m2=mfaab+mfccb;
               m1+=m2;
               m0+=m1;
               oMdrho+=m0;
               m0=mfabb+mfcbb;
               m1=mfbab+mfbcb;
               m2=mfbba+mfbbc;
               m0+=m1+m2;
               m0+=mfbbb; //hat gefehlt
               oMdrho = 1. - (oMdrho + m0);
			   //oMdrho = rho - (oMdrho + m0);


               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal wadjust;
               LBMReal qudricLimit = 0.01;
               ////////////////////////////////////////////////////////////////////////////////////
               //Hin
               ////////////////////////////////////////////////////////////////////////////////////
               // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m2    = mfaaa + mfaac;
               m1    = mfaac - mfaaa;
               m0    = m2          + mfaab;
               mfaaa = m0;
               m0   += c1o36 * oMdrho;   
               mfaab = m1 -        m0 * uz;
               mfaac = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfabc;
               m1    = mfabc  - mfaba;
               m0    = m2          + mfabb;
               mfaba = m0;
               m0   += c1o9 * oMdrho;
               mfabb = m1 -        m0 * uz;
               mfabc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfacc;
               m1    = mfacc  - mfaca;
               m0    = m2          + mfacb;
               mfaca = m0;
               m0   += c1o36 * oMdrho;
               mfacb = m1 -        m0 * uz;
               mfacc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbac;
               m1    = mfbac - mfbaa;
               m0    = m2          + mfbab;
               mfbaa = m0;
               m0   += c1o9 * oMdrho;
               mfbab = m1 -        m0 * uz;
               mfbac = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbba  + mfbbc;
               m1    = mfbbc  - mfbba;
               m0    = m2          + mfbbb;
               mfbba = m0;
               m0   += c4o9 * oMdrho;
               mfbbb = m1 -        m0 * uz;
               mfbbc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbca  + mfbcc;
               m1    = mfbcc  - mfbca;
               m0    = m2          + mfbcb;
               mfbca = m0;
               m0   += c1o9 * oMdrho;
               mfbcb = m1 -        m0 * uz;
               mfbcc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcac;
               m1    = mfcac - mfcaa;
               m0    = m2          + mfcab;
               mfcaa = m0;
               m0   += c1o36 * oMdrho;
               mfcab = m1 -        m0 * uz;
               mfcac = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcba  + mfcbc;
               m1    = mfcbc  - mfcba;
               m0    = m2          + mfcbb;
               mfcba = m0;
               m0   += c1o9 * oMdrho;
               mfcbb = m1 -        m0 * uz;
               mfcbc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcca  + mfccc;
               m1    = mfccc  - mfcca;
               m0    = m2          + mfccb;
               mfcca = m0;
               m0   += c1o36 * oMdrho;
               mfccb = m1 -        m0 * uz;
               mfccc = m2 - 2. *   m1 * uz + uz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m2    = mfaaa + mfaca;
               m1    = mfaca - mfaaa;
               m0    = m2          + mfaba;
               mfaaa = m0;
               m0   += c1o6 * oMdrho;
               mfaba = m1 -        m0 * uy;
               mfaca = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab  + mfacb;
               m1    = mfacb  - mfaab;
               m0    = m2          + mfabb;
               mfaab = m0;
               mfabb = m1 -        m0 * uy;
               mfacb = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac  + mfacc;
               m1    = mfacc  - mfaac;
               m0    = m2          + mfabc;
               mfaac = m0;
               m0   += c1o18 * oMdrho;
               mfabc = m1 -        m0 * uy;
               mfacc = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbca;
               m1    = mfbca - mfbaa;
               m0    = m2          + mfbba;
               mfbaa = m0;
               m0   += c2o3 * oMdrho;
               mfbba = m1 -        m0 * uy;
               mfbca = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbab  + mfbcb;
               m1    = mfbcb  - mfbab;
               m0    = m2          + mfbbb;
               mfbab = m0;
               mfbbb = m1 -        m0 * uy;
               mfbcb = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbac  + mfbcc;
               m1    = mfbcc  - mfbac;
               m0    = m2          + mfbbc;
               mfbac = m0;
               m0   += c2o9 * oMdrho;
               mfbbc = m1 -        m0 * uy;
               mfbcc = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcca;
               m1    = mfcca - mfcaa;
               m0    = m2          + mfcba;
               mfcaa = m0;
               m0   += c1o6 * oMdrho;
               mfcba = m1 -        m0 * uy;
               mfcca = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcab  + mfccb;
               m1    = mfccb  - mfcab;
               m0    = m2          + mfcbb;
               mfcab = m0;
               mfcbb = m1 -        m0 * uy;
               mfccb = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcac  + mfccc;
               m1    = mfccc  - mfcac;
               m0    = m2          + mfcbc;
               mfcac = m0;
               m0   += c1o18 * oMdrho;
               mfcbc = m1 -        m0 * uy;
               mfccc = m2 - 2. *   m1 * uy + uy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m2    = mfaaa + mfcaa;
               m1    = mfcaa - mfaaa;
               m0    = m2          + mfbaa;
               mfaaa = m0;
               m0   += 1. * oMdrho;
               mfbaa = m1 -        m0 * ux;
               mfcaa = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfcba;
               m1    = mfcba  - mfaba;
               m0    = m2          + mfbba;
               mfaba = m0;
               mfbba = m1 -        m0 * ux;
               mfcba = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfcca;
               m1    = mfcca  - mfaca;
               m0    = m2          + mfbca;
               mfaca = m0;
               m0   += c1o3 * oMdrho;
               mfbca = m1 -        m0 * ux;
               mfcca = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab + mfcab;
               m1    = mfcab - mfaab;
               m0    = m2          + mfbab;
               mfaab = m0;
               mfbab = m1 -        m0 * ux;
               mfcab = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabb  + mfcbb;
               m1    = mfcbb  - mfabb;
               m0    = m2          + mfbbb;
               mfabb = m0;
               mfbbb = m1 -        m0 * ux;
               mfcbb = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacb  + mfccb;
               m1    = mfccb  - mfacb;
               m0    = m2          + mfbcb;
               mfacb = m0;
               mfbcb = m1 -        m0 * ux;
               mfccb = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac + mfcac;
               m1    = mfcac - mfaac;
               m0    = m2          + mfbac;
               mfaac = m0;
               m0   += c1o3 * oMdrho;
               mfbac = m1 -        m0 * ux;
               mfcac = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabc  + mfcbc;
               m1    = mfcbc  - mfabc;
               m0    = m2          + mfbbc;
               mfabc = m0;
               mfbbc = m1 -        m0 * ux;
               mfcbc = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacc  + mfccc;
               m1    = mfccc  - mfacc;
               m0    = m2          + mfbcc;
               mfacc = m0;
               m0   += c1o9 * oMdrho;
               mfbcc = m1 -        m0 * ux;
               mfccc = m2 - 2. *   m1 * ux + ux2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal OxxPyyPzz = 1.; //omega2 or bulk viscosity
               LBMReal OxyyPxzz  = 1.;//-s9;//2+s9;//
               //LBMReal OxyyMxzz  = 1.;//2+s9;//
               LBMReal O4        = 1.;
               LBMReal O5        = 1.;
               LBMReal O6        = 1.;

               //Cum 4.
               //LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
               //LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
               //LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

               LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 ) * mfabb + 2. * mfbba * mfbab);
               LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 ) * mfbab + 2. * mfbba * mfabb);
               LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 ) * mfbba + 2. * mfbab * mfabb);

			   LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho);
			   LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho);
			   LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho);

			   //LBMReal CUMcca = mfcca - ((mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(-p1/c1o3)*oMdrho);
			   //LBMReal CUMcac = mfcac - ((mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(-p1/c1o3)*oMdrho);
			   //LBMReal CUMacc = mfacc - ((mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(-p1/c1o3)*oMdrho);

               //Cum 5.
               LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
               LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
               LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

               //Cum 6.
               LBMReal CUMccc = mfccc  +((-4. *  mfbbb * mfbbb 
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                  -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
                  +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                  +  2. * (mfcaa * mfaca * mfaac)
                  + 16. *  mfbba * mfbab * mfabb)
                  - c1o3* (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
                  - c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(1.-2.* oMdrho)- c1o27* oMdrho * oMdrho*(-2.* oMdrho)
                  +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                  +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;

               //2.
               // linear combinations
               LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               LBMReal mxxMyy    = mfcaa - mfaca;
               LBMReal mxxMzz         = mfcaa - mfaac;

               LBMReal dxux = -c1o2 * collFactorM *(mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz*(mfaaa - mxxPyyPzz);
               LBMReal dyuy = dxux + collFactorM * c3o2 * mxxMyy;
               LBMReal dzuz = dxux + collFactorM * c3o2 * mxxMzz;

			   /*LBMReal Dxy =-three*collFactorM*mfbba;
			   LBMReal Dxz =-three*collFactorM*mfbab;
			   LBMReal Dyz =-three*collFactorM*mfabb;

			   LBMReal strainMag = sqrt(2*(dxux*dxux + dyuy*dyuy + dzuz*dzuz) + Dxy*Dxy + Dxz*Dxz + Dyz*Dyz);
			   LBMReal intVis = 3*abs(denom - 1e-9)*strainMag;
			   LBMReal fluidVis = (1.0/collFactorM - 0.5)/3.0;
			   collFactorM = 1.0/((fluidVis + intVis)*3.0 + 0.5);*/
			   (*divU)(x1,x2,x3) = dxux + dyuy + dzuz;
			   
			   
			   //relax
               mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (ux2 * dxux + uy2 * dyuy + uz2 * dzuz);
               mxxMyy    += collFactorM * (-mxxMyy) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uy2 * dyuy);
               mxxMzz    += collFactorM * (-mxxMzz) - 3. * (1. - c1o2 * collFactorM) * (ux2 * dxux - uz2 * dzuz);

               mfabb     += collFactorM * (-mfabb);
               mfbab     += collFactorM * (-mfbab);
               mfbba     += collFactorM * (-mfbba);

               // linear combinations back
               mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);

               //3.
               // linear combinations
               LBMReal mxxyPyzz = mfcba + mfabc;
               LBMReal mxxyMyzz = mfcba - mfabc;

               LBMReal mxxzPyyz = mfcab + mfacb;
               LBMReal mxxzMyyz = mfcab - mfacb;

               LBMReal mxyyPxzz = mfbca + mfbac;
               LBMReal mxyyMxzz = mfbca - mfbac;

               //relax
               wadjust    = OxyyMxzz+(1.-OxyyMxzz)*fabs(mfbbb)/(fabs(mfbbb)+qudricLimit);
               mfbbb     += wadjust * (-mfbbb);
               wadjust    = OxyyPxzz+(1.-OxyyPxzz)*fabs(mxxyPyzz)/(fabs(mxxyPyzz)+qudricLimit);
               mxxyPyzz  += wadjust * (-mxxyPyzz);
               wadjust    = OxyyMxzz+(1.-OxyyMxzz)*fabs(mxxyMyzz)/(fabs(mxxyMyzz)+qudricLimit);
               mxxyMyzz  += wadjust * (-mxxyMyzz);
               wadjust    = OxyyPxzz+(1.-OxyyPxzz)*fabs(mxxzPyyz)/(fabs(mxxzPyyz)+qudricLimit);
               mxxzPyyz  += wadjust * (-mxxzPyyz);
               wadjust    = OxyyMxzz+(1.-OxyyMxzz)*fabs(mxxzMyyz)/(fabs(mxxzMyyz)+qudricLimit);
               mxxzMyyz  += wadjust * (-mxxzMyyz);
               wadjust    = OxyyPxzz+(1.-OxyyPxzz)*fabs(mxyyPxzz)/(fabs(mxyyPxzz)+qudricLimit);
               mxyyPxzz  += wadjust * (-mxyyPxzz);
               wadjust    = OxyyMxzz+(1.-OxyyMxzz)*fabs(mxyyMxzz)/(fabs(mxyyMxzz)+qudricLimit);
               mxyyMxzz  += wadjust * (-mxyyMxzz);

               // linear combinations back
               mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
               mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
               mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
               mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
               mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
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

               mfcbb = CUMcbb + ((mfcaa + c1o3 ) * mfabb + 2. * mfbba * mfbab);
               mfbcb = CUMbcb + ((mfaca + c1o3 ) * mfbab + 2. * mfbba * mfabb);
               mfbbc = CUMbbc + ((mfaac + c1o3 ) * mfbba + 2. * mfbab * mfabb);

               mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

			   //mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(-p1/c1o3)*oMdrho;
			   //mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(-p1/c1o3)*oMdrho;
			   //mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(-p1/c1o3)*oMdrho;

               //5.
               mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
               mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
               mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

               //6.
               mfccc = CUMccc  -((-4. *  mfbbb * mfbbb 
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
                  -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
                  +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                  +  2. * (mfcaa * mfaca * mfaac)
                  + 16. *  mfbba * mfbab * mfabb)
                  - c1o3* (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
                  - c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(1.-2.* oMdrho)- c1o27* oMdrho * oMdrho*(-2.* oMdrho)
                  +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                  +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) -c1o27*oMdrho;

               ////////////////////////////////////////////////////////////////////////////////////
               //forcing
               mfbaa=-mfbaa;
               mfaba=-mfaba;
               mfaab=-mfaab;
               //////////////////////////////////////////////////////////////////////////////////////

               ////////////////////////////////////////////////////////////////////////////////////
               //back
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m0 =  mfaac * c1o2 +      mfaab * (uz - c1o2) + (mfaaa + 1. * oMdrho) * (     uz2 - uz) * c1o2;
               m1 = -mfaac        - 2. * mfaab *  uz         +  mfaaa                * (1. - uz2)              - 1. * oMdrho * uz2;
               m2 =  mfaac * c1o2 +      mfaab * (uz + c1o2) + (mfaaa + 1. * oMdrho) * (     uz2 + uz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfabc * c1o2 +      mfabb * (uz - c1o2) + mfaba * (     uz2 - uz) * c1o2;
               m1 = -mfabc        - 2. * mfabb *  uz         + mfaba * (1. - uz2);
               m2 =  mfabc * c1o2 +      mfabb * (uz + c1o2) + mfaba * (     uz2 + uz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfacb * (uz - c1o2) + (mfaca + c1o3 * oMdrho) * (     uz2 - uz) * c1o2;
               m1 = -mfacc        - 2. * mfacb *  uz         +  mfaca                  * (1. - uz2)              - c1o3 * oMdrho * uz2;
               m2 =  mfacc * c1o2 +      mfacb * (uz + c1o2) + (mfaca + c1o3 * oMdrho) * (     uz2 + uz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbac * c1o2 +      mfbab * (uz - c1o2) + mfbaa * (     uz2 - uz) * c1o2;
               m1 = -mfbac        - 2. * mfbab *  uz         + mfbaa * (1. - uz2);
               m2 =  mfbac * c1o2 +      mfbab * (uz + c1o2) + mfbaa * (     uz2 + uz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbbc * c1o2 +      mfbbb * (uz - c1o2) + mfbba * (     uz2 - uz) * c1o2;
               m1 = -mfbbc        - 2. * mfbbb *  uz         + mfbba * (1. - uz2);
               m2 =  mfbbc * c1o2 +      mfbbb * (uz + c1o2) + mfbba * (     uz2 + uz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbcb * (uz - c1o2) + mfbca * (     uz2 - uz) * c1o2;
               m1 = -mfbcc        - 2. * mfbcb *  uz         + mfbca * (1. - uz2);
               m2 =  mfbcc * c1o2 +      mfbcb * (uz + c1o2) + mfbca * (     uz2 + uz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfcab * (uz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     uz2 - uz) * c1o2;
               m1 = -mfcac        - 2. * mfcab *  uz         +  mfcaa                  * (1. - uz2)              - c1o3 * oMdrho * uz2;
               m2 =  mfcac * c1o2 +      mfcab * (uz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     uz2 + uz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfcbb * (uz - c1o2) + mfcba * (     uz2 - uz) * c1o2;
               m1 = -mfcbc        - 2. * mfcbb *  uz         + mfcba * (1. - uz2);
               m2 =  mfcbc * c1o2 +      mfcbb * (uz + c1o2) + mfcba * (     uz2 + uz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfccb * (uz - c1o2) + (mfcca + c1o9 * oMdrho) * (     uz2 - uz) * c1o2;
               m1 = -mfccc        - 2. * mfccb *  uz         +  mfcca                  * (1. - uz2)              - c1o9 * oMdrho * uz2;
               m2 =  mfccc * c1o2 +      mfccb * (uz + c1o2) + (mfcca + c1o9 * oMdrho) * (     uz2 + uz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 =  mfaca * c1o2 +      mfaba * (uy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfaca        - 2. * mfaba *  uy         +  mfaaa                  * (1. - uy2)              - c1o6 * oMdrho * uy2;
               m2 =  mfaca * c1o2 +      mfaba * (uy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     uy2 + uy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacb * c1o2 +      mfabb * (uy - c1o2) + (mfaab + c2o3 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfacb        - 2. * mfabb *  uy         +  mfaab                  * (1. - uy2)              - c2o3 * oMdrho * uy2;
               m2 =  mfacb * c1o2 +      mfabb * (uy + c1o2) + (mfaab + c2o3 * oMdrho) * (     uy2 + uy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfabc * (uy - c1o2) + (mfaac + c1o6 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfacc        - 2. * mfabc *  uy         +  mfaac                  * (1. - uy2)              - c1o6 * oMdrho * uy2;
               m2 =  mfacc * c1o2 +      mfabc * (uy + c1o2) + (mfaac + c1o6 * oMdrho) * (     uy2 + uy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbca * c1o2 +      mfbba * (uy - c1o2) + mfbaa * (     uy2 - uy) * c1o2;
               m1 = -mfbca        - 2. * mfbba *  uy         + mfbaa * (1. - uy2);
               m2 =  mfbca * c1o2 +      mfbba * (uy + c1o2) + mfbaa * (     uy2 + uy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcb * c1o2 +      mfbbb * (uy - c1o2) + mfbab * (     uy2 - uy) * c1o2;
               m1 = -mfbcb        - 2. * mfbbb *  uy         + mfbab * (1. - uy2);
               m2 =  mfbcb * c1o2 +      mfbbb * (uy + c1o2) + mfbab * (     uy2 + uy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbbc * (uy - c1o2) + mfbac * (     uy2 - uy) * c1o2;
               m1 = -mfbcc        - 2. * mfbbc *  uy         + mfbac * (1. - uy2);
               m2 =  mfbcc * c1o2 +      mfbbc * (uy + c1o2) + mfbac * (     uy2 + uy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfcba * (uy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfcca        - 2. * mfcba *  uy         +  mfcaa                   * (1. - uy2)              - c1o18 * oMdrho * uy2;
               m2 =  mfcca * c1o2 +      mfcba * (uy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     uy2 + uy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfcbb * (uy - c1o2) + (mfcab + c2o9 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfccb        - 2. * mfcbb *  uy         +  mfcab                  * (1. - uy2)              - c2o9 * oMdrho * uy2;
               m2 =  mfccb * c1o2 +      mfcbb * (uy + c1o2) + (mfcab + c2o9 * oMdrho) * (     uy2 + uy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfcbc * (uy - c1o2) + (mfcac + c1o18 * oMdrho) * (     uy2 - uy) * c1o2;
               m1 = -mfccc        - 2. * mfcbc *  uy         +  mfcac                   * (1. - uy2)              - c1o18 * oMdrho * uy2;
               m2 =  mfccc * c1o2 +      mfcbc * (uy + c1o2) + (mfcac + c1o18 * oMdrho) * (     uy2 + uy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 =  mfcaa * c1o2 +      mfbaa * (ux - c1o2) + (mfaaa + c1o36 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcaa        - 2. * mfbaa *  ux         +  mfaaa                   * (1. - ux2)              - c1o36 * oMdrho * ux2;
               m2 =  mfcaa * c1o2 +      mfbaa * (ux + c1o2) + (mfaaa + c1o36 * oMdrho) * (     ux2 + ux) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcba * c1o2 +      mfbba * (ux - c1o2) + (mfaba + c1o9 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcba        - 2. * mfbba *  ux         +  mfaba                  * (1. - ux2)              - c1o9 * oMdrho * ux2;
               m2 =  mfcba * c1o2 +      mfbba * (ux + c1o2) + (mfaba + c1o9 * oMdrho) * (     ux2 + ux) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfbca * (ux - c1o2) + (mfaca + c1o36 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcca        - 2. * mfbca *  ux         +  mfaca                   * (1. - ux2)              - c1o36 * oMdrho * ux2;
               m2 =  mfcca * c1o2 +      mfbca * (ux + c1o2) + (mfaca + c1o36 * oMdrho) * (     ux2 + ux) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcab * c1o2 +      mfbab * (ux - c1o2) + (mfaab + c1o9 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcab        - 2. * mfbab *  ux         +  mfaab                  * (1. - ux2)              - c1o9 * oMdrho * ux2;
               m2 =  mfcab * c1o2 +      mfbab * (ux + c1o2) + (mfaab + c1o9 * oMdrho) * (     ux2 + ux) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfcbb * c1o2 +      mfbbb * (ux - c1o2) + (mfabb + c4o9 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcbb        - 2. * mfbbb *  ux         +  mfabb                  * (1. - ux2)              - c4o9 * oMdrho * ux2;
               m2 =  mfcbb * c1o2 +      mfbbb * (ux + c1o2) + (mfabb + c4o9 * oMdrho) * (     ux2 + ux) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfbcb * (ux - c1o2) + (mfacb + c1o9 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfccb        - 2. * mfbcb *  ux         +  mfacb                  * (1. - ux2)              - c1o9 * oMdrho * ux2;
               m2 =  mfccb * c1o2 +      mfbcb * (ux + c1o2) + (mfacb + c1o9 * oMdrho) * (     ux2 + ux) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfbac * (ux - c1o2) + (mfaac + c1o36 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcac        - 2. * mfbac *  ux         +  mfaac                   * (1. - ux2)              - c1o36 * oMdrho * ux2;
               m2 =  mfcac * c1o2 +      mfbac * (ux + c1o2) + (mfaac + c1o36 * oMdrho) * (     ux2 + ux) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfbbc * (ux - c1o2) + (mfabc + c1o9 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfcbc        - 2. * mfbbc *  ux         +  mfabc                  * (1. - ux2)              - c1o9 * oMdrho * ux2;
               m2 =  mfcbc * c1o2 +      mfbbc * (ux + c1o2) + (mfabc + c1o9 * oMdrho) * (     ux2 + ux) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfbcc * (ux - c1o2) + (mfacc + c1o36 * oMdrho) * (     ux2 - ux) * c1o2;
               m1 = -mfccc        - 2. * mfbcc *  ux         +  mfacc                   * (1. - ux2)              - c1o36 * oMdrho * ux2;
               m2 =  mfccc * c1o2 +      mfbcc * (ux + c1o2) + (mfacc + c1o36 * oMdrho) * (     ux2 + ux) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;



///////////////////////////////////////////////////////////////////////////

/*
			    m000 = f1[ZERO] + f1[E] + f1[W] + f1[N] + f1[S] + f1[T] + f1[B] + f1[NE] + f1[SW] + f1[SE] + f1[NW] + f1[TE] + f1[BW] + f1[BE] + f1[TW] + f1[TN] + f1[BS] + f1[BN] + 
					   f1[TS] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW] + f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW];
				
				m100 = f1[BE] + f1[BNE] - f1[BNW] + f1[BSE] - f1[BSW] - f1[BW] + f1[E] + f1[NE] - f1[NW] + f1[SE] - f1[SW] + f1[TE] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW] - f1[TW] - f1[W];
				m010 = f1[BN] + f1[BNE] + f1[BNW] - f1[BS] - f1[BSE] - f1[BSW] + f1[N] + f1[NE] + f1[NW] - f1[S] - f1[SE] - f1[SW] + f1[TN] + f1[TNE] + f1[TNW] - f1[TS] - f1[TSE] - f1[TSW];
				m001 = f1[T] + f1[TE] + f1[TN] + f1[TNE] + f1[TNW] + f1[TS] + f1[TSE] + f1[TSW] + f1[TW] - f1[B] - f1[BE] - f1[BN] - f1[BNE] - f1[BNW] - f1[BS] - f1[BSE] - f1[BSW] - f1[BW];

				m110 =  f1[BNE] - f1[BNW] - f1[BSE] + f1[BSW] + f1[NE] - f1[NW] - f1[SE] + f1[SW] + f1[TNE] - f1[TNW] - f1[TSE] + f1[TSW];
				m101 = -f1[BE] - f1[BNE] + f1[BNW] - f1[BSE] + f1[BSW] + f1[BW] + f1[TE] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW] - f1[TW];
				m011 = -f1[BN] - f1[BNE] - f1[BNW] + f1[BS] + f1[BSE] + f1[BSW] + f1[TN] + f1[TNE] + f1[TNW] - f1[TS] - f1[TSE] - f1[TSW];
				m200 =  f1[BE] + f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW] + f1[BW] + f1[E] + f1[NE] + f1[NW] + f1[SE] + f1[SW] + f1[TE] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW] + f1[TW] + f1[W];
				m020 =  f1[BN] + f1[BNE] + f1[BNW] + f1[BS] + f1[BSE] + f1[BSW] + f1[N] + f1[NE] + f1[NW] + f1[S] + f1[SE] + f1[SW] + f1[TN] + f1[TNE] + f1[TNW] + f1[TS] + f1[TSE] + f1[TSW];
				m002 =  f1[B] + f1[BE] + f1[BN] + f1[BNE] + f1[BNW] + f1[BS] + f1[BSE] + f1[BSW] + f1[BW] + f1[T] + f1[TE] + f1[TN] + f1[TNE] + f1[TNW] + f1[TS] + f1[TSE] + f1[TSW] + f1[TW];
				m120 =  f1[BNE] - f1[BNW] + f1[BSE] - f1[BSW] + f1[NE] - f1[NW] + f1[SE] - f1[SW] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW];
				m102 =  f1[BE] + f1[BNE] - f1[BNW] + f1[BSE] - f1[BSW] - f1[BW] + f1[TE] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW] - f1[TW];
				m210 =  f1[BNE] + f1[BNW] - f1[BSE] - f1[BSW] + f1[NE] + f1[NW] - f1[SE] - f1[SW] + f1[TNE] + f1[TNW] - f1[TSE] - f1[TSW];
				m012 =  f1[BN] + f1[BNE] + f1[BNW] - f1[BS] - f1[BSE] - f1[BSW] + f1[TN] + f1[TNE] + f1[TNW] - f1[TS] - f1[TSE] - f1[TSW];
				m201 = -f1[BE] - f1[BNE] - f1[BNW] - f1[BSE] - f1[BSW] - f1[BW] + f1[TE] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW] + f1[TW];
				m021 = -f1[BN] - f1[BNE] - f1[BNW] - f1[BS] - f1[BSE] - f1[BSW] + f1[TN] + f1[TNE] + f1[TNW] + f1[TS] + f1[TSE] + f1[TSW];
				m111 = -f1[BNE] + f1[BNW] + f1[BSE] - f1[BSW] + f1[TNE] - f1[TNW] - f1[TSE] + f1[TSW];
				m220 =  f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW] + f1[NE] + f1[NW] + f1[SE] + f1[SW] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW];
				m202 =  f1[BE] + f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW] + f1[BW] + f1[TE] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW] + f1[TW];
				m022 =  f1[BN] + f1[BNE] + f1[BNW] + f1[BS] + f1[BSE] + f1[BSW] + f1[TN] + f1[TNE] + f1[TNW] + f1[TS] + f1[TSE] + f1[TSW];
				m211 = -f1[BNE] - f1[BNW] + f1[BSE] + f1[BSW] + f1[TNE] + f1[TNW] - f1[TSE] - f1[TSW];
				m121 = -f1[BNE] + f1[BNW] - f1[BSE] + f1[BSW] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW];
				m112 =  f1[BNE] - f1[BNW] - f1[BSE] + f1[BSW] + f1[TNE] - f1[TNW] - f1[TSE] + f1[TSW];
				m221 = -f1[BNE] - f1[BNW] - f1[BSE] - f1[BSW] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW];
				m212 =  f1[BNE] + f1[BNW] - f1[BSE] - f1[BSW] + f1[TNE] + f1[TNW] - f1[TSE] - f1[TSW];
				m122 =  f1[BNE] - f1[BNW] + f1[BSE] - f1[BSW] + f1[TNE] - f1[TNW] + f1[TSE] - f1[TSW];
				m222 =  f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW];

				k200 = m200 - m000*ux2;
				k020 = m020 - m000*uy2;
				k002 = m002 - m000*uz2; 
				k110 = m110 - m000*ux_uy;          
				k101 = m101 - m000*ux_uz;
				k011 = m011 - m000*uy_uz;

				k021 = m021 - (m020*uz + 2.0*uy*k011);          
				k012 = m012 - (uy*m002 + 2.0*uz*k011);
				k102 = m102 - (ux*m002 + 2.0*uz*k101);
				k201 = m201 - (m200*uz + 2.0*ux*k101 );
				k210 = m210 - (m200*uy + 2.0*ux*k110 );
				k120 = m120 - (ux*m020 + 2.0*uy*k110);
				k111 = m111 - (ux_uy_uz + ux*k011 + uy*k101 + uz*k110);

				k121 = m121 - (uz*m120+ 2.0*ux_uy*k011+ux*k021+uy2*k101+2.0*uy*k111); 
				k112 = m112 - (uy*m102+ 2.0*ux_uz*k011+ux*k012+uz2*k110+2.0*uz*k111);
				k211 = m211 - (uz*m210+ 2.0*ux_uy*k101+uy*k201+ux2*k011+2.0*ux*k111); 


				k220 = m220 - ( ux2*m020+ 4.0*ux_uy*k110 + 2.0*ux*k120 + uy2*k200 + 2.0*uy*k210); 
				k022 = m022 - ( uy2*m002+ 4.0*uy_uz*k011 + 2.0*uy*k012 + uz2*k020 + 2.0*uz*k021);
				k202 = m202 - ( ux2*m002+ 4.0*ux_uz*k101 + 2.0*ux*k102 + uz2*k200 + 2.0*uz*k201); 

				k221 = m221 - (uz*m220+
					2.0*ux2*uy*k011  + ux2*k021 + 2.0*ux*uy2*k101  +4.0*ux_uy*k111 + 
					2.0*ux*k121 + uy2*k201 + 2.0*uy*k211);
				k122 = m122 - (ux*m022 
					+ 2.0*uy2*uz*k101 + uy2*k102 + 2.0*uy*uz2*k110 + 4.0*uy_uz*k111 +
					2.0*uy*k112 + uz2*k120 + 2.0*uz*k121);
				k212 = m212 - (uy*m202
					+ 2.0*ux2*uz*k011 +ux2*k012  + 2.0*ux*uz2*k110 + 4.0*ux_uz*k111 +
					2.0*ux*k112 + uz2*k210 + 2.0*uz*k211);

				k222 = m222 - (ux2*m022 
					+ 4.0* ux*uy2*uz*k101   + 
					2.0* ux*uy2*k102 + 4.0* ux_uy*uz2*k110+ 
					8.0* ux_uy_uz*k111 + 4.0* ux_uy*k112 + 
					2.0* ux*uz2*k120 + 4.0* ux_uz*k121 + 2.0*ux*k122 + 
					uy2*uz2*k200 + 2.0* uy2*uz*k201  + 
					uy2*k202     + 2.0* uy*uz2*k210  + 4.0*uy_uz*k211 + 
					2.0* uy*k212 + uz2*k220 +2.0* uz*k221);

				//////////////// Central moments to Cumulants \\\\\\\\\\\\\\\\\\\\\\\\\

				c200 = k200;
				c020 = k020;
				c002 = k002;
				c110 = k110;
				c101 = k101;
				c011 = k011;

				c021 = k021;
				c012 = k012;
				c102 = k102;
				c201 = k201;
				c210 = k210;
				c120 = k120;
				c111 = k111;

				c121 = k121 - (k020*k101 + 2*k011*k110)/m000;
				c211 = k211 - (k200*k011 + 2*k110*k101)/m000;
				c112 = k112 - (k002*k110 + 2*k011*k101)/m000;

				c220 = k220 - (k200*k020 + 2*k110*k110)/m000;
				c202 = k202 - (k200*k002 + 2*k101*k101)/m000;
				c022 = k022 - (k020*k002 + 2*k011*k011)/m000;

				c122 = k122 - (k002*k120 + k020*k102 + 4*k011*k111 + 2*(k101*k021 + k110*k012))/m000;
				c212 = k212 - (k200*k012 + k002*k210 + 4*k101*k111 + 2*(k110*k102 + k011*k201))/m000;
				c221 = k221 - (k200*k021 + k020*k201 + 4*k110*k111 + 2*(k101*k120 + k011*k210))/m000;

				c222 = k222 - (4*k111*k111 + k200*k022 + k020*k202 + k002*k220 + 4*(k011*k211 + k101*k121 + k110*k112)
					+ 2*(k120*k102 + k210*k012 + k201*k021))/m000
					+ (16*k110*k101*k011 + 4*(k101*k101*k020 + k011*k011*k200 + k110*k110*k002) + 2*k200*k020*k002)/(m000*m000);

				c200_pl_c020_pl_c002 = c200 + c020 + c002;
				c200_mi_c020 = c200 - c020;
				c200_mi_c002 = c200 - c002;

				c210_pl_c012  =  c210+c012;
				c210_mi_c012  =  c210-c012;
				c201_pl_c021  =  c201+c021;
				c201_mi_c021  =  c201-c021;
				c120_pl_c102  =  c120+c102;
				c120_mi_c102  =  c120-c102;

				c220_pl_c202_pl_c022  = c220 + c202 + c022;
				c220_mi2_c202_pl_c022 = c220 - 2.0*c202 + c022;
				c220_pl_c202_mi2_c022 = c220 + c202 - 2.0*c022;

				/////////////////////////// Relaxation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

				c110 = c110*(1-w1);
				c101 = c101*(1-w1);
				c011 = c011*(1-w1);

				c200_mi_c020 = c200_mi_c020*(1-w1); 
				c200_mi_c002 = c200_mi_c002*(1-w1); 

				c200_pl_c020_pl_c002 = c200_pl_c020_pl_c002*(1-w2) + m000*w2;

				c120_pl_c102 = c120_pl_c102*(1-w3);
				c210_pl_c012 = c210_pl_c012*(1-w3);
				c201_pl_c021 = c201_pl_c021*(1-w3); 


				c120_mi_c102 = c120_mi_c102*(1-w4); 
				c210_mi_c012 = c210_mi_c012*(1-w4);
				c201_mi_c021 = c201_mi_c021*(1-w4);

				c111 = c111*(1-w5);

				c220_mi2_c202_pl_c022 = c220_mi2_c202_pl_c022*(1-w6);
				c220_pl_c202_mi2_c022 = c220_pl_c202_mi2_c022*(1-w6);

				c220_pl_c202_pl_c022 =  c220_pl_c202_pl_c022*(1-w7);

				c211 = c211*(1-w8);
				c121 = c121*(1-w8);
				c112 = c112*(1-w8);


				c221 = c221*(1-w9);
				c212 = c212*(1-w9);
				c122 = c122*(1-w9);

				c222 = c222*(1-w10);


				c200 = c1o3 *c200_mi_c020 + c1o3 *c200_mi_c002 +  c1o3* c200_pl_c020_pl_c002; 
				c020 = -2*c1o3* c200_mi_c020+ c1o3* c200_mi_c002 +  c1o3* c200_pl_c020_pl_c002; 
				c002 = c1o3   * c200_mi_c020 -2*c1o3 *c200_mi_c002 +  c1o3* c200_pl_c020_pl_c002; 

				c210 = (c210_mi_c012 + c210_pl_c012)*0.5; 
				c012 = 0.5*(-c210_mi_c012 + c210_pl_c012); 
				c120 =(c120_mi_c102 + c120_pl_c102)*0.5; 
				c102 = 0.5*(-c120_mi_c102 + c120_pl_c102); 
				c201 = (c201_mi_c021 + c201_pl_c021)*0.5; 
				c021 = 0.5*(-c201_mi_c021 + c201_pl_c021);

				c220 =  c1o3* c220_mi2_c202_pl_c022 + c1o3* c220_pl_c202_mi2_c022 +  c1o3*c220_pl_c202_pl_c022; 
				c202 = -c1o3* c220_mi2_c202_pl_c022 + c1o3* c220_pl_c202_pl_c022; 
				c022 = -c1o3* c220_pl_c202_mi2_c022 + c1o3* c220_pl_c202_pl_c022;


				////////////////////// Cumulants to Central moments   \\\\\\\\\\\\\\\\\\\\\\\\\

				k200 = c200;
				k020 = c020;
				k002 = c002;
				k110 = c110;
				k101 = c101;
				k011 = c011;

				k021 = c021;
				k012 = c012;
				k102 = c102;
				k201 = c201;
				k210 = c210;
				k120 = c120;
				k111 = c111;

				k121 = c121 + (k020*k101 + 2*k011*k110)/m000;
				k211 = c211 + (k200*k011 + 2*k110*k101)/m000;
				k112 = c112 + (k002*k110 + 2*k011*k101)/m000;

				k220 = c220 + (k200*k020 + 2*k110*k110)/m000;
				k202 = c202 + (k200*k002 + 2*k101*k101)/m000;
				k022 = c022 + (k020*k002 + 2*k011*k011)/m000;

				k122 = c122 + (k002*k120 + k020*k102 + 4*k011*k111 + 2*(k101*k021 + k110*k012))/m000;
				k212 = c212 + (k200*k012 + k002*k210 + 4*k101*k111 + 2*(k110*k102 + k011*k201))/m000;
				k221 = c221 + (k200*k021 + k020*k201 + 4*k110*k111 + 2*(k101*k120 + k011*k210))/m000;

				k222 = c222 + (4*k111*k111 + k200*k022 + k020*k202 + k002*k220 + 4*(k011*k211 + k101*k121 + k110*k112)
					+ 2*(k120*k102 + k210*k012 + k201*k021))/m000
					- (16*k110*k101*k011 + 4*(k101*k101*k020 + k011*k011*k200 + k110*k110*k002) + 2*k200*k020*k002)/(m000*m000);

				///////////////////////////////////////////////////////////////////////////////


				m200 = k200 + m000*ux2;
				m020 = k020 + m000*uy2;
				m002 = k002 + m000*uz2; 
				m110 = k110 + m000*ux_uy;          
				m101 = k101 + m000*ux_uz;
				m011 = k011 + m000*uy_uz;

				m021 = m020*uz + 2.0*uy*k011  + k021;          
				m012 = uy*m002 + 2.0*uz*k011 + k012;
				m102 = ux*m002 + 2.0*uz*k101 + k102;
				m112 = uy*m102 + 2.0*ux_uz*k011+ux*k012+uz2*k110+2.0*uz*k111+k112;

				m201 = m200*uz + 2.0*ux*k101  + k201;
				m210 = m200*uy + 2.0*ux*k110  + k210;
				m211 = uz*m210 + 2.0*ux_uy*k101+uy*k201+ux2*k011+2.0*ux*k111+k211; 

				m120 = ux*m020 + 2.0*uy*k110 + k120;
				m121 = uz*m120 + 2.0*ux_uy*k011+ux*k021+uy2*k101+2.0*uy*k111+ k121; 

				m111 = ux_uy_uz + ux*k011 + uy*k101 + uz*k110+k111;

				m220 = ux2*m020 + 4.0*ux_uy*k110 + 2.0*ux*k120  + uy2*k200 + 2.0*uy*k210 +  k220; 

				m221 = uz*m220 + 2.0*ux2*uy*k011  + ux2*k021 + 2.0*ux*uy2*k101  +4.0*ux_uy*k111+ 
					2.0*ux*k121  +    uy2*k201  + 2.0*uy*k211  + k221;

				m022 = uy2*m002 + 4.0*uy_uz*k011  +2.0*uy*k012+ uz2*k020+ 2.0*uz*k021 +  k022;

				m122 = ux*m022 + 2.0*uy2*uz*k101 + uy2*k102 + 
					2.0*uy*uz2*k110 + 4.0*uy_uz*k111 + 2.0*uy*k112 +uz2*k120 + 2.0*uz*k121 + k122;

				m202 = ux2*m002 + 4.0*ux_uz*k101+ 2.0*ux*k102 + uz2*k200+ 2.0*uz*k201 + k202; 

				m212 = uy*m202 + 2.0*ux2*uz*k011 +ux2*k012  + 2.0*ux*uz2*k110 + 4.0*ux_uz*k111 
					+ 2.0*ux*k112  + uz2*k210 +  2.0*uz*k211 + k212;

				m222 = ux2*m022 
					+ 4.0* ux*uy2*uz*k101   + 
					2.0* ux*uy2*k102 + 4.0* ux_uy*uz2*k110+ 
					8.0* ux_uy_uz*k111 + 4.0* ux_uy*k112 + 
					2.0* ux*uz2*k120 + 4.0* ux_uz*k121 + 2.0*ux*k122 + 
					uy2*uz2*k200 + 2.0* uy2*uz*k201  + 
					uy2*k202     + 2.0* uy*uz2*k210  + 4.0*uy_uz*k211 + 
					2.0* uy*k212 + uz2*k220 +2.0* uz*k221 + k222;

				f1[ZERO] = (-m200 + m220 - m222 + m202 - m020 + m022 - m002 + m000);
				f1[E] = 0.5* (m200 -  m220 + m222 - m202 - m120 + m122 - m102 +m100);
				f1[W] = 0.5* (m200 - m220 + m222 - m202 + m120 - m122 + m102 -m100);
				f1[N] = 0.5* (-m210 - m220 + m222 + m212 + m020 - m022 - m012 +m010);
				f1[S] = 0.5* (m210 -  m220 + m222 - m212 + m020 - m022 + m012 -m010);
				f1[T] = 0.5* (m221 +  m222 - m201 - m202 - m021 - m022 + m002 +m001);
				f1[B] = 0.5* (-m221 + m222 + m201  - m202 + m021 - m022 + m002-m001);

				f1[NE] = 0.25*( m210  + m220- m222 - m212 + m110+ m120- m122 -m112); 
				f1[SW] = 0.25*(-m210 + m220- m222 + m212 + m110- m120+ m122 -m112); 
				f1[SE] = 0.25*(-m210 + m220- m222 + m212 - m110+ m120- m122 +m112); 
				f1[NW] = 0.25*( m210  + m220- m222 - m212 - m110- m120+ m122 + m112); 
				f1[TE] = 0.25*(-m221 - m222 + m201 + m202 - m121 - m122 + m101 + m102); 
				f1[BW] = 0.25*( m221  -m222 - m201 + m202 - m121 + m122 + m101 - m102);
				f1[BE] = 0.25*(m221 - m222 - m201 + m202 + m121 - m122 - m101 +m102);
				f1[TW] = 0.25*(-m221 - m222 + m201 + m202 + m121 + m122 - m101 -m102); 
				f1[TN] = 0.25*(-m221 - m222 - m211 - m212 + m021 + m022 + m011+m012);
				f1[BS] = 0.25*( m221 - m222 - m211 + m212 - m021 + m022 + m011 - m012);
				f1[BN] = 0.25*( m221 - m222 + m211 - m212 - m021 + m022 - m011 + m012);
				f1[TS] = 0.25*(-m221 - m222 + m211 + m212 + m021 + m022 - m011 -m012); 

				f1[TNE]=0.125*( m221 + m222 + m211 + m212 + m121 + m122 + m111 + m112); 
				f1[BNE]=0.125*(-m221 + m222 -m211 + m212 -m121 + m122 -m111 + m112);
				f1[TSE]=0.125*( m221 + m222 - m211 - m212 + m121 + m122 - m111 - m112); 
				f1[BSE]=0.125*(-m221 + m222 +m211 - m212 -m121 + m122 +m111 - m112); 
				f1[TNW]=0.125*( m221 + m222 + m211 + m212 - m121 - m122 - m111 - m112); 
				f1[BNW]=0.125*(-m221 + m222 -m211 + m212 +m121 - m122 +m111 - m112); 
				f1[TSW]=0.125*( m221 + m222 - m211 - m212 - m121 - m122 + m111 + m112); 
				f1[BSW]=0.125*(-m221 + m222+m211 - m212+m121 - m122-m111 + m112);
			   
*/
			   //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
			   LBMReal rho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
			   +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
			   +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;
				
				/*LBMReal rho_post = f1[ZERO] + f1[E] + f1[W] + f1[N] + f1[S] + f1[T] + f1[B] + f1[NE] + f1[SW] + f1[SE] + f1[NW] + f1[TE] + f1[BW] + 
					f1[BE] + f1[TW] + f1[TN] + f1[BS] + f1[BN] + f1[TS] + f1[TNE] + f1[TNW] + f1[TSE] + f1[TSW] + f1[BNE] + f1[BNW] + f1[BSE] + f1[BSW]; */
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = rho1 - rho_post;
#ifdef SINGLEPRECISION
               if(dif > 10.0E-7 || dif < -10.0E-7)
#else
               if(dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
				   UB_THROW(UbException(UB_EXARGS,"rho="+UbSystem::toString(rho)+", rho_post="+UbSystem::toString(rho_post)
                     +" dif="+UbSystem::toString(dif)
                     +" rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)));
                  //UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
                  //exit(EXIT_FAILURE);
               }
#endif
               			   
			   
			   mfcbb = rho*c1o3*(mfcbb ) + 0.5*forcingTerm[E  ] ;
			   mfbcb = rho*c1o3*(mfbcb ) + 0.5*forcingTerm[N  ] ;
			   mfbbc = rho*c1o3*(mfbbc ) + 0.5*forcingTerm[T  ] ;
			   mfccb = rho*c1o3*(mfccb ) + 0.5*forcingTerm[NE ] ;
			   mfacb = rho*c1o3*(mfacb ) + 0.5*forcingTerm[NW ] ;
			   mfcbc = rho*c1o3*(mfcbc ) + 0.5*forcingTerm[TE ] ;
			   mfabc = rho*c1o3*(mfabc ) + 0.5*forcingTerm[TW ] ;
			   mfbcc = rho*c1o3*(mfbcc ) + 0.5*forcingTerm[TN ] ;
			   mfbac = rho*c1o3*(mfbac ) + 0.5*forcingTerm[TS ] ;
			   mfccc = rho*c1o3*(mfccc ) + 0.5*forcingTerm[TNE] ;
			   mfacc = rho*c1o3*(mfacc ) + 0.5*forcingTerm[TNW] ;
			   mfcac = rho*c1o3*(mfcac ) + 0.5*forcingTerm[TSE] ;
			   mfaac = rho*c1o3*(mfaac ) + 0.5*forcingTerm[TSW] ;
			   mfabb = rho*c1o3*(mfabb ) + 0.5*forcingTerm[W  ] ;
			   mfbab = rho*c1o3*(mfbab ) + 0.5*forcingTerm[S  ] ;
			   mfbba = rho*c1o3*(mfbba ) + 0.5*forcingTerm[B  ] ;
			   mfaab = rho*c1o3*(mfaab ) + 0.5*forcingTerm[SW ] ;
			   mfcab = rho*c1o3*(mfcab ) + 0.5*forcingTerm[SE ] ;
			   mfaba = rho*c1o3*(mfaba ) + 0.5*forcingTerm[BW ] ;
			   mfcba = rho*c1o3*(mfcba ) + 0.5*forcingTerm[BE ] ;
			   mfbaa = rho*c1o3*(mfbaa ) + 0.5*forcingTerm[BS ] ;
			   mfbca = rho*c1o3*(mfbca ) + 0.5*forcingTerm[BN ] ;
			   mfaaa = rho*c1o3*(mfaaa ) + 0.5*forcingTerm[BSW] ;
			   mfcaa = rho*c1o3*(mfcaa ) + 0.5*forcingTerm[BSE] ;
			   mfaca = rho*c1o3*(mfaca ) + 0.5*forcingTerm[BNW] ;
			   mfcca = rho*c1o3*(mfcca ) + 0.5*forcingTerm[BNE] ;
			   mfbbb = rho*c1o3*(mfbbb ) + 0.5*forcingTerm[ZERO];
			   
			   //////////////////////////////////////////////////////////////////////////
			   //write distribution for F
			   //////////////////////////////////////////////////////////////////////////

			   (*this->localDistributionsF)(D3Q27System::ET_E,   x1,  x2,  x3) = mfabb;
			   (*this->localDistributionsF)(D3Q27System::ET_N,   x1,  x2,  x3) = mfbab;
			   (*this->localDistributionsF)(D3Q27System::ET_T,   x1,  x2,  x3) = mfbba;
			   (*this->localDistributionsF)(D3Q27System::ET_NE,  x1,  x2,  x3) = mfaab;
			   (*this->localDistributionsF)(D3Q27System::ET_NW,  x1p, x2,  x3) = mfcab;
			   (*this->localDistributionsF)(D3Q27System::ET_TE,  x1,  x2,  x3) = mfaba;
			   (*this->localDistributionsF)(D3Q27System::ET_TW,  x1p, x2,  x3) = mfcba;
			   (*this->localDistributionsF)(D3Q27System::ET_TN,  x1,  x2,  x3) = mfbaa;
			   (*this->localDistributionsF)(D3Q27System::ET_TS,  x1,  x2p, x3) = mfbca;
			   (*this->localDistributionsF)(D3Q27System::ET_TNE, x1,  x2,  x3) = mfaaa;
			   (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2,  x3) = mfcaa;
			   (*this->localDistributionsF)(D3Q27System::ET_TSE, x1,  x2p, x3) = mfaca;
			   (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

			   (*this->nonLocalDistributionsF)(D3Q27System::ET_W,   x1p, x2,  x3 ) = mfcbb;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_S,   x1,  x2p, x3 ) = mfbcb;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_B,   x1,  x2,  x3p) = mfbbc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = mfccb;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = mfacb;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BW,  x1p, x2,  x3p) = mfcbc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BE,  x1,  x2,  x3p) = mfabc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BS,  x1,  x2p, x3p) = mfbcc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BN,  x1,  x2,  x3p) = mfbac;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1,  x2p, x3p) = mfacc;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2,  x3p) = mfcac;
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1,  x2,  x3p) = mfaac;

			   (*this->zeroDistributionsF)(x1,x2,x3) = mfbbb;


			   
/////////////////////  P H A S E - F I E L D   S O L V E R /////////////////////////////////////////			   
/*			   
			   
			   mfcbb = (*this->localDistributionsH)(D3Q27System::ET_E, x1,x2,x3);
			   mfbcb = (*this->localDistributionsH)(D3Q27System::ET_N,x1,x2,x3); 
			   mfbbc = (*this->localDistributionsH)(D3Q27System::ET_T,x1,x2,x3);
			   mfccb = (*this->localDistributionsH)(D3Q27System::ET_NE,x1,x2,x3);
			   mfacb = (*this->localDistributionsH)(D3Q27System::ET_NW,x1p,x2,x3);
			   mfcbc = (*this->localDistributionsH)(D3Q27System::ET_TE,x1,x2,x3);
			   mfabc = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p,x2,x3);
			   mfbcc = (*this->localDistributionsH)(D3Q27System::ET_TN,x1,x2,x3);
			   mfbac = (*this->localDistributionsH)(D3Q27System::ET_TS,x1,x2p,x3);
			   mfccc = (*this->localDistributionsH)(D3Q27System::ET_TNE,x1,x2,x3);
			   mfacc = (*this->localDistributionsH)(D3Q27System::ET_TNW,x1p,x2,x3);
			   mfcac = (*this->localDistributionsH)(D3Q27System::ET_TSE,x1,x2p,x3);
			   mfaac = (*this->localDistributionsH)(D3Q27System::ET_TSW,x1p,x2p,x3);
			   mfabb = (*this->nonLocalDistributionsH)(D3Q27System::ET_W,x1p,x2,x3  );
			   mfbab = (*this->nonLocalDistributionsH)(D3Q27System::ET_S,x1,x2p,x3  );
			   mfbba = (*this->nonLocalDistributionsH)(D3Q27System::ET_B,x1,x2,x3p  );
			   mfaab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,x1p,x2p,x3 );
			   mfcab = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,x1,x2p,x3 );
			   mfaba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,x1p,x2,x3p );
			   mfcba = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,x1,x2,x3p );
			   mfbaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,x1,x2p,x3p );
			   mfbca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,x1,x2,x3p );
			   mfaaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW,x1p,x2p,x3p);
			   mfcaa = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE,x1,x2p,x3p);
			   mfaca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW,x1p,x2,x3p);
			   mfcca = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE,x1,x2,x3p);
			   mfbbb = (*this->zeroDistributionsH)(x1,x2,x3);

			   
			   LBMReal hSource = ( (tauH - 0.5)*(1.0 - phi[ZERO])*(phi[ZERO])/denom ); // + phi[ZERO]*(dxux + dyuy + dzuz);


			   LBMReal drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
				   (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				   ((mfabb + mfcbb) + (mfbab + mfbcb)) + (mfbba + mfbbc)) + mfbbb;

			   
			   LBMReal collFactorPhi = 1.0 / tauH;
			   oMdrho = one; // comp special
			   ////////////////////////////////////////////////////////////////////////////////////

			   // 						LBMReal wadjust;
			   // 						LBMReal qudricLimitP = 0.01f;// * 0.0001f;
			   // 						LBMReal qudricLimitM = 0.01f;// * 0.0001f;
			   // 						LBMReal qudricLimitD = 0.01f;// * 0.001f;
			   //LBMReal s9 = minusomega;
			   //test
			   //s9 = 0.;
			   ////////////////////////////////////////////////////////////////////////////////////
			   //Hin
			   ////////////////////////////////////////////////////////////////////////////////////
			   // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   m2 = mfaaa + mfaac; // 2nd raw moment
			   m1 = mfaac - mfaaa; // 1st raw moment
			   m0 = m2 + mfaab;    // zeroth raw moment
			   mfaaa = m0;
			   m0 += c1o36 * oMdrho;
			   mfaab = m1 - m0 * uz;                   // this corresponds to a (central moment of order 1 = first raw moment  - velocity * zeroth moment) 
			   mfaac = m2 - two*	m1 * uz + uz2 * m0; // this corresponds to a (central moment of order 2 = second raw moment - 2* velocity * 1st raw moment + velocity^2 * zeroth moment)
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaba + mfabc;
			   m1 = mfabc - mfaba;
			   m0 = m2 + mfabb;
			   mfaba = m0;
			   m0 += c1o9 * oMdrho;
			   mfabb = m1 - m0 * uz;
			   mfabc = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaca + mfacc;
			   m1 = mfacc - mfaca;
			   m0 = m2 + mfacb;
			   mfaca = m0;
			   m0 += c1o36 * oMdrho;
			   mfacb = m1 - m0 * uz;
			   mfacc = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbaa + mfbac;
			   m1 = mfbac - mfbaa;
			   m0 = m2 + mfbab;
			   mfbaa = m0;
			   m0 += c1o9 * oMdrho;
			   mfbab = m1 - m0 * uz;
			   mfbac = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbba + mfbbc;
			   m1 = mfbbc - mfbba;
			   m0 = m2 + mfbbb;
			   mfbba = m0;
			   m0 += c4o9 * oMdrho;
			   mfbbb = m1 - m0 * uz;
			   mfbbc = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbca + mfbcc;
			   m1 = mfbcc - mfbca;
			   m0 = m2 + mfbcb;
			   mfbca = m0;
			   m0 += c1o9 * oMdrho;
			   mfbcb = m1 - m0 * uz;
			   mfbcc = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcaa + mfcac;
			   m1 = mfcac - mfcaa;
			   m0 = m2 + mfcab;
			   mfcaa = m0;
			   m0 += c1o36 * oMdrho;
			   mfcab = m1 - m0 * uz;
			   mfcac = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcba + mfcbc;
			   m1 = mfcbc - mfcba;
			   m0 = m2 + mfcbb;
			   mfcba = m0;
			   m0 += c1o9 * oMdrho;
			   mfcbb = m1 - m0 * uz;
			   mfcbc = m2 - two*	m1 * uz + uz2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcca + mfccc;
			   m1 = mfccc - mfcca;
			   m0 = m2 + mfccb;
			   mfcca = m0;
			   m0 += c1o36 * oMdrho;
			   mfccb = m1 - m0 * uz;
			   mfccc = m2 - two*	m1 * uz + uz2 * m0;
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
			   mfaba = m1 - m0 * uy;
			   mfaca = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaab + mfacb;
			   m1 = mfacb - mfaab;
			   m0 = m2 + mfabb;
			   mfaab = m0;
			   mfabb = m1 - m0 * uy;
			   mfacb = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaac + mfacc;
			   m1 = mfacc - mfaac;
			   m0 = m2 + mfabc;
			   mfaac = m0;
			   m0 += c1o18 * oMdrho;
			   mfabc = m1 - m0 * uy;
			   mfacc = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbaa + mfbca;
			   m1 = mfbca - mfbaa;
			   m0 = m2 + mfbba;
			   mfbaa = m0;
			   m0 += c2o3 * oMdrho;
			   mfbba = m1 - m0 * uy;
			   mfbca = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbab + mfbcb;
			   m1 = mfbcb - mfbab;
			   m0 = m2 + mfbbb;
			   mfbab = m0;
			   mfbbb = m1 - m0 * uy;
			   mfbcb = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfbac + mfbcc;
			   m1 = mfbcc - mfbac;
			   m0 = m2 + mfbbc;
			   mfbac = m0;
			   m0 += c2o9 * oMdrho;
			   mfbbc = m1 - m0 * uy;
			   mfbcc = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcaa + mfcca;
			   m1 = mfcca - mfcaa;
			   m0 = m2 + mfcba;
			   mfcaa = m0;
			   m0 += c1o6 * oMdrho;
			   mfcba = m1 - m0 * uy;
			   mfcca = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcab + mfccb;
			   m1 = mfccb - mfcab;
			   m0 = m2 + mfcbb;
			   mfcab = m0;
			   mfcbb = m1 - m0 * uy;
			   mfccb = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfcac + mfccc;
			   m1 = mfccc - mfcac;
			   m0 = m2 + mfcbc;
			   mfcac = m0;
			   m0 += c1o18 * oMdrho;
			   mfcbc = m1 - m0 * uy;
			   mfccc = m2 - two*	m1 * uy + uy2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   m2 = mfaaa + mfcaa;
			   m1 = mfcaa - mfaaa;
			   m0 = m2 + mfbaa;
			   mfaaa = m0;
			   m0 += one* oMdrho;
			   mfbaa = m1 - m0 * ux;
			   mfcaa = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaba + mfcba;
			   m1 = mfcba - mfaba;
			   m0 = m2 + mfbba;
			   mfaba = m0;
			   mfbba = m1 - m0 * ux;
			   mfcba = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaca + mfcca;
			   m1 = mfcca - mfaca;
			   m0 = m2 + mfbca;
			   mfaca = m0;
			   m0 += c1o3 * oMdrho;
			   mfbca = m1 - m0 * ux;
			   mfcca = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaab + mfcab;
			   m1 = mfcab - mfaab;
			   m0 = m2 + mfbab;
			   mfaab = m0;
			   mfbab = m1 - m0 * ux;
			   mfcab = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfabb + mfcbb;
			   m1 = mfcbb - mfabb;
			   m0 = m2 + mfbbb;
			   mfabb = m0;
			   mfbbb = m1 - m0 * ux;
			   mfcbb = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfacb + mfccb;
			   m1 = mfccb - mfacb;
			   m0 = m2 + mfbcb;
			   mfacb = m0;
			   mfbcb = m1 - m0 * ux;
			   mfccb = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfaac + mfcac;
			   m1 = mfcac - mfaac;
			   m0 = m2 + mfbac;
			   mfaac = m0;
			   m0 += c1o3 * oMdrho;
			   mfbac = m1 - m0 * ux;
			   mfcac = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfabc + mfcbc;
			   m1 = mfcbc - mfabc;
			   m0 = m2 + mfbbc;
			   mfabc = m0;
			   mfbbc = m1 - m0 * ux;
			   mfcbc = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m2 = mfacc + mfccc;
			   m1 = mfccc - mfacc;
			   m0 = m2 + mfbcc;
			   mfacc = m0;
			   m0 += c1o9 * oMdrho;
			   mfbcc = m1 - m0 * ux;
			   mfccc = m2 - two*	m1 * ux + ux2 * m0;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////


			   ////////////////////////////////////////////////////////////////////////////////////
			   // Collision
			   ////////////////////////////////////////////////////////////////////////////////////
			   LBMReal m000 = mfaaa;
			   LBMReal m100 = mfbaa;
			   LBMReal m010 = mfaba;
			   LBMReal m001 = mfaab;
			   LBMReal m110 = mfbba;
			   LBMReal m101 = mfbab;
			   LBMReal m011 = mfabb;
			   LBMReal m111 = mfbbb;
			   LBMReal m200 = mfcaa - m000 / 3.0;
			   LBMReal m020 = mfaca - m000 / 3.0;
			   LBMReal m002 = mfaac - m000 / 3.0;
			   LBMReal m210 = mfcba - m010 / 3.0;
			   LBMReal m012 = mfabc - m010 / 3.0;
			   LBMReal m201 = mfcab - m001 / 3.0;
			   LBMReal m021 = mfacb - m001 / 3.0;
			   LBMReal m120 = mfbca - m100 / 3.0;
			   LBMReal m102 = mfbac - m100 / 3.0;
			   LBMReal m220 = mfcca - m000 / 9.0;
			   LBMReal m202 = mfcac - m000 / 9.0;
			   LBMReal m022 = mfacc - m000 / 9.0;
			   LBMReal m221 = mfccb - m001 / 9.0;
			   LBMReal m212 = mfcbc - m010 / 9.0;
			   LBMReal m122 = mfbcc - m100 / 9.0;
			   LBMReal m211 = mfcbb - m011 / 3.0;
			   LBMReal m121 = mfbcb - m101 / 3.0;
			   LBMReal m112 = mfbbc - m110 / 3.0;
			   LBMReal m222 = mfccc - m000 / 27.0;


			   m100 = (1.0 - collFactorPhi)*m100 + hSource*dX1_phi*collFactorPhi / 3.0;
			   m010 = (1.0 - collFactorPhi)*m010 + hSource*dX2_phi*collFactorPhi / 3.0;
			   m001 = (1.0 - collFactorPhi)*m001 + hSource*dX3_phi*collFactorPhi / 3.0;

			   m110 = 0.0;
			   m101 = 0.0;
			   m011 = 0.0;

			   m111 = 0.0;

			   //(200)//
			   m200 = m000 / 3.0;
			   m020 = m000 / 3.0;
			   m002 = m000 / 3.0;
			   ////

			   //(210)//
			   m210 = m010 / 3.0;
			   m201 = m001 / 3.0;
			   m120 = m100 / 3.0;


			   m102 = m100 / 3.0;
			   m021 = m001 / 3.0;
			   m012 = m010 / 3.0;
			   ////


			   //(220)//
			   m220 = m000 / 9.0;
			   m202 = m000 / 9.0;
			   m022 = m000 / 9.0;
			   ////

			   //(221)//
			   m221 = m001 / 9.0;
			   m212 = m010 / 9.0;
			   m122 = m100 / 9.0;
			   ////

			   //(211)//
			   m211 = m011 / 3.0;
			   m121 = m101 / 3.0;
			   m112 = m110 / 3.0;
			   ////

			   //(222)//
			   m222 = m000 / 27.0;
			   ////
			   mfaaa = m000 ;
			   mfbaa = m100;
			   mfaba = m010;
			   mfaab = m001;
			   mfbba = m110;
			   mfbab = m101;
			   mfabb = m011;
			   mfbbb = m111;
			   mfcaa = m200;
			   mfaca = m020;
			   mfaac = m002;
			   mfcba = m210;
			   mfabc = m012;
			   mfcab = m201;
			   mfacb = m021;
			   mfbca = m120;
			   mfbac = m102;
			   mfcca = m220;
			   mfcac = m202;
			   mfacc = m022;
			   mfccb = m221;
			   mfcbc = m212;
			   mfbcc = m122;
			   mfcbb = m211;
			   mfbcb = m121;
			   mfbbc = m112;
			   mfccc = m222;
			   ////////////////////////////////////////////////////////////////////////////////////
			   //back
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Z - Dir
			   m0 = mfaac * c1o2 + mfaab * (uz - c1o2) + (mfaaa + one* oMdrho) * (uz2 - uz) * c1o2;
			   m1 = -mfaac - two* mfaab *  uz + mfaaa                * (one - uz2) - one* oMdrho * uz2;
			   m2 = mfaac * c1o2 + mfaab * (uz + c1o2) + (mfaaa + one* oMdrho) * (uz2 + uz) * c1o2;
			   mfaaa = m0;
			   mfaab = m1;
			   mfaac = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfabc * c1o2 + mfabb * (uz - c1o2) + mfaba * (uz2 - uz) * c1o2;
			   m1 = -mfabc - two* mfabb *  uz + mfaba * (one - uz2);
			   m2 = mfabc * c1o2 + mfabb * (uz + c1o2) + mfaba * (uz2 + uz) * c1o2;
			   mfaba = m0;
			   mfabb = m1;
			   mfabc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacc * c1o2 + mfacb * (uz - c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
			   m1 = -mfacc - two* mfacb *  uz + mfaca                  * (one - uz2) - c1o3 * oMdrho * uz2;
			   m2 = mfacc * c1o2 + mfacb * (uz + c1o2) + (mfaca + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
			   mfaca = m0;
			   mfacb = m1;
			   mfacc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfbac * c1o2 + mfbab * (uz - c1o2) + mfbaa * (uz2 - uz) * c1o2;
			   m1 = -mfbac - two* mfbab *  uz + mfbaa * (one - uz2);
			   m2 = mfbac * c1o2 + mfbab * (uz + c1o2) + mfbaa * (uz2 + uz) * c1o2;
			   mfbaa = m0;
			   mfbab = m1;
			   mfbac = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbbc * c1o2 + mfbbb * (uz - c1o2) + mfbba * (uz2 - uz) * c1o2;
			   m1 = -mfbbc - two* mfbbb *  uz + mfbba * (one - uz2);
			   m2 = mfbbc * c1o2 + mfbbb * (uz + c1o2) + mfbba * (uz2 + uz) * c1o2;
			   mfbba = m0;
			   mfbbb = m1;
			   mfbbc = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcc * c1o2 + mfbcb * (uz - c1o2) + mfbca * (uz2 - uz) * c1o2;
			   m1 = -mfbcc - two* mfbcb *  uz + mfbca * (one - uz2);
			   m2 = mfbcc * c1o2 + mfbcb * (uz + c1o2) + mfbca * (uz2 + uz) * c1o2;
			   mfbca = m0;
			   mfbcb = m1;
			   mfbcc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcac * c1o2 + mfcab * (uz - c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 - uz) * c1o2;
			   m1 = -mfcac - two* mfcab *  uz + mfcaa                  * (one - uz2) - c1o3 * oMdrho * uz2;
			   m2 = mfcac * c1o2 + mfcab * (uz + c1o2) + (mfcaa + c1o3 * oMdrho) * (uz2 + uz) * c1o2;
			   mfcaa = m0;
			   mfcab = m1;
			   mfcac = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfcbc * c1o2 + mfcbb * (uz - c1o2) + mfcba * (uz2 - uz) * c1o2;
			   m1 = -mfcbc - two* mfcbb *  uz + mfcba * (one - uz2);
			   m2 = mfcbc * c1o2 + mfcbb * (uz + c1o2) + mfcba * (uz2 + uz) * c1o2;
			   mfcba = m0;
			   mfcbb = m1;
			   mfcbc = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfccb * (uz - c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 - uz) * c1o2;
			   m1 = -mfccc - two* mfccb *  uz + mfcca                  * (one - uz2) - c1o9 * oMdrho * uz2;
			   m2 = mfccc * c1o2 + mfccb * (uz + c1o2) + (mfcca + c1o9 * oMdrho) * (uz2 + uz) * c1o2;
			   mfcca = m0;
			   mfccb = m1;
			   mfccc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // Y - Dir
			   m0 = mfaca * c1o2 + mfaba * (uy - c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfaca - two* mfaba *  uy + mfaaa                  * (one - uy2) - c1o6 * oMdrho * uy2;
			   m2 = mfaca * c1o2 + mfaba * (uy + c1o2) + (mfaaa + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
			   mfaaa = m0;
			   mfaba = m1;
			   mfaca = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacb * c1o2 + mfabb * (uy - c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfacb - two* mfabb *  uy + mfaab                  * (one - uy2) - c2o3 * oMdrho * uy2;
			   m2 = mfacb * c1o2 + mfabb * (uy + c1o2) + (mfaab + c2o3 * oMdrho) * (uy2 + uy) * c1o2;
			   mfaab = m0;
			   mfabb = m1;
			   mfacb = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfacc * c1o2 + mfabc * (uy - c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfacc - two* mfabc *  uy + mfaac                  * (one - uy2) - c1o6 * oMdrho * uy2;
			   m2 = mfacc * c1o2 + mfabc * (uy + c1o2) + (mfaac + c1o6 * oMdrho) * (uy2 + uy) * c1o2;
			   mfaac = m0;
			   mfabc = m1;
			   mfacc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfbca * c1o2 + mfbba * (uy - c1o2) + mfbaa * (uy2 - uy) * c1o2;
			   m1 = -mfbca - two* mfbba *  uy + mfbaa * (one - uy2);
			   m2 = mfbca * c1o2 + mfbba * (uy + c1o2) + mfbaa * (uy2 + uy) * c1o2;
			   mfbaa = m0;
			   mfbba = m1;
			   mfbca = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcb * c1o2 + mfbbb * (uy - c1o2) + mfbab * (uy2 - uy) * c1o2;
			   m1 = -mfbcb - two* mfbbb *  uy + mfbab * (one - uy2);
			   m2 = mfbcb * c1o2 + mfbbb * (uy + c1o2) + mfbab * (uy2 + uy) * c1o2;
			   mfbab = m0;
			   mfbbb = m1;
			   mfbcb = m2;
			   /////////b//////////////////////////////////////////////////////////////////////////
			   m0 = mfbcc * c1o2 + mfbbc * (uy - c1o2) + mfbac * (uy2 - uy) * c1o2;
			   m1 = -mfbcc - two* mfbbc *  uy + mfbac * (one - uy2);
			   m2 = mfbcc * c1o2 + mfbbc * (uy + c1o2) + mfbac * (uy2 + uy) * c1o2;
			   mfbac = m0;
			   mfbbc = m1;
			   mfbcc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcca * c1o2 + mfcba * (uy - c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfcca - two* mfcba *  uy + mfcaa                   * (one - uy2) - c1o18 * oMdrho * uy2;
			   m2 = mfcca * c1o2 + mfcba * (uy + c1o2) + (mfcaa + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
			   mfcaa = m0;
			   mfcba = m1;
			   mfcca = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccb * c1o2 + mfcbb * (uy - c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfccb - two* mfcbb *  uy + mfcab                  * (one - uy2) - c2o9 * oMdrho * uy2;
			   m2 = mfccb * c1o2 + mfcbb * (uy + c1o2) + (mfcab + c2o9 * oMdrho) * (uy2 + uy) * c1o2;
			   mfcab = m0;
			   mfcbb = m1;
			   mfccb = m2;
			   /////////c//////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfcbc * (uy - c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 - uy) * c1o2;
			   m1 = -mfccc - two* mfcbc *  uy + mfcac                   * (one - uy2) - c1o18 * oMdrho * uy2;
			   m2 = mfccc * c1o2 + mfcbc * (uy + c1o2) + (mfcac + c1o18 * oMdrho) * (uy2 + uy) * c1o2;
			   mfcac = m0;
			   mfcbc = m1;
			   mfccc = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
			   ////////////////////////////////////////////////////////////////////////////////////
			   // X - Dir
			   m0 = mfcaa * c1o2 + mfbaa * (ux - c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcaa - two* mfbaa *  ux + mfaaa                   * (one - ux2) - c1o36 * oMdrho * ux2;
			   m2 = mfcaa * c1o2 + mfbaa * (ux + c1o2) + (mfaaa + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
			   mfaaa = m0;
			   mfbaa = m1;
			   mfcaa = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcba * c1o2 + mfbba * (ux - c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcba - two* mfbba *  ux + mfaba                  * (one - ux2) - c1o9 * oMdrho * ux2;
			   m2 = mfcba * c1o2 + mfbba * (ux + c1o2) + (mfaba + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
			   mfaba = m0;
			   mfbba = m1;
			   mfcba = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcca * c1o2 + mfbca * (ux - c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcca - two* mfbca *  ux + mfaca                   * (one - ux2) - c1o36 * oMdrho * ux2;
			   m2 = mfcca * c1o2 + mfbca * (ux + c1o2) + (mfaca + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
			   mfaca = m0;
			   mfbca = m1;
			   mfcca = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcab * c1o2 + mfbab * (ux - c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcab - two* mfbab *  ux + mfaab                  * (one - ux2) - c1o9 * oMdrho * ux2;
			   m2 = mfcab * c1o2 + mfbab * (ux + c1o2) + (mfaab + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
			   mfaab = m0;
			   mfbab = m1;
			   mfcab = m2;
			   ///////////b////////////////////////////////////////////////////////////////////////
			   m0 = mfcbb * c1o2 + mfbbb * (ux - c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcbb - two* mfbbb *  ux + mfabb                  * (one - ux2) - c4o9 * oMdrho * ux2;
			   m2 = mfcbb * c1o2 + mfbbb * (ux + c1o2) + (mfabb + c4o9 * oMdrho) * (ux2 + ux) * c1o2;
			   mfabb = m0;
			   mfbbb = m1;
			   mfcbb = m2;
			   ///////////b////////////////////////////////////////////////////////////////////////
			   m0 = mfccb * c1o2 + mfbcb * (ux - c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfccb - two* mfbcb *  ux + mfacb                  * (one - ux2) - c1o9 * oMdrho * ux2;
			   m2 = mfccb * c1o2 + mfbcb * (ux + c1o2) + (mfacb + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
			   mfacb = m0;
			   mfbcb = m1;
			   mfccb = m2;
			   ////////////////////////////////////////////////////////////////////////////////////
			   ////////////////////////////////////////////////////////////////////////////////////
			   m0 = mfcac * c1o2 + mfbac * (ux - c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcac - two* mfbac *  ux + mfaac                   * (one - ux2) - c1o36 * oMdrho * ux2;
			   m2 = mfcac * c1o2 + mfbac * (ux + c1o2) + (mfaac + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
			   mfaac = m0;
			   mfbac = m1;
			   mfcac = m2;
			   ///////////c////////////////////////////////////////////////////////////////////////
			   m0 = mfcbc * c1o2 + mfbbc * (ux - c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfcbc - two* mfbbc *  ux + mfabc                  * (one - ux2) - c1o9 * oMdrho * ux2;
			   m2 = mfcbc * c1o2 + mfbbc * (ux + c1o2) + (mfabc + c1o9 * oMdrho) * (ux2 + ux) * c1o2;
			   mfabc = m0;
			   mfbbc = m1;
			   mfcbc = m2;
			   ///////////c////////////////////////////////////////////////////////////////////////
			   m0 = mfccc * c1o2 + mfbcc * (ux - c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 - ux) * c1o2;
			   m1 = -mfccc - two* mfbcc *  ux + mfacc                   * (one - ux2) - c1o36 * oMdrho * ux2;
			   m2 = mfccc * c1o2 + mfbcc * (ux + c1o2) + (mfacc + c1o36 * oMdrho) * (ux2 + ux) * c1o2;
			   mfacc = m0;
			   mfbcc = m1;
			   mfccc = m2;

			   ////////////////////////////////////////////////////////////////////////////////////

			   //////////////////////////////////////////////////////////////////////////
			   //proof correctness
			   //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
			   LBMReal drho_post = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca)
				   + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc)
				   + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;
			   //LBMReal dif = fabs(rho - rho_post);
			   dif = drho - drho_post;
#ifdef SINGLEPRECISION
			   if (dif > 10.0E-7 || dif < -10.0E-7)
#else
			   if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
			   {
				   UB_THROW(UbException(UB_EXARGS, "rho=" + UbSystem::toString(drho) + ", rho_post=" + UbSystem::toString(drho_post)
					   + " dif=" + UbSystem::toString(dif)
					   + " rho is not correct for node " + UbSystem::toString(x1) + "," + UbSystem::toString(x2) + "," + UbSystem::toString(x3)));
				   //UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
				   //exit(EXIT_FAILURE);
			   }
#endif			   
			   
			   
			   //////////////////////////////////////////////////////////////////////////
			   //write distribution for H
			   //////////////////////////////////////////////////////////////////////////

			   //(*this->localDistributionsH)(D3Q27System::ET_E,   x1,  x2,  x3) = mfabb;
			   //(*this->localDistributionsH)(D3Q27System::ET_N,   x1,  x2,  x3) = mfbab;
			   //(*this->localDistributionsH)(D3Q27System::ET_T,   x1,  x2,  x3) = mfbba;
			   //(*this->localDistributionsH)(D3Q27System::ET_NE,  x1,  x2,  x3) = mfaab;
			   //(*this->localDistributionsH)(D3Q27System::ET_NW,  x1p, x2,  x3) = mfcab;
			   //(*this->localDistributionsH)(D3Q27System::ET_TE,  x1,  x2,  x3) = mfaba;
			   //(*this->localDistributionsH)(D3Q27System::ET_TW,  x1p, x2,  x3) = mfcba;
			   //(*this->localDistributionsH)(D3Q27System::ET_TN,  x1,  x2,  x3) = mfbaa;
			   //(*this->localDistributionsH)(D3Q27System::ET_TS,  x1,  x2p, x3) = mfbca;
			   //(*this->localDistributionsH)(D3Q27System::ET_TNE, x1,  x2,  x3) = mfaaa;
			   //(*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2,  x3) = mfcaa;
			   //(*this->localDistributionsH)(D3Q27System::ET_TSE, x1,  x2p, x3) = mfaca;
			   //(*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = mfcca;

			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_W,   x1p, x2,  x3 ) = mfcbb;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_S,   x1,  x2p, x3 ) = mfbcb;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_B,   x1,  x2,  x3p) = mfbbc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = mfccb;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = mfacb;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BW,  x1p, x2,  x3p) = mfcbc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BE,  x1,  x2,  x3p) = mfabc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BS,  x1,  x2p, x3p) = mfbcc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BN,  x1,  x2,  x3p) = mfbac;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1,  x2p, x3p) = mfacc;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2,  x3p) = mfcac;
			   //(*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1,  x2,  x3p) = mfaac;

			   //(*this->zeroDistributionsH)(x1,x2,x3) = mfbbb;
			   
*/
			   
			   
/////////////////////   PHASE-FIELD BGK SOLVER ///////////////////////////////

			   h[E  ] = (*this->localDistributionsH)(D3Q27System::ET_E, x1,x2,x3);
			   h[N  ] = (*this->localDistributionsH)(D3Q27System::ET_N,x1,x2,x3); 
			   h[T  ] = (*this->localDistributionsH)(D3Q27System::ET_T,x1,x2,x3);
			   h[NE ] = (*this->localDistributionsH)(D3Q27System::ET_NE,x1,x2,x3);
			   h[NW ] = (*this->localDistributionsH)(D3Q27System::ET_NW,x1p,x2,x3);
			   h[TE ] = (*this->localDistributionsH)(D3Q27System::ET_TE,x1,x2,x3);
			   h[TW ] = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p,x2,x3);
			   h[TN ] = (*this->localDistributionsH)(D3Q27System::ET_TN,x1,x2,x3);
			   h[TS ] = (*this->localDistributionsH)(D3Q27System::ET_TS,x1,x2p,x3);
			   h[TNE] = (*this->localDistributionsH)(D3Q27System::ET_TNE,x1,x2,x3);
			   h[TNW] = (*this->localDistributionsH)(D3Q27System::ET_TNW,x1p,x2,x3);
			   h[TSE] = (*this->localDistributionsH)(D3Q27System::ET_TSE,x1,x2p,x3);
			   h[TSW] = (*this->localDistributionsH)(D3Q27System::ET_TSW,x1p,x2p,x3);

			   h[W  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_W,x1p,x2,x3  );
			   h[S  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_S,x1,x2p,x3  );
			   h[B  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_B,x1,x2,x3p  );
			   h[SW ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,x1p,x2p,x3 );
			   h[SE ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,x1,x2p,x3 );
			   h[BW ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,x1p,x2,x3p );
			   h[BE ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,x1,x2,x3p );
			   h[BS ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,x1,x2p,x3p );
			   h[BN ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,x1,x2,x3p );
			   h[BSW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW,x1p,x2p,x3p);
			   h[BSE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE,x1,x2p,x3p);
			   h[BNW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW,x1p,x2,x3p);
			   h[BNE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE,x1,x2,x3p);

			   h[ZERO] = (*this->zeroDistributionsH)(x1,x2,x3);
			   
			   
			   //LBMReal denom = sqrt(dX1_phi*dX1_phi + dX2_phi*dX2_phi + dX3_phi*dX3_phi) + 1e-15;
			   //LBMReal di = sqrt(8*kappa/beta);
			   LBMReal tauH1 = 3.0*mob + 0.5;
			   for (int dir = STARTF; dir < (ENDF+1); dir++)
			   {
				   LBMReal velProd = DX1[dir]*ux + DX2[dir]*uy + DX3[dir]*uz;
				   LBMReal velSq1 = velProd*velProd;
				   LBMReal hEq, gEq;
				   
				   if (dir != ZERO)
				   {
					   //LBMReal dirGrad_phi = DX1[dir]*dX1_phi+DX2[dir]*dX2_phi+DX3[dir]*dX3_phi;
					   LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]])/2.0;
					   LBMReal hSource = (tauH - 0.5)*(1.0 - phi[ZERO])*(phi[ZERO])*(dirGrad_phi)/denom; // + phi[ZERO]*(dxux + dyuy + dzuz);
						   
					   //LBMReal hSource =((phi[ZERO]>phiH || phi[ZERO]<phiL) ? 0.1 : 1.0) * 3.0*mob*(-4.0)/di*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)*(dirGrad_phi)/denom;
					   //LBMReal hSource = 3.0*mob*(-4.0)/di*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)*(dirGrad_phi)/denom;
					   hEq = phi[ZERO]*WEIGTH[dir]*(1.0 + 3.0*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2)) + hSource*WEIGTH[dir];
					   //gEq = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx2+vy2+vz2))*c1o3 + (p1-rho*c1o3)*WEIGTH[dir];
					   //h[dir] = hEq; //h[dir] - (h[dir] - hEq)/(tauH + 0.5));  /// This corresponds with the collision factor of 1.0 which equals (tauH + 0.5). 
					   h[dir] = h[dir] - (h[dir] - hEq)/(tauH); // + WEIGTH[dir]*phi[ZERO]*(dxux + dyuy + dzuz);
					   //h[dir] = h[dir] - (h[dir] - hEq)/(tauH1);
					   //g[dir] = g[dir] - collFactorM*(g[dir]-gEq) + 0.5*forcingTerm[dir];

				   } 
				   else
				   {
					   hEq = phi[ZERO]*WEIGTH[ZERO]*(1.0 - 1.5*(ux2+uy2+uz2));
					   //gEq = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx2+vy2+vz2))*c1o3 + (p1-rho*c1o3)*WEIGTH[dir];
					   //h[dir] = hEq;
					   h[ZERO] = h[ZERO] - (h[ZERO] - hEq)/(tauH); // + WEIGTH[ZERO]*phi[ZERO]*(dxux + dyuy + dzuz);
					   //g[dir] = g[dir] - collFactorM*(g[dir]-gEq) + 0.5*forcingTerm[dir];
				   }
			   }
			   
			   
			   (*this->localDistributionsH)(D3Q27System::ET_E,   x1,  x2,  x3) = h[D3Q27System::INV_E];
			   (*this->localDistributionsH)(D3Q27System::ET_N,   x1,  x2,  x3) = h[D3Q27System::INV_N];
			   (*this->localDistributionsH)(D3Q27System::ET_T,   x1,  x2,  x3) = h[D3Q27System::INV_T];
			   (*this->localDistributionsH)(D3Q27System::ET_NE,  x1,  x2,  x3) = h[D3Q27System::INV_NE];
			   (*this->localDistributionsH)(D3Q27System::ET_NW,  x1p, x2,  x3) = h[D3Q27System::INV_NW];
			   (*this->localDistributionsH)(D3Q27System::ET_TE,  x1,  x2,  x3) = h[D3Q27System::INV_TE];
			   (*this->localDistributionsH)(D3Q27System::ET_TW,  x1p, x2,  x3) = h[D3Q27System::INV_TW];
			   (*this->localDistributionsH)(D3Q27System::ET_TN,  x1,  x2,  x3) = h[D3Q27System::INV_TN];
			   (*this->localDistributionsH)(D3Q27System::ET_TS,  x1,  x2p, x3) = h[D3Q27System::INV_TS];
			   (*this->localDistributionsH)(D3Q27System::ET_TNE, x1,  x2,  x3) = h[D3Q27System::INV_TNE];
			   (*this->localDistributionsH)(D3Q27System::ET_TNW, x1p, x2,  x3) = h[D3Q27System::INV_TNW];
			   (*this->localDistributionsH)(D3Q27System::ET_TSE, x1,  x2p, x3) = h[D3Q27System::INV_TSE];
			   (*this->localDistributionsH)(D3Q27System::ET_TSW, x1p, x2p, x3) = h[D3Q27System::INV_TSW];
			   
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_W,   x1p, x2,  x3 ) = h[D3Q27System::INV_W ];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_S,   x1,  x2p, x3 ) = h[D3Q27System::INV_S ];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_B,   x1,  x2,  x3p) = h[D3Q27System::INV_B ];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = h[D3Q27System::INV_SW];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = h[D3Q27System::INV_SE];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,  x1p, x2,  x3p) = h[D3Q27System::INV_BW];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,  x1,  x2,  x3p) = h[D3Q27System::INV_BE];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,  x1,  x2p, x3p) = h[D3Q27System::INV_BS];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,  x1,  x2,  x3p) = h[D3Q27System::INV_BN];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW, x1p, x2p, x3p) = h[D3Q27System::INV_BSW];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE, x1,  x2p, x3p) = h[D3Q27System::INV_BSE];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW, x1p, x2,  x3p) = h[D3Q27System::INV_BNW];
			   (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE, x1,  x2,  x3p) = h[D3Q27System::INV_BNE];
			   
			   (*this->zeroDistributionsH)(x1,x2,x3) = h[D3Q27System::ZERO];			   
			   
			   
/////////////////////   END OF OLD BGK SOLVER ///////////////////////////////
			   





               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
   dataSet->setPhaseField(divU);
   
   }
}
//////////////////////////////////////////////////////////////////////////
double MultiphaseCumulantLBMKernel::getCallculationTime()
{
   //return timer.getDuration();
   return timer.getTotalTime();
}
//////////////////////////////////////////////////////////////////////////

LBMReal MultiphaseCumulantLBMKernel::gradX1_phi()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX1[k]*phi[k];
	}
	return 3.0*sum;
}

LBMReal MultiphaseCumulantLBMKernel::gradX2_phi()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX2[k]*phi[k];
	}
	return 3.0*sum;
}

LBMReal MultiphaseCumulantLBMKernel::gradX3_phi()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX3[k]*phi[k];
	}
	return 3.0*sum;
}

LBMReal MultiphaseCumulantLBMKernel::nabla2_phi()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*(phi[k] - phi[ZERO]);
	}
	return 6.0*sum;
}
///// Commnets neeeded ////////

void MultiphaseCumulantLBMKernel::computePhasefield()
{
	using namespace D3Q27System;
	DistributionArray3DPtr distributionsH = dataSet->getHdistributions();

	//const int bcArrayMaxX1 = (int)distributionsH->getNX1();
	//const int bcArrayMaxX2 = (int)distributionsH->getNX2();
	//const int bcArrayMaxX3 = (int)distributionsH->getNX3();

	int minX1 = ghostLayerWidth;
	int minX2 = ghostLayerWidth;
	int minX3 = ghostLayerWidth;
	int maxX1 = (int)distributionsH->getNX1() - ghostLayerWidth;
	int maxX2 = (int)distributionsH->getNX2() - ghostLayerWidth;
	int maxX3 = (int)distributionsH->getNX3() - ghostLayerWidth;

	//------------- Computing the phase-field ------------------
	for(int x3 = minX3; x3 < maxX3; x3++)
	{
		for(int x2 = minX2; x2 < maxX2; x2++)
		{
			for(int x1 = minX1; x1 < maxX1; x1++)
			{
				//if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
				{
					int x1p = x1 + 1;
					int x2p = x2 + 1;
					int x3p = x3 + 1;

					h[E  ] = (*this->localDistributionsH)(D3Q27System::ET_E, x1,x2,x3);
					h[N  ] = (*this->localDistributionsH)(D3Q27System::ET_N,x1,x2,x3); 
					h[T  ] = (*this->localDistributionsH)(D3Q27System::ET_T,x1,x2,x3);
					h[NE ] = (*this->localDistributionsH)(D3Q27System::ET_NE,x1,x2,x3);
					h[NW ] = (*this->localDistributionsH)(D3Q27System::ET_NW,x1p,x2,x3);
					h[TE ] = (*this->localDistributionsH)(D3Q27System::ET_TE,x1,x2,x3);
					h[TW ] = (*this->localDistributionsH)(D3Q27System::ET_TW, x1p,x2,x3);
					h[TN ] = (*this->localDistributionsH)(D3Q27System::ET_TN,x1,x2,x3);
					h[TS ] = (*this->localDistributionsH)(D3Q27System::ET_TS,x1,x2p,x3);
					h[TNE] = (*this->localDistributionsH)(D3Q27System::ET_TNE,x1,x2,x3);
					h[TNW] = (*this->localDistributionsH)(D3Q27System::ET_TNW,x1p,x2,x3);
					h[TSE] = (*this->localDistributionsH)(D3Q27System::ET_TSE,x1,x2p,x3);
					h[TSW] = (*this->localDistributionsH)(D3Q27System::ET_TSW,x1p,x2p,x3);

					h[W  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_W,x1p,x2,x3  );
					h[S  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_S,x1,x2p,x3  );
					h[B  ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_B,x1,x2,x3p  );
					h[SW ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_SW,x1p,x2p,x3 );
					h[SE ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_SE,x1,x2p,x3 );
					h[BW ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BW,x1p,x2,x3p );
					h[BE ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BE,x1,x2,x3p );
					h[BS ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BS,x1,x2p,x3p );
					h[BN ] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BN,x1,x2,x3p );
					h[BSW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSW,x1p,x2p,x3p);
					h[BSE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BSE,x1,x2p,x3p);
					h[BNW] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNW,x1p,x2,x3p);
					h[BNE] = (*this->nonLocalDistributionsH)(D3Q27System::ET_BNE,x1,x2,x3p);

					h[ZERO] = (*this->zeroDistributionsH)(x1,x2,x3);

					/*(*this->phaseField)(x1,x2,x3) = h[ZERO] + h[E] + h[W] + h[N] + h[S] + h[T] + h[B] + h[NE] + h[SW] + h[SE] + h[NW] + h[TE] + h[BW] + 
						h[BE] + h[TW] + h[TN] + h[BS] + h[BN] + h[TS] + h[TNE] + h[TNW] + h[TSE] + h[TSW] + h[BNE] + h[BNW] + h[BSE] + h[BSW];*/

				}
			}
		}
	}
	//----------------------------------------------------------
	
/*
	/////// Filling ghost nodes for FD computations //////////
	for(int x1 = minX1; x1 < maxX1; x1++)
	{
		for(int x2 = minX2; x2 < maxX2; x2++)
		{
			int x3 = 0;
			(*phaseField)(x1, x2, x3) = (*phaseField)(x1, x2, maxX3-1);
			x3 = maxX3;
			(*phaseField)(x1, x2, x3) = (*phaseField)(x1, x2, minX3);
		}
	}
	for(int x2 = minX2; x2 < maxX2; x2++)
	{
		for(int x3 = minX3; x3 < maxX3; x3++)
		{
			int x1 = 0;
			(*phaseField)(x1, x2, x3) = (*phaseField)(maxX1-1, x2, x3);
			x1 = maxX1;
			(*phaseField)(x1, x2, x3) = (*phaseField)(minX1, x2, x3);
		}
	}
	for(int x1 = minX1; x1 < maxX1; x1++)
	{
		for(int x3 = minX3; x3 < maxX3; x3++)
		{
			int x2 = 0;
			(*phaseField)(x1, x2, x3) = (*phaseField)(x1, maxX2-1, x3);
			x2 = maxX2;
			(*phaseField)(x1, x2, x3) = (*phaseField)(x1, minX2, x3);
		}
	}
	(*phaseField)(0, 0,     0    ) = (*phaseField)(maxX1-1, maxX2-1, maxX3-1);
	(*phaseField)(0, 0,     maxX3) = (*phaseField)(maxX1-1, maxX2-1, minX3  );
	(*phaseField)(0, maxX2, 0    ) = (*phaseField)(maxX1-1, minX2, maxX3-1  );
	(*phaseField)(0, maxX2, maxX3) = (*phaseField)(maxX1-1, minX2, minX3    );

	(*phaseField)(maxX1, 0,     0    ) = (*phaseField)(minX1, maxX2-1, maxX3-1);
	(*phaseField)(maxX1, 0,     maxX3) = (*phaseField)(minX1, maxX2-1, minX3  );
	(*phaseField)(maxX1, maxX2, 0    ) = (*phaseField)(minX1, minX2, maxX3-1  );
	(*phaseField)(maxX1, maxX2, maxX3) = (*phaseField)(minX1, minX2, minX3    );

	///////////////////////////////////////////////////////// 
*/
}

void MultiphaseCumulantLBMKernel::findNeighbors(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr ph, int x1, int x2, int x3)
{
	using namespace D3Q27System;
	
	BCArray3DPtr bcArray = this->getBCProcessor()->getBCArray();

	phi[ZERO] = (*ph)(x1,x2,x3);
	
	LBMReal a = -0.5*sqrt(2*beta/kappa)*cos(contactAngle*PI/180);
	LBMReal a1 = 1 + a;
	
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		
		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]))
		{
			phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
		} 
		else
		{
			/*
			if (phi[ZERO] < 1e-2)
			{
				phi[k] = (*ph)(x1 + DX1[INVDIR[k]], x2 + DX2[INVDIR[k]], x3 + DX3[INVDIR[k]]);
			}
			else
			{
				LBMReal phi_f = (*ph)(x1 + DX1[k], x2, x3 + DX3[k]);
				phi[k] = (a1 - sqrt(a1*a1 - 4*a*phi_f) )/a - phi_f;
			}
			*/
			
			phi[k] = (*ph)(x1, x2, x3);

			//if (bcArray->isSolid(x1 + DX1[k], x2, x3))
			//{
			//	phi[k] = (*ph)(x1, x2, x3);
			//	//if (!bcArray->isSolid(x1 , x2 + DX2[k], x3 + DX3[k]))
			//	//{
			//	//	//phi[k] = (*ph)(x1 , x2 + DX2[k], x3 + DX3[k]);
			//	//	LBMReal phi_f = (*ph)(x1 , x2 + DX2[k], x3 + DX3[k]);
			//	//	phi[k] = (a1 - sqrt(a1*a1 - 4*a*phi_f) )/a - phi_f;
			//	//} 
			//	//else
			//	//{
			//	//	phi[k] = (*ph)(x1, x2, x3);
			//	//}
			//}
			//
			//if (bcArray->isSolid(x1 , x2 + DX2[k], x3))
			//{
			//	phi[k] = (*ph)(x1, x2, x3);
			//	//if (!bcArray->isSolid(x1 + DX1[k], x2 , x3 + DX3[k]))
			//	//{
			//	//	//phi[k] = (*ph)(x1 + DX1[k], x2 , x3 + DX3[k]);
			//	//	LBMReal phi_f = (*ph)(x1 + DX1[k], x2 , x3 + DX3[k]);
			//	//	phi[k] = (a1 - sqrt(a1*a1 - 4*a*phi_f) )/a - phi_f;
			//	//} 
			//	//else
			//	//{
			//	//	phi[k] = (*ph)(x1, x2, x3);
			//	//}
			//}


			//if (bcArray->isSolid(x1 , x2, x3+ DX3[k]))
			//{
			//	if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3))
			//	{
			//		//phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3);
			//		LBMReal phi_f = (*ph)(x1 + DX1[k], x2 + DX2[k], x3);
			//		phi[k] = (a1 - sqrt(a1*a1 - 4*a*phi_f) )/a - phi_f;
			//	} 
			//	else
			//	{
			//		phi[k] = (*ph)(x1, x2, x3);
			//	}
			//}


			/*if (bcArray->isSolid(x1 + DX1[k], x2, x3)) phi[k] = (*ph)(x1 , x2 + DX2[k], x3 + DX3[k]);
			if (bcArray->isSolid(x1 , x2 + DX2[k], x3)) phi[k] = (*ph)(x1 + DX1[k], x2 , x3 + DX3[k]);
			if (bcArray->isSolid(x1 , x2, x3+ DX3[k])) phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 );*/

			/*if (phi[ZERO] < 0.00001)
			{
			phi[k] = 0.0;
			} 
			else
			{
			phi[k] = 0.5;
			}*/
			
			//phi[k] = 0.5;
			//phi[k] = (*ph)(x1, x2, x3);
			//phi[k] = (*ph)(x1 + DX1[INVDIR[k]], x2 + DX2[INVDIR[k]], x3 + DX3[INVDIR[k]]);
			
			
		}
	}

	/*
	phi[E  ] = (*ph)(x1 + DX1[E  ], x2 + DX2[E  ], x3 + DX3[E  ]);
	phi[N  ] = (*ph)(x1 + DX1[N  ], x2 + DX2[N  ], x3 + DX3[N  ]);
	phi[T  ] = (*ph)(x1 + DX1[T  ], x2 + DX2[T  ], x3 + DX3[T  ]);
	phi[W  ] = (*ph)(x1 + DX1[W  ], x2 + DX2[W  ], x3 + DX3[W  ]);
	phi[S  ] = (*ph)(x1 + DX1[S  ], x2 + DX2[S  ], x3 + DX3[S  ]);
	phi[B  ] = (*ph)(x1 + DX1[B  ], x2 + DX2[B  ], x3 + DX3[B  ]);
	phi[NE ] = (*ph)(x1 + DX1[NE ], x2 + DX2[NE ], x3 + DX3[NE ]);
	phi[NW ] = (*ph)(x1 + DX1[NW ], x2 + DX2[NW ], x3 + DX3[NW ]);
	phi[TE ] = (*ph)(x1 + DX1[TE ], x2 + DX2[TE ], x3 + DX3[TE ]);
	phi[TW ] = (*ph)(x1 + DX1[TW ], x2 + DX2[TW ], x3 + DX3[TW ]);
	phi[TN ] = (*ph)(x1 + DX1[TN ], x2 + DX2[TN ], x3 + DX3[TN ]);
	phi[TS ] = (*ph)(x1 + DX1[TS ], x2 + DX2[TS ], x3 + DX3[TS ]);
	phi[SW ] = (*ph)(x1 + DX1[SW ], x2 + DX2[SW ], x3 + DX3[SW ]);
	phi[SE ] = (*ph)(x1 + DX1[SE ], x2 + DX2[SE ], x3 + DX3[SE ]);
	phi[BW ] = (*ph)(x1 + DX1[BW ], x2 + DX2[BW ], x3 + DX3[BW ]);
	phi[BE ] = (*ph)(x1 + DX1[BE ], x2 + DX2[BE ], x3 + DX3[BE ]);
	phi[BS ] = (*ph)(x1 + DX1[BS ], x2 + DX2[BS ], x3 + DX3[BS ]);
	phi[BN ] = (*ph)(x1 + DX1[BN ], x2 + DX2[BN ], x3 + DX3[BN ]);
	phi[BSW] = (*ph)(x1 + DX1[BSW], x2 + DX2[BSW], x3 + DX3[BSW]);
	phi[BSE] = (*ph)(x1 + DX1[BSE], x2 + DX2[BSE], x3 + DX3[BSE]);
	phi[BNW] = (*ph)(x1 + DX1[BNW], x2 + DX2[BNW], x3 + DX3[BNW]);
	phi[BNE] = (*ph)(x1 + DX1[BNE], x2 + DX2[BNE], x3 + DX3[BNE]);
	phi[TNE] = (*ph)(x1 + DX1[TNE], x2 + DX2[TNE], x3 + DX3[TNE]);
	phi[TNW] = (*ph)(x1 + DX1[TNW], x2 + DX2[TNW], x3 + DX3[TNW]);
	phi[TSE] = (*ph)(x1 + DX1[TSE], x2 + DX2[TSE], x3 + DX3[TSE]);
	phi[TSW] = (*ph)(x1 + DX1[TSW], x2 + DX2[TSW], x3 + DX3[TSW]);
	*/
}

void MultiphaseCumulantLBMKernel::swapDistributions()
{
   dataSet->getFdistributions()->swap();
   dataSet->getHdistributions()->swap();
   //computePhasefield();
}