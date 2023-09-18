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
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr p1Field(new CbArray3D<LBMReal,IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, -999.0));
   CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr p1Field_filtered(new CbArray3D<LBMReal,IndexerX3X2X1>(bcArrayMaxX1, bcArrayMaxX2, bcArrayMaxX3, -999.0));
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
				if ((*phaseField)(x1,x2,x3) > 1.0) (*phaseField)(x1,x2,x3)=1.0;

				mfcbb = (*this->localDistributionsF)(D3Q27System::ET_E, x1,x2,x3);
				mfbcb = (*this->localDistributionsF)(D3Q27System::ET_N,x1,x2,x3); 
				mfbbc = (*this->localDistributionsF)(D3Q27System::ET_T,x1,x2,x3);
				mfccb = (*this->localDistributionsF)(D3Q27System::ET_NE,x1,x2,x3);
				mfacb = (*this->localDistributionsF)(D3Q27System::ET_NW,x1p,x2,x3);
				mfcbc = (*this->localDistributionsF)(D3Q27System::ET_TE,x1,x2,x3);
				mfabc = (*this->localDistributionsF)(D3Q27System::ET_TW, x1p,x2,x3);
				mfbcc = (*this->localDistributionsF)(D3Q27System::ET_TN,x1,x2,x3);
				mfbac = (*this->localDistributionsF)(D3Q27System::ET_TS,x1,x2p,x3);
				mfccc = (*this->localDistributionsF)(D3Q27System::ET_TNE,x1,x2,x3);
				mfacc = (*this->localDistributionsF)(D3Q27System::ET_TNW,x1p,x2,x3);
				mfcac = (*this->localDistributionsF)(D3Q27System::ET_TSE,x1,x2p,x3);
				mfaac = (*this->localDistributionsF)(D3Q27System::ET_TSW,x1p,x2p,x3);
				mfabb = (*this->nonLocalDistributionsF)(D3Q27System::ET_W,x1p,x2,x3  );
				mfbab = (*this->nonLocalDistributionsF)(D3Q27System::ET_S,x1,x2p,x3  );
				mfbba = (*this->nonLocalDistributionsF)(D3Q27System::ET_B,x1,x2,x3p  );
				mfaab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SW,x1p,x2p,x3 );
				mfcab = (*this->nonLocalDistributionsF)(D3Q27System::ET_SE,x1,x2p,x3 );
				mfaba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BW,x1p,x2,x3p );
				mfcba = (*this->nonLocalDistributionsF)(D3Q27System::ET_BE,x1,x2,x3p );
				mfbaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BS,x1,x2p,x3p );
				mfbca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BN,x1,x2,x3p );
				mfaaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW,x1p,x2p,x3p);
				mfcaa = (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE,x1,x2p,x3p);
				mfaca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW,x1p,x2,x3p);
				mfcca = (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE,x1,x2,x3p);
				mfbbb = (*this->zeroDistributionsF)(x1,x2,x3);
				
				(*p1Field)(x1,x2,x3) = ((mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
					+(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
					+(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb);
				(*p1Field_filtered)(x1,x2,x3) = (*p1Field)(x1,x2,x3);
			}
		 }
	  }
   }
   
   LBMReal sum1, AA;
   for(int x3 = minX3; x3 < maxX3; x3++)
   {
	   for(int x2 = minX2; x2 < maxX2; x2++)
	   {
		   for(int x1 = minX1; x1 < maxX1; x1++)
		   {
			   if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
			   {
				   int cnum = 0;
				   sum1 = 0.0;
				   for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
				   //for (int k = FSTARTDIR ; k <= 5 ; k++)
				   {
					   if(!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]))
					   {
						   //cnum++;
						   sum1 += (*p1Field)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])*WEIGTH[k];

					   }
					   else
					   {
						   //cnum++;
						   //sum1 += (*p1Field)(x1 + DX1[INVDIR[k]], x2 + DX2[INVDIR[k]], x3 + DX3[INVDIR[k]])*WEIGTH[k];
						   //std::cout << x1 << '  ' << x2  << '  ' << x3 << '  '<< k << std::endl;
						   sum1 += ((*p1Field)(x1, x2, x3))*WEIGTH[k];
					   }
				   }
				   //LBMReal av = (sum1+(*p1Field)(x1, x2, x3))/7.0;
				   //(*p1Field_filtered)(x1, x2, x3) = av;
				   
				   /*LBMReal d0 = 2.0/9.0;
				   LBMReal d1 = -1.0/9.6;
				   LBMReal wf = 0.25;
				   AA = 3*d0*(*p1Field)(x1, x2, x3) + d1*((*p1Field)(x1+1, x2, x3)+(*p1Field)(x1-1, x2, x3)) 
					   + d1*((*p1Field)(x1, x2+1, x3)+(*p1Field)(x1, x2-1, x3)) + d1*((*p1Field)(x1, x2, x3+1)+(*p1Field)(x1, x2, x3-1));
				   (*p1Field_filtered)(x1, x2, x3) = (*p1Field)(x1, x2, x3) - AA*wf;*/
				   
				   (*p1Field_filtered)(x1, x2, x3) = ((*p1Field)(x1, x2, x3))*WEIGTH[ZERO] + sum1;
				   //(*p1Field)(x1, x2, x3) = ((*p1Field)(x1, x2, x3))*WEIGTH[ZERO] + sum1;

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
			   findNeighbors(phaseField, p1Field_filtered, x1, x2, x3);
			   //findNeighbors(phaseField, p1Field, x1, x2, x3);

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

			   g[E  ]  = mfcbb;
			   g[N  ]  = mfbcb;
			   g[T  ]  = mfbbc;
			   g[NE ]  = mfccb;
			   g[NW ]  = mfacb;
			   g[TE ]  = mfcbc;
			   g[TW ]  = mfabc;
			   g[TN ]  = mfbcc;
			   g[TS ]  = mfbac;
			   g[TNE]  = mfccc;
			   g[TNW]  = mfacc;
			   g[TSE]  = mfcac;
			   g[TSW]  = mfaac;
			   g[W  ]  = mfabb;
			   g[S  ]  = mfbab;
			   g[B  ]  = mfbba;
			   g[SW ]  = mfaab;
			   g[SE ]  = mfcab;
			   g[BW ]  = mfaba;
			   g[BE ]  = mfcba;
			   g[BS ]  = mfbaa;
			   g[BN ]  = mfbca;
			   g[BSW]  = mfaaa;
			   g[BSE]  = mfcaa;
			   g[BNW]  = mfaca;
			   g[BNE]  = mfcca;
			   g[ZERO] = mfbbb;

			   LBMReal rhoH = 997.0;
			   LBMReal rhoL = rhoH/densityRatio;

			   //LBMReal rhoToPhi = (1.0 - 1.0/densityRatio);
			   LBMReal rhoToPhi = (rhoH - rhoL)/(phiH - phiL);

			   //collFactorM = phi[ZERO]*collFactorL + (1-phi[ZERO])*collFactorG;
			   //collFactorM = phi[ZERO]*collFactorG + (1-phi[ZERO])*collFactorL;
			   
			   //LBMReal tauH = 1.0;
			   LBMReal di = sqrt(8*kappa/beta);
			   
			   LBMReal dX1_phi = gradX1_phi();
			   LBMReal dX2_phi = gradX2_phi();
			   LBMReal dX3_phi = gradX3_phi();

			   LBMReal dX1_pr1 = gradX1_pr1();
			   LBMReal dX2_pr1 = gradX2_pr1();
			   LBMReal dX3_pr1 = gradX3_pr1();
			   
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
			   //mu = mu*10;
			   
			   //LBMReal rhoToPhi = (1.0/densityRatio - 1.0);
			   
			   			   

			   //----------- Calculating Macroscopic Values -------------

			   //LBMReal rho = phi[ZERO] + (1.0 - phi[ZERO])*1.0/densityRatio;
			   LBMReal rho = rhoH + rhoToPhi*(phi[ZERO] - phiH);
			   if (phi[ZERO] > 1.0) rho = rhoH;
			   //LBMReal rho = phi[ZERO]*1.0/densityRatio + (1.0 - phi[ZERO]);

			   if (withForcing)
			   {
				   //muX1 = static_cast<double>(x1-1+ix1*maxX1);
				   //muX2 = static_cast<double>(x2-1+ix2*maxX2);
				   //muX3 = static_cast<double>(x3-1+ix3*maxX3);

				   forcingX1 = muForcingX1.Eval();
				   forcingX2 = muForcingX2.Eval();
				   forcingX3 = muForcingX3.Eval();

				   //LBMReal rho_m = 1.0/densityRatio;
				   //LBMReal rho_m = (rhoL + rhoH)/2.0;
				   LBMReal rho_m = rhoH;
				   forcingX1 = forcingX1*(rho-rho_m);
				   forcingX2 = forcingX2*(rho-rho_m);
				   forcingX3 = forcingX3*(rho-rho_m);

				   //ux += forcingX1*deltaT*0.5; // X
				   //uy += forcingX2*deltaT*0.5; // Y
				   //uz += forcingX3*deltaT*0.5; // Z
			   }

			   /*LBMReal p1 = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
				   +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
				   +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;
			   
			   p1 = p1*rho*c1o3;*/
			   
			   LBMReal rho2 = rho*rho;

			   /*LBMReal Sx = rhoToPhi*dX1_phi*p1/rho2;
			   LBMReal Sy = rhoToPhi*dX2_phi*p1/rho2;
			   LBMReal Sz = rhoToPhi*dX3_phi*p1/rho2;*/
			   
			   LBMReal Sx = -1.0*dX1_pr1/rho;
			   LBMReal Sy = -1.0*dX2_pr1/rho;
			   LBMReal Sz = -1.0*dX3_pr1/rho;

			   
			   /*if ((phi[ZERO] < 0.1)||(phi[ZERO] > 0.9))
			   {
			   Sx = 0.0;
			   Sy = 0.0;
			   Sz = 0.0;
			   }*/

			   LBMReal ux = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
				   (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
				   (mfcbb-mfabb))  + (mu*dX1_phi + forcingX1)/(2*rho) + 0.5*Sx;

			   LBMReal uy = ((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
				   (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				   (mfbcb-mfbab))  + (mu*dX2_phi + forcingX2)/(2*rho) + 0.5*Sy;

			   LBMReal uz = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
				   (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				   (mfbbc-mfbba))  + (mu*dX3_phi + forcingX3)/(2*rho) + 0.5*Sz;

			   
			   //+ (ux*rhoToPhi*dX1_phi*c1o3 + uy*rhoToPhi*dX2_phi*c1o3 + uz*rhoToPhi*dX3_phi*c1o3)/2.0;
			   
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
				   //LBMReal gamma = WEIGTH[dir]*(1.0 + 3*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2));
				   LBMReal gamma = WEIGTH[dir]*(3*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2));
				   
				   //forcingTerm[dir] = (DX1[dir] - ux)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX1_phi + gamma*mu*dX1_phi) + 
					//   (DX2[dir] - uy)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX2_phi + gamma*mu*dX2_phi) + 
					//   (DX3[dir] - uz)*((gamma - WEIGTH[dir])*c1o3*rhoToPhi*dX3_phi + gamma*mu*dX3_phi);
				   
				   LBMReal fac1 = (gamma - WEIGTH[dir])*c1o3*rhoToPhi;
				   LBMReal dirGrad_phi = (phi[dir] - phi[INVDIR[dir]])/2.0;
				   //LBMReal dirGrad_phi = DX1[dir]*dX1_phi + DX2[dir]*dX2_phi + DX3[dir]*dX3_phi;
				   
				   /*forcingTerm[dir] =  (- (ux)*(fac1*dX1_phi + gamma*mu*dX1_phi) - 
				   (uy)*(fac1*dX2_phi + gamma*mu*dX2_phi) - 
				   (uz)*(fac1*dX3_phi + gamma*mu*dX3_phi)) + (fac1*dirGrad_phi + gamma*mu*dirGrad_phi + DX1[dir]*forcingX1 + DX2[dir]*forcingX2 + DX3[dir]*forcingX3);*/
				   
				   
				   forcingTerm[dir] =  (-ux)*     (gamma*(mu*dX1_phi/rho + forcingX1/rho + Sx)) +
					   				   (-uy)*     (gamma*(mu*dX2_phi/rho + forcingX2/rho + Sy)) +
					   				   (-uz)*     (gamma*(mu*dX3_phi/rho + forcingX3/rho + Sz)) +
									   (DX1[dir])*(gamma*(mu*dX1_phi/rho + forcingX1/rho + Sx)) +
									   (DX2[dir])*(gamma*(mu*dX2_phi/rho + forcingX2/rho + Sy)) +
									   (DX3[dir])*(gamma*(mu*dX3_phi/rho + forcingX3/rho + Sz));

			   }

			   LBMReal gamma = WEIGTH[ZERO]*(-1.5*(ux2+uy2+uz2));
			   /*forcingTerm[ZERO] = -(ux)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX1_phi + gamma*mu*dX1_phi) - 
			   (uy)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX2_phi + gamma*mu*dX2_phi) - 
			   (uz)*((gamma - WEIGTH[ZERO])*c1o3*rhoToPhi*dX3_phi + gamma*mu*dX3_phi);*/
			   LBMReal fac1 = (gamma - WEIGTH[ZERO])*c1o3*rhoToPhi;
			   forcingTerm[ZERO] = (-ux)*(gamma*(mu*dX1_phi/rho + forcingX1/rho + Sx)) +
								   (-uy)*(gamma*(mu*dX2_phi/rho + forcingX2/rho + Sy)) +
								   (-uz)*(gamma*(mu*dX3_phi/rho + forcingX3/rho + Sz));

			   //--------------------------------------------------------

			   g[E  ]  += 0.5*forcingTerm[E  ] ;
			   g[N  ]  += 0.5*forcingTerm[N  ] ;
			   g[T  ]  += 0.5*forcingTerm[T  ] ;
			   g[NE ]  += 0.5*forcingTerm[NE ] ;
			   g[NW ]  += 0.5*forcingTerm[NW ] ;
			   g[TE ]  += 0.5*forcingTerm[TE ] ;
			   g[TW ]  += 0.5*forcingTerm[TW ] ;
			   g[TN ]  += 0.5*forcingTerm[TN ] ;
			   g[TS ]  += 0.5*forcingTerm[TS ] ;
			   g[TNE]  += 0.5*forcingTerm[TNE] ;
			   g[TNW]  += 0.5*forcingTerm[TNW] ;
			   g[TSE]  += 0.5*forcingTerm[TSE] ;
			   g[TSW]  += 0.5*forcingTerm[TSW] ;
			   g[W  ]  += 0.5*forcingTerm[W  ] ;
			   g[S  ]  += 0.5*forcingTerm[S  ] ;
			   g[B  ]  += 0.5*forcingTerm[B  ] ;
			   g[SW ]  += 0.5*forcingTerm[SW ] ;
			   g[SE ]  += 0.5*forcingTerm[SE ] ;
			   g[BW ]  += 0.5*forcingTerm[BW ] ;
			   g[BE ]  += 0.5*forcingTerm[BE ] ;
			   g[BS ]  += 0.5*forcingTerm[BS ] ;
			   g[BN ]  += 0.5*forcingTerm[BN ] ;
			   g[BSW]  += 0.5*forcingTerm[BSW] ;
			   g[BSE]  += 0.5*forcingTerm[BSE] ;
			   g[BNW]  += 0.5*forcingTerm[BNW] ;
			   g[BNE]  += 0.5*forcingTerm[BNE] ;
			   g[ZERO] += 0.5*forcingTerm[ZERO];
			   
			   for (int dir = STARTF; dir < (ENDF+1); dir++)
			   {
				   LBMReal velProd = DX1[dir]*ux + DX2[dir]*uy + DX3[dir]*uz;
				   LBMReal velSq1 = velProd*velProd;
				   LBMReal gamma = WEIGTH[dir]*(3*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2));
				   LBMReal hEq, gEq;

				   if (dir != ZERO)
				   {

					   //gEq = p1*WEIGTH[dir]/(rho*c1o3) + gamma;
					   gEq = gamma;

					   g[dir] = g[dir] - collFactorM*(g[dir]-gEq) + 0.5*forcingTerm[dir];

				   } 
				   else
				   {
					   //gEq = p1*WEIGTH[dir]/(rho*c1o3) + gamma;
					   gEq = pr1[ZERO]*WEIGTH1[dir] + gamma;

					   g[dir] = g[dir] - collFactorM*(g[dir]-gEq) + 0.5*forcingTerm[dir];
				   }
			   }



			   (*this->localDistributionsF)(D3Q27System::ET_E,   x1,  x2,  x3) = g[D3Q27System::INV_E];
			   (*this->localDistributionsF)(D3Q27System::ET_N,   x1,  x2,  x3) = g[D3Q27System::INV_N];
			   (*this->localDistributionsF)(D3Q27System::ET_T,   x1,  x2,  x3) = g[D3Q27System::INV_T];
			   (*this->localDistributionsF)(D3Q27System::ET_NE,  x1,  x2,  x3) = g[D3Q27System::INV_NE];
			   (*this->localDistributionsF)(D3Q27System::ET_NW,  x1p, x2,  x3) = g[D3Q27System::INV_NW];
			   (*this->localDistributionsF)(D3Q27System::ET_TE,  x1,  x2,  x3) = g[D3Q27System::INV_TE];
			   (*this->localDistributionsF)(D3Q27System::ET_TW,  x1p, x2,  x3) = g[D3Q27System::INV_TW];
			   (*this->localDistributionsF)(D3Q27System::ET_TN,  x1,  x2,  x3) = g[D3Q27System::INV_TN];
			   (*this->localDistributionsF)(D3Q27System::ET_TS,  x1,  x2p, x3) = g[D3Q27System::INV_TS];
			   (*this->localDistributionsF)(D3Q27System::ET_TNE, x1,  x2,  x3) = g[D3Q27System::INV_TNE];
			   (*this->localDistributionsF)(D3Q27System::ET_TNW, x1p, x2,  x3) = g[D3Q27System::INV_TNW];
			   (*this->localDistributionsF)(D3Q27System::ET_TSE, x1,  x2p, x3) = g[D3Q27System::INV_TSE];
			   (*this->localDistributionsF)(D3Q27System::ET_TSW, x1p, x2p, x3) = g[D3Q27System::INV_TSW];

			   (*this->nonLocalDistributionsF)(D3Q27System::ET_W,   x1p, x2,  x3 ) = g[D3Q27System::INV_W ];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_S,   x1,  x2p, x3 ) = g[D3Q27System::INV_S ];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_B,   x1,  x2,  x3p) = g[D3Q27System::INV_B ];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SW,  x1p, x2p, x3 ) = g[D3Q27System::INV_SW];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_SE,  x1,  x2p, x3 ) = g[D3Q27System::INV_SE];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BW,  x1p, x2,  x3p) = g[D3Q27System::INV_BW];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BE,  x1,  x2,  x3p) = g[D3Q27System::INV_BE];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BS,  x1,  x2p, x3p) = g[D3Q27System::INV_BS];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BN,  x1,  x2,  x3p) = g[D3Q27System::INV_BN];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSW, x1p, x2p, x3p) = g[D3Q27System::INV_BSW];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BSE, x1,  x2p, x3p) = g[D3Q27System::INV_BSE];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNW, x1p, x2,  x3p) = g[D3Q27System::INV_BNW];
			   (*this->nonLocalDistributionsF)(D3Q27System::ET_BNE, x1,  x2,  x3p) = g[D3Q27System::INV_BNE];

			   (*this->zeroDistributionsF)(x1,x2,x3) = g[D3Q27System::ZERO];

			   







			   
/////////////////////  P H A S E - F I E L D   S O L V E R /////////////////////////////////////////			   

			   
			   
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
					   //LBMReal hSource = (tauH - 0.5)*(1.0 - phi[ZERO])*(phi[ZERO])*(dirGrad_phi)/denom; // + phi[ZERO]*(dxux + dyuy + dzuz);
						   
					   //LBMReal hSource =((phi[ZERO]>phiH || phi[ZERO]<phiL) ? 0.1 : 1.0) * 3.0*mob*(-4.0)/di*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)*(dirGrad_phi)/denom;
					   LBMReal hSource = 3.0*mob*(-4.0)/di*(phi[ZERO] - phiL)*(phi[ZERO] - phiH)*(dirGrad_phi)/denom;
					   hEq = phi[ZERO]*WEIGTH[dir]*(1.0 + 3.0*velProd + 4.5*velSq1 - 1.5*(ux2+uy2+uz2)) + hSource*WEIGTH[dir];
					   //gEq = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx2+vy2+vz2))*c1o3 + (p1-rho*c1o3)*WEIGTH[dir];
					   //h[dir] = hEq; //h[dir] - (h[dir] - hEq)/(tauH + 0.5));  /// This corresponds with the collision factor of 1.0 which equals (tauH + 0.5). 
					   h[dir] = h[dir] - (h[dir] - hEq)/(tauH1); // + WEIGTH[dir]*phi[ZERO]*(dxux + dyuy + dzuz);
					   //h[dir] = h[dir] - (h[dir] - hEq)/(tauH1);
					   //g[dir] = g[dir] - collFactorM*(g[dir]-gEq) + 0.5*forcingTerm[dir];

				   } 
				   else
				   {
					   hEq = phi[ZERO]*WEIGTH[ZERO]*(1.0 - 1.5*(ux2+uy2+uz2));
					   //gEq = rho*WEIGTH[dir]*(1 + 3*velProd + 4.5*velSq1 - 1.5*(vx2+vy2+vz2))*c1o3 + (p1-rho*c1o3)*WEIGTH[dir];
					   //h[dir] = hEq;
					   h[ZERO] = h[ZERO] - (h[ZERO] - hEq)/(tauH1); // + WEIGTH[ZERO]*phi[ZERO]*(dxux + dyuy + dzuz);
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



LBMReal MultiphaseCumulantLBMKernel::gradX1_pr1()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX1[k]*pr1[k];
	}
	return 3.0*sum;
	//return 0.5*(pr1[E] - pr1[W]);
}

LBMReal MultiphaseCumulantLBMKernel::gradX2_pr1()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX2[k]*pr1[k];
	}
	return 3.0*sum;
	//return 0.5*(pr1[N] - pr1[S]);
}

LBMReal MultiphaseCumulantLBMKernel::gradX3_pr1()
{
	using namespace D3Q27System;
	LBMReal sum = 0.0;
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		sum += WEIGTH[k]*DX3[k]*pr1[k];
	}
	return 3.0*sum;
	//return 0.5*(pr1[T] - pr1[B]);
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

void MultiphaseCumulantLBMKernel::findNeighbors(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr ph, CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf, int x1, int x2, int x3)
{
	using namespace D3Q27System;
	
	BCArray3DPtr bcArray = this->getBCProcessor()->getBCArray();

	phi[ZERO] = (*ph)(x1,x2,x3);
	pr1[ZERO] = (*pf)(x1,x2,x3);

	LBMReal a = -0.5*sqrt(2*beta/kappa)*cos(contactAngle*PI/180);
	LBMReal a1 = 1 + a;
	
	for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
	{
		
		if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]))
		{
			phi[k] = (*ph)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
			pr1[k] = (*pf)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]);
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
			pr1[k] = (*pf)(x1, x2, x3);

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

void MultiphaseCumulantLBMKernel::pressureFiltering(CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf /*Pressure-Field*/, CbArray3D<LBMReal,IndexerX3X2X1>::CbArray3DPtr pf_filtered /*Pressure-Field*/)
{
	using namespace D3Q27System;
	
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

	for(int x3 = 0; x3 <= maxX3; x3++)
	{
		for(int x2 = 0; x2 <= maxX2; x2++)
		{
			for(int x1 = 0; x1 <= maxX1; x1++)
			{
				if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
				{
					int cnum = 0;
					LBMReal sum = 0.0;
					for (int k = FSTARTDIR ; k <= FENDDIR ; k++)
					{
						//if (!bcArray->isSolid(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k]))
						{
							cnum++;
							sum += (*pf)(x1 + DX1[k], x2 + DX2[k], x3 + DX3[k])*WEIGTH[k];
							 
						}
					}
					LBMReal av = sum/cnum;
					(*pf_filtered)(x1, x2, x3) = ((*pf)(x1, x2, x3))*WEIGTH[ZERO] + sum;

				}
			}
		}
	}


	
	
}



void MultiphaseCumulantLBMKernel::swapDistributions()
{
   dataSet->getFdistributions()->swap();
   dataSet->getHdistributions()->swap();
   //computePhasefield();
}