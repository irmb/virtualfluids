#include "LBMKernelETD3Q27CCLB_Geier.h"
#include "D3Q27System.h"
#include "D3Q27NoSlipBCAdapter.h"
#include "D3Q27DensityBCAdapter.h"
#include "D3Q27VelocityBCAdapter.h"
#include "SimulationParameters.h"
#include "D3Q27InterpolationProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <math.h>

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27CCLB_Geier::LBMKernelETD3Q27CCLB_Geier()
{

}
//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27CCLB_Geier::LBMKernelETD3Q27CCLB_Geier(int nx1, int nx2, int nx3, int option) 
   : LBMKernelETD3Q27(nx1, nx2, nx3),
     option(option)
{
   this->compressible = false;
}
//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27CCLB_Geier::~LBMKernelETD3Q27CCLB_Geier(void)
{

}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27CCLB_Geier::init()
{
   //DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+ghostLayerWitdh*2, nx2+ghostLayerWitdh*2, nx3+ghostLayerWitdh*2, -999.0));
   DistributionArray3DPtr d(new D3Q27EsoTwist3DSplittedVector(nx1+2, nx2+2, nx3+2, -999.0));
   dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
LBMKernel3DPtr LBMKernelETD3Q27CCLB_Geier::clone()
{
   LBMKernel3DPtr kernel(new LBMKernelETD3Q27CCLB_Geier(nx1, nx2, nx3, option));
   boost::dynamic_pointer_cast<LBMKernelETD3Q27CCLB_Geier>(kernel)->init();
   kernel->setCollisionFactor(this->collFactor);
   kernel->setBCProcessor(bcProcessor->clone(kernel));
   kernel->setWithForcing(withForcing);
   kernel->setForcingX1(muForcingX1);
   kernel->setForcingX2(muForcingX2);
   kernel->setForcingX3(muForcingX3);
   kernel->setIndex(ix1, ix2, ix3);
   kernel->setDeltaT(deltaT);
   if (option == 0)
   {
      boost::dynamic_pointer_cast<LBMKernelETD3Q27CCLB_Geier>(kernel)->OxyyMxzz = 1.0;
   }
   else if (option == 1)
   {
      boost::dynamic_pointer_cast<LBMKernelETD3Q27CCLB_Geier>(kernel)->OxyyMxzz = 2.0 +(-collFactor);
   }
   
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27CCLB_Geier::calculate()
{
   timer.resetAndStart();
   collideAll2();
   timer.stop();
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27CCLB_Geier::collideAll()
{
   using namespace D3Q27System;

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);

      muDeltaT = deltaT;

      muForcingX1.DefineVar("dx",&muDeltaT);
      muForcingX2.DefineVar("dx",&muDeltaT);
      muForcingX3.DefineVar("dx",&muDeltaT);

      muNue = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

      muForcingX1.DefineVar("nue",&muNue);
      muForcingX2.DefineVar("nue",&muNue);
      muForcingX3.DefineVar("nue",&muNue);

      LBMReal forcingX1 = 0;
      LBMReal forcingX2 = 0;
      LBMReal forcingX3 = 0;
   }
   /////////////////////////////////////

   s9 = - collFactor;
   c1o27=1.0/27.0;
   c2o3=2.0/3.0;
   w2=-1.0; //MXXpMYYpMZZ bulk viscosity
   w7=-1.0;//s9; //ORDER 4 Isotropic
   w9=-1.0;
   w10=-1.0;//s9;//-1.0; // ORDER 6 Isotropic
   w1=s9;
   // wenn es mal an den Ecken nicht gut aussieht -2.0-s9 probieren
   w3=-1.0;//-2.0-s9;//-1.0;//MXXYpMYZZ
   w4=-1.0;//-2.0-s9;//-1.0;//MXXYmMYZZ
   w5=-1.0;//-2.0-s9;//-1.0;//MYXZ
   w6=-1.0; //MXXYYpm2p
   w8=-1.0; //M_zXXYZ 


   localDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(this->getBCProcessor())->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray.getNX1();
   const int bcArrayMaxX2 = (int)bcArray.getNX2();
   const int bcArrayMaxX3 = (int)bcArray.getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

   for(int x3 = minX3; x3 < maxX3; x3++)
   {
      for(int x2 = minX2; x2 < maxX2; x2++)
      {
         for(int x1 = minX1; x1 < maxX1; x1++)
         {
            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3))
            {
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////

               //E   N  T
               //c   c  c
               //////////
               //W   S  B
               //a   a  a
               
               //Rest ist b

               LBMReal mfaaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1); 
               LBMReal mfaab = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
               LBMReal mfaac = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
               LBMReal mfaba = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
               LBMReal mfabb = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
               LBMReal mfabc = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
               LBMReal mfbaa = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
               LBMReal mfbab = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
               LBMReal mfbac = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
               LBMReal mfbba = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
               LBMReal mfbbb = (*this->zeroDistributions)(x1,x2,x3);
               LBMReal mfbbc = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               LBMReal mfaca = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
               LBMReal mfacb = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
               LBMReal mfacc = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
               LBMReal mfcaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
               LBMReal mfcab = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
               LBMReal mfcac = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
               LBMReal mfcca = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);
               LBMReal mfccb = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               LBMReal mfccc = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               LBMReal mfbca = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
               LBMReal mfbcb = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3); 
               LBMReal mfbcc = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               LBMReal mfcba = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
               LBMReal mfcbb = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               LBMReal mfcbc = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               LBMReal m0, m1, m2;

               LBMReal rho=(mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;

               mfaaa -=c1o216;
               mfaab -=c1o54;
               mfaac -=c1o216;
               mfaba -=c1o54;
               mfabb -=c2o27;
               mfabc -=c1o54;
               mfbaa -=c1o54;
               mfbab -=c2o27;
               mfbac -=c1o54;
               mfbba -=c2o27;
               mfbbb -=c8o27;
               mfbbc -=c2o27;
               mfaca -=c1o216;
               mfacb -=c1o54;
               mfacc -=c1o216;
               mfcaa -=c1o216;
               mfcab -=c1o54;
               mfcac -=c1o216;
               mfcca -=c1o216;
               mfccb -=c1o54;
               mfccc -=c1o216;
               mfbca -=c1o54;
               mfbcb -=c2o27;
               mfbcc -=c1o54;
               mfcba -=c1o54;
               mfcbb -=c2o27;
               mfcbc -=c1o54;

               LBMReal vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
                  (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                  (mfcbb-mfabb))/rho;

               LBMReal vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
                  (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                  (mfbcb-mfbab))/rho;

               LBMReal vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
                  (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                  (mfbbc-mfbba))/rho; 

               LBMReal vx2=vvx*vvx;
               LBMReal vy2=vvy*vvy;
               LBMReal vz2=vvz*vvz;



               ////////////////////////////////////////////////////////////////////////////////////
               //Hin
               ////////////////////////////////////////////////////////////////////////////////////
               // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m2    = mfaaa	+ mfaac;
               m1    = mfaac	- mfaaa;
               m0    = m2		+ mfaab;
               mfaaa = m0;
               m0   += c1o36;	
               mfaab = m1 -		m0 * vvz;
               mfaac = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfabc;
               m1    = mfabc  - mfaba;
               m0    = m2		+ mfabb;
               mfaba = m0;
               m0   += c1o9;
               mfabb = m1 -		m0 * vvz;
               mfabc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfacc;
               m1    = mfacc  - mfaca;
               m0    = m2		+ mfacb;
               mfaca = m0;
               m0   += c1o36;
               mfacb = m1 -		m0 * vvz;
               mfacc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa	+ mfbac;
               m1    = mfbac	- mfbaa;
               m0    = m2		+ mfbab;
               mfbaa = m0;
               m0   += c1o9;
               mfbab = m1 -		m0 * vvz;
               mfbac = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbba  + mfbbc;
               m1    = mfbbc  - mfbba;
               m0    = m2		+ mfbbb;
               mfbba = m0;
               m0   += c4o9;
               mfbbb = m1 -		m0 * vvz;
               mfbbc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbca  + mfbcc;
               m1    = mfbcc  - mfbca;
               m0    = m2		+ mfbcb;
               mfbca = m0;
               m0   += c1o9;
               mfbcb = m1 -		m0 * vvz;
               mfbcc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa	+ mfcac;
               m1    = mfcac	- mfcaa;
               m0    = m2		+ mfcab;
               mfcaa = m0;
               m0   += c1o36;
               mfcab = m1 -		m0 * vvz;
               mfcac = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcba  + mfcbc;
               m1    = mfcbc  - mfcba;
               m0    = m2		+ mfcbb;
               mfcba = m0;
               m0   += c1o9;
               mfcbb = m1 -		m0 * vvz;
               mfcbc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcca  + mfccc;
               m1    = mfccc  - mfcca;
               m0    = m2		+ mfccb;
               mfcca = m0;
               m0   += c1o36;
               mfccb = m1 -		m0 * vvz;
               mfccc = m2 - 2. *	m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m2    = mfaaa	+ mfaca;
               m1    = mfaca	- mfaaa;
               m0    = m2		+ mfaba;
               mfaaa = m0;
               m0   += c1o6;
               mfaba = m1 -		m0 * vvy;
               mfaca = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab  + mfacb;
               m1    = mfacb  - mfaab;
               m0    = m2		+ mfabb;
               mfaab = m0;
               mfabb = m1 -		m0 * vvy;
               mfacb = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac  + mfacc;
               m1    = mfacc  - mfaac;
               m0    = m2		+ mfabc;
               mfaac = m0;
               m0   += c1o18;
               mfabc = m1 -		m0 * vvy;
               mfacc = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa	+ mfbca;
               m1    = mfbca	- mfbaa;
               m0    = m2		+ mfbba;
               mfbaa = m0;
               m0   += c2o3;
               mfbba = m1 -		m0 * vvy;
               mfbca = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbab  + mfbcb;
               m1    = mfbcb  - mfbab;
               m0    = m2		+ mfbbb;
               mfbab = m0;
               mfbbb = m1 -		m0 * vvy;
               mfbcb = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbac  + mfbcc;
               m1    = mfbcc  - mfbac;
               m0    = m2		+ mfbbc;
               mfbac = m0;
               m0   += c2o9;
               mfbbc = m1 -		m0 * vvy;
               mfbcc = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa	+ mfcca;
               m1    = mfcca	- mfcaa;
               m0    = m2		+ mfcba;
               mfcaa = m0;
               m0   += c1o6;
               mfcba = m1 -		m0 * vvy;
               mfcca = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcab  + mfccb;
               m1    = mfccb  - mfcab;
               m0    = m2		+ mfcbb;
               mfcab = m0;
               mfcbb = m1 -		m0 * vvy;
               mfccb = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcac  + mfccc;
               m1    = mfccc  - mfcac;
               m0    = m2		+ mfcbc;
               mfcac = m0;
               m0   += c1o18;
               mfcbc = m1 -		m0 * vvy;
               mfccc = m2 - 2. *	m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m2    = mfaaa	+ mfcaa;
               m1    = mfcaa	- mfaaa;
               m0    = m2		+ mfbaa;
               mfaaa = m0;
               m0   += 1.;
               mfbaa = m1 -		m0 * vvx;
               mfcaa = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfcba;
               m1    = mfcba  - mfaba;
               m0    = m2		+ mfbba;
               mfaba = m0;
               mfbba = m1 -		m0 * vvx;
               mfcba = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfcca;
               m1    = mfcca  - mfaca;
               m0    = m2		+ mfbca;
               mfaca = m0;
               m0   += c1o3;
               mfbca = m1 -		m0 * vvx;
               mfcca = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab	+ mfcab;
               m1    = mfcab	- mfaab;
               m0    = m2		+ mfbab;
               mfaab = m0;
               mfbab = m1 -		m0 * vvx;
               mfcab = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabb  + mfcbb;
               m1    = mfcbb  - mfabb;
               m0    = m2		+ mfbbb;
               mfabb = m0;
               mfbbb = m1 -		m0 * vvx;
               mfcbb = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacb  + mfccb;
               m1    = mfccb  - mfacb;
               m0    = m2		+ mfbcb;
               mfacb = m0;
               mfbcb = m1 -		m0 * vvx;
               mfccb = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac	+ mfcac;
               m1    = mfcac	- mfaac;
               m0    = m2		+ mfbac;
               mfaac = m0;
               m0   += c1o3;
               mfbac = m1 -		m0 * vvx;
               mfcac = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabc  + mfcbc;
               m1    = mfcbc  - mfabc;
               m0    = m2		+ mfbbc;
               mfabc = m0;
               mfbbc = m1 -		m0 * vvx;
               mfcbc = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacc  + mfccc;
               m1    = mfccc  - mfacc;
               m0    = m2		+ mfbcc;
               mfacc = m0;
               m0   += c1o9;
               mfbcc = m1 -		m0 * vvx;
               mfccc = m2 - 2. *	m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////


               //////////////////////////////////////////////////////////////////////////////////////
               //// BGK
               //////////////////////////////////////////////////////////////////////////////////////
               ////2.
               //mfabb += -s9 * (-mfabb);
               //mfbab += -s9 * (-mfbab);
               //mfbba += -s9 * (-mfbba);
               //
               //mfcaa += -s9 * (c1o3 * mfaaa - mfcaa);
               //mfaca += -s9 * (c1o3 * mfaaa - mfaca);
               //mfaac += -s9 * (c1o3 * mfaaa - mfaac);
               //
               ////3.
               //mfabc += -s9 * (-mfabc);
               //mfbac += -s9 * (-mfbac);
               //
               //mfacb += -s9 * (-mfacb);
               //mfbca += -s9 * (-mfbca);

               //mfcab += -s9 * (-mfcab);
               //mfcba += -s9 * (-mfcba);

               //mfbbb += -s9 * (-mfbbb);

               ////4.
               //mfacc += -s9 * (c1o9 * mfaaa - mfacc);
               //mfcac += -s9 * (c1o9 * mfaaa - mfcac);
               //mfcca += -s9 * (c1o9 * mfaaa - mfcca);

               //mfbbc += -s9 * (-mfbbc);
               //mfbcb += -s9 * (-mfbcb);
               //mfcbb += -s9 * (-mfcbb);

               ////5.
               //mfbcc += -s9 * (-mfbcc);
               //mfcbc += -s9 * (-mfcbc);
               //mfccb += -s9 * (-mfccb);

               ////6.
               //mfccc += -s9 * (c1o27 * mfaaa - mfccc);
               ////////////////////////////////////////////////////////////////////////////////////
               //Central Moments Style
               ////////////////////////////////////////////////////////////////////////////////////
               //Cum 4.
               LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * rho) * mfabb + 2. * mfbba * mfbab) / rho;
               LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * rho) * mfbab + 2. * mfbba * mfabb) / rho;
               LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * rho) * mfbba + 2. * mfbab * mfabb) / rho; 

               LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba) / rho - c1o3 * (mfcaa + mfaca);
               LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab) / rho - c1o3 * (mfcaa + mfaac);
               LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb) / rho - c1o3 * (mfaac + mfaca);

               //Cum 5.
               LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac);
               LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc);
               LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab);

               //Cum 6.
               LBMReal CUMccc = mfccc +(-4. *  mfbbb * mfbbb  
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
                  -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
                  +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                  +  2. * (mfcaa * mfaca * mfaac)
                  + 16. *  mfbba * mfbab * mfabb) / (rho * rho)
                  - c1o3* (mfacc + mfcac + mfcca)
                  + c1o9* (mfcaa + mfaca + mfaac)
                  +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                  +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho;

               //4.
               CUMacc += 1.0 * (-CUMacc); 
               CUMcac += 1.0 * (-CUMcac); 
               CUMcca += 1.0 * (-CUMcca); 

               CUMbbc += 1.0 * (-CUMbbc); 
               CUMbcb += 1.0 * (-CUMbcb); 
               CUMcbb += 1.0 * (-CUMcbb); 

               //5.
               CUMbcc += 1.0 * (-CUMbcc);
               CUMcbc += 1.0 * (-CUMcbc);
               CUMccb += 1.0 * (-CUMccb);

               //6.
               CUMccc += 1.0 * (-CUMccc);

     
               //
               //mfabb += -s9 * (-mfabb);
               //mfbab += -s9 * (-mfbab);
               //mfbba += -s9 * (-mfbba);
               //
               //mfcaa += -s9 * (c1o3 * mfaaa - mfcaa);
               //mfaca += -s9 * (c1o3 * mfaaa - mfaca);
               //mfaac += -s9 * (c1o3 * mfaaa - mfaac);
               LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               LBMReal mxxMyy    = mfcaa - mfaca;
               LBMReal mxxMzz	   = mfcaa - mfaac;

               //relax
               mxxPyyPzz += 1.0*(mfaaa-mxxPyyPzz);
               mxxMyy    += -s9 * (-mxxMyy);
               mxxMzz    += -s9 * (-mxxMzz);
               mfabb     += -s9 * (-mfabb);
               mfbab     += -s9 * (-mfbab);
               mfbba     += -s9 * (-mfbba);

               mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);


               //3.
               //mfabc += 1.0 * (-mfabc);
               //mfbac += 1.0 * (-mfbac);
               //
               //mfacb += 1.0 * (-mfacb);
               //mfbca += 1.0 * (-mfbca);

               //mfcab += 1.0 * (-mfcab);
               //mfcba += 1.0 * (-mfcba);

               //mfbbb += 1.0 * (-mfbbb);

               //3.
               // linear combinations
               LBMReal mxxyPyzz = mfcba + mfabc;
               LBMReal mxxyMyzz = mfcba - mfabc;

               LBMReal mxxzPyyz = mfcab + mfacb;
               LBMReal mxxzMyyz = mfcab - mfacb;

               LBMReal mxyyPxzz = mfbca + mfbac;
               LBMReal mxyyMxzz = mfbca - mfbac;

               //relax


               mfbbb     +=  1.0* (-mfbbb);
               mxxyPyzz  +=  1.0 * (-mxxyPyzz);
               mxxyMyzz  += 1.0 * (-mxxyMyzz);
               mxxzPyyz  += 1.0 * (-mxxzPyyz);
               mxxzMyyz  += 1.0 * (-mxxzMyyz);
               mxyyPxzz  += 1.0 * (-mxyyPxzz);
               mxyyMxzz  += 1.0 * (-mxyyMxzz);
 
               // linear combinations back
               mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
               mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
               mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
               mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
               mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
               mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
               //4.
               mfacc += 1.0 * (c1o9 * mfaaa - mfacc);
               mfcac += 1.0 * (c1o9 * mfaaa - mfcac);
               mfcca += 1.0 * (c1o9 * mfaaa - mfcca);

               mfbbc += 1.0 * (-mfbbc);
               mfbcb += 1.0 * (-mfbcb);
               mfcbb += 1.0 * (-mfcbb);

               //5.
               mfbcc += 1.0 * (-mfbcc);
               mfcbc += 1.0 * (-mfcbc);
               mfccb += 1.0 * (-mfccb);

               ////6.
               mfccc += 1.0 * (c1o27 * mfaaa - mfccc);
 
               //back cumulants to central moments
               //4.
               mfcbb = CUMcbb + ((mfcaa + c1o3 ) * mfabb + 2. * mfbba * mfbab) / rho;
               mfbcb = CUMbcb + ((mfaca + c1o3 ) * mfbab + 2. * mfbba * mfabb) / rho;
               mfbbc = CUMbbc + ((mfaac + c1o3 ) * mfbba + 2. * mfbab * mfabb) / rho; 

               //here is problem
               mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) / rho + c1o3 * (mfcaa + mfaca)/rho-(1.0-1.0/rho)*c1o9;
               mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) / rho + c1o3 * (mfcaa + mfaac)/rho-(1.0-1.0/rho)*c1o9;
               mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) / rho + c1o3 * (mfaac + mfaca)/rho-(1.0-1.0/rho)*c1o9;

               ////5.
               mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) / rho + c1o3 * (mfbca + mfbac)/rho;
               mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) / rho + c1o3 * (mfcba + mfabc)/rho;
               mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) / rho + c1o3 * (mfacb + mfcab)/rho;

               ////6.
               mfccc = CUMccc  -((-4. *  mfbbb * mfbbb  
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
                  -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
                  +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                  +  2. * (mfcaa * mfaca * mfaac)
                  + 16. *  mfbba * mfbab * mfabb) / (rho * rho)
                  - c1o3* (mfacc + mfcac + mfcca)/rho
                  - c1o9* (mfcaa + mfaca + mfaac)/rho
                  +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                  +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho) -(1.0-1.0/rho)*c1o27;
               //////////////////////////////////////////////////////////////////////////////////////


               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               //////////////////////////////////////////////////////////////////////////////////////
               //LBMReal OxxPyyPzz = 1.;
               //LBMReal OxyyPxzz  = 1.0;
               //LBMReal OxyyMxzz  = 1.0;
               //LBMReal O4        = 1.;
               //LBMReal O5        = 1.;
               //LBMReal O6        = 1.;

               ////LBMReal OxxPyyPzz = -s9;
               ////LBMReal OxyyPxzz  = -s9;
               ////LBMReal OxyyMxzz  = -s9;
               ////LBMReal O4        = -s9;
               ////LBMReal O5        = -s9;
               ////LBMReal O6        = -s9;


               ////Cum 4.
               //LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * rho) * mfabb + 2. * mfbba * mfbab) / rho;
               //LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * rho) * mfbab + 2. * mfbba * mfabb) / rho;
               //LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * rho) * mfbba + 2. * mfbab * mfabb) / rho; 

               //LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba) / rho - c1o3 * (mfcaa + mfaca);
               //LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab) / rho - c1o3 * (mfcaa + mfaac);
               //LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb) / rho - c1o3 * (mfaac + mfaca);

               ////Cum 5.
               //LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac);
               //LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc);
               //LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab);

               ////Cum 6.
               //LBMReal CUMccc = mfccc +(-4. *  mfbbb * mfbbb  
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
               //   +  2. * (mfcaa * mfaca * mfaac)
               //   + 16. *  mfbba * mfbab * mfabb) / (rho * rho)
               //   - c1o3* (mfacc + mfcac + mfcca)
               //   + c1o9* (mfcaa + mfaca + mfaac)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho;

               //////////
               ////2.

               //// linear combinations
               ////LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               ////LBMReal mxxMyy    = mfcaa - mfaca;
               ////LBMReal mxxMzz    = mfcaa - mfaac;
               ////{
               ////   LBMReal dxux = c1o2 * (s9 *(mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
               ////   LBMReal dyuy = dxux - s9 * c3o2 * mxxMyy;
               ////   LBMReal dzuz = dxux - s9 * c3o2 * mxxMzz;
               ////   //relax
               ////   mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
               ////   mxxMyy    += -s9 * (-mxxMyy) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vy2 * dyuy);
               ////   mxxMzz    += -s9 * (-mxxMzz) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vz2 * dzuz);
               ////}
               ////mfabb     += -s9 * (-mfabb);
               ////mfbab     += -s9 * (-mfbab);
               ////mfbba     += -s9 * (-mfbba);
               ////// linear combinations back
               ////mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               ////mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               ////mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);
               ////////////

               //////// simple 2nd moments without high order Galilean invariance
               ////2.
               //// linear combinations
               //LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               //LBMReal mxxMyy    = mfcaa - mfaca;
               //LBMReal mxxMzz	   = mfcaa - mfaac;

               ////relax
               //mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);
               //mxxMyy    += -s9 * (-mxxMyy);
               //mxxMzz    += -s9 * (-mxxMzz);
               //mfabb     += -s9 * (-mfabb);
               //mfbab     += -s9 * (-mfbab);
               //mfbba     += -s9 * (-mfbba);

               //// linear combinations back
               //mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);

               ////3.
               //// linear combinations
               //LBMReal mxxyPyzz = mfcba + mfabc;
               //LBMReal mxxyMyzz = mfcba - mfabc;

               //LBMReal mxxzPyyz = mfcab + mfacb;
               //LBMReal mxxzMyyz = mfcab - mfacb;

               //LBMReal mxyyPxzz = mfbca + mfbac;
               //LBMReal mxyyMxzz = mfbca - mfbac;

               ////relax


               //mfbbb     +=  OxyyMxzz* (-mfbbb);
               //mxxyPyzz  +=  OxyyPxzz * (-mxxyPyzz);
               //mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
               //mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
               //mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
               //mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
               //mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
               ////wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
               ////mfbbb     += wadjust * (-mfbbb);
               ////wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
               ////mxxyPyzz  += wadjust * (-mxxyPyzz);
               ////wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
               ////mxxyMyzz  += wadjust * (-mxxyMyzz);
               ////wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
               ////mxxzPyyz  += wadjust * (-mxxzPyyz);
               ////wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
               ////mxxzMyyz  += wadjust * (-mxxzMyyz);
               ////wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
               ////mxyyPxzz  += wadjust * (-mxyyPxzz);
               ////wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
               ////mxyyMxzz  += wadjust * (-mxyyMxzz);

               //// linear combinations back
               //mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
               //mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
               //mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
               //mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
               //mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
               //mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

               ////4.
               //CUMacc += O4 * (-CUMacc); 
               //CUMcac += O4 * (-CUMcac); 
               //CUMcca += O4 * (-CUMcca); 

               //CUMbbc += O4 * (-CUMbbc); 
               //CUMbcb += O4 * (-CUMbcb); 
               //CUMcbb += O4 * (-CUMcbb); 

               ////5.
               //CUMbcc += O5 * (-CUMbcc);
               //CUMcbc += O5 * (-CUMcbc);
               //CUMccb += O5 * (-CUMccb);

               ////6.
               //CUMccc += O6 * (-CUMccc);

               ////back cumulants to central moments
               ////4.
               //mfcbb = CUMcbb + ((mfcaa + c1o3 * rho) * mfabb + 2. * mfbba * mfbab) / rho;
               //mfbcb = CUMbcb + ((mfaca + c1o3 * rho) * mfbab + 2. * mfbba * mfabb) / rho;
               //mfbbc = CUMbbc + ((mfaac + c1o3 * rho) * mfbba + 2. * mfbab * mfabb) / rho; 

               //mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) / rho + c1o3 * (mfcaa + mfaca);
               //mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) / rho + c1o3 * (mfcaa + mfaac);
               //mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) / rho + c1o3 * (mfaac + mfaca);

               ////5.
               //mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) / rho + c1o3 * (mfbca + mfbac);
               //mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) / rho + c1o3 * (mfcba + mfabc);
               //mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) / rho + c1o3 * (mfacb + mfcab);

               ////6.
               //mfccc = CUMccc  -((-4. *  mfbbb * mfbbb  
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
               //   +  2. * (mfcaa * mfaca * mfaac)
               //   + 16. *  mfbba * mfbab * mfabb) / (rho * rho)
               //   - c1o3* (mfacc + mfcac + mfcca)
               //   + c1o9* (mfcaa + mfaca + mfaac)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3 / rho);
               //////////////////////////////////////////////////////////////////////////////////////


               ////////////////////////////////////////////////////////////////////////////////////
               //back
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + 1.) * (     vz2 - vvz) * c1o2; 
               m1 = -mfaac        - 2. * mfaab *  vvz         +  mfaaa       * (1. - vz2)              - 1. * vz2; 
               m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + 1.) * (     vz2 + vvz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
               m1 = -mfabc        - 2. * mfabb *  vvz         + mfaba * (1. - vz2); 
               m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3) * (     vz2 - vvz) * c1o2; 
               m1 = -mfacc        - 2. * mfacb *  vvz         +  mfaca         * (1. - vz2)              - c1o3 * vz2; 
               m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3) * (     vz2 + vvz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
               m1 = -mfbac        - 2. * mfbab *  vvz         + mfbaa * (1. - vz2); 
               m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
               m1 = -mfbbc        - 2. * mfbbb *  vvz         + mfbba * (1. - vz2); 
               m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
               m1 = -mfbcc        - 2. * mfbcb *  vvz         + mfbca * (1. - vz2); 
               m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3) * (     vz2 - vvz) * c1o2; 
               m1 = -mfcac        - 2. * mfcab *  vvz         +  mfcaa         * (1. - vz2)              - c1o3 * vz2; 
               m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3) * (     vz2 + vvz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
               m1 = -mfcbc        - 2. * mfcbb *  vvz         + mfcba * (1. - vz2); 
               m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9) * (     vz2 - vvz) * c1o2; 
               m1 = -mfccc        - 2. * mfccb *  vvz         +  mfcca         * (1. - vz2)              - c1o9 * vz2; 
               m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9) * (     vz2 + vvz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6) * (     vy2 - vvy) * c1o2; 
               m1 = -mfaca        - 2. * mfaba *  vvy         +  mfaaa         * (1. - vy2)              - c1o6 * vy2; 
               m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6) * (     vy2 + vvy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3) * (     vy2 - vvy) * c1o2; 
               m1 = -mfacb        - 2. * mfabb *  vvy         +  mfaab         * (1. - vy2)              - c2o3 * vy2; 
               m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3) * (     vy2 + vvy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6) * (     vy2 - vvy) * c1o2; 
               m1 = -mfacc        - 2. * mfabc *  vvy         +  mfaac         * (1. - vy2)              - c1o6 * vy2; 
               m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6) * (     vy2 + vvy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
               m1 = -mfbca        - 2. * mfbba *  vvy         + mfbaa * (1. - vy2); 
               m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
               m1 = -mfbcb        - 2. * mfbbb *  vvy         + mfbab * (1. - vy2); 
               m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
               m1 = -mfbcc        - 2. * mfbbc *  vvy         + mfbac * (1. - vy2); 
               m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18) * (     vy2 - vvy) * c1o2; 
               m1 = -mfcca        - 2. * mfcba *  vvy         +  mfcaa          * (1. - vy2)              - c1o18 * vy2; 
               m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18) * (     vy2 + vvy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9) * (     vy2 - vvy) * c1o2; 
               m1 = -mfccb        - 2. * mfcbb *  vvy         +  mfcab         * (1. - vy2)              - c2o9 * vy2; 
               m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9) * (     vy2 + vvy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18) * (     vy2 - vvy) * c1o2; 
               m1 = -mfccc        - 2. * mfcbc *  vvy         +  mfcac          * (1. - vy2)              - c1o18 * vy2; 
               m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18) * (     vy2 + vvy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcaa        - 2. * mfbaa *  vvx         +  mfaaa          * (1. - vx2)              - c1o36 * vx2; 
               m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36) * (     vx2 + vvx) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcba        - 2. * mfbba *  vvx         +  mfaba         * (1. - vx2)              - c1o9 * vx2; 
               m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9) * (     vx2 + vvx) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcca        - 2. * mfbca *  vvx         +  mfaca          * (1. - vx2)              - c1o36 * vx2; 
               m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36) * (     vx2 + vvx) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcab        - 2. * mfbab *  vvx         +  mfaab         * (1. - vx2)              - c1o9 * vx2; 
               m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9) * (     vx2 + vvx) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcbb        - 2. * mfbbb *  vvx         +  mfabb         * (1. - vx2)              - c4o9 * vx2; 
               m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9) * (     vx2 + vvx) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9) * (     vx2 - vvx) * c1o2; 
               m1 = -mfccb        - 2. * mfbcb *  vvx         +  mfacb         * (1. - vx2)              - c1o9 * vx2; 
               m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9) * (     vx2 + vvx) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcac        - 2. * mfbac *  vvx         +  mfaac          * (1. - vx2)              - c1o36 * vx2; 
               m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36) * (     vx2 + vvx) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9) * (     vx2 - vvx) * c1o2; 
               m1 = -mfcbc        - 2. * mfbbc *  vvx         +  mfabc         * (1. - vx2)              - c1o9 * vx2; 
               m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9) * (     vx2 + vvx) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36) * (     vx2 - vvx) * c1o2; 
               m1 = -mfccc        - 2. * mfbcc *  vvx         +  mfacc          * (1. - vx2)              - c1o36 * vx2; 
               m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36) * (     vx2 + vvx) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               mfaaa +=c1o216;
               mfaab +=c1o54;
               mfaac +=c1o216;
               mfaba +=c1o54;
               mfabb +=c2o27;
               mfabc +=c1o54;
               mfbaa +=c1o54;
               mfbab +=c2o27;
               mfbac +=c1o54;
               mfbba +=c2o27;
               mfbbb +=c8o27;
               mfbbc +=c2o27;
               mfaca +=c1o216;
               mfacb +=c1o54;
               mfacc +=c1o216;
               mfcaa +=c1o216;
               mfcab +=c1o54;
               mfcac +=c1o216;
               mfcca +=c1o216;
               mfccb +=c1o54;
               mfccc +=c1o216;
               mfbca +=c1o54;
               mfbcb +=c2o27;
               mfbcc +=c1o54;
               mfcba +=c1o54;
               mfcbb +=c2o27;
               mfcbc +=c1o54;
          
               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               LBMReal rho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                     +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                     +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb; 
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = rho - rho_post;
#ifdef SINGLEPRECISION
               if(dif > 10.0E-7 || dif < -10.0E-7)
#else
               if(dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
                  UB_THROW(UbException(UB_EXARGS,"rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)));
                  //UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
                  //exit(EXIT_FAILURE);
               }
#endif
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3)    = mfabb;
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3)    = mfbab;
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3)    = mfbba;
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3)   = mfaab;
               (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3)   = mfcab;
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3)   = mfaba;
               (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3)   = mfcba;
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3)   = mfbaa;
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3)   = mfbca;
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3)  = mfaaa;
               (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3)  = mfcaa;
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3)  = mfaca;
               (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3)  = mfcca;

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = mfcbb;
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = mfbcb;
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = mfbbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = mfccb;
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = mfacb;
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = mfcbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = mfabc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = mfbcc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = mfbac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = mfccc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = mfacc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = mfcac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = mfaac;

               (*this->zeroDistributions)(x1,x2,x3) = mfbbb;
               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27CCLB_Geier::collideAll2()
{
   using namespace D3Q27System;

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);

      muDeltaT = deltaT;

      muForcingX1.DefineVar("dx",&muDeltaT);
      muForcingX2.DefineVar("dx",&muDeltaT);
      muForcingX3.DefineVar("dx",&muDeltaT);

      muNue = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

      muForcingX1.DefineVar("nue",&muNue);
      muForcingX2.DefineVar("nue",&muNue);
      muForcingX3.DefineVar("nue",&muNue);

      LBMReal forcingX1 = 0;
      LBMReal forcingX2 = 0;
      LBMReal forcingX3 = 0;
   }
   /////////////////////////////////////

   s9 = - collFactor;
   c1o27=1.0/27.0;
   c2o3=2.0/3.0;
   w2=-1.0; //MXXpMYYpMZZ bulk viscosity
   w7=-1.0;//s9; //ORDER 4 Isotropic
   w9=-1.0;
   w10=-1.0;//s9;//-1.0; // ORDER 6 Isotropic
   w1=s9;
   // wenn es mal an den Ecken nicht gut aussieht -2.0-s9 probieren
   w3=-1.0;//-2.0-s9;//-1.0;//MXXYpMYZZ
   w4=-1.0;//-2.0-s9;//-1.0;//MXXYmMYZZ
   w5=-1.0;//-2.0-s9;//-1.0;//MYXZ
   w6=-1.0; //MXXYYpm2p
   w8=-1.0; //M_zXXYZ 


   localDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(this->getBCProcessor())->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray.getNX1();
   const int bcArrayMaxX2 = (int)bcArray.getNX2();
   const int bcArrayMaxX3 = (int)bcArray.getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

   for(int x3 = minX3; x3 < maxX3; x3++)
   {
      for(int x2 = minX2; x2 < maxX2; x2++)
      {
         for(int x1 = minX1; x1 < maxX1; x1++)
         {
            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3))
            {
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////

               //E   N  T
               //c   c  c
               //////////
               //W   S  B
               //a   a  a

               //Rest ist b

               LBMReal mfaaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1); 
               LBMReal mfaab = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
               LBMReal mfaac = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
               LBMReal mfaba = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
               LBMReal mfabb = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
               LBMReal mfabc = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
               LBMReal mfbaa = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
               LBMReal mfbab = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
               LBMReal mfbac = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
               LBMReal mfbba = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
               LBMReal mfbbb = (*this->zeroDistributions)(x1,x2,x3);
               LBMReal mfbbc = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               LBMReal mfaca = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
               LBMReal mfacb = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
               LBMReal mfacc = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
               LBMReal mfcaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
               LBMReal mfcab = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
               LBMReal mfcac = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
               LBMReal mfcca = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);
               LBMReal mfccb = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               LBMReal mfccc = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               LBMReal mfbca = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
               LBMReal mfbcb = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3); 
               LBMReal mfbcc = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               LBMReal mfcba = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
               LBMReal mfcbb = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               LBMReal mfcbc = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               LBMReal m0, m1, m2;


               LBMReal rho=(mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;

               //mfaaa -=c1o216;
               //mfaab -=c1o54;
               //mfaac -=c1o216;
               //mfaba -=c1o54;
               //mfabb -=c2o27;
               //mfabc -=c1o54;
               //mfbaa -=c1o54;
               //mfbab -=c2o27;
               //mfbac -=c1o54;
               //mfbba -=c2o27;
               //mfbbb -=c8o27;
               //mfbbc -=c2o27;
               //mfaca -=c1o216;
               //mfacb -=c1o54;
               //mfacc -=c1o216;
               //mfcaa -=c1o216;
               //mfcab -=c1o54;
               //mfcac -=c1o216;
               //mfcca -=c1o216;
               //mfccb -=c1o54;
               //mfccc -=c1o216;
               //mfbca -=c1o54;
               //mfbcb -=c2o27;
               //mfbcc -=c1o54;
               //mfcba -=c1o54;
               //mfcbb -=c2o27;
               //mfcbc -=c1o54;
               
               LBMReal vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
                  (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                  (mfcbb-mfabb));
               LBMReal vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
                  (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                  (mfbcb-mfbab));
               LBMReal vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
                  (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                  (mfbbc-mfbba));

               //forcing 
               ///////////////////////////////////////////////////////////////////////////////////////////
               if (withForcing)
               {
                  muX1 = static_cast<double>(x1+ix1*bcArrayMaxX1);
                  muX2 = static_cast<double>(x2+ix2*bcArrayMaxX2);
                  muX3 = static_cast<double>(x3+ix3*bcArrayMaxX3);

                  forcingX1 = muForcingX1.Eval();
                  forcingX2 = muForcingX2.Eval();
                  forcingX3 = muForcingX3.Eval();

                  vvx += forcingX1*0.5; // X
                  vvy += forcingX2*0.5; // Y
                  vvz += forcingX3*0.5; // Z
               }
               ///////////////////////////////////////////////////////////////////////////////////////////               
               
               ////////////////////////////////////////////////////////////////////////////////////
               //fast
               //LBMReal oMdrho = 1. - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca +
               //                                   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
               //                                   mfabb+mfcbb + mfbab+mfbcb  +  mfbba+mfbbc);//fehlt mfbbb
               //LBMReal vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac +
               //                         mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
               //                         mfcbb-mfabb;
               //LBMReal vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac +
               //                         mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
               //                         mfbcb-mfbab;
               //LBMReal vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca +
               //                         mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
               //                         mfbbc-mfbba;
               ////////////////////////////////////////////////////////////////////////////////////
               // oMdrho assembler style -------> faaaaaastaaaa
               //LBMReal m0, m1, m2;
               LBMReal oMdrho;
               {
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
               }
               //LBMReal vvx;
               LBMReal vx2;
               //{
               //     vvx = mfccc-mfaaa;
               //     m0  = mfcac-mfaca;
               //     m1  = mfcaa-mfacc;
               //     m2  = mfcca-mfaac;
               //     vvx+= m0;
               //     m1 += m2;
               //     vvx+= m1;
               //     vx2 = mfcba-mfabc;
               //     m0  = mfcbc-mfaba;
               //     m1  = mfcab-mfacb;
               //     m2  = mfccb-mfaab;
               //     vx2+= m0;
               //     m1 += m2;
               //     vx2+= m1;
               //     vvx+= vx2;
               //     vx2 = mfcbb-mfabb;
               //     vvx+= vx2;
               //}
               //LBMReal vvy;
               LBMReal vy2;
               //{
               //     vvy = mfccc-mfaaa;
               //     m0  = mfaca-mfcac;
               //     m1  = mfacc-mfcaa;
               //     m2  = mfcca-mfaac;
               //     vvy+= m0;
               //     m1 += m2;
               //     vvy+= m1;
               //     vy2 = mfbca-mfbac;
               //     m0  = mfbcc-mfbaa;
               //     m1  = mfacb-mfcab;
               //     m2  = mfccb-mfaab;
               //     vy2+= m0;
               //     m1 += m2;
               //     vy2+= m1;
               //     vvy+= vy2;
               //     vy2 = mfbcb-mfbab;
               //     vvy+= vy2;
               //}
               //LBMReal vvz;
               LBMReal vz2;
               //{
               //     vvz = mfccc-mfaaa;
               //     m0  = mfcac-mfaca;
               //     m1  = mfacc-mfcaa;
               //     m2  = mfaac-mfcca;
               //     vvz+= m0;
               //     m1 += m2;
               //     vvz+= m1;
               //     vz2 = mfbac-mfbca;
               //     m0  = mfbcc-mfbaa;
               //     m1  = mfabc-mfcba;
               //     m2  = mfcbc-mfaba;
               //     vz2+= m0;
               //     m1 += m2;
               //     vz2+= m1;
               //     vvz+= vz2;
               //     vz2 = mfbbc-mfbba;
               //     vvz+= vz2;
               //}
               vx2=vvx*vvx;
               vy2=vvy*vvy;
               vz2=vvz*vvz;
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal wadjust;
               LBMReal qudricLimit = 0.01;
               //LBMReal s9 = minusomega;
               //test
               //s9 = 0.;
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
               mfaab = m1 -        m0 * vvz;
               mfaac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfabc;
               m1    = mfabc  - mfaba;
               m0    = m2          + mfabb;
               mfaba = m0;
               m0   += c1o9 * oMdrho;
               mfabb = m1 -        m0 * vvz;
               mfabc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfacc;
               m1    = mfacc  - mfaca;
               m0    = m2          + mfacb;
               mfaca = m0;
               m0   += c1o36 * oMdrho;
               mfacb = m1 -        m0 * vvz;
               mfacc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbac;
               m1    = mfbac - mfbaa;
               m0    = m2          + mfbab;
               mfbaa = m0;
               m0   += c1o9 * oMdrho;
               mfbab = m1 -        m0 * vvz;
               mfbac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbba  + mfbbc;
               m1    = mfbbc  - mfbba;
               m0    = m2          + mfbbb;
               mfbba = m0;
               m0   += c4o9 * oMdrho;
               mfbbb = m1 -        m0 * vvz;
               mfbbc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbca  + mfbcc;
               m1    = mfbcc  - mfbca;
               m0    = m2          + mfbcb;
               mfbca = m0;
               m0   += c1o9 * oMdrho;
               mfbcb = m1 -        m0 * vvz;
               mfbcc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcac;
               m1    = mfcac - mfcaa;
               m0    = m2          + mfcab;
               mfcaa = m0;
               m0   += c1o36 * oMdrho;
               mfcab = m1 -        m0 * vvz;
               mfcac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcba  + mfcbc;
               m1    = mfcbc  - mfcba;
               m0    = m2          + mfcbb;
               mfcba = m0;
               m0   += c1o9 * oMdrho;
               mfcbb = m1 -        m0 * vvz;
               mfcbc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcca  + mfccc;
               m1    = mfccc  - mfcca;
               m0    = m2          + mfccb;
               mfcca = m0;
               m0   += c1o36 * oMdrho;
               mfccb = m1 -        m0 * vvz;
               mfccc = m2 - 2. *   m1 * vvz + vz2 * m0;
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
               mfaba = m1 -        m0 * vvy;
               mfaca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab  + mfacb;
               m1    = mfacb  - mfaab;
               m0    = m2          + mfabb;
               mfaab = m0;
               mfabb = m1 -        m0 * vvy;
               mfacb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac  + mfacc;
               m1    = mfacc  - mfaac;
               m0    = m2          + mfabc;
               mfaac = m0;
               m0   += c1o18 * oMdrho;
               mfabc = m1 -        m0 * vvy;
               mfacc = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbca;
               m1    = mfbca - mfbaa;
               m0    = m2          + mfbba;
               mfbaa = m0;
               m0   += c2o3 * oMdrho;
               mfbba = m1 -        m0 * vvy;
               mfbca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbab  + mfbcb;
               m1    = mfbcb  - mfbab;
               m0    = m2          + mfbbb;
               mfbab = m0;
               mfbbb = m1 -        m0 * vvy;
               mfbcb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbac  + mfbcc;
               m1    = mfbcc  - mfbac;
               m0    = m2          + mfbbc;
               mfbac = m0;
               m0   += c2o9 * oMdrho;
               mfbbc = m1 -        m0 * vvy;
               mfbcc = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcca;
               m1    = mfcca - mfcaa;
               m0    = m2          + mfcba;
               mfcaa = m0;
               m0   += c1o6 * oMdrho;
               mfcba = m1 -        m0 * vvy;
               mfcca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcab  + mfccb;
               m1    = mfccb  - mfcab;
               m0    = m2          + mfcbb;
               mfcab = m0;
               mfcbb = m1 -        m0 * vvy;
               mfccb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcac  + mfccc;
               m1    = mfccc  - mfcac;
               m0    = m2          + mfcbc;
               mfcac = m0;
               m0   += c1o18 * oMdrho;
               mfcbc = m1 -        m0 * vvy;
               mfccc = m2 - 2. *   m1 * vvy + vy2 * m0;
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
               mfbaa = m1 -        m0 * vvx;
               mfcaa = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfcba;
               m1    = mfcba  - mfaba;
               m0    = m2          + mfbba;
               mfaba = m0;
               mfbba = m1 -        m0 * vvx;
               mfcba = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfcca;
               m1    = mfcca  - mfaca;
               m0    = m2          + mfbca;
               mfaca = m0;
               m0   += c1o3 * oMdrho;
               mfbca = m1 -        m0 * vvx;
               mfcca = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab + mfcab;
               m1    = mfcab - mfaab;
               m0    = m2          + mfbab;
               mfaab = m0;
               mfbab = m1 -        m0 * vvx;
               mfcab = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabb  + mfcbb;
               m1    = mfcbb  - mfabb;
               m0    = m2          + mfbbb;
               mfabb = m0;
               mfbbb = m1 -        m0 * vvx;
               mfcbb = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacb  + mfccb;
               m1    = mfccb  - mfacb;
               m0    = m2          + mfbcb;
               mfacb = m0;
               mfbcb = m1 -        m0 * vvx;
               mfccb = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac + mfcac;
               m1    = mfcac - mfaac;
               m0    = m2          + mfbac;
               mfaac = m0;
               m0   += c1o3 * oMdrho;
               mfbac = m1 -        m0 * vvx;
               mfcac = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabc  + mfcbc;
               m1    = mfcbc  - mfabc;
               m0    = m2          + mfbbc;
               mfabc = m0;
               mfbbc = m1 -        m0 * vvx;
               mfcbc = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacc  + mfccc;
               m1    = mfccc  - mfacc;
               m0    = m2          + mfbcc;
               mfacc = m0;
               m0   += c1o9 * oMdrho;
               mfbcc = m1 -        m0 * vvx;
               mfccc = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////////////////
               //// BGK
               //////////////////////////////////////////////////////////////////////////////////////
               ////2.
               //mfabb += -s9 * (-mfabb);
               //mfbab += -s9 * (-mfbab);
               //mfbba += -s9 * (-mfbba);
               //
               //mfcaa += -s9 * (c1o3 * mfaaa - mfcaa);
               //mfaca += -s9 * (c1o3 * mfaaa - mfaca);
               //mfaac += -s9 * (c1o3 * mfaaa - mfaac);
               //
               ////3.
               //mfabc += -s9 * (-mfabc);
               //mfbac += -s9 * (-mfbac);
               //
               //mfacb += -s9 * (-mfacb);
               //mfbca += -s9 * (-mfbca);
               //mfcab += -s9 * (-mfcab);
               //mfcba += -s9 * (-mfcba);
               //mfbbb += -s9 * (-mfbbb);
               ////4.
               //mfacc += -s9 * (c1o9 * mfaaa - mfacc);
               //mfcac += -s9 * (c1o9 * mfaaa - mfcac);
               //mfcca += -s9 * (c1o9 * mfaaa - mfcca);
               //mfbbc += -s9 * (-mfbbc);
               //mfbcb += -s9 * (-mfbcb);
               //mfcbb += -s9 * (-mfcbb);
               ////5.
               //mfbcc += -s9 * (-mfbcc);
               //mfcbc += -s9 * (-mfcbc);
               //mfccb += -s9 * (-mfccb);
               ////6.
               //mfccc += -s9 * (c1o27 * mfaaa - mfccc);
               //////////////////////////////////////////////////////////////////////////////////////

               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal OxxPyyPzz = 1.;
               LBMReal OxyyPxzz  = 1.;//-s9;//2+s9;//
               //LBMReal OxyyMxzz  = 1.;//2+s9;//
               LBMReal O4        = 1.;
               LBMReal O5        = 1.;
               LBMReal O6        = 1.;

               //Cum 4.
               LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab);
               LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb);
               LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb);

               LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

               //Cum 5.
               LBMReal CUMbcc = (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
               LBMReal CUMcbc = (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
               LBMReal CUMccb = (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

               //Cum 6.
               LBMReal CUMccc = mfccc  +((-4. *  mfbbb * mfbbb 
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
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

               {
                  LBMReal dxux = c1o2 * (s9 *(mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
                  LBMReal dyuy = dxux - s9 * c3o2 * mxxMyy;
                  LBMReal dzuz = dxux - s9 * c3o2 * mxxMzz;

                  //relax
                  mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
                  mxxMyy    += -s9 * (-mxxMyy) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vy2 * dyuy);
                  mxxMzz    += -s9 * (-mxxMzz) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vz2 * dzuz);
               }
               mfabb     += -s9 * (-mfabb);
               mfbab     += -s9 * (-mfbab);
               mfbba     += -s9 * (-mfbba);

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
               mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab);
               mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb);
               mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb);

               mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

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
               //// Cumulants
               //////////////////////////////////////////////////////////////////////////////////////
               //LBMReal OxxPyyPzz = 1.;
               //LBMReal OxyyPxzz  = 2+s9;//1.;
               //LBMReal OxyyMxzz  = 2+s9;//1.;
               //LBMReal O4        = 1.;
               //LBMReal O5        = 1.;
               //LBMReal O6        = 1.;
               ////Cum 4.
               //LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
               //LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
               //LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
               //LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba)- c1o3 * (mfcaa + mfaca);
               //LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab)- c1o3 * (mfcaa + mfaac);
               //LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb)- c1o3 * (mfaac + mfaca);
               ////Cum 5.
               //LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5)
               //   - c1o3 * (mfbca + mfbac); //O(eps^3)
               //LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
               //   - c1o3 * (mfcba + mfabc); //O(eps^3)
               //LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
               //   - c1o3 * (mfacb + mfcab);//O(eps^3)
               ////Cum 6.
               //LBMReal CUMccc = mfccc +(-4. *  mfbbb * mfbbb  //O(eps^6)
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) // O(eps^4)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) // O(eps^6)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) // O(esp^6)
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) //O(eps^6)
               //   +  2. * (mfcaa * mfaca * mfaac) //O(eps^6)
               //   + 16. *  mfbba * mfbab * mfabb) //O(eps^6)
               //   - c1o3* (mfacc + mfcac + mfcca) //O(eps^2)
               //   + c1o9* (mfcaa + mfaca + mfaac) //O(eps^2)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)//O(eps^4)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3;//O(eps^4)
               ////2.
               //// linear combinations
               //LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               //LBMReal mxxMyy    = mfcaa - mfaca;
               //LBMReal mxxMzz         = mfcaa - mfaac;
               //{
               //   LBMReal dxux = c1o2 * (s9 *(mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
               //   LBMReal dyuy = dxux - s9 * c3o2 * mxxMyy;
               //   LBMReal dzuz = dxux - s9 * c3o2 * mxxMzz;
               //   //relax
               //   mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
               //   mxxMyy    += -s9 * (-mxxMyy) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vy2 * dyuy);
               //   mxxMzz    += -s9 * (-mxxMzz) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vz2 * dzuz);
               //}
               //mfabb     += -s9 * (-mfabb);
               //mfbab     += -s9 * (-mfbab);
               //mfbba     += -s9 * (-mfbba);
               //// linear combinations back
               //mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);
               ////3.
               //// linear combinations
               //LBMReal mxxyPyzz = mfcba + mfabc;
               //LBMReal mxxyMyzz = mfcba - mfabc;
               //LBMReal mxxzPyyz = mfcab + mfacb;
               //LBMReal mxxzMyyz = mfcab - mfacb;
               //LBMReal mxyyPxzz = mfbca + mfbac;
               //LBMReal mxyyMxzz = mfbca - mfbac;
               ////relax
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
               //mfbbb     += wadjust * (-mfbbb);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
               //mxxyPyzz  += wadjust * (-mxxyPyzz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
               //mxxyMyzz  += wadjust * (-mxxyMyzz);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
               //mxxzPyyz  += wadjust * (-mxxzPyyz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
               //mxxzMyyz  += wadjust * (-mxxzMyyz);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
               //mxyyPxzz  += wadjust * (-mxyyPxzz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
               //mxyyMxzz  += wadjust * (-mxyyMxzz);
               //// linear combinations back
               //mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
               //mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
               //mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
               //mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
               //mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
               //mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
               ////4.
               //CUMacc += O4 * (-CUMacc);
               //CUMcac += O4 * (-CUMcac);
               //CUMcca += O4 * (-CUMcca);
               //CUMbbc += O4 * (-CUMbbc);
               //CUMbcb += O4 * (-CUMbcb);
               //CUMcbb += O4 * (-CUMcbb);
               ////5.
               //CUMbcc += O5 * (-CUMbcc);
               //CUMcbc += O5 * (-CUMcbc);
               //CUMccb += O5 * (-CUMccb);
               ////6.
               //CUMccc += O6 * (-CUMccc);
               ////back cumulants to central moments
               ////4.
               //mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
               //mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
               //mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
               //mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca);
               //mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac);
               //mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca);
               ////5.
               //mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac);
               //mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc);
               //mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab);
               ////6.
               //mfccc = CUMccc  -((-4. *  mfbbb * mfbbb 
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
               //   +  2. * (mfcaa * mfaca * mfaac)
               //   + 16. *  mfbba * mfbab * mfabb)
               //   - c1o3* (mfacc + mfcac + mfcca)
               //   + c1o9* (mfcaa + mfaca + mfaac)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3);
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //back
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + 1. * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfaac        - 2. * mfaab *  vvz         +  mfaaa                * (1. - vz2)              - 1. * oMdrho * vz2;
               m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + 1. * oMdrho) * (     vz2 + vvz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2;
               m1 = -mfabc        - 2. * mfabb *  vvz         + mfaba * (1. - vz2);
               m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfacc        - 2. * mfacb *  vvz         +  mfaca                  * (1. - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2;
               m1 = -mfbac        - 2. * mfbab *  vvz         + mfbaa * (1. - vz2);
               m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2;
               m1 = -mfbbc        - 2. * mfbbb *  vvz         + mfbba * (1. - vz2);
               m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2;
               m1 = -mfbcc        - 2. * mfbcb *  vvz         + mfbca * (1. - vz2);
               m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfcac        - 2. * mfcab *  vvz         +  mfcaa                  * (1. - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2;
               m1 = -mfcbc        - 2. * mfcbb *  vvz         + mfcba * (1. - vz2);
               m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfccc        - 2. * mfccb *  vvz         +  mfcca                  * (1. - vz2)              - c1o9 * oMdrho * vz2;
               m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfaca        - 2. * mfaba *  vvy         +  mfaaa                  * (1. - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfacb        - 2. * mfabb *  vvy         +  mfaab                  * (1. - vy2)              - c2o3 * oMdrho * vy2;
               m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfacc        - 2. * mfabc *  vvy         +  mfaac                  * (1. - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2;
               m1 = -mfbca        - 2. * mfbba *  vvy         + mfbaa * (1. - vy2);
               m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2;
               m1 = -mfbcb        - 2. * mfbbb *  vvy         + mfbab * (1. - vy2);
               m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2;
               m1 = -mfbcc        - 2. * mfbbc *  vvy         + mfbac * (1. - vy2);
               m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfcca        - 2. * mfcba *  vvy         +  mfcaa                   * (1. - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfccb        - 2. * mfcbb *  vvy         +  mfcab                  * (1. - vy2)              - c2o9 * oMdrho * vy2;
               m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfccc        - 2. * mfcbc *  vvy         +  mfcac                   * (1. - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcaa        - 2. * mfbaa *  vvx         +  mfaaa                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcba        - 2. * mfbba *  vvx         +  mfaba                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcca        - 2. * mfbca *  vvx         +  mfaca                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcab        - 2. * mfbab *  vvx         +  mfaab                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcbb        - 2. * mfbbb *  vvx         +  mfabb                  * (1. - vx2)              - c4o9 * oMdrho * vx2;
               m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfccb        - 2. * mfbcb *  vvx         +  mfacb                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcac        - 2. * mfbac *  vvx         +  mfaac                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcbc        - 2. * mfbbc *  vvx         +  mfabc                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfccc        - 2. * mfbcc *  vvx         +  mfacc                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;

               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mfaaa +=c1o216;
               //mfaab +=c1o54;
               //mfaac +=c1o216;
               //mfaba +=c1o54;
               //mfabb +=c2o27;
               //mfabc +=c1o54;
               //mfbaa +=c1o54;
               //mfbab +=c2o27;
               //mfbac +=c1o54;
               //mfbba +=c2o27;
               //mfbbb +=c8o27;
               //mfbbc +=c2o27;
               //mfaca +=c1o216;
               //mfacb +=c1o54;
               //mfacc +=c1o216;
               //mfcaa +=c1o216;
               //mfcab +=c1o54;
               //mfcac +=c1o216;
               //mfcca +=c1o216;
               //mfccb +=c1o54;
               //mfccc +=c1o216;
               //mfbca +=c1o54;
               //mfbcb +=c2o27;
               //mfbcc +=c1o54;
               //mfcba +=c1o54;
               //mfcbb +=c2o27;
               //mfcbc +=c1o54;

               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               LBMReal rho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb; 
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = rho - rho_post;
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
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3)    = mfabb;
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3)    = mfbab;
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3)    = mfbba;
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3)   = mfaab;
               (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3)   = mfcab;
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3)   = mfaba;
               (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3)   = mfcba;
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3)   = mfbaa;
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3)   = mfbca;
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3)  = mfaaa;
               (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3)  = mfcaa;
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3)  = mfaca;
               (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3)  = mfcca;

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = mfcbb;
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = mfbcb;
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = mfbbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = mfccb;
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = mfacb;
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = mfcbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = mfabc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = mfbcc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = mfbac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = mfccc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = mfacc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = mfcac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = mfaac;

               (*this->zeroDistributions)(x1,x2,x3) = mfbbb;
               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27CCLB_Geier::collideAll3()
{
   using namespace D3Q27System;

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);

      muDeltaT = deltaT;

      muForcingX1.DefineVar("dx",&muDeltaT);
      muForcingX2.DefineVar("dx",&muDeltaT);
      muForcingX3.DefineVar("dx",&muDeltaT);

      muNue = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

      muForcingX1.DefineVar("nue",&muNue);
      muForcingX2.DefineVar("nue",&muNue);
      muForcingX3.DefineVar("nue",&muNue);

      LBMReal forcingX1 = 0;
      LBMReal forcingX2 = 0;
      LBMReal forcingX3 = 0;
   }
   /////////////////////////////////////

   s9 = - collFactor;
   c1o27=1.0/27.0;
   c2o3=2.0/3.0;
   w2=-1.0; //MXXpMYYpMZZ bulk viscosity
   w7=-1.0;//s9; //ORDER 4 Isotropic
   w9=-1.0;
   w10=-1.0;//s9;//-1.0; // ORDER 6 Isotropic
   w1=s9;
   // wenn es mal an den Ecken nicht gut aussieht -2.0-s9 probieren
   w3=-1.0;//-2.0-s9;//-1.0;//MXXYpMYZZ
   w4=-1.0;//-2.0-s9;//-1.0;//MXXYmMYZZ
   w5=-1.0;//-2.0-s9;//-1.0;//MYXZ
   w6=-1.0; //MXXYYpm2p
   w8=-1.0; //M_zXXYZ 


   localDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = boost::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(this->getBCProcessor())->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray.getNX1();
   const int bcArrayMaxX2 = (int)bcArray.getNX2();
   const int bcArrayMaxX3 = (int)bcArray.getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

   for(int x3 = minX3; x3 < maxX3; x3++)
   {
      for(int x2 = minX2; x2 < maxX2; x2++)
      {
         for(int x1 = minX1; x1 < maxX1; x1++)
         {
            if(!bcArray.isSolid(x1,x2,x3) && !bcArray.isUndefined(x1,x2,x3))
            {
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////

               //E   N  T
               //c   c  c
               //////////
               //W   S  B
               //a   a  a

               //Rest ist b

               LBMReal mfaaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1); 
               LBMReal mfaab = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3 );
               LBMReal mfaac = (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3);
               LBMReal mfaba = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,x3+1 );
               LBMReal mfabb = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,x3  );
               LBMReal mfabc = (*this->localDistributions)(D3Q27System::ET_TW, x1+1,x2,x3);
               LBMReal mfbaa = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2+1,x3+1 );
               LBMReal mfbab = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2+1,x3  );
               LBMReal mfbac = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2+1,x3);
               LBMReal mfbba = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3+1  );
               LBMReal mfbbb = (*this->zeroDistributions)(x1,x2,x3);
               LBMReal mfbbc = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               LBMReal mfaca = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,x3+1);
               LBMReal mfacb = (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,x3);
               LBMReal mfacc = (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,x3);
               LBMReal mfcaa = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2+1,x3+1);
               LBMReal mfcab = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2+1,x3 );
               LBMReal mfcac = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2+1,x3);
               LBMReal mfcca = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3+1);
               LBMReal mfccb = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               LBMReal mfccc = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               LBMReal mfbca = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3+1 );
               LBMReal mfbcb = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3); 
               LBMReal mfbcc = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               LBMReal mfcba = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3+1 );
               LBMReal mfcbb = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               LBMReal mfcbc = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               LBMReal m0, m1, m2;


               LBMReal rho=(mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;

               //mfaaa -=c1o216;
               //mfaab -=c1o54;
               //mfaac -=c1o216;
               //mfaba -=c1o54;
               //mfabb -=c2o27;
               //mfabc -=c1o54;
               //mfbaa -=c1o54;
               //mfbab -=c2o27;
               //mfbac -=c1o54;
               //mfbba -=c2o27;
               //mfbbb -=c8o27;
               //mfbbc -=c2o27;
               //mfaca -=c1o216;
               //mfacb -=c1o54;
               //mfacc -=c1o216;
               //mfcaa -=c1o216;
               //mfcab -=c1o54;
               //mfcac -=c1o216;
               //mfcca -=c1o216;
               //mfccb -=c1o54;
               //mfccc -=c1o216;
               //mfbca -=c1o54;
               //mfbcb -=c2o27;
               //mfbcc -=c1o54;
               //mfcba -=c1o54;
               //mfcbb -=c2o27;
               //mfcbc -=c1o54;

               LBMReal vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
                  (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                  (mfcbb-mfabb));
               LBMReal vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
                  (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                  (mfbcb-mfbab));
               LBMReal vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
                  (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                  (mfbbc-mfbba));
               ////////////////////////////////////////////////////////////////////////////////////
               //fast
               //LBMReal oMdrho = 1. - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca +
               //                                   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
               //                                   mfabb+mfcbb + mfbab+mfbcb  +  mfbba+mfbbc);//fehlt mfbbb
               //LBMReal vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac +
               //                         mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
               //                         mfcbb-mfabb;
               //LBMReal vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac +
               //                         mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
               //                         mfbcb-mfbab;
               //LBMReal vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca +
               //                         mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
               //                         mfbbc-mfbba;
               ////////////////////////////////////////////////////////////////////////////////////
               // oMdrho assembler style -------> faaaaaastaaaa
               //LBMReal m0, m1, m2;
               LBMReal oMdrho;
               {
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
               }
               //LBMReal vvx;
               LBMReal vx2;
               //{
               //     vvx = mfccc-mfaaa;
               //     m0  = mfcac-mfaca;
               //     m1  = mfcaa-mfacc;
               //     m2  = mfcca-mfaac;
               //     vvx+= m0;
               //     m1 += m2;
               //     vvx+= m1;
               //     vx2 = mfcba-mfabc;
               //     m0  = mfcbc-mfaba;
               //     m1  = mfcab-mfacb;
               //     m2  = mfccb-mfaab;
               //     vx2+= m0;
               //     m1 += m2;
               //     vx2+= m1;
               //     vvx+= vx2;
               //     vx2 = mfcbb-mfabb;
               //     vvx+= vx2;
               //}
               //LBMReal vvy;
               LBMReal vy2;
               //{
               //     vvy = mfccc-mfaaa;
               //     m0  = mfaca-mfcac;
               //     m1  = mfacc-mfcaa;
               //     m2  = mfcca-mfaac;
               //     vvy+= m0;
               //     m1 += m2;
               //     vvy+= m1;
               //     vy2 = mfbca-mfbac;
               //     m0  = mfbcc-mfbaa;
               //     m1  = mfacb-mfcab;
               //     m2  = mfccb-mfaab;
               //     vy2+= m0;
               //     m1 += m2;
               //     vy2+= m1;
               //     vvy+= vy2;
               //     vy2 = mfbcb-mfbab;
               //     vvy+= vy2;
               //}
               //LBMReal vvz;
               LBMReal vz2;
               //{
               //     vvz = mfccc-mfaaa;
               //     m0  = mfcac-mfaca;
               //     m1  = mfacc-mfcaa;
               //     m2  = mfaac-mfcca;
               //     vvz+= m0;
               //     m1 += m2;
               //     vvz+= m1;
               //     vz2 = mfbac-mfbca;
               //     m0  = mfbcc-mfbaa;
               //     m1  = mfabc-mfcba;
               //     m2  = mfcbc-mfaba;
               //     vz2+= m0;
               //     m1 += m2;
               //     vz2+= m1;
               //     vvz+= vz2;
               //     vz2 = mfbbc-mfbba;
               //     vvz+= vz2;
               //}
               vx2=vvx*vvx;
               vy2=vvy*vvy;
               vz2=vvz*vvz;
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal wadjust;
               LBMReal qudricLimit = 0.01;
               //LBMReal s9 = minusomega;
               //test
               //s9 = 0.;
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
               mfaab = m1 -        m0 * vvz;
               mfaac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfabc;
               m1    = mfabc  - mfaba;
               m0    = m2          + mfabb;
               mfaba = m0;
               m0   += c1o9 * oMdrho;
               mfabb = m1 -        m0 * vvz;
               mfabc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfacc;
               m1    = mfacc  - mfaca;
               m0    = m2          + mfacb;
               mfaca = m0;
               m0   += c1o36 * oMdrho;
               mfacb = m1 -        m0 * vvz;
               mfacc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbac;
               m1    = mfbac - mfbaa;
               m0    = m2          + mfbab;
               mfbaa = m0;
               m0   += c1o9 * oMdrho;
               mfbab = m1 -        m0 * vvz;
               mfbac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbba  + mfbbc;
               m1    = mfbbc  - mfbba;
               m0    = m2          + mfbbb;
               mfbba = m0;
               m0   += c4o9 * oMdrho;
               mfbbb = m1 -        m0 * vvz;
               mfbbc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbca  + mfbcc;
               m1    = mfbcc  - mfbca;
               m0    = m2          + mfbcb;
               mfbca = m0;
               m0   += c1o9 * oMdrho;
               mfbcb = m1 -        m0 * vvz;
               mfbcc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcac;
               m1    = mfcac - mfcaa;
               m0    = m2          + mfcab;
               mfcaa = m0;
               m0   += c1o36 * oMdrho;
               mfcab = m1 -        m0 * vvz;
               mfcac = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcba  + mfcbc;
               m1    = mfcbc  - mfcba;
               m0    = m2          + mfcbb;
               mfcba = m0;
               m0   += c1o9 * oMdrho;
               mfcbb = m1 -        m0 * vvz;
               mfcbc = m2 - 2. *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcca  + mfccc;
               m1    = mfccc  - mfcca;
               m0    = m2          + mfccb;
               mfcca = m0;
               m0   += c1o36 * oMdrho;
               mfccb = m1 -        m0 * vvz;
               mfccc = m2 - 2. *   m1 * vvz + vz2 * m0;
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
               mfaba = m1 -        m0 * vvy;
               mfaca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab  + mfacb;
               m1    = mfacb  - mfaab;
               m0    = m2          + mfabb;
               mfaab = m0;
               mfabb = m1 -        m0 * vvy;
               mfacb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac  + mfacc;
               m1    = mfacc  - mfaac;
               m0    = m2          + mfabc;
               mfaac = m0;
               m0   += c1o18 * oMdrho;
               mfabc = m1 -        m0 * vvy;
               mfacc = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbca;
               m1    = mfbca - mfbaa;
               m0    = m2          + mfbba;
               mfbaa = m0;
               m0   += c2o3 * oMdrho;
               mfbba = m1 -        m0 * vvy;
               mfbca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbab  + mfbcb;
               m1    = mfbcb  - mfbab;
               m0    = m2          + mfbbb;
               mfbab = m0;
               mfbbb = m1 -        m0 * vvy;
               mfbcb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbac  + mfbcc;
               m1    = mfbcc  - mfbac;
               m0    = m2          + mfbbc;
               mfbac = m0;
               m0   += c2o9 * oMdrho;
               mfbbc = m1 -        m0 * vvy;
               mfbcc = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcca;
               m1    = mfcca - mfcaa;
               m0    = m2          + mfcba;
               mfcaa = m0;
               m0   += c1o6 * oMdrho;
               mfcba = m1 -        m0 * vvy;
               mfcca = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcab  + mfccb;
               m1    = mfccb  - mfcab;
               m0    = m2          + mfcbb;
               mfcab = m0;
               mfcbb = m1 -        m0 * vvy;
               mfccb = m2 - 2. *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcac  + mfccc;
               m1    = mfccc  - mfcac;
               m0    = m2          + mfcbc;
               mfcac = m0;
               m0   += c1o18 * oMdrho;
               mfcbc = m1 -        m0 * vvy;
               mfccc = m2 - 2. *   m1 * vvy + vy2 * m0;
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
               mfbaa = m1 -        m0 * vvx;
               mfcaa = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfcba;
               m1    = mfcba  - mfaba;
               m0    = m2          + mfbba;
               mfaba = m0;
               mfbba = m1 -        m0 * vvx;
               mfcba = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfcca;
               m1    = mfcca  - mfaca;
               m0    = m2          + mfbca;
               mfaca = m0;
               m0   += c1o3 * oMdrho;
               mfbca = m1 -        m0 * vvx;
               mfcca = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab + mfcab;
               m1    = mfcab - mfaab;
               m0    = m2          + mfbab;
               mfaab = m0;
               mfbab = m1 -        m0 * vvx;
               mfcab = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabb  + mfcbb;
               m1    = mfcbb  - mfabb;
               m0    = m2          + mfbbb;
               mfabb = m0;
               mfbbb = m1 -        m0 * vvx;
               mfcbb = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacb  + mfccb;
               m1    = mfccb  - mfacb;
               m0    = m2          + mfbcb;
               mfacb = m0;
               mfbcb = m1 -        m0 * vvx;
               mfccb = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac + mfcac;
               m1    = mfcac - mfaac;
               m0    = m2          + mfbac;
               mfaac = m0;
               m0   += c1o3 * oMdrho;
               mfbac = m1 -        m0 * vvx;
               mfcac = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabc  + mfcbc;
               m1    = mfcbc  - mfabc;
               m0    = m2          + mfbbc;
               mfabc = m0;
               mfbbc = m1 -        m0 * vvx;
               mfcbc = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacc  + mfccc;
               m1    = mfccc  - mfacc;
               m0    = m2          + mfbcc;
               mfacc = m0;
               m0   += c1o9 * oMdrho;
               mfbcc = m1 -        m0 * vvx;
               mfccc = m2 - 2. *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////////////////
               //// BGK
               //////////////////////////////////////////////////////////////////////////////////////
               ////2.
               //mfabb += -s9 * (-mfabb);
               //mfbab += -s9 * (-mfbab);
               //mfbba += -s9 * (-mfbba);
               //
               //mfcaa += -s9 * (c1o3 * mfaaa - mfcaa);
               //mfaca += -s9 * (c1o3 * mfaaa - mfaca);
               //mfaac += -s9 * (c1o3 * mfaaa - mfaac);
               //
               ////3.
               //mfabc += -s9 * (-mfabc);
               //mfbac += -s9 * (-mfbac);
               //
               //mfacb += -s9 * (-mfacb);
               //mfbca += -s9 * (-mfbca);
               //mfcab += -s9 * (-mfcab);
               //mfcba += -s9 * (-mfcba);
               //mfbbb += -s9 * (-mfbbb);
               ////4.
               //mfacc += -s9 * (c1o9 * mfaaa - mfacc);
               //mfcac += -s9 * (c1o9 * mfaaa - mfcac);
               //mfcca += -s9 * (c1o9 * mfaaa - mfcca);
               //mfbbc += -s9 * (-mfbbc);
               //mfbcb += -s9 * (-mfbcb);
               //mfcbb += -s9 * (-mfcbb);
               ////5.
               //mfbcc += -s9 * (-mfbcc);
               //mfcbc += -s9 * (-mfcbc);
               //mfccb += -s9 * (-mfccb);
               ////6.
               //mfccc += -s9 * (c1o27 * mfaaa - mfccc);
               //////////////////////////////////////////////////////////////////////////////////////

               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal OxxPyyPzz = 1.;
               LBMReal OxyyPxzz  = 1.;//-s9;//2+s9;//
               LBMReal OxyyMxzz  = 2+s9;//
               LBMReal O4        = 1.;
               LBMReal O5        = 1.;
               LBMReal O6        = 1.;

               //Cum 4.
               LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab);
               LBMReal CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb);
               LBMReal CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb);

               LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

               //Cum 5.
               LBMReal CUMbcc = (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
               LBMReal CUMcbc = (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
               LBMReal CUMccb = (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

               //Cum 6.
               LBMReal CUMccc = mfccc  +((-4. *  mfbbb * mfbbb 
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
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

               {
                  LBMReal dxux = c1o2 * (s9 *(mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
                  LBMReal dyuy = dxux - s9 * c3o2 * mxxMyy;
                  LBMReal dzuz = dxux - s9 * c3o2 * mxxMzz;

                  //relax
                  mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
                  mxxMyy    += -s9 * (-mxxMyy) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vy2 * dyuy);
                  mxxMzz    += -s9 * (-mxxMzz) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vz2 * dzuz);
               }
               mfabb     += -s9 * (-mfabb);
               mfbab     += -s9 * (-mfbab);
               mfbba     += -s9 * (-mfbba);

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
               mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab);
               mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb);
               mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb);

               mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

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

               //////////////////////////////////////////////////////////////////////////////////////
               //// Cumulants
               //////////////////////////////////////////////////////////////////////////////////////
               //LBMReal OxxPyyPzz = 1.;
               //LBMReal OxyyPxzz  = 2+s9;//1.;
               //LBMReal OxyyMxzz  = 2+s9;//1.;
               //LBMReal O4        = 1.;
               //LBMReal O5        = 1.;
               //LBMReal O6        = 1.;
               ////Cum 4.
               //LBMReal CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
               //LBMReal CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
               //LBMReal CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
               //LBMReal CUMcca = mfcca - (mfcaa * mfaca + 2. * mfbba * mfbba)- c1o3 * (mfcaa + mfaca);
               //LBMReal CUMcac = mfcac - (mfcaa * mfaac + 2. * mfbab * mfbab)- c1o3 * (mfcaa + mfaac);
               //LBMReal CUMacc = mfacc - (mfaac * mfaca + 2. * mfabb * mfabb)- c1o3 * (mfaac + mfaca);
               ////Cum 5.
               //LBMReal CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) //O(eps^5)
               //   - c1o3 * (mfbca + mfbac); //O(eps^3)
               //LBMReal CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) //O(eps^5)
               //   - c1o3 * (mfcba + mfabc); //O(eps^3)
               //LBMReal CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) //O(eps^5)
               //   - c1o3 * (mfacb + mfcab);//O(eps^3)
               ////Cum 6.
               //LBMReal CUMccc = mfccc +(-4. *  mfbbb * mfbbb  //O(eps^6)
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca) // O(eps^4)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc) // O(eps^6)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) // O(esp^6)
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac) //O(eps^6)
               //   +  2. * (mfcaa * mfaca * mfaac) //O(eps^6)
               //   + 16. *  mfbba * mfbab * mfabb) //O(eps^6)
               //   - c1o3* (mfacc + mfcac + mfcca) //O(eps^2)
               //   + c1o9* (mfcaa + mfaca + mfaac) //O(eps^2)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)//O(eps^4)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3;//O(eps^4)
               ////2.
               //// linear combinations
               //LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
               //LBMReal mxxMyy    = mfcaa - mfaca;
               //LBMReal mxxMzz         = mfcaa - mfaac;
               //{
               //   LBMReal dxux = c1o2 * (s9 *(mxxMyy + mxxMzz) + (mfaaa - mxxPyyPzz));
               //   LBMReal dyuy = dxux - s9 * c3o2 * mxxMyy;
               //   LBMReal dzuz = dxux - s9 * c3o2 * mxxMzz;
               //   //relax
               //   mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- 3. * (1. - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
               //   mxxMyy    += -s9 * (-mxxMyy) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vy2 * dyuy);
               //   mxxMzz    += -s9 * (-mxxMzz) - 3. * (1. + c1o2 * s9) * (vx2 * dxux + vz2 * dzuz);
               //}
               //mfabb     += -s9 * (-mfabb);
               //mfbab     += -s9 * (-mfbab);
               //mfbba     += -s9 * (-mfbba);
               //// linear combinations back
               //mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaca = c1o3 * (-2. *  mxxMyy +      mxxMzz + mxxPyyPzz);
               //mfaac = c1o3 * (       mxxMyy - 2. * mxxMzz + mxxPyyPzz);
               ////3.
               //// linear combinations
               //LBMReal mxxyPyzz = mfcba + mfabc;
               //LBMReal mxxyMyzz = mfcba - mfabc;
               //LBMReal mxxzPyyz = mfcab + mfacb;
               //LBMReal mxxzMyyz = mfcab - mfacb;
               //LBMReal mxyyPxzz = mfbca + mfbac;
               //LBMReal mxyyMxzz = mfbca - mfbac;
               ////relax
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mfbbb)/(abs(mfbbb)+qudricLimit);
               //mfbbb     += wadjust * (-mfbbb);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimit);
               //mxxyPyzz  += wadjust * (-mxxyPyzz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimit);
               //mxxyMyzz  += wadjust * (-mxxyMyzz);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimit);
               //mxxzPyyz  += wadjust * (-mxxzPyyz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimit);
               //mxxzMyyz  += wadjust * (-mxxzMyyz);
               //wadjust    = OxyyPxzz+(1.-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimit);
               //mxyyPxzz  += wadjust * (-mxyyPxzz);
               //wadjust    = OxyyMxzz+(1.-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimit);
               //mxyyMxzz  += wadjust * (-mxyyMxzz);
               //// linear combinations back
               //mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
               //mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
               //mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
               //mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
               //mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
               //mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
               ////4.
               //CUMacc += O4 * (-CUMacc);
               //CUMcac += O4 * (-CUMcac);
               //CUMcca += O4 * (-CUMcca);
               //CUMbbc += O4 * (-CUMbbc);
               //CUMbcb += O4 * (-CUMbcb);
               //CUMcbb += O4 * (-CUMcbb);
               ////5.
               //CUMbcc += O5 * (-CUMbcc);
               //CUMcbc += O5 * (-CUMcbc);
               //CUMccb += O5 * (-CUMccb);
               ////6.
               //CUMccc += O6 * (-CUMccc);
               ////back cumulants to central moments
               ////4.
               //mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + 2. * mfbba * mfbab);
               //mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + 2. * mfbba * mfabb);
               //mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + 2. * mfbab * mfabb);
               //mfcca = CUMcca + (mfcaa * mfaca + 2. * mfbba * mfbba) + c1o3 * (mfcaa + mfaca);
               //mfcac = CUMcac + (mfcaa * mfaac + 2. * mfbab * mfbab) + c1o3 * (mfcaa + mfaac);
               //mfacc = CUMacc + (mfaac * mfaca + 2. * mfabb * mfabb) + c1o3 * (mfaac + mfaca);
               ////5.
               //mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + 4. * mfabb * mfbbb + 2. * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac);
               //mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + 4. * mfbab * mfbbb + 2. * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc);
               //mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + 4. * mfbba * mfbbb + 2. * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab);
               ////6.
               //mfccc = CUMccc  -((-4. *  mfbbb * mfbbb 
               //   -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
               //   -  4. * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
               //   -  2. * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
               //   +( 4. * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
               //   +  2. * (mfcaa * mfaca * mfaac)
               //   + 16. *  mfbba * mfbab * mfabb)
               //   - c1o3* (mfacc + mfcac + mfcca)
               //   + c1o9* (mfcaa + mfaca + mfaac)
               //   +( 2. * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
               //   +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3);
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //back
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + 1. * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfaac        - 2. * mfaab *  vvz         +  mfaaa                * (1. - vz2)              - 1. * oMdrho * vz2;
               m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + 1. * oMdrho) * (     vz2 + vvz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2;
               m1 = -mfabc        - 2. * mfabb *  vvz         + mfaba * (1. - vz2);
               m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfacc        - 2. * mfacb *  vvz         +  mfaca                  * (1. - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2;
               m1 = -mfbac        - 2. * mfbab *  vvz         + mfbaa * (1. - vz2);
               m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2;
               m1 = -mfbbc        - 2. * mfbbb *  vvz         + mfbba * (1. - vz2);
               m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2;
               m1 = -mfbcc        - 2. * mfbcb *  vvz         + mfbca * (1. - vz2);
               m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfcac        - 2. * mfcab *  vvz         +  mfcaa                  * (1. - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2;
               m1 = -mfcbc        - 2. * mfcbb *  vvz         + mfcba * (1. - vz2);
               m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2;
               m1 = -mfccc        - 2. * mfccb *  vvz         +  mfcca                  * (1. - vz2)              - c1o9 * oMdrho * vz2;
               m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfaca        - 2. * mfaba *  vvy         +  mfaaa                  * (1. - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfacb        - 2. * mfabb *  vvy         +  mfaab                  * (1. - vy2)              - c2o3 * oMdrho * vy2;
               m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfacc        - 2. * mfabc *  vvy         +  mfaac                  * (1. - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2;
               m1 = -mfbca        - 2. * mfbba *  vvy         + mfbaa * (1. - vy2);
               m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2;
               m1 = -mfbcb        - 2. * mfbbb *  vvy         + mfbab * (1. - vy2);
               m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2;
               m1 = -mfbcc        - 2. * mfbbc *  vvy         + mfbac * (1. - vy2);
               m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfcca        - 2. * mfcba *  vvy         +  mfcaa                   * (1. - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfccb        - 2. * mfcbb *  vvy         +  mfcab                  * (1. - vy2)              - c2o9 * oMdrho * vy2;
               m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2;
               m1 = -mfccc        - 2. * mfcbc *  vvy         +  mfcac                   * (1. - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcaa        - 2. * mfbaa *  vvx         +  mfaaa                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcba        - 2. * mfbba *  vvx         +  mfaba                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcca        - 2. * mfbca *  vvx         +  mfaca                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcab        - 2. * mfbab *  vvx         +  mfaab                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcbb        - 2. * mfbbb *  vvx         +  mfabb                  * (1. - vx2)              - c4o9 * oMdrho * vx2;
               m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfccb        - 2. * mfbcb *  vvx         +  mfacb                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcac        - 2. * mfbac *  vvx         +  mfaac                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfcbc        - 2. * mfbbc *  vvx         +  mfabc                  * (1. - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2;
               m1 = -mfccc        - 2. * mfbcc *  vvx         +  mfacc                   * (1. - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;

               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mfaaa +=c1o216;
               //mfaab +=c1o54;
               //mfaac +=c1o216;
               //mfaba +=c1o54;
               //mfabb +=c2o27;
               //mfabc +=c1o54;
               //mfbaa +=c1o54;
               //mfbab +=c2o27;
               //mfbac +=c1o54;
               //mfbba +=c2o27;
               //mfbbb +=c8o27;
               //mfbbc +=c2o27;
               //mfaca +=c1o216;
               //mfacb +=c1o54;
               //mfacc +=c1o216;
               //mfcaa +=c1o216;
               //mfcab +=c1o54;
               //mfcac +=c1o216;
               //mfcca +=c1o216;
               //mfccb +=c1o54;
               //mfccc +=c1o216;
               //mfbca +=c1o54;
               //mfbcb +=c2o27;
               //mfbcc +=c1o54;
               //mfcba +=c1o54;
               //mfcbb +=c2o27;
               //mfcbc +=c1o54;

               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               LBMReal rho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb; 
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = rho - rho_post;
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
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3)    = mfabb;
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3)    = mfbab;
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3)    = mfbba;
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3)   = mfaab;
               (*this->localDistributions)(D3Q27System::ET_NW,x1+1,x2,  x3)   = mfcab;
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3)   = mfaba;
               (*this->localDistributions)(D3Q27System::ET_TW,x1+1,x2,  x3)   = mfcba;
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3)   = mfbaa;
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2+1,x3)   = mfbca;
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3)  = mfaaa;
               (*this->localDistributions)(D3Q27System::ET_TNW,x1+1,x2,  x3)  = mfcaa;
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2+1,x3)  = mfaca;
               (*this->localDistributions)(D3Q27System::ET_TSW,x1+1,x2+1,x3)  = mfcca;

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1+1,x2,  x3    ) = mfcbb;
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2+1,x3    ) = mfbcb;
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3+1  ) = mfbbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1+1,x2+1,x3   ) = mfccb;
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2+1,x3   ) = mfacb;
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1+1,x2,  x3+1 ) = mfcbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3+1 ) = mfabc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2+1,x3+1 ) = mfbcc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3+1 ) = mfbac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1+1,x2+1,x3+1) = mfccc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2+1,x3+1) = mfacc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1+1,x2,  x3+1) = mfcac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3+1) = mfaac;

               (*this->zeroDistributions)(x1,x2,x3) = mfbbb;
               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
double LBMKernelETD3Q27CCLB_Geier::getCallculationTime()
{
   //return timer.getDuration();
   return timer.getTotalTime();
}
