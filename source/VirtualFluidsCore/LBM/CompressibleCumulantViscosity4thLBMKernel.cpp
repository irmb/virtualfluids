#include "CompressibleCumulantViscosity4thLBMKernel.h"
#include "D3Q27System.h"
#include "InterpolationProcessor.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include <math.h>
#include "DataSet3D.h"
#include "LBMKernel.h"

#define PROOF_CORRECTNESS

//////////////////////////////////////////////////////////////////////////
CompressibleCumulantViscosity4thLBMKernel::CompressibleCumulantViscosity4thLBMKernel()
{
   this->nx1 = 0;
   this->nx2 = 0;
   this->nx3 = 0;
   this->parameter = NORMAL;
   this->OxyyMxzz = 1.0;
   this->compressible = true;
   this->bulkOmegaToOmega = false;
   this->OxxPyyPzz = 1.0;
}
//////////////////////////////////////////////////////////////////////////
CompressibleCumulantViscosity4thLBMKernel::CompressibleCumulantViscosity4thLBMKernel(int nx1, int nx2, int nx3, Parameter p)
{
   this->nx1 = nx1;
   this->nx2 = nx2;
   this->nx3 = nx3;
   this->parameter = p;
   this->OxyyMxzz = 1.0;
   this->compressible = true;
   this->bulkOmegaToOmega = false;
   this->OxxPyyPzz = 1.0;
}
//////////////////////////////////////////////////////////////////////////
CompressibleCumulantViscosity4thLBMKernel::~CompressibleCumulantViscosity4thLBMKernel(void)
{

}
//////////////////////////////////////////////////////////////////////////
void CompressibleCumulantViscosity4thLBMKernel::init()
{
   SPtr<DistributionArray3D> d(new D3Q27EsoTwist3DSplittedVector(nx1+2, nx2+2, nx3+2, 0.0));
   dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> CompressibleCumulantViscosity4thLBMKernel::clone()
{
   SPtr<LBMKernel> kernel(new CompressibleCumulantViscosity4thLBMKernel(nx1, nx2, nx3, parameter));
   dynamicPointerCast<CompressibleCumulantViscosity4thLBMKernel>(kernel)->init();
   kernel->setCollisionFactor(this->collFactor);
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
      dynamicPointerCast<CompressibleCumulantViscosity4thLBMKernel>(kernel)->OxyyMxzz = 1.0;
      break;
   case MAGIC:
      dynamicPointerCast<CompressibleCumulantViscosity4thLBMKernel>(kernel)->OxyyMxzz = 2.0 +(-collFactor);
      break;
   }

   if (bulkOmegaToOmega)
   {
      dynamicPointerCast<CompressibleCumulantViscosity4thLBMKernel>(kernel)->OxxPyyPzz = collFactor;
   }
   else
   {
      dynamicPointerCast<CompressibleCumulantViscosity4thLBMKernel>(kernel)->OxxPyyPzz = one;
   }
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleCumulantViscosity4thLBMKernel::calculate()
{
   timer.resetAndStart();
   collideAll();
   timer.stop();
}
//////////////////////////////////////////////////////////////////////////
void CompressibleCumulantViscosity4thLBMKernel::collideAll()
{
   using namespace D3Q27System;
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

      muNu = (1.0/3.0)*(1.0/collFactor - 1.0/2.0);

      muForcingX1.DefineVar("nu", &muNu);
      muForcingX2.DefineVar("nu", &muNu);
      muForcingX3.DefineVar("nu", &muNu);

      LBMReal forcingX1 = 0;
      LBMReal forcingX2 = 0;
      LBMReal forcingX3 = 0;
   }
   /////////////////////////////////////

   localDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray->getNX1();
   const int bcArrayMaxX2 = (int)bcArray->getNX2();
   const int bcArrayMaxX3 = (int)bcArray->getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

   LBMReal omega = collFactor;
   LBMReal OxyyPxzz  = eight*(-two+omega)*(one+two*omega)/(-eight-fourteen*omega+seven*omega*omega);//one;
   LBMReal OxyyMxzz  = eight*(-two+omega)*(-seven+four*omega)/(fiftysix-fifty*omega+nine*omega*omega);//one;
   LBMReal Oxyz      = twentyfour*(-two+omega)*(-two-seven*omega+three*omega*omega)/(fourtyeight+c152*omega-c130*omega*omega+twentynine*omega*omega*omega);


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

               LBMReal mfbbb = (*this->zeroDistributions)(x1, x2, x3);

               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal drho = ((((mfccc+mfaaa)+(mfaca+mfcac))+((mfacc+mfcaa)+(mfaac+mfcca)))+
                  (((mfbac+mfbca)+(mfbaa+mfbcc))+((mfabc+mfcba)+(mfaba+mfcbc))+((mfacb+mfcab)+(mfaab+mfccb)))+
                  ((mfabb+mfcbb)+(mfbab+mfbcb))+(mfbba+mfbbc))+mfbbb;

               LBMReal rho = one+drho;
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal vvx = ((((mfccc-mfaaa)+(mfcac-mfaca))+((mfcaa-mfacc)+(mfcca-mfaac)))+
                  (((mfcba-mfabc)+(mfcbc-mfaba))+((mfcab-mfacb)+(mfccb-mfaab)))+
                  (mfcbb-mfabb))/rho;
               LBMReal vvy = ((((mfccc-mfaaa)+(mfaca-mfcac))+((mfacc-mfcaa)+(mfcca-mfaac)))+
                  (((mfbca-mfbac)+(mfbcc-mfbaa))+((mfacb-mfcab)+(mfccb-mfaab)))+
                  (mfbcb-mfbab))/rho;
               LBMReal vvz = ((((mfccc-mfaaa)+(mfcac-mfaca))+((mfacc-mfcaa)+(mfaac-mfcca)))+
                  (((mfbac-mfbca)+(mfbcc-mfbaa))+((mfabc-mfcba)+(mfcbc-mfaba)))+
                  (mfbbc-mfbba))/rho;
               ////////////////////////////////////////////////////////////////////////////////////

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
         ////////////////////////////////////////////////////////////////////////////////////
               LBMReal oMdrho = one; // comp special
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal m0, m1, m2;
               LBMReal vx2;
               LBMReal vy2;
               LBMReal vz2;
               vx2 = vvx*vvx;
               vy2 = vvy*vvy;
               vz2 = vvz*vvz;
               ////////////////////////////////////////////////////////////////////////////////////
               LBMReal wadjust;
               LBMReal qudricLimitP = 0.01;// * 0.0001f;
               LBMReal qudricLimitM = 0.01;// * 0.0001f;
               LBMReal qudricLimitD = 0.01;// * 0.001f;
               //LBMReal s9 = minusomega;
               //test
               //s9 = 0.;
               ////////////////////////////////////////////////////////////////////////////////////
               //Hin
               ////////////////////////////////////////////////////////////////////////////////////
               // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Z - Dir
               m2 = mfaaa+mfaac;
               m1 = mfaac-mfaaa;
               m0 = m2+mfaab;
               mfaaa = m0;
               m0 += c1o36 * oMdrho;
               mfaab = m1-m0 * vvz;
               mfaac = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaba+mfabc;
               m1 = mfabc-mfaba;
               m0 = m2+mfabb;
               mfaba = m0;
               m0 += c1o9 * oMdrho;
               mfabb = m1-m0 * vvz;
               mfabc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaca+mfacc;
               m1 = mfacc-mfaca;
               m0 = m2+mfacb;
               mfaca = m0;
               m0 += c1o36 * oMdrho;
               mfacb = m1-m0 * vvz;
               mfacc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbaa+mfbac;
               m1 = mfbac-mfbaa;
               m0 = m2+mfbab;
               mfbaa = m0;
               m0 += c1o9 * oMdrho;
               mfbab = m1-m0 * vvz;
               mfbac = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbba+mfbbc;
               m1 = mfbbc-mfbba;
               m0 = m2+mfbbb;
               mfbba = m0;
               m0 += c4o9 * oMdrho;
               mfbbb = m1-m0 * vvz;
               mfbbc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbca+mfbcc;
               m1 = mfbcc-mfbca;
               m0 = m2+mfbcb;
               mfbca = m0;
               m0 += c1o9 * oMdrho;
               mfbcb = m1-m0 * vvz;
               mfbcc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcaa+mfcac;
               m1 = mfcac-mfcaa;
               m0 = m2+mfcab;
               mfcaa = m0;
               m0 += c1o36 * oMdrho;
               mfcab = m1-m0 * vvz;
               mfcac = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcba+mfcbc;
               m1 = mfcbc-mfcba;
               m0 = m2+mfcbb;
               mfcba = m0;
               m0 += c1o9 * oMdrho;
               mfcbb = m1-m0 * vvz;
               mfcbc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcca+mfccc;
               m1 = mfccc-mfcca;
               m0 = m2+mfccb;
               mfcca = m0;
               m0 += c1o36 * oMdrho;
               mfccb = m1-m0 * vvz;
               mfccc = m2-two*	m1 * vvz+vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m2 = mfaaa+mfaca;
               m1 = mfaca-mfaaa;
               m0 = m2+mfaba;
               mfaaa = m0;
               m0 += c1o6 * oMdrho;
               mfaba = m1-m0 * vvy;
               mfaca = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaab+mfacb;
               m1 = mfacb-mfaab;
               m0 = m2+mfabb;
               mfaab = m0;
               mfabb = m1-m0 * vvy;
               mfacb = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaac+mfacc;
               m1 = mfacc-mfaac;
               m0 = m2+mfabc;
               mfaac = m0;
               m0 += c1o18 * oMdrho;
               mfabc = m1-m0 * vvy;
               mfacc = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbaa+mfbca;
               m1 = mfbca-mfbaa;
               m0 = m2+mfbba;
               mfbaa = m0;
               m0 += c2o3 * oMdrho;
               mfbba = m1-m0 * vvy;
               mfbca = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbab+mfbcb;
               m1 = mfbcb-mfbab;
               m0 = m2+mfbbb;
               mfbab = m0;
               mfbbb = m1-m0 * vvy;
               mfbcb = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfbac+mfbcc;
               m1 = mfbcc-mfbac;
               m0 = m2+mfbbc;
               mfbac = m0;
               m0 += c2o9 * oMdrho;
               mfbbc = m1-m0 * vvy;
               mfbcc = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcaa+mfcca;
               m1 = mfcca-mfcaa;
               m0 = m2+mfcba;
               mfcaa = m0;
               m0 += c1o6 * oMdrho;
               mfcba = m1-m0 * vvy;
               mfcca = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcab+mfccb;
               m1 = mfccb-mfcab;
               m0 = m2+mfcbb;
               mfcab = m0;
               mfcbb = m1-m0 * vvy;
               mfccb = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfcac+mfccc;
               m1 = mfccc-mfcac;
               m0 = m2+mfcbc;
               mfcac = m0;
               m0 += c1o18 * oMdrho;
               mfcbc = m1-m0 * vvy;
               mfccc = m2-two*	m1 * vvy+vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m2 = mfaaa+mfcaa;
               m1 = mfcaa-mfaaa;
               m0 = m2+mfbaa;
               mfaaa = m0;
               m0 += one* oMdrho;
               mfbaa = m1-m0 * vvx;
               mfcaa = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaba+mfcba;
               m1 = mfcba-mfaba;
               m0 = m2+mfbba;
               mfaba = m0;
               mfbba = m1-m0 * vvx;
               mfcba = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaca+mfcca;
               m1 = mfcca-mfaca;
               m0 = m2+mfbca;
               mfaca = m0;
               m0 += c1o3 * oMdrho;
               mfbca = m1-m0 * vvx;
               mfcca = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaab+mfcab;
               m1 = mfcab-mfaab;
               m0 = m2+mfbab;
               mfaab = m0;
               mfbab = m1-m0 * vvx;
               mfcab = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfabb+mfcbb;
               m1 = mfcbb-mfabb;
               m0 = m2+mfbbb;
               mfabb = m0;
               mfbbb = m1-m0 * vvx;
               mfcbb = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfacb+mfccb;
               m1 = mfccb-mfacb;
               m0 = m2+mfbcb;
               mfacb = m0;
               mfbcb = m1-m0 * vvx;
               mfccb = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfaac+mfcac;
               m1 = mfcac-mfaac;
               m0 = m2+mfbac;
               mfaac = m0;
               m0 += c1o3 * oMdrho;
               mfbac = m1-m0 * vvx;
               mfcac = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfabc+mfcbc;
               m1 = mfcbc-mfabc;
               m0 = m2+mfbbc;
               mfabc = m0;
               mfbbc = m1-m0 * vvx;
               mfcbc = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2 = mfacc+mfccc;
               m1 = mfccc-mfacc;
               m0 = m2+mfbcc;
               mfacc = m0;
               m0 += c1o9 * oMdrho;
               mfbcc = m1-m0 * vvx;
               mfccc = m2-two*	m1 * vvx+vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////


               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               ////////////////////////////////////////////////////////////////////////////////////

               //LBMReal OxxPyyPzz = one; // bulk viscosity

               ////////////////////////////////////////////////////////////
               //3.
               //////////////////////////////
               //LBMReal OxyyPxzz = one;//three  * (two - omega) / (three  - omega);//
               //LBMReal OxyyMxzz = one;//six    * (two - omega) / (six    - omega);//
               //LBMReal Oxyz = one;//twelve * (two - omega) / (twelve + omega);//

               //LBMReal OxyyPxzz  = eight*(-two+omega)*(one+two*omega)/(-eight-fourteen*omega+seven*omega*omega);//one;

               //LBMReal OxyyMxzz  = eight*(-two+omega)*(-seven+four*omega)/(fiftysix-fifty*omega+nine*omega*omega);//one;

               //LBMReal Oxyz      = twentyfour*(-two+omega)*(-two-seven*omega+three*omega*omega)/(fourtyeight+c152*omega-c130*omega*omega+twentynine*omega*omega*omega);
               //////////////////////////////
               //LBMReal OxyyPxzz  = two-omega;//
               //LBMReal OxyyMxzz  = two-omega;//
               //////////////////////////////
               //LBMReal OxyyPxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
               //LBMReal OxyyMxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
               //////////////////////////////
               //LBMReal OxyyPxzz  = omega;//BGK
               //LBMReal OxyyMxzz  = omega;//BGK
               //////////////////////////////
               //LBMReal OxyyPxzz  = (one + omega) / two;//1P5
               //LBMReal OxyyMxzz  = (one + omega) / two;//1P5
               //////////////////////////////
               //LBMReal OxyyPxzz  = (three - omega) / two;//0P5
               //LBMReal OxyyMxzz  = (three - omega) / two;//0P5
               //////////////////////////////
               //LBMReal OxyyPxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
               //LBMReal OxyyMxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
               ////////////////////////////////////////////////////////////
               //4.
               //////////////////////////////
               LBMReal O4 = one;
               //////////////////////////////
               //LBMReal O4        = omega;//TRT
               ////////////////////////////////////////////////////////////
               //5.
               //////////////////////////////
               LBMReal O5 = one;
               ////////////////////////////////////////////////////////////
               //6.
               //////////////////////////////
               LBMReal O6 = one;
               ////////////////////////////////////////////////////////////


               //central moments to cumulants
               //4.
               LBMReal CUMcbb = mfcbb-((mfcaa+c1o3) * mfabb+two * mfbba * mfbab)/rho;	//ab 15.05.2015 verwendet
               LBMReal CUMbcb = mfbcb-((mfaca+c1o3) * mfbab+two * mfbba * mfabb)/rho; //ab 15.05.2015 verwendet
               LBMReal CUMbbc = mfbbc-((mfaac+c1o3) * mfbba+two * mfbab * mfabb)/rho; //ab 15.05.2015 verwendet

               LBMReal CUMcca = mfcca-(((mfcaa * mfaca+two * mfbba * mfbba)+c1o3 * (mfcaa+mfaca))/rho-c1o9*(drho/rho));
               LBMReal CUMcac = mfcac-(((mfcaa * mfaac+two * mfbab * mfbab)+c1o3 * (mfcaa+mfaac))/rho-c1o9*(drho/rho));
               LBMReal CUMacc = mfacc-(((mfaac * mfaca+two * mfabb * mfabb)+c1o3 * (mfaac+mfaca))/rho-c1o9*(drho/rho));

               //5.
               LBMReal CUMbcc = mfbcc-((mfaac * mfbca+mfaca * mfbac+four * mfabb * mfbbb+two * (mfbab * mfacb+mfbba * mfabc))+c1o3 * (mfbca+mfbac))/rho;
               LBMReal CUMcbc = mfcbc-((mfaac * mfcba+mfcaa * mfabc+four * mfbab * mfbbb+two * (mfabb * mfcab+mfbba * mfbac))+c1o3 * (mfcba+mfabc))/rho;
               LBMReal CUMccb = mfccb-((mfcaa * mfacb+mfaca * mfcab+four * mfbba * mfbbb+two * (mfbab * mfbca+mfabb * mfcba))+c1o3 * (mfacb+mfcab))/rho;

               //6.

               LBMReal CUMccc = mfccc+((-four *  mfbbb * mfbbb
                  -(mfcaa * mfacc+mfaca * mfcac+mfaac * mfcca)
                  -four * (mfabb * mfcbb+mfbab * mfbcb+mfbba * mfbbc)
                  -two * (mfbca * mfbac+mfcba * mfabc+mfcab * mfacb))/rho
                  +(four * (mfbab * mfbab * mfaca+mfabb * mfabb * mfcaa+mfbba * mfbba * mfaac)
                     +two * (mfcaa * mfaca * mfaac)
                     +sixteen *  mfbba * mfbab * mfabb)/(rho * rho)
                  -c1o3 * (mfacc+mfcac+mfcca)/rho
                  -c1o9 * (mfcaa+mfaca+mfaac)/rho
                  +(two * (mfbab * mfbab+mfabb * mfabb+mfbba * mfbba)
                     +(mfaac * mfaca+mfaac * mfcaa+mfaca * mfcaa)+c1o3 *(mfaac+mfaca+mfcaa))/(rho * rho) * c2o3
                  +c1o27*((drho * drho-drho)/(rho*rho)));
               //+ c1o27*(one -three/rho +two/(rho*rho)));




   //2.
   // linear combinations
               LBMReal mxxPyyPzz = mfcaa+mfaca+mfaac;
               LBMReal mxxMyy = mfcaa-mfaca;
               LBMReal mxxMzz = mfcaa-mfaac;

               //////////////////////////////////////////////////////////////////////////
      // 			LBMReal magicBulk=(CUMacc+CUMcac+CUMcca)*(one/OxxPyyPzz-c1o2)*c3o2*8.;

               //////////////////////////////////////////////////////////////////////////
               //limiter-Scheise Teil 1
               //LBMReal oxxyy,oxxzz,oxy,oxz,oyz;
               //LBMReal smag=0.001;
               //oxxyy    = omega+(one-omega)*abs(mxxMyy)/(abs(mxxMyy)+smag);
               //oxxzz    = omega+(one-omega)*abs(mxxMzz)/(abs(mxxMzz)+smag);
               //oxy      = omega+(one-omega)*abs(mfbba)/(abs(mfbba)+smag);
               //oxz      = omega+(one-omega)*abs(mfbab)/(abs(mfbab)+smag);
               //oyz      = omega+(one-omega)*abs(mfabb)/(abs(mfabb)+smag);

               ////////////////////////////////////////////////////////////////////////////
               ////Teil 1b
               //LBMReal constante = 1000.0;
               //LBMReal nuEddi = constante * abs(mxxPyyPzz);
               //LBMReal omegaLimit = one / (one / omega + three * nuEddi);

               //{
               //	LBMReal dxux = c1o2 * (-omegaLimit) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
               //	LBMReal dyuy = dxux + omegaLimit * c3o2 * mxxMyy;
               //	LBMReal dzuz = dxux + omegaLimit * c3o2 * mxxMzz;

                  ////relax
                  //mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
                  //mxxMyy    += omegaLimit * (-mxxMyy) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vy2 * dyuy);
                  //mxxMzz    += omegaLimit * (-mxxMzz) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vz2 * dzuz);

               //}
               //mfabb     += omegaLimit * (-mfabb);
               //mfbab     += omegaLimit * (-mfbab);
               //mfbba     += omegaLimit * (-mfbba);
               ////////////////////////////////////////////////////////////////////////////

               ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)

               LBMReal dxux = c1o2 * (-omega) *(mxxMyy+mxxMzz)+c1o2 *  OxxPyyPzz * (mfaaa-mxxPyyPzz);
               LBMReal dyuy = dxux+omega * c3o2 * mxxMyy;
               LBMReal dzuz = dxux+omega * c3o2 * mxxMzz;

               LBMReal Dxy =-three*omega*mfbba;
               LBMReal Dxz =-three*omega*mfbab;
               LBMReal Dyz =-three*omega*mfabb;



               //relax
               mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz)-three * (one-c1o2 * OxxPyyPzz) * (vx2 * dxux+vy2 * dyuy+vz2 * dzuz);//-magicBulk*OxxPyyPzz;
               mxxMyy += omega * (-mxxMyy)-three * (one+c1o2 * (-omega)) * (vx2 * dxux-vy2 * dyuy);
               mxxMzz += omega * (-mxxMzz)-three * (one+c1o2 * (-omega)) * (vx2 * dxux-vz2 * dzuz);

               //////////////////////////////////////////////////////////////////////////
               //limiter-Scheise Teil 2
               //mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
               //mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
               //////////////////////////////////////////////////////////////////////////


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////no correction
            //mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
            //mxxMyy    += -(-omega) * (-mxxMyy);
            //mxxMzz    += -(-omega) * (-mxxMzz);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               mfabb += omega * (-mfabb);
               mfbab += omega * (-mfbab);
               mfbba += omega * (-mfbba);

               //////////////////////////////////////////////////////////////////////////
               //limiter-Scheise Teil 3
               //mfabb     += oyz * (-mfabb);
               //mfbab     += oxz * (-mfbab);
               //mfbba     += oxy * (-mfbba);
               //////////////////////////////////////////////////////////////////////////

               // linear combinations back
               mfcaa = c1o3 * (mxxMyy+mxxMzz+mxxPyyPzz);
               mfaca = c1o3 * (-two*  mxxMyy+mxxMzz+mxxPyyPzz);
               mfaac = c1o3 * (mxxMyy-two* mxxMzz+mxxPyyPzz);

               //3.
               // linear combinations

               LBMReal mxxyPyzz = mfcba+mfabc;
               LBMReal mxxyMyzz = mfcba-mfabc;

               LBMReal mxxzPyyz = mfcab+mfacb;
               LBMReal mxxzMyyz = mfcab-mfacb;

               LBMReal mxyyPxzz = mfbca+mfbac;
               LBMReal mxyyMxzz = mfbca-mfbac;

               //relax
               //////////////////////////////////////////////////////////////////////////
               //das ist der limiter
               wadjust = Oxyz+(one-Oxyz)*abs(mfbbb)/(abs(mfbbb)+qudricLimitD);
               mfbbb += wadjust * (-mfbbb);
               wadjust = OxyyPxzz+(one-OxyyPxzz)*abs(mxxyPyzz)/(abs(mxxyPyzz)+qudricLimitP);
               mxxyPyzz += wadjust * (-mxxyPyzz);
               wadjust = OxyyMxzz+(one-OxyyMxzz)*abs(mxxyMyzz)/(abs(mxxyMyzz)+qudricLimitM);
               mxxyMyzz += wadjust * (-mxxyMyzz);
               wadjust = OxyyPxzz+(one-OxyyPxzz)*abs(mxxzPyyz)/(abs(mxxzPyyz)+qudricLimitP);
               mxxzPyyz += wadjust * (-mxxzPyyz);
               wadjust = OxyyMxzz+(one-OxyyMxzz)*abs(mxxzMyyz)/(abs(mxxzMyyz)+qudricLimitM);
               mxxzMyyz += wadjust * (-mxxzMyyz);
               wadjust = OxyyPxzz+(one-OxyyPxzz)*abs(mxyyPxzz)/(abs(mxyyPxzz)+qudricLimitP);
               mxyyPxzz += wadjust * (-mxyyPxzz);
               wadjust = OxyyMxzz+(one-OxyyMxzz)*abs(mxyyMxzz)/(abs(mxyyMxzz)+qudricLimitM);
               mxyyMxzz += wadjust * (-mxyyMxzz);
               //////////////////////////////////////////////////////////////////////////
               //ohne limiter
               //mfbbb     += OxyyMxzz * (-mfbbb);
               //mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
               //mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
               //mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
               //mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
               //mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
               //mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
               //////////////////////////////////////////////////////////////////////////

               //// linear combinations back
               mfcba = (mxxyMyzz+mxxyPyzz) * c1o2;
               mfabc = (-mxxyMyzz+mxxyPyzz) * c1o2;
               mfcab = (mxxzMyyz+mxxzPyyz) * c1o2;
               mfacb = (-mxxzMyyz+mxxzPyyz) * c1o2;
               mfbca = (mxyyMxzz+mxyyPxzz) * c1o2;
               mfbac = (-mxyyMxzz+mxyyPxzz) * c1o2;

               //4.
               //////////////////////////////////////////////////////////////////////////
               //mit limiter
            //	wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
               //CUMacc    += wadjust * (-CUMacc);
            //	wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
               //CUMcac    += wadjust * (-CUMcac); 
            //	wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
               //CUMcca    += wadjust * (-CUMcca); 

            //	wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
               //CUMbbc    += wadjust * (-CUMbbc); 
            //	wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
               //CUMbcb    += wadjust * (-CUMbcb); 
            //	wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
               //CUMcbb    += wadjust * (-CUMcbb); 
               //////////////////////////////////////////////////////////////////////////
               //////////////////////////////////////////////////////////////////////////
               LBMReal A = (four + two*omega - three*omega*omega) / (two - seven*omega + five*omega*omega);
               LBMReal B = (four + twentyeight*omega - fourteen*omega*omega) / (six - twentyone*omega + fiveteen*omega*omega);
               //////////////////////////////////////////////////////////////////////////


               //ohne limiter
               //CUMacc += O4 * (-CUMacc);
               //CUMcac += O4 * (-CUMcac);
               //CUMcca += O4 * (-CUMcca);

               //CUMbbc += O4 * (-CUMbbc);
               //CUMbcb += O4 * (-CUMbcb);
               //CUMcbb += O4 * (-CUMcbb);
               CUMacc = -O4*(one / omega - c1o2) * (dyuy + dzuz) * c2o3 * A + (one - O4) * (CUMacc);
               CUMcac = -O4*(one / omega - c1o2) * (dxux + dzuz) * c2o3 * A + (one - O4) * (CUMcac);
               CUMcca = -O4*(one / omega - c1o2) * (dyuy + dxux) * c2o3 * A + (one - O4) * (CUMcca);
               CUMbbc = -O4*(one / omega - c1o2) * Dxy           * c1o3 * B + (one - O4) * (CUMbbc);
               CUMbcb = -O4*(one / omega - c1o2) * Dxz           * c1o3 * B + (one - O4) * (CUMbcb);
               CUMcbb = -O4*(one / omega - c1o2) * Dyz           * c1o3 * B + (one - O4) * (CUMcbb);
               //////////////////////////////////////////////////////////////////////////


               //5.
               CUMbcc += O5 * (-CUMbcc);
               CUMcbc += O5 * (-CUMcbc);
               CUMccb += O5 * (-CUMccb);

               //6.
               CUMccc += O6 * (-CUMccc);



               //back cumulants to central moments
               //4.
               mfcbb = CUMcbb+((mfcaa+c1o3) * mfabb+two * mfbba * mfbab)/rho;
               mfbcb = CUMbcb+((mfaca+c1o3) * mfbab+two * mfbba * mfabb)/rho;
               mfbbc = CUMbbc+((mfaac+c1o3) * mfbba+two * mfbab * mfabb)/rho;

               mfcca = CUMcca+(((mfcaa * mfaca+two * mfbba * mfbba)+c1o3 * (mfcaa+mfaca))/rho-c1o9*(drho/rho));//(one/rho-one));
               mfcac = CUMcac+(((mfcaa * mfaac+two * mfbab * mfbab)+c1o3 * (mfcaa+mfaac))/rho-c1o9*(drho/rho));//(one/rho-one));
               mfacc = CUMacc+(((mfaac * mfaca+two * mfabb * mfabb)+c1o3 * (mfaac+mfaca))/rho-c1o9*(drho/rho));//(one/rho-one));

               //5.
               mfbcc = CUMbcc+((mfaac * mfbca+mfaca * mfbac+four * mfabb * mfbbb+two * (mfbab * mfacb+mfbba * mfabc))+c1o3 * (mfbca+mfbac))/rho;
               mfcbc = CUMcbc+((mfaac * mfcba+mfcaa * mfabc+four * mfbab * mfbbb+two * (mfabb * mfcab+mfbba * mfbac))+c1o3 * (mfcba+mfabc))/rho;
               mfccb = CUMccb+((mfcaa * mfacb+mfaca * mfcab+four * mfbba * mfbbb+two * (mfbab * mfbca+mfabb * mfcba))+c1o3 * (mfacb+mfcab))/rho;

               //6.

               mfccc = CUMccc-((-four *  mfbbb * mfbbb
                  -(mfcaa * mfacc+mfaca * mfcac+mfaac * mfcca)
                  -four * (mfabb * mfcbb+mfbab * mfbcb+mfbba * mfbbc)
                  -two * (mfbca * mfbac+mfcba * mfabc+mfcab * mfacb))/rho
                  +(four * (mfbab * mfbab * mfaca+mfabb * mfabb * mfcaa+mfbba * mfbba * mfaac)
                     +two * (mfcaa * mfaca * mfaac)
                     +sixteen *  mfbba * mfbab * mfabb)/(rho * rho)
                  -c1o3 * (mfacc+mfcac+mfcca)/rho
                  -c1o9 * (mfcaa+mfaca+mfaac)/rho
                  +(two * (mfbab * mfbab+mfabb * mfabb+mfbba * mfbba)
                     +(mfaac * mfaca+mfaac * mfcaa+mfaca * mfcaa)+c1o3 *(mfaac+mfaca+mfcaa))/(rho * rho) * c2o3
                  +c1o27*((drho * drho-drho)/(rho*rho)));
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
               m0 = mfaac * c1o2+mfaab * (vvz-c1o2)+(mfaaa+one* oMdrho) * (vz2-vvz) * c1o2;
               m1 = -mfaac-two* mfaab *  vvz+mfaaa                * (one-vz2)-one* oMdrho * vz2;
               m2 = mfaac * c1o2+mfaab * (vvz+c1o2)+(mfaaa+one* oMdrho) * (vz2+vvz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfabc * c1o2+mfabb * (vvz-c1o2)+mfaba * (vz2-vvz) * c1o2;
               m1 = -mfabc-two* mfabb *  vvz+mfaba * (one-vz2);
               m2 = mfabc * c1o2+mfabb * (vvz+c1o2)+mfaba * (vz2+vvz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfacc * c1o2+mfacb * (vvz-c1o2)+(mfaca+c1o3 * oMdrho) * (vz2-vvz) * c1o2;
               m1 = -mfacc-two* mfacb *  vvz+mfaca                  * (one-vz2)-c1o3 * oMdrho * vz2;
               m2 = mfacc * c1o2+mfacb * (vvz+c1o2)+(mfaca+c1o3 * oMdrho) * (vz2+vvz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfbac * c1o2+mfbab * (vvz-c1o2)+mfbaa * (vz2-vvz) * c1o2;
               m1 = -mfbac-two* mfbab *  vvz+mfbaa * (one-vz2);
               m2 = mfbac * c1o2+mfbab * (vvz+c1o2)+mfbaa * (vz2+vvz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 = mfbbc * c1o2+mfbbb * (vvz-c1o2)+mfbba * (vz2-vvz) * c1o2;
               m1 = -mfbbc-two* mfbbb *  vvz+mfbba * (one-vz2);
               m2 = mfbbc * c1o2+mfbbb * (vvz+c1o2)+mfbba * (vz2+vvz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 = mfbcc * c1o2+mfbcb * (vvz-c1o2)+mfbca * (vz2-vvz) * c1o2;
               m1 = -mfbcc-two* mfbcb *  vvz+mfbca * (one-vz2);
               m2 = mfbcc * c1o2+mfbcb * (vvz+c1o2)+mfbca * (vz2+vvz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcac * c1o2+mfcab * (vvz-c1o2)+(mfcaa+c1o3 * oMdrho) * (vz2-vvz) * c1o2;
               m1 = -mfcac-two* mfcab *  vvz+mfcaa                  * (one-vz2)-c1o3 * oMdrho * vz2;
               m2 = mfcac * c1o2+mfcab * (vvz+c1o2)+(mfcaa+c1o3 * oMdrho) * (vz2+vvz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 = mfcbc * c1o2+mfcbb * (vvz-c1o2)+mfcba * (vz2-vvz) * c1o2;
               m1 = -mfcbc-two* mfcbb *  vvz+mfcba * (one-vz2);
               m2 = mfcbc * c1o2+mfcbb * (vvz+c1o2)+mfcba * (vz2+vvz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 = mfccc * c1o2+mfccb * (vvz-c1o2)+(mfcca+c1o9 * oMdrho) * (vz2-vvz) * c1o2;
               m1 = -mfccc-two* mfccb *  vvz+mfcca                  * (one-vz2)-c1o9 * oMdrho * vz2;
               m2 = mfccc * c1o2+mfccb * (vvz+c1o2)+(mfcca+c1o9 * oMdrho) * (vz2+vvz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 = mfaca * c1o2+mfaba * (vvy-c1o2)+(mfaaa+c1o6 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfaca-two* mfaba *  vvy+mfaaa                  * (one-vy2)-c1o6 * oMdrho * vy2;
               m2 = mfaca * c1o2+mfaba * (vvy+c1o2)+(mfaaa+c1o6 * oMdrho) * (vy2+vvy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfacb * c1o2+mfabb * (vvy-c1o2)+(mfaab+c2o3 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfacb-two* mfabb *  vvy+mfaab                  * (one-vy2)-c2o3 * oMdrho * vy2;
               m2 = mfacb * c1o2+mfabb * (vvy+c1o2)+(mfaab+c2o3 * oMdrho) * (vy2+vvy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfacc * c1o2+mfabc * (vvy-c1o2)+(mfaac+c1o6 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfacc-two* mfabc *  vvy+mfaac                  * (one-vy2)-c1o6 * oMdrho * vy2;
               m2 = mfacc * c1o2+mfabc * (vvy+c1o2)+(mfaac+c1o6 * oMdrho) * (vy2+vvy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfbca * c1o2+mfbba * (vvy-c1o2)+mfbaa * (vy2-vvy) * c1o2;
               m1 = -mfbca-two* mfbba *  vvy+mfbaa * (one-vy2);
               m2 = mfbca * c1o2+mfbba * (vvy+c1o2)+mfbaa * (vy2+vvy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 = mfbcb * c1o2+mfbbb * (vvy-c1o2)+mfbab * (vy2-vvy) * c1o2;
               m1 = -mfbcb-two* mfbbb *  vvy+mfbab * (one-vy2);
               m2 = mfbcb * c1o2+mfbbb * (vvy+c1o2)+mfbab * (vy2+vvy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 = mfbcc * c1o2+mfbbc * (vvy-c1o2)+mfbac * (vy2-vvy) * c1o2;
               m1 = -mfbcc-two* mfbbc *  vvy+mfbac * (one-vy2);
               m2 = mfbcc * c1o2+mfbbc * (vvy+c1o2)+mfbac * (vy2+vvy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcca * c1o2+mfcba * (vvy-c1o2)+(mfcaa+c1o18 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfcca-two* mfcba *  vvy+mfcaa                   * (one-vy2)-c1o18 * oMdrho * vy2;
               m2 = mfcca * c1o2+mfcba * (vvy+c1o2)+(mfcaa+c1o18 * oMdrho) * (vy2+vvy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 = mfccb * c1o2+mfcbb * (vvy-c1o2)+(mfcab+c2o9 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfccb-two* mfcbb *  vvy+mfcab                  * (one-vy2)-c2o9 * oMdrho * vy2;
               m2 = mfccb * c1o2+mfcbb * (vvy+c1o2)+(mfcab+c2o9 * oMdrho) * (vy2+vvy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 = mfccc * c1o2+mfcbc * (vvy-c1o2)+(mfcac+c1o18 * oMdrho) * (vy2-vvy) * c1o2;
               m1 = -mfccc-two* mfcbc *  vvy+mfcac                   * (one-vy2)-c1o18 * oMdrho * vy2;
               m2 = mfccc * c1o2+mfcbc * (vvy+c1o2)+(mfcac+c1o18 * oMdrho) * (vy2+vvy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 = mfcaa * c1o2+mfbaa * (vvx-c1o2)+(mfaaa+c1o36 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcaa-two* mfbaa *  vvx+mfaaa                   * (one-vx2)-c1o36 * oMdrho * vx2;
               m2 = mfcaa * c1o2+mfbaa * (vvx+c1o2)+(mfaaa+c1o36 * oMdrho) * (vx2+vvx) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcba * c1o2+mfbba * (vvx-c1o2)+(mfaba+c1o9 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcba-two* mfbba *  vvx+mfaba                  * (one-vx2)-c1o9 * oMdrho * vx2;
               m2 = mfcba * c1o2+mfbba * (vvx+c1o2)+(mfaba+c1o9 * oMdrho) * (vx2+vvx) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcca * c1o2+mfbca * (vvx-c1o2)+(mfaca+c1o36 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcca-two* mfbca *  vvx+mfaca                   * (one-vx2)-c1o36 * oMdrho * vx2;
               m2 = mfcca * c1o2+mfbca * (vvx+c1o2)+(mfaca+c1o36 * oMdrho) * (vx2+vvx) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcab * c1o2+mfbab * (vvx-c1o2)+(mfaab+c1o9 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcab-two* mfbab *  vvx+mfaab                  * (one-vx2)-c1o9 * oMdrho * vx2;
               m2 = mfcab * c1o2+mfbab * (vvx+c1o2)+(mfaab+c1o9 * oMdrho) * (vx2+vvx) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 = mfcbb * c1o2+mfbbb * (vvx-c1o2)+(mfabb+c4o9 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcbb-two* mfbbb *  vvx+mfabb                  * (one-vx2)-c4o9 * oMdrho * vx2;
               m2 = mfcbb * c1o2+mfbbb * (vvx+c1o2)+(mfabb+c4o9 * oMdrho) * (vx2+vvx) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 = mfccb * c1o2+mfbcb * (vvx-c1o2)+(mfacb+c1o9 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfccb-two* mfbcb *  vvx+mfacb                  * (one-vx2)-c1o9 * oMdrho * vx2;
               m2 = mfccb * c1o2+mfbcb * (vvx+c1o2)+(mfacb+c1o9 * oMdrho) * (vx2+vvx) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 = mfcac * c1o2+mfbac * (vvx-c1o2)+(mfaac+c1o36 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcac-two* mfbac *  vvx+mfaac                   * (one-vx2)-c1o36 * oMdrho * vx2;
               m2 = mfcac * c1o2+mfbac * (vvx+c1o2)+(mfaac+c1o36 * oMdrho) * (vx2+vvx) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 = mfcbc * c1o2+mfbbc * (vvx-c1o2)+(mfabc+c1o9 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfcbc-two* mfbbc *  vvx+mfabc                  * (one-vx2)-c1o9 * oMdrho * vx2;
               m2 = mfcbc * c1o2+mfbbc * (vvx+c1o2)+(mfabc+c1o9 * oMdrho) * (vx2+vvx) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 = mfccc * c1o2+mfbcc * (vvx-c1o2)+(mfacc+c1o36 * oMdrho) * (vx2-vvx) * c1o2;
               m1 = -mfccc-two* mfbcc *  vvx+mfacc                   * (one-vx2)-c1o36 * oMdrho * vx2;
               m2 = mfccc * c1o2+mfbcc * (vvx+c1o2)+(mfacc+c1o36 * oMdrho) * (vx2+vvx) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////

               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               LBMReal drho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;
               //LBMReal dif = fabs(rho - rho_post);
               LBMReal dif = drho - drho_post;
#ifdef SINGLEPRECISION
               if (dif > 10.0E-7 || dif < -10.0E-7)
#else
               if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
                  UB_THROW(UbException(UB_EXARGS, "rho="+UbSystem::toString(drho)+", rho_post="+UbSystem::toString(drho_post)
                     +" dif="+UbSystem::toString(dif)
                     +" rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)));
                  //UBLOG(logERROR,"LBMKernelETD3Q27CCLB::collideAll(): rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3));
                  //exit(EXIT_FAILURE);
               }
#endif
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3)    = mfabb;
               (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3)    = mfbab;
               (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3)    = mfbba;
               (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3)   = mfaab;
               (*this->localDistributions)(D3Q27System::ET_NW, x1p, x2, x3)   = mfcab;
               (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3)   = mfaba;
               (*this->localDistributions)(D3Q27System::ET_TW, x1p, x2, x3)   = mfcba;
               (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3)   = mfbaa;
               (*this->localDistributions)(D3Q27System::ET_TS, x1, x2p, x3)   = mfbca;
               (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3)  = mfaaa;
               (*this->localDistributions)(D3Q27System::ET_TNW, x1p, x2, x3)  = mfcaa;
               (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2p, x3)  = mfaca;
               (*this->localDistributions)(D3Q27System::ET_TSW, x1p, x2p, x3)  = mfcca;

               (*this->nonLocalDistributions)(D3Q27System::ET_W, x1p, x2, x3) = mfcbb;
               (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2p, x3) = mfbcb;
               (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3p) = mfbbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1p, x2p, x3) = mfccb;
               (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2p, x3) = mfacb;
               (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1p, x2, x3p) = mfcbc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3p) = mfabc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2p, x3p) = mfbcc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3p) = mfbac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1p, x2p, x3p) = mfccc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2p, x3p) = mfacc;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1p, x2, x3p) = mfcac;
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3p) = mfaac;

               (*this->zeroDistributions)(x1, x2, x3) = mfbbb;
               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
double CompressibleCumulantViscosity4thLBMKernel::getCalculationTime()
{
   //return timer.getDuration();
   return timer.getTotalTime();
}
//////////////////////////////////////////////////////////////////////////
void CompressibleCumulantViscosity4thLBMKernel::setBulkOmegaToOmega(bool value)
{
   bulkOmegaToOmega = value;
}
