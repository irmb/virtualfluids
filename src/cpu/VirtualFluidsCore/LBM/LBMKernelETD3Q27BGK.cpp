#include "LBMKernelETD3Q27BGK.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27EsoTwist3DSoA.h"
#include "DataSet3D.h"
#include "BCSet.h"
#include "BCArray3D.h"
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;
//using namespace UbMath;

//#define PROOF_CORRECTNESS


//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27BGK::LBMKernelETD3Q27BGK() 
{
   this->compressible = false;
}
//////////////////////////////////////////////////////////////////////////
LBMKernelETD3Q27BGK::~LBMKernelETD3Q27BGK(void)
= default;
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27BGK::initDataSet()
{
   SPtr<DistributionArray3D> d(new D3Q27EsoTwist3DSplittedVector(nx[0]+2, nx[1]+2, nx[2]+2, -999.9));
   dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> LBMKernelETD3Q27BGK::clone()
{
   SPtr<LBMKernel> kernel(new LBMKernelETD3Q27BGK());
   std::dynamic_pointer_cast<LBMKernelETD3Q27BGK>(kernel)->initDataSet();
   kernel->setCollisionFactor(this->collFactor);
   kernel->setBCSet(bcSet->clone(kernel));
   kernel->setWithForcing(withForcing);
   kernel->setForcingX1(muForcingX1);
   kernel->setForcingX2(muForcingX2);
   kernel->setForcingX3(muForcingX3);
   kernel->setIndex(ix1, ix2, ix3);
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void LBMKernelETD3Q27BGK::calculate(int  /*step*/)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);
      forcingX1 = c0o1;
      forcingX2 = c0o1;
      forcingX3 = c0o1;
   }
   /////////////////////////////////////

   localDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();
   real f[D3Q27System::ENDF+1];
   real feq[D3Q27System::ENDF+1];
   real drho,vx1,vx2,vx3;
   const int bcArrayMaxX1 = (int)bcArray->getNX1();
   const int bcArrayMaxX2 = (int)bcArray->getNX2();
   const int bcArrayMaxX3 = (int)bcArray->getNX3();

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
            if(!bcArray->isSolid(x1,x2,x3) && !bcArray->isUndefined(x1,x2,x3))
            {
               int x1p = x1 + 1;
               int x2p = x2 + 1;
               int x3p = x3 + 1;
               //////////////////////////////////////////////////////////////////////////
               //read distribution
               ////////////////////////////////////////////////////////////////////////////
               f[DIR_000] = (*this->zeroDistributions)(x1,x2,x3);

               f[DIR_P00] = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               f[DIR_0P0] = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3);
               f[DIR_00P] = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               f[DIR_PP0] = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               f[DIR_MP0] = (*this->localDistributions)(D3Q27System::ET_NW,x1p,x2,x3);
               f[DIR_P0P] = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               f[DIR_M0P] = (*this->localDistributions)(D3Q27System::ET_TW, x1p,x2,x3);
               f[DIR_0PP] = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               f[DIR_0MP] = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2p,x3);
               f[DIR_PPP] = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               f[DIR_MPP] = (*this->localDistributions)(D3Q27System::ET_TNW,x1p,x2,x3);
               f[DIR_PMP] = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2p,x3);
               f[DIR_MMP] = (*this->localDistributions)(D3Q27System::ET_TSW,x1p,x2p,x3);

               f[DIR_M00] = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1p,x2,x3  );
               f[DIR_0M0] = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2p,x3  );
               f[DIR_00M] = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3p  );
               f[DIR_MM0] = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1p,x2p,x3 );
               f[DIR_PM0] = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2p,x3 );
               f[DIR_M0M] = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1p,x2,x3p );
               f[DIR_P0M] = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3p );
               f[DIR_0MM] = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2p,x3p );
               f[DIR_0PM] = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3p );
               f[DIR_MMM] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1p,x2p,x3p);
               f[DIR_PMM] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2p,x3p);
               f[DIR_MPM] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1p,x2,x3p);
               f[DIR_PPM] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3p);
               //////////////////////////////////////////////////////////////////////////

               drho = f[DIR_000] + f[DIR_P00] + f[DIR_M00] + f[DIR_0P0] + f[DIR_0M0] + f[DIR_00P] + f[DIR_00M]
               + f[DIR_PP0] + f[DIR_MM0] + f[DIR_PM0] + f[DIR_MP0] + f[DIR_P0P] + f[DIR_M0M] + f[DIR_P0M]
               + f[DIR_M0P] + f[DIR_0PP] + f[DIR_0MM] + f[DIR_0PM] + f[DIR_0MP] + f[DIR_PPP] + f[DIR_MMP]
               + f[DIR_PMP] + f[DIR_MPP] + f[DIR_PPM] + f[DIR_MMM] + f[DIR_PMM] + f[DIR_MPM];

               vx1 = f[DIR_P00] - f[DIR_M00] + f[DIR_PP0] - f[DIR_MM0] + f[DIR_PM0] - f[DIR_MP0] + f[DIR_P0P] - f[DIR_M0M]
               + f[DIR_P0M] - f[DIR_M0P] + f[DIR_PPP] - f[DIR_MMP] + f[DIR_PMP] - f[DIR_MPP] + f[DIR_PPM] - f[DIR_MMM]
               + f[DIR_PMM] - f[DIR_MPM]; 

               vx2 = f[DIR_0P0] - f[DIR_0M0] + f[DIR_PP0] - f[DIR_MM0] - f[DIR_PM0] + f[DIR_MP0] + f[DIR_0PP] - f[DIR_0MM] + f[DIR_0PM]
               - f[DIR_0MP] + f[DIR_PPP] - f[DIR_MMP] - f[DIR_PMP] + f[DIR_MPP] + f[DIR_PPM] - f[DIR_MMM] - f[DIR_PMM] 
               + f[DIR_MPM]; 

               vx3 = f[DIR_00P] - f[DIR_00M] + f[DIR_P0P] - f[DIR_M0M] - f[DIR_P0M] + f[DIR_M0P] + f[DIR_0PP] - f[DIR_0MM] - f[DIR_0PM]
               + f[DIR_0MP] + f[DIR_PPP] + f[DIR_MMP] + f[DIR_PMP] + f[DIR_MPP] - f[DIR_PPM] - f[DIR_MMM] - f[DIR_PMM] 
               - f[DIR_MPM];

               real cu_sq= c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

               feq[DIR_000] =  c8o27*(drho-cu_sq);
               feq[DIR_P00] =  c2o27*(drho+c3o1*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq);
               feq[DIR_M00] =  c2o27*(drho+c3o1*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq);
               feq[DIR_0P0] =  c2o27*(drho+c3o1*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq);
               feq[DIR_0M0] =  c2o27*(drho+c3o1*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq);
               feq[DIR_00P] =  c2o27*(drho+c3o1*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq);
               feq[DIR_00M] =  c2o27*(drho+c3o1*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq);
               feq[DIR_PP0] = c1o54*(drho+c3o1*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq);
               feq[DIR_MM0] = c1o54*(drho+c3o1*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
               feq[DIR_PM0] = c1o54*(drho+c3o1*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq);
               feq[DIR_MP0] = c1o54*(drho+c3o1*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
               feq[DIR_P0P] = c1o54*(drho+c3o1*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq);
               feq[DIR_M0M] = c1o54*(drho+c3o1*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
               feq[DIR_P0M] = c1o54*(drho+c3o1*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq);
               feq[DIR_M0P] = c1o54*(drho+c3o1*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
               feq[DIR_0PP] = c1o54*(drho+c3o1*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq);
               feq[DIR_0MM] = c1o54*(drho+c3o1*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
               feq[DIR_0PM] = c1o54*(drho+c3o1*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq);
               feq[DIR_0MP] = c1o54*(drho+c3o1*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
               feq[DIR_PPP]= c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
               feq[DIR_MMM]= c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
               feq[DIR_PPM]= c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
               feq[DIR_MMP]= c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
               feq[DIR_PMP]= c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
               feq[DIR_MPM]= c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
               feq[DIR_PMM]= c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
               feq[DIR_MPP]= c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

               //Relaxation
               f[DIR_000] += (feq[DIR_000]-f[DIR_000])*collFactor;
               f[DIR_P00] += (feq[DIR_P00]-f[DIR_P00])*collFactor;
               f[DIR_M00] += (feq[DIR_M00]-f[DIR_M00])*collFactor;
               f[DIR_0P0] += (feq[DIR_0P0]-f[DIR_0P0])*collFactor;
               f[DIR_0M0] += (feq[DIR_0M0]-f[DIR_0M0])*collFactor;
               f[DIR_00P] += (feq[DIR_00P]-f[DIR_00P])*collFactor;
               f[DIR_00M] += (feq[DIR_00M]-f[DIR_00M])*collFactor;
               f[DIR_PP0] += (feq[DIR_PP0]-f[DIR_PP0])*collFactor;
               f[DIR_MM0] += (feq[DIR_MM0]-f[DIR_MM0])*collFactor;
               f[DIR_PM0] += (feq[DIR_PM0]-f[DIR_PM0])*collFactor;
               f[DIR_MP0] += (feq[DIR_MP0]-f[DIR_MP0])*collFactor;
               f[DIR_P0P] += (feq[DIR_P0P]-f[DIR_P0P])*collFactor;
               f[DIR_M0M] += (feq[DIR_M0M]-f[DIR_M0M])*collFactor;
               f[DIR_P0M] += (feq[DIR_P0M]-f[DIR_P0M])*collFactor;
               f[DIR_M0P] += (feq[DIR_M0P]-f[DIR_M0P])*collFactor;
               f[DIR_0PP] += (feq[DIR_0PP]-f[DIR_0PP])*collFactor;
               f[DIR_0MM] += (feq[DIR_0MM]-f[DIR_0MM])*collFactor;
               f[DIR_0PM] += (feq[DIR_0PM]-f[DIR_0PM])*collFactor;
               f[DIR_0MP] += (feq[DIR_0MP]-f[DIR_0MP])*collFactor;

               f[DIR_PPP] += (feq[DIR_PPP]-f[DIR_PPP])*collFactor;
               f[DIR_MMM] += (feq[DIR_MMM]-f[DIR_MMM])*collFactor;
               f[DIR_PPM] += (feq[DIR_PPM]-f[DIR_PPM])*collFactor;
               f[DIR_MMP] += (feq[DIR_MMP]-f[DIR_MMP])*collFactor;
               f[DIR_PMP] += (feq[DIR_PMP]-f[DIR_PMP])*collFactor;
               f[DIR_MPM] += (feq[DIR_MPM]-f[DIR_MPM])*collFactor;
               f[DIR_PMM] += (feq[DIR_PMM]-f[DIR_PMM])*collFactor;
               f[DIR_MPP] += (feq[DIR_MPP]-f[DIR_MPP])*collFactor;

               //////////////////////////////////////////////////////////////////////////
               //forcing
               if (withForcing)
               {
                  muX1 = x1+ix1*bcArrayMaxX1;
                  muX2 = x2+ix2*bcArrayMaxX2;
                  muX3 = x3+ix3*bcArrayMaxX3;

                  forcingX1 = muForcingX1.Eval();
                  forcingX2 = muForcingX2.Eval();
                  forcingX3 = muForcingX3.Eval();

                  f[DIR_000] += c0o1;
                  f[DIR_P00] +=  c3o1*c2o27  *  (forcingX1)                    ;
                  f[DIR_M00] +=  c3o1*c2o27  *  (-forcingX1)                   ;
                  f[DIR_0P0] +=  c3o1*c2o27  *             (forcingX2)         ;
                  f[DIR_0M0] +=  c3o1*c2o27  *             (-forcingX2)        ;
                  f[DIR_00P] +=  c3o1*c2o27  *                     (forcingX3) ;
                  f[DIR_00M] +=  c3o1*c2o27  *                     (-forcingX3);
                  f[DIR_PP0] +=  c3o1*c1o54 * ( forcingX1+forcingX2          ) ;
                  f[DIR_MM0 ] +=  c3o1*c1o54 * (-forcingX1-forcingX2          ) ;
                  f[DIR_PM0 ] +=  c3o1*c1o54 * ( forcingX1-forcingX2          ) ;
                  f[DIR_MP0 ] +=  c3o1*c1o54 * (-forcingX1+forcingX2          ) ;
                  f[DIR_P0P ] +=  c3o1*c1o54 * ( forcingX1          +forcingX3) ;
                  f[DIR_M0M ] +=  c3o1*c1o54 * (-forcingX1          -forcingX3) ;
                  f[DIR_P0M ] +=  c3o1*c1o54 * ( forcingX1          -forcingX3) ;
                  f[DIR_M0P ] +=  c3o1*c1o54 * (-forcingX1          +forcingX3) ;
                  f[DIR_0PP ] +=  c3o1*c1o54 * (           forcingX2+forcingX3) ;
                  f[DIR_0MM ] +=  c3o1*c1o54 * (          -forcingX2-forcingX3) ;
                  f[DIR_0PM ] +=  c3o1*c1o54 * (           forcingX2-forcingX3) ;
                  f[DIR_0MP ] +=  c3o1*c1o54 * (          -forcingX2+forcingX3) ;
                  f[DIR_PPP] +=  c3o1*c1o216* ( forcingX1+forcingX2+forcingX3) ;
                  f[DIR_MMM] +=  c3o1*c1o216* (-forcingX1-forcingX2-forcingX3) ;
                  f[DIR_PPM] +=  c3o1*c1o216* ( forcingX1+forcingX2-forcingX3) ;
                  f[DIR_MMP] +=  c3o1*c1o216* (-forcingX1-forcingX2+forcingX3) ;
                  f[DIR_PMP] +=  c3o1*c1o216* ( forcingX1-forcingX2+forcingX3) ;
                  f[DIR_MPM] +=  c3o1*c1o216* (-forcingX1+forcingX2-forcingX3) ;
                  f[DIR_PMM] +=  c3o1*c1o216* ( forcingX1-forcingX2-forcingX3) ;
                  f[DIR_MPP] +=  c3o1*c1o216* (-forcingX1+forcingX2+forcingX3) ;
               }
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               real rho_post = f[REST] + f[DIR_P00] + f[W] + f[N] + f[S] + f[T] + f[B] 
               + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
               + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
               + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];
               real dif = drho - rho_post;
#ifdef SINGLEPRECISION
               if(dif > 10.0E-7 || dif < -10.0E-7)
#else
               if(dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
                  UB_THROW(UbException(UB_EXARGS,"rho is not correct"));
               }
#endif
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[INV_P00];
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[INV_0P0];
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[INV_00P];
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[INV_PP0];
               (*this->localDistributions)(D3Q27System::ET_NW,x1p,x2,  x3) = f[INV_MP0];
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[INV_P0P];
               (*this->localDistributions)(D3Q27System::ET_TW,x1p,x2,  x3) = f[INV_M0P];
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[INV_0PP];
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2p,x3) = f[INV_0MP];
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[INV_PPP];
               (*this->localDistributions)(D3Q27System::ET_TNW,x1p,x2,  x3) = f[INV_MPP];
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2p,x3) = f[INV_PMP];
               (*this->localDistributions)(D3Q27System::ET_TSW,x1p,x2p,x3) = f[INV_MMP];

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1p,x2,  x3    ) = f[INV_M00 ];
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2p,x3    ) = f[INV_0M0 ];
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3p  ) = f[INV_00M ];
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1p,x2p,x3   ) = f[INV_MM0];
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2p,x3   ) = f[INV_PM0];
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1p,x2,  x3p ) = f[INV_M0M];
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3p ) = f[INV_P0M];
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2p,x3p ) = f[INV_0MM];
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3p ) = f[INV_0PM];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1p,x2p,x3p) = f[INV_MMM];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2p,x3p) = f[INV_PMM];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1p,x2,  x3p) = f[INV_MPM];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3p) = f[INV_PPM];

               (*this->zeroDistributions)(x1,x2,x3) = f[DIR_000];
               //////////////////////////////////////////////////////////////////////////


            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
real LBMKernelETD3Q27BGK::getCalculationTime()
{
   return c0o1;
}
