#include "LBMKernelETD3Q27BGK.h"
#include "D3Q27System.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "D3Q27EsoTwist3DSoA.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "BCArray3D.h"

//#define PROOF_CORRECTNESS

using namespace UbMath;

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
   kernel->setBCProcessor(bcProcessor->clone(kernel));
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

   //initializing of forcing stuff 
   if (withForcing)
   {
      muForcingX1.DefineVar("x1",&muX1); muForcingX1.DefineVar("x2",&muX2); muForcingX1.DefineVar("x3",&muX3);
      muForcingX2.DefineVar("x1",&muX1); muForcingX2.DefineVar("x2",&muX2); muForcingX2.DefineVar("x3",&muX3);
      muForcingX3.DefineVar("x1",&muX1); muForcingX3.DefineVar("x2",&muX2); muForcingX3.DefineVar("x3",&muX3);
      forcingX1 = 0;
      forcingX2 = 0;
      forcingX3 = 0;
   }
   /////////////////////////////////////

   localDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = std::dynamic_pointer_cast<D3Q27EsoTwist3DSplittedVector>(dataSet->getFdistributions())->getZeroDistributions();

   SPtr<BCArray3D> bcArray = this->getBCProcessor()->getBCArray();
   LBMReal f[D3Q27System::ENDF+1];
   LBMReal feq[D3Q27System::ENDF+1];
   LBMReal drho,vx1,vx2,vx3;
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

               f[E] = (*this->localDistributions)(D3Q27System::ET_E, x1,x2,x3);
               f[N] = (*this->localDistributions)(D3Q27System::ET_N,x1,x2,x3);  
               f[T] = (*this->localDistributions)(D3Q27System::ET_T,x1,x2,x3);
               f[NE] = (*this->localDistributions)(D3Q27System::ET_NE,x1,x2,x3);
               f[NW] = (*this->localDistributions)(D3Q27System::ET_NW,x1p,x2,x3);
               f[TE] = (*this->localDistributions)(D3Q27System::ET_TE,x1,x2,x3);
               f[TW] = (*this->localDistributions)(D3Q27System::ET_TW, x1p,x2,x3);
               f[TN] = (*this->localDistributions)(D3Q27System::ET_TN,x1,x2,x3);
               f[TS] = (*this->localDistributions)(D3Q27System::ET_TS,x1,x2p,x3);
               f[TNE] = (*this->localDistributions)(D3Q27System::ET_TNE,x1,x2,x3);
               f[TNW] = (*this->localDistributions)(D3Q27System::ET_TNW,x1p,x2,x3);
               f[TSE] = (*this->localDistributions)(D3Q27System::ET_TSE,x1,x2p,x3);
               f[TSW] = (*this->localDistributions)(D3Q27System::ET_TSW,x1p,x2p,x3);

               f[W ] = (*this->nonLocalDistributions)(D3Q27System::ET_W,x1p,x2,x3  );
               f[S ] = (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,x2p,x3  );
               f[B ] = (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,x2,x3p  );
               f[SW] = (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1p,x2p,x3 );
               f[SE] = (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,x2p,x3 );
               f[BW] = (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1p,x2,x3p );
               f[BE] = (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,x2,x3p );
               f[BS] = (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,x2p,x3p );
               f[BN] = (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,x2,x3p );
               f[BSW] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1p,x2p,x3p);
               f[BSE] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,x2p,x3p);
               f[BNW] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1p,x2,x3p);
               f[BNE] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,x2,x3p);
               //////////////////////////////////////////////////////////////////////////

               drho = f[DIR_000] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
               + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
               + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
               + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];

               vx1 = f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
               + f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
               + f[BSE] - f[BNW]; 

               vx2 = f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
               - f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
               + f[BNW]; 

               vx3 = f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
               + f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
               - f[BNW];

               LBMReal cu_sq=1.5*(vx1*vx1+vx2*vx2+vx3*vx3);

               feq[DIR_000] =  c8o27*(drho-cu_sq);
               feq[E] =  c2o27*(drho+3.0*( vx1   )+c9o2*( vx1   )*( vx1   )-cu_sq);
               feq[W] =  c2o27*(drho+3.0*(-vx1   )+c9o2*(-vx1   )*(-vx1   )-cu_sq);
               feq[N] =  c2o27*(drho+3.0*(    vx2)+c9o2*(    vx2)*(    vx2)-cu_sq);
               feq[S] =  c2o27*(drho+3.0*(   -vx2)+c9o2*(   -vx2)*(   -vx2)-cu_sq);
               feq[T] =  c2o27*(drho+3.0*( vx3   )+c9o2*(    vx3)*(    vx3)-cu_sq);
               feq[B] =  c2o27*(drho+3.0*(   -vx3)+c9o2*(   -vx3)*(   -vx3)-cu_sq);
               feq[NE] = c1o54*(drho+3.0*( vx1+vx2)+c9o2*( vx1+vx2)*( vx1+vx2)-cu_sq);
               feq[SW] = c1o54*(drho+3.0*(-vx1-vx2)+c9o2*(-vx1-vx2)*(-vx1-vx2)-cu_sq);
               feq[SE] = c1o54*(drho+3.0*( vx1-vx2)+c9o2*( vx1-vx2)*( vx1-vx2)-cu_sq);
               feq[NW] = c1o54*(drho+3.0*(-vx1+vx2)+c9o2*(-vx1+vx2)*(-vx1+vx2)-cu_sq);
               feq[TE] = c1o54*(drho+3.0*( vx1+vx3)+c9o2*( vx1+vx3)*( vx1+vx3)-cu_sq);
               feq[BW] = c1o54*(drho+3.0*(-vx1-vx3)+c9o2*(-vx1-vx3)*(-vx1-vx3)-cu_sq);
               feq[BE] = c1o54*(drho+3.0*( vx1-vx3)+c9o2*( vx1-vx3)*( vx1-vx3)-cu_sq);
               feq[TW] = c1o54*(drho+3.0*(-vx1+vx3)+c9o2*(-vx1+vx3)*(-vx1+vx3)-cu_sq);
               feq[TN] = c1o54*(drho+3.0*( vx2+vx3)+c9o2*( vx2+vx3)*( vx2+vx3)-cu_sq);
               feq[BS] = c1o54*(drho+3.0*(-vx2-vx3)+c9o2*(-vx2-vx3)*(-vx2-vx3)-cu_sq);
               feq[BN] = c1o54*(drho+3.0*( vx2-vx3)+c9o2*( vx2-vx3)*( vx2-vx3)-cu_sq);
               feq[TS] = c1o54*(drho+3.0*(-vx2+vx3)+c9o2*(-vx2+vx3)*(-vx2+vx3)-cu_sq);
               feq[TNE]= c1o216*(drho+3.0*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
               feq[BSW]= c1o216*(drho+3.0*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
               feq[BNE]= c1o216*(drho+3.0*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
               feq[TSW]= c1o216*(drho+3.0*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
               feq[TSE]= c1o216*(drho+3.0*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
               feq[BNW]= c1o216*(drho+3.0*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
               feq[BSE]= c1o216*(drho+3.0*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
               feq[TNW]= c1o216*(drho+3.0*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);

               //Relaxation
               f[DIR_000] += (feq[DIR_000]-f[DIR_000])*collFactor;
               f[E] += (feq[E]-f[E])*collFactor;
               f[W] += (feq[W]-f[W])*collFactor;
               f[N] += (feq[N]-f[N])*collFactor;
               f[S] += (feq[S]-f[S])*collFactor;
               f[T] += (feq[T]-f[T])*collFactor;
               f[B] += (feq[B]-f[B])*collFactor;
               f[NE] += (feq[NE]-f[NE])*collFactor;
               f[SW] += (feq[SW]-f[SW])*collFactor;
               f[SE] += (feq[SE]-f[SE])*collFactor;
               f[NW] += (feq[NW]-f[NW])*collFactor;
               f[TE] += (feq[TE]-f[TE])*collFactor;
               f[BW] += (feq[BW]-f[BW])*collFactor;
               f[BE] += (feq[BE]-f[BE])*collFactor;
               f[TW] += (feq[TW]-f[TW])*collFactor;
               f[TN] += (feq[TN]-f[TN])*collFactor;
               f[BS] += (feq[BS]-f[BS])*collFactor;
               f[BN] += (feq[BN]-f[BN])*collFactor;
               f[TS] += (feq[TS]-f[TS])*collFactor;

               f[TNE] += (feq[TNE]-f[TNE])*collFactor;
               f[BSW] += (feq[BSW]-f[BSW])*collFactor;
               f[BNE] += (feq[BNE]-f[BNE])*collFactor;
               f[TSW] += (feq[TSW]-f[TSW])*collFactor;
               f[TSE] += (feq[TSE]-f[TSE])*collFactor;
               f[BNW] += (feq[BNW]-f[BNW])*collFactor;
               f[BSE] += (feq[BSE]-f[BSE])*collFactor;
               f[TNW] += (feq[TNW]-f[TNW])*collFactor;

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

                  f[DIR_000] +=                   0.0                        ;
                  f[E  ] +=  3.0*c2o27  *  (forcingX1)                    ;
                  f[W  ] +=  3.0*c2o27  *  (-forcingX1)                   ;
                  f[N  ] +=  3.0*c2o27  *             (forcingX2)         ;
                  f[S  ] +=  3.0*c2o27  *             (-forcingX2)        ;
                  f[T  ] +=  3.0*c2o27  *                     (forcingX3) ;
                  f[B  ] +=  3.0*c2o27  *                     (-forcingX3);
                  f[NE ] +=  3.0*c1o54 * ( forcingX1+forcingX2          ) ;
                  f[SW ] +=  3.0*c1o54 * (-forcingX1-forcingX2          ) ;
                  f[SE ] +=  3.0*c1o54 * ( forcingX1-forcingX2          ) ;
                  f[NW ] +=  3.0*c1o54 * (-forcingX1+forcingX2          ) ;
                  f[TE ] +=  3.0*c1o54 * ( forcingX1          +forcingX3) ;
                  f[BW ] +=  3.0*c1o54 * (-forcingX1          -forcingX3) ;
                  f[BE ] +=  3.0*c1o54 * ( forcingX1          -forcingX3) ;
                  f[TW ] +=  3.0*c1o54 * (-forcingX1          +forcingX3) ;
                  f[TN ] +=  3.0*c1o54 * (           forcingX2+forcingX3) ;
                  f[BS ] +=  3.0*c1o54 * (          -forcingX2-forcingX3) ;
                  f[BN ] +=  3.0*c1o54 * (           forcingX2-forcingX3) ;
                  f[TS ] +=  3.0*c1o54 * (          -forcingX2+forcingX3) ;
                  f[TNE] +=  3.0*c1o216* ( forcingX1+forcingX2+forcingX3) ;
                  f[BSW] +=  3.0*c1o216* (-forcingX1-forcingX2-forcingX3) ;
                  f[BNE] +=  3.0*c1o216* ( forcingX1+forcingX2-forcingX3) ;
                  f[TSW] +=  3.0*c1o216* (-forcingX1-forcingX2+forcingX3) ;
                  f[TSE] +=  3.0*c1o216* ( forcingX1-forcingX2+forcingX3) ;
                  f[BNW] +=  3.0*c1o216* (-forcingX1+forcingX2-forcingX3) ;
                  f[BSE] +=  3.0*c1o216* ( forcingX1-forcingX2-forcingX3) ;
                  f[TNW] +=  3.0*c1o216* (-forcingX1+forcingX2+forcingX3) ;
               }
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               LBMReal rho_post = f[REST] + f[E] + f[W] + f[N] + f[S] + f[T] + f[B] 
               + f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE]
               + f[TW] + f[TN] + f[BS] + f[BN] + f[TS] + f[TNE] + f[TSW]
               + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW];
               LBMReal dif = drho - rho_post;
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
               (*this->localDistributions)(D3Q27System::ET_E,x1,  x2,  x3) = f[D3Q27System::INV_E];
               (*this->localDistributions)(D3Q27System::ET_N,x1,  x2,  x3) = f[D3Q27System::INV_N];
               (*this->localDistributions)(D3Q27System::ET_T,x1,  x2,  x3) = f[D3Q27System::INV_T];
               (*this->localDistributions)(D3Q27System::ET_NE,x1,  x2,  x3) = f[D3Q27System::INV_NE];
               (*this->localDistributions)(D3Q27System::ET_NW,x1p,x2,  x3) = f[D3Q27System::INV_NW];
               (*this->localDistributions)(D3Q27System::ET_TE,x1,  x2,  x3) = f[D3Q27System::INV_TE];
               (*this->localDistributions)(D3Q27System::ET_TW,x1p,x2,  x3) = f[D3Q27System::INV_TW];
               (*this->localDistributions)(D3Q27System::ET_TN,x1,  x2,  x3) = f[D3Q27System::INV_TN];
               (*this->localDistributions)(D3Q27System::ET_TS,x1,  x2p,x3) = f[D3Q27System::INV_TS];
               (*this->localDistributions)(D3Q27System::ET_TNE,x1,  x2,  x3) = f[D3Q27System::INV_TNE];
               (*this->localDistributions)(D3Q27System::ET_TNW,x1p,x2,  x3) = f[D3Q27System::INV_TNW];
               (*this->localDistributions)(D3Q27System::ET_TSE,x1,  x2p,x3) = f[D3Q27System::INV_TSE];
               (*this->localDistributions)(D3Q27System::ET_TSW,x1p,x2p,x3) = f[D3Q27System::INV_TSW];

               (*this->nonLocalDistributions)(D3Q27System::ET_W,x1p,x2,  x3    ) = f[D3Q27System::INV_W ];
               (*this->nonLocalDistributions)(D3Q27System::ET_S,x1,  x2p,x3    ) = f[D3Q27System::INV_S ];
               (*this->nonLocalDistributions)(D3Q27System::ET_B,x1,  x2,  x3p  ) = f[D3Q27System::INV_B ];
               (*this->nonLocalDistributions)(D3Q27System::ET_SW,x1p,x2p,x3   ) = f[D3Q27System::INV_SW];
               (*this->nonLocalDistributions)(D3Q27System::ET_SE,x1,  x2p,x3   ) = f[D3Q27System::INV_SE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BW,x1p,x2,  x3p ) = f[D3Q27System::INV_BW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BE,x1,  x2,  x3p ) = f[D3Q27System::INV_BE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BS,x1,  x2p,x3p ) = f[D3Q27System::INV_BS];
               (*this->nonLocalDistributions)(D3Q27System::ET_BN,x1,  x2,  x3p ) = f[D3Q27System::INV_BN];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSW,x1p,x2p,x3p) = f[D3Q27System::INV_BSW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BSE,x1,  x2p,x3p) = f[D3Q27System::INV_BSE];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNW,x1p,x2,  x3p) = f[D3Q27System::INV_BNW];
               (*this->nonLocalDistributions)(D3Q27System::ET_BNE,x1,  x2,  x3p) = f[D3Q27System::INV_BNE];

               (*this->zeroDistributions)(x1,x2,x3) = f[D3Q27System::DIR_000];
               //////////////////////////////////////////////////////////////////////////


            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
double LBMKernelETD3Q27BGK::getCalculationTime()
{
   return 0.0;
}
