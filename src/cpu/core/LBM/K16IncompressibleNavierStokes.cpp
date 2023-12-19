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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#include "K16IncompressibleNavierStokes.h"
#include "D3Q27System.h"
#include "Interpolator.h"
#include "EsoSplit.h"
#include "DataSet3D.h"
#include <cmath>
#include "Block3D.h"

#define PROOF_CORRECTNESS

//using namespace UbMath;
using namespace vf::basics::constant;

//////////////////////////////////////////////////////////////////////////
K16IncompressibleNavierStokes::K16IncompressibleNavierStokes()
{
   this->parameter = NORMAL;
   this->OxyyMxzz = c1o1;
   this->compressible = false;
}
//////////////////////////////////////////////////////////////////////////
K16IncompressibleNavierStokes::~K16IncompressibleNavierStokes(void)
= default;
//////////////////////////////////////////////////////////////////////////
void K16IncompressibleNavierStokes::initDataSet()
{
   SPtr<DistributionArray3D> d(new EsoSplit(nx[0]+2, nx[1]+2, nx[2]+2, -999.9));
   dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> K16IncompressibleNavierStokes::clone()
{
   SPtr<K16IncompressibleNavierStokes> kernel(new K16IncompressibleNavierStokes());
   kernel->setNX(nx);
   kernel->initDataSet();
   kernel->setCollisionFactor(this->collFactor);
   kernel->setBCSet(bcSet->clone(kernel));
   kernel->setWithForcing(withForcing);
   kernel->setForcingX1(muForcingX1);
   kernel->setForcingX2(muForcingX2);
   kernel->setForcingX3(muForcingX3);
   kernel->setIndex(ix1, ix2, ix3);
   kernel->setDeltaT(deltaT);
   kernel->setBlock(block.lock());

   switch (parameter)
   {
   case NORMAL:
      kernel->OxyyMxzz = c1o1;
      break;
   case MAGIC:
      kernel->OxyyMxzz = c2o1 +(-collFactor);
      break;
   }
   return kernel;
}
//////////////////////////////////////////////////////////////////////////
void K16IncompressibleNavierStokes::calculate(int step)
{
   using namespace D3Q27System;
   using namespace std;
   using namespace vf::lbm::dir;

   //timer.start();

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

      muNu = (c1o1/c3o1)*(c1o1/collFactor - c1o1/c2o1);

      muForcingX1.DefineVar("nu", &muNu);
      muForcingX2.DefineVar("nu", &muNu);
      muForcingX3.DefineVar("nu", &muNu);

//      real forcingX1 = 0;
//      real forcingX2 = 0;
//      real forcingX3 = 0;
   }
   /////////////////////////////////////

   localDistributions = dynamicPointerCast<EsoSplit>(dataSet->getFdistributions())->getLocalDistributions();
   nonLocalDistributions = dynamicPointerCast<EsoSplit>(dataSet->getFdistributions())->getNonLocalDistributions();
   zeroDistributions = dynamicPointerCast<EsoSplit>(dataSet->getFdistributions())->getZeroDistributions();

   SPtr<BCArray3D> bcArray = this->getBCSet()->getBCArray();

   const int bcArrayMaxX1 = (int)bcArray->getNX1();
   const int bcArrayMaxX2 = (int)bcArray->getNX2();
   const int bcArrayMaxX3 = (int)bcArray->getNX3();

   int minX1 = ghostLayerWidth;
   int minX2 = ghostLayerWidth;
   int minX3 = ghostLayerWidth;
   int maxX1 = bcArrayMaxX1-ghostLayerWidth;
   int maxX2 = bcArrayMaxX2-ghostLayerWidth;
   int maxX3 = bcArrayMaxX3-ghostLayerWidth;

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

               real mfcbb = (*this->localDistributions)(eP00, x1, x2, x3);
               real mfbcb = (*this->localDistributions)(e0P0, x1, x2, x3);
               real mfbbc = (*this->localDistributions)(e00P, x1, x2, x3);
               real mfccb = (*this->localDistributions)(ePP0, x1, x2, x3);
               real mfacb = (*this->localDistributions)(eMP0, x1p, x2, x3);
               real mfcbc = (*this->localDistributions)(eP0P, x1, x2, x3);
               real mfabc = (*this->localDistributions)(eM0P, x1p, x2, x3);
               real mfbcc = (*this->localDistributions)(e0PP, x1, x2, x3);
               real mfbac = (*this->localDistributions)(e0MP, x1, x2p, x3);
               real mfccc = (*this->localDistributions)(ePPP, x1, x2, x3);
               real mfacc = (*this->localDistributions)(eMPP, x1p, x2, x3);
               real mfcac = (*this->localDistributions)(ePMP, x1, x2p, x3);
               real mfaac = (*this->localDistributions)(eMMP, x1p, x2p, x3);

               real mfabb = (*this->nonLocalDistributions)(eM00, x1p, x2, x3);
               real mfbab = (*this->nonLocalDistributions)(e0M0, x1, x2p, x3);
               real mfbba = (*this->nonLocalDistributions)(e00M, x1, x2, x3p);
               real mfaab = (*this->nonLocalDistributions)(eMM0, x1p, x2p, x3);
               real mfcab = (*this->nonLocalDistributions)(ePM0, x1, x2p, x3);
               real mfaba = (*this->nonLocalDistributions)(eM0M, x1p, x2, x3p);
               real mfcba = (*this->nonLocalDistributions)(eP0M, x1, x2, x3p);
               real mfbaa = (*this->nonLocalDistributions)(e0MM, x1, x2p, x3p);
               real mfbca = (*this->nonLocalDistributions)(e0PM, x1, x2, x3p);
               real mfaaa = (*this->nonLocalDistributions)(eMMM, x1p, x2p, x3p);
               real mfcaa = (*this->nonLocalDistributions)(ePMM, x1, x2p, x3p);
               real mfaca = (*this->nonLocalDistributions)(eMPM, x1p, x2, x3p);
               real mfcca = (*this->nonLocalDistributions)(ePPM, x1, x2, x3p);

               real mfbbb = (*this->zeroDistributions)(x1, x2, x3);

               real m0, m1, m2;

               real rho=(mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;

               real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
                  (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                  (mfcbb-mfabb));
               real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
                  (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                  (mfbcb-mfbab));
               real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
                  (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                  (mfbbc-mfbba));

               //forcing 
               ///////////////////////////////////////////////////////////////////////////////////////////
               if (withForcing)
               {
                  muX1 = static_cast<real>(x1-1+ix1*maxX1);
                  muX2 = static_cast<real>(x2-1+ix2*maxX2);
                  muX3 = static_cast<real>(x3-1+ix3*maxX3);

                  forcingX1 = muForcingX1.Eval();
                  forcingX2 = muForcingX2.Eval();
                  forcingX3 = muForcingX3.Eval();

                  vvx += forcingX1*deltaT*c1o2; // X
                  vvy += forcingX2*deltaT*c1o2; // Y
                  vvz += forcingX3*deltaT*c1o2; // Z
               }
               ///////////////////////////////////////////////////////////////////////////////////////////               
               real oMdrho;

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
               m0+=mfbbb; 
               oMdrho = c1o1 - (oMdrho + m0);

               real vx2;
               real vy2;
               real vz2;
               vx2=vvx*vvx;
               vy2=vvy*vvy;
               vz2=vvz*vvz;
               ////////////////////////////////////////////////////////////////////////////////////
               real wadjust;
               real qudricLimit = c1o100;
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
               mfaac = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfabc;
               m1    = mfabc  - mfaba;
               m0    = m2          + mfabb;
               mfaba = m0;
               m0   += c1o9 * oMdrho;
               mfabb = m1 -        m0 * vvz;
               mfabc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfacc;
               m1    = mfacc  - mfaca;
               m0    = m2          + mfacb;
               mfaca = m0;
               m0   += c1o36 * oMdrho;
               mfacb = m1 -        m0 * vvz;
               mfacc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbac;
               m1    = mfbac - mfbaa;
               m0    = m2          + mfbab;
               mfbaa = m0;
               m0   += c1o9 * oMdrho;
               mfbab = m1 -        m0 * vvz;
               mfbac = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbba  + mfbbc;
               m1    = mfbbc  - mfbba;
               m0    = m2          + mfbbb;
               mfbba = m0;
               m0   += c4o9 * oMdrho;
               mfbbb = m1 -        m0 * vvz;
               mfbbc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbca  + mfbcc;
               m1    = mfbcc  - mfbca;
               m0    = m2          + mfbcb;
               mfbca = m0;
               m0   += c1o9 * oMdrho;
               mfbcb = m1 -        m0 * vvz;
               mfbcc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcac;
               m1    = mfcac - mfcaa;
               m0    = m2          + mfcab;
               mfcaa = m0;
               m0   += c1o36 * oMdrho;
               mfcab = m1 -        m0 * vvz;
               mfcac = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcba  + mfcbc;
               m1    = mfcbc  - mfcba;
               m0    = m2          + mfcbb;
               mfcba = m0;
               m0   += c1o9 * oMdrho;
               mfcbb = m1 -        m0 * vvz;
               mfcbc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcca  + mfccc;
               m1    = mfccc  - mfcca;
               m0    = m2          + mfccb;
               mfcca = m0;
               m0   += c1o36 * oMdrho;
               mfccb = m1 -        m0 * vvz;
               mfccc = m2 - c2o1 *   m1 * vvz + vz2 * m0;
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
               mfaca = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab  + mfacb;
               m1    = mfacb  - mfaab;
               m0    = m2          + mfabb;
               mfaab = m0;
               mfabb = m1 -        m0 * vvy;
               mfacb = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac  + mfacc;
               m1    = mfacc  - mfaac;
               m0    = m2          + mfabc;
               mfaac = m0;
               m0   += c1o18 * oMdrho;
               mfabc = m1 -        m0 * vvy;
               mfacc = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbaa + mfbca;
               m1    = mfbca - mfbaa;
               m0    = m2          + mfbba;
               mfbaa = m0;
               m0   += c2o3 * oMdrho;
               mfbba = m1 -        m0 * vvy;
               mfbca = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbab  + mfbcb;
               m1    = mfbcb  - mfbab;
               m0    = m2          + mfbbb;
               mfbab = m0;
               mfbbb = m1 -        m0 * vvy;
               mfbcb = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfbac  + mfbcc;
               m1    = mfbcc  - mfbac;
               m0    = m2          + mfbbc;
               mfbac = m0;
               m0   += c2o9 * oMdrho;
               mfbbc = m1 -        m0 * vvy;
               mfbcc = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcaa + mfcca;
               m1    = mfcca - mfcaa;
               m0    = m2          + mfcba;
               mfcaa = m0;
               m0   += c1o6 * oMdrho;
               mfcba = m1 -        m0 * vvy;
               mfcca = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcab  + mfccb;
               m1    = mfccb  - mfcab;
               m0    = m2          + mfcbb;
               mfcab = m0;
               mfcbb = m1 -        m0 * vvy;
               mfccb = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfcac  + mfccc;
               m1    = mfccc  - mfcac;
               m0    = m2          + mfcbc;
               mfcac = m0;
               m0   += c1o18 * oMdrho;
               mfcbc = m1 -        m0 * vvy;
               mfccc = m2 - c2o1 *   m1 * vvy + vy2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m2    = mfaaa + mfcaa;
               m1    = mfcaa - mfaaa;
               m0    = m2          + mfbaa;
               mfaaa = m0;
               m0   += c1o1 * oMdrho;
               mfbaa = m1 -        m0 * vvx;
               mfcaa = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaba  + mfcba;
               m1    = mfcba  - mfaba;
               m0    = m2          + mfbba;
               mfaba = m0;
               mfbba = m1 -        m0 * vvx;
               mfcba = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaca  + mfcca;
               m1    = mfcca  - mfaca;
               m0    = m2          + mfbca;
               mfaca = m0;
               m0   += c1o3 * oMdrho;
               mfbca = m1 -        m0 * vvx;
               mfcca = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaab + mfcab;
               m1    = mfcab - mfaab;
               m0    = m2          + mfbab;
               mfaab = m0;
               mfbab = m1 -        m0 * vvx;
               mfcab = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabb  + mfcbb;
               m1    = mfcbb  - mfabb;
               m0    = m2          + mfbbb;
               mfabb = m0;
               mfbbb = m1 -        m0 * vvx;
               mfcbb = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacb  + mfccb;
               m1    = mfccb  - mfacb;
               m0    = m2          + mfbcb;
               mfacb = m0;
               mfbcb = m1 -        m0 * vvx;
               mfccb = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfaac + mfcac;
               m1    = mfcac - mfaac;
               m0    = m2          + mfbac;
               mfaac = m0;
               m0   += c1o3 * oMdrho;
               mfbac = m1 -        m0 * vvx;
               mfcac = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfabc  + mfcbc;
               m1    = mfcbc  - mfabc;
               m0    = m2          + mfbbc;
               mfabc = m0;
               mfbbc = m1 -        m0 * vvx;
               mfcbc = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               m2    = mfacc  + mfccc;
               m1    = mfccc  - mfacc;
               m0    = m2          + mfbcc;
               mfacc = m0;
               m0   += c1o9 * oMdrho;
               mfbcc = m1 -        m0 * vvx;
               mfccc = m2 - c2o1 *   m1 * vvx + vx2 * m0;
               ////////////////////////////////////////////////////////////////////////////////////
               // Cumulants
               ////////////////////////////////////////////////////////////////////////////////////
               real OxxPyyPzz = c1o1; //omega2 or bulk viscosity
               real OxyyPxzz  = c1o1;//-s9;//2+s9;//
               //real OxyyMxzz  = c1o1;//2+s9;//
               real O4        = c1o1;
               real O5        = c1o1;
               real O6        = c1o1;

               //Cum 4.
               //real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + 2. * mfbba * mfbab); // till 18.05.2015
               //real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + 2. * mfbba * mfabb); // till 18.05.2015
               //real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + 2. * mfbab * mfabb); // till 18.05.2015

               real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab);
               real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb);
               real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb);

               real CUMcca = mfcca - ((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
               real CUMcac = mfcac - ((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
               real CUMacc = mfacc - ((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);

               //Cum 5.
               real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
               real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
               real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

               //Cum 6.
               real CUMccc = mfccc  +((-c4o1 *  mfbbb * mfbbb
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                  -  c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
                  +(c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                     +  c2o1 * (mfcaa * mfaca * mfaac)
                     + c16o1 *  mfbba * mfbab * mfabb)
                  - c1o3* (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
                  - c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
                  +(c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                     +       (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;

               //2.
               // linear combinations
               real mxxPyyPzz = mfcaa + mfaca + mfaac;
               real mxxMyy    = mfcaa - mfaca;
               real mxxMzz         = mfcaa - mfaac;

               real dxux = -c1o2 * collFactor *(mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz*(mfaaa - mxxPyyPzz);
               real dyuy = dxux + collFactor * c3o2 * mxxMyy;
               real dzuz = dxux + collFactor * c3o2 * mxxMzz;

               //relax
               mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
               mxxMyy    += collFactor * (-mxxMyy) - c3o1 * (c1o1 - c1o2 * collFactor) * (vx2 * dxux - vy2 * dyuy);
               mxxMzz    += collFactor * (-mxxMzz) - c3o1 * (c1o1 - c1o2 * collFactor) * (vx2 * dxux - vz2 * dzuz);

               mfabb     += collFactor * (-mfabb);
               mfbab     += collFactor * (-mfbab);
               mfbba     += collFactor * (-mfbba);

               // linear combinations back
               mfcaa = c1o3 * (mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaca = c1o3 * (-c2o1 *  mxxMyy +      mxxMzz + mxxPyyPzz);
               mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

               //3.
               // linear combinations
               real mxxyPyzz = mfcba + mfabc;
               real mxxyMyzz = mfcba - mfabc;

               real mxxzPyyz = mfcab + mfacb;
               real mxxzMyyz = mfcab - mfacb;

               real mxyyPxzz = mfbca + mfbac;
               real mxyyMxzz = mfbca - mfbac;

               //relax
               wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*fabs(mfbbb)/(fabs(mfbbb)+qudricLimit);
               mfbbb     += wadjust * (-mfbbb);
               wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*fabs(mxxyPyzz)/(fabs(mxxyPyzz)+qudricLimit);
               mxxyPyzz  += wadjust * (-mxxyPyzz);
               wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*fabs(mxxyMyzz)/(fabs(mxxyMyzz)+qudricLimit);
               mxxyMyzz  += wadjust * (-mxxyMyzz);
               wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*fabs(mxxzPyyz)/(fabs(mxxzPyyz)+qudricLimit);
               mxxzPyyz  += wadjust * (-mxxzPyyz);
               wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*fabs(mxxzMyyz)/(fabs(mxxzMyyz)+qudricLimit);
               mxxzMyyz  += wadjust * (-mxxzMyyz);
               wadjust    = OxyyPxzz+(c1o1-OxyyPxzz)*fabs(mxyyPxzz)/(fabs(mxyyPxzz)+qudricLimit);
               mxyyPxzz  += wadjust * (-mxyyPxzz);
               wadjust    = OxyyMxzz+(c1o1-OxyyMxzz)*fabs(mxyyMxzz)/(fabs(mxyyMxzz)+qudricLimit);
               mxyyMxzz  += wadjust * (-mxyyMxzz);

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

               mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab);
               mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb);
               mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb);

               mfcca = CUMcca + (mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfcac = CUMcac + (mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-1)*oMdrho;
               mfacc = CUMacc + (mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-1)*oMdrho;

               //5.
               mfbcc = CUMbcc + (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) * oMdrho;
               mfcbc = CUMcbc + (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) * oMdrho;
               mfccb = CUMccb + (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) * oMdrho;

               //6.
               mfccc = CUMccc  -((-c4o1 *  mfbbb * mfbbb
                  -       (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                  -  c4o1 * (mfabb * mfcbb + mfbac * mfbca + mfbba * mfbbc)
                  -  c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
                  +(c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                     +  c2o1 * (mfcaa * mfaca * mfaac)
                     + c16o1 *  mfbba * mfbab * mfabb)
                  - c1o3* (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
                  - c1o9* (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
                  +(c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
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
               m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + c1o1 * oMdrho) * (vz2 - vvz) * c1o2;
               m1 = -mfaac        - c2o1 * mfaab *  vvz         +  mfaaa                * (c1o1 - vz2)              - c1o1 * oMdrho * vz2;
               m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + c1o1 * oMdrho) * (vz2 + vvz) * c1o2;
               mfaaa = m0;
               mfaab = m1;
               mfaac = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
               m1 = -mfabc        - c2o1 * mfabb *  vvz         + mfaba * (c1o1 - vz2);
               m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
               mfaba = m0;
               mfabb = m1;
               mfabc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
               m1 = -mfacc        - c2o1 * mfacb *  vvz         +  mfaca                  * (c1o1 - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
               mfaca = m0;
               mfacb = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
               m1 = -mfbac        - c2o1 * mfbab *  vvz         + mfbaa * (c1o1 - vz2);
               m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
               mfbaa = m0;
               mfbab = m1;
               mfbac = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
               m1 = -mfbbc        - c2o1 * mfbbb *  vvz         + mfbba * (c1o1 - vz2);
               m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
               mfbba = m0;
               mfbbb = m1;
               mfbbc = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
               m1 = -mfbcc        - c2o1 * mfbcb *  vvz         + mfbca * (c1o1 - vz2);
               m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
               mfbca = m0;
               mfbcb = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
               m1 = -mfcac        - c2o1 * mfcab *  vvz         +  mfcaa                  * (c1o1 - vz2)              - c1o3 * oMdrho * vz2;
               m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
               mfcaa = m0;
               mfcab = m1;
               mfcac = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
               m1 = -mfcbc        - c2o1 * mfcbb *  vvz         + mfcba * (c1o1 - vz2);
               m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
               mfcba = m0;
               mfcbb = m1;
               mfcbc = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
               m1 = -mfccc        - c2o1 * mfccb *  vvz         +  mfcca                  * (c1o1 - vz2)              - c1o9 * oMdrho * vz2;
               m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
               mfcca = m0;
               mfccb = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // Y - Dir
               m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfaca        - c2o1 * mfaba *  vvy         +  mfaaa                  * (c1o1 - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
               mfaaa = m0;
               mfaba = m1;
               mfaca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfacb        - c2o1 * mfabb *  vvy         +  mfaab                  * (c1o1 - vy2)              - c2o3 * oMdrho * vy2;
               m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
               mfaab = m0;
               mfabb = m1;
               mfacb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfacc        - c2o1 * mfabc *  vvy         +  mfaac                  * (c1o1 - vy2)              - c1o6 * oMdrho * vy2;
               m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
               mfaac = m0;
               mfabc = m1;
               mfacc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
               m1 = -mfbca        - c2o1 * mfbba *  vvy         + mfbaa * (c1o1 - vy2);
               m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
               mfbaa = m0;
               mfbba = m1;
               mfbca = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
               m1 = -mfbcb        - c2o1 * mfbbb *  vvy         + mfbab * (c1o1 - vy2);
               m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
               mfbab = m0;
               mfbbb = m1;
               mfbcb = m2;
               /////////b//////////////////////////////////////////////////////////////////////////
               m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
               m1 = -mfbcc        - c2o1 * mfbbc *  vvy         + mfbac * (c1o1 - vy2);
               m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
               mfbac = m0;
               mfbbc = m1;
               mfbcc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfcca        - c2o1 * mfcba *  vvy         +  mfcaa                   * (c1o1 - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
               mfcaa = m0;
               mfcba = m1;
               mfcca = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfccb        - c2o1 * mfcbb *  vvy         +  mfcab                  * (c1o1 - vy2)              - c2o9 * oMdrho * vy2;
               m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
               mfcab = m0;
               mfcbb = m1;
               mfccb = m2;
               /////////c//////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
               m1 = -mfccc        - c2o1 * mfcbc *  vvy         +  mfcac                   * (c1o1 - vy2)              - c1o18 * oMdrho * vy2;
               m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
               mfcac = m0;
               mfcbc = m1;
               mfccc = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
               ////////////////////////////////////////////////////////////////////////////////////
               // X - Dir
               m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcaa        - c2o1 * mfbaa *  vvx         +  mfaaa                   * (c1o1 - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
               mfaaa = m0;
               mfbaa = m1;
               mfcaa = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcba        - c2o1 * mfbba *  vvx         +  mfaba                  * (c1o1 - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
               mfaba = m0;
               mfbba = m1;
               mfcba = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcca        - c2o1 * mfbca *  vvx         +  mfaca                   * (c1o1 - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
               mfaca = m0;
               mfbca = m1;
               mfcca = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcab        - c2o1 * mfbab *  vvx         +  mfaab                  * (c1o1 - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
               mfaab = m0;
               mfbab = m1;
               mfcab = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcbb        - c2o1 * mfbbb *  vvx         +  mfabb                  * (c1o1 - vx2)              - c4o9 * oMdrho * vx2;
               m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
               mfabb = m0;
               mfbbb = m1;
               mfcbb = m2;
               ///////////b////////////////////////////////////////////////////////////////////////
               m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfccb        - c2o1 * mfbcb *  vvx         +  mfacb                  * (c1o1 - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
               mfacb = m0;
               mfbcb = m1;
               mfccb = m2;
               ////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////////////////////////////////////////////
               m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcac        - c2o1 * mfbac *  vvx         +  mfaac                   * (c1o1 - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
               mfaac = m0;
               mfbac = m1;
               mfcac = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfcbc        - c2o1 * mfbbc *  vvx         +  mfabc                  * (c1o1 - vx2)              - c1o9 * oMdrho * vx2;
               m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
               mfabc = m0;
               mfbbc = m1;
               mfcbc = m2;
               ///////////c////////////////////////////////////////////////////////////////////////
               m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
               m1 = -mfccc        - c2o1 * mfbcc *  vvx         +  mfacc                   * (c1o1 - vx2)              - c1o36 * oMdrho * vx2;
               m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
               mfacc = m0;
               mfbcc = m1;
               mfccc = m2;

               //////////////////////////////////////////////////////////////////////////
               //proof correctness
               //////////////////////////////////////////////////////////////////////////
#ifdef  PROOF_CORRECTNESS
               real rho_post = (mfaaa+mfaac+mfaca+mfcaa+mfacc+mfcac+mfccc+mfcca)
                  +(mfaab+mfacb+mfcab+mfccb)+(mfaba+mfabc+mfcba+mfcbc)+(mfbaa+mfbac+mfbca+mfbcc)
                  +(mfabb+mfcbb)+(mfbab+mfbcb)+(mfbba+mfbbc)+mfbbb;
               //real dif = fabs(rho - rho_post);
               real dif = rho - rho_post;
#ifdef SINGLEPRECISION
               if (dif > 10.0E-7 || dif < -10.0E-7)
#else
               if (dif > 10.0E-15 || dif < -10.0E-15)
#endif
               {
                  UB_THROW(UbException(UB_EXARGS, "rho="+UbSystem::toString(rho)+", rho_post="+UbSystem::toString(rho_post)
                     +" dif="+UbSystem::toString(dif)
                     +" rho is not correct for node "+UbSystem::toString(x1)+","+UbSystem::toString(x2)+","+UbSystem::toString(x3)
                     +" in " + block.lock()->toString()+" step = "+UbSystem::toString(step)));
               }
#endif
               //////////////////////////////////////////////////////////////////////////
               //write distribution
               //////////////////////////////////////////////////////////////////////////
               (*this->localDistributions)(eP00, x1, x2, x3)    = mfabb;
               (*this->localDistributions)(e0P0, x1, x2, x3)    = mfbab;
               (*this->localDistributions)(e00P, x1, x2, x3)    = mfbba;
               (*this->localDistributions)(ePP0, x1, x2, x3)   = mfaab;
               (*this->localDistributions)(eMP0, x1p, x2, x3)   = mfcab;
               (*this->localDistributions)(eP0P, x1, x2, x3)   = mfaba;
               (*this->localDistributions)(eM0P, x1p, x2, x3)   = mfcba;
               (*this->localDistributions)(e0PP, x1, x2, x3)   = mfbaa;
               (*this->localDistributions)(e0MP, x1, x2p, x3)   = mfbca;
               (*this->localDistributions)(ePPP, x1, x2, x3)  = mfaaa;
               (*this->localDistributions)(eMPP, x1p, x2, x3)  = mfcaa;
               (*this->localDistributions)(ePMP, x1, x2p, x3)  = mfaca;
               (*this->localDistributions)(eMMP, x1p, x2p, x3)  = mfcca;

               (*this->nonLocalDistributions)(eM00, x1p, x2, x3) = mfcbb;
               (*this->nonLocalDistributions)(e0M0, x1, x2p, x3) = mfbcb;
               (*this->nonLocalDistributions)(e00M, x1, x2, x3p) = mfbbc;
               (*this->nonLocalDistributions)(eMM0, x1p, x2p, x3) = mfccb;
               (*this->nonLocalDistributions)(ePM0, x1, x2p, x3) = mfacb;
               (*this->nonLocalDistributions)(eM0M, x1p, x2, x3p) = mfcbc;
               (*this->nonLocalDistributions)(eP0M, x1, x2, x3p) = mfabc;
               (*this->nonLocalDistributions)(e0MM, x1, x2p, x3p) = mfbcc;
               (*this->nonLocalDistributions)(e0PM, x1, x2, x3p) = mfbac;
               (*this->nonLocalDistributions)(eMMM, x1p, x2p, x3p) = mfccc;
               (*this->nonLocalDistributions)(ePMM, x1, x2p, x3p) = mfacc;
               (*this->nonLocalDistributions)(eMPM, x1p, x2, x3p) = mfcac;
               (*this->nonLocalDistributions)(ePPM, x1, x2, x3p) = mfaac;

               (*this->zeroDistributions)(x1, x2, x3) = mfbbb;
               //////////////////////////////////////////////////////////////////////////

            }
         }
      }
   }
   //timer.stop();
}
//////////////////////////////////////////////////////////////////////////
real K16IncompressibleNavierStokes::getCalculationTime()
{
   return timer.getTimeInSeconds();
}
//////////////////////////////////////////////////////////////////////////
void K16IncompressibleNavierStokes::setRelaxationParameter(Parameter p)
{
   parameter = p;
}


//! \}
