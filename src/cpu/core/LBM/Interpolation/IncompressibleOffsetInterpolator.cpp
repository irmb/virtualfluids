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
//! \file IncompressibleOffsetInterpolator.cpp
//! \ingroup Interpolation
//! \author Konstantin Kutscher
//=======================================================================================

#include "IncompressibleOffsetInterpolator.h"
#include "D3Q27System.h"



//////////////////////////////////////////////////////////////////////////
IncompressibleOffsetInterpolator::IncompressibleOffsetInterpolator(real omegaC, real omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{

}

//////////////////////////////////////////////////////////////////////////
InterpolationProcessorPtr IncompressibleOffsetInterpolator::clone()
{
   InterpolationProcessorPtr iproc = InterpolationProcessorPtr (new IncompressibleOffsetInterpolator(this->omegaC, this->omegaF));
   //dynamicPointerCast<D3Q27IncompressibleOffsetInterpolationProcessor>(iproc)->forcingC = forcingC;
   //dynamicPointerCast<D3Q27IncompressibleOffsetInterpolationProcessor>(iproc)->forcingF = forcingF;
   return iproc;
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::setOmegas( real omegaC, real omegaF )
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::setOffsets(real xoff, real yoff, real zoff)
{
   this->xoff = xoff;
   this->yoff = yoff;
   this->zoff = zoff;     
   this->xoff_sq = xoff * xoff;
   this->yoff_sq = yoff * yoff;
   this->zoff_sq = zoff * zoff;
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff)
{
    using namespace vf::basics::constant;

   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellC, omegaC, c1o2);
   calcInterpolatedNode(icellF.BSW, omegaF, -c1o4, -c1o4, -c1o4, calcPressBSW(), -c1o1, -c1o1, -c1o1);
   calcInterpolatedNode(icellF.BNE, omegaF,  c1o4,  c1o4, -c1o4, calcPressBNE(),  c1o1,  c1o1, -c1o1);
   calcInterpolatedNode(icellF.TNW, omegaF, -c1o4,  c1o4,  c1o4, calcPressTNW(), -c1o1,  c1o1,  c1o1);
   calcInterpolatedNode(icellF.TSE, omegaF,  c1o4, -c1o4,  c1o4, calcPressTSE(),  c1o1, -c1o1,  c1o1);
   calcInterpolatedNode(icellF.BNW, omegaF, -c1o4,  c1o4, -c1o4, calcPressBNW(), -c1o1,  c1o1, -c1o1);
   calcInterpolatedNode(icellF.BSE, omegaF,  c1o4, -c1o4, -c1o4, calcPressBSE(),  c1o1, -c1o1, -c1o1);
   calcInterpolatedNode(icellF.TSW, omegaF, -c1o4, -c1o4,  c1o4, calcPressTSW(), -c1o1, -c1o1,  c1o1);
   calcInterpolatedNode(icellF.TNE, omegaF,  c1o4,  c1o4,  c1o4, calcPressTNE(),  c1o1,  c1o1,  c1o1);
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellF, omegaF, vf::basics::constant::c2o1);
   calcInterpolatedNodeFC(icellC, omegaC);
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::calcMoments(const real* const f, real omega, real& press, real& vx1, real& vx2, real& vx3, 
                                                    real& kxy, real& kyz, real& kxz, real& kxxMyy, real& kxxMzz)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   //UBLOG(logINFO,"D3Q27System::dM0M  = " << D3Q27System::dM0M);
   //UBLOG(logINFO,"BW  = " << BW);;

   real rho = c0o1;
   D3Q27System::calcIncompMacroscopicValues(f,rho,vx1,vx2,vx3);
   
   //////////////////////////////////////////////////////////////////////////
   //DRAFT
   //if (omega == omegaC)
   //{
   //   vx1 += forcingC*0.5;
   //} 
   //else
   //{
   //   vx1 += forcingF*0.5;
   //}
   //////////////////////////////////////////////////////////////////////////

   //press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);
   press = rho; //interpolate rho!

   kxy   = -c3o1*omega*((((f[dMMP]+f[dPPM])-(f[dMPP]+f[dPMM]))+((f[dMMM]+f[dPPP])-(f[dMPM]+f[dPMP])))+((f[dMM0]+f[dPP0])-(f[dMP0]+f[dPM0]))-(vx1*vx2));// might not be optimal MG 25.2.13
   kyz   = -c3o1*omega*((((f[dMMM]+f[dPPP])-(f[dPMP]+f[dMPM]))+((f[dPMM]+f[dMPP])-(f[dMMP]+f[dPPM])))+((f[d0MM]+f[d0PP])-(f[d0MP]+f[d0PM]))-(vx2*vx3));
   kxz   = -c3o1*omega*((((f[dMPM]+f[dPMP])-(f[dMMP]+f[dPPM]))+((f[dMMM]+f[dPPP])-(f[dPMM]+f[dMPP])))+((f[dM0M]+f[dP0P])-(f[dM0P]+f[dP0M]))-(vx1*vx3));
   kxxMyy = -c3o1/c2o1*omega*((((f[dM0M]+f[dP0P])-(f[d0MM]+f[d0PP]))+((f[dM0P]+f[dP0M])-(f[d0MP]+f[d0PM])))+((f[dM00]+f[dP00])-(f[d0M0]+f[d0P0]))-(vx1*vx1-vx2*vx2));
   kxxMzz = -c3o1/c2o1*omega*((((f[dMP0]+f[dPM0])-(f[d0MM]+f[d0PP]))+((f[dMM0]+f[dPP0])-(f[d0MP]+f[d0PM])))+((f[dM00]+f[dP00])-(f[d00M]+f[d00P]))-(vx1*vx1-vx3*vx3));
   //kxxMzz = -c3o1/c2o1*omega*(((((f[NW]+f[SE])-(f[BS]+f[TN]))+((f[SW]+f[NE])-(f[17]+f[BN])))+((f[W]+f[dP00])-(f[B]+f[T])))-(vx1*vx1-vx3*vx3));

   //UBLOG(logINFO, "t1 = "<<(((f[NW]+f[SE])-(f[BS]+f[TN]))+((f[SW]+f[NE])-(f[17]+f[BN])))+((f[W]+f[dP00])-(f[B]+f[T])));
   //UBLOG(logINFO, "kxxMzz = "<<kxxMzz);

   //UBLOG(logINFO,"f[BW]  = " << f[BW] << " BW  = " << BW);
   //UBLOG(logINFO,"f[BE]  = " << f[BE] << " BE  = " << BE);
   //UBLOG(logINFO,"f[NW]  = " << f[NW] << " NW  = " << NW);
   //UBLOG(logINFO,"f[SE]  = " << f[SE] << " SE  = " << SE);
   //UBLOG(logINFO,"f[BS]  = " << f[BS] << " BS  = " << BS);
   //UBLOG(logINFO,"f[TN]  = " << f[TN] << " TN  = " << TN);
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::calcInterpolatedCoefficiets(const D3Q27ICell& icell, real omega, real eps_new)
{
    using namespace vf::basics::constant;

   real        vx1_SWT,vx2_SWT,vx3_SWT;
   real        vx1_NWT,vx2_NWT,vx3_NWT;
   real        vx1_NET,vx2_NET,vx3_NET;
   real        vx1_SET,vx2_SET,vx3_SET;
   real        vx1_SWB,vx2_SWB,vx3_SWB;
   real        vx1_NWB,vx2_NWB,vx3_NWB;
   real        vx1_NEB,vx2_NEB,vx3_NEB;
   real        vx1_SEB,vx2_SEB,vx3_SEB;

   real        kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT;
   real        kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT;
   real        kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET;
   real        kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET;
   real        kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB;
   real        kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB;
   real        kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB;
   real        kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB;

   calcMoments(icell.TSW,omega,press_SWT,vx1_SWT,vx2_SWT,vx3_SWT, kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT);
   calcMoments(icell.TNW,omega,press_NWT,vx1_NWT,vx2_NWT,vx3_NWT, kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT);
   calcMoments(icell.TNE,omega,press_NET,vx1_NET,vx2_NET,vx3_NET, kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET);
   calcMoments(icell.TSE,omega,press_SET,vx1_SET,vx2_SET,vx3_SET, kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET);
   calcMoments(icell.BSW,omega,press_SWB,vx1_SWB,vx2_SWB,vx3_SWB, kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB);
   calcMoments(icell.BNW,omega,press_NWB,vx1_NWB,vx2_NWB,vx3_NWB, kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB);
   calcMoments(icell.BNE,omega,press_NEB,vx1_NEB,vx2_NEB,vx3_NEB, kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB);
   calcMoments(icell.BSE,omega,press_SEB,vx1_SEB,vx2_SEB,vx3_SEB, kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB);

   //real dxRho=c1o4*((press_NET-press_SWB)+(press_SET-press_NWB)+(press_NEB-press_SWT)+(press_SEB-press_NWT));
   //real dyRho=c1o4*((press_NET-press_SWB)-(press_SET-press_NWB)+(press_NEB-press_SWT)-(press_SEB-press_NWT));
   //real dzRho=c1o4*((press_NET-press_SWB)+(press_SET-press_NWB)-(press_NEB-press_SWT)-(press_SEB-press_NWT));

   //   kxyFromfcNEQ_SWT+=vx1_SWT*dyRho+vx2_SWT*dxRho;
   //   kxyFromfcNEQ_NWT+=vx1_NWT*dyRho+vx2_NWT*dxRho;
   //   kxyFromfcNEQ_NET+=vx1_NET*dyRho+vx2_NET*dxRho;
   //   kxyFromfcNEQ_SET+=vx1_SET*dyRho+vx2_SET*dxRho;
   //   kxyFromfcNEQ_SWB+=vx1_SWB*dyRho+vx2_SWB*dxRho;
   //   kxyFromfcNEQ_NWB+=vx1_NWB*dyRho+vx2_NWB*dxRho;
   //   kxyFromfcNEQ_NEB+=vx1_NEB*dyRho+vx2_NEB*dxRho;
   //   kxyFromfcNEQ_SEB+=vx1_SEB*dyRho+vx2_SEB*dxRho;

   //   kyzFromfcNEQ_SWT+=vx3_SWT*dyRho+vx2_SWT*dzRho;
   //   kyzFromfcNEQ_NWT+=vx3_NWT*dyRho+vx2_NWT*dzRho;
   //   kyzFromfcNEQ_NET+=vx3_NET*dyRho+vx2_NET*dzRho;
   //   kyzFromfcNEQ_SET+=vx3_SET*dyRho+vx2_SET*dzRho;
   //   kyzFromfcNEQ_SWB+=vx3_SWB*dyRho+vx2_SWB*dzRho;
   //   kyzFromfcNEQ_NWB+=vx3_NWB*dyRho+vx2_NWB*dzRho;
   //   kyzFromfcNEQ_NEB+=vx3_NEB*dyRho+vx2_NEB*dzRho;
   //   kyzFromfcNEQ_SEB+=vx3_SEB*dyRho+vx2_SEB*dzRho;

   //   kxzFromfcNEQ_SWT+=vx1_SWT*dzRho+vx3_SWT*dxRho;
   //   kxzFromfcNEQ_NWT+=vx1_NWT*dzRho+vx3_NWT*dxRho;
   //   kxzFromfcNEQ_NET+=vx1_NET*dzRho+vx3_NET*dxRho;
   //   kxzFromfcNEQ_SET+=vx1_SET*dzRho+vx3_SET*dxRho;
   //   kxzFromfcNEQ_SWB+=vx1_SWB*dzRho+vx3_SWB*dxRho;
   //   kxzFromfcNEQ_NWB+=vx1_NWB*dzRho+vx3_NWB*dxRho;
   //   kxzFromfcNEQ_NEB+=vx1_NEB*dzRho+vx3_NEB*dxRho;
   //   kxzFromfcNEQ_SEB+=vx1_SEB*dzRho+vx3_SEB*dxRho;

   //   kxxMyyFromfcNEQ_SWT+=vx1_SWT*dxRho-vx2_SWT*dyRho;
   //   kxxMyyFromfcNEQ_NWT+=vx1_NWT*dxRho-vx2_NWT*dyRho;
   //   kxxMyyFromfcNEQ_NET+=vx1_NET*dxRho-vx2_NET*dyRho;
   //   kxxMyyFromfcNEQ_SET+=vx1_SET*dxRho-vx2_SET*dyRho;
   //   kxxMyyFromfcNEQ_SWB+=vx1_SWB*dxRho-vx2_SWB*dyRho;
   //   kxxMyyFromfcNEQ_NWB+=vx1_NWB*dxRho-vx2_NWB*dyRho;
   //   kxxMyyFromfcNEQ_NEB+=vx1_NEB*dxRho-vx2_NEB*dyRho;
   //   kxxMyyFromfcNEQ_SEB+=vx1_SEB*dxRho-vx2_SEB*dyRho;

   //   kxxMzzFromfcNEQ_SWT+=vx1_SWT*dxRho-vx3_SWT*dzRho;
   //   kxxMzzFromfcNEQ_NWT+=vx1_NWT*dxRho-vx3_NWT*dzRho;
   //   kxxMzzFromfcNEQ_NET+=vx1_NET*dxRho-vx3_NET*dzRho;
   //   kxxMzzFromfcNEQ_SET+=vx1_SET*dxRho-vx3_SET*dzRho;
   //   kxxMzzFromfcNEQ_SWB+=vx1_SWB*dxRho-vx3_SWB*dzRho;
   //   kxxMzzFromfcNEQ_NWB+=vx1_NWB*dxRho-vx3_NWB*dzRho;
   //   kxxMzzFromfcNEQ_NEB+=vx1_NEB*dxRho-vx3_NEB*dzRho;
   //   kxxMzzFromfcNEQ_SEB+=vx1_SEB*dxRho-vx3_SEB*dzRho;


      //kxxMzzFromfcNEQ_SWT=0.0;
      //kxxMzzFromfcNEQ_NWT=0.0;
      //kxxMzzFromfcNEQ_NET=0.0;
      //kxxMzzFromfcNEQ_SET=0.0;
      //kxxMzzFromfcNEQ_SWB=0.0;
      //kxxMzzFromfcNEQ_NWB=0.0;
      //kxxMzzFromfcNEQ_NEB=0.0;
      //kxxMzzFromfcNEQ_SEB=0.0;





   a0 = (-kxxMyyFromfcNEQ_NEB - kxxMyyFromfcNEQ_NET + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_NWT -
      kxxMyyFromfcNEQ_SEB - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_SWT -
      kxxMzzFromfcNEQ_NEB - kxxMzzFromfcNEQ_NET + kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_NWT -
      kxxMzzFromfcNEQ_SEB - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_SWT -
      c2o1*kxyFromfcNEQ_NEB - c2o1*kxyFromfcNEQ_NET - c2o1*kxyFromfcNEQ_NWB - c2o1*kxyFromfcNEQ_NWT +
      c2o1*kxyFromfcNEQ_SEB + c2o1*kxyFromfcNEQ_SET + c2o1*kxyFromfcNEQ_SWB + c2o1*kxyFromfcNEQ_SWT +
      c2o1*kxzFromfcNEQ_NEB - c2o1*kxzFromfcNEQ_NET + c2o1*kxzFromfcNEQ_NWB - c2o1*kxzFromfcNEQ_NWT +
      c2o1*kxzFromfcNEQ_SEB - c2o1*kxzFromfcNEQ_SET + c2o1*kxzFromfcNEQ_SWB - c2o1*kxzFromfcNEQ_SWT +
      c8o1*vx1_NEB + c8o1*vx1_NET + c8o1*vx1_NWB + c8o1*vx1_NWT + c8o1*vx1_SEB +
      c8o1*vx1_SET + c8o1*vx1_SWB + c8o1*vx1_SWT + c2o1*vx2_NEB + c2o1*vx2_NET -
      c2o1*vx2_NWB - c2o1*vx2_NWT - c2o1*vx2_SEB - c2o1*vx2_SET + c2o1*vx2_SWB +
      c2o1*vx2_SWT - c2o1*vx3_NEB + c2o1*vx3_NET + c2o1*vx3_NWB - c2o1*vx3_NWT -
      c2o1*vx3_SEB + c2o1*vx3_SET + c2o1*vx3_SWB - c2o1*vx3_SWT)/c64o1;
   b0 = (c2o1*kxxMyyFromfcNEQ_NEB + c2o1*kxxMyyFromfcNEQ_NET + c2o1*kxxMyyFromfcNEQ_NWB + c2o1*kxxMyyFromfcNEQ_NWT -
      c2o1*kxxMyyFromfcNEQ_SEB - c2o1*kxxMyyFromfcNEQ_SET - c2o1*kxxMyyFromfcNEQ_SWB - c2o1*kxxMyyFromfcNEQ_SWT -
      kxxMzzFromfcNEQ_NEB - kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_SWT -
      c2o1*kxyFromfcNEQ_NEB - c2o1*kxyFromfcNEQ_NET + c2o1*kxyFromfcNEQ_NWB + c2o1*kxyFromfcNEQ_NWT -
      c2o1*kxyFromfcNEQ_SEB - c2o1*kxyFromfcNEQ_SET + c2o1*kxyFromfcNEQ_SWB + c2o1*kxyFromfcNEQ_SWT +
      c2o1*kyzFromfcNEQ_NEB - c2o1*kyzFromfcNEQ_NET + c2o1*kyzFromfcNEQ_NWB - c2o1*kyzFromfcNEQ_NWT +
      c2o1*kyzFromfcNEQ_SEB - c2o1*kyzFromfcNEQ_SET + c2o1*kyzFromfcNEQ_SWB - c2o1*kyzFromfcNEQ_SWT +
      c2o1*vx1_NEB + c2o1*vx1_NET - c2o1*vx1_NWB - c2o1*vx1_NWT -
      c2o1*vx1_SEB - c2o1*vx1_SET + c2o1*vx1_SWB + c2o1*vx1_SWT +
      c8o1*vx2_NEB + c8o1*vx2_NET + c8o1*vx2_NWB + c8o1*vx2_NWT +
      c8o1*vx2_SEB + c8o1*vx2_SET + c8o1*vx2_SWB + c8o1*vx2_SWT -
      c2o1*vx3_NEB + c2o1*vx3_NET - c2o1*vx3_NWB + c2o1*vx3_NWT +
      c2o1*vx3_SEB - c2o1*vx3_SET + c2o1*vx3_SWB - c2o1*vx3_SWT)/c64o1;
   c0 = (kxxMyyFromfcNEQ_NEB - kxxMyyFromfcNEQ_NET + kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT -
      c2o1*kxxMzzFromfcNEQ_NEB + c2o1*kxxMzzFromfcNEQ_NET - c2o1*kxxMzzFromfcNEQ_NWB + c2o1*kxxMzzFromfcNEQ_NWT -
      c2o1*kxxMzzFromfcNEQ_SEB + c2o1*kxxMzzFromfcNEQ_SET - c2o1*kxxMzzFromfcNEQ_SWB + c2o1*kxxMzzFromfcNEQ_SWT -
      c2o1*kxzFromfcNEQ_NEB - c2o1*kxzFromfcNEQ_NET + c2o1*kxzFromfcNEQ_NWB + c2o1*kxzFromfcNEQ_NWT -
      c2o1*kxzFromfcNEQ_SEB - c2o1*kxzFromfcNEQ_SET + c2o1*kxzFromfcNEQ_SWB + c2o1*kxzFromfcNEQ_SWT -
      c2o1*kyzFromfcNEQ_NEB - c2o1*kyzFromfcNEQ_NET - c2o1*kyzFromfcNEQ_NWB - c2o1*kyzFromfcNEQ_NWT +
      c2o1*kyzFromfcNEQ_SEB + c2o1*kyzFromfcNEQ_SET + c2o1*kyzFromfcNEQ_SWB + c2o1*kyzFromfcNEQ_SWT -
      c2o1*vx1_NEB + c2o1*vx1_NET + c2o1*vx1_NWB - c2o1*vx1_NWT -
      c2o1*vx1_SEB + c2o1*vx1_SET + c2o1*vx1_SWB - c2o1*vx1_SWT -
      c2o1*vx2_NEB + c2o1*vx2_NET - c2o1*vx2_NWB + c2o1*vx2_NWT +
      c2o1*vx2_SEB - c2o1*vx2_SET + c2o1*vx2_SWB - c2o1*vx2_SWT +
      c8o1*vx3_NEB + c8o1*vx3_NET + c8o1*vx3_NWB + c8o1*vx3_NWT +
      c8o1*vx3_SEB + c8o1*vx3_SET + c8o1*vx3_SWB + c8o1*vx3_SWT)/c64o1;
   ax = (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT + vx1_SEB + vx1_SET - vx1_SWB - vx1_SWT)/c4o1;
   bx = (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT + vx2_SEB + vx2_SET - vx2_SWB - vx2_SWT)/c4o1;
   cx = (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT + vx3_SEB + vx3_SET - vx3_SWB - vx3_SWT)/c4o1;
   axx= (kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      c2o1*vx2_NEB + c2o1*vx2_NET - c2o1*vx2_NWB - c2o1*vx2_NWT -
      c2o1*vx2_SEB - c2o1*vx2_SET + c2o1*vx2_SWB + c2o1*vx2_SWT -
      c2o1*vx3_NEB + c2o1*vx3_NET + c2o1*vx3_NWB - c2o1*vx3_NWT -
      c2o1*vx3_SEB + c2o1*vx3_SET + c2o1*vx3_SWB - c2o1*vx3_SWT)/c16o1;
   bxx= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET - kxyFromfcNEQ_NWB - kxyFromfcNEQ_NWT +
      kxyFromfcNEQ_SEB + kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      c2o1*vx1_NEB - c2o1*vx1_NET + c2o1*vx1_NWB + c2o1*vx1_NWT +
      c2o1*vx1_SEB + c2o1*vx1_SET - c2o1*vx1_SWB - c2o1*vx1_SWT)/c8o1;
   cxx= (kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB - kxzFromfcNEQ_NWT +
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB - kxzFromfcNEQ_SWT +
      c2o1*vx1_NEB - c2o1*vx1_NET - c2o1*vx1_NWB + c2o1*vx1_NWT +
      c2o1*vx1_SEB - c2o1*vx1_SET - c2o1*vx1_SWB + c2o1*vx1_SWT)/c8o1;
   ay = (vx1_NEB + vx1_NET + vx1_NWB + vx1_NWT - vx1_SEB - vx1_SET - vx1_SWB - vx1_SWT)/c4o1;
   by = (vx2_NEB + vx2_NET + vx2_NWB + vx2_NWT - vx2_SEB - vx2_SET - vx2_SWB - vx2_SWT)/c4o1;
   cy = (vx3_NEB + vx3_NET + vx3_NWB + vx3_NWT - vx3_SEB - vx3_SET - vx3_SWB - vx3_SWT)/c4o1;
   ayy= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET + kxyFromfcNEQ_NWB + kxyFromfcNEQ_NWT -
      kxyFromfcNEQ_SEB - kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      c2o1*vx2_NEB - c2o1*vx2_NET + c2o1*vx2_NWB + c2o1*vx2_NWT +
      c2o1*vx2_SEB + c2o1*vx2_SET - c2o1*vx2_SWB - c2o1*vx2_SWT)/c8o1;
   byy= (-c2o1*kxxMyyFromfcNEQ_NEB - c2o1*kxxMyyFromfcNEQ_NET - c2o1*kxxMyyFromfcNEQ_NWB - c2o1*kxxMyyFromfcNEQ_NWT +
      c2o1*kxxMyyFromfcNEQ_SEB + c2o1*kxxMyyFromfcNEQ_SET + c2o1*kxxMyyFromfcNEQ_SWB + c2o1*kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET + kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_NWT -
      kxxMzzFromfcNEQ_SEB - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      c2o1*vx1_NEB + c2o1*vx1_NET - c2o1*vx1_NWB - c2o1*vx1_NWT -
      c2o1*vx1_SEB - c2o1*vx1_SET + c2o1*vx1_SWB + c2o1*vx1_SWT -
      c2o1*vx3_NEB + c2o1*vx3_NET - c2o1*vx3_NWB + c2o1*vx3_NWT +
      c2o1*vx3_SEB - c2o1*vx3_SET + c2o1*vx3_SWB - c2o1*vx3_SWT)/c16o1;
   cyy= (kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET + kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB - kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB - kyzFromfcNEQ_SWT +
      c2o1*vx2_NEB - c2o1*vx2_NET + c2o1*vx2_NWB - c2o1*vx2_NWT -
      c2o1*vx2_SEB + c2o1*vx2_SET - c2o1*vx2_SWB + c2o1*vx2_SWT)/c8o1;
   az = (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT - vx1_SEB + vx1_SET - vx1_SWB + vx1_SWT)/c4o1;
   bz = (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT - vx2_SEB + vx2_SET - vx2_SWB + vx2_SWT)/c4o1;
   cz = (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT - vx3_SEB + vx3_SET - vx3_SWB + vx3_SWT)/c4o1;
   azz= (-kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB + kxzFromfcNEQ_NWT -
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB + kxzFromfcNEQ_SWT +
      c2o1*vx3_NEB - c2o1*vx3_NET - c2o1*vx3_NWB + c2o1*vx3_NWT +
      c2o1*vx3_SEB - c2o1*vx3_SET - c2o1*vx3_SWB + c2o1*vx3_SWT)/c8o1;
   bzz= (-kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET - kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB + kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB + kyzFromfcNEQ_SWT +
      c2o1*vx3_NEB - c2o1*vx3_NET + c2o1*vx3_NWB - c2o1*vx3_NWT -
      c2o1*vx3_SEB + c2o1*vx3_SET - c2o1*vx3_SWB + c2o1*vx3_SWT)/c8o1;
   czz= (-kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_NWT -
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_SWT +
      c2o1*kxxMzzFromfcNEQ_NEB - c2o1*kxxMzzFromfcNEQ_NET + c2o1*kxxMzzFromfcNEQ_NWB - c2o1*kxxMzzFromfcNEQ_NWT +
      c2o1*kxxMzzFromfcNEQ_SEB - c2o1*kxxMzzFromfcNEQ_SET + c2o1*kxxMzzFromfcNEQ_SWB - c2o1*kxxMzzFromfcNEQ_SWT -
      c2o1*vx1_NEB + c2o1*vx1_NET + c2o1*vx1_NWB - c2o1*vx1_NWT -
      c2o1*vx1_SEB + c2o1*vx1_SET + c2o1*vx1_SWB - c2o1*vx1_SWT -
      c2o1*vx2_NEB + c2o1*vx2_NET - c2o1*vx2_NWB + c2o1*vx2_NWT +
      c2o1*vx2_SEB - c2o1*vx2_SET + c2o1*vx2_SWB - c2o1*vx2_SWT)/c16o1;
   axy= (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT - vx1_SEB - vx1_SET + vx1_SWB + vx1_SWT)/c2o1;
   bxy= (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT - vx2_SEB - vx2_SET + vx2_SWB + vx2_SWT)/c2o1;
   cxy= (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT - vx3_SEB - vx3_SET + vx3_SWB + vx3_SWT)/c2o1;
   axz= (-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT - vx1_SEB + vx1_SET + vx1_SWB - vx1_SWT)/c2o1;
   bxz= (-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT - vx2_SEB + vx2_SET + vx2_SWB - vx2_SWT)/c2o1;
   cxz= (-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT - vx3_SEB + vx3_SET + vx3_SWB - vx3_SWT)/c2o1;
   ayz= (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT + vx1_SEB - vx1_SET + vx1_SWB - vx1_SWT)/c2o1;
   byz= (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT + vx2_SEB - vx2_SET + vx2_SWB - vx2_SWT)/c2o1;
   cyz= (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT + vx3_SEB - vx3_SET + vx3_SWB - vx3_SWT)/c2o1;
   axyz=-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT + vx1_SEB - vx1_SET - vx1_SWB + vx1_SWT;
   bxyz=-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT + vx2_SEB - vx2_SET - vx2_SWB + vx2_SWT;
   cxyz=-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT + vx3_SEB - vx3_SET - vx3_SWB + vx3_SWT;

   /////////////////////B�SE!!!
   //axx=0;   ayy=0;   azz=0;
   //bxx=0;   byy=0;   bzz=0;
   //cxx=0;   cyy=0;   czz=0;
   ////////////////////!!!B�SE

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //kxyAverage       =(kxyFromfcNEQ_SWB+
   //                   kxyFromfcNEQ_SWT+
   //                   kxyFromfcNEQ_SET+
   //                   kxyFromfcNEQ_SEB+
   //                   kxyFromfcNEQ_NWB+
   //                   kxyFromfcNEQ_NWT+
   //                   kxyFromfcNEQ_NET+
   //                   kxyFromfcNEQ_NEB)*c1o8-(ay+bx);
   //kyzAverage       =(kyzFromfcNEQ_SWB+
   //                   kyzFromfcNEQ_SWT+
   //                   kyzFromfcNEQ_SET+
   //                   kyzFromfcNEQ_SEB+
   //                   kyzFromfcNEQ_NWB+
   //                   kyzFromfcNEQ_NWT+
   //                   kyzFromfcNEQ_NET+
   //                   kyzFromfcNEQ_NEB)*c1o8-(bz+cy);
   //kxzAverage       =(kxzFromfcNEQ_SWB+
   //                   kxzFromfcNEQ_SWT+
   //                   kxzFromfcNEQ_SET+
   //                   kxzFromfcNEQ_SEB+
   //                   kxzFromfcNEQ_NWB+
   //                   kxzFromfcNEQ_NWT+
   //                   kxzFromfcNEQ_NET+
   //                   kxzFromfcNEQ_NEB)*c1o8-(az+cx);
   //kxxMyyAverage    =(kxxMyyFromfcNEQ_SWB+
   //                   kxxMyyFromfcNEQ_SWT+
   //                   kxxMyyFromfcNEQ_SET+
   //                   kxxMyyFromfcNEQ_SEB+
   //                   kxxMyyFromfcNEQ_NWB+
   //                   kxxMyyFromfcNEQ_NWT+
   //                   kxxMyyFromfcNEQ_NET+
   //                   kxxMyyFromfcNEQ_NEB)*c1o8-(ax-by);
   //kxxMzzAverage    =(kxxMzzFromfcNEQ_SWB+
   //                  kxxMzzFromfcNEQ_SWT+
   //                  kxxMzzFromfcNEQ_SET+
   //                  kxxMzzFromfcNEQ_SEB+
   //                  kxxMzzFromfcNEQ_NWB+
   //                  kxxMzzFromfcNEQ_NWT+
   //                  kxxMzzFromfcNEQ_NET+
   //                  kxxMzzFromfcNEQ_NEB)*c1o8-(ax-cz);
   kxyAverage       = c0o1;//(kxyFromfcNEQ_SWB+
                    //kxyFromfcNEQ_SWT+
                    //kxyFromfcNEQ_SET+
                    //kxyFromfcNEQ_SEB+
                    //kxyFromfcNEQ_NWB+
                    //kxyFromfcNEQ_NWT+
                    //kxyFromfcNEQ_NET+
                    //kxyFromfcNEQ_NEB)*c1o8-(ay+bx);
   kyzAverage       = c0o1;//(kyzFromfcNEQ_SWB+
                       //kyzFromfcNEQ_SWT+
                       //kyzFromfcNEQ_SET+
                       //kyzFromfcNEQ_SEB+
                       //kyzFromfcNEQ_NWB+
                       //kyzFromfcNEQ_NWT+
                       //kyzFromfcNEQ_NET+
                       //kyzFromfcNEQ_NEB)*c1o8-(bz+cy);
   kxzAverage       = c0o1;//(kxzFromfcNEQ_SWB+
                       //kxzFromfcNEQ_SWT+
                       //kxzFromfcNEQ_SET+
                       //kxzFromfcNEQ_SEB+
                       //kxzFromfcNEQ_NWB+
                       //kxzFromfcNEQ_NWT+
                       //kxzFromfcNEQ_NET+
                       //kxzFromfcNEQ_NEB)*c1o8-(az+cx);
   kxxMyyAverage    = c0o1;//(kxxMyyFromfcNEQ_SWB+
                       //kxxMyyFromfcNEQ_SWT+
                       //kxxMyyFromfcNEQ_SET+
                       //kxxMyyFromfcNEQ_SEB+
                       //kxxMyyFromfcNEQ_NWB+
                       //kxxMyyFromfcNEQ_NWT+
                       //kxxMyyFromfcNEQ_NET+
                       //kxxMyyFromfcNEQ_NEB)*c1o8-(ax-by);
   kxxMzzAverage    = c0o1;//(kxxMzzFromfcNEQ_SWB+
                       //kxxMzzFromfcNEQ_SWT+
                       //kxxMzzFromfcNEQ_SET+
                       //kxxMzzFromfcNEQ_SEB+
                       //kxxMzzFromfcNEQ_NWB+
                       //kxxMzzFromfcNEQ_NWT+
                       //kxxMzzFromfcNEQ_NET+
                       //kxxMzzFromfcNEQ_NEB)*c1o8-(ax-cz);
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //
   // Bernd das Brot
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   a0 = a0 + xoff * ax + yoff * ay + zoff * az + xoff_sq * axx + yoff_sq * ayy + zoff_sq * azz + xoff*yoff*axy + xoff*zoff*axz + yoff*zoff*ayz + xoff*yoff*zoff*axyz ;
   ax = ax + c2o1 * xoff * axx + yoff * axy + zoff * axz + yoff*zoff*axyz;
   ay = ay + c2o1 * yoff * ayy + xoff * axy + zoff * ayz + xoff*zoff*axyz;
   az = az + c2o1 * zoff * azz + xoff * axz + yoff * ayz + xoff*yoff*axyz;
   b0 = b0 + xoff * bx + yoff * by + zoff * bz + xoff_sq * bxx + yoff_sq * byy + zoff_sq * bzz + xoff*yoff*bxy + xoff*zoff*bxz + yoff*zoff*byz + xoff*yoff*zoff*bxyz;
   bx = bx + c2o1 * xoff * bxx + yoff * bxy + zoff * bxz + yoff*zoff*bxyz;
   by = by + c2o1 * yoff * byy + xoff * bxy + zoff * byz + xoff*zoff*bxyz;
   bz = bz + c2o1 * zoff * bzz + xoff * bxz + yoff * byz + xoff*yoff*bxyz;
   c0 = c0 + xoff * cx + yoff * cy + zoff * cz + xoff_sq * cxx + yoff_sq * cyy + zoff_sq * czz + xoff*yoff*cxy + xoff*zoff*cxz + yoff*zoff*cyz + xoff*yoff*zoff*cxyz;
   cx = cx + c2o1 * xoff * cxx + yoff * cxy + zoff * cxz + yoff*zoff*cxyz;
   cy = cy + c2o1 * yoff * cyy + xoff * cxy + zoff * cyz + xoff*zoff*cxyz;
   cz = cz + c2o1 * zoff * czz + xoff * cxz + yoff * cyz + xoff*yoff*cxyz;
   axy= axy + zoff*axyz;
   axz= axz + yoff*axyz;
   ayz= ayz + xoff*axyz;
   bxy= bxy + zoff*bxyz;
   bxz= bxz + yoff*bxyz;
   byz= byz + xoff*bxyz;
   cxy= cxy + zoff*cxyz;
   cxz= cxz + yoff*cxyz;
   cyz= cyz + xoff*cxyz;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   const real o = omega;

   f_E = eps_new*((c2o1*(-c2o1*ax + by + cz-kxxMzzAverage-kxxMyyAverage))/(c27o1*o));
   f_N = eps_new*((c2o1*(ax - c2o1*by + cz+c2o1*kxxMyyAverage-kxxMzzAverage))/(c27o1*o));
   f_T = eps_new*((c2o1*(ax + by - c2o1*cz-kxxMyyAverage+c2o1*kxxMzzAverage))/(c27o1*o));
   f_NE = eps_new*(-(ax + c3o1*ay + c3o1*bx + by - c2o1*cz+c2o1*kxxMyyAverage-kxxMyyAverage+c3o1*kxyAverage)/(c54o1*o));
   f_SE = eps_new*(-(ax - c3o1*ay - c3o1*bx + by - c2o1*cz+c2o1*kxxMyyAverage-kxxMyyAverage-c3o1*kxyAverage)/(c54o1*o));
   f_TE = eps_new*(-(ax + c3o1*az - c2o1*by + c3o1*cx + cz+c2o1*kxxMyyAverage-kxxMzzAverage+c3o1*kxzAverage)/(c54o1*o));
   f_BE = eps_new*(-(ax - c3o1*az - c2o1*by - c3o1*cx + cz+c2o1*kxxMyyAverage-kxxMzzAverage-c3o1*kxzAverage)/(c54o1*o));
   f_TN = eps_new*(-(-c2o1*ax + by + c3o1*bz + c3o1*cy + cz-kxxMyyAverage-kxxMzzAverage+c3o1*kyzAverage)/(c54o1*o));
   f_BN = eps_new*(-(-c2o1*ax + by - c3o1*bz - c3o1*cy + cz-kxxMyyAverage-kxxMzzAverage-c3o1*kyzAverage)/(c54o1*o));
   f_ZERO = c0o1;
   f_TNE = eps_new*(-(ay + az + bx + bz + cx + cy+kxyAverage+kxzAverage+kyzAverage)/(c72o1*o));
   f_TSW = eps_new*((-ay + az - bx + bz + cx + cy-kxyAverage+kxzAverage+kyzAverage)/(c72o1*o));
   f_TSE = eps_new*((ay - az + bx + bz - cx + cy+kxyAverage-kxzAverage+kyzAverage)/(c72o1*o));
   f_TNW = eps_new*((ay + az + bx - bz + cx - cy+kxyAverage+kxzAverage-kyzAverage)/(c72o1*o));

   x_E = c1o4*eps_new*((c2o1*(-c4o1*axx + bxy + cxz))/(c27o1*o));
   x_N = c1o4*eps_new*((c2o1*(c2o1*axx - c2o1*bxy + cxz))/(c27o1*o));
   x_T = c1o4*eps_new*((c2o1*(c2o1*axx + bxy - c2o1*cxz))/(c27o1*o));
   x_NE = c1o4*eps_new*(-((c2o1*axx + c3o1*axy + c6o1*bxx + bxy - c2o1*cxz))/(c54o1*o));
   x_SE = c1o4*eps_new*(-((c2o1*axx - c3o1*axy - c6o1*bxx + bxy - c2o1*cxz))/(c54o1*o));
   x_TE = c1o4*eps_new*(-((c2o1*axx + c3o1*axz - c2o1*bxy + c6o1*cxx + cxz))/(c54o1*o));
   x_BE = c1o4*eps_new*(-((c2o1*axx - c3o1*axz - c2o1*bxy - c6o1*cxx + cxz))/(c54o1*o));
   x_TN = c1o4*eps_new*(-((-c4o1*axx + bxy + c3o1*bxz + c3o1*cxy + cxz))/(c54o1*o));
   x_BN = c1o4*eps_new*(-((-c4o1*axx + bxy - c3o1*bxz - c3o1*cxy + cxz))/(c54o1*o));
   x_ZERO = c0o1;
   x_TNE = c1o4*eps_new*(-((axy + axz + c2o1*bxx + bxz + c2o1*cxx + cxy))/(c72o1*o));
   x_TSW = c1o4*eps_new*(((-axy + axz - c2o1*bxx + bxz + c2o1*cxx + cxy))/(c72o1*o));
   x_TSE = c1o4*eps_new*(((axy - axz + c2o1*bxx + bxz - c2o1*cxx + cxy))/(c72o1*o));
   x_TNW = c1o4*eps_new*(((axy + axz + c2o1*bxx - bxz + c2o1*cxx - cxy))/(c72o1*o));

   y_E = c1o4*eps_new*(c2o1*(-c2o1*axy + c2o1*byy + cyz))/(c27o1*o);
   y_N = c1o4*eps_new*(c2o1*(axy - c4o1*byy + cyz))/(c27o1*o);
   y_T = c1o4*eps_new*(c2o1*(axy + c2o1*byy - c2o1*cyz))/(c27o1*o);
   y_NE = c1o4*eps_new*(-((axy + c6o1*ayy + c3o1*bxy + c2o1*byy - c2o1*cyz))/(c54o1*o));
   y_SE = c1o4*eps_new*(-((axy - c6o1*ayy - c3o1*bxy + c2o1*byy - c2o1*cyz))/(c54o1*o));
   y_TE = c1o4*eps_new*(-((axy + c3o1*ayz - c4o1*byy + c3o1*cxy + cyz))/(c54o1*o));
   y_BE = c1o4*eps_new*(-((axy - c3o1*ayz - c4o1*byy - c3o1*cxy + cyz))/(c54o1*o));
   y_TN = c1o4*eps_new*(-((-c2o1*axy + c2o1*byy + c3o1*byz + c6o1*cyy + cyz))/(c54o1*o));
   y_BN = c1o4*eps_new*(-((-c2o1*axy + c2o1*byy - c3o1*byz - c6o1*cyy + cyz))/(c54o1*o));
   y_ZERO = c0o1;
   y_TNE = c1o4*eps_new*(-((c2o1*ayy + ayz + bxy + byz + cxy + c2o1*cyy))/(c72o1*o));
   y_TSW = c1o4*eps_new*(((-c2o1*ayy + ayz - bxy + byz + cxy + c2o1*cyy))/(c72o1*o));
   y_TSE = c1o4*eps_new*(((c2o1*ayy - ayz + bxy + byz - cxy + c2o1*cyy))/(c72o1*o));
   y_TNW = c1o4*eps_new*(((c2o1*ayy + ayz + bxy - byz + cxy - c2o1*cyy))/(c72o1*o));

   z_E = c1o4*eps_new*((c2o1*(-c2o1*axz + byz + c2o1*czz))/(c27o1*o));
   z_N = c1o4*eps_new*((c2o1*(axz - c2o1*byz + c2o1*czz))/(c27o1*o));
   z_T = c1o4*eps_new*((c2o1*(axz + byz - c4o1*czz))/(c27o1*o));
   z_NE = c1o4*eps_new*(-((axz + c3o1*ayz + c3o1*bxz + byz - c4o1*czz))/(c54o1*o));
   z_SE = c1o4*eps_new*(-((axz - c3o1*ayz - c3o1*bxz + byz - c4o1*czz))/(c54o1*o));
   z_TE = c1o4*eps_new*(-((axz + c6o1*azz - c2o1*byz + c3o1*cxz + c2o1*czz))/(c54o1*o));
   z_BE = c1o4*eps_new*(-((axz - c6o1*azz - c2o1*byz - c3o1*cxz + c2o1*czz))/(c54o1*o));
   z_TN = c1o4*eps_new*(-((-c2o1*axz + byz + c6o1*bzz + c3o1*cyz + c2o1*czz))/(c54o1*o));
   z_BN = c1o4*eps_new*(-((-c2o1*axz + byz - c6o1*bzz - c3o1*cyz + c2o1*czz))/(c54o1*o));
   z_ZERO = c0o1;
   z_TNE = c1o4*eps_new*(-((ayz + c2o1*azz + bxz + c2o1*bzz + cxz + cyz))/(c72o1*o));
   z_TSW = c1o4*eps_new*(((-ayz + c2o1*azz - bxz + c2o1*bzz + cxz + cyz))/(c72o1*o));
   z_TSE = c1o4*eps_new*(((ayz - c2o1*azz + bxz + c2o1*bzz - cxz + cyz))/(c72o1*o));
   z_TNW = c1o4*eps_new*(((ayz + c2o1*azz + bxz - c2o1*bzz + cxz - cyz))/(c72o1*o));

   xy_E   =   c1o16*eps_new *((                       c2o1*cxyz)/(c27o1*o));
   xy_N   =   c1o16*eps_new *((                       c2o1*cxyz)/(c27o1*o));
   xy_T   = -(c1o16*eps_new *((                       c4o1*cxyz)/(c27o1*o)));
   xy_NE  =   c1o16*eps_new *(                            cxyz /(c27o1*o));
   xy_SE  =   c1o16*eps_new *(                            cxyz /(c27o1*o));
   xy_TE  = -(c1o16*eps_new *(( c3o1*axyz            +     cxyz)/(c54o1*o)));
   xy_BE  = -(c1o16*eps_new *((-c3o1*axyz            +     cxyz)/(c54o1*o)));
   xy_TN  = -(c1o16*eps_new *((            c3o1*bxyz +     cxyz)/(c54o1*o)));
   xy_BN  = -(c1o16*eps_new *((          - c3o1*bxyz +     cxyz)/(c54o1*o)));
   //xy_ZERO=   c1o16*eps_new;
   xy_TNE = -(c1o16*eps_new *((     axyz +     bxyz           )/(c72o1*o)));
   xy_TSW =   c1o16*eps_new *((     axyz +     bxyz           )/(c72o1*o));
   xy_TSE =   c1o16*eps_new *((-    axyz +     bxyz           )/(c72o1*o));
   xy_TNW =   c1o16*eps_new *((     axyz -     bxyz           )/(c72o1*o));

   xz_E   =   c1o16*eps_new *((            c2o1*bxyz           )/(c27o1*o));
   xz_N   = -(c1o16*eps_new *((            c4o1*bxyz           )/(c27o1*o)));
   xz_T   =   c1o16*eps_new *((            c2o1*bxyz           )/(c27o1*o));
   xz_NE  = -(c1o16*eps_new *(( c3o1*axyz +     bxyz           )/(c54o1*o)));
   xz_SE  = -(c1o16*eps_new *((-c3o1*axyz +     bxyz           )/(c54o1*o)));
   xz_TE  =   c1o16*eps_new *((                bxyz           )/(c27o1*o));
   xz_BE  =   c1o16*eps_new *((                bxyz           )/(c27o1*o));
   xz_TN  = -(c1o16*eps_new *((                bxyz + c3o1*cxyz)/(c54o1*o)));
   xz_BN  = -(c1o16*eps_new *((                bxyz - c3o1*cxyz)/(c54o1*o)));
   //xz_ZERO=   c1o16*eps_new;
   xz_TNE = -(c1o16*eps_new *((     axyz            +     cxyz)/(c72o1*o)));
   xz_TSW =   c1o16*eps_new *((-    axyz            +     cxyz)/(c72o1*o));
   xz_TSE =   c1o16*eps_new *((     axyz            +     cxyz)/(c72o1*o));
   xz_TNW =   c1o16*eps_new *((     axyz            -     cxyz)/(c72o1*o));

   yz_E   = -(c1o16*eps_new *(( c4o1*axyz                      )/(c27o1*o)));
   yz_N   =   c1o16*eps_new *(( c2o1*axyz                      )/(c27o1*o));
   yz_T   =   c1o16*eps_new *(( c2o1*axyz                      )/(c27o1*o));
   yz_NE  = -(c1o16*eps_new *((     axyz + c3o1*bxyz           )/(c54o1*o)));
   yz_SE  = -(c1o16*eps_new *((     axyz - c3o1*bxyz           )/(c54o1*o)));
   yz_TE  = -(c1o16*eps_new *((     axyz            + c3o1*cxyz)/(c54o1*o)));
   yz_BE  = -(c1o16*eps_new *((     axyz            - c3o1*cxyz)/(c54o1*o)));
   yz_TN  =   c1o16*eps_new *((     axyz                      )/(c27o1*o));
   yz_BN  =   c1o16*eps_new *((     axyz                      )/(c27o1*o));
   //yz_ZERO=   c1o16*eps_new;
   yz_TNE = -(c1o16*eps_new *((                bxyz +     cxyz)/(c72o1*o)));
   yz_TSW =   c1o16*eps_new *((          -     bxyz +     cxyz)/(c72o1*o));
   yz_TSE =   c1o16*eps_new *((                bxyz -     cxyz)/(c72o1*o));
   yz_TNW =   c1o16*eps_new *((                bxyz +     cxyz)/(c72o1*o));
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::calcInterpolatedNode(real* f, real  /*omega*/, real  /*x*/, real  /*y*/, real  /*z*/, real press, real xs, real ys, real zs)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   real rho  = press ;//+ (2.*axx*x+axy*y+axz*z+axyz*y*z+ax + 2.*byy*y+bxy*x+byz*z+bxyz*x*z+by + 2.*czz*z+cxz*x+cyz*y+cxyz*x*y+cz)/3.;
   real vx1  = a0 + c1o4*( xs*ax + ys*ay + zs*az) + c1o16*(axx + xs*ys*axy + xs*zs*axz + ayy + ys*zs*ayz + azz) + c1o64*(xs*ys*zs*axyz);
   real vx2  = b0 + c1o4*( xs*bx + ys*by + zs*bz) + c1o16*(bxx + xs*ys*bxy + xs*zs*bxz + byy + ys*zs*byz + bzz) + c1o64*(xs*ys*zs*bxyz);
   real vx3  = c0 + c1o4*( xs*cx + ys*cy + zs*cz) + c1o16*(cxx + xs*ys*cxy + xs*zs*cxz + cyy + ys*zs*cyz + czz) + c1o64*(xs*ys*zs*cxyz);

   //////////////////////////////////////////////////////////////////////////
   //DRAFT
   //vx1 -= forcingF*0.5;
   //////////////////////////////////////////////////////////////////////////

   real feq[ENDF+1];
   D3Q27System::calcIncompFeq(feq,rho,vx1,vx2,vx3);

   f[dP00]    = f_E    + xs*x_E    + ys*y_E    + zs*z_E    + xs*ys*xy_E    + xs*zs*xz_E    + ys*zs*yz_E    + feq[dP00];
   f[dM00]    = f_E    + xs*x_E    + ys*y_E    + zs*z_E    + xs*ys*xy_E    + xs*zs*xz_E    + ys*zs*yz_E    + feq[dM00];
   f[d0P0]    = f_N    + xs*x_N    + ys*y_N    + zs*z_N    + xs*ys*xy_N    + xs*zs*xz_N    + ys*zs*yz_N    + feq[d0P0];
   f[d0M0]    = f_N    + xs*x_N    + ys*y_N    + zs*z_N    + xs*ys*xy_N    + xs*zs*xz_N    + ys*zs*yz_N    + feq[d0M0];
   f[d00P]    = f_T    + xs*x_T    + ys*y_T    + zs*z_T    + xs*ys*xy_T    + xs*zs*xz_T    + ys*zs*yz_T    + feq[d00P];
   f[d00M]    = f_T    + xs*x_T    + ys*y_T    + zs*z_T    + xs*ys*xy_T    + xs*zs*xz_T    + ys*zs*yz_T    + feq[d00M];
   f[dPP0]   = f_NE   + xs*x_NE   + ys*y_NE   + zs*z_NE   + xs*ys*xy_NE   + xs*zs*xz_NE   + ys*zs*yz_NE   + feq[dPP0];
   f[dMM0]   = f_NE   + xs*x_NE   + ys*y_NE   + zs*z_NE   + xs*ys*xy_NE   + xs*zs*xz_NE   + ys*zs*yz_NE   + feq[dMM0];
   f[dPM0]   = f_SE   + xs*x_SE   + ys*y_SE   + zs*z_SE   + xs*ys*xy_SE   + xs*zs*xz_SE   + ys*zs*yz_SE   + feq[dPM0];
   f[dMP0]   = f_SE   + xs*x_SE   + ys*y_SE   + zs*z_SE   + xs*ys*xy_SE   + xs*zs*xz_SE   + ys*zs*yz_SE   + feq[dMP0];
   f[dP0P]   = f_TE   + xs*x_TE   + ys*y_TE   + zs*z_TE   + xs*ys*xy_TE   + xs*zs*xz_TE   + ys*zs*yz_TE   + feq[dP0P];
   f[dM0M]   = f_TE   + xs*x_TE   + ys*y_TE   + zs*z_TE   + xs*ys*xy_TE   + xs*zs*xz_TE   + ys*zs*yz_TE   + feq[dM0M];
   f[dP0M]   = f_BE   + xs*x_BE   + ys*y_BE   + zs*z_BE   + xs*ys*xy_BE   + xs*zs*xz_BE   + ys*zs*yz_BE   + feq[dP0M];
   f[dM0P]   = f_BE   + xs*x_BE   + ys*y_BE   + zs*z_BE   + xs*ys*xy_BE   + xs*zs*xz_BE   + ys*zs*yz_BE   + feq[dM0P];
   f[d0PP]   = f_TN   + xs*x_TN   + ys*y_TN   + zs*z_TN   + xs*ys*xy_TN   + xs*zs*xz_TN   + ys*zs*yz_TN   + feq[d0PP];
   f[d0MM]   = f_TN   + xs*x_TN   + ys*y_TN   + zs*z_TN   + xs*ys*xy_TN   + xs*zs*xz_TN   + ys*zs*yz_TN   + feq[d0MM];
   f[d0PM]   = f_BN   + xs*x_BN   + ys*y_BN   + zs*z_BN   + xs*ys*xy_BN   + xs*zs*xz_BN   + ys*zs*yz_BN   + feq[d0PM];
   f[d0MP]   = f_BN   + xs*x_BN   + ys*y_BN   + zs*z_BN   + xs*ys*xy_BN   + xs*zs*xz_BN   + ys*zs*yz_BN   + feq[d0MP];
   f[dPPP]  = f_TNE  + xs*x_TNE  + ys*y_TNE  + zs*z_TNE  + xs*ys*xy_TNE  + xs*zs*xz_TNE  + ys*zs*yz_TNE  + feq[dPPP];
   f[dMMP]  = f_TSW  + xs*x_TSW  + ys*y_TSW  + zs*z_TSW  + xs*ys*xy_TSW  + xs*zs*xz_TSW  + ys*zs*yz_TSW  + feq[dMMP];
   f[dPMP]  = f_TSE  + xs*x_TSE  + ys*y_TSE  + zs*z_TSE  + xs*ys*xy_TSE  + xs*zs*xz_TSE  + ys*zs*yz_TSE  + feq[dPMP];
   f[dMPP]  = f_TNW  + xs*x_TNW  + ys*y_TNW  + zs*z_TNW  + xs*ys*xy_TNW  + xs*zs*xz_TNW  + ys*zs*yz_TNW  + feq[dMPP];
   f[dPPM]  = f_TSW  + xs*x_TSW  + ys*y_TSW  + zs*z_TSW  + xs*ys*xy_TSW  + xs*zs*xz_TSW  + ys*zs*yz_TSW  + feq[dPPM];
   f[dMMM]  = f_TNE  + xs*x_TNE  + ys*y_TNE  + zs*z_TNE  + xs*ys*xy_TNE  + xs*zs*xz_TNE  + ys*zs*yz_TNE  + feq[dMMM];
   f[dPMM]  = f_TNW  + xs*x_TNW  + ys*y_TNW  + zs*z_TNW  + xs*ys*xy_TNW  + xs*zs*xz_TNW  + ys*zs*yz_TNW  + feq[dPMM];
   f[dMPM]  = f_TSE  + xs*x_TSE  + ys*y_TSE  + zs*z_TSE  + xs*ys*xy_TSE  + xs*zs*xz_TSE  + ys*zs*yz_TSE  + feq[dMPM];
   f[d000] = f_ZERO + xs*x_ZERO + ys*y_ZERO + zs*z_ZERO                                                 + feq[d000];
}
//////////////////////////////////////////////////////////////////////////
//Position SWB -0.25, -0.25, -0.25
real IncompressibleOffsetInterpolator::calcPressBSW()
{
    using namespace vf::basics::constant;

   return   press_SWT * (c9o64 + c3o16 * xoff + c3o16 * yoff - c9o16 * zoff) +
      press_NWT * (c3o64 + c1o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c3o16 * xoff + c1o16 * yoff - c3o16 * zoff) +
      press_NET * (c1o64 - c1o16 * xoff - c1o16 * yoff - c1o16 * zoff) +
      press_NEB * (c3o64 - c3o16 * xoff - c3o16 * yoff + c1o16 * zoff) +
      press_NWB * (c9o64 + c3o16 * xoff - c9o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c9o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c27o64 + c9o16 * xoff + c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SWT -0.25, -0.25, 0.25
real IncompressibleOffsetInterpolator::calcPressTSW()
{
    using namespace vf::basics::constant;

   return   press_SWT * (c27o64 + c9o16 * xoff + c9o16 * yoff - c9o16 * zoff) +
      press_NWT * (c9o64 + c3o16 * xoff - c9o16 * yoff - c3o16 * zoff) +
      press_SET * (c9o64 - c9o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_NET * (c3o64 - c3o16 * xoff - c3o16 * yoff - c1o16 * zoff) +
      press_NEB * (c1o64 - c1o16 * xoff - c1o16 * yoff + c1o16 * zoff) +
      press_NWB * (c3o64 + c1o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c3o64 - c3o16 * xoff + c1o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c3o16 * xoff + c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SET 0.25, -0.25, 0.25
real IncompressibleOffsetInterpolator::calcPressTSE()
{
    using namespace vf::basics::constant;

   return   press_SET * (c27o64 - c9o16 * xoff + c9o16 * yoff - c9o16 * zoff) +
      press_NET * (c9o64 - c3o16 * xoff - c9o16 * yoff - c3o16 * zoff) +
      press_SWT * (c9o64 + c9o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_NWT * (c3o64 + c3o16 * xoff - c3o16 * yoff - c1o16 * zoff) +
      press_NWB * (c1o64 + c1o16 * xoff - c1o16 * yoff + c1o16 * zoff) +
      press_NEB * (c3o64 - c1o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c3o64 + c3o16 * xoff + c1o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c3o16 * xoff + c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SEB 0.25, -0.25, -0.25
real IncompressibleOffsetInterpolator::calcPressBSE()
{
    using namespace vf::basics::constant;

   return   press_SET * (c9o64 - c3o16 * xoff + c3o16 * yoff - c9o16 * zoff) +
      press_NET * (c3o64 - c1o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c3o16 * xoff + c1o16 * yoff - c3o16 * zoff) +
      press_NWT * (c1o64 + c1o16 * xoff - c1o16 * yoff - c1o16 * zoff) +
      press_NWB * (c3o64 + c3o16 * xoff - c3o16 * yoff + c1o16 * zoff) +
      press_NEB * (c9o64 - c3o16 * xoff - c9o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c9o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c27o64 - c9o16 * xoff + c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWB -0.25, 0.25, -0.25
real IncompressibleOffsetInterpolator::calcPressBNW()
{
    using namespace vf::basics::constant;

   return   press_NWT * (c9o64 + c3o16 * xoff - c3o16 * yoff - c9o16 * zoff) +
      press_NET * (c3o64 - c3o16 * xoff - c1o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c1o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c1o64 - c1o16 * xoff + c1o16 * yoff - c1o16 * zoff) +
      press_SEB * (c3o64 - c3o16 * xoff + c3o16 * yoff + c1o16 * zoff) +
      press_NEB * (c9o64 - c9o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SWB * (c9o64 + c3o16 * xoff + c9o16 * yoff + c3o16 * zoff) +
      press_NWB * (c27o64 + c9o16 * xoff - c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWT -0.25, 0.25, 0.25
real IncompressibleOffsetInterpolator::calcPressTNW()
{
    using namespace vf::basics::constant;

   return   press_NWT * (c27o64 + c9o16 * xoff - c9o16 * yoff - c9o16 * zoff) +
      press_NET * (c9o64 - c9o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c9o64 + c3o16 * xoff + c9o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c3o16 * xoff + c3o16 * yoff - c1o16 * zoff) +
      press_SEB * (c1o64 - c1o16 * xoff + c1o16 * yoff + c1o16 * zoff) +
      press_NEB * (c3o64 - c3o16 * xoff - c1o16 * yoff + c3o16 * zoff) +
      press_SWB * (c3o64 + c1o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_NWB * (c9o64 + c3o16 * xoff - c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NET 0.25, 0.25, 0.25
real IncompressibleOffsetInterpolator::calcPressTNE()
{
    using namespace vf::basics::constant;

   return   press_NET * (c27o64 - c9o16 * xoff - c9o16 * yoff - c9o16 * zoff) +
      press_NWT * (c9o64 + c9o16 * xoff - c3o16 * yoff - c3o16 * zoff) +
      press_SET * (c9o64 - c3o16 * xoff + c9o16 * yoff - c3o16 * zoff) +
      press_SWT * (c3o64 + c3o16 * xoff + c3o16 * yoff - c1o16 * zoff) +
      press_SWB * (c1o64 + c1o16 * xoff + c1o16 * yoff + c1o16 * zoff) +
      press_NWB * (c3o64 + c3o16 * xoff - c1o16 * yoff + c3o16 * zoff) +
      press_SEB * (c3o64 - c1o16 * xoff + c3o16 * yoff + c3o16 * zoff) +
      press_NEB * (c9o64 - c3o16 * xoff - c3o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NEB 0.25, 0.25, -0.25
real IncompressibleOffsetInterpolator::calcPressBNE()
{
    using namespace vf::basics::constant;

   return   press_NET * (c9o64 - c3o16 * xoff - c3o16 * yoff - c9o16 * zoff) +
      press_NWT * (c3o64 + c3o16 * xoff - c1o16 * yoff - c3o16 * zoff) +
      press_SET * (c3o64 - c1o16 * xoff + c3o16 * yoff - c3o16 * zoff) +
      press_SWT * (c1o64 + c1o16 * xoff + c1o16 * yoff - c1o16 * zoff) +
      press_SWB * (c3o64 + c3o16 * xoff + c3o16 * yoff + c1o16 * zoff) +
      press_NWB * (c9o64 + c9o16 * xoff - c3o16 * yoff + c3o16 * zoff) +
      press_SEB * (c9o64 - c3o16 * xoff + c9o16 * yoff + c3o16 * zoff) +
      press_NEB * (c27o64 - c9o16 * xoff - c9o16 * yoff + c9o16 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position C 0.0, 0.0, 0.0
void IncompressibleOffsetInterpolator::calcInterpolatedNodeFC(real* f, real omega)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;
   using namespace vf::basics::constant;

   real press  =  press_NET * (c4o32 - c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) +
      press_NWT * (c4o32 + c1o4 * xoff - c1o4 * yoff - c1o4 * zoff) +
      press_SET * (c4o32 - c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) +
      press_SWT * (c4o32 + c1o4 * xoff + c1o4 * yoff - c1o4 * zoff) +
      press_NEB * (c4o32 - c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) +
      press_NWB * (c4o32 + c1o4 * xoff - c1o4 * yoff + c1o4 * zoff) +
      press_SEB * (c4o32 - c1o4 * xoff + c1o4 * yoff + c1o4 * zoff) +
      press_SWB * (c4o32 + c1o4 * xoff + c1o4 * yoff + c1o4 * zoff);
   real vx1  = a0;
   real vx2  = b0;
   real vx3  = c0;

   real rho = press ;//+ (ax+by+cz)/3.;

   //////////////////////////////////////////////////////////////////////////
   //DRAFT
   //vx1 -= forcingC*0.5;
   //////////////////////////////////////////////////////////////////////////

   real feq[ENDF+1];
   D3Q27System::calcIncompFeq(feq,rho,vx1,vx2,vx3);

   real eps_new = c2o1;
   real o  = omega;
//   real op = 1.;

   //f_E    = eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op)/(54.*o*op));
   //f_N    = f_E + eps_new *((2.*(ax - by))/(9.*o));
   //f_T    = f_E + eps_new *((2.*(ax - cz))/(9.*o));
   //f_NE   = eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op))/(54.*o*op));
   //f_SE   = f_NE + eps_new *((  ay + bx )/(9.*o));
   //f_TE   = eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op))/(54.*o*op));
   //f_BE   = f_TE + eps_new *((  az + cx )/(9.*o));
   //f_TN   = eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op)/(54.*o*op));
   //f_BN   = f_TN + eps_new *((  bz + cy )/(9.*o));
   //f_ZERO = eps_new *((5.*(ax + by + cz))/(9.*op));
   //f_TNE  = eps_new *(-(ay + az + bx + bz + cx + cy)/(72.*o));
   //f_TSW  = - eps_new *((ay + bx)/(36.*o)) - f_TNE;
   //f_TSE  = - eps_new *((az + cx)/(36.*o)) - f_TNE;
   //f_TNW  = - eps_new *((bz + cy)/(36.*o)) - f_TNE;

   f_E = eps_new*((c2o1*(-c2o1*ax + by + cz-kxxMzzAverage-kxxMyyAverage))/(c27o1*o));
   f_N = eps_new*((c2o1*(ax - c2o1*by + cz+c2o1*kxxMyyAverage-kxxMzzAverage))/(c27o1*o));
   f_T = eps_new*((c2o1*(ax + by - c2o1*cz-kxxMyyAverage+c2o1*kxxMzzAverage))/(c27o1*o));
   f_NE = eps_new*(-(ax + c3o1*ay + c3o1*bx + by - c2o1*cz+c2o1*kxxMyyAverage-kxxMyyAverage+c3o1*kxyAverage)/(c54o1*o));
   f_SE = eps_new*(-(ax - c3o1*ay - c3o1*bx + by - c2o1*cz+c2o1*kxxMyyAverage-kxxMyyAverage-c3o1*kxyAverage)/(c54o1*o));
   f_TE = eps_new*(-(ax + c3o1*az - c2o1*by + c3o1*cx + cz+c2o1*kxxMyyAverage-kxxMzzAverage+c3o1*kxzAverage)/(c54o1*o));
   f_BE = eps_new*(-(ax - c3o1*az - c2o1*by - c3o1*cx + cz+c2o1*kxxMyyAverage-kxxMzzAverage-c3o1*kxzAverage)/(c54o1*o));
   f_TN = eps_new*(-(-c2o1*ax + by + c3o1*bz + c3o1*cy + cz-kxxMyyAverage-kxxMzzAverage+c3o1*kyzAverage)/(c54o1*o));
   f_BN = eps_new*(-(-c2o1*ax + by - c3o1*bz - c3o1*cy + cz-kxxMyyAverage-kxxMzzAverage-c3o1*kyzAverage)/(c54o1*o));
   f_ZERO = c0o1;
   f_TNE = eps_new*(-(ay + az + bx + bz + cx + cy+kxyAverage+kxzAverage+kyzAverage)/(c72o1*o));
   f_TSW = eps_new*((-ay + az - bx + bz + cx + cy-kxyAverage+kxzAverage+kyzAverage)/(c72o1*o));
   f_TSE = eps_new*((ay - az + bx + bz - cx + cy+kxyAverage-kxzAverage+kyzAverage)/(c72o1*o));
   f_TNW = eps_new*((ay + az + bx - bz + cx - cy+kxyAverage+kxzAverage-kyzAverage)/(c72o1*o));

   f[dP00]    = f_E    + feq[dP00];
   f[dM00]    = f_E    + feq[dM00];
   f[d0P0]    = f_N    + feq[d0P0];
   f[d0M0]    = f_N    + feq[d0M0];
   f[d00P]    = f_T    + feq[d00P];
   f[d00M]    = f_T    + feq[d00M];
   f[dPP0]   = f_NE   + feq[dPP0];
   f[dMM0]   = f_NE   + feq[dMM0];
   f[dPM0]   = f_SE   + feq[dPM0];
   f[dMP0]   = f_SE   + feq[dMP0];
   f[dP0P]   = f_TE   + feq[dP0P];
   f[dM0M]   = f_TE   + feq[dM0M];
   f[dP0M]   = f_BE   + feq[dP0M];
   f[dM0P]   = f_BE   + feq[dM0P];
   f[d0PP]   = f_TN   + feq[d0PP];
   f[d0MM]   = f_TN   + feq[d0MM];
   f[d0PM]   = f_BN   + feq[d0PM];
   f[d0MP]   = f_BN   + feq[d0MP];
   f[dPPP]  = f_TNE  + feq[dPPP];
   f[dMPP]  = f_TNW  + feq[dMPP];
   f[dPMP]  = f_TSE  + feq[dPMP];
   f[dMMP]  = f_TSW  + feq[dMMP];
   f[dPPM]  = f_TSW  + feq[dPPM];
   f[dMPM]  = f_TSE  + feq[dMPM];
   f[dPMM]  = f_TNW  + feq[dPMM];
   f[dMMM]  = f_TNE  + feq[dMMM];
   f[d000] = f_ZERO + feq[d000];
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::calcInterpolatedVelocity(real x, real y, real z, real& vx1, real& vx2, real& vx3)
{
    vx1  = a0 + ax*x + ay*y + az*z + axx*x*x + ayy*y*y + azz*z*z + axy*x*y + axz*x*z + ayz*y*z+axyz*x*y*z;
    vx2  = b0 + bx*x + by*y + bz*z + bxx*x*x + byy*y*y + bzz*z*z + bxy*x*y + bxz*x*z + byz*y*z+bxyz*x*y*z;
    vx3  = c0 + cx*x + cy*y + cz*z + cxx*x*x + cyy*y*y + czz*z*z + cxy*x*y + cxz*x*z + cyz*y*z+cxyz*x*y*z;
}
//////////////////////////////////////////////////////////////////////////
void IncompressibleOffsetInterpolator::calcInterpolatedShearStress(real x, real y, real z,real& tauxx, real& tauyy, real& tauzz,real& tauxy, real& tauxz, real& tauyz)
{
    using namespace vf::basics::constant;

    tauxx=ax+c2o1*axx*x+axy*y+axz*z+axyz*y*z;
    tauyy=by+c2o1*byy*y+bxy*x+byz*z+bxyz*x*z;
    tauzz=cz+c2o1*czz*z+cxz*x+cyz*y+cxyz*x*y;
    tauxy=c1o2*((ay+c2o1*ayy*y+axy*x+ayz*z+axyz*x*z)+(bx+c2o1*bxx*x+bxy*y+bxz*z+bxyz*y*z));
    tauxz=c1o2*((az+c2o1*azz*z+axz*x+ayz*y+axyz*x*y)+(cx+c2o1*cxx*x+cxy*y+cxz*z+cxyz*y*z));
    tauyz=c1o2*((bz+c2o1*bzz*z+bxz*x+byz*y+bxyz*x*y)+(cy+c2o1*cyy*y+cxy*x+cyz*z+cxyz*x*z));
}
