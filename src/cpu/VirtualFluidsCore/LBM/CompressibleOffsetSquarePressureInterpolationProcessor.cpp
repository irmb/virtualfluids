#include "CompressibleOffsetSquarePressureInterpolationProcessor.h"
#include "D3Q27System.h"

//using namespace UbMath;
using namespace vf::lbm::constant;

CompressibleOffsetSquarePressureInterpolationProcessor::CompressibleOffsetSquarePressureInterpolationProcessor()
    
{
   this->bulkOmegaToOmega = false;
   this->OxxPyyPzzC = c1o1;
   this->OxxPyyPzzF = c1o1;
}
//////////////////////////////////////////////////////////////////////////
CompressibleOffsetSquarePressureInterpolationProcessor::CompressibleOffsetSquarePressureInterpolationProcessor(real omegaC, real omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{
   this->bulkOmegaToOmega = false;
   this->OxxPyyPzzC = c1o1;
   this->OxxPyyPzzF = c1o1;
}
//////////////////////////////////////////////////////////////////////////
CompressibleOffsetSquarePressureInterpolationProcessor::~CompressibleOffsetSquarePressureInterpolationProcessor()
= default;
//////////////////////////////////////////////////////////////////////////
InterpolationProcessorPtr CompressibleOffsetSquarePressureInterpolationProcessor::clone()
{
   InterpolationProcessorPtr iproc = InterpolationProcessorPtr (new CompressibleOffsetSquarePressureInterpolationProcessor(this->omegaC, this->omegaF));
   if (bulkOmegaToOmega)
   {
      dynamicPointerCast<CompressibleOffsetSquarePressureInterpolationProcessor>(iproc)->OxxPyyPzzC = omegaC;
      dynamicPointerCast<CompressibleOffsetSquarePressureInterpolationProcessor>(iproc)->OxxPyyPzzF = omegaF;
   }
   else
   {
      dynamicPointerCast<CompressibleOffsetSquarePressureInterpolationProcessor>(iproc)->OxxPyyPzzC = c1o1;
      dynamicPointerCast<CompressibleOffsetSquarePressureInterpolationProcessor>(iproc)->OxxPyyPzzF = c1o1;
   }
   return iproc;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::setOmegas( real omegaC, real omegaF )
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::setOffsets(real xoff, real yoff, real zoff)
{
   this->xoff = xoff;
   this->yoff = yoff;
   this->zoff = zoff;     
   this->xoff_sq = xoff * xoff;
   this->yoff_sq = yoff * yoff;
   this->zoff_sq = zoff * zoff;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellC, omegaC, 0.5);
   calcInterpolatedNodeCF(icellF.BSW, omegaF, -0.25, -0.25, -0.25, calcPressBSW(), -1, -1, -1);
   calcInterpolatedNodeCF(icellF.BNE, omegaF,  0.25,  0.25, -0.25, calcPressBNE(),  1,  1, -1);
   calcInterpolatedNodeCF(icellF.TNW, omegaF, -0.25,  0.25,  0.25, calcPressTNW(), -1,  1,  1);
   calcInterpolatedNodeCF(icellF.TSE, omegaF,  0.25, -0.25,  0.25, calcPressTSE(),  1, -1,  1);
   calcInterpolatedNodeCF(icellF.BNW, omegaF, -0.25,  0.25, -0.25, calcPressBNW(), -1,  1, -1);
   calcInterpolatedNodeCF(icellF.BSE, omegaF,  0.25, -0.25, -0.25, calcPressBSE(),  1, -1, -1);
   calcInterpolatedNodeCF(icellF.TSW, omegaF, -0.25, -0.25,  0.25, calcPressTSW(), -1, -1,  1);
   calcInterpolatedNodeCF(icellF.TNE, omegaF,  0.25,  0.25,  0.25, calcPressTNE(),  1,  1,  1);
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellF, omegaF, 2.0);
   calcInterpolatedNodeFC(icellC, omegaC);
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::calcMoments(const real* const f, real omega, real& press, real& vx1, real& vx2, real& vx3, 
                                                    real& kxy, real& kyz, real& kxz, real& kxxMyy, real& kxxMzz)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;

   real drho = 0.0;
   D3Q27System::calcCompMacroscopicValues(f,drho,vx1,vx2,vx3);
   
   press = drho; //interpolate rho!

   kxy   = -3.*omega*((((f[DIR_MMP]+f[DIR_PPM])-(f[DIR_MPP]+f[DIR_PMM]))+((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_MPM]+f[DIR_PMP])))+((f[DIR_MM0]+f[DIR_PP0])-(f[DIR_MP0]+f[DIR_PM0]))/(c1o1 + drho)-(vx1*vx2));// might not be optimal MG 25.2.13
   kyz   = -3.*omega*((((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_PMP]+f[DIR_MPM]))+((f[DIR_PMM]+f[DIR_MPP])-(f[DIR_MMP]+f[DIR_PPM])))+((f[DIR_0MM]+f[DIR_0PP])-(f[DIR_0MP]+f[DIR_0PM]))/(c1o1 + drho)-(vx2*vx3));
   kxz   = -3.*omega*((((f[DIR_MPM]+f[DIR_PMP])-(f[DIR_MMP]+f[DIR_PPM]))+((f[DIR_MMM]+f[DIR_PPP])-(f[DIR_PMM]+f[DIR_MPP])))+((f[DIR_M0M]+f[DIR_P0P])-(f[DIR_M0P]+f[DIR_P0M]))/(c1o1 + drho)-(vx1*vx3));
   kxxMyy = -3./2.*omega*((((f[DIR_M0M]+f[DIR_P0P])-(f[DIR_0MM]+f[DIR_0PP]))+((f[DIR_M0P]+f[DIR_P0M])-(f[DIR_0MP]+f[DIR_0PM])))+((f[DIR_M00]+f[DIR_P00])-(f[DIR_0M0]+f[DIR_0P0]))/(c1o1 + drho)-(vx1*vx1-vx2*vx2));
   kxxMzz = -3./2.*omega*((((f[DIR_MP0]+f[DIR_PM0])-(f[DIR_0MM]+f[DIR_0PP]))+((f[DIR_MM0]+f[DIR_PP0])-(f[DIR_0MP]+f[DIR_0PM])))+((f[DIR_M00]+f[DIR_P00])-(f[DIR_00M]+f[DIR_00P]))/(c1o1 + drho)-(vx1*vx1-vx3*vx3));
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::calcInterpolatedCoefficiets(const D3Q27ICell& icell, real omega, real eps_new)
{
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

   //LBMReal dxRho=c1o4*((press_NET-press_SWB)+(press_SET-press_NWB)+(press_NEB-press_SWT)+(press_SEB-press_NWT));
   //LBMReal dyRho=c1o4*((press_NET-press_SWB)-(press_SET-press_NWB)+(press_NEB-press_SWT)-(press_SEB-press_NWT));
   //LBMReal dzRho=c1o4*((press_NET-press_SWB)+(press_SET-press_NWB)-(press_NEB-press_SWT)-(press_SEB-press_NWT));

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
      2.*kxyFromfcNEQ_NEB - 2.*kxyFromfcNEQ_NET - 2.*kxyFromfcNEQ_NWB - 2.*kxyFromfcNEQ_NWT +
      2.*kxyFromfcNEQ_SEB + 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_SWT +
      2.*kxzFromfcNEQ_NEB - 2.*kxzFromfcNEQ_NET + 2.*kxzFromfcNEQ_NWB - 2.*kxzFromfcNEQ_NWT +
      2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_SWB - 2.*kxzFromfcNEQ_SWT +
      8.*vx1_NEB + 8.*vx1_NET + 8.*vx1_NWB + 8.*vx1_NWT + 8.*vx1_SEB +
      8.*vx1_SET + 8.*vx1_SWB + 8.*vx1_SWT + 2.*vx2_NEB + 2.*vx2_NET -
      2.*vx2_NWB - 2.*vx2_NWT - 2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB +
      2.*vx2_SWT - 2.*vx3_NEB + 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/64.;
   b0 = (2.*kxxMyyFromfcNEQ_NEB + 2.*kxxMyyFromfcNEQ_NET + 2.*kxxMyyFromfcNEQ_NWB + 2.*kxxMyyFromfcNEQ_NWT -
      2.*kxxMyyFromfcNEQ_SEB - 2.*kxxMyyFromfcNEQ_SET - 2.*kxxMyyFromfcNEQ_SWB - 2.*kxxMyyFromfcNEQ_SWT -
      kxxMzzFromfcNEQ_NEB - kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_SWT -
      2.*kxyFromfcNEQ_NEB - 2.*kxyFromfcNEQ_NET + 2.*kxyFromfcNEQ_NWB + 2.*kxyFromfcNEQ_NWT -
      2.*kxyFromfcNEQ_SEB - 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_SWT +
      2.*kyzFromfcNEQ_NEB - 2.*kyzFromfcNEQ_NET + 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_NWT +
      2.*kyzFromfcNEQ_SEB - 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_SWB - 2.*kyzFromfcNEQ_SWT +
      2.*vx1_NEB + 2.*vx1_NET - 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB - 2.*vx1_SET + 2.*vx1_SWB + 2.*vx1_SWT +
      8.*vx2_NEB + 8.*vx2_NET + 8.*vx2_NWB + 8.*vx2_NWT +
      8.*vx2_SEB + 8.*vx2_SET + 8.*vx2_SWB + 8.*vx2_SWT -
      2.*vx3_NEB + 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/64.;
   c0 = (kxxMyyFromfcNEQ_NEB - kxxMyyFromfcNEQ_NET + kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT -
      2.*kxxMzzFromfcNEQ_NEB + 2.*kxxMzzFromfcNEQ_NET - 2.*kxxMzzFromfcNEQ_NWB + 2.*kxxMzzFromfcNEQ_NWT -
      2.*kxxMzzFromfcNEQ_SEB + 2.*kxxMzzFromfcNEQ_SET - 2.*kxxMzzFromfcNEQ_SWB + 2.*kxxMzzFromfcNEQ_SWT -
      2.*kxzFromfcNEQ_NEB - 2.*kxzFromfcNEQ_NET + 2.*kxzFromfcNEQ_NWB + 2.*kxzFromfcNEQ_NWT -
      2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_SWB + 2.*kxzFromfcNEQ_SWT -
      2.*kyzFromfcNEQ_NEB - 2.*kyzFromfcNEQ_NET - 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_NWT +
      2.*kyzFromfcNEQ_SEB + 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_SWB + 2.*kyzFromfcNEQ_SWT -
      2.*vx1_NEB + 2.*vx1_NET + 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB + 2.*vx1_SET + 2.*vx1_SWB - 2.*vx1_SWT -
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB - 2.*vx2_SWT +
      8.*vx3_NEB + 8.*vx3_NET + 8.*vx3_NWB + 8.*vx3_NWT +
      8.*vx3_SEB + 8.*vx3_SET + 8.*vx3_SWB + 8.*vx3_SWT)/64.;
   ax = (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT + vx1_SEB + vx1_SET - vx1_SWB - vx1_SWT)/4.;
   bx = (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT + vx2_SEB + vx2_SET - vx2_SWB - vx2_SWT)/4.;
   cx = (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT + vx3_SEB + vx3_SET - vx3_SWB - vx3_SWT)/4.;
   axx= (kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB - kxxMyyFromfcNEQ_NWT +
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET - kxxMzzFromfcNEQ_NWB - kxxMzzFromfcNEQ_NWT +
      kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB - 2.*vx2_NWT -
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB + 2.*vx2_SWT -
      2.*vx3_NEB + 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/16.;
   bxx= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET - kxyFromfcNEQ_NWB - kxyFromfcNEQ_NWT +
      kxyFromfcNEQ_SEB + kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      2.*vx1_NEB - 2.*vx1_NET + 2.*vx1_NWB + 2.*vx1_NWT +
      2.*vx1_SEB + 2.*vx1_SET - 2.*vx1_SWB - 2.*vx1_SWT)/8.;
   cxx= (kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB - kxzFromfcNEQ_NWT +
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB - kxzFromfcNEQ_SWT +
      2.*vx1_NEB - 2.*vx1_NET - 2.*vx1_NWB + 2.*vx1_NWT +
      2.*vx1_SEB - 2.*vx1_SET - 2.*vx1_SWB + 2.*vx1_SWT)/8.;
   ay = (vx1_NEB + vx1_NET + vx1_NWB + vx1_NWT - vx1_SEB - vx1_SET - vx1_SWB - vx1_SWT)/4.;
   by = (vx2_NEB + vx2_NET + vx2_NWB + vx2_NWT - vx2_SEB - vx2_SET - vx2_SWB - vx2_SWT)/4.;
   cy = (vx3_NEB + vx3_NET + vx3_NWB + vx3_NWT - vx3_SEB - vx3_SET - vx3_SWB - vx3_SWT)/4.;
   ayy= (kxyFromfcNEQ_NEB + kxyFromfcNEQ_NET + kxyFromfcNEQ_NWB + kxyFromfcNEQ_NWT -
      kxyFromfcNEQ_SEB - kxyFromfcNEQ_SET - kxyFromfcNEQ_SWB - kxyFromfcNEQ_SWT -
      2.*vx2_NEB - 2.*vx2_NET + 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB + 2.*vx2_SET - 2.*vx2_SWB - 2.*vx2_SWT)/8.;
   byy= (-2.*kxxMyyFromfcNEQ_NEB - 2.*kxxMyyFromfcNEQ_NET - 2.*kxxMyyFromfcNEQ_NWB - 2.*kxxMyyFromfcNEQ_NWT +
      2.*kxxMyyFromfcNEQ_SEB + 2.*kxxMyyFromfcNEQ_SET + 2.*kxxMyyFromfcNEQ_SWB + 2.*kxxMyyFromfcNEQ_SWT +
      kxxMzzFromfcNEQ_NEB + kxxMzzFromfcNEQ_NET + kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_NWT -
      kxxMzzFromfcNEQ_SEB - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_SWT +
      2.*vx1_NEB + 2.*vx1_NET - 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB - 2.*vx1_SET + 2.*vx1_SWB + 2.*vx1_SWT -
      2.*vx3_NEB + 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET + 2.*vx3_SWB - 2.*vx3_SWT)/16.;
   cyy= (kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET + kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB - kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB - kyzFromfcNEQ_SWT +
      2.*vx2_NEB - 2.*vx2_NET + 2.*vx2_NWB - 2.*vx2_NWT -
      2.*vx2_SEB + 2.*vx2_SET - 2.*vx2_SWB + 2.*vx2_SWT)/8.;
   az = (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT - vx1_SEB + vx1_SET - vx1_SWB + vx1_SWT)/4.;
   bz = (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT - vx2_SEB + vx2_SET - vx2_SWB + vx2_SWT)/4.;
   cz = (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT - vx3_SEB + vx3_SET - vx3_SWB + vx3_SWT)/4.;
   azz= (-kxzFromfcNEQ_NEB + kxzFromfcNEQ_NET - kxzFromfcNEQ_NWB + kxzFromfcNEQ_NWT -
      kxzFromfcNEQ_SEB + kxzFromfcNEQ_SET - kxzFromfcNEQ_SWB + kxzFromfcNEQ_SWT +
      2.*vx3_NEB - 2.*vx3_NET - 2.*vx3_NWB + 2.*vx3_NWT +
      2.*vx3_SEB - 2.*vx3_SET - 2.*vx3_SWB + 2.*vx3_SWT)/8.;
   bzz= (-kyzFromfcNEQ_NEB + kyzFromfcNEQ_NET - kyzFromfcNEQ_NWB + kyzFromfcNEQ_NWT -
      kyzFromfcNEQ_SEB + kyzFromfcNEQ_SET - kyzFromfcNEQ_SWB + kyzFromfcNEQ_SWT +
      2.*vx3_NEB - 2.*vx3_NET + 2.*vx3_NWB - 2.*vx3_NWT -
      2.*vx3_SEB + 2.*vx3_SET - 2.*vx3_SWB + 2.*vx3_SWT)/8.;
   czz= (-kxxMyyFromfcNEQ_NEB + kxxMyyFromfcNEQ_NET - kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_NWT -
      kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_SWT +
      2.*kxxMzzFromfcNEQ_NEB - 2.*kxxMzzFromfcNEQ_NET + 2.*kxxMzzFromfcNEQ_NWB - 2.*kxxMzzFromfcNEQ_NWT +
      2.*kxxMzzFromfcNEQ_SEB - 2.*kxxMzzFromfcNEQ_SET + 2.*kxxMzzFromfcNEQ_SWB - 2.*kxxMzzFromfcNEQ_SWT -
      2.*vx1_NEB + 2.*vx1_NET + 2.*vx1_NWB - 2.*vx1_NWT -
      2.*vx1_SEB + 2.*vx1_SET + 2.*vx1_SWB - 2.*vx1_SWT -
      2.*vx2_NEB + 2.*vx2_NET - 2.*vx2_NWB + 2.*vx2_NWT +
      2.*vx2_SEB - 2.*vx2_SET + 2.*vx2_SWB - 2.*vx2_SWT)/16.;
   axy= (vx1_NEB + vx1_NET - vx1_NWB - vx1_NWT - vx1_SEB - vx1_SET + vx1_SWB + vx1_SWT)/2.;
   bxy= (vx2_NEB + vx2_NET - vx2_NWB - vx2_NWT - vx2_SEB - vx2_SET + vx2_SWB + vx2_SWT)/2.;
   cxy= (vx3_NEB + vx3_NET - vx3_NWB - vx3_NWT - vx3_SEB - vx3_SET + vx3_SWB + vx3_SWT)/2.;
   axz= (-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT - vx1_SEB + vx1_SET + vx1_SWB - vx1_SWT)/2.;
   bxz= (-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT - vx2_SEB + vx2_SET + vx2_SWB - vx2_SWT)/2.;
   cxz= (-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT - vx3_SEB + vx3_SET + vx3_SWB - vx3_SWT)/2.;
   ayz= (-vx1_NEB + vx1_NET - vx1_NWB + vx1_NWT + vx1_SEB - vx1_SET + vx1_SWB - vx1_SWT)/2.;
   byz= (-vx2_NEB + vx2_NET - vx2_NWB + vx2_NWT + vx2_SEB - vx2_SET + vx2_SWB - vx2_SWT)/2.;
   cyz= (-vx3_NEB + vx3_NET - vx3_NWB + vx3_NWT + vx3_SEB - vx3_SET + vx3_SWB - vx3_SWT)/2.;
   axyz=-vx1_NEB + vx1_NET + vx1_NWB - vx1_NWT + vx1_SEB - vx1_SET - vx1_SWB + vx1_SWT;
   bxyz=-vx2_NEB + vx2_NET + vx2_NWB - vx2_NWT + vx2_SEB - vx2_SET - vx2_SWB + vx2_SWT;
   cxyz=-vx3_NEB + vx3_NET + vx3_NWB - vx3_NWT + vx3_SEB - vx3_SET - vx3_SWB + vx3_SWT;


   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   kxyAverage       =0;//(kxyFromfcNEQ_SWB+
                       //kxyFromfcNEQ_SWT+
                       //kxyFromfcNEQ_SET+
                       //kxyFromfcNEQ_SEB+
                       //kxyFromfcNEQ_NWB+
                       //kxyFromfcNEQ_NWT+
                       //kxyFromfcNEQ_NET+
                       //kxyFromfcNEQ_NEB)*c1o8-(ay+bx);
   kyzAverage       =0;//(kyzFromfcNEQ_SWB+
                       //kyzFromfcNEQ_SWT+
                       //kyzFromfcNEQ_SET+
                       //kyzFromfcNEQ_SEB+
                       //kyzFromfcNEQ_NWB+
                       //kyzFromfcNEQ_NWT+
                       //kyzFromfcNEQ_NET+
                       //kyzFromfcNEQ_NEB)*c1o8-(bz+cy);
   kxzAverage       =0;//(kxzFromfcNEQ_SWB+
                       //kxzFromfcNEQ_SWT+
                       //kxzFromfcNEQ_SET+
                       //kxzFromfcNEQ_SEB+
                       //kxzFromfcNEQ_NWB+
                       //kxzFromfcNEQ_NWT+
                       //kxzFromfcNEQ_NET+
                       //kxzFromfcNEQ_NEB)*c1o8-(az+cx);
   kxxMyyAverage    =0;//(kxxMyyFromfcNEQ_SWB+
                       //kxxMyyFromfcNEQ_SWT+
                       //kxxMyyFromfcNEQ_SET+
                       //kxxMyyFromfcNEQ_SEB+
                       //kxxMyyFromfcNEQ_NWB+
                       //kxxMyyFromfcNEQ_NWT+
                       //kxxMyyFromfcNEQ_NET+
                       //kxxMyyFromfcNEQ_NEB)*c1o8-(ax-by);
   kxxMzzAverage    =0;//(kxxMzzFromfcNEQ_SWB+
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
   ax = ax + 2. * xoff * axx + yoff * axy + zoff * axz + yoff*zoff*axyz;
   ay = ay + 2. * yoff * ayy + xoff * axy + zoff * ayz + xoff*zoff*axyz;
   az = az + 2. * zoff * azz + xoff * axz + yoff * ayz + xoff*yoff*axyz;
   b0 = b0 + xoff * bx + yoff * by + zoff * bz + xoff_sq * bxx + yoff_sq * byy + zoff_sq * bzz + xoff*yoff*bxy + xoff*zoff*bxz + yoff*zoff*byz + xoff*yoff*zoff*bxyz;
   bx = bx + 2. * xoff * bxx + yoff * bxy + zoff * bxz + yoff*zoff*bxyz;
   by = by + 2. * yoff * byy + xoff * bxy + zoff * byz + xoff*zoff*bxyz;
   bz = bz + 2. * zoff * bzz + xoff * bxz + yoff * byz + xoff*yoff*bxyz;
   c0 = c0 + xoff * cx + yoff * cy + zoff * cz + xoff_sq * cxx + yoff_sq * cyy + zoff_sq * czz + xoff*yoff*cxy + xoff*zoff*cxz + yoff*zoff*cyz + xoff*yoff*zoff*cxyz;
   cx = cx + 2. * xoff * cxx + yoff * cxy + zoff * cxz + yoff*zoff*cxyz;
   cy = cy + 2. * yoff * cyy + xoff * cxy + zoff * cyz + xoff*zoff*cxyz;
   cz = cz + 2. * zoff * czz + xoff * cxz + yoff * cyz + xoff*yoff*cxyz;
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

   f_E = eps_new*((2*(-2*ax + by + cz-kxxMzzAverage-kxxMyyAverage))/(27.*o));
   f_N = eps_new*((2*(ax - 2*by + cz+2*kxxMyyAverage-kxxMzzAverage))/(27.*o));
   f_T = eps_new*((2*(ax + by - 2*cz-kxxMyyAverage+2*kxxMzzAverage))/(27.*o));
   f_NE = eps_new*(-(ax + 3*ay + 3*bx + by - 2*cz+2*kxxMyyAverage-kxxMyyAverage+3*kxyAverage)/(54.*o));
   f_SE = eps_new*(-(ax - 3*ay - 3*bx + by - 2*cz+2*kxxMyyAverage-kxxMyyAverage-3*kxyAverage)/(54.*o));
   f_TE = eps_new*(-(ax + 3*az - 2*by + 3*cx + cz+2*kxxMyyAverage-kxxMzzAverage+3*kxzAverage)/(54.*o));
   f_BE = eps_new*(-(ax - 3*az - 2*by - 3*cx + cz+2*kxxMyyAverage-kxxMzzAverage-3*kxzAverage)/(54.*o));
   f_TN = eps_new*(-(-2*ax + by + 3*bz + 3*cy + cz-kxxMyyAverage-kxxMzzAverage+3*kyzAverage)/(54.*o));
   f_BN = eps_new*(-(-2*ax + by - 3*bz - 3*cy + cz-kxxMyyAverage-kxxMzzAverage-3*kyzAverage)/(54.*o));
   f_ZERO = 0.;
   f_TNE = eps_new*(-(ay + az + bx + bz + cx + cy+kxyAverage+kxzAverage+kyzAverage)/(72.*o));
   f_TSW = eps_new*((-ay + az - bx + bz + cx + cy-kxyAverage+kxzAverage+kyzAverage)/(72.*o));
   f_TSE = eps_new*((ay - az + bx + bz - cx + cy+kxyAverage-kxzAverage+kyzAverage)/(72.*o));
   f_TNW = eps_new*((ay + az + bx - bz + cx - cy+kxyAverage+kxzAverage-kyzAverage)/(72.*o));

   x_E = 0.25*eps_new*((2*(-4*axx + bxy + cxz))/(27.*o));
   x_N = 0.25*eps_new*((2*(2*axx - 2*bxy + cxz))/(27.*o));
   x_T = 0.25*eps_new*((2*(2*axx + bxy - 2*cxz))/(27.*o));
   x_NE = 0.25*eps_new*(-((2*axx + 3*axy + 6*bxx + bxy - 2*cxz))/(54.*o));
   x_SE = 0.25*eps_new*(-((2*axx - 3*axy - 6*bxx + bxy - 2*cxz))/(54.*o));
   x_TE = 0.25*eps_new*(-((2*axx + 3*axz - 2*bxy + 6*cxx + cxz))/(54.*o));
   x_BE = 0.25*eps_new*(-((2*axx - 3*axz - 2*bxy - 6*cxx + cxz))/(54.*o));
   x_TN = 0.25*eps_new*(-((-4*axx + bxy + 3*bxz + 3*cxy + cxz))/(54.*o));
   x_BN = 0.25*eps_new*(-((-4*axx + bxy - 3*bxz - 3*cxy + cxz))/(54.*o));
   x_ZERO = 0.;
   x_TNE = 0.25*eps_new*(-((axy + axz + 2*bxx + bxz + 2*cxx + cxy))/(72.*o));
   x_TSW = 0.25*eps_new*(((-axy + axz - 2*bxx + bxz + 2*cxx + cxy))/(72.*o));
   x_TSE = 0.25*eps_new*(((axy - axz + 2*bxx + bxz - 2*cxx + cxy))/(72.*o));
   x_TNW = 0.25*eps_new*(((axy + axz + 2*bxx - bxz + 2*cxx - cxy))/(72.*o));

   y_E = 0.25*eps_new*(2*(-2*axy + 2*byy + cyz))/(27.*o);
   y_N = 0.25*eps_new*(2*(axy - 4*byy + cyz))/(27.*o);
   y_T = 0.25*eps_new*(2*(axy + 2*byy - 2*cyz))/(27.*o);
   y_NE = 0.25*eps_new*(-((axy + 6*ayy + 3*bxy + 2*byy - 2*cyz))/(54.*o));
   y_SE = 0.25*eps_new*(-((axy - 6*ayy - 3*bxy + 2*byy - 2*cyz))/(54.*o));
   y_TE = 0.25*eps_new*(-((axy + 3*ayz - 4*byy + 3*cxy + cyz))/(54.*o));
   y_BE = 0.25*eps_new*(-((axy - 3*ayz - 4*byy - 3*cxy + cyz))/(54.*o));
   y_TN = 0.25*eps_new*(-((-2*axy + 2*byy + 3*byz + 6*cyy + cyz))/(54.*o));
   y_BN = 0.25*eps_new*(-((-2*axy + 2*byy - 3*byz - 6*cyy + cyz))/(54.*o));
   y_ZERO = 0.;
   y_TNE = 0.25*eps_new*(-((2*ayy + ayz + bxy + byz + cxy + 2*cyy))/(72.*o));
   y_TSW = 0.25*eps_new*(((-2*ayy + ayz - bxy + byz + cxy + 2*cyy))/(72.*o));
   y_TSE = 0.25*eps_new*(((2*ayy - ayz + bxy + byz - cxy + 2*cyy))/(72.*o));
   y_TNW = 0.25*eps_new*(((2*ayy + ayz + bxy - byz + cxy - 2*cyy))/(72.*o));

   z_E = 0.25*eps_new*((2*(-2*axz + byz + 2*czz))/(27.*o));
   z_N = 0.25*eps_new*((2*(axz - 2*byz + 2*czz))/(27.*o));
   z_T = 0.25*eps_new*((2*(axz + byz - 4*czz))/(27.*o));
   z_NE = 0.25*eps_new*(-((axz + 3*ayz + 3*bxz + byz - 4*czz))/(54.*o));
   z_SE = 0.25*eps_new*(-((axz - 3*ayz - 3*bxz + byz - 4*czz))/(54.*o));
   z_TE = 0.25*eps_new*(-((axz + 6*azz - 2*byz + 3*cxz + 2*czz))/(54.*o));
   z_BE = 0.25*eps_new*(-((axz - 6*azz - 2*byz - 3*cxz + 2*czz))/(54.*o));
   z_TN = 0.25*eps_new*(-((-2*axz + byz + 6*bzz + 3*cyz + 2*czz))/(54.*o));
   z_BN = 0.25*eps_new*(-((-2*axz + byz - 6*bzz - 3*cyz + 2*czz))/(54.*o));
   z_ZERO = 0.;
   z_TNE = 0.25*eps_new*(-((ayz + 2*azz + bxz + 2*bzz + cxz + cyz))/(72.*o));
   z_TSW = 0.25*eps_new*(((-ayz + 2*azz - bxz + 2*bzz + cxz + cyz))/(72.*o));
   z_TSE = 0.25*eps_new*(((ayz - 2*azz + bxz + 2*bzz - cxz + cyz))/(72.*o));
   z_TNW = 0.25*eps_new*(((ayz + 2*azz + bxz - 2*bzz + cxz - cyz))/(72.*o));

   xy_E   =   0.0625*eps_new *((                       2.*cxyz)/(27.*o));
   xy_N   =   0.0625*eps_new *((                       2.*cxyz)/(27.*o));
   xy_T   = -(0.0625*eps_new *((                       4.*cxyz)/(27.*o)));
   xy_NE  =   0.0625*eps_new *(                            cxyz /(27.*o));
   xy_SE  =   0.0625*eps_new *(                            cxyz /(27.*o));
   xy_TE  = -(0.0625*eps_new *(( 3.*axyz            +     cxyz)/(54.*o)));
   xy_BE  = -(0.0625*eps_new *((-3.*axyz            +     cxyz)/(54.*o)));
   xy_TN  = -(0.0625*eps_new *((            3.*bxyz +     cxyz)/(54.*o)));
   xy_BN  = -(0.0625*eps_new *((          - 3.*bxyz +     cxyz)/(54.*o)));
   //xy_ZERO=   0.0625*eps_new;
   xy_TNE = -(0.0625*eps_new *((     axyz +     bxyz           )/(72.*o)));
   xy_TSW =   0.0625*eps_new *((     axyz +     bxyz           )/(72.*o));
   xy_TSE =   0.0625*eps_new *((-    axyz +     bxyz           )/(72.*o));
   xy_TNW =   0.0625*eps_new *((     axyz -     bxyz           )/(72.*o));

   xz_E   =   0.0625*eps_new *((            2.*bxyz           )/(27.*o));
   xz_N   = -(0.0625*eps_new *((            4.*bxyz           )/(27.*o)));
   xz_T   =   0.0625*eps_new *((            2.*bxyz           )/(27.*o));
   xz_NE  = -(0.0625*eps_new *(( 3.*axyz +     bxyz           )/(54.*o)));
   xz_SE  = -(0.0625*eps_new *((-3.*axyz +     bxyz           )/(54.*o)));
   xz_TE  =   0.0625*eps_new *((                bxyz           )/(27.*o));
   xz_BE  =   0.0625*eps_new *((                bxyz           )/(27.*o));
   xz_TN  = -(0.0625*eps_new *((                bxyz + 3.*cxyz)/(54.*o)));
   xz_BN  = -(0.0625*eps_new *((                bxyz - 3.*cxyz)/(54.*o)));
   //xz_ZERO=   0.0625*eps_new;
   xz_TNE = -(0.0625*eps_new *((     axyz            +     cxyz)/(72.*o)));
   xz_TSW =   0.0625*eps_new *((-    axyz            +     cxyz)/(72.*o));
   xz_TSE =   0.0625*eps_new *((     axyz            +     cxyz)/(72.*o));
   xz_TNW =   0.0625*eps_new *((     axyz            -     cxyz)/(72.*o));

   yz_E   = -(0.0625*eps_new *(( 4.*axyz                      )/(27.*o)));
   yz_N   =   0.0625*eps_new *(( 2.*axyz                      )/(27.*o));
   yz_T   =   0.0625*eps_new *(( 2.*axyz                      )/(27.*o));
   yz_NE  = -(0.0625*eps_new *((     axyz + 3.*bxyz           )/(54.*o)));
   yz_SE  = -(0.0625*eps_new *((     axyz - 3.*bxyz           )/(54.*o)));
   yz_TE  = -(0.0625*eps_new *((     axyz            + 3.*cxyz)/(54.*o)));
   yz_BE  = -(0.0625*eps_new *((     axyz            - 3.*cxyz)/(54.*o)));
   yz_TN  =   0.0625*eps_new *((     axyz                      )/(27.*o));
   yz_BN  =   0.0625*eps_new *((     axyz                      )/(27.*o));
   //yz_ZERO=   0.0625*eps_new;
   yz_TNE = -(0.0625*eps_new *((                bxyz +     cxyz)/(72.*o)));
   yz_TSW =   0.0625*eps_new *((          -     bxyz +     cxyz)/(72.*o));
   yz_TSE =   0.0625*eps_new *((                bxyz -     cxyz)/(72.*o));
   yz_TNW =   0.0625*eps_new *((                bxyz +     cxyz)/(72.*o));
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::calcInterpolatedNodeCF(real* f, real omega, real x, real y, real z, real press, real xs, real ys, real zs)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;

   real eps_new = 0.5;
   real o = omega;
   //bulk viscosity
   real oP = OxxPyyPzzF;

   real rho  = press ;//+ (2.*axx*x+axy*y+axz*z+axyz*y*z+ax + 2.*byy*y+bxy*x+byz*z+bxyz*x*z+by + 2.*czz*z+cxz*x+cyz*y+cxyz*x*y+cz)/3.;

   real laplaceRho = (xoff!=0.0 || yoff!=0.0 || zoff!= 0.0) ? 0.0 :(-3.0*(by*by+ax*ax+cz*cz)-6.0*(ay*bx+bz*cy+az*cx))*(1.0+rho);

   rho=rho+laplaceRho*(3.0/16.0);

   real vx1  = a0 + 0.25*( xs*ax + ys*ay + zs*az) + 0.0625*(axx + xs*ys*axy + xs*zs*axz + ayy + ys*zs*ayz + azz) + 0.015625*(xs*ys*zs*axyz);
   real vx2  = b0 + 0.25*( xs*bx + ys*by + zs*bz) + 0.0625*(bxx + xs*ys*bxy + xs*zs*bxz + byy + ys*zs*byz + bzz) + 0.015625*(xs*ys*zs*bxyz);
   real vx3  = c0 + 0.25*( xs*cx + ys*cy + zs*cz) + 0.0625*(cxx + xs*ys*cxy + xs*zs*cxz + cyy + ys*zs*cyz + czz) + 0.015625*(xs*ys*zs*cxyz);

   real mfcbb = c0o1;
   real mfabb = c0o1;
   real mfbcb = c0o1;
   real mfbab = c0o1;
   real mfbbc = c0o1;
   real mfbba = c0o1;
   real mfccb = c0o1;
   real mfaab = c0o1;
   real mfcab = c0o1;
   real mfacb = c0o1;
   real mfcbc = c0o1;
   real mfaba = c0o1;
   real mfcba = c0o1;
   real mfabc = c0o1;
   real mfbcc = c0o1;
   real mfbaa = c0o1;
   real mfbca = c0o1;
   real mfbac = c0o1;
   real mfbbb = c0o1;
   real mfccc = c0o1;
   real mfaac = c0o1;
   real mfcac = c0o1;
   real mfacc = c0o1;
   real mfcca = c0o1;
   real mfaaa = c0o1;
   real mfcaa = c0o1;
   real mfaca = c0o1;

   mfaaa = rho; // if drho is interpolated directly

   real vx1Sq = vx1*vx1;
   real vx2Sq = vx2*vx2;
   real vx3Sq = vx3*vx3;
   real oMdrho = c1o1;

   //2.f

   // linear combinations
   real mxxPyyPzz = mfaaa - c2o3*(ax + by + c2o1*axx*x + bxy*x + axy*y + c2o1*byy*y + axz*z + byz*z + bxyz*x*z + axyz*y*z + cz - cxz*x + cyz*y + cxyz*x*y + c2o1*czz*z)*eps_new / oP* (c1o1 + press);
   real mxxMyy    = -c2o3*(ax - by + kxxMyyAverage + c2o1*axx*x - bxy*x + axy*y - c2o1*byy*y + axz*z - byz*z - bxyz*x*z + axyz*y*z)*eps_new/o * (c1o1 + press);
   real mxxMzz    = -c2o3*(ax - cz + kxxMzzAverage + c2o1*axx*x - cxz*x + axy*y - cyz*y - cxyz*x*y + axz*z - c2o1*czz*z + axyz*y*z)*eps_new/o * (c1o1 + press);

   mfabb     = -c1o3 * (bz + cy + kyzAverage + bxz*x + cxy*x + byz*y + c2o1*cyy*y + bxyz*x*y + c2o1*bzz*z + cyz*z + cxyz*x*z)*eps_new/o * (c1o1 + press);
   mfbab     = -c1o3 * (az + cx + kxzAverage + axz*x + c2o1*cxx*x + ayz*y + cxy*y + axyz*x*y + c2o1*azz*z + cxz*z + cxyz*y*z)*eps_new/o * (c1o1 + press);
   mfbba     = -c1o3 * (ay + bx + kxyAverage + axy*x + c2o1*bxx*x + c2o1*ayy*y + bxy*y + ayz*z + bxz*z + axyz*x*z + bxyz*y*z)*eps_new/o * (c1o1 + press);

   // linear combinations back
   mfcaa = c1o3 * (mxxMyy +       mxxMzz + mxxPyyPzz) ;
   mfaca = c1o3 * (-c2o1 * mxxMyy +       mxxMzz + mxxPyyPzz) ;
   mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) ;

   //three
   mfbbb = c0o1;
   real mxxyPyzz = c0o1;
   real mxxyMyzz = c0o1;
   real mxxzPyyz = c0o1;
   real mxxzMyyz = c0o1;
   real mxyyPxzz =  c0o1;
   real mxyyMxzz = c0o1;

   // linear combinations back
   mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
   mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
   mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
   mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
   mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
   mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

   //4.f
   mfacc = mfaaa*c1o9;
   mfcac = mfacc;
   mfcca = mfacc;

   //5.

   //6.

   mfccc = mfaaa*c1o27;
   ////////////////////////////////////////////////////////////////////////////////////
   //back
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Z - Dir
   real m0 =  mfaac * c1o2 +      mfaab * (vx3 - c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq - vx3) * c1o2;
   real m1 = -mfaac        - c2o1 * mfaab *  vx3         +  mfaaa                * (c1o1 - vx3Sq)              - c1o1 * oMdrho * vx3Sq;
   real m2 =  mfaac * c1o2 +      mfaab * (vx3 + c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaaa = m0;
   mfaab = m1;
   mfaac = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfabc * c1o2 +      mfabb * (vx3 - c1o2) + mfaba * (vx3Sq - vx3) * c1o2;
   m1 = -mfabc        - c2o1 * mfabb *  vx3         + mfaba * (c1o1 - vx3Sq);
   m2 =  mfabc * c1o2 +      mfabb * (vx3 + c1o2) + mfaba * (vx3Sq + vx3) * c1o2;
   mfaba = m0;
   mfabb = m1;
   mfabc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfacb * (vx3 - c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfacc        - c2o1 * mfacb *  vx3         +  mfaca                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfacc * c1o2 +      mfacb * (vx3 + c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaca = m0;
   mfacb = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbac * c1o2 +      mfbab * (vx3 - c1o2) + mfbaa * (vx3Sq - vx3) * c1o2;
   m1 = -mfbac        - c2o1 * mfbab *  vx3         + mfbaa * (c1o1 - vx3Sq);
   m2 =  mfbac * c1o2 +      mfbab * (vx3 + c1o2) + mfbaa * (vx3Sq + vx3) * c1o2;
   mfbaa = m0;
   mfbab = m1;
   mfbac = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbbc * c1o2 +      mfbbb * (vx3 - c1o2) + mfbba * (vx3Sq - vx3) * c1o2;
   m1 = -mfbbc        - c2o1 * mfbbb *  vx3         + mfbba * (c1o1 - vx3Sq);
   m2 =  mfbbc * c1o2 +      mfbbb * (vx3 + c1o2) + mfbba * (vx3Sq + vx3) * c1o2;
   mfbba = m0;
   mfbbb = m1;
   mfbbc = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbcb * (vx3 - c1o2) + mfbca * (vx3Sq - vx3) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbcb *  vx3         + mfbca * (c1o1 - vx3Sq);
   m2 =  mfbcc * c1o2 +      mfbcb * (vx3 + c1o2) + mfbca * (vx3Sq + vx3) * c1o2;
   mfbca = m0;
   mfbcb = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfcab * (vx3 - c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfcac        - c2o1 * mfcab *  vx3         +  mfcaa                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfcac * c1o2 +      mfcab * (vx3 + c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcaa = m0;
   mfcab = m1;
   mfcac = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfcbb * (vx3 - c1o2) + mfcba * (vx3Sq - vx3) * c1o2;
   m1 = -mfcbc        - c2o1 * mfcbb *  vx3         + mfcba * (c1o1 - vx3Sq);
   m2 =  mfcbc * c1o2 +      mfcbb * (vx3 + c1o2) + mfcba * (vx3Sq + vx3) * c1o2;
   mfcba = m0;
   mfcbb = m1;
   mfcbc = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfccb * (vx3 - c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfccc        - c2o1 * mfccb *  vx3         +  mfcca                  * (c1o1 - vx3Sq)              - c1o9 * oMdrho * vx3Sq;
   m2 =  mfccc * c1o2 +      mfccb * (vx3 + c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcca = m0;
   mfccb = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Y - Dir
   m0 =  mfaca * c1o2 +      mfaba * (vx2 - c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfaca        - c2o1 * mfaba *  vx2         +  mfaaa                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfaca * c1o2 +      mfaba * (vx2 + c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaaa = m0;
   mfaba = m1;
   mfaca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacb * c1o2 +      mfabb * (vx2 - c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacb        - c2o1 * mfabb *  vx2         +  mfaab                  * (c1o1 - vx2Sq)              - c2o3 * oMdrho * vx2Sq;
   m2 =  mfacb * c1o2 +      mfabb * (vx2 + c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaab = m0;
   mfabb = m1;
   mfacb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfabc * (vx2 - c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacc        - c2o1 * mfabc *  vx2         +  mfaac                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfacc * c1o2 +      mfabc * (vx2 + c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaac = m0;
   mfabc = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbca * c1o2 +      mfbba * (vx2 - c1o2) + mfbaa * (vx2Sq - vx2) * c1o2;
   m1 = -mfbca        - c2o1 * mfbba *  vx2         + mfbaa * (c1o1 - vx2Sq);
   m2 =  mfbca * c1o2 +      mfbba * (vx2 + c1o2) + mfbaa * (vx2Sq + vx2) * c1o2;
   mfbaa = m0;
   mfbba = m1;
   mfbca = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcb * c1o2 +      mfbbb * (vx2 - c1o2) + mfbab * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcb        - c2o1 * mfbbb *  vx2         + mfbab * (c1o1 - vx2Sq);
   m2 =  mfbcb * c1o2 +      mfbbb * (vx2 + c1o2) + mfbab * (vx2Sq + vx2) * c1o2;
   mfbab = m0;
   mfbbb = m1;
   mfbcb = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbbc * (vx2 - c1o2) + mfbac * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbbc *  vx2         + mfbac * (c1o1 - vx2Sq);
   m2 =  mfbcc * c1o2 +      mfbbc * (vx2 + c1o2) + mfbac * (vx2Sq + vx2) * c1o2;
   mfbac = m0;
   mfbbc = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfcba * (vx2 - c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfcca        - c2o1 * mfcba *  vx2         +  mfcaa                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfcca * c1o2 +      mfcba * (vx2 + c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcaa = m0;
   mfcba = m1;
   mfcca = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfcbb * (vx2 - c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccb        - c2o1 * mfcbb *  vx2         +  mfcab                  * (c1o1 - vx2Sq)              - c2o9 * oMdrho * vx2Sq;
   m2 =  mfccb * c1o2 +      mfcbb * (vx2 + c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcab = m0;
   mfcbb = m1;
   mfccb = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfcbc * (vx2 - c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccc        - c2o1 * mfcbc *  vx2         +  mfcac                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfccc * c1o2 +      mfcbc * (vx2 + c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcac = m0;
   mfcbc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // X - Dir
   m0 =  mfcaa * c1o2 +      mfbaa * (vx1 - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcaa        - c2o1 * mfbaa *  vx1         +  mfaaa                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcaa * c1o2 +      mfbaa * (vx1 + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaaa = m0;
   mfbaa = m1;
   mfcaa = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcba * c1o2 +      mfbba * (vx1 - c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcba        - c2o1 * mfbba *  vx1         +  mfaba                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcba * c1o2 +      mfbba * (vx1 + c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaba = m0;
   mfbba = m1;
   mfcba = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfbca * (vx1 - c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcca        - c2o1 * mfbca *  vx1         +  mfaca                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcca * c1o2 +      mfbca * (vx1 + c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaca = m0;
   mfbca = m1;
   mfcca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcab * c1o2 +      mfbab * (vx1 - c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcab        - c2o1 * mfbab *  vx1         +  mfaab                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcab * c1o2 +      mfbab * (vx1 + c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaab = m0;
   mfbab = m1;
   mfcab = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfcbb * c1o2 +      mfbbb * (vx1 - c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbb        - c2o1 * mfbbb *  vx1         +  mfabb                  * (c1o1 - vx1Sq)              - c4o9 * oMdrho * vx1Sq;
   m2 =  mfcbb * c1o2 +      mfbbb * (vx1 + c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabb = m0;
   mfbbb = m1;
   mfcbb = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfbcb * (vx1 - c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccb        - c2o1 * mfbcb *  vx1         +  mfacb                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfccb * c1o2 +      mfbcb * (vx1 + c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacb = m0;
   mfbcb = m1;
   mfccb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfbac * (vx1 - c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcac        - c2o1 * mfbac *  vx1         +  mfaac                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcac * c1o2 +      mfbac * (vx1 + c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaac = m0;
   mfbac = m1;
   mfcac = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfbbc * (vx1 - c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbc        - c2o1 * mfbbc *  vx1         +  mfabc                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcbc * c1o2 +      mfbbc * (vx1 + c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabc = m0;
   mfbbc = m1;
   mfcbc = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfbcc * (vx1 - c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccc        - c2o1 * mfbcc *  vx1         +  mfacc                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfccc * c1o2 +      mfbcc * (vx1 + c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacc = m0;
   mfbcc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////

   f[DIR_P00]    = mfcbb;
   f[DIR_M00]    = mfabb;
   f[DIR_0P0]    = mfbcb;
   f[DIR_0M0]    = mfbab;
   f[DIR_00P]    = mfbbc;
   f[DIR_00M]    = mfbba;
   f[DIR_PP0]   = mfccb;
   f[DIR_MM0]   = mfaab;
   f[DIR_PM0]   = mfcab;
   f[DIR_MP0]   = mfacb;
   f[DIR_P0P]   = mfcbc;
   f[DIR_M0M]   = mfaba;
   f[DIR_P0M]   = mfcba;
   f[DIR_M0P]   = mfabc;
   f[DIR_0PP]   = mfbcc;
   f[DIR_0MM]   = mfbaa;
   f[DIR_0PM]   = mfbca;
   f[DIR_0MP]   = mfbac;
   f[DIR_000] = mfbbb;
   f[DIR_PPP]  = mfccc;
   f[DIR_PMP]  = mfcac;
   f[DIR_PPM]  = mfcca;
   f[DIR_PMM]  = mfcaa;
   f[DIR_MPP]  = mfacc;
   f[DIR_MMP]  = mfaac;
   f[DIR_MPM]  = mfaca;
   f[DIR_MMM]  = mfaaa;
}
//////////////////////////////////////////////////////////////////////////
//Position SWB -0.25, -0.25, -0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressBSW()
{
   return   press_SWT * (0.140625 + 0.1875 * xoff + 0.1875 * yoff - 0.5625 * zoff) +
      press_NWT * (0.046875 + 0.0625 * xoff - 0.1875 * yoff - 0.1875 * zoff) +
      press_SET * (0.046875 - 0.1875 * xoff + 0.0625 * yoff - 0.1875 * zoff) +
      press_NET * (0.015625 - 0.0625 * xoff - 0.0625 * yoff - 0.0625 * zoff) +
      press_NEB * (0.046875 - 0.1875 * xoff - 0.1875 * yoff + 0.0625 * zoff) +
      press_NWB * (0.140625 + 0.1875 * xoff - 0.5625 * yoff + 0.1875 * zoff) +
      press_SEB * (0.140625 - 0.5625 * xoff + 0.1875 * yoff + 0.1875 * zoff) +
      press_SWB * (0.421875 + 0.5625 * xoff + 0.5625 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SWT -0.25, -0.25, 0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressTSW()
{
   return   press_SWT * (0.421875 + 0.5625 * xoff + 0.5625 * yoff - 0.5625 * zoff) +
      press_NWT * (0.140625 + 0.1875 * xoff - 0.5625 * yoff - 0.1875 * zoff) +
      press_SET * (0.140625 - 0.5625 * xoff + 0.1875 * yoff - 0.1875 * zoff) +
      press_NET * (0.046875 - 0.1875 * xoff - 0.1875 * yoff - 0.0625 * zoff) +
      press_NEB * (0.015625 - 0.0625 * xoff - 0.0625 * yoff + 0.0625 * zoff) +
      press_NWB * (0.046875 + 0.0625 * xoff - 0.1875 * yoff + 0.1875 * zoff) +
      press_SEB * (0.046875 - 0.1875 * xoff + 0.0625 * yoff + 0.1875 * zoff) +
      press_SWB * (0.140625 + 0.1875 * xoff + 0.1875 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SET 0.25, -0.25, 0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressTSE()
{
   return   press_SET * (0.421875 - 0.5625 * xoff + 0.5625 * yoff - 0.5625 * zoff) +
      press_NET * (0.140625 - 0.1875 * xoff - 0.5625 * yoff - 0.1875 * zoff) +
      press_SWT * (0.140625 + 0.5625 * xoff + 0.1875 * yoff - 0.1875 * zoff) +
      press_NWT * (0.046875 + 0.1875 * xoff - 0.1875 * yoff - 0.0625 * zoff) +
      press_NWB * (0.015625 + 0.0625 * xoff - 0.0625 * yoff + 0.0625 * zoff) +
      press_NEB * (0.046875 - 0.0625 * xoff - 0.1875 * yoff + 0.1875 * zoff) +
      press_SWB * (0.046875 + 0.1875 * xoff + 0.0625 * yoff + 0.1875 * zoff) +
      press_SEB * (0.140625 - 0.1875 * xoff + 0.1875 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position SEB 0.25, -0.25, -0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressBSE()
{
   return   press_SET * (0.140625 - 0.1875 * xoff + 0.1875 * yoff - 0.5625 * zoff) +
      press_NET * (0.046875 - 0.0625 * xoff - 0.1875 * yoff - 0.1875 * zoff) +
      press_SWT * (0.046875 + 0.1875 * xoff + 0.0625 * yoff - 0.1875 * zoff) +
      press_NWT * (0.015625 + 0.0625 * xoff - 0.0625 * yoff - 0.0625 * zoff) +
      press_NWB * (0.046875 + 0.1875 * xoff - 0.1875 * yoff + 0.0625 * zoff) +
      press_NEB * (0.140625 - 0.1875 * xoff - 0.5625 * yoff + 0.1875 * zoff) +
      press_SWB * (0.140625 + 0.5625 * xoff + 0.1875 * yoff + 0.1875 * zoff) +
      press_SEB * (0.421875 - 0.5625 * xoff + 0.5625 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWB -0.25, 0.25, -0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressBNW()
{
   return   press_NWT * (0.140625 + 0.1875 * xoff - 0.1875 * yoff - 0.5625 * zoff) +
      press_NET * (0.046875 - 0.1875 * xoff - 0.0625 * yoff - 0.1875 * zoff) +
      press_SWT * (0.046875 + 0.0625 * xoff + 0.1875 * yoff - 0.1875 * zoff) +
      press_SET * (0.015625 - 0.0625 * xoff + 0.0625 * yoff - 0.0625 * zoff) +
      press_SEB * (0.046875 - 0.1875 * xoff + 0.1875 * yoff + 0.0625 * zoff) +
      press_NEB * (0.140625 - 0.5625 * xoff - 0.1875 * yoff + 0.1875 * zoff) +
      press_SWB * (0.140625 + 0.1875 * xoff + 0.5625 * yoff + 0.1875 * zoff) +
      press_NWB * (0.421875 + 0.5625 * xoff - 0.5625 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NWT -0.25, 0.25, 0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressTNW()
{
   return   press_NWT * (0.421875 + 0.5625 * xoff - 0.5625 * yoff - 0.5625 * zoff) +
      press_NET * (0.140625 - 0.5625 * xoff - 0.1875 * yoff - 0.1875 * zoff) +
      press_SWT * (0.140625 + 0.1875 * xoff + 0.5625 * yoff - 0.1875 * zoff) +
      press_SET * (0.046875 - 0.1875 * xoff + 0.1875 * yoff - 0.0625 * zoff) +
      press_SEB * (0.015625 - 0.0625 * xoff + 0.0625 * yoff + 0.0625 * zoff) +
      press_NEB * (0.046875 - 0.1875 * xoff - 0.0625 * yoff + 0.1875 * zoff) +
      press_SWB * (0.046875 + 0.0625 * xoff + 0.1875 * yoff + 0.1875 * zoff) +
      press_NWB * (0.140625 + 0.1875 * xoff - 0.1875 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NET 0.25, 0.25, 0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressTNE()
{
   return   press_NET * (0.421875 - 0.5625 * xoff - 0.5625 * yoff - 0.5625 * zoff) +
      press_NWT * (0.140625 + 0.5625 * xoff - 0.1875 * yoff - 0.1875 * zoff) +
      press_SET * (0.140625 - 0.1875 * xoff + 0.5625 * yoff - 0.1875 * zoff) +
      press_SWT * (0.046875 + 0.1875 * xoff + 0.1875 * yoff - 0.0625 * zoff) +
      press_SWB * (0.015625 + 0.0625 * xoff + 0.0625 * yoff + 0.0625 * zoff) +
      press_NWB * (0.046875 + 0.1875 * xoff - 0.0625 * yoff + 0.1875 * zoff) +
      press_SEB * (0.046875 - 0.0625 * xoff + 0.1875 * yoff + 0.1875 * zoff) +
      press_NEB * (0.140625 - 0.1875 * xoff - 0.1875 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position NEB 0.25, 0.25, -0.25
real CompressibleOffsetSquarePressureInterpolationProcessor::calcPressBNE()
{
   return   press_NET * (0.140625 - 0.1875 * xoff - 0.1875 * yoff - 0.5625 * zoff) +
      press_NWT * (0.046875 + 0.1875 * xoff - 0.0625 * yoff - 0.1875 * zoff) +
      press_SET * (0.046875 - 0.0625 * xoff + 0.1875 * yoff - 0.1875 * zoff) +
      press_SWT * (0.015625 + 0.0625 * xoff + 0.0625 * yoff - 0.0625 * zoff) +
      press_SWB * (0.046875 + 0.1875 * xoff + 0.1875 * yoff + 0.0625 * zoff) +
      press_NWB * (0.140625 + 0.5625 * xoff - 0.1875 * yoff + 0.1875 * zoff) +
      press_SEB * (0.140625 - 0.1875 * xoff + 0.5625 * yoff + 0.1875 * zoff) +
      press_NEB * (0.421875 - 0.5625 * xoff - 0.5625 * yoff + 0.5625 * zoff);
}
//////////////////////////////////////////////////////////////////////////
//Position C 0.0, 0.0, 0.0
void CompressibleOffsetSquarePressureInterpolationProcessor::calcInterpolatedNodeFC(real* f, real omega)
{
   using namespace D3Q27System;
   using namespace vf::lbm::dir;

   real press  =  press_NET * (0.125 - 0.25 * xoff - 0.25 * yoff - 0.25 * zoff) +
      press_NWT * (0.125 + 0.25 * xoff - 0.25 * yoff - 0.25 * zoff) +
      press_SET * (0.125 - 0.25 * xoff + 0.25 * yoff - 0.25 * zoff) +
      press_SWT * (0.125 + 0.25 * xoff + 0.25 * yoff - 0.25 * zoff) +
      press_NEB * (0.125 - 0.25 * xoff - 0.25 * yoff + 0.25 * zoff) +
      press_NWB * (0.125 + 0.25 * xoff - 0.25 * yoff + 0.25 * zoff) +
      press_SEB * (0.125 - 0.25 * xoff + 0.25 * yoff + 0.25 * zoff) +
      press_SWB * (0.125 + 0.25 * xoff + 0.25 * yoff + 0.25 * zoff);
   real vx1  = a0;
   real vx2  = b0;
   real vx3  = c0;
  
   
   real rho = press ;//+ (ax+by+cz)/3.;

   real laplaceRho = (xoff!=0.0 || yoff!=0.0 || zoff!= 0.0) ? 0.0 :(-3.0*(by*by+ax*ax+cz*cz)-6.0*(ay*bx+bz*cy+az*cx))*(1.0+rho);

   rho=rho-laplaceRho*0.25;

   real eps_new = 2.0;
   real o  = omega;
   //bulk viscosity
   real oP = OxxPyyPzzC;

   real mfcbb = c0o1;
   real mfabb = c0o1;
   real mfbcb = c0o1;
   real mfbab = c0o1;
   real mfbbc = c0o1;
   real mfbba = c0o1;
   real mfccb = c0o1;
   real mfaab = c0o1;
   real mfcab = c0o1;
   real mfacb = c0o1;
   real mfcbc = c0o1;
   real mfaba = c0o1;
   real mfcba = c0o1;
   real mfabc = c0o1;
   real mfbcc = c0o1;
   real mfbaa = c0o1;
   real mfbca = c0o1;
   real mfbac = c0o1;
   real mfbbb = c0o1;
   real mfccc = c0o1;
   real mfaac = c0o1;
   real mfcac = c0o1;
   real mfacc = c0o1;
   real mfcca = c0o1;
   real mfaaa = c0o1;
   real mfcaa = c0o1;
   real mfaca = c0o1;

   mfaaa = rho; // if drho is interpolated directly

   real vx1Sq = vx1*vx1;
   real vx2Sq = vx2*vx2;
   real vx3Sq = vx3*vx3;
   real oMdrho = c1o1;
   //oMdrho = c1o1 - mfaaa;

   //2.f
   // linear combinations

/////////////////////////
   real mxxPyyPzz = mfaaa    -c2o3*(ax+by+cz)*eps_new/oP*(c1o1+press);

   real mxxMyy    = -c2o3*((ax - by)+kxxMyyAverage)*eps_new/o * (c1o1 + press);
   real mxxMzz    = -c2o3*((ax - cz)+kxxMzzAverage)*eps_new/o * (c1o1 + press);

   mfabb     = -c1o3 * ((bz + cy)+kyzAverage)*eps_new/o * (c1o1 + press);
   mfbab     = -c1o3 * ((az + cx)+kxzAverage)*eps_new/o * (c1o1 + press);
   mfbba     = -c1o3 * ((ay + bx)+kxyAverage)*eps_new/o * (c1o1 + press);

   ////////////////////////
   // linear combinations back
   mfcaa = c1o3 * (mxxMyy +       mxxMzz + mxxPyyPzz);
   mfaca = c1o3 * (-c2o1 * mxxMyy +       mxxMzz + mxxPyyPzz);
   mfaac = c1o3 * (mxxMyy - c2o1 * mxxMzz + mxxPyyPzz);

   //three
   mfbbb = c0o1;

   real mxxyPyzz = c0o1;
   real mxxyMyzz = c0o1;
   real mxxzPyyz = c0o1;
   real mxxzMyyz = c0o1;
   real mxyyPxzz =  c0o1;
   real mxyyMxzz = c0o1;

   // linear combinations back
   mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
   mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
   mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
   mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
   mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
   mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

   //4.f
   mfacc = mfaaa*c1o9;
   mfcac = mfacc;
   mfcca = mfacc;
   //5.

   //6.
   mfccc = mfaaa*c1o27;
   ////////////////////////////////////////////////////////////////////////////////////
   //back
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Z - Dir
   real m0 =  mfaac * c1o2 +      mfaab * (vx3 - c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq - vx3) * c1o2;
   real m1 = -mfaac        - c2o1 * mfaab *  vx3         +  mfaaa                * (c1o1 - vx3Sq)              - c1o1 * oMdrho * vx3Sq;
   real m2 =  mfaac * c1o2 +      mfaab * (vx3 + c1o2) + (mfaaa + c1o1 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaaa = m0;
   mfaab = m1;
   mfaac = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfabc * c1o2 +      mfabb * (vx3 - c1o2) + mfaba * (vx3Sq - vx3) * c1o2;
   m1 = -mfabc        - c2o1 * mfabb *  vx3         + mfaba * (c1o1 - vx3Sq);
   m2 =  mfabc * c1o2 +      mfabb * (vx3 + c1o2) + mfaba * (vx3Sq + vx3) * c1o2;
   mfaba = m0;
   mfabb = m1;
   mfabc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfacb * (vx3 - c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfacc        - c2o1 * mfacb *  vx3         +  mfaca                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfacc * c1o2 +      mfacb * (vx3 + c1o2) + (mfaca + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfaca = m0;
   mfacb = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbac * c1o2 +      mfbab * (vx3 - c1o2) + mfbaa * (vx3Sq - vx3) * c1o2;
   m1 = -mfbac        - c2o1 * mfbab *  vx3         + mfbaa * (c1o1 - vx3Sq);
   m2 =  mfbac * c1o2 +      mfbab * (vx3 + c1o2) + mfbaa * (vx3Sq + vx3) * c1o2;
   mfbaa = m0;
   mfbab = m1;
   mfbac = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbbc * c1o2 +      mfbbb * (vx3 - c1o2) + mfbba * (vx3Sq - vx3) * c1o2;
   m1 = -mfbbc        - c2o1 * mfbbb *  vx3         + mfbba * (c1o1 - vx3Sq);
   m2 =  mfbbc * c1o2 +      mfbbb * (vx3 + c1o2) + mfbba * (vx3Sq + vx3) * c1o2;
   mfbba = m0;
   mfbbb = m1;
   mfbbc = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbcb * (vx3 - c1o2) + mfbca * (vx3Sq - vx3) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbcb *  vx3         + mfbca * (c1o1 - vx3Sq);
   m2 =  mfbcc * c1o2 +      mfbcb * (vx3 + c1o2) + mfbca * (vx3Sq + vx3) * c1o2;
   mfbca = m0;
   mfbcb = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfcab * (vx3 - c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfcac        - c2o1 * mfcab *  vx3         +  mfcaa                  * (c1o1 - vx3Sq)              - c1o3 * oMdrho * vx3Sq;
   m2 =  mfcac * c1o2 +      mfcab * (vx3 + c1o2) + (mfcaa + c1o3 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcaa = m0;
   mfcab = m1;
   mfcac = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfcbb * (vx3 - c1o2) + mfcba * (vx3Sq - vx3) * c1o2;
   m1 = -mfcbc        - c2o1 * mfcbb *  vx3         + mfcba * (c1o1 - vx3Sq);
   m2 =  mfcbc * c1o2 +      mfcbb * (vx3 + c1o2) + mfcba * (vx3Sq + vx3) * c1o2;
   mfcba = m0;
   mfcbb = m1;
   mfcbc = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfccb * (vx3 - c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq - vx3) * c1o2;
   m1 = -mfccc        - c2o1 * mfccb *  vx3         +  mfcca                  * (c1o1 - vx3Sq)              - c1o9 * oMdrho * vx3Sq;
   m2 =  mfccc * c1o2 +      mfccb * (vx3 + c1o2) + (mfcca + c1o9 * oMdrho) * (vx3Sq + vx3) * c1o2;
   mfcca = m0;
   mfccb = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // Y - Dir
   m0 =  mfaca * c1o2 +      mfaba * (vx2 - c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfaca        - c2o1 * mfaba *  vx2         +  mfaaa                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfaca * c1o2 +      mfaba * (vx2 + c1o2) + (mfaaa + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaaa = m0;
   mfaba = m1;
   mfaca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacb * c1o2 +      mfabb * (vx2 - c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacb        - c2o1 * mfabb *  vx2         +  mfaab                  * (c1o1 - vx2Sq)              - c2o3 * oMdrho * vx2Sq;
   m2 =  mfacb * c1o2 +      mfabb * (vx2 + c1o2) + (mfaab + c2o3 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaab = m0;
   mfabb = m1;
   mfacb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfacc * c1o2 +      mfabc * (vx2 - c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfacc        - c2o1 * mfabc *  vx2         +  mfaac                  * (c1o1 - vx2Sq)              - c1o6 * oMdrho * vx2Sq;
   m2 =  mfacc * c1o2 +      mfabc * (vx2 + c1o2) + (mfaac + c1o6 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfaac = m0;
   mfabc = m1;
   mfacc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfbca * c1o2 +      mfbba * (vx2 - c1o2) + mfbaa * (vx2Sq - vx2) * c1o2;
   m1 = -mfbca        - c2o1 * mfbba *  vx2         + mfbaa * (c1o1 - vx2Sq);
   m2 =  mfbca * c1o2 +      mfbba * (vx2 + c1o2) + mfbaa * (vx2Sq + vx2) * c1o2;
   mfbaa = m0;
   mfbba = m1;
   mfbca = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcb * c1o2 +      mfbbb * (vx2 - c1o2) + mfbab * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcb        - c2o1 * mfbbb *  vx2         + mfbab * (c1o1 - vx2Sq);
   m2 =  mfbcb * c1o2 +      mfbbb * (vx2 + c1o2) + mfbab * (vx2Sq + vx2) * c1o2;
   mfbab = m0;
   mfbbb = m1;
   mfbcb = m2;
   /////////b//////////////////////////////////////////////////////////////////////////
   m0 =  mfbcc * c1o2 +      mfbbc * (vx2 - c1o2) + mfbac * (vx2Sq - vx2) * c1o2;
   m1 = -mfbcc        - c2o1 * mfbbc *  vx2         + mfbac * (c1o1 - vx2Sq);
   m2 =  mfbcc * c1o2 +      mfbbc * (vx2 + c1o2) + mfbac * (vx2Sq + vx2) * c1o2;
   mfbac = m0;
   mfbbc = m1;
   mfbcc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfcba * (vx2 - c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfcca        - c2o1 * mfcba *  vx2         +  mfcaa                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfcca * c1o2 +      mfcba * (vx2 + c1o2) + (mfcaa + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcaa = m0;
   mfcba = m1;
   mfcca = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfcbb * (vx2 - c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccb        - c2o1 * mfcbb *  vx2         +  mfcab                  * (c1o1 - vx2Sq)              - c2o9 * oMdrho * vx2Sq;
   m2 =  mfccb * c1o2 +      mfcbb * (vx2 + c1o2) + (mfcab + c2o9 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcab = m0;
   mfcbb = m1;
   mfccb = m2;
   /////////c//////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfcbc * (vx2 - c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq - vx2) * c1o2;
   m1 = -mfccc        - c2o1 * mfcbc *  vx2         +  mfcac                   * (c1o1 - vx2Sq)              - c1o18 * oMdrho * vx2Sq;
   m2 =  mfccc * c1o2 +      mfcbc * (vx2 + c1o2) + (mfcac + c1o18 * oMdrho) * (vx2Sq + vx2) * c1o2;
   mfcac = m0;
   mfcbc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
   ////////////////////////////////////////////////////////////////////////////////////
   // X - Dir
   m0 =  mfcaa * c1o2 +      mfbaa * (vx1 - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcaa        - c2o1 * mfbaa *  vx1         +  mfaaa                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcaa * c1o2 +      mfbaa * (vx1 + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaaa = m0;
   mfbaa = m1;
   mfcaa = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcba * c1o2 +      mfbba * (vx1 - c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcba        - c2o1 * mfbba *  vx1         +  mfaba                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcba * c1o2 +      mfbba * (vx1 + c1o2) + (mfaba + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaba = m0;
   mfbba = m1;
   mfcba = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcca * c1o2 +      mfbca * (vx1 - c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcca        - c2o1 * mfbca *  vx1         +  mfaca                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcca * c1o2 +      mfbca * (vx1 + c1o2) + (mfaca + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaca = m0;
   mfbca = m1;
   mfcca = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcab * c1o2 +      mfbab * (vx1 - c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcab        - c2o1 * mfbab *  vx1         +  mfaab                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcab * c1o2 +      mfbab * (vx1 + c1o2) + (mfaab + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaab = m0;
   mfbab = m1;
   mfcab = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfcbb * c1o2 +      mfbbb * (vx1 - c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbb        - c2o1 * mfbbb *  vx1         +  mfabb                  * (c1o1 - vx1Sq)              - c4o9 * oMdrho * vx1Sq;
   m2 =  mfcbb * c1o2 +      mfbbb * (vx1 + c1o2) + (mfabb + c4o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabb = m0;
   mfbbb = m1;
   mfcbb = m2;
   ///////////b////////////////////////////////////////////////////////////////////////
   m0 =  mfccb * c1o2 +      mfbcb * (vx1 - c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccb        - c2o1 * mfbcb *  vx1         +  mfacb                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfccb * c1o2 +      mfbcb * (vx1 + c1o2) + (mfacb + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacb = m0;
   mfbcb = m1;
   mfccb = m2;
   ////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////
   m0 =  mfcac * c1o2 +      mfbac * (vx1 - c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcac        - c2o1 * mfbac *  vx1         +  mfaac                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfcac * c1o2 +      mfbac * (vx1 + c1o2) + (mfaac + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfaac = m0;
   mfbac = m1;
   mfcac = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfcbc * c1o2 +      mfbbc * (vx1 - c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfcbc        - c2o1 * mfbbc *  vx1         +  mfabc                  * (c1o1 - vx1Sq)              - c1o9 * oMdrho * vx1Sq;
   m2 =  mfcbc * c1o2 +      mfbbc * (vx1 + c1o2) + (mfabc + c1o9 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfabc = m0;
   mfbbc = m1;
   mfcbc = m2;
   ///////////c////////////////////////////////////////////////////////////////////////
   m0 =  mfccc * c1o2 +      mfbcc * (vx1 - c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq - vx1) * c1o2;
   m1 = -mfccc        - c2o1 * mfbcc *  vx1         +  mfacc                   * (c1o1 - vx1Sq)              - c1o36 * oMdrho * vx1Sq;
   m2 =  mfccc * c1o2 +      mfbcc * (vx1 + c1o2) + (mfacc + c1o36 * oMdrho) * (vx1Sq + vx1) * c1o2;
   mfacc = m0;
   mfbcc = m1;
   mfccc = m2;
   ////////////////////////////////////////////////////////////////////////////////////

   f[DIR_P00]    = mfcbb;
   f[DIR_M00]    = mfabb;
   f[DIR_0P0]    = mfbcb;
   f[DIR_0M0]    = mfbab;
   f[DIR_00P]    = mfbbc;
   f[DIR_00M]    = mfbba;
   f[DIR_PP0]   = mfccb;
   f[DIR_MM0]   = mfaab;
   f[DIR_PM0]   = mfcab;
   f[DIR_MP0]   = mfacb;
   f[DIR_P0P]   = mfcbc;
   f[DIR_M0M]   = mfaba;
   f[DIR_P0M]   = mfcba;
   f[DIR_M0P]   = mfabc;
   f[DIR_0PP]   = mfbcc;
   f[DIR_0MM]   = mfbaa;
   f[DIR_0PM]   = mfbca;
   f[DIR_0MP]   = mfbac;
   f[DIR_000] = mfbbb;
   f[DIR_PPP]  = mfccc;
   f[DIR_PMP]  = mfcac;
   f[DIR_PPM]  = mfcca;
   f[DIR_PMM]  = mfcaa;
   f[DIR_MPP]  = mfacc;
   f[DIR_MMP]  = mfaac;
   f[DIR_MPM]  = mfaca;
   f[DIR_MMM]  = mfaaa;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::calcInterpolatedVelocity(real x, real y, real z, real& vx1, real& vx2, real& vx3)
{
	vx1  = a0 + ax*x + ay*y + az*z + axx*x*x + ayy*y*y + azz*z*z + axy*x*y + axz*x*z + ayz*y*z+axyz*x*y*z;
	vx2  = b0 + bx*x + by*y + bz*z + bxx*x*x + byy*y*y + bzz*z*z + bxy*x*y + bxz*x*z + byz*y*z+bxyz*x*y*z;
	vx3  = c0 + cx*x + cy*y + cz*z + cxx*x*x + cyy*y*y + czz*z*z + cxy*x*y + cxz*x*z + cyz*y*z+cxyz*x*y*z;
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::calcInterpolatedShearStress(real x, real y, real z,real& tauxx, real& tauyy, real& tauzz,real& tauxy, real& tauxz, real& tauyz)
{
	tauxx=ax+2*axx*x+axy*y+axz*z+axyz*y*z;
	tauyy=by+2*byy*y+bxy*x+byz*z+bxyz*x*z;
	tauzz=cz+2*czz*z+cxz*x+cyz*y+cxyz*x*y;
	tauxy=0.5*((ay+2.0*ayy*y+axy*x+ayz*z+axyz*x*z)+(bx+2.0*bxx*x+bxy*y+bxz*z+bxyz*y*z));
	tauxz=0.5*((az+2.0*azz*z+axz*x+ayz*y+axyz*x*y)+(cx+2.0*cxx*x+cxy*y+cxz*z+cxyz*y*z));
	tauyz=0.5*((bz+2.0*bzz*z+bxz*x+byz*y+bxyz*x*y)+(cy+2.0*cyy*y+cxy*x+cyz*z+cxyz*x*z));
}
//////////////////////////////////////////////////////////////////////////
void CompressibleOffsetSquarePressureInterpolationProcessor::setBulkOmegaToOmega(bool value)
{
   bulkOmegaToOmega = value;
}

