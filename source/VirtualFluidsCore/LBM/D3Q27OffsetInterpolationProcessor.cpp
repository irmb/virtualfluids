#include "D3Q27OffsetInterpolationProcessor.h"

#include <boost/foreach.hpp>

D3Q27OffsetInterpolationProcessor::D3Q27OffsetInterpolationProcessor()
   : omegaC(0.0), omegaF(0.0)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
D3Q27OffsetInterpolationProcessor::D3Q27OffsetInterpolationProcessor(LBMReal omegaC, LBMReal omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
D3Q27OffsetInterpolationProcessor::~D3Q27OffsetInterpolationProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
D3Q27InterpolationProcessorPtr D3Q27OffsetInterpolationProcessor::clone()
{
   return D3Q27InterpolationProcessorPtr (new D3Q27OffsetInterpolationProcessor(this->omegaC, this->omegaF));
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::init()
{
   calcFeqsForDirFct = NULL;
   calcMacrosFct     = NULL;
   calcFeqFct        = NULL;
   //SimulationParametersPtr param = SimulationParameters::getInstanz();

   //if(param->isCompressibleModel())
   //{
      calcMacrosFct     = &D3Q27System::calcCompMacroscopicValues;
      calcFeqFct        = &D3Q27System::calcCompFeq; 
   //}
   //else
   //{
   //   calcMacrosFct     = &D3Q27System::calcIncompMacroscopicValues;
   //   calcFeqFct        = &D3Q27System::calcIncompFeq; 
   //}
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::setOmegas( LBMReal omegaC, LBMReal omegaF )
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::setOffsets(LBMReal xoff, LBMReal yoff, LBMReal zoff)
{
   this->xoff = xoff;
   this->yoff = yoff;
   this->zoff = zoff;     
   this->xoff_sq = xoff * xoff;
   this->yoff_sq = yoff * yoff;
   this->zoff_sq = zoff * zoff;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellC, omegaC, 0.5);
   calcInterpolatedNode(icellF.BSW, omegaF, -0.25, -0.25, -0.25, calcPressBSW(), -1, -1, -1);
   calcInterpolatedNode(icellF.BNE, omegaF,  0.25,  0.25, -0.25, calcPressBNE(),  1,  1, -1);
   calcInterpolatedNode(icellF.TNW, omegaF, -0.25,  0.25,  0.25, calcPressTNW(), -1,  1,  1);
   calcInterpolatedNode(icellF.TSE, omegaF,  0.25, -0.25,  0.25, calcPressTSE(),  1, -1,  1);
   calcInterpolatedNode(icellF.BNW, omegaF, -0.25,  0.25, -0.25, calcPressBNW(), -1,  1, -1);
   calcInterpolatedNode(icellF.BSE, omegaF,  0.25, -0.25, -0.25, calcPressBSE(),  1, -1, -1);
   calcInterpolatedNode(icellF.TSW, omegaF, -0.25, -0.25,  0.25, calcPressTSW(), -1, -1,  1);
   calcInterpolatedNode(icellF.TNE, omegaF,  0.25,  0.25,  0.25, calcPressTNE(),  1,  1,  1);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff)
{
   setOffsets(xoff, yoff, zoff);
   calcInterpolatedCoefficiets(icellF, omegaF, 2.0);
   calcInterpolatedNodeFC(icellC, omegaC);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::calcMoments(const LBMReal* const f, LBMReal omega, LBMReal& press, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3, 
   LBMReal& kxy, LBMReal& kyz, LBMReal& kxz, LBMReal& kxxMyy, LBMReal& kxxMzz)
{
   using namespace D3Q27System;
   
   LBMReal rho = 0.0;
   calcMacrosFct(f,rho,vx1,vx2,vx3);
   //press = D3Q27System::calcPress(f,rho,vx1,vx2,vx3);
   press = rho; //interpolate rho!

   kxy    = -3.*omega*(f[SW]+f[TSW]+f[BSW]-f[NW]-f[TNW]-f[BNW]-f[SE]-f[TSE]-f[BSE]+f[NE]+f[TNE]+f[BNE]-(vx1*vx2));
   kyz    = -3.*omega*(f[BS]+f[BSW]+f[BSE]-f[TS]-f[TSW]-f[TSE]-f[BN]-f[BNW]-f[BNE]+f[TN]+f[TNW]+f[TNE]-(vx2*vx3));
   kxz    = -3.*omega*(f[BW]+f[BNW]+f[BSW]-f[TW]-f[TNW]-f[TSW]-f[BE]-f[BNE]-f[BSE]+f[TE]+f[TNE]+f[TSE]-(vx1*vx3));
   kxxMyy = -3./2.*omega*(f[BW]+f[W]+f[TW]-f[BS]-f[S]-f[TS]-f[BN]-f[N]-f[TN]+f[BE]+f[E]+f[TE]-(vx1*vx1-vx2*vx2));
   kxxMzz = -3./2.*omega*(f[SW]+f[W]+f[NW]-f[BS]-f[TS]-f[B]-f[T]-f[BN]-f[TN]+f[SE]+f[E]+f[NE]-(vx1*vx1-vx3*vx3));
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::calcInterpolatedCoefficiets(const D3Q27ICell& icell, LBMReal omega, LBMReal eps_new)
{
   LBMReal        vx1_SWT,vx2_SWT,vx3_SWT;
   LBMReal        vx1_NWT,vx2_NWT,vx3_NWT;
   LBMReal        vx1_NET,vx2_NET,vx3_NET;
   LBMReal        vx1_SET,vx2_SET,vx3_SET;
   LBMReal        vx1_SWB,vx2_SWB,vx3_SWB;
   LBMReal        vx1_NWB,vx2_NWB,vx3_NWB;
   LBMReal        vx1_NEB,vx2_NEB,vx3_NEB;
   LBMReal        vx1_SEB,vx2_SEB,vx3_SEB;

   LBMReal        kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT;
   LBMReal        kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT;
   LBMReal        kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET;
   LBMReal        kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET;
   LBMReal        kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB;
   LBMReal        kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB;
   LBMReal        kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB;
   LBMReal        kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB;

   calcMoments(icell.TSW,omega,press_SWT,vx1_SWT,vx2_SWT,vx3_SWT, kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT);
   calcMoments(icell.TNW,omega,press_NWT,vx1_NWT,vx2_NWT,vx3_NWT, kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT);
   calcMoments(icell.TNE,omega,press_NET,vx1_NET,vx2_NET,vx3_NET, kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET);
   calcMoments(icell.TSE,omega,press_SET,vx1_SET,vx2_SET,vx3_SET, kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET);
   calcMoments(icell.BSW,omega,press_SWB,vx1_SWB,vx2_SWB,vx3_SWB, kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB);
   calcMoments(icell.BNW,omega,press_NWB,vx1_NWB,vx2_NWB,vx3_NWB, kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB);
   calcMoments(icell.BNE,omega,press_NEB,vx1_NEB,vx2_NEB,vx3_NEB, kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB);
   calcMoments(icell.BSE,omega,press_SEB,vx1_SEB,vx2_SEB,vx3_SEB, kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB);

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

   const LBMReal o = omega;

   f_E = eps_new*((2*(-2*ax + by + cz))/(27.*o));
   f_N = eps_new*((2*(ax - 2*by + cz))/(27.*o));
   f_T = eps_new*((2*(ax + by - 2*cz))/(27.*o));
   f_NE = eps_new*(-(ax + 3*ay + 3*bx + by - 2*cz)/(54.*o));
   f_SE = eps_new*(-(ax - 3*ay - 3*bx + by - 2*cz)/(54.*o));
   f_TE = eps_new*(-(ax + 3*az - 2*by + 3*cx + cz)/(54.*o));
   f_BE = eps_new*(-(ax - 3*az - 2*by - 3*cx + cz)/(54.*o));
   f_TN = eps_new*(-(-2*ax + by + 3*bz + 3*cy + cz)/(54.*o));
   f_BN = eps_new*(-(-2*ax + by - 3*bz - 3*cy + cz)/(54.*o));
   f_ZERO = 0.;
   f_TNE = eps_new*(-(ay + az + bx + bz + cx + cy)/(72.*o));
   f_TSW = eps_new*((-ay + az - bx + bz + cx + cy)/(72.*o));
   f_TSE = eps_new*((ay - az + bx + bz - cx + cy)/(72.*o));
   f_TNW = eps_new*((ay + az + bx - bz + cx - cy)/(72.*o));

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
void D3Q27OffsetInterpolationProcessor::calcInterpolatedNode(LBMReal* f, LBMReal omega, LBMReal x, LBMReal y, LBMReal z, LBMReal press, LBMReal xs, LBMReal ys, LBMReal zs)
{
   using namespace D3Q27System;

   LBMReal rho  = press ;//+ (2.*axx*x+axy*y+axz*z+axyz*y*z+ax + 2.*byy*y+bxy*x+byz*z+bxyz*x*z+by + 2.*czz*z+cxz*x+cyz*y+cxyz*x*y+cz)/3.;
   LBMReal vx1  = a0 + 0.25*( xs*ax + ys*ay + zs*az) + 0.0625*(axx + xs*ys*axy + xs*zs*axz + ayy + ys*zs*ayz + azz) + 0.015625*(xs*ys*zs*axyz);
   LBMReal vx2  = b0 + 0.25*( xs*bx + ys*by + zs*bz) + 0.0625*(bxx + xs*ys*bxy + xs*zs*bxz + byy + ys*zs*byz + bzz) + 0.015625*(xs*ys*zs*bxyz);
   LBMReal vx3  = c0 + 0.25*( xs*cx + ys*cy + zs*cz) + 0.0625*(cxx + xs*ys*cxy + xs*zs*cxz + cyy + ys*zs*cyz + czz) + 0.015625*(xs*ys*zs*cxyz);

   LBMReal feq[ENDF+1];
   calcFeqFct(feq,rho,vx1,vx2,vx3);

   f[E]    = f_E    + xs*x_E    + ys*y_E    + zs*z_E    + xs*ys*xy_E    + xs*zs*xz_E    + ys*zs*yz_E    + feq[E];
   f[W]    = f_E    + xs*x_E    + ys*y_E    + zs*z_E    + xs*ys*xy_E    + xs*zs*xz_E    + ys*zs*yz_E    + feq[W];
   f[N]    = f_N    + xs*x_N    + ys*y_N    + zs*z_N    + xs*ys*xy_N    + xs*zs*xz_N    + ys*zs*yz_N    + feq[N];
   f[S]    = f_N    + xs*x_N    + ys*y_N    + zs*z_N    + xs*ys*xy_N    + xs*zs*xz_N    + ys*zs*yz_N    + feq[S];
   f[T]    = f_T    + xs*x_T    + ys*y_T    + zs*z_T    + xs*ys*xy_T    + xs*zs*xz_T    + ys*zs*yz_T    + feq[T];
   f[B]    = f_T    + xs*x_T    + ys*y_T    + zs*z_T    + xs*ys*xy_T    + xs*zs*xz_T    + ys*zs*yz_T    + feq[B];
   f[NE]   = f_NE   + xs*x_NE   + ys*y_NE   + zs*z_NE   + xs*ys*xy_NE   + xs*zs*xz_NE   + ys*zs*yz_NE   + feq[NE];
   f[SW]   = f_NE   + xs*x_NE   + ys*y_NE   + zs*z_NE   + xs*ys*xy_NE   + xs*zs*xz_NE   + ys*zs*yz_NE   + feq[SW];
   f[SE]   = f_SE   + xs*x_SE   + ys*y_SE   + zs*z_SE   + xs*ys*xy_SE   + xs*zs*xz_SE   + ys*zs*yz_SE   + feq[SE];
   f[NW]   = f_SE   + xs*x_SE   + ys*y_SE   + zs*z_SE   + xs*ys*xy_SE   + xs*zs*xz_SE   + ys*zs*yz_SE   + feq[NW];
   f[TE]   = f_TE   + xs*x_TE   + ys*y_TE   + zs*z_TE   + xs*ys*xy_TE   + xs*zs*xz_TE   + ys*zs*yz_TE   + feq[TE];
   f[BW]   = f_TE   + xs*x_TE   + ys*y_TE   + zs*z_TE   + xs*ys*xy_TE   + xs*zs*xz_TE   + ys*zs*yz_TE   + feq[BW];
   f[BE]   = f_BE   + xs*x_BE   + ys*y_BE   + zs*z_BE   + xs*ys*xy_BE   + xs*zs*xz_BE   + ys*zs*yz_BE   + feq[BE];
   f[TW]   = f_BE   + xs*x_BE   + ys*y_BE   + zs*z_BE   + xs*ys*xy_BE   + xs*zs*xz_BE   + ys*zs*yz_BE   + feq[TW];
   f[TN]   = f_TN   + xs*x_TN   + ys*y_TN   + zs*z_TN   + xs*ys*xy_TN   + xs*zs*xz_TN   + ys*zs*yz_TN   + feq[TN];
   f[BS]   = f_TN   + xs*x_TN   + ys*y_TN   + zs*z_TN   + xs*ys*xy_TN   + xs*zs*xz_TN   + ys*zs*yz_TN   + feq[BS];
   f[BN]   = f_BN   + xs*x_BN   + ys*y_BN   + zs*z_BN   + xs*ys*xy_BN   + xs*zs*xz_BN   + ys*zs*yz_BN   + feq[BN];
   f[TS]   = f_BN   + xs*x_BN   + ys*y_BN   + zs*z_BN   + xs*ys*xy_BN   + xs*zs*xz_BN   + ys*zs*yz_BN   + feq[TS];
   f[TNE]  = f_TNE  + xs*x_TNE  + ys*y_TNE  + zs*z_TNE  + xs*ys*xy_TNE  + xs*zs*xz_TNE  + ys*zs*yz_TNE  + feq[TNE];
   f[TSW]  = f_TSW  + xs*x_TSW  + ys*y_TSW  + zs*z_TSW  + xs*ys*xy_TSW  + xs*zs*xz_TSW  + ys*zs*yz_TSW  + feq[TSW];
   f[TSE]  = f_TSE  + xs*x_TSE  + ys*y_TSE  + zs*z_TSE  + xs*ys*xy_TSE  + xs*zs*xz_TSE  + ys*zs*yz_TSE  + feq[TSE];
   f[TNW]  = f_TNW  + xs*x_TNW  + ys*y_TNW  + zs*z_TNW  + xs*ys*xy_TNW  + xs*zs*xz_TNW  + ys*zs*yz_TNW  + feq[TNW];
   f[BNE]  = f_TSW  + xs*x_TSW  + ys*y_TSW  + zs*z_TSW  + xs*ys*xy_TSW  + xs*zs*xz_TSW  + ys*zs*yz_TSW  + feq[BNE];
   f[BSW]  = f_TNE  + xs*x_TNE  + ys*y_TNE  + zs*z_TNE  + xs*ys*xy_TNE  + xs*zs*xz_TNE  + ys*zs*yz_TNE  + feq[BSW];
   f[BSE]  = f_TNW  + xs*x_TNW  + ys*y_TNW  + zs*z_TNW  + xs*ys*xy_TNW  + xs*zs*xz_TNW  + ys*zs*yz_TNW  + feq[BSE];
   f[BNW]  = f_TSE  + xs*x_TSE  + ys*y_TSE  + zs*z_TSE  + xs*ys*xy_TSE  + xs*zs*xz_TSE  + ys*zs*yz_TSE  + feq[BNW];
   f[ZERO] = f_ZERO + xs*x_ZERO + ys*y_ZERO + zs*z_ZERO                                                 + feq[ZERO];
}
//////////////////////////////////////////////////////////////////////////
//Position SWB -0.25, -0.25, -0.25
LBMReal D3Q27OffsetInterpolationProcessor::calcPressBSW()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressTSW()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressTSE()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressBSE()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressBNW()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressTNW()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressTNE()
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
LBMReal D3Q27OffsetInterpolationProcessor::calcPressBNE()
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
void D3Q27OffsetInterpolationProcessor::calcInterpolatedNodeFC(LBMReal* f, LBMReal omega)
{
   using namespace D3Q27System;

   LBMReal press  =  press_NET * (0.125 - 0.25 * xoff - 0.25 * yoff - 0.25 * zoff) +
                     press_NWT * (0.125 + 0.25 * xoff - 0.25 * yoff - 0.25 * zoff) +
                     press_SET * (0.125 - 0.25 * xoff + 0.25 * yoff - 0.25 * zoff) +
                     press_SWT * (0.125 + 0.25 * xoff + 0.25 * yoff - 0.25 * zoff) +
                     press_NEB * (0.125 - 0.25 * xoff - 0.25 * yoff + 0.25 * zoff) +
                     press_NWB * (0.125 + 0.25 * xoff - 0.25 * yoff + 0.25 * zoff) +
                     press_SEB * (0.125 - 0.25 * xoff + 0.25 * yoff + 0.25 * zoff) +
                     press_SWB * (0.125 + 0.25 * xoff + 0.25 * yoff + 0.25 * zoff);
   LBMReal vx1  = a0;
   LBMReal vx2  = b0;
   LBMReal vx3  = c0;

   LBMReal rho = press ;//+ (ax+by+cz)/3.;

   LBMReal feq[ENDF+1];
   calcFeqFct(feq,rho,vx1,vx2,vx3);

   LBMReal eps_new = 2.;
   LBMReal o  = omega;
   LBMReal op = 1.;

   f_E    = eps_new *((5.*ax*o + 5.*by*o + 5.*cz*o - 8.*ax*op + 4.*by*op + 4.*cz*op)/(54.*o*op));
   f_N    = f_E + eps_new *((2.*(ax - by))/(9.*o));
   f_T    = f_E + eps_new *((2.*(ax - cz))/(9.*o));
   f_NE   = eps_new *(-(5.*cz*o + 3.*(ay + bx)*op - 2.*cz*op + ax*(5.*o + op) + by*(5.*o + op))/(54.*o*op));
   f_SE   = f_NE + eps_new *((  ay + bx )/(9.*o));
   f_TE   = eps_new *(-(5.*cz*o + by*(5.*o - 2.*op) + 3.*(az + cx)*op + cz*op + ax*(5.*o + op))/(54.*o*op));
   f_BE   = f_TE + eps_new *((  az + cx )/(9.*o));
   f_TN   = eps_new *(-(5.*ax*o + 5.*by*o + 5.*cz*o - 2.*ax*op + by*op + 3.*bz*op + 3.*cy*op + cz*op)/(54.*o*op));
   f_BN   = f_TN + eps_new *((  bz + cy )/(9.*o));
   f_ZERO = eps_new *((5.*(ax + by + cz))/(9.*op));
   f_TNE  = eps_new *(-(ay + az + bx + bz + cx + cy)/(72.*o));
   f_TSW  = - eps_new *((ay + bx)/(36.*o)) - f_TNE;
   f_TSE  = - eps_new *((az + cx)/(36.*o)) - f_TNE;
   f_TNW  = - eps_new *((bz + cy)/(36.*o)) - f_TNE;


   f[E]    = f_E    + feq[E];
   f[W]    = f_E    + feq[W];
   f[N]    = f_N    + feq[N];
   f[S]    = f_N    + feq[S];
   f[T]    = f_T    + feq[T];
   f[B]    = f_T    + feq[B];
   f[NE]   = f_NE   + feq[NE];
   f[SW]   = f_NE   + feq[SW];
   f[SE]   = f_SE   + feq[SE];
   f[NW]   = f_SE   + feq[NW];
   f[TE]   = f_TE   + feq[TE];
   f[BW]   = f_TE   + feq[BW];
   f[BE]   = f_BE   + feq[BE];
   f[TW]   = f_BE   + feq[TW];
   f[TN]   = f_TN   + feq[TN];
   f[BS]   = f_TN   + feq[BS];
   f[BN]   = f_BN   + feq[BN];
   f[TS]   = f_BN   + feq[TS];
   f[TNE]  = f_TNE  + feq[TNE];
   f[TNW]  = f_TNW  + feq[TNW];
   f[TSE]  = f_TSE  + feq[TSE];
   f[TSW]  = f_TSW  + feq[TSW];
   f[BNE]  = f_TSW  + feq[BNE];
   f[BNW]  = f_TSE  + feq[BNW];
   f[BSE]  = f_TNW  + feq[BSE];
   f[BSW]  = f_TNE  + feq[BSW];
   f[ZERO] = f_ZERO + feq[ZERO];
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::calcInterpolatedVelocity(LBMReal x, LBMReal y, LBMReal z, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3)
{
   vx1  = a0 + ax*x + ay*y + az*z + axx*x*x + ayy*y*y + azz*z*z + axy*x*y + axz*x*z + ayz*y*z+axyz*x*y*z;
   vx2  = b0 + bx*x + by*y + bz*z + bxx*x*x + byy*y*y + bzz*z*z + bxy*x*y + bxz*x*z + byz*y*z+bxyz*x*y*z;
   vx3  = c0 + cx*x + cy*y + cz*z + cxx*x*x + cyy*y*y + czz*z*z + cxy*x*y + cxz*x*z + cyz*y*z+cxyz*x*y*z;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27OffsetInterpolationProcessor::calcInterpolatedShearStress(LBMReal x, LBMReal y, LBMReal z,LBMReal& tauxx, LBMReal& tauyy, LBMReal& tauzz,LBMReal& tauxy, LBMReal& tauxz, LBMReal& tauyz)
{
   tauxx=ax+2*axx*x+axy*y+axz*z+axyz*y*z;
   tauyy=by+2*byy*y+bxy*x+byz*z+bxyz*x*z;
   tauzz=cz+2*czz*z+cxz*x+cyz*y+cxyz*x*y;
   tauxy=0.5*((ay+2.0*ayy*y+axy*x+ayz*z+axyz*x*z)+(bx+2.0*bxx*x+bxy*y+bxz*z+bxyz*y*z));
   tauxz=0.5*((az+2.0*azz*z+axz*x+ayz*y+axyz*x*y)+(cx+2.0*cxx*x+cxy*y+cxz*z+cxyz*y*z));
   tauyz=0.5*((bz+2.0*bzz*z+bxz*x+byz*y+bxyz*x*y)+(cy+2.0*cyy*y+cxy*x+cyz*z+cxyz*x*z));
}
//////////////////////////////////////////////////////////////////////////

