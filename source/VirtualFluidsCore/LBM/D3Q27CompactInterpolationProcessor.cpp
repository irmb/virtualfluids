#include "D3Q27CompactInterpolationProcessor.h"
#include "SimulationParameters.h"

#include <boost/foreach.hpp>

D3Q27CompactInterpolationProcessor::D3Q27CompactInterpolationProcessor()
                              : omegaC(0.0), omegaF(0.0)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
D3Q27CompactInterpolationProcessor::D3Q27CompactInterpolationProcessor(LBMReal omegaC, LBMReal omegaF)
                                                          : omegaC(omegaC), omegaF(omegaF)
{
   init();
}
//////////////////////////////////////////////////////////////////////////
D3Q27CompactInterpolationProcessor::~D3Q27CompactInterpolationProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
D3Q27InterpolationProcessorPtr D3Q27CompactInterpolationProcessor::clone()
{
   return D3Q27InterpolationProcessorPtr (new D3Q27CompactInterpolationProcessor(this->omegaC, this->omegaF));
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::init()
{
   calcFeqsForDirFct = NULL;
   calcMacrosFct     = NULL;
   calcFeqFct        = NULL;
   SimulationParametersPtr param = SimulationParameters::getInstanz();

   if(param->isCompressibleModel())
   {
      calcMacrosFct     = &D3Q27System::calcCompMacroscopicValues;
      calcFeqFct        = &D3Q27System::calcCompFeq; 
   }
   else
   {
      calcMacrosFct     = &D3Q27System::calcIncompMacroscopicValues;
      calcFeqFct        = &D3Q27System::calcIncompFeq; 
   }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::setOmegas( LBMReal omegaC, LBMReal omegaF )
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF)
{
   calcInterpolatedCoefficiets(icellC, omegaC);
   calcInterpolatedNode(icellF.BSW, omegaF, -0.25, -0.25, -0.25);
   calcInterpolatedNode(icellF.BNE, omegaF,  0.25,  0.25, -0.25);
   calcInterpolatedNode(icellF.TNW, omegaF, -0.25,  0.25,  0.25);
   calcInterpolatedNode(icellF.TSE, omegaF,  0.25, -0.25,  0.25);
   calcInterpolatedNode(icellF.BNW, omegaF, -0.25,  0.25, -0.25);
   calcInterpolatedNode(icellF.BSE, omegaF,  0.25, -0.25, -0.25);
   calcInterpolatedNode(icellF.TSW, omegaF, -0.25, -0.25,  0.25);
   calcInterpolatedNode(icellF.TNE, omegaF,  0.25,  0.25,  0.25);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC)
{
   calcInterpolatedCoefficiets(icellF, omegaF);
   calcInterpolatedNode(icellC, omegaC, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::interpolate8to1(D3Q27ICell& icellF, LBMReal* icellC, double x1, double x2, double x3, LBMReal omega)
{
   calcInterpolatedCoefficiets(icellF, omega);
   calcInterpolatedNode(icellC, omega, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::calcMoments(const LBMReal* const f, LBMReal omega, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3, 
                                              LBMReal& kxy, LBMReal& kyz, LBMReal& kxz, LBMReal& kxxMyy, LBMReal& kxxMzz)
{
   using namespace D3Q27System;

   LBMReal fd[ENDF+1];
   LBMReal feq[ENDF+1];

   calcMacrosFct(f,rho,vx1,vx2,vx3);
   calcFeqFct(feq,rho,vx1,vx2,vx3);

   fd[E]    = f[E] - feq[E];
   fd[W]    = f[W] - feq[W];
   fd[N]    = f[N] - feq[N];
   fd[S]    = f[S] - feq[S];
   fd[T]    = f[T] - feq[T];
   fd[B]    = f[B] - feq[B];
   fd[NE]   = f[NE] - feq[NE];
   fd[SW]   = f[SW] - feq[SW];
   fd[SE]   = f[SE] - feq[SE];
   fd[NW]   = f[NW] - feq[NW];
   fd[TE]   = f[TE] - feq[TE];
   fd[BW]   = f[BW] - feq[BW];
   fd[BE]   = f[BE] - feq[BE];
   fd[TW]   = f[TW] - feq[TW];
   fd[TN]   = f[TN] - feq[TN];
   fd[BS]   = f[BS] - feq[BS];
   fd[BN]   = f[BN] - feq[BN];
   fd[TS]   = f[TS] - feq[TS];
   fd[TNE]  = f[TNE] - feq[TNE];
   fd[TNW]  = f[TNW] - feq[TNW];
   fd[TSE]  = f[TSE] - feq[TSE];
   fd[TSW]  = f[TSW] - feq[TSW];
   fd[BNE]  = f[BNE] - feq[BNE];
   fd[BNW]  = f[BNW] - feq[BNW];
   fd[BSE]  = f[BSE] - feq[BSE];
   fd[BSW]  = f[BSW] - feq[BSW];
   fd[ZERO] = f[ZERO] - feq[ZERO];

   kxy    = -3.*omega*(fd[SW]+fd[TSW]+fd[BSW]-fd[NW]-fd[TNW]-fd[BNW]-fd[SE]-fd[TSE]-fd[BSE]+fd[NE]+fd[TNE]+fd[BNE]);
   kyz    = -3.*omega*(fd[BS]+fd[BSW]+fd[BSE]-fd[TS]-fd[TSW]-fd[TSE]-fd[BN]-fd[BNW]-fd[BNE]+fd[TN]+fd[TNW]+fd[TNE]);
   kxz    = -3.*omega*(fd[BW]+fd[BNW]+fd[BSW]-fd[TW]-fd[TNW]-fd[TSW]-fd[BE]-fd[BNE]-fd[BSE]+fd[TE]+fd[TNE]+fd[TSE]);
   kxxMyy = -3./2.*omega*(fd[BW]+fd[W]+fd[TW]-fd[BS]-fd[S]-fd[TS]-fd[BN]-fd[N]-fd[TN]+fd[BE]+fd[E]+fd[TE]);
   kxxMzz = -3./2.*omega*(fd[SW]+fd[W]+fd[NW]-fd[BS]-fd[TS]-fd[B]-fd[T]-fd[BN]-fd[TN]+fd[SE]+fd[E]+fd[NE]);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::calcInterpolatedCoefficiets(const D3Q27ICell& icell, LBMReal omega)
{
   LBMReal epsylon = 1.0;
   
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
   
   calcMoments(icell.TSW,omega,drho_SWT,vx1_SWT,vx2_SWT,vx3_SWT, kxyFromfcNEQ_SWT, kyzFromfcNEQ_SWT, kxzFromfcNEQ_SWT, kxxMyyFromfcNEQ_SWT, kxxMzzFromfcNEQ_SWT);
   calcMoments(icell.TNW,omega,drho_NWT,vx1_NWT,vx2_NWT,vx3_NWT, kxyFromfcNEQ_NWT, kyzFromfcNEQ_NWT, kxzFromfcNEQ_NWT, kxxMyyFromfcNEQ_NWT, kxxMzzFromfcNEQ_NWT);
   calcMoments(icell.TNE,omega,drho_NET,vx1_NET,vx2_NET,vx3_NET, kxyFromfcNEQ_NET, kyzFromfcNEQ_NET, kxzFromfcNEQ_NET, kxxMyyFromfcNEQ_NET, kxxMzzFromfcNEQ_NET);
   calcMoments(icell.TSE,omega,drho_SET,vx1_SET,vx2_SET,vx3_SET, kxyFromfcNEQ_SET, kyzFromfcNEQ_SET, kxzFromfcNEQ_SET, kxxMyyFromfcNEQ_SET, kxxMzzFromfcNEQ_SET);
   calcMoments(icell.BSW,omega,drho_SWB,vx1_SWB,vx2_SWB,vx3_SWB, kxyFromfcNEQ_SWB, kyzFromfcNEQ_SWB, kxzFromfcNEQ_SWB, kxxMyyFromfcNEQ_SWB, kxxMzzFromfcNEQ_SWB);
   calcMoments(icell.BNW,omega,drho_NWB,vx1_NWB,vx2_NWB,vx3_NWB, kxyFromfcNEQ_NWB, kyzFromfcNEQ_NWB, kxzFromfcNEQ_NWB, kxxMyyFromfcNEQ_NWB, kxxMzzFromfcNEQ_NWB);
   calcMoments(icell.BNE,omega,drho_NEB,vx1_NEB,vx2_NEB,vx3_NEB, kxyFromfcNEQ_NEB, kyzFromfcNEQ_NEB, kxzFromfcNEQ_NEB, kxxMyyFromfcNEQ_NEB, kxxMzzFromfcNEQ_NEB);
   calcMoments(icell.BSE,omega,drho_SEB,vx1_SEB,vx2_SEB,vx3_SEB, kxyFromfcNEQ_SEB, kyzFromfcNEQ_SEB, kxzFromfcNEQ_SEB, kxxMyyFromfcNEQ_SEB, kxxMzzFromfcNEQ_SEB);

   a0 = 8.*vx1_SWB + 8.*vx1_NWT + 8.*vx1_SET + 8.*vx1_NEB + 2.*vx2_SWB - 2.*vx2_NWT - 2.*vx2_SET + 2.*vx2_NEB + 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET - 2.*vx3_NEB + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + 2.*kxyFromfcNEQ_SWB - 2.*kxyFromfcNEQ_NWT + 2.*kxyFromfcNEQ_SET - 2.*kxyFromfcNEQ_NEB + 2.*kxzFromfcNEQ_SWB - 2.*kxzFromfcNEQ_NWT - 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB;
   ax = (-4.*vx1_SWB - 4.*vx1_NWT + 4.*vx1_SET + 4.*vx1_NEB - kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB - kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT - kyzFromfcNEQ_SET - kyzFromfcNEQ_NEB)/epsylon;
   ay = (-2.*vx1_SWB + 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB + 2.*vx2_SWB + 2.*vx2_NWT - 2.*vx2_SET - 2.*vx2_NEB - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB)/epsylon;
   az = (-2.*vx1_SWB + 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB + 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET - 2.*vx3_NEB - kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB)/epsylon;
   axx= (2.*vx2_SWB - 2.*vx2_NWT - 2.*vx2_SET + 2.*vx2_NEB + 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET - 2.*vx3_NEB - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   ayy= (-2.*vx2_SWB + 2.*vx2_NWT + 2.*vx2_SET - 2.*vx2_NEB - 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET + 2.*vx3_NEB - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB - 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_NWT - 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT - kyzFromfcNEQ_SET - kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   azz= (-2.*vx2_SWB + 2.*vx2_NWT + 2.*vx2_SET - 2.*vx2_NEB - 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET + 2.*vx3_NEB + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB - 2.*kxzFromfcNEQ_SWB + 2.*kxzFromfcNEQ_NWT + 2.*kxzFromfcNEQ_SET - 2.*kxzFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT - kyzFromfcNEQ_SET - kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   axy= (2.*vx1_SWB - 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB + 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET - 2.*vx3_NEB - kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB)/(epsylon*epsylon);
   axz= (2.*vx1_SWB - 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB + 2.*vx2_SWB + 2.*vx2_NWT - 2.*vx2_SET - 2.*vx2_NEB - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB)/(epsylon*epsylon);
   ayz= (-kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB - kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT - kyzFromfcNEQ_SET - kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   b0 = 2.*vx1_SWB - 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB + 8.*vx2_SWB + 8.*vx2_NWT + 8.*vx2_SET + 8.*vx2_NEB + 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET - 2.*vx3_NEB - 2.*kxxMyyFromfcNEQ_SWB + 2.*kxxMyyFromfcNEQ_NWT - 2.*kxxMyyFromfcNEQ_SET + 2.*kxxMyyFromfcNEQ_NEB + kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + 2.*kxyFromfcNEQ_SWB + 2.*kxyFromfcNEQ_NWT - 2.*kxyFromfcNEQ_SET - 2.*kxyFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB + 2.*kyzFromfcNEQ_SWB - 2.*kyzFromfcNEQ_NWT - 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_NEB;
   bx = (2.*vx1_SWB - 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB - 2.*vx2_SWB - 2.*vx2_NWT + 2.*vx2_SET + 2.*vx2_NEB + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB)/epsylon;
   by = (-4.*vx2_SWB + 4.*vx2_NWT - 4.*vx2_SET + 4.*vx2_NEB - kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB + kxzFromfcNEQ_SWB - kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET - kxzFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/epsylon;
   bz = (-2.*vx2_SWB + 2.*vx2_NWT + 2.*vx2_SET - 2.*vx2_NEB + 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET - 2.*vx3_NEB + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/epsylon;
   bxx= (-2.*vx1_SWB + 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB - 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET + 2.*vx3_NEB + kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB - 2.*kxyFromfcNEQ_SWB - 2.*kxyFromfcNEQ_NWT + 2.*kxyFromfcNEQ_SET + 2.*kxyFromfcNEQ_NEB - kxzFromfcNEQ_SWB - kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET - kxzFromfcNEQ_NEB)/(epsylon*epsylon);
   byy= (2.*vx1_SWB - 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB + 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET - 2.*vx3_NEB + 2.*kxxMyyFromfcNEQ_SWB - 2.*kxxMyyFromfcNEQ_NWT + 2.*kxxMyyFromfcNEQ_SET - 2.*kxxMyyFromfcNEQ_NEB - kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB)/(epsylon*epsylon);
   bzz= (-2.*vx1_SWB + 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB - 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET + 2.*vx3_NEB - kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB - kxzFromfcNEQ_SWB - kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET - kxzFromfcNEQ_NEB - 2.*kyzFromfcNEQ_SWB + 2.*kyzFromfcNEQ_NWT + 2.*kyzFromfcNEQ_SET - 2.*kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   bxy= (2.*vx2_SWB - 2.*vx2_NWT - 2.*vx2_SET + 2.*vx2_NEB + 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET - 2.*vx3_NEB + kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB - kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET + kxxMzzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   bxz= (-kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB + kxzFromfcNEQ_SWB - kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET - kxzFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   byz= (2.*vx1_SWB - 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB + 2.*vx2_SWB + 2.*vx2_NWT - 2.*vx2_SET - 2.*vx2_NEB + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB)/(epsylon*epsylon);
   c0 = 2.*vx1_SWB - 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB + 2.*vx2_SWB + 2.*vx2_NWT - 2.*vx2_SET - 2.*vx2_NEB + 8.*vx3_SWB + 8.*vx3_NWT + 8.*vx3_SET + 8.*vx3_NEB + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB - 2.*kxxMzzFromfcNEQ_SWB + 2.*kxxMzzFromfcNEQ_NWT + 2.*kxxMzzFromfcNEQ_SET - 2.*kxxMzzFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB + 2.*kxzFromfcNEQ_SWB + 2.*kxzFromfcNEQ_NWT - 2.*kxzFromfcNEQ_SET - 2.*kxzFromfcNEQ_NEB + 2.*kyzFromfcNEQ_SWB - 2.*kyzFromfcNEQ_NWT + 2.*kyzFromfcNEQ_SET - 2.*kyzFromfcNEQ_NEB;
   cx = (2.*vx1_SWB - 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB - 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET + 2.*vx3_NEB + kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB)/epsylon;
   cy = (2.*vx2_SWB - 2.*vx2_NWT - 2.*vx2_SET + 2.*vx2_NEB - 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET + 2.*vx3_NEB - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/epsylon;
   cz = (-4.*vx3_SWB + 4.*vx3_NWT + 4.*vx3_SET - 4.*vx3_NEB + kxyFromfcNEQ_SWB - kxyFromfcNEQ_NWT - kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB - kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/epsylon;
   cxx= (-2.*vx1_SWB + 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB - 2.*vx2_SWB - 2.*vx2_NWT + 2.*vx2_SET + 2.*vx2_NEB + kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT - kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB - kxyFromfcNEQ_SWB - kxyFromfcNEQ_NWT - kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB - 2.*kxzFromfcNEQ_SWB - 2.*kxzFromfcNEQ_NWT + 2.*kxzFromfcNEQ_SET + 2.*kxzFromfcNEQ_NEB)/(epsylon*epsylon);
   cyy= (-2.*vx1_SWB + 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB - 2.*vx2_SWB - 2.*vx2_NWT + 2.*vx2_SET + 2.*vx2_NEB - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB - kxyFromfcNEQ_SWB - kxyFromfcNEQ_NWT - kxyFromfcNEQ_SET - kxyFromfcNEQ_NEB - 2.*kyzFromfcNEQ_SWB + 2.*kyzFromfcNEQ_NWT - 2.*kyzFromfcNEQ_SET + 2.*kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   czz= (2.*vx1_SWB - 2.*vx1_NWT + 2.*vx1_SET - 2.*vx1_NEB + 2.*vx2_SWB + 2.*vx2_NWT - 2.*vx2_SET - 2.*vx2_NEB - kxxMyyFromfcNEQ_SWB + kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET - kxxMyyFromfcNEQ_NEB + 2.*kxxMzzFromfcNEQ_SWB - 2.*kxxMzzFromfcNEQ_NWT - 2.*kxxMzzFromfcNEQ_SET + 2.*kxxMzzFromfcNEQ_NEB + kxyFromfcNEQ_SWB + kxyFromfcNEQ_NWT + kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB)/(epsylon*epsylon);
   cxy= (kxyFromfcNEQ_SWB - kxyFromfcNEQ_NWT - kxyFromfcNEQ_SET + kxyFromfcNEQ_NEB - kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT - kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB - kyzFromfcNEQ_SWB - kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   cxz= (2.*vx2_SWB - 2.*vx2_NWT - 2.*vx2_SET + 2.*vx2_NEB + 2.*vx3_SWB - 2.*vx3_NWT + 2.*vx3_SET - 2.*vx3_NEB - kxxMyyFromfcNEQ_SWB - kxxMyyFromfcNEQ_NWT + kxxMyyFromfcNEQ_SET + kxxMyyFromfcNEQ_NEB + kxxMzzFromfcNEQ_SWB + kxxMzzFromfcNEQ_NWT - kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + kyzFromfcNEQ_SWB + kyzFromfcNEQ_NWT + kyzFromfcNEQ_SET + kyzFromfcNEQ_NEB)/(epsylon*epsylon);
   cyz= (2.*vx1_SWB - 2.*vx1_NWT - 2.*vx1_SET + 2.*vx1_NEB + 2.*vx3_SWB + 2.*vx3_NWT - 2.*vx3_SET - 2.*vx3_NEB + kxxMzzFromfcNEQ_SWB - kxxMzzFromfcNEQ_NWT + kxxMzzFromfcNEQ_SET - kxxMzzFromfcNEQ_NEB + kxzFromfcNEQ_SWB + kxzFromfcNEQ_NWT + kxzFromfcNEQ_SET + kxzFromfcNEQ_NEB)/(epsylon*epsylon);
   //////////////////////////////////////////////////////////////////////////
   //another 4 points
   a0 += 8.*vx1_SWT + 8.*vx1_NWB + 8.*vx1_SEB + 8.*vx1_NET + 2.*vx2_SWT - 2.*vx2_NWB - 2.*vx2_SEB + 2.*vx2_NET - 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB + 2.*vx3_NET + kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET + kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET + 2.*kxyFromfcNEQ_SWT - 2.*kxyFromfcNEQ_NWB + 2.*kxyFromfcNEQ_SEB - 2.*kxyFromfcNEQ_NET - 2.*kxzFromfcNEQ_SWT + 2.*kxzFromfcNEQ_NWB + 2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET;
   ax += (-4.*vx1_SWT - 4.*vx1_NWB + 4.*vx1_SEB + 4.*vx1_NET - 1.*kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET + kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/epsylon;
   ay += (-2.*vx1_SWT + 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET + 2.*vx2_SWT + 2.*vx2_NWB - 2.*vx2_SEB - 2.*vx2_NET - 1.*kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET + kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET)/epsylon;
   az += (2.*vx1_SWT - 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET + 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB - 2.*vx3_NET + kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET + kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET)/epsylon;
   axx += (2.*vx2_SWT - 2.*vx2_NWB - 2.*vx2_SEB + 2.*vx2_NET - 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB + 2.*vx3_NET - 1.*kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET - 1.*kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   ayy += (-2.*vx2_SWT + 2.*vx2_NWB + 2.*vx2_SEB - 2.*vx2_NET + 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB - 2.*vx3_NET - 1.*kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET + kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET - 2.*kxyFromfcNEQ_SWT + 2.*kxyFromfcNEQ_NWB - 2.*kxyFromfcNEQ_SEB + 2.*kxyFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/(epsylon*epsylon);
   azz += (-2.*vx2_SWT + 2.*vx2_NWB + 2.*vx2_SEB - 2.*vx2_NET + 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB - 2.*vx3_NET + kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET - 1.*kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET + 2.*kxzFromfcNEQ_SWT - 2.*kxzFromfcNEQ_NWB - 2.*kxzFromfcNEQ_SEB + 2.*kxzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/(epsylon*epsylon);
   axy += (2.*vx1_SWT - 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB + 2.*vx3_NET - 1.*kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET)/(epsylon*epsylon);
   axz += (-2.*vx1_SWT + 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx2_SWT - 2.*vx2_NWB + 2.*vx2_SEB + 2.*vx2_NET + kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET - 1.*kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET)/(epsylon*epsylon);
   ayz += (kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   b0 += 2.*vx1_SWT - 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET + 8.*vx2_SWT + 8.*vx2_NWB + 8.*vx2_SEB + 8.*vx2_NET - 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB + 2.*vx3_NET - 2.*kxxMyyFromfcNEQ_SWT + 2.*kxxMyyFromfcNEQ_NWB - 2.*kxxMyyFromfcNEQ_SEB + 2.*kxxMyyFromfcNEQ_NET + kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET + 2.*kxyFromfcNEQ_SWT + 2.*kxyFromfcNEQ_NWB - 2.*kxyFromfcNEQ_SEB - 2.*kxyFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET - 2.*kyzFromfcNEQ_SWT + 2.*kyzFromfcNEQ_NWB + 2.*kyzFromfcNEQ_SEB - 2.*kyzFromfcNEQ_NET;
   bx += (2.*vx1_SWT - 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET - 2.*vx2_SWT - 2.*vx2_NWB + 2.*vx2_SEB + 2.*vx2_NET + kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET + kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET)/epsylon;
   by += (-4.*vx2_SWT + 4.*vx2_NWB - 4.*vx2_SEB + 4.*vx2_NET - 1.*kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/epsylon;
   bz += (2.*vx2_SWT - 2.*vx2_NWB - 2.*vx2_SEB + 2.*vx2_NET + 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB - 2.*vx3_NET - 1.*kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET + kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/epsylon;
   bxx += (-2.*vx1_SWT + 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET + 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB - 2.*vx3_NET + kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET - 2.*kxyFromfcNEQ_SWT - 2.*kxyFromfcNEQ_NWB + 2.*kxyFromfcNEQ_SEB + 2.*kxyFromfcNEQ_NET + kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET)/(epsylon*epsylon);
   byy += (2.*vx1_SWT - 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB + 2.*vx3_NET + 2.*kxxMyyFromfcNEQ_SWT - 2.*kxxMyyFromfcNEQ_NWB + 2.*kxxMyyFromfcNEQ_SEB - 2.*kxxMyyFromfcNEQ_NET - 1.*kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET)/(epsylon*epsylon);
   bzz += (-2.*vx1_SWT + 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET + 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB - 2.*vx3_NET - 1.*kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET + kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET + 2.*kyzFromfcNEQ_SWT - 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_SEB + 2.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   bxy += (2.*vx2_SWT - 2.*vx2_NWB - 2.*vx2_SEB + 2.*vx2_NET - 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB + 2.*vx3_NET + kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET - 1.*kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   bxz += (kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET + kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/(epsylon*epsylon);
   byz += (-2.*vx1_SWT + 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx2_SWT - 2.*vx2_NWB + 2.*vx2_SEB + 2.*vx2_NET - 1.*kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET - 1.*kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET)/(epsylon*epsylon);
   c0 += -2.*vx1_SWT + 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx2_SWT - 2.*vx2_NWB + 2.*vx2_SEB + 2.*vx2_NET + 8.*vx3_SWT + 8.*vx3_NWB + 8.*vx3_SEB + 8.*vx3_NET - 1.*kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET + 2.*kxxMzzFromfcNEQ_SWT - 2.*kxxMzzFromfcNEQ_NWB - 2.*kxxMzzFromfcNEQ_SEB + 2.*kxxMzzFromfcNEQ_NET - 1.*kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET + 2.*kxzFromfcNEQ_SWT + 2.*kxzFromfcNEQ_NWB - 2.*kxzFromfcNEQ_SEB - 2.*kxzFromfcNEQ_NET + 2.*kyzFromfcNEQ_SWT - 2.*kyzFromfcNEQ_NWB + 2.*kyzFromfcNEQ_SEB - 2.*kyzFromfcNEQ_NET;
   cx += (-2.*vx1_SWT + 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET - 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB + 2.*vx3_NET - 1.*kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET + kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET)/epsylon;
   cy += (-2.*vx2_SWT + 2.*vx2_NWB + 2.*vx2_SEB - 2.*vx2_NET - 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB + 2.*vx3_NET + kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET - 1.*kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB + kxxMzzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/epsylon;
   cz += (4.*vx3_SWT - 4.*vx3_NWB - 4.*vx3_SEB + 4.*vx3_NET + kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET + kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB + kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET + kyzFromfcNEQ_SWT + kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/epsylon;
   cxx += (2.*vx1_SWT - 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET + 2.*vx2_SWT + 2.*vx2_NWB - 2.*vx2_SEB - 2.*vx2_NET - 1.*kxxMyyFromfcNEQ_SWT + kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB - 1.*kxxMyyFromfcNEQ_NET + kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET - 2.*kxzFromfcNEQ_SWT - 2.*kxzFromfcNEQ_NWB + 2.*kxzFromfcNEQ_SEB + 2.*kxzFromfcNEQ_NET)/(epsylon*epsylon);
   cyy += (2.*vx1_SWT - 2.*vx1_NWB + 2.*vx1_SEB - 2.*vx1_NET + 2.*vx2_SWT + 2.*vx2_NWB - 2.*vx2_SEB - 2.*vx2_NET + kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET + kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB + kxyFromfcNEQ_NET - 2.*kyzFromfcNEQ_SWT + 2.*kyzFromfcNEQ_NWB - 2.*kyzFromfcNEQ_SEB + 2.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   czz += (-2.*vx1_SWT + 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx2_SWT - 2.*vx2_NWB + 2.*vx2_SEB + 2.*vx2_NET + kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB - 1.*kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET - 2.*kxxMzzFromfcNEQ_SWT + 2.*kxxMzzFromfcNEQ_NWB + 2.*kxxMzzFromfcNEQ_SEB - 2.*kxxMzzFromfcNEQ_NET - 1.*kxyFromfcNEQ_SWT - 1.*kxyFromfcNEQ_NWB - 1.*kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET)/(epsylon*epsylon);
   cxy += (-1.*kxyFromfcNEQ_SWT + kxyFromfcNEQ_NWB + kxyFromfcNEQ_SEB - 1.*kxyFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT + kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB + kxzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB + kyzFromfcNEQ_SEB + kyzFromfcNEQ_NET)/(epsylon*epsylon);
   cxz += (2.*vx2_SWT - 2.*vx2_NWB - 2.*vx2_SEB + 2.*vx2_NET - 2.*vx3_SWT + 2.*vx3_NWB - 2.*vx3_SEB + 2.*vx3_NET - 1.*kxxMyyFromfcNEQ_SWT - 1.*kxxMyyFromfcNEQ_NWB + kxxMyyFromfcNEQ_SEB + kxxMyyFromfcNEQ_NET + kxxMzzFromfcNEQ_SWT + kxxMzzFromfcNEQ_NWB - 1.*kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET - 1.*kyzFromfcNEQ_SWT - 1.*kyzFromfcNEQ_NWB - 1.*kyzFromfcNEQ_SEB - 1.*kyzFromfcNEQ_NET)/(epsylon*epsylon);
   cyz += (2.*vx1_SWT - 2.*vx1_NWB - 2.*vx1_SEB + 2.*vx1_NET - 2.*vx3_SWT - 2.*vx3_NWB + 2.*vx3_SEB + 2.*vx3_NET + kxxMzzFromfcNEQ_SWT - 1.*kxxMzzFromfcNEQ_NWB + kxxMzzFromfcNEQ_SEB - 1.*kxxMzzFromfcNEQ_NET - 1.*kxzFromfcNEQ_SWT - 1.*kxzFromfcNEQ_NWB - 1.*kxzFromfcNEQ_SEB - 1.*kxzFromfcNEQ_NET)/(epsylon*epsylon);
   //////////////////////////////////////////////////////////////////////////
   //normalization
   a0  /= 2.;
   ax  /= 2.;
   ay  /= 2.;
   az  /= 2.;
   axx /= 2.;
   ayy /= 2.;
   azz /= 2.;
   axy /= 2.;
   axz /= 2.;
   ayz /= 2.;
   b0  /= 2.;
   bx  /= 2.;
   by  /= 2.;
   bz  /= 2.;
   bxx /= 2.;
   byy /= 2.;
   bzz /= 2.;
   bxy /= 2.;
   bxz /= 2.;
   byz /= 2.;
   c0  /= 2.;
   cx  /= 2.;
   cy  /= 2.;
   cz  /= 2.;
   cxx /= 2.;
   cyy /= 2.;
   czz /= 2.;
   cxy /= 2.;
   cxz /= 2.;
   cyz /= 2.;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   a0  /= 32.;
   ax  /= 8.;
   ay  /= 8.;
   az  /= 8.;
   axx /= 8.;
   ayy /= 8.;
   azz /= 8.;
   axy /= 4.;
   axz /= 4.;
   ayz /= 4.;
   b0  /= 32.;
   bx  /= 8.;
   by  /= 8.;
   bz  /= 8.;
   bxx /= 8.;
   byy /= 8.;
   bzz /= 8.;
   bxy /= 4.;
   bxz /= 4.;
   byz /= 4.;
   c0  /= 32.;
   cx  /= 8.;
   cy  /= 8.;
   cz  /= 8.;
   cxx /= 8.;
   cyy /= 8.;
   czz /= 8.;
   cxy /= 4.;
   cxz /= 4.;
   cyz /= 4.;
}
//////////////////////////////////////////////////////////////////////////
void D3Q27CompactInterpolationProcessor::calcInterpolatedNode(LBMReal* f, LBMReal omega, LBMReal x, LBMReal y, LBMReal z)
{
   using namespace D3Q27System;
   
   LBMReal eps_new = 0.5;
   LBMReal op = 1.0;
   const LBMReal o = omega;


   LBMReal drho = (drho_NEB*(1. + 2.*x + 2.*y + 4.*x*y - 2.*z - 4.*x*z - 4.*y*z - 8.*x*y*z) + 
                     drho_SET*(1. + 2.*x - 2.*y - 4.*x*y + 2.*z + 4.*x*z - 4.*y*z - 8.*x*y*z) + 
                     drho_NWT*(1. - 2.*x + 2.*y - 4.*x*y + 2.*z - 4.*x*z + 4.*y*z - 8.*x*y*z) + 
                     drho_SWB*(1. - 2.*x - 2.*y + 4.*x*y - 2.*z + 4.*x*z + 4.*y*z - 8.*x*y*z) + 
                     drho_SWT*(1. - 2.*x - 2.*y + 4.*x*y + 2.*z - 4.*x*z - 4.*y*z + 8.*x*y*z) + 
                     drho_NWB*(1. - 2.*x + 2.*y - 4.*x*y - 2.*z + 4.*x*z - 4.*y*z + 8.*x*y*z) + 
                     drho_SEB*(1. + 2.*x - 2.*y - 4.*x*y - 2.*z - 4.*x*z + 4.*y*z + 8.*x*y*z) + 
                     drho_NET*(1. + 2.*x + 2.*y + 4.*x*y + 2.*z + 4.*x*z + 4.*y*z + 8.*x*y*z))/8.;
   LBMReal vx1  = a0 + ax*x + ay*y + az*z + axx*x*x + ayy*y*y + azz*z*z + axy*x*y + axz*x*z + ayz*y*z;
   LBMReal vx2  = b0 + bx*x + by*y + bz*z + bxx*x*x + byy*y*y + bzz*z*z + bxy*x*y + bxz*x*z + byz*y*z;
   LBMReal vx3  = c0 + cx*x + cy*y + cz*z + cxx*x*x + cyy*y*y + czz*z*z + cxy*x*y + cxz*x*z + cyz*y*z;

   LBMReal feq[ENDF+1];
   calcFeqFct(feq,drho,vx1,vx2,vx3);

   f[E]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-by - cz - bxy*x - cxz*x - 2*byy*y - cyz*y - byz*z - 2*czz*z + 2*(ax + 2*axx*x + axy*y + axz*z)))/(9.*o))/2.;
   f[W]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-by - cz - bxy*x - cxz*x - 2*byy*y - cyz*y - byz*z - 2*czz*z + 2*(ax + 2*axx*x + axy*y + axz*z)))/(9.*o))/2.;
   f[N]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-ax - cz - 2*axx*x - cxz*x - axy*y - cyz*y - axz*z - 2*czz*z + 2*(by + bxy*x + 2*byy*y + byz*z)))/(9.*o))/2.;
   f[S]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-ax - cz - 2*axx*x - cxz*x - axy*y - cyz*y - axz*z - 2*czz*z + 2*(by + bxy*x + 2*byy*y + byz*z)))/(9.*o))/2.;
   f[T]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-ax - by - 2*axx*x - bxy*x - axy*y - 2*byy*y - axz*z - byz*z + 2*(cz + cxz*x + cyz*y + 2*czz*z)))/(9.*o))/2.;
   f[B]    = eps_new *((-5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(9.*op) - 
            (2*(-ax - by - 2*axx*x - bxy*x - axy*y - 2*byy*y - axz*z - byz*z + 2*(cz + cxz*x + cyz*y + 2*czz*z)))/(9.*o))/2.;
   f[NE]   = eps_new *(-(ay + bx + axy*x + 2*bxx*x + 2*ayy*y + bxy*y + ayz*z + bxz*z)/(12.*o));
   f[SW]   = eps_new *(-(ay + bx + axy*x + 2*bxx*x + 2*ayy*y + bxy*y + ayz*z + bxz*z)/(12.*o));
   f[SE]   = eps_new *(ay + bx + axy*x + 2*bxx*x + 2*ayy*y + bxy*y + ayz*z + bxz*z)/(12.*o);
   f[NW]   = eps_new *(ay + bx + axy*x + 2*bxx*x + 2*ayy*y + bxy*y + ayz*z + bxz*z)/(12.*o);
   f[TE]   = eps_new *(-(az + cx + axz*x + 2*cxx*x + ayz*y + cxy*y + 2*azz*z + cxz*z)/(12.*o));
   f[BW]   = eps_new *(-(az + cx + axz*x + 2*cxx*x + ayz*y + cxy*y + 2*azz*z + cxz*z)/(12.*o));
   f[BE]   = eps_new *(az + cx + axz*x + 2*cxx*x + ayz*y + cxy*y + 2*azz*z + cxz*z)/(12.*o);
   f[TW]   = eps_new *(az + cx + axz*x + 2*cxx*x + ayz*y + cxy*y + 2*azz*z + cxz*z)/(12.*o);
   f[TN]   = eps_new *(-(bz + cy + bxz*x + cxy*x + byz*y + 2*cyy*y + 2*bzz*z + cyz*z)/(12.*o));
   f[BS]   = eps_new *(-(bz + cy + bxz*x + cxy*x + byz*y + 2*cyy*y + 2*bzz*z + cyz*z)/(12.*o));
   f[BN]   = eps_new *(bz + cy + bxz*x + cxy*x + byz*y + 2*cyy*y + 2*bzz*z + cyz*z)/(12.*o);
   f[TS]   = eps_new *(bz + cy + bxz*x + cxy*x + byz*y + 2*cyy*y + 2*bzz*z + cyz*z)/(12.*o);
   f[ZERO] = eps_new *(5*(ax + by + cz + 2*axx*x + bxy*x + cxz*x + axy*y + 2*byy*y + cyz*y + axz*z + byz*z + 2*czz*z))/(3.*op) + 
            (2*(-by - cz - bxy*x - cxz*x - 2*byy*y - cyz*y - byz*z - 2*czz*z + 2*(ax + 2*axx*x + axy*y + axz*z)))/(9.*o) + 
            (2*(-ax - cz - 2*axx*x - cxz*x - axy*y - cyz*y - axz*z - 2*czz*z + 2*(by + bxy*x + 2*byy*y + byz*z)))/(9.*o) + 
            (2*(-ax - by - 2*axx*x - bxy*x - axy*y - 2*byy*y - axz*z - byz*z + 2*(cz + cxz*x + cyz*y + 2*czz*z)))/(9.*o);
   f[TNE]  = 0.0;
   f[TNW]  = 0.0;
   f[TSE]  = 0.0;
   f[TSW]  = 0.0;
   f[BNE]  = 0.0;
   f[BNW]  = 0.0;
   f[BSE]  = 0.0;
   f[BSW]  = 0.0;

   f[E]    += feq[E];
   f[W]    += feq[W];
   f[N]    += feq[N];
   f[S]    += feq[S];
   f[T]    += feq[T];
   f[B]    += feq[B];
   f[NE]   += feq[NE];
   f[SW]   += feq[SW];
   f[SE]   += feq[SE];
   f[NW]   += feq[NW];
   f[TE]   += feq[TE];
   f[BW]   += feq[BW];
   f[BE]   += feq[BE];
   f[TW]   += feq[TW];
   f[TN]   += feq[TN];
   f[BS]   += feq[BS];
   f[BN]   += feq[BN];
   f[TS]   += feq[TS];
   f[TNE]  += feq[TNE];
   f[TNW]  += feq[TNW];
   f[TSE]  += feq[TSE];
   f[TSW]  += feq[TSW];
   f[BNE]  += feq[BNE];
   f[BNW]  += feq[BNW];
   f[BSE]  += feq[BSE];
   f[BSW]  += feq[BSW];
   f[ZERO] += feq[ZERO];
}
//////////////////////////////////////////////////////////////////////////

