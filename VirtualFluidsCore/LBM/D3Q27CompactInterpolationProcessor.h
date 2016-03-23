#ifndef D3Q27CompactInterpolationProcessor_H_
#define D3Q27CompactInterpolationProcessor_H_

#include "D3Q27InterpolationProcessor.h"
#include "D3Q27System.h"


//////////////////////////////////////////////////////////////////////////
//it work only for cascaded LBM
//super compact interpolation method by Martin Geier
//////////////////////////////////////////////////////////////////////////

class D3Q27CompactInterpolationProcessor;
typedef boost::shared_ptr<D3Q27CompactInterpolationProcessor> D3Q27CompactInterpolationProcessorPtr;

class D3Q27CompactInterpolationProcessor : public D3Q27InterpolationProcessor
{
public:
   D3Q27CompactInterpolationProcessor();
   D3Q27CompactInterpolationProcessor(LBMReal omegaC, LBMReal omegaF);
   virtual ~D3Q27CompactInterpolationProcessor();
   D3Q27InterpolationProcessorPtr clone();
   void setOmegas(LBMReal omegaC, LBMReal omegaF);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff);
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC); 
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff); 
   void interpolate8to1(D3Q27ICell& icellF, LBMReal* icellC, double x1, double x2, double x3, LBMReal omega);
protected:
private:
   typedef void      (*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   typedef LBMReal (*CalcFeqForDirFct)(const int& /*direction*/,const LBMReal& /*(d)rho*/,const LBMReal& /*vx1*/,const LBMReal& /*vx2*/,const LBMReal& /*vx3*/);
   typedef  void     (*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/,const LBMReal& /*rho*/,const LBMReal& /*vx1*/,const LBMReal& /*vx2*/,const LBMReal& /*vx3*/);
   
   CalcFeqForDirFct calcFeqsForDirFct;
   CalcMacrosFct    calcMacrosFct;
   CalcFeqFct       calcFeqFct;
   LBMReal omegaC, omegaF;
   LBMReal a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz, b0, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz, c0, cx, cy, cz, cxx, cyy, czz, cxy, cxz, cyz;
   LBMReal drho_SWT, drho_NWT, drho_NET, drho_SET, drho_SWB, drho_NWB, drho_NEB, drho_SEB;

   void init();
   void calcMoments(const LBMReal* const f, LBMReal omega, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3, 
                           LBMReal& kxy, LBMReal& kyz, LBMReal& kxz, LBMReal& kxxMyy, LBMReal& kxxMzz);
   void calcInterpolatedCoefficiets(const D3Q27ICell& icell, LBMReal omega);
   void calcInterpolatedNode(LBMReal* f, LBMReal omega, LBMReal x, LBMReal y, LBMReal z);
};
//////////////////////////////////////////////////////////////////////////
inline void D3Q27CompactInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff)
{
   this->interpolateCoarseToFine(icellC, icellF);
}
//////////////////////////////////////////////////////////////////////////
inline void D3Q27CompactInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff)
{
   this->interpolateFineToCoarse(icellF, icellC);
}

#endif
