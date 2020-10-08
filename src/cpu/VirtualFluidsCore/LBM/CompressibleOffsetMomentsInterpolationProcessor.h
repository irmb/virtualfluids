#ifndef CompressibleOffsetMomentsInterpolationProcessor_H_
#define CompressibleOffsetMomentsInterpolationProcessor_H_

#include "InterpolationProcessor.h"
#include "D3Q27System.h"

//////////////////////////////////////////////////////////////////////////
//it works only for cascaded LBM
//super compact interpolation method by Martin Geier
//////////////////////////////////////////////////////////////////////////

class CompressibleOffsetMomentsInterpolationProcessor;

class CompressibleOffsetMomentsInterpolationProcessor : public InterpolationProcessor
{
public:
   CompressibleOffsetMomentsInterpolationProcessor();
   CompressibleOffsetMomentsInterpolationProcessor(LBMReal omegaC, LBMReal omegaF);
   ~CompressibleOffsetMomentsInterpolationProcessor() override;
   InterpolationProcessorPtr clone() override;
   void setOmegas(LBMReal omegaC, LBMReal omegaF) override;
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF) override;
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff) override;
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC) override; 
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff) override; 
   void setBulkViscosity(LBMReal shearViscosity, LBMReal bulkViscosity);
protected:   
private:
   LBMReal omegaC, omegaF;
   LBMReal a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz, b0, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz, c0, cx, cy, cz, cxx, cyy, czz, cxy, cxz, cyz, axyz, bxyz, cxyz;
   LBMReal xoff,    yoff,    zoff;
   LBMReal xoff_sq, yoff_sq, zoff_sq;
   LBMReal press_SWT, press_NWT, press_NET, press_SET, press_SWB, press_NWB, press_NEB, press_SEB;

   LBMReal  f_E,  f_N,  f_T,  f_NE,  f_SE,  f_BE,  f_TE,  f_TN,  f_BN,  f_TNE,  f_TNW,  f_TSE,  f_TSW,  f_ZERO;
   LBMReal  x_E,  x_N,  x_T,  x_NE,  x_SE,  x_BE,  x_TE,  x_TN,  x_BN,  x_TNE,  x_TNW,  x_TSE,  x_TSW,  x_ZERO;
   LBMReal  y_E,  y_N,  y_T,  y_NE,  y_SE,  y_BE,  y_TE,  y_TN,  y_BN,  y_TNE,  y_TNW,  y_TSE,  y_TSW,  y_ZERO;
   LBMReal  z_E,  z_N,  z_T,  z_NE,  z_SE,  z_BE,  z_TE,  z_TN,  z_BN,  z_TNE,  z_TNW,  z_TSE,  z_TSW,  z_ZERO;
   LBMReal xy_E, xy_N, xy_T, xy_NE, xy_SE, xy_BE, xy_TE, xy_TN, xy_BN, xy_TNE, xy_TNW, xy_TSE, xy_TSW/*, xy_ZERO*/;
   LBMReal xz_E, xz_N, xz_T, xz_NE, xz_SE, xz_BE, xz_TE, xz_TN, xz_BN, xz_TNE, xz_TNW, xz_TSE, xz_TSW/*, xz_ZERO*/;
   LBMReal yz_E, yz_N, yz_T, yz_NE, yz_SE, yz_BE, yz_TE, yz_TN, yz_BN, yz_TNE, yz_TNW, yz_TSE, yz_TSW/*, yz_ZERO*/;

   LBMReal kxyAverage, kyzAverage, kxzAverage, kxxMyyAverage, kxxMzzAverage; 

//   LBMReal a,b,c;

   // bulk viscosity
   LBMReal shearViscosity;
   LBMReal bulkViscosity;
   LBMReal OxxPyyPzzC;
   LBMReal OxxPyyPzzF;

   void setOffsets(LBMReal xoff, LBMReal yoff, LBMReal zoff) override;
   void calcMoments(const LBMReal* const f, LBMReal omega, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3, 
      LBMReal& kxy, LBMReal& kyz, LBMReal& kxz, LBMReal& kxxMyy, LBMReal& kxxMzz);
   void calcInterpolatedCoefficiets(const D3Q27ICell& icell, LBMReal omega, LBMReal eps_new) override;
   void calcInterpolatedNodeCF(LBMReal* f, LBMReal omega, LBMReal x, LBMReal y, LBMReal z, LBMReal press, LBMReal xs, LBMReal ys, LBMReal zs);
   LBMReal calcPressBSW();
   LBMReal calcPressTSW();
   LBMReal calcPressTSE();
   LBMReal calcPressBSE();
   LBMReal calcPressBNW();
   LBMReal calcPressTNW();
   LBMReal calcPressTNE();
   LBMReal calcPressBNE();
   void calcInterpolatedNodeFC(LBMReal* f, LBMReal omega) override;
   void calcInterpolatedVelocity(LBMReal x, LBMReal y, LBMReal z,LBMReal& vx1, LBMReal& vx2, LBMReal& vx3) override;
   void calcInterpolatedShearStress(LBMReal x, LBMReal y, LBMReal z,LBMReal& tauxx, LBMReal& tauyy, LBMReal& tauzz,LBMReal& tauxy, LBMReal& tauxz, LBMReal& tauyz) override;
};

//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF)
{
   this->interpolateCoarseToFine(icellC, icellF, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC)
{
   this->interpolateFineToCoarse(icellF, icellC, 0.0, 0.0, 0.0);
}

#endif
