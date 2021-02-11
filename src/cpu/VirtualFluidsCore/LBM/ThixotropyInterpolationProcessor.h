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
//! \file ThixotropyInterpolationProcessor.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef ThixotropyInterpolationProcessor_H_
#define ThixotropyInterpolationProcessor_H_

#include "InterpolationProcessor.h"
#include "D3Q27System.h"

//! \brief A class implements an interpolation function of grid refinement for thixotropic fluid.

class ThixotropyInterpolationProcessor : public InterpolationProcessor
{
public:
   ThixotropyInterpolationProcessor();
   ThixotropyInterpolationProcessor(LBMReal omegaC, LBMReal omegaF, LBMReal omegaMin);
   virtual ~ThixotropyInterpolationProcessor();
   InterpolationProcessorPtr clone();
   void setOmegas(LBMReal omegaC, LBMReal omegaF);
   void setOmegaMin(LBMReal omegaMin);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, LBMReal xoff, LBMReal yoff, LBMReal zoff);
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC); 
   void interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC, LBMReal xoff, LBMReal yoff, LBMReal zoff); 
   //LBMReal forcingC, forcingF;
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

   LBMReal rho;
   LBMReal shearRate;

   LBMReal omegaMin;

   void setOffsets(LBMReal xoff, LBMReal yoff, LBMReal zoff);
   void calcMoments(const LBMReal* const f, LBMReal omegaInf, LBMReal& rho, LBMReal& vx1, LBMReal& vx2, LBMReal& vx3,
      LBMReal& kxy, LBMReal& kyz, LBMReal& kxz, LBMReal& kxxMyy, LBMReal& kxxMzz);
   void calcInterpolatedCoefficiets_intern(const D3Q27ICell& icell, LBMReal omega, LBMReal eps_new, LBMReal x, LBMReal y, LBMReal z, LBMReal xs, LBMReal ys, LBMReal zs);
   void calcInterpolatedNode(LBMReal* f, /*LBMReal omega,*/ LBMReal x, LBMReal y, LBMReal z, LBMReal press, LBMReal xs, LBMReal ys, LBMReal zs);
   LBMReal calcPressBSW();
   LBMReal calcPressTSW();
   LBMReal calcPressTSE();
   LBMReal calcPressBSE();
   LBMReal calcPressBNW();
   LBMReal calcPressTNW();
   LBMReal calcPressTNE();
   LBMReal calcPressBNE();
   void calcInterpolatedNodeFC(LBMReal* f, LBMReal omega);
   void calcInterpolatedVelocity(LBMReal x, LBMReal y, LBMReal z,LBMReal& vx1, LBMReal& vx2, LBMReal& vx3);
   void calcInterpolatedShearStress(LBMReal x, LBMReal y, LBMReal z,LBMReal& tauxx, LBMReal& tauyy, LBMReal& tauzz,LBMReal& tauxy, LBMReal& tauxz, LBMReal& tauyz);
};

//////////////////////////////////////////////////////////////////////////
inline void ThixotropyInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF)
{
   this->interpolateCoarseToFine(icellC, icellF, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
inline void ThixotropyInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, LBMReal* icellC)
{
   this->interpolateFineToCoarse(icellF, icellC, 0.0, 0.0, 0.0);
}

#endif
