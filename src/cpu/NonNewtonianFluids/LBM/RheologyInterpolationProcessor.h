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
//! \file RheologyInterpolationProcessor.h
//! \ingroup LBM
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RheologyInterpolationProcessor_H_
#define RheologyInterpolationProcessor_H_

#include "Interpolator.h"
#include "D3Q27System.h"

//! \brief A class implements an interpolation function of grid refinement for thixotropic fluid.

class RheologyInterpolationProcessor : public Interpolator
{
public:
   RheologyInterpolationProcessor();
   RheologyInterpolationProcessor(real omegaC, real omegaF, real omegaMin);
   virtual ~RheologyInterpolationProcessor();
   InterpolationProcessorPtr clone();
   void setOmegas(real omegaC, real omegaF);
   void setOmegaMin(real omegaMin);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF);
   void interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff);
   void interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC); 
   void interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff); 
   //real forcingC, forcingF;
protected:   
private:
   real omegaC, omegaF;
   real a0, ax, ay, az, axx, ayy, azz, axy, axz, ayz, b0, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz, c0, cx, cy, cz, cxx, cyy, czz, cxy, cxz, cyz, axyz, bxyz, cxyz;
   real xoff,    yoff,    zoff;
   real xoff_sq, yoff_sq, zoff_sq;
   real press_SWT, press_NWT, press_NET, press_SET, press_SWB, press_NWB, press_NEB, press_SEB;

   real  f_E,  f_N,  f_T,  f_NE,  f_SE,  f_BE,  f_TE,  f_TN,  f_BN,  f_TNE,  f_TNW,  f_TSE,  f_TSW,  f_ZERO;
   real  x_E,  x_N,  x_T,  x_NE,  x_SE,  x_BE,  x_TE,  x_TN,  x_BN,  x_TNE,  x_TNW,  x_TSE,  x_TSW,  x_ZERO;
   real  y_E,  y_N,  y_T,  y_NE,  y_SE,  y_BE,  y_TE,  y_TN,  y_BN,  y_TNE,  y_TNW,  y_TSE,  y_TSW,  y_ZERO;
   real  z_E,  z_N,  z_T,  z_NE,  z_SE,  z_BE,  z_TE,  z_TN,  z_BN,  z_TNE,  z_TNW,  z_TSE,  z_TSW,  z_ZERO;
   real xy_E, xy_N, xy_T, xy_NE, xy_SE, xy_BE, xy_TE, xy_TN, xy_BN, xy_TNE, xy_TNW, xy_TSE, xy_TSW/*, xy_ZERO*/;
   real xz_E, xz_N, xz_T, xz_NE, xz_SE, xz_BE, xz_TE, xz_TN, xz_BN, xz_TNE, xz_TNW, xz_TSE, xz_TSW/*, xz_ZERO*/;
   real yz_E, yz_N, yz_T, yz_NE, yz_SE, yz_BE, yz_TE, yz_TN, yz_BN, yz_TNE, yz_TNW, yz_TSE, yz_TSW/*, yz_ZERO*/;

   real kxyAverage, kyzAverage, kxzAverage, kxxMyyAverage, kxxMzzAverage; 

   real rho;
   real shearRate;

   real omegaMin;

   void setOffsets(real xoff, real yoff, real zoff);
   void calcMoments(const real* const f, real omegaInf, real& rho, real& vx1, real& vx2, real& vx3,
      real& kxy, real& kyz, real& kxz, real& kxxMyy, real& kxxMzz);
   void calcInterpolatedCoefficiets_intern(const D3Q27ICell& icell, real omega, real eps_new, real x, real y, real z, real xs, real ys, real zs);
   void calcInterpolatedNode(real* f, /*real omega,*/ real x, real y, real z, real press, real xs, real ys, real zs);
   real calcPressBSW();
   real calcPressTSW();
   real calcPressTSE();
   real calcPressBSE();
   real calcPressBNW();
   real calcPressTNW();
   real calcPressTNE();
   real calcPressBNE();
   void calcInterpolatedNodeFC(real* f, real omega);
   void calcInterpolatedVelocity(real x, real y, real z,real& vx1, real& vx2, real& vx3);
   void calcInterpolatedShearStress(real x, real y, real z,real& tauxx, real& tauyy, real& tauzz,real& tauxy, real& tauxz, real& tauyz);
};

//////////////////////////////////////////////////////////////////////////
inline void RheologyInterpolationProcessor::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF)
{
   this->interpolateCoarseToFine(icellC, icellF, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
inline void RheologyInterpolationProcessor::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC)
{
   this->interpolateFineToCoarse(icellF, icellC, 0.0, 0.0, 0.0);
}

#endif
