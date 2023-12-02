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
//! \author Henrik Asmuth, Martin Schönherr
//! iMEM approach (see, Asmuth et. al (2021), https://doi.org/10.1063/5.0065701)
//! Note, that the iMEM function is currently only implemented for straight walls with z-normal and q=0.5.
//! Other wall models could be implemented in the iMEM by replacing the formulations from Monin-Obukhov similarity theory (MOST)
//! with other formulations, e.g., for smooth walls.
//! iMEM so far most extensively tested with StressBounceBackCompressible_Device, but StressCompressible_Device also seems to be stable and working.
//=======================================================================================
#ifndef iMEM_H
#define iMEM_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__ void iMEM(
    uint k, uint kN,
    real* _wallNormalX, real* _wallNormalY, real* _wallNormalZ,
    real* vx, real* vy, real* vz,
    real* vx_el,      real* vy_el,      real* vz_el,      //!>mean (temporally filtered) velocities at exchange location
    real* vx_w_mean,  real* vy_w_mean,  real* vz_w_mean,  //!>mean (temporally filtered) velocities at wall-adjactent node
    real  vx_w_inst,  real  vy_w_inst,  real  vz_w_inst,  //!>instantaneous velocities at wall-adjactent node
    real  rho,
    int* samplingOffset,
    real q,
    real forceFactor,                                     //!>e.g., 1.0 for simple-bounce back, or (1+q) for interpolated single-node bounce-back as in Geier et al (2015)
    real eps,                                             //!>filter constant in temporal averaging
    real* z0,                                             //!>aerodynamic roughness length
    bool  hasWallModelMonitor,
    real* u_star_monitor,
    real wallMomentumX, real wallMomentumY, real wallMomentumZ,
    real& wallVelocityX, real& wallVelocityY, real&wallVelocityZ)
{
      real wallNormalX = _wallNormalX[k];
      real wallNormalY = _wallNormalY[k];
      real wallNormalZ = _wallNormalZ[k];

      //Sample velocity at exchange location and filter temporally
      real _vx_el = eps*vx[kN]+(1.0-eps)*vx_el[k];
      real _vy_el = eps*vy[kN]+(1.0-eps)*vy_el[k];
      real _vz_el = eps*vz[kN]+(1.0-eps)*vz_el[k];
      vx_el[k] = _vx_el;
      vy_el[k] = _vy_el;
      vz_el[k] = _vz_el;

      //filter velocity at wall-adjacent node
      real _vx_w_mean = eps*vx_w_inst+(1.0-eps)*vx_w_mean[k];
      real _vy_w_mean = eps*vy_w_inst+(1.0-eps)*vy_w_mean[k];
      real _vz_w_mean = eps*vz_w_inst+(1.0-eps)*vz_w_mean[k];
      vx_w_mean[k] = _vx_w_mean;
      vy_w_mean[k] = _vy_w_mean;
      vz_w_mean[k] = _vz_w_mean;

      //Subtract wall-normal velocity components
      real vDotN_el = _vx_el*wallNormalX + _vy_el*wallNormalY + _vz_el*wallNormalZ;
      _vx_el -= vDotN_el*wallNormalX;
      _vy_el -= vDotN_el*wallNormalY;
      _vz_el -= vDotN_el*wallNormalZ;
      real vMag_el = sqrt( _vx_el*_vx_el + _vy_el*_vy_el + _vz_el*_vz_el );

      real vDotN_w_mean = _vx_w_mean*wallNormalX + _vy_w_mean*wallNormalY + _vz_w_mean*wallNormalZ;
      _vx_w_mean -= vDotN_w_mean*wallNormalX;
      _vy_w_mean -= vDotN_w_mean*wallNormalY;
      _vz_w_mean -= vDotN_w_mean*wallNormalZ;
      real vMag_w_mean = sqrt( _vx_w_mean*_vx_w_mean + _vy_w_mean*_vy_w_mean + _vz_w_mean*_vz_w_mean );

      real vDotN_w = vx_w_inst*wallNormalX + vy_w_inst*wallNormalY + vz_w_inst*wallNormalZ;
      real _vx_w = vx_w_inst-vDotN_w*wallNormalX;
      real _vy_w = vy_w_inst-vDotN_w*wallNormalY;
      real _vz_w = vz_w_inst-vDotN_w*wallNormalZ;

      //Compute wall shear stress tau_w via MOST
      real z = (real)samplingOffset[k] + q; //assuming q=0.5, could be replaced by wall distance via wall normal
      real kappa = 0.4;
      real u_star = vMag_el*kappa/(log(z/z0[k]));
      if(hasWallModelMonitor) u_star_monitor[k] = u_star;
      real tau_w = u_star*u_star;                  //Note: this is actually tau_w/rho
      real A = 1.0;                                //wall area (obviously 1 for grid aligned walls, can come from grid builder later for complex geometries)

      //Scale wall shear stress with near wall velocity, i.e., Schumann-Grötzbach (SG) approach
      real F_w_x = (tau_w*A) * (_vx_w/vMag_w_mean);//(_vx_el/vMag_el)
      real F_w_y = (tau_w*A) * (_vy_w/vMag_w_mean);//(_vy_el/vMag_el)
      real F_w_z = (tau_w*A) * (_vz_w/vMag_w_mean);//(_vz_el/vMag_el)
      //                                                ^^^^^^^^^^^^--- old alternative: do not scale SG-like but only set direction via velocity at exchange location

      //Momentum to be applied via wall velocity
      real wallMomDotN = wallMomentumX*wallNormalX+wallMomentumY*wallNormalY+wallMomentumZ*wallNormalZ;
      real F_x =  F_w_x - ( wallMomentumX - wallMomDotN*wallNormalX )/rho;
      real F_y =  F_w_y - ( wallMomentumY - wallMomDotN*wallNormalY )/rho;
      real F_z =  F_w_z - ( wallMomentumZ - wallMomDotN*wallNormalZ )/rho;

      //Compute  wall velocity and clip (clipping only necessary for initial boundary layer development)
      real clipWallVelo = 2.0;
      real clipVx = clipWallVelo*_vx_el;
      real clipVy = clipWallVelo*_vy_el;
      real clipVz = clipWallVelo*_vz_el;

      wallVelocityX = clipVx > -clipVx? min(clipVx, max(-clipVx, -3.0*F_x*forceFactor)): max(clipVx, min(-clipVx, -3.0*F_x*forceFactor));
      wallVelocityY = clipVy > -clipVy? min(clipVy, max(-clipVy, -3.0*F_y*forceFactor)): max(clipVy, min(-clipVy, -3.0*F_y*forceFactor));
      wallVelocityZ = clipVz > -clipVz? min(clipVz, max(-clipVz, -3.0*F_z*forceFactor)): max(clipVz, min(-clipVz, -3.0*F_z*forceFactor));
}

#endif
