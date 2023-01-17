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
//! \file StressBcs27.cu
//! \author Henrik Asmuth
//! \date 16/05/2022
//! \brief Kernels for StressBC using the iMEM approach
//!
//! Both kernels prescribe a wall shear stress using the iMEM apprach (see, Asmuth et. al (2021), https://doi.org/10.1063/5.0065701)
//! QStressDeviceComp27 couples the iMEM to the single-node interpolated bounce-back.
//! BBStressDevice27 couples the iMEM to a simple bounce-back.
//! Note, that the iMEM function is currently only implemented for straight walls with z-normal and q=0.5.
//! Other wall models could be implemented in the iMEM by replacing the formulations from Monin-Obukhov similarity theory (MOST)
//! with other formulations, e.g., for smooth walls.
//! iMEM so far most extensively tested with BBStressDevice27, but QStressDeviceComp27 also seems to be stable and working.
//=======================================================================================

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "Kernel/Utilities/DistributionHelper.cuh"
#include <lbm/constants/NumericConstants.h>
#include "KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

//////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__ void iMEM(uint k, uint kN,
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


//////////////////////////////////////////////////////////////////////////////
__global__ void QStressDeviceComp27(real* DD,
											   int* k_Q,
                                    int* k_N,
											   real* QQ,
                                    unsigned int numberOfBCnodes,
                                    real om1,
                                    real* turbViscosity,
                                    real* vx,
                                    real* vy,
                                    real* vz,
                                    real* normalX,
                                    real* normalY,
                                    real* normalZ,
                                    real* vx_el,
                                    real* vy_el,
                                    real* vz_el,
                                    real* vx_w_mean,
                                    real* vy_w_mean,
                                    real* vz_w_mean,
                                    int* samplingOffset,
                                    real* z0,
                                    bool  hasWallModelMonitor,
                                    real* u_star_monitor,
                                    real* Fx_monitor,
                                    real* Fy_monitor,
                                    real* Fz_monitor,
											   unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned long long numberOfLBnodes,
                                    bool isEvenTimestep)
{

   Distributions27 D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index
   const unsigned  y = blockIdx.x;   // Globaler y-Index
   const unsigned  z = blockIdx.y;   // Globaler z-Index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k< numberOfBCnodes/*numberOfBCnodes*/)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB,
            *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
            *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
            *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
            *q_dirBSE, *q_dirBNW;
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;      //get right adress of post-coll f's
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[DIR_P00])[ke   ];     //post-coll f's
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho);


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho);

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho);

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

      real om_turb = om1 / (c1o1 + c3o1*om1*max(c0o1, turbViscosity[k_Q[k]]));
      //////////////////////////////////////////////////////////////////////////

      D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, !isEvenTimestep);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Compute incoming f's with zero wall velocity
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // incoming f's from bounce back
      real f_E_in = 0.0,  f_W_in = 0.0,  f_N_in = 0.0,  f_S_in = 0.0,  f_T_in = 0.0,  f_B_in = 0.0,   f_NE_in = 0.0,  f_SW_in = 0.0,  f_SE_in = 0.0,  f_NW_in = 0.0,  f_TE_in = 0.0,  f_BW_in = 0.0,  f_BE_in = 0.0, f_TW_in = 0.0, f_TN_in = 0.0, f_BS_in = 0.0, f_BN_in = 0.0, f_TS_in = 0.0, f_TNE_in = 0.0, f_TSW_in = 0.0, f_TSE_in = 0.0, f_TNW_in = 0.0, f_BNE_in = 0.0, f_BSW_in = 0.0, f_BSE_in = 0.0, f_BNW_in = 0.0;
      // momentum exchanged with wall at rest
      real wallMomentumX = 0.0, wallMomentumY = 0.0, wallMomentumZ = 0.0;
      real velocityLB = 0.0;
      
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_W_in = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, om_turb);
         wallMomentumX += f_E+f_W_in;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_E_in = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, om_turb);
         wallMomentumX -= f_W+f_E_in;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_S_in = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, om_turb);
         wallMomentumY += f_N+f_S_in;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_N_in = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, om_turb);
         wallMomentumY -= f_S+f_N_in;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_B_in = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, om_turb);
         wallMomentumZ += f_T+f_B_in;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         f_T_in = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, om_turb);
         wallMomentumZ -= f_B+f_T_in;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_SW_in = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, om_turb);
         wallMomentumX += f_NE+f_SW_in;
         wallMomentumY += f_NE+f_SW_in;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_NE_in = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, om_turb);
         wallMomentumX -= f_SW+f_NE_in;
         wallMomentumY -= f_SW+f_NE_in;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_NW_in = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, om_turb);
         wallMomentumX += f_SE+f_NW_in;
         wallMomentumY -= f_SE+f_NW_in;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_SE_in = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, om_turb);
         wallMomentumX -= f_NW+f_SE_in;
         wallMomentumY += f_NW+f_SE_in;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_BW_in = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, om_turb);
         wallMomentumX += f_TE+f_BW_in;
         wallMomentumZ += f_TE+f_BW_in;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_TE_in = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, om_turb);
         wallMomentumX -= f_BW+f_TE_in;
         wallMomentumZ -= f_BW+f_TE_in;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_TW_in = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, om_turb);
         wallMomentumX += f_BE+f_TW_in;
         wallMomentumZ -= f_BE+f_TW_in;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_BE_in = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, om_turb);
         wallMomentumX -= f_TW+f_BE_in;
         wallMomentumZ += f_TW+f_BE_in;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_BS_in = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, om_turb);
         wallMomentumY += f_TN+f_BS_in;
         wallMomentumZ += f_TN+f_BS_in;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_TN_in = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, om_turb);
         wallMomentumY -= f_BS+f_TN_in;
         wallMomentumZ -= f_BS+f_TN_in;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_TS_in = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, om_turb);
         wallMomentumY += f_BN+f_TS_in;
         wallMomentumZ -= f_BN+f_TS_in;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         f_BN_in = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, om_turb);
         wallMomentumY -= f_TS+f_BN_in;
         wallMomentumZ += f_TS+f_BN_in;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_BSW_in = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, om_turb);
         wallMomentumX += f_TNE+f_BSW_in;
         wallMomentumY += f_TNE+f_BSW_in;
         wallMomentumZ += f_TNE+f_BSW_in;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_TNE_in = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, om_turb);
         wallMomentumX -= f_BSW+f_TNE_in;
         wallMomentumY -= f_BSW+f_TNE_in;
         wallMomentumZ -= f_BSW+f_TNE_in;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_TSW_in = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, om_turb);
         wallMomentumX += f_BNE+f_TSW_in;
         wallMomentumY += f_BNE+f_TSW_in;
         wallMomentumZ -= f_BNE+f_TSW_in;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_BNE_in = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, om_turb);
         wallMomentumX -= f_TSW+f_BNE_in;
         wallMomentumY -= f_TSW+f_BNE_in;
         wallMomentumZ += f_TSW+f_BNE_in;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_BNW_in = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, om_turb);
         wallMomentumX += f_TSE+f_BNW_in;
         wallMomentumY -= f_TSE+f_BNW_in;
         wallMomentumZ += f_TSE+f_BNW_in;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_TSE_in = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, om_turb);
         wallMomentumX -= f_BNW+f_TSE_in;
         wallMomentumY += f_BNW+f_TSE_in;
         wallMomentumZ -= f_BNW+f_TSE_in;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_TNW_in = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, om_turb);
         wallMomentumX += f_BSE+f_TNW_in;
         wallMomentumY -= f_BSE+f_TNW_in;
         wallMomentumZ -= f_BSE+f_TNW_in;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         f_BSE_in = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, om_turb);
         wallMomentumX -= f_TNW+f_BSE_in;
         wallMomentumY += f_TNW+f_BSE_in;
         wallMomentumZ += f_TNW+f_BSE_in;
      }

      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Compute wall velocity
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real VeloX=0.0, VeloY=0.0, VeloZ=0.0;

      q = 0.5f;
      real eps = 0.001f;

      iMEM( k, k_N[k],
            normalX, normalY, normalZ,
            vx, vy, vz,
            vx_el,      vy_el,      vz_el,
            vx_w_mean,  vy_w_mean,  vz_w_mean,
            vx1,        vx2,        vx3,
            c1o1+drho,
            samplingOffset,
            q,
            1.0+q,
            eps,
            z0,
            hasWallModelMonitor,
            u_star_monitor,
            wallMomentumX, wallMomentumY, wallMomentumZ,
            VeloX, VeloY, VeloZ);

      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Add wall velocity and write f's
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M00])[kw] = f_W_in - (c6o1*c2o27*( VeloX     ))/(c1o1+q);
         wallMomentumX += -(c6o1*c2o27*( VeloX     ))/(c1o1+q);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P00])[ke] = f_E_in - (c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c2o27*(-VeloX     ))/(c1o1+q);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0M0])[ks] = f_S_in - (c6o1*c2o27*( VeloY     ))/(c1o1+q);
         wallMomentumY += - (c6o1*c2o27*( VeloY     ))/(c1o1+q);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0P0])[kn] = f_N_in - (c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         wallMomentumY -=  -(c6o1*c2o27*(-VeloY     ))/(c1o1+q);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00M])[kb] = f_B_in - (c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c2o27*( VeloZ     ))/(c1o1+q);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00P])[kt] = f_T_in - (c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         wallMomentumZ -= -(c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MM0])[ksw] = f_SW_in - (c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         wallMomentumX +=  -(c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         wallMomentumY +=  -(c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PP0])[kne] = f_NE_in - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MP0])[knw] = f_NW_in - (c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         wallMomentumX += -(c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         wallMomentumY -= -(c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PM0])[kse] = f_SE_in - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0M])[kbw] = f_BW_in - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0P])[kte] = f_TE_in - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0P])[ktw] = f_TW_in - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0M])[kbe] = f_BE_in - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MM])[kbs] = f_BS_in - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PP])[ktn] = f_TN_in - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MP])[kts] = f_TS_in - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PM])[kbn] = f_BN_in - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMM])[kbsw] = f_BSW_in - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPP])[ktne] = f_TNE_in - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMP])[ktsw] = f_TSW_in - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPM])[kbne] = f_BNE_in - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPM])[kbnw] = f_BNW_in - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMP])[ktse] = f_TSE_in - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPP])[ktnw] = f_TNW_in - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMM])[kbse] = f_BSE_in - (c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      if(hasWallModelMonitor)
      {
         Fx_monitor[k] = wallMomentumX;
         Fy_monitor[k] = wallMomentumY;
         Fz_monitor[k] = wallMomentumZ;
      }

   }
}

//////////////////////////////////////////////////////////////////////////////
__global__ void BBStressDevice27( real* DD,
											            int* k_Q,
                                             int* k_N,
                                             real* QQ,
                                             unsigned int  numberOfBCnodes,
                                             real* vx,
                                             real* vy,
                                             real* vz,
                                             real* normalX,
                                             real* normalY,
                                             real* normalZ,
                                             real* vx_el,
                                             real* vy_el,
                                             real* vz_el,
                                             real* vx_w_mean,
                                             real* vy_w_mean,
                                             real* vz_w_mean,
                                             int* samplingOffset,
                                             real* z0,
                                             bool  hasWallModelMonitor,
                                             real* u_star_monitor,
                                             real* Fx_monitor,
                                             real* Fy_monitor,
                                             real* Fz_monitor,
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned long long numberOfLBnodes,
                                             bool isEvenTimestep)
{

   Distributions27 D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index
   const unsigned  y = blockIdx.x;   // Globaler y-Index
   const unsigned  z = blockIdx.y;   // Globaler z-Index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k< numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB,
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW;
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho);


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho);

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////

      D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, !isEvenTimestep);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E_in,  f_W_in,  f_N_in,  f_S_in,  f_T_in,  f_B_in,   f_NE_in,  f_SW_in,  f_SE_in,  f_NW_in,  f_TE_in,  f_BW_in,  f_BE_in,
         f_TW_in, f_TN_in, f_BS_in, f_BN_in, f_TS_in, f_TNE_in, f_TSW_in, f_TSE_in, f_TNW_in, f_BNE_in, f_BSW_in, f_BSE_in, f_BNW_in;

      // momentum exchanged with wall at rest
      real wallMomentumX = 0.0, wallMomentumY = 0.0, wallMomentumZ = 0.0;

      real q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_W_in=f_E;
         wallMomentumX += f_E+f_W_in;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_E_in=f_W;
          wallMomentumX -= f_W+f_E_in;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_S_in=f_N;
         wallMomentumY += f_N+f_S_in;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_N_in=f_S;
         wallMomentumY -= f_S+f_N_in;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_B_in=f_T;
         wallMomentumZ += f_T+f_B_in;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_T_in=f_B;
         wallMomentumZ -= f_B+f_T_in;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_SW_in=f_NE;
         wallMomentumX += f_NE+f_SW_in;
         wallMomentumY += f_NE+f_SW_in;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_NE_in=f_SW;
         wallMomentumX -= f_SW+f_NE_in;
         wallMomentumY -= f_SW+f_NE_in;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_NW_in=f_SE;
         wallMomentumX += f_SE+f_NW_in;
         wallMomentumY -= f_SE+f_NW_in;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_SE_in=f_NW;
         wallMomentumX -= f_NW+f_SE_in;
         wallMomentumY += f_NW+f_SE_in;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BW_in=f_TE;
         wallMomentumX += f_TE+f_BW_in;
         wallMomentumZ += f_TE+f_BW_in;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TE_in=f_BW;
         wallMomentumX -= f_BW+f_TE_in;
         wallMomentumZ -= f_BW+f_TE_in;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TW_in=f_BE;
         wallMomentumX += f_BE+f_TW_in;
         wallMomentumZ -= f_BE+f_TW_in;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BE_in=f_TW;
         wallMomentumX -= f_TW+f_BE_in;
         wallMomentumZ += f_TW+f_BE_in;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BS_in=f_TN;
         wallMomentumY += f_TN+f_BS_in;
         wallMomentumZ += f_TN+f_BS_in;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TN_in=f_BS;
         wallMomentumY -= f_BS+f_TN_in;
         wallMomentumZ -= f_BS+f_TN_in;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TS_in=f_BN;
         wallMomentumY += f_BN+f_TS_in;
         wallMomentumZ -= f_BN+f_TS_in;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BN_in=f_TS;
         wallMomentumY -= f_TS+f_BN_in;
         wallMomentumZ += f_TS+f_BN_in;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BSW_in=f_TNE;
         wallMomentumX += f_TNE+f_BSW_in;
         wallMomentumY += f_TNE+f_BSW_in;
         wallMomentumZ += f_TNE+f_BSW_in;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TNE_in=f_BSW;
         wallMomentumX -= f_BSW+f_TNE_in;
         wallMomentumY -= f_BSW+f_TNE_in;
         wallMomentumZ -= f_BSW+f_TNE_in;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TSW_in=f_BNE;
         wallMomentumX += f_BNE+f_TSW_in;
         wallMomentumY += f_BNE+f_TSW_in;
         wallMomentumZ -= f_BNE+f_TSW_in;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BNE_in=f_TSW;
         wallMomentumX -= f_TSW+f_BNE_in;
         wallMomentumY -= f_TSW+f_BNE_in;
         wallMomentumZ += f_TSW+f_BNE_in;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BNW_in=f_TSE;
         wallMomentumX += f_TSE+f_BNW_in;
         wallMomentumY -= f_TSE+f_BNW_in;
         wallMomentumZ += f_TSE+f_BNW_in;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TSE_in=f_BNW;
         wallMomentumX -= f_BNW+f_TSE_in;
         wallMomentumY += f_BNW+f_TSE_in;
         wallMomentumZ -= f_BNW+f_TSE_in;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TNW_in=f_BSE;
         wallMomentumX += f_BSE+f_TNW_in;
         wallMomentumY -= f_BSE+f_TNW_in;
         wallMomentumZ -= f_BSE+f_TNW_in;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BSE_in=f_TNW;
         wallMomentumX -= f_TNW+f_BSE_in;
         wallMomentumY += f_TNW+f_BSE_in;
         wallMomentumZ += f_TNW+f_BSE_in;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Compute wall velocity
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real VeloX=0.0, VeloY=0.0, VeloZ=0.0;

      q = q_dirB[k];
      real eps = 0.001f;

      iMEM( k, k_N[k],
         normalX, normalY, normalZ,
         vx, vy, vz,
         vx_el,      vy_el,      vz_el,
         vx_w_mean,  vy_w_mean,  vz_w_mean,
         vx1,        vx2,        vx3,
         c1o1+drho,
         samplingOffset,
         q,
         1.0,
         eps,
         z0,
         hasWallModelMonitor,
         u_star_monitor,
         wallMomentumX, wallMomentumY, wallMomentumZ,
         VeloX, VeloY, VeloZ);

      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Add wall velocity and write f's
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M00])[kw] = f_W_in - (c6o1*c2o27*( VeloX     ));
         wallMomentumX += -(c6o1*c2o27*( VeloX     ));
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P00])[ke] = f_E_in - (c6o1*c2o27*(-VeloX     ));
         wallMomentumX -= - (c6o1*c2o27*(-VeloX     ));
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0M0])[ks] = f_S_in - (c6o1*c2o27*( VeloY     ));
         wallMomentumY += - (c6o1*c2o27*( VeloY     ));
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0P0])[kn] = f_N_in - (c6o1*c2o27*(-VeloY     ));
         wallMomentumY -=  -(c6o1*c2o27*(-VeloY     ));
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00M])[kb] = f_B_in - (c6o1*c2o27*( VeloZ     ));
         wallMomentumZ += - (c6o1*c2o27*( VeloZ     ));
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00P])[kt] = f_T_in - (c6o1*c2o27*(-VeloZ     ));
         wallMomentumZ -= -(c6o1*c2o27*(-VeloZ     ));
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MM0])[ksw] = f_SW_in - (c6o1*c1o54*(VeloX+VeloY));
         wallMomentumX +=  -(c6o1*c1o54*(VeloX+VeloY));
         wallMomentumY +=  -(c6o1*c1o54*(VeloX+VeloY));
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PP0])[kne] = f_NE_in - (c6o1*c1o54*(-VeloX-VeloY));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloY));
         wallMomentumY -= - (c6o1*c1o54*(-VeloX-VeloY));
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MP0])[knw] = f_NW_in - (c6o1*c1o54*( VeloX-VeloY));
         wallMomentumX += -(c6o1*c1o54*( VeloX-VeloY));
         wallMomentumY -= -(c6o1*c1o54*( VeloX-VeloY));
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PM0])[kse] = f_SE_in - (c6o1*c1o54*(-VeloX+VeloY));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloY));
         wallMomentumY += - (c6o1*c1o54*(-VeloX+VeloY));
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0M])[kbw] = f_BW_in - (c6o1*c1o54*( VeloX+VeloZ));
         wallMomentumX += - (c6o1*c1o54*( VeloX+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( VeloX+VeloZ));
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0P])[kte] = f_TE_in - (c6o1*c1o54*(-VeloX-VeloZ));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*(-VeloX-VeloZ));
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0P])[ktw] = f_TW_in - (c6o1*c1o54*( VeloX-VeloZ));
         wallMomentumX += - (c6o1*c1o54*( VeloX-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( VeloX-VeloZ));
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0M])[kbe] = f_BE_in - (c6o1*c1o54*(-VeloX+VeloZ));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*(-VeloX+VeloZ));
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MM])[kbs] = f_BS_in - (c6o1*c1o54*( VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o54*( VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( VeloY+VeloZ));
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PP])[ktn] = f_TN_in - (c6o1*c1o54*( -VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o54*( -VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( -VeloY-VeloZ));
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MP])[kts] = f_TS_in - (c6o1*c1o54*( VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o54*( VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( VeloY-VeloZ));
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PM])[kbn] = f_BN_in - (c6o1*c1o54*( -VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o54*( -VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( -VeloY+VeloZ));
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMM])[kbsw] = f_BSW_in - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPP])[ktne] = f_TNE_in - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMP])[ktsw] = f_TSW_in - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPM])[kbne] = f_BNE_in - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPM])[kbnw] = f_BNW_in - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMP])[ktse] = f_TSE_in - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPP])[ktnw] = f_TNW_in - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMM])[kbse] = f_BSE_in - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
      }

      if(hasWallModelMonitor)
      {
         Fx_monitor[k] = wallMomentumX;
         Fy_monitor[k] = wallMomentumY;
         Fz_monitor[k] = wallMomentumZ;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
__global__ void BBStressPressureDevice27( real* DD,
											            int* k_Q,
                                             int* k_N,
                                             real* QQ,
                                             unsigned int  numberOfBCnodes,
                                             real* vx,
                                             real* vy,
                                             real* vz,
                                             real* normalX,
                                             real* normalY,
                                             real* normalZ,
                                             real* vx_el,
                                             real* vy_el,
                                             real* vz_el,
                                             real* vx_w_mean,
                                             real* vy_w_mean,
                                             real* vz_w_mean,
                                             int* samplingOffset,
                                             real* z0,
                                             bool  hasWallModelMonitor,
                                             real* u_star_monitor,
                                             real* Fx_monitor,
                                             real* Fy_monitor,
                                             real* Fz_monitor,
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned long long numberOfLBnodes,
                                             bool isEvenTimestep)
{
   Distributions27 D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);

   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index
   const unsigned  y = blockIdx.x;   // Globaler y-Index
   const unsigned  z = blockIdx.y;   // Globaler z-Index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k< numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB,
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW;
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];

      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho);


      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho);

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho);

      //////////////////////////////////////////////////////////////////////////
      D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E_in,  f_W_in,  f_N_in,  f_S_in,  f_T_in,  f_B_in,   f_NE_in,  f_SW_in,  f_SE_in,  f_NW_in,  f_TE_in,  f_BW_in,  f_BE_in,
         f_TW_in, f_TN_in, f_BS_in, f_BN_in, f_TS_in, f_TNE_in, f_TSW_in, f_TSE_in, f_TNW_in, f_BNE_in, f_BSW_in, f_BSE_in, f_BNW_in;

      // momentum exchanged with wall at rest
      real wallMomentumX = 0.0, wallMomentumY = 0.0, wallMomentumZ = 0.0;

      real q;
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_W_in=f_E - c2o27 * drho;
         wallMomentumX += f_E+f_W_in;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_E_in=f_W - c2o27 * drho;
          wallMomentumX -= f_W+f_E_in;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_S_in=f_N - c2o27 * drho;
         wallMomentumY += f_N+f_S_in;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_N_in=f_S - c2o27 * drho;
         wallMomentumY -= f_S+f_N_in;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_B_in=f_T - c2o27 * drho;
         wallMomentumZ += f_T+f_B_in;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_T_in=f_B - c2o27 * drho;
         wallMomentumZ -= f_B+f_T_in;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_SW_in=f_NE - c1o54 * drho;
         wallMomentumX += f_NE+f_SW_in;
         wallMomentumY += f_NE+f_SW_in;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_NE_in=f_SW - c1o54 * drho;
         wallMomentumX -= f_SW+f_NE_in;
         wallMomentumY -= f_SW+f_NE_in;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_NW_in=f_SE - c1o54 * drho;
         wallMomentumX += f_SE+f_NW_in;
         wallMomentumY -= f_SE+f_NW_in;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_SE_in=f_NW - c1o54 * drho;
         wallMomentumX -= f_NW+f_SE_in;
         wallMomentumY += f_NW+f_SE_in;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BW_in=f_TE - c1o54 * drho;
         wallMomentumX += f_TE+f_BW_in;
         wallMomentumZ += f_TE+f_BW_in;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TE_in=f_BW - c1o54 * drho;
         wallMomentumX -= f_BW+f_TE_in;
         wallMomentumZ -= f_BW+f_TE_in;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TW_in=f_BE - c1o54 * drho;
         wallMomentumX += f_BE+f_TW_in;
         wallMomentumZ -= f_BE+f_TW_in;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BE_in=f_TW - c1o54 * drho;
         wallMomentumX -= f_TW+f_BE_in;
         wallMomentumZ += f_TW+f_BE_in;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BS_in=f_TN - c1o54 * drho;
         wallMomentumY += f_TN+f_BS_in;
         wallMomentumZ += f_TN+f_BS_in;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TN_in=f_BS - c1o54 * drho;
         wallMomentumY -= f_BS+f_TN_in;
         wallMomentumZ -= f_BS+f_TN_in;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TS_in=f_BN - c1o54 * drho;
         wallMomentumY += f_BN+f_TS_in;
         wallMomentumZ -= f_BN+f_TS_in;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BN_in=f_TS - c1o54 * drho;
         wallMomentumY -= f_TS+f_BN_in;
         wallMomentumZ += f_TS+f_BN_in;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BSW_in=f_TNE - c1o216 * drho;
         wallMomentumX += f_TNE+f_BSW_in;
         wallMomentumY += f_TNE+f_BSW_in;
         wallMomentumZ += f_TNE+f_BSW_in;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TNE_in=f_BSW - c1o216 * drho;
         wallMomentumX -= f_BSW+f_TNE_in;
         wallMomentumY -= f_BSW+f_TNE_in;
         wallMomentumZ -= f_BSW+f_TNE_in;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TSW_in=f_BNE - c1o216 * drho;
         wallMomentumX += f_BNE+f_TSW_in;
         wallMomentumY += f_BNE+f_TSW_in;
         wallMomentumZ -= f_BNE+f_TSW_in;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BNE_in=f_TSW - c1o216 * drho;
         wallMomentumX -= f_TSW+f_BNE_in;
         wallMomentumY -= f_TSW+f_BNE_in;
         wallMomentumZ += f_TSW+f_BNE_in;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BNW_in=f_TSE - c1o216 * drho;
         wallMomentumX += f_TSE+f_BNW_in;
         wallMomentumY -= f_TSE+f_BNW_in;
         wallMomentumZ += f_TSE+f_BNW_in;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TSE_in=f_BNW - c1o216 * drho;
         wallMomentumX -= f_BNW+f_TSE_in;
         wallMomentumY += f_BNW+f_TSE_in;
         wallMomentumZ -= f_BNW+f_TSE_in;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_TNW_in=f_BSE - c1o216 * drho;
         wallMomentumX += f_BSE+f_TNW_in;
         wallMomentumY -= f_BSE+f_TNW_in;
         wallMomentumZ -= f_BSE+f_TNW_in;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         f_BSE_in=f_TNW - c1o216 * drho;
         wallMomentumX -= f_TNW+f_BSE_in;
         wallMomentumY += f_TNW+f_BSE_in;
         wallMomentumZ += f_TNW+f_BSE_in;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Compute wall velocity
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real VeloX=0.0, VeloY=0.0, VeloZ=0.0;

      q = q_dirB[k];
      real eps = 0.001f;

      iMEM( k, k_N[k],
         normalX, normalY, normalZ,
         vx, vy, vz,
         vx_el,      vy_el,      vz_el,
         vx_w_mean,  vy_w_mean,  vz_w_mean,
         vx1,        vx2,        vx3,
         c1o1+drho,
         samplingOffset,
         q,
         1.0,
         eps,
         z0,
         hasWallModelMonitor,
         u_star_monitor,
         wallMomentumX, wallMomentumY, wallMomentumZ,
         VeloX, VeloY, VeloZ);

      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // //Add wall velocity and write f's
      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M00])[kw] = f_W_in - (c6o1*c2o27*( VeloX     ));
         wallMomentumX += -(c6o1*c2o27*( VeloX     ));
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P00])[ke] = f_E_in - (c6o1*c2o27*(-VeloX     ));
         wallMomentumX -= - (c6o1*c2o27*(-VeloX     ));
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0M0])[ks] = f_S_in - (c6o1*c2o27*( VeloY     ));
         wallMomentumY += - (c6o1*c2o27*( VeloY     ));
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0P0])[kn] = f_N_in - (c6o1*c2o27*(-VeloY     ));
         wallMomentumY -=  -(c6o1*c2o27*(-VeloY     ));
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00M])[kb] = f_B_in - (c6o1*c2o27*( VeloZ     ));
         wallMomentumZ += - (c6o1*c2o27*( VeloZ     ));
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00P])[kt] = f_T_in - (c6o1*c2o27*(-VeloZ     ));
         wallMomentumZ -= -(c6o1*c2o27*(-VeloZ     ));
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MM0])[ksw] = f_SW_in - (c6o1*c1o54*(VeloX+VeloY));
         wallMomentumX +=  -(c6o1*c1o54*(VeloX+VeloY));
         wallMomentumY +=  -(c6o1*c1o54*(VeloX+VeloY));
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PP0])[kne] = f_NE_in - (c6o1*c1o54*(-VeloX-VeloY));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloY));
         wallMomentumY -= - (c6o1*c1o54*(-VeloX-VeloY));
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MP0])[knw] = f_NW_in - (c6o1*c1o54*( VeloX-VeloY));
         wallMomentumX += -(c6o1*c1o54*( VeloX-VeloY));
         wallMomentumY -= -(c6o1*c1o54*( VeloX-VeloY));
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PM0])[kse] = f_SE_in - (c6o1*c1o54*(-VeloX+VeloY));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloY));
         wallMomentumY += - (c6o1*c1o54*(-VeloX+VeloY));
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0M])[kbw] = f_BW_in - (c6o1*c1o54*( VeloX+VeloZ));
         wallMomentumX += - (c6o1*c1o54*( VeloX+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( VeloX+VeloZ));
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0P])[kte] = f_TE_in - (c6o1*c1o54*(-VeloX-VeloZ));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*(-VeloX-VeloZ));
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0P])[ktw] = f_TW_in - (c6o1*c1o54*( VeloX-VeloZ));
         wallMomentumX += - (c6o1*c1o54*( VeloX-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( VeloX-VeloZ));
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0M])[kbe] = f_BE_in - (c6o1*c1o54*(-VeloX+VeloZ));
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*(-VeloX+VeloZ));
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MM])[kbs] = f_BS_in - (c6o1*c1o54*( VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o54*( VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( VeloY+VeloZ));
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PP])[ktn] = f_TN_in - (c6o1*c1o54*( -VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o54*( -VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( -VeloY-VeloZ));
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MP])[kts] = f_TS_in - (c6o1*c1o54*( VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o54*( VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o54*( VeloY-VeloZ));
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PM])[kbn] = f_BN_in - (c6o1*c1o54*( -VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o54*( -VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o54*( -VeloY+VeloZ));
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMM])[kbsw] = f_BSW_in - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*( VeloX+VeloY+VeloZ));
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPP])[ktne] = f_TNE_in - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ));
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMP])[ktsw] = f_TSW_in - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*( VeloX+VeloY-VeloZ));
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPM])[kbne] = f_BNE_in - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*(-VeloX-VeloY+VeloZ));
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPM])[kbnw] = f_BNW_in - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*( VeloX-VeloY+VeloZ));
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMP])[ktse] = f_TSE_in - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ));
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPP])[ktnw] = f_TNW_in - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
         wallMomentumZ -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ));
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMM])[kbse] = f_BSE_in - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
         wallMomentumZ += - (c6o1*c1o216*(-VeloX+VeloY+VeloZ));
      }

      if(hasWallModelMonitor)
      {
         Fx_monitor[k] = wallMomentumX;
         Fy_monitor[k] = wallMomentumY;
         Fz_monitor[k] = wallMomentumZ;
      }

   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////