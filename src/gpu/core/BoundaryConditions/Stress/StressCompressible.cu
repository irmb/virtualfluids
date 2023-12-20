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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Henrik Asmuth, Martin Schönherr
//! \brief Kernel for StressBC using the iMEM approach
//!
//! kernel prescribe a wall shear stress using the iMEM apprach (see, Asmuth et. al (2021), https://doi.org/10.1063/5.0065701)
//! StressCompressible_Device couples the iMEM to the single-node interpolated bounce-back.
//! StressBounceBackCompressible_Device couples the iMEM to a simple bounce-back.
//! Note, that the iMEM function is currently only implemented for straight walls with z-normal and q=0.5.
//! Other wall models could be implemented in the iMEM by replacing the formulations from Monin-Obukhov similarity theory (MOST)
//! with other formulations, e.g., for smooth walls.
//! iMEM so far most extensively tested with StressBounceBackCompressible_Device, but StressCompressible_Device also seems to be stable and working.
//=======================================================================================

#include "BoundaryConditions/Stress/iMEM.cuh"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void StressCompressible_Device(
    real* DD,
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
      q_dirE   = &QQ[dP00 * numberOfBCnodes];
      q_dirW   = &QQ[dM00 * numberOfBCnodes];
      q_dirN   = &QQ[d0P0 * numberOfBCnodes];
      q_dirS   = &QQ[d0M0 * numberOfBCnodes];
      q_dirT   = &QQ[d00P * numberOfBCnodes];
      q_dirB   = &QQ[d00M * numberOfBCnodes];
      q_dirNE  = &QQ[dPP0 * numberOfBCnodes];
      q_dirSW  = &QQ[dMM0 * numberOfBCnodes];
      q_dirSE  = &QQ[dPM0 * numberOfBCnodes];
      q_dirNW  = &QQ[dMP0 * numberOfBCnodes];
      q_dirTE  = &QQ[dP0P * numberOfBCnodes];
      q_dirBW  = &QQ[dM0M * numberOfBCnodes];
      q_dirBE  = &QQ[dP0M * numberOfBCnodes];
      q_dirTW  = &QQ[dM0P * numberOfBCnodes];
      q_dirTN  = &QQ[d0PP * numberOfBCnodes];
      q_dirBS  = &QQ[d0MM * numberOfBCnodes];
      q_dirBN  = &QQ[d0PM * numberOfBCnodes];
      q_dirTS  = &QQ[d0MP * numberOfBCnodes];
      q_dirTNE = &QQ[dPPP * numberOfBCnodes];
      q_dirTSW = &QQ[dMMP * numberOfBCnodes];
      q_dirTSE = &QQ[dPMP * numberOfBCnodes];
      q_dirTNW = &QQ[dMPP * numberOfBCnodes];
      q_dirBNE = &QQ[dPPM * numberOfBCnodes];
      q_dirBSW = &QQ[dMMM * numberOfBCnodes];
      q_dirBSE = &QQ[dPMM * numberOfBCnodes];
      q_dirBNW = &QQ[dMPM * numberOfBCnodes];
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

      f_W    = (D.f[dP00])[ke   ];     //post-coll f's
      f_E    = (D.f[dM00])[kw   ];
      f_S    = (D.f[d0P0])[kn   ];
      f_N    = (D.f[d0M0])[ks   ];
      f_B    = (D.f[d00P])[kt   ];
      f_T    = (D.f[d00M])[kb   ];
      f_SW   = (D.f[dPP0])[kne  ];
      f_NE   = (D.f[dMM0])[ksw  ];
      f_NW   = (D.f[dPM0])[kse  ];
      f_SE   = (D.f[dMP0])[knw  ];
      f_BW   = (D.f[dP0P])[kte  ];
      f_TE   = (D.f[dM0M])[kbw  ];
      f_TW   = (D.f[dP0M])[kbe  ];
      f_BE   = (D.f[dM0P])[ktw  ];
      f_BS   = (D.f[d0PP])[ktn  ];
      f_TN   = (D.f[d0MM])[kbs  ];
      f_TS   = (D.f[d0PM])[kbn  ];
      f_BN   = (D.f[d0MP])[kts  ];
      f_BSW  = (D.f[dPPP])[ktne ];
      f_BNE  = (D.f[dMMP])[ktsw ];
      f_BNW  = (D.f[dPMP])[ktse ];
      f_BSE  = (D.f[dMPP])[ktnw ];
      f_TSW  = (D.f[dPPM])[kbne ];
      f_TNE  = (D.f[dMMM])[kbsw ];
      f_TNW  = (D.f[dPMM])[kbse ];
      f_TSE  = (D.f[dMPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);

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
         // f_W_in = getInterpolatedDistributionForNoSlipBC(q, f_E, f_W, feq, om_turb);
         f_W_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_E, f_W, feq, om_turb, drho, c2o27);
         wallMomentumX += f_E+f_W_in;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         // f_E_in = getInterpolatedDistributionForNoSlipBC(q, f_W, f_E, feq, om_turb);
         f_E_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_W, f_E, feq, om_turb, drho, c2o27);
         wallMomentumX -= f_W+f_E_in;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         // f_S_in = getInterpolatedDistributionForNoSlipBC(q, f_N, f_S, feq, om_turb);
         f_S_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_N, f_S, feq, om_turb, drho, c2o27);
         wallMomentumY += f_N+f_S_in;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         // f_N_in = getInterpolatedDistributionForNoSlipBC(q, f_S, f_N, feq, om_turb);
         f_N_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_S, f_N, feq, om_turb, drho, c2o27);
         wallMomentumY -= f_S+f_N_in;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         // f_B_in = getInterpolatedDistributionForNoSlipBC(q, f_T, f_B, feq, om_turb);
         f_B_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_T, f_B, feq, om_turb, drho, c2o27);
         wallMomentumZ += f_T+f_B_in;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         // f_T_in = getInterpolatedDistributionForNoSlipBC(q, f_B, f_T, feq, om_turb);
         f_T_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_B, f_T, feq, om_turb, drho, c2o27);
         wallMomentumZ -= f_B+f_T_in;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_SW_in = getInterpolatedDistributionForNoSlipBC(q, f_NE, f_SW, feq, om_turb);
         f_SW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_NE, f_SW, feq, om_turb, drho, c2o27);
         wallMomentumX += f_NE+f_SW_in;
         wallMomentumY += f_NE+f_SW_in;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_NE_in = getInterpolatedDistributionForNoSlipBC(q, f_SW, f_NE, feq, om_turb);
         f_NE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_SW, f_NE, feq, om_turb, drho, c1o54);
         wallMomentumX -= f_SW+f_NE_in;
         wallMomentumY -= f_SW+f_NE_in;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_NW_in = getInterpolatedDistributionForNoSlipBC(q, f_SE, f_NW, feq, om_turb);
         f_NW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_SE, f_NW, feq, om_turb, drho, c1o54);
         wallMomentumX += f_SE+f_NW_in;
         wallMomentumY -= f_SE+f_NW_in;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_SE_in = getInterpolatedDistributionForNoSlipBC(q, f_NW, f_SE, feq, om_turb);
         f_SE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_NW, f_SE, feq, om_turb, drho, c1o54);
         wallMomentumX -= f_NW+f_SE_in;
         wallMomentumY += f_NW+f_SE_in;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_BW_in = getInterpolatedDistributionForNoSlipBC(q, f_TE, f_BW, feq, om_turb);
         f_BW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TE, f_BW, feq, om_turb, drho, c1o54);
         wallMomentumX += f_TE+f_BW_in;
         wallMomentumZ += f_TE+f_BW_in;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_TE_in = getInterpolatedDistributionForNoSlipBC(q, f_BW, f_TE, feq, om_turb);
         f_TE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BW, f_TE, feq, om_turb, drho, c1o54);
         wallMomentumX -= f_BW+f_TE_in;
         wallMomentumZ -= f_BW+f_TE_in;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_TW_in = getInterpolatedDistributionForNoSlipBC(q, f_BE, f_TW, feq, om_turb);
         f_TW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BE, f_TW, feq, om_turb, drho, c1o54);
         wallMomentumX += f_BE+f_TW_in;
         wallMomentumZ -= f_BE+f_TW_in;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_BE_in = getInterpolatedDistributionForNoSlipBC(q, f_TW, f_BE, feq, om_turb);
         f_BE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TW, f_BE, feq, om_turb, drho, c1o54);
         wallMomentumX -= f_TW+f_BE_in;
         wallMomentumZ += f_TW+f_BE_in;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_BS_in = getInterpolatedDistributionForNoSlipBC(q, f_TN, f_BS, feq, om_turb);
         f_BS_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TN, f_BS, feq, om_turb, drho, c1o54);
         wallMomentumY += f_TN+f_BS_in;
         wallMomentumZ += f_TN+f_BS_in;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_TN_in = getInterpolatedDistributionForNoSlipBC(q, f_BS, f_TN, feq, om_turb);
         f_TN_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BS, f_TN, feq, om_turb, drho, c1o54);
         wallMomentumY -= f_BS+f_TN_in;
         wallMomentumZ -= f_BS+f_TN_in;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_TS_in = getInterpolatedDistributionForNoSlipBC(q, f_BN, f_TS, feq, om_turb);
         f_TS_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BN, f_TS, feq, om_turb, drho, c1o54);
         wallMomentumY += f_BN+f_TS_in;
         wallMomentumZ -= f_BN+f_TS_in;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         // f_BN_in = getInterpolatedDistributionForNoSlipBC(q, f_TS, f_BN, feq, om_turb);
         f_BN_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TS, f_BN, feq, om_turb, drho, c1o54);
         wallMomentumY -= f_TS+f_BN_in;
         wallMomentumZ += f_TS+f_BN_in;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_BSW_in = getInterpolatedDistributionForNoSlipBC(q, f_TNE, f_BSW, feq, om_turb);
         f_BSW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TNE, f_BSW, feq, om_turb, drho, c1o216);
         wallMomentumX += f_TNE+f_BSW_in;
         wallMomentumY += f_TNE+f_BSW_in;
         wallMomentumZ += f_TNE+f_BSW_in;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_TNE_in = getInterpolatedDistributionForNoSlipBC(q, f_BSW, f_TNE, feq, om_turb);
         f_TNE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BSW, f_TNE, feq, om_turb, drho, c1o216);
         wallMomentumX -= f_BSW+f_TNE_in;
         wallMomentumY -= f_BSW+f_TNE_in;
         wallMomentumZ -= f_BSW+f_TNE_in;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_TSW_in = getInterpolatedDistributionForNoSlipBC(q, f_BNE, f_TSW, feq, om_turb);
         f_TSW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BNE, f_TSW, feq, om_turb, drho, c1o216);
         wallMomentumX += f_BNE+f_TSW_in;
         wallMomentumY += f_BNE+f_TSW_in;
         wallMomentumZ -= f_BNE+f_TSW_in;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_BNE_in = getInterpolatedDistributionForNoSlipBC(q, f_TSW, f_BNE, feq, om_turb);
         f_BNE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TSW, f_BNE, feq, om_turb, drho, c1o216);
         wallMomentumX -= f_TSW+f_BNE_in;
         wallMomentumY -= f_TSW+f_BNE_in;
         wallMomentumZ += f_TSW+f_BNE_in;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_BNW_in = getInterpolatedDistributionForNoSlipBC(q, f_TSE, f_BNW, feq, om_turb);
         f_BNW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TSE, f_BNW, feq, om_turb, drho, c1o216);
         wallMomentumX += f_TSE+f_BNW_in;
         wallMomentumY -= f_TSE+f_BNW_in;
         wallMomentumZ += f_TSE+f_BNW_in;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_TSE_in = getInterpolatedDistributionForNoSlipBC(q, f_BNW, f_TSE, feq, om_turb);
         f_TSE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BNW, f_TSE, feq, om_turb, drho, c1o216);
         wallMomentumX -= f_BNW+f_TSE_in;
         wallMomentumY += f_BNW+f_TSE_in;
         wallMomentumZ -= f_BNW+f_TSE_in;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_TNW_in = getInterpolatedDistributionForNoSlipBC(q, f_BSE, f_TNW, feq, om_turb);
         f_TNW_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_BSE, f_TNW, feq, om_turb, drho, c1o216);
         wallMomentumX += f_BSE+f_TNW_in;
         wallMomentumY -= f_BSE+f_TNW_in;
         wallMomentumZ -= f_BSE+f_TNW_in;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         // f_BSE_in = getInterpolatedDistributionForNoSlipBC(q, f_TNW, f_BSE, feq, om_turb);
         f_BSE_in = getInterpolatedDistributionForNoSlipWithPressureBC(q, f_TNW, f_BSE, feq, om_turb, drho, c1o216);
         wallMomentumX -= f_TNW+f_BSE_in;
         wallMomentumY += f_TNW+f_BSE_in;
         wallMomentumZ += f_TNW+f_BSE_in;
      }

      // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
         (D.f[dM00])[kw] = f_W_in - (c6o1*c2o27*( VeloX     ))/(c1o1+q);
         wallMomentumX += -(c6o1*c2o27*( VeloX     ))/(c1o1+q);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP00])[ke] = f_E_in - (c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c2o27*(-VeloX     ))/(c1o1+q);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0M0])[ks] = f_S_in - (c6o1*c2o27*( VeloY     ))/(c1o1+q);
         wallMomentumY += - (c6o1*c2o27*( VeloY     ))/(c1o1+q);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0P0])[kn] = f_N_in - (c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         wallMomentumY -=  -(c6o1*c2o27*(-VeloY     ))/(c1o1+q);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00M])[kb] = f_B_in - (c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c2o27*( VeloZ     ))/(c1o1+q);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00P])[kt] = f_T_in - (c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         wallMomentumZ -= -(c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMM0])[ksw] = f_SW_in - (c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         wallMomentumX +=  -(c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         wallMomentumY +=  -(c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPP0])[kne] = f_NE_in - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMP0])[knw] = f_NW_in - (c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         wallMomentumX += -(c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         wallMomentumY -= -(c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPM0])[kse] = f_SE_in - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0M])[kbw] = f_BW_in - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0P])[kte] = f_TE_in - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0P])[ktw] = f_TW_in - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0M])[kbe] = f_BE_in - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MM])[kbs] = f_BS_in - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PP])[ktn] = f_TN_in - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MP])[kts] = f_TS_in - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PM])[kbn] = f_BN_in - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMM])[kbsw] = f_BSW_in - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPP])[ktne] = f_TNE_in - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMP])[ktsw] = f_TSW_in - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPM])[kbne] = f_BNE_in - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPM])[kbnw] = f_BNW_in - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         wallMomentumZ += - (c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMP])[ktse] = f_TSE_in - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumX -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumY += - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPP])[ktnw] = f_TNW_in - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumX += - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumY -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         wallMomentumZ -= - (c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMM])[kbse] = f_BSE_in - (c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
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

//! \}
