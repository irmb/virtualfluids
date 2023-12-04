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
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////////
__global__ void SlipCompressible_Device(
    real* distributions, 
    int* subgridDistanceIndices, 
    real* subgridDistances,
    unsigned int numberOfBCnodes,
    real omega, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   //! The slip boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   if(nodeIndex < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[nodeIndex];
      unsigned int kzero= indexOfBCnode;
      unsigned int ke   = indexOfBCnode;
      unsigned int kw   = neighborX[indexOfBCnode];
      unsigned int kn   = indexOfBCnode;
      unsigned int ks   = neighborY[indexOfBCnode];
      unsigned int kt   = indexOfBCnode;
      unsigned int kb   = neighborZ[indexOfBCnode];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = indexOfBCnode;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = indexOfBCnode;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = indexOfBCnode;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = indexOfBCnode;
      unsigned int kbsw = neighborZ[ksw];
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real f_W    = (dist.f[dP00])[ke   ];
      real f_E    = (dist.f[dM00])[kw   ];
      real f_S    = (dist.f[d0P0])[kn   ];
      real f_N    = (dist.f[d0M0])[ks   ];
      real f_B    = (dist.f[d00P])[kt   ];
      real f_T    = (dist.f[d00M])[kb   ];
      real f_SW   = (dist.f[dPP0])[kne  ];
      real f_NE   = (dist.f[dMM0])[ksw  ];
      real f_NW   = (dist.f[dPM0])[kse  ];
      real f_SE   = (dist.f[dMP0])[knw  ];
      real f_BW   = (dist.f[dP0P])[kte  ];
      real f_TE   = (dist.f[dM0M])[kbw  ];
      real f_TW   = (dist.f[dP0M])[kbe  ];
      real f_BE   = (dist.f[dM0P])[ktw  ];
      real f_BS   = (dist.f[d0PP])[ktn  ];
      real f_TN   = (dist.f[d0MM])[kbs  ];
      real f_TS   = (dist.f[d0PM])[kbn  ];
      real f_BN   = (dist.f[d0MP])[kts  ];
      real f_BSW  = (dist.f[dPPP])[ktne ];
      real f_BNE  = (dist.f[dMMP])[ktsw ];
      real f_BNW  = (dist.f[dPMP])[ktse ];
      real f_BSE  = (dist.f[dMPP])[ktnw ];
      real f_TSW  = (dist.f[dPPM])[kbne ];
      real f_TNE  = (dist.f[dMMM])[kbsw ];
      real f_TNW  = (dist.f[dPMM])[kbse ];
      real f_TSE  = (dist.f[dMPM])[kbnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[d0M0])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[d0P0])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[d00M])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[d00P])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[dPP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[dMM0])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dMM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[dPP0])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dPM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[dMP0])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dMP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[dPM0])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[dM0M])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[dP0P])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[dM0P])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[dP0M])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[d0MM])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[d0PP])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, omega, velocityBC, c1o54);
      }


      q = (subgridD.q[d0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[d0MP])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[d0PM])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dPPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[dMMM])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[dPPP])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, omega, velocityBC, c1o216);
      }


      q = (subgridD.q[dPPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[dMMP])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[dPPM])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[dMPM])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[dPMP])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[dMPP])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[dPMM])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, omega, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
