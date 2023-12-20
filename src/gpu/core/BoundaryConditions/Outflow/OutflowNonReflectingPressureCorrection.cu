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
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "lbm/MacroscopicQuantities.h"
#include "Utilities/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;


__host__ __device__ real computeOutflowDistribution(const real* const &f, const real* const &f1, const int dir, const real rhoCorrection, const real cs, const real weight)
{
   return f1[dir  ] * cs + (c1o1 - cs) * f[dir  ] - weight *rhoCorrection;
}

__global__ void OutflowNonReflectingPressureCorrection_Device(
    real* rhoBC,
    real* distributions,
    int* k_Q,
    int* k_N,
    int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    int direction,
    real densityCorrectionFactor)
{
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get the node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////

   if( nodeIndex >= numberOfBCnodes ) return;

   ////////////////////////////////////////////////////////////////////////////////
   //index

   uint k_000 = k_Q[nodeIndex];
   uint k_M00 = neighborX[k_000];
   uint k_0M0 = neighborY[k_000];
   uint k_00M = neighborZ[k_000];
   uint k_MM0 = neighborY[k_M00];
   uint k_M0M = neighborZ[k_M00];
   uint k_0MM = neighborZ[k_0M0];
   uint k_MMM = neighborZ[k_MM0];

   ////////////////////////////////////////////////////////////////////////////////
   //index of neighbor
   uint kN_000 = k_N[nodeIndex];
   uint kN_M00 = neighborX[k_000];
   uint kN_0M0 = neighborY[k_000];
   uint kN_00M = neighborZ[k_000];
   uint kN_MM0 = neighborY[k_M00];
   uint kN_M0M = neighborZ[k_M00];
   uint kN_0MM = neighborZ[k_0M0];
   uint kN_MMM = neighborZ[k_MM0];
   ////////////////////////////////////////////////////////////////////////////////
   Distributions27 dist;
   getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
   real f[27], fN[27];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   f[d000] = (dist.f[d000])[k_000];
   f[dP00] = (dist.f[dP00])[k_000];
   f[dM00] = (dist.f[dM00])[k_M00];
   f[d0P0] = (dist.f[d0P0])[k_000];
   f[d0M0] = (dist.f[d0M0])[k_0M0];
   f[d00P] = (dist.f[d00P])[k_000];
   f[d00M] = (dist.f[d00M])[k_00M];
   f[dPP0] = (dist.f[dPP0])[k_000];
   f[dMM0] = (dist.f[dMM0])[k_MM0];
   f[dPM0] = (dist.f[dPM0])[k_0M0];
   f[dMP0] = (dist.f[dMP0])[k_M00];
   f[dP0P] = (dist.f[dP0P])[k_000];
   f[dM0M] = (dist.f[dM0M])[k_M0M];
   f[dP0M] = (dist.f[dP0M])[k_00M];
   f[dM0P] = (dist.f[dM0P])[k_M00];
   f[d0PP] = (dist.f[d0PP])[k_000];
   f[d0MM] = (dist.f[d0MM])[k_0MM];
   f[d0PM] = (dist.f[d0PM])[k_00M];
   f[d0MP] = (dist.f[d0MP])[k_0M0];
   f[dPPP] = (dist.f[dPPP])[k_000];
   f[dMPP] = (dist.f[dMPP])[k_M00];
   f[dPMP] = (dist.f[dPMP])[k_0M0];
   f[dMMP] = (dist.f[dMMP])[k_MM0];
   f[dPPM] = (dist.f[dPPM])[k_00M];
   f[dMPM] = (dist.f[dMPM])[k_M0M];
   f[dPMM] = (dist.f[dPMM])[k_0MM];
   f[dMMM] = (dist.f[dMMM])[k_MMM];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   fN[d000] = (dist.f[d000])[kN_000];
   fN[dP00] = (dist.f[dP00])[kN_000];
   fN[dM00] = (dist.f[dM00])[kN_M00];
   fN[d0P0] = (dist.f[d0P0])[kN_000];
   fN[d0M0] = (dist.f[d0M0])[kN_0M0];
   fN[d00P] = (dist.f[d00P])[kN_000];
   fN[d00M] = (dist.f[d00M])[kN_00M];
   fN[dPP0] = (dist.f[dPP0])[kN_000];
   fN[dMM0] = (dist.f[dMM0])[kN_MM0];
   fN[dPM0] = (dist.f[dPM0])[kN_0M0];
   fN[dMP0] = (dist.f[dMP0])[kN_M00];
   fN[dP0P] = (dist.f[dP0P])[kN_000];
   fN[dM0M] = (dist.f[dM0M])[kN_M0M];
   fN[dP0M] = (dist.f[dP0M])[kN_00M];
   fN[dM0P] = (dist.f[dM0P])[kN_M00];
   fN[d0PP] = (dist.f[d0PP])[kN_000];
   fN[d0MM] = (dist.f[d0MM])[kN_0MM];
   fN[d0PM] = (dist.f[d0PM])[kN_00M];
   fN[d0MP] = (dist.f[d0MP])[kN_0M0];
   fN[dPPP] = (dist.f[dPPP])[kN_000];
   fN[dMPP] = (dist.f[dMPP])[kN_M00];
   fN[dPMP] = (dist.f[dPMP])[kN_0M0];
   fN[dMMP] = (dist.f[dMMP])[kN_MM0];
   fN[dPPM] = (dist.f[dPPM])[kN_00M];
   fN[dMPM] = (dist.f[dMPM])[kN_M0M];
   fN[dPMM] = (dist.f[dPMM])[kN_0MM];
   fN[dMMM] = (dist.f[dMMM])[kN_MMM];
   //////////////////////////////////////////////////////////////////////////
   real drho = vf::lbm::getDensity(f);

   real rhoCorrection = densityCorrectionFactor*drho;

   real cs = c1o1 / sqrtf(c3o1);

   getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

   switch(direction)
   {
      case dM00:
         (dist.f[dP00])[k_000] = computeOutflowDistribution(f, fN, dP00  , rhoCorrection, cs, c2o27);
         (dist.f[dPM0])[k_0M0] = computeOutflowDistribution(f, fN, dPM0, rhoCorrection, cs, c1o54);
         (dist.f[dPP0])[k_000] = computeOutflowDistribution(f, fN, dPP0, rhoCorrection, cs, c1o54);
         (dist.f[dP0M])[k_00M] = computeOutflowDistribution(f, fN, dP0M, rhoCorrection, cs, c1o54);
         (dist.f[dP0P])[k_000] = computeOutflowDistribution(f, fN, dP0P, rhoCorrection, cs, c1o54);
         (dist.f[dPMP])[k_0M0] = computeOutflowDistribution(f, fN, dPMP, rhoCorrection, cs, c1o216);
         (dist.f[dPPP])[k_000] = computeOutflowDistribution(f, fN, dPPP, rhoCorrection, cs, c1o216);
         (dist.f[dPMM])[k_0MM] = computeOutflowDistribution(f, fN, dPMM, rhoCorrection, cs, c1o216);
         (dist.f[dPPM])[k_00M] = computeOutflowDistribution(f, fN, dPPM, rhoCorrection, cs, c1o216);
         break;

      case dP00:
         (dist.f[dM00])[k_M00] = computeOutflowDistribution(f, fN, dM00, rhoCorrection, cs, c2o27);
         (dist.f[dMM0])[k_MM0] = computeOutflowDistribution(f, fN, dMM0, rhoCorrection, cs, c1o54);
         (dist.f[dMP0])[k_M00] = computeOutflowDistribution(f, fN, dMP0, rhoCorrection, cs, c1o54);
         (dist.f[dM0M])[k_M0M] = computeOutflowDistribution(f, fN, dM0M, rhoCorrection, cs, c1o54);
         (dist.f[dM0P])[k_M00] = computeOutflowDistribution(f, fN, dM0P, rhoCorrection, cs, c1o54);
         (dist.f[dMMP])[k_MM0] = computeOutflowDistribution(f, fN, dMMP, rhoCorrection, cs, c1o216);
         (dist.f[dMPP])[k_M00] = computeOutflowDistribution(f, fN, dMPP, rhoCorrection, cs, c1o216);
         (dist.f[dMMM])[k_MMM] = computeOutflowDistribution(f, fN, dMMM, rhoCorrection, cs, c1o216);
         (dist.f[dMPM])[k_M0M] = computeOutflowDistribution(f, fN, dMPM, rhoCorrection, cs, c1o216);
         break;

      case d0M0:
         (dist.f[d0P0])[k_000] = computeOutflowDistribution(f, fN, d0P0, rhoCorrection, cs, c2o27);
         (dist.f[dPP0])[k_000] = computeOutflowDistribution(f, fN, dPP0, rhoCorrection, cs, c1o54);
         (dist.f[dMP0])[k_M00] = computeOutflowDistribution(f, fN, dMP0, rhoCorrection, cs, c1o54);
         (dist.f[d0PP])[k_000] = computeOutflowDistribution(f, fN, d0PP, rhoCorrection, cs, c1o54);
         (dist.f[d0PM])[k_00M] = computeOutflowDistribution(f, fN, d0PM, rhoCorrection, cs, c1o54);
         (dist.f[dPPP])[k_000] = computeOutflowDistribution(f, fN, dPPP, rhoCorrection, cs, c1o216);
         (dist.f[dMPP])[k_M00] = computeOutflowDistribution(f, fN, dMPP, rhoCorrection, cs, c1o216);
         (dist.f[dPPM])[k_00M] = computeOutflowDistribution(f, fN, dPPM, rhoCorrection, cs, c1o216);
         (dist.f[dMPM])[k_M0M] = computeOutflowDistribution(f, fN, dMPM, rhoCorrection, cs, c1o216);
         break;

      case d0P0:
         (dist.f[d0M0])[k_0M0] =computeOutflowDistribution(f, fN, d0M0, rhoCorrection, cs, c2o27);
         (dist.f[dPM0])[k_0M0] =computeOutflowDistribution(f, fN, dPM0, rhoCorrection, cs, c1o54);
         (dist.f[dMM0])[k_MM0] =computeOutflowDistribution(f, fN, dMM0, rhoCorrection, cs, c1o54);
         (dist.f[d0MP])[k_0M0] =computeOutflowDistribution(f, fN, d0MP, rhoCorrection, cs, c1o54);
         (dist.f[d0MM])[k_0MM] =computeOutflowDistribution(f, fN, d0MM, rhoCorrection, cs, c1o54);
         (dist.f[dPMP])[k_0M0] =computeOutflowDistribution(f, fN, dPMP, rhoCorrection, cs, c1o216);
         (dist.f[dMMP])[k_MM0] =computeOutflowDistribution(f, fN, dMMP, rhoCorrection, cs, c1o216);
         (dist.f[dPMM])[k_0MM] =computeOutflowDistribution(f, fN, dPMM, rhoCorrection, cs, c1o216);
         (dist.f[dMMM])[k_MMM] =computeOutflowDistribution(f, fN, dMMM, rhoCorrection, cs, c1o216);
         break;

      case d00M:
         (dist.f[d00P])[k_000] = computeOutflowDistribution(f, fN, d00P, rhoCorrection, cs, c2o27);
         (dist.f[dP0P])[k_000] = computeOutflowDistribution(f, fN, dP0P, rhoCorrection, cs, c1o54);
         (dist.f[dM0P])[k_M00] = computeOutflowDistribution(f, fN, dM0P, rhoCorrection, cs, c1o54);
         (dist.f[d0PP])[k_000] = computeOutflowDistribution(f, fN, d0PP, rhoCorrection, cs, c1o54);
         (dist.f[d0MP])[k_0M0] = computeOutflowDistribution(f, fN, d0MP, rhoCorrection, cs, c1o54);
         (dist.f[dPPP])[k_000] = computeOutflowDistribution(f, fN, dPPP, rhoCorrection, cs, c1o216);
         (dist.f[dMPP])[k_M00] = computeOutflowDistribution(f, fN, dMPP, rhoCorrection, cs, c1o216);
         (dist.f[dPMP])[k_0M0] = computeOutflowDistribution(f, fN, dPMP, rhoCorrection, cs, c1o216);
         (dist.f[dMMP])[k_MM0] = computeOutflowDistribution(f, fN, dMMP, rhoCorrection, cs, c1o216);
         break;

      case d00P:
         (dist.f[d00M])[k_00M] = computeOutflowDistribution(f, fN, d00M, rhoCorrection, cs, c2o27);
         (dist.f[dP0M])[k_00M] = computeOutflowDistribution(f, fN, dP0M, rhoCorrection, cs, c1o54);
         (dist.f[dM0M])[k_M0M] = computeOutflowDistribution(f, fN, dM0M, rhoCorrection, cs, c1o54);
         (dist.f[d0PM])[k_00M] = computeOutflowDistribution(f, fN, d0PM, rhoCorrection, cs, c1o54);
         (dist.f[d0MM])[k_0MM] = computeOutflowDistribution(f, fN, d0MM, rhoCorrection, cs, c1o54);
         (dist.f[dPPM])[k_00M] = computeOutflowDistribution(f, fN, dPPM, rhoCorrection, cs, c1o216);
         (dist.f[dMPM])[k_M0M] = computeOutflowDistribution(f, fN, dMPM, rhoCorrection, cs, c1o216);
         (dist.f[dPMM])[k_0MM] = computeOutflowDistribution(f, fN, dPMM, rhoCorrection, cs, c1o216);
         (dist.f[dMMM])[k_MMM] = computeOutflowDistribution(f, fN, dMMM, rhoCorrection, cs, c1o216);
         break;
      default:
         break;
   }
}


//! \}
