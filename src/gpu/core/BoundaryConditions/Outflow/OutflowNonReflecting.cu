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
#include "Utilities/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__host__ __device__ real computeOutflowDistribution(const real* const &f, const real* const &f1, const int dir, const real cs)
{
   return f1[dir] * cs + (c1o1 - cs) * f[dir];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void OutflowNonReflecting_Device(
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
    int direction)
{
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get the node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////

   if(nodeIndex >= numberOfBCnodes) return;

   ////////////////////////////////////////////////////////////////////////////////
   //index
   unsigned int KQK  = k_Q[nodeIndex];
   // unsigned int kzero= KQK;
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
   //index1
   unsigned int K1QK  = k_N[nodeIndex];
   //unsigned int k1zero= K1QK;
   unsigned int k1e   = K1QK;
   unsigned int k1w   = neighborX[K1QK];
   unsigned int k1n   = K1QK;
   unsigned int k1s   = neighborY[K1QK];
   unsigned int k1t   = K1QK;
   unsigned int k1b   = neighborZ[K1QK];
   unsigned int k1sw  = neighborY[k1w];
   unsigned int k1ne  = K1QK;
   unsigned int k1se  = k1s;
   unsigned int k1nw  = k1w;
   unsigned int k1bw  = neighborZ[k1w];
   unsigned int k1te  = K1QK;
   unsigned int k1be  = k1b;
   unsigned int k1tw  = k1w;
   unsigned int k1bs  = neighborZ[k1s];
   unsigned int k1tn  = K1QK;
   unsigned int k1bn  = k1b;
   unsigned int k1ts  = k1s;
   unsigned int k1tse = k1s;
   unsigned int k1bnw = k1bw;
   unsigned int k1tnw = k1w;
   unsigned int k1bse = k1bs;
   unsigned int k1tsw = k1sw;
   unsigned int k1bne = k1b;
   unsigned int k1tne = K1QK;
   unsigned int k1bsw = neighborZ[k1sw];
   ////////////////////////////////////////////////////////////////////////////////
   Distributions27 dist;
   getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
   real f[27], f1[27];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   f1[dP00] = (dist.f[dP00])[k1e   ];
   f1[dM00] = (dist.f[dM00])[k1w   ];
   f1[d0P0] = (dist.f[d0P0])[k1n   ];
   f1[d0M0] = (dist.f[d0M0])[k1s   ];
   f1[d00P] = (dist.f[d00P])[k1t   ];
   f1[d00M] = (dist.f[d00M])[k1b   ];
   f1[dPP0] = (dist.f[dPP0])[k1ne  ];
   f1[dMM0] = (dist.f[dMM0])[k1sw  ];
   f1[dPM0] = (dist.f[dPM0])[k1se  ];
   f1[dMP0] = (dist.f[dMP0])[k1nw  ];
   f1[dP0P] = (dist.f[dP0P])[k1te  ];
   f1[dM0M] = (dist.f[dM0M])[k1bw  ];
   f1[dP0M] = (dist.f[dP0M])[k1be  ];
   f1[dM0P] = (dist.f[dM0P])[k1tw  ];
   f1[d0PP] = (dist.f[d0PP])[k1tn  ];
   f1[d0MM] = (dist.f[d0MM])[k1bs  ];
   f1[d0PM] = (dist.f[d0PM])[k1bn  ];
   f1[d0MP] = (dist.f[d0MP])[k1ts  ];
   // f1[d000] = (dist.f[d000])[k1zero];
   f1[dPPP] = (dist.f[dPPP])[k1tne ];
   f1[dMMP] = (dist.f[dMMP])[k1tsw ];
   f1[dPMP] = (dist.f[dPMP])[k1tse ];
   f1[dMPP] = (dist.f[dMPP])[k1tnw ];
   f1[dPPM] = (dist.f[dPPM])[k1bne ];
   f1[dMMM] = (dist.f[dMMM])[k1bsw ];
   f1[dPMM] = (dist.f[dPMM])[k1bse ];
   f1[dMPM] = (dist.f[dMPM])[k1bnw ];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   f[dP00] = (dist.f[dP00])[ke   ];
   f[dM00] = (dist.f[dM00])[kw   ];
   f[d0P0] = (dist.f[d0P0])[kn   ];
   f[d0M0] = (dist.f[d0M0])[ks   ];
   f[d00P] = (dist.f[d00P])[kt   ];
   f[d00M] = (dist.f[d00M])[kb   ];
   f[dPP0] = (dist.f[dPP0])[kne  ];
   f[dMM0] = (dist.f[dMM0])[ksw  ];
   f[dPM0] = (dist.f[dPM0])[kse  ];
   f[dMP0] = (dist.f[dMP0])[knw  ];
   f[dP0P] = (dist.f[dP0P])[kte  ];
   f[dM0M] = (dist.f[dM0M])[kbw  ];
   f[dP0M] = (dist.f[dP0M])[kbe  ];
   f[dM0P] = (dist.f[dM0P])[ktw  ];
   f[d0PP] = (dist.f[d0PP])[ktn  ];
   f[d0MM] = (dist.f[d0MM])[kbs  ];
   f[d0PM] = (dist.f[d0PM])[kbn  ];
   f[d0MP] = (dist.f[d0MP])[kts  ];
   // f[d000] = (dist.f[d000])[kzero];
   f[dPPP] = (dist.f[dPPP])[ktne ];
   f[dMMP] = (dist.f[dMMP])[ktsw ];
   f[dPMP] = (dist.f[dPMP])[ktse ];
   f[dMPP] = (dist.f[dMPP])[ktnw ];
   f[dPPM] = (dist.f[dPPM])[kbne ];
   f[dMMM] = (dist.f[dMMM])[kbsw ];
   f[dPMM] = (dist.f[dPMM])[kbse ];
   f[dMPM] = (dist.f[dMPM])[kbnw ];
   //////////////////////////////////////////////////////////////////////////


   real cs = c1o1 / sqrtf(c3o1);

   //////////////////////////////////////////////////////////////////////////
   getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);
   switch(direction)
   {
      case dM00:
         (dist.f[dP00])[ke   ] = computeOutflowDistribution(f, f1, dP00, cs);
         (dist.f[dPM0])[kse  ] = computeOutflowDistribution(f, f1, dPM0, cs);
         (dist.f[dPP0])[kne  ] = computeOutflowDistribution(f, f1, dPP0, cs);
         (dist.f[dP0M])[kbe  ] = computeOutflowDistribution(f, f1, dP0M, cs);
         (dist.f[dP0P])[kte  ] = computeOutflowDistribution(f, f1, dP0P, cs);
         (dist.f[dPMP])[ktse ] = computeOutflowDistribution(f, f1, dPMP, cs);
         (dist.f[dPPP])[ktne ] = computeOutflowDistribution(f, f1, dPPP, cs);
         (dist.f[dPMM])[kbse ] = computeOutflowDistribution(f, f1, dPMM, cs);
         (dist.f[dPPM])[kbne ] = computeOutflowDistribution(f, f1, dPPM, cs);
         break;

      case dP00:
         (dist.f[dM00])[kw   ] = computeOutflowDistribution(f, f1, dM00, cs);
         (dist.f[dMM0])[ksw  ] = computeOutflowDistribution(f, f1, dMM0, cs);
         (dist.f[dMP0])[knw  ] = computeOutflowDistribution(f, f1, dMP0, cs);
         (dist.f[dM0M])[kbw  ] = computeOutflowDistribution(f, f1, dM0M, cs);
         (dist.f[dM0P])[ktw  ] = computeOutflowDistribution(f, f1, dM0P, cs);
         (dist.f[dMMP])[ktsw ] = computeOutflowDistribution(f, f1, dMMP, cs);
         (dist.f[dMPP])[ktnw ] = computeOutflowDistribution(f, f1, dMPP, cs);
         (dist.f[dMMM])[kbsw ] = computeOutflowDistribution(f, f1, dMMM, cs);
         (dist.f[dMPM])[kbnw ] = computeOutflowDistribution(f, f1, dMPM, cs);
         break;

      case d0M0:
         (dist.f[d0P0])[kn   ] = computeOutflowDistribution(f, f1, d0P0, cs);
         (dist.f[dPP0])[kne  ] = computeOutflowDistribution(f, f1, dPP0, cs);
         (dist.f[dMP0])[knw  ] = computeOutflowDistribution(f, f1, dMP0, cs);
         (dist.f[d0PP])[ktn  ] = computeOutflowDistribution(f, f1, d0PP, cs);
         (dist.f[d0PM])[kbn  ] = computeOutflowDistribution(f, f1, d0PM, cs);
         (dist.f[dPPP])[ktne ] = computeOutflowDistribution(f, f1, dPPP, cs);
         (dist.f[dMPP])[ktnw ] = computeOutflowDistribution(f, f1, dMPP, cs);
         (dist.f[dPPM])[kbne ] = computeOutflowDistribution(f, f1, dPPM, cs);
         (dist.f[dMPM])[kbnw ] = computeOutflowDistribution(f, f1, dMPM, cs);
         break;

      case d0P0:
         (dist.f[d0M0])[ks   ] = computeOutflowDistribution(f, f1, d0M0, cs);
         (dist.f[dPM0])[kse  ] = computeOutflowDistribution(f, f1, dPM0, cs);
         (dist.f[dMM0])[ksw  ] = computeOutflowDistribution(f, f1, dMM0, cs);
         (dist.f[d0MP])[kts  ] = computeOutflowDistribution(f, f1, d0MP, cs);
         (dist.f[d0MM])[kbs  ] = computeOutflowDistribution(f, f1, d0MM, cs);
         (dist.f[dPMP])[ktse ] = computeOutflowDistribution(f, f1, dPMP, cs);
         (dist.f[dMMP])[ktsw ] = computeOutflowDistribution(f, f1, dMMP, cs);
         (dist.f[dPMM])[kbse ] = computeOutflowDistribution(f, f1, dPMM, cs);
         (dist.f[dMMM])[kbsw ] = computeOutflowDistribution(f, f1, dMMM, cs);
         break;

      case d00M:
         (dist.f[d00P])[kt   ] = computeOutflowDistribution(f, f1, d00P, cs);
         (dist.f[dP0P])[kte  ] = computeOutflowDistribution(f, f1, dP0P, cs);
         (dist.f[dM0P])[ktw  ] = computeOutflowDistribution(f, f1, dM0P, cs);
         (dist.f[d0PP])[ktn  ] = computeOutflowDistribution(f, f1, d0PP, cs);
         (dist.f[d0MP])[kts  ] = computeOutflowDistribution(f, f1, d0MP, cs);
         (dist.f[dPPP])[ktne ] = computeOutflowDistribution(f, f1, dPPP, cs);
         (dist.f[dMPP])[ktnw ] = computeOutflowDistribution(f, f1, dMPP, cs);
         (dist.f[dPMP])[ktse ] = computeOutflowDistribution(f, f1, dPMP, cs);
         (dist.f[dMMP])[ktsw ] = computeOutflowDistribution(f, f1, dMMP, cs);
         break;

      case d00P:
         (dist.f[d00M])[kb   ] = computeOutflowDistribution(f, f1, d00M, cs);
         (dist.f[dP0M])[kbe  ] = computeOutflowDistribution(f, f1, dP0M, cs);
         (dist.f[dM0M])[kbw  ] = computeOutflowDistribution(f, f1, dM0M, cs);
         (dist.f[d0PM])[kbn  ] = computeOutflowDistribution(f, f1, d0PM, cs);
         (dist.f[d0MM])[kbs  ] = computeOutflowDistribution(f, f1, d0MM, cs);
         (dist.f[dPPM])[kbne ] = computeOutflowDistribution(f, f1, dPPM, cs);
         (dist.f[dMPM])[kbnw ] = computeOutflowDistribution(f, f1, dMPM, cs);
         (dist.f[dPMM])[kbse ] = computeOutflowDistribution(f, f1, dPMM, cs);
         (dist.f[dMMM])[kbsw ] = computeOutflowDistribution(f, f1, dMMM, cs);
         break;
      default:
         break;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//! \}
