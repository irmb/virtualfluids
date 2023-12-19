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
#include <basics/constants/NumericConstants.h>
#include "Utilities/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void NoSlipBounceBack_Device(
    real* distributions, 
    int* subgridDistanceIndices, 
    real* subgridDistances,
    unsigned int numberOfBCnodes, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The no-slip boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   //////////////////////////////////////////////////////////////////////////
   // run for all indices in size of boundary condition (numberOfBCnodes)
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
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - rewrite distributions if there is a sub-grid distance (q) in same direction
      real q;
      q = (subgridD.q[dP00])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dM00])[kw  ]=f_E  ;
      q = (subgridD.q[dM00])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dP00])[ke  ]=f_W  ;
      q = (subgridD.q[d0P0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0M0])[ks  ]=f_N  ;
      q = (subgridD.q[d0M0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0P0])[kn  ]=f_S  ;
      q = (subgridD.q[d00P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d00M])[kb  ]=f_T  ;
      q = (subgridD.q[d00M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d00P])[kt  ]=f_B  ;
      q = (subgridD.q[dPP0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMM0])[ksw ]=f_NE ;
      q = (subgridD.q[dMM0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPP0])[kne ]=f_SW ;
      q = (subgridD.q[dPM0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMP0])[knw ]=f_SE ;
      q = (subgridD.q[dMP0])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPM0])[kse ]=f_NW ;
      q = (subgridD.q[dP0P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dM0M])[kbw ]=f_TE ;
      q = (subgridD.q[dM0M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dP0P])[kte ]=f_BW ;
      q = (subgridD.q[dP0M])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dM0P])[ktw ]=f_BE ;
      q = (subgridD.q[dM0P])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dP0M])[kbe ]=f_TW ;
      q = (subgridD.q[d0PP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0MM])[kbs ]=f_TN ;
      q = (subgridD.q[d0MM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0PP])[ktn ]=f_BS ;
      q = (subgridD.q[d0PM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0MP])[kts ]=f_BN ;
      q = (subgridD.q[d0MP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[d0PM])[kbn ]=f_TS ;
      q = (subgridD.q[dPPP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMMM])[kbsw]=f_TNE;
      q = (subgridD.q[dMMM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPPP])[ktne]=f_BSW;
      q = (subgridD.q[dPPM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMMP])[ktsw]=f_BNE;
      q = (subgridD.q[dMMP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPPM])[kbne]=f_TSW;
      q = (subgridD.q[dPMP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMPM])[kbnw]=f_TSE;
      q = (subgridD.q[dMPM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPMP])[ktse]=f_BNW;
      q = (subgridD.q[dPMM])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dMPP])[ktnw]=f_BSE;
      q = (subgridD.q[dMPP])[nodeIndex];   if (q>=c0o1 && q<=c1o1)    (dist.f[dPMM])[kbse]=f_TNW;
   }
}


//! \}
