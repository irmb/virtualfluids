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
//! \addtogroup gpu_PreProcessor PreProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "Calculation/Calculation.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm;
using namespace vf::lbm::dir;


__global__ void InitAdvectionDiffusionCompressible_Device(
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* typeOfGridNode,
    real* concentration,
    real* velocityX,
    real* velocityY,
    real* velocityZ,
    unsigned long long numberOfLBnodes,
    real* distributionsAD,
    bool isEvenTimestep)
{
    //////////////////////////////////////////////////////////////////////////
    //! The initialization is executed in the following steps
    //!
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned  x = threadIdx.x;  // Globaler x-Index
    const unsigned  y = blockIdx.x;   // Globaler y-Index
    const unsigned  z = blockIdx.y;   // Globaler z-Index

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx*(ny*z + y) + x;

    //////////////////////////////////////////////////////////////////////////
    // run for all indices in size_Mat and fluid nodes
    if ((k < numberOfLBnodes) && (typeOfGridNode[k] == GEO_FLUID))
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
        //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 distAD;
        if (isEvenTimestep)
        {
            distAD.f[dP00] = &distributionsAD[dP00 * numberOfLBnodes];
            distAD.f[dM00] = &distributionsAD[dM00 * numberOfLBnodes];
            distAD.f[d0P0] = &distributionsAD[d0P0 * numberOfLBnodes];
            distAD.f[d0M0] = &distributionsAD[d0M0 * numberOfLBnodes];
            distAD.f[d00P] = &distributionsAD[d00P * numberOfLBnodes];
            distAD.f[d00M] = &distributionsAD[d00M * numberOfLBnodes];
            distAD.f[dPP0] = &distributionsAD[dPP0 * numberOfLBnodes];
            distAD.f[dMM0] = &distributionsAD[dMM0 * numberOfLBnodes];
            distAD.f[dPM0] = &distributionsAD[dPM0 * numberOfLBnodes];
            distAD.f[dMP0] = &distributionsAD[dMP0 * numberOfLBnodes];
            distAD.f[dP0P] = &distributionsAD[dP0P * numberOfLBnodes];
            distAD.f[dM0M] = &distributionsAD[dM0M * numberOfLBnodes];
            distAD.f[dP0M] = &distributionsAD[dP0M * numberOfLBnodes];
            distAD.f[dM0P] = &distributionsAD[dM0P * numberOfLBnodes];
            distAD.f[d0PP] = &distributionsAD[d0PP * numberOfLBnodes];
            distAD.f[d0MM] = &distributionsAD[d0MM * numberOfLBnodes];
            distAD.f[d0PM] = &distributionsAD[d0PM * numberOfLBnodes];
            distAD.f[d0MP] = &distributionsAD[d0MP * numberOfLBnodes];
            distAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
            distAD.f[dPPP] = &distributionsAD[dPPP * numberOfLBnodes];
            distAD.f[dMMP] = &distributionsAD[dMMP * numberOfLBnodes];
            distAD.f[dPMP] = &distributionsAD[dPMP * numberOfLBnodes];
            distAD.f[dMPP] = &distributionsAD[dMPP * numberOfLBnodes];
            distAD.f[dPPM] = &distributionsAD[dPPM * numberOfLBnodes];
            distAD.f[dMMM] = &distributionsAD[dMMM * numberOfLBnodes];
            distAD.f[dPMM] = &distributionsAD[dPMM * numberOfLBnodes];
            distAD.f[dMPM] = &distributionsAD[dMPM * numberOfLBnodes];
        }
        else
        {
            distAD.f[dM00] = &distributionsAD[dP00 * numberOfLBnodes];
            distAD.f[dP00] = &distributionsAD[dM00 * numberOfLBnodes];
            distAD.f[d0M0] = &distributionsAD[d0P0 * numberOfLBnodes];
            distAD.f[d0P0] = &distributionsAD[d0M0 * numberOfLBnodes];
            distAD.f[d00M] = &distributionsAD[d00P * numberOfLBnodes];
            distAD.f[d00P] = &distributionsAD[d00M * numberOfLBnodes];
            distAD.f[dMM0] = &distributionsAD[dPP0 * numberOfLBnodes];
            distAD.f[dPP0] = &distributionsAD[dMM0 * numberOfLBnodes];
            distAD.f[dMP0] = &distributionsAD[dPM0 * numberOfLBnodes];
            distAD.f[dPM0] = &distributionsAD[dMP0 * numberOfLBnodes];
            distAD.f[dM0M] = &distributionsAD[dP0P * numberOfLBnodes];
            distAD.f[dP0P] = &distributionsAD[dM0M * numberOfLBnodes];
            distAD.f[dM0P] = &distributionsAD[dP0M * numberOfLBnodes];
            distAD.f[dP0M] = &distributionsAD[dM0P * numberOfLBnodes];
            distAD.f[d0MM] = &distributionsAD[d0PP * numberOfLBnodes];
            distAD.f[d0PP] = &distributionsAD[d0MM * numberOfLBnodes];
            distAD.f[d0MP] = &distributionsAD[d0PM * numberOfLBnodes];
            distAD.f[d0PM] = &distributionsAD[d0MP * numberOfLBnodes];
            distAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
            distAD.f[dMMM] = &distributionsAD[dPPP * numberOfLBnodes];
            distAD.f[dPPM] = &distributionsAD[dMMP * numberOfLBnodes];
            distAD.f[dMPM] = &distributionsAD[dPMP * numberOfLBnodes];
            distAD.f[dPMM] = &distributionsAD[dMPP * numberOfLBnodes];
            distAD.f[dMMP] = &distributionsAD[dPPM * numberOfLBnodes];
            distAD.f[dPPP] = &distributionsAD[dMMM * numberOfLBnodes];
            distAD.f[dMPP] = &distributionsAD[dPMM * numberOfLBnodes];
            distAD.f[dPMP] = &distributionsAD[dMPM * numberOfLBnodes];
        }
        //////////////////////////////////////////////////////////////////////////
        //! - Set local velocities and concetration
        //!
        real conc = concentration[k];
        real  vx1 = velocityX[k];
        real  vx2 = velocityY[k];
        real  vx3 = velocityZ[k];
        //////////////////////////////////////////////////////////////////////////
        //! - Set neighbor indices (necessary for indirect addressing)
        //!
        uint kzero = k;
        uint ke    = k;
        uint kw    = neighborX[k];
        uint kn    = k;
        uint ks    = neighborY[k];
        uint kt    = k;
        uint kb    = neighborZ[k];
        uint ksw   = neighborY[kw];
        uint kne   = k;
        uint kse   = ks;
        uint knw   = kw;
        uint kbw   = neighborZ[kw];
        uint kte   = k;
        uint kbe   = kb;
        uint ktw   = kw;
        uint kbs   = neighborZ[ks];
        uint ktn   = k;
        uint kbn   = kb;
        uint kts   = ks;
        uint ktse  = ks;
        uint kbnw  = kbw;
        uint ktnw  = kw;
        uint kbse  = kbs;
        uint ktsw  = ksw;
        uint kbne  = kb;
        uint ktne  = k;
        uint kbsw  = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        //! - Calculate the equilibrium and set the distributions
        //!
        real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);

        (distAD.f[d000])[kzero] = c8o27  * conc * (c1o1 - cu_sq);
        (distAD.f[dP00])[ke   ] = c2o27  * conc * (c1o1 + c3o1 * ( vx1            ) + c9o2 * ( vx1            ) * ( vx1            ) - cu_sq);
        (distAD.f[dM00])[kw   ] = c2o27  * conc * (c1o1 + c3o1 * (-vx1            ) + c9o2 * (-vx1            ) * (-vx1            ) - cu_sq);
        (distAD.f[d0P0])[kn   ] = c2o27  * conc * (c1o1 + c3o1 * (       vx2      ) + c9o2 * (       vx2      ) * (       vx2      ) - cu_sq);
        (distAD.f[d0M0])[ks   ] = c2o27  * conc * (c1o1 + c3o1 * (     - vx2      ) + c9o2 * (     - vx2      ) * (     - vx2      ) - cu_sq);
        (distAD.f[d00P])[kt   ] = c2o27  * conc * (c1o1 + c3o1 * (             vx3) + c9o2 * (             vx3) * (             vx3) - cu_sq);
        (distAD.f[d00M])[kb   ] = c2o27  * conc * (c1o1 + c3o1 * (           - vx3) + c9o2 * (           - vx3) * (           - vx3) - cu_sq);
        (distAD.f[dPP0])[kne  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 + vx2      ) + c9o2 * ( vx1 + vx2      ) * ( vx1 + vx2      ) - cu_sq);
        (distAD.f[dMM0])[ksw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 - vx2      ) + c9o2 * (-vx1 - vx2      ) * (-vx1 - vx2      ) - cu_sq);
        (distAD.f[dPM0])[kse  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 - vx2      ) + c9o2 * ( vx1 - vx2      ) * ( vx1 - vx2      ) - cu_sq);
        (distAD.f[dMP0])[knw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 + vx2      ) + c9o2 * (-vx1 + vx2      ) * (-vx1 + vx2      ) - cu_sq);
        (distAD.f[dP0P])[kte  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       + vx3) + c9o2 * ( vx1       + vx3) * ( vx1       + vx3) - cu_sq);
        (distAD.f[dM0M])[kbw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       - vx3) + c9o2 * (-vx1       - vx3) * (-vx1       - vx3) - cu_sq);
        (distAD.f[dP0M])[kbe  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       - vx3) + c9o2 * ( vx1       - vx3) * ( vx1       - vx3) - cu_sq);
        (distAD.f[dM0P])[ktw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       + vx3) + c9o2 * (-vx1       + vx3) * (-vx1       + vx3) - cu_sq);
        (distAD.f[d0PP])[ktn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 + vx3) + c9o2 * (       vx2 + vx3) * (       vx2 + vx3) - cu_sq);
        (distAD.f[d0MM])[kbs  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 - vx3) + c9o2 * (     - vx2 - vx3) * (     - vx2 - vx3) - cu_sq);
        (distAD.f[d0PM])[kbn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 - vx3) + c9o2 * (       vx2 - vx3) * (       vx2 - vx3) - cu_sq);
        (distAD.f[d0MP])[kts  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 + vx3) + c9o2 * (     - vx2 + vx3) * (     - vx2 + vx3) - cu_sq);
        (distAD.f[dPPP])[ktne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 + vx3) + c9o2 * ( vx1 + vx2 + vx3) * ( vx1 + vx2 + vx3) - cu_sq);
        (distAD.f[dMMM])[kbsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
        (distAD.f[dPPM])[kbne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 - vx3) + c9o2 * ( vx1 + vx2 - vx3) * ( vx1 + vx2 - vx3) - cu_sq);
        (distAD.f[dMMP])[ktsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
        (distAD.f[dPMP])[ktse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 + vx3) + c9o2 * ( vx1 - vx2 + vx3) * ( vx1 - vx2 + vx3) - cu_sq);
        (distAD.f[dMPM])[kbnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
        (distAD.f[dPMM])[kbse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 - vx3) + c9o2 * ( vx1 - vx2 - vx3) * ( vx1 - vx2 - vx3) - cu_sq);
        (distAD.f[dMPP])[ktnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
    }
}

//! \}
