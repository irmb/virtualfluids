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
//! \author Martin Schoenherr
//=======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"


__global__ void InitK18K20NavierStokesCompressible_Device(unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int* geoD,
    real* rho,
    real* ux,
    real* uy,
    real* uz,
    unsigned int size_Mat,
    real* G6,
    bool EvenOrOdd)
{
    ////////////////////////////////////////////////////////////////////////////////
    const unsigned  x = threadIdx.x;  // Globaler x-Index 
    const unsigned  y = blockIdx.x;   // Globaler y-Index 
    const unsigned  z = blockIdx.y;   // Globaler z-Index 

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx*(ny*z + y) + x;
    //////////////////////////////////////////////////////////////////////////

    if (k < size_Mat)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = geoD[k];

        if (BC != GEO_SOLID &&  BC != GEO_VOID)
        {
            Distributions6 D;
            if (EvenOrOdd == true)
            {
                D.g[dP00] = &G6[dP00   *size_Mat];
                D.g[dM00] = &G6[dM00   *size_Mat];
                D.g[d0P0] = &G6[d0P0   *size_Mat];
                D.g[d0M0] = &G6[d0M0   *size_Mat];
                D.g[d00P] = &G6[d00P   *size_Mat];
                D.g[d00M] = &G6[d00M   *size_Mat];
            }
            else
            {
                D.g[dM00] = &G6[dP00   *size_Mat];
                D.g[dP00] = &G6[dM00   *size_Mat];
                D.g[d0M0] = &G6[d0P0   *size_Mat];
                D.g[d0P0] = &G6[d0M0   *size_Mat];
                D.g[d00M] = &G6[d00P   *size_Mat];
                D.g[d00P] = &G6[d00M   *size_Mat];
            }
            //////////////////////////////////////////////////////////////////////////
            //index
            //////////////////////////////////////////////////////////////////////////
            // unsigned int kzero = k;
            unsigned int ke = k;
            unsigned int kw = neighborX[k];
            unsigned int kn = k;
            unsigned int ks = neighborY[k];
            unsigned int kt = k;
            unsigned int kb = neighborZ[k];
            //////////////////////////////////////////////////////////////////////////

            (D.g[dP00])[ke] = 0.0f;
            (D.g[dM00])[kw] = 0.0f;
            (D.g[d0P0])[kn] = 0.0f;
            (D.g[d0M0])[ks] = 0.0f;
            (D.g[d00P])[kt] = 0.0f;
            (D.g[d00M])[kb] = 0.0f;
        }
    }
}