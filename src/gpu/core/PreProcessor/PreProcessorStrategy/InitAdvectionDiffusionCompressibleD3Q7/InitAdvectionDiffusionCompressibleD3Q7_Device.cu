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

__global__ void InitAdvectionDiffusionCompressibleD3Q7_Device(unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int* geoD,
    real* Conc,
    real* ux,
    real* uy,
    real* uz,
    unsigned int size_Mat,
    real* DD7,
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

    if (k<size_Mat)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = geoD[k];

        if (BC != GEO_SOLID && BC != GEO_VOID)
        {
            Distributions7 D7;
            if (EvenOrOdd == true)
            {
                D7.f[0] = &DD7[0 * size_Mat];
                D7.f[1] = &DD7[1 * size_Mat];
                D7.f[2] = &DD7[2 * size_Mat];
                D7.f[3] = &DD7[3 * size_Mat];
                D7.f[4] = &DD7[4 * size_Mat];
                D7.f[5] = &DD7[5 * size_Mat];
                D7.f[6] = &DD7[6 * size_Mat];
            }
            else
            {
                D7.f[0] = &DD7[0 * size_Mat];
                D7.f[2] = &DD7[1 * size_Mat];
                D7.f[1] = &DD7[2 * size_Mat];
                D7.f[4] = &DD7[3 * size_Mat];
                D7.f[3] = &DD7[4 * size_Mat];
                D7.f[6] = &DD7[5 * size_Mat];
                D7.f[5] = &DD7[6 * size_Mat];
            }
            //////////////////////////////////////////////////////////////////////////
            real ConcD = Conc[k];
            real   vx1 = ux[k];
            real   vx2 = uy[k];
            real   vx3 = uz[k];
            real lambdaD = -c3o1 + sqrt(c3o1);
            real Diffusivity = c1o20;
            real Lam = -(c1o2 + c1o1 / lambdaD);
            real nue_d = Lam / c3o1;
            real ae = Diffusivity / nue_d - c1o1;
            real ux_sq = vx1 * vx1;
            real uy_sq = vx2 * vx2;
            real uz_sq = vx3 * vx3;
            //////////////////////////////////////////////////////////////////////////
            //index
            //////////////////////////////////////////////////////////////////////////
            unsigned int kzero = k;
            unsigned int ke = k;
            unsigned int kw = neighborX[k];
            unsigned int kn = k;
            unsigned int ks = neighborY[k];
            unsigned int kt = k;
            unsigned int kb = neighborZ[k];
            //////////////////////////////////////////////////////////////////////////

            (D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-c3o1)) - (ux_sq + uy_sq + uz_sq));
            (D7.f[1])[ke] = ConcD*(c1o6*(ae + c1o1) + c1o2*(ux_sq)+vx1*c1o2);
            (D7.f[2])[kw] = ConcD*(c1o6*(ae + c1o1) + c1o2*(ux_sq)-vx1*c1o2);
            (D7.f[3])[kn] = ConcD*(c1o6*(ae + c1o1) + c1o2*(uy_sq)+vx2*c1o2);
            (D7.f[4])[ks] = ConcD*(c1o6*(ae + c1o1) + c1o2*(uy_sq)-vx2*c1o2);
            (D7.f[5])[kt] = ConcD*(c1o6*(ae + c1o1) + c1o2*(uz_sq)+vx3*c1o2);
            (D7.f[6])[kb] = ConcD*(c1o6*(ae + c1o1) + c1o2*(uz_sq)-vx3*c1o2);
        }
    }
}