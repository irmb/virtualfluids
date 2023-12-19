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
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "CP27.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "Calculation/Calculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__ void CalcCP27(real* DD, 
                                    int* cpIndex, 
                                    int nonCp, 
                                    double *cpPress,
                                    unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned long long numberOfLBnodes, 
                                    bool isEvenTimestep)
{
    Distributions27 D;
    if (isEvenTimestep==true)
    {
        D.f[dP00] = &DD[dP00 * numberOfLBnodes];
        D.f[dM00] = &DD[dM00 * numberOfLBnodes];
        D.f[d0P0] = &DD[d0P0 * numberOfLBnodes];
        D.f[d0M0] = &DD[d0M0 * numberOfLBnodes];
        D.f[d00P] = &DD[d00P * numberOfLBnodes];
        D.f[d00M] = &DD[d00M * numberOfLBnodes];
        D.f[dPP0] = &DD[dPP0 * numberOfLBnodes];
        D.f[dMM0] = &DD[dMM0 * numberOfLBnodes];
        D.f[dPM0] = &DD[dPM0 * numberOfLBnodes];
        D.f[dMP0] = &DD[dMP0 * numberOfLBnodes];
        D.f[dP0P] = &DD[dP0P * numberOfLBnodes];
        D.f[dM0M] = &DD[dM0M * numberOfLBnodes];
        D.f[dP0M] = &DD[dP0M * numberOfLBnodes];
        D.f[dM0P] = &DD[dM0P * numberOfLBnodes];
        D.f[d0PP] = &DD[d0PP * numberOfLBnodes];
        D.f[d0MM] = &DD[d0MM * numberOfLBnodes];
        D.f[d0PM] = &DD[d0PM * numberOfLBnodes];
        D.f[d0MP] = &DD[d0MP * numberOfLBnodes];
        D.f[d000] = &DD[d000 * numberOfLBnodes];
        D.f[dPPP] = &DD[dPPP * numberOfLBnodes];
        D.f[dMMP] = &DD[dMMP * numberOfLBnodes];
        D.f[dPMP] = &DD[dPMP * numberOfLBnodes];
        D.f[dMPP] = &DD[dMPP * numberOfLBnodes];
        D.f[dPPM] = &DD[dPPM * numberOfLBnodes];
        D.f[dMMM] = &DD[dMMM * numberOfLBnodes];
        D.f[dPMM] = &DD[dPMM * numberOfLBnodes];
        D.f[dMPM] = &DD[dMPM * numberOfLBnodes];
    } 
    else
    {
        D.f[dM00] = &DD[dP00 * numberOfLBnodes];
        D.f[dP00] = &DD[dM00 * numberOfLBnodes];
        D.f[d0M0] = &DD[d0P0 * numberOfLBnodes];
        D.f[d0P0] = &DD[d0M0 * numberOfLBnodes];
        D.f[d00M] = &DD[d00P * numberOfLBnodes];
        D.f[d00P] = &DD[d00M * numberOfLBnodes];
        D.f[dMM0] = &DD[dPP0 * numberOfLBnodes];
        D.f[dPP0] = &DD[dMM0 * numberOfLBnodes];
        D.f[dMP0] = &DD[dPM0 * numberOfLBnodes];
        D.f[dPM0] = &DD[dMP0 * numberOfLBnodes];
        D.f[dM0M] = &DD[dP0P * numberOfLBnodes];
        D.f[dP0P] = &DD[dM0M * numberOfLBnodes];
        D.f[dM0P] = &DD[dP0M * numberOfLBnodes];
        D.f[dP0M] = &DD[dM0P * numberOfLBnodes];
        D.f[d0MM] = &DD[d0PP * numberOfLBnodes];
        D.f[d0PP] = &DD[d0MM * numberOfLBnodes];
        D.f[d0MP] = &DD[d0PM * numberOfLBnodes];
        D.f[d0PM] = &DD[d0MP * numberOfLBnodes];
        D.f[d000] = &DD[d000 * numberOfLBnodes];
        D.f[dPPP] = &DD[dMMM * numberOfLBnodes];
        D.f[dMMP] = &DD[dPPM * numberOfLBnodes];
        D.f[dPMP] = &DD[dMPM * numberOfLBnodes];
        D.f[dMPP] = &DD[dPMM * numberOfLBnodes];
        D.f[dPPM] = &DD[dMMP * numberOfLBnodes];
        D.f[dMMM] = &DD[dPPP * numberOfLBnodes];
        D.f[dPMM] = &DD[dMPP * numberOfLBnodes];
        D.f[dMPM] = &DD[dPMP * numberOfLBnodes];
    }
    ////////////////////////////////////////////////////////////////////////////////
    const unsigned  x = threadIdx.x;  // Globaler x-Index 
    const unsigned  y = blockIdx.x;   // Globaler y-Index 
    const unsigned  z = blockIdx.y;   // Globaler z-Index 

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx*(ny*z + y) + x;
    //////////////////////////////////////////////////////////////////////////

    if(k<nonCp)
    {
        ////////////////////////////////////////////////////////////////////////////////
        //index
        unsigned int KQK  = cpIndex[k];
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
        double PressCP;

        PressCP  =   (D.f[dP00])[ke  ]+ (D.f[dM00])[kw  ]+ 
                     (D.f[d0P0])[kn  ]+ (D.f[d0M0])[ks  ]+
                     (D.f[d00P])[kt  ]+ (D.f[d00M])[kb  ]+
                     (D.f[dPP0])[kne ]+ (D.f[dMM0])[ksw ]+
                     (D.f[dPM0])[kse ]+ (D.f[dMP0])[knw ]+
                     (D.f[dP0P])[kte ]+ (D.f[dM0M])[kbw ]+
                     (D.f[dP0M])[kbe ]+ (D.f[dM0P])[ktw ]+
                     (D.f[d0PP])[ktn ]+ (D.f[d0MM])[kbs ]+
                     (D.f[d0PM])[kbn ]+ (D.f[d0MP])[kts ]+
                     (D.f[d000])[kzero]+ 
                     (D.f[dPPP])[ktne]+ (D.f[dMMP])[ktsw]+ 
                     (D.f[dPMP])[ktse]+ (D.f[dMPP])[ktnw]+ 
                     (D.f[dPPM])[kbne]+ (D.f[dMMM])[kbsw]+ 
                     (D.f[dPMM])[kbse]+ (D.f[dMPM])[kbnw];
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cpPress[k] = PressCP;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}

void CalcCPtop27(real* DD, int* cpIndex, int nonCp, double* cpPress, unsigned int* neighborX, unsigned int* neighborY,
                 unsigned int* neighborZ, unsigned long long numberOfLBnodes, bool isEvenTimestep,
                 unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

    CalcCP27<<<grid.grid, grid.threads>>>(DD, cpIndex, nonCp, cpPress, neighborX, neighborY, neighborZ, numberOfLBnodes,
                                          isEvenTimestep);
    getLastCudaError("CalcCP27 execution failed");
}

void CalcCPbottom27(real* DD, int* cpIndex, int nonCp, double* cpPress, unsigned int* neighborX, unsigned int* neighborY,
                    unsigned int* neighborZ, unsigned long long numberOfLBnodes, bool isEvenTimestep,
                    unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, nonCp);

    CalcCP27<<<grid.grid, grid.threads>>>(DD, cpIndex, nonCp, cpPress, neighborX, neighborY, neighborZ, numberOfLBnodes,
                                          isEvenTimestep);
    getLastCudaError("CalcCP27 execution failed");
}

//! \}
