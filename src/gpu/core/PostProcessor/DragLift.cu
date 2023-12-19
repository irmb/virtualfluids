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

#include "DragLift.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "Calculation/Calculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__ void DragLiftPost27(  real* DD, 
                                            int* k_Q, 
                                            real* QQ,
                                            int numberOfBCnodes, 
                                            double *DragX,
                                            double *DragY,
                                            double *DragZ,
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

    if(k<numberOfBCnodes)
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
        //unsigned int kzero= KQK;
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

        f_W    = (D.f[dP00])[ke   ];
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
        double    OnE   = c0o1, OnW   = c0o1, OnN   = c0o1, OnS   = c0o1, OnT = c0o1, OnB = c0o1, 
                OnNE  = c0o1, OnSW  = c0o1, OnSE  = c0o1, OnNW  = c0o1, 
                OnTE  = c0o1, OnBW  = c0o1, OnBE  = c0o1, OnTW  = c0o1,
                OnTN  = c0o1, OnBS  = c0o1, OnBN  = c0o1, OnTS  = c0o1, 
                OnTNE = c0o1, OnTSW = c0o1, OnTSE = c0o1, OnTNW = c0o1, 
                OnBNE = c0o1, OnBSW = c0o1, OnBSE = c0o1, OnBNW = c0o1;
        ////////////////////////////////////////////////////////////////////////////////
        real q;
        q = q_dirE[k];        if (q>=c0o1 && q<=c1o1) OnE   = c1o1;
        q = q_dirW[k];        if (q>=c0o1 && q<=c1o1) OnW   = c1o1;
        q = q_dirN[k];        if (q>=c0o1 && q<=c1o1) OnN   = c1o1;
        q = q_dirS[k];        if (q>=c0o1 && q<=c1o1) OnS   = c1o1;
        q = q_dirT[k];        if (q>=c0o1 && q<=c1o1) OnT   = c1o1;
        q = q_dirB[k];        if (q>=c0o1 && q<=c1o1) OnB   = c1o1;
        q = q_dirNE[k];        if (q>=c0o1 && q<=c1o1) OnNE  = c1o1;
        q = q_dirSW[k];        if (q>=c0o1 && q<=c1o1) OnSW  = c1o1;
        q = q_dirSE[k];        if (q>=c0o1 && q<=c1o1) OnSE  = c1o1;
        q = q_dirNW[k];        if (q>=c0o1 && q<=c1o1) OnNW  = c1o1;
        q = q_dirTE[k];        if (q>=c0o1 && q<=c1o1) OnTE  = c1o1;
        q = q_dirBW[k];        if (q>=c0o1 && q<=c1o1) OnBW  = c1o1;
        q = q_dirBE[k];        if (q>=c0o1 && q<=c1o1) OnBE  = c1o1;
        q = q_dirTW[k];        if (q>=c0o1 && q<=c1o1) OnTW  = c1o1;
        q = q_dirTN[k];        if (q>=c0o1 && q<=c1o1) OnTN  = c1o1;
        q = q_dirBS[k];        if (q>=c0o1 && q<=c1o1) OnBS  = c1o1;
        q = q_dirBN[k];        if (q>=c0o1 && q<=c1o1) OnBN  = c1o1;
        q = q_dirTS[k];        if (q>=c0o1 && q<=c1o1) OnTS  = c1o1;
        q = q_dirTNE[k];    if (q>=c0o1 && q<=c1o1) OnTNE = c1o1;
        q = q_dirBSW[k];    if (q>=c0o1 && q<=c1o1) OnBSW = c1o1;
        q = q_dirBNE[k];    if (q>=c0o1 && q<=c1o1) OnBNE = c1o1;
        q = q_dirTSW[k];    if (q>=c0o1 && q<=c1o1) OnTSW = c1o1;
        q = q_dirTSE[k];    if (q>=c0o1 && q<=c1o1) OnTSE = c1o1;
        q = q_dirBNW[k];    if (q>=c0o1 && q<=c1o1) OnBNW = c1o1;
        q = q_dirBSE[k];    if (q>=c0o1 && q<=c1o1) OnBSE = c1o1;
        q = q_dirTNW[k];    if (q>=c0o1 && q<=c1o1) OnTNW = c1o1;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double dragX, dragY, dragZ;

        dragX = ((((f_BSE * OnBSE) - (f_TNW * OnTNW))   + 
                  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
                  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_BE  * OnBE ) - (f_TW  * OnTW ))   + 
                  ((f_TE  * OnTE ) - (f_BW  * OnBW ))   + 
                  ((f_SE  * OnSE ) - (f_NW  * OnNW ))   + 
                  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
                  ((f_E   * OnE  ) - (f_W   * OnW  ));

        dragY = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
                  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
                  ((f_BNW * OnBNW) - (f_TSE * OnTSE))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_BN  * OnBN ) - (f_TS  * OnTS ))   + 
                  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
                  ((f_NW  * OnNW ) - (f_SE  * OnSE ))   + 
                  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
                  ((f_N   * OnN  ) - (f_S   * OnS  ));

        dragZ = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
                  ((f_TSW * OnTSW) - (f_BNE * OnBNE))   + 
                  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_TS  * OnTS ) - (f_BN  * OnBN ))   + 
                  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
                  ((f_TW  * OnTW ) - (f_BE  * OnBE ))   + 
                  ((f_TE  * OnTE ) - (f_BW  * OnBW )))) + 
                  ((f_T   * OnT  ) - (f_B   * OnB  ));
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        DragX[k] = -dragX;
        DragY[k] = -dragY;
        DragZ[k] = -dragZ;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}

__global__ void DragLiftPre27(   real* DD, 
                                            int* k_Q, 
                                            real* QQ,
                                            int numberOfBCnodes, 
                                            double *DragX,
                                            double *DragY,
                                            double *DragZ,
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

    if(k<numberOfBCnodes)
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
        //unsigned int kzero= KQK;
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

        f_E   = (D.f[dP00])[ke   ];
        f_W   = (D.f[dM00])[kw   ];
        f_N   = (D.f[d0P0])[kn   ];
        f_S   = (D.f[d0M0])[ks   ];
        f_T   = (D.f[d00P])[kt   ];
        f_B   = (D.f[d00M])[kb   ];
        f_NE  = (D.f[dPP0])[kne  ];
        f_SW  = (D.f[dMM0])[ksw  ];
        f_SE  = (D.f[dPM0])[kse  ];
        f_NW  = (D.f[dMP0])[knw  ];
        f_TE  = (D.f[dP0P])[kte  ];
        f_BW  = (D.f[dM0M])[kbw  ];
        f_BE  = (D.f[dP0M])[kbe  ];
        f_TW  = (D.f[dM0P])[ktw  ];
        f_TN  = (D.f[d0PP])[ktn  ];
        f_BS  = (D.f[d0MM])[kbs  ];
        f_BN  = (D.f[d0PM])[kbn  ];
        f_TS  = (D.f[d0MP])[kts  ];
        f_TNE = (D.f[dPPP])[ktne ];
        f_TSW = (D.f[dMMP])[ktsw ];
        f_TSE = (D.f[dPMP])[ktse ];
        f_TNW = (D.f[dMPP])[ktnw ];
        f_BNE = (D.f[dPPM])[kbne ];
        f_BSW = (D.f[dMMM])[kbsw ];
        f_BSE = (D.f[dPMM])[kbse ];
        f_BNW = (D.f[dMPM])[kbnw ];
         ////////////////////////////////////////////////////////////////////////////////
        double    OnE   = c0o1, OnW   = c0o1, OnN   = c0o1, OnS   = c0o1, OnT = c0o1, OnB = c0o1, 
                OnNE  = c0o1, OnSW  = c0o1, OnSE  = c0o1, OnNW  = c0o1, 
                OnTE  = c0o1, OnBW  = c0o1, OnBE  = c0o1, OnTW  = c0o1,
                OnTN  = c0o1, OnBS  = c0o1, OnBN  = c0o1, OnTS  = c0o1, 
                OnTNE = c0o1, OnTSW = c0o1, OnTSE = c0o1, OnTNW = c0o1, 
                OnBNE = c0o1, OnBSW = c0o1, OnBSE = c0o1, OnBNW = c0o1;
        ////////////////////////////////////////////////////////////////////////////////
        real q;
        q = q_dirE[k];        if (q>=c0o1 && q<=c1o1) OnW   = c1o1;
        q = q_dirW[k];        if (q>=c0o1 && q<=c1o1) OnE   = c1o1;
        q = q_dirN[k];        if (q>=c0o1 && q<=c1o1) OnS   = c1o1;
        q = q_dirS[k];        if (q>=c0o1 && q<=c1o1) OnN   = c1o1;
        q = q_dirT[k];        if (q>=c0o1 && q<=c1o1) OnB   = c1o1;
        q = q_dirB[k];        if (q>=c0o1 && q<=c1o1) OnT   = c1o1;
        q = q_dirNE[k];        if (q>=c0o1 && q<=c1o1) OnSW  = c1o1;
        q = q_dirSW[k];        if (q>=c0o1 && q<=c1o1) OnNE  = c1o1;
        q = q_dirSE[k];        if (q>=c0o1 && q<=c1o1) OnNW  = c1o1;
        q = q_dirNW[k];        if (q>=c0o1 && q<=c1o1) OnSE  = c1o1;
        q = q_dirTE[k];        if (q>=c0o1 && q<=c1o1) OnBW  = c1o1;
        q = q_dirBW[k];        if (q>=c0o1 && q<=c1o1) OnTE  = c1o1;
        q = q_dirBE[k];        if (q>=c0o1 && q<=c1o1) OnTW  = c1o1;
        q = q_dirTW[k];        if (q>=c0o1 && q<=c1o1) OnBE  = c1o1;
        q = q_dirTN[k];        if (q>=c0o1 && q<=c1o1) OnBS  = c1o1;
        q = q_dirBS[k];        if (q>=c0o1 && q<=c1o1) OnTN  = c1o1;
        q = q_dirBN[k];        if (q>=c0o1 && q<=c1o1) OnTS  = c1o1;
        q = q_dirTS[k];        if (q>=c0o1 && q<=c1o1) OnBN  = c1o1;
        q = q_dirTNE[k];    if (q>=c0o1 && q<=c1o1) OnBSW = c1o1;
        q = q_dirBSW[k];    if (q>=c0o1 && q<=c1o1) OnTNE = c1o1;
        q = q_dirBNE[k];    if (q>=c0o1 && q<=c1o1) OnTSW = c1o1;
        q = q_dirTSW[k];    if (q>=c0o1 && q<=c1o1) OnBNE = c1o1;
        q = q_dirTSE[k];    if (q>=c0o1 && q<=c1o1) OnBNW = c1o1;
        q = q_dirBNW[k];    if (q>=c0o1 && q<=c1o1) OnTSE = c1o1;
        q = q_dirBSE[k];    if (q>=c0o1 && q<=c1o1) OnTNW = c1o1;
        q = q_dirTNW[k];    if (q>=c0o1 && q<=c1o1) OnBSE = c1o1;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double dragX, dragY, dragZ;

        dragX = ((((f_BSE * OnBSE) - (f_TNW * OnTNW))   + 
                  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
                  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_BE  * OnBE ) - (f_TW  * OnTW ))   + 
                  ((f_TE  * OnTE ) - (f_BW  * OnBW ))   + 
                  ((f_SE  * OnSE ) - (f_NW  * OnNW ))   + 
                  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
                  ((f_E   * OnE  ) - (f_W   * OnW  ));

        dragY = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
                  ((f_BNE * OnBNE) - (f_TSW * OnTSW))   + 
                  ((f_BNW * OnBNW) - (f_TSE * OnTSE))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_BN  * OnBN ) - (f_TS  * OnTS ))   + 
                  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
                  ((f_NW  * OnNW ) - (f_SE  * OnSE ))   + 
                  ((f_NE  * OnNE ) - (f_SW  * OnSW )))) + 
                  ((f_N   * OnN  ) - (f_S   * OnS  ));

        dragZ = ((((f_TNW * OnTNW) - (f_BSE * OnBSE))   + 
                  ((f_TSW * OnTSW) - (f_BNE * OnBNE))   + 
                  ((f_TSE * OnTSE) - (f_BNW * OnBNW))   + 
                  ((f_TNE * OnTNE) - (f_BSW * OnBSW)))  + 
                 (((f_TS  * OnTS ) - (f_BN  * OnBN ))   + 
                  ((f_TN  * OnTN ) - (f_BS  * OnBS ))   + 
                  ((f_TW  * OnTW ) - (f_BE  * OnBE ))   + 
                  ((f_TE  * OnTE ) - (f_BW  * OnBW )))) + 
                  ((f_T   * OnT  ) - (f_B   * OnB  ));
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        DragX[k] = -dragX;
        DragY[k] = -dragY;
        DragZ[k] = -dragZ;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}

void DragLiftPostD27(real* DD, int* k_Q, real* QQ, int numberOfBCnodes, double* DragX, double* DragY, double* DragZ,
                     unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                     unsigned long long numberOfLBnodes, bool isEvenTimestep, unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    DragLiftPost27<<<grid.grid, grid.threads>>>(DD, k_Q, QQ, numberOfBCnodes, DragX, DragY, DragZ, neighborX, neighborY,
                                                neighborZ, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("DragLiftPost27 execution failed");
}

void DragLiftPreD27(real* DD, int* k_Q, real* QQ, int numberOfBCnodes, double* DragX, double* DragY, double* DragZ,
                    unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                    unsigned long long numberOfLBnodes, bool isEvenTimestep, unsigned int numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfBCnodes);

    DragLiftPre27<<<grid.grid, grid.threads>>>(DD, k_Q, QQ, numberOfBCnodes, DragX, DragY, DragZ, neighborX, neighborY,
                                               neighborZ, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("DragLiftPre27 execution failed");
}

//! \}
