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
//! \author Martin Schoenherr
//=======================================================================================
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"

#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__ void AdvectionDiffusionBounceBack_Device(
    real* DD, 
    real* DD27, 
    real* temp,
    real diffusivity,
    int* k_Q, 
    real* QQ,
    unsigned int numberOfBCnodes, 
    real om1, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[dP00] = &DD27[dP00 * numberOfLBnodes];
      D27.f[dM00] = &DD27[dM00 * numberOfLBnodes];
      D27.f[d0P0] = &DD27[d0P0 * numberOfLBnodes];
      D27.f[d0M0] = &DD27[d0M0 * numberOfLBnodes];
      D27.f[d00P] = &DD27[d00P * numberOfLBnodes];
      D27.f[d00M] = &DD27[d00M * numberOfLBnodes];
      D27.f[dPP0] = &DD27[dPP0 * numberOfLBnodes];
      D27.f[dMM0] = &DD27[dMM0 * numberOfLBnodes];
      D27.f[dPM0] = &DD27[dPM0 * numberOfLBnodes];
      D27.f[dMP0] = &DD27[dMP0 * numberOfLBnodes];
      D27.f[dP0P] = &DD27[dP0P * numberOfLBnodes];
      D27.f[dM0M] = &DD27[dM0M * numberOfLBnodes];
      D27.f[dP0M] = &DD27[dP0M * numberOfLBnodes];
      D27.f[dM0P] = &DD27[dM0P * numberOfLBnodes];
      D27.f[d0PP] = &DD27[d0PP * numberOfLBnodes];
      D27.f[d0MM] = &DD27[d0MM * numberOfLBnodes];
      D27.f[d0PM] = &DD27[d0PM * numberOfLBnodes];
      D27.f[d0MP] = &DD27[d0MP * numberOfLBnodes];
      D27.f[d000] = &DD27[d000 * numberOfLBnodes];
      D27.f[dPPP] = &DD27[dPPP * numberOfLBnodes];
      D27.f[dMMP] = &DD27[dMMP * numberOfLBnodes];
      D27.f[dPMP] = &DD27[dPMP * numberOfLBnodes];
      D27.f[dMPP] = &DD27[dMPP * numberOfLBnodes];
      D27.f[dPPM] = &DD27[dPPM * numberOfLBnodes];
      D27.f[dMMM] = &DD27[dMMM * numberOfLBnodes];
      D27.f[dPMM] = &DD27[dPMM * numberOfLBnodes];
      D27.f[dMPM] = &DD27[dMPM * numberOfLBnodes];
   } 
   else
   {
      D27.f[dM00] = &DD27[dP00 * numberOfLBnodes];
      D27.f[dP00] = &DD27[dM00 * numberOfLBnodes];
      D27.f[d0M0] = &DD27[d0P0 * numberOfLBnodes];
      D27.f[d0P0] = &DD27[d0M0 * numberOfLBnodes];
      D27.f[d00M] = &DD27[d00P * numberOfLBnodes];
      D27.f[d00P] = &DD27[d00M * numberOfLBnodes];
      D27.f[dMM0] = &DD27[dPP0 * numberOfLBnodes];
      D27.f[dPP0] = &DD27[dMM0 * numberOfLBnodes];
      D27.f[dMP0] = &DD27[dPM0 * numberOfLBnodes];
      D27.f[dPM0] = &DD27[dMP0 * numberOfLBnodes];
      D27.f[dM0M] = &DD27[dP0P * numberOfLBnodes];
      D27.f[dP0P] = &DD27[dM0M * numberOfLBnodes];
      D27.f[dM0P] = &DD27[dP0M * numberOfLBnodes];
      D27.f[dP0M] = &DD27[dM0P * numberOfLBnodes];
      D27.f[d0MM] = &DD27[d0PP * numberOfLBnodes];
      D27.f[d0PP] = &DD27[d0MM * numberOfLBnodes];
      D27.f[d0MP] = &DD27[d0PM * numberOfLBnodes];
      D27.f[d0PM] = &DD27[d0MP * numberOfLBnodes];
      D27.f[d000] = &DD27[d000 * numberOfLBnodes];
      D27.f[dPPP] = &DD27[dMMM * numberOfLBnodes];
      D27.f[dMMP] = &DD27[dPPM * numberOfLBnodes];
      D27.f[dPMP] = &DD27[dMPM * numberOfLBnodes];
      D27.f[dMPP] = &DD27[dPMM * numberOfLBnodes];
      D27.f[dPPM] = &DD27[dMMP * numberOfLBnodes];
      D27.f[dMMM] = &DD27[dPPP * numberOfLBnodes];
      D27.f[dPMM] = &DD27[dMPP * numberOfLBnodes];
      D27.f[dMPM] = &DD27[dPMP * numberOfLBnodes];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;
   const unsigned  y = blockIdx.x; 
   const unsigned  z = blockIdx.y; 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      real  *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
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
      real f27_W    = (D27.f[dP00])[ke   ];
      real f27_E    = (D27.f[dM00])[kw   ];
      real f27_S    = (D27.f[d0P0])[kn   ];
      real f27_N    = (D27.f[d0M0])[ks   ];
      real f27_B    = (D27.f[d00P])[kt   ];
      real f27_T    = (D27.f[d00M])[kb   ];
      real f27_SW   = (D27.f[dPP0])[kne  ];
      real f27_NE   = (D27.f[dMM0])[ksw  ];
      real f27_NW   = (D27.f[dPM0])[kse  ];
      real f27_SE   = (D27.f[dMP0])[knw  ];
      real f27_BW   = (D27.f[dP0P])[kte  ];
      real f27_TE   = (D27.f[dM0M])[kbw  ];
      real f27_TW   = (D27.f[dP0M])[kbe  ];
      real f27_BE   = (D27.f[dM0P])[ktw  ];
      real f27_BS   = (D27.f[d0PP])[ktn  ];
      real f27_TN   = (D27.f[d0MM])[kbs  ];
      real f27_TS   = (D27.f[d0PM])[kbn  ];
      real f27_BN   = (D27.f[d0MP])[kts  ];
      real f27_BSW  = (D27.f[dPPP])[ktne ];
      real f27_BNE  = (D27.f[dMMP])[ktsw ];
      real f27_BNW  = (D27.f[dPMP])[ktse ];
      real f27_BSE  = (D27.f[dMPP])[ktnw ];
      real f27_TSW  = (D27.f[dPPM])[kbne ];
      real f27_TNE  = (D27.f[dMMM])[kbsw ];
      real f27_TNW  = (D27.f[dPMM])[kbse ];
      real f27_TSE  = (D27.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
      {
         D27.f[dP00] = &DD27[dP00 * numberOfLBnodes];
         D27.f[dM00] = &DD27[dM00 * numberOfLBnodes];
         D27.f[d0P0] = &DD27[d0P0 * numberOfLBnodes];
         D27.f[d0M0] = &DD27[d0M0 * numberOfLBnodes];
         D27.f[d00P] = &DD27[d00P * numberOfLBnodes];
         D27.f[d00M] = &DD27[d00M * numberOfLBnodes];
         D27.f[dPP0] = &DD27[dPP0 * numberOfLBnodes];
         D27.f[dMM0] = &DD27[dMM0 * numberOfLBnodes];
         D27.f[dPM0] = &DD27[dPM0 * numberOfLBnodes];
         D27.f[dMP0] = &DD27[dMP0 * numberOfLBnodes];
         D27.f[dP0P] = &DD27[dP0P * numberOfLBnodes];
         D27.f[dM0M] = &DD27[dM0M * numberOfLBnodes];
         D27.f[dP0M] = &DD27[dP0M * numberOfLBnodes];
         D27.f[dM0P] = &DD27[dM0P * numberOfLBnodes];
         D27.f[d0PP] = &DD27[d0PP * numberOfLBnodes];
         D27.f[d0MM] = &DD27[d0MM * numberOfLBnodes];
         D27.f[d0PM] = &DD27[d0PM * numberOfLBnodes];
         D27.f[d0MP] = &DD27[d0MP * numberOfLBnodes];
         D27.f[d000] = &DD27[d000 * numberOfLBnodes];
         D27.f[dPPP] = &DD27[dPPP * numberOfLBnodes];
         D27.f[dMMP] = &DD27[dMMP * numberOfLBnodes];
         D27.f[dPMP] = &DD27[dPMP * numberOfLBnodes];
         D27.f[dMPP] = &DD27[dMPP * numberOfLBnodes];
         D27.f[dPPM] = &DD27[dPPM * numberOfLBnodes];
         D27.f[dMMM] = &DD27[dMMM * numberOfLBnodes];
         D27.f[dPMM] = &DD27[dPMM * numberOfLBnodes];
         D27.f[dMPM] = &DD27[dMPM * numberOfLBnodes];
      } 
      else
      {
         D27.f[dM00] = &DD27[dP00 * numberOfLBnodes];
         D27.f[dP00] = &DD27[dM00 * numberOfLBnodes];
         D27.f[d0M0] = &DD27[d0P0 * numberOfLBnodes];
         D27.f[d0P0] = &DD27[d0M0 * numberOfLBnodes];
         D27.f[d00M] = &DD27[d00P * numberOfLBnodes];
         D27.f[d00P] = &DD27[d00M * numberOfLBnodes];
         D27.f[dMM0] = &DD27[dPP0 * numberOfLBnodes];
         D27.f[dPP0] = &DD27[dMM0 * numberOfLBnodes];
         D27.f[dMP0] = &DD27[dPM0 * numberOfLBnodes];
         D27.f[dPM0] = &DD27[dMP0 * numberOfLBnodes];
         D27.f[dM0M] = &DD27[dP0P * numberOfLBnodes];
         D27.f[dP0P] = &DD27[dM0M * numberOfLBnodes];
         D27.f[dM0P] = &DD27[dP0M * numberOfLBnodes];
         D27.f[dP0M] = &DD27[dM0P * numberOfLBnodes];
         D27.f[d0MM] = &DD27[d0PP * numberOfLBnodes];
         D27.f[d0PP] = &DD27[d0MM * numberOfLBnodes];
         D27.f[d0MP] = &DD27[d0PM * numberOfLBnodes];
         D27.f[d0PM] = &DD27[d0MP * numberOfLBnodes];
         D27.f[d000] = &DD27[d000 * numberOfLBnodes];
         D27.f[dPPP] = &DD27[dMMM * numberOfLBnodes];
         D27.f[dMMP] = &DD27[dPPM * numberOfLBnodes];
         D27.f[dPMP] = &DD27[dMPM * numberOfLBnodes];
         D27.f[dMPP] = &DD27[dPMM * numberOfLBnodes];
         D27.f[dPPM] = &DD27[dMMP * numberOfLBnodes];
         D27.f[dMMM] = &DD27[dPPP * numberOfLBnodes];
         D27.f[dPMM] = &DD27[dMPP * numberOfLBnodes];
         D27.f[dMPM] = &DD27[dPMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real q;
      q = q_dirE[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dM00])[kw  ]=f27_E  ;
      q = q_dirW[k];   if (q>=c0o1 && q<=c1o1) (D27.f[dP00])[ke  ]=f27_W  ;
      q = q_dirN[k];   if (q>=c0o1 && q<=c1o1) (D27.f[d0M0])[ks  ]=f27_N  ;
      q = q_dirS[k];   if (q>=c0o1 && q<=c1o1) (D27.f[d0P0])[kn  ]=f27_S  ;
      q = q_dirT[k];   if (q>=c0o1 && q<=c1o1) (D27.f[d00M])[kb  ]=f27_T  ;
      q = q_dirB[k];   if (q>=c0o1 && q<=c1o1) (D27.f[d00P])[kt  ]=f27_B  ;
      q = q_dirNE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dMM0])[ksw ]=f27_NE ;
      q = q_dirSW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dPP0])[kne ]=f27_SW ;
      q = q_dirSE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dMP0])[knw ]=f27_SE ;
      q = q_dirNW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dPM0])[kse ]=f27_NW ;
      q = q_dirTE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dM0M])[kbw ]=f27_TE ;
      q = q_dirBW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dP0P])[kte ]=f27_BW ;
      q = q_dirBE[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dM0P])[ktw ]=f27_BE ;
      q = q_dirTW[k];  if (q>=c0o1 && q<=c1o1) (D27.f[dP0M])[kbe ]=f27_TW ;
      q = q_dirTN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[d0MM])[kbs ]=f27_TN ;
      q = q_dirBS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[d0PP])[ktn ]=f27_BS ;
      q = q_dirBN[k];  if (q>=c0o1 && q<=c1o1) (D27.f[d0MP])[kts ]=f27_BN ;
      q = q_dirTS[k];  if (q>=c0o1 && q<=c1o1) (D27.f[d0PM])[kbn ]=f27_TS ;
      q = q_dirTNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dMMM])[kbsw]=f27_TNE;
      q = q_dirBSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dPPP])[ktne]=f27_BSW;
      q = q_dirBNE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dMMP])[ktsw]=f27_BNE;
      q = q_dirTSW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dPPM])[kbne]=f27_TSW;
      q = q_dirTSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dMPM])[kbnw]=f27_TSE;
      q = q_dirBNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dPMP])[ktse]=f27_BNW;
      q = q_dirBSE[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dMPP])[ktnw]=f27_BSE;
      q = q_dirTNW[k]; if (q>=c0o1 && q<=c1o1) (D27.f[dPMM])[kbse]=f27_TNW;
   }
}

//! \}
