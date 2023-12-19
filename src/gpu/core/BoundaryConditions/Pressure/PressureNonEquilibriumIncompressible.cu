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


__global__ void PressureNonEquilibriumIncompressible_Device(
    real* rhoBC,
    real* DD,
    int* k_Q,
    int* k_N,
    int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep,
    size_t direction)
{
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
      //index
      unsigned int KQK  = k_Q[k];
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
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
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
      Distributions27 D;
      if (isEvenTimestep==true) //// ACHTUNG PREColl !!!!!!!!!!!!!!
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
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dP00])[k1e   ];
      f1_E    = (D.f[dM00])[k1w   ];
      f1_S    = (D.f[d0P0])[k1n   ];
      f1_N    = (D.f[d0M0])[k1s   ];
      f1_B    = (D.f[d00P])[k1t   ];
      f1_T    = (D.f[d00M])[k1b   ];
      f1_SW   = (D.f[dPP0])[k1ne  ];
      f1_NE   = (D.f[dMM0])[k1sw  ];
      f1_NW   = (D.f[dPM0])[k1se  ];
      f1_SE   = (D.f[dMP0])[k1nw  ];
      f1_BW   = (D.f[dP0P])[k1te  ];
      f1_TE   = (D.f[dM0M])[k1bw  ];
      f1_TW   = (D.f[dP0M])[k1be  ];
      f1_BE   = (D.f[dM0P])[k1tw  ];
      f1_BS   = (D.f[d0PP])[k1tn  ];
      f1_TN   = (D.f[d0MM])[k1bs  ];
      f1_TS   = (D.f[d0PM])[k1bn  ];
      f1_BN   = (D.f[d0MP])[k1ts  ];
      f1_ZERO = (D.f[d000])[k1zero];
      f1_BSW  = (D.f[dPPP])[k1tne ];
      f1_BNE  = (D.f[dMMP])[k1tsw ];
      f1_BNW  = (D.f[dPMP])[k1tse ];
      f1_BSE  = (D.f[dMPP])[k1tnw ];
      f1_TSW  = (D.f[dPPM])[k1bne ];
      f1_TNE  = (D.f[dMMM])[k1bsw ];
      f1_TNW  = (D.f[dPMM])[k1bse ];
      f1_TSE  = (D.f[dMPM])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      real vx1      =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                    ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                    (f1_E - f1_W);


      real vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                   (f1_N - f1_S);

      real vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                   (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                   (f1_T - f1_B);

      real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      f1_ZERO  -= c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

     drho1 = rhoBC[k];

      f1_ZERO  += c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));


      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes (only for the relevant directions)
      //!

      switch (direction)
      {
         case dM00:
            (D.f[dP00])[ke   ] = f1_W   ;
            (D.f[d0P0])[kn   ] = f1_S   ;
            (D.f[d0M0])[ks   ] = f1_N   ;
            (D.f[d00P])[kt   ] = f1_B   ;
            (D.f[d00M])[kb   ] = f1_T   ;
            (D.f[dPP0])[kne  ] = f1_SW  ;
            (D.f[dPM0])[kse  ] = f1_NW  ;
            (D.f[dP0P])[kte  ] = f1_BW  ;
            (D.f[dP0M])[kbe  ] = f1_TW  ;
            (D.f[d0PP])[ktn  ] = f1_BS  ;
            (D.f[d0MM])[kbs  ] = f1_TN  ;
            (D.f[d0PM])[kbn  ] = f1_TS  ;
            (D.f[d0MP])[kts  ] = f1_BN  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dPPP])[ktne ] = f1_BSW ;
            (D.f[dPMP])[ktse ] = f1_BNW ;
            (D.f[dPPM])[kbne ] = f1_TSW ;
            (D.f[dPMM])[kbse ] = f1_TNW ;
            break;
         case dP00:
            (D.f[dM00])[kw   ] = f1_E   ;
            (D.f[d0P0])[kn   ] = f1_S   ;
            (D.f[d0M0])[ks   ] = f1_N   ;
            (D.f[d00P])[kt   ] = f1_B   ;
            (D.f[d00M])[kb   ] = f1_T   ;
            (D.f[dMM0])[ksw  ] = f1_NE  ;
            (D.f[dMP0])[knw  ] = f1_SE  ;
            (D.f[dM0M])[kbw  ] = f1_TE  ;
            (D.f[dM0P])[ktw  ] = f1_BE  ;
            (D.f[d0PP])[ktn  ] = f1_BS  ;
            (D.f[d0MM])[kbs  ] = f1_TN  ;
            (D.f[d0PM])[kbn  ] = f1_TS  ;
            (D.f[d0MP])[kts  ] = f1_BN  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dMMP])[ktsw ] = f1_BNE ;
            (D.f[dMPP])[ktnw ] = f1_BSE ;
            (D.f[dMMM])[kbsw ] = f1_TNE ;
            (D.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         case d0M0:
            (D.f[dP00])[ke   ] = f1_W   ;
            (D.f[dM00])[kw   ] = f1_E   ;
            (D.f[d0P0])[kn   ] = f1_S   ;
            (D.f[d00P])[kt   ] = f1_B   ;
            (D.f[d00M])[kb   ] = f1_T   ;
            (D.f[dPP0])[kne  ] = f1_SW  ;
            (D.f[dMP0])[knw  ] = f1_SE  ;
            (D.f[dP0P])[kte  ] = f1_BW  ;
            (D.f[dM0M])[kbw  ] = f1_TE  ;
            (D.f[dP0M])[kbe  ] = f1_TW  ;
            (D.f[dM0P])[ktw  ] = f1_BE  ;
            (D.f[d0PP])[ktn  ] = f1_BS  ;
            (D.f[d0PM])[kbn  ] = f1_TS  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dPPP])[ktne ] = f1_BSW ;
            (D.f[dMPP])[ktnw ] = f1_BSE ;
            (D.f[dPPM])[kbne ] = f1_TSW ;
            (D.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         case d0P0:
            (D.f[dP00])[ke   ] = f1_W   ;
            (D.f[dM00])[kw   ] = f1_E   ;
            (D.f[d0M0])[ks   ] = f1_N   ;
            (D.f[d00P])[kt   ] = f1_B   ;
            (D.f[d00M])[kb   ] = f1_T   ;
            (D.f[dMM0])[ksw  ] = f1_NE  ;
            (D.f[dPM0])[kse  ] = f1_NW  ;
            (D.f[dP0P])[kte  ] = f1_BW  ;
            (D.f[dM0M])[kbw  ] = f1_TE  ;
            (D.f[dP0M])[kbe  ] = f1_TW  ;
            (D.f[dM0P])[ktw  ] = f1_BE  ;
            (D.f[d0MM])[kbs  ] = f1_TN  ;
            (D.f[d0MP])[kts  ] = f1_BN  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dMMP])[ktsw ] = f1_BNE ;
            (D.f[dPMP])[ktse ] = f1_BNW ;
            (D.f[dMMM])[kbsw ] = f1_TNE ;
            (D.f[dPMM])[kbse ] = f1_TNW ;
            break;
         case d00M:
            (D.f[dP00])[ke   ] = f1_W   ;
            (D.f[dM00])[kw   ] = f1_E   ;
            (D.f[d0P0])[kn   ] = f1_S   ;
            (D.f[d0M0])[ks   ] = f1_N   ;
            (D.f[d00P])[kt   ] = f1_B   ;
            (D.f[dPP0])[kne  ] = f1_SW  ;
            (D.f[dMM0])[ksw  ] = f1_NE  ;
            (D.f[dPM0])[kse  ] = f1_NW  ;
            (D.f[dMP0])[knw  ] = f1_SE  ;
            (D.f[dP0P])[kte  ] = f1_BW  ;
            (D.f[dM0P])[ktw  ] = f1_BE  ;
            (D.f[d0PP])[ktn  ] = f1_BS  ;
            (D.f[d0MP])[kts  ] = f1_BN  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dPPP])[ktne ] = f1_BSW ;
            (D.f[dMMP])[ktsw ] = f1_BNE ;
            (D.f[dPMP])[ktse ] = f1_BNW ;
            (D.f[dMPP])[ktnw ] = f1_BSE ;
            break;
         case d00P:
            (D.f[dP00])[ke   ] = f1_W   ;
            (D.f[dM00])[kw   ] = f1_E   ;
            (D.f[d0P0])[kn   ] = f1_S   ;
            (D.f[d0M0])[ks   ] = f1_N   ;
            (D.f[d00M])[kb   ] = f1_T   ;
            (D.f[dPP0])[kne  ] = f1_SW  ;
            (D.f[dMM0])[ksw  ] = f1_NE  ;
            (D.f[dPM0])[kse  ] = f1_NW  ;
            (D.f[dMP0])[knw  ] = f1_SE  ;
            (D.f[dM0M])[kbw  ] = f1_TE  ;
            (D.f[dP0M])[kbe  ] = f1_TW  ;
            (D.f[d0MM])[kbs  ] = f1_TN  ;
            (D.f[d0PM])[kbn  ] = f1_TS  ;
            (D.f[d000])[kzero] = f1_ZERO;
            (D.f[dPPM])[kbne ] = f1_TSW ;
            (D.f[dMMM])[kbsw ] = f1_TNE ;
            (D.f[dPMM])[kbse ] = f1_TNW ;
            (D.f[dMPM])[kbnw ] = f1_TSE ;
            break;
         default:
            break; 
      }
   }
}


//! \}
