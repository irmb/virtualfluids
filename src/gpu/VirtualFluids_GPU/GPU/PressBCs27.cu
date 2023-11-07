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
//! \file PressBCs27.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "lbm/MacroscopicQuantities.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QInflowScaleByPressDevice27(
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
    bool isEvenTimestep)
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
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f1_E    = (D.f[dP00])[k1e   ];
      real f1_W    = (D.f[dM00])[k1w   ];
      real f1_N    = (D.f[d0P0])[k1n   ];
      real f1_S    = (D.f[d0M0])[k1s   ];
      real f1_T    = (D.f[d00P])[k1t   ];
      real f1_B    = (D.f[d00M])[k1b   ];
      real f1_NE   = (D.f[dPP0])[k1ne  ];
      real f1_SW   = (D.f[dMM0])[k1sw  ];
      real f1_SE   = (D.f[dPM0])[k1se  ];
      real f1_NW   = (D.f[dMP0])[k1nw  ];
      real f1_TE   = (D.f[dP0P])[k1te  ];
      real f1_BW   = (D.f[dM0M])[k1bw  ];
      real f1_BE   = (D.f[dP0M])[k1be  ];
      real f1_TW   = (D.f[dM0P])[k1tw  ];
      real f1_TN   = (D.f[d0PP])[k1tn  ];
      real f1_BS   = (D.f[d0MM])[k1bs  ];
      real f1_BN   = (D.f[d0PM])[k1bn  ];
      real f1_TS   = (D.f[d0MP])[k1ts  ];
      //real f1_ZERO = (D.f[d000])[k1zero];
      real f1_TNE  = (D.f[dPPP])[k1tne ];
      real f1_TSW  = (D.f[dMMP])[k1tsw ];
      real f1_TSE  = (D.f[dPMP])[k1tse ];
      real f1_TNW  = (D.f[dMPP])[k1tnw ];
      real f1_BNE  = (D.f[dPPM])[k1bne ];
      real f1_BSW  = (D.f[dMMM])[k1bsw ];
      real f1_BSE  = (D.f[dPMM])[k1bse ];
      real f1_BNW  = (D.f[dMPM])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E    = (D.f[dP00])[ke   ];
      real f_W    = (D.f[dM00])[kw   ];
      real f_N    = (D.f[d0P0])[kn   ];
      real f_S    = (D.f[d0M0])[ks   ];
      real f_T    = (D.f[d00P])[kt   ];
      real f_B    = (D.f[d00M])[kb   ];
      real f_NE   = (D.f[dPP0])[kne  ];
      real f_SW   = (D.f[dMM0])[ksw  ];
      real f_SE   = (D.f[dPM0])[kse  ];
      real f_NW   = (D.f[dMP0])[knw  ];
      real f_TE   = (D.f[dP0P])[kte  ];
      real f_BW   = (D.f[dM0M])[kbw  ];
      real f_BE   = (D.f[dP0M])[kbe  ];
      real f_TW   = (D.f[dM0P])[ktw  ];
      real f_TN   = (D.f[d0PP])[ktn  ];
      real f_BS   = (D.f[d0MM])[kbs  ];
      real f_BN   = (D.f[d0PM])[kbn  ];
      real f_TS   = (D.f[d0MP])[kts  ];
      //real f_ZERO = (D.f[d000])[kzero];
      real f_TNE  = (D.f[dPPP])[ktne ];
      real f_TSW  = (D.f[dMMP])[ktsw ];
      real f_TSE  = (D.f[dPMP])[ktse ];
      real f_TNW  = (D.f[dMPP])[ktnw ];
      real f_BNE  = (D.f[dPPM])[kbne ];
      real f_BSW  = (D.f[dMMM])[kbsw ];
      real f_BSE  = (D.f[dPMM])[kbse ];
      real f_BNW  = (D.f[dMPM])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      // real vx1, vx2, vx3;
      real drho, drho1;
      //////////////////////////////////////////////////////////////////////////
     //Dichte
      drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW +
                f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[d000])[k1zero]);
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);
      //////////////////////////////////////////////////////////////////////////
     //Schallgeschwindigkeit
     real cs = c1o1 / sqrtf(c3o1);
      //////////////////////////////////////////////////////////////////////////
     real rhoInterpol = drho1 * cs + (c1o1 - cs) * drho;
     //real diffRho = (rhoBC[k] + one) / (rhoInterpol + one);
     real diffRhoToAdd = rhoBC[k] - rhoInterpol;
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //no velocity
     //////////////////////////////////////////
      f_E    = f1_E   * cs + (c1o1 - cs) * f_E   ;
      f_W    = f1_W   * cs + (c1o1 - cs) * f_W   ;
      f_N    = f1_N   * cs + (c1o1 - cs) * f_N   ;
      f_S    = f1_S   * cs + (c1o1 - cs) * f_S   ;
      f_T    = f1_T   * cs + (c1o1 - cs) * f_T   ;
      f_B    = f1_B   * cs + (c1o1 - cs) * f_B   ;
      f_NE   = f1_NE  * cs + (c1o1 - cs) * f_NE  ;
      f_SW   = f1_SW  * cs + (c1o1 - cs) * f_SW  ;
      f_SE   = f1_SE  * cs + (c1o1 - cs) * f_SE  ;
      f_NW   = f1_NW  * cs + (c1o1 - cs) * f_NW  ;
      f_TE   = f1_TE  * cs + (c1o1 - cs) * f_TE  ;
      f_BW   = f1_BW  * cs + (c1o1 - cs) * f_BW  ;
      f_BE   = f1_BE  * cs + (c1o1 - cs) * f_BE  ;
      f_TW   = f1_TW  * cs + (c1o1 - cs) * f_TW  ;
      f_TN   = f1_TN  * cs + (c1o1 - cs) * f_TN  ;
      f_BS   = f1_BS  * cs + (c1o1 - cs) * f_BS  ;
      f_BN   = f1_BN  * cs + (c1o1 - cs) * f_BN  ;
      f_TS   = f1_TS  * cs + (c1o1 - cs) * f_TS  ;
      f_TNE  = f1_TNE * cs + (c1o1 - cs) * f_TNE ;
      f_TSW  = f1_TSW * cs + (c1o1 - cs) * f_TSW ;
      f_TSE  = f1_TSE * cs + (c1o1 - cs) * f_TSE ;
      f_TNW  = f1_TNW * cs + (c1o1 - cs) * f_TNW ;
      f_BNE  = f1_BNE * cs + (c1o1 - cs) * f_BNE ;
      f_BSW  = f1_BSW * cs + (c1o1 - cs) * f_BSW ;
      f_BSE  = f1_BSE * cs + (c1o1 - cs) * f_BSE ;
      f_BNW  = f1_BNW * cs + (c1o1 - cs) * f_BNW ;
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //scale by press
     //////////////////////////////////////////
     //f_E    = (f_E   + c2over27 ) * diffRho - c2over27 ;
   //   f_W    = (f_W   + c2over27 ) * diffRho - c2over27 ;
   //   f_N    = (f_N   + c2over27 ) * diffRho - c2over27 ;
   //   f_S    = (f_S   + c2over27 ) * diffRho - c2over27 ;
   //   f_T    = (f_T   + c2over27 ) * diffRho - c2over27 ;
   //   f_B    = (f_B   + c2over27 ) * diffRho - c2over27 ;
     //f_NE   = (f_NE  + c1over54 ) * diffRho - c1over54 ;
   //   f_SW   = (f_SW  + c1over54 ) * diffRho - c1over54 ;
   //   f_SE   = (f_SE  + c1over54 ) * diffRho - c1over54 ;
   //   f_NW   = (f_NW  + c1over54 ) * diffRho - c1over54 ;
   //   f_TE   = (f_TE  + c1over54 ) * diffRho - c1over54 ;
   //   f_BW   = (f_BW  + c1over54 ) * diffRho - c1over54 ;
   //   f_BE   = (f_BE  + c1over54 ) * diffRho - c1over54 ;
   //   f_TW   = (f_TW  + c1over54 ) * diffRho - c1over54 ;
   //   f_TN   = (f_TN  + c1over54 ) * diffRho - c1over54 ;
   //   f_BS   = (f_BS  + c1over54 ) * diffRho - c1over54 ;
   //   f_BN   = (f_BN  + c1over54 ) * diffRho - c1over54 ;
   //   f_TS   = (f_TS  + c1over54 ) * diffRho - c1over54 ;
   //   f_TNE  = (f_TNE + c1over216) * diffRho - c1over216;
   //   f_TSW  = (f_TSW + c1over216) * diffRho - c1over216;
   //   f_TSE  = (f_TSE + c1over216) * diffRho - c1over216;
   //   f_TNW  = (f_TNW + c1over216) * diffRho - c1over216;
   //   f_BNE  = (f_BNE + c1over216) * diffRho - c1over216;
   //   f_BSW  = (f_BSW + c1over216) * diffRho - c1over216;
   //   f_BSE  = (f_BSE + c1over216) * diffRho - c1over216;
   //   f_BNW  = (f_BNW + c1over216) * diffRho - c1over216;
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // add press
     //////////////////////////////////////////
     f_E    = (f_E   + c2o27  * diffRhoToAdd);
      f_W    = (f_W   + c2o27  * diffRhoToAdd);
      f_N    = (f_N   + c2o27  * diffRhoToAdd);
      f_S    = (f_S   + c2o27  * diffRhoToAdd);
      f_T    = (f_T   + c2o27  * diffRhoToAdd);
      f_B    = (f_B   + c2o27  * diffRhoToAdd);
     f_NE   = (f_NE  + c1o54  * diffRhoToAdd);
      f_SW   = (f_SW  + c1o54  * diffRhoToAdd);
      f_SE   = (f_SE  + c1o54  * diffRhoToAdd);
      f_NW   = (f_NW  + c1o54  * diffRhoToAdd);
      f_TE   = (f_TE  + c1o54  * diffRhoToAdd);
      f_BW   = (f_BW  + c1o54  * diffRhoToAdd);
      f_BE   = (f_BE  + c1o54  * diffRhoToAdd);
      f_TW   = (f_TW  + c1o54  * diffRhoToAdd);
      f_TN   = (f_TN  + c1o54  * diffRhoToAdd);
      f_BS   = (f_BS  + c1o54  * diffRhoToAdd);
      f_BN   = (f_BN  + c1o54  * diffRhoToAdd);
      f_TS   = (f_TS  + c1o54  * diffRhoToAdd);
      f_TNE  = (f_TNE + c1o216 * diffRhoToAdd);
      f_TSW  = (f_TSW + c1o216 * diffRhoToAdd);
      f_TSE  = (f_TSE + c1o216 * diffRhoToAdd);
      f_TNW  = (f_TNW + c1o216 * diffRhoToAdd);
      f_BNE  = (f_BNE + c1o216 * diffRhoToAdd);
      f_BSW  = (f_BSW + c1o216 * diffRhoToAdd);
      f_BSE  = (f_BSE + c1o216 * diffRhoToAdd);
      f_BNW  = (f_BNW + c1o216 * diffRhoToAdd);
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
     // -X
     //(D.f[dP00])[ke   ] = f_E   ;
     //(D.f[dPM0])[kse  ] = f_SE  ;
     //(D.f[dPP0])[kne  ] = f_NE  ;
     //(D.f[dP0M])[kbe  ] = f_BE  ;
     //(D.f[dP0P])[kte  ] = f_TE  ;
     //(D.f[dPMP])[ktse ] = f_TSE ;
     //(D.f[dPPP])[ktne ] = f_TNE ;
     //(D.f[dPMM])[kbse ] = f_BSE ;
     //(D.f[dPPM])[kbne ] = f_BNE ;
     // X
     (D.f[dM00])[kw   ] = f_W   ;
     (D.f[dMM0])[ksw  ] = f_SW  ;
     (D.f[dMP0])[knw  ] = f_NW  ;
     (D.f[dM0M])[kbw  ] = f_BW  ;
     (D.f[dM0P])[ktw  ] = f_TW  ;
     (D.f[dMMP])[ktsw ] = f_TSW ;
     (D.f[dMPP])[ktnw ] = f_TNW ;
     (D.f[dMMM])[kbsw ] = f_BSW ;
     (D.f[dMPM])[kbnw ] = f_BNW ;
     // Y
     //(D.f[d0M0])[ks   ] = f_S   ;
     //(D.f[dPM0])[kse  ] = f_SE  ;
     //(D.f[dMM0])[ksw  ] = f_SW  ;
     //(D.f[d0MP])[kts  ] = f_TS  ;
     //(D.f[d0MM])[kbs  ] = f_BS  ;
     //(D.f[dPMP])[ktse ] = f_TSE ;
     //(D.f[dMMP])[ktsw ] = f_TSW ;
     //(D.f[dPMM])[kbse ] = f_BSE ;
     //(D.f[dMMM])[kbsw ] = f_BSW ;
     // Z
     //(D.f[d00M])[kb   ] = f_B   ;
     //(D.f[dP0M])[kbe  ] = f_BE  ;
     //(D.f[dM0M])[kbw  ] = f_BW  ;
     //(D.f[d0PM])[kbn  ] = f_BN  ;
     //(D.f[d0MM])[kbs  ] = f_BS  ;
     //(D.f[dPPM])[kbne ] = f_BNE ;
     //(D.f[dMPM])[kbnw ] = f_BNW ;
     //(D.f[dPMM])[kbse ] = f_BSE ;
     //(D.f[dMMM])[kbsw ] = f_BSW ;
      //////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceIncompNEQ27(
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
    bool isEvenTimestep)
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

     //if(vx1 < zero){
       // vx1 *= 0.9;
     //}
     //if(vx2 < zero){
       // vx2 *= c1o10;//0.9;
     //}

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

     //drho1 = (drho1 + rhoBC[k])/2.f;
     //drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[dP00])[ke   ] = f1_W   ;
      (D.f[dM00])[kw   ] = f1_E   ;
      (D.f[d0P0])[kn   ] = f1_S   ;
      (D.f[d0M0])[ks   ] = f1_N   ;
      (D.f[d00P])[kt   ] = f1_B   ;
      (D.f[d00M])[kb   ] = f1_T   ;
      (D.f[dPP0])[kne  ] = f1_SW  ;
      (D.f[dMM0])[ksw  ] = f1_NE  ;
      (D.f[dPM0])[kse  ] = f1_NW  ;
      (D.f[dMP0])[knw  ] = f1_SE  ;
      (D.f[dP0P])[kte  ] = f1_BW  ;
      (D.f[dM0M])[kbw  ] = f1_TE  ;
      (D.f[dP0M])[kbe  ] = f1_TW  ;
      (D.f[dM0P])[ktw  ] = f1_BE  ;
      (D.f[d0PP])[ktn  ] = f1_BS  ;
      (D.f[d0MM])[kbs  ] = f1_TN  ;
      (D.f[d0PM])[kbn  ] = f1_TS  ;
      (D.f[d0MP])[kts  ] = f1_BN  ;
      (D.f[d000])[kzero] = f1_ZERO;
      (D.f[dPPP])[ktne ] = f1_BSW ;
      (D.f[dMMP])[ktsw ] = f1_BNE ;
      (D.f[dPMP])[ktse ] = f1_BNW ;
      (D.f[dMPP])[ktnw ] = f1_BSE ;
      (D.f[dPPM])[kbne ] = f1_TSW ;
      (D.f[dMMM])[kbsw ] = f1_TNE ;
      (D.f[dPMM])[kbse ] = f1_TNW ;
      (D.f[dMPM])[kbnw ] = f1_TSE ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceNEQ27(
    real* rhoBC,
    real* distributions,
    int* bcNodeIndices,
    int* bcNeighborIndices,
    int numberOfBCnodes,
    real omega1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   //! The pressure boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   ////////////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(nodeIndex < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local pressure
      //!
      real rhoBClocal = rhoBC[nodeIndex];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int KQK  = bcNodeIndices[nodeIndex];
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
      //! - Set neighbor indices (necessary for indirect addressing) for neighboring node
      //!
      unsigned int K1QK  = bcNeighborIndices[nodeIndex];
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
      //! - Set local distributions for neighboring node
      //!
      real f1_W    = (dist.f[dP00])[k1e   ];
      real f1_E    = (dist.f[dM00])[k1w   ];
      real f1_S    = (dist.f[d0P0])[k1n   ];
      real f1_N    = (dist.f[d0M0])[k1s   ];
      real f1_B    = (dist.f[d00P])[k1t   ];
      real f1_T    = (dist.f[d00M])[k1b   ];
      real f1_SW   = (dist.f[dPP0])[k1ne  ];
      real f1_NE   = (dist.f[dMM0])[k1sw  ];
      real f1_NW   = (dist.f[dPM0])[k1se  ];
      real f1_SE   = (dist.f[dMP0])[k1nw  ];
      real f1_BW   = (dist.f[dP0P])[k1te  ];
      real f1_TE   = (dist.f[dM0M])[k1bw  ];
      real f1_TW   = (dist.f[dP0M])[k1be  ];
      real f1_BE   = (dist.f[dM0P])[k1tw  ];
      real f1_BS   = (dist.f[d0PP])[k1tn  ];
      real f1_TN   = (dist.f[d0MM])[k1bs  ];
      real f1_TS   = (dist.f[d0PM])[k1bn  ];
      real f1_BN   = (dist.f[d0MP])[k1ts  ];
      real f1_ZERO = (dist.f[d000])[k1zero];
      real f1_BSW  = (dist.f[dPPP])[k1tne ];
      real f1_BNE  = (dist.f[dMMP])[k1tsw ];
      real f1_BNW  = (dist.f[dPMP])[k1tse ];
      real f1_BSE  = (dist.f[dMPP])[k1tnw ];
      real f1_TSW  = (dist.f[dPPM])[k1bne ];
      real f1_TNE  = (dist.f[dMMM])[k1bsw ];
      real f1_TNW  = (dist.f[dPMM])[k1bse ];
      real f1_TSE  = (dist.f[dMPM])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities (for neighboring node)
      //!
      real drho1 = f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                   f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW +
                   f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((dist.f[d000])[kzero]);

      real vx1  = (((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                   (f1_E - f1_W)) / (c1o1 + drho1);

      real vx2  = ((-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                   (f1_N - f1_S)) / (c1o1 + drho1);

      real vx3  = (((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                   (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                   (f1_T - f1_B)) / (c1o1 + drho1);

      real cusq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      //! subtract the equilibrium (eq) to obtain the non-equilibrium (neq) (for neighboring node)
      //!
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

      ////////////////////////////////////////////////////////////////////////////////
      //! redefine drho1 with rhoBClocal
      //!
      drho1 = rhoBClocal;

      ////////////////////////////////////////////////////////////////////////////////
      //! add the equilibrium (eq), which is calculated with rhoBClocal (for neighboring node)
      //!
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

      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes
      //!
      (dist.f[dP00])[ke   ] = f1_W   ;
      (dist.f[dM00])[kw   ] = f1_E   ;
      (dist.f[d0P0])[kn   ] = f1_S   ;
      (dist.f[d0M0])[ks   ] = f1_N   ;
      (dist.f[d00P])[kt   ] = f1_B   ;
      (dist.f[d00M])[kb   ] = f1_T   ;
      (dist.f[dPP0])[kne  ] = f1_SW  ;
      (dist.f[dMM0])[ksw  ] = f1_NE  ;
      (dist.f[dPM0])[kse  ] = f1_NW  ;
      (dist.f[dMP0])[knw  ] = f1_SE  ;
      (dist.f[dP0P])[kte  ] = f1_BW  ;
      (dist.f[dM0M])[kbw  ] = f1_TE  ;
      (dist.f[dP0M])[kbe  ] = f1_TW  ;
      (dist.f[dM0P])[ktw  ] = f1_BE  ;
      (dist.f[d0PP])[ktn  ] = f1_BS  ;
      (dist.f[d0MM])[kbs  ] = f1_TN  ;
      (dist.f[d0PM])[kbn  ] = f1_TS  ;
      (dist.f[d0MP])[kts  ] = f1_BN  ;
      (dist.f[d000])[kzero] = f1_ZERO;
      (dist.f[dPPP])[ktne ] = f1_BSW ;
      (dist.f[dMMP])[ktsw ] = f1_BNE ;
      (dist.f[dPMP])[ktse ] = f1_BNW ;
      (dist.f[dMPP])[ktnw ] = f1_BSE ;
      (dist.f[dPPM])[kbne ] = f1_TSW ;
      (dist.f[dMMM])[kbsw ] = f1_TNE ;
      (dist.f[dPMM])[kbse ] = f1_TNW ;
      (dist.f[dMPM])[kbnw ] = f1_TSE ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
__global__ void LB_BC_Press_East27(
    int nx,
    int ny,
    int tz,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DD,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
   //thread-index
   int ty = blockIdx.x;
   int tx = threadIdx.x;

   int  k, k1, nxny;                   // Zugriff auf arrays im device

   int  x = tx + STARTOFFX;  // Globaler x-Index
   int  y = ty + STARTOFFY;  // Globaler y-Index
   int  z = tz + STARTOFFZ;  // Globaler z-Index

   k = nx*(ny*z + y) + x;
   nxny = nx*ny;
   k1 = k-nxny;

   if( bcMatD[k] == GEO_PRESS && bcMatD[k1] == GEO_FLUID)
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
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = k;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = k;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = k;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = k;
      unsigned int kbsw = neighborZ[ksw];
      //unsigned int kzero= k;
      //unsigned int ke   = k;
      //unsigned int kw   = k + 1;
      //unsigned int kn   = k;
      //unsigned int ks   = k + nx;
      //unsigned int kt   = k;
      //unsigned int kb   = k + nxny;
      //unsigned int ksw  = k + nx + 1;
      //unsigned int kne  = k;
      //unsigned int kse  = k + nx;
      //unsigned int knw  = k + 1;
      //unsigned int kbw  = k + nxny + 1;
      //unsigned int kte  = k;
      //unsigned int kbe  = k + nxny;
      //unsigned int ktw  = k + 1;
      //unsigned int kbs  = k + nxny + nx;
      //unsigned int ktn  = k;
      //unsigned int kbn  = k + nxny;
      //unsigned int kts  = k + nx;
      //unsigned int ktse = k + nx;
      //unsigned int kbnw = k + nxny + 1;
      //unsigned int ktnw = k + 1;
      //unsigned int kbse = k + nxny + nx;
      //unsigned int ktsw = k + nx + 1;
      //unsigned int kbne = k + nxny;
      //unsigned int ktne = k;
      //unsigned int kbsw = k + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int k1zero= k1;
      unsigned int k1e   = k1;
      unsigned int k1w   = neighborX[k1];
      unsigned int k1n   = k1;
      unsigned int k1s   = neighborY[k1];
      unsigned int k1t   = k1;
      unsigned int k1b   = neighborZ[k1];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = k1;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = k1;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = k1;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = k1;
      unsigned int k1bsw = neighborZ[k1sw];
      //unsigned int k1zero= k1;
      //unsigned int k1e   = k1;
      //unsigned int k1w   = k1 + 1;
      //unsigned int k1n   = k1;
      //unsigned int k1s   = k1 + nx;
      //unsigned int k1t   = k1;
      //unsigned int k1b   = k1 + nxny;
      //unsigned int k1sw  = k1 + nx + 1;
      //unsigned int k1ne  = k1;
      //unsigned int k1se  = k1 + nx;
      //unsigned int k1nw  = k1 + 1;
      //unsigned int k1bw  = k1 + nxny + 1;
      //unsigned int k1te  = k1;
      //unsigned int k1be  = k1 + nxny;
      //unsigned int k1tw  = k1 + 1;
      //unsigned int k1bs  = k1 + nxny + nx;
      //unsigned int k1tn  = k1;
      //unsigned int k1bn  = k1 + nxny;
      //unsigned int k1ts  = k1 + nx;
      //unsigned int k1tse = k1 + nx;
      //unsigned int k1bnw = k1 + nxny + 1;
      //unsigned int k1tnw = k1 + 1;
      //unsigned int k1bse = k1 + nxny + nx;
      //unsigned int k1tsw = k1 + nx + 1;
      //unsigned int k1bne = k1 + nxny;
      //unsigned int k1tne = k1;
      //unsigned int k1bsw = k1 + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
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

      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                        f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      (D.f[dP00])[ke   ] = f1_W   -c2o27*drho1;
      (D.f[dM00])[kw   ] = f1_E   -c2o27*drho1;
      (D.f[d0P0])[kn   ] = f1_S   -c2o27*drho1;
      (D.f[d0M0])[ks   ] = f1_N   -c2o27*drho1;
      (D.f[d00P])[kt   ] = f1_B   -c2o27*drho1;
      (D.f[d00M])[kb   ] = f1_T   -c2o27*drho1;
      (D.f[dPP0])[kne  ] = f1_SW  -c1o54*drho1;
      (D.f[dMM0])[ksw  ] = f1_NE  -c1o54*drho1;
      (D.f[dPM0])[kse  ] = f1_NW  -c1o54*drho1;
      (D.f[dMP0])[knw  ] = f1_SE  -c1o54*drho1;
      (D.f[dP0P])[kte  ] = f1_BW  -c1o54*drho1;
      (D.f[dM0M])[kbw  ] = f1_TE  -c1o54*drho1;
      (D.f[dP0M])[kbe  ] = f1_TW  -c1o54*drho1;
      (D.f[dM0P])[ktw  ] = f1_BE  -c1o54*drho1;
      (D.f[d0PP])[ktn  ] = f1_BS  -c1o54*drho1;
      (D.f[d0MM])[kbs  ] = f1_TN  -c1o54*drho1;
      (D.f[d0PM])[kbn  ] = f1_TS  -c1o54*drho1;
      (D.f[d0MP])[kts  ] = f1_BN  -c1o54*drho1;
      (D.f[d000])[kzero] = f1_ZERO-c8o27*drho1;
      (D.f[dPPP])[ktne ] = f1_BSW -c1o216*drho1;
      (D.f[dMMP])[ktsw ] = f1_BNE -c1o216*drho1;
      (D.f[dPMP])[ktse ] = f1_BNW -c1o216*drho1;
      (D.f[dMPP])[ktnw ] = f1_BSE -c1o216*drho1;
      (D.f[dPPM])[kbne ] = f1_TSW -c1o216*drho1;
      (D.f[dMMM])[kbsw ] = f1_TNE -c1o216*drho1;
      (D.f[dPMM])[kbse ] = f1_TNW -c1o216*drho1;
      (D.f[dMPM])[kbnw ] = f1_TSE -c1o216*drho1;
   }
   __syncthreads();
}
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QPressDevice27(
    real* rhoBC,
    real* DD,
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
      real q, vx1, vx2, vx3, drho;
      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                  (f_E - f_W);


      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S);

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B);

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      drho = rhoBC[k];
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM00])[kw]=c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         //(D.f[dP00])[ke]=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP00])[ke]=c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         //(D.f[dM00])[kw]=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0M0])[ks]=c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         //(D.f[d0P0])[kn]=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0P0])[kn]=c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         //(D.f[d0M0])[ks]=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00M])[kb]=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         //(D.f[d00P])[kt]=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00P])[kt]=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         //(D.f[d00M])[kb]=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMM0])[ksw]=c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         //(D.f[dPP0])[kne]=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPP0])[kne]=c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         //(D.f[dMM0])[ksw]=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMP0])[knw]=c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         //(D.f[dPM0])[kse]=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPM0])[kse]=c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         //(D.f[dMP0])[knw]=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0M])[kbw]=c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         //(D.f[dP0P])[kte]=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0P])[kte]=c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         //(D.f[dM0M])[kbw]=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0P])[ktw]=c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         //(D.f[dP0M])[kbe]=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0M])[kbe]=c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         //(D.f[dM0P])[ktw]=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MM])[kbs]=c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         //(D.f[d0PP])[ktn]=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PP])[ktn]=c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         //(D.f[d0MM])[kbs]=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MP])[kts]=c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         //(D.f[d0PM])[kbn]=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PM])[kbn]=c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         //(D.f[d0MP])[kts]=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMM])[kbsw]=c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         //(D.f[dPPP])[ktne]=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPP])[ktne]=c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         //(D.f[dMMM])[kbsw]=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMP])[ktsw]=c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         //(D.f[dPPM])[kbne]=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPM])[kbne]=c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         //(D.f[dMMP])[ktsw]=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPM])[kbnw]=c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         //(D.f[dPMP])[ktse]=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMP])[ktse]=c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         //(D.f[dMPM])[kbnw]=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPP])[ktnw]=c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         //(D.f[dPMM])[kbse]=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMM])[kbse]=c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         //(D.f[dMPP])[ktnw]=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceAntiBB27(
    real* rhoBC,
    real* vx,
    real* vy,
    real* vz,
    real* DD,
    int* k_Q,
    real* QQ,
    int numberOfBCnodes,
    real om1,
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
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

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
      f_ZERO = (D.f[d000])[kzero];
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1, vx2, vx3, drho;
      //vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //            (f_E - f_W);


      //vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //            (f_N - f_S);

      //vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //            (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //            (f_T - f_B);

      //real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      real drho    = f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+
                  f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      drho = drho - rhoBC[k];
     drho *= 0.01f;
      ////////////////////////////////////////////////////////////////////////////////
     real q;
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM00])[kw]=f_W-c2o27*drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP00])[ke]=f_E-c2o27*drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0M0])[ks]=f_S-c2o27*drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0P0])[kn]=f_N-c2o27*drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00M])[kb]=f_B-c2o27*drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d00P])[kt]=f_T-c2o27*drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMM0])[ksw]=f_SW-c1o54*drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPP0])[kne]=f_NE-c1o54*drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMP0])[knw]=f_NW-c1o54*drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPM0])[kse]=f_SE-c1o54*drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0M])[kbw]=f_BW-c1o54*drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0P])[kte]=f_TE-c1o54*drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dM0P])[ktw]=f_TW-c1o54*drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dP0M])[kbe]=f_BE-c1o54*drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MM])[kbs]=f_BS-c1o54*drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PP])[ktn]=f_TN-c1o54*drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0MP])[kts]=f_TS-c1o54*drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[d0PM])[kbn]=f_BN-c1o54*drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMM])[kbsw]=f_BSW-c1o216*drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPP])[ktne]=f_TNE-c1o216*drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMMP])[ktsw]=f_TSW-c1o216*drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPPM])[kbne]=f_BNE-c1o216*drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPM])[kbnw]=f_BNW-c1o216*drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMP])[ktse]=f_TSE-c1o216*drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dMPP])[ktnw]=f_TNW-c1o216*drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[dPMM])[kbse]=f_BSE-c1o216*drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceFixBackflow27(
    real* rhoBC,
    real* DD,
    int* k_Q,
    int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
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
      real deltaRho;
      ////////////////////////////////////////////////////////////////////////////////
      deltaRho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         (D.f[dM00])[kw]       = c2o27  * deltaRho;
         (D.f[dP00])[ke]       = c2o27  * deltaRho;
         (D.f[d0M0])[ks]       = c2o27  * deltaRho;
         (D.f[d0P0])[kn]       = c2o27  * deltaRho;
         (D.f[d00M])[kb]       = c2o27  * deltaRho;
         (D.f[d00P])[kt]       = c2o27  * deltaRho;
         (D.f[dMM0])[ksw]     = c1o54  * deltaRho;
         (D.f[dPP0])[kne]     = c1o54  * deltaRho;
         (D.f[dMP0])[knw]     = c1o54  * deltaRho;
         (D.f[dPM0])[kse]     = c1o54  * deltaRho;
         (D.f[dM0M])[kbw]     = c1o54  * deltaRho;
         (D.f[dP0P])[kte]     = c1o54  * deltaRho;
         (D.f[dM0P])[ktw]     = c1o54  * deltaRho;
         (D.f[dP0M])[kbe]     = c1o54  * deltaRho;
         (D.f[d0MM])[kbs]     = c1o54  * deltaRho;
         (D.f[d0PP])[ktn]     = c1o54  * deltaRho;
         (D.f[d0MP])[kts]     = c1o54  * deltaRho;
         (D.f[d0PM])[kbn]     = c1o54  * deltaRho;
         (D.f[dMMM])[kbsw]   = c1o216 * deltaRho;
         (D.f[dPPP])[ktne]   = c1o216 * deltaRho;
         (D.f[dMMP])[ktsw]   = c1o216 * deltaRho;
         (D.f[dPPM])[kbne]   = c1o216 * deltaRho;
         (D.f[dMPM])[kbnw]   = c1o216 * deltaRho;
         (D.f[dPMP])[ktse]   = c1o216 * deltaRho;
         (D.f[dMPP])[ktnw]   = c1o216 * deltaRho;
         (D.f[dPMM])[kbse]   = c1o216 * deltaRho;
         (D.f[d000])[kzero] = c8o27  * deltaRho;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceDirDepBot27(
    real* rhoBC,
    real* DD,
    int* k_Q,
    int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
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
      real rho;
      ////////////////////////////////////////////////////////////////////////////////
      rho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E,f_W,f_N,f_S,f_T,f_NE,f_SW,f_SE,f_NW,f_TE,f_TW,f_TN,f_TS,f_ZERO,f_TNE,f_TSW,f_TSE,f_TNW;//,
            //f_B,f_BW,f_BE,f_BS,f_BN,f_BSW,f_BNE,f_BNW,f_BSE;

      f_E    = (D.f[dP00])[ke   ];
      f_W    = (D.f[dM00])[kw   ];
      f_N    = (D.f[d0P0])[kn   ];
      f_S    = (D.f[d0M0])[ks   ];
      f_T    = (D.f[d00P])[kt   ];
      f_NE   = (D.f[dPP0])[kne  ];
      f_SW   = (D.f[dMM0])[ksw  ];
      f_SE   = (D.f[dPM0])[kse  ];
      f_NW   = (D.f[dMP0])[knw  ];
      f_TE   = (D.f[dP0P])[kte  ];
      f_TW   = (D.f[dM0P])[ktw  ];
      f_TN   = (D.f[d0PP])[ktn  ];
      f_TS   = (D.f[d0MP])[kts  ];
      f_ZERO = (D.f[d000])[kzero];
      f_TNE  = (D.f[dPPP])[ktne ];
      f_TSW  = (D.f[dMMP])[ktsw ];
      f_TSE  = (D.f[dPMP])[ktse ];
      f_TNW  = (D.f[dMPP])[ktnw ];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //f_B   = (four*rho- four*f_SW-     eight*f_TSW-four*f_W-   eight*f_TW- four*f_NW-     eight*f_TNW-four*f_S-   eight*f_TS-four*f_ZERO+     f_T-four*f_N-   eight*f_TN- four*f_SE-     eight*f_TSE-four*f_E-   eight*f_TE- four*f_NE-     eight*f_TNE)/nine;
      //f_BW  = ( two*rho+      f_SW-      four*f_TSW+     f_W-    four*f_TW+      f_NW-      four*f_TNW- two*f_S-    four*f_TS- two*f_ZERO-four*f_T- two*f_N-    four*f_TN- five*f_SE-      four*f_TSE-five*f_E+fourteen*f_TE- five*f_NE-      four*f_TNE)/eighteen;
      //f_BE  = ( two*rho- five*f_SW-      four*f_TSW-five*f_W+fourteen*f_TW- five*f_NW-      four*f_TNW- two*f_S-    four*f_TS- two*f_ZERO-four*f_T- two*f_N-    four*f_TN+      f_SE-      four*f_TSE+     f_E-    four*f_TE+      f_NE-      four*f_TNE)/eighteen;
      //f_BS  = ( two*rho+      f_SW-      four*f_TSW- two*f_W-    four*f_TW- five*f_NW-      four*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N+fourteen*f_TN+      f_SE-      four*f_TSE- two*f_E-    four*f_TE- five*f_NE-      four*f_TNE)/eighteen;
      //f_BN  = ( two*rho- five*f_SW-      four*f_TSW- two*f_W-    four*f_TW+      f_NW-      four*f_TNW-five*f_S+fourteen*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN- five*f_SE-      four*f_TSE- two*f_E-    four*f_TE+      f_NE-      four*f_TNE)/eighteen;
      //f_BSW = ( two*rho+ four*f_SW-      four*f_TSW+     f_W-    four*f_TW-  two*f_NW-      four*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N-    four*f_TN-  two*f_SE-      four*f_TSE-five*f_E-    four*f_TE-eight*f_NE+sixtyeight*f_TNE)/seventytwo;
      //f_BNE = ( two*rho-eight*f_SW+sixtyeight*f_TSW-five*f_W-    four*f_TW-  two*f_NW-      four*f_TNW-five*f_S-    four*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN-  two*f_SE-      four*f_TSE+     f_E-    four*f_TE+ four*f_NE-      four*f_TNE)/seventytwo;
      //f_BNW = ( two*rho-  two*f_SW-      four*f_TSW+     f_W-    four*f_TW+ four*f_NW-      four*f_TNW-five*f_S-    four*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN-eight*f_SE+sixtyeight*f_TSE-five*f_E-    four*f_TE-  two*f_NE-      four*f_TNE)/seventytwo;
      //f_BSE = ( two*rho-  two*f_SW-      four*f_TSW-five*f_W-    four*f_TW-eight*f_NW+sixtyeight*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N-    four*f_TN+ four*f_SE-      four*f_TSE+     f_E-    four*f_TE-  two*f_NE-      four*f_TNE)/seventytwo;

      //real drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      //real vx1     =  (f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW);
      //real vx2     =  (f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW);
      //real vx3     =  (f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW);

      //real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //(D.f[d000])[kzero] = c8over27*  (drho-cusq);
      //(D.f[dP00])[ke]    = c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
      //(D.f[dM00])[kw]    = c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
      //(D.f[d0P0])[kn]     = c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
      //(D.f[d0M0])[ks]    = c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
      //(D.f[d00P])[kt]    = c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
      //(D.f[d00M])[kb]    = c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
      //(D.f[dPP0])[kne]   = c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
      //(D.f[dMM0])[ksw]   = c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
      //(D.f[dPM0])[kse]   =  c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
      //(D.f[dMP0])[knw]   =  c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
      //(D.f[dP0P])[kte]   =  c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
      //(D.f[dM0M])[kbw]   =  c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
      //(D.f[dP0M])[kbe]   =  c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
      //(D.f[dM0P])[ktw]   =  c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
      //(D.f[d0PP])[ktn]   =  c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
      //(D.f[d0MM])[kbs]   =  c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
      //(D.f[d0PM])[kbn]   =  c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
      //(D.f[d0MP])[kts]   =  c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
      //(D.f[dPPP])[ktne]  =  c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
      //(D.f[dMMM])[kbsw]  =  c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
      //(D.f[dPPM])[kbne]  =  c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
      //(D.f[dMMP])[ktsw]  =  c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
      //(D.f[dPMP])[ktse]  =  c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
      //(D.f[dMPM])[kbnw]  =  c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
      //(D.f[dPMM])[kbse]  =  c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
      //(D.f[dMPP])[ktnw]  =  c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
      real drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_NE+f_SW+f_SE+f_NW+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      real dTop   =    f_T+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      (D.f[d00M])[kb]     = (f_T+c2o27)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c2o27;
      (D.f[dM0M])[kbw]   = (f_TW+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[dP0M])[kbe]   = (f_TE+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[d0MM])[kbs]   = (f_TS+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[d0PM])[kbn]   = (f_TN+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[dMMM])[kbsw] = (f_TSW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[dPPM])[kbne] = (f_TNE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[dMPM])[kbnw] = (f_TNW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[dPMM])[kbse] = (f_TSE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




































__host__ __device__ real computeOutflowDistribution(const real* const &f, const real* const &f1, const int dir, const real cs)
{
   return f1[dir] * cs + (c1o1 - cs) * f[dir];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressNoRhoDevice27(
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
      case MZZ:
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

      case PZZ:
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

      case ZMZ:
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

      case ZPZ:
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

      case ZZM:
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

      case ZZP:
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



































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ real computeOutflowDistribution(const real* const &f, const real* const &f1, const int dir, const real rhoCorrection, const real cs, const real weight)
{
   return f1[dir  ] * cs + (c1o1 - cs) * f[dir  ] - weight *rhoCorrection;
}

__global__ void QPressZeroRhoOutflowDevice27(
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
      case MZZ:
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

      case PZZ:
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

      case ZMZ:
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

      case ZPZ:
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

      case ZZM:
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

      case ZZP:
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceOld27(
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
    bool isEvenTimestep)
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
      if (isEvenTimestep==false)
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

     //drho1 = (drho1 + rhoBC[k])/2.f;
     drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[dP00])[ke   ] = f1_W   -c2o27*drho1;   //  c1o100;  // zero;  //
      (D.f[dM00])[kw   ] = f1_E   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[d0P0])[kn   ] = f1_S   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[d0M0])[ks   ] = f1_N   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[d00P])[kt   ] = f1_B   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[d00M])[kb   ] = f1_T   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[dPP0])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dMM0])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dPM0])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dMP0])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dP0P])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dM0M])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dP0M])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dM0P])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0PP])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0MM])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0PM])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0MP])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d000])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[dPPP])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMMP])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPMP])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMPP])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPPM])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMMM])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPMM])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMPM])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceEQZ27(
    real* rhoBC,
    real* DD,
    int* k_Q,
    int* k_N,
    real* kTestRE,
    int numberOfBCnodes,
    real om1,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
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
    //   Distributions27 kDistTest;
    //      kDistTest.f[dP00] = &kTestRE[dP00 * numberOfBCnodes];
    //      kDistTest.f[dM00] = &kTestRE[dM00 * numberOfBCnodes];
    //      kDistTest.f[d0P0] = &kTestRE[d0P0 * numberOfBCnodes];
    //      kDistTest.f[d0M0] = &kTestRE[d0M0 * numberOfBCnodes];
    //      kDistTest.f[d00P] = &kTestRE[d00P * numberOfBCnodes];
    //      kDistTest.f[d00M] = &kTestRE[d00M * numberOfBCnodes];
    //      kDistTest.f[dPP0] = &kTestRE[dPP0 * numberOfBCnodes];
    //      kDistTest.f[dMM0] = &kTestRE[dMM0 * numberOfBCnodes];
    //      kDistTest.f[dPM0] = &kTestRE[dPM0 * numberOfBCnodes];
    //      kDistTest.f[dMP0] = &kTestRE[dMP0 * numberOfBCnodes];
    //      kDistTest.f[dP0P] = &kTestRE[dP0P * numberOfBCnodes];
    //      kDistTest.f[dM0M] = &kTestRE[dM0M * numberOfBCnodes];
    //      kDistTest.f[dP0M] = &kTestRE[dP0M * numberOfBCnodes];
    //      kDistTest.f[dM0P] = &kTestRE[dM0P * numberOfBCnodes];
    //      kDistTest.f[d0PP] = &kTestRE[d0PP * numberOfBCnodes];
    //      kDistTest.f[d0MM] = &kTestRE[d0MM * numberOfBCnodes];
    //      kDistTest.f[d0PM] = &kTestRE[d0PM * numberOfBCnodes];
    //      kDistTest.f[d0MP] = &kTestRE[d0MP * numberOfBCnodes];
    //      kDistTest.f[d000] = &kTestRE[d000 * numberOfBCnodes];
    //      kDistTest.f[dPPP] = &kTestRE[dPPP * numberOfBCnodes];
    //      kDistTest.f[dMMP] = &kTestRE[dMMP * numberOfBCnodes];
    //      kDistTest.f[dPMP] = &kTestRE[dPMP * numberOfBCnodes];
    //      kDistTest.f[dMPP] = &kTestRE[dMPP * numberOfBCnodes];
    //      kDistTest.f[dPPM] = &kTestRE[dPPM * numberOfBCnodes];
    //      kDistTest.f[dMMM] = &kTestRE[dMMM * numberOfBCnodes];
    //      kDistTest.f[dPMM] = &kTestRE[dPMM * numberOfBCnodes];
    //      kDistTest.f[dMPM] = &kTestRE[dMPM * numberOfBCnodes];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   //f1_W    = (D.f[dP00])[k1e   ];
   //   //f1_E    = (D.f[dM00])[k1w   ];
   //   //f1_S    = (D.f[d0P0])[k1n   ];
   //   //f1_N    = (D.f[d0M0])[k1s   ];
   //   //f1_B    = (D.f[d00P])[k1t   ];
   //   //f1_T    = (D.f[d00M])[k1b   ];
   //   //f1_SW   = (D.f[dPP0])[k1ne  ];
   //   //f1_NE   = (D.f[dMM0])[k1sw  ];
   //   //f1_NW   = (D.f[dPM0])[k1se  ];
   //   //f1_SE   = (D.f[dMP0])[k1nw  ];
   //   //f1_BW   = (D.f[dP0P])[k1te  ];
   //   //f1_TE   = (D.f[dM0M])[k1bw  ];
   //   //f1_TW   = (D.f[dP0M])[k1be  ];
   //   //f1_BE   = (D.f[dM0P])[k1tw  ];
   //   //f1_BS   = (D.f[d0PP])[k1tn  ];
   //   //f1_TN   = (D.f[d0MM])[k1bs  ];
   //   //f1_TS   = (D.f[d0PM])[k1bn  ];
   //   //f1_BN   = (D.f[d0MP])[k1ts  ];
   //   //f1_ZERO = (D.f[d000])[k1zero];
   //   //f1_BSW  = (D.f[dPPP])[k1tne ];
   //   //f1_BNE  = (D.f[dMMP])[k1tsw ];
   //   //f1_BNW  = (D.f[dPMP])[k1tse ];
   //   //f1_BSE  = (D.f[dMPP])[k1tnw ];
   //   //f1_TSW  = (D.f[dPPM])[k1bne ];
   //   //f1_TNE  = (D.f[dMMM])[k1bsw ];
   //   //f1_TNW  = (D.f[dPMM])[k1bse ];
   //   //f1_TSE  = (D.f[dMPM])[k1bnw ];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   f1_E    = (D.f[dP00])[k1e   ];
   //   f1_W    = (D.f[dM00])[k1w   ];
   //   f1_N    = (D.f[d0P0])[k1n   ];
   //   f1_S    = (D.f[d0M0])[k1s   ];
   //   f1_T    = (D.f[d00P])[k1t   ];
   //   f1_B    = (D.f[d00M])[k1b   ];
   //   f1_NE   = (D.f[dPP0])[k1ne  ];
   //   f1_SW   = (D.f[dMM0])[k1sw  ];
   //   f1_SE   = (D.f[dPM0])[k1se  ];
   //   f1_NW   = (D.f[dMP0])[k1nw  ];
   //   f1_TE   = (D.f[dP0P])[k1te  ];
   //   f1_BW   = (D.f[dM0M])[k1bw  ];
   //   f1_BE   = (D.f[dP0M])[k1be  ];
   //   f1_TW   = (D.f[dM0P])[k1tw  ];
   //   f1_TN   = (D.f[d0PP])[k1tn  ];
   //   f1_BS   = (D.f[d0MM])[k1bs  ];
   //   f1_BN   = (D.f[d0PM])[k1bn  ];
   //   f1_TS   = (D.f[d0MP])[k1ts  ];
   //   f1_ZERO = (D.f[d000])[k1zero];
   //   f1_TNE  = (D.f[dPPP])[k1tne ];
   //   f1_TSW  = (D.f[dMMP])[k1tsw ];
   //   f1_TSE  = (D.f[dPMP])[k1tse ];
   //   f1_TNW  = (D.f[dMPP])[k1tnw ];
   //   f1_BNE  = (D.f[dPPM])[k1bne ];
   //   f1_BSW  = (D.f[dMMM])[k1bsw ];
   //   f1_BSE  = (D.f[dPMM])[k1bse ];
   //   f1_BNW  = (D.f[dMPM])[k1bnw ];
   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////
   //   real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+ f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;
     //real vx1      = (((f1_TNE-f1_BSW)+(f1_BSE-f1_TNW)+(f1_BNE-f1_TSW)+(f1_TSE-f1_BNW)) + (((f1_NE-f1_SW)+(f1_TE-f1_BW))+((f1_SE-f1_NW)+(f1_BE-f1_TW))) + (f1_E-f1_W)) / (one + drho1);
     //real vx2      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_BNE-f1_TSW)+(f1_BNW-f1_TSE)) + (((f1_NE-f1_SW)+(f1_TN-f1_BS))+((f1_BN-f1_TS)+(f1_NW-f1_SE))) + (f1_N-f1_S)) / (one + drho1);
     //real vx3      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_TSW-f1_BNE)+(f1_TSE-f1_BNW)) + (((f1_TE-f1_BW)+(f1_TN-f1_BS))+((f1_TW-f1_BE)+(f1_TS-f1_BN))) + (f1_T-f1_B)) / (one + drho1);
   //   //////////////////////////////////////////////////////////////////////////
     ////real omega = om1;
   //   real cusq  = c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
   //   //////////////////////////////////////////////////////////////////////////
     ////T�st MK
     ////if(vx1 < zero) vx1 = zero;
   //   //////////////////////////////////////////////////////////////////////////
     ////becomes higher with neighbor source and lower with local source
   //   //real fZERO = c8over27*  (rhoBC[k]-(one + rhoBC[k])*(cusq))                                                           ;
   //   //real fE    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //   //real fW    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //   //real fN    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //   //real fS    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //   //real fT    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //   //real fB    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //   //real fNE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //   //real fSW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //   //real fSE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //   //real fNW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //   //real fTE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //   //real fBW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //   //real fBE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //   //real fTW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //   //real fTN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //   //real fBS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //   //real fBN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //   //real fTS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //   //real fTNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //   //real fBSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //   //real fBNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //   //real fTSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //   //real fTSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //   //real fBNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //   //real fBSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //   //real fTNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
   //   //////////////////////////////////////////////////////////////////////////
     //// based on VirtualFluids (kucher + fard)
   //   real fZERO = c8over27  * rhoBC[k] * (one                                                                      - cusq);
   //   real fE    = c2over27  * rhoBC[k] * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   //   real fW    = c2over27  * rhoBC[k] * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   //   real fN    = c2over27  * rhoBC[k] * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   //   real fS    = c2over27  * rhoBC[k] * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   //   real fT    = c2over27  * rhoBC[k] * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   //   real fB    = c2over27  * rhoBC[k] * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   //   real fNE   = c1over54  * rhoBC[k] * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   //   real fSW   = c1over54  * rhoBC[k] * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   //   real fSE   = c1over54  * rhoBC[k] * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   //   real fNW   = c1over54  * rhoBC[k] * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   //   real fTE   = c1over54  * rhoBC[k] * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   //   real fBW   = c1over54  * rhoBC[k] * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   //   real fBE   = c1over54  * rhoBC[k] * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   //   real fTW   = c1over54  * rhoBC[k] * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   //   real fTN   = c1over54  * rhoBC[k] * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   //   real fBS   = c1over54  * rhoBC[k] * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   //   real fBN   = c1over54  * rhoBC[k] * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   //   real fTS   = c1over54  * rhoBC[k] * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   //   real fTNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   //   real fBSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   //   real fBNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   //   real fTSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   //   real fTSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   //   real fBNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   //   real fBSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   //   real fTNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);
   ////   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //////test
   ////   real fZERO = c8over27  * ((drho1 + rhoBC[k]) / two) * (one                                                                      - cusq);
   ////   real fE    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   ////   real fW    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   ////   real fN    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   ////   real fS    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   ////   real fT    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   ////   real fB    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   ////   real fNE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   ////   real fSW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   ////   real fSE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   ////   real fNW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   ////   real fTE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   ////   real fBW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   ////   real fBE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   ////   real fTW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   ////   real fTN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   ////   real fBS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   ////   real fBN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   ////   real fTS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   ////   real fTNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   ////   real fBSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   ////   real fBNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   ////   real fTSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   ////   real fTSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   ////   real fBNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   ////   real fBSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   ////   real fTNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);

         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //double mfabb = (D.f[dP00])[k1e   ];
         //double mfcbb = (D.f[dM00])[k1w   ];
         //double mfbab = (D.f[d0P0])[k1n   ];
         //double mfbcb = (D.f[d0M0])[k1s   ];
         //double mfbba = (D.f[d00P])[k1t   ];
         //double mfbbc = (D.f[d00M])[k1b   ];
         //double mfaab = (D.f[dPP0])[k1ne  ];
         //double mfccb = (D.f[dMM0])[k1sw  ];
         //double mfacb = (D.f[dPM0])[k1se  ];
         //double mfcab = (D.f[dMP0])[k1nw  ];
         //double mfaba = (D.f[dP0P])[k1te  ];
         //double mfcbc = (D.f[dM0M])[k1bw  ];
         //double mfabc = (D.f[dP0M])[k1be  ];
         //double mfcba = (D.f[dM0P])[k1tw  ];
         //double mfbaa = (D.f[d0PP])[k1tn  ];
         //double mfbcc = (D.f[d0MM])[k1bs  ];
         //double mfbac = (D.f[d0PM])[k1bn  ];
         //double mfbca = (D.f[d0MP])[k1ts  ];
         //double mfbbb = (D.f[d000])[k1zero];
         //double mfaaa = (D.f[dPPP])[k1tne ];
         //double mfcca = (D.f[dMMP])[k1tsw ];
         //double mfaca = (D.f[dPMP])[k1tse ];
         //double mfcaa = (D.f[dMPP])[k1tnw ];
         //double mfaac = (D.f[dPPM])[k1bne ];
         //double mfccc = (D.f[dMMM])[k1bsw ];
         //double mfacc = (D.f[dPMM])[k1bse ];
         //double mfcac = (D.f[dMPM])[k1bnw ];
         real mfabb = (D.f[dP00])[k1e   ];
         real mfcbb = (D.f[dM00])[k1w   ];
         real mfbab = (D.f[d0P0])[k1n   ];
         real mfbcb = (D.f[d0M0])[k1s   ];
         real mfbba = (D.f[d00P])[k1t   ];
         real mfbbc = (D.f[d00M])[k1b   ];
         real mfaab = (D.f[dPP0])[k1ne  ];
         real mfccb = (D.f[dMM0])[k1sw  ];
         real mfacb = (D.f[dPM0])[k1se  ];
         real mfcab = (D.f[dMP0])[k1nw  ];
         real mfaba = (D.f[dP0P])[k1te  ];
         real mfcbc = (D.f[dM0M])[k1bw  ];
         real mfabc = (D.f[dP0M])[k1be  ];
         real mfcba = (D.f[dM0P])[k1tw  ];
         real mfbaa = (D.f[d0PP])[k1tn  ];
         real mfbcc = (D.f[d0MM])[k1bs  ];
         real mfbac = (D.f[d0PM])[k1bn  ];
         real mfbca = (D.f[d0MP])[k1ts  ];
         real mfbbb = (D.f[d000])[k1zero];
         real mfaaa = (D.f[dPPP])[k1tne ];
         real mfcca = (D.f[dMMP])[k1tsw ];
         real mfaca = (D.f[dPMP])[k1tse ];
         real mfcaa = (D.f[dMPP])[k1tnw ];
         real mfaac = (D.f[dPPM])[k1bne ];
         real mfccc = (D.f[dMMM])[k1bsw ];
         real mfacc = (D.f[dPMM])[k1bse ];
         real mfcac = (D.f[dMPM])[k1bnw ];

         //real mfcbb = (D.f[dP00])[ke   ];
         //real mfabb = (D.f[dM00])[kw   ];
         //real mfbcb = (D.f[d0P0])[kn   ];
         //real mfbab = (D.f[d0M0])[ks   ];
         //real mfbbc = (D.f[d00P])[kt   ];
         //real mfbba = (D.f[d00M])[kb   ];
         //real mfccb = (D.f[dPP0])[kne  ];
         //real mfaab = (D.f[dMM0])[ksw  ];
         //real mfcab = (D.f[dPM0])[kse  ];
         //real mfacb = (D.f[dMP0])[knw  ];
         //real mfcbc = (D.f[dP0P])[kte  ];
         //real mfaba = (D.f[dM0M])[kbw  ];
         //real mfcba = (D.f[dP0M])[kbe  ];
         //real mfabc = (D.f[dM0P])[ktw  ];
         //real mfbcc = (D.f[d0PP])[ktn  ];
         //real mfbaa = (D.f[d0MM])[kbs  ];
         //real mfbca = (D.f[d0PM])[kbn  ];
         //real mfbac = (D.f[d0MP])[kts  ];
         //real mfbbb = (D.f[d000])[kzero];
         //real mfccc = (D.f[dPPP])[ktne ];
         //real mfaac = (D.f[dMMP])[ktsw ];
         //real mfcac = (D.f[dPMP])[ktse ];
         //real mfacc = (D.f[dMPP])[ktnw ];
         //real mfcca = (D.f[dPPM])[kbne ];
         //real mfaaa = (D.f[dMMM])[kbsw ];
         //real mfcaa = (D.f[dPMM])[kbse ];
         //real mfaca = (D.f[dMPM])[kbnw ];
         ////////////////////////////////////////////////////////////////////////////////////
         //real rho   = (((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) +
         //				(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
         //				((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb) + one;//!!!!Achtung + one
         ////////////////////////////////////////////////////////////////////////////////////
         real rho = rhoBC[k];
         ////////////////////////////////////////////////////////////////////////////////////
         real OoRho = c1o1 / (rho * 1.5f);
         ////////////////////////////////////////////////////////////////////////////////////
         real vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) +
                       (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                         (mfcbb-mfabb)) * OoRho;
         real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) +
                         (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                           (mfbcb-mfbab)) * OoRho;
         real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) +
                         (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                           (mfbbc-mfbba)) * OoRho;
         /////////////////////////
         //Test Values
         //double vvx    = 0.016;
         //double vvy    = zero;
         //double vvz    = zero;
         ////////////////////////////////////////////////////////////////////////////////////////
         ////round off error test
         //if(vvx!=zero){
         //	(kDistTest.f[dP00])[k] = mfabb;
         //	(kDistTest.f[dM00])[k] = mfcbb;
         //	(kDistTest.f[d0P0])[k] = mfbab;
         //	(kDistTest.f[d0M0])[k] = mfbcb;
         //	(kDistTest.f[d00P])[k] = mfbba;
         //	(kDistTest.f[d00M])[k] = mfbbc;
         //	(kDistTest.f[dPP0])[k] = mfaab;
         //	(kDistTest.f[dMM0])[k] = mfccb;
         //	(kDistTest.f[dPM0])[k] = mfacb;
         //	(kDistTest.f[dMP0])[k] = mfcab;
         //	(kDistTest.f[dP0P])[k] = mfaba;
         //	(kDistTest.f[dM0M])[k] = mfcbc;
         //	(kDistTest.f[dP0M])[k] = mfabc;
         //	(kDistTest.f[dM0P])[k] = mfcba;
         //	(kDistTest.f[d0PP])[k] = mfbaa;
         //	(kDistTest.f[d0MM])[k] = mfbcc;
         //	(kDistTest.f[d0PM])[k] = mfbac;
         //	(kDistTest.f[d0MP])[k] = mfbca;
         //	(kDistTest.f[d000])[k] = KQK;
         //	(kDistTest.f[dPPP])[k] = mfaaa;
         //	(kDistTest.f[dMMP])[k] = mfcca;
         //	(kDistTest.f[dPMP])[k] = mfaca;
         //	(kDistTest.f[dMPP])[k] = mfcaa;
         //	(kDistTest.f[dPPM])[k] = mfaac;
         //	(kDistTest.f[dMMM])[k] = mfccc;
         //	(kDistTest.f[dPMM])[k] = mfacc;
         //	(kDistTest.f[dMPM])[k] = mfcac;
         //}else{
         //	(kDistTest.f[dP00])[k] = zero;
         //	(kDistTest.f[dM00])[k] = zero;
         //	(kDistTest.f[d0P0])[k] = zero;
         //	(kDistTest.f[d0M0])[k] = zero;
         //	(kDistTest.f[d00P])[k] = zero;
         //	(kDistTest.f[d00M])[k] = zero;
         //	(kDistTest.f[dPP0])[k] = zero;
         //	(kDistTest.f[dMM0])[k] = zero;
         //	(kDistTest.f[dPM0])[k] = zero;
         //	(kDistTest.f[dMP0])[k] = zero;
         //	(kDistTest.f[dP0P])[k] = zero;
         //	(kDistTest.f[dM0M])[k] = zero;
         //	(kDistTest.f[dP0M])[k] = zero;
         //	(kDistTest.f[dM0P])[k] = zero;
         //	(kDistTest.f[d0PP])[k] = zero;
         //	(kDistTest.f[d0MM])[k] = zero;
         //	(kDistTest.f[d0PM])[k] = zero;
         //	(kDistTest.f[d0MP])[k] = zero;
         //	(kDistTest.f[d000])[k] = zero;
         //	(kDistTest.f[dPPP])[k] = zero;
         //	(kDistTest.f[dMMP])[k] = zero;
         //	(kDistTest.f[dPMP])[k] = zero;
         //	(kDistTest.f[dMPP])[k] = zero;
         //	(kDistTest.f[dPPM])[k] = zero;
         //	(kDistTest.f[dMMM])[k] = zero;
         //	(kDistTest.f[dPMM])[k] = zero;
         //	(kDistTest.f[dMPM])[k] = zero;
         //}

         //////////////////////////////////////////////////////////////////////////////////////
         //// first bad fix for negative x velocity
         ////if(vvx > zero) vvx = zero;
         //////////////////////////////////////////////////////////////////////////////////////
         ////// second bad fix for negative x velocity
         ////if(vvx > zero){
         ////	vvx = -vvx;
         ////	vvy = -vvy;
         ////	vvz = -vvz;
         ////}
         ////////////////////////////////////////////////////////////////////////////////////
         double vx2    = vvx * vvx;
         double vy2    = vvy * vvy;
         double vz2    = vvz * vvz;
         //////////////////////////////////////////////////////////////////////////////////
         //original
            real XXb    = -c2o3 + vx2;
            real XXc    = -c1o2 * (XXb + c1o1 + vvx);
            real XXa    = XXc + vvx;
            real YYb    = -c2o3 + vy2;
            real YYc    = -c1o2 * (YYb + c1o1 + vvy);
            real YYa    = YYc + vvy;
            real ZZb    = -c2o3 + vz2;
            real ZZc    = -c1o2 * (ZZb + c1o1 + vvz);
            real ZZa    = ZZc + vvz;
         //////////////////////////////////////////////////////////////////////////////////
         //unkonditioniert
            mfcbb = -(rhoBC[k] + c1o1) * XXc * YYb * ZZb - c2o27;
         mfabb = -(rhoBC[k] + c1o1) * XXa * YYb * ZZb - c2o27;
         mfbcb = -(rhoBC[k] + c1o1) * XXb * YYc * ZZb - c2o27;
         mfbab = -(rhoBC[k] + c1o1) * XXb * YYa * ZZb - c2o27;
         mfbbc = -(rhoBC[k] + c1o1) * XXb * YYb * ZZc - c2o27;
         mfbba = -(rhoBC[k] + c1o1) * XXb * YYb * ZZa - c2o27;
         mfccb = -(rhoBC[k] + c1o1) * XXc * YYc * ZZb - c1o54;
         mfaab = -(rhoBC[k] + c1o1) * XXa * YYa * ZZb - c1o54;
         mfcab = -(rhoBC[k] + c1o1) * XXc * YYa * ZZb - c1o54;
         mfacb = -(rhoBC[k] + c1o1) * XXa * YYc * ZZb - c1o54;
         mfcbc = -(rhoBC[k] + c1o1) * XXc * YYb * ZZc - c1o54;
         mfaba = -(rhoBC[k] + c1o1) * XXa * YYb * ZZa - c1o54;
         mfcba = -(rhoBC[k] + c1o1) * XXc * YYb * ZZa - c1o54;
         mfabc = -(rhoBC[k] + c1o1) * XXa * YYb * ZZc - c1o54;
         mfbcc = -(rhoBC[k] + c1o1) * XXb * YYc * ZZc - c1o54;
         mfbaa = -(rhoBC[k] + c1o1) * XXb * YYa * ZZa - c1o54;
         mfbca = -(rhoBC[k] + c1o1) * XXb * YYc * ZZa - c1o54;
         mfbac = -(rhoBC[k] + c1o1) * XXb * YYa * ZZc - c1o54;
         mfbbb = -(rhoBC[k] + c1o1) * XXb * YYb * ZZb - c8o27;
         mfccc = -(rhoBC[k] + c1o1) * XXc * YYc * ZZc - c1o216;
         mfaac = -(rhoBC[k] + c1o1) * XXa * YYa * ZZc - c1o216;
         mfcac = -(rhoBC[k] + c1o1) * XXc * YYa * ZZc - c1o216;
         mfacc = -(rhoBC[k] + c1o1) * XXa * YYc * ZZc - c1o216;
         mfcca = -(rhoBC[k] + c1o1) * XXc * YYc * ZZa - c1o216;
         mfaaa = -(rhoBC[k] + c1o1) * XXa * YYa * ZZa - c1o216;
         mfcaa = -(rhoBC[k] + c1o1) * XXc * YYa * ZZa - c1o216;
         mfaca = -(rhoBC[k] + c1o1) * XXa * YYc * ZZa - c1o216;
         //////////////////////////////////////////////////////////
         ////konditioniert
         //double OneOver216RhoPlusOne = c1over216*(rhoBC[k]+one);
         //double OnoOver216Rho        = c1over216*rhoBC[k];
         //mfcbb = OnoOver216Rho*sixteen + OneOver216RhoPlusOne*twelve*(-(two*vy2) - two*vz2 + three*vy2*vz2 + vvx*(-two + three*vy2)*(-two + three*vz2) + vx2*(-two + three*vy2)*(-two + three*vz2));
         //mfabb = OnoOver216Rho*sixteen - OneOver216RhoPlusOne*twelve*(two*vy2 + two*vz2 - three*vy2*vz2 + vvx*(-two + three*vy2)*(-two + three*vz2) + vx2*(-four + six*vy2 + six*vz2 - nine*vy2*vz2));
         //mfbcb = four*(-(four*OneOver216RhoPlusOne) + four*OnoOver216Rho + OneOver216RhoPlusOne*(-two + three*vx2)*(one + three*vvy + three*vy2)*(-two + three*vz2));
         //mfbab = four*(four*OnoOver216Rho - OneOver216RhoPlusOne*three*(vvy*(-two + three*vx2)*(-two + three*vz2) - one*vx2*(one + three*vy2)*(-two + three*vz2) + two*(-(two*vy2) + vz2 + three*vy2*vz2)));
         //mfbbc = four*(-(four*OneOver216RhoPlusOne) + four*OnoOver216Rho + OneOver216RhoPlusOne*(-two + three*vx2)*(-two + three*vy2)*(one + three*vvz + three*vz2));
         //mfbba = four*(four*OnoOver216Rho - OneOver216RhoPlusOne*three*(vvz*(-two + three*vx2)*(-two + three*vy2) - one*vx2*(-two + three*vy2)*(one + three*vz2) + two*(vy2 - two*vz2 + three*vy2*vz2)));
         //mfccb = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) - two*vy2 - six*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(one + three*vvy + three*vy2)*(-two + three*vz2))));
         //mfaab = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) - two*vy2 - six*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-two + three*vz2))));
         //mfcab = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 + two*vy2 + six*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-two + three*vz2)));
         //mfacb = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 + two*vy2 + six*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(one + three*vvy + three*vy2)*(-two + three*vz2)));
         //mfcbc = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) + vy2 + three*vx2*vy2 + vvz*(one + three*vx2)*(-two + three*vy2) - two*vz2 - six*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(one + three*vvz + three*vz2))));
         //mfaba = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) + vy2 + three*vx2*vy2 - one*vvz*(one + three*vx2)*(-two + three*vy2) - two*vz2 - six*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(-one + three*vvz - three*vz2))));
         //mfcba = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 - one*vy2 - three*vx2*vy2 + vvz*(one + three*vx2)*(-two + three*vy2) + two*vz2 + six*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(-one + three*vvz - three*vz2)));
         //mfabc = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 - one*vy2 - three*vx2*vy2 - one*vvz*(one + three*vx2)*(-two + three*vy2) + two*vz2 + six*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(one + three*vvz + three*vz2)));
         //mfbcc = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(vx2 - two*vy2 + three*vx2*vy2 + vvz*(-two + three*vx2)*(one + three*vy2) - two*vz2 + three*vx2*vz2 - six*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(one + three*vvz + three*vz2))));
         //mfbaa = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(vx2 - two*vy2 + three*vx2*vy2 - one*vvz*(-two + three*vx2)*(one + three*vy2) - two*vz2 + three*vx2*vz2 - six*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(-one + three*vvz - three*vz2))));
         //mfbca = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(-(one*vx2) + two*vy2 - three*vx2*vy2 + vvz*(-two + three*vx2)*(one + three*vy2) + two*vz2 - three*vx2*vz2 + six*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(-one + three*vvz - three*vz2)));
         //mfbac = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(-(one*vx2) + two*vy2 - three*vx2*vy2 - one*vvz*(-two + three*vx2)*(one + three*vy2) + two*vz2 - three*vx2*vz2 + six*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(one + three*vvz + three*vz2)));
         //mfbbb = eight*(eight*OnoOver216Rho + OneOver216RhoPlusOne*three*(four*vy2 + four*vz2 - six*vy2*vz2 + vx2*(-two + three*vy2)*(-two + three*vz2)));
         //mfccc = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(one + three*vvz + three*vz2) + vvx*(one + three*vvy + three*vy2)*(one + three*vvz + three*vz2));
         //mfaac = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(one + three*vvz + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(one + three*vvz + three*vz2));
         //mfcac = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(one + three*vvz + three*vz2) - one*vvx*(-one + three*vvy - three*vy2)*(one + three*vvz + three*vz2));
         //mfacc = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(one + three*vvz + three*vz2) - one*vvx*(one + three*vvy + three*vy2)*(one + three*vvz + three*vz2));
         //mfcca = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) - one*vvx*(one + three*vvy + three*vy2)*(-one + three*vvz - three*vz2));
         //mfaaa = OnoOver216Rho - OneOver216RhoPlusOne*three*(vvz - one*vx2 + three*vvz*vx2 - one*vy2 + three*vvz*vy2 - three*vx2*vy2 + nine*vvz*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-one + three*vvz - three*vz2));
         //mfcaa = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-one + three*vvz - three*vz2));
         //mfaca = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(one + three*vvy + three*vy2)*(-one + three*vvz - three*vz2));

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //if (isEvenTimestep==true)
      //{
      //   D.f[dP00] = &DD[dP00 * size_Mat];
      //   D.f[dM00] = &DD[dM00 * size_Mat];
      //   D.f[d0P0] = &DD[d0P0 * size_Mat];
      //   D.f[d0M0] = &DD[d0M0 * size_Mat];
      //   D.f[d00P] = &DD[d00P * size_Mat];
      //   D.f[d00M] = &DD[d00M * size_Mat];
      //   D.f[dPP0] = &DD[dPP0 * size_Mat];
      //   D.f[dMM0] = &DD[dMM0 * size_Mat];
      //   D.f[dPM0] = &DD[dPM0 * size_Mat];
      //   D.f[dMP0] = &DD[dMP0 * size_Mat];
      //   D.f[dP0P] = &DD[dP0P * size_Mat];
      //   D.f[dM0M] = &DD[dM0M * size_Mat];
      //   D.f[dP0M] = &DD[dP0M * size_Mat];
      //   D.f[dM0P] = &DD[dM0P * size_Mat];
      //   D.f[d0PP] = &DD[d0PP * size_Mat];
      //   D.f[d0MM] = &DD[d0MM * size_Mat];
      //   D.f[d0PM] = &DD[d0PM * size_Mat];
      //   D.f[d0MP] = &DD[d0MP * size_Mat];
      //   D.f[d000] = &DD[d000 * size_Mat];
      //   D.f[dPPP] = &DD[dPPP * size_Mat];
      //   D.f[dMMP] = &DD[dMMP * size_Mat];
      //   D.f[dPMP] = &DD[dPMP * size_Mat];
      //   D.f[dMPP] = &DD[dMPP * size_Mat];
      //   D.f[dPPM] = &DD[dPPM * size_Mat];
      //   D.f[dMMM] = &DD[dMMM * size_Mat];
      //   D.f[dPMM] = &DD[dPMM * size_Mat];
      //   D.f[dMPM] = &DD[dMPM * size_Mat];
      //}
      //else
      //{
      //   D.f[dM00] = &DD[dP00 * size_Mat];
      //   D.f[dP00] = &DD[dM00 * size_Mat];
      //   D.f[d0M0] = &DD[d0P0 * size_Mat];
      //   D.f[d0P0] = &DD[d0M0 * size_Mat];
      //   D.f[d00M] = &DD[d00P * size_Mat];
      //   D.f[d00P] = &DD[d00M * size_Mat];
      //   D.f[dMM0] = &DD[dPP0 * size_Mat];
      //   D.f[dPP0] = &DD[dMM0 * size_Mat];
      //   D.f[dMP0] = &DD[dPM0 * size_Mat];
      //   D.f[dPM0] = &DD[dMP0 * size_Mat];
      //   D.f[dM0M] = &DD[dP0P * size_Mat];
      //   D.f[dP0P] = &DD[dM0M * size_Mat];
      //   D.f[dM0P] = &DD[dP0M * size_Mat];
      //   D.f[dP0M] = &DD[dM0P * size_Mat];
      //   D.f[d0MM] = &DD[d0PP * size_Mat];
      //   D.f[d0PP] = &DD[d0MM * size_Mat];
      //   D.f[d0MP] = &DD[d0PM * size_Mat];
      //   D.f[d0PM] = &DD[d0MP * size_Mat];
      //   D.f[d000] = &DD[d000 * size_Mat];
      //   D.f[dPPP] = &DD[dMMM * size_Mat];
      //   D.f[dMMP] = &DD[dPPM * size_Mat];
      //   D.f[dPMP] = &DD[dMPM * size_Mat];
      //   D.f[dMPP] = &DD[dPMM * size_Mat];
      //   D.f[dPPM] = &DD[dMMP * size_Mat];
      //   D.f[dMMM] = &DD[dPPP * size_Mat];
      //   D.f[dPMM] = &DD[dMPP * size_Mat];
      //   D.f[dMPM] = &DD[dPMP * size_Mat];
      //}
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();

         (D.f[dP00])[ke   ] = mfabb;//mfcbb;
         (D.f[dM00])[kw   ] = mfcbb;//mfabb;
         (D.f[d0P0])[kn   ] = mfbab;//mfbcb;
         (D.f[d0M0])[ks   ] = mfbcb;//mfbab;
         (D.f[d00P])[kt   ] = mfbba;//mfbbc;
         (D.f[d00M])[kb   ] = mfbbc;//mfbba;
         (D.f[dPP0])[kne  ] = mfaab;//mfccb;
         (D.f[dMM0])[ksw  ] = mfccb;//mfaab;
         (D.f[dPM0])[kse  ] = mfacb;//mfcab;
         (D.f[dMP0])[knw  ] = mfcab;//mfacb;
         (D.f[dP0P])[kte  ] = mfaba;//mfcbc;
         (D.f[dM0M])[kbw  ] = mfcbc;//mfaba;
         (D.f[dP0M])[kbe  ] = mfabc;//mfcba;
         (D.f[dM0P])[ktw  ] = mfcba;//mfabc;
         (D.f[d0PP])[ktn  ] = mfbaa;//mfbcc;
         (D.f[d0MM])[kbs  ] = mfbcc;//mfbaa;
         (D.f[d0PM])[kbn  ] = mfbac;//mfbca;
         (D.f[d0MP])[kts  ] = mfbca;//mfbac;
         (D.f[d000])[kzero] = mfbbb;//mfbbb;
         (D.f[dPPP])[ktne ] = mfaaa;//mfccc;
         (D.f[dMMP])[ktsw ] = mfcca;//mfaac;
         (D.f[dPMP])[ktse ] = mfaca;//mfcac;
         (D.f[dMPP])[ktnw ] = mfcaa;//mfacc;
         (D.f[dPPM])[kbne ] = mfaac;//mfcca;
         (D.f[dMMM])[kbsw ] = mfccc;//mfaaa;
         (D.f[dPMM])[kbse ] = mfacc;//mfcaa;
         (D.f[dMPM])[kbnw ] = mfcac;//mfaca;
         //(D.f[dP00])[ke   ] = mfcbb;
         //(D.f[dM00])[kw   ] = mfabb;
         //(D.f[d0P0])[kn   ] = mfbcb;
         //(D.f[d0M0])[ks   ] = mfbab;
         //(D.f[d00P])[kt   ] = mfbbc;
         //(D.f[d00M])[kb   ] = mfbba;
         //(D.f[dPP0])[kne  ] = mfccb;
         //(D.f[dMM0])[ksw  ] = mfaab;
         //(D.f[dPM0])[kse  ] = mfcab;
         //(D.f[dMP0])[knw  ] = mfacb;
         //(D.f[dP0P])[kte  ] = mfcbc;
         //(D.f[dM0M])[kbw  ] = mfaba;
         //(D.f[dP0M])[kbe  ] = mfcba;
         //(D.f[dM0P])[ktw  ] = mfabc;
         //(D.f[d0PP])[ktn  ] = mfbcc;
         //(D.f[d0MM])[kbs  ] = mfbaa;
         //(D.f[d0PM])[kbn  ] = mfbca;
         //(D.f[d0MP])[kts  ] = mfbac;
         //(D.f[d000])[kzero] = mfbbb;
         //(D.f[dPPP])[ktne ] = mfccc;
         //(D.f[dMMP])[ktsw ] = mfaac;
         //(D.f[dPMP])[ktse ] = mfcac;
         //(D.f[dMPP])[ktnw ] = mfacc;
         //(D.f[dPPM])[kbne ] = mfcca;
         //(D.f[dMMM])[kbsw ] = mfaaa;
         //(D.f[dPMM])[kbse ] = mfcaa;
         //(D.f[dMPM])[kbnw ] = mfaca;

      //(D.f[dP00])[ke   ] = fE ;  //f1_E ;   //fW;    //fE ;
      //(D.f[dM00])[kw   ] = fW ;  //f1_W ;   //fE;    //fW ;
      //(D.f[d0P0])[kn   ] = fN ;  //f1_N ;   //fS;    //fN ;
      //(D.f[d0M0])[ks   ] = fS ;  //f1_S ;   //fN;    //fS ;
      //(D.f[d00P])[kt   ] = fT ;  //f1_T ;   //fB;    //fT ;
      //(D.f[d00M])[kb   ] = fB ;  //f1_B ;   //fT;    //fB ;
      //(D.f[dPP0])[kne  ] = fNE;  //f1_NE;   //fSW;   //fNE;
      //(D.f[dMM0])[ksw  ] = fSW;  //f1_SW;   //fNE;   //fSW;
      //(D.f[dPM0])[kse  ] = fSE;  //f1_SE;   //fNW;   //fSE;
      //(D.f[dMP0])[knw  ] = fNW;  //f1_NW;   //fSE;   //fNW;
      //(D.f[dP0P])[kte  ] = fTE;  //f1_TE;   //fBW;   //fTE;
      //(D.f[dM0M])[kbw  ] = fBW;  //f1_BW;   //fTE;   //fBW;
      //(D.f[dP0M])[kbe  ] = fBE;  //f1_BE;   //fTW;   //fBE;
      //(D.f[dM0P])[ktw  ] = fTW;  //f1_TW;   //fBE;   //fTW;
      //(D.f[d0PP])[ktn  ] = fTN;  //f1_TN;   //fBS;   //fTN;
      //(D.f[d0MM])[kbs  ] = fBS;  //f1_BS;   //fTN;   //fBS;
      //(D.f[d0PM])[kbn  ] = fBN;  //f1_BN;   //fTS;   //fBN;
      //(D.f[d0MP])[kts  ] = fTS;  //f1_TS;   //fBN;   //fTS;
      //(D.f[d000])[kzero] = fZERO;//f1_ZERO; //fZERO; //fZERO;
      //(D.f[dPPP])[ktne ] = fTNE; //f1_TNE;  //fBSW;  //fTNE;
      //(D.f[dMMM])[kbsw ] = fBSW; //f1_BSW;  //fTNE;  //fBSW;
      //(D.f[dPPM])[kbne ] = fBNE; //f1_BNE;  //fTSW;  //fBNE;
      //(D.f[dMMP])[ktsw ] = fTSW; //f1_TSW;  //fBNE;  //fTSW;
      //(D.f[dPMP])[ktse ] = fTSE; //f1_TSE;  //fBNW;  //fTSE;
      //(D.f[dMPM])[kbnw ] = fBNW; //f1_BNW;  //fTSE;  //fBNW;
      //(D.f[dPMM])[kbse ] = fBSE; //f1_BSE;  //fTNW;  //fBSE;
      //(D.f[dMPP])[ktnw ] = fTNW; //f1_TNW;  //fBSE;  //fTNW;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceZero27(
    real* DD,
    int* k_Q,
    unsigned int numberOfBCnodes,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
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
      Distributions27 D;
      if (isEvenTimestep==false)
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
      //__syncthreads();
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      (D.f[dP00])[ke   ] =c0o1;
      (D.f[dM00])[kw   ] =c0o1;
      (D.f[d0P0])[kn   ] =c0o1;
      (D.f[d0M0])[ks   ] =c0o1;
      (D.f[d00P])[kt   ] =c0o1;
      (D.f[d00M])[kb   ] =c0o1;
      (D.f[dPP0])[kne  ] =c0o1;
      (D.f[dMM0])[ksw  ] =c0o1;
      (D.f[dPM0])[kse  ] =c0o1;
      (D.f[dMP0])[knw  ] =c0o1;
      (D.f[dP0P])[kte  ] =c0o1;
      (D.f[dM0M])[kbw  ] =c0o1;
      (D.f[dP0M])[kbe  ] =c0o1;
      (D.f[dM0P])[ktw  ] =c0o1;
      (D.f[d0PP])[ktn  ] =c0o1;
      (D.f[d0MM])[kbs  ] =c0o1;
      (D.f[d0PM])[kbn  ] =c0o1;
      (D.f[d0MP])[kts  ] =c0o1;
      (D.f[d000])[kzero] =c0o1;
      (D.f[dPPP])[ktne ] =c0o1;
      (D.f[dMMP])[ktsw ] =c0o1;
      (D.f[dPMP])[ktse ] =c0o1;
      (D.f[dMPP])[ktnw ] =c0o1;
      (D.f[dPPM])[kbne ] =c0o1;
      (D.f[dMMM])[kbsw ] =c0o1;
      (D.f[dPMM])[kbse ] =c0o1;
      (D.f[dMPM])[kbnw ] =c0o1;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceFake27(
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
    bool isEvenTimestep)
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
      if (isEvenTimestep==false)
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

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3;
      vx1    =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                  (f1_E - f1_W);


      vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                  (f1_N - f1_S);

      vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                  (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                  (f1_T - f1_B);

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

     //drho1 = (drho1 + rhoBC[k])/2.f;
     drho1 = drho1 - rhoBC[k];

      __syncthreads();

      (D.f[dP00])[ke   ] = c2o27* (rhoBC[k]+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      (D.f[dM00])[kw   ] = c2o27* (rhoBC[k]+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      (D.f[d0P0])[kn   ] = c2o27* (rhoBC[k]+c3o1*(    -vx2    )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      (D.f[d0M0])[ks   ] = c2o27* (rhoBC[k]+c3o1*(     vx2    )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      (D.f[d00P])[kt   ] = c2o27* (rhoBC[k]+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      (D.f[d00M])[kb   ] = c2o27* (rhoBC[k]+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      (D.f[dPP0])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dMM0])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dPM0])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dMP0])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dP0P])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dM0M])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dP0M])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[dM0P])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0PP])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0MM])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0PM])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d0MP])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[d000])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[dPPP])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMMP])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPMP])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMPP])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPPM])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMMM])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dPMM])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[dMPM])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////
__global__ void QPressDevice27_IntBB(
    real* rho,
    real* DD,
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

   if(k < numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //real VeloX = vx[k];
      //real VeloY = vy[k];
      //real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
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
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
         f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);

      vx1    = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
         ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
         (f_E - f_W))/(c1o1+drho);


      vx2    =  ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
         ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
         (f_N - f_S))/(c1o1+drho);

      vx3    =  (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
         (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
         (f_T - f_B))/(c1o1+drho);

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
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
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[d000])[k]=c1o10;
      real rhoDiff = drho - rho[k];
      real VeloX = vx1;
      real VeloY = vx2;
      real VeloZ = vx3;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c2o27*(rhoDiff + c6o1*( VeloX     )))/(c1o1+q);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c2o27*(rhoDiff + c6o1*(-VeloX     )))/(c1o1+q);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[d0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c2o27*(rhoDiff + c6o1*( VeloY     )))/(c1o1+q);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[d0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c2o27*(rhoDiff + c6o1*(-VeloY     )))/(c1o1+q);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D.f[d00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c2o27*(rhoDiff + c6o1*( VeloZ     )))/(c1o1+q);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[d00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c2o27*(rhoDiff + c6o1*(-VeloZ     )))/(c1o1+q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[dMM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c1o54*(rhoDiff + c6o1*(VeloX+VeloY)))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[dPP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloY)))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[dMP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloY)))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[dPM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloY)))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[dM0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c1o54*(rhoDiff + c6o1*( VeloX+VeloZ)))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[dP0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloZ)))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[dM0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloZ)))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[dP0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloZ)))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[d0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c1o54*(rhoDiff + c6o1*( VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[d0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c1o54*(rhoDiff + c6o1*( -VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[d0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c1o54*(rhoDiff + c6o1*( VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[d0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c1o54*(rhoDiff + c6o1*( -VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[dMMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[dPPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[dMMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[dPPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[dMPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[dPMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[dMPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         (D.f[dPMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY+VeloZ)))/(c1o1+q);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


