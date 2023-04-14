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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f1_E    = (D.f[DIR_P00])[k1e   ];
      real f1_W    = (D.f[DIR_M00])[k1w   ];
      real f1_N    = (D.f[DIR_0P0])[k1n   ];
      real f1_S    = (D.f[DIR_0M0])[k1s   ];
      real f1_T    = (D.f[DIR_00P])[k1t   ];
      real f1_B    = (D.f[DIR_00M])[k1b   ];
      real f1_NE   = (D.f[DIR_PP0])[k1ne  ];
      real f1_SW   = (D.f[DIR_MM0])[k1sw  ];
      real f1_SE   = (D.f[DIR_PM0])[k1se  ];
      real f1_NW   = (D.f[DIR_MP0])[k1nw  ];
      real f1_TE   = (D.f[DIR_P0P])[k1te  ];
      real f1_BW   = (D.f[DIR_M0M])[k1bw  ];
      real f1_BE   = (D.f[DIR_P0M])[k1be  ];
      real f1_TW   = (D.f[DIR_M0P])[k1tw  ];
      real f1_TN   = (D.f[DIR_0PP])[k1tn  ];
      real f1_BS   = (D.f[DIR_0MM])[k1bs  ];
      real f1_BN   = (D.f[DIR_0PM])[k1bn  ];
      real f1_TS   = (D.f[DIR_0MP])[k1ts  ];
      //real f1_ZERO = (D.f[DIR_000])[k1zero];
      real f1_TNE  = (D.f[DIR_PPP])[k1tne ];
      real f1_TSW  = (D.f[DIR_MMP])[k1tsw ];
      real f1_TSE  = (D.f[DIR_PMP])[k1tse ];
      real f1_TNW  = (D.f[DIR_MPP])[k1tnw ];
      real f1_BNE  = (D.f[DIR_PPM])[k1bne ];
      real f1_BSW  = (D.f[DIR_MMM])[k1bsw ];
      real f1_BSE  = (D.f[DIR_PMM])[k1bse ];
      real f1_BNW  = (D.f[DIR_MPM])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E    = (D.f[DIR_P00])[ke   ];
      real f_W    = (D.f[DIR_M00])[kw   ];
      real f_N    = (D.f[DIR_0P0])[kn   ];
      real f_S    = (D.f[DIR_0M0])[ks   ];
      real f_T    = (D.f[DIR_00P])[kt   ];
      real f_B    = (D.f[DIR_00M])[kb   ];
      real f_NE   = (D.f[DIR_PP0])[kne  ];
      real f_SW   = (D.f[DIR_MM0])[ksw  ];
      real f_SE   = (D.f[DIR_PM0])[kse  ];
      real f_NW   = (D.f[DIR_MP0])[knw  ];
      real f_TE   = (D.f[DIR_P0P])[kte  ];
      real f_BW   = (D.f[DIR_M0M])[kbw  ];
      real f_BE   = (D.f[DIR_P0M])[kbe  ];
      real f_TW   = (D.f[DIR_M0P])[ktw  ];
      real f_TN   = (D.f[DIR_0PP])[ktn  ];
      real f_BS   = (D.f[DIR_0MM])[kbs  ];
      real f_BN   = (D.f[DIR_0PM])[kbn  ];
      real f_TS   = (D.f[DIR_0MP])[kts  ];
      //real f_ZERO = (D.f[DIR_000])[kzero];
      real f_TNE  = (D.f[DIR_PPP])[ktne ];
      real f_TSW  = (D.f[DIR_MMP])[ktsw ];
      real f_TSE  = (D.f[DIR_PMP])[ktse ];
      real f_TNW  = (D.f[DIR_MPP])[ktnw ];
      real f_BNE  = (D.f[DIR_PPM])[kbne ];
      real f_BSW  = (D.f[DIR_MMM])[kbsw ];
      real f_BSE  = (D.f[DIR_PMM])[kbse ];
      real f_BNW  = (D.f[DIR_MPM])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      // real vx1, vx2, vx3;
      real drho, drho1;
      //////////////////////////////////////////////////////////////////////////
     //Dichte
      drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW +
                f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[DIR_000])[k1zero]);
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
     // -X
     //(D.f[DIR_P00])[ke   ] = f_E   ;
     //(D.f[DIR_PM0])[kse  ] = f_SE  ;
     //(D.f[DIR_PP0])[kne  ] = f_NE  ;
     //(D.f[DIR_P0M])[kbe  ] = f_BE  ;
     //(D.f[DIR_P0P])[kte  ] = f_TE  ;
     //(D.f[DIR_PMP])[ktse ] = f_TSE ;
     //(D.f[DIR_PPP])[ktne ] = f_TNE ;
     //(D.f[DIR_PMM])[kbse ] = f_BSE ;
     //(D.f[DIR_PPM])[kbne ] = f_BNE ;
     // X
     (D.f[DIR_M00])[kw   ] = f_W   ;
     (D.f[DIR_MM0])[ksw  ] = f_SW  ;
     (D.f[DIR_MP0])[knw  ] = f_NW  ;
     (D.f[DIR_M0M])[kbw  ] = f_BW  ;
     (D.f[DIR_M0P])[ktw  ] = f_TW  ;
     (D.f[DIR_MMP])[ktsw ] = f_TSW ;
     (D.f[DIR_MPP])[ktnw ] = f_TNW ;
     (D.f[DIR_MMM])[kbsw ] = f_BSW ;
     (D.f[DIR_MPM])[kbnw ] = f_BNW ;
     // Y
     //(D.f[DIR_0M0])[ks   ] = f_S   ;
     //(D.f[DIR_PM0])[kse  ] = f_SE  ;
     //(D.f[DIR_MM0])[ksw  ] = f_SW  ;
     //(D.f[DIR_0MP])[kts  ] = f_TS  ;
     //(D.f[DIR_0MM])[kbs  ] = f_BS  ;
     //(D.f[DIR_PMP])[ktse ] = f_TSE ;
     //(D.f[DIR_MMP])[ktsw ] = f_TSW ;
     //(D.f[DIR_PMM])[kbse ] = f_BSE ;
     //(D.f[DIR_MMM])[kbsw ] = f_BSW ;
     // Z
     //(D.f[DIR_00M])[kb   ] = f_B   ;
     //(D.f[DIR_P0M])[kbe  ] = f_BE  ;
     //(D.f[DIR_M0M])[kbw  ] = f_BW  ;
     //(D.f[DIR_0PM])[kbn  ] = f_BN  ;
     //(D.f[DIR_0MM])[kbs  ] = f_BS  ;
     //(D.f[DIR_PPM])[kbne ] = f_BNE ;
     //(D.f[DIR_MPM])[kbnw ] = f_BNW ;
     //(D.f[DIR_PMM])[kbse ] = f_BSE ;
     //(D.f[DIR_MMM])[kbsw ] = f_BSW ;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[DIR_P00])[k1e   ];
      f1_E    = (D.f[DIR_M00])[k1w   ];
      f1_S    = (D.f[DIR_0P0])[k1n   ];
      f1_N    = (D.f[DIR_0M0])[k1s   ];
      f1_B    = (D.f[DIR_00P])[k1t   ];
      f1_T    = (D.f[DIR_00M])[k1b   ];
      f1_SW   = (D.f[DIR_PP0])[k1ne  ];
      f1_NE   = (D.f[DIR_MM0])[k1sw  ];
      f1_NW   = (D.f[DIR_PM0])[k1se  ];
      f1_SE   = (D.f[DIR_MP0])[k1nw  ];
      f1_BW   = (D.f[DIR_P0P])[k1te  ];
      f1_TE   = (D.f[DIR_M0M])[k1bw  ];
      f1_TW   = (D.f[DIR_P0M])[k1be  ];
      f1_BE   = (D.f[DIR_M0P])[k1tw  ];
      f1_BS   = (D.f[DIR_0PP])[k1tn  ];
      f1_TN   = (D.f[DIR_0MM])[k1bs  ];
      f1_TS   = (D.f[DIR_0PM])[k1bn  ];
      f1_BN   = (D.f[DIR_0MP])[k1ts  ];
      f1_ZERO = (D.f[DIR_000])[k1zero];
      f1_BSW  = (D.f[DIR_PPP])[k1tne ];
      f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
      f1_BNW  = (D.f[DIR_PMP])[k1tse ];
      f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
      f1_TSW  = (D.f[DIR_PPM])[k1bne ];
      f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
      f1_TNW  = (D.f[DIR_PMM])[k1bse ];
      f1_TSE  = (D.f[DIR_MPM])[k1bnw ];

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

      (D.f[DIR_P00])[ke   ] = f1_W   ;
      (D.f[DIR_M00])[kw   ] = f1_E   ;
      (D.f[DIR_0P0])[kn   ] = f1_S   ;
      (D.f[DIR_0M0])[ks   ] = f1_N   ;
      (D.f[DIR_00P])[kt   ] = f1_B   ;
      (D.f[DIR_00M])[kb   ] = f1_T   ;
      (D.f[DIR_PP0])[kne  ] = f1_SW  ;
      (D.f[DIR_MM0])[ksw  ] = f1_NE  ;
      (D.f[DIR_PM0])[kse  ] = f1_NW  ;
      (D.f[DIR_MP0])[knw  ] = f1_SE  ;
      (D.f[DIR_P0P])[kte  ] = f1_BW  ;
      (D.f[DIR_M0M])[kbw  ] = f1_TE  ;
      (D.f[DIR_P0M])[kbe  ] = f1_TW  ;
      (D.f[DIR_M0P])[ktw  ] = f1_BE  ;
      (D.f[DIR_0PP])[ktn  ] = f1_BS  ;
      (D.f[DIR_0MM])[kbs  ] = f1_TN  ;
      (D.f[DIR_0PM])[kbn  ] = f1_TS  ;
      (D.f[DIR_0MP])[kts  ] = f1_BN  ;
      (D.f[DIR_000])[kzero] = f1_ZERO;
      (D.f[DIR_PPP])[ktne ] = f1_BSW ;
      (D.f[DIR_MMP])[ktsw ] = f1_BNE ;
      (D.f[DIR_PMP])[ktse ] = f1_BNW ;
      (D.f[DIR_MPP])[ktnw ] = f1_BSE ;
      (D.f[DIR_PPM])[kbne ] = f1_TSW ;
      (D.f[DIR_MMM])[kbsw ] = f1_TNE ;
      (D.f[DIR_PMM])[kbse ] = f1_TNW ;
      (D.f[DIR_MPM])[kbnw ] = f1_TSE ;
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
      real f1_W    = (dist.f[DIR_P00])[k1e   ];
      real f1_E    = (dist.f[DIR_M00])[k1w   ];
      real f1_S    = (dist.f[DIR_0P0])[k1n   ];
      real f1_N    = (dist.f[DIR_0M0])[k1s   ];
      real f1_B    = (dist.f[DIR_00P])[k1t   ];
      real f1_T    = (dist.f[DIR_00M])[k1b   ];
      real f1_SW   = (dist.f[DIR_PP0])[k1ne  ];
      real f1_NE   = (dist.f[DIR_MM0])[k1sw  ];
      real f1_NW   = (dist.f[DIR_PM0])[k1se  ];
      real f1_SE   = (dist.f[DIR_MP0])[k1nw  ];
      real f1_BW   = (dist.f[DIR_P0P])[k1te  ];
      real f1_TE   = (dist.f[DIR_M0M])[k1bw  ];
      real f1_TW   = (dist.f[DIR_P0M])[k1be  ];
      real f1_BE   = (dist.f[DIR_M0P])[k1tw  ];
      real f1_BS   = (dist.f[DIR_0PP])[k1tn  ];
      real f1_TN   = (dist.f[DIR_0MM])[k1bs  ];
      real f1_TS   = (dist.f[DIR_0PM])[k1bn  ];
      real f1_BN   = (dist.f[DIR_0MP])[k1ts  ];
      real f1_ZERO = (dist.f[DIR_000])[k1zero];
      real f1_BSW  = (dist.f[DIR_PPP])[k1tne ];
      real f1_BNE  = (dist.f[DIR_MMP])[k1tsw ];
      real f1_BNW  = (dist.f[DIR_PMP])[k1tse ];
      real f1_BSE  = (dist.f[DIR_MPP])[k1tnw ];
      real f1_TSW  = (dist.f[DIR_PPM])[k1bne ];
      real f1_TNE  = (dist.f[DIR_MMM])[k1bsw ];
      real f1_TNW  = (dist.f[DIR_PMM])[k1bse ];
      real f1_TSE  = (dist.f[DIR_MPM])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities (for neighboring node)
      //!
      real drho1 = f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                   f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW +
                   f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((dist.f[DIR_000])[kzero]);

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
      (dist.f[DIR_P00])[ke   ] = f1_W   ;
      (dist.f[DIR_M00])[kw   ] = f1_E   ;
      (dist.f[DIR_0P0])[kn   ] = f1_S   ;
      (dist.f[DIR_0M0])[ks   ] = f1_N   ;
      (dist.f[DIR_00P])[kt   ] = f1_B   ;
      (dist.f[DIR_00M])[kb   ] = f1_T   ;
      (dist.f[DIR_PP0])[kne  ] = f1_SW  ;
      (dist.f[DIR_MM0])[ksw  ] = f1_NE  ;
      (dist.f[DIR_PM0])[kse  ] = f1_NW  ;
      (dist.f[DIR_MP0])[knw  ] = f1_SE  ;
      (dist.f[DIR_P0P])[kte  ] = f1_BW  ;
      (dist.f[DIR_M0M])[kbw  ] = f1_TE  ;
      (dist.f[DIR_P0M])[kbe  ] = f1_TW  ;
      (dist.f[DIR_M0P])[ktw  ] = f1_BE  ;
      (dist.f[DIR_0PP])[ktn  ] = f1_BS  ;
      (dist.f[DIR_0MM])[kbs  ] = f1_TN  ;
      (dist.f[DIR_0PM])[kbn  ] = f1_TS  ;
      (dist.f[DIR_0MP])[kts  ] = f1_BN  ;
      (dist.f[DIR_000])[kzero] = f1_ZERO;
      (dist.f[DIR_PPP])[ktne ] = f1_BSW ;
      (dist.f[DIR_MMP])[ktsw ] = f1_BNE ;
      (dist.f[DIR_PMP])[ktse ] = f1_BNW ;
      (dist.f[DIR_MPP])[ktnw ] = f1_BSE ;
      (dist.f[DIR_PPM])[kbne ] = f1_TSW ;
      (dist.f[DIR_MMM])[kbsw ] = f1_TNE ;
      (dist.f[DIR_PMM])[kbse ] = f1_TNW ;
      (dist.f[DIR_MPM])[kbnw ] = f1_TSE ;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
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

      f1_W    = (D.f[DIR_P00])[k1e   ];
      f1_E    = (D.f[DIR_M00])[k1w   ];
      f1_S    = (D.f[DIR_0P0])[k1n   ];
      f1_N    = (D.f[DIR_0M0])[k1s   ];
      f1_B    = (D.f[DIR_00P])[k1t   ];
      f1_T    = (D.f[DIR_00M])[k1b   ];
      f1_SW   = (D.f[DIR_PP0])[k1ne  ];
      f1_NE   = (D.f[DIR_MM0])[k1sw  ];
      f1_NW   = (D.f[DIR_PM0])[k1se  ];
      f1_SE   = (D.f[DIR_MP0])[k1nw  ];
      f1_BW   = (D.f[DIR_P0P])[k1te  ];
      f1_TE   = (D.f[DIR_M0M])[k1bw  ];
      f1_TW   = (D.f[DIR_P0M])[k1be  ];
      f1_BE   = (D.f[DIR_M0P])[k1tw  ];
      f1_BS   = (D.f[DIR_0PP])[k1tn  ];
      f1_TN   = (D.f[DIR_0MM])[k1bs  ];
      f1_TS   = (D.f[DIR_0PM])[k1bn  ];
      f1_BN   = (D.f[DIR_0MP])[k1ts  ];
      f1_ZERO = (D.f[DIR_000])[k1zero];
      f1_BSW  = (D.f[DIR_PPP])[k1tne ];
      f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
      f1_BNW  = (D.f[DIR_PMP])[k1tse ];
      f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
      f1_TSW  = (D.f[DIR_PPM])[k1bne ];
      f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
      f1_TNW  = (D.f[DIR_PMM])[k1bse ];
      f1_TSE  = (D.f[DIR_MPM])[k1bnw ];

      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                        f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      (D.f[DIR_P00])[ke   ] = f1_W   -c2o27*drho1;
      (D.f[DIR_M00])[kw   ] = f1_E   -c2o27*drho1;
      (D.f[DIR_0P0])[kn   ] = f1_S   -c2o27*drho1;
      (D.f[DIR_0M0])[ks   ] = f1_N   -c2o27*drho1;
      (D.f[DIR_00P])[kt   ] = f1_B   -c2o27*drho1;
      (D.f[DIR_00M])[kb   ] = f1_T   -c2o27*drho1;
      (D.f[DIR_PP0])[kne  ] = f1_SW  -c1o54*drho1;
      (D.f[DIR_MM0])[ksw  ] = f1_NE  -c1o54*drho1;
      (D.f[DIR_PM0])[kse  ] = f1_NW  -c1o54*drho1;
      (D.f[DIR_MP0])[knw  ] = f1_SE  -c1o54*drho1;
      (D.f[DIR_P0P])[kte  ] = f1_BW  -c1o54*drho1;
      (D.f[DIR_M0M])[kbw  ] = f1_TE  -c1o54*drho1;
      (D.f[DIR_P0M])[kbe  ] = f1_TW  -c1o54*drho1;
      (D.f[DIR_M0P])[ktw  ] = f1_BE  -c1o54*drho1;
      (D.f[DIR_0PP])[ktn  ] = f1_BS  -c1o54*drho1;
      (D.f[DIR_0MM])[kbs  ] = f1_TN  -c1o54*drho1;
      (D.f[DIR_0PM])[kbn  ] = f1_TS  -c1o54*drho1;
      (D.f[DIR_0MP])[kts  ] = f1_BN  -c1o54*drho1;
      (D.f[DIR_000])[kzero] = f1_ZERO-c8o27*drho1;
      (D.f[DIR_PPP])[ktne ] = f1_BSW -c1o216*drho1;
      (D.f[DIR_MMP])[ktsw ] = f1_BNE -c1o216*drho1;
      (D.f[DIR_PMP])[ktse ] = f1_BNW -c1o216*drho1;
      (D.f[DIR_MPP])[ktnw ] = f1_BSE -c1o216*drho1;
      (D.f[DIR_PPM])[kbne ] = f1_TSW -c1o216*drho1;
      (D.f[DIR_MMM])[kbsw ] = f1_TNE -c1o216*drho1;
      (D.f[DIR_PMM])[kbse ] = f1_TNW -c1o216*drho1;
      (D.f[DIR_MPM])[kbnw ] = f1_TSE -c1o216*drho1;
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
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   }
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M00])[kw]=c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         //(D.f[DIR_P00])[ke]=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P00])[ke]=c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         //(D.f[DIR_M00])[kw]=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0M0])[ks]=c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         //(D.f[DIR_0P0])[kn]=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0P0])[kn]=c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         //(D.f[DIR_0M0])[ks]=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00M])[kb]=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         //(D.f[DIR_00P])[kt]=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00P])[kt]=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         //(D.f[DIR_00M])[kb]=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MM0])[ksw]=c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         //(D.f[DIR_PP0])[kne]=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PP0])[kne]=c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         //(D.f[DIR_MM0])[ksw]=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MP0])[knw]=c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         //(D.f[DIR_PM0])[kse]=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PM0])[kse]=c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         //(D.f[DIR_MP0])[knw]=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0M])[kbw]=c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         //(D.f[DIR_P0P])[kte]=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0P])[kte]=c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         //(D.f[DIR_M0M])[kbw]=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0P])[ktw]=c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         //(D.f[DIR_P0M])[kbe]=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0M])[kbe]=c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         //(D.f[DIR_M0P])[ktw]=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MM])[kbs]=c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         //(D.f[DIR_0PP])[ktn]=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PP])[ktn]=c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         //(D.f[DIR_0MM])[kbs]=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MP])[kts]=c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         //(D.f[DIR_0PM])[kbn]=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PM])[kbn]=c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         //(D.f[DIR_0MP])[kts]=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMM])[kbsw]=c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         //(D.f[DIR_PPP])[ktne]=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPP])[ktne]=c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         //(D.f[DIR_MMM])[kbsw]=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMP])[ktsw]=c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         //(D.f[DIR_PPM])[kbne]=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPM])[kbne]=c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         //(D.f[DIR_MMP])[ktsw]=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPM])[kbnw]=c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         //(D.f[DIR_PMP])[ktse]=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMP])[ktse]=c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         //(D.f[DIR_MPM])[kbnw]=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPP])[ktnw]=c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         //(D.f[DIR_PMM])[kbse]=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMM])[kbse]=c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         //(D.f[DIR_MPP])[ktnw]=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
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
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   }
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
      f_ZERO = (D.f[DIR_000])[kzero];
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M00])[kw]=f_W-c2o27*drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P00])[ke]=f_E-c2o27*drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0M0])[ks]=f_S-c2o27*drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0P0])[kn]=f_N-c2o27*drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00M])[kb]=f_B-c2o27*drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_00P])[kt]=f_T-c2o27*drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MM0])[ksw]=f_SW-c1o54*drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PP0])[kne]=f_NE-c1o54*drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MP0])[knw]=f_NW-c1o54*drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PM0])[kse]=f_SE-c1o54*drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0M])[kbw]=f_BW-c1o54*drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0P])[kte]=f_TE-c1o54*drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_M0P])[ktw]=f_TW-c1o54*drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_P0M])[kbe]=f_BE-c1o54*drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MM])[kbs]=f_BS-c1o54*drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PP])[ktn]=f_TN-c1o54*drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0MP])[kts]=f_TS-c1o54*drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_0PM])[kbn]=f_BN-c1o54*drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMM])[kbsw]=f_BSW-c1o216*drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPP])[ktne]=f_TNE-c1o216*drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MMP])[ktsw]=f_TSW-c1o216*drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PPM])[kbne]=f_BNE-c1o216*drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPM])[kbnw]=f_BNW-c1o216*drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMP])[ktse]=f_TSE-c1o216*drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_MPP])[ktnw]=f_TNW-c1o216*drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[DIR_PMM])[kbse]=f_BSE-c1o216*drho;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         (D.f[DIR_M00])[kw]       = c2o27  * deltaRho;
         (D.f[DIR_P00])[ke]       = c2o27  * deltaRho;
         (D.f[DIR_0M0])[ks]       = c2o27  * deltaRho;
         (D.f[DIR_0P0])[kn]       = c2o27  * deltaRho;
         (D.f[DIR_00M])[kb]       = c2o27  * deltaRho;
         (D.f[DIR_00P])[kt]       = c2o27  * deltaRho;
         (D.f[DIR_MM0])[ksw]     = c1o54  * deltaRho;
         (D.f[DIR_PP0])[kne]     = c1o54  * deltaRho;
         (D.f[DIR_MP0])[knw]     = c1o54  * deltaRho;
         (D.f[DIR_PM0])[kse]     = c1o54  * deltaRho;
         (D.f[DIR_M0M])[kbw]     = c1o54  * deltaRho;
         (D.f[DIR_P0P])[kte]     = c1o54  * deltaRho;
         (D.f[DIR_M0P])[ktw]     = c1o54  * deltaRho;
         (D.f[DIR_P0M])[kbe]     = c1o54  * deltaRho;
         (D.f[DIR_0MM])[kbs]     = c1o54  * deltaRho;
         (D.f[DIR_0PP])[ktn]     = c1o54  * deltaRho;
         (D.f[DIR_0MP])[kts]     = c1o54  * deltaRho;
         (D.f[DIR_0PM])[kbn]     = c1o54  * deltaRho;
         (D.f[DIR_MMM])[kbsw]   = c1o216 * deltaRho;
         (D.f[DIR_PPP])[ktne]   = c1o216 * deltaRho;
         (D.f[DIR_MMP])[ktsw]   = c1o216 * deltaRho;
         (D.f[DIR_PPM])[kbne]   = c1o216 * deltaRho;
         (D.f[DIR_MPM])[kbnw]   = c1o216 * deltaRho;
         (D.f[DIR_PMP])[ktse]   = c1o216 * deltaRho;
         (D.f[DIR_MPP])[ktnw]   = c1o216 * deltaRho;
         (D.f[DIR_PMM])[kbse]   = c1o216 * deltaRho;
         (D.f[DIR_000])[kzero] = c8o27  * deltaRho;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E,f_W,f_N,f_S,f_T,f_NE,f_SW,f_SE,f_NW,f_TE,f_TW,f_TN,f_TS,f_ZERO,f_TNE,f_TSW,f_TSE,f_TNW;//,
            //f_B,f_BW,f_BE,f_BS,f_BN,f_BSW,f_BNE,f_BNW,f_BSE;

      f_E    = (D.f[DIR_P00])[ke   ];
      f_W    = (D.f[DIR_M00])[kw   ];
      f_N    = (D.f[DIR_0P0])[kn   ];
      f_S    = (D.f[DIR_0M0])[ks   ];
      f_T    = (D.f[DIR_00P])[kt   ];
      f_NE   = (D.f[DIR_PP0])[kne  ];
      f_SW   = (D.f[DIR_MM0])[ksw  ];
      f_SE   = (D.f[DIR_PM0])[kse  ];
      f_NW   = (D.f[DIR_MP0])[knw  ];
      f_TE   = (D.f[DIR_P0P])[kte  ];
      f_TW   = (D.f[DIR_M0P])[ktw  ];
      f_TN   = (D.f[DIR_0PP])[ktn  ];
      f_TS   = (D.f[DIR_0MP])[kts  ];
      f_ZERO = (D.f[DIR_000])[kzero];
      f_TNE  = (D.f[DIR_PPP])[ktne ];
      f_TSW  = (D.f[DIR_MMP])[ktsw ];
      f_TSE  = (D.f[DIR_PMP])[ktse ];
      f_TNW  = (D.f[DIR_MPP])[ktnw ];
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

      //(D.f[DIR_000])[kzero] = c8over27*  (drho-cusq);
      //(D.f[DIR_P00])[ke]    = c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
      //(D.f[DIR_M00])[kw]    = c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
      //(D.f[DIR_0P0])[kn]     = c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
      //(D.f[DIR_0M0])[ks]    = c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
      //(D.f[DIR_00P])[kt]    = c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
      //(D.f[DIR_00M])[kb]    = c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
      //(D.f[DIR_PP0])[kne]   = c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
      //(D.f[DIR_MM0])[ksw]   = c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
      //(D.f[DIR_PM0])[kse]   =  c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
      //(D.f[DIR_MP0])[knw]   =  c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
      //(D.f[DIR_P0P])[kte]   =  c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
      //(D.f[DIR_M0M])[kbw]   =  c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
      //(D.f[DIR_P0M])[kbe]   =  c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
      //(D.f[DIR_M0P])[ktw]   =  c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
      //(D.f[DIR_0PP])[ktn]   =  c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
      //(D.f[DIR_0MM])[kbs]   =  c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
      //(D.f[DIR_0PM])[kbn]   =  c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
      //(D.f[DIR_0MP])[kts]   =  c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
      //(D.f[DIR_PPP])[ktne]  =  c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
      //(D.f[DIR_MMM])[kbsw]  =  c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
      //(D.f[DIR_PPM])[kbne]  =  c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
      //(D.f[DIR_MMP])[ktsw]  =  c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
      //(D.f[DIR_PMP])[ktse]  =  c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
      //(D.f[DIR_MPM])[kbnw]  =  c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
      //(D.f[DIR_PMM])[kbse]  =  c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
      //(D.f[DIR_MPP])[ktnw]  =  c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
      real drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_NE+f_SW+f_SE+f_NW+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      real dTop   =    f_T+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      (D.f[DIR_00M])[kb]     = (f_T+c2o27)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c2o27;
      (D.f[DIR_M0M])[kbw]   = (f_TW+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[DIR_P0M])[kbe]   = (f_TE+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[DIR_0MM])[kbs]   = (f_TS+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[DIR_0PM])[kbn]   = (f_TN+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[DIR_MMM])[kbsw] = (f_TSW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[DIR_PPM])[kbne] = (f_TNE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[DIR_MPM])[kbnw] = (f_TNW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[DIR_PMM])[kbse] = (f_TSE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
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
   f1[DIR_P00] = (dist.f[DIR_P00])[k1e   ];
   f1[DIR_M00] = (dist.f[DIR_M00])[k1w   ];
   f1[DIR_0P0] = (dist.f[DIR_0P0])[k1n   ];
   f1[DIR_0M0] = (dist.f[DIR_0M0])[k1s   ];
   f1[DIR_00P] = (dist.f[DIR_00P])[k1t   ];
   f1[DIR_00M] = (dist.f[DIR_00M])[k1b   ];
   f1[DIR_PP0] = (dist.f[DIR_PP0])[k1ne  ];
   f1[DIR_MM0] = (dist.f[DIR_MM0])[k1sw  ];
   f1[DIR_PM0] = (dist.f[DIR_PM0])[k1se  ];
   f1[DIR_MP0] = (dist.f[DIR_MP0])[k1nw  ];
   f1[DIR_P0P] = (dist.f[DIR_P0P])[k1te  ];
   f1[DIR_M0M] = (dist.f[DIR_M0M])[k1bw  ];
   f1[DIR_P0M] = (dist.f[DIR_P0M])[k1be  ];
   f1[DIR_M0P] = (dist.f[DIR_M0P])[k1tw  ];
   f1[DIR_0PP] = (dist.f[DIR_0PP])[k1tn  ];
   f1[DIR_0MM] = (dist.f[DIR_0MM])[k1bs  ];
   f1[DIR_0PM] = (dist.f[DIR_0PM])[k1bn  ];
   f1[DIR_0MP] = (dist.f[DIR_0MP])[k1ts  ];
   // f1[DIR_000] = (dist.f[DIR_000])[k1zero];
   f1[DIR_PPP] = (dist.f[DIR_PPP])[k1tne ];
   f1[DIR_MMP] = (dist.f[DIR_MMP])[k1tsw ];
   f1[DIR_PMP] = (dist.f[DIR_PMP])[k1tse ];
   f1[DIR_MPP] = (dist.f[DIR_MPP])[k1tnw ];
   f1[DIR_PPM] = (dist.f[DIR_PPM])[k1bne ];
   f1[DIR_MMM] = (dist.f[DIR_MMM])[k1bsw ];
   f1[DIR_PMM] = (dist.f[DIR_PMM])[k1bse ];
   f1[DIR_MPM] = (dist.f[DIR_MPM])[k1bnw ];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   f[DIR_P00] = (dist.f[DIR_P00])[ke   ];
   f[DIR_M00] = (dist.f[DIR_M00])[kw   ];
   f[DIR_0P0] = (dist.f[DIR_0P0])[kn   ];
   f[DIR_0M0] = (dist.f[DIR_0M0])[ks   ];
   f[DIR_00P] = (dist.f[DIR_00P])[kt   ];
   f[DIR_00M] = (dist.f[DIR_00M])[kb   ];
   f[DIR_PP0] = (dist.f[DIR_PP0])[kne  ];
   f[DIR_MM0] = (dist.f[DIR_MM0])[ksw  ];
   f[DIR_PM0] = (dist.f[DIR_PM0])[kse  ];
   f[DIR_MP0] = (dist.f[DIR_MP0])[knw  ];
   f[DIR_P0P] = (dist.f[DIR_P0P])[kte  ];
   f[DIR_M0M] = (dist.f[DIR_M0M])[kbw  ];
   f[DIR_P0M] = (dist.f[DIR_P0M])[kbe  ];
   f[DIR_M0P] = (dist.f[DIR_M0P])[ktw  ];
   f[DIR_0PP] = (dist.f[DIR_0PP])[ktn  ];
   f[DIR_0MM] = (dist.f[DIR_0MM])[kbs  ];
   f[DIR_0PM] = (dist.f[DIR_0PM])[kbn  ];
   f[DIR_0MP] = (dist.f[DIR_0MP])[kts  ];
   // f[DIR_000] = (dist.f[DIR_000])[kzero];
   f[DIR_PPP] = (dist.f[DIR_PPP])[ktne ];
   f[DIR_MMP] = (dist.f[DIR_MMP])[ktsw ];
   f[DIR_PMP] = (dist.f[DIR_PMP])[ktse ];
   f[DIR_MPP] = (dist.f[DIR_MPP])[ktnw ];
   f[DIR_PPM] = (dist.f[DIR_PPM])[kbne ];
   f[DIR_MMM] = (dist.f[DIR_MMM])[kbsw ];
   f[DIR_PMM] = (dist.f[DIR_PMM])[kbse ];
   f[DIR_MPM] = (dist.f[DIR_MPM])[kbnw ];
   //////////////////////////////////////////////////////////////////////////


   real cs = c1o1 / sqrtf(c3o1);

   //////////////////////////////////////////////////////////////////////////
   getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);
   switch(direction)
   {
      case MZZ:
         (dist.f[DIR_P00])[ke   ] = computeOutflowDistribution(f, f1, DIR_P00, cs);
         (dist.f[DIR_PM0])[kse  ] = computeOutflowDistribution(f, f1, DIR_PM0, cs);
         (dist.f[DIR_PP0])[kne  ] = computeOutflowDistribution(f, f1, DIR_PP0, cs);
         (dist.f[DIR_P0M])[kbe  ] = computeOutflowDistribution(f, f1, DIR_P0M, cs);
         (dist.f[DIR_P0P])[kte  ] = computeOutflowDistribution(f, f1, DIR_P0P, cs);
         (dist.f[DIR_PMP])[ktse ] = computeOutflowDistribution(f, f1, DIR_PMP, cs);
         (dist.f[DIR_PPP])[ktne ] = computeOutflowDistribution(f, f1, DIR_PPP, cs);
         (dist.f[DIR_PMM])[kbse ] = computeOutflowDistribution(f, f1, DIR_PMM, cs);
         (dist.f[DIR_PPM])[kbne ] = computeOutflowDistribution(f, f1, DIR_PPM, cs);
         break;

      case PZZ:
         (dist.f[DIR_M00])[kw   ] = computeOutflowDistribution(f, f1, DIR_M00, cs);
         (dist.f[DIR_MM0])[ksw  ] = computeOutflowDistribution(f, f1, DIR_MM0, cs);
         (dist.f[DIR_MP0])[knw  ] = computeOutflowDistribution(f, f1, DIR_MP0, cs);
         (dist.f[DIR_M0M])[kbw  ] = computeOutflowDistribution(f, f1, DIR_M0M, cs);
         (dist.f[DIR_M0P])[ktw  ] = computeOutflowDistribution(f, f1, DIR_M0P, cs);
         (dist.f[DIR_MMP])[ktsw ] = computeOutflowDistribution(f, f1, DIR_MMP, cs);
         (dist.f[DIR_MPP])[ktnw ] = computeOutflowDistribution(f, f1, DIR_MPP, cs);
         (dist.f[DIR_MMM])[kbsw ] = computeOutflowDistribution(f, f1, DIR_MMM, cs);
         (dist.f[DIR_MPM])[kbnw ] = computeOutflowDistribution(f, f1, DIR_MPM, cs);
         break;

      case ZMZ:
         (dist.f[DIR_0P0])[kn   ] = computeOutflowDistribution(f, f1, DIR_0P0, cs);
         (dist.f[DIR_PP0])[kne  ] = computeOutflowDistribution(f, f1, DIR_PP0, cs);
         (dist.f[DIR_MP0])[knw  ] = computeOutflowDistribution(f, f1, DIR_MP0, cs);
         (dist.f[DIR_0PP])[ktn  ] = computeOutflowDistribution(f, f1, DIR_0PP, cs);
         (dist.f[DIR_0PM])[kbn  ] = computeOutflowDistribution(f, f1, DIR_0PM, cs);
         (dist.f[DIR_PPP])[ktne ] = computeOutflowDistribution(f, f1, DIR_PPP, cs);
         (dist.f[DIR_MPP])[ktnw ] = computeOutflowDistribution(f, f1, DIR_MPP, cs);
         (dist.f[DIR_PPM])[kbne ] = computeOutflowDistribution(f, f1, DIR_PPM, cs);
         (dist.f[DIR_MPM])[kbnw ] = computeOutflowDistribution(f, f1, DIR_MPM, cs);
         break;

      case ZPZ:
         (dist.f[DIR_0M0])[ks   ] = computeOutflowDistribution(f, f1, DIR_0M0, cs);
         (dist.f[DIR_PM0])[kse  ] = computeOutflowDistribution(f, f1, DIR_PM0, cs);
         (dist.f[DIR_MM0])[ksw  ] = computeOutflowDistribution(f, f1, DIR_MM0, cs);
         (dist.f[DIR_0MP])[kts  ] = computeOutflowDistribution(f, f1, DIR_0MP, cs);
         (dist.f[DIR_0MM])[kbs  ] = computeOutflowDistribution(f, f1, DIR_0MM, cs);
         (dist.f[DIR_PMP])[ktse ] = computeOutflowDistribution(f, f1, DIR_PMP, cs);
         (dist.f[DIR_MMP])[ktsw ] = computeOutflowDistribution(f, f1, DIR_MMP, cs);
         (dist.f[DIR_PMM])[kbse ] = computeOutflowDistribution(f, f1, DIR_PMM, cs);
         (dist.f[DIR_MMM])[kbsw ] = computeOutflowDistribution(f, f1, DIR_MMM, cs);
         break;

      case ZZM:
         (dist.f[DIR_00P])[kt   ] = computeOutflowDistribution(f, f1, DIR_00P, cs);
         (dist.f[DIR_P0P])[kte  ] = computeOutflowDistribution(f, f1, DIR_P0P, cs);
         (dist.f[DIR_M0P])[ktw  ] = computeOutflowDistribution(f, f1, DIR_M0P, cs);
         (dist.f[DIR_0PP])[ktn  ] = computeOutflowDistribution(f, f1, DIR_0PP, cs);
         (dist.f[DIR_0MP])[kts  ] = computeOutflowDistribution(f, f1, DIR_0MP, cs);
         (dist.f[DIR_PPP])[ktne ] = computeOutflowDistribution(f, f1, DIR_PPP, cs);
         (dist.f[DIR_MPP])[ktnw ] = computeOutflowDistribution(f, f1, DIR_MPP, cs);
         (dist.f[DIR_PMP])[ktse ] = computeOutflowDistribution(f, f1, DIR_PMP, cs);
         (dist.f[DIR_MMP])[ktsw ] = computeOutflowDistribution(f, f1, DIR_MMP, cs);
         break;

      case ZZP:
         (dist.f[DIR_00M])[kb   ] = computeOutflowDistribution(f, f1, DIR_00M, cs);
         (dist.f[DIR_P0M])[kbe  ] = computeOutflowDistribution(f, f1, DIR_P0M, cs);
         (dist.f[DIR_M0M])[kbw  ] = computeOutflowDistribution(f, f1, DIR_M0M, cs);
         (dist.f[DIR_0PM])[kbn  ] = computeOutflowDistribution(f, f1, DIR_0PM, cs);
         (dist.f[DIR_0MM])[kbs  ] = computeOutflowDistribution(f, f1, DIR_0MM, cs);
         (dist.f[DIR_PPM])[kbne ] = computeOutflowDistribution(f, f1, DIR_PPM, cs);
         (dist.f[DIR_MPM])[kbnw ] = computeOutflowDistribution(f, f1, DIR_MPM, cs);
         (dist.f[DIR_PMM])[kbse ] = computeOutflowDistribution(f, f1, DIR_PMM, cs);
         (dist.f[DIR_MMM])[kbsw ] = computeOutflowDistribution(f, f1, DIR_MMM, cs);
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
   f[DIR_000] = (dist.f[DIR_000])[k_000];
   f[DIR_P00] = (dist.f[DIR_P00])[k_000];
   f[DIR_M00] = (dist.f[DIR_M00])[k_M00];
   f[DIR_0P0] = (dist.f[DIR_0P0])[k_000];
   f[DIR_0M0] = (dist.f[DIR_0M0])[k_0M0];
   f[DIR_00P] = (dist.f[DIR_00P])[k_000];
   f[DIR_00M] = (dist.f[DIR_00M])[k_00M];
   f[DIR_PP0] = (dist.f[DIR_PP0])[k_000];
   f[DIR_MM0] = (dist.f[DIR_MM0])[k_MM0];
   f[DIR_PM0] = (dist.f[DIR_PM0])[k_0M0];
   f[DIR_MP0] = (dist.f[DIR_MP0])[k_M00];
   f[DIR_P0P] = (dist.f[DIR_P0P])[k_000];
   f[DIR_M0M] = (dist.f[DIR_M0M])[k_M0M];
   f[DIR_P0M] = (dist.f[DIR_P0M])[k_00M];
   f[DIR_M0P] = (dist.f[DIR_M0P])[k_M00];
   f[DIR_0PP] = (dist.f[DIR_0PP])[k_000];
   f[DIR_0MM] = (dist.f[DIR_0MM])[k_0MM];
   f[DIR_0PM] = (dist.f[DIR_0PM])[k_00M];
   f[DIR_0MP] = (dist.f[DIR_0MP])[k_0M0];
   f[DIR_PPP] = (dist.f[DIR_PPP])[k_000];
   f[DIR_MPP] = (dist.f[DIR_MPP])[k_M00];
   f[DIR_PMP] = (dist.f[DIR_PMP])[k_0M0];
   f[DIR_MMP] = (dist.f[DIR_MMP])[k_MM0];
   f[DIR_PPM] = (dist.f[DIR_PPM])[k_00M];
   f[DIR_MPM] = (dist.f[DIR_MPM])[k_M0M];
   f[DIR_PMM] = (dist.f[DIR_PMM])[k_0MM];
   f[DIR_MMM] = (dist.f[DIR_MMM])[k_MMM];
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   fN[DIR_000] = (dist.f[DIR_000])[kN_000];
   fN[DIR_P00] = (dist.f[DIR_P00])[kN_000];
   fN[DIR_M00] = (dist.f[DIR_M00])[kN_M00];
   fN[DIR_0P0] = (dist.f[DIR_0P0])[kN_000];
   fN[DIR_0M0] = (dist.f[DIR_0M0])[kN_0M0];
   fN[DIR_00P] = (dist.f[DIR_00P])[kN_000];
   fN[DIR_00M] = (dist.f[DIR_00M])[kN_00M];
   fN[DIR_PP0] = (dist.f[DIR_PP0])[kN_000];
   fN[DIR_MM0] = (dist.f[DIR_MM0])[kN_MM0];
   fN[DIR_PM0] = (dist.f[DIR_PM0])[kN_0M0];
   fN[DIR_MP0] = (dist.f[DIR_MP0])[kN_M00];
   fN[DIR_P0P] = (dist.f[DIR_P0P])[kN_000];
   fN[DIR_M0M] = (dist.f[DIR_M0M])[kN_M0M];
   fN[DIR_P0M] = (dist.f[DIR_P0M])[kN_00M];
   fN[DIR_M0P] = (dist.f[DIR_M0P])[kN_M00];
   fN[DIR_0PP] = (dist.f[DIR_0PP])[kN_000];
   fN[DIR_0MM] = (dist.f[DIR_0MM])[kN_0MM];
   fN[DIR_0PM] = (dist.f[DIR_0PM])[kN_00M];
   fN[DIR_0MP] = (dist.f[DIR_0MP])[kN_0M0];
   fN[DIR_PPP] = (dist.f[DIR_PPP])[kN_000];
   fN[DIR_MPP] = (dist.f[DIR_MPP])[kN_M00];
   fN[DIR_PMP] = (dist.f[DIR_PMP])[kN_0M0];
   fN[DIR_MMP] = (dist.f[DIR_MMP])[kN_MM0];
   fN[DIR_PPM] = (dist.f[DIR_PPM])[kN_00M];
   fN[DIR_MPM] = (dist.f[DIR_MPM])[kN_M0M];
   fN[DIR_PMM] = (dist.f[DIR_PMM])[kN_0MM];
   fN[DIR_MMM] = (dist.f[DIR_MMM])[kN_MMM];
   //////////////////////////////////////////////////////////////////////////
   real drho = vf::lbm::getDensity(f);

   real rhoCorrection = densityCorrectionFactor*drho;

   real cs = c1o1 / sqrtf(c3o1);

   getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

   switch(direction)
   {
      case MZZ:
         (dist.f[DIR_P00])[k_000] = computeOutflowDistribution(f, fN, DIR_P00  , rhoCorrection, cs, c2o27);
         (dist.f[DIR_PM0])[k_0M0] = computeOutflowDistribution(f, fN, DIR_PM0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PP0])[k_000] = computeOutflowDistribution(f, fN, DIR_PP0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_P0M])[k_00M] = computeOutflowDistribution(f, fN, DIR_P0M, rhoCorrection, cs, c1o54);
         (dist.f[DIR_P0P])[k_000] = computeOutflowDistribution(f, fN, DIR_P0P, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PMP])[k_0M0] = computeOutflowDistribution(f, fN, DIR_PMP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PPP])[k_000] = computeOutflowDistribution(f, fN, DIR_PPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PMM])[k_0MM] = computeOutflowDistribution(f, fN, DIR_PMM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PPM])[k_00M] = computeOutflowDistribution(f, fN, DIR_PPM, rhoCorrection, cs, c1o216);
         break;

      case PZZ:
         (dist.f[DIR_M00])[k_M00] = computeOutflowDistribution(f, fN, DIR_M00, rhoCorrection, cs, c2o27);
         (dist.f[DIR_MM0])[k_MM0] = computeOutflowDistribution(f, fN, DIR_MM0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_MP0])[k_M00] = computeOutflowDistribution(f, fN, DIR_MP0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_M0M])[k_M0M] = computeOutflowDistribution(f, fN, DIR_M0M, rhoCorrection, cs, c1o54);
         (dist.f[DIR_M0P])[k_M00] = computeOutflowDistribution(f, fN, DIR_M0P, rhoCorrection, cs, c1o54);
         (dist.f[DIR_MMP])[k_MM0] = computeOutflowDistribution(f, fN, DIR_MMP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPP])[k_M00] = computeOutflowDistribution(f, fN, DIR_MPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MMM])[k_MMM] = computeOutflowDistribution(f, fN, DIR_MMM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPM])[k_M0M] = computeOutflowDistribution(f, fN, DIR_MPM, rhoCorrection, cs, c1o216);
         break;

      case ZMZ:
         (dist.f[DIR_0P0])[k_000] = computeOutflowDistribution(f, fN, DIR_0P0, rhoCorrection, cs, c2o27);
         (dist.f[DIR_PP0])[k_000] = computeOutflowDistribution(f, fN, DIR_PP0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_MP0])[k_M00] = computeOutflowDistribution(f, fN, DIR_MP0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0PP])[k_000] = computeOutflowDistribution(f, fN, DIR_0PP, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0PM])[k_00M] = computeOutflowDistribution(f, fN, DIR_0PM, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PPP])[k_000] = computeOutflowDistribution(f, fN, DIR_PPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPP])[k_M00] = computeOutflowDistribution(f, fN, DIR_MPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PPM])[k_00M] = computeOutflowDistribution(f, fN, DIR_PPM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPM])[k_M0M] = computeOutflowDistribution(f, fN, DIR_MPM, rhoCorrection, cs, c1o216);
         break;

      case ZPZ:
         (dist.f[DIR_0M0])[k_0M0] =computeOutflowDistribution(f, fN, DIR_0M0, rhoCorrection, cs, c2o27);
         (dist.f[DIR_PM0])[k_0M0] =computeOutflowDistribution(f, fN, DIR_PM0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_MM0])[k_MM0] =computeOutflowDistribution(f, fN, DIR_MM0, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0MP])[k_0M0] =computeOutflowDistribution(f, fN, DIR_0MP, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0MM])[k_0MM] =computeOutflowDistribution(f, fN, DIR_0MM, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PMP])[k_0M0] =computeOutflowDistribution(f, fN, DIR_PMP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MMP])[k_MM0] =computeOutflowDistribution(f, fN, DIR_MMP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PMM])[k_0MM] =computeOutflowDistribution(f, fN, DIR_PMM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MMM])[k_MMM] =computeOutflowDistribution(f, fN, DIR_MMM, rhoCorrection, cs, c1o216);
         break;

      case ZZM:
         (dist.f[DIR_00P])[k_000] = computeOutflowDistribution(f, fN, DIR_00P, rhoCorrection, cs, c2o27);
         (dist.f[DIR_P0P])[k_000] = computeOutflowDistribution(f, fN, DIR_P0P, rhoCorrection, cs, c1o54);
         (dist.f[DIR_M0P])[k_M00] = computeOutflowDistribution(f, fN, DIR_M0P, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0PP])[k_000] = computeOutflowDistribution(f, fN, DIR_0PP, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0MP])[k_0M0] = computeOutflowDistribution(f, fN, DIR_0MP, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PPP])[k_000] = computeOutflowDistribution(f, fN, DIR_PPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPP])[k_M00] = computeOutflowDistribution(f, fN, DIR_MPP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PMP])[k_0M0] = computeOutflowDistribution(f, fN, DIR_PMP, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MMP])[k_MM0] = computeOutflowDistribution(f, fN, DIR_MMP, rhoCorrection, cs, c1o216);
         break;

      case ZZP:
         (dist.f[DIR_00M])[k_00M] = computeOutflowDistribution(f, fN, DIR_00M, rhoCorrection, cs, c2o27);
         (dist.f[DIR_P0M])[k_00M] = computeOutflowDistribution(f, fN, DIR_P0M, rhoCorrection, cs, c1o54);
         (dist.f[DIR_M0M])[k_M0M] = computeOutflowDistribution(f, fN, DIR_M0M, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0PM])[k_00M] = computeOutflowDistribution(f, fN, DIR_0PM, rhoCorrection, cs, c1o54);
         (dist.f[DIR_0MM])[k_0MM] = computeOutflowDistribution(f, fN, DIR_0MM, rhoCorrection, cs, c1o54);
         (dist.f[DIR_PPM])[k_00M] = computeOutflowDistribution(f, fN, DIR_PPM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MPM])[k_M0M] = computeOutflowDistribution(f, fN, DIR_MPM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_PMM])[k_0MM] = computeOutflowDistribution(f, fN, DIR_PMM, rhoCorrection, cs, c1o216);
         (dist.f[DIR_MMM])[k_MMM] = computeOutflowDistribution(f, fN, DIR_MMM, rhoCorrection, cs, c1o216);
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[DIR_P00])[k1e   ];
      f1_E    = (D.f[DIR_M00])[k1w   ];
      f1_S    = (D.f[DIR_0P0])[k1n   ];
      f1_N    = (D.f[DIR_0M0])[k1s   ];
      f1_B    = (D.f[DIR_00P])[k1t   ];
      f1_T    = (D.f[DIR_00M])[k1b   ];
      f1_SW   = (D.f[DIR_PP0])[k1ne  ];
      f1_NE   = (D.f[DIR_MM0])[k1sw  ];
      f1_NW   = (D.f[DIR_PM0])[k1se  ];
      f1_SE   = (D.f[DIR_MP0])[k1nw  ];
      f1_BW   = (D.f[DIR_P0P])[k1te  ];
      f1_TE   = (D.f[DIR_M0M])[k1bw  ];
      f1_TW   = (D.f[DIR_P0M])[k1be  ];
      f1_BE   = (D.f[DIR_M0P])[k1tw  ];
      f1_BS   = (D.f[DIR_0PP])[k1tn  ];
      f1_TN   = (D.f[DIR_0MM])[k1bs  ];
      f1_TS   = (D.f[DIR_0PM])[k1bn  ];
      f1_BN   = (D.f[DIR_0MP])[k1ts  ];
      f1_ZERO = (D.f[DIR_000])[k1zero];
      f1_BSW  = (D.f[DIR_PPP])[k1tne ];
      f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
      f1_BNW  = (D.f[DIR_PMP])[k1tse ];
      f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
      f1_TSW  = (D.f[DIR_PPM])[k1bne ];
      f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
      f1_TNW  = (D.f[DIR_PMM])[k1bse ];
      f1_TSE  = (D.f[DIR_MPM])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

     //drho1 = (drho1 + rhoBC[k])/2.f;
     drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[DIR_P00])[ke   ] = f1_W   -c2o27*drho1;   //  c1o100;  // zero;  //
      (D.f[DIR_M00])[kw   ] = f1_E   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0P0])[kn   ] = f1_S   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0M0])[ks   ] = f1_N   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_00P])[kt   ] = f1_B   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_00M])[kb   ] = f1_T   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PP0])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MM0])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PM0])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MP0])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_P0P])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_M0M])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_P0M])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_M0P])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0PP])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0MM])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0PM])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0MP])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_000])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PPP])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MMP])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PMP])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MPP])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PPM])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MMM])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PMM])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MPM])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////
    //   Distributions27 kDistTest;
    //      kDistTest.f[DIR_P00] = &kTestRE[DIR_P00 * numberOfBCnodes];
    //      kDistTest.f[DIR_M00] = &kTestRE[DIR_M00 * numberOfBCnodes];
    //      kDistTest.f[DIR_0P0] = &kTestRE[DIR_0P0 * numberOfBCnodes];
    //      kDistTest.f[DIR_0M0] = &kTestRE[DIR_0M0 * numberOfBCnodes];
    //      kDistTest.f[DIR_00P] = &kTestRE[DIR_00P * numberOfBCnodes];
    //      kDistTest.f[DIR_00M] = &kTestRE[DIR_00M * numberOfBCnodes];
    //      kDistTest.f[DIR_PP0] = &kTestRE[DIR_PP0 * numberOfBCnodes];
    //      kDistTest.f[DIR_MM0] = &kTestRE[DIR_MM0 * numberOfBCnodes];
    //      kDistTest.f[DIR_PM0] = &kTestRE[DIR_PM0 * numberOfBCnodes];
    //      kDistTest.f[DIR_MP0] = &kTestRE[DIR_MP0 * numberOfBCnodes];
    //      kDistTest.f[DIR_P0P] = &kTestRE[DIR_P0P * numberOfBCnodes];
    //      kDistTest.f[DIR_M0M] = &kTestRE[DIR_M0M * numberOfBCnodes];
    //      kDistTest.f[DIR_P0M] = &kTestRE[DIR_P0M * numberOfBCnodes];
    //      kDistTest.f[DIR_M0P] = &kTestRE[DIR_M0P * numberOfBCnodes];
    //      kDistTest.f[DIR_0PP] = &kTestRE[DIR_0PP * numberOfBCnodes];
    //      kDistTest.f[DIR_0MM] = &kTestRE[DIR_0MM * numberOfBCnodes];
    //      kDistTest.f[DIR_0PM] = &kTestRE[DIR_0PM * numberOfBCnodes];
    //      kDistTest.f[DIR_0MP] = &kTestRE[DIR_0MP * numberOfBCnodes];
    //      kDistTest.f[DIR_000] = &kTestRE[DIR_000 * numberOfBCnodes];
    //      kDistTest.f[DIR_PPP] = &kTestRE[DIR_PPP * numberOfBCnodes];
    //      kDistTest.f[DIR_MMP] = &kTestRE[DIR_MMP * numberOfBCnodes];
    //      kDistTest.f[DIR_PMP] = &kTestRE[DIR_PMP * numberOfBCnodes];
    //      kDistTest.f[DIR_MPP] = &kTestRE[DIR_MPP * numberOfBCnodes];
    //      kDistTest.f[DIR_PPM] = &kTestRE[DIR_PPM * numberOfBCnodes];
    //      kDistTest.f[DIR_MMM] = &kTestRE[DIR_MMM * numberOfBCnodes];
    //      kDistTest.f[DIR_PMM] = &kTestRE[DIR_PMM * numberOfBCnodes];
    //      kDistTest.f[DIR_MPM] = &kTestRE[DIR_MPM * numberOfBCnodes];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   //f1_W    = (D.f[DIR_P00])[k1e   ];
   //   //f1_E    = (D.f[DIR_M00])[k1w   ];
   //   //f1_S    = (D.f[DIR_0P0])[k1n   ];
   //   //f1_N    = (D.f[DIR_0M0])[k1s   ];
   //   //f1_B    = (D.f[DIR_00P])[k1t   ];
   //   //f1_T    = (D.f[DIR_00M])[k1b   ];
   //   //f1_SW   = (D.f[DIR_PP0])[k1ne  ];
   //   //f1_NE   = (D.f[DIR_MM0])[k1sw  ];
   //   //f1_NW   = (D.f[DIR_PM0])[k1se  ];
   //   //f1_SE   = (D.f[DIR_MP0])[k1nw  ];
   //   //f1_BW   = (D.f[DIR_P0P])[k1te  ];
   //   //f1_TE   = (D.f[DIR_M0M])[k1bw  ];
   //   //f1_TW   = (D.f[DIR_P0M])[k1be  ];
   //   //f1_BE   = (D.f[DIR_M0P])[k1tw  ];
   //   //f1_BS   = (D.f[DIR_0PP])[k1tn  ];
   //   //f1_TN   = (D.f[DIR_0MM])[k1bs  ];
   //   //f1_TS   = (D.f[DIR_0PM])[k1bn  ];
   //   //f1_BN   = (D.f[DIR_0MP])[k1ts  ];
   //   //f1_ZERO = (D.f[DIR_000])[k1zero];
   //   //f1_BSW  = (D.f[DIR_PPP])[k1tne ];
   //   //f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
   //   //f1_BNW  = (D.f[DIR_PMP])[k1tse ];
   //   //f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
   //   //f1_TSW  = (D.f[DIR_PPM])[k1bne ];
   //   //f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
   //   //f1_TNW  = (D.f[DIR_PMM])[k1bse ];
   //   //f1_TSE  = (D.f[DIR_MPM])[k1bnw ];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   f1_E    = (D.f[DIR_P00])[k1e   ];
   //   f1_W    = (D.f[DIR_M00])[k1w   ];
   //   f1_N    = (D.f[DIR_0P0])[k1n   ];
   //   f1_S    = (D.f[DIR_0M0])[k1s   ];
   //   f1_T    = (D.f[DIR_00P])[k1t   ];
   //   f1_B    = (D.f[DIR_00M])[k1b   ];
   //   f1_NE   = (D.f[DIR_PP0])[k1ne  ];
   //   f1_SW   = (D.f[DIR_MM0])[k1sw  ];
   //   f1_SE   = (D.f[DIR_PM0])[k1se  ];
   //   f1_NW   = (D.f[DIR_MP0])[k1nw  ];
   //   f1_TE   = (D.f[DIR_P0P])[k1te  ];
   //   f1_BW   = (D.f[DIR_M0M])[k1bw  ];
   //   f1_BE   = (D.f[DIR_P0M])[k1be  ];
   //   f1_TW   = (D.f[DIR_M0P])[k1tw  ];
   //   f1_TN   = (D.f[DIR_0PP])[k1tn  ];
   //   f1_BS   = (D.f[DIR_0MM])[k1bs  ];
   //   f1_BN   = (D.f[DIR_0PM])[k1bn  ];
   //   f1_TS   = (D.f[DIR_0MP])[k1ts  ];
   //   f1_ZERO = (D.f[DIR_000])[k1zero];
   //   f1_TNE  = (D.f[DIR_PPP])[k1tne ];
   //   f1_TSW  = (D.f[DIR_MMP])[k1tsw ];
   //   f1_TSE  = (D.f[DIR_PMP])[k1tse ];
   //   f1_TNW  = (D.f[DIR_MPP])[k1tnw ];
   //   f1_BNE  = (D.f[DIR_PPM])[k1bne ];
   //   f1_BSW  = (D.f[DIR_MMM])[k1bsw ];
   //   f1_BSE  = (D.f[DIR_PMM])[k1bse ];
   //   f1_BNW  = (D.f[DIR_MPM])[k1bnw ];
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
         //double mfabb = (D.f[DIR_P00])[k1e   ];
         //double mfcbb = (D.f[DIR_M00])[k1w   ];
         //double mfbab = (D.f[DIR_0P0])[k1n   ];
         //double mfbcb = (D.f[DIR_0M0])[k1s   ];
         //double mfbba = (D.f[DIR_00P])[k1t   ];
         //double mfbbc = (D.f[DIR_00M])[k1b   ];
         //double mfaab = (D.f[DIR_PP0])[k1ne  ];
         //double mfccb = (D.f[DIR_MM0])[k1sw  ];
         //double mfacb = (D.f[DIR_PM0])[k1se  ];
         //double mfcab = (D.f[DIR_MP0])[k1nw  ];
         //double mfaba = (D.f[DIR_P0P])[k1te  ];
         //double mfcbc = (D.f[DIR_M0M])[k1bw  ];
         //double mfabc = (D.f[DIR_P0M])[k1be  ];
         //double mfcba = (D.f[DIR_M0P])[k1tw  ];
         //double mfbaa = (D.f[DIR_0PP])[k1tn  ];
         //double mfbcc = (D.f[DIR_0MM])[k1bs  ];
         //double mfbac = (D.f[DIR_0PM])[k1bn  ];
         //double mfbca = (D.f[DIR_0MP])[k1ts  ];
         //double mfbbb = (D.f[DIR_000])[k1zero];
         //double mfaaa = (D.f[DIR_PPP])[k1tne ];
         //double mfcca = (D.f[DIR_MMP])[k1tsw ];
         //double mfaca = (D.f[DIR_PMP])[k1tse ];
         //double mfcaa = (D.f[DIR_MPP])[k1tnw ];
         //double mfaac = (D.f[DIR_PPM])[k1bne ];
         //double mfccc = (D.f[DIR_MMM])[k1bsw ];
         //double mfacc = (D.f[DIR_PMM])[k1bse ];
         //double mfcac = (D.f[DIR_MPM])[k1bnw ];
         real mfabb = (D.f[DIR_P00])[k1e   ];
         real mfcbb = (D.f[DIR_M00])[k1w   ];
         real mfbab = (D.f[DIR_0P0])[k1n   ];
         real mfbcb = (D.f[DIR_0M0])[k1s   ];
         real mfbba = (D.f[DIR_00P])[k1t   ];
         real mfbbc = (D.f[DIR_00M])[k1b   ];
         real mfaab = (D.f[DIR_PP0])[k1ne  ];
         real mfccb = (D.f[DIR_MM0])[k1sw  ];
         real mfacb = (D.f[DIR_PM0])[k1se  ];
         real mfcab = (D.f[DIR_MP0])[k1nw  ];
         real mfaba = (D.f[DIR_P0P])[k1te  ];
         real mfcbc = (D.f[DIR_M0M])[k1bw  ];
         real mfabc = (D.f[DIR_P0M])[k1be  ];
         real mfcba = (D.f[DIR_M0P])[k1tw  ];
         real mfbaa = (D.f[DIR_0PP])[k1tn  ];
         real mfbcc = (D.f[DIR_0MM])[k1bs  ];
         real mfbac = (D.f[DIR_0PM])[k1bn  ];
         real mfbca = (D.f[DIR_0MP])[k1ts  ];
         real mfbbb = (D.f[DIR_000])[k1zero];
         real mfaaa = (D.f[DIR_PPP])[k1tne ];
         real mfcca = (D.f[DIR_MMP])[k1tsw ];
         real mfaca = (D.f[DIR_PMP])[k1tse ];
         real mfcaa = (D.f[DIR_MPP])[k1tnw ];
         real mfaac = (D.f[DIR_PPM])[k1bne ];
         real mfccc = (D.f[DIR_MMM])[k1bsw ];
         real mfacc = (D.f[DIR_PMM])[k1bse ];
         real mfcac = (D.f[DIR_MPM])[k1bnw ];

         //real mfcbb = (D.f[DIR_P00])[ke   ];
         //real mfabb = (D.f[DIR_M00])[kw   ];
         //real mfbcb = (D.f[DIR_0P0])[kn   ];
         //real mfbab = (D.f[DIR_0M0])[ks   ];
         //real mfbbc = (D.f[DIR_00P])[kt   ];
         //real mfbba = (D.f[DIR_00M])[kb   ];
         //real mfccb = (D.f[DIR_PP0])[kne  ];
         //real mfaab = (D.f[DIR_MM0])[ksw  ];
         //real mfcab = (D.f[DIR_PM0])[kse  ];
         //real mfacb = (D.f[DIR_MP0])[knw  ];
         //real mfcbc = (D.f[DIR_P0P])[kte  ];
         //real mfaba = (D.f[DIR_M0M])[kbw  ];
         //real mfcba = (D.f[DIR_P0M])[kbe  ];
         //real mfabc = (D.f[DIR_M0P])[ktw  ];
         //real mfbcc = (D.f[DIR_0PP])[ktn  ];
         //real mfbaa = (D.f[DIR_0MM])[kbs  ];
         //real mfbca = (D.f[DIR_0PM])[kbn  ];
         //real mfbac = (D.f[DIR_0MP])[kts  ];
         //real mfbbb = (D.f[DIR_000])[kzero];
         //real mfccc = (D.f[DIR_PPP])[ktne ];
         //real mfaac = (D.f[DIR_MMP])[ktsw ];
         //real mfcac = (D.f[DIR_PMP])[ktse ];
         //real mfacc = (D.f[DIR_MPP])[ktnw ];
         //real mfcca = (D.f[DIR_PPM])[kbne ];
         //real mfaaa = (D.f[DIR_MMM])[kbsw ];
         //real mfcaa = (D.f[DIR_PMM])[kbse ];
         //real mfaca = (D.f[DIR_MPM])[kbnw ];
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
         //	(kDistTest.f[DIR_P00])[k] = mfabb;
         //	(kDistTest.f[DIR_M00])[k] = mfcbb;
         //	(kDistTest.f[DIR_0P0])[k] = mfbab;
         //	(kDistTest.f[DIR_0M0])[k] = mfbcb;
         //	(kDistTest.f[DIR_00P])[k] = mfbba;
         //	(kDistTest.f[DIR_00M])[k] = mfbbc;
         //	(kDistTest.f[DIR_PP0])[k] = mfaab;
         //	(kDistTest.f[DIR_MM0])[k] = mfccb;
         //	(kDistTest.f[DIR_PM0])[k] = mfacb;
         //	(kDistTest.f[DIR_MP0])[k] = mfcab;
         //	(kDistTest.f[DIR_P0P])[k] = mfaba;
         //	(kDistTest.f[DIR_M0M])[k] = mfcbc;
         //	(kDistTest.f[DIR_P0M])[k] = mfabc;
         //	(kDistTest.f[DIR_M0P])[k] = mfcba;
         //	(kDistTest.f[DIR_0PP])[k] = mfbaa;
         //	(kDistTest.f[DIR_0MM])[k] = mfbcc;
         //	(kDistTest.f[DIR_0PM])[k] = mfbac;
         //	(kDistTest.f[DIR_0MP])[k] = mfbca;
         //	(kDistTest.f[DIR_000])[k] = KQK;
         //	(kDistTest.f[DIR_PPP])[k] = mfaaa;
         //	(kDistTest.f[DIR_MMP])[k] = mfcca;
         //	(kDistTest.f[DIR_PMP])[k] = mfaca;
         //	(kDistTest.f[DIR_MPP])[k] = mfcaa;
         //	(kDistTest.f[DIR_PPM])[k] = mfaac;
         //	(kDistTest.f[DIR_MMM])[k] = mfccc;
         //	(kDistTest.f[DIR_PMM])[k] = mfacc;
         //	(kDistTest.f[DIR_MPM])[k] = mfcac;
         //}else{
         //	(kDistTest.f[DIR_P00])[k] = zero;
         //	(kDistTest.f[DIR_M00])[k] = zero;
         //	(kDistTest.f[DIR_0P0])[k] = zero;
         //	(kDistTest.f[DIR_0M0])[k] = zero;
         //	(kDistTest.f[DIR_00P])[k] = zero;
         //	(kDistTest.f[DIR_00M])[k] = zero;
         //	(kDistTest.f[DIR_PP0])[k] = zero;
         //	(kDistTest.f[DIR_MM0])[k] = zero;
         //	(kDistTest.f[DIR_PM0])[k] = zero;
         //	(kDistTest.f[DIR_MP0])[k] = zero;
         //	(kDistTest.f[DIR_P0P])[k] = zero;
         //	(kDistTest.f[DIR_M0M])[k] = zero;
         //	(kDistTest.f[DIR_P0M])[k] = zero;
         //	(kDistTest.f[DIR_M0P])[k] = zero;
         //	(kDistTest.f[DIR_0PP])[k] = zero;
         //	(kDistTest.f[DIR_0MM])[k] = zero;
         //	(kDistTest.f[DIR_0PM])[k] = zero;
         //	(kDistTest.f[DIR_0MP])[k] = zero;
         //	(kDistTest.f[DIR_000])[k] = zero;
         //	(kDistTest.f[DIR_PPP])[k] = zero;
         //	(kDistTest.f[DIR_MMP])[k] = zero;
         //	(kDistTest.f[DIR_PMP])[k] = zero;
         //	(kDistTest.f[DIR_MPP])[k] = zero;
         //	(kDistTest.f[DIR_PPM])[k] = zero;
         //	(kDistTest.f[DIR_MMM])[k] = zero;
         //	(kDistTest.f[DIR_PMM])[k] = zero;
         //	(kDistTest.f[DIR_MPM])[k] = zero;
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
      //   D.f[DIR_P00] = &DD[DIR_P00 * size_Mat];
      //   D.f[DIR_M00] = &DD[DIR_M00 * size_Mat];
      //   D.f[DIR_0P0] = &DD[DIR_0P0 * size_Mat];
      //   D.f[DIR_0M0] = &DD[DIR_0M0 * size_Mat];
      //   D.f[DIR_00P] = &DD[DIR_00P * size_Mat];
      //   D.f[DIR_00M] = &DD[DIR_00M * size_Mat];
      //   D.f[DIR_PP0] = &DD[DIR_PP0 * size_Mat];
      //   D.f[DIR_MM0] = &DD[DIR_MM0 * size_Mat];
      //   D.f[DIR_PM0] = &DD[DIR_PM0 * size_Mat];
      //   D.f[DIR_MP0] = &DD[DIR_MP0 * size_Mat];
      //   D.f[DIR_P0P] = &DD[DIR_P0P * size_Mat];
      //   D.f[DIR_M0M] = &DD[DIR_M0M * size_Mat];
      //   D.f[DIR_P0M] = &DD[DIR_P0M * size_Mat];
      //   D.f[DIR_M0P] = &DD[DIR_M0P * size_Mat];
      //   D.f[DIR_0PP] = &DD[DIR_0PP * size_Mat];
      //   D.f[DIR_0MM] = &DD[DIR_0MM * size_Mat];
      //   D.f[DIR_0PM] = &DD[DIR_0PM * size_Mat];
      //   D.f[DIR_0MP] = &DD[DIR_0MP * size_Mat];
      //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
      //   D.f[DIR_PPP] = &DD[DIR_PPP * size_Mat];
      //   D.f[DIR_MMP] = &DD[DIR_MMP * size_Mat];
      //   D.f[DIR_PMP] = &DD[DIR_PMP * size_Mat];
      //   D.f[DIR_MPP] = &DD[DIR_MPP * size_Mat];
      //   D.f[DIR_PPM] = &DD[DIR_PPM * size_Mat];
      //   D.f[DIR_MMM] = &DD[DIR_MMM * size_Mat];
      //   D.f[DIR_PMM] = &DD[DIR_PMM * size_Mat];
      //   D.f[DIR_MPM] = &DD[DIR_MPM * size_Mat];
      //}
      //else
      //{
      //   D.f[DIR_M00] = &DD[DIR_P00 * size_Mat];
      //   D.f[DIR_P00] = &DD[DIR_M00 * size_Mat];
      //   D.f[DIR_0M0] = &DD[DIR_0P0 * size_Mat];
      //   D.f[DIR_0P0] = &DD[DIR_0M0 * size_Mat];
      //   D.f[DIR_00M] = &DD[DIR_00P * size_Mat];
      //   D.f[DIR_00P] = &DD[DIR_00M * size_Mat];
      //   D.f[DIR_MM0] = &DD[DIR_PP0 * size_Mat];
      //   D.f[DIR_PP0] = &DD[DIR_MM0 * size_Mat];
      //   D.f[DIR_MP0] = &DD[DIR_PM0 * size_Mat];
      //   D.f[DIR_PM0] = &DD[DIR_MP0 * size_Mat];
      //   D.f[DIR_M0M] = &DD[DIR_P0P * size_Mat];
      //   D.f[DIR_P0P] = &DD[DIR_M0M * size_Mat];
      //   D.f[DIR_M0P] = &DD[DIR_P0M * size_Mat];
      //   D.f[DIR_P0M] = &DD[DIR_M0P * size_Mat];
      //   D.f[DIR_0MM] = &DD[DIR_0PP * size_Mat];
      //   D.f[DIR_0PP] = &DD[DIR_0MM * size_Mat];
      //   D.f[DIR_0MP] = &DD[DIR_0PM * size_Mat];
      //   D.f[DIR_0PM] = &DD[DIR_0MP * size_Mat];
      //   D.f[DIR_000] = &DD[DIR_000 * size_Mat];
      //   D.f[DIR_PPP] = &DD[DIR_MMM * size_Mat];
      //   D.f[DIR_MMP] = &DD[DIR_PPM * size_Mat];
      //   D.f[DIR_PMP] = &DD[DIR_MPM * size_Mat];
      //   D.f[DIR_MPP] = &DD[DIR_PMM * size_Mat];
      //   D.f[DIR_PPM] = &DD[DIR_MMP * size_Mat];
      //   D.f[DIR_MMM] = &DD[DIR_PPP * size_Mat];
      //   D.f[DIR_PMM] = &DD[DIR_MPP * size_Mat];
      //   D.f[DIR_MPM] = &DD[DIR_PMP * size_Mat];
      //}
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();

         (D.f[DIR_P00])[ke   ] = mfabb;//mfcbb;
         (D.f[DIR_M00])[kw   ] = mfcbb;//mfabb;
         (D.f[DIR_0P0])[kn   ] = mfbab;//mfbcb;
         (D.f[DIR_0M0])[ks   ] = mfbcb;//mfbab;
         (D.f[DIR_00P])[kt   ] = mfbba;//mfbbc;
         (D.f[DIR_00M])[kb   ] = mfbbc;//mfbba;
         (D.f[DIR_PP0])[kne  ] = mfaab;//mfccb;
         (D.f[DIR_MM0])[ksw  ] = mfccb;//mfaab;
         (D.f[DIR_PM0])[kse  ] = mfacb;//mfcab;
         (D.f[DIR_MP0])[knw  ] = mfcab;//mfacb;
         (D.f[DIR_P0P])[kte  ] = mfaba;//mfcbc;
         (D.f[DIR_M0M])[kbw  ] = mfcbc;//mfaba;
         (D.f[DIR_P0M])[kbe  ] = mfabc;//mfcba;
         (D.f[DIR_M0P])[ktw  ] = mfcba;//mfabc;
         (D.f[DIR_0PP])[ktn  ] = mfbaa;//mfbcc;
         (D.f[DIR_0MM])[kbs  ] = mfbcc;//mfbaa;
         (D.f[DIR_0PM])[kbn  ] = mfbac;//mfbca;
         (D.f[DIR_0MP])[kts  ] = mfbca;//mfbac;
         (D.f[DIR_000])[kzero] = mfbbb;//mfbbb;
         (D.f[DIR_PPP])[ktne ] = mfaaa;//mfccc;
         (D.f[DIR_MMP])[ktsw ] = mfcca;//mfaac;
         (D.f[DIR_PMP])[ktse ] = mfaca;//mfcac;
         (D.f[DIR_MPP])[ktnw ] = mfcaa;//mfacc;
         (D.f[DIR_PPM])[kbne ] = mfaac;//mfcca;
         (D.f[DIR_MMM])[kbsw ] = mfccc;//mfaaa;
         (D.f[DIR_PMM])[kbse ] = mfacc;//mfcaa;
         (D.f[DIR_MPM])[kbnw ] = mfcac;//mfaca;
         //(D.f[DIR_P00])[ke   ] = mfcbb;
         //(D.f[DIR_M00])[kw   ] = mfabb;
         //(D.f[DIR_0P0])[kn   ] = mfbcb;
         //(D.f[DIR_0M0])[ks   ] = mfbab;
         //(D.f[DIR_00P])[kt   ] = mfbbc;
         //(D.f[DIR_00M])[kb   ] = mfbba;
         //(D.f[DIR_PP0])[kne  ] = mfccb;
         //(D.f[DIR_MM0])[ksw  ] = mfaab;
         //(D.f[DIR_PM0])[kse  ] = mfcab;
         //(D.f[DIR_MP0])[knw  ] = mfacb;
         //(D.f[DIR_P0P])[kte  ] = mfcbc;
         //(D.f[DIR_M0M])[kbw  ] = mfaba;
         //(D.f[DIR_P0M])[kbe  ] = mfcba;
         //(D.f[DIR_M0P])[ktw  ] = mfabc;
         //(D.f[DIR_0PP])[ktn  ] = mfbcc;
         //(D.f[DIR_0MM])[kbs  ] = mfbaa;
         //(D.f[DIR_0PM])[kbn  ] = mfbca;
         //(D.f[DIR_0MP])[kts  ] = mfbac;
         //(D.f[DIR_000])[kzero] = mfbbb;
         //(D.f[DIR_PPP])[ktne ] = mfccc;
         //(D.f[DIR_MMP])[ktsw ] = mfaac;
         //(D.f[DIR_PMP])[ktse ] = mfcac;
         //(D.f[DIR_MPP])[ktnw ] = mfacc;
         //(D.f[DIR_PPM])[kbne ] = mfcca;
         //(D.f[DIR_MMM])[kbsw ] = mfaaa;
         //(D.f[DIR_PMM])[kbse ] = mfcaa;
         //(D.f[DIR_MPM])[kbnw ] = mfaca;

      //(D.f[DIR_P00])[ke   ] = fE ;  //f1_E ;   //fW;    //fE ;
      //(D.f[DIR_M00])[kw   ] = fW ;  //f1_W ;   //fE;    //fW ;
      //(D.f[DIR_0P0])[kn   ] = fN ;  //f1_N ;   //fS;    //fN ;
      //(D.f[DIR_0M0])[ks   ] = fS ;  //f1_S ;   //fN;    //fS ;
      //(D.f[DIR_00P])[kt   ] = fT ;  //f1_T ;   //fB;    //fT ;
      //(D.f[DIR_00M])[kb   ] = fB ;  //f1_B ;   //fT;    //fB ;
      //(D.f[DIR_PP0])[kne  ] = fNE;  //f1_NE;   //fSW;   //fNE;
      //(D.f[DIR_MM0])[ksw  ] = fSW;  //f1_SW;   //fNE;   //fSW;
      //(D.f[DIR_PM0])[kse  ] = fSE;  //f1_SE;   //fNW;   //fSE;
      //(D.f[DIR_MP0])[knw  ] = fNW;  //f1_NW;   //fSE;   //fNW;
      //(D.f[DIR_P0P])[kte  ] = fTE;  //f1_TE;   //fBW;   //fTE;
      //(D.f[DIR_M0M])[kbw  ] = fBW;  //f1_BW;   //fTE;   //fBW;
      //(D.f[DIR_P0M])[kbe  ] = fBE;  //f1_BE;   //fTW;   //fBE;
      //(D.f[DIR_M0P])[ktw  ] = fTW;  //f1_TW;   //fBE;   //fTW;
      //(D.f[DIR_0PP])[ktn  ] = fTN;  //f1_TN;   //fBS;   //fTN;
      //(D.f[DIR_0MM])[kbs  ] = fBS;  //f1_BS;   //fTN;   //fBS;
      //(D.f[DIR_0PM])[kbn  ] = fBN;  //f1_BN;   //fTS;   //fBN;
      //(D.f[DIR_0MP])[kts  ] = fTS;  //f1_TS;   //fBN;   //fTS;
      //(D.f[DIR_000])[kzero] = fZERO;//f1_ZERO; //fZERO; //fZERO;
      //(D.f[DIR_PPP])[ktne ] = fTNE; //f1_TNE;  //fBSW;  //fTNE;
      //(D.f[DIR_MMM])[kbsw ] = fBSW; //f1_BSW;  //fTNE;  //fBSW;
      //(D.f[DIR_PPM])[kbne ] = fBNE; //f1_BNE;  //fTSW;  //fBNE;
      //(D.f[DIR_MMP])[ktsw ] = fTSW; //f1_TSW;  //fBNE;  //fTSW;
      //(D.f[DIR_PMP])[ktse ] = fTSE; //f1_TSE;  //fBNW;  //fTSE;
      //(D.f[DIR_MPM])[kbnw ] = fBNW; //f1_BNW;  //fTSE;  //fBNW;
      //(D.f[DIR_PMM])[kbse ] = fBSE; //f1_BSE;  //fTNW;  //fBSE;
      //(D.f[DIR_MPP])[ktnw ] = fTNW; //f1_TNW;  //fBSE;  //fTNW;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      (D.f[DIR_P00])[ke   ] =c0o1;
      (D.f[DIR_M00])[kw   ] =c0o1;
      (D.f[DIR_0P0])[kn   ] =c0o1;
      (D.f[DIR_0M0])[ks   ] =c0o1;
      (D.f[DIR_00P])[kt   ] =c0o1;
      (D.f[DIR_00M])[kb   ] =c0o1;
      (D.f[DIR_PP0])[kne  ] =c0o1;
      (D.f[DIR_MM0])[ksw  ] =c0o1;
      (D.f[DIR_PM0])[kse  ] =c0o1;
      (D.f[DIR_MP0])[knw  ] =c0o1;
      (D.f[DIR_P0P])[kte  ] =c0o1;
      (D.f[DIR_M0M])[kbw  ] =c0o1;
      (D.f[DIR_P0M])[kbe  ] =c0o1;
      (D.f[DIR_M0P])[ktw  ] =c0o1;
      (D.f[DIR_0PP])[ktn  ] =c0o1;
      (D.f[DIR_0MM])[kbs  ] =c0o1;
      (D.f[DIR_0PM])[kbn  ] =c0o1;
      (D.f[DIR_0MP])[kts  ] =c0o1;
      (D.f[DIR_000])[kzero] =c0o1;
      (D.f[DIR_PPP])[ktne ] =c0o1;
      (D.f[DIR_MMP])[ktsw ] =c0o1;
      (D.f[DIR_PMP])[ktse ] =c0o1;
      (D.f[DIR_MPP])[ktnw ] =c0o1;
      (D.f[DIR_PPM])[kbne ] =c0o1;
      (D.f[DIR_MMM])[kbsw ] =c0o1;
      (D.f[DIR_PMM])[kbse ] =c0o1;
      (D.f[DIR_MPM])[kbnw ] =c0o1;
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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[DIR_P00])[k1e   ];
      f1_E    = (D.f[DIR_M00])[k1w   ];
      f1_S    = (D.f[DIR_0P0])[k1n   ];
      f1_N    = (D.f[DIR_0M0])[k1s   ];
      f1_B    = (D.f[DIR_00P])[k1t   ];
      f1_T    = (D.f[DIR_00M])[k1b   ];
      f1_SW   = (D.f[DIR_PP0])[k1ne  ];
      f1_NE   = (D.f[DIR_MM0])[k1sw  ];
      f1_NW   = (D.f[DIR_PM0])[k1se  ];
      f1_SE   = (D.f[DIR_MP0])[k1nw  ];
      f1_BW   = (D.f[DIR_P0P])[k1te  ];
      f1_TE   = (D.f[DIR_M0M])[k1bw  ];
      f1_TW   = (D.f[DIR_P0M])[k1be  ];
      f1_BE   = (D.f[DIR_M0P])[k1tw  ];
      f1_BS   = (D.f[DIR_0PP])[k1tn  ];
      f1_TN   = (D.f[DIR_0MM])[k1bs  ];
      f1_TS   = (D.f[DIR_0PM])[k1bn  ];
      f1_BN   = (D.f[DIR_0MP])[k1ts  ];
      f1_ZERO = (D.f[DIR_000])[k1zero];
      f1_BSW  = (D.f[DIR_PPP])[k1tne ];
      f1_BNE  = (D.f[DIR_MMP])[k1tsw ];
      f1_BNW  = (D.f[DIR_PMP])[k1tse ];
      f1_BSE  = (D.f[DIR_MPP])[k1tnw ];
      f1_TSW  = (D.f[DIR_PPM])[k1bne ];
      f1_TNE  = (D.f[DIR_MMM])[k1bsw ];
      f1_TNW  = (D.f[DIR_PMM])[k1bse ];
      f1_TSE  = (D.f[DIR_MPM])[k1bnw ];

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

      (D.f[DIR_P00])[ke   ] = c2o27* (rhoBC[k]+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      (D.f[DIR_M00])[kw   ] = c2o27* (rhoBC[k]+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      (D.f[DIR_0P0])[kn   ] = c2o27* (rhoBC[k]+c3o1*(    -vx2    )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      (D.f[DIR_0M0])[ks   ] = c2o27* (rhoBC[k]+c3o1*(     vx2    )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      (D.f[DIR_00P])[kt   ] = c2o27* (rhoBC[k]+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      (D.f[DIR_00M])[kb   ] = c2o27* (rhoBC[k]+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      (D.f[DIR_PP0])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MM0])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PM0])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MP0])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_P0P])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_M0M])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_P0M])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_M0P])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0PP])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0MM])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0PM])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_0MP])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_000])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PPP])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MMP])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PMP])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MPP])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PPM])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MMM])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_PMM])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[DIR_MPM])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //
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
      D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
   }
   else
   {
      D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
      D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
      D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
      D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
      D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
      D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
      D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
      D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
      D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
      D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
      D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
      D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
      D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
      D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
      D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
      D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
      D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
      D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
      D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
      D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
      D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
      D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
      D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
      D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
      D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
      D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
      D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
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
      q_dirE   = &QQ[DIR_P00 * numberOfBCnodes];
      q_dirW   = &QQ[DIR_M00 * numberOfBCnodes];
      q_dirN   = &QQ[DIR_0P0 * numberOfBCnodes];
      q_dirS   = &QQ[DIR_0M0 * numberOfBCnodes];
      q_dirT   = &QQ[DIR_00P * numberOfBCnodes];
      q_dirB   = &QQ[DIR_00M * numberOfBCnodes];
      q_dirNE  = &QQ[DIR_PP0 * numberOfBCnodes];
      q_dirSW  = &QQ[DIR_MM0 * numberOfBCnodes];
      q_dirSE  = &QQ[DIR_PM0 * numberOfBCnodes];
      q_dirNW  = &QQ[DIR_MP0 * numberOfBCnodes];
      q_dirTE  = &QQ[DIR_P0P * numberOfBCnodes];
      q_dirBW  = &QQ[DIR_M0M * numberOfBCnodes];
      q_dirBE  = &QQ[DIR_P0M * numberOfBCnodes];
      q_dirTW  = &QQ[DIR_M0P * numberOfBCnodes];
      q_dirTN  = &QQ[DIR_0PP * numberOfBCnodes];
      q_dirBS  = &QQ[DIR_0MM * numberOfBCnodes];
      q_dirBN  = &QQ[DIR_0PM * numberOfBCnodes];
      q_dirTS  = &QQ[DIR_0MP * numberOfBCnodes];
      q_dirTNE = &QQ[DIR_PPP * numberOfBCnodes];
      q_dirTSW = &QQ[DIR_MMP * numberOfBCnodes];
      q_dirTSE = &QQ[DIR_PMP * numberOfBCnodes];
      q_dirTNW = &QQ[DIR_MPP * numberOfBCnodes];
      q_dirBNE = &QQ[DIR_PPM * numberOfBCnodes];
      q_dirBSW = &QQ[DIR_MMM * numberOfBCnodes];
      q_dirBSE = &QQ[DIR_PMM * numberOfBCnodes];
      q_dirBNW = &QQ[DIR_MPM * numberOfBCnodes];
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

      f_W    = (D.f[DIR_P00])[ke   ];
      f_E    = (D.f[DIR_M00])[kw   ];
      f_S    = (D.f[DIR_0P0])[kn   ];
      f_N    = (D.f[DIR_0M0])[ks   ];
      f_B    = (D.f[DIR_00P])[kt   ];
      f_T    = (D.f[DIR_00M])[kb   ];
      f_SW   = (D.f[DIR_PP0])[kne  ];
      f_NE   = (D.f[DIR_MM0])[ksw  ];
      f_NW   = (D.f[DIR_PM0])[kse  ];
      f_SE   = (D.f[DIR_MP0])[knw  ];
      f_BW   = (D.f[DIR_P0P])[kte  ];
      f_TE   = (D.f[DIR_M0M])[kbw  ];
      f_TW   = (D.f[DIR_P0M])[kbe  ];
      f_BE   = (D.f[DIR_M0P])[ktw  ];
      f_BS   = (D.f[DIR_0PP])[ktn  ];
      f_TN   = (D.f[DIR_0MM])[kbs  ];
      f_TS   = (D.f[DIR_0PM])[kbn  ];
      f_BN   = (D.f[DIR_0MP])[kts  ];
      f_BSW  = (D.f[DIR_PPP])[ktne ];
      f_BNE  = (D.f[DIR_MMP])[ktsw ];
      f_BNW  = (D.f[DIR_PMP])[ktse ];
      f_BSE  = (D.f[DIR_MPP])[ktnw ];
      f_TSW  = (D.f[DIR_PPM])[kbne ];
      f_TNE  = (D.f[DIR_MMM])[kbsw ];
      f_TNW  = (D.f[DIR_PMM])[kbse ];
      f_TSE  = (D.f[DIR_MPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
         f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
         f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[DIR_000])[kzero]);

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
         D.f[DIR_P00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_M00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      }
      else
      {
         D.f[DIR_M00] = &DD[DIR_P00 * numberOfLBnodes];
         D.f[DIR_P00] = &DD[DIR_M00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[DIR_000] = &DD[DIR_000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[DIR_000])[k]=c1o10;
      real rhoDiff = drho - rho[k];
      real VeloX = vx1;
      real VeloY = vx2;
      real VeloZ = vx3;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D.f[DIR_M00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c2o27*(rhoDiff + c6o1*( VeloX     )))/(c1o1+q);
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[DIR_P00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c2o27*(rhoDiff + c6o1*(-VeloX     )))/(c1o1+q);
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[DIR_0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c2o27*(rhoDiff + c6o1*( VeloY     )))/(c1o1+q);
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[DIR_0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c2o27*(rhoDiff + c6o1*(-VeloY     )))/(c1o1+q);
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D.f[DIR_00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c2o27*(rhoDiff + c6o1*( VeloZ     )))/(c1o1+q);
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c2o27* (drho+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[DIR_00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c2o27*(rhoDiff + c6o1*(-VeloZ     )))/(c1o1+q);
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[DIR_MM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c1o54*(rhoDiff + c6o1*(VeloX+VeloY)))/(c1o1+q);
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[DIR_PP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloY)))/(c1o1+q);
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[DIR_MP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloY)))/(c1o1+q);
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[DIR_PM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloY)))/(c1o1+q);
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[DIR_M0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c1o54*(rhoDiff + c6o1*( VeloX+VeloZ)))/(c1o1+q);
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[DIR_P0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloZ)))/(c1o1+q);
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[DIR_M0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloZ)))/(c1o1+q);
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[DIR_P0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloZ)))/(c1o1+q);
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[DIR_0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c1o54*(rhoDiff + c6o1*( VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[DIR_0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c1o54*(rhoDiff + c6o1*( -VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[DIR_0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c1o54*(rhoDiff + c6o1*( VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o54* (drho+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[DIR_0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c1o54*(rhoDiff + c6o1*( -VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[DIR_MMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[DIR_PPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[DIR_MMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[DIR_PPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[DIR_MPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY+VeloZ)))/(c1o1+q);
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[DIR_PMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[DIR_MPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY-VeloZ)))/(c1o1+q);
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         feq=c1o216*(drho+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
         (D.f[DIR_PMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY+VeloZ)))/(c1o1+q);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


