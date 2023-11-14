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
//! \file SlipBCs27.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDevice27(
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
   Distributions27 D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);

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
      real f_W    = (D.f[dP00])[ke   ];
      real f_E    = (D.f[dM00])[kw   ];
      real f_S    = (D.f[d0P0])[kn   ];
      real f_N    = (D.f[d0M0])[ks   ];
      real f_B    = (D.f[d00P])[kt   ];
      real f_T    = (D.f[d00M])[kb   ];
      real f_SW   = (D.f[dPP0])[kne  ];
      real f_NE   = (D.f[dMM0])[ksw  ];
      real f_NW   = (D.f[dPM0])[kse  ];
      real f_SE   = (D.f[dMP0])[knw  ];
      real f_BW   = (D.f[dP0P])[kte  ];
      real f_TE   = (D.f[dM0M])[kbw  ];
      real f_TW   = (D.f[dP0M])[kbe  ];
      real f_BE   = (D.f[dM0P])[ktw  ];
      real f_BS   = (D.f[d0PP])[ktn  ];
      real f_TN   = (D.f[d0MM])[kbs  ];
      real f_TS   = (D.f[d0PM])[kbn  ];
      real f_BN   = (D.f[d0MP])[kts  ];
      real f_BSW  = (D.f[dPPP])[ktne ];
      real f_BNE  = (D.f[dMMP])[ktsw ];
      real f_BNW  = (D.f[dPMP])[ktse ];
      real f_BSE  = (D.f[dMPP])[ktnw ];
      real f_TSW  = (D.f[dPPM])[kbne ];
      real f_TNE  = (D.f[dMMM])[kbsw ];
      real f_TNW  = (D.f[dPMM])[kbse ];
      real f_TSE  = (D.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

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

      D = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes, !isEvenTimestep);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Test
      //(D.f[d000])[k]=c1o10;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real fac = c1o1;//c99o100;
      real VeloX = fac*vx1;
      real VeloY = fac*vx2;
      real VeloZ = fac*vx3;
      bool x = false;
      bool y = false;
      bool z = false;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         x = true;
         feq=c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq); 
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-feq*om1)/(c1o1-om1)+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q);
         //(D.f[dM00])[kw]=zero;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         x = true;
         feq=c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-feq*om1)/(c1o1-om1)+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q);
         //(D.f[dP00])[ke]=zero;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = c0o1;
         VeloZ = fac*vx3;
         y = true;
         feq=c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         (D.f[d0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-feq*om1)/(c1o1-om1)+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q);
         //(D.f[d0M0])[ks]=zero;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = c0o1;
         VeloZ = fac*vx3;
         y = true;
         feq=c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         (D.f[d0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-feq*om1)/(c1o1-om1)+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q);
         //(D.f[d0P0])[kn]=zero;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = c0o1;
         z = true;
         feq=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); 
         (D.f[d00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-feq*om1)/(c1o1-om1)+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q);
         //(D.f[d00M])[kb]=one;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = c0o1;
         z = true;
         feq=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         (D.f[d00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-feq*om1)/(c1o1-om1)+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q);
         //(D.f[d00P])[kt]=zero;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         (D.f[dMM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-feq*om1)/(c1o1-om1)+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q);
         //(D.f[dMM0])[ksw]=zero;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         (D.f[dPP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-feq*om1)/(c1o1-om1)+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q);
         //(D.f[dPP0])[kne]=zero;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         (D.f[dMP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-feq*om1)/(c1o1-om1)+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q);
         //(D.f[dMP0])[knw]=zero;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         (D.f[dPM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-feq*om1)/(c1o1-om1)+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q);
         //(D.f[dPM0])[kse]=zero;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         (D.f[dM0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-feq*om1)/(c1o1-om1)+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q);
         //(D.f[dM0M])[kbw]=zero;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         (D.f[dP0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-feq*om1)/(c1o1-om1)+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q);
         //(D.f[dP0P])[kte]=zero;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         (D.f[dM0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-feq*om1)/(c1o1-om1)+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q);
         //(D.f[dM0P])[ktw]=zero;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         (D.f[dP0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-feq*om1)/(c1o1-om1)+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q);
         //(D.f[dP0M])[kbe]=zero;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         (D.f[d0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-feq*om1)/(c1o1-om1)+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q);
         //(D.f[d0MM])[kbs]=zero;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         (D.f[d0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-feq*om1)/(c1o1-om1)+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q);
         //(D.f[d0PP])[ktn]=zero;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         (D.f[d0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-feq*om1)/(c1o1-om1)+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q);
         //(D.f[d0MP])[kts]=zero;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         (D.f[d0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-feq*om1)/(c1o1-om1)+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q);
         //(D.f[d0PM])[kbn]=zero;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         (D.f[dMMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-feq*om1)/(c1o1-om1)+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[dMMM])[kbsw]=zero;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         (D.f[dPPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-feq*om1)/(c1o1-om1)+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[dPPP])[ktne]=zero;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         (D.f[dMMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-feq*om1)/(c1o1-om1)+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[dMMP])[ktsw]=zero;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         (D.f[dPPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-feq*om1)/(c1o1-om1)+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[dPPM])[kbne]=zero;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         (D.f[dMPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-feq*om1)/(c1o1-om1)+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q);
         //(D.f[dMPM])[kbnw]=zero;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         (D.f[dPMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-feq*om1)/(c1o1-om1)+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q);
         //(D.f[dPMP])[ktse]=zero;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         (D.f[dMPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-feq*om1)/(c1o1-om1)+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q);
         //(D.f[dMPP])[ktnw]=zero;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = fac*vx1;
         VeloY = fac*vx2;
         VeloZ = fac*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         feq=c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         (D.f[dPMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-feq*om1)/(c1o1-om1)+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q);
         //(D.f[dPMM])[kbse]=zero;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDeviceComp27(
    real* distributions, 
    int* subgridDistanceIndices, 
    real* subgridDistances,
    unsigned int numberOfBCnodes,
    real omega, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   //! The slip boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

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
      unsigned int kzero= indexOfBCnode;
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
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[d0M0])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[d0P0])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[d00M])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[d00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[d00P])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, omega, velocityBC, c2o27);
      }

      q = (subgridD.q[dPP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[dMM0])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dMM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[dPP0])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dPM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[dMP0])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dMP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[dPM0])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[dM0M])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloZ;
         (dist.f[dP0P])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[dM0P])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[dP0M])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[d0MM])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[d0PP])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, omega, velocityBC, c1o54);
      }


      q = (subgridD.q[d0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[d0MP])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[d0PM])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, omega, velocityBC, c1o54);
      }

      q = (subgridD.q[dPPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[dMMM])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[dPPP])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, omega, velocityBC, c1o216);
      }


      q = (subgridD.q[dPPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[dMMP])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[dPPM])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[dMPM])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[dPMP])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[dMPP])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, omega, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[dPMM])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, omega, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



























//////////////////////////////////////////////////////////////////////////////
__global__ void BBSlipDeviceComp27(
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
   //! The slip boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

   if(nodeIndex < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist = vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local subgrid distances (q's)
      //!
      SubgridDistances27 subgridD;
      getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);
      
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int indexOfBCnode  = subgridDistanceIndices[nodeIndex];
      unsigned int kzero= indexOfBCnode;
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
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      // real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!

      dist = vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, !isEvenTimestep);
      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real q, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getBounceBackDistributionForVeloBC(f_W, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getBounceBackDistributionForVeloBC(f_E, velocityBC, c2o27);
      }

      q = (subgridD.q[d0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityBC = VeloY;
         (dist.f[d0M0])[ks] = getBounceBackDistributionForVeloBC(f_S, velocityBC, c2o27);
      }

      q = (subgridD.q[d0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityBC = -VeloY;
         (dist.f[d0P0])[kn] = getBounceBackDistributionForVeloBC(f_N, velocityBC, c2o27);
      }

      q = (subgridD.q[d00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityBC = VeloZ;
         (dist.f[d00M])[kb] = getBounceBackDistributionForVeloBC(f_B, velocityBC, c2o27);
      }

      q = (subgridD.q[d00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityBC = -VeloZ;
         (dist.f[d00P])[kt] = getBounceBackDistributionForVeloBC(f_T, velocityBC, c2o27);
      }

      q = (subgridD.q[dPP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityBC = VeloX + VeloY;
         (dist.f[dMM0])[ksw] = getBounceBackDistributionForVeloBC(f_SW, velocityBC, c1o54);
      }

      q = (subgridD.q[dMM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityBC = -VeloX - VeloY;
         (dist.f[dPP0])[kne] = getBounceBackDistributionForVeloBC(f_NE, velocityBC, c1o54);
      }

      q = (subgridD.q[dPM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityBC = VeloX - VeloY;
         (dist.f[dMP0])[knw] = getBounceBackDistributionForVeloBC(f_NW, velocityBC, c1o54);
      }

      q = (subgridD.q[dMP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityBC = -VeloX + VeloY;
         (dist.f[dPM0])[kse] = getBounceBackDistributionForVeloBC(f_SE, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX + VeloZ;
         (dist.f[dM0M])[kbw] = getBounceBackDistributionForVeloBC(f_BW, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

        velocityBC = -VeloX - VeloZ;
        (dist.f[dP0P])[kte] = getBounceBackDistributionForVeloBC(f_TE, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX - VeloZ;
         (dist.f[dM0P])[ktw] = getBounceBackDistributionForVeloBC(f_TW, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloX + VeloZ;
         (dist.f[dP0M])[kbe] = getBounceBackDistributionForVeloBC(f_BE, velocityBC, c1o54);
      }

      q = (subgridD.q[d0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloY + VeloZ;
         (dist.f[d0MM])[kbs] = getBounceBackDistributionForVeloBC(f_BS, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloY - VeloZ;
         (dist.f[d0PP])[ktn] = getBounceBackDistributionForVeloBC(f_TN, velocityBC, c1o54);
      }


      q = (subgridD.q[d0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloY - VeloZ;
         (dist.f[d0MP])[kts] = getBounceBackDistributionForVeloBC(f_TS, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloY + VeloZ;
         (dist.f[d0PM])[kbn] = getBounceBackDistributionForVeloBC(f_BN, velocityBC, c1o54);
      }

      q = (subgridD.q[dPPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[dMMM])[kbsw] = getBounceBackDistributionForVeloBC(f_TNE, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[dPPP])[ktne] = getBounceBackDistributionForVeloBC(f_TNE, velocityBC, c1o216);
      }


      q = (subgridD.q[dPPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[dMMP])[ktsw] = getBounceBackDistributionForVeloBC(f_TSW, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[dPPM])[kbne] = getBounceBackDistributionForVeloBC(f_BNE, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[dMPM])[kbnw] = getBounceBackDistributionForVeloBC(f_BNW, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[dPMP])[ktse] = getBounceBackDistributionForVeloBC(f_TSE, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[dMPP])[ktnw] = getBounceBackDistributionForVeloBC(f_TNW, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[dPMM])[kbse] = getBounceBackDistributionForVeloBC(f_BSE, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




























////////////////////////////////////////////////////////////////////////////
__global__ void QSlipDeviceComp27TurbViscosity(
    real* distributions, 
    int* subgridDistanceIndices, 
    real* subgridDistances,
    unsigned int numberOfBCnodes,
    real omega, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* turbViscosity,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   //! The slip boundary condition is executed in the following steps
   //!

   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

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
      unsigned int kzero= indexOfBCnode;
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
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - compute local relaxation rate
      //!
      real om_turb = omega / (c1o1 + c3o1* omega* max(c0o1, turbViscosity[indexOfBCnode]) );

      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloBC(q, f_E, f_W, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloBC(q, f_W, f_E, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[d0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[d0M0])[ks] = getInterpolatedDistributionForVeloBC(q, f_N, f_S, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[d0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[d0P0])[kn] = getInterpolatedDistributionForVeloBC(q, f_S, f_N, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[d00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[d00M])[kb] = getInterpolatedDistributionForVeloBC(q, f_T, f_B, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[d00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[d00P])[kt] = getInterpolatedDistributionForVeloBC(q, f_B, f_T, feq, om_turb, velocityBC, c2o27);
      }

      q = (subgridD.q[dPP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[dMM0])[ksw] = getInterpolatedDistributionForVeloBC(q, f_NE, f_SW, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dMM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[dPP0])[kne] = getInterpolatedDistributionForVeloBC(q, f_SW, f_NE, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dPM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[dMP0])[knw] = getInterpolatedDistributionForVeloBC(q, f_SE, f_NW, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dMP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[dPM0])[kse] = getInterpolatedDistributionForVeloBC(q, f_NW, f_SE, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[dM0M])[kbw] = getInterpolatedDistributionForVeloBC(q, f_TE, f_BW, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

        velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloZ;
        (dist.f[dP0P])[kte] = getInterpolatedDistributionForVeloBC(q, f_BW, f_TE, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[dM0P])[ktw] = getInterpolatedDistributionForVeloBC(q, f_BE, f_TW, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[dP0M])[kbe] = getInterpolatedDistributionForVeloBC(q, f_TW, f_BE, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[d0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[d0MM])[kbs] = getInterpolatedDistributionForVeloBC(q, f_TN, f_BS, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[d0PP])[ktn] = getInterpolatedDistributionForVeloBC(q, f_BS, f_TN, feq, om_turb, velocityBC, c1o54);
      }


      q = (subgridD.q[d0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[d0MP])[kts] = getInterpolatedDistributionForVeloBC(q, f_BN, f_TS, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[d0PM])[kbn] = getInterpolatedDistributionForVeloBC(q, f_TS, f_BN, feq, om_turb, velocityBC, c1o54);
      }

      q = (subgridD.q[dPPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[dMMM])[kbsw] = getInterpolatedDistributionForVeloBC(q, f_TNE, f_BSW, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[dPPP])[ktne] = getInterpolatedDistributionForVeloBC(q, f_BSW, f_TNE, feq, om_turb, velocityBC, c1o216);
      }


      q = (subgridD.q[dPPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[dMMP])[ktsw] = getInterpolatedDistributionForVeloBC(q, f_BNE, f_TSW, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[dPPM])[kbne] = getInterpolatedDistributionForVeloBC(q, f_TSW, f_BNE, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[dMPM])[kbnw] = getInterpolatedDistributionForVeloBC(q, f_TSE, f_BNW, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[dPMP])[ktse] = getInterpolatedDistributionForVeloBC(q, f_BNW, f_TSE, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[dMPP])[ktnw] = getInterpolatedDistributionForVeloBC(q, f_BSE, f_TNW, feq, om_turb, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[dPMM])[kbse] = getInterpolatedDistributionForVeloBC(q, f_TNW, f_BSE, feq, om_turb, velocityBC, c1o216);
      }
   }
}
////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////
__global__ void QSlipPressureDeviceComp27TurbViscosity(
    real* distributions, 
    int* subgridDistanceIndices, 
    real* subgridDistances,
    unsigned int numberOfBCnodes,
    real omega, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* turbViscosity,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep)
{
   //! The slip boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned nodeIndex = getNodeIndex();

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
      unsigned int kzero= indexOfBCnode;
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
      //! - Calculate macroscopic quantities
      //!
      real drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                  f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                  f_T + f_B + f_N + f_S + f_E + f_W + ((dist.f[d000])[kzero]); 

      real vx1  = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                   (f_E - f_W)) / (c1o1 + drho);

      real vx2  = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                   ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                   (f_N - f_S)) / (c1o1 + drho);

      real vx3  = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                   (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                   (f_T - f_B)) / (c1o1 + drho);

      real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

      ////////////////////////////////////////////////////////////////////////////////
      //! - change the pointer to write the results in the correct array
      //!
      getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - compute local relaxation rate
      //!
      real om_turb = omega / (c1o1 + c3o1* omega* max(c0o1, turbViscosity[indexOfBCnode]) );

      ////////////////////////////////////////////////////////////////////////////////
      //! - Multiply the local velocities by the slipLength
      //!
      real slipLength = c1o1;
      real VeloX = slipLength*vx1;
      real VeloY = slipLength*vx2;
      real VeloZ = slipLength*vx3;

      ////////////////////////////////////////////////////////////////////////////////
      //! - Update distributions with subgrid distance (q) between zero and one
      //!
      real feq, q, velocityLB, velocityBC;

      bool x = false;
      bool y = false;
      bool z = false;

      q = (subgridD.q[dP00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)  // only update distribution for q between zero and one
      {
         VeloX = c0o1;
         x = true;

         velocityLB = vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloX;
         (dist.f[dM00])[kw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_E, f_W, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[dM00])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = c0o1;
         x = true;

         velocityLB = -vx1;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloX;
         (dist.f[dP00])[ke] = getInterpolatedDistributionForVeloWithPressureBC(q, f_W, f_E, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[d0P0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloY;
         (dist.f[d0M0])[ks] = getInterpolatedDistributionForVeloWithPressureBC(q, f_N, f_S, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[d0M0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = c0o1;
         y = true;

         velocityLB = -vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloY;
         (dist.f[d0P0])[kn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_S, f_N, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[d00P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = VeloZ;
         (dist.f[d00M])[kb] = getInterpolatedDistributionForVeloWithPressureBC(q, f_T, f_B, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[d00M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = c0o1;
         z = true;

         velocityLB = -vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
         velocityBC = -VeloZ;
         (dist.f[d00P])[kt] = getInterpolatedDistributionForVeloWithPressureBC(q, f_B, f_T, feq, om_turb, drho, velocityBC, c2o27);
      }

      q = (subgridD.q[dPP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloY;
         (dist.f[dMM0])[ksw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NE, f_SW, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dMM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX - VeloY;
         (dist.f[dPP0])[kne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SW, f_NE, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dPM0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = vx1 - vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloY;
         (dist.f[dMP0])[knw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_SE, f_NW, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dMP0])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;

         velocityLB = -vx1 + vx2;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloY;
         (dist.f[dPM0])[kse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_NW, f_SE, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX + VeloZ;
         (dist.f[dM0M])[kbw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TE, f_BW, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
        VeloX = slipLength*vx1;
        VeloZ = slipLength*vx3;
        if (x == true) VeloX = c0o1;
        if (z == true) VeloZ = c0o1;

        velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloZ;
        (dist.f[dP0P])[kte] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BW, f_TE, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dP0M])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx1 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloX - VeloZ;
         (dist.f[dM0P])[ktw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BE, f_TW, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dM0P])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx1 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloX + VeloZ;
         (dist.f[dP0M])[kbe] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TW, f_BE, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[d0PP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY + VeloZ;
         (dist.f[d0MM])[kbs] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TN, f_BS, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY - VeloZ;
         (dist.f[d0PP])[ktn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BS, f_TN, feq, om_turb, drho, velocityBC, c1o54);
      }


      q = (subgridD.q[d0PM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = VeloY - VeloZ;
         (dist.f[d0MP])[kts] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BN, f_TS, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[d0MP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;

         velocityLB = -vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
         velocityBC = -VeloY + VeloZ;
         (dist.f[d0PM])[kbn] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TS, f_BN, feq, om_turb, drho, velocityBC, c1o54);
      }

      q = (subgridD.q[dPPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY + VeloZ;
         (dist.f[dMMM])[kbsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNE, f_BSW, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY - VeloZ;
         (dist.f[dPPP])[ktne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSW, f_TNE, feq, om_turb, drho, velocityBC, c1o216);
      }


      q = (subgridD.q[dPPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX + VeloY - VeloZ;
         (dist.f[dMMP])[ktsw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNE, f_TSW, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dMMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX - VeloY + VeloZ;
         (dist.f[dPPM])[kbne] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSW, f_BNE, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY + VeloZ;
         (dist.f[dMPM])[kbnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TSE, f_BNW, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY - VeloZ;
         (dist.f[dPMP])[ktse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BNW, f_TSE, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dPMM])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = vx1 - vx2 - vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = VeloX - VeloY - VeloZ;
         (dist.f[dMPP])[ktnw] = getInterpolatedDistributionForVeloWithPressureBC(q, f_BSE, f_TNW, feq, om_turb, drho, velocityBC, c1o216);
      }

      q = (subgridD.q[dMPP])[nodeIndex];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = slipLength*vx1;
         VeloY = slipLength*vx2;
         VeloZ = slipLength*vx3;
         if (x == true) VeloX = c0o1;
         if (y == true) VeloY = c0o1;
         if (z == true) VeloZ = c0o1;
         velocityLB = -vx1 + vx2 + vx3;
         feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
         velocityBC = -VeloX + VeloY + VeloZ;
         (dist.f[dPMM])[kbse] = getInterpolatedDistributionForVeloWithPressureBC(q, f_TNW, f_BSE, feq, om_turb, drho, velocityBC, c1o216);
      }
   }
}

// __global__ void QSlipDeviceComp27TurbViscosity(real* DD, 
//                                              int* k_Q, 
//                                              real* QQ,
//                                              unsigned int numberOfBCnodes,
//                                              real om1, 
//                                              unsigned int* neighborX,
//                                              unsigned int* neighborY,
//                                              unsigned int* neighborZ,
//                                   real* turbViscosity,
//                                              unsigned int size_Mat, 
//                                              bool isEvenTimestep)
// {
//    Distributions27 D;
//    if (isEvenTimestep==true)
//    {
//       D.f[dP00] = &DD[dP00 * size_Mat];
//       D.f[dM00] = &DD[dM00 * size_Mat];
//       D.f[d0P0] = &DD[d0P0 * size_Mat];
//       D.f[d0M0] = &DD[d0M0 * size_Mat];
//       D.f[d00P] = &DD[d00P * size_Mat];
//       D.f[d00M] = &DD[d00M * size_Mat];
//       D.f[dPP0] = &DD[dPP0 * size_Mat];
//       D.f[dMM0] = &DD[dMM0 * size_Mat];
//       D.f[dPM0] = &DD[dPM0 * size_Mat];
//       D.f[dMP0] = &DD[dMP0 * size_Mat];
//       D.f[dP0P] = &DD[dP0P * size_Mat];
//       D.f[dM0M] = &DD[dM0M * size_Mat];
//       D.f[dP0M] = &DD[dP0M * size_Mat];
//       D.f[dM0P] = &DD[dM0P * size_Mat];
//       D.f[d0PP] = &DD[d0PP * size_Mat];
//       D.f[d0MM] = &DD[d0MM * size_Mat];
//       D.f[d0PM] = &DD[d0PM * size_Mat];
//       D.f[d0MP] = &DD[d0MP * size_Mat];
//       D.f[d000] = &DD[d000 * size_Mat];
//       D.f[dPPP] = &DD[dPPP * size_Mat];
//       D.f[dMMP] = &DD[dMMP * size_Mat];
//       D.f[dPMP] = &DD[dPMP * size_Mat];
//       D.f[dMPP] = &DD[dMPP * size_Mat];
//       D.f[dPPM] = &DD[dPPM * size_Mat];
//       D.f[dMMM] = &DD[dMMM * size_Mat];
//       D.f[dPMM] = &DD[dPMM * size_Mat];
//       D.f[dMPM] = &DD[dMPM * size_Mat];
//    } 
//    else
//    {
//       D.f[dM00] = &DD[dP00 * size_Mat];
//       D.f[dP00] = &DD[dM00 * size_Mat];
//       D.f[d0M0] = &DD[d0P0 * size_Mat];
//       D.f[d0P0] = &DD[d0M0 * size_Mat];
//       D.f[d00M] = &DD[d00P * size_Mat];
//       D.f[d00P] = &DD[d00M * size_Mat];
//       D.f[dMM0] = &DD[dPP0 * size_Mat];
//       D.f[dPP0] = &DD[dMM0 * size_Mat];
//       D.f[dMP0] = &DD[dPM0 * size_Mat];
//       D.f[dPM0] = &DD[dMP0 * size_Mat];
//       D.f[dM0M] = &DD[dP0P * size_Mat];
//       D.f[dP0P] = &DD[dM0M * size_Mat];
//       D.f[dM0P] = &DD[dP0M * size_Mat];
//       D.f[dP0M] = &DD[dM0P * size_Mat];
//       D.f[d0MM] = &DD[d0PP * size_Mat];
//       D.f[d0PP] = &DD[d0MM * size_Mat];
//       D.f[d0MP] = &DD[d0PM * size_Mat];
//       D.f[d0PM] = &DD[d0MP * size_Mat];
//       D.f[d000] = &DD[d000 * size_Mat];
//       D.f[dPPP] = &DD[dMMM * size_Mat];
//       D.f[dMMP] = &DD[dPPM * size_Mat];
//       D.f[dPMP] = &DD[dMPM * size_Mat];
//       D.f[dMPP] = &DD[dPMM * size_Mat];
//       D.f[dPPM] = &DD[dMMP * size_Mat];
//       D.f[dMMM] = &DD[dPPP * size_Mat];
//       D.f[dPMM] = &DD[dMPP * size_Mat];
//       D.f[dMPM] = &DD[dPMP * size_Mat];
//    }
//    ////////////////////////////////////////////////////////////////////////////////
//    const unsigned  x = threadIdx.x;  // Globaler x-Index 
//    const unsigned  y = blockIdx.x;   // Globaler y-Index 
//    const unsigned  z = blockIdx.y;   // Globaler z-Index 

//    const unsigned nx = blockDim.x;
//    const unsigned ny = gridDim.x;

//    const unsigned k = nx*(ny*z + y) + x;
//    //////////////////////////////////////////////////////////////////////////

//    if(k<numberOfBCnodes)
//    {
//       ////////////////////////////////////////////////////////////////////////////////
//       real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
//             *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
//             *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
//             *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
//             *q_dirBSE, *q_dirBNW; 
//       q_dirE   = &QQ[dP00 * numberOfBCnodes];
//       q_dirW   = &QQ[dM00 * numberOfBCnodes];
//       q_dirN   = &QQ[d0P0 * numberOfBCnodes];
//       q_dirS   = &QQ[d0M0 * numberOfBCnodes];
//       q_dirT   = &QQ[d00P * numberOfBCnodes];
//       q_dirB   = &QQ[d00M * numberOfBCnodes];
//       q_dirNE  = &QQ[dPP0 * numberOfBCnodes];
//       q_dirSW  = &QQ[dMM0 * numberOfBCnodes];
//       q_dirSE  = &QQ[dPM0 * numberOfBCnodes];
//       q_dirNW  = &QQ[dMP0 * numberOfBCnodes];
//       q_dirTE  = &QQ[dP0P * numberOfBCnodes];
//       q_dirBW  = &QQ[dM0M * numberOfBCnodes];
//       q_dirBE  = &QQ[dP0M * numberOfBCnodes];
//       q_dirTW  = &QQ[dM0P * numberOfBCnodes];
//       q_dirTN  = &QQ[d0PP * numberOfBCnodes];
//       q_dirBS  = &QQ[d0MM * numberOfBCnodes];
//       q_dirBN  = &QQ[d0PM * numberOfBCnodes];
//       q_dirTS  = &QQ[d0MP * numberOfBCnodes];
//       q_dirTNE = &QQ[dPPP * numberOfBCnodes];
//       q_dirTSW = &QQ[dMMP * numberOfBCnodes];
//       q_dirTSE = &QQ[dPMP * numberOfBCnodes];
//       q_dirTNW = &QQ[dMPP * numberOfBCnodes];
//       q_dirBNE = &QQ[dPPM * numberOfBCnodes];
//       q_dirBSW = &QQ[dMMM * numberOfBCnodes];
//       q_dirBSE = &QQ[dPMM * numberOfBCnodes];
//       q_dirBNW = &QQ[dMPM * numberOfBCnodes];
//       ////////////////////////////////////////////////////////////////////////////////
//       //index
//       unsigned int KQK  = k_Q[k];
//       unsigned int kzero= KQK;
//       unsigned int ke   = KQK;
//       unsigned int kw   = neighborX[KQK];
//       unsigned int kn   = KQK;
//       unsigned int ks   = neighborY[KQK];
//       unsigned int kt   = KQK;
//       unsigned int kb   = neighborZ[KQK];
//       unsigned int ksw  = neighborY[kw];
//       unsigned int kne  = KQK;
//       unsigned int kse  = ks;
//       unsigned int knw  = kw;
//       unsigned int kbw  = neighborZ[kw];
//       unsigned int kte  = KQK;
//       unsigned int kbe  = kb;
//       unsigned int ktw  = kw;
//       unsigned int kbs  = neighborZ[ks];
//       unsigned int ktn  = KQK;
//       unsigned int kbn  = kb;
//       unsigned int kts  = ks;
//       unsigned int ktse = ks;
//       unsigned int kbnw = kbw;
//       unsigned int ktnw = kw;
//       unsigned int kbse = kbs;
//       unsigned int ktsw = ksw;
//       unsigned int kbne = kb;
//       unsigned int ktne = KQK;
//       unsigned int kbsw = neighborZ[ksw];
      
//       ////////////////////////////////////////////////////////////////////////////////
//       real f_W    = (D.f[dP00])[ke   ];
//       real f_E    = (D.f[dM00])[kw   ];
//       real f_S    = (D.f[d0P0])[kn   ];
//       real f_N    = (D.f[d0M0])[ks   ];
//       real f_B    = (D.f[d00P])[kt   ];
//       real f_T    = (D.f[d00M])[kb   ];
//       real f_SW   = (D.f[dPP0])[kne  ];
//       real f_NE   = (D.f[dMM0])[ksw  ];
//       real f_NW   = (D.f[dPM0])[kse  ];
//       real f_SE   = (D.f[dMP0])[knw  ];
//       real f_BW   = (D.f[dP0P])[kte  ];
//       real f_TE   = (D.f[dM0M])[kbw  ];
//       real f_TW   = (D.f[dP0M])[kbe  ];
//       real f_BE   = (D.f[dM0P])[ktw  ];
//       real f_BS   = (D.f[d0PP])[ktn  ];
//       real f_TN   = (D.f[d0MM])[kbs  ];
//       real f_TS   = (D.f[d0PM])[kbn  ];
//       real f_BN   = (D.f[d0MP])[kts  ];
//       real f_BSW  = (D.f[dPPP])[ktne ];
//       real f_BNE  = (D.f[dMMP])[ktsw ];
//       real f_BNW  = (D.f[dPMP])[ktse ];
//       real f_BSE  = (D.f[dMPP])[ktnw ];
//       real f_TSW  = (D.f[dPPM])[kbne ];
//       real f_TNE  = (D.f[dMMM])[kbsw ];
//       real f_TNW  = (D.f[dPMM])[kbse ];
//       real f_TSE  = (D.f[dMPM])[kbnw ];
//       ////////////////////////////////////////////////////////////////////////////////
//       real vx1, vx2, vx3, drho, feq, q;
//       drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
//                 f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
//                 f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

//       vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
//                 ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
//                 (f_E - f_W)) / (c1o1 + drho); 
         

//       vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
//                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
//                  (f_N - f_S)) / (c1o1 + drho); 

//       vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
//                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
//                  (f_T - f_B)) / (c1o1 + drho); 

//       real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

//       //////////////////////////////////////////////////////////////////////////
//       if (isEvenTimestep==false)
//       {
//          D.f[dP00] = &DD[dP00 * size_Mat];
//          D.f[dM00] = &DD[dM00 * size_Mat];
//          D.f[d0P0] = &DD[d0P0 * size_Mat];
//          D.f[d0M0] = &DD[d0M0 * size_Mat];
//          D.f[d00P] = &DD[d00P * size_Mat];
//          D.f[d00M] = &DD[d00M * size_Mat];
//          D.f[dPP0] = &DD[dPP0 * size_Mat];
//          D.f[dMM0] = &DD[dMM0 * size_Mat];
//          D.f[dPM0] = &DD[dPM0 * size_Mat];
//          D.f[dMP0] = &DD[dMP0 * size_Mat];
//          D.f[dP0P] = &DD[dP0P * size_Mat];
//          D.f[dM0M] = &DD[dM0M * size_Mat];
//          D.f[dP0M] = &DD[dP0M * size_Mat];
//          D.f[dM0P] = &DD[dM0P * size_Mat];
//          D.f[d0PP] = &DD[d0PP * size_Mat];
//          D.f[d0MM] = &DD[d0MM * size_Mat];
//          D.f[d0PM] = &DD[d0PM * size_Mat];
//          D.f[d0MP] = &DD[d0MP * size_Mat];
//          D.f[d000] = &DD[d000 * size_Mat];
//          D.f[dPPP] = &DD[dPPP * size_Mat];
//          D.f[dMMP] = &DD[dMMP * size_Mat];
//          D.f[dPMP] = &DD[dPMP * size_Mat];
//          D.f[dMPP] = &DD[dMPP * size_Mat];
//          D.f[dPPM] = &DD[dPPM * size_Mat];
//          D.f[dMMM] = &DD[dMMM * size_Mat];
//          D.f[dPMM] = &DD[dPMM * size_Mat];
//          D.f[dMPM] = &DD[dMPM * size_Mat];
//       } 
//       else
//       {
//          D.f[dM00] = &DD[dP00 * size_Mat];
//          D.f[dP00] = &DD[dM00 * size_Mat];
//          D.f[d0M0] = &DD[d0P0 * size_Mat];
//          D.f[d0P0] = &DD[d0M0 * size_Mat];
//          D.f[d00M] = &DD[d00P * size_Mat];
//          D.f[d00P] = &DD[d00M * size_Mat];
//          D.f[dMM0] = &DD[dPP0 * size_Mat];
//          D.f[dPP0] = &DD[dMM0 * size_Mat];
//          D.f[dMP0] = &DD[dPM0 * size_Mat];
//          D.f[dPM0] = &DD[dMP0 * size_Mat];
//          D.f[dM0M] = &DD[dP0P * size_Mat];
//          D.f[dP0P] = &DD[dM0M * size_Mat];
//          D.f[dM0P] = &DD[dP0M * size_Mat];
//          D.f[dP0M] = &DD[dM0P * size_Mat];
//          D.f[d0MM] = &DD[d0PP * size_Mat];
//          D.f[d0PP] = &DD[d0MM * size_Mat];
//          D.f[d0MP] = &DD[d0PM * size_Mat];
//          D.f[d0PM] = &DD[d0MP * size_Mat];
//          D.f[d000] = &DD[d000 * size_Mat];
//          D.f[dPPP] = &DD[dMMM * size_Mat];
//          D.f[dMMP] = &DD[dPPM * size_Mat];
//          D.f[dPMP] = &DD[dMPM * size_Mat];
//          D.f[dMPP] = &DD[dPMM * size_Mat];
//          D.f[dPPM] = &DD[dMMP * size_Mat];
//          D.f[dMMM] = &DD[dPPP * size_Mat];
//          D.f[dPMM] = &DD[dMPP * size_Mat];
//          D.f[dMPM] = &DD[dPMP * size_Mat];
//       }
//       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       //Test
//       //(D.f[d000])[k]=c1o10;
//       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       real om_turb = om1 / (c1o1 + c3o1*om1*max(c0o1, turbViscosity[k_Q[k]]));
     
//      real fac = c1o1;//c99o100;
//       real VeloX = fac*vx1;
//       real VeloY = fac*vx2;
//       real VeloZ = fac*vx3;
//       bool x = false;
//       bool y = false;
//       bool z = false;

//       q = q_dirE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = c0o1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          x = true;
//          feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
//          //(D.f[dM00])[kw]=(one-q)/(one+q)*(f_E-feq*om1)/(one-om1)+(q*(f_E+f_W)-six*c2over27*( VeloX     ))/(one+q);
//          //(D.f[dM00])[kw]=zero;
//       }

//       q = q_dirW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = c0o1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          x = true;
//          feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
//          //(D.f[dP00])[ke]=(one-q)/(one+q)*(f_W-feq*om_turb)/(one-om_turb)+(q*(f_W+f_E)-six*c2over27*(-VeloX     ))/(one+q);
//          //(D.f[dP00])[ke]=zero;
//       }

//       q = q_dirN[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = c0o1;
//          VeloZ = fac*vx3;
//          y = true;
//          feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
//          //(D.f[d0M0])[ks]=(one-q)/(one+q)*(f_N-feq*om_turb)/(one-om_turb)+(q*(f_N+f_S)-six*c2over27*( VeloY     ))/(one+q);
//          //(D.f[d0M0])[ks]=zero;
//       }

//       q = q_dirS[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = c0o1;
//          VeloZ = fac*vx3;
//          y = true;
//          feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
//          //(D.f[d0P0])[kn]=(one-q)/(one+q)*(f_S-feq*om_turb)/(one-om_turb)+(q*(f_S+f_N)-six*c2over27*(-VeloY     ))/(one+q);
//          //(D.f[d0P0])[kn]=zero;
//       }

//       q = q_dirT[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = c0o1;
//          z = true;
//          feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
//          //(D.f[d00M])[kb]=(one-q)/(one+q)*(f_T-feq*om_turb)/(one-om_turb)+(q*(f_T+f_B)-six*c2over27*( VeloZ     ))/(one+q);
//          //(D.f[d00M])[kb]=one;
//       }

//       q = q_dirB[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = c0o1;
//          z = true;
//          feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
//          //feq=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
//          //(D.f[d00P])[kt]=(one-q)/(one+q)*(f_B-feq*om_turb)/(one-om_turb)+(q*(f_B+f_T)-six*c2over27*(-VeloZ     ))/(one+q);
//          //(D.f[d00P])[kt]=zero;
//       }

//       q = q_dirNE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
//          //(D.f[dMM0])[ksw]=(one-q)/(one+q)*(f_NE-feq*om_turb)/(one-om_turb)+(q*(f_NE+f_SW)-six*c1over54*(VeloX+VeloY))/(one+q);
//          //(D.f[dMM0])[ksw]=zero;
//       }

//       q = q_dirSW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
//          //(D.f[dPP0])[kne]=(one-q)/(one+q)*(f_SW-feq*om_turb)/(one-om_turb)+(q*(f_SW+f_NE)-six*c1over54*(-VeloX-VeloY))/(one+q);
//          //(D.f[dPP0])[kne]=zero;
//       }

//       q = q_dirSE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
//          //(D.f[dMP0])[knw]=(one-q)/(one+q)*(f_SE-feq*om_turb)/(one-om_turb)+(q*(f_SE+f_NW)-six*c1over54*( VeloX-VeloY))/(one+q);
//          //(D.f[dMP0])[knw]=zero;
//       }

//       q = q_dirNW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
//          //(D.f[dPM0])[kse]=(one-q)/(one+q)*(f_NW-feq*om_turb)/(one-om_turb)+(q*(f_NW+f_SE)-six*c1over54*(-VeloX+VeloY))/(one+q);
//          //(D.f[dPM0])[kse]=zero;
//       }

//       q = q_dirTE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (z == true) VeloZ = c0o1;
//       //  if (k==10000) printf("AFTER x: %u \t  y: %u \t z: %u \n  VeloX: %f \t VeloY: %f \t VeloZ: %f \n\n", x,y,z, VeloX,VeloY,VeloZ);
//          feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dM0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
//          //(D.f[dM0M])[kbw]=(one-q)/(one+q)*(f_TE-feq*om_turb)/(one-om_turb)+(q*(f_TE+f_BW)-six*c1over54*( VeloX+VeloZ))/(one+q);
//          //(D.f[dM0M])[kbw]=zero;
//       }

//       q = q_dirBW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dP0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
//          //(D.f[dP0P])[kte]=(one-q)/(one+q)*(f_BW-feq*om_turb)/(one-om_turb)+(q*(f_BW+f_TE)-six*c1over54*(-VeloX-VeloZ))/(one+q);
//          //(D.f[dP0P])[kte]=zero;
//       }

//       q = q_dirBE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dM0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
//          //(D.f[dM0P])[ktw]=(one-q)/(one+q)*(f_BE-feq*om_turb)/(one-om_turb)+(q*(f_BE+f_TW)-six*c1over54*( VeloX-VeloZ))/(one+q);
//          //(D.f[dM0P])[ktw]=zero;
//       }

//       q = q_dirTW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dP0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
//          //(D.f[dP0M])[kbe]=(one-q)/(one+q)*(f_TW-feq*om_turb)/(one-om_turb)+(q*(f_TW+f_BE)-six*c1over54*(-VeloX+VeloZ))/(one+q);
//          //(D.f[dP0M])[kbe]=zero;
//       }

//       q = q_dirTN[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
//          //(D.f[d0MM])[kbs]=(one-q)/(one+q)*(f_TN-feq*om_turb)/(one-om_turb)+(q*(f_TN+f_BS)-six*c1over54*( VeloY+VeloZ))/(one+q);
//          //(D.f[d0MM])[kbs]=zero;
//       }

//       q = q_dirBS[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
//          //(D.f[d0PP])[ktn]=(one-q)/(one+q)*(f_BS-feq*om_turb)/(one-om_turb)+(q*(f_BS+f_TN)-six*c1over54*( -VeloY-VeloZ))/(one+q);
//          //(D.f[d0PP])[ktn]=zero;
//       }

//       q = q_dirBN[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
//          //(D.f[d0MP])[kts]=(one-q)/(one+q)*(f_BN-feq*om_turb)/(one-om_turb)+(q*(f_BN+f_TS)-six*c1over54*( VeloY-VeloZ))/(one+q);
//          //(D.f[d0MP])[kts]=zero;
//       }

//       q = q_dirTS[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[d0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
//          //feq=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
//          //(D.f[d0PM])[kbn]=(one-q)/(one+q)*(f_TS-feq*om_turb)/(one-om_turb)+(q*(f_TS+f_BN)-six*c1over54*( -VeloY+VeloZ))/(one+q);
//          //(D.f[d0PM])[kbn]=zero;
//       }

//       q = q_dirTNE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
//          //(D.f[dMMM])[kbsw]=(one-q)/(one+q)*(f_TNE-feq*om_turb)/(one-om_turb)+(q*(f_TNE+f_BSW)-six*c1over216*( VeloX+VeloY+VeloZ))/(one+q);
//          //(D.f[dMMM])[kbsw]=zero;
//       }

//       q = q_dirBSW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
//          //(D.f[dPPP])[ktne]=(one-q)/(one+q)*(f_BSW-feq*om_turb)/(one-om_turb)+(q*(f_BSW+f_TNE)-six*c1over216*(-VeloX-VeloY-VeloZ))/(one+q);
//          //(D.f[dPPP])[ktne]=zero;
//       }

//       q = q_dirBNE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
//          //(D.f[dMMP])[ktsw]=(one-q)/(one+q)*(f_BNE-feq*om_turb)/(one-om_turb)+(q*(f_BNE+f_TSW)-six*c1over216*( VeloX+VeloY-VeloZ))/(one+q);
//          //(D.f[dMMP])[ktsw]=zero;
//       }

//       q = q_dirTSW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
//          //(D.f[dPPM])[kbne]=(one-q)/(one+q)*(f_TSW-feq*om_turb)/(one-om_turb)+(q*(f_TSW+f_BNE)-six*c1over216*(-VeloX-VeloY+VeloZ))/(one+q);
//          //(D.f[dPPM])[kbne]=zero;
//       }

//       q = q_dirTSE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
//          //(D.f[dMPM])[kbnw]=(one-q)/(one+q)*(f_TSE-feq*om_turb)/(one-om_turb)+(q*(f_TSE+f_BNW)-six*c1over216*( VeloX-VeloY+VeloZ))/(one+q);
//          //(D.f[dMPM])[kbnw]=zero;
//       }

//       q = q_dirBNW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
//          //(D.f[dPMP])[ktse]=(one-q)/(one+q)*(f_BNW-feq*om_turb)/(one-om_turb)+(q*(f_BNW+f_TSE)-six*c1over216*(-VeloX+VeloY-VeloZ))/(one+q);
//          //(D.f[dPMP])[ktse]=zero;
//       }

//       q = q_dirBSE[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dMPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
//          //(D.f[dMPP])[ktnw]=(one-q)/(one+q)*(f_BSE-feq*om_turb)/(one-om_turb)+(q*(f_BSE+f_TNW)-six*c1over216*( VeloX-VeloY-VeloZ))/(one+q);
//          //(D.f[dMPP])[ktnw]=zero;
//       }

//       q = q_dirTNW[k];
//       if (q>=c0o1 && q<=c1o1)
//       {
//          VeloX = fac*vx1;
//          VeloY = fac*vx2;
//          VeloZ = fac*vx3;
//          if (x == true) VeloX = c0o1;
//          if (y == true) VeloY = c0o1;
//          if (z == true) VeloZ = c0o1;
//          feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
//          (D.f[dPMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om_turb)/(c1o1-om_turb))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
//          //feq=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
//          //(D.f[dPMM])[kbse]=(one-q)/(one+q)*(f_TNW-feq*om_turb)/(one-om_turb)+(q*(f_TNW+f_BSE)-six*c1over216*(-VeloX+VeloY+VeloZ))/(one+q);
//          //(D.f[dPMM])[kbse]=zero;
//       }
//    }
// }






































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipGeomDeviceComp27(
    real* DD, 
    int* k_Q, 
    real* QQ,
    unsigned int  numberOfBCnodes,
    real om1, 
    real* NormalX,
    real* NormalY,
    real* NormalZ,
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

   if(k< numberOfBCnodes)
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
      real *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[dP00 * numberOfBCnodes];
      nx_dirW   = &NormalX[dM00 * numberOfBCnodes];
      nx_dirN   = &NormalX[d0P0 * numberOfBCnodes];
      nx_dirS   = &NormalX[d0M0 * numberOfBCnodes];
      nx_dirT   = &NormalX[d00P * numberOfBCnodes];
      nx_dirB   = &NormalX[d00M * numberOfBCnodes];
      nx_dirNE  = &NormalX[dPP0 * numberOfBCnodes];
      nx_dirSW  = &NormalX[dMM0 * numberOfBCnodes];
      nx_dirSE  = &NormalX[dPM0 * numberOfBCnodes];
      nx_dirNW  = &NormalX[dMP0 * numberOfBCnodes];
      nx_dirTE  = &NormalX[dP0P * numberOfBCnodes];
      nx_dirBW  = &NormalX[dM0M * numberOfBCnodes];
      nx_dirBE  = &NormalX[dP0M * numberOfBCnodes];
      nx_dirTW  = &NormalX[dM0P * numberOfBCnodes];
      nx_dirTN  = &NormalX[d0PP * numberOfBCnodes];
      nx_dirBS  = &NormalX[d0MM * numberOfBCnodes];
      nx_dirBN  = &NormalX[d0PM * numberOfBCnodes];
      nx_dirTS  = &NormalX[d0MP * numberOfBCnodes];
      nx_dirTNE = &NormalX[dPPP * numberOfBCnodes];
      nx_dirTSW = &NormalX[dMMP * numberOfBCnodes];
      nx_dirTSE = &NormalX[dPMP * numberOfBCnodes];
      nx_dirTNW = &NormalX[dMPP * numberOfBCnodes];
      nx_dirBNE = &NormalX[dPPM * numberOfBCnodes];
      nx_dirBSW = &NormalX[dMMM * numberOfBCnodes];
      nx_dirBSE = &NormalX[dPMM * numberOfBCnodes];
      nx_dirBNW = &NormalX[dMPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[dP00 * numberOfBCnodes];
      ny_dirW   = &NormalY[dM00 * numberOfBCnodes];
      ny_dirN   = &NormalY[d0P0 * numberOfBCnodes];
      ny_dirS   = &NormalY[d0M0 * numberOfBCnodes];
      ny_dirT   = &NormalY[d00P * numberOfBCnodes];
      ny_dirB   = &NormalY[d00M * numberOfBCnodes];
      ny_dirNE  = &NormalY[dPP0 * numberOfBCnodes];
      ny_dirSW  = &NormalY[dMM0 * numberOfBCnodes];
      ny_dirSE  = &NormalY[dPM0 * numberOfBCnodes];
      ny_dirNW  = &NormalY[dMP0 * numberOfBCnodes];
      ny_dirTE  = &NormalY[dP0P * numberOfBCnodes];
      ny_dirBW  = &NormalY[dM0M * numberOfBCnodes];
      ny_dirBE  = &NormalY[dP0M * numberOfBCnodes];
      ny_dirTW  = &NormalY[dM0P * numberOfBCnodes];
      ny_dirTN  = &NormalY[d0PP * numberOfBCnodes];
      ny_dirBS  = &NormalY[d0MM * numberOfBCnodes];
      ny_dirBN  = &NormalY[d0PM * numberOfBCnodes];
      ny_dirTS  = &NormalY[d0MP * numberOfBCnodes];
      ny_dirTNE = &NormalY[dPPP * numberOfBCnodes];
      ny_dirTSW = &NormalY[dMMP * numberOfBCnodes];
      ny_dirTSE = &NormalY[dPMP * numberOfBCnodes];
      ny_dirTNW = &NormalY[dMPP * numberOfBCnodes];
      ny_dirBNE = &NormalY[dPPM * numberOfBCnodes];
      ny_dirBSW = &NormalY[dMMM * numberOfBCnodes];
      ny_dirBSE = &NormalY[dPMM * numberOfBCnodes];
      ny_dirBNW = &NormalY[dMPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[dP00 * numberOfBCnodes];
      nz_dirW   = &NormalZ[dM00 * numberOfBCnodes];
      nz_dirN   = &NormalZ[d0P0 * numberOfBCnodes];
      nz_dirS   = &NormalZ[d0M0 * numberOfBCnodes];
      nz_dirT   = &NormalZ[d00P * numberOfBCnodes];
      nz_dirB   = &NormalZ[d00M * numberOfBCnodes];
      nz_dirNE  = &NormalZ[dPP0 * numberOfBCnodes];
      nz_dirSW  = &NormalZ[dMM0 * numberOfBCnodes];
      nz_dirSE  = &NormalZ[dPM0 * numberOfBCnodes];
      nz_dirNW  = &NormalZ[dMP0 * numberOfBCnodes];
      nz_dirTE  = &NormalZ[dP0P * numberOfBCnodes];
      nz_dirBW  = &NormalZ[dM0M * numberOfBCnodes];
      nz_dirBE  = &NormalZ[dP0M * numberOfBCnodes];
      nz_dirTW  = &NormalZ[dM0P * numberOfBCnodes];
      nz_dirTN  = &NormalZ[d0PP * numberOfBCnodes];
      nz_dirBS  = &NormalZ[d0MM * numberOfBCnodes];
      nz_dirBN  = &NormalZ[d0PM * numberOfBCnodes];
      nz_dirTS  = &NormalZ[d0MP * numberOfBCnodes];
      nz_dirTNE = &NormalZ[dPPP * numberOfBCnodes];
      nz_dirTSW = &NormalZ[dMMP * numberOfBCnodes];
      nz_dirTSE = &NormalZ[dPMP * numberOfBCnodes];
      nz_dirTNW = &NormalZ[dMPP * numberOfBCnodes];
      nz_dirBNE = &NormalZ[dPPM * numberOfBCnodes];
      nz_dirBSW = &NormalZ[dMMM * numberOfBCnodes];
      nz_dirBSE = &NormalZ[dPMM * numberOfBCnodes];
      nz_dirBNW = &NormalZ[dMPM * numberOfBCnodes];
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
      real f_W    = (D.f[dP00])[ke   ];
      real f_E    = (D.f[dM00])[kw   ];
      real f_S    = (D.f[d0P0])[kn   ];
      real f_N    = (D.f[d0M0])[ks   ];
      real f_B    = (D.f[d00P])[kt   ];
      real f_T    = (D.f[d00M])[kb   ];
      real f_SW   = (D.f[dPP0])[kne  ];
      real f_NE   = (D.f[dMM0])[ksw  ];
      real f_NW   = (D.f[dPM0])[kse  ];
      real f_SE   = (D.f[dMP0])[knw  ];
      real f_BW   = (D.f[dP0P])[kte  ];
      real f_TE   = (D.f[dM0M])[kbw  ];
      real f_TW   = (D.f[dP0M])[kbe  ];
      real f_BE   = (D.f[dM0P])[ktw  ];
      real f_BS   = (D.f[d0PP])[ktn  ];
      real f_TN   = (D.f[d0MM])[kbs  ];
      real f_TS   = (D.f[d0PM])[kbn  ];
      real f_BN   = (D.f[d0MP])[kts  ];
      real f_BSW  = (D.f[dPPP])[ktne ];
      real f_BNE  = (D.f[dMMP])[ktsw ];
      real f_BNW  = (D.f[dPMP])[ktse ];
      real f_BSE  = (D.f[dMPP])[ktnw ];
      real f_TSW  = (D.f[dPPM])[kbne ];
      real f_TNE  = (D.f[dMMM])[kbsw ];
      real f_TNW  = (D.f[dPMM])[kbse ];
      real f_TSE  = (D.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

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
      real VeloX = vx1;
      real VeloY = vx2;
      real VeloZ = vx3;
      real fac = c0o1;//0.5;
       real phi = c0o1;
      //real alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real kxyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (c1o1 + drho) - ((vx1*vx2)));
      real kyzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (c1o1 + drho) - ((vx2*vx3)));
      real kxzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (c1o1 + drho) - ((vx1*vx3)));

      real kxxFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx1*vx1)));
      real kyyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx2*vx2)));
      real kzzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx3*vx3)));

      real magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

      //fac = fac * magS / (c1o3 * (one / om1 - c1o2));
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real *facAst = &QQ[d000 * numberOfBCnodes];

      //fac = fac * alpha + facAst[k] * (one - alpha);
      //facAst[k] = fac;
      //(&QQ[d000 * numberOfBCnodes])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////real uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real phi = expf(magS/0.01f) - one;
      //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real C = five;
      //real kappa = 0.41f;
      //real phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
      //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real sum = zero, count = zero;
   //   q = q_dirE   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirE[k] ); count += one;
   //   q = q_dirW   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirW[k])); count += one;
   //   q = q_dirN   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirN[k] ); count += one;
   //   q = q_dirS   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirS[k])); count += one;
   //   q = q_dirT   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirT[k] ); count += one;
   //   q = q_dirB   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirB[k])); count += one;
   //   q = q_dirNE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirNE[k]  + ny_dirNE[k])/(sqrtf(two))); count += one;
   //   q = q_dirSW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirSW[k]) - ny_dirSW[k])/(sqrtf(two))); count += one;
   //   q = q_dirSE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirSE[k]  - ny_dirSE[k])/(sqrtf(two))); count += one;
   //   q = q_dirNW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirNW[k]) + ny_dirNW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTE[k]  + nz_dirTE[k])/(sqrtf(two))); count += one;
   //   q = q_dirBW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBW[k]) - nz_dirBW[k])/(sqrtf(two))); count += one;
   //   q = q_dirBE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBE[k]  - nz_dirBE[k])/(sqrtf(two))); count += one;
   //   q = q_dirTW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTW[k]) + nz_dirTW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirTN[k]  + nz_dirTN[k])/(sqrtf(two))); count += one;
   //   q = q_dirBS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirBS[k]) - nz_dirBS[k])/(sqrtf(two))); count += one;
   //   q = q_dirBN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirBN[k]  - nz_dirBN[k])/(sqrtf(two))); count += one;
   //   q = q_dirTS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirTS[k]) + nz_dirTS[k])/(sqrtf(two))); count += one;
   //   q = q_dirTNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTSW[k])- ny_dirTSW[k] + nz_dirTSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTNW[k])+ ny_dirTNW[k] + nz_dirTNW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBSW[k])- ny_dirBSW[k] - nz_dirBSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBNW[k])+ ny_dirBNW[k] - nz_dirBNW[k])/(sqrtf(three))); count += one;
      //real qMed = sum/count;
      //real phi = fac / (qMed + fac);
      //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real testQ = c2o1;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirE[k]) + fac);
         VeloX *= phi;
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c6o1*c2o27*( VeloX     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirW[k]) + fac);
         VeloX *= phi;
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c6o1*c2o27*(-VeloX     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( ny_dirN[k]) + fac);
         VeloY *= phi;
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[d0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c6o1*c2o27*( VeloY     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-ny_dirS[k]) + fac);
         VeloY *= phi;
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[d0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c6o1*c2o27*(-VeloY     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nz_dirT[k]) + fac);
         VeloZ *= phi;
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c6o1*c2o27*( VeloZ     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nz_dirB[k]) + fac);
         VeloZ *= phi;
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c6o1*c2o27*(-VeloZ     ))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k];
         VeloY = vx2 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirNE[k] + ny_dirNE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dMM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c6o1*c1o54*(VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k];
         VeloY = vx2 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirSW[k] - ny_dirSW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dPP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c6o1*c1o54*(-VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k];
         VeloY = vx2 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirSE[k] - ny_dirSE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dMP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c6o1*c1o54*( VeloX-VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k];
         VeloY = vx2 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirNW[k] + ny_dirNW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dPM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c6o1*c1o54*(-VeloX+VeloY))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k];
         VeloZ = vx3 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirTE[k] + nz_dirTE[k]) + fac);
         VeloX *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dM0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c6o1*c1o54*( VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k];
         VeloZ = vx3 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirBW[k] - nz_dirBW[k]) + fac);
         VeloX *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dP0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c6o1*c1o54*(-VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k];
         VeloZ = vx3 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirBE[k] - nz_dirBE[k]) + fac);
         VeloX *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dM0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c6o1*c1o54*( VeloX-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k];
         VeloZ = vx3 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirTW[k] + nz_dirTW[k]) + fac);
         VeloX *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dP0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c6o1*c1o54*(-VeloX+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k];
         VeloZ = vx3 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( ny_dirTN[k] + nz_dirTN[k]) + fac);
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c6o1*c1o54*( VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k];
         VeloZ = vx3 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-ny_dirBS[k] - nz_dirBS[k]) + fac);
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c6o1*c1o54*( -VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k];
         VeloZ = vx3 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( ny_dirBN[k] - nz_dirBN[k]) + fac);
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c6o1*c1o54*( VeloY-VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k];
         VeloZ = vx3 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-ny_dirTS[k] + nz_dirTS[k]) + fac);
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c6o1*c1o54*( -VeloY+VeloZ))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k];
         VeloY = vx2 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k];
         VeloZ = vx3 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c6o1*c1o216*( VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k];
         VeloY = vx2 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k];
         VeloZ = vx3 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirBSW[k] - ny_dirBSW[k] - nz_dirBSW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c6o1*c1o216*(-VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k];
         VeloY = vx2 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k];
         VeloZ = vx3 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c6o1*c1o216*( VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k];
         VeloY = vx2 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k];
         VeloZ = vx3 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirTSW[k] - ny_dirTSW[k] + nz_dirTSW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c6o1*c1o216*(-VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k];
         VeloY = vx2 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k];
         VeloZ = vx3 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c6o1*c1o216*( VeloX-VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k];
         VeloY = vx2 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k];
         VeloZ = vx3 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirBNW[k] + ny_dirBNW[k] - nz_dirBNW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c6o1*c1o216*(-VeloX+VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k];
         VeloY = vx2 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k];
         VeloZ = vx3 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = (phi > one) ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs( nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c6o1*c1o216*( VeloX-VeloY-VeloZ))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k];
         VeloY = vx2 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k];
         VeloZ = vx3 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k];
         //phi = fac * (one + magS / (c0p0000002+uk) * (one-q));
         //phi = phi > one ? one:phi;
         //phi = fac; //Test
         q = testQ; //AAAAHHHHHH bitte wieder raus nehmen!!!!
         phi = fac / (q * fabs(-nx_dirTNW[k] + ny_dirTNW[k] + nz_dirTNW[k]) + fac);
         VeloX *= phi;
         VeloY *= phi;
         VeloZ *= phi;
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c6o1*c1o216*(-VeloX+VeloY+VeloZ))/(c1o1+q) - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
__global__ void QSlipNormDeviceComp27(
    real* DD, 
    int* k_Q, 
    real* QQ,
    unsigned int  numberOfBCnodes,
    real om1, 
    real* NormalX,
    real* NormalY,
    real* NormalZ,
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

   if(k< numberOfBCnodes)
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
      real *nx_dirE,   *nx_dirW,   *nx_dirN,   *nx_dirS,   *nx_dirT,   *nx_dirB, 
              *nx_dirNE,  *nx_dirSW,  *nx_dirSE,  *nx_dirNW,  *nx_dirTE,  *nx_dirBW,
              *nx_dirBE,  *nx_dirTW,  *nx_dirTN,  *nx_dirBS,  *nx_dirBN,  *nx_dirTS,
              *nx_dirTNE, *nx_dirTSW, *nx_dirTSE, *nx_dirTNW, *nx_dirBNE, *nx_dirBSW,
              *nx_dirBSE, *nx_dirBNW; 
      nx_dirE   = &NormalX[dP00 * numberOfBCnodes];
      nx_dirW   = &NormalX[dM00 * numberOfBCnodes];
      nx_dirN   = &NormalX[d0P0 * numberOfBCnodes];
      nx_dirS   = &NormalX[d0M0 * numberOfBCnodes];
      nx_dirT   = &NormalX[d00P * numberOfBCnodes];
      nx_dirB   = &NormalX[d00M * numberOfBCnodes];
      nx_dirNE  = &NormalX[dPP0 * numberOfBCnodes];
      nx_dirSW  = &NormalX[dMM0 * numberOfBCnodes];
      nx_dirSE  = &NormalX[dPM0 * numberOfBCnodes];
      nx_dirNW  = &NormalX[dMP0 * numberOfBCnodes];
      nx_dirTE  = &NormalX[dP0P * numberOfBCnodes];
      nx_dirBW  = &NormalX[dM0M * numberOfBCnodes];
      nx_dirBE  = &NormalX[dP0M * numberOfBCnodes];
      nx_dirTW  = &NormalX[dM0P * numberOfBCnodes];
      nx_dirTN  = &NormalX[d0PP * numberOfBCnodes];
      nx_dirBS  = &NormalX[d0MM * numberOfBCnodes];
      nx_dirBN  = &NormalX[d0PM * numberOfBCnodes];
      nx_dirTS  = &NormalX[d0MP * numberOfBCnodes];
      nx_dirTNE = &NormalX[dPPP * numberOfBCnodes];
      nx_dirTSW = &NormalX[dMMP * numberOfBCnodes];
      nx_dirTSE = &NormalX[dPMP * numberOfBCnodes];
      nx_dirTNW = &NormalX[dMPP * numberOfBCnodes];
      nx_dirBNE = &NormalX[dPPM * numberOfBCnodes];
      nx_dirBSW = &NormalX[dMMM * numberOfBCnodes];
      nx_dirBSE = &NormalX[dPMM * numberOfBCnodes];
      nx_dirBNW = &NormalX[dMPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *ny_dirE,   *ny_dirW,   *ny_dirN,   *ny_dirS,   *ny_dirT,   *ny_dirB, 
              *ny_dirNE,  *ny_dirSW,  *ny_dirSE,  *ny_dirNW,  *ny_dirTE,  *ny_dirBW,
              *ny_dirBE,  *ny_dirTW,  *ny_dirTN,  *ny_dirBS,  *ny_dirBN,  *ny_dirTS,
              *ny_dirTNE, *ny_dirTSW, *ny_dirTSE, *ny_dirTNW, *ny_dirBNE, *ny_dirBSW,
              *ny_dirBSE, *ny_dirBNW; 
      ny_dirE   = &NormalY[dP00 * numberOfBCnodes];
      ny_dirW   = &NormalY[dM00 * numberOfBCnodes];
      ny_dirN   = &NormalY[d0P0 * numberOfBCnodes];
      ny_dirS   = &NormalY[d0M0 * numberOfBCnodes];
      ny_dirT   = &NormalY[d00P * numberOfBCnodes];
      ny_dirB   = &NormalY[d00M * numberOfBCnodes];
      ny_dirNE  = &NormalY[dPP0 * numberOfBCnodes];
      ny_dirSW  = &NormalY[dMM0 * numberOfBCnodes];
      ny_dirSE  = &NormalY[dPM0 * numberOfBCnodes];
      ny_dirNW  = &NormalY[dMP0 * numberOfBCnodes];
      ny_dirTE  = &NormalY[dP0P * numberOfBCnodes];
      ny_dirBW  = &NormalY[dM0M * numberOfBCnodes];
      ny_dirBE  = &NormalY[dP0M * numberOfBCnodes];
      ny_dirTW  = &NormalY[dM0P * numberOfBCnodes];
      ny_dirTN  = &NormalY[d0PP * numberOfBCnodes];
      ny_dirBS  = &NormalY[d0MM * numberOfBCnodes];
      ny_dirBN  = &NormalY[d0PM * numberOfBCnodes];
      ny_dirTS  = &NormalY[d0MP * numberOfBCnodes];
      ny_dirTNE = &NormalY[dPPP * numberOfBCnodes];
      ny_dirTSW = &NormalY[dMMP * numberOfBCnodes];
      ny_dirTSE = &NormalY[dPMP * numberOfBCnodes];
      ny_dirTNW = &NormalY[dMPP * numberOfBCnodes];
      ny_dirBNE = &NormalY[dPPM * numberOfBCnodes];
      ny_dirBSW = &NormalY[dMMM * numberOfBCnodes];
      ny_dirBSE = &NormalY[dPMM * numberOfBCnodes];
      ny_dirBNW = &NormalY[dMPM * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      real *nz_dirE,   *nz_dirW,   *nz_dirN,   *nz_dirS,   *nz_dirT,   *nz_dirB, 
              *nz_dirNE,  *nz_dirSW,  *nz_dirSE,  *nz_dirNW,  *nz_dirTE,  *nz_dirBW,
              *nz_dirBE,  *nz_dirTW,  *nz_dirTN,  *nz_dirBS,  *nz_dirBN,  *nz_dirTS,
              *nz_dirTNE, *nz_dirTSW, *nz_dirTSE, *nz_dirTNW, *nz_dirBNE, *nz_dirBSW,
              *nz_dirBSE, *nz_dirBNW; 
      nz_dirE   = &NormalZ[dP00 * numberOfBCnodes];
      nz_dirW   = &NormalZ[dM00 * numberOfBCnodes];
      nz_dirN   = &NormalZ[d0P0 * numberOfBCnodes];
      nz_dirS   = &NormalZ[d0M0 * numberOfBCnodes];
      nz_dirT   = &NormalZ[d00P * numberOfBCnodes];
      nz_dirB   = &NormalZ[d00M * numberOfBCnodes];
      nz_dirNE  = &NormalZ[dPP0 * numberOfBCnodes];
      nz_dirSW  = &NormalZ[dMM0 * numberOfBCnodes];
      nz_dirSE  = &NormalZ[dPM0 * numberOfBCnodes];
      nz_dirNW  = &NormalZ[dMP0 * numberOfBCnodes];
      nz_dirTE  = &NormalZ[dP0P * numberOfBCnodes];
      nz_dirBW  = &NormalZ[dM0M * numberOfBCnodes];
      nz_dirBE  = &NormalZ[dP0M * numberOfBCnodes];
      nz_dirTW  = &NormalZ[dM0P * numberOfBCnodes];
      nz_dirTN  = &NormalZ[d0PP * numberOfBCnodes];
      nz_dirBS  = &NormalZ[d0MM * numberOfBCnodes];
      nz_dirBN  = &NormalZ[d0PM * numberOfBCnodes];
      nz_dirTS  = &NormalZ[d0MP * numberOfBCnodes];
      nz_dirTNE = &NormalZ[dPPP * numberOfBCnodes];
      nz_dirTSW = &NormalZ[dMMP * numberOfBCnodes];
      nz_dirTSE = &NormalZ[dPMP * numberOfBCnodes];
      nz_dirTNW = &NormalZ[dMPP * numberOfBCnodes];
      nz_dirBNE = &NormalZ[dPPM * numberOfBCnodes];
      nz_dirBSW = &NormalZ[dMMM * numberOfBCnodes];
      nz_dirBSE = &NormalZ[dPMM * numberOfBCnodes];
      nz_dirBNW = &NormalZ[dMPM * numberOfBCnodes];
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
      real f_W    = (D.f[dP00])[ke   ];
      real f_E    = (D.f[dM00])[kw   ];
      real f_S    = (D.f[d0P0])[kn   ];
      real f_N    = (D.f[d0M0])[ks   ];
      real f_B    = (D.f[d00P])[kt   ];
      real f_T    = (D.f[d00M])[kb   ];
      real f_SW   = (D.f[dPP0])[kne  ];
      real f_NE   = (D.f[dMM0])[ksw  ];
      real f_NW   = (D.f[dPM0])[kse  ];
      real f_SE   = (D.f[dMP0])[knw  ];
      real f_BW   = (D.f[dP0P])[kte  ];
      real f_TE   = (D.f[dM0M])[kbw  ];
      real f_TW   = (D.f[dP0M])[kbe  ];
      real f_BE   = (D.f[dM0P])[ktw  ];
      real f_BS   = (D.f[d0PP])[ktn  ];
      real f_TN   = (D.f[d0MM])[kbs  ];
      real f_TS   = (D.f[d0PM])[kbn  ];
      real f_BN   = (D.f[d0MP])[kts  ];
      real f_BSW  = (D.f[dPPP])[ktne ];
      real f_BNE  = (D.f[dMMP])[ktsw ];
      real f_BNW  = (D.f[dPMP])[ktse ];
      real f_BSE  = (D.f[dMPP])[ktnw ];
      real f_TSW  = (D.f[dPPM])[kbne ];
      real f_TNE  = (D.f[dMMM])[kbsw ];
      real f_TNW  = (D.f[dPMM])[kbse ];
      real f_TSE  = (D.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3, drho, feq, q;
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]); 

      vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                (f_E - f_W)) / (c1o1 + drho); 
         

      vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                 ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                 (f_N - f_S)) / (c1o1 + drho); 

      vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                 (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                 (f_T - f_B)) / (c1o1 + drho); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);

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
      real VeloX = vx1;
      real VeloY = vx2;
      real VeloZ = vx3;
      real fac = c1o100;//0.5;
       //real phi = c0o1;
      real alpha = c1o100;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real kxyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE ) / (c1o1 + drho) - ((vx1*vx2)));
      real kyzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW ) / (c1o1 + drho) - ((vx2*vx3)));
      real kxzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE ) / (c1o1 + drho) - ((vx1*vx3)));

      real kxxFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_E+f_NE+f_SE+f_TE+f_BE+f_W+f_NW+f_SW+f_TW+f_BW+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx1*vx1)));
      real kyyFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_N+f_NE+f_NW+f_TN+f_BN+f_S+f_SE+f_SW+f_TS+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx2*vx2)));
      real kzzFromfcNEQ = -(c3o1 * om1 / (c1o1-om1))*((f_T+f_TE+f_TW+f_TN+f_BS+f_B+f_BE+f_BW+f_BN+f_BS+ f_TNE+f_TSE+f_BNE+f_TNE+ f_TNW+f_TSW+f_BNW+f_TNW ) / (c1o1 + drho) - ((c1o3*drho + vx3*vx3)));

      real magS = sqrtf(kxyFromfcNEQ*kxyFromfcNEQ + kyzFromfcNEQ*kyzFromfcNEQ + kxzFromfcNEQ*kxzFromfcNEQ + kxxFromfcNEQ*kxxFromfcNEQ + kyyFromfcNEQ*kyyFromfcNEQ + kzzFromfcNEQ*kzzFromfcNEQ);

      fac = fac * magS / (c1o3 * (c1o1 / om1 - c1o2));
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real *facAst = &QQ[d000 * numberOfBCnodes];

      fac = fac * alpha + facAst[k] * (c1o1 - alpha);
      facAst[k] = fac;
      //(&QQ[d000 * numberOfBCnodes])[KQK] = fac;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////real uk = sqrtf(vx1*vx1 + vx2*vx2 + vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real phi = expf(magS/0.01f) - one;
      //phi = (phi > one) ? one:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real C = five;
      //real kappa = 0.41f;
      //real phi = (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))) - one) / (C * kappa * c1o2 * logf(magS / (c1o3 * (one / om1 - c1o2))));
      //phi = (phi < zero) ? zero:phi;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //real sum = zero, count = zero;
   //   q = q_dirE   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirE[k] ); count += one;
   //   q = q_dirW   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirW[k])); count += one;
   //   q = q_dirN   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirN[k] ); count += one;
   //   q = q_dirS   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirS[k])); count += one;
   //   q = q_dirT   [k]; if (q>=zero && q<=one) sum += (q *   nx_dirT[k] ); count += one;
   //   q = q_dirB   [k]; if (q>=zero && q<=one) sum += (q * (-nx_dirB[k])); count += one;
   //   q = q_dirNE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirNE[k]  + ny_dirNE[k])/(sqrtf(two))); count += one;
   //   q = q_dirSW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirSW[k]) - ny_dirSW[k])/(sqrtf(two))); count += one;
   //   q = q_dirSE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirSE[k]  - ny_dirSE[k])/(sqrtf(two))); count += one;
   //   q = q_dirNW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirNW[k]) + ny_dirNW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTE[k]  + nz_dirTE[k])/(sqrtf(two))); count += one;
   //   q = q_dirBW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBW[k]) - nz_dirBW[k])/(sqrtf(two))); count += one;
   //   q = q_dirBE  [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBE[k]  - nz_dirBE[k])/(sqrtf(two))); count += one;
   //   q = q_dirTW  [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTW[k]) + nz_dirTW[k])/(sqrtf(two))); count += one;
   //   q = q_dirTN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirTN[k]  + nz_dirTN[k])/(sqrtf(two))); count += one;
   //   q = q_dirBS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirBS[k]) - nz_dirBS[k])/(sqrtf(two))); count += one;
   //   q = q_dirBN  [k]; if (q>=zero && q<=one) sum += (q * (  ny_dirBN[k]  - nz_dirBN[k])/(sqrtf(two))); count += one;
   //   q = q_dirTS  [k]; if (q>=zero && q<=one) sum += (q * ((-ny_dirTS[k]) + nz_dirTS[k])/(sqrtf(two))); count += one;
   //   q = q_dirTNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTSW[k])- ny_dirTSW[k] + nz_dirTSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirTSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirTNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirTNW[k])+ ny_dirTNW[k] + nz_dirTNW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBSW[k])- ny_dirBSW[k] - nz_dirBSW[k])/(sqrtf(three))); count += one;
   //   q = q_dirBSE [k]; if (q>=zero && q<=one) sum += (q * (  nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k])/(sqrtf(three))); count += one;
   //   q = q_dirBNW [k]; if (q>=zero && q<=one) sum += (q * ((-nx_dirBNW[k])+ ny_dirBNW[k] - nz_dirBNW[k])/(sqrtf(three))); count += one;
      //real qMed = sum/count;
      //real phi = fac / (qMed + fac);
      //phi = (phi > one) ? one:one;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real sliplength = 0.9f;//c1o2;
      real qSlip = c0o1;
      real un = c0o1;
      real ut = c0o1;
      real tangential = c0o1;
      //real cSmallSingle = c0p0000002;

      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k];
         un = fabs((vx1 * nx_dirE[k] + vx2 * ny_dirE[k] + vx3 * nz_dirE[k]) * nx_dirE[k]);
         ut = fabs(VeloX);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*( vx1        )*/+c9o2*( vx1        )*( vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dM00])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k];
         un = fabs(-(vx1 * nx_dirW[k] + vx2 * ny_dirW[k] + vx3 * nz_dirW[k]) * nx_dirW[k]);
         ut = fabs(-VeloX);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*(-vx1        )*/+c9o2*(-vx1        )*(-vx1        ) * (c1o1 + drho)-cu_sq); 
         (D.f[dP00])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k];
         un = fabs( (vx1 * nx_dirN[k] + vx2 * ny_dirN[k] + vx3 * nz_dirN[k]) * ny_dirN[k]);
         ut = fabs( VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( ny_dirN[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*(    vx2     )*/+c9o2*(     vx2    )*(     vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[d0M0])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k];
         un = fabs(-(vx1 * nx_dirS[k] + vx2 * ny_dirS[k] + vx3 * nz_dirS[k]) * ny_dirS[k]);
         ut = fabs(-VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-ny_dirS[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*(   -vx2     )*/+c9o2*(    -vx2    )*(    -vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[d0P0])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = vx3 - (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k];
         un = fabs( (vx1 * nx_dirT[k] + vx2 * ny_dirT[k] + vx3 * nz_dirT[k]) * nz_dirT[k]);
         ut = fabs( VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nz_dirT[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*(         vx3)*/+c9o2*(         vx3)*(         vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d00M])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloZ = vx3 - (vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k];
         un = fabs(-(vx1 * nx_dirB[k] + vx2 * ny_dirB[k] + vx3 * nz_dirB[k]) * nz_dirB[k]);
         ut = fabs(-VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nz_dirB[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c2o27* (drho/*+three*(        -vx3)*/+c9o2*(        -vx3)*(        -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d00P])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T))/(c1o1+q) - c2o27 * drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k];
         VeloY = vx2 - (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k];
         un = fabs( (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * nx_dirNE[k] + (vx1 * nx_dirNE[k] + vx2 * ny_dirNE[k] + vx3 * nz_dirNE[k]) * ny_dirNE[k]);
         ut = fabs( VeloX + VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirNE[k]+ny_dirNE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*( vx1+vx2    )*/+c9o2*( vx1+vx2    )*( vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dMM0])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k];
         VeloY = vx2 - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k];
         un = fabs(-(vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * nx_dirSW[k] - (vx1 * nx_dirSW[k] + vx2 * ny_dirSW[k] + vx3 * nz_dirSW[k]) * ny_dirSW[k]);
         ut = fabs(-VeloX - VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirSW[k]-ny_dirSW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1-vx2    )*/+c9o2*(-vx1-vx2    )*(-vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dPP0])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k];
         VeloY = vx2 - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k];
         un = fabs( (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * nx_dirSE[k] - (vx1 * nx_dirSE[k] + vx2 * ny_dirSE[k] + vx3 * nz_dirSE[k]) * ny_dirSE[k]);
         ut = fabs( VeloX - VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirSE[k]-ny_dirSE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*( vx1-vx2    )*/+c9o2*( vx1-vx2    )*( vx1-vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dMP0])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k];
         VeloY = vx2 - (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k];
         un = fabs(-(vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * nx_dirNW[k] + (vx1 * nx_dirNW[k] + vx2 * ny_dirNW[k] + vx3 * nz_dirNW[k]) * ny_dirNW[k]);
         ut = fabs(-VeloX + VeloY);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirNW[k]+ny_dirNW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1+vx2    )*/+c9o2*(-vx1+vx2    )*(-vx1+vx2    ) * (c1o1 + drho)-cu_sq); 
         (D.f[dPM0])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k];
         VeloZ = vx3 - (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k];
         un = fabs( (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nx_dirTE[k] + (vx1 * nx_dirTE[k] + vx2 * ny_dirTE[k] + vx3 * nz_dirTE[k]) * nz_dirTE[k]);
         ut = fabs( VeloX + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirTE[k]+nz_dirTE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*( vx1    +vx3)*/+c9o2*( vx1    +vx3)*( vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dM0M])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k];
         VeloZ = vx3 - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k];
         un = fabs(-(vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nx_dirBW[k] - (vx1 * nx_dirBW[k] + vx2 * ny_dirBW[k] + vx3 * nz_dirBW[k]) * nz_dirBW[k]);
         ut = fabs(-VeloX - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirBW[k]-nz_dirBW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1    -vx3)*/+c9o2*(-vx1    -vx3)*(-vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dP0P])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k];
         VeloZ = vx3 - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k];
         un = fabs( (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nx_dirBE[k] - (vx1 * nx_dirBE[k] + vx2 * ny_dirBE[k] + vx3 * nz_dirBE[k]) * nz_dirBE[k]);
         ut = fabs( VeloX - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirBE[k]-nz_dirBE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*( vx1    -vx3)*/+c9o2*( vx1    -vx3)*( vx1    -vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dM0P])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k];
         VeloZ = vx3 - (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k];
         un = fabs(-(vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nx_dirTW[k] + (vx1 * nx_dirTW[k] + vx2 * ny_dirTW[k] + vx3 * nz_dirTW[k]) * nz_dirTW[k]);
         ut = fabs(-VeloX + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirTW[k]+nz_dirTW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(-vx1    +vx3)*/+c9o2*(-vx1    +vx3)*(-vx1    +vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dP0M])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k];
         VeloZ = vx3 - (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k];
         un = fabs( (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * ny_dirTN[k] + (vx1 * nx_dirTN[k] + vx2 * ny_dirTN[k] + vx3 * nz_dirTN[k]) * nz_dirTN[k]);
         ut = fabs( VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( ny_dirTN[k]+nz_dirTN[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(     vx2+vx3)*/+c9o2*(     vx2+vx3)*(     vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0MM])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k];
         VeloZ = vx3 - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k];
         un = fabs(-(vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * ny_dirBS[k] - (vx1 * nx_dirBS[k] + vx2 * ny_dirBS[k] + vx3 * nz_dirBS[k]) * nz_dirBS[k]);
         ut = fabs(-VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-ny_dirBS[k]-nz_dirBS[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(    -vx2-vx3)*/+c9o2*(    -vx2-vx3)*(    -vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0PP])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k];
         VeloZ = vx3 - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k];
         un = fabs( (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * ny_dirBN[k] - (vx1 * nx_dirBN[k] + vx2 * ny_dirBN[k] + vx3 * nz_dirBN[k]) * nz_dirBN[k]);
         ut = fabs( VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( ny_dirBN[k]-nz_dirBN[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(     vx2-vx3)*/+c9o2*(     vx2-vx3)*(     vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0MP])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloY = vx2 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k];
         VeloZ = vx3 - (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k];
         un = fabs(-(vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * ny_dirTS[k] + (vx1 * nx_dirTS[k] + vx2 * ny_dirTS[k] + vx3 * nz_dirTS[k]) * nz_dirTS[k]);
         ut = fabs(-VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-ny_dirTS[k]+nz_dirTS[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o54* (drho/*+three*(    -vx2+vx3)*/+c9o2*(    -vx2+vx3)*(    -vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[d0PM])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN))/(c1o1+q) - c1o54 * drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k];
         VeloY = vx2 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k];
         VeloZ = vx3 - (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k];
         un = fabs( (vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nx_dirTNE[k] 
                   +(vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * ny_dirTNE[k] 
                   +(vx1 * nx_dirTNE[k] + vx2 * ny_dirTNE[k] + vx3 * nz_dirTNE[k]) * nz_dirTNE[k]);
         ut = fabs( VeloX + VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirTNE[k] + ny_dirTNE[k] + nz_dirTNE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*( vx1+vx2+vx3)*/+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMMM])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k];
         VeloY = vx2 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k];
         VeloZ = vx3 - (vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k];
         un = fabs(-(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nx_dirBSW[k] 
                   -(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * ny_dirBSW[k] 
                   -(vx1 * nx_dirBSW[k] + vx2 * ny_dirBSW[k] + vx3 * nz_dirBSW[k]) * nz_dirBSW[k]);
         ut = fabs(-VeloX - VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirBSW[k] - ny_dirBSW[k] - nz_dirBSW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1-vx2-vx3)*/+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPPP])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k];
         VeloY = vx2 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k];
         VeloZ = vx3 - (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k];
         un = fabs( (vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nx_dirBNE[k] 
                   +(vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * ny_dirBNE[k] 
                   -(vx1 * nx_dirBNE[k] + vx2 * ny_dirBNE[k] + vx3 * nz_dirBNE[k]) * nz_dirBNE[k]);
         ut = fabs( VeloX + VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirBNE[k] + ny_dirBNE[k] - nz_dirBNE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*( vx1+vx2-vx3)*/+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMMP])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k];
         VeloY = vx2 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k];
         VeloZ = vx3 - (vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k];
         un = fabs(-(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nx_dirTSW[k] 
                   -(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * ny_dirTSW[k] 
                   +(vx1 * nx_dirTSW[k] + vx2 * ny_dirTSW[k] + vx3 * nz_dirTSW[k]) * nz_dirTSW[k]);
         ut = fabs(-VeloX - VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirTSW[k] - ny_dirTSW[k] + nz_dirTSW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1-vx2+vx3)*/+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPPM])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k];
         VeloY = vx2 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k];
         VeloZ = vx3 - (vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k];
         un = fabs(+(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nx_dirTSE[k] 
                   -(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * ny_dirTSE[k] 
                   +(vx1 * nx_dirTSE[k] + vx2 * ny_dirTSE[k] + vx3 * nz_dirTSE[k]) * nz_dirTSE[k]);
         ut = fabs(+VeloX - VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirTSE[k] - ny_dirTSE[k] + nz_dirTSE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*( vx1-vx2+vx3)*/+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMPM])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k];
         VeloY = vx2 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k];
         VeloZ = vx3 - (vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k];
         un = fabs(-(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nx_dirBNW[k] 
                   +(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * ny_dirBNW[k] 
                   -(vx1 * nx_dirBNW[k] + vx2 * ny_dirBNW[k] + vx3 * nz_dirBNW[k]) * nz_dirBNW[k]);
         ut = fabs(-VeloX + VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirBNW[k] + ny_dirBNW[k] - nz_dirBNW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1+vx2-vx3)*/+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPMP])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k];
         VeloY = vx2 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k];
         VeloZ = vx3 - (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k];
         un = fabs( (vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nx_dirBSE[k] 
                   -(vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * ny_dirBSE[k] 
                   -(vx1 * nx_dirBSE[k] + vx2 * ny_dirBSE[k] + vx3 * nz_dirBSE[k]) * nz_dirBSE[k]);
         ut = fabs( VeloX - VeloY - VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs( nx_dirBSE[k] - ny_dirBSE[k] - nz_dirBSE[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*( vx1-vx2-vx3)*/+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dMPP])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW))/(c1o1+q) - c1o216 * drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         VeloX = vx1 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k];
         VeloY = vx2 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k];
         VeloZ = vx3 - (vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k];
         un = fabs(-(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nx_dirTNW[k] 
                   +(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * ny_dirTNW[k] 
                   +(vx1 * nx_dirTNW[k] + vx2 * ny_dirTNW[k] + vx3 * nz_dirTNW[k]) * nz_dirTNW[k]);
         ut = fabs(-VeloX + VeloY + VeloZ);
         tangential = ut / (ut + un + cSmallSingle);
         qSlip = sliplength * fabs(-nx_dirTNW[k] + ny_dirTNW[k] + nz_dirTNW[k]);        //sliplength * e_i * n_i
         //qSlip = (qSlip < zero) ? zero:qSlip;
         //tangential = (tangential > one) ? one:tangential;
         q = (q + qSlip)/(c1o1 + qSlip * (c1o1 - tangential) / (cSmallSingle + q));
         feq=c1o216*(drho/*+three*(-vx1+vx2+vx3)*/+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3) * (c1o1 + drho)-cu_sq); 
         (D.f[dPMM])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE))/(c1o1+q) - c1o216 * drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
