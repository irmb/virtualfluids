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
#include "basics/constants/NumericConstants.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__ void AdvectionDiffusionDirichlet_Device(
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
      real f_ZERO = (D.f[d000])[kzero];
      real f_BSW  = (D.f[dPPP])[ktne ];
      real f_BNE  = (D.f[dMMP])[ktsw ];
      real f_BNW  = (D.f[dPMP])[ktse ];
      real f_BSE  = (D.f[dMPP])[ktnw ];
      real f_TSW  = (D.f[dPPM])[kbne ];
      real f_TNE  = (D.f[dMMM])[kbsw ];
      real f_TNW  = (D.f[dPMM])[kbse ];
      real f_TSE  = (D.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real rho0   =  (f_TNE+f_BSW)+(f_TSW+f_BNE)+(f_TSE+f_BNW)+(f_TNW+f_BSE)+(f_NE+f_SW)+(f_NW+f_SE)+(f_TE+f_BW)+(f_BE+f_TW)+(f_TN+f_BS)+(f_BN+f_TS)+(f_E+f_W)+(f_N+f_S)+(f_T+f_B)+f_ZERO;
      real rho    =  rho0 + c1o1;
      real OORho  =  c1o1/rho;
      real vx1    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)+(f_BSE-f_TNW) +(f_NE-f_SW)+(f_SE-f_NW)+(f_TE-f_BW)+(f_BE-f_TW)+(f_E-f_W));
      real vx2    =  OORho*((f_TNE-f_BSW)+(f_BNE-f_TSW)+(f_BNW-f_TSE)+(f_TNW-f_BSE) +(f_NE-f_SW)+(f_NW-f_SE)+(f_TN-f_BS)+(f_BN-f_TS)+(f_N-f_S));
      real vx3    =  OORho*((f_TNE-f_BSW)+(f_TSW-f_BNE)+(f_TSE-f_BNW)+(f_TNW-f_BSE) +(f_TE-f_BW)+(f_TW-f_BE)+(f_TN-f_BS)+(f_TS-f_BN)+(f_T-f_B));
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
      real f27_ZERO = (D27.f[d000])[kzero];
      real f27_BSW  = (D27.f[dPPP])[ktne ];
      real f27_BNE  = (D27.f[dMMP])[ktsw ];
      real f27_BNW  = (D27.f[dPMP])[ktse ];
      real f27_BSE  = (D27.f[dMPP])[ktnw ];
      real f27_TSW  = (D27.f[dPPM])[kbne ];
      real f27_TNE  = (D27.f[dMMM])[kbsw ];
      real f27_TNW  = (D27.f[dPMM])[kbse ];
      real f27_TSE  = (D27.f[dMPM])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      ////////////////////////////////////////////////////////////////////////////////
      real ConcD =   f27_TSE + f27_TNW + f27_TNE + f27_TSW + f27_BSE + f27_BNW + f27_BNE + f27_BSW +
         f27_BN  + f27_TS  + f27_TN  + f27_BS  + f27_BE  + f27_TW  + f27_TE  + f27_BW  + f27_SE + f27_NW + f27_NE + f27_SW + 
         f27_T   + f27_B   + f27_N   + f27_S   + f27_E   + f27_W   + f27_ZERO; 

      real feq27_E    =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feq27_W    =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feq27_N    =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feq27_S    =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feq27_T    =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feq27_B    =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feq27_NE   =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feq27_SW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feq27_SE   =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feq27_NW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feq27_TE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feq27_BW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feq27_BE   =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feq27_TW   =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feq27_TN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feq27_BS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feq27_BN   =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feq27_TS   =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feq27_TNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feq27_BSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feq27_BNE  =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feq27_TSW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feq27_TSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feq27_BNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feq27_BSE  =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feq27_TNW  =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real TempD = temp[k];

      real feqW27_E    =   c2o27* TempD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      real feqW27_W    =   c2o27* TempD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      real feqW27_N    =   c2o27* TempD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      real feqW27_S    =   c2o27* TempD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      real feqW27_T    =   c2o27* TempD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      real feqW27_B    =   c2o27* TempD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      real feqW27_NE   =   c1o54* TempD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
      real feqW27_SW   =   c1o54* TempD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
      real feqW27_SE   =   c1o54* TempD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
      real feqW27_NW   =   c1o54* TempD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
      real feqW27_TE   =   c1o54* TempD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
      real feqW27_BW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
      real feqW27_BE   =   c1o54* TempD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
      real feqW27_TW   =   c1o54* TempD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
      real feqW27_TN   =   c1o54* TempD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
      real feqW27_BS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
      real feqW27_BN   =   c1o54* TempD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
      real feqW27_TS   =   c1o54* TempD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
      real feqW27_TNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
      real feqW27_BSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
      real feqW27_BNE  =   c1o216*TempD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
      real feqW27_TSW  =   c1o216*TempD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
      real feqW27_TSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
      real feqW27_BNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
      real feqW27_BSE  =   c1o216*TempD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
      real feqW27_TNW  =   c1o216*TempD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
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
      real omegaD = c3o1 - sqrt(c3o1);
      real q;

      q = q_dirE[  ke   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dM00])[kw  ]=(c2o1*feqW27_W  -(f27_E  *(q*omegaD-c1o1)-omegaD*feq27_E  *(q-c1o1))/(omegaD-c1o1)+f27_W  *q)/(q+c1o1);
      q = q_dirW[  kw   ]; if (q>=c0o1 && q<=c1o1) (D27.f[dP00])[ke  ]=(c2o1*feqW27_E  -(f27_W  *(q*omegaD-c1o1)-omegaD*feq27_W  *(q-c1o1))/(omegaD-c1o1)+f27_E  *q)/(q+c1o1);
      q = q_dirN[  kn   ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0M0])[ks  ]=(c2o1*feqW27_S  -(f27_N  *(q*omegaD-c1o1)-omegaD*feq27_N  *(q-c1o1))/(omegaD-c1o1)+f27_S  *q)/(q+c1o1);
      q = q_dirS[  ks   ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0P0])[kn  ]=(c2o1*feqW27_N  -(f27_S  *(q*omegaD-c1o1)-omegaD*feq27_S  *(q-c1o1))/(omegaD-c1o1)+f27_N  *q)/(q+c1o1);
      q = q_dirT[  kt   ]; if (q>=c0o1 && q<=c1o1) (D27.f[d00M])[kb  ]=(c2o1*feqW27_B  -(f27_T  *(q*omegaD-c1o1)-omegaD*feq27_T  *(q-c1o1))/(omegaD-c1o1)+f27_B  *q)/(q+c1o1);
      q = q_dirB[  kb   ]; if (q>=c0o1 && q<=c1o1) (D27.f[d00P])[kt  ]=(c2o1*feqW27_T  -(f27_B  *(q*omegaD-c1o1)-omegaD*feq27_B  *(q-c1o1))/(omegaD-c1o1)+f27_T  *q)/(q+c1o1);
      q = q_dirNE[ kne  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMM0])[ksw ]=(c2o1*feqW27_SW -(f27_NE *(q*omegaD-c1o1)-omegaD*feq27_NE *(q-c1o1))/(omegaD-c1o1)+f27_SW *q)/(q+c1o1);
      q = q_dirSW[ ksw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPP0])[kne ]=(c2o1*feqW27_NE -(f27_SW *(q*omegaD-c1o1)-omegaD*feq27_SW *(q-c1o1))/(omegaD-c1o1)+f27_NE *q)/(q+c1o1);
      q = q_dirSE[ kse  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMP0])[knw ]=(c2o1*feqW27_NW -(f27_SE *(q*omegaD-c1o1)-omegaD*feq27_SE *(q-c1o1))/(omegaD-c1o1)+f27_NW *q)/(q+c1o1);
      q = q_dirNW[ knw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPM0])[kse ]=(c2o1*feqW27_SE -(f27_NW *(q*omegaD-c1o1)-omegaD*feq27_NW *(q-c1o1))/(omegaD-c1o1)+f27_SE *q)/(q+c1o1);
      q = q_dirTE[ kte  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dM0M])[kbw ]=(c2o1*feqW27_BW -(f27_TE *(q*omegaD-c1o1)-omegaD*feq27_TE *(q-c1o1))/(omegaD-c1o1)+f27_BW *q)/(q+c1o1);
      q = q_dirBW[ kbw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dP0P])[kte ]=(c2o1*feqW27_TE -(f27_BW *(q*omegaD-c1o1)-omegaD*feq27_BW *(q-c1o1))/(omegaD-c1o1)+f27_TE *q)/(q+c1o1);
      q = q_dirBE[ kbe  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dM0P])[ktw ]=(c2o1*feqW27_TW -(f27_BE *(q*omegaD-c1o1)-omegaD*feq27_BE *(q-c1o1))/(omegaD-c1o1)+f27_TW *q)/(q+c1o1);
      q = q_dirTW[ ktw  ]; if (q>=c0o1 && q<=c1o1) (D27.f[dP0M])[kbe ]=(c2o1*feqW27_BE -(f27_TW *(q*omegaD-c1o1)-omegaD*feq27_TW *(q-c1o1))/(omegaD-c1o1)+f27_BE *q)/(q+c1o1);
      q = q_dirTN[ ktn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0MM])[kbs ]=(c2o1*feqW27_BS -(f27_TN *(q*omegaD-c1o1)-omegaD*feq27_TN *(q-c1o1))/(omegaD-c1o1)+f27_BS *q)/(q+c1o1);
      q = q_dirBS[ kbs  ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0PP])[ktn ]=(c2o1*feqW27_TN -(f27_BS *(q*omegaD-c1o1)-omegaD*feq27_BS *(q-c1o1))/(omegaD-c1o1)+f27_TN *q)/(q+c1o1);
      q = q_dirBN[ kbn  ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0MP])[kts ]=(c2o1*feqW27_TS -(f27_BN *(q*omegaD-c1o1)-omegaD*feq27_BN *(q-c1o1))/(omegaD-c1o1)+f27_TS *q)/(q+c1o1);
      q = q_dirTS[ kts  ]; if (q>=c0o1 && q<=c1o1) (D27.f[d0PM])[kbn ]=(c2o1*feqW27_BN -(f27_TS *(q*omegaD-c1o1)-omegaD*feq27_TS *(q-c1o1))/(omegaD-c1o1)+f27_BN *q)/(q+c1o1);
      q = q_dirTNE[ktne ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMMM])[kbsw]=(c2o1*feqW27_BSW-(f27_TNE*(q*omegaD-c1o1)-omegaD*feq27_TNE*(q-c1o1))/(omegaD-c1o1)+f27_BSW*q)/(q+c1o1);
      q = q_dirBSW[kbsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPPP])[ktne]=(c2o1*feqW27_TNE-(f27_BSW*(q*omegaD-c1o1)-omegaD*feq27_BSW*(q-c1o1))/(omegaD-c1o1)+f27_TNE*q)/(q+c1o1);
      q = q_dirBNE[kbne ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMMP])[ktsw]=(c2o1*feqW27_TSW-(f27_BNE*(q*omegaD-c1o1)-omegaD*feq27_BNE*(q-c1o1))/(omegaD-c1o1)+f27_TSW*q)/(q+c1o1);
      q = q_dirTSW[ktsw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPPM])[kbne]=(c2o1*feqW27_BNE-(f27_TSW*(q*omegaD-c1o1)-omegaD*feq27_TSW*(q-c1o1))/(omegaD-c1o1)+f27_BNE*q)/(q+c1o1);
      q = q_dirTSE[ktse ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMPM])[kbnw]=(c2o1*feqW27_BNW-(f27_TSE*(q*omegaD-c1o1)-omegaD*feq27_TSE*(q-c1o1))/(omegaD-c1o1)+f27_BNW*q)/(q+c1o1);
      q = q_dirBNW[kbnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPMP])[ktse]=(c2o1*feqW27_TSE-(f27_BNW*(q*omegaD-c1o1)-omegaD*feq27_BNW*(q-c1o1))/(omegaD-c1o1)+f27_TSE*q)/(q+c1o1);
      q = q_dirBSE[kbse ]; if (q>=c0o1 && q<=c1o1) (D27.f[dMPP])[ktnw]=(c2o1*feqW27_TNW-(f27_BSE*(q*omegaD-c1o1)-omegaD*feq27_BSE*(q-c1o1))/(omegaD-c1o1)+f27_TNW*q)/(q+c1o1);
      q = q_dirTNW[ktnw ]; if (q>=c0o1 && q<=c1o1) (D27.f[dPMM])[kbse]=(c2o1*feqW27_BSE-(f27_TNW*(q*omegaD-c1o1)-omegaD*feq27_TNW*(q-c1o1))/(omegaD-c1o1)+f27_BSE*q)/(q+c1o1);
   }
}
