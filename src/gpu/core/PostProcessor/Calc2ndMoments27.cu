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
//======================================================================================
#include "Calc2ndMoments27.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "Calculation/Calculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc2ndMomentsIncompSP27(  real* kxyFromfcNEQ,
                                                        real* kyzFromfcNEQ,
                                                        real* kxzFromfcNEQ,
                                                        real* kxxMyyFromfcNEQ,
                                                        real* kxxMzzFromfcNEQ,
                                                        unsigned int* geoD,
                                                        unsigned int* neighborX,
                                                        unsigned int* neighborY,
                                                        unsigned int* neighborZ,
                                                        unsigned long long numberOfLBnodes,
                                                        real* DD,
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

   if(k < numberOfLBnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      //unsigned int kzero= k;
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
      //////////////////////////////////////////////////////////////////////////
      real        f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,/*f_ZERO,*/f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;
      f_E    = (D.f[dP00])[ke   ];
      f_W    = (D.f[dM00])[kw   ];
      f_N    = (D.f[d0P0])[kn   ];
      f_S    = (D.f[d0M0])[ks   ];
      f_T    = (D.f[d00P])[kt   ];
      f_B    = (D.f[d00M])[kb   ];
      f_NE   = (D.f[dPP0])[kne  ];
      f_SW   = (D.f[dMM0])[ksw  ];
      f_SE   = (D.f[dPM0])[kse  ];
      f_NW   = (D.f[dMP0])[knw  ];
      f_TE   = (D.f[dP0P])[kte  ];
      f_BW   = (D.f[dM0M])[kbw  ];
      f_BE   = (D.f[dP0M])[kbe  ];
      f_TW   = (D.f[dM0P])[ktw  ];
      f_TN   = (D.f[d0PP])[ktn  ];
      f_BS   = (D.f[d0MM])[kbs  ];
      f_BN   = (D.f[d0PM])[kbn  ];
      f_TS   = (D.f[d0MP])[kts  ];
      //f_ZERO = (D.f[d000])[kzero];
      f_TNE  = (D.f[dPPP])[ktne ];
      f_TSW  = (D.f[dMMP])[ktsw ];
      f_TSE  = (D.f[dPMP])[ktse ];
      f_TNW  = (D.f[dMPP])[ktnw ];
      f_BNE  = (D.f[dPPM])[kbne ];
      f_BSW  = (D.f[dMMM])[kbsw ];
      f_BSE  = (D.f[dPMM])[kbse ];
      f_BNW  = (D.f[dMPM])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3;
      kxyFromfcNEQ[k]       = c0o1;
      kyzFromfcNEQ[k]       = c0o1;
      kxzFromfcNEQ[k]       = c0o1;
      kxxMyyFromfcNEQ[k]    = c0o1;
      kxxMzzFromfcNEQ[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
          vx1                = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W);
          vx2                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S);
          vx3                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B);
          kxyFromfcNEQ[k]    = -c3o1 *(f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE-(vx1*vx2));
          kyzFromfcNEQ[k]    = -c3o1 *(f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW-(vx2*vx3));
          kxzFromfcNEQ[k]    = -c3o1 *(f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE-(vx1*vx3));
          kxxMyyFromfcNEQ[k] = -c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));        //all dP00+dM00 minus all d0P0+d0M0 (no combinations of xy left)
          kxxMzzFromfcNEQ[k] = -c3o2 * (f_SW+f_W+f_NW-f_BS-f_TS-f_B-f_T-f_BN-f_TN+f_SE+f_E+f_NE-(vx1*vx1-vx3*vx3));        //all dP00+dM00 minus all d00P+d00M (no combinations of xz left)
      }
   }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc2ndMomentsCompSP27(real* kxyFromfcNEQ,
                                                    real* kyzFromfcNEQ,
                                                    real* kxzFromfcNEQ,
                                                    real* kxxMyyFromfcNEQ,
                                                    real* kxxMzzFromfcNEQ,
                                                    unsigned int* geoD,
                                                    unsigned int* neighborX,
                                                    unsigned int* neighborY,
                                                    unsigned int* neighborZ,
                                                    unsigned long long numberOfLBnodes,
                                                    real* DD,
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

   if(k < numberOfLBnodes)
   {
      //////////////////////////////////////////////////////////////////////////
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
      //////////////////////////////////////////////////////////////////////////
      real f_ZERO;
      real        f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;
      f_E    = (D.f[dP00])[ke   ];
      f_W    = (D.f[dM00])[kw   ];
      f_N    = (D.f[d0P0])[kn   ];
      f_S    = (D.f[d0M0])[ks   ];
      f_T    = (D.f[d00P])[kt   ];
      f_B    = (D.f[d00M])[kb   ];
      f_NE   = (D.f[dPP0])[kne  ];
      f_SW   = (D.f[dMM0])[ksw  ];
      f_SE   = (D.f[dPM0])[kse  ];
      f_NW   = (D.f[dMP0])[knw  ];
      f_TE   = (D.f[dP0P])[kte  ];
      f_BW   = (D.f[dM0M])[kbw  ];
      f_BE   = (D.f[dP0M])[kbe  ];
      f_TW   = (D.f[dM0P])[ktw  ];
      f_TN   = (D.f[d0PP])[ktn  ];
      f_BS   = (D.f[d0MM])[kbs  ];
      f_BN   = (D.f[d0PM])[kbn  ];
      f_TS   = (D.f[d0MP])[kts  ];
      f_ZERO = (D.f[d000])[kzero];
      f_TNE  = (D.f[dPPP])[ktne ];
      f_TSW  = (D.f[dMMP])[ktsw ];
      f_TSE  = (D.f[dPMP])[ktse ];
      f_TNW  = (D.f[dMPP])[ktnw ];
      f_BNE  = (D.f[dPPM])[kbne ];
      f_BSW  = (D.f[dMMM])[kbsw ];
      f_BSE  = (D.f[dPMM])[kbse ];
      f_BNW  = (D.f[dMPM])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      real drho;
      real vx1, vx2, vx3, rho;
      kxyFromfcNEQ[k]       = c0o1;
      kyzFromfcNEQ[k]       = c0o1;
      kxzFromfcNEQ[k]       = c0o1;
      kxxMyyFromfcNEQ[k]    = c0o1;
      kxxMzzFromfcNEQ[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
          drho               = ((f_TNE+f_BSW)+(f_BSE+f_TNW)+(f_BNE+f_TSW)+(f_TSE+f_BNW)) +
                                ((f_NE+f_SW)+(f_TE+f_BW)+(f_SE+f_NW)+(f_BE+f_TW)+(f_BN+f_TS)+(f_TN+f_BS)) +
                                ((f_E-f_W) + (f_N-f_S) + (f_T-f_B)) + f_ZERO;
          rho                = drho + c1o1;
          vx1                = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W) / rho;
          vx2                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S) / rho;
          vx3                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B) / rho;
          kxyFromfcNEQ[k]    = -c3o1 *(f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE-(vx1*vx2));
          kyzFromfcNEQ[k]    = -c3o1 *(f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW-(vx2*vx3));
          kxzFromfcNEQ[k]    = -c3o1 *(f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE-(vx1*vx3));
          kxxMyyFromfcNEQ[k] = -c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));        //all dP00+dM00 minus all d0P0+d0M0 (no combinations of xy left)
          kxxMzzFromfcNEQ[k] = -c3o2 * (f_SW+f_W+f_NW-f_BS-f_TS-f_B-f_T-f_BN-f_TN+f_SE+f_E+f_NE-(vx1*vx1-vx3*vx3));        //all dP00+dM00 minus all d00P+d00M (no combinations of xz left)
      }
   }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc3rdMomentsIncompSP27(  real* CUMbbb,
                                                        real* CUMabc,
                                                        real* CUMbac,
                                                        real* CUMbca,
                                                        real* CUMcba,
                                                        real* CUMacb,
                                                        real* CUMcab,
                                                        unsigned int* bcMatD,
                                                        unsigned int* neighborX,
                                                        unsigned int* neighborY,
                                                        unsigned int* neighborZ,
                                                        real* DDStart,
                                                        unsigned long long numberOfLBnodes,
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

    if(k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
        {
            Distributions27 D;
            if (EvenOrOdd==true)
            {
                D.f[dP00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dM00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00M * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dMPM * numberOfLBnodes];
            }
            else
            {
                D.f[dM00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dP00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00M * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dMPM * numberOfLBnodes];
            }

            ////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kw   = neighborX[k];
            unsigned int ks   = neighborY[k];
            unsigned int kb   = neighborZ[k];
            unsigned int ksw  = neighborY[kw];
            unsigned int kbw  = neighborZ[kw];
            unsigned int kbs  = neighborZ[ks];
            unsigned int kbsw = neighborZ[ksw];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real mfcbb = (D.f[dP00])[k  ];
            real mfabb = (D.f[dM00])[kw ];
            real mfbcb = (D.f[d0P0])[k  ];
            real mfbab = (D.f[d0M0])[ks ];
            real mfbbc = (D.f[d00P])[k  ];
            real mfbba = (D.f[d00M])[kb ];
            real mfccb = (D.f[dPP0])[k  ];
            real mfaab = (D.f[dMM0])[ksw];
            real mfcab = (D.f[dPM0])[ks ];
            real mfacb = (D.f[dMP0])[kw ];
            real mfcbc = (D.f[dP0P])[k  ];
            real mfaba = (D.f[dM0M])[kbw];
            real mfcba = (D.f[dP0M])[kb ];
            real mfabc = (D.f[dM0P])[kw ];
            real mfbcc = (D.f[d0PP])[k  ];
            real mfbaa = (D.f[d0MM])[kbs];
            real mfbca = (D.f[d0PM])[kb ];
            real mfbac = (D.f[d0MP])[ks ];
            real mfbbb = (D.f[d000])[k  ];
            real mfccc = (D.f[dPPP])[k  ];
            real mfaac = (D.f[dMMP])[ksw];
            real mfcac = (D.f[dPMP])[ks ];
            real mfacc = (D.f[dMPP])[kw ];
            real mfcca = (D.f[dPPM])[kb ];
            real mfaaa = (D.f[dMMM])[kbsw];
            real mfcaa = (D.f[dPMM])[kbs];
            real mfaca = (D.f[dMPM])[kbw];
            ////////////////////////////////////////////////////////////////////////////////////
            real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
                             (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                               (mfcbb-mfabb));
            real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
                             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                               (mfbcb-mfbab));
            real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
                             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                               (mfbbc-mfbba));
            ////////////////////////////////////////////////////////////////////////////////////
            real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
                                   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
                                   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);
            ////////////////////////////////////////////////////////////////////////////////////
            real m0, m1, m2;    
            real vx2;
            real vy2;
            real vz2;
            vx2=vvx*vvx;
            vy2=vvy*vvy;
            vz2=vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            m2    = mfaaa    + mfaac;
            m1    = mfaac    - mfaaa;
            m0    = m2        + mfaab;
            mfaaa = m0;
            m0   += c1o36 * oMdrho;    
            mfaab = m1 -        m0 * vvz;
            mfaac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfabc;
            m1    = mfabc  - mfaba;
            m0    = m2        + mfabb;
            mfaba = m0;
            m0   += c1o9 * oMdrho;
            mfabb = m1 -        m0 * vvz;
            mfabc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfacc;
            m1    = mfacc  - mfaca;
            m0    = m2        + mfacb;
            mfaca = m0;
            m0   += c1o36 * oMdrho;
            mfacb = m1 -        m0 * vvz;
            mfacc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbac;
            m1    = mfbac    - mfbaa;
            m0    = m2        + mfbab;
            mfbaa = m0;
            m0   += c1o9 * oMdrho;
            mfbab = m1 -        m0 * vvz;
            mfbac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbba  + mfbbc;
            m1    = mfbbc  - mfbba;
            m0    = m2        + mfbbb;
            mfbba = m0;
            m0   += c4o9 * oMdrho;
            mfbbb = m1 -        m0 * vvz;
            mfbbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbca  + mfbcc;
            m1    = mfbcc  - mfbca;
            m0    = m2        + mfbcb;
            mfbca = m0;
            m0   += c1o9 * oMdrho;
            mfbcb = m1 -        m0 * vvz;
            mfbcc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcac;
            m1    = mfcac    - mfcaa;
            m0    = m2        + mfcab;
            mfcaa = m0;
            m0   += c1o36 * oMdrho;
            mfcab = m1 -        m0 * vvz;
            mfcac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcba  + mfcbc;
            m1    = mfcbc  - mfcba;
            m0    = m2        + mfcbb;
            mfcba = m0;
            m0   += c1o9 * oMdrho;
            mfcbb = m1 -        m0 * vvz;
            mfcbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcca  + mfccc;
            m1    = mfccc  - mfcca;
            m0    = m2        + mfccb;
            mfcca = m0;
            m0   += c1o36 * oMdrho;
            mfccb = m1 -        m0 * vvz;
            mfccc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m2    = mfaaa    + mfaca;
            m1    = mfaca    - mfaaa;
            m0    = m2        + mfaba;
            mfaaa = m0;
            m0   += c1o6 * oMdrho;
            mfaba = m1 -        m0 * vvy;
            mfaca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab  + mfacb;
            m1    = mfacb  - mfaab;
            m0    = m2        + mfabb;
            mfaab = m0;
            mfabb = m1 -        m0 * vvy;
            mfacb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac  + mfacc;
            m1    = mfacc  - mfaac;
            m0    = m2        + mfabc;
            mfaac = m0;
            m0   += c1o18 * oMdrho;
            mfabc = m1 -        m0 * vvy;
            mfacc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbca;
            m1    = mfbca    - mfbaa;
            m0    = m2        + mfbba;
            mfbaa = m0;
            m0   += c2o3 * oMdrho;
            mfbba = m1 -        m0 * vvy;
            mfbca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbab  + mfbcb;
            m1    = mfbcb  - mfbab;
            m0    = m2        + mfbbb;
            mfbab = m0;
            mfbbb = m1 -        m0 * vvy;
            mfbcb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbac  + mfbcc;
            m1    = mfbcc  - mfbac;
            m0    = m2        + mfbbc;
            mfbac = m0;
            m0   += c2o9 * oMdrho;
            mfbbc = m1 -        m0 * vvy;
            mfbcc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcca;
            m1    = mfcca    - mfcaa;
            m0    = m2        + mfcba;
            mfcaa = m0;
            m0   += c1o6 * oMdrho;
            mfcba = m1 -        m0 * vvy;
            mfcca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcab  + mfccb;
            m1    = mfccb  - mfcab;
            m0    = m2        + mfcbb;
            mfcab = m0;
            mfcbb = m1 -        m0 * vvy;
            mfccb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcac  + mfccc;
            m1    = mfccc  - mfcac;
            m0    = m2        + mfcbc;
            mfcac = m0;
            m0   += c1o18 * oMdrho;
            mfcbc = m1 -        m0 * vvy;
            mfccc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m2    = mfaaa    + mfcaa;
            m1    = mfcaa    - mfaaa;
            m0    = m2        + mfbaa;
            mfaaa = m0;
            m0   += c1o1* oMdrho;
            mfbaa = m1 -        m0 * vvx;
            mfcaa = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfcba;
            m1    = mfcba  - mfaba;
            m0    = m2        + mfbba;
            mfaba = m0;
            mfbba = m1 -        m0 * vvx;
            mfcba = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfcca;
            m1    = mfcca  - mfaca;
            m0    = m2        + mfbca;
            mfaca = m0;
            m0   += c1o3 * oMdrho;
            mfbca = m1 -        m0 * vvx;
            mfcca = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab    + mfcab;
            m1    = mfcab    - mfaab;
            m0    = m2        + mfbab;
            mfaab = m0;
            mfbab = m1 -        m0 * vvx;
            mfcab = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabb  + mfcbb;
            m1    = mfcbb  - mfabb;
            m0    = m2        + mfbbb;
            mfabb = m0;
            mfbbb = m1 -        m0 * vvx;
            mfcbb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacb  + mfccb;
            m1    = mfccb  - mfacb;
            m0    = m2        + mfbcb;
            mfacb = m0;
            mfbcb = m1 -        m0 * vvx;
            mfccb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac    + mfcac;
            m1    = mfcac    - mfaac;
            m0    = m2        + mfbac;
            mfaac = m0;
            m0   += c1o3 * oMdrho;
            mfbac = m1 -        m0 * vvx;
            mfcac = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabc  + mfcbc;
            m1    = mfcbc  - mfabc;
            m0    = m2        + mfbbc;
            mfabc = m0;
            mfbbc = m1 -        m0 * vvx;
            mfcbc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacc  + mfccc;
            m1    = mfccc  - mfacc;
            m0    = m2        + mfbcc;
            mfacc = m0;
            m0   += c1o9 * oMdrho;
            mfbcc = m1 -        m0 * vvx;
            mfccc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            //3.
            CUMbbb[k] = mfbbb;
            CUMabc[k] = mfabc;
            CUMbac[k] = mfbac;
            CUMbca[k] = mfbca;
            CUMcba[k] = mfcba;
            CUMacb[k] = mfacb;
            CUMcab[k] = mfcab;
            ////////////////////////////////////////////////////////////////////////////////////
        }                                                                                                                    
    }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc3rdMomentsCompSP27(real* CUMbbb,
                                                    real* CUMabc,
                                                    real* CUMbac,
                                                    real* CUMbca,
                                                    real* CUMcba,
                                                    real* CUMacb,
                                                    real* CUMcab,
                                                    unsigned int* bcMatD,
                                                    unsigned int* neighborX,
                                                    unsigned int* neighborY,
                                                    unsigned int* neighborZ,
                                                    real* DDStart,
                                                    unsigned long long numberOfLBnodes,
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

    if(k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
        {
            Distributions27 D;
            if (EvenOrOdd==true)
            {
                D.f[dP00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dM00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00M * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dMPM * numberOfLBnodes];
            }
            else
            {
                D.f[dM00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dP00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00M * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dMPM * numberOfLBnodes];
            }

            ////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kw   = neighborX[k];
            unsigned int ks   = neighborY[k];
            unsigned int kb   = neighborZ[k];
            unsigned int ksw  = neighborY[kw];
            unsigned int kbw  = neighborZ[kw];
            unsigned int kbs  = neighborZ[ks];
            unsigned int kbsw = neighborZ[ksw];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real mfcbb = (D.f[dP00])[k  ];
            real mfabb = (D.f[dM00])[kw ];
            real mfbcb = (D.f[d0P0])[k  ];
            real mfbab = (D.f[d0M0])[ks ];
            real mfbbc = (D.f[d00P])[k  ];
            real mfbba = (D.f[d00M])[kb ];
            real mfccb = (D.f[dPP0])[k  ];
            real mfaab = (D.f[dMM0])[ksw];
            real mfcab = (D.f[dPM0])[ks ];
            real mfacb = (D.f[dMP0])[kw ];
            real mfcbc = (D.f[dP0P])[k  ];
            real mfaba = (D.f[dM0M])[kbw];
            real mfcba = (D.f[dP0M])[kb ];
            real mfabc = (D.f[dM0P])[kw ];
            real mfbcc = (D.f[d0PP])[k  ];
            real mfbaa = (D.f[d0MM])[kbs];
            real mfbca = (D.f[d0PM])[kb ];
            real mfbac = (D.f[d0MP])[ks ];
            real mfbbb = (D.f[d000])[k  ];
            real mfccc = (D.f[dPPP])[k  ];
            real mfaac = (D.f[dMMP])[ksw];
            real mfcac = (D.f[dPMP])[ks ];
            real mfacc = (D.f[dMPP])[kw ];
            real mfcca = (D.f[dPPM])[kb ];
            real mfaaa = (D.f[dMMM])[kbsw];
            real mfcaa = (D.f[dPMM])[kbs];
            real mfaca = (D.f[dMPM])[kbw];
            ////////////////////////////////////////////////////////////////////////////////////
            real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
                            (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
                            ((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb;

            real rho = c1o1+drho;
            ////////////////////////////////////////////////////////////////////////////////////
            real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
                             (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                               (mfcbb-mfabb)) / rho;
            real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
                             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                               (mfbcb-mfbab)) / rho;
            real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
                             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                               (mfbbc-mfbba)) / rho;
            ////////////////////////////////////////////////////////////////////////////////////
            real oMdrho = c1o1; // comp special
            ////////////////////////////////////////////////////////////////////////////////////
            real m0, m1, m2;    
            real vx2;
            real vy2;
            real vz2;
            vx2=vvx*vvx;
            vy2=vvy*vvy;
            vz2=vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            m2    = mfaaa    + mfaac;
            m1    = mfaac    - mfaaa;
            m0    = m2        + mfaab;
            mfaaa = m0;
            m0   += c1o36 * oMdrho;    
            mfaab = m1 -        m0 * vvz;
            mfaac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfabc;
            m1    = mfabc  - mfaba;
            m0    = m2        + mfabb;
            mfaba = m0;
            m0   += c1o9 * oMdrho;
            mfabb = m1 -        m0 * vvz;
            mfabc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfacc;
            m1    = mfacc  - mfaca;
            m0    = m2        + mfacb;
            mfaca = m0;
            m0   += c1o36 * oMdrho;
            mfacb = m1 -        m0 * vvz;
            mfacc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbac;
            m1    = mfbac    - mfbaa;
            m0    = m2        + mfbab;
            mfbaa = m0;
            m0   += c1o9 * oMdrho;
            mfbab = m1 -        m0 * vvz;
            mfbac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbba  + mfbbc;
            m1    = mfbbc  - mfbba;
            m0    = m2        + mfbbb;
            mfbba = m0;
            m0   += c4o9 * oMdrho;
            mfbbb = m1 -        m0 * vvz;
            mfbbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbca  + mfbcc;
            m1    = mfbcc  - mfbca;
            m0    = m2        + mfbcb;
            mfbca = m0;
            m0   += c1o9 * oMdrho;
            mfbcb = m1 -        m0 * vvz;
            mfbcc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcac;
            m1    = mfcac    - mfcaa;
            m0    = m2        + mfcab;
            mfcaa = m0;
            m0   += c1o36 * oMdrho;
            mfcab = m1 -        m0 * vvz;
            mfcac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcba  + mfcbc;
            m1    = mfcbc  - mfcba;
            m0    = m2        + mfcbb;
            mfcba = m0;
            m0   += c1o9 * oMdrho;
            mfcbb = m1 -        m0 * vvz;
            mfcbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcca  + mfccc;
            m1    = mfccc  - mfcca;
            m0    = m2        + mfccb;
            mfcca = m0;
            m0   += c1o36 * oMdrho;
            mfccb = m1 -        m0 * vvz;
            mfccc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m2    = mfaaa    + mfaca;
            m1    = mfaca    - mfaaa;
            m0    = m2        + mfaba;
            mfaaa = m0;
            m0   += c1o6 * oMdrho;
            mfaba = m1 -        m0 * vvy;
            mfaca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab  + mfacb;
            m1    = mfacb  - mfaab;
            m0    = m2        + mfabb;
            mfaab = m0;
            mfabb = m1 -        m0 * vvy;
            mfacb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac  + mfacc;
            m1    = mfacc  - mfaac;
            m0    = m2        + mfabc;
            mfaac = m0;
            m0   += c1o18 * oMdrho;
            mfabc = m1 -        m0 * vvy;
            mfacc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbca;
            m1    = mfbca    - mfbaa;
            m0    = m2        + mfbba;
            mfbaa = m0;
            m0   += c2o3 * oMdrho;
            mfbba = m1 -        m0 * vvy;
            mfbca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbab  + mfbcb;
            m1    = mfbcb  - mfbab;
            m0    = m2        + mfbbb;
            mfbab = m0;
            mfbbb = m1 -        m0 * vvy;
            mfbcb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbac  + mfbcc;
            m1    = mfbcc  - mfbac;
            m0    = m2        + mfbbc;
            mfbac = m0;
            m0   += c2o9 * oMdrho;
            mfbbc = m1 -        m0 * vvy;
            mfbcc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcca;
            m1    = mfcca    - mfcaa;
            m0    = m2        + mfcba;
            mfcaa = m0;
            m0   += c1o6 * oMdrho;
            mfcba = m1 -        m0 * vvy;
            mfcca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcab  + mfccb;
            m1    = mfccb  - mfcab;
            m0    = m2        + mfcbb;
            mfcab = m0;
            mfcbb = m1 -        m0 * vvy;
            mfccb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcac  + mfccc;
            m1    = mfccc  - mfcac;
            m0    = m2        + mfcbc;
            mfcac = m0;
            m0   += c1o18 * oMdrho;
            mfcbc = m1 -        m0 * vvy;
            mfccc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m2    = mfaaa    + mfcaa;
            m1    = mfcaa    - mfaaa;
            m0    = m2        + mfbaa;
            mfaaa = m0;
            m0   += c1o1* oMdrho;
            mfbaa = m1 -        m0 * vvx;
            mfcaa = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfcba;
            m1    = mfcba  - mfaba;
            m0    = m2        + mfbba;
            mfaba = m0;
            mfbba = m1 -        m0 * vvx;
            mfcba = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfcca;
            m1    = mfcca  - mfaca;
            m0    = m2        + mfbca;
            mfaca = m0;
            m0   += c1o3 * oMdrho;
            mfbca = m1 -        m0 * vvx;
            mfcca = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab    + mfcab;
            m1    = mfcab    - mfaab;
            m0    = m2        + mfbab;
            mfaab = m0;
            mfbab = m1 -        m0 * vvx;
            mfcab = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabb  + mfcbb;
            m1    = mfcbb  - mfabb;
            m0    = m2        + mfbbb;
            mfabb = m0;
            mfbbb = m1 -        m0 * vvx;
            mfcbb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacb  + mfccb;
            m1    = mfccb  - mfacb;
            m0    = m2        + mfbcb;
            mfacb = m0;
            mfbcb = m1 -        m0 * vvx;
            mfccb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac    + mfcac;
            m1    = mfcac    - mfaac;
            m0    = m2        + mfbac;
            mfaac = m0;
            m0   += c1o3 * oMdrho;
            mfbac = m1 -        m0 * vvx;
            mfcac = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabc  + mfcbc;
            m1    = mfcbc  - mfabc;
            m0    = m2        + mfbbc;
            mfabc = m0;
            mfbbc = m1 -        m0 * vvx;
            mfcbc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacc  + mfccc;
            m1    = mfccc  - mfacc;
            m0    = m2        + mfbcc;
            mfacc = m0;
            m0   += c1o9 * oMdrho;
            mfbcc = m1 -        m0 * vvx;
            mfccc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            //3.
            CUMbbb[k] = mfbbb;
            CUMabc[k] = mfabc;
            CUMbac[k] = mfbac;
            CUMbca[k] = mfbca;
            CUMcba[k] = mfcba;
            CUMacb[k] = mfacb;
            CUMcab[k] = mfcab;
            ////////////////////////////////////////////////////////////////////////////////////
        }                                                                                                                    
    }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcHigherMomentsIncompSP27(   real* CUMcbb,
                                                            real* CUMbcb,
                                                            real* CUMbbc,
                                                            real* CUMcca,
                                                            real* CUMcac,
                                                            real* CUMacc,
                                                            real* CUMbcc,
                                                            real* CUMcbc,
                                                            real* CUMccb,
                                                            real* CUMccc,
                                                            unsigned int* bcMatD,
                                                            unsigned int* neighborX,
                                                            unsigned int* neighborY,
                                                            unsigned int* neighborZ,
                                                            real* DDStart,
                                                            unsigned long long numberOfLBnodes,
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

    if(k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
        {
            Distributions27 D;
            if (EvenOrOdd==true)
            {
                D.f[dP00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dM00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00M * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dMPM * numberOfLBnodes];
            }
            else
            {
                D.f[dM00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dP00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00M * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dMPM * numberOfLBnodes];
            }

            ////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kw   = neighborX[k];
            unsigned int ks   = neighborY[k];
            unsigned int kb   = neighborZ[k];
            unsigned int ksw  = neighborY[kw];
            unsigned int kbw  = neighborZ[kw];
            unsigned int kbs  = neighborZ[ks];
            unsigned int kbsw = neighborZ[ksw];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real mfcbb = (D.f[dP00])[k  ];
            real mfabb = (D.f[dM00])[kw ];
            real mfbcb = (D.f[d0P0])[k  ];
            real mfbab = (D.f[d0M0])[ks ];
            real mfbbc = (D.f[d00P])[k  ];
            real mfbba = (D.f[d00M])[kb ];
            real mfccb = (D.f[dPP0])[k  ];
            real mfaab = (D.f[dMM0])[ksw];
            real mfcab = (D.f[dPM0])[ks ];
            real mfacb = (D.f[dMP0])[kw ];
            real mfcbc = (D.f[dP0P])[k  ];
            real mfaba = (D.f[dM0M])[kbw];
            real mfcba = (D.f[dP0M])[kb ];
            real mfabc = (D.f[dM0P])[kw ];
            real mfbcc = (D.f[d0PP])[k  ];
            real mfbaa = (D.f[d0MM])[kbs];
            real mfbca = (D.f[d0PM])[kb ];
            real mfbac = (D.f[d0MP])[ks ];
            real mfbbb = (D.f[d000])[k  ];
            real mfccc = (D.f[dPPP])[k  ];
            real mfaac = (D.f[dMMP])[ksw];
            real mfcac = (D.f[dPMP])[ks ];
            real mfacc = (D.f[dMPP])[kw ];
            real mfcca = (D.f[dPPM])[kb ];
            real mfaaa = (D.f[dMMM])[kbsw];
            real mfcaa = (D.f[dPMM])[kbs];
            real mfaca = (D.f[dMPM])[kbw];
            ////////////////////////////////////////////////////////////////////////////////////
            real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
                             (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                               (mfcbb-mfabb));
            real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
                             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                               (mfbcb-mfbab));
            real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
                             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                               (mfbbc-mfbba));
            ////////////////////////////////////////////////////////////////////////////////////
            real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
                                   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
                                   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);
            ////////////////////////////////////////////////////////////////////////////////////
            real m0, m1, m2;    
            real vx2;
            real vy2;
            real vz2;
            vx2=vvx*vvx;
            vy2=vvy*vvy;
            vz2=vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            m2    = mfaaa    + mfaac;
            m1    = mfaac    - mfaaa;
            m0    = m2        + mfaab;
            mfaaa = m0;
            m0   += c1o36 * oMdrho;    
            mfaab = m1 -        m0 * vvz;
            mfaac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfabc;
            m1    = mfabc  - mfaba;
            m0    = m2        + mfabb;
            mfaba = m0;
            m0   += c1o9 * oMdrho;
            mfabb = m1 -        m0 * vvz;
            mfabc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfacc;
            m1    = mfacc  - mfaca;
            m0    = m2        + mfacb;
            mfaca = m0;
            m0   += c1o36 * oMdrho;
            mfacb = m1 -        m0 * vvz;
            mfacc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbac;
            m1    = mfbac    - mfbaa;
            m0    = m2        + mfbab;
            mfbaa = m0;
            m0   += c1o9 * oMdrho;
            mfbab = m1 -        m0 * vvz;
            mfbac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbba  + mfbbc;
            m1    = mfbbc  - mfbba;
            m0    = m2        + mfbbb;
            mfbba = m0;
            m0   += c4o9 * oMdrho;
            mfbbb = m1 -        m0 * vvz;
            mfbbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbca  + mfbcc;
            m1    = mfbcc  - mfbca;
            m0    = m2        + mfbcb;
            mfbca = m0;
            m0   += c1o9 * oMdrho;
            mfbcb = m1 -        m0 * vvz;
            mfbcc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcac;
            m1    = mfcac    - mfcaa;
            m0    = m2        + mfcab;
            mfcaa = m0;
            m0   += c1o36 * oMdrho;
            mfcab = m1 -        m0 * vvz;
            mfcac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcba  + mfcbc;
            m1    = mfcbc  - mfcba;
            m0    = m2        + mfcbb;
            mfcba = m0;
            m0   += c1o9 * oMdrho;
            mfcbb = m1 -        m0 * vvz;
            mfcbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcca  + mfccc;
            m1    = mfccc  - mfcca;
            m0    = m2        + mfccb;
            mfcca = m0;
            m0   += c1o36 * oMdrho;
            mfccb = m1 -        m0 * vvz;
            mfccc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m2    = mfaaa    + mfaca;
            m1    = mfaca    - mfaaa;
            m0    = m2        + mfaba;
            mfaaa = m0;
            m0   += c1o6 * oMdrho;
            mfaba = m1 -        m0 * vvy;
            mfaca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab  + mfacb;
            m1    = mfacb  - mfaab;
            m0    = m2        + mfabb;
            mfaab = m0;
            mfabb = m1 -        m0 * vvy;
            mfacb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac  + mfacc;
            m1    = mfacc  - mfaac;
            m0    = m2        + mfabc;
            mfaac = m0;
            m0   += c1o18 * oMdrho;
            mfabc = m1 -        m0 * vvy;
            mfacc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbca;
            m1    = mfbca    - mfbaa;
            m0    = m2        + mfbba;
            mfbaa = m0;
            m0   += c2o3 * oMdrho;
            mfbba = m1 -        m0 * vvy;
            mfbca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbab  + mfbcb;
            m1    = mfbcb  - mfbab;
            m0    = m2        + mfbbb;
            mfbab = m0;
            mfbbb = m1 -        m0 * vvy;
            mfbcb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbac  + mfbcc;
            m1    = mfbcc  - mfbac;
            m0    = m2        + mfbbc;
            mfbac = m0;
            m0   += c2o9 * oMdrho;
            mfbbc = m1 -        m0 * vvy;
            mfbcc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcca;
            m1    = mfcca    - mfcaa;
            m0    = m2        + mfcba;
            mfcaa = m0;
            m0   += c1o6 * oMdrho;
            mfcba = m1 -        m0 * vvy;
            mfcca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcab  + mfccb;
            m1    = mfccb  - mfcab;
            m0    = m2        + mfcbb;
            mfcab = m0;
            mfcbb = m1 -        m0 * vvy;
            mfccb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcac  + mfccc;
            m1    = mfccc  - mfcac;
            m0    = m2        + mfcbc;
            mfcac = m0;
            m0   += c1o18 * oMdrho;
            mfcbc = m1 -        m0 * vvy;
            mfccc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m2    = mfaaa    + mfcaa;
            m1    = mfcaa    - mfaaa;
            m0    = m2        + mfbaa;
            mfaaa = m0;
            m0   += c1o1* oMdrho;
            mfbaa = m1 -        m0 * vvx;
            mfcaa = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfcba;
            m1    = mfcba  - mfaba;
            m0    = m2        + mfbba;
            mfaba = m0;
            mfbba = m1 -        m0 * vvx;
            mfcba = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfcca;
            m1    = mfcca  - mfaca;
            m0    = m2        + mfbca;
            mfaca = m0;
            m0   += c1o3 * oMdrho;
            mfbca = m1 -        m0 * vvx;
            mfcca = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab    + mfcab;
            m1    = mfcab    - mfaab;
            m0    = m2        + mfbab;
            mfaab = m0;
            mfbab = m1 -        m0 * vvx;
            mfcab = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabb  + mfcbb;
            m1    = mfcbb  - mfabb;
            m0    = m2        + mfbbb;
            mfabb = m0;
            mfbbb = m1 -        m0 * vvx;
            mfcbb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacb  + mfccb;
            m1    = mfccb  - mfacb;
            m0    = m2        + mfbcb;
            mfacb = m0;
            mfbcb = m1 -        m0 * vvx;
            mfccb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac    + mfcac;
            m1    = mfcac    - mfaac;
            m0    = m2        + mfbac;
            mfaac = m0;
            m0   += c1o3 * oMdrho;
            mfbac = m1 -        m0 * vvx;
            mfcac = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabc  + mfcbc;
            m1    = mfcbc  - mfabc;
            m0    = m2        + mfbbc;
            mfabc = m0;
            mfbbc = m1 -        m0 * vvx;
            mfcbc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacc  + mfccc;
            m1    = mfccc  - mfacc;
            m0    = m2        + mfbcc;
            mfacc = m0;
            m0   += c1o9 * oMdrho;
            mfbcc = m1 -        m0 * vvx;
            mfccc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            //Cum 4.
            CUMcbb[k]      = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab); 
            CUMbcb[k]      = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb); 
            CUMbbc[k]      = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); 

            CUMcca[k]      = mfcca - ((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
            CUMcac[k]      = mfcac - ((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
            CUMacc[k]      = mfacc - ((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);

            //Cum 5.
            CUMbcc[k]      = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
            CUMcbc[k]      = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
            CUMccb[k]      = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

            //Cum 6.
            CUMccc[k]      = mfccc  +((-c4o1 *  mfbbb * mfbbb  
                            -           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                            -    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                            -     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
                            +(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                            +     c2o1 * (mfcaa * mfaca * mfaac)
                            + c16o1 *  mfbba * mfbab * mfabb)
                            -    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
                            -    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
                            +(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                            +           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;
            ////////////////////////////////////////////////////////////////////////////////////
        }                                                                                                                    
    }
}































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcHigherMomentsCompSP27( real* CUMcbb,
                                                        real* CUMbcb,
                                                        real* CUMbbc,
                                                        real* CUMcca,
                                                        real* CUMcac,
                                                        real* CUMacc,
                                                        real* CUMbcc,
                                                        real* CUMcbc,
                                                        real* CUMccb,
                                                        real* CUMccc,
                                                        unsigned int* bcMatD,
                                                        unsigned int* neighborX,
                                                        unsigned int* neighborY,
                                                        unsigned int* neighborZ,
                                                        real* DDStart,
                                                        unsigned long long numberOfLBnodes,
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

    if(k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
        {
            Distributions27 D;
            if (EvenOrOdd==true)
            {
                D.f[dP00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dM00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00M * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dMPM * numberOfLBnodes];
            }
            else
            {
                D.f[dM00] = &DDStart[dP00 * numberOfLBnodes];
                D.f[dP00] = &DDStart[dM00 * numberOfLBnodes];
                D.f[d0M0] = &DDStart[d0P0 * numberOfLBnodes];
                D.f[d0P0] = &DDStart[d0M0 * numberOfLBnodes];
                D.f[d00M] = &DDStart[d00P * numberOfLBnodes];
                D.f[d00P] = &DDStart[d00M * numberOfLBnodes];
                D.f[dMM0] = &DDStart[dPP0 * numberOfLBnodes];
                D.f[dPP0] = &DDStart[dMM0 * numberOfLBnodes];
                D.f[dMP0] = &DDStart[dPM0 * numberOfLBnodes];
                D.f[dPM0] = &DDStart[dMP0 * numberOfLBnodes];
                D.f[dM0M] = &DDStart[dP0P * numberOfLBnodes];
                D.f[dP0P] = &DDStart[dM0M * numberOfLBnodes];
                D.f[dM0P] = &DDStart[dP0M * numberOfLBnodes];
                D.f[dP0M] = &DDStart[dM0P * numberOfLBnodes];
                D.f[d0MM] = &DDStart[d0PP * numberOfLBnodes];
                D.f[d0PP] = &DDStart[d0MM * numberOfLBnodes];
                D.f[d0MP] = &DDStart[d0PM * numberOfLBnodes];
                D.f[d0PM] = &DDStart[d0MP * numberOfLBnodes];
                D.f[d000] = &DDStart[d000 * numberOfLBnodes];
                D.f[dMMM] = &DDStart[dPPP * numberOfLBnodes];
                D.f[dPPM] = &DDStart[dMMP * numberOfLBnodes];
                D.f[dMPM] = &DDStart[dPMP * numberOfLBnodes];
                D.f[dPMM] = &DDStart[dMPP * numberOfLBnodes];
                D.f[dMMP] = &DDStart[dPPM * numberOfLBnodes];
                D.f[dPPP] = &DDStart[dMMM * numberOfLBnodes];
                D.f[dMPP] = &DDStart[dPMM * numberOfLBnodes];
                D.f[dPMP] = &DDStart[dMPM * numberOfLBnodes];
            }

            ////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kw   = neighborX[k];
            unsigned int ks   = neighborY[k];
            unsigned int kb   = neighborZ[k];
            unsigned int ksw  = neighborY[kw];
            unsigned int kbw  = neighborZ[kw];
            unsigned int kbs  = neighborZ[ks];
            unsigned int kbsw = neighborZ[ksw];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real mfcbb = (D.f[dP00])[k  ];
            real mfabb = (D.f[dM00])[kw ];
            real mfbcb = (D.f[d0P0])[k  ];
            real mfbab = (D.f[d0M0])[ks ];
            real mfbbc = (D.f[d00P])[k  ];
            real mfbba = (D.f[d00M])[kb ];
            real mfccb = (D.f[dPP0])[k  ];
            real mfaab = (D.f[dMM0])[ksw];
            real mfcab = (D.f[dPM0])[ks ];
            real mfacb = (D.f[dMP0])[kw ];
            real mfcbc = (D.f[dP0P])[k  ];
            real mfaba = (D.f[dM0M])[kbw];
            real mfcba = (D.f[dP0M])[kb ];
            real mfabc = (D.f[dM0P])[kw ];
            real mfbcc = (D.f[d0PP])[k  ];
            real mfbaa = (D.f[d0MM])[kbs];
            real mfbca = (D.f[d0PM])[kb ];
            real mfbac = (D.f[d0MP])[ks ];
            real mfbbb = (D.f[d000])[k  ];
            real mfccc = (D.f[dPPP])[k  ];
            real mfaac = (D.f[dMMP])[ksw];
            real mfcac = (D.f[dPMP])[ks ];
            real mfacc = (D.f[dMPP])[kw ];
            real mfcca = (D.f[dPPM])[kb ];
            real mfaaa = (D.f[dMMM])[kbsw];
            real mfcaa = (D.f[dPMM])[kbs];
            real mfaca = (D.f[dMPM])[kbw];
            ////////////////////////////////////////////////////////////////////////////////////
            real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
                            (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
                            ((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb;

            real rho = c1o1+drho;
            ////////////////////////////////////////////////////////////////////////////////////
            real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
                             (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
                               (mfcbb-mfabb)) / rho;
            real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
                             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
                               (mfbcb-mfbab)) / rho;
            real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
                             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
                               (mfbbc-mfbba)) / rho;
            ////////////////////////////////////////////////////////////////////////////////////
            real oMdrho = c1o1; // comp special
            ////////////////////////////////////////////////////////////////////////////////////
            real m0, m1, m2;    
            real vx2;
            real vy2;
            real vz2;
            vx2=vvx*vvx;
            vy2=vvy*vvy;
            vz2=vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            m2    = mfaaa    + mfaac;
            m1    = mfaac    - mfaaa;
            m0    = m2        + mfaab;
            mfaaa = m0;
            m0   += c1o36 * oMdrho;    
            mfaab = m1 -        m0 * vvz;
            mfaac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfabc;
            m1    = mfabc  - mfaba;
            m0    = m2        + mfabb;
            mfaba = m0;
            m0   += c1o9 * oMdrho;
            mfabb = m1 -        m0 * vvz;
            mfabc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfacc;
            m1    = mfacc  - mfaca;
            m0    = m2        + mfacb;
            mfaca = m0;
            m0   += c1o36 * oMdrho;
            mfacb = m1 -        m0 * vvz;
            mfacc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbac;
            m1    = mfbac    - mfbaa;
            m0    = m2        + mfbab;
            mfbaa = m0;
            m0   += c1o9 * oMdrho;
            mfbab = m1 -        m0 * vvz;
            mfbac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbba  + mfbbc;
            m1    = mfbbc  - mfbba;
            m0    = m2        + mfbbb;
            mfbba = m0;
            m0   += c4o9 * oMdrho;
            mfbbb = m1 -        m0 * vvz;
            mfbbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbca  + mfbcc;
            m1    = mfbcc  - mfbca;
            m0    = m2        + mfbcb;
            mfbca = m0;
            m0   += c1o9 * oMdrho;
            mfbcb = m1 -        m0 * vvz;
            mfbcc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcac;
            m1    = mfcac    - mfcaa;
            m0    = m2        + mfcab;
            mfcaa = m0;
            m0   += c1o36 * oMdrho;
            mfcab = m1 -        m0 * vvz;
            mfcac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcba  + mfcbc;
            m1    = mfcbc  - mfcba;
            m0    = m2        + mfcbb;
            mfcba = m0;
            m0   += c1o9 * oMdrho;
            mfcbb = m1 -        m0 * vvz;
            mfcbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcca  + mfccc;
            m1    = mfccc  - mfcca;
            m0    = m2        + mfccb;
            mfcca = m0;
            m0   += c1o36 * oMdrho;
            mfccb = m1 -        m0 * vvz;
            mfccc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m2    = mfaaa    + mfaca;
            m1    = mfaca    - mfaaa;
            m0    = m2        + mfaba;
            mfaaa = m0;
            m0   += c1o6 * oMdrho;
            mfaba = m1 -        m0 * vvy;
            mfaca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab  + mfacb;
            m1    = mfacb  - mfaab;
            m0    = m2        + mfabb;
            mfaab = m0;
            mfabb = m1 -        m0 * vvy;
            mfacb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac  + mfacc;
            m1    = mfacc  - mfaac;
            m0    = m2        + mfabc;
            mfaac = m0;
            m0   += c1o18 * oMdrho;
            mfabc = m1 -        m0 * vvy;
            mfacc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbaa    + mfbca;
            m1    = mfbca    - mfbaa;
            m0    = m2        + mfbba;
            mfbaa = m0;
            m0   += c2o3 * oMdrho;
            mfbba = m1 -        m0 * vvy;
            mfbca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbab  + mfbcb;
            m1    = mfbcb  - mfbab;
            m0    = m2        + mfbbb;
            mfbab = m0;
            mfbbb = m1 -        m0 * vvy;
            mfbcb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfbac  + mfbcc;
            m1    = mfbcc  - mfbac;
            m0    = m2        + mfbbc;
            mfbac = m0;
            m0   += c2o9 * oMdrho;
            mfbbc = m1 -        m0 * vvy;
            mfbcc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcaa    + mfcca;
            m1    = mfcca    - mfcaa;
            m0    = m2        + mfcba;
            mfcaa = m0;
            m0   += c1o6 * oMdrho;
            mfcba = m1 -        m0 * vvy;
            mfcca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcab  + mfccb;
            m1    = mfccb  - mfcab;
            m0    = m2        + mfcbb;
            mfcab = m0;
            mfcbb = m1 -        m0 * vvy;
            mfccb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfcac  + mfccc;
            m1    = mfccc  - mfcac;
            m0    = m2        + mfcbc;
            mfcac = m0;
            m0   += c1o18 * oMdrho;
            mfcbc = m1 -        m0 * vvy;
            mfccc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m2    = mfaaa    + mfcaa;
            m1    = mfcaa    - mfaaa;
            m0    = m2        + mfbaa;
            mfaaa = m0;
            m0   += c1o1* oMdrho;
            mfbaa = m1 -        m0 * vvx;
            mfcaa = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaba  + mfcba;
            m1    = mfcba  - mfaba;
            m0    = m2        + mfbba;
            mfaba = m0;
            mfbba = m1 -        m0 * vvx;
            mfcba = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaca  + mfcca;
            m1    = mfcca  - mfaca;
            m0    = m2        + mfbca;
            mfaca = m0;
            m0   += c1o3 * oMdrho;
            mfbca = m1 -        m0 * vvx;
            mfcca = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaab    + mfcab;
            m1    = mfcab    - mfaab;
            m0    = m2        + mfbab;
            mfaab = m0;
            mfbab = m1 -        m0 * vvx;
            mfcab = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabb  + mfcbb;
            m1    = mfcbb  - mfabb;
            m0    = m2        + mfbbb;
            mfabb = m0;
            mfbbb = m1 -        m0 * vvx;
            mfcbb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacb  + mfccb;
            m1    = mfccb  - mfacb;
            m0    = m2        + mfbcb;
            mfacb = m0;
            mfbcb = m1 -        m0 * vvx;
            mfccb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfaac    + mfcac;
            m1    = mfcac    - mfaac;
            m0    = m2        + mfbac;
            mfaac = m0;
            m0   += c1o3 * oMdrho;
            mfbac = m1 -        m0 * vvx;
            mfcac = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfabc  + mfcbc;
            m1    = mfcbc  - mfabc;
            m0    = m2        + mfbbc;
            mfabc = m0;
            mfbbc = m1 -        m0 * vvx;
            mfcbc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2    = mfacc  + mfccc;
            m1    = mfccc  - mfacc;
            m0    = m2        + mfbcc;
            mfacc = m0;
            m0   += c1o9 * oMdrho;
            mfbcc = m1 -        m0 * vvx;
            mfccc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////

            real OxxPyyPzz = c1o1;
            real omega = c1o1 / (c3o1*0.001 + c1o2);
            real d00M = (c4o1 * omega * OxxPyyPzz * (c9o1 * omega - c16o1) - c4o1 * omega * omega - c2o1 * OxxPyyPzz * OxxPyyPzz * (c2o1 + c9o1 * omega * (omega - c2o1))) /
                (c3o1 * (omega - OxxPyyPzz) * (OxxPyyPzz * (c2o1 + c3o1 * omega) - c8o1 * omega));

            CUMbcc[k] = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)*(c1o1 + rho*c6o1*d00M / (c2o1 + c3o1 * d00M))) / rho;
            CUMcbc[k] = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)*(c1o1 + rho*c6o1*d00M / (c2o1 + c3o1 * d00M))) / rho;
            CUMccb[k] = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)*(c1o1 + rho*c6o1*d00M / (c2o1 + c3o1 * d00M))) / rho;

            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            //central moments to cumulants
            //4.
            CUMcbb[k]      = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;    
            CUMbcb[k]      = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho; 
            CUMbbc[k]      = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho; 
                     
            CUMcca[k]      = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
            CUMcac[k]      = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
            CUMacc[k]      = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

            //5.
            //CUMbcc[k]      = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
            //CUMcbc[k]      = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
            //CUMccb[k]      = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
            
            //6.
            CUMccc[k]      = mfccc + ((-c4o1 *  mfbbb * mfbbb  
                            -           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                            -    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                            -     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
                            +(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                            +     c2o1 * (mfcaa * mfaca * mfaac)
                            + c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
                            -    c1o3 * (mfacc + mfcac + mfcca) /rho 
                            -    c1o9 * (mfcaa + mfaca + mfaac) /rho 
                            +(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
                            +           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
                            + c1o27*((drho * drho - drho)/(rho*rho)));
            ////////////////////////////////////////////////////////////////////////////////////
        }                                                                                                                    
    }
}

void Calc2ndMomentsIncompSP27(real* kxyFromfcNEQ, real* kyzFromfcNEQ, real* kxzFromfcNEQ, real* kxxMyyFromfcNEQ,
                              real* kxxMzzFromfcNEQ, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                              unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                              real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc2ndMomentsIncompSP27<<<grid.grid, grid.threads>>>(kxyFromfcNEQ, kyzFromfcNEQ, kxzFromfcNEQ, kxxMyyFromfcNEQ,
                                                            kxxMzzFromfcNEQ, geoD, neighborX, neighborY, neighborZ,
                                                            numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalc2ndMomentsIncompSP27 execution failed");
}

void Calc2ndMomentsCompSP27(real* kxyFromfcNEQ, real* kyzFromfcNEQ, real* kxzFromfcNEQ, real* kxxMyyFromfcNEQ,
                            real* kxxMzzFromfcNEQ, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                            unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                            real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc2ndMomentsCompSP27<<<grid.grid, grid.threads>>>(kxyFromfcNEQ, kyzFromfcNEQ, kxzFromfcNEQ, kxxMyyFromfcNEQ,
                                                          kxxMzzFromfcNEQ, geoD, neighborX, neighborY, neighborZ,
                                                          numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalc2ndMomentsCompSP27 execution failed");
}

void Calc3rdMomentsIncompSP27(real* CUMbbb, real* CUMabc, real* CUMbac, real* CUMbca, real* CUMcba, real* CUMacb,
                              real* CUMcab, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY,
                              unsigned int* neighborZ, unsigned long long numberOfLBnodes, unsigned int numberOfThreads,
                              real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc3rdMomentsIncompSP27<<<grid.grid, grid.threads>>>(CUMbbb, CUMabc, CUMbac, CUMbca, CUMcba, CUMacb, CUMcab, geoD,
                                                            neighborX, neighborY, neighborZ, DD, numberOfLBnodes,
                                                            isEvenTimestep);
    getLastCudaError("LBCalc3rdMomentsIncompSP27 execution failed");
}

void Calc3rdMomentsCompSP27(real* CUMbbb, real* CUMabc, real* CUMbac, real* CUMbca, real* CUMcba, real* CUMacb, real* CUMcab,
                            unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                            unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalc3rdMomentsCompSP27<<<grid.grid, grid.threads>>>(CUMbbb, CUMabc, CUMbac, CUMbca, CUMcba, CUMacb, CUMcab, geoD,
                                                          neighborX, neighborY, neighborZ, DD, numberOfLBnodes,
                                                          isEvenTimestep);
    getLastCudaError("LBCalc3rdMomentsCompSP27 execution failed");
}

void CalcHigherMomentsIncompSP27(real* CUMcbb, real* CUMbcb, real* CUMbbc, real* CUMcca, real* CUMcac, real* CUMacc,
                                 real* CUMbcc, real* CUMcbc, real* CUMccb, real* CUMccc, unsigned int* geoD,
                                 unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                                 unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD,
                                 bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcHigherMomentsIncompSP27<<<grid.grid, grid.threads>>>(CUMcbb, CUMbcb, CUMbbc, CUMcca, CUMcac, CUMacc, CUMbcc,
                                                               CUMcbc, CUMccb, CUMccc, geoD, neighborX, neighborY, neighborZ,
                                                               DD, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBCalcHigherMomentsIncompSP27 execution failed");
}

void CalcHigherMomentsCompSP27(real* CUMcbb, real* CUMbcb, real* CUMbbc, real* CUMcca, real* CUMcac, real* CUMacc,
                               real* CUMbcc, real* CUMcbc, real* CUMccb, real* CUMccc, unsigned int* geoD,
                               unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                               unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD,
                               bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBCalcHigherMomentsCompSP27<<<grid.grid, grid.threads>>>(CUMcbb, CUMbcb, CUMbbc, CUMcca, CUMcac, CUMacc, CUMbcc, CUMcbc,
                                                             CUMccb, CUMccc, geoD, neighborX, neighborY, neighborZ, DD,
                                                             numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBCalcHigherMomentsCompSP27 execution failed");
}
