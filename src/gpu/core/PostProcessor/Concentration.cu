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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \
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
#include "Concentration.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "LBM/LB.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void CalcConc27(
    real* concentration,
    uint* typeOfGridNode,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributionsAD,
    bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // global x-index
   const unsigned  y = blockIdx.x;   // global y-index
   const unsigned  z = blockIdx.y;   // global z-index

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////
   // run for all indices in size_Mat and fluid nodes
   if ((k < numberOfLBnodes) && (typeOfGridNode[k] == GEO_FLUID))
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 distAD;
      if (isEvenTimestep)
      {
         distAD.f[dP00] = &distributionsAD[dP00 * numberOfLBnodes];
         distAD.f[dM00] = &distributionsAD[dM00 * numberOfLBnodes];
         distAD.f[d0P0] = &distributionsAD[d0P0 * numberOfLBnodes];
         distAD.f[d0M0] = &distributionsAD[d0M0 * numberOfLBnodes];
         distAD.f[d00P] = &distributionsAD[d00P * numberOfLBnodes];
         distAD.f[d00M] = &distributionsAD[d00M * numberOfLBnodes];
         distAD.f[dPP0] = &distributionsAD[dPP0 * numberOfLBnodes];
         distAD.f[dMM0] = &distributionsAD[dMM0 * numberOfLBnodes];
         distAD.f[dPM0] = &distributionsAD[dPM0 * numberOfLBnodes];
         distAD.f[dMP0] = &distributionsAD[dMP0 * numberOfLBnodes];
         distAD.f[dP0P] = &distributionsAD[dP0P * numberOfLBnodes];
         distAD.f[dM0M] = &distributionsAD[dM0M * numberOfLBnodes];
         distAD.f[dP0M] = &distributionsAD[dP0M * numberOfLBnodes];
         distAD.f[dM0P] = &distributionsAD[dM0P * numberOfLBnodes];
         distAD.f[d0PP] = &distributionsAD[d0PP * numberOfLBnodes];
         distAD.f[d0MM] = &distributionsAD[d0MM * numberOfLBnodes];
         distAD.f[d0PM] = &distributionsAD[d0PM * numberOfLBnodes];
         distAD.f[d0MP] = &distributionsAD[d0MP * numberOfLBnodes];
         distAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
         distAD.f[dPPP] = &distributionsAD[dPPP * numberOfLBnodes];
         distAD.f[dMMP] = &distributionsAD[dMMP * numberOfLBnodes];
         distAD.f[dPMP] = &distributionsAD[dPMP * numberOfLBnodes];
         distAD.f[dMPP] = &distributionsAD[dMPP * numberOfLBnodes];
         distAD.f[dPPM] = &distributionsAD[dPPM * numberOfLBnodes];
         distAD.f[dMMM] = &distributionsAD[dMMM * numberOfLBnodes];
         distAD.f[dPMM] = &distributionsAD[dPMM * numberOfLBnodes];
         distAD.f[dMPM] = &distributionsAD[dMPM * numberOfLBnodes];
      }
      else
      {
         distAD.f[dM00] = &distributionsAD[dP00 * numberOfLBnodes];
         distAD.f[dP00] = &distributionsAD[dM00 * numberOfLBnodes];
         distAD.f[d0M0] = &distributionsAD[d0P0 * numberOfLBnodes];
         distAD.f[d0P0] = &distributionsAD[d0M0 * numberOfLBnodes];
         distAD.f[d00M] = &distributionsAD[d00P * numberOfLBnodes];
         distAD.f[d00P] = &distributionsAD[d00M * numberOfLBnodes];
         distAD.f[dMM0] = &distributionsAD[dPP0 * numberOfLBnodes];
         distAD.f[dPP0] = &distributionsAD[dMM0 * numberOfLBnodes];
         distAD.f[dMP0] = &distributionsAD[dPM0 * numberOfLBnodes];
         distAD.f[dPM0] = &distributionsAD[dMP0 * numberOfLBnodes];
         distAD.f[dM0M] = &distributionsAD[dP0P * numberOfLBnodes];
         distAD.f[dP0P] = &distributionsAD[dM0M * numberOfLBnodes];
         distAD.f[dM0P] = &distributionsAD[dP0M * numberOfLBnodes];
         distAD.f[dP0M] = &distributionsAD[dM0P * numberOfLBnodes];
         distAD.f[d0MM] = &distributionsAD[d0PP * numberOfLBnodes];
         distAD.f[d0PP] = &distributionsAD[d0MM * numberOfLBnodes];
         distAD.f[d0MP] = &distributionsAD[d0PM * numberOfLBnodes];
         distAD.f[d0PM] = &distributionsAD[d0MP * numberOfLBnodes];
         distAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
         distAD.f[dPPP] = &distributionsAD[dMMM * numberOfLBnodes];
         distAD.f[dMMP] = &distributionsAD[dPPM * numberOfLBnodes];
         distAD.f[dPMP] = &distributionsAD[dMPM * numberOfLBnodes];
         distAD.f[dMPP] = &distributionsAD[dPMM * numberOfLBnodes];
         distAD.f[dPPM] = &distributionsAD[dMMP * numberOfLBnodes];
         distAD.f[dMMM] = &distributionsAD[dPPP * numberOfLBnodes];
         distAD.f[dPMM] = &distributionsAD[dMPP * numberOfLBnodes];
         distAD.f[dMPM] = &distributionsAD[dPMP * numberOfLBnodes];
      }
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      uint ke   = k;
      uint kw   = neighborX[k];
      uint kn   = k;
      uint ks   = neighborY[k];
      uint kt   = k;
      uint kb   = neighborZ[k];
      uint ksw  = neighborY[kw];
      uint kne  = k;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = k;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = k;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = k;
      uint kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions
      //!
      real mfcbb = (distAD.f[dP00])[ke  ];
      real mfabb = (distAD.f[dM00])[kw  ];
      real mfbcb = (distAD.f[d0P0])[kn  ];
      real mfbab = (distAD.f[d0M0])[ks  ];
      real mfbbc = (distAD.f[d00P])[kt  ];
      real mfbba = (distAD.f[d00M])[kb  ];
      real mfccb = (distAD.f[dPP0])[kne ];
      real mfaab = (distAD.f[dMM0])[ksw ];
      real mfcab = (distAD.f[dPM0])[kse ];
      real mfacb = (distAD.f[dMP0])[knw ];
      real mfcbc = (distAD.f[dP0P])[kte ];
      real mfaba = (distAD.f[dM0M])[kbw ];
      real mfcba = (distAD.f[dP0M])[kbe ];
      real mfabc = (distAD.f[dM0P])[ktw ];
      real mfbcc = (distAD.f[d0PP])[ktn ];
      real mfbaa = (distAD.f[d0MM])[kbs ];
      real mfbca = (distAD.f[d0PM])[kbn ];
      real mfbac = (distAD.f[d0MP])[kts ];
      real mfbbb = (distAD.f[d000])[k   ];
      real mfccc = (distAD.f[dPPP])[ktne];
      real mfaac = (distAD.f[dMMP])[ktsw];
      real mfcac = (distAD.f[dPMP])[ktse];
      real mfacc = (distAD.f[dMPP])[ktnw];
      real mfcca = (distAD.f[dPPM])[kbne];
      real mfaaa = (distAD.f[dMMM])[kbsw];
      real mfcaa = (distAD.f[dPMM])[kbse];
      real mfaca = (distAD.f[dMPM])[kbnw];
      //////////////////////////////////////////////////////////////////////////
      //! - Calculate concentration using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
      //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
      //!
      concentration[k] =
       ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
          (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
           ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb;

   }
}



















////////////////////////////////////////////////////////////////////////////////
__global__ void CalcConc7( real* Conc,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned long long numberOfLBnodes,
                                          real* DD7,
                                          bool isEvenTimestep)
{
   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   } 
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
   }
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
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= k;
      unsigned int ke   = k;
      unsigned int kw   = neighborX[k];
      unsigned int kn   = k;
      unsigned int ks   = neighborY[k];
      unsigned int kt   = k;
      unsigned int kb   = neighborZ[k];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D7.f[1])[ke   ]+ (D7.f[2])[kw  ]+ 
                        (D7.f[3])[kn   ]+ (D7.f[4])[ks  ]+
                        (D7.f[5])[kt   ]+ (D7.f[6])[kb  ]+
                        (D7.f[0])[kzero];  
      }
   }
}

























// DEPRECATED (2022)
//////////////////////////////////////////////////////////////////////////////////
// __global__ void LBCalcMacThS27(real* Conc,
//                                           unsigned int* geoD,
//                                           unsigned int* neighborX,
//                                           unsigned int* neighborY,
//                                           unsigned int* neighborZ,
//                                           unsigned int size_Mat,
//                                           real* DD27,
//                                           bool isEvenTimestep)
// {
//    Distributions27 D27;
//    if (isEvenTimestep==true)
//    {
//       D27.f[dP00] = &DD27[dP00 * size_Mat];
//       D27.f[dM00] = &DD27[dM00 * size_Mat];
//       D27.f[d0P0] = &DD27[d0P0 * size_Mat];
//       D27.f[d0M0] = &DD27[d0M0 * size_Mat];
//       D27.f[d00P] = &DD27[d00P * size_Mat];
//       D27.f[d00M] = &DD27[d00M * size_Mat];
//       D27.f[dPP0] = &DD27[dPP0 * size_Mat];
//       D27.f[dMM0] = &DD27[dMM0 * size_Mat];
//       D27.f[dPM0] = &DD27[dPM0 * size_Mat];
//       D27.f[dMP0] = &DD27[dMP0 * size_Mat];
//       D27.f[dP0P] = &DD27[dP0P * size_Mat];
//       D27.f[dM0M] = &DD27[dM0M * size_Mat];
//       D27.f[dP0M] = &DD27[dP0M * size_Mat];
//       D27.f[dM0P] = &DD27[dM0P * size_Mat];
//       D27.f[d0PP] = &DD27[d0PP * size_Mat];
//       D27.f[d0MM] = &DD27[d0MM * size_Mat];
//       D27.f[d0PM] = &DD27[d0PM * size_Mat];
//       D27.f[d0MP] = &DD27[d0MP * size_Mat];
//       D27.f[d000] = &DD27[d000 * size_Mat];
//       D27.f[dPPP] = &DD27[dPPP * size_Mat];
//       D27.f[dMMP] = &DD27[dMMP * size_Mat];
//       D27.f[dPMP] = &DD27[dPMP * size_Mat];
//       D27.f[dMPP] = &DD27[dMPP * size_Mat];
//       D27.f[dPPM] = &DD27[dPPM * size_Mat];
//       D27.f[dMMM] = &DD27[dMMM * size_Mat];
//       D27.f[dPMM] = &DD27[dPMM * size_Mat];
//       D27.f[dMPM] = &DD27[dMPM * size_Mat];
//    }
//    else
//    {
//       D27.f[dM00] = &DD27[dP00 * size_Mat];
//       D27.f[dP00] = &DD27[dM00 * size_Mat];
//       D27.f[d0M0] = &DD27[d0P0 * size_Mat];
//       D27.f[d0P0] = &DD27[d0M0 * size_Mat];
//       D27.f[d00M] = &DD27[d00P * size_Mat];
//       D27.f[d00P] = &DD27[d00M * size_Mat];
//       D27.f[dMM0] = &DD27[dPP0 * size_Mat];
//       D27.f[dPP0] = &DD27[dMM0 * size_Mat];
//       D27.f[dMP0] = &DD27[dPM0 * size_Mat];
//       D27.f[dPM0] = &DD27[dMP0 * size_Mat];
//       D27.f[dM0M] = &DD27[dP0P * size_Mat];
//       D27.f[dP0P] = &DD27[dM0M * size_Mat];
//       D27.f[dM0P] = &DD27[dP0M * size_Mat];
//       D27.f[dP0M] = &DD27[dM0P * size_Mat];
//       D27.f[d0MM] = &DD27[d0PP * size_Mat];
//       D27.f[d0PP] = &DD27[d0MM * size_Mat];
//       D27.f[d0MP] = &DD27[d0PM * size_Mat];
//       D27.f[d0PM] = &DD27[d0MP * size_Mat];
//       D27.f[d000] = &DD27[d000 * size_Mat];
//       D27.f[dMMM] = &DD27[dPPP * size_Mat];
//       D27.f[dPPM] = &DD27[dMMP * size_Mat];
//       D27.f[dMPM] = &DD27[dPMP * size_Mat];
//       D27.f[dPMM] = &DD27[dMPP * size_Mat];
//       D27.f[dMMP] = &DD27[dPPM * size_Mat];
//       D27.f[dPPP] = &DD27[dMMM * size_Mat];
//       D27.f[dMPP] = &DD27[dPMM * size_Mat];
//       D27.f[dPMP] = &DD27[dMPM * size_Mat];
//    }
//    ////////////////////////////////////////////////////////////////////////////////
//    const unsigned  x = threadIdx.x;  // Globaler x-Index 
//    const unsigned  y = blockIdx.x;   // Globaler y-Index 
//    const unsigned  z = blockIdx.y;   // Globaler z-Index 

//    const unsigned nx = blockDim.x;
//    const unsigned ny = gridDim.x;

//    const unsigned k = nx*(ny*z + y) + x;
//    //////////////////////////////////////////////////////////////////////////

//    if(k<size_Mat)
//    {
//       //////////////////////////////////////////////////////////////////////////
//       //index
//       unsigned int kzero= k;
//       unsigned int ke   = k;
//       unsigned int kw   = neighborX[k];
//       unsigned int kn   = k;
//       unsigned int ks   = neighborY[k];
//       unsigned int kt   = k;
//       unsigned int kb   = neighborZ[k];
//       unsigned int ksw  = neighborY[kw];
//       unsigned int kne  = k;
//       unsigned int kse  = ks;
//       unsigned int knw  = kw;
//       unsigned int kbw  = neighborZ[kw];
//       unsigned int kte  = k;
//       unsigned int kbe  = kb;
//       unsigned int ktw  = kw;
//       unsigned int kbs  = neighborZ[ks];
//       unsigned int ktn  = k;
//       unsigned int kbn  = kb;
//       unsigned int kts  = ks;
//       unsigned int ktse = ks;
//       unsigned int kbnw = kbw;
//       unsigned int ktnw = kw;
//       unsigned int kbse = kbs;
//       unsigned int ktsw = ksw;
//       unsigned int kbne = kb;
//       unsigned int ktne = k;
//       unsigned int kbsw = neighborZ[ksw];
//       //////////////////////////////////////////////////////////////////////////
//       Conc[k] = c0o1;

//       if(geoD[k] == GEO_FLUID)
//       {
//          Conc[k]    =   (D27.f[dP00])[ke  ]+ (D27.f[dM00])[kw  ]+ 
//                         (D27.f[d0P0])[kn  ]+ (D27.f[d0M0])[ks  ]+
//                         (D27.f[d00P])[kt  ]+ (D27.f[d00M])[kb  ]+
//                         (D27.f[dPP0])[kne ]+ (D27.f[dMM0])[ksw ]+
//                         (D27.f[dPM0])[kse ]+ (D27.f[dMP0])[knw ]+
//                         (D27.f[dP0P])[kte ]+ (D27.f[dM0M])[kbw ]+
//                         (D27.f[dP0M])[kbe ]+ (D27.f[dM0P])[ktw ]+
//                         (D27.f[d0PP])[ktn ]+ (D27.f[d0MM])[kbs ]+
//                         (D27.f[d0PM])[kbn ]+ (D27.f[d0MP])[kts ]+
//                         (D27.f[d000])[kzero]+ 
//                         (D27.f[dPPP])[ktne]+ (D27.f[dMMP])[ktsw]+
//                         (D27.f[dPMP])[ktse]+ (D27.f[dMPP])[ktnw]+
//                         (D27.f[dPPM])[kbne]+ (D27.f[dMMM])[kbsw]+
//                         (D27.f[dPMM])[kbse]+ (D27.f[dMPM])[kbnw];
//       }
//    }   
// }



















////////////////////////////////////////////////////////////////////////////////
__global__ void GetPlaneConc7(real* Conc,
                                            int* kPC,
                                            unsigned int numberOfPointskPC,
                                            unsigned int* geoD,
                                            unsigned int* neighborX,
                                            unsigned int* neighborY,
                                            unsigned int* neighborZ,
                                            unsigned long long numberOfLBnodes,
                                            real* DD7,
                                            bool isEvenTimestep)
{
   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[1] = &DD7[1*numberOfLBnodes];
      D7.f[2] = &DD7[2*numberOfLBnodes];
      D7.f[3] = &DD7[3*numberOfLBnodes];
      D7.f[4] = &DD7[4*numberOfLBnodes];
      D7.f[5] = &DD7[5*numberOfLBnodes];
      D7.f[6] = &DD7[6*numberOfLBnodes];
   } 
   else
   {
      D7.f[0] = &DD7[0*numberOfLBnodes];
      D7.f[2] = &DD7[1*numberOfLBnodes];
      D7.f[1] = &DD7[2*numberOfLBnodes];
      D7.f[4] = &DD7[3*numberOfLBnodes];
      D7.f[3] = &DD7[4*numberOfLBnodes];
      D7.f[6] = &DD7[5*numberOfLBnodes];
      D7.f[5] = &DD7[6*numberOfLBnodes];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfPointskPC)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kPC[k];
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D7.f[1])[ke   ]+ (D7.f[2])[kw  ]+ 
                        (D7.f[3])[kn   ]+ (D7.f[4])[ks  ]+
                        (D7.f[5])[kt   ]+ (D7.f[6])[kb  ]+
                        (D7.f[0])[kzero];  
      }
   }
}






































////////////////////////////////////////////////////////////////////////////////
__global__ void GetPlaneConc27(real* Conc,
                                             int* kPC,
                                             unsigned int numberOfPointskPC,
                                             unsigned int* geoD,
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned long long numberOfLBnodes,
                                             real* DD27,
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
      D27.f[dMMM] = &DD27[dPPP * numberOfLBnodes];
      D27.f[dPPM] = &DD27[dMMP * numberOfLBnodes];
      D27.f[dMPM] = &DD27[dPMP * numberOfLBnodes];
      D27.f[dPMM] = &DD27[dMPP * numberOfLBnodes];
      D27.f[dMMP] = &DD27[dPPM * numberOfLBnodes];
      D27.f[dPPP] = &DD27[dMMM * numberOfLBnodes];
      D27.f[dMPP] = &DD27[dPMM * numberOfLBnodes];
      D27.f[dPMP] = &DD27[dMPM * numberOfLBnodes];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfPointskPC)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      unsigned int kzero= kPC[k];
      unsigned int ke   = kzero;
      unsigned int kw   = neighborX[kzero];
      unsigned int kn   = kzero;
      unsigned int ks   = neighborY[kzero];
      unsigned int kt   = kzero;
      unsigned int kb   = neighborZ[kzero];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = kzero;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = kzero;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = kzero;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = kzero;
      unsigned int kbsw = neighborZ[ksw];
      //////////////////////////////////////////////////////////////////////////
      Conc[k] = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
         Conc[k]    =   (D27.f[dP00])[ke  ]+ (D27.f[dM00])[kw  ]+ 
                        (D27.f[d0P0])[kn  ]+ (D27.f[d0M0])[ks  ]+
                        (D27.f[d00P])[kt  ]+ (D27.f[d00M])[kb  ]+
                        (D27.f[dPP0])[kne ]+ (D27.f[dMM0])[ksw ]+
                        (D27.f[dPM0])[kse ]+ (D27.f[dMP0])[knw ]+
                        (D27.f[dP0P])[kte ]+ (D27.f[dM0M])[kbw ]+
                        (D27.f[dP0M])[kbe ]+ (D27.f[dM0P])[ktw ]+
                        (D27.f[d0PP])[ktn ]+ (D27.f[d0MM])[kbs ]+
                        (D27.f[d0PM])[kbn ]+ (D27.f[d0MP])[kts ]+
                        (D27.f[d000])[kzero]+ 
                        (D27.f[dPPP])[ktne]+ (D27.f[dMMP])[ktsw]+
                        (D27.f[dPMP])[ktse]+ (D27.f[dMPP])[ktnw]+
                        (D27.f[dPPM])[kbne]+ (D27.f[dMMM])[kbsw]+
                        (D27.f[dPMM])[kbse]+ (D27.f[dMPM])[kbnw];
      }
   }   
}

void PlaneConcThS7(real* Conc, int* kPC, unsigned int numberOfPointskPC, unsigned int* geoD, unsigned int* neighborX,
                   unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                   unsigned int numberOfThreads, real* DD7, bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

   GetPlaneConc7<<<grid.grid, grid.threads>>>(Conc, kPC, numberOfPointskPC, geoD, neighborX, neighborY, neighborZ,
                                              numberOfLBnodes, DD7, isEvenTimestep);
   getLastCudaError("GetPlaneConc7 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void PlaneConcThS27(real* Conc, int* kPC, unsigned int numberOfPointskPC, unsigned int* geoD, unsigned int* neighborX,
                    unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                    unsigned int numberOfThreads, real* DD27, bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskPC);

   GetPlaneConc27<<<grid.grid, grid.threads>>>(Conc, kPC, numberOfPointskPC, geoD, neighborX, neighborY, neighborZ,
                                               numberOfLBnodes, DD27, isEvenTimestep);
   getLastCudaError("GetPlaneConc27 execution failed");
}

void CalcConcentration27(unsigned int numberOfThreads, real* Conc, unsigned int* geoD, unsigned int* neighborX,
                         unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes, real* DD27,
                         bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

   CalcConc27<<<grid.grid, grid.threads>>>(Conc, geoD, neighborX, neighborY, neighborZ, numberOfLBnodes, DD27,
                                           isEvenTimestep);
   getLastCudaError("CalcConc27 execution failed");
}

void CalcConcThS7(real* Conc, unsigned int* geoD, unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                  unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD7, bool isEvenTimestep)
{
   vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

   CalcConc7<<<grid.grid, grid.threads>>>(Conc, geoD, neighborX, neighborY, neighborZ, numberOfLBnodes, DD7, isEvenTimestep);
   getLastCudaError("CalcConc7 execution failed");
}
