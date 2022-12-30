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
//! \file CalcConc27.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
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
         distAD.f[DIR_P00   ] = &distributionsAD[DIR_P00   *numberOfLBnodes];
         distAD.f[DIR_M00   ] = &distributionsAD[DIR_M00   *numberOfLBnodes];
         distAD.f[DIR_0P0   ] = &distributionsAD[DIR_0P0   *numberOfLBnodes];
         distAD.f[DIR_0M0   ] = &distributionsAD[DIR_0M0   *numberOfLBnodes];
         distAD.f[DIR_00P   ] = &distributionsAD[DIR_00P   *numberOfLBnodes];
         distAD.f[DIR_00M   ] = &distributionsAD[DIR_00M   *numberOfLBnodes];
         distAD.f[DIR_PP0  ] = &distributionsAD[DIR_PP0  *numberOfLBnodes];
         distAD.f[DIR_MM0  ] = &distributionsAD[DIR_MM0  *numberOfLBnodes];
         distAD.f[DIR_PM0  ] = &distributionsAD[DIR_PM0  *numberOfLBnodes];
         distAD.f[DIR_MP0  ] = &distributionsAD[DIR_MP0  *numberOfLBnodes];
         distAD.f[DIR_P0P  ] = &distributionsAD[DIR_P0P  *numberOfLBnodes];
         distAD.f[DIR_M0M  ] = &distributionsAD[DIR_M0M  *numberOfLBnodes];
         distAD.f[DIR_P0M  ] = &distributionsAD[DIR_P0M  *numberOfLBnodes];
         distAD.f[DIR_M0P  ] = &distributionsAD[DIR_M0P  *numberOfLBnodes];
         distAD.f[DIR_0PP  ] = &distributionsAD[DIR_0PP  *numberOfLBnodes];
         distAD.f[DIR_0MM  ] = &distributionsAD[DIR_0MM  *numberOfLBnodes];
         distAD.f[DIR_0PM  ] = &distributionsAD[DIR_0PM  *numberOfLBnodes];
         distAD.f[DIR_0MP  ] = &distributionsAD[DIR_0MP  *numberOfLBnodes];
         distAD.f[DIR_000] = &distributionsAD[DIR_000*numberOfLBnodes];
         distAD.f[DIR_PPP ] = &distributionsAD[DIR_PPP *numberOfLBnodes];
         distAD.f[DIR_MMP ] = &distributionsAD[DIR_MMP *numberOfLBnodes];
         distAD.f[DIR_PMP ] = &distributionsAD[DIR_PMP *numberOfLBnodes];
         distAD.f[DIR_MPP ] = &distributionsAD[DIR_MPP *numberOfLBnodes];
         distAD.f[DIR_PPM ] = &distributionsAD[DIR_PPM *numberOfLBnodes];
         distAD.f[DIR_MMM ] = &distributionsAD[DIR_MMM *numberOfLBnodes];
         distAD.f[DIR_PMM ] = &distributionsAD[DIR_PMM *numberOfLBnodes];
         distAD.f[DIR_MPM ] = &distributionsAD[DIR_MPM *numberOfLBnodes];
      }
      else
      {
         distAD.f[DIR_M00   ] = &distributionsAD[DIR_P00   *numberOfLBnodes];
         distAD.f[DIR_P00   ] = &distributionsAD[DIR_M00   *numberOfLBnodes];
         distAD.f[DIR_0M0   ] = &distributionsAD[DIR_0P0   *numberOfLBnodes];
         distAD.f[DIR_0P0   ] = &distributionsAD[DIR_0M0   *numberOfLBnodes];
         distAD.f[DIR_00M   ] = &distributionsAD[DIR_00P   *numberOfLBnodes];
         distAD.f[DIR_00P   ] = &distributionsAD[DIR_00M   *numberOfLBnodes];
         distAD.f[DIR_MM0  ] = &distributionsAD[DIR_PP0  *numberOfLBnodes];
         distAD.f[DIR_PP0  ] = &distributionsAD[DIR_MM0  *numberOfLBnodes];
         distAD.f[DIR_MP0  ] = &distributionsAD[DIR_PM0  *numberOfLBnodes];
         distAD.f[DIR_PM0  ] = &distributionsAD[DIR_MP0  *numberOfLBnodes];
         distAD.f[DIR_M0M  ] = &distributionsAD[DIR_P0P  *numberOfLBnodes];
         distAD.f[DIR_P0P  ] = &distributionsAD[DIR_M0M  *numberOfLBnodes];
         distAD.f[DIR_M0P  ] = &distributionsAD[DIR_P0M  *numberOfLBnodes];
         distAD.f[DIR_P0M  ] = &distributionsAD[DIR_M0P  *numberOfLBnodes];
         distAD.f[DIR_0MM  ] = &distributionsAD[DIR_0PP  *numberOfLBnodes];
         distAD.f[DIR_0PP  ] = &distributionsAD[DIR_0MM  *numberOfLBnodes];
         distAD.f[DIR_0MP  ] = &distributionsAD[DIR_0PM  *numberOfLBnodes];
         distAD.f[DIR_0PM  ] = &distributionsAD[DIR_0MP  *numberOfLBnodes];
         distAD.f[DIR_000] = &distributionsAD[DIR_000*numberOfLBnodes];
         distAD.f[DIR_PPP ] = &distributionsAD[DIR_MMM *numberOfLBnodes];
         distAD.f[DIR_MMP ] = &distributionsAD[DIR_PPM *numberOfLBnodes];
         distAD.f[DIR_PMP ] = &distributionsAD[DIR_MPM *numberOfLBnodes];
         distAD.f[DIR_MPP ] = &distributionsAD[DIR_PMM *numberOfLBnodes];
         distAD.f[DIR_PPM ] = &distributionsAD[DIR_MMP *numberOfLBnodes];
         distAD.f[DIR_MMM ] = &distributionsAD[DIR_PPP *numberOfLBnodes];
         distAD.f[DIR_PMM ] = &distributionsAD[DIR_MPP *numberOfLBnodes];
         distAD.f[DIR_MPM ] = &distributionsAD[DIR_PMP *numberOfLBnodes];
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
	  real mfcbb = (distAD.f[DIR_P00   ])[ke  ];
	  real mfabb = (distAD.f[DIR_M00   ])[kw  ];
	  real mfbcb = (distAD.f[DIR_0P0   ])[kn  ];
	  real mfbab = (distAD.f[DIR_0M0   ])[ks  ];
	  real mfbbc = (distAD.f[DIR_00P   ])[kt  ];
	  real mfbba = (distAD.f[DIR_00M   ])[kb  ];
	  real mfccb = (distAD.f[DIR_PP0  ])[kne ];
	  real mfaab = (distAD.f[DIR_MM0  ])[ksw ];
	  real mfcab = (distAD.f[DIR_PM0  ])[kse ];
	  real mfacb = (distAD.f[DIR_MP0  ])[knw ];
	  real mfcbc = (distAD.f[DIR_P0P  ])[kte ];
	  real mfaba = (distAD.f[DIR_M0M  ])[kbw ];
	  real mfcba = (distAD.f[DIR_P0M  ])[kbe ];
	  real mfabc = (distAD.f[DIR_M0P  ])[ktw ];
	  real mfbcc = (distAD.f[DIR_0PP  ])[ktn ];
	  real mfbaa = (distAD.f[DIR_0MM  ])[kbs ];
	  real mfbca = (distAD.f[DIR_0PM  ])[kbn ];
	  real mfbac = (distAD.f[DIR_0MP  ])[kts ];
	  real mfbbb = (distAD.f[DIR_000])[k   ];
	  real mfccc = (distAD.f[DIR_PPP ])[ktne];
	  real mfaac = (distAD.f[DIR_MMP ])[ktsw];
	  real mfcac = (distAD.f[DIR_PMP ])[ktse];
	  real mfacc = (distAD.f[DIR_MPP ])[ktnw];
	  real mfcca = (distAD.f[DIR_PPM ])[kbne];
	  real mfaaa = (distAD.f[DIR_MMM ])[kbsw];
	  real mfcaa = (distAD.f[DIR_PMM ])[kbse];
	  real mfaca = (distAD.f[DIR_MPM ])[kbnw];
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
//       D27.f[DIR_P00   ] = &DD27[DIR_P00   *size_Mat];
//       D27.f[DIR_M00   ] = &DD27[DIR_M00   *size_Mat];
//       D27.f[DIR_0P0   ] = &DD27[DIR_0P0   *size_Mat];
//       D27.f[DIR_0M0   ] = &DD27[DIR_0M0   *size_Mat];
//       D27.f[DIR_00P   ] = &DD27[DIR_00P   *size_Mat];
//       D27.f[DIR_00M   ] = &DD27[DIR_00M   *size_Mat];
//       D27.f[DIR_PP0  ] = &DD27[DIR_PP0  *size_Mat];
//       D27.f[DIR_MM0  ] = &DD27[DIR_MM0  *size_Mat];
//       D27.f[DIR_PM0  ] = &DD27[DIR_PM0  *size_Mat];
//       D27.f[DIR_MP0  ] = &DD27[DIR_MP0  *size_Mat];
//       D27.f[DIR_P0P  ] = &DD27[DIR_P0P  *size_Mat];
//       D27.f[DIR_M0M  ] = &DD27[DIR_M0M  *size_Mat];
//       D27.f[DIR_P0M  ] = &DD27[DIR_P0M  *size_Mat];
//       D27.f[DIR_M0P  ] = &DD27[DIR_M0P  *size_Mat];
//       D27.f[DIR_0PP  ] = &DD27[DIR_0PP  *size_Mat];
//       D27.f[DIR_0MM  ] = &DD27[DIR_0MM  *size_Mat];
//       D27.f[DIR_0PM  ] = &DD27[DIR_0PM  *size_Mat];
//       D27.f[DIR_0MP  ] = &DD27[DIR_0MP  *size_Mat];
//       D27.f[DIR_000] = &DD27[DIR_000*size_Mat];
//       D27.f[DIR_PPP ] = &DD27[DIR_PPP *size_Mat];
//       D27.f[DIR_MMP ] = &DD27[DIR_MMP *size_Mat];
//       D27.f[DIR_PMP ] = &DD27[DIR_PMP *size_Mat];
//       D27.f[DIR_MPP ] = &DD27[DIR_MPP *size_Mat];
//       D27.f[DIR_PPM ] = &DD27[DIR_PPM *size_Mat];
//       D27.f[DIR_MMM ] = &DD27[DIR_MMM *size_Mat];
//       D27.f[DIR_PMM ] = &DD27[DIR_PMM *size_Mat];
//       D27.f[DIR_MPM ] = &DD27[DIR_MPM *size_Mat];
//    }
//    else
//    {
//       D27.f[DIR_M00   ] = &DD27[DIR_P00   *size_Mat];
//       D27.f[DIR_P00   ] = &DD27[DIR_M00   *size_Mat];
//       D27.f[DIR_0M0   ] = &DD27[DIR_0P0   *size_Mat];
//       D27.f[DIR_0P0   ] = &DD27[DIR_0M0   *size_Mat];
//       D27.f[DIR_00M   ] = &DD27[DIR_00P   *size_Mat];
//       D27.f[DIR_00P   ] = &DD27[DIR_00M   *size_Mat];
//       D27.f[DIR_MM0  ] = &DD27[DIR_PP0  *size_Mat];
//       D27.f[DIR_PP0  ] = &DD27[DIR_MM0  *size_Mat];
//       D27.f[DIR_MP0  ] = &DD27[DIR_PM0  *size_Mat];
//       D27.f[DIR_PM0  ] = &DD27[DIR_MP0  *size_Mat];
//       D27.f[DIR_M0M  ] = &DD27[DIR_P0P  *size_Mat];
//       D27.f[DIR_P0P  ] = &DD27[DIR_M0M  *size_Mat];
//       D27.f[DIR_M0P  ] = &DD27[DIR_P0M  *size_Mat];
//       D27.f[DIR_P0M  ] = &DD27[DIR_M0P  *size_Mat];
//       D27.f[DIR_0MM  ] = &DD27[DIR_0PP  *size_Mat];
//       D27.f[DIR_0PP  ] = &DD27[DIR_0MM  *size_Mat];
//       D27.f[DIR_0MP  ] = &DD27[DIR_0PM  *size_Mat];
//       D27.f[DIR_0PM  ] = &DD27[DIR_0MP  *size_Mat];
//       D27.f[DIR_000] = &DD27[DIR_000*size_Mat];
//       D27.f[DIR_MMM ] = &DD27[DIR_PPP *size_Mat];
//       D27.f[DIR_PPM ] = &DD27[DIR_MMP *size_Mat];
//       D27.f[DIR_MPM ] = &DD27[DIR_PMP *size_Mat];
//       D27.f[DIR_PMM ] = &DD27[DIR_MPP *size_Mat];
//       D27.f[DIR_MMP ] = &DD27[DIR_PPM *size_Mat];
//       D27.f[DIR_PPP ] = &DD27[DIR_MMM *size_Mat];
//       D27.f[DIR_MPP ] = &DD27[DIR_PMM *size_Mat];
//       D27.f[DIR_PMP ] = &DD27[DIR_MPM *size_Mat];
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
//          Conc[k]    =   (D27.f[DIR_P00   ])[ke  ]+ (D27.f[DIR_M00   ])[kw  ]+ 
//                         (D27.f[DIR_0P0   ])[kn  ]+ (D27.f[DIR_0M0   ])[ks  ]+
//                         (D27.f[DIR_00P   ])[kt  ]+ (D27.f[DIR_00M   ])[kb  ]+
//                         (D27.f[DIR_PP0  ])[kne ]+ (D27.f[DIR_MM0  ])[ksw ]+
//                         (D27.f[DIR_PM0  ])[kse ]+ (D27.f[DIR_MP0  ])[knw ]+
//                         (D27.f[DIR_P0P  ])[kte ]+ (D27.f[DIR_M0M  ])[kbw ]+
//                         (D27.f[DIR_P0M  ])[kbe ]+ (D27.f[DIR_M0P  ])[ktw ]+
//                         (D27.f[DIR_0PP  ])[ktn ]+ (D27.f[DIR_0MM  ])[kbs ]+
//                         (D27.f[DIR_0PM  ])[kbn ]+ (D27.f[DIR_0MP  ])[kts ]+
//                         (D27.f[DIR_000])[kzero]+ 
//                         (D27.f[DIR_PPP ])[ktne]+ (D27.f[DIR_MMP ])[ktsw]+
//                         (D27.f[DIR_PMP ])[ktse]+ (D27.f[DIR_MPP ])[ktnw]+
//                         (D27.f[DIR_PPM ])[kbne]+ (D27.f[DIR_MMM ])[kbsw]+
//                         (D27.f[DIR_PMM ])[kbse]+ (D27.f[DIR_MPM ])[kbnw];
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
      D27.f[DIR_P00   ] = &DD27[DIR_P00   *numberOfLBnodes];
      D27.f[DIR_M00   ] = &DD27[DIR_M00   *numberOfLBnodes];
      D27.f[DIR_0P0   ] = &DD27[DIR_0P0   *numberOfLBnodes];
      D27.f[DIR_0M0   ] = &DD27[DIR_0M0   *numberOfLBnodes];
      D27.f[DIR_00P   ] = &DD27[DIR_00P   *numberOfLBnodes];
      D27.f[DIR_00M   ] = &DD27[DIR_00M   *numberOfLBnodes];
      D27.f[DIR_PP0  ] = &DD27[DIR_PP0  *numberOfLBnodes];
      D27.f[DIR_MM0  ] = &DD27[DIR_MM0  *numberOfLBnodes];
      D27.f[DIR_PM0  ] = &DD27[DIR_PM0  *numberOfLBnodes];
      D27.f[DIR_MP0  ] = &DD27[DIR_MP0  *numberOfLBnodes];
      D27.f[DIR_P0P  ] = &DD27[DIR_P0P  *numberOfLBnodes];
      D27.f[DIR_M0M  ] = &DD27[DIR_M0M  *numberOfLBnodes];
      D27.f[DIR_P0M  ] = &DD27[DIR_P0M  *numberOfLBnodes];
      D27.f[DIR_M0P  ] = &DD27[DIR_M0P  *numberOfLBnodes];
      D27.f[DIR_0PP  ] = &DD27[DIR_0PP  *numberOfLBnodes];
      D27.f[DIR_0MM  ] = &DD27[DIR_0MM  *numberOfLBnodes];
      D27.f[DIR_0PM  ] = &DD27[DIR_0PM  *numberOfLBnodes];
      D27.f[DIR_0MP  ] = &DD27[DIR_0MP  *numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000*numberOfLBnodes];
      D27.f[DIR_PPP ] = &DD27[DIR_PPP *numberOfLBnodes];
      D27.f[DIR_MMP ] = &DD27[DIR_MMP *numberOfLBnodes];
      D27.f[DIR_PMP ] = &DD27[DIR_PMP *numberOfLBnodes];
      D27.f[DIR_MPP ] = &DD27[DIR_MPP *numberOfLBnodes];
      D27.f[DIR_PPM ] = &DD27[DIR_PPM *numberOfLBnodes];
      D27.f[DIR_MMM ] = &DD27[DIR_MMM *numberOfLBnodes];
      D27.f[DIR_PMM ] = &DD27[DIR_PMM *numberOfLBnodes];
      D27.f[DIR_MPM ] = &DD27[DIR_MPM *numberOfLBnodes];
   }
   else
   {
      D27.f[DIR_M00   ] = &DD27[DIR_P00   *numberOfLBnodes];
      D27.f[DIR_P00   ] = &DD27[DIR_M00   *numberOfLBnodes];
      D27.f[DIR_0M0   ] = &DD27[DIR_0P0   *numberOfLBnodes];
      D27.f[DIR_0P0   ] = &DD27[DIR_0M0   *numberOfLBnodes];
      D27.f[DIR_00M   ] = &DD27[DIR_00P   *numberOfLBnodes];
      D27.f[DIR_00P   ] = &DD27[DIR_00M   *numberOfLBnodes];
      D27.f[DIR_MM0  ] = &DD27[DIR_PP0  *numberOfLBnodes];
      D27.f[DIR_PP0  ] = &DD27[DIR_MM0  *numberOfLBnodes];
      D27.f[DIR_MP0  ] = &DD27[DIR_PM0  *numberOfLBnodes];
      D27.f[DIR_PM0  ] = &DD27[DIR_MP0  *numberOfLBnodes];
      D27.f[DIR_M0M  ] = &DD27[DIR_P0P  *numberOfLBnodes];
      D27.f[DIR_P0P  ] = &DD27[DIR_M0M  *numberOfLBnodes];
      D27.f[DIR_M0P  ] = &DD27[DIR_P0M  *numberOfLBnodes];
      D27.f[DIR_P0M  ] = &DD27[DIR_M0P  *numberOfLBnodes];
      D27.f[DIR_0MM  ] = &DD27[DIR_0PP  *numberOfLBnodes];
      D27.f[DIR_0PP  ] = &DD27[DIR_0MM  *numberOfLBnodes];
      D27.f[DIR_0MP  ] = &DD27[DIR_0PM  *numberOfLBnodes];
      D27.f[DIR_0PM  ] = &DD27[DIR_0MP  *numberOfLBnodes];
      D27.f[DIR_000] = &DD27[DIR_000*numberOfLBnodes];
      D27.f[DIR_MMM ] = &DD27[DIR_PPP *numberOfLBnodes];
      D27.f[DIR_PPM ] = &DD27[DIR_MMP *numberOfLBnodes];
      D27.f[DIR_MPM ] = &DD27[DIR_PMP *numberOfLBnodes];
      D27.f[DIR_PMM ] = &DD27[DIR_MPP *numberOfLBnodes];
      D27.f[DIR_MMP ] = &DD27[DIR_PPM *numberOfLBnodes];
      D27.f[DIR_PPP ] = &DD27[DIR_MMM *numberOfLBnodes];
      D27.f[DIR_MPP ] = &DD27[DIR_PMM *numberOfLBnodes];
      D27.f[DIR_PMP ] = &DD27[DIR_MPM *numberOfLBnodes];
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
         Conc[k]    =   (D27.f[DIR_P00   ])[ke  ]+ (D27.f[DIR_M00   ])[kw  ]+ 
                        (D27.f[DIR_0P0   ])[kn  ]+ (D27.f[DIR_0M0   ])[ks  ]+
                        (D27.f[DIR_00P   ])[kt  ]+ (D27.f[DIR_00M   ])[kb  ]+
                        (D27.f[DIR_PP0  ])[kne ]+ (D27.f[DIR_MM0  ])[ksw ]+
                        (D27.f[DIR_PM0  ])[kse ]+ (D27.f[DIR_MP0  ])[knw ]+
                        (D27.f[DIR_P0P  ])[kte ]+ (D27.f[DIR_M0M  ])[kbw ]+
                        (D27.f[DIR_P0M  ])[kbe ]+ (D27.f[DIR_M0P  ])[ktw ]+
                        (D27.f[DIR_0PP  ])[ktn ]+ (D27.f[DIR_0MM  ])[kbs ]+
                        (D27.f[DIR_0PM  ])[kbn ]+ (D27.f[DIR_0MP  ])[kts ]+
                        (D27.f[DIR_000])[kzero]+ 
                        (D27.f[DIR_PPP ])[ktne]+ (D27.f[DIR_MMP ])[ktsw]+
                        (D27.f[DIR_PMP ])[ktse]+ (D27.f[DIR_MPP ])[ktnw]+
                        (D27.f[DIR_PPM ])[kbne]+ (D27.f[DIR_MMM ])[kbsw]+
                        (D27.f[DIR_PMM ])[kbse]+ (D27.f[DIR_MPM ])[kbnw];
      }
   }   
}