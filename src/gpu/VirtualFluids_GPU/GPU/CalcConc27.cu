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
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void CalcConc27(
	real* concentration,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat,
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
   if ((k < size_Mat) && (typeOfGridNode[k] == GEO_FLUID))
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 distAD;
      if (isEvenTimestep)
      {
         distAD.f[E   ] = &distributionsAD[E   *size_Mat];
         distAD.f[W   ] = &distributionsAD[W   *size_Mat];
         distAD.f[N   ] = &distributionsAD[N   *size_Mat];
         distAD.f[S   ] = &distributionsAD[S   *size_Mat];
         distAD.f[T   ] = &distributionsAD[T   *size_Mat];
         distAD.f[B   ] = &distributionsAD[B   *size_Mat];
         distAD.f[NE  ] = &distributionsAD[NE  *size_Mat];
         distAD.f[SW  ] = &distributionsAD[SW  *size_Mat];
         distAD.f[SE  ] = &distributionsAD[SE  *size_Mat];
         distAD.f[NW  ] = &distributionsAD[NW  *size_Mat];
         distAD.f[TE  ] = &distributionsAD[TE  *size_Mat];
         distAD.f[BW  ] = &distributionsAD[BW  *size_Mat];
         distAD.f[BE  ] = &distributionsAD[BE  *size_Mat];
         distAD.f[TW  ] = &distributionsAD[TW  *size_Mat];
         distAD.f[TN  ] = &distributionsAD[TN  *size_Mat];
         distAD.f[BS  ] = &distributionsAD[BS  *size_Mat];
         distAD.f[BN  ] = &distributionsAD[BN  *size_Mat];
         distAD.f[TS  ] = &distributionsAD[TS  *size_Mat];
         distAD.f[dirREST] = &distributionsAD[dirREST*size_Mat];
         distAD.f[TNE ] = &distributionsAD[TNE *size_Mat];
         distAD.f[TSW ] = &distributionsAD[TSW *size_Mat];
         distAD.f[TSE ] = &distributionsAD[TSE *size_Mat];
         distAD.f[TNW ] = &distributionsAD[TNW *size_Mat];
         distAD.f[BNE ] = &distributionsAD[BNE *size_Mat];
         distAD.f[BSW ] = &distributionsAD[BSW *size_Mat];
         distAD.f[BSE ] = &distributionsAD[BSE *size_Mat];
         distAD.f[BNW ] = &distributionsAD[BNW *size_Mat];
      }
      else
      {
         distAD.f[W   ] = &distributionsAD[E   *size_Mat];
         distAD.f[E   ] = &distributionsAD[W   *size_Mat];
         distAD.f[S   ] = &distributionsAD[N   *size_Mat];
         distAD.f[N   ] = &distributionsAD[S   *size_Mat];
         distAD.f[B   ] = &distributionsAD[T   *size_Mat];
         distAD.f[T   ] = &distributionsAD[B   *size_Mat];
         distAD.f[SW  ] = &distributionsAD[NE  *size_Mat];
         distAD.f[NE  ] = &distributionsAD[SW  *size_Mat];
         distAD.f[NW  ] = &distributionsAD[SE  *size_Mat];
         distAD.f[SE  ] = &distributionsAD[NW  *size_Mat];
         distAD.f[BW  ] = &distributionsAD[TE  *size_Mat];
         distAD.f[TE  ] = &distributionsAD[BW  *size_Mat];
         distAD.f[TW  ] = &distributionsAD[BE  *size_Mat];
         distAD.f[BE  ] = &distributionsAD[TW  *size_Mat];
         distAD.f[BS  ] = &distributionsAD[TN  *size_Mat];
         distAD.f[TN  ] = &distributionsAD[BS  *size_Mat];
         distAD.f[TS  ] = &distributionsAD[BN  *size_Mat];
         distAD.f[BN  ] = &distributionsAD[TS  *size_Mat];
         distAD.f[dirREST] = &distributionsAD[dirREST*size_Mat];
         distAD.f[TNE ] = &distributionsAD[BSW *size_Mat];
         distAD.f[TSW ] = &distributionsAD[BNE *size_Mat];
         distAD.f[TSE ] = &distributionsAD[BNW *size_Mat];
         distAD.f[TNW ] = &distributionsAD[BSE *size_Mat];
         distAD.f[BNE ] = &distributionsAD[TSW *size_Mat];
         distAD.f[BSW ] = &distributionsAD[TNE *size_Mat];
         distAD.f[BSE ] = &distributionsAD[TNW *size_Mat];
         distAD.f[BNW ] = &distributionsAD[TSE *size_Mat];
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
	  real mfcbb = (distAD.f[E   ])[ke  ];
	  real mfabb = (distAD.f[W   ])[kw  ];
	  real mfbcb = (distAD.f[N   ])[kn  ];
	  real mfbab = (distAD.f[S   ])[ks  ];
	  real mfbbc = (distAD.f[T   ])[kt  ];
	  real mfbba = (distAD.f[B   ])[kb  ];
	  real mfccb = (distAD.f[NE  ])[kne ];
	  real mfaab = (distAD.f[SW  ])[ksw ];
	  real mfcab = (distAD.f[SE  ])[kse ];
	  real mfacb = (distAD.f[NW  ])[knw ];
	  real mfcbc = (distAD.f[TE  ])[kte ];
	  real mfaba = (distAD.f[BW  ])[kbw ];
	  real mfcba = (distAD.f[BE  ])[kbe ];
	  real mfabc = (distAD.f[TW  ])[ktw ];
	  real mfbcc = (distAD.f[TN  ])[ktn ];
	  real mfbaa = (distAD.f[BS  ])[kbs ];
	  real mfbca = (distAD.f[BN  ])[kbn ];
	  real mfbac = (distAD.f[TS  ])[kts ];
	  real mfbbb = (distAD.f[dirREST])[k   ];
	  real mfccc = (distAD.f[TNE ])[ktne];
	  real mfaac = (distAD.f[TSW ])[ktsw];
	  real mfcac = (distAD.f[TSE ])[ktse];
	  real mfacc = (distAD.f[TNW ])[ktnw];
	  real mfcca = (distAD.f[BNE ])[kbne];
	  real mfaaa = (distAD.f[BSW ])[kbsw];
	  real mfcaa = (distAD.f[BSE ])[kbse];
	  real mfaca = (distAD.f[BNW ])[kbnw];
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
extern "C" __global__ void CalcConc7( real* Conc,
                                          unsigned int* geoD,
                                          unsigned int* neighborX,
                                          unsigned int* neighborY,
                                          unsigned int* neighborZ,
                                          unsigned int size_Mat,
                                          real* DD7,
                                          bool isEvenTimestep)
{
   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   } 
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
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
// extern "C" __global__ void LBCalcMacThS27(real* Conc,
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
//       D27.f[E   ] = &DD27[E   *size_Mat];
//       D27.f[W   ] = &DD27[W   *size_Mat];
//       D27.f[N   ] = &DD27[N   *size_Mat];
//       D27.f[S   ] = &DD27[S   *size_Mat];
//       D27.f[T   ] = &DD27[T   *size_Mat];
//       D27.f[B   ] = &DD27[B   *size_Mat];
//       D27.f[NE  ] = &DD27[NE  *size_Mat];
//       D27.f[SW  ] = &DD27[SW  *size_Mat];
//       D27.f[SE  ] = &DD27[SE  *size_Mat];
//       D27.f[NW  ] = &DD27[NW  *size_Mat];
//       D27.f[TE  ] = &DD27[TE  *size_Mat];
//       D27.f[BW  ] = &DD27[BW  *size_Mat];
//       D27.f[BE  ] = &DD27[BE  *size_Mat];
//       D27.f[TW  ] = &DD27[TW  *size_Mat];
//       D27.f[TN  ] = &DD27[TN  *size_Mat];
//       D27.f[BS  ] = &DD27[BS  *size_Mat];
//       D27.f[BN  ] = &DD27[BN  *size_Mat];
//       D27.f[TS  ] = &DD27[TS  *size_Mat];
//       D27.f[dirREST] = &DD27[dirREST*size_Mat];
//       D27.f[TNE ] = &DD27[TNE *size_Mat];
//       D27.f[TSW ] = &DD27[TSW *size_Mat];
//       D27.f[TSE ] = &DD27[TSE *size_Mat];
//       D27.f[TNW ] = &DD27[TNW *size_Mat];
//       D27.f[BNE ] = &DD27[BNE *size_Mat];
//       D27.f[BSW ] = &DD27[BSW *size_Mat];
//       D27.f[BSE ] = &DD27[BSE *size_Mat];
//       D27.f[BNW ] = &DD27[BNW *size_Mat];
//    }
//    else
//    {
//       D27.f[W   ] = &DD27[E   *size_Mat];
//       D27.f[E   ] = &DD27[W   *size_Mat];
//       D27.f[S   ] = &DD27[N   *size_Mat];
//       D27.f[N   ] = &DD27[S   *size_Mat];
//       D27.f[B   ] = &DD27[T   *size_Mat];
//       D27.f[T   ] = &DD27[B   *size_Mat];
//       D27.f[SW  ] = &DD27[NE  *size_Mat];
//       D27.f[NE  ] = &DD27[SW  *size_Mat];
//       D27.f[NW  ] = &DD27[SE  *size_Mat];
//       D27.f[SE  ] = &DD27[NW  *size_Mat];
//       D27.f[BW  ] = &DD27[TE  *size_Mat];
//       D27.f[TE  ] = &DD27[BW  *size_Mat];
//       D27.f[TW  ] = &DD27[BE  *size_Mat];
//       D27.f[BE  ] = &DD27[TW  *size_Mat];
//       D27.f[BS  ] = &DD27[TN  *size_Mat];
//       D27.f[TN  ] = &DD27[BS  *size_Mat];
//       D27.f[TS  ] = &DD27[BN  *size_Mat];
//       D27.f[BN  ] = &DD27[TS  *size_Mat];
//       D27.f[dirREST] = &DD27[dirREST*size_Mat];
//       D27.f[BSW ] = &DD27[TNE *size_Mat];
//       D27.f[BNE ] = &DD27[TSW *size_Mat];
//       D27.f[BNW ] = &DD27[TSE *size_Mat];
//       D27.f[BSE ] = &DD27[TNW *size_Mat];
//       D27.f[TSW ] = &DD27[BNE *size_Mat];
//       D27.f[TNE ] = &DD27[BSW *size_Mat];
//       D27.f[TNW ] = &DD27[BSE *size_Mat];
//       D27.f[TSE ] = &DD27[BNW *size_Mat];
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
//          Conc[k]    =   (D27.f[E   ])[ke  ]+ (D27.f[W   ])[kw  ]+ 
//                         (D27.f[N   ])[kn  ]+ (D27.f[S   ])[ks  ]+
//                         (D27.f[T   ])[kt  ]+ (D27.f[B   ])[kb  ]+
//                         (D27.f[NE  ])[kne ]+ (D27.f[SW  ])[ksw ]+
//                         (D27.f[SE  ])[kse ]+ (D27.f[NW  ])[knw ]+
//                         (D27.f[TE  ])[kte ]+ (D27.f[BW  ])[kbw ]+
//                         (D27.f[BE  ])[kbe ]+ (D27.f[TW  ])[ktw ]+
//                         (D27.f[TN  ])[ktn ]+ (D27.f[BS  ])[kbs ]+
//                         (D27.f[BN  ])[kbn ]+ (D27.f[TS  ])[kts ]+
//                         (D27.f[dirREST])[kzero]+ 
//                         (D27.f[TNE ])[ktne]+ (D27.f[TSW ])[ktsw]+
//                         (D27.f[TSE ])[ktse]+ (D27.f[TNW ])[ktnw]+
//                         (D27.f[BNE ])[kbne]+ (D27.f[BSW ])[kbsw]+
//                         (D27.f[BSE ])[kbse]+ (D27.f[BNW ])[kbnw];
//       }
//    }   
// }



















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void GetPlaneConc7(real* Conc,
								            int* kPC,
								            unsigned int numberOfPointskPC,
											unsigned int* geoD,
											unsigned int* neighborX,
											unsigned int* neighborY,
											unsigned int* neighborZ,
											unsigned int size_Mat,
											real* DD7,
											bool isEvenTimestep)
{
   Distributions7 D7;
   if (isEvenTimestep==true)
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[1] = &DD7[1*size_Mat];
      D7.f[2] = &DD7[2*size_Mat];
      D7.f[3] = &DD7[3*size_Mat];
      D7.f[4] = &DD7[4*size_Mat];
      D7.f[5] = &DD7[5*size_Mat];
      D7.f[6] = &DD7[6*size_Mat];
   } 
   else
   {
      D7.f[0] = &DD7[0*size_Mat];
      D7.f[2] = &DD7[1*size_Mat];
      D7.f[1] = &DD7[2*size_Mat];
      D7.f[4] = &DD7[3*size_Mat];
      D7.f[3] = &DD7[4*size_Mat];
      D7.f[6] = &DD7[5*size_Mat];
      D7.f[5] = &DD7[6*size_Mat];
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
extern "C" __global__ void GetPlaneConc27(real* Conc,
								             int* kPC,
								             unsigned int numberOfPointskPC,
											 unsigned int* geoD,
											 unsigned int* neighborX,
											 unsigned int* neighborY,
											 unsigned int* neighborZ,
											 unsigned int size_Mat,
											 real* DD27,
											 bool isEvenTimestep)
{
   Distributions27 D27;
   if (isEvenTimestep==true)
   {
      D27.f[E   ] = &DD27[E   *size_Mat];
      D27.f[W   ] = &DD27[W   *size_Mat];
      D27.f[N   ] = &DD27[N   *size_Mat];
      D27.f[S   ] = &DD27[S   *size_Mat];
      D27.f[T   ] = &DD27[T   *size_Mat];
      D27.f[B   ] = &DD27[B   *size_Mat];
      D27.f[NE  ] = &DD27[NE  *size_Mat];
      D27.f[SW  ] = &DD27[SW  *size_Mat];
      D27.f[SE  ] = &DD27[SE  *size_Mat];
      D27.f[NW  ] = &DD27[NW  *size_Mat];
      D27.f[TE  ] = &DD27[TE  *size_Mat];
      D27.f[BW  ] = &DD27[BW  *size_Mat];
      D27.f[BE  ] = &DD27[BE  *size_Mat];
      D27.f[TW  ] = &DD27[TW  *size_Mat];
      D27.f[TN  ] = &DD27[TN  *size_Mat];
      D27.f[BS  ] = &DD27[BS  *size_Mat];
      D27.f[BN  ] = &DD27[BN  *size_Mat];
      D27.f[TS  ] = &DD27[TS  *size_Mat];
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
      D27.f[TNE ] = &DD27[TNE *size_Mat];
      D27.f[TSW ] = &DD27[TSW *size_Mat];
      D27.f[TSE ] = &DD27[TSE *size_Mat];
      D27.f[TNW ] = &DD27[TNW *size_Mat];
      D27.f[BNE ] = &DD27[BNE *size_Mat];
      D27.f[BSW ] = &DD27[BSW *size_Mat];
      D27.f[BSE ] = &DD27[BSE *size_Mat];
      D27.f[BNW ] = &DD27[BNW *size_Mat];
   }
   else
   {
      D27.f[W   ] = &DD27[E   *size_Mat];
      D27.f[E   ] = &DD27[W   *size_Mat];
      D27.f[S   ] = &DD27[N   *size_Mat];
      D27.f[N   ] = &DD27[S   *size_Mat];
      D27.f[B   ] = &DD27[T   *size_Mat];
      D27.f[T   ] = &DD27[B   *size_Mat];
      D27.f[SW  ] = &DD27[NE  *size_Mat];
      D27.f[NE  ] = &DD27[SW  *size_Mat];
      D27.f[NW  ] = &DD27[SE  *size_Mat];
      D27.f[SE  ] = &DD27[NW  *size_Mat];
      D27.f[BW  ] = &DD27[TE  *size_Mat];
      D27.f[TE  ] = &DD27[BW  *size_Mat];
      D27.f[TW  ] = &DD27[BE  *size_Mat];
      D27.f[BE  ] = &DD27[TW  *size_Mat];
      D27.f[BS  ] = &DD27[TN  *size_Mat];
      D27.f[TN  ] = &DD27[BS  *size_Mat];
      D27.f[TS  ] = &DD27[BN  *size_Mat];
      D27.f[BN  ] = &DD27[TS  *size_Mat];
      D27.f[dirREST] = &DD27[dirREST*size_Mat];
      D27.f[BSW ] = &DD27[TNE *size_Mat];
      D27.f[BNE ] = &DD27[TSW *size_Mat];
      D27.f[BNW ] = &DD27[TSE *size_Mat];
      D27.f[BSE ] = &DD27[TNW *size_Mat];
      D27.f[TSW ] = &DD27[BNE *size_Mat];
      D27.f[TNE ] = &DD27[BSW *size_Mat];
      D27.f[TNW ] = &DD27[BSE *size_Mat];
      D27.f[TSE ] = &DD27[BNW *size_Mat];
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
         Conc[k]    =   (D27.f[E   ])[ke  ]+ (D27.f[W   ])[kw  ]+ 
                        (D27.f[N   ])[kn  ]+ (D27.f[S   ])[ks  ]+
                        (D27.f[T   ])[kt  ]+ (D27.f[B   ])[kb  ]+
                        (D27.f[NE  ])[kne ]+ (D27.f[SW  ])[ksw ]+
                        (D27.f[SE  ])[kse ]+ (D27.f[NW  ])[knw ]+
                        (D27.f[TE  ])[kte ]+ (D27.f[BW  ])[kbw ]+
                        (D27.f[BE  ])[kbe ]+ (D27.f[TW  ])[ktw ]+
                        (D27.f[TN  ])[ktn ]+ (D27.f[BS  ])[kbs ]+
                        (D27.f[BN  ])[kbn ]+ (D27.f[TS  ])[kts ]+
                        (D27.f[dirREST])[kzero]+ 
                        (D27.f[TNE ])[ktne]+ (D27.f[TSW ])[ktsw]+
                        (D27.f[TSE ])[ktse]+ (D27.f[TNW ])[ktnw]+
                        (D27.f[BNE ])[kbne]+ (D27.f[BSW ])[kbsw]+
                        (D27.f[BSE ])[kbse]+ (D27.f[BNW ])[kbnw];
      }
   }   
}