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
   //! - Get node index coordinates from thredIdx, blockIdx, blockDim and gridDim.
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
         distAD.f[dirE   ] = &distributionsAD[dirE   *size_Mat];
         distAD.f[dirW   ] = &distributionsAD[dirW   *size_Mat];
         distAD.f[dirN   ] = &distributionsAD[dirN   *size_Mat];
         distAD.f[dirS   ] = &distributionsAD[dirS   *size_Mat];
         distAD.f[dirT   ] = &distributionsAD[dirT   *size_Mat];
         distAD.f[dirB   ] = &distributionsAD[dirB   *size_Mat];
         distAD.f[dirNE  ] = &distributionsAD[dirNE  *size_Mat];
         distAD.f[dirSW  ] = &distributionsAD[dirSW  *size_Mat];
         distAD.f[dirSE  ] = &distributionsAD[dirSE  *size_Mat];
         distAD.f[dirNW  ] = &distributionsAD[dirNW  *size_Mat];
         distAD.f[dirTE  ] = &distributionsAD[dirTE  *size_Mat];
         distAD.f[dirBW  ] = &distributionsAD[dirBW  *size_Mat];
         distAD.f[dirBE  ] = &distributionsAD[dirBE  *size_Mat];
         distAD.f[dirTW  ] = &distributionsAD[dirTW  *size_Mat];
         distAD.f[dirTN  ] = &distributionsAD[dirTN  *size_Mat];
         distAD.f[dirBS  ] = &distributionsAD[dirBS  *size_Mat];
         distAD.f[dirBN  ] = &distributionsAD[dirBN  *size_Mat];
         distAD.f[dirTS  ] = &distributionsAD[dirTS  *size_Mat];
         distAD.f[dirREST] = &distributionsAD[dirREST*size_Mat];
         distAD.f[dirTNE ] = &distributionsAD[dirTNE *size_Mat];
         distAD.f[dirTSW ] = &distributionsAD[dirTSW *size_Mat];
         distAD.f[dirTSE ] = &distributionsAD[dirTSE *size_Mat];
         distAD.f[dirTNW ] = &distributionsAD[dirTNW *size_Mat];
         distAD.f[dirBNE ] = &distributionsAD[dirBNE *size_Mat];
         distAD.f[dirBSW ] = &distributionsAD[dirBSW *size_Mat];
         distAD.f[dirBSE ] = &distributionsAD[dirBSE *size_Mat];
         distAD.f[dirBNW ] = &distributionsAD[dirBNW *size_Mat];
      }
      else
      {
         distAD.f[dirW   ] = &distributionsAD[dirE   *size_Mat];
         distAD.f[dirE   ] = &distributionsAD[dirW   *size_Mat];
         distAD.f[dirS   ] = &distributionsAD[dirN   *size_Mat];
         distAD.f[dirN   ] = &distributionsAD[dirS   *size_Mat];
         distAD.f[dirB   ] = &distributionsAD[dirT   *size_Mat];
         distAD.f[dirT   ] = &distributionsAD[dirB   *size_Mat];
         distAD.f[dirSW  ] = &distributionsAD[dirNE  *size_Mat];
         distAD.f[dirNE  ] = &distributionsAD[dirSW  *size_Mat];
         distAD.f[dirNW  ] = &distributionsAD[dirSE  *size_Mat];
         distAD.f[dirSE  ] = &distributionsAD[dirNW  *size_Mat];
         distAD.f[dirBW  ] = &distributionsAD[dirTE  *size_Mat];
         distAD.f[dirTE  ] = &distributionsAD[dirBW  *size_Mat];
         distAD.f[dirTW  ] = &distributionsAD[dirBE  *size_Mat];
         distAD.f[dirBE  ] = &distributionsAD[dirTW  *size_Mat];
         distAD.f[dirBS  ] = &distributionsAD[dirTN  *size_Mat];
         distAD.f[dirTN  ] = &distributionsAD[dirBS  *size_Mat];
         distAD.f[dirTS  ] = &distributionsAD[dirBN  *size_Mat];
         distAD.f[dirBN  ] = &distributionsAD[dirTS  *size_Mat];
         distAD.f[dirREST] = &distributionsAD[dirREST*size_Mat];
         distAD.f[dirTNE ] = &distributionsAD[dirBSW *size_Mat];
         distAD.f[dirTSW ] = &distributionsAD[dirBNE *size_Mat];
         distAD.f[dirTSE ] = &distributionsAD[dirBNW *size_Mat];
         distAD.f[dirTNW ] = &distributionsAD[dirBSE *size_Mat];
         distAD.f[dirBNE ] = &distributionsAD[dirTSW *size_Mat];
         distAD.f[dirBSW ] = &distributionsAD[dirTNE *size_Mat];
         distAD.f[dirBSE ] = &distributionsAD[dirTNW *size_Mat];
         distAD.f[dirBNW ] = &distributionsAD[dirTSE *size_Mat];
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
	  real mfcbb = (distAD.f[dirE   ])[ke  ];
	  real mfabb = (distAD.f[dirW   ])[kw  ];
	  real mfbcb = (distAD.f[dirN   ])[kn  ];
	  real mfbab = (distAD.f[dirS   ])[ks  ];
	  real mfbbc = (distAD.f[dirT   ])[kt  ];
	  real mfbba = (distAD.f[dirB   ])[kb  ];
	  real mfccb = (distAD.f[dirNE  ])[kne ];
	  real mfaab = (distAD.f[dirSW  ])[ksw ];
	  real mfcab = (distAD.f[dirSE  ])[kse ];
	  real mfacb = (distAD.f[dirNW  ])[knw ];
	  real mfcbc = (distAD.f[dirTE  ])[kte ];
	  real mfaba = (distAD.f[dirBW  ])[kbw ];
	  real mfcba = (distAD.f[dirBE  ])[kbe ];
	  real mfabc = (distAD.f[dirTW  ])[ktw ];
	  real mfbcc = (distAD.f[dirTN  ])[ktn ];
	  real mfbaa = (distAD.f[dirBS  ])[kbs ];
	  real mfbca = (distAD.f[dirBN  ])[kbn ];
	  real mfbac = (distAD.f[dirTS  ])[kts ];
	  real mfbbb = (distAD.f[dirREST])[k   ];
	  real mfccc = (distAD.f[dirTNE ])[ktne];
	  real mfaac = (distAD.f[dirTSW ])[ktsw];
	  real mfcac = (distAD.f[dirTSE ])[ktse];
	  real mfacc = (distAD.f[dirTNW ])[ktnw];
	  real mfcca = (distAD.f[dirBNE ])[kbne];
	  real mfaaa = (distAD.f[dirBSW ])[kbsw];
	  real mfcaa = (distAD.f[dirBSE ])[kbse];
	  real mfaca = (distAD.f[dirBNW ])[kbnw];
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
//       D27.f[dirE   ] = &DD27[dirE   *size_Mat];
//       D27.f[dirW   ] = &DD27[dirW   *size_Mat];
//       D27.f[dirN   ] = &DD27[dirN   *size_Mat];
//       D27.f[dirS   ] = &DD27[dirS   *size_Mat];
//       D27.f[dirT   ] = &DD27[dirT   *size_Mat];
//       D27.f[dirB   ] = &DD27[dirB   *size_Mat];
//       D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
//       D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
//       D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
//       D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
//       D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
//       D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
//       D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
//       D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
//       D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
//       D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
//       D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
//       D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
//       D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
//       D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
//       D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
//       D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
//       D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
//       D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
//       D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
//       D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
//       D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
//    }
//    else
//    {
//       D27.f[dirW   ] = &DD27[dirE   *size_Mat];
//       D27.f[dirE   ] = &DD27[dirW   *size_Mat];
//       D27.f[dirS   ] = &DD27[dirN   *size_Mat];
//       D27.f[dirN   ] = &DD27[dirS   *size_Mat];
//       D27.f[dirB   ] = &DD27[dirT   *size_Mat];
//       D27.f[dirT   ] = &DD27[dirB   *size_Mat];
//       D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
//       D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
//       D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
//       D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
//       D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
//       D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
//       D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
//       D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
//       D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
//       D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
//       D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
//       D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
//       D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
//       D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
//       D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
//       D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
//       D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
//       D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
//       D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
//       D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
//       D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
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
//          Conc[k]    =   (D27.f[dirE   ])[ke  ]+ (D27.f[dirW   ])[kw  ]+ 
//                         (D27.f[dirN   ])[kn  ]+ (D27.f[dirS   ])[ks  ]+
//                         (D27.f[dirT   ])[kt  ]+ (D27.f[dirB   ])[kb  ]+
//                         (D27.f[dirNE  ])[kne ]+ (D27.f[dirSW  ])[ksw ]+
//                         (D27.f[dirSE  ])[kse ]+ (D27.f[dirNW  ])[knw ]+
//                         (D27.f[dirTE  ])[kte ]+ (D27.f[dirBW  ])[kbw ]+
//                         (D27.f[dirBE  ])[kbe ]+ (D27.f[dirTW  ])[ktw ]+
//                         (D27.f[dirTN  ])[ktn ]+ (D27.f[dirBS  ])[kbs ]+
//                         (D27.f[dirBN  ])[kbn ]+ (D27.f[dirTS  ])[kts ]+
//                         (D27.f[dirZERO])[kzero]+ 
//                         (D27.f[dirTNE ])[ktne]+ (D27.f[dirTSW ])[ktsw]+
//                         (D27.f[dirTSE ])[ktse]+ (D27.f[dirTNW ])[ktnw]+
//                         (D27.f[dirBNE ])[kbne]+ (D27.f[dirBSW ])[kbsw]+
//                         (D27.f[dirBSE ])[kbse]+ (D27.f[dirBNW ])[kbnw];
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
      D27.f[dirE   ] = &DD27[dirE   *size_Mat];
      D27.f[dirW   ] = &DD27[dirW   *size_Mat];
      D27.f[dirN   ] = &DD27[dirN   *size_Mat];
      D27.f[dirS   ] = &DD27[dirS   *size_Mat];
      D27.f[dirT   ] = &DD27[dirT   *size_Mat];
      D27.f[dirB   ] = &DD27[dirB   *size_Mat];
      D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
      D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
      D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
      D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
      D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
      D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
      D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
      D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
   }
   else
   {
      D27.f[dirW   ] = &DD27[dirE   *size_Mat];
      D27.f[dirE   ] = &DD27[dirW   *size_Mat];
      D27.f[dirS   ] = &DD27[dirN   *size_Mat];
      D27.f[dirN   ] = &DD27[dirS   *size_Mat];
      D27.f[dirB   ] = &DD27[dirT   *size_Mat];
      D27.f[dirT   ] = &DD27[dirB   *size_Mat];
      D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
      D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
      D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
      D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
      D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
      D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
      D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
      D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
      D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
      D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
      D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
      D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
      D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
      D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
      D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
      D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
      D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
      D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
      D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
      D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
      D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
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
         Conc[k]    =   (D27.f[dirE   ])[ke  ]+ (D27.f[dirW   ])[kw  ]+ 
                        (D27.f[dirN   ])[kn  ]+ (D27.f[dirS   ])[ks  ]+
                        (D27.f[dirT   ])[kt  ]+ (D27.f[dirB   ])[kb  ]+
                        (D27.f[dirNE  ])[kne ]+ (D27.f[dirSW  ])[ksw ]+
                        (D27.f[dirSE  ])[kse ]+ (D27.f[dirNW  ])[knw ]+
                        (D27.f[dirTE  ])[kte ]+ (D27.f[dirBW  ])[kbw ]+
                        (D27.f[dirBE  ])[kbe ]+ (D27.f[dirTW  ])[ktw ]+
                        (D27.f[dirTN  ])[ktn ]+ (D27.f[dirBS  ])[kbs ]+
                        (D27.f[dirBN  ])[kbn ]+ (D27.f[dirTS  ])[kts ]+
                        (D27.f[dirZERO])[kzero]+ 
                        (D27.f[dirTNE ])[ktne]+ (D27.f[dirTSW ])[ktsw]+
                        (D27.f[dirTSE ])[ktse]+ (D27.f[dirTNW ])[ktnw]+
                        (D27.f[dirBNE ])[kbne]+ (D27.f[dirBSW ])[kbsw]+
                        (D27.f[dirBSE ])[kbse]+ (D27.f[dirBNW ])[kbnw];
      }
   }   
}