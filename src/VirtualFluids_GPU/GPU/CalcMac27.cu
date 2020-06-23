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
//! \file CalcMac27.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LBCalcMacCompSP27( 
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	real* rho,
	real* pressure,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat,
	real* distributions,
	bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
   //! The velocity boundary condition is executed in the following steps
   //!
   ////////////////////////////////////////////////////////////////////////////////
   //! - Get node index coordinates from thredIdx, blockIdx, blockDim and gridDim.
   //!
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

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
      Distributions27 dist;
      if (isEvenTimestep)
      {
         dist.f[dirE   ] = &distributions[dirE   *size_Mat];
         dist.f[dirW   ] = &distributions[dirW   *size_Mat];
         dist.f[dirN   ] = &distributions[dirN   *size_Mat];
         dist.f[dirS   ] = &distributions[dirS   *size_Mat];
         dist.f[dirT   ] = &distributions[dirT   *size_Mat];
         dist.f[dirB   ] = &distributions[dirB   *size_Mat];
         dist.f[dirNE  ] = &distributions[dirNE  *size_Mat];
         dist.f[dirSW  ] = &distributions[dirSW  *size_Mat];
         dist.f[dirSE  ] = &distributions[dirSE  *size_Mat];
         dist.f[dirNW  ] = &distributions[dirNW  *size_Mat];
         dist.f[dirTE  ] = &distributions[dirTE  *size_Mat];
         dist.f[dirBW  ] = &distributions[dirBW  *size_Mat];
         dist.f[dirBE  ] = &distributions[dirBE  *size_Mat];
         dist.f[dirTW  ] = &distributions[dirTW  *size_Mat];
         dist.f[dirTN  ] = &distributions[dirTN  *size_Mat];
         dist.f[dirBS  ] = &distributions[dirBS  *size_Mat];
         dist.f[dirBN  ] = &distributions[dirBN  *size_Mat];
         dist.f[dirTS  ] = &distributions[dirTS  *size_Mat];
         dist.f[dirREST] = &distributions[dirREST*size_Mat];
         dist.f[dirTNE ] = &distributions[dirTNE *size_Mat];
         dist.f[dirTSW ] = &distributions[dirTSW *size_Mat];
         dist.f[dirTSE ] = &distributions[dirTSE *size_Mat];
         dist.f[dirTNW ] = &distributions[dirTNW *size_Mat];
         dist.f[dirBNE ] = &distributions[dirBNE *size_Mat];
         dist.f[dirBSW ] = &distributions[dirBSW *size_Mat];
         dist.f[dirBSE ] = &distributions[dirBSE *size_Mat];
         dist.f[dirBNW ] = &distributions[dirBNW *size_Mat];
      } 
      else
      {
         dist.f[dirW   ] = &distributions[dirE   *size_Mat];
         dist.f[dirE   ] = &distributions[dirW   *size_Mat];
         dist.f[dirS   ] = &distributions[dirN   *size_Mat];
         dist.f[dirN   ] = &distributions[dirS   *size_Mat];
         dist.f[dirB   ] = &distributions[dirT   *size_Mat];
         dist.f[dirT   ] = &distributions[dirB   *size_Mat];
         dist.f[dirSW  ] = &distributions[dirNE  *size_Mat];
         dist.f[dirNE  ] = &distributions[dirSW  *size_Mat];
         dist.f[dirNW  ] = &distributions[dirSE  *size_Mat];
         dist.f[dirSE  ] = &distributions[dirNW  *size_Mat];
         dist.f[dirBW  ] = &distributions[dirTE  *size_Mat];
         dist.f[dirTE  ] = &distributions[dirBW  *size_Mat];
         dist.f[dirTW  ] = &distributions[dirBE  *size_Mat];
         dist.f[dirBE  ] = &distributions[dirTW  *size_Mat];
         dist.f[dirBS  ] = &distributions[dirTN  *size_Mat];
         dist.f[dirTN  ] = &distributions[dirBS  *size_Mat];
         dist.f[dirTS  ] = &distributions[dirBN  *size_Mat];
         dist.f[dirBN  ] = &distributions[dirTS  *size_Mat];
         dist.f[dirREST] = &distributions[dirREST*size_Mat];
         dist.f[dirTNE ] = &distributions[dirBSW *size_Mat];
         dist.f[dirTSW ] = &distributions[dirBNE *size_Mat];
         dist.f[dirTSE ] = &distributions[dirBNW *size_Mat];
         dist.f[dirTNW ] = &distributions[dirBSE *size_Mat];
         dist.f[dirBNE ] = &distributions[dirTSW *size_Mat];
         dist.f[dirBSW ] = &distributions[dirTNE *size_Mat];
         dist.f[dirBSE ] = &distributions[dirTNW *size_Mat];
         dist.f[dirBNW ] = &distributions[dirTSE *size_Mat];
      }
	  ////////////////////////////////////////////////////////////////////////////////
	  //! - Set neighbor indices (necessary for indirect addressing)  
	  //!
	  uint ke = k;
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
	  real mfcbb = (dist.f[dirE   ])[k   ];
	  real mfabb = (dist.f[dirW   ])[kw  ];
	  real mfbcb = (dist.f[dirN   ])[k   ];
	  real mfbab = (dist.f[dirS   ])[ks  ];
	  real mfbbc = (dist.f[dirT   ])[k   ];
	  real mfbba = (dist.f[dirB   ])[kb  ];
	  real mfccb = (dist.f[dirNE  ])[k   ];
	  real mfaab = (dist.f[dirSW  ])[ksw ];
	  real mfcab = (dist.f[dirSE  ])[ks  ];
	  real mfacb = (dist.f[dirNW  ])[kw  ];
	  real mfcbc = (dist.f[dirTE  ])[k   ];
	  real mfaba = (dist.f[dirBW  ])[kbw ];
	  real mfcba = (dist.f[dirBE  ])[kb  ];
	  real mfabc = (dist.f[dirTW  ])[kw  ];
	  real mfbcc = (dist.f[dirTN  ])[k   ];
	  real mfbaa = (dist.f[dirBS  ])[kbs ];
	  real mfbca = (dist.f[dirBN  ])[kb  ];
	  real mfbac = (dist.f[dirTS  ])[ks  ];
	  real mfbbb = (dist.f[dirREST])[k   ];
	  real mfccc = (dist.f[dirTNE ])[k   ];
	  real mfaac = (dist.f[dirTSW ])[ksw ];
	  real mfcac = (dist.f[dirTSE ])[ks  ];
	  real mfacc = (dist.f[dirTNW ])[kw  ];
	  real mfcca = (dist.f[dirBNE ])[kb  ];
	  real mfaaa = (dist.f[dirBSW ])[kbsw];
	  real mfcaa = (dist.f[dirBSE ])[kbs ];
	  real mfaca = (dist.f[dirBNW ])[kbw ];
	  ////////////////////////////////////////////////////////////////////////////////
	  //! - Set pressure, density and velocity to \f$ 1.0 \f$  
	  //!
	  pressure[k]  = c0o1;
	  rho[k]       = c0o1;
	  velocityX[k] = c0o1;
	  velocityY[k] = c0o1;
	  velocityZ[k] = c0o1;

      //////////////////////////////////////////////////////////////////////////
	  //! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
	  //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
	  //!
	  real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
      	(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
      	((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
      
      real rhoDev = c1o1 + drho;
      ////////////////////////////////////////////////////////////////////////////////////
      real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
      	(((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
      	(mfcbb - mfabb)) / rhoDev;
      real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
      	(((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
      	(mfbcb - mfbab)) / rhoDev;
      real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
      	(((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
      	(mfbbc - mfbba)) / rhoDev;
      
      rho[k]        = drho;
      velocityX[k]  = vvx;
      velocityY[k]  = vvy;
      velocityZ[k]  = vvz;
      
	  //////////////////////////////////////////////////////////////////////////
	  //! - Calculate pressure
	  //!
	  real OxxPyyPzz = c1o1;
      pressure[k] =((dist.f[dirE   ])[ke  ]+ (dist.f[dirW   ])[kw  ]+ 
                    (dist.f[dirN   ])[kn  ]+ (dist.f[dirS   ])[ks  ]+
                    (dist.f[dirT   ])[kt  ]+ (dist.f[dirB   ])[kb  ]+
                    c2o1*(
                    (dist.f[dirNE  ])[kne ]+ (dist.f[dirSW  ])[ksw ]+
                    (dist.f[dirSE  ])[kse ]+ (dist.f[dirNW  ])[knw ]+
                    (dist.f[dirTE  ])[kte ]+ (dist.f[dirBW  ])[kbw ]+
                    (dist.f[dirBE  ])[kbe ]+ (dist.f[dirTW  ])[ktw ]+
                    (dist.f[dirTN  ])[ktn ]+ (dist.f[dirBS  ])[kbs ]+
                    (dist.f[dirBN  ])[kbn ]+ (dist.f[dirTS  ])[kts ])+
                    c3o1*(
                    (dist.f[dirTNE ])[ktne]+ (dist.f[dirTSW ])[ktsw]+ 
                    (dist.f[dirTSE ])[ktse]+ (dist.f[dirTNW ])[ktnw]+ 
                    (dist.f[dirBNE ])[kbne]+ (dist.f[dirBSW ])[kbsw]+ 
                    (dist.f[dirBSE ])[kbse]+ (dist.f[dirBNW ])[kbnw])-
                    rho[k]-(velocityX[k] * velocityX[k] + velocityY[k] * velocityY[k] + velocityZ[k] * velocityZ[k]) * (c1o1+rho[k])) * (c1o1 / OxxPyyPzz - c1o2) +rho[k];
   }
}





