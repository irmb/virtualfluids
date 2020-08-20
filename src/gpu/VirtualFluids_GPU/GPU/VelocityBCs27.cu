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
//! \file VelocityBCs27.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"

//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QVelDevPlainBB27(
	real* vx,
	real* vy,
	real* vz,
	real* distributions,
	int* k_Q, 
	real* QQ,
	uint sizeQ,
	int kQ, 
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint size_Mat, 
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
   // run for all indices in size of boundary condition (kQ)
   if(k<kQ)
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
	  //! - Set local velocities
	  //!
	  real VeloX = vx[k];
	  real VeloY = vy[k];
	  real VeloZ = vz[k];
      ////////////////////////////////////////////////////////////////////////////////
	  //! - Set local subgrid distances (q's)
	  //!
      real   *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			 *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			 *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			 *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			 *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
      ////////////////////////////////////////////////////////////////////////////////
	  //! - Set neighbor indices (necessary for indirect addressing)  
	  //!
	  uint KQK = k_Q[k];
      uint ke   = KQK;
      uint kw   = neighborX[KQK];
      uint kn   = KQK;
      uint ks   = neighborY[KQK];
      uint kt   = KQK;
      uint kb   = neighborZ[KQK];
      uint ksw  = neighborY[kw];
      uint kne  = KQK;
      uint kse  = ks;
      uint knw  = kw;
      uint kbw  = neighborZ[kw];
      uint kte  = KQK;
      uint kbe  = kb;
      uint ktw  = kw;
      uint kbs  = neighborZ[ks];
      uint ktn  = KQK;
      uint kbn  = kb;
      uint kts  = ks;
      uint ktse = ks;
      uint kbnw = kbw;
      uint ktnw = kw;
      uint kbse = kbs;
      uint ktsw = ksw;
      uint kbne = kb;
      uint ktne = KQK;
      uint kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
     
      ////////////////////////////////////////////////////////////////////////////////
	  //! - Set local distributions  
	  //!
	  real f_W = (dist.f[dirE])[ke];
      real f_E    = (dist.f[dirW   ])[kw   ];
      real f_S    = (dist.f[dirN   ])[kn   ];
      real f_N    = (dist.f[dirS   ])[ks   ];
      real f_B    = (dist.f[dirT   ])[kt   ];
      real f_T    = (dist.f[dirB   ])[kb   ];
      real f_SW   = (dist.f[dirNE  ])[kne  ];
      real f_NE   = (dist.f[dirSW  ])[ksw  ];
      real f_NW   = (dist.f[dirSE  ])[kse  ];
      real f_SE   = (dist.f[dirNW  ])[knw  ];
      real f_BW   = (dist.f[dirTE  ])[kte  ];
      real f_TE   = (dist.f[dirBW  ])[kbw  ];
      real f_TW   = (dist.f[dirBE  ])[kbe  ];
      real f_BE   = (dist.f[dirTW  ])[ktw  ];
      real f_BS   = (dist.f[dirTN  ])[ktn  ];
      real f_TN   = (dist.f[dirBS  ])[kbs  ];
      real f_TS   = (dist.f[dirBN  ])[kbn  ];
      real f_BN   = (dist.f[dirTS  ])[kts  ];
      real f_BSW  = (dist.f[dirTNE ])[ktne ];
      real f_BNE  = (dist.f[dirTSW ])[ktsw ];
      real f_BNW  = (dist.f[dirTSE ])[ktse ];
      real f_BSE  = (dist.f[dirTNW ])[ktnw ];
      real f_TSW  = (dist.f[dirBNE ])[kbne ];
      real f_TNE  = (dist.f[dirBSW ])[kbsw ];
      real f_TNW  = (dist.f[dirBSE ])[kbse ];
      real f_TSE  = (dist.f[dirBNW ])[kbnw ];
	  ////////////////////////////////////////////////////////////////////////////////

	  ////////////////////////////////////////////////////////////////////////////////
	  //! - change the pointer to write the results in the correct array  
	  //!
	  if (!isEvenTimestep)
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
	  //! - rewrite distributions if there is a sub-grid distance (q) in same direction
	  real q;
      q = q_dirE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirW  ])[kw  ]=f_E   + c4o9  * (-VeloX);	
      q = q_dirW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirE  ])[ke  ]=f_W   + c4o9  * ( VeloX);	
      q = q_dirN[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirS  ])[ks  ]=f_N   + c4o9  * (-VeloY);	
      q = q_dirS[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirN  ])[kn  ]=f_S   + c4o9  * ( VeloY);	
      q = q_dirT[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirB  ])[kb  ]=f_T   + c4o9  * (-VeloZ);
      q = q_dirB[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirT  ])[kt  ]=f_B   + c4o9  * ( VeloZ);
      q = q_dirNE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirSW ])[ksw ]=f_NE  + c1o9  * (-VeloX - VeloY);
	  q = q_dirSW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirNE ])[kne ]=f_SW  + c1o9  * ( VeloX + VeloY);
	  q = q_dirSE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirNW ])[knw ]=f_SE  + c1o9  * (-VeloX + VeloY);
	  q = q_dirNW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirSE ])[kse ]=f_NW  + c1o9  * ( VeloX - VeloY);
	  q = q_dirTE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBW ])[kbw ]=f_TE  + c1o9  * (-VeloX - VeloZ);
	  q = q_dirBW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTE ])[kte ]=f_BW  + c1o9  * ( VeloX + VeloZ);
	  q = q_dirBE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTW ])[ktw ]=f_BE  + c1o9  * (-VeloX + VeloZ);
	  q = q_dirTW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBE ])[kbe ]=f_TW  + c1o9  * ( VeloX - VeloZ);
	  q = q_dirTN[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBS ])[kbs ]=f_TN  + c1o9  * (-VeloY - VeloZ);
	  q = q_dirBS[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTN ])[ktn ]=f_BS  + c1o9  * ( VeloY + VeloZ);
	  q = q_dirBN[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTS ])[kts ]=f_BN  + c1o9  * (-VeloY + VeloZ);
	  q = q_dirTS[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBN ])[kbn ]=f_TS  + c1o9  * ( VeloY - VeloZ);
      q = q_dirTNE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBSW])[kbsw]=f_TNE + c1o36 * (-VeloX - VeloY - VeloZ);
      q = q_dirBSW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTNE])[ktne]=f_BSW + c1o36 * ( VeloX + VeloY + VeloZ);
      q = q_dirBNE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTSW])[ktsw]=f_BNE + c1o36 * (-VeloX - VeloY + VeloZ);
      q = q_dirTSW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBNE])[kbne]=f_TSW + c1o36 * ( VeloX + VeloY - VeloZ);
      q = q_dirTSE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBNW])[kbnw]=f_TSE + c1o36 * (-VeloX + VeloY - VeloZ);
      q = q_dirBNW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTSE])[ktse]=f_BNW + c1o36 * ( VeloX - VeloY + VeloZ);
      q = q_dirBSE[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirTNW])[ktnw]=f_BSE + c1o36 * (-VeloX + VeloY + VeloZ);
      q = q_dirTNW[k];	if (q>=c0o1 && q<=c1o1)	(dist.f[dirBSE])[kbse]=f_TNW + c1o36 * ( VeloX - VeloY - VeloZ);
   }
}



