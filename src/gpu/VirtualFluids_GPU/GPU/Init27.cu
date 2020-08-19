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
//! \file Init27.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"

extern "C" __global__ void LBInit(
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* typeOfGridNode,
	real* rho,
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	uint size_Mat,
	real* distributions,
	bool isEvenTimestep)
{
	//////////////////////////////////////////////////////////////////////////
	//! The initialization is executed in the following steps
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
			dist.f[dirBSW ] = &distributions[dirTNE *size_Mat];
			dist.f[dirBNE ] = &distributions[dirTSW *size_Mat];
			dist.f[dirBNW ] = &distributions[dirTSE *size_Mat];
			dist.f[dirBSE ] = &distributions[dirTNW *size_Mat];
			dist.f[dirTSW ] = &distributions[dirBNE *size_Mat];
			dist.f[dirTNE ] = &distributions[dirBSW *size_Mat];
			dist.f[dirTNW ] = &distributions[dirBSE *size_Mat];
			dist.f[dirTSE ] = &distributions[dirBNW *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////
		//! - Set local velocities and density
		//!
		real drho = rho[k];
		real  vx1 = velocityX[k]; 
		real  vx2 = velocityY[k]; 
		real  vx3 = velocityZ[k]; 
		//////////////////////////////////////////////////////////////////////////
		//! - Set neighbor indices (necessary for indirect addressing)  
		//!
		uint kzero = k;
		uint ke    = k;
		uint kw    = neighborX[k];
		uint kn    = k;
		uint ks    = neighborY[k];
		uint kt    = k;
		uint kb    = neighborZ[k];
		uint ksw   = neighborY[kw];
		uint kne   = k;
		uint kse   = ks;
		uint knw   = kw;
		uint kbw   = neighborZ[kw];
		uint kte   = k;
		uint kbe   = kb;
		uint ktw   = kw;
		uint kbs   = neighborZ[ks];
		uint ktn   = k;
		uint kbn   = kb;
		uint kts   = ks;
		uint ktse  = ks;
		uint kbnw  = kbw;
		uint ktnw  = kw;
		uint kbse  = kbs;
		uint ktsw  = ksw;
		uint kbne  = kb;
		uint ktne  = k;
		uint kbsw  = neighborZ[ksw];
		//////////////////////////////////////////////////////////////////////////
		//! - Calculate the equilibrium and set the distributions
		//!
		real cu_sq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);

		(dist.f[dirREST])[kzero] = c8o27* (drho - cu_sq);
		(dist.f[dirE   ])[ke   ] = c2o27* (drho + c3o1*( vx1            ) + c9o2*( vx1            )*( vx1            ) - cu_sq);
		(dist.f[dirW   ])[kw   ] = c2o27* (drho + c3o1*(-vx1            ) + c9o2*(-vx1            )*(-vx1            ) - cu_sq);
		(dist.f[dirN   ])[kn   ] = c2o27* (drho + c3o1*(       vx2      ) + c9o2*(       vx2      )*(       vx2      ) - cu_sq);
		(dist.f[dirS   ])[ks   ] = c2o27* (drho + c3o1*(     - vx2      ) + c9o2*(     - vx2      )*(     - vx2      ) - cu_sq);
		(dist.f[dirT   ])[kt   ] = c2o27* (drho + c3o1*(             vx3) + c9o2*(             vx3)*(             vx3) - cu_sq);
		(dist.f[dirB   ])[kb   ] = c2o27* (drho + c3o1*(            -vx3) + c9o2*(            -vx3)*(            -vx3) - cu_sq);
		(dist.f[dirNE  ])[kne  ] = c1o54* (drho + c3o1*( vx1 + vx2      ) + c9o2*( vx1 + vx2      )*( vx1 + vx2      ) - cu_sq);
		(dist.f[dirSW  ])[ksw  ] = c1o54* (drho + c3o1*(-vx1 - vx2      ) + c9o2*(-vx1 - vx2      )*(-vx1 - vx2      ) - cu_sq);
		(dist.f[dirSE  ])[kse  ] = c1o54* (drho + c3o1*( vx1 - vx2      ) + c9o2*( vx1 - vx2      )*( vx1 - vx2      ) - cu_sq);
		(dist.f[dirNW  ])[knw  ] = c1o54* (drho + c3o1*(-vx1 + vx2      ) + c9o2*(-vx1 + vx2      )*(-vx1 + vx2      ) - cu_sq);
		(dist.f[dirTE  ])[kte  ] = c1o54* (drho + c3o1*( vx1       + vx3) + c9o2*( vx1       + vx3)*( vx1       + vx3) - cu_sq);
		(dist.f[dirBW  ])[kbw  ] = c1o54* (drho + c3o1*(-vx1       - vx3) + c9o2*(-vx1       - vx3)*(-vx1       - vx3) - cu_sq);
		(dist.f[dirBE  ])[kbe  ] = c1o54* (drho + c3o1*( vx1       - vx3) + c9o2*( vx1       - vx3)*( vx1       - vx3) - cu_sq);
		(dist.f[dirTW  ])[ktw  ] = c1o54* (drho + c3o1*(-vx1       + vx3) + c9o2*(-vx1       + vx3)*(-vx1       + vx3) - cu_sq);
		(dist.f[dirTN  ])[ktn  ] = c1o54* (drho + c3o1*(       vx2 + vx3) + c9o2*(       vx2 + vx3)*(       vx2 + vx3) - cu_sq);
		(dist.f[dirBS  ])[kbs  ] = c1o54* (drho + c3o1*(     - vx2 - vx3) + c9o2*(     - vx2 - vx3)*(     - vx2 - vx3) - cu_sq);
		(dist.f[dirBN  ])[kbn  ] = c1o54* (drho + c3o1*(       vx2 - vx3) + c9o2*(       vx2 - vx3)*(       vx2 - vx3) - cu_sq);
		(dist.f[dirTS  ])[kts  ] = c1o54* (drho + c3o1*(     - vx2 + vx3) + c9o2*(     - vx2 + vx3)*(     - vx2 + vx3) - cu_sq);
		(dist.f[dirTNE ])[ktne ] = c1o216*(drho + c3o1*( vx1 + vx2 + vx3) + c9o2*( vx1 + vx2 + vx3)*( vx1 + vx2 + vx3) - cu_sq);
		(dist.f[dirBSW ])[kbsw ] = c1o216*(drho + c3o1*(-vx1 - vx2 - vx3) + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cu_sq);
		(dist.f[dirBNE ])[kbne ] = c1o216*(drho + c3o1*( vx1 + vx2 - vx3) + c9o2*( vx1 + vx2 - vx3)*( vx1 + vx2 - vx3) - cu_sq);
		(dist.f[dirTSW ])[ktsw ] = c1o216*(drho + c3o1*(-vx1 - vx2 + vx3) + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cu_sq);
		(dist.f[dirTSE ])[ktse ] = c1o216*(drho + c3o1*( vx1 - vx2 + vx3) + c9o2*( vx1 - vx2 + vx3)*( vx1 - vx2 + vx3) - cu_sq);
		(dist.f[dirBNW ])[kbnw ] = c1o216*(drho + c3o1*(-vx1 + vx2 - vx3) + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cu_sq);
		(dist.f[dirBSE ])[kbse ] = c1o216*(drho + c3o1*( vx1 - vx2 - vx3) + c9o2*( vx1 - vx2 - vx3)*( vx1 - vx2 - vx3) - cu_sq);
		(dist.f[dirTNW ])[ktnw ] = c1o216*(drho + c3o1*(-vx1 + vx2 + vx3) + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cu_sq);
	}
}

