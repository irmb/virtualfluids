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
//! \file InitAdvectionDiffusion.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;


extern "C" __global__ void InitAD27(
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	uint* typeOfGridNode,
	real* concentration,
	real* velocityX,
	real* velocityY,
	real* velocityZ,
	uint size_Mat,
	real* distributionsAD,
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
			distAD.f[dirBSW ] = &distributionsAD[dirTNE *size_Mat];
			distAD.f[dirBNE ] = &distributionsAD[dirTSW *size_Mat];
			distAD.f[dirBNW ] = &distributionsAD[dirTSE *size_Mat];
			distAD.f[dirBSE ] = &distributionsAD[dirTNW *size_Mat];
			distAD.f[dirTSW ] = &distributionsAD[dirBNE *size_Mat];
			distAD.f[dirTNE ] = &distributionsAD[dirBSW *size_Mat];
			distAD.f[dirTNW ] = &distributionsAD[dirBSE *size_Mat];
			distAD.f[dirTSE ] = &distributionsAD[dirBNW *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////
		//! - Set local velocities and concetration
		//!
		real conc = concentration[k];
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

		(distAD.f[dirREST])[kzero] = c8o27  * conc * (c1o1 - cu_sq);
		(distAD.f[dirE   ])[ke   ] = c2o27  * conc * (c1o1 + c3o1 * ( vx1            ) + c9o2 * ( vx1            ) * ( vx1            ) - cu_sq);
		(distAD.f[dirW   ])[kw   ] = c2o27  * conc * (c1o1 + c3o1 * (-vx1            ) + c9o2 * (-vx1            ) * (-vx1            ) - cu_sq);
		(distAD.f[dirN   ])[kn   ] = c2o27  * conc * (c1o1 + c3o1 * (       vx2      ) + c9o2 * (       vx2      ) * (       vx2      ) - cu_sq);
		(distAD.f[dirS   ])[ks   ] = c2o27  * conc * (c1o1 + c3o1 * (     - vx2      ) + c9o2 * (     - vx2      ) * (     - vx2      ) - cu_sq);
		(distAD.f[dirT   ])[kt   ] = c2o27  * conc * (c1o1 + c3o1 * (             vx3) + c9o2 * (             vx3) * (             vx3) - cu_sq);
		(distAD.f[dirB   ])[kb   ] = c2o27  * conc * (c1o1 + c3o1 * (           - vx3) + c9o2 * (           - vx3) * (           - vx3) - cu_sq);
		(distAD.f[dirNE  ])[kne  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 + vx2      ) + c9o2 * ( vx1 + vx2      ) * ( vx1 + vx2      ) - cu_sq);
		(distAD.f[dirSW  ])[ksw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 - vx2      ) + c9o2 * (-vx1 - vx2      ) * (-vx1 - vx2      ) - cu_sq);
		(distAD.f[dirSE  ])[kse  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 - vx2      ) + c9o2 * ( vx1 - vx2      ) * ( vx1 - vx2      ) - cu_sq);
		(distAD.f[dirNW  ])[knw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 + vx2      ) + c9o2 * (-vx1 + vx2      ) * (-vx1 + vx2      ) - cu_sq);
		(distAD.f[dirTE  ])[kte  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       + vx3) + c9o2 * ( vx1       + vx3) * ( vx1       + vx3) - cu_sq);
		(distAD.f[dirBW  ])[kbw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       - vx3) + c9o2 * (-vx1       - vx3) * (-vx1       - vx3) - cu_sq);
		(distAD.f[dirBE  ])[kbe  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       - vx3) + c9o2 * ( vx1       - vx3) * ( vx1       - vx3) - cu_sq);
		(distAD.f[dirTW  ])[ktw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       + vx3) + c9o2 * (-vx1       + vx3) * (-vx1       + vx3) - cu_sq);
		(distAD.f[dirTN  ])[ktn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 + vx3) + c9o2 * (       vx2 + vx3) * (       vx2 + vx3) - cu_sq);
		(distAD.f[dirBS  ])[kbs  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 - vx3) + c9o2 * (     - vx2 - vx3) * (     - vx2 - vx3) - cu_sq);
		(distAD.f[dirBN  ])[kbn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 - vx3) + c9o2 * (       vx2 - vx3) * (       vx2 - vx3) - cu_sq);
		(distAD.f[dirTS  ])[kts  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 + vx3) + c9o2 * (     - vx2 + vx3) * (     - vx2 + vx3) - cu_sq);
		(distAD.f[dirTNE ])[ktne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 + vx3) + c9o2 * ( vx1 + vx2 + vx3) * ( vx1 + vx2 + vx3) - cu_sq);
		(distAD.f[dirBSW ])[kbsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[dirBNE ])[kbne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 - vx3) + c9o2 * ( vx1 + vx2 - vx3) * ( vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[dirTSW ])[ktsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[dirTSE ])[ktse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 + vx3) + c9o2 * ( vx1 - vx2 + vx3) * ( vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[dirBNW ])[kbnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[dirBSE ])[kbse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 - vx3) + c9o2 * ( vx1 - vx2 - vx3) * ( vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[dirTNW ])[ktnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
	}
}




















// DEPRECATED (2022)

// ////////////////////////////////////////////////////////////////////////////////
// extern "C" __global__ void InitAD27(unsigned int* neighborX,
//                                        unsigned int* neighborY,
//                                        unsigned int* neighborZ,
//                                        unsigned int* geoD,
//                                        real* Conc,
//                                        real* ux,
//                                        real* uy,
//                                        real* uz,
//                                        unsigned int size_Mat,
//                                        real* DD27,
//                                        bool EvenOrOdd)
// {
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
//       ////////////////////////////////////////////////////////////////////////////////
//       unsigned int BC;
//       BC        =   geoD[k];

//       if( BC != GEO_SOLID && BC != GEO_VOID)
//       {
//          Distributions27 D27;
//          if (EvenOrOdd==true)
//          {
//             D27.f[dirE   ] = &DD27[dirE   *size_Mat];
//             D27.f[dirW   ] = &DD27[dirW   *size_Mat];
//             D27.f[dirN   ] = &DD27[dirN   *size_Mat];
//             D27.f[dirS   ] = &DD27[dirS   *size_Mat];
//             D27.f[dirT   ] = &DD27[dirT   *size_Mat];
//             D27.f[dirB   ] = &DD27[dirB   *size_Mat];
//             D27.f[dirNE  ] = &DD27[dirNE  *size_Mat];
//             D27.f[dirSW  ] = &DD27[dirSW  *size_Mat];
//             D27.f[dirSE  ] = &DD27[dirSE  *size_Mat];
//             D27.f[dirNW  ] = &DD27[dirNW  *size_Mat];
//             D27.f[dirTE  ] = &DD27[dirTE  *size_Mat];
//             D27.f[dirBW  ] = &DD27[dirBW  *size_Mat];
//             D27.f[dirBE  ] = &DD27[dirBE  *size_Mat];
//             D27.f[dirTW  ] = &DD27[dirTW  *size_Mat];
//             D27.f[dirTN  ] = &DD27[dirTN  *size_Mat];
//             D27.f[dirBS  ] = &DD27[dirBS  *size_Mat];
//             D27.f[dirBN  ] = &DD27[dirBN  *size_Mat];
//             D27.f[dirTS  ] = &DD27[dirTS  *size_Mat];
//             D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
//             D27.f[dirTNE ] = &DD27[dirTNE *size_Mat];
//             D27.f[dirTSW ] = &DD27[dirTSW *size_Mat];
//             D27.f[dirTSE ] = &DD27[dirTSE *size_Mat];
//             D27.f[dirTNW ] = &DD27[dirTNW *size_Mat];
//             D27.f[dirBNE ] = &DD27[dirBNE *size_Mat];
//             D27.f[dirBSW ] = &DD27[dirBSW *size_Mat];
//             D27.f[dirBSE ] = &DD27[dirBSE *size_Mat];
//             D27.f[dirBNW ] = &DD27[dirBNW *size_Mat];
//          }
//          else
//          {
//             D27.f[dirW   ] = &DD27[dirE   *size_Mat];
//             D27.f[dirE   ] = &DD27[dirW   *size_Mat];
//             D27.f[dirS   ] = &DD27[dirN   *size_Mat];
//             D27.f[dirN   ] = &DD27[dirS   *size_Mat];
//             D27.f[dirB   ] = &DD27[dirT   *size_Mat];
//             D27.f[dirT   ] = &DD27[dirB   *size_Mat];
//             D27.f[dirSW  ] = &DD27[dirNE  *size_Mat];
//             D27.f[dirNE  ] = &DD27[dirSW  *size_Mat];
//             D27.f[dirNW  ] = &DD27[dirSE  *size_Mat];
//             D27.f[dirSE  ] = &DD27[dirNW  *size_Mat];
//             D27.f[dirBW  ] = &DD27[dirTE  *size_Mat];
//             D27.f[dirTE  ] = &DD27[dirBW  *size_Mat];
//             D27.f[dirTW  ] = &DD27[dirBE  *size_Mat];
//             D27.f[dirBE  ] = &DD27[dirTW  *size_Mat];
//             D27.f[dirBS  ] = &DD27[dirTN  *size_Mat];
//             D27.f[dirTN  ] = &DD27[dirBS  *size_Mat];
//             D27.f[dirTS  ] = &DD27[dirBN  *size_Mat];
//             D27.f[dirBN  ] = &DD27[dirTS  *size_Mat];
//             D27.f[dirZERO] = &DD27[dirZERO*size_Mat];
//             D27.f[dirBSW ] = &DD27[dirTNE *size_Mat];
//             D27.f[dirBNE ] = &DD27[dirTSW *size_Mat];
//             D27.f[dirBNW ] = &DD27[dirTSE *size_Mat];
//             D27.f[dirBSE ] = &DD27[dirTNW *size_Mat];
//             D27.f[dirTSW ] = &DD27[dirBNE *size_Mat];
//             D27.f[dirTNE ] = &DD27[dirBSW *size_Mat];
//             D27.f[dirTNW ] = &DD27[dirBSE *size_Mat];
//             D27.f[dirTSE ] = &DD27[dirBNW *size_Mat];
//          }
//          //////////////////////////////////////////////////////////////////////////
//          real ConcD = Conc[k];
//          real   vx1 = ux[k];
//          real   vx2 = uy[k];
//          real   vx3 = uz[k];
//          //real lambdaD     = -three + sqrt(three);
//          //real Diffusivity = c1o20;
//          //real Lam         = -(c1o2+one/lambdaD);
//          //real nue_d       = Lam/three;
//          //real ae          = Diffusivity/nue_d - one;
//          //real ux_sq       = vx1 * vx1;
//          //real uy_sq       = vx2 * vx2;
//          //real uz_sq       = vx3 * vx3;
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //D3Q7
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //index
//          //unsigned int kzero= k;
//          //unsigned int ke   = k;
//          //unsigned int kw   = neighborX[k];
//          //unsigned int kn   = k;
//          //unsigned int ks   = neighborY[k];
//          //unsigned int kt   = k;
//          //unsigned int kb   = neighborZ[k];
//          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //(D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-three))-(ux_sq+uy_sq+uz_sq));
//          //(D7.f[1])[ke   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)+vx1*c1o2);
//          //(D7.f[2])[kw   ] = ConcD*(c1o6*(ae+one)+c1o2*(ux_sq)-vx1*c1o2);
//          //(D7.f[3])[kn   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)+vx2*c1o2);
//          //(D7.f[4])[ks   ] = ConcD*(c1o6*(ae+one)+c1o2*(uy_sq)-vx2*c1o2);
//          //(D7.f[5])[kt   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)+vx3*c1o2);
//          //(D7.f[6])[kb   ] = ConcD*(c1o6*(ae+one)+c1o2*(uz_sq)-vx3*c1o2);
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //D3Q27
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //index
//          unsigned int kzero= k;
//          unsigned int ke   = k;
//          unsigned int kw   = neighborX[k];
//          unsigned int kn   = k;
//          unsigned int ks   = neighborY[k];
//          unsigned int kt   = k;
//          unsigned int kb   = neighborZ[k];
//          unsigned int ksw  = neighborY[kw];
//          unsigned int kne  = k;
//          unsigned int kse  = ks;
//          unsigned int knw  = kw;
//          unsigned int kbw  = neighborZ[kw];
//          unsigned int kte  = k;
//          unsigned int kbe  = kb;
//          unsigned int ktw  = kw;
//          unsigned int kbs  = neighborZ[ks];
//          unsigned int ktn  = k;
//          unsigned int kbn  = kb;
//          unsigned int kts  = ks;
//          unsigned int ktse = ks;
//          unsigned int kbnw = kbw;
//          unsigned int ktnw = kw;
//          unsigned int kbse = kbs;
//          unsigned int ktsw = ksw;
//          unsigned int kbne = kb;
//          unsigned int ktne = k;
//          unsigned int kbsw = neighborZ[ksw];
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

//          (D27.f[dirZERO])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
//          (D27.f[dirE   ])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
//          (D27.f[dirW   ])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
//          (D27.f[dirN   ])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
//          (D27.f[dirS   ])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
//          (D27.f[dirT   ])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
//          (D27.f[dirB   ])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
//          (D27.f[dirNE  ])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
//          (D27.f[dirSW  ])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
//          (D27.f[dirSE  ])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
//          (D27.f[dirNW  ])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
//          (D27.f[dirTE  ])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
//          (D27.f[dirBW  ])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
//          (D27.f[dirBE  ])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
//          (D27.f[dirTW  ])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
//          (D27.f[dirTN  ])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
//          (D27.f[dirBS  ])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
//          (D27.f[dirBN  ])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
//          (D27.f[dirTS  ])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
//          (D27.f[dirTNE ])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
//          (D27.f[dirBSW ])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
//          (D27.f[dirBNE ])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
//          (D27.f[dirTSW ])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
//          (D27.f[dirTSE ])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
//          (D27.f[dirBNW ])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
//          (D27.f[dirBSE ])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
//          (D27.f[dirTNW ])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       }
//    }
// }


















////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void InitAD7( unsigned int* neighborX,
                                    unsigned int* neighborY,
                                    unsigned int* neighborZ,
                                    unsigned int* geoD,
                                    real* Conc,
                                    real* ux,
                                    real* uy,
                                    real* uz,
                                    unsigned int size_Mat,
                                    real* DD7,
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

   if(k<size_Mat)
   {
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions7 D7;
         if (EvenOrOdd==true)
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
         //////////////////////////////////////////////////////////////////////////
         real ConcD = Conc[k];
         real   vx1 = ux[k];
         real   vx2 = uy[k];
         real   vx3 = uz[k];
         real lambdaD     = -c3o1 + sqrt(c3o1);
         real Diffusivity = c1o20;
         real Lam         = -(c1o2+c1o1/lambdaD);
         real nue_d       = Lam/c3o1;
         real ae          = Diffusivity/nue_d - c1o1;
         real ux_sq       = vx1 * vx1;
         real uy_sq       = vx2 * vx2;
         real uz_sq       = vx3 * vx3;
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
         unsigned int kzero= k;
         unsigned int ke   = k;
         unsigned int kw   = neighborX[k];
         unsigned int kn   = k;
         unsigned int ks   = neighborY[k];
         unsigned int kt   = k;
         unsigned int kb   = neighborZ[k];
         //////////////////////////////////////////////////////////////////////////

         (D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-c3o1))-(ux_sq+uy_sq+uz_sq));
         (D7.f[1])[ke   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)+vx1*c1o2);
         (D7.f[2])[kw   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)-vx1*c1o2);
         (D7.f[3])[kn   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)+vx2*c1o2);
         (D7.f[4])[ks   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)-vx2*c1o2);
         (D7.f[5])[kt   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)+vx3*c1o2);
         (D7.f[6])[kb   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)-vx3*c1o2);
      }
   }
}