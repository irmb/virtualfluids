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
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

__global__ void InitAD27(
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
	//! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
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
			distAD.f[REST] = &distributionsAD[REST*size_Mat];
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
			distAD.f[REST] = &distributionsAD[REST*size_Mat];
			distAD.f[BSW ] = &distributionsAD[TNE *size_Mat];
			distAD.f[BNE ] = &distributionsAD[TSW *size_Mat];
			distAD.f[BNW ] = &distributionsAD[TSE *size_Mat];
			distAD.f[BSE ] = &distributionsAD[TNW *size_Mat];
			distAD.f[TSW ] = &distributionsAD[BNE *size_Mat];
			distAD.f[TNE ] = &distributionsAD[BSW *size_Mat];
			distAD.f[TNW ] = &distributionsAD[BSE *size_Mat];
			distAD.f[TSE ] = &distributionsAD[BNW *size_Mat];
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

		(distAD.f[REST])[kzero] = c8o27  * conc * (c1o1 - cu_sq);
		(distAD.f[E   ])[ke   ] = c2o27  * conc * (c1o1 + c3o1 * ( vx1            ) + c9o2 * ( vx1            ) * ( vx1            ) - cu_sq);
		(distAD.f[W   ])[kw   ] = c2o27  * conc * (c1o1 + c3o1 * (-vx1            ) + c9o2 * (-vx1            ) * (-vx1            ) - cu_sq);
		(distAD.f[N   ])[kn   ] = c2o27  * conc * (c1o1 + c3o1 * (       vx2      ) + c9o2 * (       vx2      ) * (       vx2      ) - cu_sq);
		(distAD.f[S   ])[ks   ] = c2o27  * conc * (c1o1 + c3o1 * (     - vx2      ) + c9o2 * (     - vx2      ) * (     - vx2      ) - cu_sq);
		(distAD.f[T   ])[kt   ] = c2o27  * conc * (c1o1 + c3o1 * (             vx3) + c9o2 * (             vx3) * (             vx3) - cu_sq);
		(distAD.f[B   ])[kb   ] = c2o27  * conc * (c1o1 + c3o1 * (           - vx3) + c9o2 * (           - vx3) * (           - vx3) - cu_sq);
		(distAD.f[NE  ])[kne  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 + vx2      ) + c9o2 * ( vx1 + vx2      ) * ( vx1 + vx2      ) - cu_sq);
		(distAD.f[SW  ])[ksw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 - vx2      ) + c9o2 * (-vx1 - vx2      ) * (-vx1 - vx2      ) - cu_sq);
		(distAD.f[SE  ])[kse  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 - vx2      ) + c9o2 * ( vx1 - vx2      ) * ( vx1 - vx2      ) - cu_sq);
		(distAD.f[NW  ])[knw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 + vx2      ) + c9o2 * (-vx1 + vx2      ) * (-vx1 + vx2      ) - cu_sq);
		(distAD.f[TE  ])[kte  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       + vx3) + c9o2 * ( vx1       + vx3) * ( vx1       + vx3) - cu_sq);
		(distAD.f[BW  ])[kbw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       - vx3) + c9o2 * (-vx1       - vx3) * (-vx1       - vx3) - cu_sq);
		(distAD.f[BE  ])[kbe  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       - vx3) + c9o2 * ( vx1       - vx3) * ( vx1       - vx3) - cu_sq);
		(distAD.f[TW  ])[ktw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       + vx3) + c9o2 * (-vx1       + vx3) * (-vx1       + vx3) - cu_sq);
		(distAD.f[TN  ])[ktn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 + vx3) + c9o2 * (       vx2 + vx3) * (       vx2 + vx3) - cu_sq);
		(distAD.f[BS  ])[kbs  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 - vx3) + c9o2 * (     - vx2 - vx3) * (     - vx2 - vx3) - cu_sq);
		(distAD.f[BN  ])[kbn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 - vx3) + c9o2 * (       vx2 - vx3) * (       vx2 - vx3) - cu_sq);
		(distAD.f[TS  ])[kts  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 + vx3) + c9o2 * (     - vx2 + vx3) * (     - vx2 + vx3) - cu_sq);
		(distAD.f[TNE ])[ktne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 + vx3) + c9o2 * ( vx1 + vx2 + vx3) * ( vx1 + vx2 + vx3) - cu_sq);
		(distAD.f[BSW ])[kbsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[BNE ])[kbne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 - vx3) + c9o2 * ( vx1 + vx2 - vx3) * ( vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[TSW ])[ktsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[TSE ])[ktse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 + vx3) + c9o2 * ( vx1 - vx2 + vx3) * ( vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[BNW ])[kbnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[BSE ])[kbse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 - vx3) + c9o2 * ( vx1 - vx2 - vx3) * ( vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[TNW ])[ktnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
	}
}




















// DEPRECATED (2022)

// ////////////////////////////////////////////////////////////////////////////////
// __global__ void InitAD27(unsigned int* neighborX,
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
//             D27.f[E   ] = &DD27[E   *size_Mat];
//             D27.f[W   ] = &DD27[W   *size_Mat];
//             D27.f[N   ] = &DD27[N   *size_Mat];
//             D27.f[S   ] = &DD27[S   *size_Mat];
//             D27.f[T   ] = &DD27[T   *size_Mat];
//             D27.f[B   ] = &DD27[B   *size_Mat];
//             D27.f[NE  ] = &DD27[NE  *size_Mat];
//             D27.f[SW  ] = &DD27[SW  *size_Mat];
//             D27.f[SE  ] = &DD27[SE  *size_Mat];
//             D27.f[NW  ] = &DD27[NW  *size_Mat];
//             D27.f[TE  ] = &DD27[TE  *size_Mat];
//             D27.f[BW  ] = &DD27[BW  *size_Mat];
//             D27.f[BE  ] = &DD27[BE  *size_Mat];
//             D27.f[TW  ] = &DD27[TW  *size_Mat];
//             D27.f[TN  ] = &DD27[TN  *size_Mat];
//             D27.f[BS  ] = &DD27[BS  *size_Mat];
//             D27.f[BN  ] = &DD27[BN  *size_Mat];
//             D27.f[TS  ] = &DD27[TS  *size_Mat];
//             D27.f[REST] = &DD27[REST*size_Mat];
//             D27.f[TNE ] = &DD27[TNE *size_Mat];
//             D27.f[TSW ] = &DD27[TSW *size_Mat];
//             D27.f[TSE ] = &DD27[TSE *size_Mat];
//             D27.f[TNW ] = &DD27[TNW *size_Mat];
//             D27.f[BNE ] = &DD27[BNE *size_Mat];
//             D27.f[BSW ] = &DD27[BSW *size_Mat];
//             D27.f[BSE ] = &DD27[BSE *size_Mat];
//             D27.f[BNW ] = &DD27[BNW *size_Mat];
//          }
//          else
//          {
//             D27.f[W   ] = &DD27[E   *size_Mat];
//             D27.f[E   ] = &DD27[W   *size_Mat];
//             D27.f[S   ] = &DD27[N   *size_Mat];
//             D27.f[N   ] = &DD27[S   *size_Mat];
//             D27.f[B   ] = &DD27[T   *size_Mat];
//             D27.f[T   ] = &DD27[B   *size_Mat];
//             D27.f[SW  ] = &DD27[NE  *size_Mat];
//             D27.f[NE  ] = &DD27[SW  *size_Mat];
//             D27.f[NW  ] = &DD27[SE  *size_Mat];
//             D27.f[SE  ] = &DD27[NW  *size_Mat];
//             D27.f[BW  ] = &DD27[TE  *size_Mat];
//             D27.f[TE  ] = &DD27[BW  *size_Mat];
//             D27.f[TW  ] = &DD27[BE  *size_Mat];
//             D27.f[BE  ] = &DD27[TW  *size_Mat];
//             D27.f[BS  ] = &DD27[TN  *size_Mat];
//             D27.f[TN  ] = &DD27[BS  *size_Mat];
//             D27.f[TS  ] = &DD27[BN  *size_Mat];
//             D27.f[BN  ] = &DD27[TS  *size_Mat];
//             D27.f[REST] = &DD27[REST*size_Mat];
//             D27.f[BSW ] = &DD27[TNE *size_Mat];
//             D27.f[BNE ] = &DD27[TSW *size_Mat];
//             D27.f[BNW ] = &DD27[TSE *size_Mat];
//             D27.f[BSE ] = &DD27[TNW *size_Mat];
//             D27.f[TSW ] = &DD27[BNE *size_Mat];
//             D27.f[TNE ] = &DD27[BSW *size_Mat];
//             D27.f[TNW ] = &DD27[BSE *size_Mat];
//             D27.f[TSE ] = &DD27[BNW *size_Mat];
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

//          (D27.f[REST])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
//          (D27.f[E   ])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
//          (D27.f[W   ])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
//          (D27.f[N   ])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
//          (D27.f[S   ])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
//          (D27.f[T   ])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
//          (D27.f[B   ])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
//          (D27.f[NE  ])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
//          (D27.f[SW  ])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
//          (D27.f[SE  ])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
//          (D27.f[NW  ])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
//          (D27.f[TE  ])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
//          (D27.f[BW  ])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
//          (D27.f[BE  ])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
//          (D27.f[TW  ])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
//          (D27.f[TN  ])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
//          (D27.f[BS  ])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
//          (D27.f[BN  ])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
//          (D27.f[TS  ])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
//          (D27.f[TNE ])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
//          (D27.f[BSW ])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
//          (D27.f[BNE ])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
//          (D27.f[TSW ])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
//          (D27.f[TSE ])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
//          (D27.f[BNW ])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
//          (D27.f[BSE ])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
//          (D27.f[TNW ])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
//          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//       }
//    }
// }


















////////////////////////////////////////////////////////////////////////////////
__global__ void InitAD7( unsigned int* neighborX,
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