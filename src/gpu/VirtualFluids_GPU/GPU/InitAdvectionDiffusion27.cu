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
			distAD.f[DIR_P00   ] = &distributionsAD[DIR_P00   *size_Mat];
			distAD.f[DIR_M00   ] = &distributionsAD[DIR_M00   *size_Mat];
			distAD.f[DIR_0P0   ] = &distributionsAD[DIR_0P0   *size_Mat];
			distAD.f[DIR_0M0   ] = &distributionsAD[DIR_0M0   *size_Mat];
			distAD.f[DIR_00P   ] = &distributionsAD[DIR_00P   *size_Mat];
			distAD.f[DIR_00M   ] = &distributionsAD[DIR_00M   *size_Mat];
			distAD.f[DIR_PP0  ] = &distributionsAD[DIR_PP0  *size_Mat];
			distAD.f[DIR_MM0  ] = &distributionsAD[DIR_MM0  *size_Mat];
			distAD.f[DIR_PM0  ] = &distributionsAD[DIR_PM0  *size_Mat];
			distAD.f[DIR_MP0  ] = &distributionsAD[DIR_MP0  *size_Mat];
			distAD.f[DIR_P0P  ] = &distributionsAD[DIR_P0P  *size_Mat];
			distAD.f[DIR_M0M  ] = &distributionsAD[DIR_M0M  *size_Mat];
			distAD.f[DIR_P0M  ] = &distributionsAD[DIR_P0M  *size_Mat];
			distAD.f[DIR_M0P  ] = &distributionsAD[DIR_M0P  *size_Mat];
			distAD.f[DIR_0PP  ] = &distributionsAD[DIR_0PP  *size_Mat];
			distAD.f[DIR_0MM  ] = &distributionsAD[DIR_0MM  *size_Mat];
			distAD.f[DIR_0PM  ] = &distributionsAD[DIR_0PM  *size_Mat];
			distAD.f[DIR_0MP  ] = &distributionsAD[DIR_0MP  *size_Mat];
			distAD.f[DIR_000] = &distributionsAD[DIR_000*size_Mat];
			distAD.f[DIR_PPP ] = &distributionsAD[DIR_PPP *size_Mat];
			distAD.f[DIR_MMP ] = &distributionsAD[DIR_MMP *size_Mat];
			distAD.f[DIR_PMP ] = &distributionsAD[DIR_PMP *size_Mat];
			distAD.f[DIR_MPP ] = &distributionsAD[DIR_MPP *size_Mat];
			distAD.f[DIR_PPM ] = &distributionsAD[DIR_PPM *size_Mat];
			distAD.f[DIR_MMM ] = &distributionsAD[DIR_MMM *size_Mat];
			distAD.f[DIR_PMM ] = &distributionsAD[DIR_PMM *size_Mat];
			distAD.f[DIR_MPM ] = &distributionsAD[DIR_MPM *size_Mat];
		}
		else
		{
			distAD.f[DIR_M00   ] = &distributionsAD[DIR_P00   *size_Mat];
			distAD.f[DIR_P00   ] = &distributionsAD[DIR_M00   *size_Mat];
			distAD.f[DIR_0M0   ] = &distributionsAD[DIR_0P0   *size_Mat];
			distAD.f[DIR_0P0   ] = &distributionsAD[DIR_0M0   *size_Mat];
			distAD.f[DIR_00M   ] = &distributionsAD[DIR_00P   *size_Mat];
			distAD.f[DIR_00P   ] = &distributionsAD[DIR_00M   *size_Mat];
			distAD.f[DIR_MM0  ] = &distributionsAD[DIR_PP0  *size_Mat];
			distAD.f[DIR_PP0  ] = &distributionsAD[DIR_MM0  *size_Mat];
			distAD.f[DIR_MP0  ] = &distributionsAD[DIR_PM0  *size_Mat];
			distAD.f[DIR_PM0  ] = &distributionsAD[DIR_MP0  *size_Mat];
			distAD.f[DIR_M0M  ] = &distributionsAD[DIR_P0P  *size_Mat];
			distAD.f[DIR_P0P  ] = &distributionsAD[DIR_M0M  *size_Mat];
			distAD.f[DIR_M0P  ] = &distributionsAD[DIR_P0M  *size_Mat];
			distAD.f[DIR_P0M  ] = &distributionsAD[DIR_M0P  *size_Mat];
			distAD.f[DIR_0MM  ] = &distributionsAD[DIR_0PP  *size_Mat];
			distAD.f[DIR_0PP  ] = &distributionsAD[DIR_0MM  *size_Mat];
			distAD.f[DIR_0MP  ] = &distributionsAD[DIR_0PM  *size_Mat];
			distAD.f[DIR_0PM  ] = &distributionsAD[DIR_0MP  *size_Mat];
			distAD.f[DIR_000] = &distributionsAD[DIR_000*size_Mat];
			distAD.f[DIR_MMM ] = &distributionsAD[DIR_PPP *size_Mat];
			distAD.f[DIR_PPM ] = &distributionsAD[DIR_MMP *size_Mat];
			distAD.f[DIR_MPM ] = &distributionsAD[DIR_PMP *size_Mat];
			distAD.f[DIR_PMM ] = &distributionsAD[DIR_MPP *size_Mat];
			distAD.f[DIR_MMP ] = &distributionsAD[DIR_PPM *size_Mat];
			distAD.f[DIR_PPP ] = &distributionsAD[DIR_MMM *size_Mat];
			distAD.f[DIR_MPP ] = &distributionsAD[DIR_PMM *size_Mat];
			distAD.f[DIR_PMP ] = &distributionsAD[DIR_MPM *size_Mat];
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

		(distAD.f[DIR_000])[kzero] = c8o27  * conc * (c1o1 - cu_sq);
		(distAD.f[DIR_P00   ])[ke   ] = c2o27  * conc * (c1o1 + c3o1 * ( vx1            ) + c9o2 * ( vx1            ) * ( vx1            ) - cu_sq);
		(distAD.f[DIR_M00   ])[kw   ] = c2o27  * conc * (c1o1 + c3o1 * (-vx1            ) + c9o2 * (-vx1            ) * (-vx1            ) - cu_sq);
		(distAD.f[DIR_0P0   ])[kn   ] = c2o27  * conc * (c1o1 + c3o1 * (       vx2      ) + c9o2 * (       vx2      ) * (       vx2      ) - cu_sq);
		(distAD.f[DIR_0M0   ])[ks   ] = c2o27  * conc * (c1o1 + c3o1 * (     - vx2      ) + c9o2 * (     - vx2      ) * (     - vx2      ) - cu_sq);
		(distAD.f[DIR_00P   ])[kt   ] = c2o27  * conc * (c1o1 + c3o1 * (             vx3) + c9o2 * (             vx3) * (             vx3) - cu_sq);
		(distAD.f[DIR_00M   ])[kb   ] = c2o27  * conc * (c1o1 + c3o1 * (           - vx3) + c9o2 * (           - vx3) * (           - vx3) - cu_sq);
		(distAD.f[DIR_PP0  ])[kne  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 + vx2      ) + c9o2 * ( vx1 + vx2      ) * ( vx1 + vx2      ) - cu_sq);
		(distAD.f[DIR_MM0  ])[ksw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 - vx2      ) + c9o2 * (-vx1 - vx2      ) * (-vx1 - vx2      ) - cu_sq);
		(distAD.f[DIR_PM0  ])[kse  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 - vx2      ) + c9o2 * ( vx1 - vx2      ) * ( vx1 - vx2      ) - cu_sq);
		(distAD.f[DIR_MP0  ])[knw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 + vx2      ) + c9o2 * (-vx1 + vx2      ) * (-vx1 + vx2      ) - cu_sq);
		(distAD.f[DIR_P0P  ])[kte  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       + vx3) + c9o2 * ( vx1       + vx3) * ( vx1       + vx3) - cu_sq);
		(distAD.f[DIR_M0M  ])[kbw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       - vx3) + c9o2 * (-vx1       - vx3) * (-vx1       - vx3) - cu_sq);
		(distAD.f[DIR_P0M  ])[kbe  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       - vx3) + c9o2 * ( vx1       - vx3) * ( vx1       - vx3) - cu_sq);
		(distAD.f[DIR_M0P  ])[ktw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       + vx3) + c9o2 * (-vx1       + vx3) * (-vx1       + vx3) - cu_sq);
		(distAD.f[DIR_0PP  ])[ktn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 + vx3) + c9o2 * (       vx2 + vx3) * (       vx2 + vx3) - cu_sq);
		(distAD.f[DIR_0MM  ])[kbs  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 - vx3) + c9o2 * (     - vx2 - vx3) * (     - vx2 - vx3) - cu_sq);
		(distAD.f[DIR_0PM  ])[kbn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 - vx3) + c9o2 * (       vx2 - vx3) * (       vx2 - vx3) - cu_sq);
		(distAD.f[DIR_0MP  ])[kts  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 + vx3) + c9o2 * (     - vx2 + vx3) * (     - vx2 + vx3) - cu_sq);
		(distAD.f[DIR_PPP ])[ktne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 + vx3) + c9o2 * ( vx1 + vx2 + vx3) * ( vx1 + vx2 + vx3) - cu_sq);
		(distAD.f[DIR_MMM ])[kbsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[DIR_PPM ])[kbne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 - vx3) + c9o2 * ( vx1 + vx2 - vx3) * ( vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[DIR_MMP ])[ktsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[DIR_PMP ])[ktse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 + vx3) + c9o2 * ( vx1 - vx2 + vx3) * ( vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[DIR_MPM ])[kbnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[DIR_PMM ])[kbse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 - vx3) + c9o2 * ( vx1 - vx2 - vx3) * ( vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[DIR_MPP ])[ktnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
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
//             D27.f[DIR_P00   ] = &DD27[DIR_P00   *size_Mat];
//             D27.f[DIR_M00   ] = &DD27[DIR_M00   *size_Mat];
//             D27.f[DIR_0P0   ] = &DD27[DIR_0P0   *size_Mat];
//             D27.f[DIR_0M0   ] = &DD27[DIR_0M0   *size_Mat];
//             D27.f[DIR_00P   ] = &DD27[DIR_00P   *size_Mat];
//             D27.f[DIR_00M   ] = &DD27[DIR_00M   *size_Mat];
//             D27.f[DIR_PP0  ] = &DD27[DIR_PP0  *size_Mat];
//             D27.f[DIR_MM0  ] = &DD27[DIR_MM0  *size_Mat];
//             D27.f[DIR_PM0  ] = &DD27[DIR_PM0  *size_Mat];
//             D27.f[DIR_MP0  ] = &DD27[DIR_MP0  *size_Mat];
//             D27.f[DIR_P0P  ] = &DD27[DIR_P0P  *size_Mat];
//             D27.f[DIR_M0M  ] = &DD27[DIR_M0M  *size_Mat];
//             D27.f[DIR_P0M  ] = &DD27[DIR_P0M  *size_Mat];
//             D27.f[DIR_M0P  ] = &DD27[DIR_M0P  *size_Mat];
//             D27.f[DIR_0PP  ] = &DD27[DIR_0PP  *size_Mat];
//             D27.f[DIR_0MM  ] = &DD27[DIR_0MM  *size_Mat];
//             D27.f[DIR_0PM  ] = &DD27[DIR_0PM  *size_Mat];
//             D27.f[DIR_0MP  ] = &DD27[DIR_0MP  *size_Mat];
//             D27.f[DIR_000] = &DD27[DIR_000*size_Mat];
//             D27.f[DIR_PPP ] = &DD27[DIR_PPP *size_Mat];
//             D27.f[DIR_MMP ] = &DD27[DIR_MMP *size_Mat];
//             D27.f[DIR_PMP ] = &DD27[DIR_PMP *size_Mat];
//             D27.f[DIR_MPP ] = &DD27[DIR_MPP *size_Mat];
//             D27.f[DIR_PPM ] = &DD27[DIR_PPM *size_Mat];
//             D27.f[DIR_MMM ] = &DD27[DIR_MMM *size_Mat];
//             D27.f[DIR_PMM ] = &DD27[DIR_PMM *size_Mat];
//             D27.f[DIR_MPM ] = &DD27[DIR_MPM *size_Mat];
//          }
//          else
//          {
//             D27.f[DIR_M00   ] = &DD27[DIR_P00   *size_Mat];
//             D27.f[DIR_P00   ] = &DD27[DIR_M00   *size_Mat];
//             D27.f[DIR_0M0   ] = &DD27[DIR_0P0   *size_Mat];
//             D27.f[DIR_0P0   ] = &DD27[DIR_0M0   *size_Mat];
//             D27.f[DIR_00M   ] = &DD27[DIR_00P   *size_Mat];
//             D27.f[DIR_00P   ] = &DD27[DIR_00M   *size_Mat];
//             D27.f[DIR_MM0  ] = &DD27[DIR_PP0  *size_Mat];
//             D27.f[DIR_PP0  ] = &DD27[DIR_MM0  *size_Mat];
//             D27.f[DIR_MP0  ] = &DD27[DIR_PM0  *size_Mat];
//             D27.f[DIR_PM0  ] = &DD27[DIR_MP0  *size_Mat];
//             D27.f[DIR_M0M  ] = &DD27[DIR_P0P  *size_Mat];
//             D27.f[DIR_P0P  ] = &DD27[DIR_M0M  *size_Mat];
//             D27.f[DIR_M0P  ] = &DD27[DIR_P0M  *size_Mat];
//             D27.f[DIR_P0M  ] = &DD27[DIR_M0P  *size_Mat];
//             D27.f[DIR_0MM  ] = &DD27[DIR_0PP  *size_Mat];
//             D27.f[DIR_0PP  ] = &DD27[DIR_0MM  *size_Mat];
//             D27.f[DIR_0MP  ] = &DD27[DIR_0PM  *size_Mat];
//             D27.f[DIR_0PM  ] = &DD27[DIR_0MP  *size_Mat];
//             D27.f[DIR_000] = &DD27[DIR_000*size_Mat];
//             D27.f[DIR_MMM ] = &DD27[DIR_PPP *size_Mat];
//             D27.f[DIR_PPM ] = &DD27[DIR_MMP *size_Mat];
//             D27.f[DIR_MPM ] = &DD27[DIR_PMP *size_Mat];
//             D27.f[DIR_PMM ] = &DD27[DIR_MPP *size_Mat];
//             D27.f[DIR_MMP ] = &DD27[DIR_PPM *size_Mat];
//             D27.f[DIR_PPP ] = &DD27[DIR_MMM *size_Mat];
//             D27.f[DIR_MPP ] = &DD27[DIR_PMM *size_Mat];
//             D27.f[DIR_PMP ] = &DD27[DIR_MPM *size_Mat];
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

//          (D27.f[DIR_000])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
//          (D27.f[DIR_P00   ])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
//          (D27.f[DIR_M00   ])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
//          (D27.f[DIR_0P0   ])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
//          (D27.f[DIR_0M0   ])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
//          (D27.f[DIR_00P   ])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
//          (D27.f[DIR_00M   ])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
//          (D27.f[DIR_PP0  ])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
//          (D27.f[DIR_MM0  ])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
//          (D27.f[DIR_PM0  ])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
//          (D27.f[DIR_MP0  ])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
//          (D27.f[DIR_P0P  ])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
//          (D27.f[DIR_M0M  ])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
//          (D27.f[DIR_P0M  ])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
//          (D27.f[DIR_M0P  ])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
//          (D27.f[DIR_0PP  ])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
//          (D27.f[DIR_0MM  ])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
//          (D27.f[DIR_0PM  ])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
//          (D27.f[DIR_0MP  ])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
//          (D27.f[DIR_PPP ])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
//          (D27.f[DIR_MMM ])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
//          (D27.f[DIR_PPM ])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
//          (D27.f[DIR_MMP ])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
//          (D27.f[DIR_PMP ])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
//          (D27.f[DIR_MPM ])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
//          (D27.f[DIR_PMM ])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
//          (D27.f[DIR_MPP ])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
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