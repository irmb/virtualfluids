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
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
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
	unsigned long long numberOfLBnodes,
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
			distAD.f[dMMM] = &distributionsAD[dPPP * numberOfLBnodes];
			distAD.f[dPPM] = &distributionsAD[dMMP * numberOfLBnodes];
			distAD.f[dMPM] = &distributionsAD[dPMP * numberOfLBnodes];
			distAD.f[dPMM] = &distributionsAD[dMPP * numberOfLBnodes];
			distAD.f[dMMP] = &distributionsAD[dPPM * numberOfLBnodes];
			distAD.f[dPPP] = &distributionsAD[dMMM * numberOfLBnodes];
			distAD.f[dMPP] = &distributionsAD[dPMM * numberOfLBnodes];
			distAD.f[dPMP] = &distributionsAD[dMPM * numberOfLBnodes];
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

		(distAD.f[d000])[kzero] = c8o27  * conc * (c1o1 - cu_sq);
		(distAD.f[dP00])[ke   ] = c2o27  * conc * (c1o1 + c3o1 * ( vx1            ) + c9o2 * ( vx1            ) * ( vx1            ) - cu_sq);
		(distAD.f[dM00])[kw   ] = c2o27  * conc * (c1o1 + c3o1 * (-vx1            ) + c9o2 * (-vx1            ) * (-vx1            ) - cu_sq);
		(distAD.f[d0P0])[kn   ] = c2o27  * conc * (c1o1 + c3o1 * (       vx2      ) + c9o2 * (       vx2      ) * (       vx2      ) - cu_sq);
		(distAD.f[d0M0])[ks   ] = c2o27  * conc * (c1o1 + c3o1 * (     - vx2      ) + c9o2 * (     - vx2      ) * (     - vx2      ) - cu_sq);
		(distAD.f[d00P])[kt   ] = c2o27  * conc * (c1o1 + c3o1 * (             vx3) + c9o2 * (             vx3) * (             vx3) - cu_sq);
		(distAD.f[d00M])[kb   ] = c2o27  * conc * (c1o1 + c3o1 * (           - vx3) + c9o2 * (           - vx3) * (           - vx3) - cu_sq);
		(distAD.f[dPP0])[kne  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 + vx2      ) + c9o2 * ( vx1 + vx2      ) * ( vx1 + vx2      ) - cu_sq);
		(distAD.f[dMM0])[ksw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 - vx2      ) + c9o2 * (-vx1 - vx2      ) * (-vx1 - vx2      ) - cu_sq);
		(distAD.f[dPM0])[kse  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1 - vx2      ) + c9o2 * ( vx1 - vx2      ) * ( vx1 - vx2      ) - cu_sq);
		(distAD.f[dMP0])[knw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1 + vx2      ) + c9o2 * (-vx1 + vx2      ) * (-vx1 + vx2      ) - cu_sq);
		(distAD.f[dP0P])[kte  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       + vx3) + c9o2 * ( vx1       + vx3) * ( vx1       + vx3) - cu_sq);
		(distAD.f[dM0M])[kbw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       - vx3) + c9o2 * (-vx1       - vx3) * (-vx1       - vx3) - cu_sq);
		(distAD.f[dP0M])[kbe  ] = c1o54  * conc * (c1o1 + c3o1 * ( vx1       - vx3) + c9o2 * ( vx1       - vx3) * ( vx1       - vx3) - cu_sq);
		(distAD.f[dM0P])[ktw  ] = c1o54  * conc * (c1o1 + c3o1 * (-vx1       + vx3) + c9o2 * (-vx1       + vx3) * (-vx1       + vx3) - cu_sq);
		(distAD.f[d0PP])[ktn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 + vx3) + c9o2 * (       vx2 + vx3) * (       vx2 + vx3) - cu_sq);
		(distAD.f[d0MM])[kbs  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 - vx3) + c9o2 * (     - vx2 - vx3) * (     - vx2 - vx3) - cu_sq);
		(distAD.f[d0PM])[kbn  ] = c1o54  * conc * (c1o1 + c3o1 * (       vx2 - vx3) + c9o2 * (       vx2 - vx3) * (       vx2 - vx3) - cu_sq);
		(distAD.f[d0MP])[kts  ] = c1o54  * conc * (c1o1 + c3o1 * (     - vx2 + vx3) + c9o2 * (     - vx2 + vx3) * (     - vx2 + vx3) - cu_sq);
		(distAD.f[dPPP])[ktne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 + vx3) + c9o2 * ( vx1 + vx2 + vx3) * ( vx1 + vx2 + vx3) - cu_sq);
		(distAD.f[dMMM])[kbsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[dPPM])[kbne ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 + vx2 - vx3) + c9o2 * ( vx1 + vx2 - vx3) * ( vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[dMMP])[ktsw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[dPMP])[ktse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 + vx3) + c9o2 * ( vx1 - vx2 + vx3) * ( vx1 - vx2 + vx3) - cu_sq);
		(distAD.f[dMPM])[kbnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
		(distAD.f[dPMM])[kbse ] = c1o216 * conc * (c1o1 + c3o1 * ( vx1 - vx2 - vx3) + c9o2 * ( vx1 - vx2 - vx3) * ( vx1 - vx2 - vx3) - cu_sq);
		(distAD.f[dMPP])[ktnw ] = c1o216 * conc * (c1o1 + c3o1 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
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
//             D27.f[dP00] = &DD27[dP00 * size_Mat];
//             D27.f[dM00] = &DD27[dM00 * size_Mat];
//             D27.f[d0P0] = &DD27[d0P0 * size_Mat];
//             D27.f[d0M0] = &DD27[d0M0 * size_Mat];
//             D27.f[d00P] = &DD27[d00P * size_Mat];
//             D27.f[d00M] = &DD27[d00M * size_Mat];
//             D27.f[dPP0] = &DD27[dPP0 * size_Mat];
//             D27.f[dMM0] = &DD27[dMM0 * size_Mat];
//             D27.f[dPM0] = &DD27[dPM0 * size_Mat];
//             D27.f[dMP0] = &DD27[dMP0 * size_Mat];
//             D27.f[dP0P] = &DD27[dP0P * size_Mat];
//             D27.f[dM0M] = &DD27[dM0M * size_Mat];
//             D27.f[dP0M] = &DD27[dP0M * size_Mat];
//             D27.f[dM0P] = &DD27[dM0P * size_Mat];
//             D27.f[d0PP] = &DD27[d0PP * size_Mat];
//             D27.f[d0MM] = &DD27[d0MM * size_Mat];
//             D27.f[d0PM] = &DD27[d0PM * size_Mat];
//             D27.f[d0MP] = &DD27[d0MP * size_Mat];
//             D27.f[d000] = &DD27[d000 * size_Mat];
//             D27.f[dPPP] = &DD27[dPPP * size_Mat];
//             D27.f[dMMP] = &DD27[dMMP * size_Mat];
//             D27.f[dPMP] = &DD27[dPMP * size_Mat];
//             D27.f[dMPP] = &DD27[dMPP * size_Mat];
//             D27.f[dPPM] = &DD27[dPPM * size_Mat];
//             D27.f[dMMM] = &DD27[dMMM * size_Mat];
//             D27.f[dPMM] = &DD27[dPMM * size_Mat];
//             D27.f[dMPM] = &DD27[dMPM * size_Mat];
//          }
//          else
//          {
//             D27.f[dM00] = &DD27[dP00 * size_Mat];
//             D27.f[dP00] = &DD27[dM00 * size_Mat];
//             D27.f[d0M0] = &DD27[d0P0 * size_Mat];
//             D27.f[d0P0] = &DD27[d0M0 * size_Mat];
//             D27.f[d00M] = &DD27[d00P * size_Mat];
//             D27.f[d00P] = &DD27[d00M * size_Mat];
//             D27.f[dMM0] = &DD27[dPP0 * size_Mat];
//             D27.f[dPP0] = &DD27[dMM0 * size_Mat];
//             D27.f[dMP0] = &DD27[dPM0 * size_Mat];
//             D27.f[dPM0] = &DD27[dMP0 * size_Mat];
//             D27.f[dM0M] = &DD27[dP0P * size_Mat];
//             D27.f[dP0P] = &DD27[dM0M * size_Mat];
//             D27.f[dM0P] = &DD27[dP0M * size_Mat];
//             D27.f[dP0M] = &DD27[dM0P * size_Mat];
//             D27.f[d0MM] = &DD27[d0PP * size_Mat];
//             D27.f[d0PP] = &DD27[d0MM * size_Mat];
//             D27.f[d0MP] = &DD27[d0PM * size_Mat];
//             D27.f[d0PM] = &DD27[d0MP * size_Mat];
//             D27.f[d000] = &DD27[d000 * size_Mat];
//             D27.f[dMMM] = &DD27[dPPP * size_Mat];
//             D27.f[dPPM] = &DD27[dMMP * size_Mat];
//             D27.f[dMPM] = &DD27[dPMP * size_Mat];
//             D27.f[dPMM] = &DD27[dMPP * size_Mat];
//             D27.f[dMMP] = &DD27[dPPM * size_Mat];
//             D27.f[dPPP] = &DD27[dMMM * size_Mat];
//             D27.f[dMPP] = &DD27[dPMM * size_Mat];
//             D27.f[dPMP] = &DD27[dMPM * size_Mat];
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

//          (D27.f[d000])[kzero] =   c8o27* ConcD*(c1o1-cu_sq);
//          (D27.f[dP00])[ke   ] =   c2o27* ConcD*(c1o1+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
//          (D27.f[dM00])[kw   ] =   c2o27* ConcD*(c1o1+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
//          (D27.f[d0P0])[kn   ] =   c2o27* ConcD*(c1o1+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
//          (D27.f[d0M0])[ks   ] =   c2o27* ConcD*(c1o1+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
//          (D27.f[d00P])[kt   ] =   c2o27* ConcD*(c1o1+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
//          (D27.f[d00M])[kb   ] =   c2o27* ConcD*(c1o1+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
//          (D27.f[dPP0])[kne  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
//          (D27.f[dMM0])[ksw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
//          (D27.f[dPM0])[kse  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
//          (D27.f[dMP0])[knw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
//          (D27.f[dP0P])[kte  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
//          (D27.f[dM0M])[kbw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
//          (D27.f[dP0M])[kbe  ] =   c1o54* ConcD*(c1o1+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
//          (D27.f[dM0P])[ktw  ] =   c1o54* ConcD*(c1o1+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
//          (D27.f[d0PP])[ktn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
//          (D27.f[d0MM])[kbs  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
//          (D27.f[d0PM])[kbn  ] =   c1o54* ConcD*(c1o1+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
//          (D27.f[d0MP])[kts  ] =   c1o54* ConcD*(c1o1+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
//          (D27.f[dPPP])[ktne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
//          (D27.f[dMMM])[kbsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
//          (D27.f[dPPM])[kbne ] =   c1o216*ConcD*(c1o1+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
//          (D27.f[dMMP])[ktsw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
//          (D27.f[dPMP])[ktse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
//          (D27.f[dMPM])[kbnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
//          (D27.f[dPMM])[kbse ] =   c1o216*ConcD*(c1o1+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
//          (D27.f[dMPP])[ktnw ] =   c1o216*ConcD*(c1o1+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
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
                                    unsigned long long numberOfLBnodes,
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

   if(k<numberOfLBnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions7 D7;
         if (EvenOrOdd==true)
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