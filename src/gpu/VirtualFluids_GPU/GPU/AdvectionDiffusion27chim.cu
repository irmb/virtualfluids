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
//! \file AdvectionDiffusion27chim.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;


////////////////////////////////////////////////////////////////////////////////
//! \brief forward chimera transformation \ref forwardChimera
//! - Chimera transform from distributions to central moments as defined in Eq. (43)-(45) in \ref
//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
inline __device__ void forwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2) {
	real m1 = (mfa + mfc) + mfb;
	real m2 = mfc - mfa;
	mfc     = (mfc + mfa) + (v2*m1 - c2o1*vv*m2);
	mfb     = m2 - vv*m1;
	mfa     = m1;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief backward chimera transformation \ref backwardChimera
//! - Chimera transform from  central moments to distributions as defined in Eq. (88)-(96) in \ref
//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
inline __device__ void backwardChimera(real &mfa, real &mfb, real &mfc, real vv, real v2) {
	real ma = (mfc + mfa*(v2 - vv))*c1o2 + mfb*(vv - c1o2);
	real mb = ((mfa - mfc) - mfa*v2) - c2o1*mfb*vv;
	mfc     = (mfc + mfa*(v2 + vv))*c1o2 + mfb*(vv + c1o2);
	mfb     = mb;
	mfa     = ma;
}


////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void Factorized_Central_Moments_Advection_Diffusion_Device_Kernel(
	real omegaDiffusivity,
	uint* typeOfGridNode,
	uint* neighborX,
	uint* neighborY,
	uint* neighborZ,
	real* distributions,
	real* distributionsAD,
	int size_Mat,
	real* forces,
	bool isEvenTimestep)
{
	//////////////////////////////////////////////////////////////////////////
	//! Cumulant K17 Kernel is based on \ref
	//! <a href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040 ]</b></a>
	//! and \ref
	//! <a href="https://doi.org/10.1016/j.jcp.2017.07.004"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.07.004 ]</b></a>
	//!
	//! The cumulant kernel is executed in the following steps
	//!
	////////////////////////////////////////////////////////////////////////////////
	//! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
	//!
	const unsigned  x = threadIdx.x;
	const unsigned  y = blockIdx.x;
	const unsigned  z = blockIdx.y;

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
		////////////////////////////////////////////////////////////////////////////////
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
		////////////////////////////////////////////////////////////////////////////////
		//! - Set neighbor indices (necessary for indirect addressing)
		uint kw   = neighborX[k];
		uint ks   = neighborY[k];
		uint kb   = neighborZ[k];
		uint ksw  = neighborY[kw];
		uint kbw  = neighborZ[kw];
		uint kbs  = neighborZ[ks];
		uint kbsw = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////////
		//! - Set local distributions Fluid
		//!
		real fcbb = (dist.f[dirE   ])[k];
		real fabb = (dist.f[dirW   ])[kw];
		real fbcb = (dist.f[dirN   ])[k];
		real fbab = (dist.f[dirS   ])[ks];
		real fbbc = (dist.f[dirT   ])[k];
		real fbba = (dist.f[dirB   ])[kb];
		real fccb = (dist.f[dirNE  ])[k];
		real faab = (dist.f[dirSW  ])[ksw];
		real fcab = (dist.f[dirSE  ])[ks];
		real facb = (dist.f[dirNW  ])[kw];
		real fcbc = (dist.f[dirTE  ])[k];
		real faba = (dist.f[dirBW  ])[kbw];
		real fcba = (dist.f[dirBE  ])[kb];
		real fabc = (dist.f[dirTW  ])[kw];
		real fbcc = (dist.f[dirTN  ])[k];
		real fbaa = (dist.f[dirBS  ])[kbs];
		real fbca = (dist.f[dirBN  ])[kb];
		real fbac = (dist.f[dirTS  ])[ks];
		real fbbb = (dist.f[dirREST])[k];
		real fccc = (dist.f[dirTNE ])[k];
		real faac = (dist.f[dirTSW ])[ksw];
		real fcac = (dist.f[dirTSE ])[ks];
		real facc = (dist.f[dirTNW ])[kw];
		real fcca = (dist.f[dirBNE ])[kb];
		real faaa = (dist.f[dirBSW ])[kbsw];
		real fcaa = (dist.f[dirBSE ])[kbs];
		real faca = (dist.f[dirBNW ])[kbw];
		////////////////////////////////////////////////////////////////////////////////////
		//! - Set local distributions Advection Diffusion
		//!
		real mfcbb = (distAD.f[dirE   ])[k];
		real mfabb = (distAD.f[dirW   ])[kw];
		real mfbcb = (distAD.f[dirN   ])[k];
		real mfbab = (distAD.f[dirS   ])[ks];
		real mfbbc = (distAD.f[dirT   ])[k];
		real mfbba = (distAD.f[dirB   ])[kb];
		real mfccb = (distAD.f[dirNE  ])[k];
		real mfaab = (distAD.f[dirSW  ])[ksw];
		real mfcab = (distAD.f[dirSE  ])[ks];
		real mfacb = (distAD.f[dirNW  ])[kw];
		real mfcbc = (distAD.f[dirTE  ])[k];
		real mfaba = (distAD.f[dirBW  ])[kbw];
		real mfcba = (distAD.f[dirBE  ])[kb];
		real mfabc = (distAD.f[dirTW  ])[kw];
		real mfbcc = (distAD.f[dirTN  ])[k];
		real mfbaa = (distAD.f[dirBS  ])[kbs];
		real mfbca = (distAD.f[dirBN  ])[kb];
		real mfbac = (distAD.f[dirTS  ])[ks];
		real mfbbb = (distAD.f[dirREST])[k];
		real mfccc = (distAD.f[dirTNE ])[k];
		real mfaac = (distAD.f[dirTSW ])[ksw];
		real mfcac = (distAD.f[dirTSE ])[ks];
		real mfacc = (distAD.f[dirTNW ])[kw];
		real mfcca = (distAD.f[dirBNE ])[kb];
		real mfaaa = (distAD.f[dirBSW ])[kbsw];
		real mfcaa = (distAD.f[dirBSE ])[kbs];
		real mfaca = (distAD.f[dirBNW ])[kbw];
		////////////////////////////////////////////////////////////////////////////////////
		//! - Calculate density and velocity using pyramid summation for low round-off errors as in Eq. (J1)-(J3) \ref
		//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
		//!
		////////////////////////////////////////////////////////////////////////////////////
		// fluid component
		real drhoFluid =
			((((fccc + faaa) + (faca + fcac)) + ((facc + fcaa) + (faac + fcca))) +
			(((fbac + fbca) + (fbaa + fbcc)) + ((fabc + fcba) + (faba + fcbc)) + ((facb + fcab) + (faab + fccb))) +
			((fabb + fcbb) + (fbab + fbcb) + (fbba + fbbc))) + fbbb;

		real rhoFluid = c1o1 + drhoFluid;
		real OOrhoFluid = c1o1 / rhoFluid;

        real vvx =
			((((fccc - faaa) + (fcac - faca)) + ((fcaa - facc) + (fcca - faac))) +
			(((fcba - fabc) + (fcbc - faba)) + ((fcab - facb) + (fccb - faab))) +
			(fcbb - fabb)) * OOrhoFluid;
		real vvy =
			((((fccc - faaa) + (faca - fcac)) + ((facc - fcaa) + (fcca - faac))) +
			(((fbca - fbac) + (fbcc - fbaa)) + ((facb - fcab) + (fccb - faab))) +
			(fbcb - fbab)) * OOrhoFluid;
		real vvz =
			((((fccc - faaa) + (fcac - faca)) + ((facc - fcaa) + (faac - fcca))) +
			(((fbac - fbca) + (fbcc - fbaa)) + ((fabc - fcba) + (fcbc - faba))) +
			(fbbc - fbba)) * OOrhoFluid;
		////////////////////////////////////////////////////////////////////////////////////
		// second component
		real rho =
			((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
			(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
				((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Add half of the acceleration (body force) to the velocity as in Eq. (42) \ref
        //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
        //!
        real fx = forces[0];
        real fy = forces[1];
        real fz = -rho*forces[2];
        vvx += fx*c1o2;
        vvy += fy*c1o2;
        vvz += fz*c1o2;
        ////////////////////////////////////////////////////////////////////////////////////
		// calculate the square of velocities for this lattice node
		real vx2 = vvx*vvx;
		real vy2 = vvy*vvy;
		real vz2 = vvz*vvz;
		////////////////////////////////////////////////////////////////////////////////////
		//real omegaDiffusivity = c2o1 / (c6o1 * diffusivity + c1o1);
		////////////////////////////////////////////////////////////////////////////////////
		//! - Chimera transform from distributions to central moments as defined in Eq. (43)-(45) in \ref
		//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
		//!
		////////////////////////////////////////////////////////////////////////////////////
		// Z - Dir
		forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
		forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
		forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
		forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
		forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
		forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
		forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
		forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
		forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

		////////////////////////////////////////////////////////////////////////////////////
		// Y - Dir
		forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
		forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
		forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
		forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
		forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
		forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
		forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
		forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
		forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

		////////////////////////////////////////////////////////////////////////////////////
		// X - Dir
		forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
		forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
		forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
		forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
		forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
		forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
		forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
		forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
		forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

		////////////////////////////////////////////////////////////////////////////////////
		//! - Factorized central moments for Advection Diffusion Equation - Eq. (15)-(16) in \ref
		//! <a href="https://doi.org/10.1016/j.advwatres.2015.09.015"><b>[ X. Yang et al. (2016), DOI: 10.1016/j.advwatres.2015.09.015]</b></a>
		//!

		// linearized orthogonalization of 3rd order central moments
		real Mabc = mfabc - mfaba*c1o3;
		real Mbca = mfbca - mfbaa*c1o3;
		real Macb = mfacb - mfaab*c1o3;
		real Mcba = mfcba - mfaba*c1o3;
		real Mcab = mfcab - mfaab*c1o3;
		real Mbac = mfbac - mfbaa*c1o3;
		// linearized orthogonalization of 5th order central moments
		real Mcbc = mfcbc - mfaba*c1o9;
		real Mbcc = mfbcc - mfbaa*c1o9;
		real Mccb = mfccb - mfaab*c1o9;

		// collision of 1st order moments
		mfbaa *= c1o1 - omegaDiffusivity;
		mfaba *= c1o1 - omegaDiffusivity;
		mfaab *= c1o1 - omegaDiffusivity;

		// equilibration of 3rd order moments
		Mabc = c0o1;
		Mbca = c0o1;
		Macb = c0o1;
		Mcba = c0o1;
		Mcab = c0o1;
		Mbac = c0o1;
		mfbbb = c0o1;

		// equilibration of 5th order moments
		Mcbc = c0o1;
		Mbcc = c0o1;
		Mccb = c0o1;

		// equilibration of 2nd order moments
		mfbba = c0o1;
		mfbab = c0o1;
		mfabb = c0o1;

		mfcaa = c1o3 * rho;
		mfaca = c1o3 * rho;
		mfaac = c1o3 * rho;

		// equilibration of 4th order moments
		mfacc = c1o9 * rho;
		mfcac = c1o9 * rho;
		mfcca = c1o9 * rho;

		mfcbb = c0o1;
		mfbcb = c0o1;
		mfbbc = c0o1;

		// equilibration of 6th order moment
		mfccc = c1o27 * rho;

		// from linearized orthogonalization 3rd order central moments to central moments
		mfabc = Mabc + mfaba*c1o3;
		mfbca = Mbca + mfbaa*c1o3;
		mfacb = Macb + mfaab*c1o3;
		mfcba = Mcba + mfaba*c1o3;
		mfcab = Mcab + mfaab*c1o3;
		mfbac = Mbac + mfbaa*c1o3;

		// from linearized orthogonalization 5th order central moments to central moments
		mfcbc = Mcbc + mfaba*c1o9;
		mfbcc = Mbcc + mfbaa*c1o9;
		mfccb = Mccb + mfaab*c1o9;

		////////////////////////////////////////////////////////////////////////////////////
		//! - Chimera transform from  central moments to distributions as defined in Eq. (88)-(96) in \ref
		//! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015), DOI:10.1016/j.camwa.2015.05.001 ]</b></a>
		//!
		////////////////////////////////////////////////////////////////////////////////////
		// X - Dir
		backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
		backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
		backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
		backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
		backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
		backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
		backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
		backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
		backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

		////////////////////////////////////////////////////////////////////////////////////
		// Y - Dir
		backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
		backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
		backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
		backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
		backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
		backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
		backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
		backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
		backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

		////////////////////////////////////////////////////////////////////////////////////
		// Z - Dir
		backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
		backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
		backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
		backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
		backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
		backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
		backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
		backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
		backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

		////////////////////////////////////////////////////////////////////////////////////
		//! - Write distributions: style of reading and writing the distributions from/to
		//! stored arrays dependent on timestep is based on the esoteric twist algorithm
		//! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
		//!
		(distAD.f[dirE   ])[k   ] = mfabb;
		(distAD.f[dirW   ])[kw  ] = mfcbb;
		(distAD.f[dirN   ])[k   ] = mfbab;
		(distAD.f[dirS   ])[ks  ] = mfbcb;
		(distAD.f[dirT   ])[k   ] = mfbba;
		(distAD.f[dirB   ])[kb  ] = mfbbc;
		(distAD.f[dirNE  ])[k   ] = mfaab;
		(distAD.f[dirSW  ])[ksw ] = mfccb;
		(distAD.f[dirSE  ])[ks  ] = mfacb;
		(distAD.f[dirNW  ])[kw  ] = mfcab;
		(distAD.f[dirTE  ])[k   ] = mfaba;
		(distAD.f[dirBW  ])[kbw ] = mfcbc;
		(distAD.f[dirBE  ])[kb  ] = mfabc;
		(distAD.f[dirTW  ])[kw  ] = mfcba;
		(distAD.f[dirTN  ])[k   ] = mfbaa;
		(distAD.f[dirBS  ])[kbs ] = mfbcc;
		(distAD.f[dirBN  ])[kb  ] = mfbac;
		(distAD.f[dirTS  ])[ks  ] = mfbca;
		(distAD.f[dirREST])[k   ] = mfbbb;
		(distAD.f[dirTNE ])[k   ] = mfaaa;
		(distAD.f[dirTSE ])[ks  ] = mfaca;
		(distAD.f[dirBNE ])[kb  ] = mfaac;
		(distAD.f[dirBSE ])[kbs ] = mfacc;
		(distAD.f[dirTNW ])[kw  ] = mfcaa;
		(distAD.f[dirTSW ])[ksw ] = mfcca;
		(distAD.f[dirBNW ])[kbw ] = mfcac;
		(distAD.f[dirBSW ])[kbsw] = mfccc;
	}
}
////////////////////////////////////////////////////////////////////////////////

