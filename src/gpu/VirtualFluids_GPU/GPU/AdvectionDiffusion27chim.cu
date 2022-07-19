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
#include "lbm/constants/D3Q27.h"

#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

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
__global__ void Factorized_Central_Moments_Advection_Diffusion_Device_Kernel(
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
			dist.f[E   ] = &distributions[E   *size_Mat];
			dist.f[W   ] = &distributions[W   *size_Mat];
			dist.f[N   ] = &distributions[N   *size_Mat];
			dist.f[S   ] = &distributions[S   *size_Mat];
			dist.f[T   ] = &distributions[T   *size_Mat];
			dist.f[B   ] = &distributions[B   *size_Mat];
			dist.f[NE  ] = &distributions[NE  *size_Mat];
			dist.f[SW  ] = &distributions[SW  *size_Mat];
			dist.f[SE  ] = &distributions[SE  *size_Mat];
			dist.f[NW  ] = &distributions[NW  *size_Mat];
			dist.f[TE  ] = &distributions[TE  *size_Mat];
			dist.f[BW  ] = &distributions[BW  *size_Mat];
			dist.f[BE  ] = &distributions[BE  *size_Mat];
			dist.f[TW  ] = &distributions[TW  *size_Mat];
			dist.f[TN  ] = &distributions[TN  *size_Mat];
			dist.f[BS  ] = &distributions[BS  *size_Mat];
			dist.f[BN  ] = &distributions[BN  *size_Mat];
			dist.f[TS  ] = &distributions[TS  *size_Mat];
			dist.f[REST] = &distributions[REST*size_Mat];
			dist.f[TNE ] = &distributions[TNE *size_Mat];
			dist.f[TSW ] = &distributions[TSW *size_Mat];
			dist.f[TSE ] = &distributions[TSE *size_Mat];
			dist.f[TNW ] = &distributions[TNW *size_Mat];
			dist.f[BNE ] = &distributions[BNE *size_Mat];
			dist.f[BSW ] = &distributions[BSW *size_Mat];
			dist.f[BSE ] = &distributions[BSE *size_Mat];
			dist.f[BNW ] = &distributions[BNW *size_Mat];
		}
		else
		{
			dist.f[W   ] = &distributions[E   *size_Mat];
			dist.f[E   ] = &distributions[W   *size_Mat];
			dist.f[S   ] = &distributions[N   *size_Mat];
			dist.f[N   ] = &distributions[S   *size_Mat];
			dist.f[B   ] = &distributions[T   *size_Mat];
			dist.f[T   ] = &distributions[B   *size_Mat];
			dist.f[SW  ] = &distributions[NE  *size_Mat];
			dist.f[NE  ] = &distributions[SW  *size_Mat];
			dist.f[NW  ] = &distributions[SE  *size_Mat];
			dist.f[SE  ] = &distributions[NW  *size_Mat];
			dist.f[BW  ] = &distributions[TE  *size_Mat];
			dist.f[TE  ] = &distributions[BW  *size_Mat];
			dist.f[TW  ] = &distributions[BE  *size_Mat];
			dist.f[BE  ] = &distributions[TW  *size_Mat];
			dist.f[BS  ] = &distributions[TN  *size_Mat];
			dist.f[TN  ] = &distributions[BS  *size_Mat];
			dist.f[TS  ] = &distributions[BN  *size_Mat];
			dist.f[BN  ] = &distributions[TS  *size_Mat];
			dist.f[REST] = &distributions[REST*size_Mat];
			dist.f[BSW ] = &distributions[TNE *size_Mat];
			dist.f[BNE ] = &distributions[TSW *size_Mat];
			dist.f[BNW ] = &distributions[TSE *size_Mat];
			dist.f[BSE ] = &distributions[TNW *size_Mat];
			dist.f[TSW ] = &distributions[BNE *size_Mat];
			dist.f[TNE ] = &distributions[BSW *size_Mat];
			dist.f[TNW ] = &distributions[BSE *size_Mat];
			dist.f[TSE ] = &distributions[BNW *size_Mat];
		}
		////////////////////////////////////////////////////////////////////////////////
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
		real fcbb = (dist.f[E   ])[k];
		real fabb = (dist.f[W   ])[kw];
		real fbcb = (dist.f[N   ])[k];
		real fbab = (dist.f[S   ])[ks];
		real fbbc = (dist.f[T   ])[k];
		real fbba = (dist.f[B   ])[kb];
		real fccb = (dist.f[NE  ])[k];
		real faab = (dist.f[SW  ])[ksw];
		real fcab = (dist.f[SE  ])[ks];
		real facb = (dist.f[NW  ])[kw];
		real fcbc = (dist.f[TE  ])[k];
		real faba = (dist.f[BW  ])[kbw];
		real fcba = (dist.f[BE  ])[kb];
		real fabc = (dist.f[TW  ])[kw];
		real fbcc = (dist.f[TN  ])[k];
		real fbaa = (dist.f[BS  ])[kbs];
		real fbca = (dist.f[BN  ])[kb];
		real fbac = (dist.f[TS  ])[ks];
		real fbbb = (dist.f[REST])[k];
		real fccc = (dist.f[TNE ])[k];
		real faac = (dist.f[TSW ])[ksw];
		real fcac = (dist.f[TSE ])[ks];
		real facc = (dist.f[TNW ])[kw];
		real fcca = (dist.f[BNE ])[kb];
		real faaa = (dist.f[BSW ])[kbsw];
		real fcaa = (dist.f[BSE ])[kbs];
		real faca = (dist.f[BNW ])[kbw];
		////////////////////////////////////////////////////////////////////////////////////
		//! - Set local distributions Advection Diffusion
		//!
		real mfcbb = (distAD.f[E   ])[k];
		real mfabb = (distAD.f[W   ])[kw];
		real mfbcb = (distAD.f[N   ])[k];
		real mfbab = (distAD.f[S   ])[ks];
		real mfbbc = (distAD.f[T   ])[k];
		real mfbba = (distAD.f[B   ])[kb];
		real mfccb = (distAD.f[NE  ])[k];
		real mfaab = (distAD.f[SW  ])[ksw];
		real mfcab = (distAD.f[SE  ])[ks];
		real mfacb = (distAD.f[NW  ])[kw];
		real mfcbc = (distAD.f[TE  ])[k];
		real mfaba = (distAD.f[BW  ])[kbw];
		real mfcba = (distAD.f[BE  ])[kb];
		real mfabc = (distAD.f[TW  ])[kw];
		real mfbcc = (distAD.f[TN  ])[k];
		real mfbaa = (distAD.f[BS  ])[kbs];
		real mfbca = (distAD.f[BN  ])[kb];
		real mfbac = (distAD.f[TS  ])[ks];
		real mfbbb = (distAD.f[REST])[k];
		real mfccc = (distAD.f[TNE ])[k];
		real mfaac = (distAD.f[TSW ])[ksw];
		real mfcac = (distAD.f[TSE ])[ks];
		real mfacc = (distAD.f[TNW ])[kw];
		real mfcca = (distAD.f[BNE ])[kb];
		real mfaaa = (distAD.f[BSW ])[kbsw];
		real mfcaa = (distAD.f[BSE ])[kbs];
		real mfaca = (distAD.f[BNW ])[kbw];
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
		(distAD.f[E   ])[k   ] = mfabb;
		(distAD.f[W   ])[kw  ] = mfcbb;
		(distAD.f[N   ])[k   ] = mfbab;
		(distAD.f[S   ])[ks  ] = mfbcb;
		(distAD.f[T   ])[k   ] = mfbba;
		(distAD.f[B   ])[kb  ] = mfbbc;
		(distAD.f[NE  ])[k   ] = mfaab;
		(distAD.f[SW  ])[ksw ] = mfccb;
		(distAD.f[SE  ])[ks  ] = mfacb;
		(distAD.f[NW  ])[kw  ] = mfcab;
		(distAD.f[TE  ])[k   ] = mfaba;
		(distAD.f[BW  ])[kbw ] = mfcbc;
		(distAD.f[BE  ])[kb  ] = mfabc;
		(distAD.f[TW  ])[kw  ] = mfcba;
		(distAD.f[TN  ])[k   ] = mfbaa;
		(distAD.f[BS  ])[kbs ] = mfbcc;
		(distAD.f[BN  ])[kb  ] = mfbac;
		(distAD.f[TS  ])[ks  ] = mfbca;
		(distAD.f[REST])[k   ] = mfbbb;
		(distAD.f[TNE ])[k   ] = mfaaa;
		(distAD.f[TSE ])[ks  ] = mfaca;
		(distAD.f[BNE ])[kb  ] = mfaac;
		(distAD.f[BSE ])[kbs ] = mfacc;
		(distAD.f[TNW ])[kw  ] = mfcaa;
		(distAD.f[TSW ])[ksw ] = mfcca;
		(distAD.f[BNW ])[kbw ] = mfcac;
		(distAD.f[BSW ])[kbsw] = mfccc;
	}
}
////////////////////////////////////////////////////////////////////////////////

