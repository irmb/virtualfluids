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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
#include "Calculation/Calculation.h" 
#include <lbm/constants/D3Q27.h>
#include <lbm/ChimeraTransformation.h>
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void F16CompressibleAdvectionDiffusion_Device(
    real omegaDiffusivity,
    uint* typeOfGridNode,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    real* distributions,
    real* distributionsAD,
    unsigned long long numberOfLBnodes,
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
    if ((k < numberOfLBnodes) && (typeOfGridNode[k] == GEO_FLUID))
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
        //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        if (isEvenTimestep)
        {
            dist.f[dP00] = &distributions[dP00 * numberOfLBnodes];
            dist.f[dM00] = &distributions[dM00 * numberOfLBnodes];
            dist.f[d0P0] = &distributions[d0P0 * numberOfLBnodes];
            dist.f[d0M0] = &distributions[d0M0 * numberOfLBnodes];
            dist.f[d00P] = &distributions[d00P * numberOfLBnodes];
            dist.f[d00M] = &distributions[d00M * numberOfLBnodes];
            dist.f[dPP0] = &distributions[dPP0 * numberOfLBnodes];
            dist.f[dMM0] = &distributions[dMM0 * numberOfLBnodes];
            dist.f[dPM0] = &distributions[dPM0 * numberOfLBnodes];
            dist.f[dMP0] = &distributions[dMP0 * numberOfLBnodes];
            dist.f[dP0P] = &distributions[dP0P * numberOfLBnodes];
            dist.f[dM0M] = &distributions[dM0M * numberOfLBnodes];
            dist.f[dP0M] = &distributions[dP0M * numberOfLBnodes];
            dist.f[dM0P] = &distributions[dM0P * numberOfLBnodes];
            dist.f[d0PP] = &distributions[d0PP * numberOfLBnodes];
            dist.f[d0MM] = &distributions[d0MM * numberOfLBnodes];
            dist.f[d0PM] = &distributions[d0PM * numberOfLBnodes];
            dist.f[d0MP] = &distributions[d0MP * numberOfLBnodes];
            dist.f[d000] = &distributions[d000 * numberOfLBnodes];
            dist.f[dPPP] = &distributions[dPPP * numberOfLBnodes];
            dist.f[dMMP] = &distributions[dMMP * numberOfLBnodes];
            dist.f[dPMP] = &distributions[dPMP * numberOfLBnodes];
            dist.f[dMPP] = &distributions[dMPP * numberOfLBnodes];
            dist.f[dPPM] = &distributions[dPPM * numberOfLBnodes];
            dist.f[dMMM] = &distributions[dMMM * numberOfLBnodes];
            dist.f[dPMM] = &distributions[dPMM * numberOfLBnodes];
            dist.f[dMPM] = &distributions[dMPM * numberOfLBnodes];
        }
        else
        {
            dist.f[dM00] = &distributions[dP00 * numberOfLBnodes];
            dist.f[dP00] = &distributions[dM00 * numberOfLBnodes];
            dist.f[d0M0] = &distributions[d0P0 * numberOfLBnodes];
            dist.f[d0P0] = &distributions[d0M0 * numberOfLBnodes];
            dist.f[d00M] = &distributions[d00P * numberOfLBnodes];
            dist.f[d00P] = &distributions[d00M * numberOfLBnodes];
            dist.f[dMM0] = &distributions[dPP0 * numberOfLBnodes];
            dist.f[dPP0] = &distributions[dMM0 * numberOfLBnodes];
            dist.f[dMP0] = &distributions[dPM0 * numberOfLBnodes];
            dist.f[dPM0] = &distributions[dMP0 * numberOfLBnodes];
            dist.f[dM0M] = &distributions[dP0P * numberOfLBnodes];
            dist.f[dP0P] = &distributions[dM0M * numberOfLBnodes];
            dist.f[dM0P] = &distributions[dP0M * numberOfLBnodes];
            dist.f[dP0M] = &distributions[dM0P * numberOfLBnodes];
            dist.f[d0MM] = &distributions[d0PP * numberOfLBnodes];
            dist.f[d0PP] = &distributions[d0MM * numberOfLBnodes];
            dist.f[d0MP] = &distributions[d0PM * numberOfLBnodes];
            dist.f[d0PM] = &distributions[d0MP * numberOfLBnodes];
            dist.f[d000] = &distributions[d000 * numberOfLBnodes];
            dist.f[dMMM] = &distributions[dPPP * numberOfLBnodes];
            dist.f[dPPM] = &distributions[dMMP * numberOfLBnodes];
            dist.f[dMPM] = &distributions[dPMP * numberOfLBnodes];
            dist.f[dPMM] = &distributions[dMPP * numberOfLBnodes];
            dist.f[dMMP] = &distributions[dPPM * numberOfLBnodes];
            dist.f[dPPP] = &distributions[dMMM * numberOfLBnodes];
            dist.f[dMPP] = &distributions[dPMM * numberOfLBnodes];
            dist.f[dPMP] = &distributions[dMPM * numberOfLBnodes];
        }
        ////////////////////////////////////////////////////////////////////////////////
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
        real fcbb = (dist.f[dP00])[k];
        real fabb = (dist.f[dM00])[kw];
        real fbcb = (dist.f[d0P0])[k];
        real fbab = (dist.f[d0M0])[ks];
        real fbbc = (dist.f[d00P])[k];
        real fbba = (dist.f[d00M])[kb];
        real fccb = (dist.f[dPP0])[k];
        real faab = (dist.f[dMM0])[ksw];
        real fcab = (dist.f[dPM0])[ks];
        real facb = (dist.f[dMP0])[kw];
        real fcbc = (dist.f[dP0P])[k];
        real faba = (dist.f[dM0M])[kbw];
        real fcba = (dist.f[dP0M])[kb];
        real fabc = (dist.f[dM0P])[kw];
        real fbcc = (dist.f[d0PP])[k];
        real fbaa = (dist.f[d0MM])[kbs];
        real fbca = (dist.f[d0PM])[kb];
        real fbac = (dist.f[d0MP])[ks];
        real fbbb = (dist.f[d000])[k];
        real fccc = (dist.f[dPPP])[k];
        real faac = (dist.f[dMMP])[ksw];
        real fcac = (dist.f[dPMP])[ks];
        real facc = (dist.f[dMPP])[kw];
        real fcca = (dist.f[dPPM])[kb];
        real faaa = (dist.f[dMMM])[kbsw];
        real fcaa = (dist.f[dPMM])[kbs];
        real faca = (dist.f[dMPM])[kbw];
        ////////////////////////////////////////////////////////////////////////////////////
        //! - Set local distributions Advection Diffusion
        //!
        real mfcbb = (distAD.f[dP00])[k];
        real mfabb = (distAD.f[dM00])[kw];
        real mfbcb = (distAD.f[d0P0])[k];
        real mfbab = (distAD.f[d0M0])[ks];
        real mfbbc = (distAD.f[d00P])[k];
        real mfbba = (distAD.f[d00M])[kb];
        real mfccb = (distAD.f[dPP0])[k];
        real mfaab = (distAD.f[dMM0])[ksw];
        real mfcab = (distAD.f[dPM0])[ks];
        real mfacb = (distAD.f[dMP0])[kw];
        real mfcbc = (distAD.f[dP0P])[k];
        real mfaba = (distAD.f[dM0M])[kbw];
        real mfcba = (distAD.f[dP0M])[kb];
        real mfabc = (distAD.f[dM0P])[kw];
        real mfbcc = (distAD.f[d0PP])[k];
        real mfbaa = (distAD.f[d0MM])[kbs];
        real mfbca = (distAD.f[d0PM])[kb];
        real mfbac = (distAD.f[d0MP])[ks];
        real mfbbb = (distAD.f[d000])[k];
        real mfccc = (distAD.f[dPPP])[k];
        real mfaac = (distAD.f[dMMP])[ksw];
        real mfcac = (distAD.f[dPMP])[ks];
        real mfacc = (distAD.f[dMPP])[kw];
        real mfcca = (distAD.f[dPPM])[kb];
        real mfaaa = (distAD.f[dMMM])[kbsw];
        real mfcaa = (distAD.f[dPMM])[kbs];
        real mfaca = (distAD.f[dMPM])[kbw];
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
        (distAD.f[dP00])[k   ] = mfabb;
        (distAD.f[dM00])[kw  ] = mfcbb;
        (distAD.f[d0P0])[k   ] = mfbab;
        (distAD.f[d0M0])[ks  ] = mfbcb;
        (distAD.f[d00P])[k   ] = mfbba;
        (distAD.f[d00M])[kb  ] = mfbbc;
        (distAD.f[dPP0])[k   ] = mfaab;
        (distAD.f[dMM0])[ksw ] = mfccb;
        (distAD.f[dPM0])[ks  ] = mfacb;
        (distAD.f[dMP0])[kw  ] = mfcab;
        (distAD.f[dP0P])[k   ] = mfaba;
        (distAD.f[dM0M])[kbw ] = mfcbc;
        (distAD.f[dP0M])[kb  ] = mfabc;
        (distAD.f[dM0P])[kw  ] = mfcba;
        (distAD.f[d0PP])[k   ] = mfbaa;
        (distAD.f[d0MM])[kbs ] = mfbcc;
        (distAD.f[d0PM])[kb  ] = mfbac;
        (distAD.f[d0MP])[ks  ] = mfbca;
        (distAD.f[d000])[k   ] = mfbbb;
        (distAD.f[dPPP])[k   ] = mfaaa;
        (distAD.f[dPMP])[ks  ] = mfaca;
        (distAD.f[dPPM])[kb  ] = mfaac;
        (distAD.f[dPMM])[kbs ] = mfacc;
        (distAD.f[dMPP])[kw  ] = mfcaa;
        (distAD.f[dMMP])[ksw ] = mfcca;
        (distAD.f[dMPM])[kbw ] = mfcac;
        (distAD.f[dMMM])[kbsw] = mfccc;
    }
}
////////////////////////////////////////////////////////////////////////////////


//! \}
