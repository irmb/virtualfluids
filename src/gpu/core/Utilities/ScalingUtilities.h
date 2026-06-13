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
//! \addtogroup gpu_Utilities utilities
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
#ifndef SCALING_HELPER_FUNCTIONS_H
#define SCALING_HELPER_FUNCTIONS_H

#include "Utilities/KernelUtilities.h"

#include <lbm/constants/D3Q27.h>
#include <basics/DataTypes.h>

#include <lbm/interpolation/InterpolationCoefficients.h>
#include <lbm/interpolation/InterpolationCoefficientsAdvectionDiffusion.h>
#include <lbm/collision/TurbulentViscosity.h>

namespace vf::gpu
{

template<bool hasTurbulentViscosity> __device__ void calculateMomentSet(
    vf::lbm::MomentsOnSourceNodeSet& momentsSet,
    const unsigned nodeIndex,
    real *distribution,
    unsigned int *neighborX,
    unsigned int *neighborY,
    unsigned int *neighborZ,
    unsigned int *indices_MMM,
    real* turbulentViscosity,
    unsigned long long numberOfLBnodes,
    const real omega,
    bool isEvenTimestep
)
{
    real omega_ = omega;
    Distributions27 populationReferences;
    getPointersToDistributions(populationReferences, distribution, numberOfLBnodes, isEvenTimestep);

    ListIndices indices;

    //////////////////////////////////////////////////////////////////////////
    //! - Calculate moments for each source node 
    //!
    //////////////////////////////////////////////////////////////////////////
    // source node BSW = MMM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    unsigned int k_base_000 = indices_MMM[nodeIndex];
    unsigned int k_base_M00 = neighborX [k_base_000];
    unsigned int k_base_0M0 = neighborY [k_base_000];
    unsigned int k_base_00M = neighborZ [k_base_000];
    unsigned int k_base_MM0 = neighborY [k_base_M00];
    unsigned int k_base_M0M = neighborZ [k_base_M00];
    unsigned int k_base_0MM = neighborZ [k_base_0M0];
    unsigned int k_base_MMM = neighborZ [k_base_MM0];
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    real populations[27];

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculateMMM(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZ[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZ[indices.k_M0M];
    indices.k_0MM = neighborZ[indices.k_0MM];
    indices.k_MMM = neighborZ[indices.k_MMM];

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculateMMP(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborX[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborX[indices.k_MM0];
    indices.k_M0M = neighborX[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborX[indices.k_MMM];

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculatePMP(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM 
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborX[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborX[k_base_MM0];

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculatePMM(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node BNW = MPM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors --> indices of all source nodes
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborY[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborY[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborY[k_base_0MM];
    k_base_MMM = neighborY[k_base_MMM];
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculateMPM(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZ[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZ[indices.k_M0M];
    indices.k_0MM = neighborZ[indices.k_0MM];
    indices.k_MMM = neighborZ[indices.k_MMM];

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculateMPP(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborX[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborX[indices.k_MM0];
    indices.k_M0M = neighborX[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborX[indices.k_MMM];

    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculatePPP(populations, omega_);

    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborX[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborX[k_base_MM0];
    
    omega_ = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulentViscosity[indices.k_000]) : omega;

    getPreCollisionDistribution(populations, populationReferences, indices);
    momentsSet.calculatePPM(populations, omega_);
}

}


namespace vf::gpu::ad
{

template <bool hasTurbulentDiffusivity>
__device__ void calculateMomentSet(
    vf::lbm::ad::MomentsOnSourceNodeSet& momentsSet, 
    const unsigned nodeIndex, 
    real* distribution,
    real* distributionAD,
    const uint* neighborX, 
    const uint* neighborY, 
    const uint* neighborZ,
    const uint* indices_MMM, 
    const real* turbulentDiffusivity, 
    const unsigned long long numberOfLBnodes,
    const real omegaDiffusivity, 
    const bool isEvenTimestep
)
{
    real omega = omegaDiffusivity;
    const Distributions27 populationReferences =
        getDistributionReferences27(distribution, numberOfLBnodes, isEvenTimestep);
    const Distributions27 populationReferencesAD =
        getDistributionReferences27(distributionAD, numberOfLBnodes, isEvenTimestep);

    ListIndices indices;
    ////////////////////////////////////////////////////////////////////////////////
    //! - Calculate moments for each source node
    //!
    ////////////////////////////////////////////////////////////////////////////////
    // source node BSW = MMM
    ////////////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    uint k_base_000 = indices_MMM[nodeIndex];
    uint k_base_M00 = neighborX[k_base_000];
    uint k_base_0M0 = neighborY[k_base_000];
    uint k_base_00M = neighborZ[k_base_000];
    uint k_base_MM0 = neighborY[k_base_M00];
    uint k_base_M0M = neighborZ[k_base_M00];
    uint k_base_0MM = neighborZ[k_base_0M0];
    uint k_base_MMM = neighborZ[k_base_MM0];
    ////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    real populations[27];
    real populationAD[27];

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculateMMM(populations, populationAD, omega);

    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZ[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZ[indices.k_M0M];
    indices.k_0MM = neighborZ[indices.k_0MM];
    indices.k_MMM = neighborZ[indices.k_MMM];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculateMMP(populations, populationAD, omega);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborX[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborX[indices.k_MM0];
    indices.k_M0M = neighborX[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborX[indices.k_MMM];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculatePMP(populations, populationAD, omega);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborX[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborX[k_base_MM0];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculatePMM(populations, populationAD, omega);

    //////////////////////////////////////////////////////////////////////////
    // source node BNW = MPM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors --> indices of all source nodes
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborY[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborY[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborY[k_base_0MM];
    k_base_MMM = neighborY[k_base_MMM];
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculateMPM(populations, populationAD, omega);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZ[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZ[indices.k_M0M];
    indices.k_0MM = neighborZ[indices.k_0MM];
    indices.k_MMM = neighborZ[indices.k_MMM];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculateMPP(populations, populationAD, omega);
    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborX[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborX[indices.k_MM0];
    indices.k_M0M = neighborX[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborX[indices.k_MMM];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculatePPP(populations, populationAD, omega);
    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborX[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborX[k_base_MM0];

    omega = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(
                                           omegaDiffusivity, turbulentDiffusivity[indices.k_000])
                                     : omegaDiffusivity;

    getPreCollisionDistribution(populations, populationReferences, indices);
    getPreCollisionDistribution(populationAD, populationReferencesAD, indices);
    momentsSet.calculatePPM(populations, populationAD, omega);




}

}


#endif

//! \}
