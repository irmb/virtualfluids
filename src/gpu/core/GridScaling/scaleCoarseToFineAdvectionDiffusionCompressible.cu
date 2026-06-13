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
//! \addtogroup gpu_GridScaling GridScaling
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include <basics/DataTypes.h>
#include "Utilities/KernelUtilities.h"
#include "cuda_helper/CudaIndexCalculation.h"
#include "Utilities/ScalingUtilities.h"

#include <lbm/refinement/InterpolationAdvectionDiffusionCF.h>
#include <lbm/interpolation/InterpolationCoefficientsAdvectionDiffusion.h>
#include <lbm/collision/TurbulentViscosity.h>

namespace vf::gpu {

template <bool hasTurbulentDiffusivity> __device__ void interpolate(
    vf::lbm::ad::InterpolationCoefficients& coefficients,
    const uint nodeIndex,
    real* distributionsFine, 
    real* distributionsFineAD, 
    const uint* neighborXfine,
    const uint* neighborYfine,
    const uint* neighborZfine,
    const unsigned long long numberOfLBnodesFine,
    const uint* indicesFineMMM,
    const real omegaDiffusivityFine,
    const real* turbulentDiffusivityFine)
{
    using namespace vf::basics::constant;
    const Distributions27 distFine = getDistributionReferences27(distributionsFine, numberOfLBnodesFine, true);
    Distributions27 distFineAD = getDistributionReferences27(distributionsFineAD, numberOfLBnodesFine, true);
    real populationFine[27];
    real populationFineAD[27];
    const real epsilonNew = c1o2; // ratio of grid resolutions

     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BSW = MMM: -0.25, -0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real x = -c1o4;
    real y = -c1o4;
    real z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    uint k_base_000 = indicesFineMMM[nodeIndex];
    uint k_base_M00 = neighborXfine [k_base_000];
    uint k_base_0M0 = neighborYfine [k_base_000];
    uint k_base_00M = neighborZfine [k_base_000];
    uint k_base_MM0 = neighborYfine [k_base_M00];
    uint k_base_M0M = neighborZfine [k_base_M00];
    uint k_base_0MM = neighborZfine [k_base_0M0];
    uint k_base_MMM = neighborZfine [k_base_MM0];
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    ListIndices indices;
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;
    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (zeroth to sixth order) on destination node
    //!
    real omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TSW = MMP: -0.25, -0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y = -c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TSE = PMP: 0.25, -0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y = -c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BSE = PMM: 0.25, -0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y = -c1o4;
    z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BNW = MPM: -0.25, 0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y =  c1o4;
    z = -c1o4;
    
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborYfine[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborYfine[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborYfine[k_base_0MM];
    k_base_MMM = neighborYfine[k_base_MMM];

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

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TNW = MPP: -0.25, 0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y =  c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TNE = PPP: 0.25, 0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = c1o4;
    y = c1o4;
    z = c1o4;
    ////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //Position BNE = PPM: 0.25, 0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y =  c1o4;
    z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];

    omegaF = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivityFine, turbulentDiffusivityFine[indices.k_000]) : omegaDiffusivityFine;

    getPreCollisionDistribution(populationFine, distFine, indices);

    vf::lbm::ad::interpolateAdvectionDiffusionCF(populationFine, populationFineAD, omegaDiffusivityFine, epsilonNew, coefficients, x, y, z);

    setPreCollisionDistribution(distFineAD, indices, populationFineAD);
}

//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from coarse to fine nodes
//! \details This scaling function is designed for advection diffusion F16 kernel.
//!
//! The function is executed in the following steps:
//!
template<bool hasTurbulentDiffusivity> 
__global__ void scaleCoarseToFineAdvectionDiffusionCompressible_Device(
    real* distributionsCoarse, 
    real* distributionsFine, 
    real* distributionsCoarseAD, 
    real* distributionsFineAD, 
    const uint* neighborXcoarse,
    const uint* neighborYcoarse,
    const uint* neighborZcoarse,
    const uint* neighborXfine,
    const uint* neighborYfine,
    const uint* neighborZfine,
    const unsigned long long numberOfLBnodesCoarse, 
    const unsigned long long numberOfLBnodesFine, 
    const bool isEvenTimestep,
    const uint* indicesCoarseMMM, 
    const uint* indicesFineMMM, 
    const uint numberOfInterfaceNodes, 
    const real omegaDiffusivityCoarse, 
    const real omegaDiffusivityFine, 
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborCoarseToFine)
{
    const unsigned nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfInterfaceNodes)
        return;

    // 1.calculate moments
    vf::lbm::ad::MomentsOnSourceNodeSet momentsSet;

    ad::calculateMomentSet<hasTurbulentDiffusivity>(
        momentsSet, 
        nodeIndex, 
        distributionsCoarse, 
        distributionsCoarseAD,
        neighborXcoarse, 
        neighborYcoarse, 
        neighborZcoarse, 
        indicesCoarseMMM,
        turbulentDiffusivityCoarse, 
        numberOfLBnodesCoarse, 
        omegaDiffusivityCoarse, 
        isEvenTimestep);

    // 2.calculate coefficients
    vf::lbm::ad::InterpolationCoefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, neighborCoarseToFine.x[nodeIndex], neighborCoarseToFine.y[nodeIndex], neighborCoarseToFine.z[nodeIndex]);

    // 3. interpolate coarse to fine
    interpolate<hasTurbulentDiffusivity>(
        coefficients,
        nodeIndex,
        distributionsFine, 
        distributionsFineAD, 
        neighborXfine,
        neighborYfine,
        neighborZfine,
        numberOfLBnodesFine,
        indicesFineMMM,
        omegaDiffusivityFine,
        turbulentDiffusivityFine);
}

template __global__ void scaleCoarseToFineAdvectionDiffusionCompressible_Device<true>(
    real* distributionsCoarse, 
    real* distributionsFine, 
    real* distributionsCoarseAD, 
    real* distributionsFineAD, 
    const uint* neighborXcoarse,
    const uint* neighborYcoarse,
    const uint* neighborZcoarse,
    const uint* neighborXfine,
    const uint* neighborYfine,
    const uint* neighborZfine,
    const unsigned long long numberOfLBnodesCoarse, 
    const unsigned long long numberOfLBnodesFine, 
    const bool isEvenTimestep,
    const uint* indicesCoarseMMM, 
    const uint* indicesFineMMM, 
    const uint numberOfInterfaceNodes, 
    const real omegaDiffusivityCoarse, 
    const real omegaDiffusivityFine, 
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborCoarseToFine);

template __global__ void scaleCoarseToFineAdvectionDiffusionCompressible_Device<false>(
    real* distributionsCoarse, 
    real* distributionsFine, 
    real* distributionsCoarseAD, 
    real* distributionsFineAD, 
    const uint* neighborXcoarse,
    const uint* neighborYcoarse,
    const uint* neighborZcoarse,
    const uint* neighborXfine,
    const uint* neighborYfine,
    const uint* neighborZfine,
    const unsigned long long numberOfLBnodesCoarse, 
    const unsigned long long numberOfLBnodesFine, 
    const bool isEvenTimestep,
    const uint* indicesCoarseMMM, 
    const uint* indicesFineMMM, 
    const uint numberOfInterfaceNodes, 
    const real omegaDiffusivityCoarse, 
    const real omegaDiffusivityFine, 
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborCoarseToFine);

}
//! \}
