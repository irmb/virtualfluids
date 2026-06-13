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

#include "cuda_helper/CudaIndexCalculation.h"
#include "Utilities/KernelUtilities.h"
#include "Utilities/ScalingUtilities.h"

#include <lbm/refinement/InterpolationAdvectionDiffusionFC.h>
#include <lbm/interpolation/InterpolationCoefficientsAdvectionDiffusion.h>
#include <lbm/collision/TurbulentViscosity.h>

namespace vf::gpu {

template <bool hasTurbulentDiffusivity> __device__ void interpolate(
    vf::lbm::ad::InterpolationCoefficients& coefficients,
    const uint nodeIndex,
    real* distributionsCoarse,
    real* distributionsCoarseAD,
    const uint* neighborXcoarse,
    const uint* neighborYcoarse,
    const uint* neighborZcoarse,
    const unsigned long long numberOfLBnodesCoarse,
    const uint* indicesCoarse000,
    const real omegaCoarse,
    const real* turbulentDiffusivityCoarse,
    const bool isEvenTimestep
)
{
    // Position Coarse 0., 0., 0.
    const Distributions27 distCoarse = getDistributionReferences27(distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);
    Distributions27 distCoarseAD = getDistributionReferences27(distributionsCoarseAD, numberOfLBnodesCoarse, isEvenTimestep);

    ListIndices indices(indicesCoarse000[nodeIndex], neighborXcoarse, neighborYcoarse, neighborZcoarse);

    const real epsilonNew = vf::basics::constant::c2o1; // ratio of grid resolutions
    const real omegaCoarseNew = hasTurbulentDiffusivity ? vf::lbm::calculateOmegaWithTurbulentViscosity(omegaCoarse, turbulentDiffusivityCoarse[indices.k_000]) : omegaCoarse;
    
    real populationCoarse[27];
    getPreCollisionDistribution(populationCoarse, distCoarse, indices);
    
    real populationCoarseAD[27];
    vf::lbm::ad::interpolateAdvectionDiffusionFC(
        populationCoarse, 
        populationCoarseAD, 
        omegaCoarseNew, 
        epsilonNew,
        coefficients);

    setPreCollisionDistribution(distCoarseAD, indices, populationCoarseAD);
}



//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from fine to coarse
//! \details This scaling function is designed for the advection diffusion F16 kernel
//!
//! The function is executed in the following steps:
//!
template <bool hasTurbulentDiffusivity>
__global__ void scaleFineToCoarseAdvectionDiffusionCompressible_Device(
    real *distributionsCoarse,
    real *distributionsFine,
    real *distributionsCoarseAD,
    real *distributionsFineAD,
    const uint *neighborXcoarse,
    const uint *neighborYcoarse,
    const uint *neighborZcoarse,
    const uint *neighborXfine,
    const uint *neighborYfine,
    const uint *neighborZfine,
    const unsigned long long numberOfLBnodesCoarse,
    const unsigned long long numberOfLBnodesFine,
    const bool isEvenTimestep,
    const uint *indicesCoarse000,
    const uint *indicesFineMMM,
    const uint numberOfInterfaceNodes,
    const real omegaDiffusivityCoarse,
    const real omegaDiffusivityFine,
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborFineToCoarse
)
{
    const unsigned nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= numberOfInterfaceNodes)
        return;

    // 1.calculate moments
    vf::lbm::ad::MomentsOnSourceNodeSet momentsSet;
    ad::calculateMomentSet<hasTurbulentDiffusivity>(
        momentsSet, nodeIndex, distributionsFine, distributionsFineAD, 
        neighborXfine, neighborYfine, neighborZfine, indicesFineMMM, 
        turbulentDiffusivityFine, numberOfLBnodesFine, omegaDiffusivityFine, true);

    // 2.calculate coefficients
    vf::lbm::ad::InterpolationCoefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, neighborFineToCoarse.x[nodeIndex], neighborFineToCoarse.y[nodeIndex], neighborFineToCoarse.z[nodeIndex]);

    // 3. interpolate fine to coarse
    interpolate<hasTurbulentDiffusivity>(
        coefficients,
        nodeIndex,
        distributionsCoarse, 
        distributionsCoarseAD,
        neighborXcoarse,
        neighborYcoarse,
        neighborZcoarse,
        numberOfLBnodesCoarse,
        indicesCoarse000,
        omegaDiffusivityCoarse,
        turbulentDiffusivityCoarse,
        isEvenTimestep);
}

template __global__ void scaleFineToCoarseAdvectionDiffusionCompressible_Device<true>(    
    real *distributionsCoarse,
    real *distributionsFine,
    real *distributionsCoarseAD,
    real *distributionsFineAD,
    const uint *neighborXcoarse,
    const uint *neighborYcoarse,
    const uint *neighborZcoarse,
    const uint *neighborXfine,
    const uint *neighborYfine,
    const uint *neighborZfine,
    const unsigned long long numberOfLBnodesCoarse,
    const unsigned long long numberOfLBnodesFine,
    const bool isEvenTimestep,
    const uint *indicesCoarse000,
    const uint *indicesFineMMM,
    const uint numberOfInterfaceNodes,
    const real omegaDiffusivityCoarse,
    const real omegaDiffusivityFine,
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborFineToCoarse);

template __global__ void scaleFineToCoarseAdvectionDiffusionCompressible_Device<false>(    
    real *distributionsCoarse,
    real *distributionsFine,
    real *distributionsCoarseAD,
    real *distributionsFineAD,
    const uint *neighborXcoarse,
    const uint *neighborYcoarse,
    const uint *neighborZcoarse,
    const uint *neighborXfine,
    const uint *neighborYfine,
    const uint *neighborZfine,
    const unsigned long long numberOfLBnodesCoarse,
    const unsigned long long numberOfLBnodesFine,
    const bool isEvenTimestep,
    const uint *indicesCoarse000,
    const uint *indicesFineMMM,
    const uint numberOfInterfaceNodes,
    const real omegaDiffusivityCoarse,
    const real omegaDiffusivityFine,
    const real* turbulentDiffusivityCoarse,
    const real* turbulentDiffusivityFine,
    const ICellNeigh neighborFineToCoarse);

}
//! \}
