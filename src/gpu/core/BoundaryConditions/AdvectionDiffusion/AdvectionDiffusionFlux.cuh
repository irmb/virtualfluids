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
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Henry Korb
//=======================================================================================

#include <basics/constants/NumericConstants.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/advectionDiffusion/BoundaryConditions.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/constants/D3Q27.h>

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

using FluxBC = BoundaryConditionFactory::AdvectionDiffusionFluxBC;

template <FluxBC fluxBCType>
__global__ void AdvectionDiffusionFlux_Device(real* populationsArray,
                                                      const AdvectionDiffusionFluxBoundaryConditions bcParameters,
                                                      const real* velocityX, const real* velocityY, const real* velocityZ,
                                                      const real* turbulentDiffusivity, real diffusivity,
                                                      const real omegaDiffusivity, const uint* neighborX,
                                                      const uint* neighborY, const uint* neighborZ,
                                                      const unsigned long long numberOfLBnodes, const bool isEvenTimestep)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    using namespace vf::lbm::advection_diffusion;
    using namespace vf::gpu;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    const bool withTurbulentViscosity = fluxBCType == FluxBC::FluxTurbulentViscosityCompressible;

    Distributions27 populationReferences = getDistributionReferences27(populationsArray, numberOfLBnodes, isEvenTimestep);

    const real normalX = bcParameters.normalX[nodeIndex];
    const real normalY = bcParameters.normalY[nodeIndex];
    const real normalZ = bcParameters.normalZ[nodeIndex];

    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint indexOfBCnode = bcParameters.BCNodeIndices[nodeIndex];
    const ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    const real vx1 = velocityX[indexOfBCnode];
    const real vx2 = velocityY[indexOfBCnode];
    const real vx3 = velocityZ[indexOfBCnode];

    real populations[27];
    getPostCollisionDistribution(populations, populationReferences, listIndices);
    getPointersToDistributions(populationReferences, populationsArray, numberOfLBnodes, !isEvenTimestep);

    const real concentration = vf::lbm::getDensity(populations);
    // diffusive flux
    const real diffusiveFluxX = vf::lbm::getIncompressibleVelocityX1(populations) - vx1 * concentration;
    const real diffusiveFluxY = vf::lbm::getIncompressibleVelocityX2(populations) - vx2 * concentration;
    const real diffusiveFluxZ = vf::lbm::getIncompressibleVelocityX3(populations) - vx3 * concentration;

    const real addedDiffusivity = withTurbulentViscosity ? turbulentDiffusivity[indexOfBCnode] : c0o1;
    const real effectiveDiffusivity = addedDiffusivity + diffusivity;

    const real normalFlux = diffusiveFluxX * normalX + diffusiveFluxY * normalY + diffusiveFluxZ * normalZ;

    const real neumannFlux = -bcParameters.gradient[nodeIndex] * effectiveDiffusivity;

    const real fluxX = (diffusiveFluxX - normalFlux * normalX) + neumannFlux * normalX;
    const real fluxY = (diffusiveFluxY - normalFlux * normalY) + neumannFlux * normalY;
    const real fluxZ = (diffusiveFluxZ - normalFlux * normalZ) + neumannFlux * normalZ;

    if (fluxBCType == FluxBC::FluxBounceBack) {
        forEachNonRestDirection([&](auto direction) {
            const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
            if (subgridDistance < c0o1 || subgridDistance > c1o1)
                return;
            const real population = computePopulationSimpleBounceBackWithFlux<direction>(populations, fluxX, fluxY, fluxZ);
            writeInInverseDirection<direction>(population, listIndices, populationReferences);
        });
    } else {
        const real relaxationFrequency = vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivity, addedDiffusivity);
        forEachNonRestDirection([&](auto direction) {
            const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
            if (subgridDistance < c0o1 || subgridDistance > c1o1)
                return;
            const real population = computePopulationInterpolatedBounceBackWithFlux<direction>(
                subgridDistance, populations, vx1, vx2, vx3, relaxationFrequency, concentration, fluxX, fluxY, fluxZ);
            writeInInverseDirection<direction>(population, listIndices, populationReferences);
        });
    }
}

//! \}
