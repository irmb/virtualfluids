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
//! \author Henry Korb
//=======================================================================================
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <cuda_helper/CudaIndexCalculation.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/advectionDiffusion/BoundaryConditions.h>
#include <lbm/constants/D3Q27.h>

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Utilities/KernelUtilities.h"

namespace vf::gpu {

template <BoundaryConditionFactory::AdvectionDiffusionNeumannBC bcType>
__global__ void
AdvectionDiffusionNeumann_Device(real* populationArray, const AdvectionDiffusionNeumannBoundaryConditions bcParameters,
                                 const uint* neighborX, const uint* neighborY, const uint* neighborZ, const real* velocityX,
                                 const real* velocityY, const real* velocityZ, unsigned long long numberOfLBnodes,
                                 const real relaxationFrequency, const bool isEvenTimestep)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::advection_diffusion;
    using namespace vf::lbm::dir;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint k_000 = bcParameters.BCNodeIndices[nodeIndex];
    const ListIndices listIndices(k_000, neighborX, neighborY, neighborZ);

    Distributions27 populationReferences = getDistributionReferences27(populationArray, numberOfLBnodes, isEvenTimestep);

    real populations[NUMBER_Of_DIRECTIONS];
    getPostCollisionDistribution(populations, populationReferences, listIndices);
    const real concentrationNode = vf::lbm::getDensity(populations);
    const real gradient = bcParameters.gradient[nodeIndex];
    const real concentrationWall = concentrationNode - c1o2 * gradient; // wall normal points into the fluid domain

    getPointersToDistributions(populationReferences, populationArray, numberOfLBnodes, !isEvenTimestep);

    switch (bcType) {
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip:
            forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationSimpleAntiBounceBack<direction>(
                    populations, concentrationWall, velocityX[k_000], velocityY[k_000], velocityZ[k_000]);
                writeInInverseDirection<direction>(population, listIndices, populationReferences);
            });
            break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip:
            forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationSimpleAntiBounceBack<direction>(
                    populations, concentrationWall, bcParameters.vx[nodeIndex], bcParameters.vy[nodeIndex],
                    bcParameters.vz[nodeIndex]);
                writeInInverseDirection<direction>(population, listIndices, populationReferences);
            });
            break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedSlip: {
            const real vx1 = velocityX[k_000];
            const real vx2 = velocityY[k_000];
            const real vx3 = velocityZ[k_000];
            forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationInterpolatedAntiBounceBack<direction>(
                    subgridDistance, populations, concentrationNode, concentrationWall, vx1, vx2, vx3, vx1, vx2, vx3,
                    relaxationFrequency);
                writeInInverseDirection<direction>(population, listIndices, populationReferences);
            });
        } break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedNoSlip: {
            forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationInterpolatedAntiBounceBack<direction>(
                    subgridDistance, populations, concentrationNode, concentrationWall, velocityX[k_000], velocityY[k_000],
                    velocityZ[k_000], bcParameters.vx[nodeIndex], bcParameters.vy[nodeIndex], bcParameters.vz[nodeIndex],
                    relaxationFrequency);
                writeInInverseDirection<direction>(population, listIndices, populationReferences);
            });
        } break;
    }
}

}

//! \}
