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
#include <basics/DataTypes.h>
#include <cuda_helper/CudaIndexCalculation.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/advectionDiffusion/BoundaryConditions.h>

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Utilities/KernelUtilities.h"

template <BoundaryConditionFactory::AdvectionDiffusionDirichletBC bcType>
__global__ void AdvectionDiffusionDirichlet_Device(real* populationsArray,
                                                   const AdvectionDiffusionDirichletBoundaryConditions bcParameters,
                                                   const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                   const real* velocityX, const real* velocityY, const real* velocityZ,
                                                   unsigned long long numberOfLBnodes, const real relaxationFrequency,
                                                   const bool isEvenTimestep)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::advection_diffusion;
    using namespace vf::lbm::dir;
    using namespace vf::gpu;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint k_000 = bcParameters.BCNodeIndices[nodeIndex];
    const ListIndices listIndices(k_000, neighborX, neighborY, neighborZ);

    Distributions27 populationReferences = getDistributionReferences27(populationsArray, numberOfLBnodes, isEvenTimestep);

    real populations[27];
    getPostCollisionDistribution(populations, populationReferences, listIndices);
    const real concentrationWall = bcParameters.concentration[nodeIndex];

    getPointersToDistributions(populationReferences, populationsArray, numberOfLBnodes, !isEvenTimestep);

    switch (bcType) {
        case BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackSlip:
            vf::lbm::dir::forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationSimpleAntiBounceBack<direction>(
                    populations, concentrationWall, velocityX[k_000], velocityY[k_000], velocityZ[k_000]);
                writeInInverseDirection<direction>(population, listIndices, populationReferences);
            });
            break;
        case BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip:
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
        case BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletInterpolatedSlip: {
            const real concentrationNode = vf::lbm::getDensity(populations);
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
        case BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletInterpolatedNoSlip: {
            const real concentrationNode = vf::lbm::getDensity(populations);
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