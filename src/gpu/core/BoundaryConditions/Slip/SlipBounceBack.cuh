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
//======================================================================================
#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"
#include "basics/constants/NumericConstants.h"
#include "cuda_helper/CudaIndexCalculation.h"
#include "lbm/MacroscopicQuantities.h"
#include "lbm/constants/D3Q27.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void SlipBounceBack_Device(real* populationsArray, QforBoundaryConditions bcParams,
                                      const uint* neighborX, const uint* neighborY,
                                      const uint* neighborZ, const unsigned long long numberOfLBnodes,
                                      const bool isEvenTimestep)
{

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= bcParams.numberOfBCnodes)
        return;
    Distributions27 populationReferences;
    getPointersToDistributions(populationReferences, populationsArray, numberOfLBnodes, isEvenTimestep);

    SubgridDistances27 subgridD;
    getPointersToSubgridDistances(subgridD, bcParams.q27[0], bcParams.numberOfBCnodes);

    const uint indexOfBCnode = bcParams.k[nodeIndex];
    const ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    real populations[27];
    getPostCollisionDistribution(populations, populationReferences, listIndices);

    const real drho = vf::lbm::getDensity(populations);
    const real3 velocity { vf::lbm::getCompressibleVelocityX1(populations, drho),
                           vf::lbm::getCompressibleVelocityX2(populations, drho),
                           vf::lbm::getCompressibleVelocityX3(populations, drho) };

    const real3 normal { bcParams.normalX[nodeIndex], bcParams.normalY[nodeIndex], bcParams.normalZ[nodeIndex] };

    const real3 velocityTangential = velocity - normal * dot(normal, velocity);

    // getPointersToDistributions(populationReferences, populationsArray, numberOfLBnodes, !isEvenTimestep);

    forEachNonRestDirection([&](auto dir) {
        const real subgridDistance = (subgridD.q[dir])[nodeIndex];
        if (subgridDistance > c1o1 || subgridDistance < c0o1)
            return;
        const real weight = getWeight<dir>();
        const real velocity = getVelocity<dir>(velocityTangential.x, velocityTangential.y, velocityTangential.z);
        const real population = getBounceBackDistributionForVeloBC(populations[dir], velocity, weight);
        vf::gpu::writeInInverseDirection<dir>(population, listIndices, populationReferences);
    });
}