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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef Stress_Device_H
#define Stress_Device_H

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <cuda_helper/CudaIndexCalculation.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/collision/TurbulentViscosity.h>

#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"

#include "Stress.h"
#include "inverseMomentumExchange.cuh"
#include "wallModelMoninObukhov.h"

template <BoundaryConditionFactory::StressBC stressBCType>
__global__ void StressDevice27(GridParameter gridParams, BoundaryParameter boundaryParams,
                               WallModelParameters wallModelParams)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    using StressBC = BoundaryConditionFactory::StressBC;

    const real filterFrequency = 1e-3;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= boundaryParams.numberOfBCNodes)
        return;

    //////////////////////////////////////////////////////////////////////////
    // Load inputs
    //////////////////////////////////////////////////////////////////////////

    const Distributions27 populationReferences =
        vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes, gridParams.isEvenTimestep);

    SubgridDistances27 subgridDistances;
    vf::gpu::getPointersToSubgridDistances(subgridDistances, boundaryParams.subgridDistances,
                                           boundaryParams.numberOfBCNodes);

    const uint k_000 = boundaryParams.indices[nodeIndex];
    const vf::gpu::ListIndices listIndices(k_000, gridParams.neighborX, gridParams.neighborY, gridParams.neighborZ);

    real populations[27];
    vf::gpu::getPostCollisionDistribution(populations, populationReferences, listIndices);

    real drho;
    real3 velocityNode;
    vf::lbm::getCompressibleMacroscopicValues(populations, drho, velocityNode.x, velocityNode.y,
                                              velocityNode.z);
    const real density = c1o1;

    const real3 wallNormal { boundaryParams.normalX[nodeIndex], boundaryParams.normalY[nodeIndex],
                             boundaryParams.normalZ[nodeIndex] };

    const real3 velocityNodeWallParallel = computeWallParallelVector(velocityNode, wallNormal);

    const real3 velocityNodeMean {
        smoothAndSaveMean(velocityNode.x, filterFrequency, wallModelParams.velocityNodeX[nodeIndex]),
        smoothAndSaveMean(velocityNode.y, filterFrequency, wallModelParams.velocityNodeY[nodeIndex]),
        smoothAndSaveMean(velocityNode.z, filterFrequency, wallModelParams.velocityNodeZ[nodeIndex])
    };

    const real velocityNodeMeanWallParallelMagnitude = computeWallParallelVelocityMagnitude(velocityNodeMean, wallNormal);

    const uint samplingIndex = wallModelParams.samplingIndices[nodeIndex];
    const real3 velocityExchangeLocation { gridParams.velocityX[samplingIndex], gridParams.velocityY[samplingIndex],
                                           gridParams.velocityZ[samplingIndex] };

    const real3 velocityExchangeLocationMean {
        smoothAndSaveMean(velocityExchangeLocation.x, filterFrequency, wallModelParams.velocityExchangeLocationX[nodeIndex]),
        smoothAndSaveMean(velocityExchangeLocation.y, filterFrequency, wallModelParams.velocityExchangeLocationY[nodeIndex]),
        smoothAndSaveMean(velocityExchangeLocation.z, filterFrequency, wallModelParams.velocityExchangeLocationZ[nodeIndex])
    };

    const real velocityExchangeLocationMeanWallParallelMagnitude =
        computeWallParallelVelocityMagnitude(velocityExchangeLocationMean, wallNormal);

    //////////////////////////////////////////////////////////////////////////
    // load wall model parameters and compute wall shear stress
    //////////////////////////////////////////////////////////////////////////

    const real samplingDistance = wallModelParams.samplingDistance[nodeIndex];
    const real roughnessLength = wallModelParams.roughnessLength[nodeIndex];
    const real vonKarmanConstant = wallModelParams.vonKarmanConstant[nodeIndex];
    const real frictionVelocity =
        computeFrictionVelocity(velocityExchangeLocationMeanWallParallelMagnitude, vonKarmanConstant, samplingDistance, roughnessLength, c0o1);
    const real3 wallShearStress =
        computeWallShearStress(frictionVelocity, velocityNodeWallParallel, velocityNodeMeanWallParallelMagnitude, density);

    //////////////////////////////////////////////////////////////////////////
    // apply inverse Momentum exchange
    //////////////////////////////////////////////////////////////////////////

    real populationsBouncedBack[27];
    bool linkIsCut[27];

    switch (stressBCType) {
        case StressBC::StressBounceBackCompressible:
            computeBouncedBackDistributionsBB(subgridDistances, populations, linkIsCut,
                                              populationsBouncedBack, nodeIndex);
            break;
        case StressBC::StressBounceBackPressureCompressible:
            computeBouncedBackDistributionsBBPressure(subgridDistances, drho, populations, linkIsCut,
                                                      populationsBouncedBack, nodeIndex);
            break;
        case StressBC::StressCompressible: {
            const real relaxationFrequency = vf::lbm::calculateOmegaWithTurbulentViscosity(
                gridParams.relaxationFrequency, gridParams.turbulentViscosity[k_000]);
            computeBouncedBackDistributionsInterpolated(subgridDistances, velocityNode, drho, relaxationFrequency,
                                                        populations, linkIsCut, populationsBouncedBack,
                                                        nodeIndex);
        } break;
    }

    const real3 wallMomentumBounceBack =
        computeWallMomentumBounceBack(linkIsCut, populationsBouncedBack, populations);

    const real wallArea = c1o1;
    const real subgridDistance = (subgridDistances.q[d00M])[nodeIndex];
    const real interpolationFactor = stressBCType == StressBC::StressCompressible ? c1o1 + subgridDistance : c1o1;

    const real3 fakeWallVelocity = computeFakeWallVelocity(wallNormal, velocityExchangeLocation, wallShearStress, density,
                                                           interpolationFactor, wallArea, wallMomentumBounceBack);

    const real3 wallMomentumVelocity = writeDistributionsBB(populationReferences, linkIsCut, populationsBouncedBack,
                                                          fakeWallVelocity, density, interpolationFactor, listIndices);
    // const Distributions27 outgoingDistributions =
    //     vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes, !gridParams.isEvenTimestep);

    // const real3 wallMomentumVelocity = writeDistributionsHalfWayBB(outgoingDistributions, linkIsCut, distributionsBouncedBack,
    //                                                       fakeWallVelocity, density, interpolationFactor, listIndices);

    if (wallModelParams.hasMonitor) {
        wallModelParams.frictionVelocity[nodeIndex] = frictionVelocity;
        wallModelParams.forceX[nodeIndex] = wallMomentumBounceBack.x + wallMomentumVelocity.x;
        wallModelParams.forceY[nodeIndex] = wallMomentumBounceBack.y + wallMomentumVelocity.y;
        wallModelParams.forceZ[nodeIndex] = wallMomentumBounceBack.z + wallMomentumVelocity.z;
    }
}
#endif

//! \}
