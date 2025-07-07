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
#include "lbm/constants/D3Q27.h"
#include "wallModelMoninObukhov.h"

template <BoundaryConditionFactory::StressBC stressBCType, bool delayed>
__global__ void StressDevice27(GridParameter gridParams, QforBoundaryConditions boundaryParams,
                               WallModelParameters wallModelParams)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    using StressBC = BoundaryConditionFactory::StressBC;

    const real filterFrequency = 1e-3F;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= boundaryParams.numberOfBCnodes)
        return;

    //////////////////////////////////////////////////////////////////////////
    // Load inputs
    //////////////////////////////////////////////////////////////////////////

    Distributions27 populationReferences =
        vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes, gridParams.isEvenTimestep);

    SubgridDistances27 subgridDistances;
    vf::gpu::getPointersToSubgridDistances(subgridDistances, boundaryParams.q27[0], boundaryParams.numberOfBCnodes);

    const uint k_000 = boundaryParams.k[nodeIndex];
    const vf::gpu::ListIndices listIndices(k_000, gridParams.neighborX, gridParams.neighborY, gridParams.neighborZ);

    real populations[NUMBER_Of_DIRECTIONS];
    vf::gpu::getPostCollisionDistribution(populations, populationReferences, listIndices);

    const real drho = vf::lbm::getDensity(populations);
    const real3 velocityNode = { vf::lbm::getCompressibleVelocityX1(populations, drho),
                                 vf::lbm::getCompressibleVelocityX2(populations, drho),
                                 vf::lbm::getCompressibleVelocityX3(populations, drho) };
    const real density = c1o1;

    const real3 wallNormal { boundaryParams.normalX[nodeIndex], boundaryParams.normalY[nodeIndex],
                             boundaryParams.normalZ[nodeIndex] };

    const real3 velocityNodeTangential = computeTangentialVector(velocityNode, wallNormal);

    const real velocityNodeMeanTangentialMagnitude = smoothAndSaveMean(computeMagnitude(velocityNodeTangential), filterFrequency,
                                                                       wallModelParams.velocityMagnitudeNode[nodeIndex]);

    const uint samplingIndex = wallModelParams.samplingIndices[nodeIndex];
    const real3 velocitySample { gridParams.velocityX[samplingIndex], gridParams.velocityY[samplingIndex],
                                 gridParams.velocityZ[samplingIndex] };

    const real velocitySampleTangentialMagnitude = computeMagnitude(computeTangentialVector(velocitySample, wallNormal));

    const real velocitySampleMeanTangentialMagnitude = smoothAndSaveMean(velocitySampleTangentialMagnitude, filterFrequency,
                                                                         wallModelParams.velocityMagnitudeSample[nodeIndex]);

    //////////////////////////////////////////////////////////////////////////
    // load wall model parameters and compute wall shear stress
    //////////////////////////////////////////////////////////////////////////

    const real samplingDistance = wallModelParams.samplingDistance[nodeIndex];
    const real roughnessLength = wallModelParams.roughnessLength[nodeIndex];
    const real vonKarmanConstant = wallModelParams.vonKarmanConstant[nodeIndex];
    const real frictionVelocity = computeFrictionVelocity(velocitySampleMeanTangentialMagnitude, vonKarmanConstant,
                                                          samplingDistance, roughnessLength, c0o1);
    const real3 wallShearStress =
        computeWallShearStress(frictionVelocity, velocityNodeTangential, velocityNodeMeanTangentialMagnitude, density);

    //////////////////////////////////////////////////////////////////////////
    // apply inverse Momentum exchange
    //////////////////////////////////////////////////////////////////////////

    real populationsBouncedBack[NUMBER_Of_DIRECTIONS];
    bool linkIsCut[NUMBER_Of_DIRECTIONS];

    switch (stressBCType) {
        case StressBC::StressBounceBackCompressible:
            computeBouncedBackDistributionsBB(subgridDistances, populations, linkIsCut, populationsBouncedBack, nodeIndex);
            break;
        case StressBC::StressBounceBackWithPressureCompressible:
            computeBouncedBackDistributionsBBPressure(subgridDistances, drho, populations, linkIsCut, populationsBouncedBack,
                                                      nodeIndex);
            break;
        case StressBC::StressInterpolatedCompressible: {
            const real relaxationFrequency = vf::lbm::calculateOmegaWithTurbulentViscosity(
                gridParams.relaxationFrequency, gridParams.turbulentViscosity[k_000]);
            computeBouncedBackDistributionsInterpolated(subgridDistances, velocityNode, drho, relaxationFrequency,
                                                        populations, linkIsCut, populationsBouncedBack, nodeIndex);
        } break;
    }

    real3 wallMomentum = computeWallMomentumBounceBack(linkIsCut, populationsBouncedBack, populations);

    const real wallArea = c1o1;
    const real subgridDistance = (subgridDistances.q[d00M])[nodeIndex];
    const real interpolationFactor =
        stressBCType == StressBC::StressInterpolatedCompressible ? c1o1 + subgridDistance : c1o1;

    const real3 fakeWallVelocity = computeFakeWallVelocity(wallNormal, velocitySample, wallShearStress, density,
                                                           interpolationFactor, wallArea, wallMomentum);
    if (!delayed) {
        populationReferences = vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes,
                                                                    !gridParams.isEvenTimestep);
    };

    if (stressBCType == StressBC::StressInterpolatedCompressible)
        wallMomentum +=
            writeDistributionsInterpolatedBB(populationReferences, linkIsCut, populationsBouncedBack, fakeWallVelocity,
                                             density, subgridDistances, listIndices, nodeIndex);
    else
        wallMomentum += writeDistributionsBB(populationReferences, linkIsCut, populationsBouncedBack, fakeWallVelocity,
                                             density, listIndices);

    if (wallModelParams.hasMonitor) {
        wallModelParams.frictionVelocity[nodeIndex] = frictionVelocity;
        wallModelParams.forceX[nodeIndex] = wallMomentum.x;
        wallModelParams.forceY[nodeIndex] = wallMomentum.y;
        wallModelParams.forceZ[nodeIndex] = wallMomentum.z;
    }
}
#endif

//! \}
