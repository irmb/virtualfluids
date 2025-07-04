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
//! \author Henry Korb
//=======================================================================================
#ifndef SurfaceLayer_Device_H
#define SurfaceLayer_Device_H

#include <cmath>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <cuda_helper/CudaIndexCalculation.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/advectionDiffusion/BoundaryConditions.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/constants/D3Q27.h>

#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Calculation/Calculation.h"
#include "Stress.h"
#include "SurfaceLayer.h"
#include "Utilities/KernelUtilities.h"
#include "inverseMomentumExchange.cuh"
#include "wallModelMoninObukhov.h"

template <BoundaryConditionFactory::StressBC stressBCType, BoundaryConditionFactory::SurfaceLayerBC heatFluxBCtype,
          bool delayed>
__global__ void
SurfaceLayerDevice27(GridParameter gridParams, QforBoundaryConditions boundaryParams, WallModelParameters wallModelParams,
                     TemperatureWallModelParameters temperatureWallModelParams, TemperatureParameters temperatureParams)
{
    constexpr real filterFrequency = 1e-3F;
    constexpr uint maxIter = 100;
    constexpr real convergenceCriteria = 1e-3F;
    constexpr real zero = c0o1; // need for some std functions
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    using namespace vf::gpu;
    using namespace vf::lbm::advection_diffusion;

    using StressBC = BoundaryConditionFactory::StressBC;
    using HeatFluxBC = BoundaryConditionFactory::SurfaceLayerBC;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= boundaryParams.numberOfBCnodes)
        return;

    ///////////////////////////////////////////////////////////
    // Load and compute momentum inputs
    /////////////////////////////////////////////////////////

    Distributions27 distributionReferences =
        getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes, gridParams.isEvenTimestep);
    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, boundaryParams.q27[0], boundaryParams.numberOfBCnodes);

    const uint k_000 = boundaryParams.k[nodeIndex];
    const ListIndices listIndices(k_000, gridParams.neighborX, gridParams.neighborY, gridParams.neighborZ);

    real populationsPostCollision[27];
    getPostCollisionDistribution(populationsPostCollision, distributionReferences, listIndices);

    real drho;
    real3 velocityNode;

    vf::lbm::getCompressibleMacroscopicValues(populationsPostCollision, drho, velocityNode.x, velocityNode.y,
                                              velocityNode.z);

    const real density = c1o1;

    const real3 wallNormal { boundaryParams.normalX[nodeIndex], boundaryParams.normalY[nodeIndex],
                             boundaryParams.normalZ[nodeIndex] };
    const real3 velocityNodeWallParallel = computeWallParallelVector(velocityNode, wallNormal);

    const real velocityNodeWallParallelMagnitude = computeWallParallelVelocityMagnitude(velocityNode, wallNormal);
    const real velocityNodeMeanWallParallelMagnitude =
        smoothAndSaveMean(velocityNodeWallParallelMagnitude, filterFrequency, wallModelParams.velocityNodeX[nodeIndex]);

    const uint samplingIndex = wallModelParams.samplingIndices[nodeIndex];

    const real3 velocityExchangeLocationMean = { smoothAndSaveMean(gridParams.velocityX[samplingIndex], filterFrequency,
                                                                   wallModelParams.velocityExchangeLocationX[nodeIndex]),
                                                 smoothAndSaveMean(gridParams.velocityY[samplingIndex], filterFrequency,
                                                                   wallModelParams.velocityExchangeLocationY[nodeIndex]),
                                                 smoothAndSaveMean(gridParams.velocityZ[samplingIndex], filterFrequency,
                                                                   wallModelParams.velocityExchangeLocationZ[nodeIndex]) };

    const real velocityExchangeLocationMeanWallParallelMagnitude =
        computeWallParallelVelocityMagnitude(velocityExchangeLocationMean, wallNormal);

    ///////////////////////////////////////////////////////////
    // Load and compute temperature inputs
    ///////////////////////////////////////////////////////////

    const real temperatureRelativeNode = temperatureParams.temperature[k_000];
    temperatureWallModelParams.temperatureFirstFluidNode[nodeIndex] = temperatureRelativeNode;

    const real temperatureRelativeExchangeLocation = temperatureParams.temperature[samplingIndex];
    const real temperatureRelativeExchangeLocationMean =
        smoothAndSaveMean(temperatureRelativeExchangeLocation, filterFrequency,
                          temperatureWallModelParams.temperatureExchangeLocation[nodeIndex]);

    ///////////////////////////////////////////////////////////
    // Load wall model parameters
    ///////////////////////////////////////////////////////////

    const real samplingDistance = wallModelParams.samplingDistance[nodeIndex];
    const real roughnessLength = wallModelParams.roughnessLength[nodeIndex];
    const real vonKarmanConstant = wallModelParams.vonKarmanConstant[nodeIndex];
    const real roughnessLengthTemperature = temperatureWallModelParams.roughnessLength[nodeIndex];

    const real heatingRate = temperatureWallModelParams.heatingRate[nodeIndex];
    const real surfaceTemperature = temperatureWallModelParams.surfaceTemperature[nodeIndex] + heatingRate;
    const real temperatureDifference = temperatureRelativeExchangeLocationMean - surfaceTemperature;

    ///////////////////////////////////////////////////////////
    // Compute stress and heat flux from wall model
    ///////////////////////////////////////////////////////////
    real frictionVelocity = computeFrictionVelocity(velocityExchangeLocationMeanWallParallelMagnitude, vonKarmanConstant,
                                                    samplingDistance, roughnessLength, c0o1);
    real temperatureScale = heatFluxBCtype == HeatFluxBC::SurfaceHeatFlux
                                ? -temperatureWallModelParams.surfaceHeatFlux[nodeIndex] / frictionVelocity
                                : computeFrictionVelocity(temperatureDifference, vonKarmanConstant, samplingDistance,
                                                          roughnessLengthTemperature, c0o1);
    const real stabilityParameter =
        computeStabilityParameter(samplingDistance, temperatureParams.gravity, temperatureScale, frictionVelocity,
                                  temperatureParams.referenceTemperature, vonKarmanConstant);
    if (stabilityParameter > c0o1) {
        const real stabilityParameter =
            heatFluxBCtype == HeatFluxBC::SurfaceTemperature
                ? computeStabilityParameterFromSurfaceTemperature(
                      samplingDistance, roughnessLength, roughnessLengthTemperature,
                      velocityExchangeLocationMeanWallParallelMagnitude, temperatureDifference,
                      temperatureParams.referenceTemperature, temperatureParams.gravity)
                : computeStabilityParameterFromHeatFlux(
                      samplingDistance, roughnessLength, velocityExchangeLocationMeanWallParallelMagnitude,
                      temperatureWallModelParams.surfaceHeatFlux[nodeIndex], temperatureParams.referenceTemperature,
                      temperatureParams.gravity, vonKarmanConstant);

        const auto stabilityCorrections = computeStabilityCorrectionsStable(stabilityParameter);
        frictionVelocity = computeFrictionVelocity(velocityExchangeLocationMeanWallParallelMagnitude, vonKarmanConstant,
                                                   samplingDistance, roughnessLength, stabilityCorrections.momentum);
        temperatureScale = computeFrictionVelocity(temperatureDifference, vonKarmanConstant, samplingDistance,
                                                   roughnessLengthTemperature, stabilityCorrections.temperature);
    } else if (stabilityParameter < c0o1) {
        bool converged = false;
        uint iteration = 0;
        do {
            real stabilityParameter =
                computeStabilityParameter(samplingDistance, temperatureParams.gravity, temperatureScale, frictionVelocity,
                                          temperatureParams.referenceTemperature, vonKarmanConstant);
            const auto stabilityCorrections = computeStabilityCorrectionsUnstable(stabilityParameter);

            const real newFrictionVelocity =
                std::max(zero, computeFrictionVelocity(velocityExchangeLocationMeanWallParallelMagnitude, vonKarmanConstant,
                                                       samplingDistance, roughnessLength, stabilityCorrections.momentum));
            const real newTemperatureScale =
                computeFrictionVelocity(temperatureDifference, vonKarmanConstant, samplingDistance,
                                        roughnessLengthTemperature, stabilityCorrections.temperature);

            iteration++;
            const bool velocityConverged =
                std::abs(newFrictionVelocity - frictionVelocity) < convergenceCriteria * std::abs(frictionVelocity);
            const bool temperatureConverged =
                std::abs(newTemperatureScale - temperatureScale) < convergenceCriteria * std::abs(temperatureScale);
            converged = velocityConverged && temperatureConverged;
            temperatureScale = newTemperatureScale;
            frictionVelocity = newFrictionVelocity;
        } while ((iteration < maxIter) && !converged);
    }

    const real surfaceHeatFlux = -temperatureScale * frictionVelocity;
    wallModelParams.frictionVelocity[nodeIndex] = frictionVelocity;
    temperatureWallModelParams.temperatureScale[nodeIndex] = temperatureScale;
    temperatureWallModelParams.surfaceTemperature[nodeIndex] = surfaceTemperature;

    if (heatFluxBCtype == HeatFluxBC::SurfaceTemperature)
        temperatureWallModelParams.surfaceHeatFlux[nodeIndex] = surfaceHeatFlux;

    const real3 wallShearStress =
        computeWallShearStress(frictionVelocity, velocityNodeWallParallel, velocityNodeMeanWallParallelMagnitude, density);

    ///////////////////////////////////////////////////////////
    // Apply inverse Momentum Exchange
    ///////////////////////////////////////////////////////////

    real populationsBouncedBack[27];
    bool linkIsCut[27];

    switch (stressBCType) {
        case StressBC::StressBounceBackCompressible:
            computeBouncedBackDistributionsBB(subgridDistances, populationsPostCollision, linkIsCut, populationsBouncedBack,
                                              nodeIndex);
            break;
        case StressBC::StressBounceBackWithPressureCompressible:
            computeBouncedBackDistributionsBBPressure(subgridDistances, drho, populationsPostCollision, linkIsCut,
                                                      populationsBouncedBack, nodeIndex);
            break;
        case StressBC::StressInterpolatedCompressible: {
            const real relaxationFrequency = vf::lbm::calculateOmegaWithTurbulentViscosity(
                gridParams.relaxationFrequency, gridParams.turbulentViscosity[k_000]);
            computeBouncedBackDistributionsInterpolated(subgridDistances, velocityNode, drho, relaxationFrequency,
                                                        populationsPostCollision, linkIsCut, populationsBouncedBack,
                                                        nodeIndex);
        } break;
    }

    const real3 wallMomentumBounceBack =
        computeWallMomentumBounceBack(linkIsCut, populationsBouncedBack, populationsPostCollision);

    const real wallArea = c1o1;
    const real subgridDistance = (subgridDistances.q[d00M])[nodeIndex];
    const real interpolationFactor =
        stressBCType == StressBC::StressInterpolatedCompressible ? c1o1 + subgridDistance : c1o1;

    const real3 fakeWallVelocity = computeFakeWallVelocity(wallNormal, velocityExchangeLocationMean, wallShearStress,
                                                           density, interpolationFactor, wallArea, wallMomentumBounceBack);
    if (!delayed)
        distributionReferences = vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes,
                                                                      !gridParams.isEvenTimestep);
    const real3 wallMomentumVelocity = writeDistributionsBB(distributionReferences, linkIsCut, populationsBouncedBack,
                                                            fakeWallVelocity, density, listIndices);

    if (wallModelParams.hasMonitor) {
        wallModelParams.forceX[nodeIndex] = wallMomentumBounceBack.x + wallMomentumVelocity.x;
        wallModelParams.forceY[nodeIndex] = wallMomentumBounceBack.y + wallMomentumVelocity.y;
        wallModelParams.forceZ[nodeIndex] = wallMomentumBounceBack.z + wallMomentumVelocity.z;
    }

    ///////////////////////////////////////////////////////////
    // Apply Heat Flux Boundary Condition
    ///////////////////////////////////////////////////////////
    auto populationReferencesTemperature = vf::gpu::getDistributionReferences27(
        temperatureParams.distributionsTemperature, gridParams.numberOfNodes, gridParams.isEvenTimestep);
    real populationsTemperature[27];
    getPostCollisionDistribution(populationsTemperature, populationReferencesTemperature, listIndices);

    const real3 diffusiveFlux = real3 { vf::lbm::getIncompressibleVelocityX1(populationsTemperature),
                                        vf::lbm::getIncompressibleVelocityX2(populationsTemperature),
                                        vf::lbm::getIncompressibleVelocityX3(populationsTemperature) } -
                                velocityNode * temperatureRelativeNode;
    const real normalDiffusiveFlux = dot(diffusiveFlux, wallNormal);
    const real3 wallFlux = diffusiveFlux + wallNormal * (surfaceHeatFlux - normalDiffusiveFlux);

    if (!delayed)
        populationReferencesTemperature = vf::gpu::getDistributionReferences27(
            temperatureParams.distributionsTemperature, gridParams.numberOfNodes, !gridParams.isEvenTimestep);

    switch (stressBCType) {
        case StressBC::StressBounceBackCompressible:
        case StressBC::StressBounceBackWithPressureCompressible:
            forEachNonRestDirection([&](auto direction) {
                if (!linkIsCut[direction])
                    return;
                const real population = computePopulationSimpleBounceBackWithFlux<direction>(
                    populationsTemperature, wallFlux.x, wallFlux.y, wallFlux.z);
                writeInInverseDirection<direction>(population, listIndices, populationReferencesTemperature);
            });
            break;
        case StressBC::StressInterpolatedCompressible:
            const real diffusivityNode = temperatureParams.diffusivity + temperatureParams.turbulentDiffusivity[k_000];
            const real relaxationFrequency = vf::lbm::computeRelaxationFrequency(diffusivityNode);
            forEachNonRestDirection([&](auto direction) {
                const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
                if (subgridDistance < c0o1 || subgridDistance > c1o1)
                    return;
                const real population = computePopulationInterpolatedBounceBackWithFlux<direction>(
                    subgridDistance, populationsTemperature, velocityNode.x, velocityNode.y, velocityNode.z,
                    relaxationFrequency, temperatureRelativeNode, wallFlux.x, wallFlux.y, wallFlux.z);
                writeInInverseDirection<direction>(population, listIndices, populationReferencesTemperature);
            });
            break;
    }
}

#endif