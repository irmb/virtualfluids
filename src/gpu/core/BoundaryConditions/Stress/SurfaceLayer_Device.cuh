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
          bool useDelayedBounceBack>
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
    using SurfaceLayerBC = BoundaryConditionFactory::SurfaceLayerBC;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= boundaryParams.numberOfBCnodes)
        return;

    ///////////////////////////////////////////////////////////
    // Load and compute momentum inputs
    /////////////////////////////////////////////////////////

    Distributions27 populationReferences =
        getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes, gridParams.isEvenTimestep);
    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, boundaryParams.q27[0], boundaryParams.numberOfBCnodes);

    const uint k_000 = boundaryParams.k[nodeIndex];
    const ListIndices listIndices(k_000, gridParams.neighborX, gridParams.neighborY, gridParams.neighborZ);

    real populations[NUMBER_Of_DIRECTIONS];
    getPostCollisionDistribution(populations, populationReferences, listIndices);

    const real drho = vf::lbm::getDensity(populations);
    const real3 velocityNode = { vf::lbm::getCompressibleVelocityX1(populations, drho),
                                 vf::lbm::getCompressibleVelocityX2(populations, drho),
                                 vf::lbm::getCompressibleVelocityX3(populations, drho) };

    const real density = c1o1;

    const real3 wallNormal { boundaryParams.normalX[nodeIndex], boundaryParams.normalY[nodeIndex],
                             boundaryParams.normalZ[nodeIndex] };
    const real3 velocityNodeTangential = computeTangentialVector(velocityNode, wallNormal);

    const real velocityNodeMeanTangentialMagnitude = smoothAndSaveMean(
        computeMagnitude(velocityNodeTangential), filterFrequency, wallModelParams.velocityMagnitudeNode[nodeIndex]);

    const uint samplingIndex = wallModelParams.samplingIndices[nodeIndex];

    const real3 velocitySample = { gridParams.velocityX[samplingIndex], gridParams.velocityY[samplingIndex],
                                   gridParams.velocityZ[samplingIndex] };

    const real velocitySampleTangentialMagnitude = computeMagnitude(computeTangentialVector(velocitySample, wallNormal));

    const real velocitySampleMeanTangentialMagnitude = smoothAndSaveMean(velocitySampleTangentialMagnitude, filterFrequency,
                                                                         wallModelParams.velocityMagnitudeSample[nodeIndex]);

    ///////////////////////////////////////////////////////////
    // Load and compute temperature inputs
    ///////////////////////////////////////////////////////////

    auto populationReferencesTemperature = vf::gpu::getDistributionReferences27(
        temperatureParams.distributionsTemperature, gridParams.numberOfNodes, gridParams.isEvenTimestep);
    real populationsTemperature[NUMBER_Of_DIRECTIONS];
    getPostCollisionDistribution(populationsTemperature, populationReferencesTemperature, listIndices);

    const real temperatureRelativeNode = vf::lbm::getDensity(populationsTemperature);
    temperatureWallModelParams.temperatureNode[nodeIndex] = temperatureRelativeNode;

    const real temperatureRelativeSample = temperatureParams.temperature[samplingIndex];
    const real temperatureRelativeSampleMean = smoothAndSaveMean(temperatureRelativeSample, filterFrequency,
                                                                 temperatureWallModelParams.temperatureSample[nodeIndex]);

    ///////////////////////////////////////////////////////////
    // Load wall model parameters
    ///////////////////////////////////////////////////////////

    const real samplingDistance = wallModelParams.samplingDistance[nodeIndex];
    const real roughnessLength = wallModelParams.roughnessLength[nodeIndex];
    const real vonKarmanConstant = wallModelParams.vonKarmanConstant[nodeIndex];
    const real roughnessLengthTemperature = temperatureWallModelParams.roughnessLength[nodeIndex];

    const real heatingRate = temperatureWallModelParams.heatingRate[nodeIndex];
    const real surfaceTemperature = temperatureWallModelParams.surfaceTemperature[nodeIndex] + heatingRate;
    temperatureWallModelParams.surfaceTemperature[nodeIndex] = surfaceTemperature;
    const real temperatureDifference = temperatureRelativeSampleMean - surfaceTemperature;

    ///////////////////////////////////////////////////////////
    // Compute stress and heat flux from wall model
    ///////////////////////////////////////////////////////////
    real frictionVelocity = wallModelParams.frictionVelocity[nodeIndex];
    real surfaceHeatFlux = temperatureWallModelParams.surfaceHeatFlux[nodeIndex];

    uint iteration = 0;
    bool converged = false;
    do {
        const real stabilityParameter =
            computeStabilityParameter(samplingDistance, temperatureParams.gravity, surfaceHeatFlux, frictionVelocity,
                                      temperatureParams.referenceTemperature, vonKarmanConstant);
        const real stabilityCorrection = computeStabilityCorrectionMomentum(stabilityParameter);
        const real frictionVelocityOld = frictionVelocity;
        frictionVelocity = std::max(zero, computeFrictionVelocity(velocitySampleMeanTangentialMagnitude, vonKarmanConstant,
            samplingDistance, roughnessLength, stabilityCorrection));
        if (heatFluxBCtype == SurfaceLayerBC::SurfaceTemperature)
            surfaceHeatFlux = computeSurfaceHeatFlux(temperatureDifference, frictionVelocity, vonKarmanConstant,
                                                     samplingDistance, roughnessLengthTemperature,
                                                     computeStabilityCorrectionTemperature(stabilityParameter));
        iteration++;
        converged = std::abs(frictionVelocity - frictionVelocityOld) < convergenceCriteria * std::abs(frictionVelocityOld);
    } while ((iteration < maxIter) && !converged);

    wallModelParams.frictionVelocity[nodeIndex] = frictionVelocity;
    if (heatFluxBCtype == SurfaceLayerBC::SurfaceTemperature)
        temperatureWallModelParams.surfaceHeatFlux[nodeIndex] = surfaceHeatFlux;

    const real3 wallShearStress =
        computeWallShearStress(frictionVelocity, velocityNodeTangential, velocityNodeMeanTangentialMagnitude, density);

    ///////////////////////////////////////////////////////////
    // Apply inverse Momentum Exchange
    ///////////////////////////////////////////////////////////

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
    if (!useDelayedBounceBack)
        populationReferences = vf::gpu::getDistributionReferences27(gridParams.distributions, gridParams.numberOfNodes,
                                                                    !gridParams.isEvenTimestep);
    wallMomentum += writeDistributionsBB(populationReferences, linkIsCut, populationsBouncedBack, fakeWallVelocity, density,
                                         listIndices);

    wallModelParams.forceX[nodeIndex] = wallMomentum.x;
    wallModelParams.forceY[nodeIndex] = wallMomentum.y;
    wallModelParams.forceZ[nodeIndex] = wallMomentum.z;

    ///////////////////////////////////////////////////////////
    // Apply Heat Flux Boundary Condition
    ///////////////////////////////////////////////////////////

    const real3 diffusiveFlux = real3 { vf::lbm::getIncompressibleVelocityX1(populationsTemperature),
                                        vf::lbm::getIncompressibleVelocityX2(populationsTemperature),
                                        vf::lbm::getIncompressibleVelocityX3(populationsTemperature) } -
                                velocityNode * temperatureRelativeNode;
    const real normalDiffusiveFlux = dot(diffusiveFlux, wallNormal);
    const real3 wallFlux = diffusiveFlux + wallNormal * (surfaceHeatFlux - normalDiffusiveFlux);

    if (!useDelayedBounceBack)
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