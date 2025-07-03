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
#ifndef StressFunctions_H
#define StressFunctions_H

#include <cmath>
#include <cuda_runtime.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/constants/D3Q27.h>

#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"

constexpr real smoothAndSaveMean(real instantaneous, real filterFrequency, real& mean)
{
    const real smoothed = filterFrequency * instantaneous + (vf::basics::constant::c1o1 - filterFrequency) * mean;
    mean = smoothed;
    return smoothed;
}

constexpr real3 computeWallParallelVector(real3 quantity, real3 normal)
{
    return quantity - normal * dot(quantity, normal);
}

inline __device__ real computeWallParallelVelocityMagnitude(real3 quantity, real3 normal)
{
    return std::sqrt(square(computeWallParallelVector(quantity, normal)));
}

constexpr void findCutLinks(bool* cutQs, const SubgridDistances27& subgridD, const uint nodeIndex)
{
    forEachNonRestDirection([&](auto dir) {
        const real subgridDistance = (subgridD.q[dir])[nodeIndex];
        cutQs[dir] = subgridDistance >= c0o1 && subgridDistance <= c1o1;
    });
}

constexpr void computeBouncedBackDistributionsBB(const SubgridDistances27& subgridD, const real* distributionsPostCollision,
                                                 bool* linkIsCut, real* distributionsBouncedBack, const uint nodeIndex)
{
    findCutLinks(linkIsCut, subgridD, nodeIndex);

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = vf::lbm::dir::inverseDir<dir>();
        distributionsBouncedBack[inverseDirection] = distributionsPostCollision[dir];
    });
}

constexpr void computeBouncedBackDistributionsBBPressure(const SubgridDistances27& subgridD, const real drho,
                                                         const real* distributionsPostCollision, bool* linkIsCut,
                                                         real* distributionsBouncedBack, const uint nodeIndex)
{
    findCutLinks(linkIsCut, subgridD, nodeIndex);

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = vf::lbm::dir::inverseDir<dir>();
        const real weight = vf::lbm::dir::getWeight<dir>();
        distributionsBouncedBack[inverseDirection] = distributionsPostCollision[dir] - weight * drho;
    });
}

constexpr void computeBouncedBackDistributionsInterpolated(const SubgridDistances27& subgridD, real3 velocity, real drho,
                                                           real relaxationFrequency, const real* distributionsPostCollision,
                                                           bool* linkIsCut, real* distributionsBouncedBack, uint nodeIndex)
{
    using namespace vf::lbm::dir;
    forEachNonRestDirection([&](auto dir) {
        using namespace vf::lbm::dir;
        const real subgridDistance = (subgridD.q[dir])[nodeIndex];
        if (subgridDistance < c0o1 || subgridDistance > c1o1) {
            linkIsCut[dir] = false;
            return;
        }

        linkIsCut[dir] = true;
        const real weight = getWeight<dir>();
        const size_t inverseDirection = inverseDir<dir>();
        const real cu = getVelocity<dir>(velocity.x, velocity.y, velocity.z);
        const real cu_sq =
            c3o2 * (velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z) * (c1o1 + drho);
        const real populationPostCollision = distributionsPostCollision[dir];
        const real feq = vf::gpu::getEquilibriumForBC(drho, cu, cu_sq, weight);
        const real fInverse = distributionsPostCollision[inverseDirection];

        const real fBouncedBack = vf::gpu::getInterpolatedDistributionForNoSlipWithPressureBC(
            subgridDistance, populationPostCollision, fInverse, feq, relaxationFrequency, drho, weight);

        distributionsBouncedBack[inverseDirection] = fBouncedBack;
    });
}

template <size_t direction>
constexpr void addWallMomentum(const bool* linkIsCut, const real* distributionsBouncedBack,
                               const real* distributionsPostCollision, real3& wallMomentum)
{
}

constexpr real3 computeWallMomentumBounceBack(const bool* linkIsCut, const real* distributionsBouncedBack,
                                              const real* distributionsPostCollision)
{
    using namespace vf::lbm::dir;
    real3 wallMomentum {};
    forEachNonRestDirection([&](auto dir) {
        using namespace vf::lbm::dir;
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = inverseDir<dir>();
        const real momentum = distributionsPostCollision[dir] + distributionsBouncedBack[inverseDirection];

        wallMomentum.x += momentum * getComponentX<dir>();
        wallMomentum.y += momentum * getComponentY<dir>();
        wallMomentum.z += momentum * getComponentZ<dir>();
    });

    return wallMomentum;
}

constexpr real3 computeFakeWallVelocity(real3 wallNormal, real3 velocityForClipping, real3 wallShearStress, real density,
                                        real interpolationFactor, real wallArea, real3 wallMomentum)
{

    const real3 wallModelForce = wallShearStress * wallArea;
    const real3 wallParallelMomentum = computeWallParallelVector(wallMomentum, wallNormal);

    const real3 force = wallModelForce - wallParallelMomentum;

    // Compute  wall velocity and clip (clipping only necessary for initial boundary layer development)
    constexpr real clipWallVelo = c2o1;

    const real3 clipVelocity { std::abs(clipWallVelo * velocityForClipping.x),
                               std::abs(clipWallVelo * velocityForClipping.y),
                               std::abs(clipWallVelo * velocityForClipping.z) };

    return { std::clamp(-c3o1 * force.x * interpolationFactor / density, -clipVelocity.x, clipVelocity.x),
             std::clamp(-c3o1 * force.y * interpolationFactor / density, -clipVelocity.y, clipVelocity.y),
             std::clamp(-c1o1 * force.z * interpolationFactor / density, -clipVelocity.z, clipVelocity.z) };
}

constexpr real3 writeDistributionsBB(const Distributions27& distributionReferences, const bool* linkIsCut,
                                     const real* distributionsBouncedBack, real3 velocity, real density,
                                     real interpolationFactor, const vf::gpu::ListIndices& listIndices)
{
    using namespace vf::lbm::dir;

    velocity /= interpolationFactor;

    real3 wallMomentumAdded {};

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = inverseDir<dir>();
        const real addedMomentum = -c6o1 * density * getWeight<dir>() * getVelocity<dir>(velocity.x, velocity.y, velocity.z);
        const auto writeIndex = listIndices.getIndex<inverseDirection>();

        (distributionReferences.f[inverseDirection])[writeIndex] =
            distributionsBouncedBack[inverseDirection] + addedMomentum;

        wallMomentumAdded.x += addedMomentum * getComponentX<dir>();
        wallMomentumAdded.y += addedMomentum * getComponentY<dir>();
        wallMomentumAdded.z += addedMomentum * getComponentZ<dir>();
    });
    return wallMomentumAdded;
}

#endif
