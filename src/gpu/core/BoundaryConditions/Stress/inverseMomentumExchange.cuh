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

constexpr void findCutLinks(bool* const cutQs, const SubgridDistances27& subgridD, const uint nodeIndex)
{
    forEachNonRestDirection([&](auto dir) {
        const real subgridDistance = (subgridD.q[dir])[nodeIndex];
        cutQs[dir] = subgridDistance >= c0o1 && subgridDistance <= c1o1;
    });
}

constexpr void computeBouncedBackDistributionsBB(const SubgridDistances27& subgridD, const real* populations,
                                                 bool* const linkIsCut, real* populationsBouncedBack, const uint nodeIndex)
{
    findCutLinks(linkIsCut, subgridD, nodeIndex);

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = vf::lbm::dir::inverseDir<dir>();
        populationsBouncedBack[inverseDirection] = populations[dir];
    });
}

constexpr void computeBouncedBackDistributionsBBPressure(const SubgridDistances27& subgridD, const real drho,
                                                         const real* populations, bool* const linkIsCut,
                                                         real* populationsBouncedBack, const uint nodeIndex)
{
    findCutLinks(linkIsCut, subgridD, nodeIndex);

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = vf::lbm::dir::inverseDir<dir>();
        const real weight = vf::lbm::dir::getWeight<dir>();
        populationsBouncedBack[inverseDirection] = populations[dir] - weight * drho;
    });
}

constexpr void computeBouncedBackDistributionsInterpolated(const SubgridDistances27& subgridD, real3 velocity, real drho,
                                                           real relaxationFrequency, const real* populations,
                                                           bool* const linkIsCut, real* populationsBouncedBack,
                                                           uint nodeIndex)
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
        const real cu_sq = c3o2 * square(velocity) * (c1o1 + drho);
        const real feq = vf::gpu::getEquilibriumForBC(drho, cu, cu_sq, weight);

        const real populationBouncedBack = vf::gpu::getInterpolatedDistributionForNoSlipWithPressureBC(
            subgridDistance, populations[dir], populations[inverseDirection], feq, relaxationFrequency, drho, weight);

        populationsBouncedBack[inverseDirection] = populationBouncedBack;
    });
}

template <size_t direction>
constexpr void addWallMomentum(const bool* linkIsCut, const real* populationsBouncedBack, const real* populations,
                               real3& wallMomentum)
{
}

constexpr real3 computeWallMomentumBounceBack(const bool* linkIsCut, const real* populationsBouncedBack,
                                              const real* populations)
{
    using namespace vf::lbm::dir;
    real3 wallMomentum {};
    forEachNonRestDirection([&](auto dir) {
        using namespace vf::lbm::dir;
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = inverseDir<dir>();
        const real momentum = populations[dir] + populationsBouncedBack[inverseDirection];

        wallMomentum.x += momentum * getComponentX<dir>();
        wallMomentum.y += momentum * getComponentY<dir>();
        wallMomentum.z += momentum * getComponentZ<dir>();
    });

    return wallMomentum;
}

inline __device__ real3 computeFakeWallVelocity(const real3 wallNormal, const real3 velocityForClipping,
                                                const real3 wallShearStress, const real density,
                                                const real interpolationFactor, const real wallArea,
                                                const real3 wallMomentum)
{

    const real3 wallModelForce = wallShearStress * wallArea;
    const real3 wallParallelMomentum = wallMomentum - wallNormal * dot(wallMomentum, wallNormal);

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

constexpr real3 writeDistributionsBB(const Distributions27& populationReferences, const bool* linkIsCut,
                                     const real* populationsBouncedBack, const real3 velocity, const real density,
                                     const vf::gpu::ListIndices& listIndices)
{
    using namespace vf::lbm::dir;

    real3 wallMomentumAdded {};

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = inverseDir<dir>();
        const real addedMomentum = -c6o1 * density * getWeight<dir>() * getVelocity<dir>(velocity.x, velocity.y, velocity.z);
        vf::gpu::writeInInverseDirection<dir>(populationsBouncedBack[inverseDirection] + addedMomentum, listIndices,
                                              populationReferences);

        wallMomentumAdded.x += addedMomentum * getComponentX<dir>();
        wallMomentumAdded.y += addedMomentum * getComponentY<dir>();
        wallMomentumAdded.z += addedMomentum * getComponentZ<dir>();
    });
    return wallMomentumAdded;
}

constexpr real3 writeDistributionsInterpolatedBB(const Distributions27& populationReferences, const bool* linkIsCut,
                                                 const real* populationsBouncedBack, const real3 velocity,
                                                 const real density, const SubgridDistances27& subgridDistances,
                                                 const vf::gpu::ListIndices& listIndices, const uint nodeIndex)
{
    using namespace vf::lbm::dir;

    real3 wallMomentumAdded {};

    forEachNonRestDirection([&](auto dir) {
        if (!linkIsCut[dir])
            return;
        const size_t inverseDirection = inverseDir<dir>();
        const real microVelocity = getVelocity<dir>(velocity.x, velocity.y, velocity.z);
        const real addedMomentum = -c6o1 * density * getWeight<dir>() * microVelocity / subgridDistances.q[dir][nodeIndex];
        vf::gpu::writeInInverseDirection<dir>(populationsBouncedBack[inverseDirection] + addedMomentum, listIndices,
                                              populationReferences);

        wallMomentumAdded.x += addedMomentum * getComponentX<dir>();
        wallMomentumAdded.y += addedMomentum * getComponentY<dir>();
        wallMomentumAdded.z += addedMomentum * getComponentZ<dir>();
    });
    return wallMomentumAdded;
}

#endif
