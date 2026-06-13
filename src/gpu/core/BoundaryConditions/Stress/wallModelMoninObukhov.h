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
#ifndef WallModel_MoninObukhov_H_
#define WallModel_MoninObukhov_H_

#include <algorithm>
#include <cmath>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/constants/D3Q27.h>

namespace vf::gpu {

static constexpr real stabilityFactorMomentum = 4.8F;
static constexpr real stabilityFactorTemperature = 7.8F;

constexpr real smoothAndSaveMean(real instantaneous, real filterFrequency, real& mean)
{
    const real smoothed = filterFrequency * instantaneous + (vf::basics::constant::c1o1 - filterFrequency) * mean;
    mean = smoothed;
    return smoothed;
}

constexpr real3 computeTangentialVector(real3 quantity, real3 normal)
{
    return quantity - normal * dot(quantity, normal);
}

inline __device__ real computeMagnitude(real3 vector)
{
    return std::sqrt(square(vector));
}

inline __device__ real computeFrictionVelocity(const real velocity, const real vonKarmanConstant, const real samplingDistance,
                                               const real roughnessLength, const real stabilityCorrection)
{
    return velocity * vonKarmanConstant / (std::log(samplingDistance / roughnessLength) - stabilityCorrection);
}

inline __device__ real computeSurfaceHeatFlux(const real temperatureDifference, const real frictionVelocity,
                                              const real vonKarmanConstant, const real samplingDistance, const real roughnessLength,
                                              const real stabilityCorrection)
{
    return -temperatureDifference * frictionVelocity * vonKarmanConstant /
           (std::log(samplingDistance / roughnessLength) - stabilityCorrection);
}

constexpr real3 computeWallShearStress(const real frictionVelocity, const real3 velocityTangential,
                                       const real velocityTangentialMeanMagnitude, const real density)
{
    // Scale wall shear stress with near wall velocity, i.e., Schumann-Grötzbach
    // (SG) approach
    const real wallShearStress = frictionVelocity * frictionVelocity * density;
    return velocityTangential * (wallShearStress / velocityTangentialMeanMagnitude);
}

constexpr real computeStabilityParameter(const real height, const real gravity, const real surfaceHeatFlux,
                                         const real frictionVelocity, const real referenceTemperature,
                                         const real vonKarmanConstant)
{
    using namespace vf::basics::constant;
    if (surfaceHeatFlux == c0o1)
        return c0o1;
    const real limit = c1o1;
    const real numerator = -vonKarmanConstant * gravity * surfaceHeatFlux * height;
    const real denominator = frictionVelocity * frictionVelocity * frictionVelocity * referenceTemperature;
    return std::clamp(numerator / denominator, -limit, limit);
}

// Compute stability parameters according to <a
// href="https://shop.elsevier.com/books/introduction-to-micrometeorology/holton/978-0-12-059354-5"><b>[p. 223, Arya 2001,
//! ISBN:9780120593545 ]</b></a>
inline __device__ __host__ real computeStabilityCorrectionTemperature(const real stabilityParameter)
{
    using namespace vf::basics::constant;
    if (stabilityParameter >= c0o1)
        return -stabilityFactorTemperature * stabilityParameter;
    return c2o1 * std::log(c1o2 * (c1o1 + std::sqrt(c1o1 - c15o1 * stabilityParameter)));
}
inline __device__ real computeStabilityCorrectionMomentum(const real stabilityParameter)
{
    using namespace vf::basics::constant;
    if (stabilityParameter >= c0o1)
        return -stabilityFactorMomentum * stabilityParameter;
    const real tmp = std::sqrt(std::sqrt(c1o1 - c15o1 * stabilityParameter));
    return std::log(c1o2 * (c1o1 + tmp * tmp) * c1o2 * (c1o1 + tmp) * c1o2 * (c1o1 + tmp)) - c2o1 * std::atan(tmp) +
           cPi * c1o2;
}

}

#endif // WallModel_MoninObukhov_H_