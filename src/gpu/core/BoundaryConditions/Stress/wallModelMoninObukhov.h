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

#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"

static constexpr real stabilityFactorMomentum = 4.8F;
static constexpr real stabilityFactorTemperature = 7.8F;

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

struct StabilityCorrections
{
    real momentum, temperature;
};

inline __device__ real computeFrictionVelocity(const real velocity, const real vonKarmanConstant, const real coordZ,
                                               const real roughnessLength, const real stabilityCorrection)
{
    return velocity * vonKarmanConstant / (std::log(coordZ / roughnessLength) - stabilityCorrection);
}

constexpr real3 computeWallShearStress(const real frictionVelocity, const real3 velocityWallParallel,
                                       const real velocityWallParallelMeanMagnitude, const real density)
{
    // Scale wall shear stress with near wall velocity, i.e., Schumann-Grötzbach
    // (SG) approach
    const real wallShearStress = frictionVelocity * frictionVelocity * density;
    return velocityWallParallel * (wallShearStress / velocityWallParallelMeanMagnitude);
}

constexpr real computeStabilityParameter(const real height, const real gravity, const real temperatureScale,
                                         const real frictionVelocity, const real referenceTemperature,
                                         real vonKarmanConstant)
{
    using namespace vf::basics::constant;
    if (temperatureScale == c0o1)
        return c0o1;
    const real limit = c1o1;
    const real numerator = vonKarmanConstant * gravity * temperatureScale * height;
    const real denominator = frictionVelocity * frictionVelocity * referenceTemperature;
    return std::clamp(numerator / denominator, -limit, limit);
}

// Compute stability parameters according to <a
// href="https://shop.elsevier.com/books/introduction-to-micrometeorology/holton/978-0-12-059354-5"><b>[p. 223, Arya 2001,
//! ISBN:9780120593545 ]</b></a>
inline __device__ __host__ StabilityCorrections computeStabilityCorrectionsUnstable(real stabilityParameter)
{
    using namespace vf::basics::constant;
    const real tmp = std::sqrt(std::sqrt(c1o1 - c15o1 * stabilityParameter));
    const real momentum =
        std::log(c1o2 * (c1o1 + tmp * tmp) * c1o2 * (c1o1 + tmp) * c1o2 * (c1o1 + tmp)) - c2o1 * std::atan(tmp) + cPi * c1o2;
    const real temperature = c2o1 * std::log(c1o2 * (c1o1 + tmp * tmp));
    return { momentum, temperature };
}
constexpr StabilityCorrections computeStabilityCorrectionsStable(real stabilityParameter)
{
    return { -stabilityFactorMomentum * stabilityParameter, -stabilityFactorTemperature * stabilityParameter }; // Beare 2006
}

// follows
inline __device__ __host__ real computeStabilityParameterFromSurfaceTemperature(
    const real coordZ, const real roughnessLength, const real roughnessLengthTemperature, const real velocity,
    const real temperatureDifference, const real referenceTemperature, const real gravity)
{
    const real logDistanceMomentum = std::log(coordZ / roughnessLength);
    const real logDistanceTemperature = std::log(coordZ / roughnessLengthTemperature);
    const real a1 = stabilityFactorTemperature / logDistanceTemperature;
    const real a2 = stabilityFactorMomentum / logDistanceMomentum;
    const real a3 = gravity * coordZ / (velocity * velocity * referenceTemperature) * temperatureDifference *
                    logDistanceMomentum * logDistanceMomentum / logDistanceTemperature;
    const real enumerator = (c2o1 * a2 * a3 - c1o1 + std::sqrt(c1o1 + c4o1 * a3 * (a1 - a2)));
    const real denominator = (c2o1 * (a1 - a2 * a2 * a3));
    return enumerator / denominator;
}

inline __device__ __host__ real computeStabilityParameterFromHeatFlux(const real coordZ, const real roughnessLength,
                                                                      const real velocity, const real surfaceHeatFlux,
                                                                      const real referenceTemperature, const real gravity,
                                                                      const real vonKarmanConstant)
{
    const real logDistance = std::log(coordZ / roughnessLength);
    const real b1 = -gravity * coordZ * surfaceHeatFlux /
                    (vonKarmanConstant * vonKarmanConstant * velocity * velocity * velocity * referenceTemperature) *
                    logDistance * logDistance * logDistance;
    const real b2 = stabilityFactorMomentum / logDistance;
    const real sqrtB1B2 = std::sqrt(c3o1 * b1 * b2 * b2 * b2);
    return c2o1 / sqrtB1B2 * std::cos(c1o3 * std::acos(-c3o2 / b2 * sqrtB1B2 - c2Pi * c1o3)) - c1o1 / b2;
}

#endif // WallModel_MoninObukhov_H_