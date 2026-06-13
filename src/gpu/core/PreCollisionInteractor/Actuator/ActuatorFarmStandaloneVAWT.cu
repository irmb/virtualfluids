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
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \ingroup gpu_core core
//! \{
//! \author Nils Horneff
//======================================================================================

#include "ActuatorFarmInlines.h"
#include "ActuatorFarmStandaloneVAWT.h"
#include "Cuda/CudaMemoryManager.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <basics/constants/NumericConstants.h>
#include <logger/Logger.h>

using namespace vf::basics::constant;

namespace {
real interpMonotone(real x, const std::vector<real>& xs, const std::vector<real>& ys)
{
    if (xs.empty() || ys.empty())
        return c0o1;
    if (xs.size() != ys.size())
        throw std::runtime_error("ActuatorFarmStandaloneVAWT::interpMonotone: xs and ys need same size.");
    if (xs.size() == 1)
        return ys[0];

    x = std::clamp(x, xs.front(), xs.back());
    const auto it = std::upper_bound(xs.begin(), xs.end(), x);
    const size_t i1 = static_cast<size_t>(it - xs.begin());
    if (i1 == 0)
        return ys[0];
    if (i1 >= xs.size())
        return ys.back();

    const size_t i0 = i1 - 1;
    const real x0 = xs[i0];
    const real x1 = xs[i1];
    const real denom = std::max(x1 - x0, std::numeric_limits<real>::epsilon());
    const real t = (x - x0) / denom;
    return (c1o1 - t) * ys[i0] + t * ys[i1];
}
} // namespace

using namespace vf::gpu;

void ActuatorFarmStandaloneVAWT::init()
{
    ActuatorFarm::init();
    this->updateCoordinatesVAWT(c0o1, c0o1);
    this->cudaMemoryManager->cudaCopyCoordsHtoD(this);
}

std::vector<real> ActuatorFarmStandaloneVAWT::computeBladeRadii(const real diameter, const uint numberOfNodesPerBlade)
{
    return std::vector<real>(numberOfNodesPerBlade, c1o2 * diameter);
}

std::vector<real> ActuatorFarmStandaloneVAWT::computeBladeHeights(const real rotorHeight, const uint numberOfNodesPerBlade)
{
    std::vector<real> bladeHeights(numberOfNodesPerBlade, c0o1);
    if (numberOfNodesPerBlade == 0)
        return bladeHeights;

    const real deltaZ = rotorHeight / static_cast<real>(numberOfNodesPerBlade);
    for (uint i = 0; i < numberOfNodesPerBlade; ++i)
        bladeHeights[i] = (static_cast<real>(i) + c1o2) * deltaZ;
    return bladeHeights;
}

std::vector<real> ActuatorFarmStandaloneVAWT::solveGauss(std::vector<std::vector<real>> D, std::vector<real> b)
{
    const uint n = static_cast<uint>(D.size());
    constexpr real eps = real(1e-12);

    for (uint k = 0; k < n; ++k) {
        uint pivotRow = k;
        real maxAbs = std::fabs(D[k][k]);
        for (uint row = k + 1; row < n; ++row) {
            const real value = std::fabs(D[row][k]);
            if (value > maxAbs) {
                maxAbs = value;
                pivotRow = row;
            }
        }

        if (pivotRow != k) {
            std::swap(D[pivotRow], D[k]);
            std::swap(b[pivotRow], b[k]);
        }

        real pivot = D[k][k];
        if (std::fabs(pivot) < eps) {
            pivot = pivot >= c0o1 ? eps : -eps;
            D[k][k] = pivot;
        }

        for (uint row = k + 1; row < n; ++row) {
            const real lowerValue = D[row][k] / pivot;
            D[row][k] = lowerValue;
            for (uint col = k + 1; col < n; ++col)
                D[row][col] -= lowerValue * D[k][col];
            b[row] -= lowerValue * b[k];
        }
    }

    for (int row = static_cast<int>(n) - 1; row >= 0; --row) {
        for (uint col = static_cast<uint>(row) + 1; col < n; ++col)
            b[static_cast<uint>(row)] -= D[static_cast<uint>(row)][col] * b[col];

        real diagonal = D[static_cast<uint>(row)][static_cast<uint>(row)];
        if (std::fabs(diagonal) < eps)
            diagonal = diagonal >= c0o1 ? eps : -eps;
        b[static_cast<uint>(row)] /= diagonal;
    }

    return b;
}

std::vector<real> ActuatorFarmStandaloneVAWT::computeEndEffectsDistribution(
    const real rotorHeight, const real bladeChord, const uint numberOfNodesPerBlade)
{
    const uint n = numberOfNodesPerBlade;
    std::vector<real> fEnd(n, c0o1);
    if (n == 0)
        return fEnd;

    std::vector<real> c(n, bladeChord);
    std::vector<real> angleOfAttackRad(n, real(0.1));
    std::vector<real> theta(n, c0o1);
    std::vector<real> relVelMag(n, c1o1);
    for (uint m = 0; m < n; ++m)
        theta[m] = ((static_cast<real>(m) + c1o2) / static_cast<real>(n)) * cPi;

    std::vector<std::vector<real>> D(n, std::vector<real>(n, c0o1));
    for (uint i = 0; i < n; ++i) {
        const real mode = static_cast<real>(i + 1);
        for (uint m = 0; m < n; ++m) {
            D[m][i] = ((c2o1 * rotorHeight) / (cPi * c[m])) * std::sin(mode * theta[m]) +
                      mode * std::sin(mode * theta[m]) / std::sin(theta[m]);
        }
    }

    const std::vector<real> coefficients = solveGauss(std::move(D), std::move(angleOfAttackRad));

    std::vector<real> circulation(n, c0o1);
    std::vector<real> clDistribution(n, c0o1);
    for (uint m = 0; m < n; ++m) {
        real sumA = c0o1;
        for (uint i = 0; i < n; ++i) {
            const real mode = static_cast<real>(i + 1);
            sumA += coefficients[i] * std::sin(mode * theta[m]);
        }
        circulation[m] = c2o1 * rotorHeight * relVelMag[m] * sumA;
        clDistribution[m] = circulation[m] / (c1o2 * c[m] * relVelMag[m]);
    }

    const real clMax = *std::max_element(clDistribution.begin(), clDistribution.end());
    for (uint m = 0; m < n; ++m)
        fEnd[m] = clDistribution[m] / clMax;

    return fEnd;
}

real ActuatorFarmStandaloneVAWT::get_flowCurvature(const real vrel, const real xchord, const real rotorSpeed,
                                                   const real bladeChord)
{
    return (xchord + c1o4) * bladeChord * rotorSpeed / vrel;
}

void ActuatorFarmStandaloneVAWT::updateForcesAndCoordinates(real time, real deltaT)
{
    this->updateCoordinatesVAWT(time, deltaT);
    this->updateForcesVAWT(time, deltaT);
}

void ActuatorFarmStandaloneVAWT::updateCoordinatesVAWT([[maybe_unused]] real time, const real deltaT)
{
    const uint bladesPerTurbine = this->getNumberOfBladesPerTurbine();
    const real deltaAzimuth = c2Pi / static_cast<real>(bladesPerTurbine);
    real radius = c0o1;

    for (uint turbine = 0; turbine < this->numberOfTurbines; ++turbine) {
        const real azimuthNewRaw = this->azimuths[turbine] + deltaT * this->rotorSpeeds[turbine];
        real azimuthNew = std::fmod(azimuthNewRaw, c2Pi);
        if (azimuthNew < c0o1)
            azimuthNew += c2Pi;
        this->azimuths[turbine] = azimuthNew;

        for (uint blade = 0; blade < bladesPerTurbine; ++blade) {
            const real localAzimuthRad = azimuthNew + static_cast<real>(blade) * deltaAzimuth;
            const real cosTheta = std::cos(localAzimuthRad);
            const real sinTheta = std::sin(localAzimuthRad);

            for (uint bladeNode = 0; bladeNode < this->numberOfPointsPerBlade; ++bladeNode) {
                const uint node =
                    calcPointIndexInBladeArrays({ turbine, blade, bladeNode }, this->numberOfPointsPerBlade, bladesPerTurbine);

                radius = this->bladeRadii[bladeNode];
                const real localX = -sinTheta * radius;
                const real localY = cosTheta * radius;
                const real localZ = -this->rotorHeight * c1o2 + this->bladeHeights[bladeNode];

                this->getAllBladeCoordsX()[node] = localX + this->turbinePosXH[turbine];
                this->getAllBladeCoordsY()[node] = localY + this->turbinePosYH[turbine];
                this->getAllBladeCoordsZ()[node] = localZ + this->turbinePosZH[turbine];
            }
        }
    }
}

void ActuatorFarmStandaloneVAWT::updateForcesVAWT([[maybe_unused]] real time, [[maybe_unused]] real deltaT)
{
    if (this->numberOfPointsPerBlade == 0)
        return;

    const uint bladesPerTurbine = this->getNumberOfBladesPerTurbine();
    const real deltaAzimuth = c2Pi / static_cast<real>(bladesPerTurbine);
    const real deltaZ = this->rotorHeight / static_cast<real>(this->numberOfPointsPerBlade);
    real radius = c0o1;
    const real dynamicPressureReference = c1o2 * this->para->getDensityRatio() * this->velocityInlet * this->velocityInlet;

    for (uint turbine = 0; turbine < this->numberOfTurbines; ++turbine) {
        const real rotorSpeed = this->rotorSpeeds[turbine];
        const real azimuthRad = this->azimuths[turbine];

        for (uint blade = 0; blade < bladesPerTurbine; ++blade) {
            const real localAzimuthRad = azimuthRad + static_cast<real>(blade) * deltaAzimuth;
            const real cosTheta = std::cos(localAzimuthRad);
            const real sinTheta = std::sin(localAzimuthRad);

            for (uint bladeNode = 0; bladeNode < this->numberOfPointsPerBlade; ++bladeNode) {
                const uint node =
                    calcPointIndexInBladeArrays({ turbine, blade, bladeNode }, this->numberOfPointsPerBlade, bladesPerTurbine);
                radius = this->bladeRadii[bladeNode];

                const real velocityX = this->getAllBladeVelocitiesX()[node];
                const real velocityY = this->getAllBladeVelocitiesY()[node];

                const real velocityNormal = sinTheta * velocityX - cosTheta * velocityY;
                const real velocityTangential = cosTheta * velocityX + sinTheta * velocityY + rotorSpeed * radius;
                const real flowAngleRad = std::atan2(velocityNormal, velocityTangential);
                const real vrelSq = velocityNormal * velocityNormal + velocityTangential * velocityTangential;
                const real vrel = std::max(std::sqrt(vrelSq), real(1e-6));

                real flowCurvatureCorrection = c0o1;
                if (this->flagFlowCurvature)
                    flowCurvatureCorrection =
                        get_flowCurvature(vrel, this->bladeMountingPoint, rotorSpeed, this->bladeChord);

                const real angleOfAttackRad = flowAngleRad + flowCurvatureCorrection - this->bladePitch*cPio180;
                real polarLiftCoefficient = interpMonotone(angleOfAttackRad * c180oPi, this->polarAngleOfAttackDeg, this->polarLiftCoefficient);
                const real polarDragCoefficient = interpMonotone(angleOfAttackRad * c180oPi, this->polarAngleOfAttackDeg, this->polarDragCoefficient);

                if (this->flagEndEffects)
                    polarLiftCoefficient *= this->endEffectsDistribution[bladeNode];

                const real dynamicPressure = c1o2 * this->para->getDensityRatio() * vrelSq;
                const real lift = dynamicPressure * this->bladeChord * deltaZ * polarLiftCoefficient;
                const real drag = dynamicPressure * this->bladeChord * deltaZ * polarDragCoefficient;
                const real pn = lift * std::cos(flowAngleRad) + drag * std::sin(flowAngleRad);
                const real pt = lift * std::sin(flowAngleRad) - drag * std::cos(flowAngleRad);

                this->getAllBladeForcesX()[node] = -(pn * sinTheta - pt * cosTheta);
                this->getAllBladeForcesY()[node] = -(-pn * cosTheta - pt * sinTheta);
                this->getAllBladeForcesZ()[node] = c0o1;
                if (this->useLocalSmearingWidth()) {
                    const real epsilonChord = c1o4 * this->bladeChord;
                    const real epsilonDrag  = c1o2 * polarDragCoefficient * this->bladeChord;
                    const real epsilonMesh  = c4o1 * this->para->getScaledLengthRatio(this->level);
                    this->getAllBladeLocalSmearingWidth()[node] = std::max(epsilonChord, std::max(epsilonDrag, epsilonMesh));
                }

                const real normalizedForceScale = std::max(dynamicPressureReference * radius * deltaZ, real(1e-12));
                this->forceNormal[node] = pn / normalizedForceScale;
                this->forceTangential[node] = pt / normalizedForceScale;
                this->angleOfAttackDeg[node] = angleOfAttackRad*c180oPi;

                real azimuthWrapped = std::fmod(localAzimuthRad, c2Pi);
                if (azimuthWrapped < c0o1)
                    azimuthWrapped += c2Pi;
                this->azimuthDeg[node] = azimuthWrapped*c180oPi;
            }
        }
    }
}

void ActuatorFarmStandaloneVAWT::appendOutputData(std::vector<std::string>& dataNames,
                                                  std::vector<std::vector<double>>& nodeData) const
{
    const uint totalPoints = this->getTotalNumberOfPoints();

    const auto appendField = [&](const std::string& name, const std::vector<real>& source) {
        dataNames.push_back(name);
        std::vector<double> values(totalPoints, 0.0);
        for (size_t i = 0; i < static_cast<size_t>(this->numberOfBladePoints); ++i)
            values[i] = source[i];
    
        nodeData.push_back(std::move(values));
    };

    appendField("forceNormal", this->forceNormal);
    appendField("forceTangential", this->forceTangential);
    appendField("angleOfAttackDeg", this->angleOfAttackDeg);
    appendField("azimuthDeg", this->azimuthDeg);
}

//! \}
