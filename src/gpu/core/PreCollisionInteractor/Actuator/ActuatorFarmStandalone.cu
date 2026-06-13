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
//! \author Henrik Asmuth, Henry Korb, Nils Horneff
//======================================================================================
#include "ActuatorFarmInlines.h"
#include "ActuatorFarmStandalone.h"

#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <cuda_helper/CudaGrid.h>
#include <logger/Logger.h>

using namespace vf::basics::constant;
namespace vf::gpu {

std::vector<real> ActuatorFarmStandalone::computeBladeRadii(const real diameter, const uint numberOfPointsPerBlade)
{
    const real dr = c1o2 * diameter / numberOfPointsPerBlade;
    std::vector<real> bladeRadii(numberOfPointsPerBlade);
    for (uint point = 0; point < numberOfPointsPerBlade; point++)
        bladeRadii[point] = dr * (point + c1o2);
    return bladeRadii;
}

void ActuatorFarmStandalone::updateForcesAndCoordinates([[maybe_unused]] real time, real deltaT)
{
    this->updateBladeForcesAndCoordinates(deltaT);
    if (numberOfHubPoints > 0)
        this->updateHubForces();
    if (numberOfTowerPoints > 0)
        this->updateTowerForces();
}

void ActuatorFarmStandalone::updateBladeForcesAndCoordinates(real deltaT)
{
    const real c0 = c20o1 * c1o10;
    const real deltaAzimuth = c2Pi / this->numberOfBlades;
    const real liftCoefficient = c1o1;
    const real dragCoefficient = c0o1;
    const bool overrideNormalCoefficient = this->bladeNormalCoefficients.has_value();

    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {
        const real rotorSpeed = this->rotorSpeeds[turbine];
        const real azimuthOld = this->azimuths[turbine];
        const real azimuthNew = azimuthOld + deltaT * rotorSpeed;
        this->azimuths[turbine] = azimuthNew > c2Pi ? azimuthNew - c2Pi : azimuthNew;

        for (uint blade = 0; blade < this->numberOfBlades; blade++) {
            const real localAzimuthOld = azimuthOld + blade * deltaAzimuth;
            const real localAzimuthNew = azimuthNew + blade * deltaAzimuth;

            real lastPointRadius = c0o1;
            real currentPointradius = c0o1;
            real nextPointradius = this->bladeRadii[0];

            for (uint bladePoint = 0; bladePoint < this->numberOfPointsPerBlade; bladePoint++) {
                const uint point = calcPointIndexInBladeArrays({ turbine, blade, bladePoint }, this->numberOfPointsPerBlade,
                                                               this->numberOfBlades);

                real uRel, vRel, wRel;
                rotateFromGlobalToBlade(uRel, vRel, wRel, this->getAllBladeVelocitiesX()[point],
                                        this->getAllBladeVelocitiesY()[point], this->getAllBladeVelocitiesZ()[point],
                                        localAzimuthOld);

                lastPointRadius = currentPointradius;
                currentPointradius = nextPointradius;
                nextPointradius = bladePoint < this->numberOfPointsPerBlade - 1 ? this->bladeRadii[bladePoint + c1o1]
                                                                                   : this->diameter * c1o2;

                const real dr = c1o2 * (nextPointradius - lastPointRadius);

                vRel += currentPointradius * rotorSpeed;
                const real uRelSq = uRel * uRel + vRel * vRel;
                const real phi = std::atan2(uRel, vRel);

                const real tmp = c4o1 * currentPointradius / this->diameter - c1o1;
                const real chord = c0 * std::sqrt(c1o1 - tmp * tmp);
                const real normalCoefficient = liftCoefficient * std::cos(phi) + dragCoefficient * std::sin(phi);
                const real tangentialCoefficient = liftCoefficient * std::sin(phi) - dragCoefficient * std::cos(phi);
                const real appliedNormalCoefficient = overrideNormalCoefficient
                    ? this->bladeNormalCoefficients.value()[bladePoint]
                    : normalCoefficient;
                const real appliedTangentialCoefficient = overrideNormalCoefficient
                    ? c0o1
                    : tangentialCoefficient;
                const real fx = -c1o2 * uRelSq * chord * para->getDensityRatio() * appliedNormalCoefficient * dr;
                const real fy = -c1o2 * uRelSq * chord * para->getDensityRatio() * appliedTangentialCoefficient * dr;

                rotateFromBladeToGlobal(fx, fy, c0o1, this->getAllBladeForcesX()[point], this->getAllBladeForcesY()[point],
                                        this->getAllBladeForcesZ()[point], localAzimuthNew);
                rotateFromBladeToGlobal(c0o1, c0o1, currentPointradius, this->getAllBladeCoordsX()[point],
                                        this->getAllBladeCoordsY()[point], this->getAllBladeCoordsZ()[point],
                                        localAzimuthNew);
                getAllBladeCoordsX()[point] += this->turbinePosXH[turbine];
                getAllBladeCoordsY()[point] += this->turbinePosYH[turbine];
                getAllBladeCoordsZ()[point] += this->turbinePosZH[turbine];
            }
        }
    }
}

void ActuatorFarmStandalone::updateHubForces()
{
    // Streamlined cylinder physics
    const real segmentLength = hubLength / static_cast<real>(numberOfHubPointsPerTurbine);
    const real surfaceAreaPerSegment = c2o1 * cPi * hubRadius * segmentLength; // Cylinder surface area per segment
    const real frontBackArea = cPi * hubRadius * hubRadius;                    // Cross-sectional area for pressure drag

    for (uint pointIndex = 0; pointIndex < numberOfHubPoints; pointIndex++) {
        // Get sampled velocity components
        const real velX = this->getAllHubVelocitiesX()[pointIndex];
        const real velY = this->getAllHubVelocitiesY()[pointIndex];
        const real velZ = this->getAllHubVelocitiesZ()[pointIndex];

        const real velocityMagnitude = std::hypot(velX, velY, velZ);

        // SKIN FRICTION DRAG (dominant for streamlined cylinder)
        real forceSkinX = c0o1, forceSkinY = c0o1, forceSkinZ = c0o1;
        if (velocityMagnitude > c10eM12) { // Avoid division by zero
            const real skinFrictionMagnitude = this->hubSkinFrictionCoeff.value() * c1o2 * para->getDensityRatio() *
                                               surfaceAreaPerSegment * velocityMagnitude * velocityMagnitude;
            // Skin friction acts opposite to velocity direction
            forceSkinX = -skinFrictionMagnitude * (velX / velocityMagnitude);
            forceSkinY = -skinFrictionMagnitude * (velY / velocityMagnitude);
            forceSkinZ = -skinFrictionMagnitude * (velZ / velocityMagnitude);
        }

        // SMALL PRESSURE DRAG (front/back faces contribution)
        const real pressureDragAreaPerSegment = frontBackArea / static_cast<real>(numberOfHubPointsPerTurbine);
        const real forcePressureX = -c1o2 * para->getDensityRatio() * this->hubDragCoeff.value() *
                                    pressureDragAreaPerSegment * velX * std::abs(velX);

        // Minor crossflow pressure effects
        const real crossflowFactor = c1o10; // Much smaller than main pressure drag
        const real forcePressureY = -crossflowFactor * c1o2 * para->getDensityRatio() * this->hubDragCoeff.value() *
                                    pressureDragAreaPerSegment * velY * std::abs(velY);
        const real forcePressureZ = -crossflowFactor * c1o2 * para->getDensityRatio() * this->hubDragCoeff.value() *
                                    pressureDragAreaPerSegment * velZ * std::abs(velZ);

        // TOTAL FORCE
        this->getAllHubForcesX()[pointIndex] = forceSkinX + forcePressureX;
        this->getAllHubForcesY()[pointIndex] = forceSkinY + forcePressureY;
        this->getAllHubForcesZ()[pointIndex] = forceSkinZ + forcePressureZ;
    }
}

void ActuatorFarmStandalone::updateTowerForces()
{
    for (uint turbine = 0; turbine < this->numberOfTurbines; turbine++) {
        const real segmentHeight = towerHeights[turbine] / static_cast<real>(numberOfTowerPointsPerTurbine);
        const real projectedArea = c2o1 * towerRadius * segmentHeight; // Crossflow projected area

        for (uint point = 0; point < numberOfTowerPointsPerTurbine; point++) {
            const uint pointIndex = turbine * numberOfTowerPointsPerTurbine + point;
            const real velX = this->getAllTowerVelocitiesX()[pointIndex];
            const real velY = this->getAllTowerVelocitiesY()[pointIndex];
            const real velZ = this->getAllTowerVelocitiesZ()[pointIndex];

            // MAIN EFFECT: Horizontal crossflow drag (dominates wake formation)
            const real velHorizontalMag = std::hypot(velX, velY);
            real forceX = c0o1, forceY = c0o1, forceZ = c0o1;

            if (velHorizontalMag > c10eM12) {
                const real dragMagnitude = c1o2 * para->getDensityRatio() * this->towerDragCoeff.value() * projectedArea *
                                           velHorizontalMag * velHorizontalMag;
                forceX = -dragMagnitude * (velX / velHorizontalMag);
                forceY = -dragMagnitude * (velY / velHorizontalMag);
            }

            forceZ = c0o1 * velZ; // Simplification: No vertical drag. Can be implemented if needed.

            // No viscous effects
            this->getAllTowerForcesX()[pointIndex] = forceX;
            this->getAllTowerForcesY()[pointIndex] = forceY;
            this->getAllTowerForcesZ()[pointIndex] = forceZ;
        }
    }
}
}

//! \}
