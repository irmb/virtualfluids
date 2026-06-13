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
//! \author Henrik Asmuth, Henry Korb, Anna Wellmann, Nils Horneff
//======================================================================================

#ifndef ACTUATOR_FARM_INLINES
#define ACTUATOR_FARM_INLINES

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include "Axis.h"

#include "Utilities/GeometryUtils.h"

namespace vf::gpu {


struct TurbinePointIndex {
    uint turbine;
    uint blade;
    uint bladePoint;
};

constexpr uint calcPointIndexInBladeArrays(uint bladePoint, uint numberOfPointsPerBlade, uint blade, uint numberOfBlades, uint turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    return bladePoint + numberOfPointsPerBlade * (blade + numberOfBlades * turbine);
}

constexpr uint calcPointIndexInBladeArrays(const TurbinePointIndex &turbineNodeIndex, uint numberOfPointsPerBlade, uint numberOfBlades)
{
    return calcPointIndexInBladeArrays(turbineNodeIndex.bladePoint, numberOfPointsPerBlade, turbineNodeIndex.blade, numberOfBlades, turbineNodeIndex.turbine);
}

constexpr void calcTurbineBladeAndBladePoint(uint node, uint &bladePoint, uint numberOfPointsPerBlade, uint &blade, uint numberOfBlades, uint &turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    turbine = node / (numberOfPointsPerBlade * numberOfBlades);
    const uint x_off = turbine * numberOfPointsPerBlade * numberOfBlades;
    blade = (node - x_off) / numberOfPointsPerBlade;
    const uint y_off = numberOfPointsPerBlade * blade + x_off;
    bladePoint = node - y_off;
}

constexpr TurbinePointIndex calcTurbineBladeAndBladePoint(uint node, uint numberOfPointsPerBlade, uint numberOfBlades)
{
    uint turbine = 0;
    uint blade = 0;
    uint bladePoint = 0;
    calcTurbineBladeAndBladePoint(node, bladePoint, numberOfPointsPerBlade, blade, numberOfBlades, turbine);
    return { /*.turbine = */ turbine, /*.blade = */ blade, /*.bladePoint = */ bladePoint }; // Designated initializers are a C++ 20 feature
}

__host__ __device__ __inline__ void rotateFromBladeToGlobal(real bladeCoordX_BF, real bladeCoordY_BF, real bladeCoordZ_BF,
                                                            real& bladeCoordX_GF, real& bladeCoordY_GF, real& bladeCoordZ_GF,
                                                            real azimuth)
{
    rotateAboutX3D(azimuth, bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF);
}

__host__ __device__ __inline__ void rotateFromGlobalToBlade(real& bladeCoordX_BF, real& bladeCoordY_BF, real& bladeCoordZ_BF,
                                                            real bladeCoordX_GF, real bladeCoordY_GF, real bladeCoordZ_GF,
                                                            real azimuth)
{
    invRotateAboutX3D(azimuth, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF, bladeCoordX_BF, bladeCoordY_BF,
                      bladeCoordZ_BF);
}
constexpr real distSqrd(real distX, real distY, real distZ)
{
    return distX * distX + distY * distY + distZ * distZ;
}
constexpr real getRotorBoundingVolumeRadius(real diameter, real smearingWidth, bool useVAWTVolume)
{
    using namespace vf::basics::constant;
    return c1o2 * diameter + c3o2 * smearingWidth;
}

constexpr bool inSphereVolume(real distX, real distY, real distZ, real diameter, real smearingWidth)
{
    const real boundingSphereRadius = getRotorBoundingVolumeRadius(diameter, smearingWidth, false);
    return distSqrd(distX, distY, distZ) < boundingSphereRadius * boundingSphereRadius;
}

template <Axis CylinderAxis, bool UseDonut>
__host__ __device__ __inline__ bool inCylinderVolume(real gridX, real gridY, real gridZ,
                                                     real centerX, real centerY, real centerZ,
                                                     real cylinderLength, real cylinderRadius,
                                                     real margin)
{
    using namespace vf::basics::constant;

    real axialCoord, axialCenter;
    real radial1, radial2, radialCenter1, radialCenter2;

    if constexpr (CylinderAxis == x) {
        axialCoord = gridX;
        axialCenter = centerX;
        radial1 = gridY;
        radial2 = gridZ;
        radialCenter1 = centerY;
        radialCenter2 = centerZ;
    } else if constexpr (CylinderAxis == y) {
        axialCoord = gridY;
        axialCenter = centerY;
        radial1 = gridX;
        radial2 = gridZ;
        radialCenter1 = centerX;
        radialCenter2 = centerZ;
    } else { // z
        axialCoord = gridZ;
        axialCenter = centerZ;
        radial1 = gridX;
        radial2 = gridY;
        radialCenter1 = centerX;
        radialCenter2 = centerY;
    }

    const real axialMin = axialCenter - cylinderLength * c1o2 - margin;
    const real axialMax = axialCenter + cylinderLength * c1o2 + margin;

    if (axialCoord < axialMin || axialCoord > axialMax)
        return false;

    const real d1 = radial1 - radialCenter1;
    const real d2 = radial2 - radialCenter2;
    const real squaredRadialDist = d1 * d1 + d2 * d2;

    const real maxRadialDist = cylinderRadius + margin;

    real minRadialDist;
    if constexpr (UseDonut)
        minRadialDist = std::max(cylinderRadius - margin, real(0));
    else
        minRadialDist = c0o1;

    return squaredRadialDist <= maxRadialDist * maxRadialDist &&
           squaredRadialDist >= minRadialDist * minRadialDist;
}

__host__ __device__ __inline__ real gaussianSmearing(real distX, real distY, real distZ, real smearingWidth)
{
    using namespace vf::basics::constant;

    return std::pow(smearingWidth * std::sqrt(cPi), -3) * std::exp(-distSqrd(distX, distY, distZ) / (smearingWidth * smearingWidth)) ;
}


constexpr bool isHubPoint(uint pointIndex, uint numberOfBladePoints, uint totalPoints)
{
    return pointIndex >= numberOfBladePoints && pointIndex < totalPoints;
}

constexpr uint getHubPointIndex(uint pointIndex, uint numberOfBladePoints)
{
    return pointIndex - numberOfBladePoints;
}



constexpr bool isTowerPoint(uint pointIndex, uint numberOfBladePoints, 
                                                uint numberOfHubPoints, uint totalPoints)
{
    const uint towerStartIndex = numberOfBladePoints + numberOfHubPoints;
    return pointIndex >= towerStartIndex && pointIndex < totalPoints;
}

constexpr uint getTowerPointIndex(uint pointIndex, uint numberOfBladePoints, 
                                                     uint numberOfHubPoints)
{
    return pointIndex - numberOfBladePoints - numberOfHubPoints;
}

}

#endif

//! \}
