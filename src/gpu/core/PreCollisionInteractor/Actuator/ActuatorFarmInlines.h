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
//! \author Henrik Asmuth, Henry Korb, Anna Wellmann
//======================================================================================

#ifndef ACTUATOR_FARM_INLINES
#define ACTUATOR_FARM_INLINES

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "Utilities/GeometryUtils.h"

using namespace vf::basics::constant;

struct TurbineNodeIndex {
    uint turbine;
    uint blade;
    uint bladeNode;
};

constexpr uint calcNodeIndexInBladeArrays(uint bladeNode, uint numberOfNodesPerBlade, uint blade, uint numberOfBlades, uint turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    return bladeNode + numberOfNodesPerBlade * (blade + numberOfBlades * turbine);
}

constexpr uint calcNodeIndexInBladeArrays(const TurbineNodeIndex &turbineNodeIndex, uint numberOfNodesPerBlade, uint numberOfBlades)
{
    return calcNodeIndexInBladeArrays(turbineNodeIndex.bladeNode, numberOfNodesPerBlade, turbineNodeIndex.blade, numberOfBlades, turbineNodeIndex.turbine);
}

constexpr void calcTurbineBladeAndBladeNode(uint node, uint &bladeNode, uint numberOfNodesPerBlade, uint &blade, uint numberOfBlades, uint &turbine)
{
    // see https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/merge_requests/248 for visualization
    turbine = node / (numberOfNodesPerBlade * numberOfBlades);
    const uint x_off = turbine * numberOfNodesPerBlade * numberOfBlades;
    blade = (node - x_off) / numberOfNodesPerBlade;
    const uint y_off = numberOfNodesPerBlade * blade + x_off;
    bladeNode = node - y_off;
}

constexpr TurbineNodeIndex calcTurbineBladeAndBladeNode(uint node, uint numberOfNodesPerBlade, uint numberOfBlades)
{
    uint turbine = 0;
    uint blade = 0;
    uint bladeNode = 0;
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    return { /*.turbine = */ turbine, /*.blade = */ blade, /*.bladeNode = */ bladeNode }; // Designated initializers are a C++ 20 feature
}

constexpr void rotateFromBladeToGlobal(
                            real bladeCoordX_BF, real bladeCoordY_BF, real bladeCoordZ_BF, 
                            real& bladeCoordX_GF, real& bladeCoordY_GF, real& bladeCoordZ_GF,
                            real azimuth)
{
    rotateAboutX3D(azimuth, bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF);
}

constexpr void rotateFromGlobalToBlade(
                            real& bladeCoordX_BF, real& bladeCoordY_BF, real& bladeCoordZ_BF, 
                            real bladeCoordX_GF, real bladeCoordY_GF, real bladeCoordZ_GF,
                            real azimuth)
{
    invRotateAboutX3D(azimuth, bladeCoordX_GF, bladeCoordY_GF, bladeCoordZ_GF, bladeCoordX_BF, bladeCoordY_BF, bladeCoordZ_BF);
}
constexpr real distSqrd(real distX, real distY, real distZ)
{
    return distX * distX + distY * distY + distZ * distZ;
}

constexpr real getBoundingSphereRadius(real diameter, real smearingWidth)
{
    return c1o2 * diameter + c3o2 * smearingWidth;
}

constexpr bool inBoundingSphere(real distX, real distY, real distZ, real diameter, real smearingWidth)
{
    const real boundingSphereRadius = getBoundingSphereRadius(diameter, smearingWidth);
    return distSqrd(distX, distY, distZ) < boundingSphereRadius * boundingSphereRadius;
}

constexpr real gaussianSmearing(real distX, real distY, real distZ, real epsilon, real factorGaussian)
{
    return factorGaussian * exp(-distSqrd(distX, distY, distZ) / (epsilon * epsilon));
}

__inline__ void swapArrays(real* &arr1, real* &arr2)
{
    real* tmp = arr1;
    arr1 = arr2;
    arr2 = tmp;
}

#endif

//! \}
