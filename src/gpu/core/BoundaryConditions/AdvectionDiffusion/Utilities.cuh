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
#ifndef ADVECTION_DIFFUSION_UTILITIES_CUH
#define ADVECTION_DIFFUSION_UTILITIES_CUH

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/advectionDiffusion/Equilibrium.h>
#include <lbm/constants/D3Q27.h>

#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"

constexpr real computePopulationAntiBounceBackInterpolated(const real subgridDistance, const real relaxationFrequency,
                                                           const real population, const real equilibrium,
                                                           const real inversePopulation, const real equilibriumInverseWall)
{
    using namespace vf::basics::constant;

    const real interpolated = (population * (subgridDistance * relaxationFrequency - c1o1) -
                               equilibrium * relaxationFrequency * (subgridDistance - c1o1)) /
                                  (relaxationFrequency - c1o1) +
                              inversePopulation * subgridDistance;
    return (c2o1 * equilibriumInverseWall - interpolated) / (subgridDistance + c1o1);
}

template <size_t direction>
constexpr void writePopulationAntiBounceBackInterpolated(const SubgridDistances27& subgridDistances, const real* populations,
                                                         const Distributions27& populationReferences, const uint nodeIndex,
                                                         const vf::gpu::ListIndices& listIndices,
                                                         const real concentrationNode, const real concentrationWall,
                                                         const real vx1, real vx2, const real vx3, const real velocityWallX,
                                                         const real velocityWallY, const real velocityWallZ,
                                                         const real relaxationFrequency)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::advection_diffusion;
    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance > c1o1 || subgridDistance < c0o1)
        return;

    const size_t inverseDirection = vf::lbm::dir::inverseDir<direction>();
    const real population = populations[direction];
    const real inversePopulation = populations[inverseDirection];
    const uint writeNode = listIndices.getIndex<inverseDirection>();
    const real equilibrium = computeEquilibrium<direction>(concentrationNode, vx1, vx2, vx3);
    const real equilibriumInverseWall =
        computeEquilibrium<inverseDirection>(concentrationWall, velocityWallX, velocityWallY, velocityWallZ);

    (populationReferences.f[inverseDirection])[writeNode] = computePopulationAntiBounceBackInterpolated(
        subgridDistance, relaxationFrequency, population, equilibrium, inversePopulation, equilibriumInverseWall);
}

template <size_t direction>
constexpr void writePopulationSimpleAntiBounceBack(const SubgridDistances27& subgridDistances, const real* populations,
                                                   const Distributions27& populationReferences, const uint nodeIndex,
                                                   const vf::gpu::ListIndices& listIndices, const real concentrationWall,
                                                   const real velocityWallX, const real velocityWallY,
                                                   const real velocityWallZ)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::advection_diffusion;

    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance > c1o1 || subgridDistance < c0o1)
        return;

    const real equilibriumWall =
        computeEquilibrium<direction>(concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    const size_t inverseDirection = vf::lbm::dir::inverseDir<direction>();
    const uint writeNode = listIndices.getIndex<inverseDirection>();

    (populationReferences.f[inverseDirection])[writeNode] = -populations[direction] + c2o1 * equilibriumWall;
}

constexpr void writePopulationsInterpolatedAntiBounceBack(uint nodeIndex, const SubgridDistances27& subgridDistances,
                                                          const Distributions27& populationReferences,
                                                          const vf::gpu::ListIndices& listIndices, const real* populations,
                                                          const real relaxationFrequency, const real vx1, const real vx2,
                                                          const real vx3, const real velocityWallX, const real velocityWallY,
                                                          const real velocityWallZ, const real concentrationNode,
                                                          const real concentrationWall)
{
    using namespace vf::lbm::dir;

    writePopulationAntiBounceBackInterpolated<dP00>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dM00>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0P0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0M0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d00P>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d00M>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPP0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMM0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPM0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMP0>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dP0P>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dM0M>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dP0M>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dM0P>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0PP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0MM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0PM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<d0MP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPPP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMMM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPPM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMMP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPMP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMPM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dPMM>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
    writePopulationAntiBounceBackInterpolated<dMPP>(subgridDistances, populations, populationReferences, nodeIndex,
                                                    listIndices, concentrationNode, concentrationWall, vx1, vx2, vx3,
                                                    velocityWallX, velocityWallY, velocityWallZ, relaxationFrequency);
}

constexpr void writePopulationsSimpleAntiBounceBack(const uint nodeIndex, const SubgridDistances27& subgridDistances,
                                                    const Distributions27& populationReferences,
                                                    const vf::gpu::ListIndices& listIndices, const real* populations,
                                                    const real velocityWallX, const real velocityWallY,
                                                    const real velocityWallZ, const real concentrationWall)
{
    using namespace vf::lbm::dir;

    writePopulationSimpleAntiBounceBack<dP00>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dM00>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0P0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0M0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d00P>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d00M>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPP0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMM0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPM0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMP0>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dP0P>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dM0M>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dP0M>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dM0P>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0PP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0MM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0PM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<d0MP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPPP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMPP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPMP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMMP>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPPM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMPM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dPMM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writePopulationSimpleAntiBounceBack<dMMM>(subgridDistances, populations, populationReferences, nodeIndex, listIndices,
                                              concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
}
#endif