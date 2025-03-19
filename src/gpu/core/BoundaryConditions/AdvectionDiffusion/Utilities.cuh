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

constexpr real computeDistributionInterpolated(real subgridDistance, real relaxationFrequency, real distribution,
                                               real equilibrium, real inverseDistribution, real equilibriumInverseWall)
{
    using namespace vf::basics::constant;

    const real wallContribution = c2o1 * equilibriumInverseWall;
    const real interpolated = (distribution * (subgridDistance * relaxationFrequency - c1o1) -
                               equilibrium * relaxationFrequency * (subgridDistance - c1o1)) /
                                  (relaxationFrequency - c1o1) +
                              inverseDistribution * subgridDistance;
    return (wallContribution - interpolated) / (subgridDistance + c1o1);
}

template <size_t direction>
constexpr void writeDistributionInterpolated(const SubgridDistances27& subgridDistances, const real* distributions,
                                             const Distributions27& distributionsGlobal, uint nodeIndex,
                                             const vf::gpu::ListIndices& listIndices, real concentrationNode,
                                             real concentrationWall, real vx1, real vx2, real vx3, real velocityWallX,
                                             real velocityWallY, real velocityWallZ, real relaxationFrequency)
{
    using namespace vf::basics::constant;

    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance > c1o1 || subgridDistance < c0o1)
        return;

    const size_t inverseDirection = vf::lbm::dir::inverseDir<direction>();
    const real distributionPostColl = distributions[direction];
    const real distributionInverse = distributions[inverseDirection];
    const uint writeNode = listIndices.getIndex<inverseDirection>();
    const real equilibriumNode =
        vf::lbm::advection_diffusion::computeEquilibrium<direction>(concentrationNode, vx1, vx2, vx3);
    const real equilibriumInverseWall = vf::lbm::advection_diffusion::computeEquilibrium<inverseDirection>(
        concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    const real dist = computeDistributionInterpolated(subgridDistance, relaxationFrequency, distributionPostColl,
                                                      equilibriumNode, distributionInverse, equilibriumInverseWall);

    (distributionsGlobal.f[inverseDirection])[writeNode] = dist;
}

template <size_t direction>
constexpr void writeDistributionSimpleAntiBounceBack(const SubgridDistances27& subgridDistances, const real* distributions,
                                                     const Distributions27& distributionsGlobal, uint nodeIndex,
                                                     const vf::gpu::ListIndices& listIndices, real concentrationWall,
                                                     real velocityWallX, real velocityWallY, real velocityWallZ)
{
    using namespace vf::basics::constant;

    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance > c1o1 || subgridDistance < c0o1)
        return;

    const real distributionBouncedBack = -distributions[direction];
    const real equilibriumWall = vf::lbm::advection_diffusion::computeEquilibrium<direction>(concentrationWall, velocityWallX,
                                                                                            velocityWallY, velocityWallZ);
    const size_t inverseDirection = vf::lbm::dir::inverseDir<direction>();
    const uint writeNode = listIndices.getIndex<inverseDirection>();

    (distributionsGlobal.f[inverseDirection])[writeNode] = distributionBouncedBack + c2o1 * equilibriumWall;
}

constexpr void writeTemperatureDistributionsInterpolatedAntiBounceBack(
    uint nodeIndex, const SubgridDistances27& subgridDistances, const DistributionReferences27& distributionReferences,
    const vf::gpu::ListIndices& listIndices, const real* distributions, real relaxationFrequency, real vx1, real vx2,
    real vx3, real velocityWallX, real velocityWallY, real velocityWallZ, real concentrationNode, real concentrationWall)
{
    using namespace vf::lbm::dir;

    writeDistributionInterpolated<dP00>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dM00>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0P0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0M0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d00P>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d00M>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPP0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMM0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPM0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMP0>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dP0P>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dM0M>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dP0M>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dM0P>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0PP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0MM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0PM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<d0MP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPPP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMMM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPPM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMMP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPMP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMPM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dPMM>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
    writeDistributionInterpolated<dMPP>(subgridDistances, distributions, distributionReferences, nodeIndex, listIndices,
                                        concentrationNode, concentrationWall, vx1, vx2, vx3, velocityWallX, velocityWallY,
                                        velocityWallZ, relaxationFrequency);
}

constexpr void writeTemperatureDistributionsSimpleAntiBounceBack(uint nodeIndex, const SubgridDistances27& subgridDistances,
                                                                 const DistributionReferences27& distributionReferences,
                                                                 const vf::gpu::ListIndices& listIndices,
                                                                 const real* distributions, real velocityWallX,
                                                                 real velocityWallY, real velocityWallZ,
                                                                 real concentrationWall)
{
    using namespace vf::lbm::dir;

    writeDistributionSimpleAntiBounceBack<dP00>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dM00>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0P0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0M0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d00P>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d00M>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPP0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMM0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPM0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMP0>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dP0P>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dM0M>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dP0M>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dM0P>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0PP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0MM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0PM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<d0MP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPPP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMPP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPMP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMMP>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPPM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMPM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dPMM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    writeDistributionSimpleAntiBounceBack<dMMM>(subgridDistances, distributions, distributionReferences, nodeIndex,
                                                listIndices, concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
}
#endif