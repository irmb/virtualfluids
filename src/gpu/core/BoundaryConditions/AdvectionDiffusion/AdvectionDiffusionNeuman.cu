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
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb
//=======================================================================================
#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Calculation/Calculation.h"
#include "Utilities.cuh"
#include "Utilities/KernelUtilities.h"
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <cuda_helper/CudaIndexCalculation.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/constants/D3Q27.h>

template <BoundaryConditionFactory::AdvectionDiffusionNeumannBC bcType>
__global__ void
AdvectionDiffusionNeumann_Device(real* distributionsConcentration, AdvectionDiffusionNeumannBoundaryConditions bcParameters,
                                 const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                 const real* velocityX, const real* velocityY, const real* velocityZ,
                                 unsigned long long numberOfLBnodes, real relaxationFrequency, bool isEvenTimestep)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    SubgridDistances27 subgridDistances;
    vf::gpu::getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint k_000 = bcParameters.BCNodeIndices[nodeIndex];
    const vf::gpu::ListIndices listIndices(k_000, neighborX, neighborY, neighborZ);

    Distributions27 distributionReferences =
        vf::gpu::getDistributionReferences27(distributionsConcentration, numberOfLBnodes, isEvenTimestep);

    real distributions[27];
    vf::gpu::getPostCollisionDistribution(distributions, distributionReferences, listIndices);
    const real concentrationNode = vf::lbm::getDensity(distributions);
    const real gradient = bcParameters.gradients[nodeIndex];
    const real concentrationWall = concentrationNode + c1o2 * gradient;

    vf::gpu::getPointersToDistributions(distributionReferences, distributionsConcentration, numberOfLBnodes,
                                        !isEvenTimestep);
    switch (bcType) {
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip:
            writeTemperatureDistributionsSimpleAntiBounceBack(nodeIndex, subgridDistances, distributionReferences,
                                                              listIndices, distributions, velocityX[k_000], velocityY[k_000],
                                                              velocityZ[k_000], concentrationWall);
            break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip:
            writeTemperatureDistributionsSimpleAntiBounceBack(
                nodeIndex, subgridDistances, distributionReferences, listIndices, distributions, bcParameters.vx[nodeIndex],
                bcParameters.vy[nodeIndex], bcParameters.vz[nodeIndex], concentrationWall);
            break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedSlip: {
            const real vx1 = velocityX[k_000];
            const real vx2 = velocityY[k_000];
            const real vx3 = velocityZ[k_000];
            writeTemperatureDistributionsInterpolatedAntiBounceBack(
                nodeIndex, subgridDistances, distributionReferences, listIndices, distributions, relaxationFrequency, vx1,
                vx2, vx3, vx1, vx2, vx3, concentrationNode, concentrationWall);
        } break;
        case BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedNoSlip: {
            writeTemperatureDistributionsInterpolatedAntiBounceBack(
                nodeIndex, subgridDistances, distributionReferences, listIndices, distributions, relaxationFrequency,
                velocityX[k_000], velocityY[k_000], velocityZ[k_000], bcParameters.vx[nodeIndex],
                bcParameters.vy[nodeIndex], bcParameters.vz[nodeIndex], concentrationNode, concentrationWall);
        } break;
    }
}

//! \}
