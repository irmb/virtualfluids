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
//! \author Martin Schoenherr, Henry Korb
//=======================================================================================
#include <basics/constants/NumericConstants.h>

#include <lbm/constants/D3Q27.h>

#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

__global__ void AdvectionDiffusionNoFluxBounceBack_Device(real* distributions,
                                                    AdvectionDiffusionNoFluxBoundaryConditions bcParameters,
                                                    const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                    unsigned long long numberOfLBnodes, bool isEvenTimestep)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    using namespace vf::gpu;

    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    Distributions27 populationReferences = getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);
    const uint indexOfBCnode = bcParameters.BCNodeIndices[nodeIndex];
    const ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    SubgridDistances27 subgridDistances;
    getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    real populations[27];
    getPostCollisionDistribution(populations, populationReferences, listIndices);

    getPointersToDistributions(populationReferences, distributions, numberOfLBnodes, !isEvenTimestep);

    forEachNonRestDirection([&](auto direction) {
        const real subgridDistance = (subgridDistances.q[direction])[nodeIndex];
        if (subgridDistance < c0o1 || subgridDistance > c1o1)
            return;
        writeInInverseDirection<direction>(populations[direction], listIndices, populationReferences);
    });
}

//! \}
