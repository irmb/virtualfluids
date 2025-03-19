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
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"

#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

template <size_t direction>
constexpr void setDistributions(uint nodeIndex, const SubgridDistances27& subgridDistanceReferences,
                                Distributions27& distributionReferences, const vf::gpu::ListIndices& listIndices,
                                const real* populations)
{
    const real subgridDistance = (subgridDistanceReferences.q[direction])[nodeIndex];
    if (subgridDistance < c0o1 || subgridDistance > c1o1)
        return;
    const size_t inverseDir = vf::lbm::dir::inverseDir<direction>();
    const uint writeIndex = listIndices.getIndex<inverseDir>();
    (distributionReferences.f[inverseDir])[writeIndex] = populations[direction];
}

__global__ void AdvectionDiffusionBounceBack_Device(real* distributions, AdvectionDiffusionNoSlipBoundaryConditions bcParameters, const uint* neighborX, const uint* neighborY,
                                                    const uint* neighborZ, unsigned long long numberOfLBnodes,
                                                    bool isEvenTimestep)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    Distributions27 distributionReferences =
        vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);
    const uint indexOfBCnode = bcParameters.BCNodeIndices[nodeIndex];
    const vf::gpu::ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    SubgridDistances27 subgridDistanceReferences;
    vf::gpu::getPointersToSubgridDistances(subgridDistanceReferences, bcParameters.q27[0], bcParameters.numberOfBCnodes);
    
    real populations[27];
    vf::gpu::getPostCollisionDistribution(populations, distributionReferences, listIndices);
    
    vf::gpu::getPointersToDistributions(distributionReferences, distributions, numberOfLBnodes, !isEvenTimestep);

    setDistributions<dM00>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dP00>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0M0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0P0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d00M>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d00P>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMM0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPP0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMP0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPM0>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dM0M>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dP0P>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dM0P>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dP0M>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0MM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0PP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0MP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<d0PM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMMM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPPP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMMP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPPM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMPM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPMP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dMPP>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
    setDistributions<dPMM>(nodeIndex, subgridDistanceReferences, distributionReferences, listIndices, populations);
}

//! \}
