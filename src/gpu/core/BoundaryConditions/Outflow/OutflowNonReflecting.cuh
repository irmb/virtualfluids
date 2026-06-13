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
//! \author Martin Schoenherr, Anna Wellmann, Henry Korb
//======================================================================================
#ifndef OUTFLOW_NON_REFLECTING_H_
#define OUTFLOW_NON_REFLECTING_H_

#include "Calculation/Calculation.h"
#include "Utilities/KernelUtilities.h"
#include "basics/constants/NumericConstants.h"
#include "cuda_helper/CudaIndexCalculation.h"
#include "lbm/MacroscopicQuantities.h"
#include "lbm/constants/D3Q27.h"

namespace vf::gpu {

template <size_t dir>
constexpr void computeOutflowDistribution(const Distributions27& populationReferences,
                                          const ListIndices& listIndices, const real* populations,
                                          const real* populationsNeighbor, const real cs, const real densityCorrection)
{
    using namespace vf::basics::constant;
    const real weight = vf::lbm::dir::getWeight<dir>();
    const real population = populationsNeighbor[dir] * cs + (c1o1 - cs) * populations[dir] - weight * densityCorrection;
    writeInSameDirection<dir>(population, listIndices, populationReferences);
}

template <bool applyPressureCorrection>
__global__ void OutflowNonReflecting_Device(real* rhoBC, real* distributions, const int* indicesBCnode,
                                            const int* indicesNeighbor, const int numberOfBCnodes, const uint* neighborX,
                                            const uint* neighborY, const uint* neighborZ,
                                            const unsigned long long numberOfLBnodes, const bool isEvenTimestep,
                                            const size_t direction, const real densityCorrectionFactor)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= numberOfBCnodes)
        return;

    const ListIndices listIndices(indicesBCnode[nodeIndex], neighborX, neighborY, neighborZ);
    const ListIndices neighborIndices(indicesNeighbor[nodeIndex], neighborX, neighborY, neighborZ);

    Distributions27 populationReferences;
    getPointersToDistributions(populationReferences, distributions, numberOfLBnodes, isEvenTimestep);
    real populations[NUMBER_Of_DIRECTIONS], populationsNeighbor[NUMBER_Of_DIRECTIONS];

    getPreCollisionDistribution(populationsNeighbor, populationReferences, neighborIndices);
    getPreCollisionDistribution(populations, populationReferences, listIndices);

    const real densityCorrection =
        applyPressureCorrection ? densityCorrectionFactor * vf::lbm::getDensity(populations) : c0o1;
    const real cs = std::sqrt(c1o3);

    getPointersToDistributions(populationReferences, distributions, numberOfLBnodes, !isEvenTimestep);

    switch (direction) {
        case dM00:
            computeOutflowDistribution<dP00>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPM0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPP0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dP0M>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dP0P>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;

        case dP00:
            computeOutflowDistribution<dM00>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMM0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMP0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dM0M>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dM0P>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;

        case d0M0:
            computeOutflowDistribution<d0P0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPP0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMP0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0PP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0PM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;

        case d0P0:
            computeOutflowDistribution<d0M0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPM0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMM0>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0MP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0MM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;

        case d00M:
            computeOutflowDistribution<d00P>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dP0P>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dM0P>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0PP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0MP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMP>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;

        case d00P:
            computeOutflowDistribution<d00M>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dP0M>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dM0M>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0PM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<d0MM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMPM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dPMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            computeOutflowDistribution<dMMM>(populationReferences, listIndices, populations, populationsNeighbor, cs, densityCorrection);
            break;
        default:
            break;
    }
}

}

//! \}
#endif