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

#include <basics/constants/NumericConstants.h>
#include <lbm/MacroscopicQuantities.h>
#include <lbm/collision/TurbulentViscosity.h>

#include "gpu/core/BoundaryConditions/AdvectionDiffusion/Utilities.cuh"
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using SlipBC = BoundaryConditionFactory::AdvectionDiffusionSlipVelocityBC;

template <SlipBC slipBCType>
__global__ void AdvectionDiffusionSlipVelocity_Device(real* populationsArray,
                                                      const AdvectionDiffusionSlipVelocityBoundaryConditions bcParameters,
                                                      const real* density, const real* velocityX, const real* velocityY,
                                                      const real* velocityZ, const real* turbulentDiffusivity,
                                                      real diffusivity, const real omegaDiffusivity, const uint* neighborX,
                                                      const uint* neighborY, const uint* neighborZ,
                                                      const unsigned long long numberOfLBnodes, const bool isEvenTimestep)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    const bool withTurbulentViscosity = slipBCType == SlipBC::SlipVelocityTurbulentViscosityCompressible;

    Distributions27 populationReferences =
        vf::gpu::getDistributionReferences27(populationsArray, numberOfLBnodes, isEvenTimestep);

    const real normalX = bcParameters.normalX[nodeIndex];
    const real normalY = bcParameters.normalY[nodeIndex];
    const real normalZ = bcParameters.normalZ[nodeIndex];

    SubgridDistances27 subgridDistances;
    vf::gpu::getPointersToSubgridDistances(subgridDistances, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint indexOfBCnode = bcParameters.BCNodeIndices[nodeIndex];
    const vf::gpu::ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    const real vx1 = velocityX[indexOfBCnode];
    const real vx2 = velocityY[indexOfBCnode];
    const real vx3 = velocityZ[indexOfBCnode];

    real populations[27];
    vf::gpu::getPostCollisionDistribution(populations, populationReferences, listIndices);
    vf::gpu::getPointersToDistributions(populationReferences, populationsArray, numberOfLBnodes, !isEvenTimestep);

    const real concentration = vf::lbm::getDensity(populations);
    // diffusive flux
    const real diffusiveFluxX = vf::lbm::getIncompressibleVelocityX1(populations) - vx1 * concentration;
    const real diffusiveFluxY = vf::lbm::getIncompressibleVelocityX2(populations) - vx2 * concentration;
    const real diffusiveFluxZ = vf::lbm::getIncompressibleVelocityX3(populations) - vx3 * concentration;

    const real addedDiffusivity = withTurbulentViscosity ? turbulentDiffusivity[indexOfBCnode] : c0o1;
    const real effectiveDiffusivity = addedDiffusivity + diffusivity;

    const real normalFlux = diffusiveFluxX * normalX + diffusiveFluxY * normalY + diffusiveFluxZ * normalZ;

    const real neumannFlux = -bcParameters.gradient[nodeIndex] * effectiveDiffusivity;

    const real fluxX = (diffusiveFluxX - normalFlux * normalX) + neumannFlux * normalX;
    const real fluxY = (diffusiveFluxY - normalFlux * normalY) + neumannFlux * normalY;
    const real fluxZ = (diffusiveFluxZ - normalFlux * normalZ) + neumannFlux * normalZ;

    if (slipBCType == SlipBC::SlipVelocityBounceBack) {
        writePopulationsBounceBackWithFlux(nodeIndex, subgridDistances, listIndices, populationReferences, populations,
                                           fluxX, fluxY, fluxZ);
    } else {
        const real drho = density[indexOfBCnode];
        const real relaxationFrequency = vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivity, addedDiffusivity);
        writePopulationsInterpolatedWithFlux(nodeIndex, subgridDistances, listIndices, populationReferences, populations,
                                             vx1, vx2, vx3, drho, relaxationFrequency, concentration, fluxX, fluxY, fluxZ);
    }
}

//! \}
