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

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using SlipBC = BoundaryConditionFactory::AdvectionDiffusionSlipVelocityBC;


template <size_t direction>
constexpr void setDistributionBC(uint nodeIndex, const SubgridDistances27& subgridDistances,
                                 const vf::gpu::ListIndices& listIndices,
                                 const DistributionReferences27 distributionReferences, const real* populationsConcentration,
                                 real fluxX, real fluxY, real fluxZ)
{
    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance < c0o1 || subgridDistance > c1o1)
        return;
    const real weight = vf::lbm::dir::getWeight<direction>();
    const real flux = vf::lbm::dir::getVelocity<direction>(fluxX, fluxY, fluxZ);
    const size_t inverseDirection = vf::lbm::dir::inverseDir<direction>();
    const uint writeIndex = listIndices.getIndex<inverseDirection>();
    const real population = populationsConcentration[direction];
    (distributionReferences.f[inverseDirection])[writeIndex] = population - c6o1 * weight * flux;
}

////////////////////////////////////////////////////////////////////////////////

constexpr real calculateDistributionInterpolated(real q, real weight, real v, real v_sq, real f, real finf,
                                                 real omegaDiffusivity, real jTangential, real concentration)
{
    const real feq = weight * concentration * (c1o1 + c3o1 * v + c9o2 * v * v * concentration - v_sq * concentration);
    return (c1o1 - q) / (c1o1 + q) * ((f - feq * omegaDiffusivity) / (c1o1 - omegaDiffusivity)) +
           (q * (f + finf) - c6o1 * weight * (jTangential)) / (c1o1 + q);
}
////////////////////////////////////////////////////////////////////////////////

template <size_t direction>
constexpr void
setDistributionBCInterpolated(uint nodeIndex, const SubgridDistances27& subgridDistances,
                              const vf::gpu::ListIndices& listIndices, const DistributionReferences27 distributionReferences,
                              const real* populationsConcentration, real velocityX, real velocityY, real velocityZ,
                              real relaxationFrequency, real cu_sq, real concentration, real fluxX, real fluxY, real fluxZ)
{
    using namespace vf::lbm::dir;
    const real subgridDistance = subgridDistances.q[direction][nodeIndex];
    if (subgridDistance < c0o1 || subgridDistance > c1o1)
        return;
    const size_t inverseDirection = inverseDir<direction>();
    const real weight = getWeight<direction>();
    const real flux = getVelocity<direction>(fluxX, fluxY, fluxZ);
    const uint writeIndex = listIndices.getIndex<inverseDirection>();
    const real population = populationsConcentration[direction];
    const real inversePopulation = populationsConcentration[inverseDirection];
    const real velocity = getVelocity<direction>(velocityX, velocityY, velocityZ);
    (distributionReferences.f[inverseDirection])[writeIndex] = calculateDistributionInterpolated(
        subgridDistance, weight, velocity, cu_sq, population, inversePopulation, relaxationFrequency, flux, concentration);
}

template <SlipBC slipBCType>
__global__ void AdvectionDiffusionSlipVelocity_Device(real* distributions, AdvectionDiffusionSlipVelocityBoundaryConditions bcParameters,
                                                      const real* density, const real* velocityX, const real* velocityY,
                                                      const real* velocityZ, const real* turbulentDiffusivity,
                                                      real diffusivity, real omegaDiffusivity,
                                                      const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                      unsigned long long numberOfLBnodes, bool isEvenTimestep)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= bcParameters.numberOfBCnodes)
        return;

    const bool withTurbulentViscosity = slipBCType == SlipBC::SlipVelocityTurbulentViscosityCompressible;

    Distributions27 distributionReferences =
        vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);

    const real normalX = bcParameters.normalX[nodeIndex];
    const real normalY = bcParameters.normalY[nodeIndex];
    const real normalZ = bcParameters.normalZ[nodeIndex];

    SubgridDistances27 subgridDistanceReferences;
    vf::gpu::getPointersToSubgridDistances(subgridDistanceReferences, bcParameters.q27[0], bcParameters.numberOfBCnodes);

    const uint indexOfBCnode = bcParameters.BCNodeIndices[nodeIndex];
    const vf::gpu::ListIndices listIndices(indexOfBCnode, neighborX, neighborY, neighborZ);

    const real vx1 = velocityX[indexOfBCnode];
    const real vx2 = velocityY[indexOfBCnode];
    const real vx3 = velocityZ[indexOfBCnode];

    real populationsConcentration[27];
    vf::gpu::getPostCollisionDistribution(populationsConcentration, distributionReferences, listIndices);
    vf::gpu::getPointersToDistributions(distributionReferences, distributions, numberOfLBnodes, !isEvenTimestep);

    const real concentration = vf::lbm::getDensity(populationsConcentration);
    // diffusive flux(?)
    const real fluxX = vf::lbm::getIncompressibleVelocityX1(populationsConcentration) - vx1 * concentration;
    const real fluxY = vf::lbm::getIncompressibleVelocityX2(populationsConcentration) - vx2 * concentration;
    const real fluxZ = vf::lbm::getIncompressibleVelocityX3(populationsConcentration) - vx3 * concentration;
    const real addedDiffusivity = withTurbulentViscosity ? turbulentDiffusivity[indexOfBCnode] : c0o1;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real effectiveDiffusivity = addedDiffusivity + diffusivity;
    const real normalFlux = fluxX * normalX + fluxY * normalY + fluxZ * normalZ;

    const real neumannFlux = -bcParameters.gradient[nodeIndex] * effectiveDiffusivity;

    const real fluxTangentialX = (fluxX - normalFlux * normalX) + neumannFlux * normalX;
    const real fluxTangentialY = (fluxY - normalFlux * normalY) + neumannFlux * normalY;
    const real fluxTangentialZ = (fluxZ - normalFlux * normalZ) + neumannFlux * normalZ;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (slipBCType == SlipBC::SlipVelocityBounceBack) {
        setDistributionBC<dP00>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dM00>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0P0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0M0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d00P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d00M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPP0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMM0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPM0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMP0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dP0P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dM0M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dP0M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dM0P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0PP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0MM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0PM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<d0MP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPPP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMPP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPMP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMMP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPPM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMPM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dPMM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);
        setDistributionBC<dMMM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                populationsConcentration, fluxTangentialX, fluxTangentialY, fluxTangentialZ);

    } else {
        const real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + density[indexOfBCnode]);
        const real relaxationFrequency = vf::lbm::calculateOmegaWithturbulentViscosity(omegaDiffusivity, addedDiffusivity);

        setDistributionBCInterpolated<dP00>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dM00>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0P0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0M0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d00P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d00M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPP0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMM0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPM0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMP0>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dP0P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dM0M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dP0M>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dM0P>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0PP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0MM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0PM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<d0MP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPPP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMPP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPMP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMMP>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPPM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMPM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dPMM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
        setDistributionBCInterpolated<dMMM>(nodeIndex, subgridDistanceReferences, listIndices, distributionReferences,
                                            populationsConcentration, vx1, vx2, vx3, relaxationFrequency, cu_sq,
                                            concentration, fluxX, fluxY, fluxZ);
    }
}

//! \}
