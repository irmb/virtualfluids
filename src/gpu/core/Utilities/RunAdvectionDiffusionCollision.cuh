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
//! \addtogroup gpu_Utilities Utilities
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb
//=======================================================================================
#include <basics/DataTypes.h>

#include <lbm/advectionDiffusion/collision/CollisionParameter.h>
#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/advectionDiffusion/TurbulentDiffusivity.h>

#include "KernelUtilities.h"
#include "core/Calculation/Calculation.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"
#include "core/Parameter/Parameter.h"

namespace vf::gpu::ad
{

struct GPUCollisionParameters
{
    unsigned long long numberOfLBnodes;
    bool isEvenTimestep;
    real relaxationFrequency;
    const uint* typeOfGridNode;
    real* distributions;
    const uint *neighborX, *neighborY, *neighborZ;
    const real *velocityX, *velocityY, *velocityZ;
    real* concentration;
    real turbulentPrandtlNumber;
    const real* turbulentViscosity;
    real* turbulentDiffusivity;
    const uint* indices;
    uint numberOfFluidNodes;
};

inline GPUCollisionParameters getCollisionParameter(LBMSimulationParameter* parD, real turbulentPrandtlNumber,
                                                        const uint* indices, uint sizeIndices)
{
    return { parD->numberOfNodes,
             parD->isEvenTimestep,
             parD->omegaDiffusivity,
             parD->typeOfGridNode,
             parD->distributionsAD.f[0],
             parD->neighborX,
             parD->neighborY,
             parD->neighborZ,
             parD->velocityX,
             parD->velocityY,
             parD->velocityZ,
             parD->concentration,
             turbulentPrandtlNumber,
             parD->turbulentViscosity,
             parD->turbulentDiffusivity,
             indices,
             sizeIndices };
}

template <typename CollisionFunctor, vf::lbm::advection_diffusion::TurbulenceModel turbulenceModel>
__global__ void runCollisionAdvectionDiffusion(CollisionFunctor collision, GPUCollisionParameters collisionParameter)
{
    using namespace vf::lbm::advection_diffusion;
    using namespace vf::basics::constant;
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();

    if (nodeIndex >= collisionParameter.numberOfFluidNodes)
        return;
        
    const uint k_000 = collisionParameter.indices[nodeIndex];
    if (collisionParameter.typeOfGridNode[k_000] != GEO_FLUID) {
        collisionParameter.concentration[k_000] = c0o1;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    Distributions27 distAD = getDistributionReferences27(
        collisionParameter.distributions, collisionParameter.numberOfLBnodes, collisionParameter.isEvenTimestep);
    const ListIndices listIndices(k_000, collisionParameter.neighborX, collisionParameter.neighborY,
                                  collisionParameter.neighborZ);

    ADCollisionParameter para;
    getPreCollisionDistribution(para.distribution, distAD, listIndices);

    switch (turbulenceModel) {
        case TurbulenceModel::None:
            para.omega = collisionParameter.relaxationFrequency;
            break;
        case TurbulenceModel::Default: {
            const real turbulentDiffusivity = calcTurbulentDiffusivityDefault(
                collisionParameter.turbulentViscosity[k_000], collisionParameter.turbulentPrandtlNumber);
            para.omega =
                vf::lbm::calculateOmegaWithTurbulentViscosity(collisionParameter.relaxationFrequency, turbulentDiffusivity);
            collisionParameter.turbulentDiffusivity[k_000] = turbulentDiffusivity;
        } break;
        case TurbulenceModel::Moeng:
        case TurbulenceModel::AMDStratified:
            para.omega = vf::lbm::calculateOmegaWithTurbulentViscosity(collisionParameter.relaxationFrequency,
                                                                       collisionParameter.turbulentDiffusivity[k_000]);
            break;
    }

    para.velocityX = collisionParameter.velocityX[k_000];
    para.velocityY = collisionParameter.velocityY[k_000];
    para.velocityZ = collisionParameter.velocityZ[k_000];

    collision(para);

    setPostCollisionDistribution(distAD, listIndices, para.distribution);

    collisionParameter.concentration[k_000] = para.concentration;
}

} // namespace vf::gpu