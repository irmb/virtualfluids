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
#include "Axis.h"
#include "DampingLayer.h"

#include <cmath>
#include <stdexcept>

#include <cuda.h>

#include <basics/DataTypes.h>

#include <cuda_helper/CudaIndexCalculation.h>

#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Utilities/KernelUtilities.h"
#include "gpu/cuda_helper/CudaGrid.h"

namespace vf::gpu {

__global__ void rayleighDampingLayerKernel(const real* coefficients, const real* velocity, real* force, const uint* indices,
                                           uint numberOfNodes)
{
    const uint localIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (localIndex >= numberOfNodes)
        return;
    const uint globalIndex = indices[localIndex];
    force[globalIndex] -= velocity[globalIndex] * coefficients[localIndex];
}

void DampingLayer::init()
{
    if (dampingLayerType == DampingLayerType::Rayleigh && !para->getIsBodyForce()) {
        throw std::runtime_error("Rayleigh damping requires body force to be enabled");
    }
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        makeDampingLayerData(level);
    }
}

void DampingLayer::makeDampingLayerData(int level)
{
    DampingLayerData& data = dampingLayerData.emplace_back();
    real* coordinate;
    std::vector<uint> indices;
    std::vector<real> dampingFactor, minimumValue;
    switch (direction) {
        case Axis::x:
            coordinate = para->getParH(level)->coordinateX;
            break;
        case Axis::y:
            coordinate = para->getParH(level)->coordinateY;
            break;
        case Axis::z:
            coordinate = para->getParH(level)->coordinateZ;
            break;
        default:
            throw std::runtime_error("Direction not implemented");
    }

    real start = startPosition, end = endPosition;
   
    for (uint index = 1; index < para->getParH(level)->numberOfNodes; index++) {
        const real coord = coordinate[index];
        if (coord < endPosition && coord > startPosition) {
            const real normalizedCoordinate = (coord - start) / (end - start);
            indices.push_back(index);
            dampingFactor.push_back(dampingFunction(normalizedCoordinate));
        }
    }

    data.numberOfNodes = static_cast<uint>(indices.size());

    cudaMemoryManager->cudaAllocDampingLayerData(this, level);
    std::copy(indices.begin(), indices.end(), data.indicesH);
    std::copy(dampingFactor.begin(), dampingFactor.end(), data.dampingCoefficientsH);
    cudaMemoryManager->cudaCopyDampingLayerDataHtoD(this, level);
}

void DampingLayer::interact(int level, uint t)
{
    auto& data = getDampingLayerData(level);
    if (data.numberOfNodes == 0)
        return;

    vf::cuda::CudaGrid grid(para->getParH(level)->numberofthreads, data.numberOfNodes);
    switch (dampingLayerType) {
        case DampingLayerType::Rayleigh:
            rayleighDampingLayerKernel<<<grid.grid, grid.threads>>>(
                data.dampingCoefficientsD, para->getParD(level)->velocityZ, para->getParD(level)->forceZ_SP, data.indicesD,
                data.numberOfNodes);
            break;
        default:
            throw std::runtime_error("DampingLayerType not implemented");
    }
}

DampingLayer::~DampingLayer()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeDampingLayerData(this, level);
    }
}

void DampingLayer::getTaggedFluidNodes(GridProvider* gridProvider)
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        auto& data = getDampingLayerData(level);
        std::vector<uint> dampingIndices(data.indicesH, data.indicesH + data.numberOfNodes);
        gridProvider->tagFluidNodeIndices(dampingIndices, CollisionTemplate::WriteMacroVars, level);
    }
}

}
