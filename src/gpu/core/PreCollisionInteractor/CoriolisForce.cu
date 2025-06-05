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
#include "Calculation/Calculation.h"
#include "CoriolisForce.h"

#include <cmath>

#include <cstdlib>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "cuda_helper/CudaGrid.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/Parameter/Parameter.h"
#include <cuda_helper/CudaIndexCalculation.h>

using namespace vf::basics::constant;

__global__ void computeCoriolis(unsigned long long numberOfNodes, const real* velocityX, const real* velocityY, real* forceX,
                                real* forceY, CoriolisForce::LevelData data, uint timestep, uint numberOfAccumulationSteps)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= numberOfNodes)
        return;
    const real coriolisForceX = -(data.velocityY - velocityY[nodeIndex]) * data.coriolisFrequency;
    const real coriolisForceY = (data.velocityX - velocityX[nodeIndex]) * data.coriolisFrequency;
    if (timestep % numberOfAccumulationSteps != 0) {
        data.forceAccumulatorX[nodeIndex] += coriolisForceX;
        data.forceAccumulatorY[nodeIndex] += coriolisForceY;
    } else {
        forceX[nodeIndex] += data.forceAccumulatorX[nodeIndex] + coriolisForceX;
        data.forceAccumulatorX[nodeIndex] = c0o1;
        forceY[nodeIndex] += data.forceAccumulatorY[nodeIndex] + coriolisForceY;
        data.forceAccumulatorY[nodeIndex] = c0o1;
    }
}

void CoriolisForce::init()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        auto& data = levelData.emplace_back();
        data.coriolisFrequency = coriolisParameter * para->getScaledTimeRatio(level);
        data.velocityX = std::cos(windDirection) * windSpeed / para->getScaledVelocityRatio(level);
        data.velocityY = std::sin(windDirection) * windSpeed / para->getScaledVelocityRatio(level);
        cudaMemoryManager->cudaAllocCoriolisForceData(this, level);
    }
}

void CoriolisForce::interact(int level, uint t)
{
    vf::cuda::CudaGrid grid(para->getParH(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    printf("%u %llu \n", para->getParH(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    auto parD = para->getParD(level);
    computeCoriolis<<<grid.grid, grid.threads>>>(parD->numberOfNodes, parD->velocityX, parD->velocityY, parD->forceX_SP,
                                                 parD->forceY_SP, levelData[level], t, numberOfAccumulationTimesteps);
}

CoriolisForce::LevelData& CoriolisForce::getLevelData(int level)
{
    return levelData[level];
}

CoriolisForce::~CoriolisForce()
{
    for (int level = 0; level <= para->getMaxLevel(); level++)
        cudaMemoryManager->cudaFreeCoriolisForceData(this, level);
}