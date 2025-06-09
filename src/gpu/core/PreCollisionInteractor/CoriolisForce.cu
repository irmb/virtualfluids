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

__global__ void computeCoriolis(unsigned long long numberOfNodes, const real* velocityX, const real* velocityY, real* forceX,
                                real* forceY, const real geostrophicWindX, const real geostrophicWindY,
                                const real coriolisParameter, CoriolisForce::LevelData levelData, const bool accumulate)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex >= numberOfNodes)
        return;
    const real newForceX = -(geostrophicWindY - velocityY[nodeIndex]) * coriolisParameter;
    const real newForceY = (geostrophicWindX - velocityX[nodeIndex]) * coriolisParameter;
    if(accumulate)
    {
        levelData.forceAccumulatorX[nodeIndex] += newForceX;
        levelData.forceAccumulatorY[nodeIndex] += newForceY;
    } else {
        forceX[nodeIndex] += levelData.forceAccumulatorX[nodeIndex] + newForceX;
        forceY[nodeIndex] += levelData.forceAccumulatorY[nodeIndex] + newForceY;
        levelData.forceAccumulatorX[nodeIndex] = c0o1;
        levelData.forceAccumulatorY[nodeIndex] = c0o1;
    }
}

void CoriolisForce::init()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        levelData.emplace_back();
        cudaMemoryManager->cudaAllocCoriolisForceData(this, level);
    }
}

void CoriolisForce::interact(int level, uint t)
{
    auto parD = para->getParD(level);
    vf::cuda::CudaGrid grid(parD->numberofthreads, parD->numberOfNodes);
    computeCoriolis<<<grid.grid, grid.threads>>>(parD->numberOfNodes, parD->velocityX, parD->velocityY, parD->forceX_SP,
                                                 parD->forceY_SP, geostrophicWindX, geostrophicWindY, coriolisParameter, levelData[level], t%nAccumulationSteps!=0);
}

CoriolisForce::~CoriolisForce()
{
    for (int level = 0; level <= para->getMaxLevel(); level++)
        cudaMemoryManager->cudaFreeCoriolisForceData(this, level);
}
