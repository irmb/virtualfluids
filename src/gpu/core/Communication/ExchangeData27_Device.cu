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
//! \addtogroup gpu_Communication Communication
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//======================================================================================

#include "ExchangeData27_Device.cuh"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "cuda_helper/CudaIndexCalculation.h"
#include "gpu/core/Utilities/KernelUtilities.h"

#include "Calculation/Calculation.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

constexpr void copyToBuffer(DistributionReferences27& buffer, const uint index, const real* populations)
{
   (buffer.f[dP00])[index] = populations[dP00];
   (buffer.f[dM00])[index] = populations[dM00];
   (buffer.f[d0P0])[index] = populations[d0P0];
   (buffer.f[d0M0])[index] = populations[d0M0];
   (buffer.f[d00P])[index] = populations[d00P];
   (buffer.f[d00M])[index] = populations[d00M];
   (buffer.f[dPP0])[index] = populations[dPP0];
   (buffer.f[dMM0])[index] = populations[dMM0];
   (buffer.f[dPM0])[index] = populations[dPM0];
   (buffer.f[dMP0])[index] = populations[dMP0];
   (buffer.f[dP0P])[index] = populations[dP0P];
   (buffer.f[dM0M])[index] = populations[dM0M];
   (buffer.f[dP0M])[index] = populations[dP0M];
   (buffer.f[dM0P])[index] = populations[dM0P];
   (buffer.f[d0PP])[index] = populations[d0PP];
   (buffer.f[d0MM])[index] = populations[d0MM];
   (buffer.f[d0PM])[index] = populations[d0PM];
   (buffer.f[d0MP])[index] = populations[d0MP];
   (buffer.f[d000])[index] = populations[d000];
   (buffer.f[dPPP])[index] = populations[dPPP];
   (buffer.f[dMMP])[index] = populations[dMMP];
   (buffer.f[dPMP])[index] = populations[dPMP];
   (buffer.f[dMPP])[index] = populations[dMPP];
   (buffer.f[dPPM])[index] = populations[dPPM];
   (buffer.f[dMMM])[index] = populations[dMMM];
   (buffer.f[dPMM])[index] = populations[dPMM];
   (buffer.f[dMPM])[index] = populations[dMPM];
}

constexpr void copyFromBuffer(const DistributionReferences27& buffer, const uint index, real* populations)
{
   populations[dP00] = (buffer.f[dP00])[index];
   populations[dM00] = (buffer.f[dM00])[index];
   populations[d0P0] = (buffer.f[d0P0])[index];
   populations[d0M0] = (buffer.f[d0M0])[index];
   populations[d00P] = (buffer.f[d00P])[index];
   populations[d00M] = (buffer.f[d00M])[index];
   populations[dPP0] = (buffer.f[dPP0])[index];
   populations[dMM0] = (buffer.f[dMM0])[index];
   populations[dPM0] = (buffer.f[dPM0])[index];
   populations[dMP0] = (buffer.f[dMP0])[index];
   populations[dP0P] = (buffer.f[dP0P])[index];
   populations[dM0M] = (buffer.f[dM0M])[index];
   populations[dP0M] = (buffer.f[dP0M])[index];
   populations[dM0P] = (buffer.f[dM0P])[index];
   populations[d0PP] = (buffer.f[d0PP])[index];
   populations[d0MM] = (buffer.f[d0MM])[index];
   populations[d0PM] = (buffer.f[d0PM])[index];
   populations[d0MP] = (buffer.f[d0MP])[index];
   populations[d000] = (buffer.f[d000])[index];
   populations[dPPP] = (buffer.f[dPPP])[index];
   populations[dMMP] = (buffer.f[dMMP])[index];
   populations[dPMP] = (buffer.f[dPMP])[index];
   populations[dMPP] = (buffer.f[dMPP])[index];
   populations[dPPM] = (buffer.f[dPPM])[index];
   populations[dMMM] = (buffer.f[dMMM])[index];
   populations[dPMM] = (buffer.f[dPMM])[index];
   populations[dMPM] = (buffer.f[dMPM])[index];
}

__global__ void getSendFsPost27(real* DD,
                                real* bufferFs,
                                const int* sendIndex,
                                const int buffmax,
                                const uint* neighborX,
                                const uint* neighborY,
                                const uint* neighborZ,
                                const unsigned long long numberOfLBnodes, 
                                const bool isEvenTimestep)
{
   const uint index = vf::cuda::get1DIndexFrom2DBlock();
   if(index >= buffmax) return;
   
   const vf::gpu::ListIndices indices(sendIndex[index], neighborX, neighborY, neighborZ);

   Distributions27 populationReferences = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes,isEvenTimestep);
   Distributions27 bufferReferences = vf::gpu::getDistributionReferences27(bufferFs, buffmax, true);

   real populations[27];
   vf::gpu::getPostCollisionDistribution(populations, populationReferences, indices);
   copyToBuffer(bufferReferences, index, populations);
}

__global__ void setRecvFsPost27(real* DD,
                                real* bufferFs,
                                const int* recvIndex,
                                const int buffmax,
                                const uint* neighborX,
                                const uint* neighborY,
                                const uint* neighborZ,
                                const unsigned long long numberOfLBnodes, 
                                const bool isEvenTimestep)
{
   const uint index = vf::cuda::get1DIndexFrom2DBlock();
   if(index >= buffmax) return;
   
   const vf::gpu::ListIndices indices(recvIndex[index], neighborX, neighborY, neighborZ);

   Distributions27 populationReferences = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes,isEvenTimestep);
   Distributions27 bufferReferences = vf::gpu::getDistributionReferences27(bufferFs, buffmax, true);
   real populations[27];

   copyFromBuffer(bufferReferences, index, populations);
   vf::gpu::setPostCollisionDistribution(populationReferences, indices, populations);
}

__global__ void getSendFsPre27(real* DD,
                               real* bufferFs,
                               const int* sendIndex,
                               const int buffmax,
                               const uint* neighborX,
                               const uint* neighborY,
                               const uint* neighborZ,
                               const unsigned long long numberOfLBnodes, 
                               const bool isEvenTimestep)
{
   const uint index = vf::cuda::get1DIndexFrom2DBlock();
   if(index >= buffmax) return;
   
   const vf::gpu::ListIndices indices(sendIndex[index], neighborX, neighborY, neighborZ);

   Distributions27 populationReferences = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes,isEvenTimestep);
   Distributions27 bufferReferences = vf::gpu::getDistributionReferences27(bufferFs, buffmax, true);
   real populations[27];
   vf::gpu::getPreCollisionDistribution(populations, populationReferences, indices);
   copyFromBuffer(bufferReferences, index, populations);
}

__global__ void setRecvFsPre27(real* DD,
                               real* bufferFs,
                               const int* recvIndex,
                               const int buffmax,
                               const uint* neighborX,
                               const uint* neighborY,
                               const uint* neighborZ,
                               const unsigned long long numberOfLBnodes, 
                               const bool isEvenTimestep)
{
   const uint index = vf::cuda::get1DIndexFrom2DBlock();
   if(index >= buffmax) return;
   
   const vf::gpu::ListIndices indices(recvIndex[index], neighborX, neighborY, neighborZ);

   Distributions27 populationReferences = vf::gpu::getDistributionReferences27(DD, numberOfLBnodes,isEvenTimestep);
   Distributions27 bufferReferences = vf::gpu::getDistributionReferences27(bufferFs, buffmax, true);
   real populations[27];

   copyFromBuffer(bufferReferences, index, populations);
   vf::gpu::setPreCollisionDistribution(populationReferences, indices, populations);
}




void GetSendFsPreDev27(
    real* DD,
    real* bufferFs,
    const int* sendIndex,
    const int buffmax,
    const uint* neighborX,
    const uint* neighborY,
    const uint* neighborZ,
    const unsigned long long numberOfLBnodes,
    const bool isEvenTimestep,
    const unsigned int numberOfThreads,
    cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    getSendFsPre27<<< grid.grid, grid.threads, 0, stream >>>(
        DD,
        bufferFs,
        sendIndex,
        buffmax,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("getSendFsPre27 execution failed");
}

void GetSendFsPostDev27(
    real* DD,
    real* bufferFs,
    const int* sendIndex,
    const int buffmax,
    const uint* neighborX,
    const uint* neighborY,
    const uint* neighborZ,
    const unsigned long long numberOfLBnodes,
    const bool isEvenTimestep,
    const unsigned int numberOfThreads,
    cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    getSendFsPost27<<< grid.grid, grid.threads, 0, stream >>>(
        DD,
        bufferFs,
        sendIndex,
        buffmax,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("getSendFsPost27 execution failed");
}

void SetRecvFsPreDev27(
    real* DD,
    real* bufferFs,
    const int* recvIndex,
    const int buffmax,
    const uint* neighborX,
    const uint* neighborY,
    const uint* neighborZ,
    const unsigned long long numberOfLBnodes,
    const bool isEvenTimestep,
    const unsigned int numberOfThreads,
    cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    setRecvFsPre27<<< grid.grid, grid.threads, 0, stream >>>(
        DD,
        bufferFs,
        recvIndex,
        buffmax,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("setRecvFsPre27 execution failed");
}

void SetRecvFsPostDev27(
    real* DD,
    real* bufferFs,
    const int* recvIndex,
    const int buffmax,
    const uint* neighborX,
    const uint* neighborY,
    const uint* neighborZ,
    const unsigned long long numberOfLBnodes,
    const bool isEvenTimestep,
    const unsigned int numberOfThreads,
    cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    setRecvFsPost27<<< grid.grid, grid.threads, 0, stream >>>(
        DD,
        bufferFs,
        recvIndex,
        buffmax,
        neighborX,
        neighborY,
        neighborZ,
        numberOfLBnodes,
        isEvenTimestep);
    getLastCudaError("setRecvFsPost27 execution failed");
}


//! \}
