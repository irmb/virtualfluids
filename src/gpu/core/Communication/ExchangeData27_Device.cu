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

namespace vf::gpu {

__global__ void getSendFsPost27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* sendIndex, const uint buffmax, const uint* neighborX,
                                const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                                const bool isEvenTimestep, const bool diffOn)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();
    if (index >= buffmax)
        return;

    const ListIndices indices(sendIndex[index], neighborX, neighborY, neighborZ);

    const Distributions27 populationReferences = getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);
    const Distributions27 bufferReferences = getDistributionReferences27(bufferFs, buffmax, true);

    forEachDirection([&](auto dir) { (bufferReferences.f[dir])[index] = readFromInverseDirection<dir>(indices, populationReferences); });

    if(diffOn)
    {
        const Distributions27 populationReferencesAD = getDistributionReferences27(populationsAD, numberOfLBnodes, isEvenTimestep);
        const Distributions27 bufferReferencesAD = getDistributionReferences27(bufferAD, buffmax, true);

        forEachDirection([&](auto dir) { (bufferReferencesAD.f[dir])[index] = readFromInverseDirection<dir>(indices, populationReferencesAD); });
    }
}

__global__ void setRecvFsPost27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* recvIndex, const uint buffmax, const uint* neighborX,
                                const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                                const bool isEvenTimestep, const bool diffOn)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();
    if (index >= buffmax)
        return;

    const ListIndices indices(recvIndex[index], neighborX, neighborY, neighborZ);

    const Distributions27 populationReferences = getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);
    const Distributions27 bufferReferences = getDistributionReferences27(bufferFs, buffmax, true);

    forEachDirection([&](auto dir) { writeInInverseDirection<dir>((bufferReferences.f[dir])[index], indices, populationReferences); });

    if(diffOn)
    {
        const Distributions27 populationReferencesAD = getDistributionReferences27(populationsAD, numberOfLBnodes, isEvenTimestep);
        const Distributions27 bufferReferencesAD = getDistributionReferences27(bufferAD, buffmax, true);

        forEachDirection([&](auto dir) { writeInInverseDirection<dir>((bufferReferencesAD.f[dir])[index], indices, populationReferencesAD); });
    }
}

__global__ void getSendFsPre27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* sendIndex, const uint buffmax, const uint* neighborX,
                               const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                               const bool isEvenTimestep, const bool diffOn)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();
    if (index >= buffmax)
        return;

    const ListIndices indices(sendIndex[index], neighborX, neighborY, neighborZ);

    const Distributions27 populationReferences = getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);
    const Distributions27 bufferReferences = getDistributionReferences27(bufferFs, buffmax, true);

    forEachDirection([&](auto dir) { (bufferReferences.f[dir])[index] = readFromSameDirection<dir>(indices, populationReferences); });
    if(diffOn)
    {
        const Distributions27 populationReferencesAD = getDistributionReferences27(populationsAD, numberOfLBnodes, isEvenTimestep);
        const Distributions27 bufferReferencesAD = getDistributionReferences27(bufferAD, buffmax, true);

        forEachDirection([&](auto dir) { (bufferReferencesAD.f[dir])[index] = readFromSameDirection<dir>(indices, populationReferencesAD); });
    }

}

__global__ void setRecvFsPre27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* recvIndex, const uint buffmax, const uint* neighborX,
                               const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                               const bool isEvenTimestep, const bool diffOn)
{
    const uint index = vf::cuda::get1DIndexFrom2DBlock();
    if (index >= buffmax)
        return;

    const ListIndices indices(recvIndex[index], neighborX, neighborY, neighborZ);

    const Distributions27 populationReferences = getDistributionReferences27(DD, numberOfLBnodes, isEvenTimestep);
    const Distributions27 bufferReferences = getDistributionReferences27(bufferFs, buffmax, true);

    forEachDirection([&](auto dir) { writeInSameDirection<dir>((bufferReferences.f[dir])[index], indices, populationReferences); });

    if(diffOn)
    {
        const Distributions27 populationReferencesAD = getDistributionReferences27(populationsAD, numberOfLBnodes, isEvenTimestep);
        const Distributions27 bufferReferencesAD = getDistributionReferences27(bufferAD, buffmax, true);

        forEachDirection([&](auto dir) { writeInSameDirection<dir>((bufferReferencesAD.f[dir])[index], indices, populationReferencesAD); });
    }
}

void GetSendFsPreDev27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* sendIndex, const uint buffmax, const uint* neighborX,
                       const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                       const bool isEvenTimestep, const bool diffOn, const unsigned int numberOfThreads, cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    getSendFsPre27<<<grid.grid, grid.threads, 0, stream>>>(DD, bufferFs, populationsAD, bufferAD, sendIndex, buffmax, neighborX, neighborY, neighborZ,
                                                           numberOfLBnodes, isEvenTimestep, diffOn);
    getLastCudaError("getSendFsPre27 execution failed");
}

void GetSendFsPostDev27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* sendIndex, const uint buffmax, const uint* neighborX,
                        const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                        const bool isEvenTimestep, const bool diffOn, const unsigned int numberOfThreads, cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    getSendFsPost27<<<grid.grid, grid.threads, 0, stream>>>(DD, bufferFs, populationsAD, bufferAD, sendIndex, buffmax, neighborX, neighborY,
                                                            neighborZ, numberOfLBnodes, isEvenTimestep, diffOn);
    getLastCudaError("getSendFsPost27 execution failed");
}

void SetRecvFsPreDev27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* recvIndex, const uint buffmax, const uint* neighborX,
                       const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                       const bool isEvenTimestep, const bool diffOn, const unsigned int numberOfThreads, cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    setRecvFsPre27<<<grid.grid, grid.threads, 0, stream>>>(DD, bufferFs, populationsAD, bufferAD, recvIndex, buffmax, neighborX, neighborY, neighborZ,
                                                           numberOfLBnodes, isEvenTimestep, diffOn);
    getLastCudaError("setRecvFsPre27 execution failed");
}

void SetRecvFsPostDev27(real* DD, real* bufferFs, real* populationsAD, real* bufferAD, const uint* recvIndex, const uint buffmax, const uint* neighborX,
                        const uint* neighborY, const uint* neighborZ, const unsigned long long numberOfLBnodes,
                        const bool isEvenTimestep, const bool diffOn, const unsigned int numberOfThreads, cudaStream_t stream)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, buffmax);

    setRecvFsPost27<<<grid.grid, grid.threads, 0, stream>>>(DD, bufferFs, populationsAD, bufferAD, recvIndex, buffmax, neighborX, neighborY,
                                                            neighborZ, numberOfLBnodes, isEvenTimestep, diffOn);
    getLastCudaError("setRecvFsPost27 execution failed");
}

}
//! \}
