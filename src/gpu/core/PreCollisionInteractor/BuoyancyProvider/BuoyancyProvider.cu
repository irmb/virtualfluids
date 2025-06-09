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
#include "BuoyancyProvider.h"

#include <algorithm>
#include <stdexcept>

#include <cub/cub.cuh>
#include <cub/device/device_segmented_reduce.cuh>
#include <helper_cuda.h>
#include <utility>

#include <thrust/execution_policy.h>
#include <thrust/reduce.h>

#include <basics/DataTypes.h>

#include <parallel/MPICommunicator.h>

#include "gpu/cuda_helper/CudaGrid.h"
#include "gpu/cuda_helper/CudaIndexCalculation.h"

#include "gpu/core/Calculation/Calculation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/Cuda/CudaStreamManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Parameter/Parameter.h"

constexpr uint getPlaneIndex(uint nodeIndex, const uint* referenceIndices)
{
    const uint index = nodeIndex + 1;
    const uint indicesPerLayer = referenceIndices[1] - 1;
    uint valueIndex = index / indicesPerLayer;

    while (index >= referenceIndices[valueIndex + 1]) {
        valueIndex++;
    }
    while (index < referenceIndices[valueIndex]) {
        valueIndex--;
    }
    return valueIndex;
}

constexpr real computeBuoyancyForce(real factor, real referenceTemperature, real temperature)
{
    return factor * (temperature - referenceTemperature);
}

__global__ void computeBuoyancyConstantValue(real* referenceTemperature, const real* temperature, real factor, real* forceZ,
                                             unsigned long long numberOfNodes)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex > numberOfNodes)
        return;

    forceZ[nodeIndex] += computeBuoyancyForce(factor, referenceTemperature[nodeIndex], temperature[nodeIndex]);
}

__global__ void computeBuoyancyPlanarAverage(BuoyancyProvider::ProfileParameters profileParameters, const real* temperature,
                                             real* referenceTemperature, const uint* numberOfNodesPerPlane, real factor,
                                             real* forceZ, unsigned long long numberOfNodes)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex > numberOfNodes)
        return;

    const uint planeIndex = getPlaneIndex(nodeIndex, profileParameters.indicesDevice);
    const real numberOfNodesInPlane = static_cast<real>(numberOfNodesPerPlane[planeIndex]);

    const real newReferenceTemperature = profileParameters.referenceTemperaturesDevice[planeIndex] / numberOfNodesInPlane;
    referenceTemperature[nodeIndex] = newReferenceTemperature;
    forceZ[nodeIndex] += computeBuoyancyForce(factor, newReferenceTemperature, temperature[nodeIndex]);
}

void findFirstIndicesInZDirection(std::vector<uint>& referenceStateIndices, std::vector<real>& referenceStateCoordsZ,
                                  Parameter* para, int level)
{
    real lastZ = para->getParH(level)->coordinateZ[1];
    referenceStateCoordsZ.push_back(lastZ);
    referenceStateIndices.push_back(c0o1);

    for (unsigned long long i = 1; i < para->getParH(level)->numberOfNodes; i++) {
        const real currentZ = para->getParH(level)->coordinateZ[i];
        if (currentZ > lastZ) {
            lastZ = currentZ;
            referenceStateCoordsZ.push_back(lastZ);
            referenceStateIndices.push_back(i);
        }
    }
    referenceStateIndices.push_back(para->getParH(level)->numberOfNodes);
}

void BuoyancyProvider::initializeProfileParameters()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        std::vector<uint> referenceStateIndices;
        std::vector<real> referenceStateCoordsZ;
        findFirstIndicesInZDirection(referenceStateIndices, referenceStateCoordsZ, para.get(), level);
        BuoyancyProvider::ProfileParameters profileParameters(static_cast<uint>(referenceStateIndices.size()) - 1);
        cudaMemoryManager->cudaAllocBuoyancyProviderProfileParameters(this, level);
        std::copy(referenceStateIndices.begin(), referenceStateIndices.end(), profileParameters.indicesHost);
        cudaMemoryManager->cudaCopyBuoyancyProviderProfileParametersHtoD(this, level);
    }
}

void countNodesPerPlane(std::vector<uint>& nodesPerPlane, LBMSimulationParameter* para, const uint* indices,
                        uint numberOfPlanes)
{
    for (uint plane = 0; plane < numberOfPlanes; plane++) {
        nodesPerPlane.push_back(std::count_if(para->typeOfGridNode + indices[plane],
                                              para->typeOfGridNode + indices[plane + 1],
                                              [](auto type) { return type == GEO_FLUID; }));
    }
}

void BuoyancyProviderConstantValue::init()
{
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderConstantValue::interact(int level, uint /**/)
{
    auto& parD = para->getParDeviceAsReference(level);

    vf::cuda::CudaGrid grid(parD.numberofthreads, uint(parD.numberOfNodes));
    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);

    computeBuoyancyConstantValue<<<grid.grid, grid.threads, 0, stream>>>(parD.localReferenceTemperature, parD.concentration,
                                                                         para->getScaledBuoyancyFactor(level),
                                                                         parD.forceZ_SP, parD.numberOfNodes);
    cudaStreamSynchronize(stream);
}

void BuoyancyProviderPlanarAverage::init()
{
    initializeProfileParameters();
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        std::vector<uint> nodesPerPlane;
        nodesPerPlane.reserve(profileParameters[level].numberOfPlanes);

        countNodesPerPlane(nodesPerPlane, para->getParH(level).get(), profileParameters[level].indicesHost,
                           profileParameters[level].numberOfPlanes);
        reductionParameters.emplace_back(static_cast<uint>(nodesPerPlane.size()));

        CubDebug(cub::DeviceSegmentedReduce::Sum(
            reductionParameters[level].temporaryMemory, reductionParameters[level].sizeOfTemporaryMemory,
            para->getParD(level)->concentration, profileParameters[level].referenceTemperaturesDevice,
            reductionParameters[level].numberOfPlanes, profileParameters[level].indicesDevice,
            profileParameters[level].indicesDevice + 1));

        cudaMemoryManager->cudaAllocBuoyancyProviderReductionParameters(this, level);
        std::copy(nodesPerPlane.begin(), nodesPerPlane.end(), reductionParameters[level].numberOfNodesPerPlaneHost);
        cudaMemoryManager->cudaCopyBuoyancyProviderReductionParametersHtoD(this, level);
    }
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderPlanarAverage::interact(int level, uint /**/)
{
    auto& reductionParams = this->reductionParameters[level]; // must not be const
    const auto& profileParams = this->profileParameters[level];
    auto& parD = para->getParDeviceAsReference(level);

    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);
    cub::DeviceSegmentedReduce::Sum(reductionParams.temporaryMemory, reductionParams.sizeOfTemporaryMemory,
                                    parD.concentration, profileParams.referenceTemperaturesDevice,
                                    static_cast<int>(reductionParams.numberOfPlanes), profileParams.indicesDevice,
                                    profileParams.indicesDevice + 1, stream);

    vf::cuda::CudaGrid grid(parD.numberofthreads, uint(parD.numberOfNodes));
    computeBuoyancyPlanarAverage<<<grid.grid, grid.threads, 0, stream>>>(
        profileParameters[level], parD.concentration, parD.localReferenceTemperature,
        reductionParams.numberOfNodesPerPlaneDevice, para->getBuoyancyFactor(), parD.forceZ_SP, parD.numberOfNodes);
    cudaStreamSynchronize(stream);
    getLastCudaError("computeBuoyancyPlanarAverage kernel failed");
}

BuoyancyProviderPlanarAverage::~BuoyancyProviderPlanarAverage()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeBuoyancyProviderProfileParameters(this, level);
        cudaMemoryManager->cudaFreeBuoyancyProviderReductionParameters(this, level);
    }
}

BuoyancyProviderPlanarAverageMultiGPU::BuoyancyProviderPlanarAverageMultiGPU(SPtr<Parameter> parameter,
                                                                             SPtr<CudaMemoryManager> cudaMemoryManager)
    : BuoyancyProvider(std::move(parameter), std::move(cudaMemoryManager))
{
    if (vf::parallel::MPICommunicator::getInstance()->getNumberOfProcesses() < 2)
        throw std::runtime_error("Only 1 process. No need to use MultiGPU.");
}

void BuoyancyProviderPlanarAverageMultiGPU::init()
{
    auto communicator = vf::parallel::MPICommunicator::getInstance();
    const int numberOfProcesses = communicator->getNumberOfProcesses();
    if (numberOfProcesses < 2)
        throw std::runtime_error("Only 1 process. No need to use MultiGPU.");

    initializeProfileParameters();

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        std::vector<uint> nodesPerPlane;
        countNodesPerPlane(nodesPerPlane, para->getParH(level).get(), profileParameters[level].indicesHost,
                           profileParameters[level].numberOfPlanes);

        std::vector<uint> nodesPerPlaneAllProcesses;
        communicator->allGather(nodesPerPlane, nodesPerPlaneAllProcesses);

        for (size_t plane = 0; plane < nodesPerPlane.size(); plane++) {
            int numberOfNodes = 0;
            for (int process = 0; process < numberOfProcesses; process++) {
                numberOfNodes += nodesPerPlaneAllProcesses[process * nodesPerPlane.size() + plane];
            }
            nodesPerPlane[plane] = numberOfNodes;
        }
        totalNumberOfNodesPerPlane.emplace_back(nodesPerPlane);

        ReductionParameters reductionParams(profileParameters.back().numberOfPlanes);
        cudaMemoryManager->cudaAllocBuoyancyProviderReductionParameters(this, level);
        // fill numberOfNodes with dummy value to reuse the same kernel as the single gpu version
        std::fill(reductionParams.numberOfNodesPerPlaneHost,
                  reductionParams.numberOfNodesPerPlaneHost + reductionParams.numberOfPlanes, 1);
        cudaMemoryManager->cudaCopyBuoyancyProviderReductionParametersHtoD(this, level);
        reductionParameters.push_back(reductionParams);
    }
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderPlanarAverageMultiGPU::interact(int level, uint /**/)
{
    auto communicator = vf::parallel::MPICommunicator::getInstance();
    auto* profileParams = &this->profileParameters[level];
    auto& nodesPerPlane = this->totalNumberOfNodesPerPlane[level];
    auto& parD = para->getParDeviceAsReference(level);

    std::vector<real> temperatureSums(profileParams->numberOfPlanes);

    for (uint plane = 0; plane < profileParams->numberOfPlanes; plane++) {
        temperatureSums.push_back(thrust::reduce(parD.concentration + profileParams->indicesHost[plane],
                                                 parD.concentration + profileParams->indicesHost[plane + 1]));
    }

    communicator->allReduceSum(temperatureSums);

    for (size_t plane = 0; plane < nodesPerPlane.size(); plane++) {
        if (nodesPerPlane[plane] == 0)
            continue;
        profileParams->referenceTemperaturesHost[plane] = temperatureSums[plane] / real(nodesPerPlane[plane]);
    }

    cudaMemoryManager->cudaCopyBuoyancyProviderProfileParametersHtoD(this, level);
    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);

    vf::cuda::CudaGrid grid(parD.numberofthreads, uint(parD.numberOfNodes));
    computeBuoyancyPlanarAverage<<<grid.grid, grid.threads, 0, stream>>>(
        profileParameters[level], parD.concentration, parD.localReferenceTemperature,
        reductionParameters[level].numberOfNodesPerPlaneDevice, para->getScaledBuoyancyFactor(level), parD.forceZ_SP,
        parD.numberOfNodes);
    cudaStreamSynchronize(stream);
    getLastCudaError("computeBuoyancyPlanarAverage kernel failed");
}

BuoyancyProviderPlanarAverageMultiGPU::~BuoyancyProviderPlanarAverageMultiGPU()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeBuoyancyProviderProfileParameters(this, level);
    }
}
