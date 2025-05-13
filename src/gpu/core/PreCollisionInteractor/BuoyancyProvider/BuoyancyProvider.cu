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

#include <stdexcept>

#include <cub/cub.cuh>
#include <cub/device/device_segmented_reduce.cuh>
#include <utility>
#include <helper_cuda.h>

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
                                             uint numberOfNodes)
{
    const uint nodeIndex = vf::cuda::get1DIndexFrom2DBlock();
    if (nodeIndex > numberOfNodes)
        return;

    forceZ[nodeIndex] += computeBuoyancyForce(factor, referenceTemperature[nodeIndex], temperature[nodeIndex]);
}

__global__ void computeBuoyancyPlanarAverage(ProfileParameters profileParameters, const real* temperature,
                                             real* referenceTemperature, const uint* numberOfNodesPerPlane, real factor,
                                             real* forceZ, uint numberOfNodes)
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


void findFirstIndicesInZDirection(std::vector<uint>& referencestateIndices, std::vector<real>& referenceStateCoordsZ,
                                  Parameter* para, int level)
{
    real lastZ = para->getParH(level)->coordinateZ[1];
    referenceStateCoordsZ.push_back(lastZ);
    referencestateIndices.push_back(c0o1);

    for (unsigned long long i = 1; i < para->getParH(level)->numberOfNodes; i++) {
        const real currentZ = para->getParH(level)->coordinateZ[i];
        if (currentZ > lastZ) {
            lastZ = currentZ;
            referenceStateCoordsZ.push_back(lastZ);
            referencestateIndices.push_back(i);
        }
    }
    referencestateIndices.push_back(para->getParH(level)->numberOfNodes);
}

ProfileParameters initializeProfileParameters(int level, Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    std::vector<uint> referenceStateIndices;
    std::vector<real> referenceStateCoordsZ;
    findFirstIndicesInZDirection(referenceStateIndices, referenceStateCoordsZ, para, level);
    ProfileParameters profileParameters(static_cast<uint>(referenceStateIndices.size()) - 1);
    cudaMemoryManager->cudaAllocBuoyancyProviderProfileParameters(&profileParameters);
    std::copy(referenceStateIndices.begin(), referenceStateIndices.end(), profileParameters.indicesHost);
    cudaMemoryManager->cudaCopyBuoyancyProviderProfileParametersHtoD(&profileParameters);
    return profileParameters;
}

void countNodesPerPlane(std::vector<uint>& nodesPerPlane, LBMSimulationParameter* para, const uint* indices,
                        uint numberOfPlanes)
{
    for (uint plane = 0; plane < numberOfPlanes; plane++) {
        uint nodes = 0;
        for (uint i = indices[plane]; i < indices[plane + 1]; i++) {
            if (para->typeOfGridNode[i] == GEO_FLUID)
                nodes++;
        }
        nodesPerPlane.push_back(nodes);
    }
}

BuoyancyProviderConstantValue::BuoyancyProviderConstantValue(SPtr<Parameter> parameter,
                                                             SPtr<CudaMemoryManager> cudaMemoryManager)
    : PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
{
    if (!parameter->getBuoyancyEnabled())
        throw std::runtime_error("BuoyancyProviderConstantValue: parameter needs to have buoyancy!");
}

void BuoyancyProviderConstantValue::init()
{
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderConstantValue::interact(int level, uint t)
{
    const real buoyancyFactor = para->getScaledBuoyancyFactor(level);
    const size_t numberOfNodes = para->getParD(level)->numberOfNodes;

    auto grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, numberOfNodes);
    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);

    computeBuoyancyConstantValue<<<grid.grid, grid.threads, 0, stream>>>(para->getParD(level)->localReferenceTemperature,
                                                                         para->getParD(level)->concentration, buoyancyFactor,
                                                                         para->getParD(level)->forceZ_SP, numberOfNodes);
    cudaStreamSynchronize(stream);
}

BuoyancyProviderPlanarAverage::BuoyancyProviderPlanarAverage(SPtr<Parameter> parameter,
                                                             SPtr<CudaMemoryManager> cudaMemoryManager)
    : PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
{
    if (!parameter->getBuoyancyEnabled())
        throw std::runtime_error("BuoyancyProviderPlanarAverage: parameter needs to have buyoancy!");
}

void BuoyancyProviderPlanarAverage::init()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        profileParameters.push_back(initializeProfileParameters(level, para.get(), cudaMemoryManager.get()));

        std::vector<uint> nodesPerPlane;

        countNodesPerPlane(nodesPerPlane, para->getParH(level).get(), profileParameters[level].indicesHost,
                           profileParameters[level].numberOfPlanes);
        reductionParameters.emplace_back(static_cast<uint>(nodesPerPlane.size()));


        CubDebug(cub::DeviceSegmentedReduce::Sum(
            reductionParameters[level].temporaryMemory, reductionParameters[level].sizeOfTemporaryMemory, para->getParD(level)->concentration,
            profileParameters[level].referenceTemperaturesDevice, reductionParameters[level].numberOfPlanes,
            profileParameters[level].indicesDevice, profileParameters[level].indicesDevice + 1));

        cudaMemoryManager->cudaAllocBuoyancyProviderReductionParameters(&reductionParameters[level]);
        std::copy(nodesPerPlane.begin(), nodesPerPlane.end(), reductionParameters[level].numberOfNodesPerPlaneHost);
        cudaMemoryManager->cudaCopyBuoyancyProviderReductionParametersHtoD(&reductionParameters[level]);
    }
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderPlanarAverage::interact(int level, uint t)
{
    auto& reductionParams = this->reductionParameters[level]; // must not be const
    const auto& profileParams = this->profileParameters[level];

    const real buoyancyFactor = para->getScaledBuoyancyFactor(level);
    const size_t numberOfNodes = para->getParD(level)->numberOfNodes;

    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);
    cub::DeviceSegmentedReduce::Sum(reductionParams.temporaryMemory, reductionParams.sizeOfTemporaryMemory,
                                    para->getParD(level)->concentration, profileParams.referenceTemperaturesDevice,
                                    static_cast<int>(reductionParams.numberOfPlanes), profileParams.indicesDevice,
                                    profileParams.indicesDevice + 1, stream);

    auto grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, numberOfNodes);
    computeBuoyancyPlanarAverage<<<grid.grid, grid.threads, 0, stream>>>(
        profileParameters[level], para->getParD(level)->concentration, para->getParD(level)->localReferenceTemperature,
        reductionParams.numberOfNodesPerPlaneDevice, buoyancyFactor, para->getParD(level)->forceZ_SP, numberOfNodes);
    cudaStreamSynchronize(stream);
    getLastCudaError("computeBuoyancyPlanarAverage kernel failed");
}

BuoyancyProviderPlanarAverage::~BuoyancyProviderPlanarAverage()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeBuoyancyProviderProfileParameters(&profileParameters[level]);
        cudaMemoryManager->cudaFreeBuoyancyProviderReductionParameters(&reductionParameters[level]);
    }
}

BuoyancyProviderPlanarAverageMultiGPU::BuoyancyProviderPlanarAverageMultiGPU(SPtr<Parameter> parameter,
                                                                             SPtr<CudaMemoryManager> cudaMemoryManager)
    : PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
{
    if (!parameter->getBuoyancyEnabled())
        throw std::runtime_error("BuoyancyProviderPlanarAverage: parameter needs to have buyoancy!");

    if (vf::parallel::MPICommunicator::getInstance()->getNumberOfProcesses() < 2)
        throw std::runtime_error("Only 1 process. No need to use MultiGPU.");
}

void BuoyancyProviderPlanarAverageMultiGPU::init()
{
    auto communicator = vf::parallel::MPICommunicator::getInstance();
    const int numberOfProcesses = communicator->getNumberOfProcesses();
    if (numberOfProcesses < 2)
        throw std::runtime_error("Only 1 process. No need to use MultiGPU.");

    for (int level = 0; level <= para->getMaxLevel(); level++) {
        profileParameters.push_back(initializeProfileParameters(level, para.get(), cudaMemoryManager.get()));
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
        cudaMemoryManager->cudaAllocBuoyancyProviderReductionParameters(&reductionParams);
        // fill numberOfNodes with dummy value to reuse the same kernel as the single gpu version
        std::fill(reductionParams.numberOfNodesPerPlaneHost,
                  reductionParams.numberOfNodesPerPlaneHost + reductionParams.numberOfPlanes, 1);
        cudaMemoryManager->cudaCopyBuoyancyProviderReductionParametersHtoD(&reductionParams);
        reductionParameters.push_back(reductionParams);
    }
    streamIndex = para->getStreamManager()->registerAndLaunchStream(CudaStreamIndex::BuoyancyProvider);
}

void BuoyancyProviderPlanarAverageMultiGPU::interact(int level, uint /**/)
{
    auto communicator = vf::parallel::MPICommunicator::getInstance();
    auto* profileParams = &this->profileParameters[level];
    auto& nodesPerPlane = this->totalNumberOfNodesPerPlane[level];

    const real buoyancyFactor = para->getScaledBuoyancyFactor(level);
    const size_t numberOfNodes = para->getParD(level)->numberOfNodes;
    const real* temperatureDevice = para->getParD(level)->concentration;

    std::vector<real> temperatureSums(profileParams->numberOfPlanes);
    std::vector<real> temperatureSumsAllProcesses;

    for (uint plane = 0; plane < profileParams->numberOfPlanes; plane++) {
        temperatureSums.push_back(thrust::reduce(thrust::device, temperatureDevice + profileParams->indicesHost[plane],
                                                 temperatureDevice + profileParams->indicesHost[plane + 1]));
    }

    communicator->allGather(temperatureSums, temperatureSumsAllProcesses);

    for (size_t plane = 0; plane < nodesPerPlane.size(); plane++) {
        if (nodesPerPlane[plane] == 0)
            continue;

        real temp = c0o1;
        for (int process = 0; process < communicator->getNumberOfProcesses(); process++) {
            temp += temperatureSumsAllProcesses[process * nodesPerPlane.size() + plane];
        }
        const real newRefTemp = temp / static_cast<real>((totalNumberOfNodesPerPlane[level])[plane]);
        profileParams->referenceTemperaturesHost[plane] = newRefTemp;
    }

    cudaMemoryManager->cudaCopyBuoyancyProviderProfileParametersHtoD(profileParams);
    auto* stream = para->getStreamManager()->getStream(CudaStreamIndex::BuoyancyProvider, streamIndex);

    auto grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, numberOfNodes);
    computeBuoyancyPlanarAverage<<<grid.grid, grid.threads, 0, stream>>>(
        profileParameters[level], temperatureDevice, para->getParD(level)->localReferenceTemperature,
        reductionParameters[level].numberOfNodesPerPlaneDevice, buoyancyFactor, para->getParD(level)->forceZ_SP,
        numberOfNodes);
    cudaStreamSynchronize(stream);
    getLastCudaError("computeBuoyancyPlanarAverage kernel failed");
}

BuoyancyProviderPlanarAverageMultiGPU::~BuoyancyProviderPlanarAverageMultiGPU()
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        cudaMemoryManager->cudaFreeBuoyancyProviderProfileParameters(&profileParameters[level]);
    }
}
