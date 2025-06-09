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
//! \file BuoyancyProvider.h
//! \author Henry Korb, Henrik Asmuth
//! \date 29/03/2023
//! \brief
//=======================================================================================

#ifndef BUOYANCY_PROVIDER_H
#define BUOYANCY_PROVIDER_H

#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>

#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/PreCollisionInteractor/PreCollisionInteractor.h"
#include "gpu/core/Parameter/Parameter.h"

class BuoyancyProvider : public PreCollisionInteractor
{
public:
    BuoyancyProvider(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager) : PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
    {
        if (!para->getBuoyancyEnabled())
            throw std::runtime_error("BuoyancyProvider: buoyancy needs to be enabled in Parameter!");
    }
    struct ProfileParameters
    {
        uint numberOfPlanes;
        real *referenceTemperaturesHost, *referenceTemperaturesDevice = nullptr;
        uint *indicesHost, *indicesDevice = nullptr;
        ProfileParameters(uint numberOfPlanes) : numberOfPlanes(numberOfPlanes)
        {
        }
    };

    struct ReductionParameters
    {
        uint numberOfPlanes;
        uint *numberOfNodesPerPlaneHost = nullptr, *numberOfNodesPerPlaneDevice = nullptr;
        size_t sizeOfTemporaryMemory = 0;
        void* temporaryMemory = nullptr;
        ReductionParameters(uint numberOfPlanes) : numberOfPlanes(numberOfPlanes)
        {
        }
    };

    ReductionParameters* getReductionParameter(int level)
    {
        return &reductionParameters[level];
    }
    ProfileParameters* getProfileParameter(int level)
    {
        return &profileParameters[level];
    }

    void initializeProfileParameters();

protected:
    std::vector<ProfileParameters> profileParameters;
    std::vector<ReductionParameters> reductionParameters;
};

class BuoyancyProviderConstantValue : public BuoyancyProvider
{

public:
    using BuoyancyProvider::BuoyancyProvider;
    ~BuoyancyProviderConstantValue() override = default;

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider*) override {};

private:
    int streamIndex;
};

class BuoyancyProviderPlanarAverage : public BuoyancyProvider
{
public:
    BuoyancyProviderPlanarAverage(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager);

    ~BuoyancyProviderPlanarAverage() override;

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* /**/) override {};

private:
    uint numberOfInitialReferenceValues;
    int streamIndex;
};

class BuoyancyProviderPlanarAverageMultiGPU : public BuoyancyProvider
{
public:
    BuoyancyProviderPlanarAverageMultiGPU(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager);

    ~BuoyancyProviderPlanarAverageMultiGPU() override;

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* /**/) override {};

private:
    uint numberOfInitialReferenceValues;
    int streamIndex;
    std::vector<std::vector<uint>> totalNumberOfNodesPerPlane;
};

#endif