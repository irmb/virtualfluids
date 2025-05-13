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

#ifndef BUOYANCY_REFERENCE_PROVIDER_H
#define BUOYANCY_REFERENCE_PROVIDER_H

#include <functional>
#include <iostream>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>

#include "PreCollisionInteractor/PreCollisionInteractor.h"

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

class BuoyancyProviderConstantValue : public PreCollisionInteractor
{
public:
    BuoyancyProviderConstantValue(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager);

    ~BuoyancyProviderConstantValue() = default;

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider*) override {};

private:
    int streamIndex;
};

class BuoyancyProviderPlanarAverage : public PreCollisionInteractor
{
public:
    BuoyancyProviderPlanarAverage(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager);

    ~BuoyancyProviderPlanarAverage();

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider*) override {};

private:
    uint numberOfInitialReferenceValues;
    int streamIndex;
    std::vector<ProfileParameters> profileParameters;
    std::vector<ReductionParameters> reductionParameters;
};

class BuoyancyProviderPlanarAverageMultiGPU : public PreCollisionInteractor
{
public:
    BuoyancyProviderPlanarAverageMultiGPU(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager);

    ~BuoyancyProviderPlanarAverageMultiGPU();

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider*) override {};

private:
    uint numberOfInitialReferenceValues;
    int streamIndex;
    std::vector<ProfileParameters> profileParameters;
    std::vector<ReductionParameters> reductionParameters;
    std::vector<std::vector<uint>> totalNumberOfNodesPerPlane;
};

#endif