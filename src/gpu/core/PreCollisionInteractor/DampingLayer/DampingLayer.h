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
//! \file DampingLayer.h
//! \author Henry Korb
//! \date 06/02/2024
//! \brief
//=======================================================================================

#ifndef DAMPING_LAYER_H
#define DAMPING_LAYER_H

#include <functional>
#include <utility>
#include <vector>

#include <basics/DataTypes.h>
#include <logger/Logger.h>

#include "Axis.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"

namespace vf::gpu {

class DampingLayer : public PreCollisionInteractor
{

public:
    struct DampingLayerData;

    enum class DampingLayerType { Rayleigh };

public:
    DampingLayer(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager, DampingLayerType dampingLayerType,
                 Axis direction, std::function<real(real)> dampingFunction, real startPosition, real endPosition)
        : dampingLayerType(dampingLayerType), direction(direction), dampingFunction(std::move(dampingFunction)),
          startPosition(startPosition), endPosition(endPosition),
          PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
    {
        switch (dampingLayerType) {
            case DampingLayerType::Rayleigh:
                VF_LOG_INFO("Using Rayleigh Damping");
                break;
        }
        switch (direction) {
            case Axis::x:
                VF_LOG_INFO("Damping from {} to {} in X direction", startPosition, endPosition);
                break;
            case Axis::y:
                VF_LOG_INFO("Damping from {} to {} in Y direction", startPosition, endPosition);
                break;
            case Axis::z:
                VF_LOG_INFO("Damping from {} to {} in Z direction", startPosition, endPosition);
                break;
        }
    };

    ~DampingLayer() override;

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override;
    DampingLayerData& getDampingLayerData(int level)
    {
        return dampingLayerData[level];
    }

private:
    void makeDampingLayerData(int level);

    DampingLayerType dampingLayerType;
    Axis direction;
    real startPosition, endPosition;
    std::function<real(real)> dampingFunction;
    std::vector<DampingLayerData> dampingLayerData;
};

struct DampingLayer::DampingLayerData
{
    uint numberOfNodes;
    uint *indicesH, *indicesD;
    real *dampingCoefficientsH, *dampingCoefficientsD;
    real *minimumValueH, *minimumValueD;
};

}

#endif