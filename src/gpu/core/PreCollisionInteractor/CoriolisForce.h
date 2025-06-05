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
//! \addtogroup gpu_PreCollisionInteractor PreCollisionInteractor
//! \{
//! \author Henry Korb
//! \date 27/02/2024
//=======================================================================================

#ifndef CORIOLIS_FORCE_H_
#define CORIOLIS_FORCE_H_

#include "PreCollisionInteractor.h"
#include "Parameter/Parameter.h"

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <logger/Logger.h>
#include <stdexcept>
#include <vector>

class CoriolisForce : public PreCollisionInteractor
{
public:
    CoriolisForce(SPtr<Parameter> parameter, SPtr<CudaMemoryManager> cudaMemoryManager, real windSpeed, real windDirection, real coriolisParameter,
                  uint numberOfAccumulationTimesteps = 100)
        : numberOfAccumulationTimesteps(numberOfAccumulationTimesteps), windSpeed(windSpeed), windDirection(windDirection), coriolisParameter(coriolisParameter),
          PreCollisionInteractor(std::move(parameter), std::move(cudaMemoryManager))
    {
        VF_LOG_INFO("using Coriolis Force with wind speed {}, wind direction {} and coriolis parameter {}", windSpeed, windDirection, coriolisParameter);
        if(!para->getIsBodyForce())
            throw std::runtime_error("Coriolis force needs body force.");
        para->setAllNodesAllFeatures(true);
    }

    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override {};
    ~CoriolisForce() override;
    struct LevelData;
    LevelData& getLevelData(int level);
private:

    const uint numberOfAccumulationTimesteps;
    const real windSpeed, windDirection, coriolisParameter;
    std::vector<LevelData> levelData;
};

struct CoriolisForce::LevelData {
    real velocityX, velocityY, coriolisFrequency;
    real* forceAccumulatorX, *forceAccumulatorY;
};


#endif //! \}