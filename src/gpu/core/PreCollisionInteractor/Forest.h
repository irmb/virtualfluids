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
//! \ingroup gpu_core core
//! \{
#ifndef FOREST_H
#define FOREST_H

#include "Parameter/Parameter.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <stdexcept>

namespace vf::gpu {
class GridProvider;

enum class PlantAreaDensityFunction {
    Flat,
    BottomHeavy,
    TopHeavy
};

class Forest : public PreCollisionInteractor
{
public:
    Forest(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaMemoryManager, real height, real dragCoefficient,
            real plantAreaIndex, int level, PlantAreaDensityFunction plantAreaDensityFunction)
        : PreCollisionInteractor(para, cudaMemoryManager),
          height(height),
          dragCoefficient(dragCoefficient),
          plantAreaIndex(plantAreaIndex),
          level(level),
          plantAreaDensityFunction(plantAreaDensityFunction),
          leafAreaDensityH(nullptr),
          leafAreaDensityD(nullptr)
    {
        // if (level > para->getMaxLevel())
        //     throw std::runtime_error("specified level does not exist.");
    }

    ~Forest() override;
    void init() override;
    void interact(int level, uint t) override;
    void getTaggedFluidNodes(GridProvider* gridProvider) override;

    uint getNumberOfIndices() const
    {
        return this->numberOfIndices;
    }

private:
    void initForestIndices();
    void initForestVelocities();
    void initLeafAreaDensity();

public:
    uint* forestIndicesH;
    uint* forestIndicesD;
    real* forestVelocitiesXPreviousTimestepD;
    real* forestVelocitiesYPreviousTimestepD;
    real* forestVelocitiesZPreviousTimestepD;
    real* leafAreaDensityH;
    real* leafAreaDensityD;

protected:
    const real height;
    const real dragCoefficient;
    const real plantAreaIndex;
    const int level;
    PlantAreaDensityFunction plantAreaDensityFunction;
    uint numberOfIndices{0};
    real forceRatio;
};

}
#endif // FOREST_H

//! \}
