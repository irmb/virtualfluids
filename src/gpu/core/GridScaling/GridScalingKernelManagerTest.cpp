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
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Anna Wellmann, Martin Schönherr
//=======================================================================================
#include <gmock/gmock.h>
#include <stdexcept>

#include "GridScalingKernelManager.h"
#include "GridScaling/GridScalingFactory.h"
#include "Parameter/Parameter.h"
#include <basics/PointerDefinitions.h>

class GridScalingKernelManagerTest_Initialization : public testing::Test
{
protected:
    GridScalingFactory scalingFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
        para->setGridX({2, 8});
        para->setGridY({2, 8});
        para->setGridZ({2, 8});
        para->setDistX({0, 0});
        para->setDistY({0, 0});
        para->setDistZ({0, 0});
    }
};

TEST_F(GridScalingKernelManagerTest_Initialization, singleLevel_noScalingFactoryProvided_doesNotThrow)
{
    // only one level --> no scaling factory needed --> no error
    para->initLBMSimulationParameter();
    EXPECT_NO_THROW(GridScalingKernelManager(para, nullptr));
}

TEST_F(GridScalingKernelManagerTest_Initialization, singleLevel_scalingFactoryProvided_doesNotThrow)
{
    // only one level --> no scaling function needed --> no error
    para->initLBMSimulationParameter();
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

TEST_F(GridScalingKernelManagerTest_Initialization, singleLevel_scalingFactoryAndFunctionProvided_doesNotThrow)
{
    // only one level, but the user provided a scaling function anyway --> no error
    para->initLBMSimulationParameter();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevels_notScalingFactoryProvided_throws)
{
    // multiple levels, but the user forgot the scaling factory --> error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    EXPECT_THROW(GridScalingKernelManager(para, nullptr), std::runtime_error);
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelWithoutInterpolationNodes_noScalingFunctionProvided_doesNotThrow)
{
    // multiple levels, but no interpolation nodes specified --> no scaling function needed --> no error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelWithoutInterpolationNodes_scalingFunctionProvided_doesNotThrow)
{
    // multiple levels and NO interpolation nodes specified, but the user provided a scaling function anyway --> no error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleRhoSq);
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelWithInterpolationNodes_noScalingFunctionProvided_throws)
{
    // multiple levels and interpolation nodes specified, but the user forgot to set the scalingFunction --> error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    para->getParD(0)->fineToCoarse.numberOfCells = 100;
    EXPECT_THROW(GridScalingKernelManager(para, &scalingFactory), std::runtime_error);
}
