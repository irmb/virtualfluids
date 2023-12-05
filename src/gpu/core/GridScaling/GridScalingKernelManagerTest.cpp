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
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);
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

