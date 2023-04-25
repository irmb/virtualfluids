#include <gmock/gmock.h>
#include <stdexcept>

#include "GridScalingKernelManager.h"
#include "Factories/GridScalingFactory.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"

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

TEST_F(GridScalingKernelManagerTest_Initialization, singleLevel)
{
    // only one level --> no scaling factory needed --> no error
    para->initLBMSimulationParameter();
    para->getParD(0)->fineToCoarse.numberOfCells = 100;
    EXPECT_NO_THROW(GridScalingKernelManager(para, nullptr));

    // only one level --> no scaling function needed --> no error
    para->initLBMSimulationParameter();
    para->getParD(0)->fineToCoarse.numberOfCells = 100;
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));

    // only one level, but the user provided a scaling function anyway --> no error
    para->initLBMSimulationParameter();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelNoScalingFunction)
{
    // multiple levels, but no interpolation nodes specified --> no scaling function needed --> no error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));

    // multiple levels and interpolation nodes specified, but the user forgot to set the scalingFunction --> error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    para->getParD(0)->fineToCoarse.numberOfCells = 100;
    EXPECT_THROW(GridScalingKernelManager(para, &scalingFactory), std::runtime_error);
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelNoScalingFactory)
{
    // multiple levels, but the user forgot the scaling factory --> error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    EXPECT_THROW(GridScalingKernelManager(para, nullptr), std::runtime_error);
}

TEST_F(GridScalingKernelManagerTest_Initialization, multipleLevelUnnecessaryScalingFunction)
{
    // multiple levels and NO interpolation nodes specified, but the user provided a scaling function anyway --> no error
    para->setMaxLevel(2);
    para->initLBMSimulationParameter();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleRhoSq);
    EXPECT_NO_THROW(GridScalingKernelManager(para, &scalingFactory));
}

