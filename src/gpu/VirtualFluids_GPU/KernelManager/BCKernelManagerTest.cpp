#include <gmock/gmock-function-mocker.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <stdexcept>

#include "BCKernelManager.h"
#include "Factories/BoundaryConditionFactory.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"

class BCKernelManagerTest_BCsNotSpecified : public testing::Test
{
protected:
    BoundaryConditionFactory bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
        para->initLBMSimulationParameter();
    }
};

TEST_F(BCKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory)); // no throw, as a default is specified
    para->getParD(0)->noSlipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
    para->getParD(0)->slipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPre_NotSpecified)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
    para->getParD(0)->pressureBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory)); // no throw, as a default is specified
    para->getParD(0)->geometryBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
    para->getParD(0)->stressBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory));
}

class BoundaryConditionFactoryMock : public BoundaryConditionFactory
{
public:
    mutable uint numberOfCalls = 0;

    [[nodiscard]] boundaryCondition getVelocityBoundaryConditionPost(bool) const override
    {
        return [this](LBMSimulationParameter *, QforBoundaryConditions *) { numberOfCalls++; };
    }
};

class BCKernelManagerTest_runBCs : public testing::Test
{
protected:
    BoundaryConditionFactoryMock bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();
    UPtr<BCKernelManager> sut;

    void SetUp() override
    {
        para->initLBMSimulationParameter();
        sut = std::make_unique<BCKernelManager>(para, &bcFactory);
    }
};

TEST_F(BCKernelManagerTest_runBCs, runVelocityBCKernelPost)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(1));

    bcFactory.numberOfCalls = 0;
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(0));
}