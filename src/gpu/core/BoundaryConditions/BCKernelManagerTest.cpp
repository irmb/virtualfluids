
#include <gmock/gmock.h>
#include <stdexcept>

#include "BoundaryConditionKernelManager.h"
#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Parameter/Parameter.h"
#include <basics/PointerDefinitions.h>

class BoundaryConditionKernelManagerTest_BCsNotSpecified : public testing::Test
{
protected:
    BoundaryConditionFactory bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
        para->initLBMSimulationParameter();
    }
};

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPostNotSpecified_withBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPreNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPreNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPostNotSpecified_withBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, precursorBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->precursorBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, precursorBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->precursorBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
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

class BoundaryConditionKernelManagerTest_runBCs : public testing::Test
{
protected:
    BoundaryConditionFactoryMock bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();
    UPtr<BoundaryConditionKernelManager> sut;

    void SetUp() override
    {
        para->initLBMSimulationParameter();
        sut = std::make_unique<BoundaryConditionKernelManager>(para, &bcFactory);
    }
};

TEST_F(BoundaryConditionKernelManagerTest_runBCs, runVelocityBCKernelPost)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(1));

    bcFactory.numberOfCalls = 0;
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(0));
}
