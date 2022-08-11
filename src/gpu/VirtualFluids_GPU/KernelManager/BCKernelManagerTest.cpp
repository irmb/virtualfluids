#include <gmock/gmock.h>
#include <stdexcept>

#include "BCKernelManager.h"
#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Parameter/Parameter.h"
#include "PointerDefinitions.h"

class BCKernelManagerTest_BCsNotSpecified : public testing::Test
{
protected:
    BoundaryConditionFactory bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
    }
};

TEST_F(BCKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPre_NotSpecified)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BCKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BCKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPost_NotSpecified)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 1;
    EXPECT_THROW(BCKernelManager(para, &bcFactory), std::runtime_error);
}
