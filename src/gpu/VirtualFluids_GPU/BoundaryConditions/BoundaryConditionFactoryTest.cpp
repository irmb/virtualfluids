#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "BoundaryConditionFactory.h"

TEST(BoundaryConditionFactoryTest, defaultVelocityBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getVelocityBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultNoSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getSlipBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultPressureBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getPressureBoundaryConditionPre();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultGeometryBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultStressBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getStressBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr, 0), std::bad_function_call);
}