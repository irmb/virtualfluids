#include <gmock/gmock.h>

#include "MacroscopicQuantities.h"
#include "constants/D3Q27.h"


/*
* given distributions, which are all one.
*/
real f[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

using namespace vf::lbm;


TEST(MacroscopicQuantitiesTest, check_density)
{
    const double density = getDensity(f);

    const double expected_density = 27.;
    ASSERT_THAT(density, testing::DoubleEq(expected_density));
}

TEST(MacroscopicQuantitiesTest, whenFsAreEqual_velocityInEachDirectionShouldBeZero)
{
    const double velocityX1 = getIncompressibleVelocityX1(f);
    const double velocityX2 = getIncompressibleVelocityX2(f);
    const double velocityX3 = getIncompressibleVelocityX3(f);

    const double expected_velocity = 0.;
    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity));
}

TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Eis2_velocityInX1ShouldBeOne)
{
    f[dir::DIR_P00] = 2.;

    const double velocityX1 = getIncompressibleVelocityX1(f);
    const double velocityX2 = getIncompressibleVelocityX2(f);
    const double velocityX3 = getIncompressibleVelocityX3(f);

    const double expected_velocity_x1 = 1.;
    const double expected_velocity_x2 = 0.;
    const double expected_velocity_x3 = 0.;

    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity_x1));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity_x2));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity_x3));
}

TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Nis2_velocityInX2ShouldBeOne)
{
    f[dir::DIR_0P0] = 2.;

    const double velocity = getIncompressibleVelocityX2(f);

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}


TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Tis2_velocityInX3ShouldBeOne)
{
    f[dir::DIR_00P] = 2.;

    const double velocity = getIncompressibleVelocityX3(f);

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}
