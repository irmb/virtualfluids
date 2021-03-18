#include <gmock/gmock.h>

#include "MacroscopicQuantities.h"
#include "D3Q27.h"


/*
* given distributions, which are all one.
*/
real f[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};



TEST(MacroscopicQuantitiesTest, check_density)
{
    const double density = VF::LBM::getDensity(f);

    const double expected_density = 27.;
    ASSERT_THAT(density, testing::DoubleEq(expected_density));
}

TEST(MacroscopicQuantitiesTest, whenFsAreEqual_velocityInEachDirectionShouldBeZero)
{
    const double velocityX1 = VF::LBM::getIncompressibleVelocityX1(f);
    const double velocityX2 = VF::LBM::getIncompressibleVelocityX2(f);
    const double velocityX3 = VF::LBM::getIncompressibleVelocityX3(f);

    const double expected_velocity = 0.;
    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity));
}

TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Eis2_velocityInXShouldBeOne)
{
    f[VF::LBM::DIR::E] = 2.;

    const double velocityX1 = VF::LBM::getIncompressibleVelocityX1(f);
    const double velocityX2 = VF::LBM::getIncompressibleVelocityX2(f);
    const double velocityX3 = VF::LBM::getIncompressibleVelocityX3(f);

    const double expected_velocity_x1 = 1.;
    const double expected_velocity_x2 = 0.;
    const double expected_velocity_x3 = 0.;

    EXPECT_THAT(velocityX1, testing::DoubleEq(expected_velocity_x1));
    EXPECT_THAT(velocityX2, testing::DoubleEq(expected_velocity_x2));
    EXPECT_THAT(velocityX3, testing::DoubleEq(expected_velocity_x3));
}

TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Nis2_velocityInX2ShouldBeOne)
{
    f[VF::LBM::DIR::N] = 2.;

    const double velocity = VF::LBM::getIncompressibleVelocityX2(f);

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}


TEST(MacroscopicQuantitiesTest, givenAllFsAreOne_when_Tis2_velocityInX3ShouldBeOne)
{
    f[VF::LBM::DIR::T] = 2.;

    const double velocity = VF::LBM::getIncompressibleVelocityX3(f);

    const double expected_velocity = 1.;
    ASSERT_THAT(velocity, testing::DoubleEq(expected_velocity));
}
