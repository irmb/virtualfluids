#include <gmock/gmock.h>

#include "VerticalCylinder.h"

class VerticalCylinderTest : public testing::Test
{
protected:
    VerticalCylinder cylinder = VerticalCylinder(0.1, 0.2, 0.3, 2.0, 8.0);
};

TEST_F(VerticalCylinderTest, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(0.1));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(0.2));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(0.3));
}

TEST_F(VerticalCylinderTest, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-4.0 + 0.3));
}

TEST_F(VerticalCylinderTest, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(4.0 + 0.3));
}