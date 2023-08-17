#include <gmock/gmock.h>

#include "VerticalCylinder.h"

class VerticalCylinderTest : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    VerticalCylinder cylinder = VerticalCylinder(center[0], center[1], center[2], radius, height);
};

TEST_F(VerticalCylinderTest, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(VerticalCylinderTest, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-0.5 * height + center[2]));
}

TEST_F(VerticalCylinderTest, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(0.5 * height + center[2]));
}

TEST_F(VerticalCylinderTest, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - radius + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + radius + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - radius + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + radius + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - 0.5 * height + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - 0.5 * height - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + 0.5 * height - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + 0.5 * height + epsilon, 0.0, 0.0));
}