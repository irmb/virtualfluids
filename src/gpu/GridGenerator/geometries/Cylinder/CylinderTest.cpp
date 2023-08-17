#include <gmock/gmock.h>

#include "Cylinder.h"


void testChangeSizeByDelta(Cylinder &cylinder)
{
    const double delta = 0.2;
    const double tolerance = 0.0000000001;
    const double oldRadius = cylinder.getRadius();
    const double oldHeight = cylinder.getHeight();
    const double oldMinimumX1 = cylinder.getX1Minimum();
    const double oldMinimumX2 = cylinder.getX2Minimum();
    const double oldMinimumX3 = cylinder.getX3Minimum();
    const double oldMaximumX1 = cylinder.getX1Maximum();
    const double oldMaximumX2 = cylinder.getX2Maximum();
    const double oldMaximumX3 = cylinder.getX3Maximum();

    cylinder.changeSizeByDelta(delta);

    EXPECT_THAT(cylinder.getRadius(), testing::Eq(oldRadius + delta));
    EXPECT_THAT(cylinder.getHeight(), testing::Eq(oldHeight + 2 * delta));
    EXPECT_THAT(cylinder.getX1Minimum(), testing::DoubleNear(oldMinimumX1 - delta, tolerance));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::DoubleNear(oldMinimumX2 - delta, tolerance));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::DoubleNear(oldMinimumX3 - delta, tolerance));
    EXPECT_THAT(cylinder.getX1Maximum(), testing::DoubleNear(oldMaximumX1 + delta, tolerance));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::DoubleNear(oldMaximumX2 + delta, tolerance));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::DoubleNear(oldMaximumX3 + delta, tolerance));
}

// CylinderTestAxisNormalToX ////////////////////////////////////////////////

class CylinderTestAxisNormalToX : public testing::Test
{
public:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder(center[0], center[1], center[2], radius, height, Cylinder::RotationalAxis::x);
};

TEST_F(CylinderTestAxisNormalToX, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToX, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-0.5 * height + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToX, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(0.5 * height + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToX, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - 0.5 * height + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - 0.5 * height - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + 0.5 * height - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + 0.5 * height + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - radius + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + radius + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - radius + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - radius - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + radius - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + radius + epsilon, 0.0, 0.0));
}

TEST_F(CylinderTestAxisNormalToX, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}

// CylinderTestAxisNormalToY ////////////////////////////////////////////////

class CylinderTestAxisNormalToY : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Cylinder::RotationalAxis::y);
};

TEST_F(CylinderTestAxisNormalToY, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToY, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-0.5 * height + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToY, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(0.5 * height + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToY, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - radius + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + radius + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - 0.5 * height + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - 0.5 * height - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + 0.5 * height - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + 0.5 * height + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - radius + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - radius - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + radius - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + radius + epsilon, 0.0, 0.0));
}

TEST_F(CylinderTestAxisNormalToY, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}

// CylinderTestAxisNormalToZ ////////////////////////////////////////////////

class CylinderTestAxisNormalToZ : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Cylinder::RotationalAxis::z);
};

TEST_F(CylinderTestAxisNormalToZ, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-0.5 * height + center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(0.5 * height + center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, isPointInObject)
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

TEST_F(CylinderTestAxisNormalToZ, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}
