#include <gmock/gmock.h>

#include "Cylinder.h"

// CylinderTestAxisNormalToX ////////////////////////////////////////////////

class CylinderTestAxisNormalToX : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder(center[0], center[1], center[2], radius, height, Cylinder::PrincipalAxis::x);
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

// CylinderTestAxisNormalToY ////////////////////////////////////////////////

class CylinderTestAxisNormalToY : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Cylinder::PrincipalAxis::y);
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

// CylinderTestAxisNormalToZ ////////////////////////////////////////////////

class CylinderTestAxisNormalToZ : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Cylinder::PrincipalAxis::z);
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