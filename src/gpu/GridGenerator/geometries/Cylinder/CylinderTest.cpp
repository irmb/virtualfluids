#include <gmock/gmock.h>

#include "Cylinder.h"

// CylinderTestAxisNormalToX ////////////////////////////////////////////////

class CylinderTestAxisNormalToX : public testing::Test
{
protected:
    Cylinder cylinder = Cylinder({ 0.1, 0.2, 0.3 }, 2.0, 8.0, Cylinder::PrincipalAxis::x);
};

TEST_F(CylinderTestAxisNormalToX, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(0.1));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(0.2));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(0.3));
}

TEST_F(CylinderTestAxisNormalToX, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-4.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-2.0 + 0.3));
}

TEST_F(CylinderTestAxisNormalToX, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(4.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(2.0 + 0.3));
}

// CylinderTestAxisNormalToY ////////////////////////////////////////////////

class CylinderTestAxisNormalToY : public testing::Test
{
protected:
    Cylinder cylinder = Cylinder({ 0.1, 0.2, 0.3 }, 2.0, 8.0, Cylinder::PrincipalAxis::y);
};

TEST_F(CylinderTestAxisNormalToY, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(0.1));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(0.2));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(0.3));
}

TEST_F(CylinderTestAxisNormalToY, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-4.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-2.0 + 0.3));
}

TEST_F(CylinderTestAxisNormalToY, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(4.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(2.0 + 0.3));
}

// CylinderTestAxisNormalToZ ////////////////////////////////////////////////

class CylinderTestAxisNormalToZ : public testing::Test
{
protected:
    Cylinder cylinder = Cylinder({ 0.1, 0.2, 0.3 }, 2.0, 8.0, Cylinder::PrincipalAxis::z);
};

TEST_F(CylinderTestAxisNormalToZ, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(0.1));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(0.2));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(0.3));
}

TEST_F(CylinderTestAxisNormalToZ, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-4.0 + 0.3));
}

TEST_F(CylinderTestAxisNormalToZ, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(2.0 + 0.1));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(2.0 + 0.2));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(4.0 + 0.3));
}