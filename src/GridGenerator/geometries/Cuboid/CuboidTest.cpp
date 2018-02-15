#include "gmock/gmock.h"


#include "Cuboid.h"


using namespace testing;

TEST(CuboidTest, createCuboid_shouldReturnCenter)
{
    Cuboid sut(1, 3, 5, 10, 11, 12);

    const double centerX1 = sut.getX1Centroid();
    const double centerX2 = sut.getX2Centroid();
    const double centerX3 = sut.getX3Centroid();

    EXPECT_THAT(centerX1, DoubleEq(5.5));
    EXPECT_THAT(centerX2, DoubleEq(7.0));
    EXPECT_THAT(centerX3, DoubleEq(8.5));
}

TEST(CuboidTest, createCuboid_shouldReturnMinimum)
{
    Cuboid sut(1, 3, 5, 10, 11, 3);

    const double minX1 = sut.getX1Minimum();
    const double minX2 = sut.getX2Minimum();
    const double minX3 = sut.getX3Minimum();

    EXPECT_THAT(minX1, DoubleEq(1));
    EXPECT_THAT(minX2, DoubleEq(3));
    EXPECT_THAT(minX3, DoubleEq(3));
}

TEST(CuboidTest, createCuboid_shouldReturnMaximum)
{
    Cuboid sut(1, 3, 5, 10, 11, 3);

    const double maxX1 = sut.getX1Maximum();
    const double maxX2 = sut.getX2Maximum();
    const double maxX3 = sut.getX3Maximum();

    EXPECT_THAT(maxX1, DoubleEq(10));
    EXPECT_THAT(maxX2, DoubleEq(11));
    EXPECT_THAT(maxX3, DoubleEq(5));
}


TEST(CuboidTest, checkIfPointIsInsideOfCuboid)
{
    Cuboid sut(1, 3, 5, 10, 11, 12);

    const double minOffset = 0;
    const double maxOffset = 0;

    EXPECT_TRUE(sut.isPointInObject(2, 4, 6, minOffset, maxOffset));
    EXPECT_FALSE(sut.isPointInObject(1, 4, 6, minOffset, maxOffset));
    EXPECT_FALSE(sut.isPointInObject(2, 11, 6, minOffset, maxOffset));
    EXPECT_FALSE(sut.isPointInObject(2, 4, 12, minOffset, maxOffset));
}

TEST(CuboidTest, checkIfPointIsInsideOfCuboid_withOffset)
{
    Cuboid sut(1, 3, 5, 10, 11, 12);

    const double minOffset = 2;
    const double maxOffset = 1;

    EXPECT_FALSE(sut.isPointInObject(3, 5, 7, minOffset, maxOffset));
    EXPECT_FALSE(sut.isPointInObject(9, 10, 11, minOffset, maxOffset));
}
