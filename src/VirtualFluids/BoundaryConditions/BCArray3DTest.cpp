#include "gmock/gmock.h"

#include "BCArray3D.h"


TEST(BCArray3DTest, checkIfCoordinatesAreInsideOfTheBlock)
{
    const int nx = 10, ny = 10, nz = 10;
    BCArray3D sut(nx, ny, nz);

    int ghostLayerWidth = 1;
    EXPECT_FALSE(sut.isInsideOfDomain(0, 5, 5, ghostLayerWidth));
    EXPECT_FALSE(sut.isInsideOfDomain(5, 0, 5, ghostLayerWidth));
    EXPECT_FALSE(sut.isInsideOfDomain(0, 5, 0, ghostLayerWidth));

    EXPECT_TRUE(sut.isInsideOfDomain(1, 1, 1, ghostLayerWidth));
    EXPECT_TRUE(sut.isInsideOfDomain(8, 8, 8, ghostLayerWidth));

    EXPECT_FALSE(sut.isInsideOfDomain(9, 5, 1, ghostLayerWidth));

    ghostLayerWidth = 2;
    EXPECT_FALSE(sut.isInsideOfDomain(8, 5, 2, ghostLayerWidth));
    EXPECT_TRUE(sut.isInsideOfDomain(7, 2, 2, ghostLayerWidth));
}
