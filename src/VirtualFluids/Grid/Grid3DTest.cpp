#include "gmock/gmock.h"

#include "Block3D.h"
#include "Grid3D.h"
#include "UbTuple.h"


TEST(BlockTest, transBlockToWorldCoordinates)
{
    const int worldPositionX1 = 0;
    const int worldPositionX2 = 0;
    const int worldPositionX3 = 0;
    const int level = 0;

    Block3DPtr block = Block3DPtr( new Block3D(worldPositionX1, worldPositionX2, worldPositionX3, level));

    int blockCoordX1 = 4;
    int blockCoordX2 = 4;
    int blockCoordX3 = 4;

    Grid3D grid;
    grid.setDeltaX(1);

    Vector3D worldCoords = grid.getNodeCoordinates(block, blockCoordX1, blockCoordX2, blockCoordX3);

    EXPECT_THAT(worldCoords[0], testing::DoubleEq(blockCoordX1 - OFFSET));
    EXPECT_THAT(worldCoords[1], testing::DoubleEq(blockCoordX1 - OFFSET));
    EXPECT_THAT(worldCoords[2], testing::DoubleEq(blockCoordX3 - OFFSET));
}

