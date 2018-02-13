#include "gmock/gmock.h"

//#include "LevelGridBuilder.h"
//#include "../Grid.cuh"

//TEST(LevelGridBuilderTest, addTwoGridsWithStartCoordinates_shouldMakeTheFistGridToLevel1)
//{
//    auto gridBuilder = LevelGridBuilder::makeShared("cpu", "D3Q27");
//    gridBuilder->addGrid(12.0, 12.0, 12.0, 14.0, 14.0, 14.0, false, false, false);
//    gridBuilder->addGrid(10.0, 10.0, 10.0, 20.0, 20.0, 20.0, true, true, true);
//
//    gridBuilder->generateGrids();
//
//    EXPECT_THAT(gridBuilder->getGrid(1)->startX, RealEq(12.25)); 
//}
