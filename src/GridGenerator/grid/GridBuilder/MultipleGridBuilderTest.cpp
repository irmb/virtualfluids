#include "gmock/gmock.h"

#include "MultipleGridBuilder.h"
#include "../GridMocks.h"

TEST(MultipleGridBuilderTest, addOneGrid_numberOfLevelsShouldBeOne)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, true, true, true);

    EXPECT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1)); 
}