#include "gmock/gmock.h"

#include "MultipleGridBuilder.h"
#include "../GridMocks.h"

TEST(MultipleGridBuilderTest, addOneGrid_numberOfLevelsShouldBeOne)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1));
}

TEST(MultipleGridBuilderTest, addTwoGrids_secondGridMustBeInsideOfFirstGrid)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    ASSERT_THROW(gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0), FinerGridBiggerThanCoarsestGridException);
}

TEST(MultipleGridBuilderTest, givenCoarseGrid_addAdditionalGrid_shouldCalculateDelta)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    const real delta = 2.0;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta);

    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);

    ASSERT_THAT(gridBuilder->getDelta(1), RealEq(delta / 2.0));
}

TEST(MultipleGridBuilderTest, firstGridMustBeCoarse)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    ASSERT_THROW(gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0), FirstGridMustBeCoarseException);
}


//TEST(MultipleGridBuilderTest, fall2)
//{
//    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
//
//    real delta = 2.0;
//    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta, true, true, true);
//    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, true, true, true);
//
//    uint level = 5;
//    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level, true, true, true);
//
//    ASSERT_THROW(gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, true, true, true), FinerGridBiggerThanCoarsestGridException);
//}
