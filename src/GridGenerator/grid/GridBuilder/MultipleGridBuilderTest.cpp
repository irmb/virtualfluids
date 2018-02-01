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

TEST(MultipleGridBuilderTest, addMultipleGrids_deltaShouldBeTheHalfOfTheLastAdded)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    const real delta = 2.0;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta);

    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);

    EXPECT_THAT(gridBuilder->getDelta(1), RealEq(delta / 2.0));
    EXPECT_THAT(gridBuilder->getDelta(2), RealEq(delta / 4.0));
    EXPECT_THAT(gridBuilder->getDelta(3), RealEq(delta / 8.0));
    EXPECT_THAT(gridBuilder->getDelta(4), RealEq(delta / 16.0));
}

TEST(MultipleGridBuilderTest, getInvalidLevel_shouldThrowException)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    ASSERT_THROW(gridBuilder->getDelta(0), InvalidLevelException);
}

TEST(MultipleGridBuilderTest, addedsecondGrid_shouldBeStaggered)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    const real givenStartX = 0.0;
    const real givenStartY = 0.0;
    const real givenStartZ = 0.0;

    const real givenEndX = 10.0;
    const real givenEndY = 10.0;
    const real givenEndZ = 10.0;
    gridBuilder->addGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ);

    const uint level = 1;
    const real staggeredOffset = 0.5 * gridBuilder->getDelta(level);

    const real expectedStartX = givenStartX + staggeredOffset;
    const real expectedStartY = givenStartY + staggeredOffset;
    const real expectedStartZ = givenStartZ + staggeredOffset;

    const real actualStartX = gridBuilder->getStartX(level);
    const real actualStartY = gridBuilder->getStartY(level);
    const real actualStartZ = gridBuilder->getStartZ(level);

    EXPECT_THAT(actualStartX, RealEq(expectedStartX));
    EXPECT_THAT(actualStartY, RealEq(expectedStartY));
    EXPECT_THAT(actualStartZ, RealEq(expectedStartZ));


    const real expectedEndX = givenEndX - staggeredOffset;
    const real expectedEndY = givenEndY - staggeredOffset;
    const real expectedEndZ = givenEndZ - staggeredOffset;


    const real actualEndX = gridBuilder->getEndX(level);
    const real actualEndY = gridBuilder->getEndY(level);
    const real actualEndZ = gridBuilder->getEndZ(level);

    EXPECT_THAT(actualEndX, RealEq(expectedEndX));
    EXPECT_THAT(actualEndY, RealEq(expectedEndY));
    EXPECT_THAT(actualEndZ, RealEq(expectedEndZ));
}

//TEST(MultipleGridBuilderTest, addsFineGridWithLevel_shouldCreateGridsBetween)
//{
//    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
//
//    real delta = 2.0;
//    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta);
//    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
//
//    uint level = 5;
//    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);
//
//    ASSERT_THROW(gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0), FinerGridBiggerThanCoarsestGridException);
//}
