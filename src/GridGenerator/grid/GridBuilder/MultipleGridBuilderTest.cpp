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


void expectStartCoordinatesAreStaggered(const real givenStartX, const real givenStartY, const real givenStartZ, const real staggeredOffset, SPtr<MultipleGridBuilder<GridDummy> > gridBuilder, const uint level)
{
    const real expectedStartX = givenStartX + staggeredOffset;
    const real expectedStartY = givenStartY + staggeredOffset;
    const real expectedStartZ = givenStartZ + staggeredOffset;

    const real actualStartX = gridBuilder->getStartX(level);
    const real actualStartY = gridBuilder->getStartY(level);
    const real actualStartZ = gridBuilder->getStartZ(level);

    EXPECT_THAT(actualStartX, RealEq(expectedStartX));
    EXPECT_THAT(actualStartY, RealEq(expectedStartY));
    EXPECT_THAT(actualStartZ, RealEq(expectedStartZ));
}

void expectEndCoordinatesAreStaggered(const real givenEndX, const real givenEndY, const real givenEndZ, const real staggeredOffset, SPtr<MultipleGridBuilder<GridDummy> > gridBuilder, const uint level)
{
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

TEST(MultipleGridBuilderTest, addedsecondGrid_shouldBeStaggered)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);

    const real givenStartX = 0.0;
    const real givenStartY = 1.0;
    const real givenStartZ = 2.0;
    const real givenEndX = 10.0;
    const real givenEndY = 11.0;
    const real givenEndZ = 12.0;
    gridBuilder->addGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ);

    const uint level = 1;
    const real staggeredOffset = 0.5 * gridBuilder->getDelta(level);
    expectStartCoordinatesAreStaggered(givenStartX,  givenStartY, givenStartZ, staggeredOffset, gridBuilder, level);
    expectEndCoordinatesAreStaggered(givenEndX, givenEndY, givenEndZ, staggeredOffset, gridBuilder, level);
}

TEST(MultipleGridBuilderTest, addsFineGridWithLevel_shouldCreateGridsBetween)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, 1.0);

    const uint level = 5;
    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(level + 1));
}

TEST(MultipleGridBuilderTest, addsFineGridWithLevelThree_shouldCalculateDeltaForLevels)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    const real startDelta = 1.0;
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, 1.0);

    uint level = 3;
    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);

    const real expectedDeltaLevel3 = startDelta / std::pow(2, level);
    ASSERT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel3));

    level--;
    const real expectedDeltaLevel2 = startDelta / std::pow(2, level);
    ASSERT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel2));

    level--;
    const real expectedDeltaLevel1 = startDelta / std::pow(2, level);
    ASSERT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel1));
}

TEST(MultipleGridBuilderTest, addsFineGridWithLevelThree_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    const real startDelta = 1.0;
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, 1.0);

    uint level = 3;
    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);

    EXPECT_THAT(gridBuilder->getStartX(level), RealEq(0.4375));
    EXPECT_THAT(gridBuilder->getStartY(level), RealEq(0.4375));
    EXPECT_THAT(gridBuilder->getStartZ(level), RealEq(0.4375));

    EXPECT_THAT(gridBuilder->getEndX(level), RealEq(9.5625));
    EXPECT_THAT(gridBuilder->getEndY(level), RealEq(9.5625));
    EXPECT_THAT(gridBuilder->getEndZ(level), RealEq(9.5625));
}

TEST(MultipleGridBuilderTest, addsFineGridWithLevelThree_shouldCreateStaggeredStartAndEndPointIntermediateGrid)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    const real startDelta = 1.0;
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, startDelta);

    uint level = 3;
    const real givenStartX = 4.0;
    const real givenStartY = 5.0;
    const real givenStartZ = 6.0;
    const real givenEndX = 10.0;
    const real givenEndY = 11.0;
    const real givenEndZ = 12.0;
    gridBuilder->addFineGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ, level);


    const real expectedStartXLevel2 = givenStartX + 0.375 - 7.0 * gridBuilder->getDelta(2);
    const real expectedStartYLevel2 = givenStartY + 0.375 - 7.0 * gridBuilder->getDelta(2);
    const real expectedStartZLevel2 = givenStartZ + 0.375 - 7.0 * gridBuilder->getDelta(2);

    const real expectedEndXLevel2 = givenEndX - 0.375 + 7.0 * gridBuilder->getDelta(2);
    const real expectedEndYLevel2 = givenEndY - 0.375 + 7.0 * gridBuilder->getDelta(2);
    const real expectedEndZLevel2 = givenEndZ - 0.375 + 7.0 * gridBuilder->getDelta(2);
    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(expectedStartXLevel2));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(expectedStartYLevel2));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(expectedStartZLevel2));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(expectedEndXLevel2));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(expectedEndYLevel2));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(expectedEndZLevel2));
}


TEST(MultipleGridBuilderTest, addsFineGridWithLevelTwoWithCoarseGridSize_shouldThrow)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    const uint level = 2;

    ASSERT_THROW(gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level), FinerGridBiggerThanCoarsestGridException);
}


TEST(MultipleGridBuilderTest, addGridAndAddFinestGrid)
{
    auto gridBuilder = MultipleGridBuilder<GridDummy>::makeShared();

    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 1.0);
    gridBuilder->addGrid(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

    gridBuilder->addFineGrid(20, 20, 20, 40, 40, 40, 3);

}
