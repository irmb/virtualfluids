#include "gmock/gmock.h"

#include "MultipleGridBuilder.h"
#include "../GridMocks.h"
#include "../GridStrategy/GridStrategyMocks.h"
#include "../GridFactory.h"

class MultipleGridBuilderAddGridTest : public testing::Test
{
    virtual void SetUp() override
    {
        gridFactory = SPtr<GridFactory>(new GridFactory());
        gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
        gridFactory->setGrid("stub");
        gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    }
public:
    SPtr<MultipleGridBuilder> gridBuilder;

private:
    SPtr<GridFactory> gridFactory;
};

TEST_F(MultipleGridBuilderAddGridTest, addOneGrid_numberOfLevelsShouldBeOne)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1));
}

TEST_F(MultipleGridBuilderAddGridTest, addTwoGridsWhereSecondGridIsBigger_GridShouldNotAdded)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);
    gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0);
    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1));
}

TEST_F(MultipleGridBuilderAddGridTest, givenCoarseGrid_addAdditionalGrid_shouldCalculateDelta)
{
    const real delta = 2.0;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta);

    gridBuilder->addGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);

    ASSERT_THAT(gridBuilder->getDelta(1), RealEq(delta * 0.5));
}

TEST_F(MultipleGridBuilderAddGridTest, addGridWithoutCoarseGrid_shouldNotbeAdded)
{
    gridBuilder->addGrid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(0));
}

TEST_F(MultipleGridBuilderAddGridTest, addMultipleGrids_deltaShouldBeTheHalfOfTheLastAdded)
{
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

TEST_F(MultipleGridBuilderAddGridTest, getInvalidLevel_shouldThrowException)
{
    ASSERT_THROW(gridBuilder->getDelta(0), std::exception);
}

TEST_F(MultipleGridBuilderAddGridTest, addFineGridWithoutCoarseGrid_ShouldNotAddingAGrid)
{
    gridBuilder->addFineGrid(0, 0, 0, 0, 0, 0, 0);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(0));
}


TEST_F(MultipleGridBuilderAddGridTest, addGridWithFloatingStartPoints_ShouldCreatedStaggeredCoordinatesToCoarseGrid)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1);

    gridBuilder->addGrid(0.1212, 0.1212, 0.221, 10.867, 10.45454, 10.12121);


    EXPECT_THAT(gridBuilder->getStartX(1), RealEq(0.25));
    EXPECT_THAT(gridBuilder->getStartY(1), RealEq(0.25));
    EXPECT_THAT(gridBuilder->getStartZ(1), RealEq(0.25));

    EXPECT_THAT(gridBuilder->getEndX(1), RealEq(9.75));
    EXPECT_THAT(gridBuilder->getEndY(1), RealEq(9.75));
    EXPECT_THAT(gridBuilder->getEndZ(1), RealEq(9.75));
}

TEST_F(MultipleGridBuilderAddGridTest, addGridWitNegativStartPoints_ShouldCreatedStaggeredCoordinatesToCoarseGrid)
{
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 15.0, 15.0, 15.0, 1);

    gridBuilder->addGrid(-20.0, -20.0, -20.0, -5.0, -5.0, -5.0);

    gridBuilder->addGrid(-15.0, -15.0, -15.0, -10.0, -10.0, -10.0);

    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(-15.625));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(-15.625));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(-15.625));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(-9.375));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(-9.375));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(-9.375));
}


void expectStartCoordinatesAreStaggered(const real givenStartX, const real givenStartY, const real givenStartZ, const real staggeredOffset, SPtr<MultipleGridBuilder> gridBuilder, const uint level)
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

void expectEndCoordinatesAreStaggered(const real givenEndX, const real givenEndY, const real givenEndZ, const real staggeredOffset, SPtr<MultipleGridBuilder> gridBuilder, const uint level)
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

TEST_F(MultipleGridBuilderAddGridTest, addedsecondGrid_shouldBeStaggered)
{
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

TEST_F(MultipleGridBuilderAddGridTest, addGridAfterCoarseGridWithFloatingStartPoint)
{
    gridBuilder->addCoarseGrid(1.2, 1.2, 1.2, 15.2, 15.2, 15.2, 1.0);

    const real givenStartX = 2.0;
    const real givenStartY = 3.0;
    const real givenStartZ = 4.0;
    const real givenEndX = 10.0;
    const real givenEndY = 11.0;
    const real givenEndZ = 12.0;
    gridBuilder->addGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ);

    const uint level = 1;

    //0.25 offset to start point
    EXPECT_THAT(gridBuilder->getStartX(level), RealEq(2.45));
    EXPECT_THAT(gridBuilder->getStartY(level), RealEq(3.45));
    EXPECT_THAT(gridBuilder->getStartZ(level), RealEq(4.45));

    //-0.25 offset to start point
    EXPECT_THAT(gridBuilder->getEndX(level), RealEq(9.55));
    EXPECT_THAT(gridBuilder->getEndY(level), RealEq(10.55));
    EXPECT_THAT(gridBuilder->getEndZ(level), RealEq(11.55));
}

TEST_F(MultipleGridBuilderAddGridTest, addedthirdGrid_shouldBeStaggered)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);

    const real givenStartX = 0.0;
    const real givenStartY = 1.0;
    const real givenStartZ = 2.0;
    const real givenEndX = 10.0;
    const real givenEndY = 11.0;
    const real givenEndZ = 12.0;

    gridBuilder->addGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ);
    gridBuilder->addGrid(3.0, 4.0, 5.0, 5.0, 6.0, 7.0);

    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(3.375));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(4.375));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(5.375));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(4.625));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(5.625));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(6.625));
}



TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevel_shouldCreateGridsBetween)
{
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, 1.0);

    const uint level = 5;
    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(level + 1));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelThree_shouldCalculateDeltaForLevels)
{
    const real startDelta = 1.0;
    gridBuilder->addCoarseGrid(-100.0, -100.0, -100.0, 100.0, 100.0, 100.0, 1.0);

    uint level = 3;
    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, level);

    const real expectedDeltaLevel3 = startDelta / std::pow(2, level);
    EXPECT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel3));

    level--;
    const real expectedDeltaLevel2 = startDelta / std::pow(2, level);
    EXPECT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel2));

    level--;
    const real expectedDeltaLevel1 = startDelta / std::pow(2, level);
    EXPECT_THAT(gridBuilder->getDelta(level), RealEq(expectedDeltaLevel1));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelThree_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, 1.0);

    uint level = 3;
    gridBuilder->addFineGrid(20.0, 20.0, 20.0, 40.0, 40.0, 40.0, level);

    EXPECT_THAT(gridBuilder->getStartX(level), RealEq(20.4375));
    EXPECT_THAT(gridBuilder->getStartY(level), RealEq(20.4375));
    EXPECT_THAT(gridBuilder->getStartZ(level), RealEq(20.4375));

    EXPECT_THAT(gridBuilder->getEndX(level), RealEq(39.5625));
    EXPECT_THAT(gridBuilder->getEndY(level), RealEq(39.5625));
    EXPECT_THAT(gridBuilder->getEndZ(level), RealEq(39.5625));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelThree_shouldCreateStaggeredStartAndEndPointIntermediateGrid)
{
    const real startDelta = 1.0;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, startDelta);

    const uint level = 3;
    const real givenStartX = 20.0;
    const real givenStartY = 21.0;
    const real givenStartZ = 22.0;
    const real givenEndX = 40.0;
    const real givenEndY = 41.0;
    const real givenEndZ = 42.0;
    gridBuilder->addFineGrid(givenStartX, givenStartY, givenStartZ, givenEndX, givenEndY, givenEndZ, level);


    const real expectedStartXLevel2 = givenStartX + 0.375 - 8.0 * gridBuilder->getDelta(2);
    const real expectedStartYLevel2 = givenStartY + 0.375 - 8.0 * gridBuilder->getDelta(2);
    const real expectedStartZLevel2 = givenStartZ + 0.375 - 8.0 * gridBuilder->getDelta(2);

    const real expectedEndXLevel2 = givenEndX - 0.375 + 8.0 * gridBuilder->getDelta(2);
    const real expectedEndYLevel2 = givenEndY - 0.375 + 8.0 * gridBuilder->getDelta(2);
    const real expectedEndZLevel2 = givenEndZ - 0.375 + 8.0 * gridBuilder->getDelta(2);
    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(expectedStartXLevel2));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(expectedStartYLevel2));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(expectedStartZLevel2));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(expectedEndXLevel2));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(expectedEndYLevel2));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(expectedEndZLevel2));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelTwoWithCoarseGridSize_gridsShouldNotBeAdded)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    gridBuilder->addFineGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 2);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelOne_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 1.0);

    const uint level = 1;
    gridBuilder->addFineGrid(5.0, 5.0, 5.0, 7.0, 7.0, 7.0, level);

    EXPECT_THAT(gridBuilder->getStartX(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartY(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartZ(1), RealEq(5.25));

    EXPECT_THAT(gridBuilder->getEndX(1), RealEq(6.75));
    EXPECT_THAT(gridBuilder->getEndY(1), RealEq(6.75));
    EXPECT_THAT(gridBuilder->getEndZ(1), RealEq(6.75));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelTwo_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    gridBuilder->addCoarseGrid(-5.0, -5.0, -5.0, 12.0, 12.0, 12.0, 1.0);

    const uint level = 2;
    gridBuilder->addFineGrid(5.0, 5.0, 5.0, 7.0, 7.0, 7.0, level);


    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(5.375));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(5.375));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(5.375));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(6.625));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(6.625));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(6.625));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelTwoAndTwoGridsBefore_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    gridBuilder->addCoarseGrid(1.0, 4.0, 2.0, 30.0, 30.0, 30.0, 1.0);
    gridBuilder->addGrid(5.0, 5.0, 5.0, 20.0, 20.0, 20.0);

    const uint level = 2;
    gridBuilder->addFineGrid(10.0, 10.0, 10.0, 12.0, 12.0, 12.0, level);

    EXPECT_THAT(gridBuilder->getStartX(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartY(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartZ(1), RealEq(5.25));

    EXPECT_THAT(gridBuilder->getEndX(1), RealEq(19.75));
    EXPECT_THAT(gridBuilder->getEndY(1), RealEq(19.75));
    EXPECT_THAT(gridBuilder->getEndZ(1), RealEq(19.75));

    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(10.375));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(10.375));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(10.375));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(11.625));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(11.625));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(11.625));
}

TEST_F(MultipleGridBuilderAddGridTest, addsFineGridWithLevelTwoAndTwoGridsBeforeAndFloatingStartingPoints_shouldCreateStaggeredStartAndEndPointFineGrid)
{
    gridBuilder->addCoarseGrid(1.2, 4.0, 2.0, 30.0, 30.0, 30.0, 1.0);
    gridBuilder->addGrid(5.0, 5.0, 5.0, 20.0, 20.0, 20.0);

    const uint level = 2;
    gridBuilder->addFineGrid(10.0, 10.0, 10.0, 12.0, 12.0, 12.0, level);


    EXPECT_THAT(gridBuilder->getStartX(2), RealEq(10.575));
    EXPECT_THAT(gridBuilder->getStartY(2), RealEq(10.375));
    EXPECT_THAT(gridBuilder->getStartZ(2), RealEq(10.375));

    EXPECT_THAT(gridBuilder->getEndX(2), RealEq(11.425));
    EXPECT_THAT(gridBuilder->getEndY(2), RealEq(11.625));
    EXPECT_THAT(gridBuilder->getEndZ(2), RealEq(11.625));
}

TEST(MultipleGridBuilderTest, everyExceptTheFinestGrid_shouldHaveAGridInterface)
{
    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
    gridFactory->setGrid("spy");
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);
    gridBuilder->addFineGrid(7.0, 7.0, 7.0, 10.0, 10.0, 10.0, 2);

    gridBuilder->buildGrids();

    auto grid0 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(0));
    auto grid1 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(1));
    EXPECT_TRUE(grid0->hasGridInterface());
    EXPECT_TRUE(grid1->hasGridInterface());
}


TEST(MultipleGridBuilderTest, afterCreatingGridInterface_FineGridsShouldNotBeHavePeriodicBoundarys)
{
    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
    gridFactory->setGrid("spy");
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);
    gridBuilder->addFineGrid(7.0, 7.0, 7.0, 10.0, 10.0, 10.0, 2);

    gridBuilder->buildGrids();

    auto grid1 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(1));
    auto grid2 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(2));

    EXPECT_TRUE(grid1->hasNoPeriodicityBoundaries());
    EXPECT_TRUE(grid2->hasNoPeriodicityBoundaries());
}

TEST(MultipleGridBuilderDifferentShapesTest, addSphereGridShouldSetMinimaMaximaToGrid)
{
    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
    gridFactory->setGrid("spy");
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);

    gridBuilder->addGrid(new Sphere(7.5, 7.5, 7.5, 2.5));

    EXPECT_THAT(gridBuilder->getStartX(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartY(1), RealEq(5.25));
    EXPECT_THAT(gridBuilder->getStartZ(1), RealEq(5.25));

    EXPECT_THAT(gridBuilder->getEndX(1), RealEq(9.75));
    EXPECT_THAT(gridBuilder->getEndY(1), RealEq(9.75));
    EXPECT_THAT(gridBuilder->getEndZ(1), RealEq(9.75));
}
