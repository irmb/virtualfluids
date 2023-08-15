#include "basics/tests/testUtilities.h"
#include "geometries/VerticalCylinder/VerticalCylinder.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"

class MultipleGridBuilderForTest: public MultipleGridBuilder{
public:
    using MultipleGridBuilder::makeGrid;
    using MultipleGridBuilder::makeRotatingGrid;
};

class MultipleGridBuilderTestFixture : public testing::Test
{
protected:
    MultipleGridBuilderForTest gridBuilder;
    SPtr<VerticalCylinder> cylinder;
    real delta = 0.1;

public:
    void SetUp() override
    {
        cylinder = std::make_shared<VerticalCylinder>(0.0, 0.0, 0.0, 2.0, 8.0);
        gridBuilder.addCoarseGrid(-10. + 0.5 * delta, -10. + 0.5 * delta, -10. + 0.5 * delta, 10. - 0.5 * delta,
                                  10. - 0.5 * delta, 10. - 0.5 * delta, delta);
    }
};

TEST_F(MultipleGridBuilderTestFixture, addCoarseGrid_addsOneGridToGridList)
{
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(1));
}

TEST(MultipleGridBuilderTest, noCoarseGrid_addGrid_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<Object> cylinder;

    testingVF::captureStdOut();
    gridBuilder.addGrid(cylinder, 1);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGrid_doesNotWarn)
{
    testingVF::captureStdOut();
    gridBuilder.addGrid(cylinder, 1);
    EXPECT_FALSE(testingVF::stdoutContainsWarning());
}

TEST(MultipleGridBuilderTest, noCoarseGrid_addGridWithoutLevel_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<Object> cylinder;

    testingVF::captureStdOut();
    gridBuilder.addGrid(cylinder);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGridWitoutLevel_doesNotWarn)
{
    testingVF::captureStdOut();
    gridBuilder.addGrid(cylinder);
    EXPECT_FALSE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGrid_addsGridToList)
{
    uint levelFine = 1;

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + 1 new level = 2 levels";
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGridAtLevelTwo_addsTwoGridsToList)
{
    uint levelFine = 2;

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3)) << "coarse level + 2 new levels = 3 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtLevelTwo_addsOneGridToList)
{
    uint levelFine = 2;

    gridBuilder.addGrid(cylinder, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3)) << "coarse level + fineLevel + 1 new levels = 3 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtLevelThree_addsTwoGridsToList)
{
    uint levelFine = 3;

    gridBuilder.addGrid(cylinder, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(4)) << "coarse level + fineLevel + 2 new levels = 4 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtSameLevel_noGridAdded)
{
    uint levelFine = 1;

    gridBuilder.addGrid(cylinder, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2))
        << "Can't add grid on same level twice --> no new grid should be added";
}

TEST_F(MultipleGridBuilderTestFixture, addGrid_hasHalfDeltaComparedToCoarseGrid)
{
    uint levelFine = 1;

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrid(levelFine)->getDelta(), testing::Eq(0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, addGridAtLevelTwo_hasHalfOfHalfDeltaComparedToCoarseGrid)
{
    uint levelFine = 2;

    gridBuilder.addGrid(cylinder, levelFine);
    EXPECT_THAT(gridBuilder.getGrid(levelFine)->getDelta(), testing::Eq(0.5 * 0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, addGridWithoutLevel_addsGridAtLevelOne)
{
    gridBuilder.addGrid(cylinder);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2));
    EXPECT_THAT(gridBuilder.getGrid(1)->getDelta(), testing::Eq(0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExists_addGridWithoutLevel_addsGridAtLevelTwo)
{
    gridBuilder.addGrid(cylinder);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2));
    EXPECT_THAT(gridBuilder.getGrid(1)->getDelta(), testing::Eq(0.5 * delta));
    gridBuilder.addGrid(cylinder);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3));
    EXPECT_THAT(gridBuilder.getGrid(2)->getDelta(), testing::Eq(0.5 * 0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, addRotatingGrid_hasCorrectDelta)
{
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(1));
    gridBuilder.addGridRotatingGrid(cylinder);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(1)) << "Rotating grid is not added to list of grids.";
    EXPECT_THAT(gridBuilder.getGrid(0)->getDelta(), RealEq(gridBuilder.getRotatingGrid()->getDelta()));

    gridBuilder.addGridRotatingGrid(cylinder);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(1)) << "Rotating grid is not added to list of grids.";
    EXPECT_THAT(gridBuilder.getGrid(0)->getDelta(), RealEq(gridBuilder.getRotatingGrid()->getDelta()));
}

TEST(MultipleGridBuilderTest, noCoarseGrid_addGridWithPredefinedDelta_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<VerticalCylinder> cylinder;

    testingVF::captureStdOut();
    gridBuilder.addGridRotatingGrid(cylinder);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, makeGrid_innerGridHasCorrectDimensions)
{
    const auto numberOfLevels = 1;
    const auto grid = gridBuilder.makeGrid(cylinder, numberOfLevels, 0);

    EXPECT_THAT(grid->getStartX(), RealNear(cylinder->getX1Minimum() - 0.25 * delta, 0.00005));
    EXPECT_THAT(grid->getStartY(), RealNear(cylinder->getX2Minimum() - 0.25 * delta, 0.00005));
    EXPECT_THAT(grid->getStartZ(), RealNear(cylinder->getX3Minimum() - 0.25 * delta, 0.00005));

    EXPECT_THAT(grid->getEndX(), RealNear(cylinder->getX1Maximum() + 0.25 * delta, 0.00005));
    EXPECT_THAT(grid->getEndY(), RealNear(cylinder->getX2Maximum() + 0.25 * delta, 0.00005));
    EXPECT_THAT(grid->getEndZ(), RealNear(cylinder->getX3Maximum() + 0.25 * delta, 0.00005));
}

TEST_F(MultipleGridBuilderTestFixture, makeRotatingGrid_innerGridHasCorrectDimensions)
{
    const auto numberOfLevels = 1;
    const auto grid = gridBuilder.makeRotatingGrid(cylinder, numberOfLevels, 0);

    EXPECT_THAT(grid->getStartX(), RealNear(cylinder->getX1Minimum() - 0.5 * delta, 0.00005));
    EXPECT_THAT(grid->getStartY(), RealNear(cylinder->getX2Minimum() - 0.5 * delta, 0.00005));
    EXPECT_THAT(grid->getStartZ(), RealNear(cylinder->getX3Minimum() - 0.5 * delta, 0.00005));

    EXPECT_THAT(grid->getEndX(), RealNear(cylinder->getX1Maximum() + 0.5 * delta, 0.00005));
    EXPECT_THAT(grid->getEndY(), RealNear(cylinder->getX2Maximum() + 0.5 * delta, 0.00005));
    EXPECT_THAT(grid->getEndZ(), RealNear(cylinder->getX3Maximum() + 0.5 * delta, 0.00005));
}
