#include "basics/tests/testUtilities.h"
#include "geometries/Sphere/Sphere.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"

class MultipleGridBuilderTestFixture : public testing::Test
{
protected:
    MultipleGridBuilder gridBuilder;
    SPtr<Object> gridShape;
    real delta = 0.1;

public:
    void SetUp() override
    {
        gridShape = std::make_shared<Sphere>(0.0, 0.0, 0.0, 0.0);
        gridBuilder.addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, delta);
    }
};

TEST_F(MultipleGridBuilderTestFixture, addCoarseGrid_addsOneGridToGridList)
{
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(1));
}

TEST(MultipleGridBuilderTest, noCoarseGrid_addGrid_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<Object> gridShape;

    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape, 1);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGrid_doesNotWarn)
{
    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape, 1);
    EXPECT_FALSE(testingVF::stdoutContainsWarning());
}

TEST(MultipleGridBuilderTest, noCoarseGrid_addGridWithoutLevel_warns)
{
    MultipleGridBuilder gridBuilder;
    SPtr<Object> gridShape;

    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape);
    EXPECT_TRUE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGridWitoutLevel_doesNotWarn)
{
    testingVF::captureStdOut();
    gridBuilder.addGrid(gridShape);
    EXPECT_FALSE(testingVF::stdoutContainsWarning());
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGrid_addsGridToList)
{
    uint levelFine = 1;

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + 1 new level = 2 levels";
}

TEST_F(MultipleGridBuilderTestFixture, coarseGridExist_addGridAtLevelTwo_addsTwoGridsToList)
{
    uint levelFine = 2;

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3)) << "coarse level + 2 new levels = 3 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtLevelTwo_addsOneGridToList)
{
    uint levelFine = 2;

    gridBuilder.addGrid(gridShape, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3)) << "coarse level + fineLevel + 1 new levels = 3 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtLevelThree_addsTwoGridsToList)
{
    uint levelFine = 3;

    gridBuilder.addGrid(gridShape, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(4)) << "coarse level + fineLevel + 2 new levels = 4 levels";
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExist_addGridAtSameLevel_noGridAdded)
{
    uint levelFine = 1;

    gridBuilder.addGrid(gridShape, 1);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2)) << "coarse level + fineLevel = 2 levels";

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2))
        << "Can't add grid on same level twice --> no new grid should be added";
}

TEST_F(MultipleGridBuilderTestFixture, addGrid_hasHalfDeltaComparedToCoarseGrid)
{
    uint levelFine = 1;

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrid(levelFine)->getDelta(), testing::Eq(0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, addGridAtLevelTwo_hasHalfOfHalfDeltaComparedToCoarseGrid)
{
    uint levelFine = 2;

    gridBuilder.addGrid(gridShape, levelFine);
    EXPECT_THAT(gridBuilder.getGrid(levelFine)->getDelta(), testing::Eq(0.5 * 0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, addGridWithoutLevel_addsGridAtLevelOne)
{
    gridBuilder.addGrid(gridShape);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2));
    EXPECT_THAT(gridBuilder.getGrid(1)->getDelta(), testing::Eq(0.5 * delta));
}

TEST_F(MultipleGridBuilderTestFixture, fineGridExists_addGridWithoutLevel_addsGridAtLevelTwo)
{
    gridBuilder.addGrid(gridShape);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(2));
    EXPECT_THAT(gridBuilder.getGrid(1)->getDelta(), testing::Eq(0.5 * delta));
    gridBuilder.addGrid(gridShape);
    EXPECT_THAT(gridBuilder.getGrids().size(), testing::Eq(3));
    EXPECT_THAT(gridBuilder.getGrid(2)->getDelta(), testing::Eq(0.5 * 0.5 * delta));
}
