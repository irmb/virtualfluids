#include "gmock/gmock.h"

#include "MultipleGridBuilder.h"
#include "../GridMocks.h"
#include "../GridStrategy/GridStrategyMocks.h"
#include "../GridFactory.h"

class MultipleGridBuilderAddGridTest : public testing::Test
{
    virtual void SetUp() override
    {
        gridFactory = GridFactory::make();
        gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
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
    gridBuilder->addGrid(new Cuboid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0));
    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(1));
}

TEST_F(MultipleGridBuilderAddGridTest, givenCoarseGrid_addAdditionalGrid_shouldCalculateDelta)
{
    const real delta = 2.0;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, delta);

    gridBuilder->addGrid(new Cuboid(4.0, 4.0, 4.0, 6.0, 6.0, 6.0));

    ASSERT_THAT(gridBuilder->getDelta(1), RealEq(delta * 0.5));
}

TEST_F(MultipleGridBuilderAddGridTest, addGridWithoutCoarseGrid_shouldNotbeAdded)
{
    gridBuilder->addGrid(new Cuboid(0.0, 0.0, 0.0, 20.0, 20.0, 20.0));

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(0));
}

TEST_F(MultipleGridBuilderAddGridTest, getInvalidLevel_shouldThrowException)
{
    ASSERT_THROW(gridBuilder->getDelta(0), std::exception);
}

TEST_F(MultipleGridBuilderAddGridTest, addGridWithoutCoarseGrid_ShouldNotAddingAGrid)
{
    gridBuilder->addGrid(new Cuboid(0, 0, 0, 0, 0, 0), 0);

    ASSERT_THAT(gridBuilder->getNumberOfLevels(), testing::Eq(0));
}


TEST_F(MultipleGridBuilderAddGridTest, addGridWithFloatingStartPoints_ShouldCreatedStaggeredCoordinatesToCoarseGrid)
{
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1);

    gridBuilder->addGrid(new Cuboid(1.76, 1.76, 1.76, 10.51, 10.51, 10.51));


    EXPECT_THAT(gridBuilder->getStartX(1), RealEq(1.75));
    EXPECT_THAT(gridBuilder->getStartY(1), RealEq(1.75));
    EXPECT_THAT(gridBuilder->getStartZ(1), RealEq(1.75));

    EXPECT_THAT(gridBuilder->getEndX(1), RealEq(10.75));
    EXPECT_THAT(gridBuilder->getEndY(1), RealEq(10.75));
    EXPECT_THAT(gridBuilder->getEndZ(1), RealEq(10.75));
}



TEST(MultipleGridBuilderTest, everyExceptTheFinestGrid_shouldHaveAGridInterface)
{
    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
    gridFactory->setGridType(TestDouble::SPY);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);
    gridBuilder->addGrid(new Cuboid(7.0, 7.0, 7.0, 10.0, 10.0, 10.0), 2);

    gridBuilder->buildGrids();

    auto grid0 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(0));
    auto grid1 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(1));
    EXPECT_TRUE(grid0->hasGridInterface());
    EXPECT_TRUE(grid1->hasGridInterface());
}


TEST(MultipleGridBuilderTest, afterCreatingGridInterface_FineGridsShouldNotBeHavePeriodicBoundarys)
{
    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(SPtr<GridStrategy>(new GridStrategyDummy()));
    gridFactory->setGridType(TestDouble::SPY);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 1.0);
    gridBuilder->addGrid(new Cuboid(7.0, 7.0, 7.0, 10.0, 10.0, 10.0), 2);

    gridBuilder->buildGrids();

    auto grid1 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(1));
    auto grid2 = std::dynamic_pointer_cast<GridSpy>(gridBuilder->getGrid(2));

    EXPECT_TRUE(grid1->hasNoPeriodicityBoundaries());
    EXPECT_TRUE(grid2->hasNoPeriodicityBoundaries());
}
