#include <array>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <memory>
#include <ostream>

#include "GridImp.h"
#include "PointerDefinitions.h"
#include "grid/Field.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/distributions/Distribution.h"

// This test is commented out because it causes a compiler error in Clang 10 --> The bug is fixed in Clang 14 (https://github.com/google/googletest/issues/2271)

// class FieldDouble : public Field
// {
// public:
//     FieldDouble() : Field(1)
//     {
//         this->allocateMemory();
//     };

//     void setToStopper(uint index)
//     {
//         this->field[index] = vf::gpu::STOPPER_SOLID;
//     }
// };

// class GridImpDouble : public GridImp
// {
// public:
//     std::array<real, 3> coordsOfTestedNode;
//     GridImpDouble(Object *object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
//                   Distribution d, uint level)
//         : GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, d, level)
//     {
//         this->neighborIndexX = new int[5];
//         this->neighborIndexY = new int[5];
//         this->neighborIndexZ = new int[5];
//     }

//     static SPtr<GridImpDouble> makeShared(Object *object, real startX, real startY, real startZ, real endX, real endY,
//                                           real endZ, real delta, Distribution d, uint level)
//     {
//         SPtr<GridImpDouble> grid(new GridImpDouble(object, startX, startY, startZ, endX, endY, endZ, delta, d, level));
//         return grid;
//     }

//     void transIndexToCoords(uint, real &x, real &y, real &z) const override
//     {
//         x = coordsOfTestedNode[0];
//         y = coordsOfTestedNode[1];
//         z = coordsOfTestedNode[2];
//     }

//     uint transCoordToIndex(const real &, const real &, const real &) const override
//     {
//         return 0;
//     }

//     void setStopperNeighborCoords(uint index) override
//     {
//         GridImp::setStopperNeighborCoords(index);
//     }

//     void setField(Field &field)
//     {
//         this->field = field;
//     }

//     MOCK_METHOD(int, getSparseIndex, (const real &x, const real &y, const real &z), (const, override));
// };

// // This is test is highly dependent on the implementation. Maybe it should be removed :(
// TEST(GridImp, setStopperNeighborCoords)
// {
//     real end = 1.0;
//     real delta = 0.1;

//     SPtr<GridImpDouble> gridImp =
//         GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, end, end, end, delta, Distribution(), 0);
//     FieldDouble field;
//     field.setToStopper(0);
//     gridImp->setField(field);

//     gridImp->coordsOfTestedNode = { end - ((real)0.5 * delta), end - ((real)0.5 * delta), end - ((real)0.5 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)0.51 * delta), end - ((real)0.51 * delta),
//                                     end - ((real)0.51 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);
//     gridImp->coordsOfTestedNode = { end - ((real)0.99 * delta), end - ((real)0.99 * delta),
//                                     end - ((real)0.99 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - delta, end - delta, end - delta };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)1.01 * delta), end - ((real)1.01 * delta),
//                                     end - ((real)1.01 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     // The grid should not be like this, so this should be fine...
//     gridImp->coordsOfTestedNode = { end, end, end };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(0);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)0.25 * delta), end - ((real)0.25 * delta),
//                                     end - ((real)0.25 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(0);
//     gridImp->setStopperNeighborCoords(0);
// }

std::array<int, 3> countInvalidNeighbors(SPtr<Grid> grid)
{
    auto countInvalidX = 0;
    auto countInvalidY = 0;
    auto countInvalidZ = 0;
    for (uint index = 0; index < grid->getSize(); index++) {
        if (grid->getNeighborsX()[index] == -1)
            countInvalidX++;
        if (grid->getNeighborsY()[index] == -1)
            countInvalidY++;
        if (grid->getNeighborsZ()[index] == -1)
            countInvalidZ++;
    }
    return { countInvalidX, countInvalidY, countInvalidZ };
}

std::array<int, 3> testFluidNodeNeighbors(SPtr<Grid> grid)
{
    auto countInvalidX = 0;
    auto countInvalidXY = 0;
    auto countInvalidXYZ = 0;
    for (uint index = 0; index < grid->getSize(); index++) {
        if (grid->getFieldEntry(index) != vf::gpu::FLUID) {
            continue;
        }

        auto neighX = grid->getNeighborsX()[index];
        if (neighX == -1) {
            countInvalidX++;
            continue;
        }

        auto neighXY = grid->getNeighborsY()[neighX];
        if (neighXY == -1) {
            countInvalidXY++;
            continue;
        }

        auto neighXYZ = grid->getNeighborsZ()[neighXY];
        if (neighXYZ == -1) {
            countInvalidXYZ++;
            continue;
        }
    }

    return { countInvalidX, countInvalidXY, countInvalidXYZ };
}

class findNeighborsIntegrationTest : public ::testing::Test
{
protected:
    SPtr<MultipleGridBuilder> gridBuilder;

    void SetUp() override
    {
        auto gridFactory = GridFactory::make();
        gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
        gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

        // init logger to avoid segmentation fault in buildGrids
        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::Level::WARNING);
        logging::Logger::timeStamp(logging::Logger::ENABLE);
        logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);
    }
};

TEST_F(findNeighborsIntegrationTest, grid1)
{
    const real dx = 0.15;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    // Only the last layer of nodes should have invalid neighbors. The grid is a cube with a side length of 9 nodes
    // -> 9 * 9 = 81 invalid nodes are expected
    auto numberOfInvalidNeighbors = countInvalidNeighbors(grid);
    auto expected = 9 * 9;
    EXPECT_THAT(numberOfInvalidNeighbors[0], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[1], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[2], testing::Eq(expected));

    // additional test: all fluid nodes should have valid neighbors
    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, grid2)
{
    const real dx = 1.0 / 64;
    gridBuilder->addCoarseGrid(-0.6, -0.6, -0.6, 0.6, 0.6, 0.6, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    // Only the last layer of nodes should have invalid neighbors. The grid is a cube with a side length of 79 nodes
    // -> 79 * 79 invalid nodes are expected
    auto numberOfInvalidNeighbors = countInvalidNeighbors(grid);
    auto expected = 79 * 79;
    EXPECT_THAT(numberOfInvalidNeighbors[0], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[1], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[2], testing::Eq(expected));

    // additional test: all fluid nodes should have valid neighbors
    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, validFluidNeighbors1)
{
    real dx = 0.17;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, validFluidNeighbors2)
{
    real dx = 0.18;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}
