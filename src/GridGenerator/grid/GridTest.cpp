#include "gmock/gmock.h"
#include "GridImp.h"

#include <vector>
#include <GridGenerator/geometries/Vertex/Vertex.h>
#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>
#include "GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "GridStrategy/GridStrategyMocks.h"


TEST(GridTest, getMinimumOnNode)
{
    auto gridStrategyDummy = SPtr<GridStrategy>(new GridStrategyDummy);

    const real startX = 0.0;
    const real startY = 0.75;
    const real startZ = 1.75;
    const real delta = 0.5;
    const auto sut = GridImp::makeShared(NULL, startX, startY, startZ, 10, 10, 10, delta, gridStrategyDummy, Distribution(), 0);

    const real exactX = 0.4;
    const real exactY = 0.8;
    const real exactZ = 2.3;
    const Vertex actual = sut->getMinimumOnNode(Vertex(exactX, exactY, exactZ));

    const real expectedX = 0.0;
    const real expectedY = 0.75;
    const real expectedZ = 2.25;
    const Vertex expected = Vertex(expectedX, expectedY, expectedZ);

    EXPECT_TRUE(actual == expected);
}

TEST(GridTest, getMaximumOnNode)
{
    auto gridStrategyDummy = SPtr<GridStrategy>(new GridStrategyDummy);

    const real startX = 0.0;
    const real startY = 0.75;
    const real startZ = 0.75;
    const real delta = 0.5;
    const auto sut = GridImp::makeShared(NULL, startX, startY, startZ, 10, 10, 10, delta, gridStrategyDummy, Distribution(), 0);

    const real exactX = 0.4;
    const real exactY = 0.8;
    const real exactZ = 1.3;
    const Vertex actual = sut->getMaximumOnNode(Vertex(exactX, exactY, exactZ));

    const real expectedX = 0.5;
    const real expectedY = 1.25;
    const real expectedZ = 1.75;
    const Vertex expected = Vertex(expectedX, expectedY, expectedZ);

    EXPECT_TRUE(actual == expected);
}

TEST(GridTest, getBoundingBoxOnNodes)
{
    auto gridStrategyDummy = SPtr<GridStrategy>(new GridStrategyDummy);

    const real startX = -1.0;
    const real startY = -1.0;
    const real startZ = -1.0;
    const real delta = 0.25;
    const auto sut = GridImp::makeShared(NULL, startX, startY, startZ, 10, 10, 10, delta, gridStrategyDummy, Distribution(), 0);

    Triangle t = Triangle(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0.0f, 0.0f, 0.0f));

    const auto actual = sut->getBoundingBoxOnNodes(t);

    EXPECT_THAT(actual.maxX, RealEq(1.25));
}




//std::vector<Vertex> getPointsInBoundingBox(Triangle t, real delta)
//{
//	int x, y, z;
//	BoundingBox box = BoundingBox::makeNodeBox(t);
//
//	std::vector<Vertex> points;
//
//	for (x = box.minX; x <= box.maxX; x++) {
//		for (y = box.minY; y <= box.maxY; y++) {
//			for (z = box.minZ; z <= box.maxZ; z++) {
//				points.push_back(Vertex((real)x, (real)y, (real)z));
//			}
//		}
//	}
//
//	return points;
//}
//
//
//TEST(GridTest, transIndexToCoordsAndCoordToIndex)
//{
//	Grid grid(nullptr, 0, 0, 0, 10, 10, 10, Distribution());
//
//	int index = 756;
//	real x, y, z;
//	grid.transIndexToCoords(index, x, y, z);
//
//	real expectedX = 6;
//	real expectedY = 5;
//	real expectedZ = 7;
//
//	EXPECT_THAT(x, expectedX);
//	EXPECT_THAT(y, expectedY);
//	EXPECT_THAT(z, expectedZ);
//
//	unsigned int newIndex = grid.transCoordToIndex(expectedX, expectedY, expectedZ);
//
//	EXPECT_THAT(newIndex, index);
//}
//
//
//TEST(GridTest, checkIfPointIsOutOfRange)
//{
//	Grid grid(nullptr, 5, 5, 5, 10, 10, 10, Distribution());
//	EXPECT_TRUE(grid.isOutOfRange(Vertex(0, 0, 0)));
//}
//
//TEST(GridTest, checkIfPointIsNotOutOfRange)
//{
//	Grid grid(nullptr, 5, 5, 5, 10, 10, 10, Distribution());
//	EXPECT_FALSE(grid.isOutOfRange(Vertex(8, 8, 5)));
//}



#include <GridGenerator/grid/distributions/Distribution.h>
#include <GridGenerator/grid/NodeValues.h>

//
//class GridStopperTest : public testing::Test
//{
//public:
//    Grid grid;
//
//    void SetUp()
//    {
//        int nx = 2;
//        int ny = 2;
//        int nz = 2;
//        int size = nx*ny*nz;
//
//        char *field = new char[size]();
//        grid = Grid(field, 0, 0, 0, nx, ny, nz, DistributionHelper::getDistribution27());
//    }
//};
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfeveryNodeBeforeIsSolid_ItShouldNotBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_FALSE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusXisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = FLUID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusYisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = FLUID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusZisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = FLUID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusXYZisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = FLUID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusXYisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = FLUID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusYZisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = FLUID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = SOLID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
//
//TEST_F(GridStopperTest, testIfNodeIsStopper_IfMinusXZisFluid_ItShouldBeAStopperNode)
//{
//    grid.field[grid.transCoordToIndex(1, 1, 1)] = SOLID; // node
//
//    grid.field[grid.transCoordToIndex(0, 1, 1)] = SOLID; // -x
//    grid.field[grid.transCoordToIndex(1, 0, 1)] = SOLID; // -y
//    grid.field[grid.transCoordToIndex(1, 1, 0)] = SOLID; // -z
//    grid.field[grid.transCoordToIndex(0, 0, 0)] = SOLID; // -xyz
//    grid.field[grid.transCoordToIndex(0, 0, 1)] = SOLID; // -xy
//    grid.field[grid.transCoordToIndex(1, 0, 0)] = SOLID; // -yz
//    grid.field[grid.transCoordToIndex(0, 1, 0)] = FLUID; // -xz
//
//    ASSERT_TRUE(grid.isStopper(grid.transCoordToIndex(1, 1, 1)));
//}
