//#include "gmock/gmock.h"
//#include "GridKernelCPU.h"
//#include "GridKernelGPU.h"
//
//#include <string>
//
//#include <stl/Triangle.h>
//#include <stl/BoundingBox.h>
//
//#include <grid/kernel/runGridKernelGPU.cuh>
//#include <grid/distributions/D3Q7.h>
//#include <utilities/Transformator.h>
//#include <utilities/io/STLReaderWriter.h>
//#include <utilities/io/GridVTKWriter.h>
//
//
//#include <grid/distributions/Distribution.h>
//#include <grid/partition/Partition.h>
//
//using namespace testing;
//
//
//#ifndef __unix__
//
//void verifyGridGPU(Grid &grid, int min, int max, int expectedNode)
//{
//    int x, z, y;
//    // -x  +x 
//    for (y = min; y <= max; y++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            x = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//            x = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//        }
//    }
//
//    // -y  +y 
//    for (x = min; x <= max; x++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            y = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//            y = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//        }
//    }
//
//    // -z  +z 
//    for (x = min; x <= max; x++)
//    {
//        for (y = min; y <= max; y++)
//        {
//            z = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//            z = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((real)x, (real)y, (real)z)), Eq(expectedNode));
//        }
//    }
//}
//
//void resetGridNodesGPU(Grid &grid, int min, int max, int newNode)
//{
//    int x, z, y;
//    // -x  +x 
//    for (y = min; y <= max; y++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            x = min;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//            x = max;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//        }
//    }
//
//    // -y  +y 
//    for (x = min; x <= max; x++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            y = min;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//            y = max;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//        }
//    }
//
//    // -z  +z 
//    for (x = min; x <= max; x++)
//    {
//        for (y = min; y <= max; y++)
//        {
//            z = min;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//            z = max;
//            grid.setFieldEntry(Vertex((real)x, (real)y, (real)z), newNode);
//        }
//    }
//}
//
//void verifyInitalClearGridFieldGPU(Grid &grid)
//{
//    unsigned int x, y, z;
//    for (unsigned int i = 0; i < grid.size; i += 10)
//    {
//        grid.transIndexToCoords(i, x, y, z);
//        if ((x == grid.nx - 1) || (y == grid.ny - 1) || (z == grid.nz - 1))
//            EXPECT_THAT(grid.field[i], Eq(0));
//        else
//            EXPECT_THAT(grid.field[i], Eq(0));
//    }
//}
//
//
//
//class GridKernelCubeIsExactlyOnNodesGPUTest : public Test {
//public:
//    GridKernelGPU* gridGPU;
//    std::vector<Triangle> triangles;
//
//    void SetUp() {
//        real length = 30.0f;
//        real width = 30.0f;
//        real high = 30.0f;
//        real delta = 1.0f;
//        Transformator trans(delta, Vertex(0.0f, 0.0f, 0.0f));
//        std::string path = PATH_TO_DATA;
//        std::string test = TESTSUITE;
//        std::string input = path + test + "STL/GridKernelCPUTest/" + "cube_ascii.stl";
//
//        int nx = (int)(length / delta);
//        int ny = (int)(width / delta);
//        int nz = (int)(high / delta);
//        triangles = STLReaderWriter::readSTL(input, trans);
//        std::vector<BoundingBox> boxes = Partition::getProcessBoxes(1, nx, nx, nz);
//
//        gridGPU = new GridKernelGPU(boxes[0], "D3Q27", Transformator());
//
//    }
//    void TearDown() {
//        delete gridGPU;
//    }
//};
//
//
//TEST_F(GridKernelCubeIsExactlyOnNodesGPUTest, meshCubeToGrid_expectSolidAndQRows_WallIsExactlyOnNodes_onGPU) {
//    gridGPU->meshGrid(&triangles[0], (int)triangles.size());
//    gridGPU->copyDataFromGPU();
//
//	std::string path = PATH_TO_DATA;
//	GridVTKWriter::writeGridToVTK(gridGPU->grid, path + "VTK_OUTPUT/testcube" , Transformator(1.0f, Vertex(0.0f, 0.0f, 0.0f)), false);
//
//    int min = 5;
//    int max = 25;
//    int expectedQ = 6;
//    int expectedSolid = 1;
//    int fluidNode = 0;
//
//    verifyGridGPU(gridGPU->grid, min - 1, max + 1, expectedQ);
//    resetGridNodesGPU(gridGPU->grid, min - 1, max + 1, fluidNode);
//
//    verifyGridGPU(gridGPU->grid, min, max, expectedSolid);
//    resetGridNodesGPU(gridGPU->grid, min, max, fluidNode);
//
//    verifyGridGPU(gridGPU->grid, min + 1, max - 1, expectedSolid);
//    resetGridNodesGPU(gridGPU->grid, min + 1, max - 1, fluidNode);
//
//    verifyInitalClearGridFieldGPU(gridGPU->grid);
//
//}
//
//
//class GridKernelCubeIsBetweenNodesGPUTest : public Test {
//public:
//    GridKernelGPU* gridGPU;
//    std::vector<Triangle> triangles;
//
//    void SetUp() {
//        real length = 30.0f;
//        real width = 30.0f;
//        real high = 30.0f;
//        real delta = 1.0f;
//        Transformator trans;
//        std::string path = PATH_TO_DATA;
//        std::string test = TESTSUITE;
//        std::string input = path + test + "STL/GridKernelCPUTest/" + "cubeBetweenNode_ascii.stl";
//
//        trans = Transformator(delta, Vertex(0.0f, 0.0f, 0.0f));
//        int nx = (int)(length / delta);
//        int ny = (int)(width / delta);
//        int nz = (int)(high / delta);
//        triangles = STLReaderWriter::readSTL(input, trans);
//        std::vector<BoundingBox> boxes = Partition::getProcessBoxes(1, nx, nx, nz);
//
//        gridGPU = new GridKernelGPU(boxes[0], "D3Q7", Transformator());
//    }
//
//    void TearDown()  {
//        delete gridGPU;
//    }
//
//
//};
//
//
//TEST_F(GridKernelCubeIsBetweenNodesGPUTest, meshCubeToGrid_expectSolidAndQRows_WallIsBetweenTwoNode_onGPU) {
//    gridGPU->meshGrid(&triangles[0], (int)triangles.size());
//    gridGPU->copyDataFromGPU();
//
//    int min = 6;
//    int max = 25;
//    int expectedQ = 6;
//    int expectedSolid = 1;
//    int fluidNode = 0;
//
//    verifyGridGPU(gridGPU->grid, min - 1, max + 1, expectedQ);
//    resetGridNodesGPU(gridGPU->grid, min - 1, max + 1, fluidNode);
//
//    verifyGridGPU(gridGPU->grid, min, max, expectedSolid);
//    resetGridNodesGPU(gridGPU->grid, min, max, fluidNode);
//
//    verifyInitalClearGridFieldGPU(gridGPU->grid);
//}
//
//
//TEST_F(GridKernelCubeIsBetweenNodesGPUTest, calculateQ_validateAllQs_shouldBe_0_5_onGPU_D3Q7){
//    gridGPU->meshGrid(&triangles[0], (int)triangles.size());
//    gridGPU->copyDataFromGPU();
//
//    std::vector<std::vector<real> > qs_ausgeduennt = DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(gridGPU->grid, gridGPU->grid.d);
//
//    for (int node = 0; node < qs_ausgeduennt.size(); node++) {
//        for (int dir = DIR_7_START; dir < DIR_7_END; dir++) {
//            real q = qs_ausgeduennt[node][dir + 1];
//            if (q != 0.0f){
//                EXPECT_THAT(q, DoubleEq(0.5));
//            }
//        }
//    }
//}
//
//class GridKernelCubeIsBetweenNodesD3Q27GPUTest : public Test {
//public:
//    GridKernelGPU* gridGPU;
//    std::vector<Triangle> triangles;
//
//    void SetUp() {
//        real length = 30.0f;
//        real width = 30.0f;
//        real high = 30.0f;
//        real delta = 1.0f;
//        Transformator trans;
//        std::string path = PATH_TO_DATA;
//        std::string test = TESTSUITE;
//        std::string input = path + test + "STL/GridKernelCPUTest/" + "cubeBetweenNode_ascii.stl";
//
//        trans = Transformator(delta, Vertex(0.0f, 0.0f, 0.0f));
//        int nx = (int)(length / delta);
//        int ny = (int)(width / delta);
//        int nz = (int)(high / delta);
//        triangles = STLReaderWriter::readSTL(input, trans);
//        std::vector<BoundingBox> boxes = Partition::getProcessBoxes(1, nx, nx, nz);
//        gridGPU = new GridKernelGPU(boxes[0], "D3Q27", Transformator());
//    }
//
//    void TearDown() {
//        delete gridGPU;
//    }
//
//};
//
//
//TEST_F(GridKernelCubeIsBetweenNodesD3Q27GPUTest, calculateQ_validateAllQs_shouldBe_0_5_onGPU_D3Q27) {
//    gridGPU->meshGrid(&triangles[0], (int)triangles.size());
//    gridGPU->copyDataFromGPU();
//
//    std::vector<std::vector<real> > qs_ausgeduennt = DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(gridGPU->grid, gridGPU->grid.d);
//
//    for (int node = 0; node < qs_ausgeduennt.size(); node++) {
//        for (int dir = DIR_7_START; dir < DIR_7_END; dir++) {
//            real q = qs_ausgeduennt[node][dir + 1];
//            if (q != 0.0f){
//                EXPECT_THAT(q, DoubleEq(0.5));
//            }
//        }
//    }
//}
//
//
//class GridKernelTest : public Test {
//public:
//    real length;
//    real width;
//    real high;
//    real delta;
//    unsigned int nx;
//    unsigned int ny;
//    unsigned int nz;
//    Grid grid;
//    Distribution d;
//    void SetUp() {
//        length = 10;
//        width = 5;
//        high = 4;
//        delta = 0.1f;
//        nx = (unsigned int)(length / delta);
//        ny = (unsigned int)(width / delta);
//        nz = (unsigned int)(high / delta);
//        d = DistributionHelper::getDistribution7();
//        grid = Grid(NULL, 0, 0, 0, nx, ny, nz, d);
//    }
//};
//
//
//
//TEST_F(GridKernelTest, transCoordToIndexAndIndexToCoordinates) {
//    int x = nx / 2;
//    int y = ny - 1;
//    int z = nz - 1;
//    unsigned int index = grid.transCoordToIndex(Vertex((real)x, (real)y, (real)z));
//
//    unsigned int newX, newY, newZ;
//    grid.transIndexToCoords(index, newX, newY, newZ);
//
//    EXPECT_THAT(x, Eq(newX));
//    EXPECT_THAT(y, Eq(newY));
//    EXPECT_THAT(z, Eq(newZ));
//}
//
//TEST_F(GridKernelTest, ifAllBetaAreSmallerThenAlpha_ItShouldreturnTrue) {
//
//    real alphaAngles[3];
//    real betaAngles[3];
//    alphaAngles[0] = 94.4f;
//    alphaAngles[1] = 92.4f;
//    alphaAngles[2] = 91.5f;
//
//    betaAngles[0] = 0.0f;
//    betaAngles[1] = 0.0f;
//    betaAngles[2] = 0.0f;
//
//	Grid grid;
//    ASSERT_TRUE(grid.isBetaSmallerThanAlpha(alphaAngles, betaAngles));
//}
//
//TEST_F(GridKernelTest, ifOneBetaIsBiggerThenAlpha_ItShouldreturnFalse) {
//
//    real alphaAngles[3];
//    real betaAngles[3];
//
//    alphaAngles[0] = 94.4f;
//    alphaAngles[1] = 92.4f;
//    alphaAngles[2] = 91.5f;
//
//    betaAngles[0] = 100.0f;
//    betaAngles[1] = 0.0f;
//    betaAngles[2] = 0.0f;
//
//	Grid grid;
//    ASSERT_FALSE(grid.isBetaSmallerThanAlpha(alphaAngles, betaAngles));
//}
//
//#endif
