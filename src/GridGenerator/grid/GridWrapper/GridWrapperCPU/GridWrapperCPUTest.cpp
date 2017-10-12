//#include "gmock/gmock.h"
//
//#include "GridKernelCPU.h"
//#include "GridKernelGPU.h"
//
//#include <string>
//
//#include <stl/Triangle.cuh>
//#include <stl/BoundingBox.cuh>
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
//void verifyGrid(Grid &grid, int min, int max, int expectedNode)
//{
//    int x, z, y;
//    // -x  +x 
//    for (y = min; y <= max; y++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            x = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//            x = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//        }
//    }
//
//    // -y  +y 
//    for (x = min; x <= max; x++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            y = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//            y = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//        }
//    }
//
//    // -z  +z 
//    for (x = min; x <= max; x++)
//    {
//        for (y = min; y <= max; y++)
//        {
//            z = min;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//            z = max;
//            EXPECT_THAT(grid.getFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z)), Eq(expectedNode));
//        }
//    }
//}
//
//void resetGridNodes(Grid &grid, int min, int max, int newNode)
//{
//    int x, z, y;
//    // -x  +x 
//    for (y = min; y <= max; y++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            x = min;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//            x = max;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//        }
//    }
//
//    // -y  +y 
//    for (x = min; x <= max; x++)
//    {
//        for (z = min; z <= max; z++)
//        {
//            y = min;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//            y = max;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//        }
//    }
//
//    // -z  +z 
//    for (x = min; x <= max; x++)
//    {
//        for (y = min; y <= max; y++)
//        {
//            z = min;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//            z = max;
//            grid.setFieldEntry(Vertex((doubflo)x, (doubflo)y, (doubflo)z), newNode);
//        }
//    }
//}
//
//void verifyInitalClearGridField(Grid &grid)
//{
//    unsigned int x, y, z;
//    for (unsigned int i = 0; i < grid.size; i += 10)
//    {
//        grid.transIndexToCoords(i, x, y, z);
//        if ((x == grid.nx - 1) || (y == grid.ny - 1) || (z == grid.nz - 1))
//            EXPECT_THAT(grid.field[i], Eq(1));
//        else
//            EXPECT_THAT(grid.field[i], Eq(0));
//    }
//}
//
//
//class GridKernelCubeIsExactlyOnNodesTest : public Test {
//public:
//    GridKernelCPU* gridCPU;
//    std::vector<Triangle> triangles;
//
//    void SetUp() {
//        doubflo length = 30.0f;
//        doubflo width = 30.0f;
//        doubflo high = 30.0f;
//        doubflo delta = 1.0f;
//        Transformator trans(delta, Vertex(0.0f, 0.0f, 0.0f));
//        std::string path = PATH_TO_DATA;
//        std::string test = TESTSUITE;
//        std::string input = path + test + "STL/GridKernelCPUTest/" + "cube_ascii.stl";
//
//        int nx = (int)(length / delta);
//        int ny = (int)(width / delta);
//        int nz = (int)(high / delta);
//        triangles = STLReaderWriter::readSTL(input, trans);
//
//        std::vector<BoundingBox> boxes = Partition::getProcessBoxes(1, nx, nx, nz);
//        gridCPU = new GridKernelCPU(boxes[0], "D3Q27", trans);
//        
//    }
//    void TearDown() {
//        delete gridCPU;
//    }
//};
//
//TEST_F(GridKernelCubeIsExactlyOnNodesTest, theInitalField_ShouldHaveZerosOnEveryEntry) {
//    verifyInitalClearGridField(gridCPU->grid);
//}
//
//TEST_F(GridKernelCubeIsExactlyOnNodesTest, meshCubeToGrid_expectSolidAndQRows_WallIsExactlyOnNodes) {
//    gridCPU->meshGrid(&triangles[0], (int)triangles.size());
//
//    int min = 5;
//    int max = 25;
//    int expectedQ = 6;
//    int expectedSolid = 1;
//    int fluidNode = 0;
//
//    verifyGrid(gridCPU->grid, min - 1, max + 1, expectedQ);
//    resetGridNodes(gridCPU->grid, min - 1, max + 1, fluidNode);
//    
//    verifyGrid(gridCPU->grid, min, max, expectedSolid);
//    resetGridNodes(gridCPU->grid, min, max, fluidNode);
//
//    verifyGrid(gridCPU->grid, min + 1, max - 1, expectedSolid);
//    resetGridNodes(gridCPU->grid, min + 1, max - 1, fluidNode);
//
//    verifyInitalClearGridField(gridCPU->grid);
//
//    /*std::vector<std::vector<doubflo> > qs_ausgeduennt = DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(*(gridCPU->grid), gridCPU->d);
//    DistributionHelper::printQs(qs_ausgeduennt, 2);*/
//}
//
//
//class GridKernelCubeIsBetweenNodesTest : public Test {
//public:
//	GridKernelCPU* gridCPU;
//    std::vector<Triangle> triangles;
//    
//    void SetUp() {
//        doubflo length = 30.0f;
//        doubflo width = 30.0f;
//        doubflo high = 30.0f;
//        doubflo delta = 1.0f;
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
//        gridCPU = new GridKernelCPU(boxes[0], "D3Q7", trans);
//    }
//
//    void TearDown()  {
//        delete gridCPU;
//    }
//
//
//};
//
//TEST_F(GridKernelCubeIsBetweenNodesTest, meshCubeToGrid_expectSolidAndQRows_WallIsBetweenTwoNode) {
//    gridCPU->meshGrid(&triangles[0], (int)triangles.size());
//
//    int min = 6;
//    int max = 25;
//    int expectedQ = 6;
//    int expectedSolid = 1;
//    int fluidNode = 0;
//
//    verifyGrid(gridCPU->grid, min - 1, max + 1, expectedQ);
//    resetGridNodes(gridCPU->grid, min - 1, max + 1, fluidNode);
//
//    verifyGrid(gridCPU->grid, min, max, expectedSolid);
//    resetGridNodes(gridCPU->grid, min, max, fluidNode);
//
//    verifyInitalClearGridField(gridCPU->grid);
//}
//
//
//TEST_F(GridKernelCubeIsBetweenNodesTest, calculateQ_validateAllQs_shouldBe_0_5_D3Q7){
//    gridCPU->meshGrid(&triangles[0], (int)triangles.size());
//
//    std::vector<std::vector<doubflo> > qs_ausgeduennt = DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(gridCPU->grid, gridCPU->grid.d);
//
//    for (int node = 0; node < qs_ausgeduennt.size(); node++) {
//        for (int dir = DIR_7_START; dir < DIR_7_END; dir++) {
//            doubflo q =  qs_ausgeduennt[node][dir+1];
//            if (q != 0.0f){
//                EXPECT_THAT(q, FloatEq(0.5));
//            }
//        }
//    }
//    //gridCPU->writeArrows("arrowsTest");
//    //GridVTKWriter::writeGridToVTK(grid, path + "cubeTestWall.vtk", trans, true);
//    //DistributionHelper::printQs(qs_ausgeduennt, 2);
//}
//
//
//class GridKernelCubeIsBetweenNodesD3Q27Test : public Test {
//public:
//	GridKernelCPU* gridCPU;
//    std::vector<Triangle> triangles;
//
//    void SetUp() {
//        doubflo length = 30.0f;
//        doubflo width = 30.0f;
//        doubflo high = 30.0f;
//        doubflo delta = 1.0f;
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
//        gridCPU = new GridKernelCPU(boxes[0], "D3Q27", trans);
//    }
//
//    void TearDown() {
//        delete gridCPU;
//    }
//
//};
//
//TEST_F(GridKernelCubeIsBetweenNodesD3Q27Test, calculateQ_validateAllQs_shouldBe_0_5_D3Q27){
//    gridCPU->meshGrid(&triangles[0], (int)triangles.size());
//
//    std::vector<std::vector<doubflo> > qs_ausgeduennt = DistributionHelper::getQsWithoutRowsWithOnlyZeroValues(gridCPU->grid, gridCPU->grid.d);
//
//    for (int node = 0; node < qs_ausgeduennt.size(); node++) {
//        for (int dir = DIR_7_START; dir < DIR_7_END; dir++) {
//            doubflo q = qs_ausgeduennt[node][dir + 1];
//            if (q != 0.0f){
//                EXPECT_THAT(q, FloatEq(0.5));
//            }
//        }
//    }
//    //gridCPU->writeArrows("arrowsTestQ27", Transformator(delta, getVertex(0.0f, 0.0f, 0.0f)));
//    //GridVTKWriter::writeGridToVTK(grid, path + "cubeTestWall.vtk", trans, true);
//    //printQs(qs_ausgeduennt, 2);
//
//}
//
//
//
//TEST_F(GridKernelCubeIsBetweenNodesTest, calculateQ_checkifQIsInFSolidNode_whenTriangleIsStraight_forD3Q7){
//    //GridkernelCPU kernel(&grid);
//
//    //Vertex p1 = getVertex(3.5, 2, 2);
//    //Vertex p2 = getVertex(3.5, 2, 9);
//    //Vertex p3 = getVertex(5.5, 9, 9);
//
//    //Vertex edge1, edge2, normal;
//    //edge1 = minus(p2, p1);
//    //edge2 = minus(p3, p2);
//    //normal = crossProduct(edge1, edge2);
//
//    //Triangle t = Triangle(p1,p2,p3,normal);
//
//    //kernel.meshGrid(&t, 1);
//    //GridVTKWriter::writeGridToVTK(grid, "firstd3Q7test.vtk", Transformator(), true);
//    ////writeArrowsToFile();
//
//    //UnstructuredGridWriter writer;
//    //doubflo v1_arr[3], v2_arr[3], v3_arr[3];
//    //convertVertexToArray(t.v1, v1_arr); convertVertexToArray(t.v2, v2_arr); convertVertexToArray(t.v3, v3_arr);
//    //writer.addTriangle(v1_arr, v2_arr, v3_arr);
//    //writer.writeUnstructuredGridToFile("triangleTest.vtu");
//
//
//
//    //VertexInteger v = getVertexInt(3, 5, 5);
//    //Distributions7 d;
//    //int dirs = 7;
//    //d.f[0] = new doubflo[grid.size * dirs];
//    ////calculateQs(grid, v, t, d.f[0]);
//}
