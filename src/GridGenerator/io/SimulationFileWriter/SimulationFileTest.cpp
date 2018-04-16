//#include "gmock/gmock.h"
//
//#include <GridGenerator/geometries/Vertex/Vertex.cuh>
//#include <GridGenerator/geometries/Triangle/Triangle.h>
//#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>
//#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>
//
//#include <GridGenerator/grid/kernel/GridKernelCPU.h>
//#include <GridGenerator/grid/partition/Partition.h>
//#include <GridGenerator/utilities/Transformator/Transformator.h>
//#include <GridGenerator/io/STLReaderWriter/STLReaderWriter.h>
//#include "SimulationFileWriter.h"
//#include "UnstructuredGridBuilder.h"
//#include "SimulationFileNames.h"
//
////
//using namespace testing;
//
//class SimulationFileTest : public Test {
//
//public:
//    std::vector<unsigned int> neighborXFoam, neighborYFoam, neighborZFoam, geoFoam;
//    std::vector<unsigned int> neighborX, neighborY, neighborZ, geo;
//
//    std::vector<real> coordXFoam, coordYFoam, coordZFoam;
//    std::vector<real> coordX, coordY, coordZ;
//
//    std::vector<std::vector<real> > inletQFoam, outletQFoam;
//    std::vector<std::vector<real> > inletQ, outletQ;
//
//    void SetUp() {
//        real length = 10.0f;
//        real width = 10.0f;
//        real high = 10.0f;
//        real delta = 1.0f;
//
//        int nx = (int)(length / delta) + 1;
//        int ny = (int)(width / delta) + 1;
//        int nz = (int)(high / delta) + 1;
//
//        Transformator trans = Transformator(delta, Vertex(0.0f, 0.0f, 0.0f));
//        std::string path = PATH_TO_DATA;
//        std::string folder = "TESTSUITE/SIMULATION_FILES/";
//        std::string input = path + folder + "gridcubeSimulation.stl";
//        std::vector<BoundingBox<int>> boxes = Partition::getProcessBoxes(1, nx, ny, nz);
//
//		GridKernelCPU gridCPU(boxes[0], "D3Q27", trans, true);
//		Geometry geom(input, boxes[0], trans);
//        if (geom.size > 0)
//            gridCPU.meshGrid(geom);
//
//        gridCPU.floodFill(Vertex(5,5,5));
//
//        UnstructuredGridBuilder builder;
//		BoundingBox<int> box;
//        builder.buildUnstructuredGrid(gridCPU.grid, box);
//
//        bool binaer = false;
//        Transformator dummy;
//        std::vector<Node> coords = builder.getNodes();
//        std::vector<std::vector<std::vector<real> > > qs = builder.getQsValues();
//        SimulationFileWriter::writeSimulationFiles(path + folder + "gridGeneration/", coords, qs, binaer, gridCPU.grid, dummy);
//
//        //open OpenFoam files
//        CoordNeighborGeoReader cngNXFoam(path + folder + "openFoam/" + simulationFileNames::neighborX, binaer, false);
//        CoordNeighborGeoReader cngNYFoam(path + folder + "openFoam/" + simulationFileNames::neighborY, binaer, false);
//        CoordNeighborGeoReader cngNZFoam(path + folder + "openFoam/" + simulationFileNames::neighborZ, binaer, false);
//        neighborXFoam = cngNXFoam.getNeighbors(0);
//        neighborYFoam = cngNYFoam.getNeighbors(0);
//        neighborZFoam = cngNZFoam.getNeighbors(0);
//        CoordNeighborGeoReader cngXFoam(path + folder + "openFoam/" + simulationFileNames::coordX, binaer, true);
//        CoordNeighborGeoReader cngYFoam(path + folder + "openFoam/" + simulationFileNames::coordY, binaer, true);
//        CoordNeighborGeoReader cngZFoam(path + folder + "openFoam/" + simulationFileNames::coordZ, binaer, true);
//        coordXFoam = cngXFoam.getCoords(0);
//        coordYFoam = cngYFoam.getCoords(0);
//        coordZFoam = cngZFoam.getCoords(0);
//        CoordNeighborGeoReader cngGeoFoam(path + folder + "openFoam/" + simulationFileNames::geoVec, binaer, false);
//        geoFoam = cngGeoFoam.getNeighbors(0);
//
//        //open grid generation files
//        CoordNeighborGeoReader cngNX(path + folder + "gridGeneration/" + simulationFileNames::neighborX, binaer, false);
//        CoordNeighborGeoReader cngNY(path + folder + "gridGeneration/" + simulationFileNames::neighborY, binaer, false);
//        CoordNeighborGeoReader cngNZ(path + folder + "gridGeneration/" + simulationFileNames::neighborZ, binaer, false);
//        neighborX = cngNX.getNeighbors(0);
//        neighborY = cngNY.getNeighbors(0);
//        neighborZ = cngNZ.getNeighbors(0);
//        CoordNeighborGeoReader cngX(path + folder + "gridGeneration/" + simulationFileNames::coordX, binaer, true);
//        CoordNeighborGeoReader cngY(path + folder + "gridGeneration/" + simulationFileNames::coordY, binaer, true);
//        CoordNeighborGeoReader cngZ(path + folder + "gridGeneration/" + simulationFileNames::coordZ, binaer, true);
//        coordX = cngX.getCoords(0);
//        coordY = cngY.getCoords(0);
//        coordZ = cngZ.getCoords(0);
//        CoordNeighborGeoReader cngGeo(path + folder + "gridGeneration/" + simulationFileNames::geoVec, binaer, false);
//        geo = cngGeo.getNeighbors(0);
//
//        BoundaryQsReader BCInletFoam(path + folder + "openFoam/" + simulationFileNames::inletBoundaryQ, false);
//        inletQFoam = BCInletFoam.getQs(0);
//
//        BoundaryQsReader BCInlet(path + folder + "gridGeneration/" + simulationFileNames::inletBoundaryQ, false);
//        inletQ = BCInlet.getQs(0);
//
//
//        BoundaryQsReader BCOutletFoam(path + folder + "openFoam/" + simulationFileNames::outletBoundaryQ, false);
//        outletQFoam = BCOutletFoam.getQs(0);
//
//        BoundaryQsReader BCOutlet(path + folder + "gridGeneration/" + simulationFileNames::outletBoundaryQ, false);
//        outletQ = BCOutlet.getQs(0);
//
//    }
//
//    void TearDown()  {
//    }
//
//}; 
//
//
//TEST_F(SimulationFileTest, testFileSizeFromOpenFoamWithGridGeneration)
//{
//    bool filesizes = neighborXFoam.size() == neighborX.size() && neighborYFoam.size() == neighborY.size() && neighborZFoam.size() == neighborZ.size() && geoFoam.size() == geo.size();
//
//	EXPECT_TRUE(filesizes);
//
//    for (int i = 1; i < neighborXFoam.size(); i++) {
//        EXPECT_THAT(coordXFoam[i] - 1.5f, DoubleEq(coordX[i]));
//        EXPECT_THAT(coordYFoam[i] - 1.5f, DoubleEq(coordY[i]));
//        EXPECT_THAT(coordZFoam[i] - 1.5f, DoubleEq(coordZ[i]));
//
//        EXPECT_THAT(neighborXFoam[i], Eq(neighborX[i]));
//        EXPECT_THAT(neighborYFoam[i], Eq(neighborY[i]));
//        EXPECT_THAT(neighborZFoam[i], Eq(neighborZ[i]));
//        EXPECT_THAT(geoFoam[i], Eq(geo[i]));
//    }
//
//    for (int i = 0; i < inletQFoam.size(); i++) {
//        bool size = inletQFoam[i].size() == inletQ[i].size();
//        EXPECT_TRUE(size);
//    }
//
//    for (int i = 0; i < inletQFoam.size(); i++) {
//        for (int j = 0; j < inletQFoam[i].size(); j++) {
//
//            EXPECT_THAT(inletQFoam[i][j], DoubleEq(inletQ[i][j]));
//        }
//    }
//
//    for (int i = 0; i < inletQFoam.size(); i++) {
//        for (int j = 0; j < inletQFoam[i].size(); j++) {
//
//            EXPECT_THAT(outletQFoam[i][j], DoubleEq(outletQ[i][j]));
//        }
//    }
//}
//
