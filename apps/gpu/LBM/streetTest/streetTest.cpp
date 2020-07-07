
#include <iostream>

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"


#include "GridGenerator/geometries/TriangularMesh/TriangularMeshStrategy.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/GridFactory.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

int main( int argc, char* argv[])
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_INTERMEDIATE);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

	real dx = 1.0;

	gridBuilder->addCoarseGrid(-256, -256, -10,
								256,  256,  40, dx);

    TriangularMesh* flatGroundSTL = TriangularMesh::make("F:/Work/Computations/NagelSchreckenberg/FlatGround.stl");

    gridBuilder->addGeometry(flatGroundSTL);

	gridBuilder->setPeriodicBoundaryCondition(true, true, false);
	
	gridBuilder->buildGrids(LBM, false); 
	
    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
	
	gridBuilder->writeGridsToVtk("F:/Work/Computations/NagelSchreckenberg/ExampleGrid");

    SimulationFileWriter::write("F:/Work/Computations/NagelSchreckenberg/grid/", gridBuilder, FILEFORMAT::BINARY);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StreetPointFinder finder;

    finder.readStreets( "F:/Work/Computations/NagelSchreckenberg/ExampleStreets.txt" );

    finder.writeVTK( "F:/Work/Computations/NagelSchreckenberg/ExampleStreets.vtk" );

    finder.findIndicesLB( gridBuilder->getGrid(0) );

    finder.writeConnectionVTK( "F:/Work/Computations/NagelSchreckenberg/ExampleStreetsConnection.vtk", gridBuilder->getGrid(0) );

    finder.writeSimulationFile("F:/Work/Computations/NagelSchreckenberg/grid/", 1.0, gridBuilder->getNumberOfLevels(), 0);

    return 0;
}
