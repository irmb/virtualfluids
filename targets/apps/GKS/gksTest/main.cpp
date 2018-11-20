//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>

#include "core/Logger/Logger.h"

#include "GridGenerator/geometries/Cuboid/Cuboid.h"

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "GksVtkAdapter/VTKInterface.h"

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"
#include "GksGpu/Initializer/Initializer.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void gksTest( std::string path )
{
    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = 1.0 / 32.0;

    gridBuilder->addCoarseGrid(-0.5, -0.5, -0.5,  
                                0.5,  0.5,  0.5, dx);

    Cuboid cube(-1.0, -1.0, 0.45, 1.0, 1.0, 0.55);

    gridBuilder->setNumberOfLayers(10,8);
    gridBuilder->addGrid( &cube, 1);

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter;

    meshAdapter.inputGrid( gridBuilder );

    meshAdapter.findQuadtreeConnectivity( gridBuilder );

    meshAdapter.findCellToCellConnectivity( gridBuilder );

    meshAdapter.countCells();

    meshAdapter.partitionCells();

    meshAdapter.generateNodes( gridBuilder );

    meshAdapter.computeCellGeometry();

    meshAdapter.generateFaces( gridBuilder );

    meshAdapter.sortFaces();

    meshAdapter.countFaces();

    meshAdapter.generateInterfaceConnectivity();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    Parameters parameters;

    parameters.force.z = -0.01;
    
    auto dataBase = std::make_shared<DataBase>( "CPU" );
    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables {
        
        real radius = cellCenter.length();

        return toConservedVariables( PrimitiveVariables( 1.0, 0.0, 0.0, 0.0, 1.0, radius ), parameters.K );

    });

    Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataHostToDevice();

    writeVtkXML( dataBase, parameters, 0, path + "grid/Test_0" );

    //////////////////////////////////////////////////////////////////////////

    for( uint iter = 0; iter < 100; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);
    }


    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );


}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/gridGenerator/" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    
    try
    {
        gksTest( path );
    }
    catch (const std::exception& e)
    {     
        *logging::out << logging::Logger::ERROR << e.what() << "\n";
    }
    catch (const std::bad_alloc& e)
    {  
        *logging::out << logging::Logger::ERROR << "Bad Alloc:" << e.what() << "\n";
    }
    catch (...)
    {
        *logging::out << logging::Logger::ERROR << "Unknown exception!\n";
    }

   return 0;
}
