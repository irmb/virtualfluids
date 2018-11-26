//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

#include "GridGenerator/geometries/Cuboid/Cuboid.h"

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "GksVtkAdapter/VTKInterface.h"

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"
#include "GksGpu/Initializer/Initializer.h"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void gksTest( std::string path )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L / 128.0;

    real Re  = 2.0e3;
    real U  = 0.1;
    real Ma = 0.1;
    
    real Pr  = 1.0;
    real K   = 0.0;

    real rho = 1.0;

    real mu = U * rho * L / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( cs * cs );

    real CFL = 0.5;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = lambda;

    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(-0.5, -0.5, -0.5,  
                                0.5,  0.5,  0.5, dx);

    //Cuboid cube(-1.0, -1.0, 0.45, 1.0, 1.0, 0.55);

    //gridBuilder->setNumberOfLayers(6,6);
    //gridBuilder->addGrid( &cube, 1);

    gridBuilder->setPeriodicBoundaryCondition(true, false, false);

    gridBuilder->buildGrids(GKS, false);

    gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    auto dataBase = std::make_shared<DataBase>( "GPU" );
    dataBase->setMesh( meshAdapter );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );

    //bcMX->findBoundaryCells( meshAdapter, [&](Vec3 center){ 
    //    return center.x < -0.5 || center.x > 0.5;
    //} );

    SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3( U, U, 0.0 ), lambda, 0.0 );

    bcPZ->findBoundaryCells( meshAdapter, [&](Vec3 center){ 
        return center.z > 0.5;
    } );

    SPtr<BoundaryCondition> bcWall = std::make_shared<IsothermalWall>( dataBase, Vec3( 0.0, 0.0, 0.0 ), lambda, 0.0 );

    bcWall->findBoundaryCells( meshAdapter, [&](Vec3 center){ 
        return center.z < 0.5;
    } );
    
    //dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPZ );
    dataBase->boundaryConditions.push_back( bcWall );

    //////////////////////////////////////////////////////////////////////////

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables {
        
        real radius = cellCenter.length();

        return toConservedVariables( PrimitiveVariables( 1.0, 0.0, 0.0, 0.0, lambda, 0.0 ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    writeVtkXML( dataBase, parameters, 0, path + "grid/Test_0" );

    //////////////////////////////////////////////////////////////////////////

    for( uint iter = 1; iter < 100000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( iter % 1000 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + "grid/Test_" + std::to_string( iter ) );
        }
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );


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
