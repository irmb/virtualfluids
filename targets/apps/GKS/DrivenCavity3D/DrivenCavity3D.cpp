//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>

#include "Core/Timer/Timer.h"
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

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void drivenCavity( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L / 64.0;

    real Re  = 1.0e2;
    real U  = 0.1;
    real Ma = 0.1;
    
    real Pr  = 1.0;
    real K   = 2.0;

    real rho = 1.0;

    real mu = U * rho * L / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( cs * cs );

    real CFL = 0.5;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

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

    //Cuboid refBox(-1.0, -1.0, 0.475, 1.0, 1.0, 0.55);
    ////Cuboid refBox(-1.0, -1.0, -1.0, 1.0, 1.0, -0.475);

    //gridBuilder->setNumberOfLayers(6,6);
    //gridBuilder->addGrid( &refBox, 1);

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + simulationName + "_Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + simulationName + "_MeshFaces.vtk" );

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcPZ   = std::make_shared<IsothermalWall>( dataBase, Vec3( U  , U  , 0.0 ), lambda, 0.0, false );
    SPtr<BoundaryCondition> bcWall = std::make_shared<IsothermalWall>( dataBase, Vec3( 0.0, 0.0, 0.0 ), lambda, 0.0, false );

    bcPZ->findBoundaryCells  ( meshAdapter, true,  [&](Vec3 center){ return center.z > 0.5; } );
    bcWall->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.z < 0.5; } );
    
    //dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPZ );
    dataBase->boundaryConditions.push_back( bcWall );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    //CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables {

        //real uLocal = U * ( cellCenter.z + 0.5 );

        //if( cellCenter.y )

        real uLocal = 0.0;

        return toConservedVariables( PrimitiveVariables( 1.0, uLocal, 0.0, 0.0, lambda, 0.0 ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );

    //////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, false, 60.0, true, 1000 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 1000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for( uint iter = 1; iter <= 10000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        //if( iter % 1000 == 0 )
        //{
        //    dataBase->copyDataDeviceToHost();

        //    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        //}

        cupsAnalyzer.run( iter );

        //convergenceAnalyzer.run( iter );
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );


}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/out/" );
    //std::string path( "out/" );
    std::string simulationName ( "DrivenCavity" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    
    try
    {
        drivenCavity( path, simulationName );
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
