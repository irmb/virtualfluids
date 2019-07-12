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
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/VerticalCylinder/VerticalCylinder.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "GksVtkAdapter/VTKInterface.h"

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"
#include "GksGpu/Initializer/Initializer.h"

#include "GksGpu/FlowStateData/FlowStateData.cuh"
#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"
#include "GksGpu/FlowStateData/ThermalDependencies.cuh"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"
#include "GksGpu/BoundaryConditions/Pressure2.h"
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"
#include "GksGpu/BoundaryConditions/HeatFlux.h"
#include "GksGpu/BoundaryConditions/CreepingMassFlux.h"
#include "GksGpu/BoundaryConditions/Open.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"
#include "GksGpu/Analyzer/PointTimeseriesAnalyzer.h"

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName, uint restartIter )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = 0.1;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real LBurner = 1.0;

    real HBurner = 0.5;

    real Pr  = 0.71;
    real K   = 5.0;
    
    real g   = 9.81;
    real rho = 1.2;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );
    setLambdaFromT( prim, 3.0 );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real mu      = 1.8e-5;
    real U       = 0.01;
    //real U       = 0.005;
    real rhoFuel = 0.5405;

    real CFL = 0.125;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( c1o1 + ( c2o1 * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";

    *logging::out << logging::Logger::INFO_HIGH << "HRR = " << U * rho * LBurner * LBurner * 800000.0 / 0.016 / 1000.0 << " kW\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.D = mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = -g;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = prim.lambda;

    parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    //parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    parameters.useHeatReleaseRateLimiter = true;
    parameters.useTemperatureLimiter     = true;
    parameters.usePassiveScalarLimiter   = true;
    parameters.useSmagorinsky            = true;

    parameters.reactionLimiter    = 1.0005;
    parameters.temperatureLimiter = 1.0e-7;

    parameters.useSpongeLayer = true;
    parameters.spongeLayerIdx = 2;

    parameters.forcingSchemeIdx = 2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(-2.1, -6.0, -0.1,  
                                2.1,  6.0,  3.1, dx);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    //TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/RoomExtended2.stl");
    //TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/RoomExtended3.stl");
#else
    //TriangularMesh* stl = TriangularMesh::make("inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("inp/RoomExtended.stl");
#endif

    gridBuilder->addGeometry(stl);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cuboid boxCoarse ( -3.0, -3.0, -0.5, 
                        3.0,  3.0,  3.5 );

    gridBuilder->addGrid( &boxCoarse, 1 );

    Cuboid boxRef ( -0.6 * LBurner, -0.6 * LBurner, -1.0, 
                     0.6 * LBurner,  0.6 * LBurner,  2.0 );
    Cuboid beamRef( -10.0, -0.15, 2.6, 10.0, 0.15, 10.0 );

    //boxRef.scale (0.5);
    //beamRef.scale(0.5);

    Conglomerate refRegion1;

    refRegion1.add( &boxRef );
    //refRegion1.add( &beamRef );

    gridBuilder->setNumberOfLayers(0,22);

    gridBuilder->addGrid( &refRegion1, 2 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    //meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(1);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcWall = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );

    bcWall->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return true; } );

    ////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcTop = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );

    bcTop->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > 3.0 || center.z < 0.0; } );

    ////////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcOpen = std::make_shared<Open>( dataBase, prim, 1.0 );

    bcOpen->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -6.0 || center.y > 2.2; } );

    ////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcPressure = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );

    bcPressure->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y > 2.2 && center.z > 3.0; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcBurner = std::make_shared<CreepingMassFlux>( dataBase, rho, U, prim.lambda );

    bcBurner->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

        return center.z < HBurner && std::sqrt(center.x*center.x) < 0.5 * LBurner - dx && std::sqrt(center.y*center.y) < 0.5 * LBurner - dx;
        //return center.z < 0.0 && std::sqrt(center.x*center.x) < 1.5 * LBurner - dx && std::sqrt(center.y*center.y) < 1.5 * LBurner - dx;
    } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcBurner );

    dataBase->boundaryConditions.push_back( bcWall );

    dataBase->boundaryConditions.push_back( bcTop );

    dataBase->boundaryConditions.push_back( bcOpen );

    dataBase->boundaryConditions.push_back( bcPressure );

    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcBurner = "   << bcBurner->numberOfCells   << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcWall = "     << bcWall->numberOfCells     << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcOpen = "     << bcOpen->numberOfCells     << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcPressure = " << bcPressure->numberOfCells << "\n";

    //////////////////////////////////////////////////////////////////////////

    auto pointTimeSeriesAnalyzer_P1 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(-1.5, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P2 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(-1.0, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P3 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(-0.5, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P4 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3( 0.0, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P5 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3( 0.5, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P6 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3( 1.0, 0.0, 2.5999), 'T' );
    auto pointTimeSeriesAnalyzer_P7 = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3( 1.5, 0.0, 2.5999), 'T' );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();
    
    if( restartIter == INVALID_INDEX )
    {
        Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {

            return toConservedVariables(prim, parameters.K);
        });

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );
    }
    else
    {
        Restart::readRestart( dataBase, path + simulationName + "_" + std::to_string( restartIter ), startIter );

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( restartIter ) + "_restart" );
    }

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, level );

    Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataDeviceToHost();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0, true, 10000 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 10000 );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 100000000; iter++ )
    {
        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        pointTimeSeriesAnalyzer_P1->run(iter, parameters);
        pointTimeSeriesAnalyzer_P2->run(iter, parameters);
        pointTimeSeriesAnalyzer_P3->run(iter, parameters);
        pointTimeSeriesAnalyzer_P4->run(iter, parameters);
        pointTimeSeriesAnalyzer_P5->run(iter, parameters);
        pointTimeSeriesAnalyzer_P6->run(iter, parameters);
        pointTimeSeriesAnalyzer_P7->run(iter, parameters);

        int crashCellIndex = dataBase->getCrashCellIndex();
        if( crashCellIndex >= 0 )
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Simulation Crashed at CellIndex = " << crashCellIndex << "\n";
            dataBase->copyDataDeviceToHost();
            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            break;
        }

        if( iter % 5000 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        if( iter % 5000 == 0 )
        {
            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
        }

        if( iter % 5000 == 0 )
        {
            pointTimeSeriesAnalyzer_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P5->writeToFile(path + simulationName + "_P5_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P6->writeToFile(path + simulationName + "_P6_TimeSeries_" + std::to_string( iter ));
            pointTimeSeriesAnalyzer_P7->writeToFile(path + simulationName + "_P7_TimeSeries_" + std::to_string( iter ));
        }

        //turbulenceAnalyzer->run( iter, parameters );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );

    //turbulenceAnalyzer->download();

    //writeTurbulenceVtkXML(dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence");
}

int main( int argc, char* argv[])
{

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/RoomFireExtended/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "RoomFire" );

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        uint restartIter = INVALID_INDEX;
        //uint restartIter = 35000;

        if( argc > 1 ) restartIter = atoi( argv[1] );

        thermalCavity( path, simulationName, restartIter );
    }
    catch (const std::exception& e)
    {     
        *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
    }
    catch (const std::bad_alloc& e)
    {  
        *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
    }
    catch (...)
    {
        *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
    }

    logFile.close();

    return 0;
}
