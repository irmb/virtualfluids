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

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName, uint restartIter )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = 0.05;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real LBurner = 1.0;

    real HBurner = 0.5;

    real Pr  = 0.71;
    real K   = 5.0;
    
    real g   = 9.81;
    real rho = 1.2;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );
    setLambdaFromT( prim, 3.0 );
    
    real mu = 1.5e-4;

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );
    real U   = 0.0125;

    real CFL = 0.125;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

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

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(-2.1, -3.5, -0.1,  
                                2.1,  3.5,  3.1, dx);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    //TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/RoomExtended.stl");
#else
    //TriangularMesh* stl = TriangularMesh::make("inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("inp/RoomExtended.stl");
#endif

    //gridBuilder->addGeometry(stl);
    
    //Cuboid box( -0.5 * LBurner, -0.5 * LBurner, -HBurner, 
    //             0.5 * LBurner,  0.5 * LBurner,  HBurner );
    //Cuboid beam( -0.15, -10.0, 2.6, 0.15, 10.0, 10.0 );

    //Conglomerate solid;

    //solid.add(&box);
    //solid.add(&beam);

    gridBuilder->addGeometry(stl);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Cuboid boxRefCoarse1 ( -0.8 * LBurner, -0.8 * LBurner, -100.0, 
    //                        0.8 * LBurner,  0.8 * LBurner,  100.0 );
    //Cuboid boxRefCoarse2 ( -0.8 * LBurner, -100,    2.3, 
    //                        0.8 * LBurner,  100,  100.0 );

    //Conglomerate refRegionCoarse;

    //refRegionCoarse.add( &boxRefCoarse1 );
    //refRegionCoarse.add( &boxRefCoarse2 );

    //gridBuilder->setNumberOfLayers(0,20);

    //gridBuilder->addGrid( &refRegionCoarse, 1 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cuboid boxRef ( -0.5 * LBurner, -0.5 * LBurner, -HBurner, 
                     0.5 * LBurner,  0.5 * LBurner,  HBurner );
    Cuboid beamRef( -10.0, -0.15, 2.6, 10.0, 0.15, 10.0 );

    boxRef.scale (0.2);
    beamRef.scale(0.2);

    Conglomerate refRegion1;

    refRegion1.add( &boxRef );
    refRegion1.add( &beamRef );

    gridBuilder->setNumberOfLayers(0,20);

    gridBuilder->addGrid( &refRegion1, 2 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    gridBuilder->writeGridsToVtk(path + "Grid_lev_");

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
    
    //SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, false );

    //bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    //bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMY;
    //SPtr<BoundaryCondition> bcPY;

    //bcMY = std::make_shared<AdiabaticWall>(dataBase, Vec3(0.0, 0.0, 0.0), false);
    //bcPY = std::make_shared<AdiabaticWall>(dataBase, Vec3(0.0, 0.0, 0.0), false);

    //bcMY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y < -0.5*L; });
    //bcPY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y >  0.5*L; });

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, true );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, true );
    
    //bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < 0.5; } );
    //bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > H  ; } );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcBurner = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), 0.5*prim.lambda,  0.0, true );
    //SPtr<BoundaryCondition> bcBurner = std::make_shared<HeatFlux>( dataBase, 100.0 );
    //SPtr<BoundaryCondition> bcBurner = std::make_shared<CreepingMassFlux>( dataBase, rho, U, prim.lambda );

    //bcBurner->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

    //    return center.z > HBurner - 0.125 * dx && center.z < HBurner && std::sqrt(center.x*center.x) < 0.5 * LBurner - dx && std::sqrt(center.y*center.y) < 0.5 * LBurner - dx;
    //} );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcSolid = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );

    //bcSolid->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

    //    return center.z > 2.5 && std::sqrt(center.x*center.x) < 0.15;
    //} );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcWindowOpen = std::make_shared<Open>( dataBase, prim, 1.0 );

    //bcWindowOpen->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

    //    return center.z > 1.0 && center.z < 2.0 && std::sqrt(center.x*center.x) > 1.5 && std::sqrt(center.y*center.y) < 1.0;
    //} );

    //SPtr<BoundaryCondition> bcWindowPressure = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );

    //bcWindowPressure->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

    //    return center.z > 2.0 && center.z < 2.8 && std::sqrt(center.x*center.x) > 1.5 && std::sqrt(center.y*center.y) < 1.0;
    //} );

    //////////////////////////////////////////////////////////////////////////

    //dataBase->boundaryConditions.push_back( bcBurner );

    //dataBase->boundaryConditions.push_back( bcMX );
    //dataBase->boundaryConditions.push_back( bcPX );
    //
    //dataBase->boundaryConditions.push_back( bcMY );
    //dataBase->boundaryConditions.push_back( bcPY );

    //dataBase->boundaryConditions.push_back( bcMZ );
    //dataBase->boundaryConditions.push_back( bcPZ );

    //dataBase->boundaryConditions.push_back( bcSolid );

    //dataBase->boundaryConditions.push_back( bcWindowOpen     );
    //dataBase->boundaryConditions.push_back( bcWindowPressure );

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

            //real rhoLocal = rho * std::exp(-(2.0 * g * H * prim.lambda) * cellCenter.z / H);

            //prim.rho = rhoLocal;

            //real r = sqrt(cellCenter.x * cellCenter.x + cellCenter.y * cellCenter.y + cellCenter.z * cellCenter.z);

            //if( r < 0.55 ) prim.S_2 = 1.0;

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

    return;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0, true, 1000 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 1000 );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 100000000; iter++ )
    {
        if( iter < 10000 )
        {
            //std::dynamic_pointer_cast<PassiveScalarDiriclet>(burner)->S_1 = 10.0 * ( real(iter) / 20000.0 );

            //parameters.mu = mu + 10.0 * mu * ( 1.0 - ( real(iter) / 10000.0 ) );

            //parameters.dt = 0.2 * dt + ( dt - 0.2 * dt ) * ( real(iter) / 40000.0 );
        }

        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( 
            //( iter >= 34920 && iter % 1 == 0 ) ||
            //( iter >= 35900 && iter % 10 == 0 ) ||
            ( iter % 1000 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        if( iter % 1000 == 0 )
        {
            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
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
