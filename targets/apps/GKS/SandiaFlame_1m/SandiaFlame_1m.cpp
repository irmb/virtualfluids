//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <exception>
#include <fstream>
#include <memory>

#include "Core/Timer/Timer.h"
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"
#include "Core/buildInfo.h"

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/VerticalCylinder/VerticalCylinder.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"

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

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"
#include "GksGpu/BoundaryConditions/Pressure.h"
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"
#include "GksGpu/BoundaryConditions/PassiveScalarDiriclet.h"
#include "GksGpu/BoundaryConditions/InflowComplete.h"
#include "GksGpu/BoundaryConditions/Open.h"
#include "GksGpu/BoundaryConditions/Inflow.h"
#include "GksGpu/BoundaryConditions/Symmetry.h"
#include "GksGpu/BoundaryConditions/Pressure2.h"
#include "GksGpu/BoundaryConditions/CreepingMassFlux.h"

#include "GksGpu/Interface/Interface.h"
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

    uint nx = 128;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 3.0;
    real H = 4.0;

    real R = 0.5;

    real dx = H / real(nx);

    real U = 0.074;

    real Pr  = 0.71;
    real K   = 2.0;
    
    real g   = 9.81;
    real rho = 1.2;
    real rhoFuel = 0.5405;

    real mu = 1.5e-5;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );

    setLambdaFromT( prim, 2.85 / T_FAKTOR );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real CFL = 0.06125;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    //real dh = 4192.0; // kJ / kmol  / T_FAKTOR
    real dh = 8000.0; // kJ / kmol  / T_FAKTOR

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";
    *logging::out << logging::Logger::INFO_HIGH << "Pr = " << Pr << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "HRR = " << U * rhoFuel * M_PI * R * R * ( dh * 100 ) / 0.016 / 1000.0 << " kW\n";

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

    parameters.rhoRef    = rho;

    parameters.heatOfReaction = dh;

    parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    //parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    parameters.useReactionLimiter      = false;
    parameters.useTemperatureLimiter   = true;
    parameters.usePassiveScalarLimiter = true;
    parameters.useSmagorinsky          = true;

    parameters.reactionLimiter = 1.0005;

    parameters.useSpongeLayer = false;
    parameters.spongeLayerIdx = 0;

    parameters.forcingSchemeIdx = 2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool threeDimensional = true;

    if( threeDimensional )
    {
        gridBuilder->addCoarseGrid(-0.5*L, -0.5*L, 0.0,
                                    0.5*L,  0.5*L, H, dx);
    }
    else
    {
        gridBuilder->addCoarseGrid(-0.5*L, -0.5*dx, 0.0,
                                    0.5*L,  0.5*dx, H, dx);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VerticalCylinder cylinder1( 0.0, 0.0, 0.0, 1.5*R, 0.25*H );
    VerticalCylinder cylinder2( 0.0, 0.0, 0.0, 1.1*R, 0.05*H );
    
    Conglomerate refRing;
    refRing.add     ( new VerticalCylinder( 0.0, 0.0, 0.0, 1.2*R, 0.1 ) );
    refRing.subtract( new VerticalCylinder( 0.0, 0.0, 0.0, 0.8*R, 1.0    ) );

    gridBuilder->setNumberOfLayers(0,10);
    
    gridBuilder->addGrid( &cylinder1 );
    //gridBuilder->addGrid( &cylinder2 );
    //gridBuilder->addGrid( &refRing );

    if( threeDimensional ) gridBuilder->setPeriodicBoundaryCondition(false, false, false);
    else                   gridBuilder->setPeriodicBoundaryCondition(false, true,  false);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    if( !threeDimensional )
        meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(1);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real openBoundaryVelocityLimiter = 1.0;

    SPtr<BoundaryCondition> bcMX = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );

    SPtr<BoundaryCondition> bcMX_2 = std::make_shared<Symmetry>( dataBase, 'x' );
    SPtr<BoundaryCondition> bcPX_2 = std::make_shared<Symmetry>( dataBase, 'x' );

    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    bcMX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L && center.z > 0.9*H; } );
    bcPX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L && center.z > 0.9*H; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMY;
    SPtr<BoundaryCondition> bcPY;

    SPtr<BoundaryCondition> bcMY_2;
    SPtr<BoundaryCondition> bcPY_2;

    if( threeDimensional )
    {
        bcMY = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
        bcPY = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );

        bcMY_2 = std::make_shared<Symmetry>( dataBase, 'y' );
        bcPY_2 = std::make_shared<Symmetry>( dataBase, 'y' );

        bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L; } );
        bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L; } );

        bcMY_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L && center.z > 0.9*H; } );
        bcPY_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L && center.z > 0.9*H; } );
    }
    else
    {
        bcMY = std::make_shared<Periodic>(dataBase);
        bcPY = std::make_shared<Periodic>(dataBase);

        bcMY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y < -0.5*dx; });
        bcPY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y >  0.5*dx; });
    }

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<InflowComplete>( dataBase, PrimitiveVariables(rho, 0.0, 0.0, 0.0, prim.lambda, 0.0, 0.0) );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<Open>( dataBase );

    SPtr<BoundaryCondition> bcPZ = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < 0.0; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > H  ; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> burner = std::make_shared<CreepingMassFlux>( dataBase, rhoFuel, U, prim.lambda );

    burner->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 
        
        if( threeDimensional )
            return center.z < 0.0 && std::sqrt(center.x*center.x + center.y*center.y) < R;
        else
            return center.z < 0.0 && std::sqrt(center.x*center.x) < R && std::sqrt(center.y*center.y) < 0.5 * dx;
    } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( burner );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX_2 );
    dataBase->boundaryConditions.push_back( bcPX_2 );

    if( threeDimensional ){
        dataBase->boundaryConditions.push_back( bcMY_2 );
        dataBase->boundaryConditions.push_back( bcPY_2 );
    }

    //////////////////////////////////////////////////////////////////////////

    auto pointTimeSeriesAnalyzer = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 0.5), 'W' );

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

            PrimitiveVariables primLocal = prim;

            return toConservedVariables(primLocal, parameters.K);
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

    auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 2000000; iter++ )
    {
        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        pointTimeSeriesAnalyzer->run(iter, parameters);

        int crashCellIndex = dataBase->getCrashCellIndex();

        if( crashCellIndex >= 0 )
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Simulation Crashed at CellIndex = " << crashCellIndex << "\n";
            dataBase->copyDataDeviceToHost();
            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            break;
        }

        if( 
            //( iter >= 39360 && iter % 1 == 0 ) || 
            ( iter % 10000 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();
            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        if( iter % 10000 == 0 /*|| iter == 39000*/)
        {
            dataBase->copyDataDeviceToHost();
            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
        }

        if( iter % 100000 == 0 )
        {
            turbulenceAnalyzer->download();

            writeTurbulenceVtkXML( dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence_" + std::to_string( iter ) );
        }

        if( iter % 10000 == 0 )
        {
            pointTimeSeriesAnalyzer->writeToFile(path + simulationName + "_TimeSeries_" + pointTimeSeriesAnalyzer->quantity + "_" + std::to_string( iter ));
        }

        turbulenceAnalyzer->run( iter, parameters );
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
    std::string path( "F:/Work/Computations/out/SandiaFlame_1m/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "Flame" );

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        uint restartIter = INVALID_INDEX;
        //uint restartIter = 50000;

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
