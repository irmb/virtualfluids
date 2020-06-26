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
#include "GksGpu/Analyzer/PointTimeSeriesAnalyzer.h"

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName, uint _gpuIndex, uint _testIndex, uint _nx, uint restartIter )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //uint nx = 64;
    uint nx = _nx;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 3.0;
    real H = 4.0;

    real R = 0.5;

    real dx = H / real(nx);

    real Pr  = 0.71;
    real K   = 2.0;
    
    real g   = 9.81;
    real rho = 1.2;

    real mu = 1.8e-5;

    real U = 0.0;
    real rhoFuel = 0.0;

    if( _testIndex == 14 ) { U = 0.074; rhoFuel = 0.5405; }    // Test 14      low flow rates
    if( _testIndex == 24 ) { U = 0.097; rhoFuel = 0.5464; }    // Test 24      medium flow rate
    if( _testIndex == 17 ) { U = 0.117; rhoFuel = 0.5641; }    // Test 17      high flow rate

    GksGpu::PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );

    GksGpu::setLambdaFromT( prim, 2.85 );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    //real CFL = 0.06125;
    real CFL = 0.125;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( c1o1 + ( c2o1 * mu ) / ( U * dx * rho ) ) ) );

    real dh = 8000.0; // kJ / kmol  / T_FAKTOR

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";
    *logging::out << logging::Logger::INFO_HIGH << "Pr = " << Pr << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "HRR = " << U * rhoFuel * M_PI * R * R * ( dh * 100 ) / 0.016 / 1000.0 << " kW\n";

    //////////////////////////////////////////////////////////////////////////

    GksGpu::Parameters parameters;

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

    parameters.viscosityModel = GksGpu::ViscosityModel::sutherlandsLaw;
    //parameters.viscosityModel = GksGpu::ViscosityModel::constant;

    parameters.enableReaction = true;

    parameters.useHeatReleaseRateLimiter = true;
    parameters.useReactionLimiter        = true;
    parameters.useTemperatureLimiter     = true;
    parameters.usePassiveScalarLimiter   = true;
    parameters.useSmagorinsky            = true;

    parameters.reactionLimiter    = 1.0005;
    parameters.temperatureLimiter = 1.0e-3;

    parameters.useSpongeLayer = true;
    parameters.spongeLayerIdx = 1;

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
    
    //gridBuilder->addGrid( &cylinder1 );
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

    GksGpu::CudaUtility::setCudaDevice(_gpuIndex);

    auto dataBase = std::make_shared<GksGpu::DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real openBoundaryVelocityLimiter = 1.0;

    SPtr<GksGpu::BoundaryCondition> bcMX = std::make_shared<GksGpu::Open>( dataBase, prim, openBoundaryVelocityLimiter );
    SPtr<GksGpu::BoundaryCondition> bcPX = std::make_shared<GksGpu::Open>( dataBase, prim, openBoundaryVelocityLimiter );

    SPtr<GksGpu::BoundaryCondition> bcMX_2 = std::make_shared<GksGpu::Symmetry>( dataBase, 'x' );
    SPtr<GksGpu::BoundaryCondition> bcPX_2 = std::make_shared<GksGpu::Symmetry>( dataBase, 'x' );

    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    bcMX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L && center.z > 0.9*H; } );
    bcPX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L && center.z > 0.9*H; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<GksGpu::BoundaryCondition> bcMY;
    SPtr<GksGpu::BoundaryCondition> bcPY;

    SPtr<GksGpu::BoundaryCondition> bcMY_2;
    SPtr<GksGpu::BoundaryCondition> bcPY_2;

    if( threeDimensional )
    {
        bcMY = std::make_shared<GksGpu::Open>( dataBase, prim, openBoundaryVelocityLimiter );
        bcPY = std::make_shared<GksGpu::Open>( dataBase, prim, openBoundaryVelocityLimiter );

        bcMY_2 = std::make_shared<GksGpu::Symmetry>( dataBase, 'y' );
        bcPY_2 = std::make_shared<GksGpu::Symmetry>( dataBase, 'y' );

        bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L; } );
        bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L; } );

        bcMY_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L && center.z > 0.9*H; } );
        bcPY_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L && center.z > 0.9*H; } );
    }
    else
    {
        bcMY = std::make_shared<GksGpu::Periodic>(dataBase);
        bcPY = std::make_shared<GksGpu::Periodic>(dataBase);

        bcMY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y < -0.5*dx; });
        bcPY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y >  0.5*dx; });
    }

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<GksGpu::BoundaryCondition> bcMZ = std::make_shared<GksGpu::AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<InflowComplete>( dataBase, PrimitiveVariables(rho, 0.0, 0.0, 0.0, prim.lambda, 0.0, 0.0) );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<Open>( dataBase );

    SPtr<GksGpu::BoundaryCondition> bcPZ = std::make_shared<GksGpu::Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < 0.0; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > H  ; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<GksGpu::BoundaryCondition> burner = std::make_shared<GksGpu::CreepingMassFlux>( dataBase, rhoFuel, U, prim.lambda );
    //SPtr<BoundaryCondition> burner = std::make_shared<Inflow>( dataBase, Vec3(0,0,U), prim.lambda, rhoFuel, 1, 0, 0, 1.0 );

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

    //dataBase->boundaryConditions.push_back( burner );

    //////////////////////////////////////////////////////////////////////////

    auto pointTimeSeriesAnalyzerU_P1 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 0.505), 'U' );
    auto pointTimeSeriesAnalyzerV_P1 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 0.505), 'V' );
    auto pointTimeSeriesAnalyzerW_P1 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 0.505), 'W' );

    auto pointTimeSeriesAnalyzerU_P2 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 0.505), 'U' );
    auto pointTimeSeriesAnalyzerV_P2 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 0.505), 'V' );
    auto pointTimeSeriesAnalyzerW_P2 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 0.505), 'W' );

    auto pointTimeSeriesAnalyzerU_P3 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 2.0), 'U' );
    auto pointTimeSeriesAnalyzerV_P3 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 2.0), 'V' );
    auto pointTimeSeriesAnalyzerW_P3 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.0, 0.0, 2.0), 'W' );

    auto pointTimeSeriesAnalyzerU_P4 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 2.0), 'U' );
    auto pointTimeSeriesAnalyzerV_P4 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 2.0), 'V' );
    auto pointTimeSeriesAnalyzerW_P4 = std::make_shared<GksGpu::PointTimeSeriesAnalyzer>( dataBase, meshAdapter, Vec3(0.5, 0.0, 2.0), 'W' );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    GksGpu::CudaUtility::printCudaMemoryUsage();
    
    if( restartIter == INVALID_INDEX )
    {
        GksGpu::Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> GksGpu::ConservedVariables {

            GksGpu::PrimitiveVariables primLocal = prim;

            return GksGpu::toConservedVariables(primLocal, parameters.K);
        });

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );
    }
    else
    {
        GksGpu::Restart::readRestart( dataBase, path + simulationName + "_" + std::to_string( restartIter ), startIter );

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( restartIter ) + "_restart" );
    }

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, level );

    GksGpu::Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataDeviceToHost();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint iterPerSecond = uint( c1o1 / parameters.dt ) + 1;

    *logging::out << logging::Logger::INFO_HIGH << "iterPerSecond = " << iterPerSecond << "\n";

    //////////////////////////////////////////////////////////////////////////

    GksGpu::CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0, true, 10000 );

    GksGpu::ConvergenceAnalyzer convergenceAnalyzer( dataBase, 10000 );

    auto turbulenceAnalyzer = std::make_shared<GksGpu::TurbulenceAnalyzer>( dataBase, 10 * iterPerSecond );

    turbulenceAnalyzer->collect_UU = true;
    turbulenceAnalyzer->collect_VV = true;
    turbulenceAnalyzer->collect_WW = true;

    turbulenceAnalyzer->allocate();

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();
    
    for( uint iter = startIter + 1; iter <= 40 * iterPerSecond; iter++ )
    {
        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        GksGpu::TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        pointTimeSeriesAnalyzerU_P1->run(iter, parameters);
        pointTimeSeriesAnalyzerV_P1->run(iter, parameters);
        pointTimeSeriesAnalyzerW_P1->run(iter, parameters);

        pointTimeSeriesAnalyzerU_P2->run(iter, parameters);
        pointTimeSeriesAnalyzerV_P2->run(iter, parameters);
        pointTimeSeriesAnalyzerW_P2->run(iter, parameters);

        pointTimeSeriesAnalyzerU_P3->run(iter, parameters);
        pointTimeSeriesAnalyzerV_P3->run(iter, parameters);
        pointTimeSeriesAnalyzerW_P3->run(iter, parameters);

        pointTimeSeriesAnalyzerU_P4->run(iter, parameters);
        pointTimeSeriesAnalyzerV_P4->run(iter, parameters);
        pointTimeSeriesAnalyzerW_P4->run(iter, parameters);

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
            GksGpu::Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
        }

        if( iter % 100000 == 0 )
        {
            turbulenceAnalyzer->download();

            writeTurbulenceVtkXML( dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence_" + std::to_string( iter ) );
        }

        if( iter % 100000 == 0 )
        {
            pointTimeSeriesAnalyzerU_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerU_P1->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerV_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerV_P1->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerW_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerW_P1->quantity + "_" + std::to_string( iter ));

            pointTimeSeriesAnalyzerU_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerU_P2->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerV_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerV_P2->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerW_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerW_P2->quantity + "_" + std::to_string( iter ));

            pointTimeSeriesAnalyzerU_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerU_P3->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerV_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerV_P3->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerW_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerW_P3->quantity + "_" + std::to_string( iter ));

            pointTimeSeriesAnalyzerU_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerU_P4->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerV_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerV_P4->quantity + "_" + std::to_string( iter ));
            pointTimeSeriesAnalyzerW_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerW_P4->quantity + "_" + std::to_string( iter ));
        }

        turbulenceAnalyzer->run( iter, parameters );
    }

    //////////////////////////////////////////////////////////////////////////

    //dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );

    turbulenceAnalyzer->download();

    writeTurbulenceVtkXML(dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence_final");

    pointTimeSeriesAnalyzerU_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerU_P1->quantity + "_final");
    pointTimeSeriesAnalyzerV_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerV_P1->quantity + "_final");
    pointTimeSeriesAnalyzerW_P1->writeToFile(path + simulationName + "_P1_TimeSeries_" + pointTimeSeriesAnalyzerW_P1->quantity + "_final");

    pointTimeSeriesAnalyzerU_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerU_P2->quantity + "_final");
    pointTimeSeriesAnalyzerV_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerV_P2->quantity + "_final");
    pointTimeSeriesAnalyzerW_P2->writeToFile(path + simulationName + "_P2_TimeSeries_" + pointTimeSeriesAnalyzerW_P2->quantity + "_final");

    pointTimeSeriesAnalyzerU_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerU_P3->quantity + "_final");
    pointTimeSeriesAnalyzerV_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerV_P3->quantity + "_final");
    pointTimeSeriesAnalyzerW_P3->writeToFile(path + simulationName + "_P3_TimeSeries_" + pointTimeSeriesAnalyzerW_P3->quantity + "_final");

    pointTimeSeriesAnalyzerU_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerU_P4->quantity + "_final");
    pointTimeSeriesAnalyzerV_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerV_P4->quantity + "_final");
    pointTimeSeriesAnalyzerW_P4->writeToFile(path + simulationName + "_P4_TimeSeries_" + pointTimeSeriesAnalyzerW_P4->quantity + "_final");
}

int main( int argc, char* argv[])
{
    uint restartIter = INVALID_INDEX;
    //uint restartIter = 90000;
        
    uint gpuIndex = 1;
    uint testIndex = 24;
    uint nx = 64;

    if( argc > 1 ) gpuIndex    = atoi( argv[1] );

    if( argc > 2 ) testIndex   = atoi( argv[2] );

    if( argc > 3 ) nx          = atoi( argv[3] );

    if( argc > 4 ) restartIter = atoi( argv[4] );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/SandiaFlame_1m/" );
#else
    std::string path( "out/" );
#endif

    path += "Test_" + std::to_string(testIndex) + "_" + std::to_string(nx) + "/";

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
        thermalCavity( path, simulationName, gpuIndex, testIndex, nx, restartIter );
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