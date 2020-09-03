
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ||          ||  ||  ||||||  |||||||| ||    ||  ||||||||  ||
//    ||        ||   ||  ||   ||    ||    ||    ||  ||    ||  ||
//     ||      ||    ||  ||||||     ||    ||    ||  ||||||||  ||
//      ||    ||     ||  ||   ||    ||     ||||||   ||    ||  ||||||    ||||||   ||   ||||||   ||||||   ||||||
//       ||  ||                                                        ||       ||   ||   ||  ||      |||    ||
//        ||||       |||||||||||||||||||||||||||||||||||||||||||||||||||||||   ||   ||||||   ||||||     |||
//                                                                    ||      ||   ||   ||  ||       ||   |||
//                    i R M B  @  T U  B r a u n s c h w e i g       ||      ||   ||   ||  ||||||   |||||||
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>

#include "Core/Timer/Timer.h"
#include "PointerDefinitions.h"
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

#include "GridGenerator/utilities/communication.h"

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
#include "GksGpu/BoundaryConditions/ConcreteHeatFlux.h"
#include "GksGpu/BoundaryConditions/Open.h"

#include "GksGpu/Communication/Communicator.h"
#include "GksGpu/Communication/MpiUtility.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"
#include "GksGpu/Analyzer/PointTimeseriesCollector.h"

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = 0.1;

    real L = 1.0;

    real Pr  = 1.0;
    real K   = 2.0;
    
    real g   = 9.81;
    real rho = 1.0;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );
    setLambdaFromT( prim, 12.0 );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real mu  = 1.0e-2;

    real cp = 0.5 * ( K + 5 ) * R_U / M_A;

    real k = mu / Pr * cp;

    real dt  = 0.000025;

    //real dt = 0.01 * dx*dx / ( 2.0 * mu );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";

    //*logging::out << logging::Logger::INFO_HIGH << "HRR = " << U * rhoFuel * LBurner * LBurner * (heatOfReaction * 100.0) / 0.016 / 1000.0 << " kW\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.D = mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = prim.lambda;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L,  
                                0.5 * L,  0.5 * L,  0.5 * L, dx);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "Grid_rank_" + std::to_string( rank ) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "MeshFaces_rank_" + std::to_string( rank ) + ".vtk" );

    //meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcWall1 = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), prim.lambda, false);

    bcWall1->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    
    SPtr<BoundaryCondition> bcWall2 = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false);

    bcWall2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x > -0.5*L; } );

    //SPtr<BoundaryCondition> bcWallHeatFlux = std::make_shared<ConcreteHeatFlux>( dataBase, 9,  0.1 * k / 1.0 / 50.0, 1.0, 50.0, 0.1, 3.0 );
    //SPtr<BoundaryCondition> bcWallHeatFlux = std::make_shared<ConcreteHeatFlux>( dataBase, 9,  1.0 * k / 1.0 / 50.0, 1.0, 50.0, 0.1, 3.0 );
    SPtr<BoundaryCondition> bcWallHeatFlux = std::make_shared<ConcreteHeatFlux>( dataBase, 9, 10.0 * k / 1.0 / 50.0, 1.0, 50.0, 0.1, 3.0 );

    bcWallHeatFlux->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x > 0.5*L && std::abs(center.y) < 0.5*L && std::abs(center.z) < 0.5*L; } );

    std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->init();

    //////////////////////////////////////////////////////////////////////////


    dataBase->boundaryConditions.push_back( bcWallHeatFlux );

    dataBase->boundaryConditions.push_back( bcWall1 );
    dataBase->boundaryConditions.push_back( bcWall2 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    dataBase->setCommunicators( meshAdapter );

    CudaUtility::printCudaMemoryUsage();
    
    Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {

        PrimitiveVariables primLocal = prim;

        return toConservedVariables(primLocal, parameters.K);
    });

    writeVtkXML(dataBase, parameters, 0, path + simulationName + "_0");

    //std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->writeVTKFile(dataBase, parameters, path + simulationName + "_Solid_0");
    writeConcreteHeatFluxVtkXML( dataBase, std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux), parameters, 0, path + simulationName + "_Solid_0" );

    //////////////////////////////////////////////////////////////////////////

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

    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "==================   S t a r t    T i m e    S t e p p i n g   =================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 10000000; iter++ )
    {
        //////////////////////////////////////////////////////////////////////////

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        //////////////////////////////////////////////////////////////////////////

        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        //////////////////////////////////////////////////////////////////////////

        int crashCellIndex = dataBase->getCrashCellIndex();
        if( crashCellIndex >= 0 )
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Simulation Crashed at CellIndex = " << crashCellIndex << "\n";
            dataBase->copyDataDeviceToHost();
            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            //std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->writeVTKFile(dataBase, parameters, path + simulationName + "_Solid_" + std::to_string( iter ));
            writeConcreteHeatFluxVtkXML( dataBase, std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux), parameters, 0, path + simulationName + "_Solid_" + std::to_string( iter ) );

            break;
        }

        if( iter % 100000 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            //std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->writeVTKFile(dataBase, parameters, path + simulationName + "_Solid_" + std::to_string( iter ));
            writeConcreteHeatFluxVtkXML( dataBase, std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux), parameters, 0, path + simulationName + "_Solid_" + std::to_string( iter ) );
        }
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );

    //turbulenceAnalyzer->download();

    //writeTurbulenceVtkXML(dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence");
}

int main( int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/ConcreteHeatFluxBCTest/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "ConcreteHeatFluxBCTest" );

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    //////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice( 0 );

    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    //////////////////////////////////////////////////////////////////////////

    try
    {
        thermalCavity( path, simulationName );
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
