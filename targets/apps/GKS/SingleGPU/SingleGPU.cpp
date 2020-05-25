//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <exception>
#include <fstream>
#include <sstream>
#include <memory>
#include <omp.h>

#include "Core/Timer/Timer.h"
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/VerticalCylinder/VerticalCylinder.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"
#include "GridGenerator/geometries/BoundingBox/BoundingBox.h"
#include "GridGenerator/utilities/communication.h"

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

#include "GksGpu/Communication/Communicator.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"
#include "GksGpu/Communication/MpiUtility.h"

//////////////////////////////////////////////////////////////////////////

real performanceTest( std::string path, std::string simulationName, uint nx )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L  = 1.0;

    real LX = L;
    real LY = L;
    real LZ = L;

    real dx = L / real(nx);

    //////////////////////////////////////////////////////////////////////////

    GksGpu::Parameters parameters;

    parameters.K  = 0;
    parameters.Pr = 1;
    parameters.mu = 0.01;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = 0.0001 * ( double(128) / double(nx) );
    parameters.dx = dx;

    parameters.lambdaRef = 1.0e-2;
    
    parameters.forcingSchemeIdx = 2;

    parameters.enableReaction = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    gridBuilder->addCoarseGrid( - 0.5*L, - 0.5*L, - 0.5*L,
                                  0.5*L,   0.5*L,   0.5*L, dx);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(true,true,true);

    gridBuilder->buildGrids(GKS, false);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto dataBase = std::make_shared<GksGpu::DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<GksGpu::BoundaryCondition> bcMX = std::make_shared<GksGpu::Periodic>( dataBase );
    SPtr<GksGpu::BoundaryCondition> bcPX = std::make_shared<GksGpu::Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<GksGpu::BoundaryCondition> bcMY = std::make_shared<GksGpu::Periodic>( dataBase );
    SPtr<GksGpu::BoundaryCondition> bcPY = std::make_shared<GksGpu::Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<GksGpu::BoundaryCondition> bcMZ = std::make_shared<GksGpu::Periodic>( dataBase );
    SPtr<GksGpu::BoundaryCondition> bcPZ = std::make_shared<GksGpu::Periodic>( dataBase );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*L; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_HIGH << "NumberOfBoundaryConditions = " << (int)dataBase->boundaryConditions.size() << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "bcMX ==> " << bcMX->numberOfCellsPerLevel[0] << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "bcPX ==> " << bcPX->numberOfCellsPerLevel[0] << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "bcMY ==> " << bcMY->numberOfCellsPerLevel[0] << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "bcPY ==> " << bcPY->numberOfCellsPerLevel[0] << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "bcMZ ==> " << bcMZ->numberOfCellsPerLevel[0] << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "bcPZ ==> " << bcPZ->numberOfCellsPerLevel[0] << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    dataBase->setCommunicators( meshAdapter );

    GksGpu::CudaUtility::printCudaMemoryUsage();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksGpu::Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> GksGpu::ConservedVariables
    {
        real U = 0.1;

        real ULocal =   0.1 + U * sin( 2.0 * M_PI * cellCenter.x ) * cos( 2.0 * M_PI * cellCenter.y ) * cos( 2.0 * M_PI * cellCenter.z );
        real VLocal =   0.1 - U * cos( 2.0 * M_PI * cellCenter.x ) * sin( 2.0 * M_PI * cellCenter.y ) * cos( 2.0 * M_PI * cellCenter.z );
        real WLocal =   0.1;

        real rho = 1.0;

        real p0 = 0.5 * rho / parameters.lambdaRef;

        real pLocal = p0 + rho * U * U / 16.0 * ( cos( 2.0 * M_PI * 2.0 * cellCenter.x ) + cos( 2.0 * M_PI * 2.0 * cellCenter.y ) ) * ( 2.0 + cos( 2.0 * M_PI * 2.0 * cellCenter.z ) );

        real rhoLocal = 2.0 * pLocal * parameters.lambdaRef;

        return GksGpu::toConservedVariables( GksGpu::PrimitiveVariables( rhoLocal, ULocal, VLocal, WLocal, parameters.lambdaRef ), parameters.K );
    });

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, level );

    GksGpu::Initializer::initializeDataUpdate(dataBase);

    //dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const uint numberOfIterations = 1000;

    GksGpu::CupsAnalyzer cupsAnalyzer( dataBase, false, 30.0, true, numberOfIterations );

    real CUPS = 0;

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= numberOfIterations; iter++ )
    {
        GksGpu::TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        CUPS = cupsAnalyzer.run( iter, parameters.dt );
    }

    //////////////////////////////////////////////////////////////////////////

    //dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + simulationName + "_final" );
    
    //////////////////////////////////////////////////////////////////////////

    int crashCellIndex = dataBase->getCrashCellIndex();
    if( crashCellIndex >= 0 )
    {
        *logging::out << logging::Logger::LOGGER_ERROR << "=================================================\n";
        *logging::out << logging::Logger::LOGGER_ERROR << "=================================================\n";
        *logging::out << logging::Logger::LOGGER_ERROR << "============= Simulation Crashed!!! =============\n";
        *logging::out << logging::Logger::LOGGER_ERROR << "=================================================\n";
        *logging::out << logging::Logger::LOGGER_ERROR << "=================================================\n";
    }

    return CUPS;
}

int main( int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/SingleGPU/" );
#else
    //std::string path( "/home/stephan/Computations/out/" );
    std::string path( "out/" );
#endif

    //////////////////////////////////////////////////////////////////////////

    try
    {
        logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
        logging::Logger::timeStamp(logging::Logger::ENABLE);

        std::string simulationName ( "SingleGPU" );

        std::ofstream file;
        file.open( path + simulationName + ".dat" );

        //std::vector<uint> nxList = {32,64,128,256};
        std::vector<uint> nxList = {128};

        for( auto nx : nxList )
        {
            logging::Logger::addStream(&std::cout);
    
            std::ofstream logFile( path + simulationName + "_nx_" + std::to_string(nx) + ".log" );
            logging::Logger::addStream(&logFile);

            GksGpu::CudaUtility::setCudaDevice( 0 );
    
            //////////////////////////////////////////////////////////////////////////

            if( sizeof(real) == 4 )
                *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
            else
                *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

            real CUPS = performanceTest( path, simulationName + "_nx_" + std::to_string(nx), nx );

            file << std::setw(5) << nx <<std::setw(20) << CUPS << std::endl;

            logFile.close();
            
            logging::Logger::resetStreams();
        }

        file.close();
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

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

   return 0;
}
