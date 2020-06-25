//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <sstream>
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

void performanceTest( std::string path, std::string simulationName, uint decompositionDimension, uint nx, bool strongScaling )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mpiWorldSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

    //CudaUtility::setCudaDevice(rank % devicesPerNode);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real H = 1.0;

    real L = 1.0;

    if( strongScaling ) L = H / double( mpiWorldSize );

    real dx = H / real(nx);

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = 0;
    parameters.Pr = 1;
    parameters.mu = 0.01;

    parameters.force.x = 0.1;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = 0.0001;
    parameters.dx = dx;

    parameters.lambdaRef = 1.0e-2;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( decompositionDimension == 1 && mpiWorldSize > 1 )
    {
        gridBuilder->addCoarseGrid( rank*L - 0.5*L - 5.0*dx, -0.5*H, -0.5*H,  
                                    rank*L + 0.5*L + 5.0*dx,  0.5*H,  0.5*H, dx);

        gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( rank*L - 0.5*L, rank*L + 0.5*L, 
                                                                         -H        ,      H,
                                                                         -H        ,      H ) );
    }else
    {
        gridBuilder->addCoarseGrid( -0.5*H, -0.5*H, -0.5*H,  
                                     0.5*H,  0.5*H,  0.5*H, dx);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //gridBuilder->setPeriodicBoundaryCondition(false, true, true);
    gridBuilder->setPeriodicBoundaryCondition(true, false, false);

    gridBuilder->buildGrids(GKS, false);

    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( decompositionDimension == 1 && mpiWorldSize > 1 )
    {
        gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
        gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, (rank + 1 + mpiWorldSize) % mpiWorldSize );

        gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
        gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, (rank - 1 + mpiWorldSize) % mpiWorldSize );
    }
    //if( decompositionDimension == 1 && mpiWorldSize > 1 && rank == 0 )
    //{
    //    gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
    //    gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, (rank + 1 + mpiWorldSize) % mpiWorldSize );
    //}
    //else
    //{
    //    gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
    //    gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, (rank - 1 + mpiWorldSize) % mpiWorldSize );
    //}

    //gridBuilder->writeGridsToVtk(path + "/Grid_rank_" + std::to_string(rank) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    meshAdapter.findPeriodicBoundaryNeighbors();

    //meshAdapter.writeMeshFaceVTK(path + "/Faces_rank_" + std::to_string(rank) + ".vtk");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.1, 0.0), false );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    //SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcMY = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    SPtr<BoundaryCondition> bcPY = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*H; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*H; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////
    if( mpiWorldSize == 1 )
    {
        dataBase->boundaryConditions.push_back( bcMX );
        dataBase->boundaryConditions.push_back( bcPX );
    }
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    //////////////////////////////////////////////////////////////////////////

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
    
    //*logging::out << logging::Logger::WARNING << int(dataBase->communicators[0].size()) << "\n";
    //*logging::out << logging::Logger::WARNING << int(dataBase->communicators[0][0].get()) << "\n";
    //*logging::out << logging::Logger::WARNING << int(dataBase->communicators[0][1].get()) << "\n";

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables
    {
        return toConservedVariables( PrimitiveVariables( 1.0, 1.0, 0.0, 0.0, parameters.lambdaRef ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, level );

    Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataDeviceToHost();

    if( rank == 0 ) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_0", mpiWorldSize );

    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" + "_rank_" + std::to_string(rank) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const uint numberOfIterations = 1000;

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0, true, numberOfIterations );

    MPI_Barrier(MPI_COMM_WORLD);

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= numberOfIterations; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);
    }

    cupsAnalyzer.run( numberOfIterations, parameters.dt );

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_final_rank_" + std::to_string(rank) );
    
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
}

int main( int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////

    int rank = 0;
    int mpiWorldSize = 1;
#ifdef USE_CUDA_AWARE_MPI
    int rank         = MpiUtility::getMpiRankBeforeInit();
    int mpiWorldSize = MpiUtility::getMpiWorldSizeBeforeInit();
#else
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
#endif

    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/MultiGPU/" );
#else
    //std::string path( "/home/stephan/Computations/out/" );
    std::string path( "out/" );
#endif

    std::string simulationName ( "MultiGPU_np_" + std::to_string(mpiWorldSize) );

    //////////////////////////////////////////////////////////////////////////

    bool strongScaling = false;
    uint nx = 128;

    if( argc > 1 ) path += argv[1]; path += "/";
    if( argc > 2 ) nx = atoi( argv[2] );
    if( argc > 3 ) strongScaling = true;

    //////////////////////////////////////////////////////////////////////////

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + "_rank_" + std::to_string(rank) + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    //////////////////////////////////////////////////////////////////////////

    // Important: for Cuda-Aware MPI the device must be set before MPI_Init()
    int deviceCount = CudaUtility::getCudaDeviceCount();

    if(deviceCount == 0)
    {
        std::stringstream msg;
        msg << "No devices devices found!" << std::endl;
        *logging::out << logging::Logger::WARNING << msg.str(); msg.str("");
    }

    CudaUtility::setCudaDevice( rank % deviceCount );

    //////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA_AWARE_MPI
    MPI_Init(&argc, &argv);
#endif
    
    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    try
    {
        performanceTest( path, simulationName, 1, nx, strongScaling );
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

    logFile.close();

    MPI_Finalize();

   return 0;
}
