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
#include "PointerDefinitions.h"
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

    int sideLengthX, sideLengthY, sideLengthZ, rankX, rankY, rankZ;

    if( decompositionDimension == 1 )
    {
        if      (mpiWorldSize == 1 ) { sideLengthX = 1 ; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 2 ) { sideLengthX = 2 ; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 4 ) { sideLengthX = 4 ; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 8 ) { sideLengthX = 8 ; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 16) { sideLengthX = 16; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 32) { sideLengthX = 32; sideLengthY = 1; sideLengthZ = 1; }

        rankX = rank;
        rankY = 0;
        rankZ = 0;
    }
    else if( decompositionDimension == 2 )
    {
        if      (mpiWorldSize == 1 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 2 ) { sideLengthX = 2; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 4 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 1; }
        else if (mpiWorldSize == 8 ) { sideLengthX = 4; sideLengthY = 2; sideLengthZ = 1; }
        else if (mpiWorldSize == 16) { sideLengthX = 4; sideLengthY = 4; sideLengthZ = 1; }
        else if (mpiWorldSize == 32) { sideLengthX = 8; sideLengthY = 4; sideLengthZ = 1; }

        rankX = rank % sideLengthX;
        rankY = rank / sideLengthX;
        rankZ = 0;
    }
    else if( decompositionDimension == 3 )
    {
        if      (mpiWorldSize == 1 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 2 ) { sideLengthX = 2; sideLengthY = 1; sideLengthZ = 1; }
        else if (mpiWorldSize == 4 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 1; }
        else if (mpiWorldSize == 8 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 2; }
        else if (mpiWorldSize == 16) { sideLengthX = 4; sideLengthY = 2; sideLengthZ = 2; }
        else if (mpiWorldSize == 32) { sideLengthX = 4; sideLengthY = 4; sideLengthZ = 2; }

        rankX =   rank %   sideLengthX;
        rankY = ( rank % ( sideLengthX * sideLengthY ) ) /   sideLengthX;
        rankZ =   rank                                   / ( sideLengthY * sideLengthX );
    }

    *logging::out << logging::Logger::INFO_HIGH << "SideLength = " << sideLengthX << " " << sideLengthY << " " << sideLengthZ << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "rank       = " << rankX << " " << rankY << " " << rankZ << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L  = 1.0;

    real LX = L;
    real LY = L;
    real LZ = L;

    real dx = L / real(nx);

    if( strongScaling )
    {
        if( decompositionDimension == 1 )
        {
            LX /= double(sideLengthX);
        }
        else if( decompositionDimension == 2 )
        {
            LX /= double(sideLengthX);
            LY /= double(sideLengthY);
        }
        else if( decompositionDimension == 3 )
        {
            LX /= double(sideLengthX);
            LY /= double(sideLengthY);
            LZ /= double(sideLengthZ);
        }
    }

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

    real xOverlap = ( sideLengthX == 1 ) ? 0.0 : 5.0*dx;
    real yOverlap = ( sideLengthY == 1 ) ? 0.0 : 5.0*dx;
    real zOverlap = ( sideLengthZ == 1 ) ? 0.0 : 5.0*dx;

    gridBuilder->addCoarseGrid(  rankX*LX    - 0.5*L - xOverlap,      rankY*LY    - 0.5*L - yOverlap,      rankZ*LZ    - 0.5*L - zOverlap,
                                (rankX*LX+1) - 0.5*L + xOverlap,     (rankY*LY+1) - 0.5*L + yOverlap,     (rankZ*LZ+1) - 0.5*L + zOverlap, dx);

    gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( rankX*LX - 0.5*L, (rankX+1)*LX - 0.5*L, 
                                                                 rankY*LY - 0.5*L, (rankY+1)*LY - 0.5*L,
                                                                 rankZ*LZ - 0.5*L, (rankZ+1)*LZ - 0.5*L  ) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(sideLengthX == 1, sideLengthY == 1, sideLengthZ == 1);

    *logging::out << logging::Logger::INFO_HIGH << "periodicity = " << (sideLengthX == 1) << " " << (sideLengthY == 1) << " " << (sideLengthZ == 1) << "\n";

    gridBuilder->buildGrids(GKS, false);

    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( mpiWorldSize > 1 )
    {
        int rankPX = ( (rankX + 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMX = ( (rankX - 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPY =    rankX                                    + ( (rankY + 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMY =    rankX                                    + ( (rankY - 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ + 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;
        int rankMZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ - 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;

        if( sideLengthX > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
        if( sideLengthX > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, rankPX);

        if( sideLengthX > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
        if( sideLengthX > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, rankMX);

        if( sideLengthY > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PY, GKS );
        if( sideLengthY > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PY, rankPY);

        if( sideLengthY > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MY, GKS );
        if( sideLengthY > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MY, rankMY);

        if( sideLengthZ > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PZ, GKS );
        if( sideLengthZ > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PZ, rankPZ);

        if( sideLengthZ > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MZ, GKS );
        if( sideLengthZ > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MZ, rankMZ);

        *logging::out << logging::Logger::INFO_HIGH << "neighborRanks = " << rankPX << " " << rankMX << " " << rankPY << " " << rankMY << " " << rankPZ << " " << rankMZ << "\n";
    }

    //gridBuilder->writeGridsToVtk(path + "/Grid_rank_" + std::to_string(rank) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto dataBase = std::make_shared<GksGpu::DataBase>("GPU");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for ( int i = 0; i < rank % GksGpu::CudaUtility::getCudaDeviceCount(); i++ ) MPI_Barrier(MPI_COMM_WORLD);

    {
        GksMeshAdapter meshAdapter(gridBuilder);

        meshAdapter.inputGrid();

        if (sideLengthX == 1 || sideLengthY == 1 || sideLengthZ == 1) meshAdapter.findPeriodicBoundaryNeighbors();

        gridBuilder->getGrid(0)->freeMemory();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SPtr<GksGpu::BoundaryCondition> bcMX = std::make_shared<GksGpu::Periodic>(dataBase);
        SPtr<GksGpu::BoundaryCondition> bcPX = std::make_shared<GksGpu::Periodic>(dataBase);

        if (sideLengthX == 1) bcMX->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.x < -0.5*L; });
        if (sideLengthX == 1) bcPX->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.x > 0.5*L; });

        //////////////////////////////////////////////////////////////////////////

        SPtr<GksGpu::BoundaryCondition> bcMY = std::make_shared<GksGpu::Periodic>(dataBase);
        SPtr<GksGpu::BoundaryCondition> bcPY = std::make_shared<GksGpu::Periodic>(dataBase);

        if (sideLengthY == 1) bcMY->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.y < -0.5*L; });
        if (sideLengthY == 1) bcPY->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.y > 0.5*L; });

        //////////////////////////////////////////////////////////////////////////

        SPtr<GksGpu::BoundaryCondition> bcMZ = std::make_shared<GksGpu::Periodic>(dataBase);
        SPtr<GksGpu::BoundaryCondition> bcPZ = std::make_shared<GksGpu::Periodic>(dataBase);

        if (sideLengthZ == 1) bcMZ->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.z < -0.5*L; });
        if (sideLengthZ == 1) bcPZ->findBoundaryCells(meshAdapter, true, [&](Vec3 center) { return center.z > 0.5*L; });

        //////////////////////////////////////////////////////////////////////////

        if (sideLengthX == 1) dataBase->boundaryConditions.push_back(bcMX);
        if (sideLengthX == 1) dataBase->boundaryConditions.push_back(bcPX);

        if (sideLengthY == 1) dataBase->boundaryConditions.push_back(bcMY);
        if (sideLengthY == 1) dataBase->boundaryConditions.push_back(bcPY);

        if (sideLengthZ == 1) dataBase->boundaryConditions.push_back(bcMZ);
        if (sideLengthZ == 1) dataBase->boundaryConditions.push_back(bcPZ);

        //////////////////////////////////////////////////////////////////////////

        *logging::out << logging::Logger::INFO_HIGH << "NumberOfBoundaryConditions = " << (int)dataBase->boundaryConditions.size() << "\n";

        if (sideLengthX == 1) *logging::out << logging::Logger::INFO_HIGH << "bcMX ==> " << bcMX->numberOfCellsPerLevel[0] << "\n";
        if (sideLengthX == 1) *logging::out << logging::Logger::INFO_HIGH << "bcPX ==> " << bcPX->numberOfCellsPerLevel[0] << "\n";

        if (sideLengthY == 1) *logging::out << logging::Logger::INFO_HIGH << "bcMY ==> " << bcMY->numberOfCellsPerLevel[0] << "\n";
        if (sideLengthY == 1) *logging::out << logging::Logger::INFO_HIGH << "bcPY ==> " << bcPY->numberOfCellsPerLevel[0] << "\n";

        if (sideLengthZ == 1) *logging::out << logging::Logger::INFO_HIGH << "bcMZ ==> " << bcMZ->numberOfCellsPerLevel[0] << "\n";
        if (sideLengthZ == 1) *logging::out << logging::Logger::INFO_HIGH << "bcPZ ==> " << bcPZ->numberOfCellsPerLevel[0] << "\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        dataBase->setMesh(meshAdapter);

        dataBase->setCommunicators(meshAdapter);

        GksGpu::CudaUtility::printCudaMemoryUsage();
    }

    for ( int i = 0; i < GksGpu::CudaUtility::getCudaDeviceCount() - rank % GksGpu::CudaUtility::getCudaDeviceCount(); i++ ) MPI_Barrier(MPI_COMM_WORLD);

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

        //ULocal = cellCenter.x;
        //VLocal = cellCenter.y;
        //WLocal = cellCenter.z;

        //rhoLocal = rank + 1;

        return GksGpu::toConservedVariables( GksGpu::PrimitiveVariables( rhoLocal, ULocal, VLocal, WLocal, parameters.lambdaRef ), parameters.K );
    });

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, level );

    GksGpu::Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataDeviceToHost();

    //if( rank == 0 ) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_0", mpiWorldSize );

    //writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" + "_rank_" + std::to_string(rank) );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const uint numberOfIterations = 1000;

    GksGpu::CupsAnalyzer cupsAnalyzer( dataBase, false, 30.0, true, numberOfIterations );

    MPI_Barrier(MPI_COMM_WORLD);

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= numberOfIterations; iter++ )
    {
        GksGpu::TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        cupsAnalyzer.run( iter, parameters.dt );
    }

    //////////////////////////////////////////////////////////////////////////

    //dataBase->copyDataDeviceToHost();

    //if( rank == 0 ) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_final", mpiWorldSize );

    //writeVtkXML( dataBase, parameters, 0, path + simulationName + "_final_rank_" + std::to_string(rank) );
    
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

    //////////////////////////////////////////////////////////////////////////

    bool strongScaling = false;
    uint nx = 128;
    uint decompositionDimension = 3;

    if( argc > 1 ) nx = atoi( argv[1] );
    if( argc > 2 ) decompositionDimension = atoi( argv[2] );
    if( argc > 3 ) strongScaling = true;

    //////////////////////////////////////////////////////////////////////////

    std::string simulationName ( "MultiGPU" );

    if( strongScaling ) simulationName += "_strongScaling";
    else                simulationName += "_weakScaling";

    simulationName += "_D_" + std::to_string(decompositionDimension);

    simulationName += "_nx_" + std::to_string(nx);

    simulationName += "_np_" + std::to_string(mpiWorldSize);

    //////////////////////////////////////////////////////////////////////////

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + "_rank_" + std::to_string(rank) + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    //////////////////////////////////////////////////////////////////////////

    // Important: for Cuda-Aware MPI the device must be set before MPI_Init()
    int deviceCount = GksGpu::CudaUtility::getCudaDeviceCount();

    if(deviceCount == 0)
    {
        std::stringstream msg;
        msg << "No devices devices found!" << std::endl;
        *logging::out << logging::Logger::WARNING << msg.str(); msg.str("");
    }

    GksGpu::CudaUtility::setCudaDevice( rank % deviceCount );

    //////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA_AWARE_MPI
    MPI_Init(&argc, &argv);
#endif
    
    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    try
    {
        performanceTest( path, simulationName, decompositionDimension, nx, strongScaling );
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
