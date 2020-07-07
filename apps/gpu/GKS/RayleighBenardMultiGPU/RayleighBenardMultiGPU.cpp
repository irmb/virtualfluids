//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>
#include <thread>

#include <mpi.h>

#include "Core/Timer/Timer.h"
#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
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

#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"
#include "GksGpu/BoundaryConditions/Pressure.h"
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"

#include "GksGpu/Communication/Communicator.h"
#include "GksGpu/Communication/MpiUtility.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"
#include "GksGpu/Analyzer/PointTimeSeriesCollector.h"
#include "GksGpu/Analyzer/HeatFluxAnalyzer.h"

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

//uint deviceMap [2] = {2,3};
uint deviceMap [2] = {0,1};

void simulation( std::string path, std::string simulationName, bool fine, uint restartIter )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mpiWorldSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

    int sideLengthX, sideLengthY, sideLengthZ, rankX, rankY, rankZ;

    if      (mpiWorldSize == 1 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 1; }
    else if (mpiWorldSize == 2 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 2; }
    else if (mpiWorldSize == 4 ) { sideLengthX = 1; sideLengthY = 2; sideLengthZ = 2; }
    else if (mpiWorldSize == 8 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 2; }
    else
    {
        throw std::runtime_error( "This number of processes is not supported for this target!" );
    }

    rankZ =   rank %   sideLengthZ;
    rankY = ( rank % ( sideLengthZ * sideLengthY ) ) /   sideLengthZ;
    rankX =   rank                                   / ( sideLengthY * sideLengthZ );

    *logging::out << logging::Logger::INFO_HIGH << "SideLength = " << sideLengthX << " " << sideLengthY << " " << sideLengthZ << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "rank       = " << rankX << " " << rankY << " " << rankZ << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 64;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L / real(nx);

    real Ra = 3.0e6;
    //real Ra = 1.0e2;

    real Ba  = 0.1;
    real eps = 0.8;
    real Pr  = 0.71;
    real K   = 2.0;
    
    real g   = 1.0;
    real rho = 1.0;

    real lambda     = Ba / ( 2.0 * g * L );
    real lambdaHot  = lambda / ( 1.0 + eps * 0.5 );
    real lambdaCold = lambda / ( 1.0 - eps * 0.5 );
    
    real mu = sqrt( Pr * eps * g * L * L * L / Ra ) * rho ;

    real cs  = sqrt( ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( 2.0 * lambda ) );
    real U   = sqrt( Ra ) * mu / ( rho * L );

    real CFL = 0.5;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( c1o1 + ( c2o1 * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " s\n";

    //////////////////////////////////////////////////////////////////////////

    GksGpu::Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = -g;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = lambda;

    parameters.viscosityModel = GksGpu::ViscosityModel::sutherlandsLaw2;
    //parameters.viscosityModel = ViscosityModel::constant;

    parameters.forcingSchemeIdx = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                M e s h    G e n e r a t i o n
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real LX = L / double(sideLengthX);
    real LY = L / double(sideLengthY);
    real LZ = L / double(sideLengthZ);

    real xOverlap = ( sideLengthX == 1 ) ? 0.0 : 5.0*dx;
    real yOverlap = ( sideLengthY == 1 ) ? 0.0 : 5.0*dx;
    real zOverlap = ( sideLengthZ == 1 ) ? 0.0 : 5.0*dx;

    real startX, endX;
    real startY, endY;
    real startZ, endZ;

    if( sideLengthX > 1 && rankX == 1 ) startX = -3.0 * dx;
    else                                startX = -0.5 * L;
    if( sideLengthX > 1 && rankX == 0 ) endX   =  3.0 * dx;
    else                                endX   =  0.5 * L;

    if( sideLengthY > 1 && rankY == 1 ) startY = -3.0 * dx;
    else                                startY = -0.5 * L;
    if( sideLengthY > 1 && rankY == 0 ) endY   =  3.0 * dx;
    else                                endY   =  0.5 * L;

    if( sideLengthZ > 1 && rankZ == 1 ) startZ = -3.0 * dx;
    else                                startZ = -0.5 * L;
    if( sideLengthZ > 1 && rankZ == 0 ) endZ   =  3.0 * dx;
    else                                endZ   =  0.5 * L;

    gridBuilder->addCoarseGrid(startX, startY, startZ,  
                               endX  , endY  , endZ  , dx);

    std::cout << __LINE__ << std::endl;

    //////////////////////////////////////////////////////////////////////////

    real refL[4] = { 0.05, 0.02, 0.025, 0.005 };

    if( fine )
    {
        refL[1] = 0.1;
        refL[2] = 0.05;
    }

    gridBuilder->setNumberOfLayers(6,6);

    //////////////////////////////////////////////////////////////////////////

    Conglomerate coarseRefLevel;

    if( sideLengthX == 1 || rankX == 0 ) coarseRefLevel.add( new Cuboid (-100.0, -100.0,           -100.0, 
                                                                          100.0, -0.5*L + refL[0],  100.0 ) );
    if( sideLengthX == 1 || rankX == 1 ) coarseRefLevel.add( new Cuboid (-100.0,  0.5*L - refL[0], -100.0, 
                                                                          100.0,  100.0,            100.0 ) );

    if( sideLengthY == 1 || rankY == 0 ) coarseRefLevel.add( new Cuboid (-100.0,           -100.0, -100.0, 
                                                                         -0.5*L + refL[0],  100.0,  100.0 ) );
    if( sideLengthY == 1 || rankY == 1 ) coarseRefLevel.add( new Cuboid ( 0.5*L - refL[0], -100.0, -100.0, 
                                                                          100.0,            100.0,  100.0  ) );

    if( sideLengthZ == 1 || rankZ == 0 ) coarseRefLevel.add( new Cuboid (-100.0, -100.0, -100.0, 
                                                                          100.0,  100.0, -0.5*L + refL[0] ) );
    if( sideLengthZ == 1 || rankZ == 1 ) coarseRefLevel.add( new Cuboid (-100.0, -100.0,  0.5*L - refL[0], 
                                                                          100.0,  100.0,  100.0           ) );

    gridBuilder->addGrid( &coarseRefLevel, 1);

    //////////////////////////////////////////////////////////////////////////

    Conglomerate firstRefLevel;

    if( sideLengthZ == 1 || rankZ == 0 ) firstRefLevel.add( new Cuboid (-100.0, -100.0, -100.0, 
                                                                         100.0,  100.0, -0.5*L + refL[1] ) );
    if( sideLengthZ == 1 || rankZ == 1 ) firstRefLevel.add( new Cuboid (-100.0, -100.0,  0.5*L - refL[1], 
                                                                         100.0,  100.0,  100.0           ) );

    gridBuilder->addGrid( &firstRefLevel, 2);

    //////////////////////////////////////////////////////////////////////////

    //Conglomerate secondRefLevel;

    //if( rank % 2 == 0 ) secondRefLevel.add( new Cuboid (-100.0,           -100.0, -100.0, 
    //                                                    -0.5*L + refL[2],  100.0,  100.0 ) );
    //else                secondRefLevel.add( new Cuboid ( 0.5*L - refL[2], -100.0, -100.0, 
    //                                                     100.0,            100.0,  100.0 ) );

    //if( rank % 2 == 0 ) secondRefLevel.add( new Cuboid (-100.0,           -100.0, -100.0,   
    //                                                    -0.5*L + refL[0],  100.0, -0.5*H + refL[2] ) );
    //else                secondRefLevel.add( new Cuboid ( 0.5*L - refL[0], -100.0, -100.0,   
    //                                                     100.0,            100.0, -0.5*H + refL[2] ) );

    //if( rank % 2 == 0 ) secondRefLevel.add( new Cuboid (-100.0,           -100.0,  0.5*H - refL[2], 
    //                                                    -0.5*L + refL[0],  100.0,  100.0   ) );
    //else                secondRefLevel.add( new Cuboid ( 0.5*L - refL[0], -100.0,  0.5*H - refL[2], 
    //                                                     100.0,            100.0,  100.0   ) );

    //gridBuilder->addGrid( &secondRefLevel, 3);

    //////////////////////////////////////////////////////////////////////////

    //Conglomerate thirdRefLevel;

    //if( rank % 2 == 0 ) thirdRefLevel.add( new Cuboid (-100.0,           -100.0, -100.0, 
    //                                                   -0.5*L + refL[3],  100.0,  100.0 ) );
    //else                thirdRefLevel.add( new Cuboid ( 0.5*L - refL[3], -100.0, -100.0, 
    //                                                    100.0,            100.0,  100.0 ) );

    //if( fine ) gridBuilder->addGrid( &thirdRefLevel, 4);

    //////////////////////////////////////////////////////////////////////////

    if( sideLengthX > 1 && rankX == 1 ) startX =    0.0;
    else                                startX = -100.0;
    if( sideLengthX > 1 && rankX == 0 ) endX   =    0.0;
    else                                endX   =  100.0;

    if( sideLengthY > 1 && rankY == 1 ) startY =    0.0;
    else                                startY = -100.0;
    if( sideLengthY > 1 && rankY == 0 ) endY   =    0.0;
    else                                endY   =  100.0;

    if( sideLengthZ > 1 && rankZ == 1 ) startZ =    0.0;
    else                                startZ = -100.0;
    if( sideLengthZ > 1 && rankZ == 0 ) endZ   =    0.0;
    else                                endZ   =  100.0;

    auto subDomainBox = std::make_shared<BoundingBox>( startX, endX, 
                                                       startY, endY, 
                                                       startZ, endZ );

    if( mpiWorldSize > 1 ) gridBuilder->setSubDomainBox( subDomainBox );

    //////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);
            
    //gridBuilder->writeGridsToVtk( path + simulationName + "_0" + "_rank_" + std::to_string(rank) + "_lev_" );

    //////////////////////////////////////////////////////////////////////////

    if( mpiWorldSize > 1 )
    {
        int rankPX = ( (rankX + 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMX = ( (rankX - 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPY =    rankX                                    + ( (rankY + 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMY =    rankX                                    + ( (rankY - 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ + 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;
        int rankMZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ - 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;

        if( sideLengthX > 1 && rankX == 0 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
        if( sideLengthX > 1 && rankX == 0 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, rankPX);

        if( sideLengthX > 1 && rankX == 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
        if( sideLengthX > 1 && rankX == 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, rankMX);

        if( sideLengthY > 1 && rankY == 0 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PY, GKS );
        if( sideLengthY > 1 && rankY == 0 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PY, rankPY);

        if( sideLengthY > 1 && rankY == 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MY, GKS );
        if( sideLengthY > 1 && rankY == 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MY, rankMY);

        if( sideLengthZ > 1 && rankZ == 0 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PZ, GKS );
        if( sideLengthZ > 1 && rankZ == 0 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PZ, rankPZ);

        if( sideLengthZ > 1 && rankZ == 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MZ, GKS );
        if( sideLengthZ > 1 && rankZ == 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MZ, rankMZ);

        *logging::out << logging::Logger::INFO_HIGH << "neighborRanks = " << rankPX << " " << rankMX << " " << rankPY << " " << rankMY << " " << rankPZ << " " << rankMZ << "\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //if( mpiWorldSize == 2 ) meshAdapter.findPeriodicBoundaryNeighbors();    

    //meshAdapter.writeMeshFaceVTK( path + simulationName + "_0" + "_rank_" + std::to_string(rank) + ".vtk" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto dataBase = std::make_shared<GksGpu::DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                 B o u n d a r y    C o n d i t i o n s
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<GksGpu::BoundaryCondition> bcMX = std::make_shared<GksGpu::AdiabaticWall>( dataBase, Vec3(0,0,0), false );
    SPtr<GksGpu::BoundaryCondition> bcPX = std::make_shared<GksGpu::AdiabaticWall>( dataBase, Vec3(0,0,0), false );

    SPtr<GksGpu::BoundaryCondition> bcMY = std::make_shared<GksGpu::AdiabaticWall>( dataBase, Vec3(0,0,0), false );
    SPtr<GksGpu::BoundaryCondition> bcPY = std::make_shared<GksGpu::AdiabaticWall>( dataBase, Vec3(0,0,0), false );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<GksGpu::BoundaryCondition> bcMZ = std::make_shared<GksGpu::IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaHot , false );
    SPtr<GksGpu::BoundaryCondition> bcPZ = std::make_shared<GksGpu::IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold, false );

    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*L; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );

    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                 I n i t i a l    C o n d i t i o n s
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    dataBase->setCommunicators( meshAdapter );

    GksGpu::CudaUtility::printCudaMemoryUsage();

    if( restartIter == INVALID_INDEX )
    {
        GksGpu::Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> GksGpu::ConservedVariables {

            //real Th = 1.0 / lambdaHot;
            //real Tc = 1.0 / lambdaCold;
            //real T = Th - (Th - Tc)*((cellCenter.x + 0.5 * L) / L);
            //real lambdaLocal = 1.0 / T;

            return GksGpu::toConservedVariables(GksGpu::PrimitiveVariables(rho, 0.0, 0.0, 0.0, lambda), parameters.K);
        });

        if (rank == 0) writeVtkXMLParallelSummaryFile(dataBase, parameters, path + simulationName + "_0", mpiWorldSize);

        writeVtkXML(dataBase, parameters, 0, path + simulationName + "_0" + "_rank_" + std::to_string(rank));
    }
    else
    {
        GksGpu::Restart::readRestart( dataBase, path + simulationName + "_" + std::to_string( restartIter ) + "_rank_" + std::to_string(rank), startIter );

        if (rank == 0) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_" + std::to_string( restartIter ) + "_restart", mpiWorldSize );

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( restartIter ) + "_restart" + "_rank_" + std::to_string(rank) );


    }

    dataBase->copyDataHostToDevice();

    GksGpu::Initializer::initializeDataUpdate(dataBase);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                  R u n
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksGpu::CupsAnalyzer cupsAnalyzer( dataBase, true, 300.0 );

    GksGpu::ConvergenceAnalyzer convergenceAnalyzer( dataBase );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 0 );
    auto turbulenceAnalyzer = std::make_shared<GksGpu::TurbulenceAnalyzer>( dataBase, 50000 );

    turbulenceAnalyzer->collect_UU = true;
    turbulenceAnalyzer->collect_VV = true;
    turbulenceAnalyzer->collect_WW = true;
    turbulenceAnalyzer->collect_UV = true;
    turbulenceAnalyzer->collect_UW = true;
    turbulenceAnalyzer->collect_VW = true;

    turbulenceAnalyzer->allocate();

    if( restartIter != INVALID_INDEX )
        turbulenceAnalyzer->readRestartFile( path + simulationName + "_Turbulence_" + std::to_string( restartIter ) + "_rank_" + std::to_string(rank) );

    //auto pointTimeSeriesCollector = std::make_shared<PointTimeSeriesCollector>();

    //for( real y = 0.5 * W; y < real( mpiWorldSize / 2 ) * W; y += W )
    //{
    //    if( subDomainBox->isInside( -0.485, y, -0.3*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( -0.485, y, -0.3*H ), 'W', 10000 );
    //    if( subDomainBox->isInside( -0.485, y, -0.1*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( -0.485, y, -0.1*H ), 'W', 10000 );
    //    if( subDomainBox->isInside( -0.485, y,  0.1*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( -0.485, y,  0.1*H ), 'W', 10000 );
    //    if( subDomainBox->isInside( -0.485, y,  0.3*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( -0.485, y,  0.3*H ), 'W', 10000 );
    //    
    //    if( subDomainBox->isInside(  0.485, y, -0.3*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3(  0.485, y, -0.3*H ), 'W', 10000 );
    //    if( subDomainBox->isInside(  0.485, y, -0.1*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3(  0.485, y, -0.1*H ), 'W', 10000 );
    //    if( subDomainBox->isInside(  0.485, y,  0.1*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3(  0.485, y,  0.1*H ), 'W', 10000 );
    //    if( subDomainBox->isInside(  0.485, y,  0.3*H ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3(  0.485, y,  0.3*H ), 'W', 10000 );
    //}

    GksGpu::HeatFluxAnalyzer heatFluxAnalyzerPZ(dataBase, bcPZ, 100, 10000, lambdaHot, lambdaCold, L);
    GksGpu::HeatFluxAnalyzer heatFluxAnalyzerMZ(dataBase, bcMZ, 100, 10000, lambdaHot, lambdaCold, L);
    //HeatFluxAnalyzer heatFluxAnalyzer(dataBase, bcPZ);

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 100000000; iter++ )
    {
        GksGpu::TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        turbulenceAnalyzer->run( iter, parameters );

        if(rankZ == 1) heatFluxAnalyzerPZ.run( iter, parameters );
        if(rankZ == 0) heatFluxAnalyzerMZ.run( iter, parameters );

        if( iter % 10000 == 0 )
        //if( iter % 25 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            if( rank == 0 ) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_" + std::to_string( iter ), mpiWorldSize );

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );

            if(rankZ == 1) heatFluxAnalyzerPZ.writeToFile( path + simulationName + "_Nu_top_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );
            if(rankZ == 0) heatFluxAnalyzerMZ.writeToFile( path + simulationName + "_Nu_bot_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );
        }

        //pointTimeSeriesCollector->run(iter, parameters);

        if( iter > 50000 && iter % 10000 == 0 )
        {
            turbulenceAnalyzer->download();
        
            if( rank == 0 ) writeTurbulenceVtkXMLParallelSummaryFile( dataBase, turbulenceAnalyzer, parameters, path + simulationName + "_Turbulence_" + std::to_string( iter ), mpiWorldSize );
        
            writeTurbulenceVtkXML( dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );
        }

        if( iter % 10000 == 0 )
        {
            GksGpu::Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ) + "_rank_" + std::to_string(rank), iter );

            turbulenceAnalyzer->writeRestartFile( path + simulationName + "_Turbulence_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );
        }

        //if( iter % 1000000 == 0 )
        //{
        //    pointTimeSeriesCollector->writeToFile(path + simulationName + "_TimeSeries_" + std::to_string( iter ) + "_rank_" + std::to_string(rank));
        //}
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();
}



int main( int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////

    bool fine = false;

    bool highAspect = true;

    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int mpiWorldSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
#else
    int rank         = GksGpu::MpiUtility::getMpiRankBeforeInit();
    int mpiWorldSize = GksGpu::MpiUtility::getMpiWorldSizeBeforeInit();
#endif

    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/RayleighBenardMultiGPU/test/" );
    //std::string path( "F:/Work/Computations/out/RayleighBenardMultiGPU/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "ThermalCavity3D" );

    if(fine) simulationName += "_fine";
    else     simulationName += "_coarse";

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

#ifndef _WIN32
    MPI_Init(&argc, &argv);
#endif

    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        uint restartIter = INVALID_INDEX;

        if( argc > 1 ) restartIter = atoi( argv[1] );

        simulation(path, simulationName, fine, restartIter);
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

    MPI_Finalize();

    return 0;
}
