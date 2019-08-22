
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

real getHRR( real t );

void thermalCavity( std::string path, std::string simulationName, uint windowIndex, uint restartIter, bool useConreteHeatFluxBC )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mpiWorldSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int sideLengthX, sideLengthY, sideLengthZ, rankX, rankY, rankZ;

    if      (mpiWorldSize == 1 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 1; }
    else if (mpiWorldSize == 2 ) { sideLengthX = 2; sideLengthY = 1; sideLengthZ = 1; }
    else if (mpiWorldSize == 4 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 1; }
    else if (mpiWorldSize == 8 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 2; }

    rankX =   rank %   sideLengthX;
    rankY = ( rank % ( sideLengthX * sideLengthY ) ) /   sideLengthX;
    rankZ =   rank                                   / ( sideLengthY * sideLengthX );

    *logging::out << logging::Logger::INFO_HIGH << "SideLength = " << sideLengthX << " " << sideLengthY << " " << sideLengthZ << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "rank       = " << rankX << " " << rankY << " " << rankZ << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dx = 0.1;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 4.0;
    real B = 3.0;

    real LBurner = 1.0;
    real HBurner = 0.5;

    real Pr  = 0.71;
    real K   = 2.0;
    
    real g   = 9.81;
    real rho = 1.2;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );
    setLambdaFromT( prim, 3.0 );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real mu      = 1.8e-5;
    //real U       = 0.025;       // 750 kW on top
    //real U       = 0.015;       // 900 kW on top
    //real U       = 0.005;       // 900 kW all around
    real rhoFuel = 0.5405;

    real heatOfReaction = real(8000.0); // J / mol 

    real specificHeatOfReaction = heatOfReaction / 0.016;

    real HRR = 750.0; // kW

    real U = HRR * 1000.0 / ( rhoFuel * LBurner * LBurner * (specificHeatOfReaction * 100.0) );

    real CFL = 0.125;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( c1o1 + ( c2o1 * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
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
    parameters.force.z = -g;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = prim.lambda;

    parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    //parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    parameters.heatOfReaction = heatOfReaction;

    parameters.useHeatReleaseRateLimiter = true;
    parameters.useTemperatureLimiter     = true;
    parameters.usePassiveScalarLimiter   = true;
    parameters.useSmagorinsky            = true;

    parameters.reactionLimiter    = 1.0005;
    parameters.temperatureLimiter = 1.0e-3;

    parameters.useSpongeLayer = true;
    parameters.spongeLayerIdx = 2;

    parameters.forcingSchemeIdx = 2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //gridBuilder->addCoarseGrid(-2.1, -1.6, -0.1,  
                                //2.1,  6.0,  5.0, dx);
    //gridBuilder->addCoarseGrid(-1.1, -1.2, -0.1,  
                                //1.1,  1.2,  2.2, dx);
    gridBuilder->addCoarseGrid(-2.1, -1.6, -0.1,  
                                2.1,  1.6,  3.1, dx);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#ifdef _WIN32
//    TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/RoomExtended7.stl");
//#else
//    //TriangularMesh* stl = TriangularMesh::make("inp/Unterzug.stl");
//    TriangularMesh* stl = TriangularMesh::make("inp/RoomExtended4.stl");
//#endif
//
//    gridBuilder->addGeometry(stl);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Conglomerate flowDomain;

    flowDomain.add( new Cuboid( -2.0, -1.5, 0.0, 2.0,  1.5, 3.0 ) );      // Room 
    flowDomain.add( new Cuboid( -2.0,  1.8, 0.0, 2.0,  5.0, 5.0 ) );      // Outside
    flowDomain.add( new Cuboid( -0.5, -1.8, 0.0, 0.5, -1.0, 2.0 ) );      // Door
    flowDomain.subtract( new Cuboid( -0.5, -0.5, -1.0, 0.5, 0.5, 0.5 ) ); // Fire
    flowDomain.subtract( new Cuboid( -3.0, -0.1,  2.6, 3.0, 0.1, 4.0 ) ); // Beam

    if( windowIndex == 0 ) flowDomain.add( new Cuboid( -1.0 ,  1.0,  1.0,    1.0 ,  3.0,  2.4 ) );      // Window large
    if( windowIndex == 1 ) flowDomain.add( new Cuboid( -0.5 ,  1.0,  1.0,    0.5 ,  3.0,  2.4 ) );      // Window medium
    if( windowIndex == 2 ) flowDomain.add( new Cuboid( -0.25,  1.0,  1.5,    0.25,  3.0,  2.0 ) );      // Window small
    if( windowIndex == 3 ) flowDomain.add( new Cuboid( -1.0 ,  1.0,  1.0,    1.0 ,  3.0,  2.0 ) );      // Window low

    Conglomerate solidDomain;

    solidDomain.add( new Cuboid(-2.2, -1.7, -0.2, 2.2,  6.1,  5.1) );
    solidDomain.subtract( &flowDomain );

    gridBuilder->addGeometry( &solidDomain );
    //gridBuilder->addGeometry( &flowDomain );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cuboid boxCoarse ( -2.0, -3.0, -0.5, 
                        3.0,  3.0,  3.5 );

    gridBuilder->addGrid( &boxCoarse, 1 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real startX = -1e99;
    real startY = -1e99;
    real startZ = -1e99;
    real endX   =  1e99;
    real endY   =  1e99;
    real endZ   =  1e99;

    if( mpiWorldSize == 2 )
    {
        if( rank == 0 ) { endX   = 0.0; }
        if( rank == 1 ) { startX = 0.0; }
    }
    if( mpiWorldSize == 4 )
    {
        if( rank == 0 ) { endX   = 0.0; endY   = 0.0; }
        if( rank == 1 ) { startX = 0.0; endY   = 0.0; }
        if( rank == 2 ) { endX   = 0.0; startY = 0.0; }
        if( rank == 3 ) { startX = 0.0; startY = 0.0; }
    }
    if( mpiWorldSize == 8 )
    {
        if( rank == 0 ) { endX   = 0.0; endY   = 0.0; endZ   = 1.9; }
        if( rank == 1 ) { startX = 0.0; endY   = 0.0; endZ   = 1.9; }
        if( rank == 2 ) { endX   = 0.0; startY = 0.0; endZ   = 1.9; }
        if( rank == 3 ) { startX = 0.0; startY = 0.0; endZ   = 1.9; }
        if( rank == 4 ) { endX   = 0.0; endY   = 0.0; startZ = 1.9; }
        if( rank == 5 ) { startX = 0.0; endY   = 0.0; startZ = 1.9; }
        if( rank == 6 ) { endX   = 0.0; startY = 0.0; startZ = 1.9; }
        if( rank == 7 ) { startX = 0.0; startY = 0.0; startZ = 1.9; }
    }

    auto subDomainBox = std::make_shared<BoundingBox>( startX, endX, 
                                                       startY, endY, 
                                                       startZ, endZ );

    gridBuilder->setSubDomainBox( subDomainBox );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cuboid roomRef( -2.1, -1.8, -1.0, 
                     2.1,  1.7, 10.0 );
    
    Cuboid windowRef( -1.1,  1.6,  0.9, 
                       1.1,  2.0,  3.0 );

    Conglomerate refRegion1;

    refRegion1.add( &roomRef );
    refRegion1.add( &windowRef );

    gridBuilder->setNumberOfLayers(0,22);

    //gridBuilder->addGrid( &refRegion1, 2 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cuboid boxRef ( -0.6 * LBurner, -0.6 * LBurner, -1.0, 
                     0.6 * LBurner,  0.6 * LBurner, 10.0 );
    Cuboid beamRef( -10.0, -0.25, 2.4, 10.0, 0.25, 10.0 );

    Conglomerate refRegion2;

    refRegion2.add( &boxRef );
    refRegion2.add( &beamRef );

    gridBuilder->setNumberOfLayers(0,22);

    //gridBuilder->addGrid( &refRegion2, 3 );

    uint maxLevel = gridBuilder->getNumberOfGridLevels() - 1;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);

    MPI_Barrier(MPI_COMM_WORLD);

    //gridBuilder->writeGridsToVtk(path + "Grid_rank_" + std::to_string( rank ) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    if( mpiWorldSize > 1 )
    {
        int rankPX = ( (rankX + 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMX = ( (rankX - 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPY =    rankX                                    + ( (rankY + 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMY =    rankX                                    + ( (rankY - 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ + 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;
        int rankMZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ - 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;

        if( sideLengthX > 1 && rankX < sideLengthX-1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
        if( sideLengthX > 1 && rankX < sideLengthX-1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, rankPX);

        if( sideLengthX > 1 && rankX > 0             ) gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
        if( sideLengthX > 1 && rankX > 0             ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, rankMX);

        if( sideLengthY > 1 && rankY < sideLengthY-1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PY, GKS );
        if( sideLengthY > 1 && rankY < sideLengthY-1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PY, rankPY);

        if( sideLengthY > 1 && rankY > 0             ) gridBuilder->findCommunicationIndices( CommunicationDirections::MY, GKS );
        if( sideLengthY > 1 && rankY > 0             ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MY, rankMY);

        if( sideLengthZ > 1 && rankZ < sideLengthZ-1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PZ, GKS );
        if( sideLengthZ > 1 && rankZ < sideLengthZ-1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PZ, rankPZ);

        if( sideLengthZ > 1 && rankZ > 0             ) gridBuilder->findCommunicationIndices( CommunicationDirections::MZ, GKS );
        if( sideLengthZ > 1 && rankZ > 0             ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MZ, rankMZ);

        *logging::out << logging::Logger::INFO_HIGH << "neighborRanks = " << rankPX << " " << rankMX << " " << rankPY << " " << rankMY << " " << rankPZ << " " << rankMZ << "\n";
    }

    //gridBuilder->writeGridsToVtk(path + "Grid_rank_" + std::to_string( rank ) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "MeshFaces_rank_" + std::to_string( rank ) + ".vtk" );

    //meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //CudaUtility::setCudaDevice( rank % CudaUtility::getCudaDeviceCount() );

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcWall = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );

    bcWall->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return true; } );
    
    SPtr<BoundaryCondition> bcWallHeatFlux = std::make_shared<ConcreteHeatFlux>( dataBase, 10, 1.0e-6, 2400.0, 880, 0.3, 3.0 );

    bcWallHeatFlux->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > 3.0 && center.y < 1.6; } );

    std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->init();

    ////////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcOpen = std::make_shared<Open>( dataBase, prim, 1.0 );

    bcOpen->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -6.0 || center.y > 1.7; } );

    ////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcPressure = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );

    bcPressure->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y > 1.7 && center.z > 5.0; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcBurner = std::make_shared<CreepingMassFlux>( dataBase, rhoFuel, U, prim.lambda );

    bcBurner->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

        return center.z < HBurner && 
            std::sqrt(center.x*center.x) < 0.5 * LBurner - dx * std::pow(0.5, maxLevel) && 
            std::sqrt(center.y*center.y) < 0.5 * LBurner - dx * std::pow(0.5, maxLevel);
        //return center.z < HBurner && std::sqrt(center.x*center.x) < 0.5 * LBurner && std::sqrt(center.y*center.y) < 0.5 * LBurner;
    } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcBurner );

    if( useConreteHeatFluxBC )
        dataBase->boundaryConditions.push_back( bcWallHeatFlux );

    dataBase->boundaryConditions.push_back( bcWall );

    dataBase->boundaryConditions.push_back( bcOpen );

    dataBase->boundaryConditions.push_back( bcPressure );

    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcBurner = "   << bcBurner->numberOfCells   << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcWall = "     << bcWall->numberOfCells     << "\n";

    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcOpen = "     << bcOpen->numberOfCells     << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "Number of cells bcPressure = " << bcPressure->numberOfCells << "\n";

    //////////////////////////////////////////////////////////////////////////

    auto pointTimeSeriesCollector = std::make_shared<PointTimeSeriesCollector>();

    for( real x = 0.0002; x < 2; x += 0.4449 )
    {
        if( subDomainBox->isInside( x, -1.4999, 2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x, -1.4999, 2.9999 ), 'T' );
        if( subDomainBox->isInside( x, -1.0,    2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x, -1.0,    2.9999 ), 'T' );
        if( subDomainBox->isInside( x, -0.5,    2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x, -0.5,    2.9999 ), 'T' );
        if( subDomainBox->isInside( x, -0.2001, 2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x, -0.2001, 2.9999 ), 'T' );

        if( subDomainBox->isInside( x, -0.2001, 2.5999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x, -0.2001, 2.5999 ), 'T' );
        if( subDomainBox->isInside( x,  0.0001, 2.5999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  0.0001, 2.5999 ), 'T' );
        if( subDomainBox->isInside( x,  0.2001, 2.5999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  0.2001, 2.5999 ), 'T' );
        
        if( subDomainBox->isInside( x,  0.2001, 2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  0.2001, 2.9999 ), 'T' );
        if( subDomainBox->isInside( x,  0.5,    2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  0.5,    2.9999 ), 'T' );
        if( subDomainBox->isInside( x,  1.0,    2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  1.0,    2.9999 ), 'T' );
        if( subDomainBox->isInside( x,  1.4999, 2.9999 ) ) pointTimeSeriesCollector->addAnalyzer( dataBase, meshAdapter, Vec3( x,  1.4999, 2.9999 ), 'T' );
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    dataBase->setCommunicators( meshAdapter );

    CudaUtility::printCudaMemoryUsage();
    
    if( restartIter == INVALID_INDEX )
    {
        Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {

            PrimitiveVariables primLocal = prim;

            //if( cellCenter.x > 0 ) primLocal.rho = 1.21;

            return toConservedVariables(primLocal, parameters.K);
        });

        if (rank == 0) writeVtkXMLParallelSummaryFile(dataBase, parameters, path + simulationName + "_0", mpiWorldSize);

        writeVtkXML(dataBase, parameters, 0, path + simulationName + "_0" + "_rank_" + std::to_string(rank));

        if( useConreteHeatFluxBC )
            std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->writeVTKFile(dataBase, parameters, path + simulationName + "_Solid_0");
    }
    else
    {
        Restart::readRestart( dataBase, path + simulationName + "_" + std::to_string( restartIter ) + "_rank_" + std::to_string(rank), startIter );

        if (rank == 0) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_" + std::to_string( restartIter ) + "_restart", mpiWorldSize );

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( restartIter ) + "_restart" + "_rank_" + std::to_string(rank) );
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

    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "==================   S t a r t    T i m e    S t e p p i n g   =================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";
    *logging::out << logging::Logger::INFO_HIGH << "================================================================================\n";

    MPI_Barrier(MPI_COMM_WORLD);

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 100000000; iter++ )
    {
        real currentHRR = getHRR( iter * parameters.dt );

        //*logging::out << logging::Logger::LOGGER_ERROR << "HRR(t=" << iter * parameters.dt << ") = " << currentHRR << "\n";

        std::dynamic_pointer_cast<CreepingMassFlux>(bcBurner)->velocity = currentHRR * 1000.0 / ( rhoFuel * LBurner * LBurner * (specificHeatOfReaction * 100.0) );

        //////////////////////////////////////////////////////////////////////////

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        //////////////////////////////////////////////////////////////////////////

        cupsAnalyzer.run( iter, parameters.dt );

        convergenceAnalyzer.run( iter );

        //////////////////////////////////////////////////////////////////////////

        pointTimeSeriesCollector->run(iter, parameters);

        int crashCellIndex = dataBase->getCrashCellIndex();
        if( crashCellIndex >= 0 )
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Simulation Crashed at CellIndex = " << crashCellIndex << "\n";
            dataBase->copyDataDeviceToHost();
            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            break;
        }

        if( iter % 10000 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            if( rank == 0 ) writeVtkXMLParallelSummaryFile( dataBase, parameters, path + simulationName + "_" + std::to_string( iter ), mpiWorldSize );

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) + "_rank_" + std::to_string(rank) );

            if( useConreteHeatFluxBC )
                std::dynamic_pointer_cast<ConcreteHeatFlux>(bcWallHeatFlux)->writeVTKFile(dataBase, parameters, path + simulationName + "_Solid_" + std::to_string( iter ));
        }

        if( iter % 10000 == 0 )
        {
            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ) + "_rank_" + std::to_string(rank), iter );
        }

        if( iter % 100000 == 0 )
        {
            pointTimeSeriesCollector->writeToFile(path + simulationName + "_TimeSeries_" + std::to_string( iter ) + "_rank_" + std::to_string(rank));
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
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int mpiWorldSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

    //////////////////////////////////////////////////////////////////////////

    uint restartIter = INVALID_INDEX;
    //uint restartIter = 140000;

    uint windowIndex = 2;

    bool useConcreteHeatFluxBC = true;

    uint defaultDevice = 1;

    if( argc > 1 ) windowIndex = atoi( argv[1] );
    if( argc > 2 ) useConcreteHeatFluxBC = true;
    if( argc > 3 ) restartIter = atoi( argv[3] );

    //////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/RoomFireExtended/" );
#else
    std::string path( "out/" );
    
    //if( argc > 1 ){
    //    path += "Window_";
    //    path += argv[1];
    //    path += "/";
    //}
#endif

    std::string simulationName ( "RoomFire" );

    if( useConcreteHeatFluxBC ) simulationName += "_heatFlux";
    else                        simulationName += "_adiabatic";

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

    if( mpiWorldSize == 1 ) CudaUtility::setCudaDevice( 0 );
    else                    CudaUtility::setCudaDevice( rank % deviceCount );

    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    //////////////////////////////////////////////////////////////////////////

    try
    {
        thermalCavity( path, simulationName, windowIndex, restartIter, useConcreteHeatFluxBC );
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






real getHRR( real t )
{
    // data from 
    real tInMin_table [] = 
    { 0.0, 
      1.2998404845645375,     
      1.8225528293767326,     
      2.3411883091040355,     
      3.690242379336123,      
      5.8053588126615825,     
      8.481044195158887,      
      9.683816463416616,      
      10.268361262016242,     
      11.2607867055371,       
      13.013838692038146,     
      14.302516331727396,     
      17.240382966469404,     
      20.679801074868717,     
      22.9733288897661 };


    real HRR_table [] = 
    { 0.0,
      658.3729425582654,
      590.0388596425946,
      480.1207528610856,
      440.4722692284047,
      414.659889148097,
      406.6507906206217,
      374.9279268493922,
      337.28487256561004,
      260.02439647836513,
      141.15465878904172,
      85.66658361941495,
      51.906257987905406,
      33.97096366089556,
      27.954675614199346 };

    uint upper = 0;

    if( t / 60.0 > tInMin_table[14] ) return HRR_table[14];

    while( tInMin_table[upper] < t / 60.0 ) upper++;

    uint lower = upper - 1;

    real HRR = HRR_table[lower] + ( ( t / 60.0 - tInMin_table[lower] )/( tInMin_table[upper] - tInMin_table[lower] ) ) * ( HRR_table[upper] - HRR_table[lower] );

    return HRR;
}