//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
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

void init( uint rank, SPtr<DataBase> dataBase, SPtr<Parameters> parameters, std::string path, std::string simulationName )
{

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(rank % 2);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 64;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L / real(nx);

    real Re = 100.0;

    real U  = 1.0;
    real Ma = 0.1;
    
    real Pr  = 0.71;
    real K   = 2.0;

    real rho = 1.0;

    //////////////////////////////////////////////////////////////////////////

    real gamma = ( K + 5 ) / ( K + 3 );

    real mu = U * rho * L / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( cs * cs );

    real CFL = 0.5;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    //////////////////////////////////////////////////////////////////////////

    parameters->K  = K;
    parameters->Pr = Pr;
    parameters->mu = mu;

    parameters->force.x = 0;
    parameters->force.y = 0;
    parameters->force.z = 0;

    parameters->dt = dt;
    parameters->dx = dx;

    parameters->lambdaRef = lambda;

    parameters->viscosityModel = ViscosityModel::sutherlandsLaw;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( rank == 0 ) gridBuilder->addCoarseGrid(-0.5*L , -0.5*L , -0.5*L ,  
                                                3.0*dx,  3.0*dx,  3.0*dx, dx);

    if( rank == 1 ) gridBuilder->addCoarseGrid(-3.0*dx, -0.5*L , -0.5*L ,  
                                                0.5*L ,  3.0*dx,  3.0*dx, dx);

    if( rank == 2 ) gridBuilder->addCoarseGrid(-0.5*L , -3.0*dx, -0.5*L ,  
                                                3.0*dx,  0.5*L ,  3.0*dx, dx);

    if( rank == 3 ) gridBuilder->addCoarseGrid(-3.0*dx, -3.0*dx, -0.5*L ,  
                                                0.5*L ,  0.5*L ,  3.0*dx, dx);

    if( rank == 4 ) gridBuilder->addCoarseGrid(-0.5*L , -0.5*L , -3.0*dx,  
                                                3.0*dx,  3.0*dx,  0.5*L , dx);

    if( rank == 5 ) gridBuilder->addCoarseGrid(-3.0*dx, -0.5*L , -3.0*dx,  
                                                0.5*L ,  3.0*dx,  0.5*L , dx);

    if( rank == 6 ) gridBuilder->addCoarseGrid(-0.5*L , -3.0*dx, -3.0*dx,  
                                                3.0*dx,  0.5*L ,  0.5*L , dx);

    if( rank == 7 ) gridBuilder->addCoarseGrid(-3.0*dx, -3.0*dx, -3.0*dx,  
                                                0.5*L ,  0.5*L ,  0.5*L , dx);

    //////////////////////////////////////////////////////////////////////////

    Cuboid cube( -0.1, -0.1, -0.1, 
                  0.1,  0.1,  0.1 );
    
    gridBuilder->setNumberOfLayers(6,6);

    gridBuilder->addGrid(&cube, 1);

    //////////////////////////////////////////////////////////////////////////

    if( rank == 0 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -1.0, 0.0, 
                                                                                 -1.0, 0.0, 
                                                                                 -1.0, 0.0 ) );

    if( rank == 1 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.0, 1.0, 
                                                                                 -1.0, 0.0, 
                                                                                 -1.0, 0.0 ) );

    if( rank == 2 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -1.0, 0.0, 
                                                                                  0.0, 1.0, 
                                                                                 -1.0, 0.0 ) );

    if( rank == 3 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.0, 1.0, 
                                                                                  0.0, 1.0, 
                                                                                 -1.0, 0.0 ) );

    if( rank == 4 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -1.0, 0.0, 
                                                                                 -1.0, 0.0, 
                                                                                  0.0, 1.0 ) );

    if( rank == 5 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.0, 1.0, 
                                                                                 -1.0, 0.0, 
                                                                                  0.0, 1.0 ) );

    if( rank == 6 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -1.0, 0.0, 
                                                                                  0.0, 1.0, 
                                                                                  0.0, 1.0 ) );

    if( rank == 7 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.0, 1.0, 
                                                                                  0.0, 1.0, 
                                                                                  0.0, 1.0 ) );

    //////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(GKS, false);
            
    if( rank == 0 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::PX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PX, 1);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PY, 2);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PZ, 4);
    }
            
    if( rank == 1 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::MX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MX, 0);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PY, 3);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PZ, 5);
    }
            
    if( rank == 2 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::PX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PX, 3);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MY, 0);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PZ, 6);
    }
            
    if( rank == 3 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::MX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MX, 2);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MY, 1);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PZ, 7);
    }
            
    if( rank == 4 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::PX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PX, 5);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PY, 6);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MZ, 0);
    }
            
    if( rank == 5 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::MX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MX, 4);

        gridBuilder->findCommunicationIndices(CommunicationDirections::PY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PY, 7);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MZ, 1);
    }
            
    if( rank == 6 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::PX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::PX, 7);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MY, 4);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MZ, 2);
    }
            
    if( rank == 7 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::MX, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MX, 6);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MY, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MY, 5);

        gridBuilder->findCommunicationIndices(CommunicationDirections::MZ, GKS);
        gridBuilder->setCommunicationProcess (CommunicationDirections::MZ, 3);
    }

    gridBuilder->writeGridsToVtk(path + "Grid_" + std::to_string( rank ) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    meshAdapter.getCommunicationIndices();

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces_" + std::to_string( threadIndex ) + ".vtk" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, false );
    SPtr<BoundaryCondition> bcPX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, false );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, false );
    SPtr<BoundaryCondition> bcPY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, false );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, false );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3(  U,   U, 0.0), lambda, false );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*L; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    dataBase->setCommunicators( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{
        return toConservedVariables( PrimitiveVariables( rho, 0.0, 0.0, 0.0, lambda ), parameters->K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    //writeVtkXML( dataBase, *parameters, 0, path + simulationName + "_" + std::to_string( threadIndex ) + "_" + std::to_string( 0 ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run( uint rank, SPtr<DataBase> dataBase, SPtr<Parameters> parameters, std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice( rank % 2 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    writeVtkXML( dataBase, *parameters, 0, path + simulationName + "_" + std::to_string( rank ) + "_" + std::to_string( 0 ) );

    CupsAnalyzer cupsAnalyzer( dataBase, true, 300.0, true, 1000 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 1000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 100000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, *parameters, 0);

        if( 
            //( iter < 10     && iter % 1     == 0 ) ||
            //( iter < 100    && iter % 10    == 0 ) ||
            //( iter < 1000   && iter % 100   == 0 ) ||
            //( iter < 10000  && iter % 1000  == 0 ) 
            ( iter < 10000000 && iter % 1000 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, *parameters, 0, path + simulationName + "_" + std::to_string( rank ) + "_" + std::to_string( iter ) );
        }

        cupsAnalyzer.run( iter );

        convergenceAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();
}



int main( int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //////////////////////////////////////////////////////////////////////////

    std::string path( "F:/Work/Computations/out/" );
    //std::string path( "out/" );
    std::string simulationName ( "DrivenCavity" );
            
    std::ofstream logFile;
            
    logFile.open( path + simulationName + "_" + std::to_string(rank) + ".log" );

    logging::Logger::addStream(&logFile);

    logging::Logger::addStream(&std::cout);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        auto dataBase = std::make_shared<DataBase>( "GPU" );

        auto parameters = std::make_shared<Parameters>();

        init( rank, dataBase, parameters, path, simulationName);

        run ( rank, dataBase, parameters, path, simulationName);
    }
    catch (const std::exception& e)
    {     
        *logging::out << logging::Logger::ERROR << e.what() << "\n";
    }
    catch (const std::bad_alloc& e)
    {  
        *logging::out << logging::Logger::ERROR << "Bad Alloc:" << e.what() << "\n";
    }
    catch (...)
    {
        *logging::out << logging::Logger::ERROR << "Unknown exception!\n";
    }

    logFile.close();

    MPI_Finalize();

    return 0;
}
