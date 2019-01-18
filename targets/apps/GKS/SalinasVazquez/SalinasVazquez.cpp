//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>
#include <thread>

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
#include "GksGpu/BoundaryConditions/SalinasVazquez.h"

#include "GksGpu/Communication/Communicator.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

//uint deviceMap [2] = {2,3};
uint deviceMap [2] = {0,1};

void init( uint threadIndex, SPtr<DataBase> dataBase, SPtr<Communicator> communicator, SPtr<Parameters> parameters, std::string path, std::string simulationName )
{

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(deviceMap[threadIndex]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 128;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;
    real H = 0.25;

    real dx = L / real(nx);

    real Ra = 1.58e9;

    real Ba  = 0.1;
    real eps = 0.132;
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

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    //////////////////////////////////////////////////////////////////////////

    parameters->K  = K;
    parameters->Pr = Pr;
    parameters->mu = mu;

    parameters->force.x = 0;
    parameters->force.y = -g;
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

    Conglomerate refRegion_1;
    Conglomerate refRegion_2;
    Conglomerate refRegion_3;
    Conglomerate refRegion_4;

    real L_1 = 0.35;
    real L_2 = 0.45;
    real L_3 = 0.475;
    real L_4 = 0.495;

    if( threadIndex == 0 ) gridBuilder->addCoarseGrid(-0.5*L , -0.5*L, -0.5*H,  
                                                       3.0*dx,  0.5*L,  0.5*H, dx);

    if( threadIndex == 1 ) gridBuilder->addCoarseGrid(-3.0*dx, -0.5*L, -0.5*H,  
                                                       0.5*L ,  0.5*L,  0.5*H, dx);

    if( threadIndex == 0 ) refRegion_1.add( new Cuboid (-1.0, -1.0, -1.0, 
                                                        -L_1,  1.0,  1.0 ) );

    if( threadIndex == 1 ) refRegion_1.add( new Cuboid ( L_1, -1.0, -1.0, 
                                                         1.0,  1.0,  1.0 ) );

    if( threadIndex == 0 ) refRegion_2.add( new Cuboid (-1.0, -1.0, -1.0, 
                                                        -L_2,  1.0,  1.0 ) );

    if( threadIndex == 1 ) refRegion_2.add( new Cuboid ( L_2, -1.0, -1.0, 
                                                         1.0,  1.0,  1.0 ) );

    if( threadIndex == 0 ) refRegion_3.add( new Cuboid (-1.0, -1.0, -1.0, 
                                                        -L_3,  1.0,  1.0 ) );

    if( threadIndex == 1 ) refRegion_3.add( new Cuboid ( L_3, -1.0, -1.0, 
                                                         1.0,  1.0,  1.0 ) );

    if( threadIndex == 0 ) refRegion_4.add( new Cuboid (-1.0, -1.0, -1.0, 
                                                        -L_4,  1.0,  1.0 ) );

    if( threadIndex == 1 ) refRegion_4.add( new Cuboid ( L_4, -1.0, -1.0, 
                                                         1.0,  1.0,  1.0 ) );

    gridBuilder->setNumberOfLayers(6,6);
    gridBuilder->addGrid( &refRegion_1, 1);
    gridBuilder->addGrid( &refRegion_2, 2);
    gridBuilder->addGrid( &refRegion_3, 3);
    //gridBuilder->addGrid( &refRegion_4, 4);

    if( threadIndex == 0 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( -1.0, 0.0, 
                                                                                        -1.0, 1.0, 
                                                                                        -1.0, 1.0 ) );

    if( threadIndex == 1 ) gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>(  0.0, 1.0, 
                                                                                        -1.0, 1.0, 
                                                                                        -1.0, 1.0 ) );

    gridBuilder->setPeriodicBoundaryCondition(false, false, true);

    gridBuilder->buildGrids(GKS, false);
            
    if( threadIndex == 0 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::PX, GKS);
    }
            
    if( threadIndex == 1 ){
        gridBuilder->findCommunicationIndices(CommunicationDirections::MX, GKS);
    }

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_" + std::to_string( threadIndex ) + "_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    meshAdapter.findPeriodicBoundaryNeighbors();    

    meshAdapter.getCommunicationIndices();

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces_" + std::to_string( threadIndex ) + ".vtk" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0) );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0) );
    SPtr<BoundaryCondition> bcMX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaHot , false );
    SPtr<BoundaryCondition> bcPX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold, false );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMY = std::make_shared<SalinasVazquez>( dataBase, lambdaHot, lambdaCold, 0.3371, -0.2642,  0.5301, -2.6438, true );
    SPtr<BoundaryCondition> bcPY = std::make_shared<SalinasVazquez>( dataBase, lambdaHot, lambdaCold, 0.6559, -0.2037, -0.5420, -2.7318, true );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0) );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0) );
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*H; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*H; } );

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

    if( threadIndex == 0 )
        communicator->initialize( meshAdapter, 1 );
    
    if( threadIndex == 1 )
        communicator->initialize( meshAdapter, 0 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        real Th = 1.0 / lambdaHot;
        real Tc = 1.0 / lambdaCold;
        real T = Th - (Th - Tc)*( (cellCenter.x + 0.5 * L) / L);
        real lambdaLocal = 1.0 / T;

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

void run( uint threadIndex, SPtr<DataBase> dataBase, SPtr<Communicator> communicator, SPtr<Parameters> parameters, std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(deviceMap[threadIndex]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    writeVtkXML( dataBase, *parameters, 0, path + simulationName + "_" + std::to_string( threadIndex ) + "_" + std::to_string( 0 ) );

    CupsAnalyzer cupsAnalyzer( dataBase, true, 300.0 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase );

    auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 1000000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, *parameters, communicator, 0);

        if( 
            //( iter < 10     && iter % 1     == 0 ) ||
            //( iter < 100    && iter % 10    == 0 ) ||
            //( iter < 1000   && iter % 100   == 0 ) ||
            //( iter < 10000  && iter % 1000  == 0 ) 
            ( iter < 10000000 && iter % 10000 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, *parameters, 0, path + simulationName + "_" + std::to_string( threadIndex ) + "_" + std::to_string( iter ) );
        }

        cupsAnalyzer.run( iter );

        convergenceAnalyzer.run( iter );

        //turbulenceAnalyzer->run( iter, *parameters );

        //if( iter % 50000 == 0 )
        //{
        //    turbulenceAnalyzer->download();

        //    writeTurbulenceVtkXML(dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence_" + std::to_string( iter ));
        //}
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();
}

int main( int argc, char* argv[])
{
    //std::string path( "F:/Work/Computations/out/" );
    std::string path( "out/" );
    std::string simulationName ( "SalinasVazquez" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        auto dataBase_0 = std::make_shared<DataBase>( "GPU" );
        auto dataBase_1 = std::make_shared<DataBase>( "GPU" );

        auto communicator_0 = std::make_shared<Communicator>( dataBase_0 );
        auto communicator_1 = std::make_shared<Communicator>( dataBase_1 );

        communicator_0->opposingCommunicator = communicator_1;
        communicator_1->opposingCommunicator = communicator_0;

        auto parameters_0 = std::make_shared<Parameters>();
        auto parameters_1 = std::make_shared<Parameters>();

        {
            std::thread thread_0(init, 0, dataBase_0, communicator_0, parameters_0, path, simulationName);
            std::thread thread_1(init, 1, dataBase_1, communicator_1, parameters_1, path, simulationName);

            thread_0.join();
            thread_1.join();
        }

        {
            std::thread thread_0(run, 0, dataBase_0, communicator_0, parameters_0, path, simulationName);
            std::thread thread_1(run, 1, dataBase_1, communicator_1, parameters_1, path, simulationName);

            thread_0.join();
            thread_1.join();
        }

        //writeVtkXML( dataBase_0, *parameters_0, 0, path + simulationName + "_" + std::to_string( 0 ) + "_" + std::to_string( 1 ) );
        //writeVtkXML( dataBase_1, *parameters_1, 0, path + simulationName + "_" + std::to_string( 1 ) + "_" + std::to_string( 1 ) );
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

   return 0;
}
