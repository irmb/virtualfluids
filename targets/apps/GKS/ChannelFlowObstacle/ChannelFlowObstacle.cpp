//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
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
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"

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

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"

#include "Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void channelFlow( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 32+1;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    real L = 1.0;
    real H = 1.0;

    real dx = H / real(nx);

    real Re  = 1.0e4;
    real U  = 0.1;
    real Ma = 0.1;
    
    real Pr  = 0.1;
    real K   = 2.0;

    real rho = 1.0;

    real mu = U * rho * H / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( cs * cs );

    real g = eight * mu * U / ( H * H );

    real p0 = c1o2 * rho / lambda;

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.force.x = g;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = lambda;

    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //gridBuilder->addCoarseGrid(-0.5*L, -0.5*H, -0.5*dx,  
                                //0.5*L,  0.5*H,  0.5*dx, dx);

    gridBuilder->addCoarseGrid(-0.5*L, -0.5*H, -0.5*H,  
                                2.5*L,  0.5*H,  0.5*H, dx);

    Cuboid cube1(-0.1, -0.1, -0.1, 0.2, 0.1, 0.1);
    Cuboid cube2(-0.1, -0.1, -0.1, 0.2, 0.1, 0.1);

    gridBuilder->setNumberOfLayers(10,6);
    gridBuilder->addGrid( &cube1, 2);
    //gridBuilder->addGrid( &cube2, 3);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TriangularMesh* CubeSTL = TriangularMesh::make("F:/Work/Computations/inp/Cube.stl");

    gridBuilder->addGeometry(CubeSTL);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setPeriodicBoundaryCondition(true, false, true);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMX = std::make_shared<Pressure>( dataBase, p0 + g * 0.5 * L );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<Pressure>( dataBase, p0 - g * 0.5 * L );
    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0, true );
    SPtr<BoundaryCondition> bcPY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0, true );
    //SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    //SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*H; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*dx; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*dx; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcCube = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), 0.8*lambda, 0.0, true );

    bcCube->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return std::fabs(center.x) <  0.1 && 
                                                                           std::fabs(center.y) <  0.1 && 
                                                                           std::fabs(center.z) <  0.1; } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );

    dataBase->boundaryConditions.push_back( bcCube );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    if( false )
    {
        Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {

            real rhoLocal = rho;// - cellCenter.x * two * lambda * g;

            //real ULocal =0.0;//8.0 * ( ( 0.25 - cellCenter.y * cellCenter.y ) * ( 0.25 - cellCenter.z * cellCenter.z ) ) * U;

            real ULocal = four * (0.25 - cellCenter.y * cellCenter.y) * U;

            return toConservedVariables(PrimitiveVariables(rhoLocal, ULocal, 0.0, 0.0, lambda, 0.0), parameters.K);
        });

        dataBase->copyDataHostToDevice();

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );
    }
    else
    {
        Restart::readRestart(dataBase, path + simulationName + "_10000.rst", startIter );

        dataBase->copyDataHostToDevice();

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string(startIter) + "_restart" );
    }

    Initializer::initializeDataUpdate(dataBase);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 2000000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( iter % 10000 == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
        }

        cupsAnalyzer.run( iter );

        convergenceAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );


}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/out/ChannelFlowObstacle/" );
    //std::string path( "out/" );
    std::string simulationName ( "ChannelFlowObstacle" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        channelFlow( path, simulationName );
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
