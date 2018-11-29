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

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "GksVtkAdapter/VTKInterface.h"

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"
#include "GksGpu/Initializer/Initializer.h"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"
#include "GksGpu/BoundaryConditions/Pressure.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/FlowStateData/AccessDeviceData.cuh"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void channelFlow( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 128;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;
    real H = 1.5;

    real dx = L / real(nx);

    real Re  = 1.0e2;
    real U  = 0.1;
    real Ma = 0.1;
    
    real Pr  = 0.1;
    real K   = 2.0;

    real rho = 1.0;

    real mu = U * rho * L / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( cs * cs );

    real g = eight * mu * U / ( H * H );

    real p0 = c1o2 * rho / lambda;

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    dt = 0.001 * ( 32.0 / real(nx) );

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

    gridBuilder->addCoarseGrid(-0.5*L, -0.5*H, -0.5*dx,  
                                0.5*L,  0.5*H,  0.5*dx, dx);

    //Cuboid cube(-1.0, -1.0, 0.45, 1.0, 1.0, 0.55);

    //gridBuilder->setNumberOfLayers(6,6);
    //gridBuilder->addGrid( &cube, 1);

    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

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
    
    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*H; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*dx; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*dx; } );

    //////////////////////////////////////////////////////////////////////////
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        real ULocal = 0.0;
        real VLocal = 0.0;
        real WLocal = U * sin( 2.0 * M_PI * cellCenter.x / L ) * cos( 4.0 / 3.0 * M_PI * cellCenter.y / L );

        return toConservedVariables( PrimitiveVariables( rho, ULocal, VLocal, WLocal, lambda, 0.0 ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    writeVtkXML( dataBase, parameters, 0, path + simulationName + ".0." + std::to_string( nx ) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 10000 * ( nx / 32 ); iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        //if( iter % ( 100 * ( nx / 32 ) ) == 0 )
        //{
        //    dataBase->copyDataDeviceToHost();

        //    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        //}

        cupsAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    writeVtkXML( dataBase, parameters, 0, path + simulationName + ".Result." + std::to_string( nx ) );


    {
        real sum1 = 0.0;
        real sum2 = 0.0;
    
        for( uint cellIdx = 0; cellIdx < dataBase->perLevelCount[0].numberOfBulkCells; cellIdx++ )
        {
            real velocity = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ] / dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];

            sum1 += velocity * U * sin( 2.0 * M_PI * dataBase->getCellCenter(cellIdx).x / L ) * cos( 4.0 / 3.0 * M_PI * dataBase->getCellCenter(cellIdx).y / L ); 
            sum2 +=            U * sin( 2.0 * M_PI * dataBase->getCellCenter(cellIdx).x / L ) * cos( 4.0 / 3.0 * M_PI * dataBase->getCellCenter(cellIdx).y / L )
                  *            U * sin( 2.0 * M_PI * dataBase->getCellCenter(cellIdx).x / L ) * cos( 4.0 / 3.0 * M_PI * dataBase->getCellCenter(cellIdx).y / L ); 
        }

        *logging::out << logging::Logger::INFO_HIGH << sum1 / sum2 << "\n";
    }
}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/out/" );
    std::string simulationName ( "ShearWave" );

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
