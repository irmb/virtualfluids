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

#include "GksGpu/Analyzer/CupsAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void gksTest( std::string path )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 256;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L*2.0*M_PI / real(nx);

    real Re  = 1.6e3;
    real U  = 1.0;
    real Ma = 0.1;
    
    real Pr  = 0.71;
    real K   = 2.0;

    real rho = 1.0;

    //////////////////////////////////////////////////////////////////////////

    real gamma = ( K + 5 ) / ( K + 3 );

    real mu = U * rho * L / Re;

    real cs = U / Ma;
    real lambda = c1o2 * ( ( K + 4.0 ) / ( K + 2.0 ) ) / ( cs * cs );

    real CFL = 0.5;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    dt = 0.0025 * ( 32.0 / real(nx) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.force.x = 0;
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

    gridBuilder->addCoarseGrid(-0.5*L*2.0*M_PI, -0.5*L*2.0*M_PI, -0.5*L*2.0*M_PI,  
                                0.5*L*2.0*M_PI,  0.5*L*2.0*M_PI,  0.5*L*2.0*M_PI, dx);

    //gridBuilder->addCoarseGrid(-2.0 * dx, -0.5*L*2.0*M_PI, -0.5*L*2.0*M_PI,  
    //                            2.0 * dx,  0.5*L*2.0*M_PI,  0.5*L*2.0*M_PI, dx);

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

    //meshAdapter.writeMeshFaceVTK( path + "out/MeshFaces.vtk" );

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    auto dataBase = std::make_shared<DataBase>( "GPU" );
    dataBase->setMesh( meshAdapter );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMX = std::make_shared<Pressure>( dataBase, p0 + 1.5 * g * L );
    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.x < -0.5*L*2.0*M_PI;
    } );
    
    //SPtr<BoundaryCondition> bcPX = std::make_shared<Pressure>( dataBase, p0 - 1.5 * g * L );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );

    bcPX->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.x > 0.5*L*2.0*M_PI;
    } );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.y < -0.5*L*2.0*M_PI;
    } );

    //SPtr<BoundaryCondition> bcPY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.y > 0.5*L*2.0*M_PI;
    } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );

    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.z < -0.5*L*2.0*M_PI;
    } );
    
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0 );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );

    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ 
        return center.z > 0.5*L*2.0*M_PI;
    } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    //////////////////////////////////////////////////////////////////////////

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        //real A =  1.0;
        //real B =  1.0;
        //real C = -2.0;
        //real a = 2.0 * M_PI;
        //real b = 2.0 * M_PI;
        //real c = 2.0 * M_PI;

        //real ULocal = U * A * cos( a * cellCenter.x ) * sin( b * cellCenter.y ) * sin( c * cellCenter.z );
        //real VLocal = U * B * sin( a * cellCenter.x ) * cos( b * cellCenter.y ) * sin( c * cellCenter.z );
        //real WLocal = U * C * sin( a * cellCenter.x ) * sin( b * cellCenter.y ) * cos( c * cellCenter.z );

        real ULocal =   U * sin( cellCenter.x / L ) * cos( cellCenter.y / L ) * cos( cellCenter.z / L );
        real VLocal = - U * cos( cellCenter.x / L ) * sin( cellCenter.y / L ) * cos( cellCenter.z / L );
        real WLocal =   0.0;

        real pLocal = 1.0 / ( Ma * gamma ) + 1.0 / 16.0 * ( cos( 2.0 * cellCenter.x / L ) + cos( 2.0 * cellCenter.y / L ) ) * ( 2.0 + cos( 2.0 * cellCenter.z / L ) );

        real rhoLocal = 2.0 * pLocal * lambda;

        return toConservedVariables( PrimitiveVariables( rhoLocal, ULocal, VLocal, WLocal, lambda, 0.0 ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    writeVtkXML( dataBase, parameters, 0, path + "out/Test_0" );

    //////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0, false, 100 );

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 8000 * ( nx / 32 ); iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( iter % ( 400 * ( nx / 32 ) ) == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + "out/Test_" + std::to_string( iter ) );
        }

        cupsAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );


}

int main( int argc, char* argv[])
{
    //std::string path( "F:/Work/Computations/" );
    std::string path( "./" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    
    try
    {
        gksTest( path );
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
