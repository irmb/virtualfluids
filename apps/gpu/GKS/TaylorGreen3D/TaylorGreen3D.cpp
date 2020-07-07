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
#include "GksGpu/Analyzer/KineticEnergyAnalyzer.h"
#include "GksGpu/Analyzer/EnstrophyAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void writeVelocityFile( SPtr<DataBase> dataBase, std::string filename );

//////////////////////////////////////////////////////////////////////////
real Re = 1.6e3;

uint dtPerL = 500;

uint nx = 64;
uint gpuIndex = 0;
//////////////////////////////////////////////////////////////////////////

void gksTest( std::string path )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice( gpuIndex );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = 2.0 * M_PI * L / real(nx);

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

    *logging::out << logging::Logger::INFO_HIGH << "dt(CFL=0.5) = " << dt << " s\n";

    dt = L / U /  dtPerL * ( 64.0 / real(nx) );

    *logging::out << logging::Logger::INFO_HIGH << "dt          = " << dt << " s\n";

    *logging::out << logging::Logger::INFO_HIGH << "mu          = " << mu << "\n";

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

    gridBuilder->addCoarseGrid(-M_PI*L, -M_PI*L, -M_PI*L,  
                                M_PI*L,  M_PI*L,  M_PI*L, dx);

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

    auto dataBase = std::make_shared<DataBase>( "GPU" );
    dataBase->setMesh( meshAdapter );

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

    bcMY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );

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

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        real ULocal =   U * sin( cellCenter.x ) * cos( cellCenter.y ) * cos( cellCenter.z );
        real VLocal = - U * cos( cellCenter.x ) * sin( cellCenter.y ) * cos( cellCenter.z );
        real WLocal =   0.0;

        real p0 = 0.5 * rho / lambda;

        real pLocal = p0 + rho * U * U / 16.0 * ( cos( 2.0 * cellCenter.x ) + cos( 2.0 * cellCenter.y ) ) * ( 2.0 + cos( 2.0 * cellCenter.z ) );

        real rhoLocal = 2.0 * pLocal * lambda;

        return toConservedVariables( PrimitiveVariables( rhoLocal, ULocal, VLocal, WLocal, lambda ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    writeVtkXML( dataBase, parameters, 0, path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_"          + std::to_string( 0 ) );
    writeVelocityFile( dataBase,          path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_Velocity_" + std::to_string( 0 ) );

    //////////////////////////////////////////////////////////////////////////

    KineticEnergyAnalyzer kineticEnergyAnalyzer( dataBase,             10, 10000 );
    EnstrophyAnalyzer     enstrophyAnalyzer    ( dataBase, parameters, 10, 10000 );

    CupsAnalyzer cupsAnalyzer( dataBase, true, 60.0, false, 100 );

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 40 * lround(L/(U*dt)); iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        kineticEnergyAnalyzer.run( iter );
        enstrophyAnalyzer.run( iter );

        if( iter % ( 5 * lround(L/(U*dt)) ) == 0 )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_"          + std::to_string( iter / lround(L/(U*dt)) ) );
            writeVelocityFile( dataBase,          path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_Velocity_" + std::to_string( iter / lround(L/(U*dt)) ) );
            kineticEnergyAnalyzer.writeToFile(    path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_EKin_"     + std::to_string( iter / lround(L/(U*dt)) ) );
            enstrophyAnalyzer.writeToFile    (    path + "TGV_3D_nx_" + std::to_string(nx) + "_dtPerL_" + std::to_string(dtPerL) + "_Enstrophy_"+ std::to_string( iter / lround(L/(U*dt)) ) );
        }

        cupsAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////
}

int main( int argc, char* argv[])
{
    if( argc > 1 ) gpuIndex = atoi( argv[1] );
    if( argc > 2 ) Re       = atof( argv[2] );
    if( argc > 3 ) nx       = atoi( argv[3] );
    if( argc > 4 ) dtPerL   = atoi( argv[4] );

    //////////////////////////////////////////////////////////////////////////

    //std::string path( "F:/Work/Computations/TaylorGreenVortex_3D/results/GKS/" );
    std::string path( "./results/GKS/" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    //////////////////////////////////////////////////////////////////////////

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "FlowStateData/AccessDeviceData.cuh"

void writeVelocityFile( SPtr<DataBase> dataBase, std::string filename )
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "writeVelocityFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + ".dat" );

    for( uint cellIndex = 0; cellIndex < dataBase->perLevelCount[0].numberOfBulkCells; cellIndex++ )
    {
        real rho = dataBase->dataHost[ RHO__(cellIndex, dataBase->numberOfCells) ];

        file << dataBase->dataHost[ RHO_U(cellIndex, dataBase->numberOfCells) ] / rho << ", ";
        file << dataBase->dataHost[ RHO_V(cellIndex, dataBase->numberOfCells) ] / rho << ", ";
        file << dataBase->dataHost[ RHO_W(cellIndex, dataBase->numberOfCells) ] / rho << std::endl;
    }

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}