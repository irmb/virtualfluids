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
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/VerticalCylinder/VerticalCylinder.h"
#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"

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
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"
#include "GksGpu/BoundaryConditions/Extrapolation.h"
#include "GksGpu/BoundaryConditions/Inflow.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 64;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real L = 3.0;
    real H = 1.0;

    real dx = H / real(nx);

    real Re  = 500.0;
    real Ba  = 0.1;
    real eps = 1.2;
    real Pr  = 0.71;
    real K   = 8.0;
    
    real U   = 1.0;

    real g   = 9.81;
    real rho = 1.2;

    real S_1 = 0.0;
    real S_2 = 0.5;

    real R_Mixture = S_1               * 8.31445984848 / 16.04e-3      // O2
				   + S_2               * 8.31445984848 / 32.00e-3      // CH4
		           + (1.0 - S_1 - S_2) * 8.31445984848 / 28.00e-3;     // N2

    real lambdaCold = 0.5 / ( R_Mixture *  300 ) * 1000.0;
    real lambdaHot  = 0.5 / ( R_Mixture * 1200 ) * 1000.0;
    
    real mu = U * rho * 0.25 * H / Re;

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * lambdaCold ) );

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt   << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "Ma = " << U/cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu   << " kg/sm\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.D  = 5.0 * mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = lambdaCold;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //gridBuilder->addCoarseGrid( 0.0, -0.5*H, -0.5*dx,  
    //                              L,  0.5*H,  0.5*dx, dx );

    gridBuilder->addCoarseGrid( 0.0, -0.5*H, -0.5*H,  
                                  L,  0.5*H,  0.5*H, dx );

    Sphere           sphere  ( 0.0, 0.0, 0.0, 0.15 );

    gridBuilder->setNumberOfLayers(0,10);

    gridBuilder->addGrid( &sphere, 2 );

    gridBuilder->setPeriodicBoundaryCondition(false, true, true);

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
    
    SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Extrapolation>( dataBase );



    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < 0.0; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >   L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*H; } );
    bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );

    //bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*dx; } );
    //bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*dx; } );
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*H; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcJetFuel   = std::make_shared<Inflow>( dataBase, Vec3(U, 0.0, 0.0), lambdaHot, rho, 1.0, 0.0, -64.0, 0.0, 0.5 );
    SPtr<BoundaryCondition> bcJetOxygen = std::make_shared<Inflow>( dataBase, Vec3(U, 0.0, 0.0), lambdaHot, rho, 1.0, 0.0, -64.0, 0.0, 0.5 );

    bcJetFuel->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 
        return center.x < 0.0 &&
               std::sqrt(center.y*center.y + center.z*center.z) < 0.125 / 4.0;
    } );

    bcJetOxygen->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 
        return center.x < 0.0 &&
               std::sqrt(center.y*center.y + center.z*center.z) < 0.125;
    } );

    //////////////////////////////////////////////////////////////////////////
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );

    dataBase->boundaryConditions.push_back( bcJetOxygen );
    dataBase->boundaryConditions.push_back( bcJetFuel );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        real rhoLocal = rho;
        real lambdaLocal = lambdaCold;
        //real lambdaLocal = lambdaHot;

        //real radius = sqrt( cellCenter.x*cellCenter.x + cellCenter.y*cellCenter.y + cellCenter.z*cellCenter.z );

        //if( radius < 0.2 )
        //{
        //    lambdaLocal = lambdaHot;
        //}

        //lambdaLocal = lambdaCold + ( lambdaHot - lambdaCold ) * exp( - 10. * ( cellCenter.x*cellCenter.x + cellCenter.y*cellCenter.y + cellCenter.z*cellCenter.z ) );

        //lambdaLocal = lambdaCold + ( lambdaHot - lambdaCold ) * ( 0.5 * M_PI + atan( - 1000.0 * ( radius - 0.1) ) ) / M_PI;

        //rhoLocal = rho * lambdaLocal / lambdaCold;

        //lambdaLocal = lambdaCold + ( lambdaHot - lambdaCold ) * exp( - 10. * ( (cellCenter.x-0.5)*(cellCenter.x-0.5) ) );

        real radius = sqrt( cellCenter.y*cellCenter.y + cellCenter.z*cellCenter.z );
        
        real factor = 0.0;
        //if( radius < 0.125 ) factor = ( 1.0 - 64.0 * radius * radius  );

        return toConservedVariables( PrimitiveVariables( rhoLocal, factor * U, 0.0, 0.0, lambdaLocal, S_1, S_2 ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    for( auto bc : dataBase->boundaryConditions ) 
        for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            bc->runBoundaryConditionKernel( dataBase, parameters, 0 );

    Initializer::initializeDataUpdate(dataBase);

    dataBase->copyDataDeviceToHost();

    writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 1000 );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 100000000; iter++ )
    {
        uint T = 20000;
        if( iter <= T )
        {
            std::dynamic_pointer_cast<Inflow>(bcJetFuel  )->lambda = lambdaCold + ( lambdaHot - lambdaCold ) * ( real(iter) / real(T) );
            std::dynamic_pointer_cast<Inflow>(bcJetOxygen)->lambda = lambdaCold + ( lambdaHot - lambdaCold ) * ( real(iter) / real(T) );
        }

        if( iter == T )
        {
            std::dynamic_pointer_cast<Inflow>(bcJetFuel)->S_1 = 1.0;
            std::dynamic_pointer_cast<Inflow>(bcJetFuel)->S_2 = 0.0;
        }

        TimeStepping::nestedTimeStep(dataBase, parameters, nullptr, 0);

        if( 
            //( iter < 10       && iter % 1      == 0 ) ||
            //( iter < 100      && iter % 10     == 0 ) ||
            //( iter < 1000     && iter % 100    == 0 ) ||
            ( iter < 20000    && iter % 1000   == 0 ) ||
            ( iter < 1000000  && iter % 10000  == 0 ) ||
            ( iter < 100000000  && iter % 100000 == 0 )
            //( iter > 18400 && iter % 10 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        cupsAnalyzer.run( iter );

        convergenceAnalyzer.run( iter );

        //turbulenceAnalyzer->run( iter, parameters );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    //writeVtkXML( dataBase, parameters, 0, path + "grid/Test_1" );

    //turbulenceAnalyzer->download();

    //writeTurbulenceVtkXML(dataBase, turbulenceAnalyzer, 0, path + simulationName + "_Turbulence");
}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/out/MethaneFlame/" );
    //std::string path( "out/" );
    std::string simulationName ( "MethaneFlame" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precison\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        thermalCavity( path, simulationName );
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
