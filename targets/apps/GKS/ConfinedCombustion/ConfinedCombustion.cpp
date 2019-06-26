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

#include "GksGpu/FlowStateData/FlowStateData.cuh"
#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"
#include "GksGpu/FlowStateData/ThermalDependencies.cuh"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"
#include "GksGpu/BoundaryConditions/Periodic.h"
#include "GksGpu/BoundaryConditions/Pressure.h"
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"

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

    real L = 1.0;

    real dx = L / real(nx);

    real U = 0.025;

    real eps = 2.0;
    real Pr  = 0.71;
    real K   = 5.0;
    
    real g   = 9.81;
    real rho = 1.2;
    
    real mu = 5.0e-1;

    PrimitiveVariables prim( rho, 0.0, 0.0, 10.0, -1.0 );

    setLambdaFromT( prim, 3.0 / T_FAKTOR );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.D  = mu;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid(-0.5*L, -0.5*dx, -0.5*dx,  
                                0.5*L,  0.5*dx,  0.5*dx, dx);

    //gridBuilder->addCoarseGrid(-0.5*L, -0.5*L, -0.5*L,  
    //                            0.5*L,  0.5*L,  0.5*L, dx);

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

    CudaUtility::setCudaDevice(1);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );

    //SPtr<BoundaryCondition> bcMX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold, false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, false );

    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> bcMY = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );
    //SPtr<BoundaryCondition> bcPY = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), false );

    //SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, false );
    //SPtr<BoundaryCondition> bcPY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, false );

    SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    //bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.125*L; } );
    //bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.125*L; } );

    bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*dx; } );
    bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*dx; } );

    //////////////////////////////////////////////////////////////////////////
    
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0.0, 0.0, 0.0), true );

    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaHot, false );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold,  0.0, true );

    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    
    //bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.125*L; } );
    //bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.125*L; } );

    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*dx; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*dx; } );

    //////////////////////////////////////////////////////////////////////////
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{

        //PrimitiveVariables primFuel = prim;

        //primFuel.S_1 = 1.0;

        //////////////////////////////////////////////////////////////////////////

        //PrimitiveVariables primAir = prim;

        //////////////////////////////////////////////////////////////////////////

        //real massFuel = 1.0; 
        //real massAir  = 2.0 * 32.0/16.0 + 2.0 * 0.767 / 0.233 * 32.0/16.0;

        //real volumeFuel = massFuel / primFuel.rho;
        //real volumeAir  = massAir  / primAir.rho;
        //
        //real volumeRatioFuel = volumeFuel / ( volumeFuel + volumeAir );

        //if(fabs(cellCenter.x) < 0.5 * volumeRatioFuel ) return toConservedVariables( primFuel, parameters.K );
        //else                                            return toConservedVariables( primAir , parameters.K );

        //////////////////////////////////////////////////////////////////////////

        //PrimitiveVariables primMix = prim;

        //primMix.S_1 = volumeRatioFuel;

        //return toConservedVariables( primMix, parameters.K );

        //////////////////////////////////////////////////////////////////////////

        if( nx == 1 )
        {
            double X_F = 0.21 / 2.21;
            double X_A = 1.0 - X_F;

            double M = X_F * M_F + X_A * M_A;

            double Y_F = X_F * M_F / M;
            double Y_A = X_A * M_A / M;

            prim.S_1 = Y_F;

            return toConservedVariables(prim, parameters.K);
        }

        //////////////////////////////////////////////////////////////////////////

        if( nx > 1 )
        {
            double X_F = 0.21 / 2.21;
            double X_A = 1.0 - X_F;

            double M = X_F * M_F + X_A * M_A;

            double Y_F = X_F * M_F / M;
            double Y_A = X_A * M_A / M;

            if (cellCenter.x < 0) prim.S_1 = 0.0;
            else                  prim.S_1 = 2.0 * Y_F;

            return toConservedVariables(prim, parameters.K);
        }

        //////////////////////////////////////////////////////////////////////////

        //if( nx > 1 )
        //{
        //    double X_F = 1.0;
        //    double X_A = 1.0 - X_F;

        //    double M = X_F * M_F + X_A * M_A;

        //    double Y_F = X_F * M_F / M;
        //    double Y_A = X_A * M_A / M;

        //    if (cellCenter.x < 0) prim.S_1 = 0.0;
        //    else                  prim.S_1 = Y_F;

        //    if (cellCenter.x < 0) prim.rho = 1.2;
        //    else                  prim.rho = 0.1;

        //    return toConservedVariables(prim, parameters.K);
        //}
    });

    //std::cout << toConservedVariables( PrimitiveVariables( rho, 0.0, 0.0, 0.0, lambdaHot, S_1, S_2 ), parameters.K ).rhoE << std::endl;

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

    ConvergenceAnalyzer convergenceAnalyzer( dataBase );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = 1; iter <= 10000; iter++ )
    {
        //if( iter < 100000 )
        //{
        //    std::dynamic_pointer_cast<IsothermalWall>(bcMX)->lambda = lambdaCold + ( lambdaHot - lambdaCold ) * ( real(iter) / 100000.0 );
        //}
        //if( iter == 100000 )
        //{
        //    //std::dynamic_pointer_cast<IsothermalWall>(bcMX)->lambda = lambdaHot;
        //    dataBase->boundaryConditions[4] = bcMX_2;
        //}

        cupsAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( 
            //( iter < 10       && iter % 1     == 0 ) ||
            //( iter < 100      && iter % 10    == 0 ) ||
            //( iter < 1000     && iter % 100   == 0 ) ||
            //( iter < 10000    && iter % 1000  == 0 ) ||
            //( iter < 1000000   && iter % 10000  == 0 )
            //( iter < 10000000 && iter % 100000 == 0 )
            //( iter <= 400000 && iter % 100 == 0 )
            ( iter % 100 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

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

#ifdef _WIN32
    std::string path( "F:/Work/Computations/out/ConfinedCombustion/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "ConfinedCombustion" );

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
