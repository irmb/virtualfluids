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
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"
#include "GksGpu/BoundaryConditions/Inflow.h"
#include "GksGpu/BoundaryConditions/Extrapolation.h"
#include "GksGpu/BoundaryConditions/Pressure.h"
#include "GksGpu/BoundaryConditions/Periodic.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 128;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;
    //real H = 0.25;
    real H = L / real(nx);

    real dx = L / real(nx);


    real Ra = 2.0e9;

    real Ba  = 0.1;
    real eps = 1.2;
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

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " s\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.force.x = 0.125 * g;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    parameters.lambdaRef = lambda;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //gridBuilder->addCoarseGrid(-0.5*L, -0.5*L, -0.5*H,  
    //                            0.5*L,  0.5*L,  0.5*H, dx);

    gridBuilder->addCoarseGrid(-    L,  0.0  , -0.5*H,  
                                4.0*L,  1.0*L,  0.5*H, dx);

    real L_1 = ( 0.5 - 0.35 ) / 2.0;
    real L_2 = ( 0.5 - 0.45 ) / 2.0;
    real L_3 = ( 0.5 - 0.475) / 2.0;
    real L_4 = ( 0.5 - 0.485) / 2.0;

    Cuboid* cubeMY_1 = new Cuboid (-2.0, -2.0, -2.0, 
                                    1.9,  L_1,  2.0 );

    Cuboid* cubeMY_2 = new Cuboid (-2.0, -2.0, -2.0, 
                                    1.8,  L_2,  2.0 );

    Cuboid* cubeMY_3 = new Cuboid (-2.0, -2.0, -2.0, 
                                    5.0,  L_3,  2.0 );

    Cuboid* cubeMY_4 = new Cuboid (-2.0, -2.0, -2.0, 
                                    2.0,  L_4,  2.0 );

    Conglomerate refRegion_1;
    refRegion_1.add(cubeMY_1);

    Conglomerate refRegion_2;
    refRegion_2.add(cubeMY_2);

    Conglomerate refRegion_3;
    refRegion_3.add(cubeMY_3);

    Conglomerate refRegion_4;
    refRegion_4.add(cubeMY_4);

    gridBuilder->setNumberOfLayers(6,6);

    gridBuilder->addGrid( &refRegion_1, 1);
    gridBuilder->addGrid( &refRegion_2, 2);
    //gridBuilder->addGrid( &refRegion_3, 3);
    //gridBuilder->addGrid( &refRegion_4, 4);

    gridBuilder->setPeriodicBoundaryCondition(false, false, true);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);  // heated
    //CudaUtility::setCudaDevice(1);  // cooled

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real inletHeight = 0.02;

    SPtr<BoundaryCondition> bcMX_1 = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0, false );
    SPtr<BoundaryCondition> bcMX_2 = std::make_shared<Inflow>( dataBase, Vec3(0.2, 0.0, 0.0), lambda, rho, 0.0, 0.0, inletHeight, -1.0 );
    //SPtr<BoundaryCondition> bcMX_3 = std::make_shared<Pressure>( dataBase, 0.5 * rho / lambda );

    SPtr<BoundaryCondition> bcPX   = std::make_shared<Pressure>( dataBase, 0.5 * rho / lambda );

    bcMX_1->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -L && center.y > inletHeight; } );
    bcMX_2->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -L && center.y < inletHeight; } );
    //bcMX_2->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.x < -L; } );
    bcPX->findBoundaryCells(   meshAdapter, true, [&](Vec3 center){ return center.x >  4.0*L; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaHot , 0.0, false );
    //SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambdaCold, 0.0, false );
    //SPtr<BoundaryCondition> bcMY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0, false );

    SPtr<BoundaryCondition> bcPY = std::make_shared<Extrapolation>( dataBase );
    //SPtr<BoundaryCondition> bcPY = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), lambda, 0.0, false );

    bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y <  0.0  ; } );
    bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  L    ; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < -0.5*H; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z >  0.5*H; } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX_1 );
    dataBase->boundaryConditions.push_back( bcMX_2 );
    //dataBase->boundaryConditions.push_back( bcMX_3 );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();

    Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables{
        
        real y = cellCenter.y;

        real factor = ( 0.0 
                      + inletHeight*y 
                      - 1.0  *y*y  ) * ( four / inletHeight / inletHeight );

        real U_local;
        if( y < inletHeight )
            U_local = 0.2 * factor;
        else
            U_local = 0.0;

        return toConservedVariables( PrimitiveVariables( rho, U_local, 0.0, 0.0, lambda/*, 0.0*/ ), parameters.K );
    });

    dataBase->copyDataHostToDevice();

    Initializer::initializeDataUpdate(dataBase);

    for( uint level = 0; level < dataBase->numberOfLevels; level++ )
    {
        for (SPtr<BoundaryCondition> bc : dataBase->boundaryConditions) {
            bc->runBoundaryConditionKernel(dataBase, parameters, level);
        }
    }

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

    for( uint iter = 1; iter <= 1000000; iter++ )
    {
        TimeStepping::nestedTimeStep(dataBase, parameters, nullptr, 0);

        if( 
            //( iter < 10     && iter % 1     == 0 ) ||
            //( iter < 100    && iter % 10    == 0 ) ||
            //( iter < 1000   && iter % 100   == 0 ) ||
            ( iter < 100000  && iter % 1000  == 0 ) ||
            ( iter < 10000000 && iter % 10000 == 0 )
          )
        {
            for( uint level = 0; level < dataBase->numberOfLevels; level++ )
            {
                for (SPtr<BoundaryCondition> bc : dataBase->boundaryConditions) {
                    bc->runBoundaryConditionKernel(dataBase, parameters, level);
                }
            }
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
    std::string path( "F:/Work/Computations/out/BoundaryJet/Heated/" );
    //std::string path( "F:/Work/Computations/out/BoundaryJet/Cooled/" );
    //std::string path( "out/" );
    std::string simulationName ( "BoundaryJet" );

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
