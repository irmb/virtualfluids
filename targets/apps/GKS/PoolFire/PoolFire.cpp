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
#include "GksGpu/BoundaryConditions/Pressure2.h"
#include "GksGpu/BoundaryConditions/AdiabaticWall.h"
#include "GksGpu/BoundaryConditions/PassiveScalarDiriclet.h"
#include "GksGpu/BoundaryConditions/InflowComplete.h"
#include "GksGpu/BoundaryConditions/Open.h"
#include "GksGpu/BoundaryConditions/Extrapolation.h"
#include "GksGpu/BoundaryConditions/Symmetry.h"
#include "GksGpu/BoundaryConditions/CreepingMassFlux.h"
#include "GksGpu/BoundaryConditions/MassCompensation.h"

#include "GksGpu/Interface/Interface.h"
#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
#include "GksGpu/Analyzer/TurbulenceAnalyzer.h"

#include "GksGpu/Restart/Restart.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

void thermalCavity( std::string path, std::string simulationName, uint restartIter )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 256;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 4.0;
    real H = 4.0;
    real W = 0.125;

    real dx = H / real(nx);

    real U = 0.025;

    real eps = 2.0;
    real Pr  = 0.71;
    real K   = 5.0;
    
    real g   = 9.81;
    real rho = 1.2;
    
    real mu = 5.0e-4;

    PrimitiveVariables prim( rho, 0.0, 0.0, 0.0, -1.0 );

    setLambdaFromT( prim, 3.0 / T_FAKTOR );

    real cs  = sqrt( ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * prim.lambda ) );

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * mu ) / ( U * dx * rho ) ) ) );

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";

    *logging::out << logging::Logger::INFO_HIGH << "F_rho = " << U * rho * dt * 1000 << " kg/m^3\n";

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

    parameters.rhoRef    = rho;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    parameters.enableReaction = true;

    *logging::out << logging::Logger::INFO_HIGH << "Pr = " << parameters.Pr << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool threeDimensional = false;

    if( threeDimensional )
    {
        gridBuilder->addCoarseGrid(-0.5*L, -0.5*L, 0.0,
                                    0.5*L, 0.5*L, H, dx);
    }
    else
    {
        gridBuilder->addCoarseGrid(-0.5*L, -0.5*dx, 0.0,
                                    0.5*L, 0.5*dx, H, dx);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
    //TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("F:/Work/Computations/inp/Ring.stl");
#else
    //TriangularMesh* stl = TriangularMesh::make("inp/Unterzug.stl");
    TriangularMesh* stl = TriangularMesh::make("inp/Ring.stl");
#endif

    //gridBuilder->addGeometry(stl);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    VerticalCylinder cylinder( 0.0, 0.0, 0.0, 1.1, 4.0   );
    
    Conglomerate refRing;

    refRing.add     ( new VerticalCylinder( 0.0, 0.0, 0.0, 0.6, 0.125 ) );
    refRing.subtract( new VerticalCylinder( 0.0, 0.0, 0.0, 0.4, 1.0    ) );
    //refRing.add     ( new VerticalCylinder( 0.0, 0.0, 0.0, 0.15, 0.125 ) );
    //refRing.subtract( new VerticalCylinder( 0.0, 0.0, 0.0, 0.05, 1.0    ) );

    gridBuilder->setNumberOfLayers(0,20);

    gridBuilder->addGrid( &cylinder, 1 );

    gridBuilder->setNumberOfLayers(10,20);

    //gridBuilder->addGrid( &refRing, 2 );
    //gridBuilder->addGrid( stl, 2 );

    if( threeDimensional ) gridBuilder->setPeriodicBoundaryCondition(false, false, false);
    else                   gridBuilder->setPeriodicBoundaryCondition(false, true,  false);

    gridBuilder->buildGrids(GKS, false);

    //gridBuilder->writeGridsToVtk(path + "grid/Grid_lev_");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    GksMeshAdapter meshAdapter( gridBuilder );

    meshAdapter.inputGrid();

    //meshAdapter.writeMeshVTK( path + "grid/Mesh.vtk" );

    //meshAdapter.writeMeshFaceVTK( path + "grid/MeshFaces.vtk" );

    if( !threeDimensional )
        meshAdapter.findPeriodicBoundaryNeighbors();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    CudaUtility::setCudaDevice(0);

    auto dataBase = std::make_shared<DataBase>( "GPU" );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    real openBoundaryVelocityLimiter = 1.0;
    
    SPtr<BoundaryCondition> bcMX = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<MassCompensation>( dataBase, rho, U, prim.lambda );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<MassCompensation>( dataBase, rho, U, prim.lambda );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, false );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, false );
    //SPtr<BoundaryCondition> bcMX = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );
    //SPtr<BoundaryCondition> bcPX = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );

    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    //SPtr<BoundaryCondition> bcMX_2 = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, false );
    //SPtr<BoundaryCondition> bcPX_2 = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, false );
    SPtr<BoundaryCondition> bcMX_2 = std::make_shared<Symmetry>( dataBase, 'x' );
    SPtr<BoundaryCondition> bcPX_2 = std::make_shared<Symmetry>( dataBase, 'x' );
    //SPtr<BoundaryCondition> bcMX_2 = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );
    //SPtr<BoundaryCondition> bcPX_2 = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );

    bcMX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L && center.z > H - 0.5; } );
    bcPX_2->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L && center.z > H - 0.5; } );

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMY;
    SPtr<BoundaryCondition> bcPY;

    if( threeDimensional )
    {
        //bcMY = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
        //bcPY = std::make_shared<Open>( dataBase, prim, openBoundaryVelocityLimiter );
        bcMY = std::make_shared<Symmetry>( dataBase, 'y' );
        bcPY = std::make_shared<Symmetry>( dataBase, 'y' );

        bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L; } );
        bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L; } );
    }
    else
    {
        bcMY = std::make_shared<Periodic>(dataBase);
        bcPY = std::make_shared<Periodic>(dataBase);

        bcMY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y < -0.5*dx; });
        bcPY->findBoundaryCells(meshAdapter, false, [&](Vec3 center) { return center.y >  0.5*dx; });
    }

    //////////////////////////////////////////////////////////////////////////
    
    SPtr<BoundaryCondition> bcMZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<IsothermalWall>( dataBase, Vec3(0, 0, 0), prim.lambda, true );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<InflowComplete>( dataBase, PrimitiveVariables(rho, 0.0, 0.0, 0.0, prim.lambda, 0.0, 0.0) );
    //SPtr<BoundaryCondition> bcMZ = std::make_shared<Open>( dataBase );

    //SPtr<BoundaryCondition> bcPZ = std::make_shared<Open>( dataBase, prim );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<Extrapolation>( dataBase );
    //SPtr<BoundaryCondition> bcPZ = std::make_shared<AdiabaticWall>( dataBase, Vec3(0, 0, 0), true );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Pressure2>( dataBase, c1o2 * prim.rho / prim.lambda );
    
    bcMZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z < 0.0; } );
    bcPZ->findBoundaryCells( meshAdapter, true, [&](Vec3 center){ return center.z > H  ; } );

    //////////////////////////////////////////////////////////////////////////

    //SPtr<BoundaryCondition> burner = std::make_shared<IsothermalWall>( dataBase, Vec3(0.0, 0.0, 0.0), 0.5*prim.lambda,  0.0, true );

    //SPtr<BoundaryCondition> burner = std::make_shared<InflowComplete>( dataBase, PrimitiveVariables(rho, 0.0, 0.0, U, prim.lambda, 1.0, 1.0) );
    SPtr<BoundaryCondition> burner = std::make_shared<CreepingMassFlux>( dataBase, rho, U, prim.lambda );

    burner->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ 

        if( threeDimensional )
            return center.z < 0.0 && std::sqrt(center.x*center.x + center.y*center.y) < 0.5;
        else
            return center.z < 0.0 && std::sqrt(center.x*center.x) < 0.5 && std::sqrt(center.y*center.y) < 0.5 * dx;
    } );

    //////////////////////////////////////////////////////////////////////////

    dataBase->boundaryConditions.push_back( bcMX );
    dataBase->boundaryConditions.push_back( bcPX );
    
    dataBase->boundaryConditions.push_back( bcMY );
    dataBase->boundaryConditions.push_back( bcPY );

    dataBase->boundaryConditions.push_back( bcMZ );
    dataBase->boundaryConditions.push_back( bcPZ );

    dataBase->boundaryConditions.push_back( bcMX_2 );
    dataBase->boundaryConditions.push_back( bcPX_2 );

    dataBase->boundaryConditions.push_back( burner );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint startIter = 0;

    dataBase->setMesh( meshAdapter );

    CudaUtility::printCudaMemoryUsage();
    
    if( restartIter == INVALID_INDEX )
    {
        Initializer::interpret(dataBase, [&](Vec3 cellCenter) -> ConservedVariables {

            PrimitiveVariables primLocal = prim;

            //primLocal.rho = rho * std::exp( - ( 2.0 * g * H * prim.lambda ) * cellCenter.z / H );

            real r = sqrt(cellCenter.x * cellCenter.x /*+ cellCenter.y * cellCenter.y*/ + cellCenter.z * cellCenter.z);

            //if( r < 0.6 ) primLocal.S_1 = 1.0 - r;

            //if( r < 0.5 ) prim.lambda /= (two - four*r*r);

            return toConservedVariables(primLocal, parameters.K);
        });

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );
    }
    else
    {
        Restart::readRestart( dataBase, path + simulationName + "_" + std::to_string( restartIter ), startIter );

        writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( restartIter ) + "_restart" );
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

    CupsAnalyzer cupsAnalyzer( dataBase, true, 30.0 );

    ConvergenceAnalyzer convergenceAnalyzer( dataBase, 1000 );

    //auto turbulenceAnalyzer = std::make_shared<TurbulenceAnalyzer>( dataBase, 50000 );

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    for( uint iter = startIter + 1; iter <= 2000000; iter++ )
    {
        uint runUpTime = 10000;

        if( iter < runUpTime )
        {
            //std::dynamic_pointer_cast<InflowComplete>(burner)->prim.S_1 =       1.0 * ( real(iter) / 20000.0 );
            //std::dynamic_pointer_cast<InflowComplete>(burner)->prim.S_2 = 1.0 - 1.0 * ( real(iter) / 20000.0 );

            //std::dynamic_pointer_cast<InflowComplete>(burner)->prim.W = U * ( real(iter) / 20000.0 );

            //std::dynamic_pointer_cast<CreepingMassFlux>(burner)->velocity = U * ( real(iter) / runUpTime );

            //parameters.mu = mu + 10.0 * mu * ( 1.0 - ( real(iter) / 20000.0 ) );

            //parameters.dt = 0.2 * dt + ( dt - 0.2 * dt ) * ( real(iter) / 40000.0 );
        }

        //if( iter == 5001 )
        //{
        //    parameters.enableReaction = false;
        //    std::dynamic_pointer_cast<CreepingMassFlux>(burner)->velocity = -1.0;
        //}

        cupsAnalyzer.run( iter );

        convergenceAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( 
            //( iter >= 20000 && iter % 1 == 0 ) || 
            ( iter % 400 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        if( iter % 10000 == 0 )
        {
            Restart::writeRestart( dataBase, path + simulationName + "_" + std::to_string( iter ), iter );
        }

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
    std::string path( "F:/Work/Computations/out/PoolFire/" );
#else
    std::string path( "out/" );
#endif

    std::string simulationName ( "PoolFire" );

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);

    if( sizeof(real) == 4 )
        *logging::out << logging::Logger::INFO_HIGH << "Using Single Precision\n";
    else
        *logging::out << logging::Logger::INFO_HIGH << "Using Double Precision\n";

    try
    {
        uint restartIter = INVALID_INDEX;
        //uint restartIter = 420000;

        if( argc > 1 ) restartIter = atoi( argv[1] );

        thermalCavity( path, simulationName, restartIter );
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
