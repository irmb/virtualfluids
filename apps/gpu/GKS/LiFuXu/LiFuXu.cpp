//#define MPI_LOGGING

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <iostream>
#include <exception>
#include <fstream>
#include <memory>

#include "Core/Timer/Timer.h"
#include "PointerDefinitions.h"
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
#include "GksGpu/Definitions/MemoryAccessPattern.h"

real solution(Vec3 point, const double U, const double V, const double D, const double time)
{
	return c1o4 * ( erf( (  0.225 - ( point.x - U * time ) ) / ( two * sqrt( D * time ) ) ) 
                  + erf( ( -0.175 + ( point.x - U * time ) ) / ( two * sqrt( D * time ) ) )
                  )
				* ( erf( (  0.225 - ( point.y - V * time ) ) / ( two * sqrt( D * time ) ) ) 
                  + erf( ( -0.175 + ( point.y - V * time ) ) / ( two * sqrt( D * time ) ) )
                  );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Vec3 cellCenter(std::shared_ptr<DataBase> dataBase, uint cellIdx)
{
	Vec3 cellCenter;

	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][0]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][0]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][1]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][1]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][2]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][2]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][3]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][3]].y;

	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][4]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][4]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][5]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][5]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][6]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][6]].y;
	cellCenter.x += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][7]].x;
	cellCenter.y += c1o8 * dataBase->nodeCoordinates[dataBase->cellToNode[cellIdx][7]].y;

	return cellCenter;
}

void printL_2Norm(const std::shared_ptr<DataBase> dataBase, const real U, const real V, const real D, const real time)
{
	dataBase->copyDataDeviceToHost();

	double l_2 = zero;
	double sum = zero;

	for (uint cellIdx = 0; cellIdx < dataBase->perLevelCount[0].numberOfBulkCells; cellIdx++)
	{
		Vec3 center = cellCenter(dataBase, cellIdx);

		//double simulatedResult = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
		double simulatedResult = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
		
        double analyticResult  = solution(center, U, V, D, time);

		double err = abs(simulatedResult - analyticResult);

		sum += analyticResult * analyticResult;

		l_2 += err * err;

		//std::cout << std::endl << err << " " << analyticResult;
	}
	l_2 = sqrt( l_2/sum );
	std::cout << std::endl << "The l2 norm is " << l_2 << std::endl;
}

void printL_MaxNorm(const std::shared_ptr<DataBase> dataBase, const real U, const real V, const real D, const real time)
{
	dataBase->copyDataDeviceToHost();

	double max = zero;

	for (uint cellIdx = 0; cellIdx < dataBase->perLevelCount[0].numberOfBulkCells; cellIdx++)
	{
		Vec3 center = cellCenter(dataBase, cellIdx);

		//double simulatedResult = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
		double simulatedResult = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
		
        double analyticResult  = solution(center, U, V, D, time);

		double err = abs(simulatedResult - analyticResult);

		if( err > max ) max = err;

		//std::cout << std::endl << err << " " << analyticResult;
	}
	std::cout << std::endl << "The max norm is " << max << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void thermalCavity( std::string path, std::string simulationName )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    uint nx = 128;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real L = 1.0;

    real dx = L / real(nx);

    real U = 100.0;

    real Ma = 0.1;

    real Pr  = 1.0;
    real K   = 2.0;
    
    real rho = 1.0;
    
    real mu = 0.01;

    real D = 1.5;

    real cs = U / Ma;
    PrimitiveVariables prim( rho, U, U, 0.0, ( ( K + 5.0 ) / ( K + 3.0 ) ) / ( 2.0 * cs * cs ) );

    real CFL = 0.25;

    real dt  = CFL * ( dx / ( ( U + cs ) * ( one + ( two * D ) / ( U * dx * rho ) ) ) );

    dt = 1.0e-6;

    *logging::out << logging::Logger::INFO_HIGH << "dt = " << dt << " s\n";
    *logging::out << logging::Logger::INFO_HIGH << "U  = " << U  << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "cs = " << cs << " m/s\n";
    *logging::out << logging::Logger::INFO_HIGH << "mu = " << mu << " kg/sm\n";

    //////////////////////////////////////////////////////////////////////////

    Parameters parameters;

    parameters.K  = K;
    parameters.Pr = Pr;
    parameters.mu = mu;

    parameters.D  = D;

    parameters.force.x = 0;
    parameters.force.y = 0;
    parameters.force.z = 0;

    parameters.dt = dt;
    parameters.dx = dx;

    //parameters.viscosityModel = ViscosityModel::sutherlandsLaw;
    parameters.viscosityModel = ViscosityModel::constant;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->addCoarseGrid( 0.0, 0.0, -0.5*dx,  
                                  L,   L,  0.5*dx, dx);

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

    SPtr<BoundaryCondition> bcMX = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPX = std::make_shared<Periodic>( dataBase );

    bcMX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x < -0.5*L; } );
    bcPX->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.x >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMY = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPY = std::make_shared<Periodic>( dataBase );

    bcMY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < -0.5*L; } );
    bcPY->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y >  0.5*L; } );

    //////////////////////////////////////////////////////////////////////////

    SPtr<BoundaryCondition> bcMZ = std::make_shared<Periodic>( dataBase );
    SPtr<BoundaryCondition> bcPZ = std::make_shared<Periodic>( dataBase );

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

        PrimitiveVariables localPrim = prim;

        //prim.S_1 = solution(cellCenter, U, U, D, 2e-3);
        prim.S_2 = solution(cellCenter, U, U, D, 2e-3);

        return toConservedVariables(localPrim, parameters.K);
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

    //////////////////////////////////////////////////////////////////////////

    cupsAnalyzer.start();

    uint maxIter = 1000;

    for( uint iter = 1; iter <= maxIter; iter++ )
    {
        cupsAnalyzer.run( iter );

        TimeStepping::nestedTimeStep(dataBase, parameters, 0);

        if( 
            ( iter % 100 == 0 )
          )
        {
            dataBase->copyDataDeviceToHost();

            writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
        }

        convergenceAnalyzer.run( iter );
    }

    //////////////////////////////////////////////////////////////////////////

    dataBase->copyDataDeviceToHost();

    std::cout << 2e-3 + maxIter * dt << std::endl;

    printL_2Norm  (dataBase, U, U, D, 2e-3 + maxIter * dt);
    printL_MaxNorm(dataBase, U, U, D, 2e-3 + maxIter * dt);
}

int main( int argc, char* argv[])
{
    std::string path( "F:/Work/Computations/out/LiFuXu/" );
    //std::string path( "out/" );
    std::string simulationName ( "LiFuXu" );

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
