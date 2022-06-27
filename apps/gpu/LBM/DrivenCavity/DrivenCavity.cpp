
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>
#include <filesystem>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////

//#include "GksMeshAdapter/GksMeshAdapter.h"

//#include "GksVtkAdapter/VTKInterface.h"
//
//#include "GksGpu/DataBase/DataBase.h"
//#include "GksGpu/Parameters/Parameters.h"
//#include "GksGpu/Initializer/Initializer.h"
//
//#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"
//
//#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
//#include "GksGpu/BoundaryConditions/IsothermalWall.h"
//
//#include "GksGpu/TimeStepping/NestedTimeStep.h"
//
//#include "GksGpu/Analyzer/CupsAnalyzer.h"
//#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"
//
//#include "GksGpu/CudaUtility/CudaUtility.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//LbmOrGks lbmOrGks = GKS;
LbmOrGks lbmOrGks = LBM;

const real L  = 1.0;

const real Re = 500.0;// 1000.0;

const real velocity  = 1.0;

const real dt = (real)1.0e-3; //0.5e-3;

const uint nx = 64;

std::string path("output/");
std::string gridPath("grid/");

std::string simulationName("DrivenCavityChim");

const uint timeStepOut = 10000;
const uint timeStepEnd = 250000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    auto gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	real dx = L / real(nx);

	//gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L,
	//							0.5 * L,  0.5 * L,  0.5 * L, dx);

	gridBuilder->addCoarseGrid(-2.0 * dx, -0.5 * L, -0.5 * L,
								2.0 * dx,  0.5 * L,  0.5 * L, dx);

    auto refBox = new Cuboid(-0.1 * L, -0.1 * L, -0.1 * L,
                              0.1 * L,  0.1 * L,  0.1 * L);

    gridBuilder->addGrid(refBox, 1);

    gridBuilder->setNumberOfLayers(0, 0);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( lbmOrGks == LBM )
    {

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        vf::basics::ConfigurationFile config;
        config.load(configPath);

        SPtr<Parameter> para = std::make_shared<Parameter>(config, communicator.getNummberOfProcess(), communicator.getPID());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        const real velocityLB = velocity * dt / dx; // LB units

	    const real vx = velocityLB / (real)sqrt(2.0); // LB units
	    const real vy = velocityLB / (real)sqrt(2.0); // LB units

        const real viscosityLB = nx * velocityLB / Re; // LB units

        VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
        VF_LOG_INFO("viscosity [dx^2/dt] = {}", viscosityLB);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		para->setDevices(std::vector<uint>{(uint)0});

        para->setOutputPath( path ); // optional, default is output/
        para ->setGridPath( gridPath );  // optional, default is grid/

        para->setOutputPrefix( simulationName );

        para->setPrintFiles(true);

        para->setMaxLevel(2);

        para->setVelocity(velocityLB);
        para->setViscosity(viscosityLB);

        para->setVelocityRatio(velocity/ velocityLB);

		//para->setMainKernel("CumulantK17CompChim");

		para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
            rho = (real)0.0;
            vx  = (real)0.0; //(6 * velocityLB * coordZ * (L - coordZ) / (L * L));
            vy  = (real)0.0;
            vz  = (real)0.0;
        });

        para->setTOut( timeStepOut );
        para->setTEnd( timeStepEnd );

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
		//gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
		//gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
	    //gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  vx,  vx, 0.0);
	    //gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        gridBuilder->writeGridsToVtk(para->getGridPath());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

        auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator);
        sim.run();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
    else
    {
     //   CudaUtility::setCudaDevice(0);
     //
     //   Parameters parameters;

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	    //const real vx = velocity / sqrt(2.0);
	    //const real vy = velocity / sqrt(2.0);

     //   parameters.K  = 2.0;
     //   parameters.Pr = 1.0;
     //
     //   const real Ma = 0.1;

     //   real rho = 1.0;

     //   real cs = velocity / Ma;
     //   real lambda = c1o2 * ( ( parameters.K + 5.0 ) / ( parameters.K + 3.0 ) ) / ( cs * cs );

     //   const real mu = velocity * L * rho / Re;

     //   *logging::out << logging::Logger::INFO_HIGH << "mu  = " << mu << " m^2/s\n";

     //   *logging::out << logging::Logger::INFO_HIGH << "CFL = " << dt * ( velocity + cs ) / dx << "\n";

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   parameters.mu = mu;

     //   parameters.dt = dt;
     //   parameters.dx = dx;

     //   parameters.lambdaRef = lambda;

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   GksMeshAdapter meshAdapter( gridBuilder );

     //   meshAdapter.inputGrid();

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   auto dataBase = std::make_shared<DataBase>( "GPU" );

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   SPtr<BoundaryCondition> bcLid  = std::make_shared<IsothermalWall>( dataBase, Vec3(  vx,  vy, 0.0 ), lambda, false );
     //   SPtr<BoundaryCondition> bcWall = std::make_shared<IsothermalWall>( dataBase, Vec3( 0.0, 0.0, 0.0 ), lambda, false );

     //   bcLid->findBoundaryCells ( meshAdapter, true,  [&](Vec3 center){ return center.z > 0.5; } );
     //   bcWall->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.z < 0.5; } );

     //   dataBase->boundaryConditions.push_back( bcLid  );
     //   dataBase->boundaryConditions.push_back( bcWall );

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   dataBase->setMesh( meshAdapter );

     //   Initializer::interpret(dataBase, [&] ( Vec3 cellCenter ) -> ConservedVariables {

     //       return toConservedVariables( PrimitiveVariables( rho, 0.0, 0.0, 0.0, lambda ), parameters.K );
     //   });

     //   dataBase->copyDataHostToDevice();

     //   Initializer::initializeDataUpdate(dataBase);

     //   writeVtkXML( dataBase, parameters, 0, path + simulationName + "_0" );

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   CupsAnalyzer cupsAnalyzer( dataBase, false, 60.0, true, 10000 );

     //   ConvergenceAnalyzer convergenceAnalyzer( dataBase, 10000 );

     //   cupsAnalyzer.start();

     //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     //   for( uint iter = 1; iter <= timeStepEnd; iter++ )
     //   {
     //       TimeStepping::nestedTimeStep(dataBase, parameters, 0);

     //       if( iter % timeStepOut == 0 )
     //       {
     //           dataBase->copyDataDeviceToHost();

     //           writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );
     //       }
     //
     //       int crashCellIndex = dataBase->getCrashCellIndex();
     //       if( crashCellIndex >= 0 )
     //       {
     //           *logging::out << logging::Logger::LOGGER_ERROR << "Simulation Crashed at CellIndex = " << crashCellIndex << "\n";
     //           dataBase->copyDataDeviceToHost();
     //           writeVtkXML( dataBase, parameters, 0, path + simulationName + "_" + std::to_string( iter ) );

     //           break;
     //       }

     //       dataBase->getCrashCellIndex();

     //       cupsAnalyzer.run( iter, parameters.dt );

     //       convergenceAnalyzer.run( iter );
     //   }
    }
}

int main( int argc, char* argv[])
{
    try
    {
        vf::logging::Logger::initalizeLogger();

        // assuming that the config files is stored parallel to this file.
        std::filesystem::path filePath = __FILE__;
        filePath.replace_filename("configDrivenCavity.txt");

        multipleLevel(filePath.string());
    }
    catch (const spdlog::spdlog_ex &ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    }
    catch (const std::bad_alloc& e)
    {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    }
    catch (const std::exception& e)
    {
        VF_LOG_CRITICAL("exception: {}", e.what());
    }
    catch (...)
    {
        VF_LOG_CRITICAL("Unknown exception!");
    }

   return 0;
}
