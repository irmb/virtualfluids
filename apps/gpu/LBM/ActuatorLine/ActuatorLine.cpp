
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Core/Input/ConfigFileReader/ConfigFileReader.h"

#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

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
#include "VirtualFluids_GPU/Visitor/ActuatorLine.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

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

const real D = 126.0; // diameter in m

const real L_x = 10*D;
const real L_y = 6*D;
const real L_z = 6*D;

const real viscosity = 1.56e-5;

const real velocity  = 9.0;

const real mach = 0.1;

const uint nodes_per_D = 32;

//std::string path("F:/Work/Computations/out/DrivenCavity/"); //LEGOLAS
//std::string path("D:/out/DrivenCavity"); //Mollok
std::string path(".");

std::string simulationName("ActuatorLine");

const uint timeStepOut = 500;
const uint timeStepEnd = 2500;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    auto gridFactory = GridFactory::make();
    gridFactory->setGridStrategy(Device::CPU);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
	vf::gpu::Communicator* comm = vf::gpu::Communicator::getInstanz();
	SPtr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance();

    std::cout << configPath << std::endl;
	SPtr<ConfigData> configData = configReader->readConfigFile(configPath.c_str());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real dx = D/real(nodes_per_D);

	gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
							   L_x,  L_y,  L_z, dx);

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

        SPtr<Parameter>    para         = Parameter::make(configData, comm);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        const real dt = dx * mach / (sqrt(3) * velocity);

        const real velocityLB = velocity * dt / dx; // LB units

        const real viscosityLB = viscosity * dt / (dx * dx); // LB units

        *logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
        *logging::out << logging::Logger::INFO_HIGH << "viscosity [10^8 dx^2/dt] = " << viscosityLB*1e8 << "\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		para->setDevices(std::vector<uint>{(uint)0});

        para->setOutputPath( path );
        para->setOutputPrefix( simulationName );

        para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

        para->setPrintFiles(true);

        para->setMaxLevel(1);

        para->setVelocity(velocityLB);
        para->setViscosity(viscosityLB);

        para->setVelocityRatio( dx / dt );

		para->setMainKernel("CumulantK17CompChim");

		para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
            rho = (real)0.0;
            vx  = velocityLB;
            vy  = (real)0.0;
            vz  = (real)0.0;
        });

        para->setTOut( timeStepOut );
        para->setTEnd( timeStepEnd );

        para->setIsBodyForce( true );


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        gridBuilder->setVelocityBoundaryCondition(SideType::MX,  velocityLB,  0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::PX,  velocityLB,  0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MY,  velocityLB,  0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::PY,  velocityLB,  0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MZ,  velocityLB,  0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  velocityLB,  0.0, 0.0);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

        SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

        real turbPos[3] = {3*D, 3*D, 3*D};
        real epsilon = 5.f;
        real density = 1.225f;
        int level = 0;

        ActuatorLine* actuator_line = new ActuatorLine((unsigned int) 3, density, (unsigned int)32, epsilon, turbPos[0], turbPos[1], turbPos[2], D, level, dt, dx);
        para->addActuator( actuator_line );

        Simulation sim;
        SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
        SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
        SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
        sim.setFactories(kernelFactory, preProcessorFactory);
        sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);        
        sim.run();
        sim.free();

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
    MPI_Init(&argc, &argv);
    std::string str, str2; 
    if ( argv != NULL )
    {
        //str = static_cast<std::string>(argv[0]);
        
        try
        {
            //////////////////////////////////////////////////////////////////////////

			// std::string targetPath;

			// targetPath = __FILE__;

			// targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);



			// std::cout << targetPath << std::endl;

			multipleLevel(path + "/configActuatorLine.txt");

            //////////////////////////////////////////////////////////////////////////
		}
        catch (const std::bad_alloc& e)
        { 
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
        }
        catch (const std::exception& e)
        {   
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
        }
        catch (...)
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
        }
    }

    MPI_Finalize();
    return 0;
}
