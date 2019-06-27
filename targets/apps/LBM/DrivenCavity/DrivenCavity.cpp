
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <exception>
#include <memory>

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Core/Input/ConfigFileReader/ConfigFileReader.h"

#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

//#include "GridGenerator/global.h"

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

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "GksVtkAdapter/VTKInterface.h"

#include "GksGpu/DataBase/DataBase.h"
#include "GksGpu/Parameters/Parameters.h"
#include "GksGpu/Initializer/Initializer.h"

#include "GksGpu/FlowStateData/FlowStateDataConversion.cuh"

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"
#include "GksGpu/BoundaryConditions/IsothermalWall.h"

#include "GksGpu/TimeStepping/NestedTimeStep.h"

#include "GksGpu/Analyzer/CupsAnalyzer.h"
#include "GksGpu/Analyzer/ConvergenceAnalyzer.h"

#include "GksGpu/CudaUtility/CudaUtility.h"

//////////////////////////////////////////////////////////////////////////

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
    
	Communicator* comm = Communicator::getInstanz();
	SPtr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance();
	SPtr<ConfigData> configData = configReader->readConfigFile(configPath);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LbmOrGks lbmOrGks = LBM;

	const real L = 1.0;

    const real Re = 1000;

    const uint nx = 64;

	const real vx = 0.05; // LB units
	const real vy = 0.05; // LB units

	const real velocity = sqrt(vx*vx + vy*vy);

    const real viscosity = nx * velocity / Re; // LB units

    *logging::out << logging::Logger::INFO_HIGH << "velocity  = " << velocity << " s\n";

    *logging::out << logging::Logger::INFO_HIGH << "viscosity = " << viscosity << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	real dx = L / real(nx);

	gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L,
								0.5 * L,  0.5 * L,  0.5 * L, dx);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(LBM, false); // buildGrids() has to be called before setting the BCs!!!!

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
	    gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  vx,  vy, 0.0);
	    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SPtr<Parameter> para = Parameter::make(configData, comm);
        //SPtr<Parameter> para = Parameter::make();

        SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

        SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        para->setVelocity(velocity);

        para->setViscosity(viscosity);

        para->setVelocityRatio(1.0 / velocity);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        CudaUtility::setCudaDevice(0);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        GksMeshAdapter meshAdapter( gridBuilder );

        meshAdapter.inputGrid();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        auto dataBase = std::make_shared<DataBase>( "GPU" );

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        real lambda = 0.1;

        SPtr<BoundaryCondition> bcLid  = std::make_shared<IsothermalWall>( dataBase, Vec3(  vx, 0.0, 0.0 ), lambda, false );
        SPtr<BoundaryCondition> bcWall = std::make_shared<IsothermalWall>( dataBase, Vec3( 0.0, 0.0, 0.0 ), lambda, false );

        bcLid->findBoundaryCells ( meshAdapter, true,  [&](Vec3 center){ return center.y > 0.5; } );
        bcWall->findBoundaryCells( meshAdapter, false, [&](Vec3 center){ return center.y < 0.5; } );

        dataBase->boundaryConditions.push_back( bcLid  );
        dataBase->boundaryConditions.push_back( bcWall );
    
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

			

			//std::stringstream targetPath; targetPath << __FILE__;
			//
			//std::cout << targetPath.str() << std::endl;

			std::string targetPath;

			targetPath = __FILE__;

#ifdef _WIN32
			targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
			targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

			std::cout << targetPath << std::endl;

			multipleLevel(targetPath + "configDrivenCavity.txt");

            //////////////////////////////////////////////////////////////////////////
		}
        catch (const std::exception& e)
        {
                
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
        }
        catch (const std::bad_alloc e)
        {
                
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
        }
        catch (...)
        {
            *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
        }
    }

   MPI_Finalize();
   return 0;
}
