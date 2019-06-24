//#define MPI_LOGGING

//Martin Branch

#include <mpi.h>
#if defined( MPI_LOGGING )
#include <mpe.h>
#endif

#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

//#include "metis.h"

#include "Core/LbmOrGks.h"
#include "Core/Input/Input.h"
#include "Core/StringUtilities/StringUtil.h"
#include "Core/Input/ConfigFileReader/ConfigFileReader.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

#include "global.h"

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

#include "geometries/Sphere/Sphere.h"
#include "geometries/VerticalCylinder/VerticalCylinder.h"
#include "geometries/Cuboid/Cuboid.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/Conglomerate/Conglomerate.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridFactory.h"

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"

#include "utilities/math/Math.h"
#include "utilities/communication.h"
#include "utilities/transformator/TransformatorImp.h"

#include "Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"


void multipleLevel(const std::string& configPath)
{
	//std::ofstream logFile( "F:/Work/Computations/gridGenerator/grid/gridGeneratorLog.txt" );
	//std::ofstream logFile("F:/Basel2019/log/gridGeneratorLog.txt");
	//logging::Logger::addStream(&logFile);

	logging::Logger::addStream(&std::cout);
	logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
	logging::Logger::timeStamp(logging::Logger::ENABLE);
	logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

	//UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

	auto gridFactory = GridFactory::make();
	gridFactory->setGridStrategy(Device::CPU);
	//gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);
	gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
	//gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE);

	auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

	SPtr<ConfigFileReader> configReader = ConfigFileReader::getNewInstance();
	SPtr<ConfigData> configData = configReader->readConfigFile(configPath);
	Communicator* comm = Communicator::getInstanz();

	SPtr<Parameter> para = Parameter::make(configData, comm);
	SPtr<CudaMemoryManager> cudaMemManager = CudaMemoryManager::make(para);
	SPtr<GridProvider> gridGenerator;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _WIN32
	//Baumbart
	std::string gridpath = "F:/Basel2019";
#else
	//Phoenix
	std::string gridpath = "/work/marschoe/Basel4GPU";
#endif // _WIN32

	
	std::ofstream logFile;
	logFile.open(gridpath + "/gridGeneratorLog.txt");

	bool useGridGenerator = false;

	if (useGridGenerator) {

		real dx = 1.0;
		real vx = 0.05;

#ifdef _WIN32
		//Baumbart
		TriangularMesh* BaselSTL = TriangularMesh::make("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl");
#else
		//Phoenix
		TriangularMesh* BaselSTL = TriangularMesh::make(gridpath + "/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl");
#endif


		gridBuilder->addCoarseGrid(-256.0, -256.0, -8.0,
			256.0, 256.0, 160.0, dx);

		gridBuilder->addGeometry(BaselSTL);

		//Merged for Wind in X Direction
		gridBuilder->setPeriodicBoundaryCondition(true, true, false);

		gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

		//////////////////////////////////////////////////////////////////////////

		gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx, 0.0, 0.0);
		gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx, 0.0, 0.0);

		gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

		//no forcing
		gridBuilder->setPressureBoundaryCondition(SideType::PY, 0.0);
		gridBuilder->setPressureBoundaryCondition(SideType::MY, 0.0);

		gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
		gridBuilder->setPressureBoundaryCondition(SideType::MX, 0.0);

		//////////////////////////////////////////////////////////////////////////
		//Merged for Wind in X Direction
		//gridBuilder->writeGridsToVtk(gridpath + "/grids/BaselUni/Basel_Grid");
		//SimulationFileWriter::write(gridpath + "/grids/BaselUni/", gridBuilder, FILEFORMAT::BINARY);
		////////////////////
		//one street closed
		gridBuilder->writeGridsToVtk(gridpath + "/grids/BaselUniOSC_1m/Basel_Grid");
		SimulationFileWriter::write(gridpath + "/grids/BaselUniOSC_1m/", gridBuilder, FILEFORMAT::BINARY);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		StreetPointFinder finder;
#ifdef _WIN32
		//Baumbart
		finder.readStreets("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Streets.txt");
#else
		//Phoenix
		finder.readStreets(gridpath + "/source/git/targets/apps/LBM/Basel/resources/Streets.txt");
#endif

		finder.writeVTK(gridpath + "/results/ExampleStreets.vtk");

		finder.findIndicesLB(gridBuilder->getGrid(0), 7.0);

		//all Streets
		//finder.writeConnectionVTK(gridpath + "/grids/BaselUni/Basel_Grid/ExampleStreetsConnection.vtk", gridBuilder->getGrid(0));
		//finder.writeSimulationFile(gridpath + "/grids/BaselUni/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		//finder.writeStreetVectorFile(gridpath + "/grids/BaselUni/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		////////////////////
		//one street closed
		finder.writeConnectionVTK(gridpath + "/grids/BaselUniOSC_1m/Basel_Grid/ExampleStreetsConnection.vtk", gridBuilder->getGrid(0));
		finder.writeSimulationFile(gridpath + "/grids/BaselUniOSC_1m/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		finder.writeStreetVectorFile(gridpath + "/grids/BaselUniOSC_1m/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		return;

		gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemManager);
		//gridGenerator = GridGenerator::make(gridBuilder, para);

	}
	else
	{
		gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemManager);
		//gridGenerator = GridReader::make(FileFormat::BINARY, para);
		//gridGenerator = GridReader::make(FileFormat::ASCII, para);
	}

	logFile.close();

	//return;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::shared_ptr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
	std::shared_ptr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();

	Simulation sim;
	SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
	sim.setFactories(kernelFactory, preProcessorFactory);
	sim.init(para, gridGenerator, fileWriter, cudaMemManager);
	sim.run();
	sim.free();
}


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	std::string str, str2;
	if (argv != NULL)
	{
		str = static_cast<std::string>(argv[0]);
		if (argc > 1)
		{
			str2 = static_cast<std::string>(argv[1]);
			try
			{
				multipleLevel(str2);
			}
			catch (const std::exception& e)
			{
				*logging::out << logging::Logger::LOGGED_ERROR << e.what() << "\n";
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (...)
			{
				std::cout << "unknown exeption" << std::endl;
			}
		}
		else
		{
			try
			{
				//multipleLevel("E:/temp/Basel2019/config/configBasel.txt"); //Tesla03
				//multipleLevel("C:/Users/schoen/Desktop/bin/ReleaseBasel/configBasel.txt"); //Baumbart 1
				multipleLevel("F:/Basel2019/configBasel.txt"); //Baumbart 2
				//multipleLevel("F:/Work/Computations/gridGenerator/inp/configTest.txt");
				//multipleLevel("C:/Users/hiwi/Desktop/configBasel.txt"); //Gamling
			}
			catch (const std::exception& e)
			{

				*logging::out << logging::Logger::LOGGED_ERROR << e.what() << "\n";
				//std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (const std::bad_alloc e)
			{

				*logging::out << logging::Logger::LOGGED_ERROR << "Bad Alloc:" << e.what() << "\n";
				//std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (...)
			{
				*logging::out << logging::Logger::LOGGED_ERROR << "Unknown exception!\n";
				//std::cout << "unknown exeption" << std::endl;
			}

			std::cout << "\nConfiguration file must be set!: lbmgm <config file>" << std::endl << std::flush;
			//MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}


	/*
	MPE_Init_log() & MPE_Finish_log() are NOT needed when
	liblmpe.a is linked with this program.  In that case,
	MPI_Init() would have called MPE_Init_log() already.
	*/
#if defined( MPI_LOGGING )
	MPE_Init_log();
#endif

#if defined( MPI_LOGGING )
	if (argv != NULL)
		MPE_Finish_log(argv[0]);
	if (str != "")
		MPE_Finish_log(str.c_str());
	else
		MPE_Finish_log("TestLog");
#endif

	MPI_Finalize();
	return 0;
}
