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

#include "Input/Input.h"
#include "StringUtilities/StringUtil.h"
#include "Input/ConfigFileReader/ConfigFileReader.h"

#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Communication/MpiCommunicator.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/Factories/BoundaryConditionFactory.h"

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
	BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
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
		auto BaselSTL = std::make_shared<TriangularMesh>("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl");
#else
		//Phoenix
		auto BaselSTL = std::make_shared<TriangularMesh>(gridpath + "/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl");
#endif


		gridBuilder->addCoarseGrid(-256.0, -256.0, -8.0,
			256.0, 256.0, 160.0, dx);

		gridBuilder->addGeometry(BaselSTL);

		//Merged for Wind in X Direction
		gridBuilder->setPeriodicBoundaryCondition(true, true, false);

		gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

		//////////////////////////////////////////////////////////////////////////

		gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx, 0.0, 0.0);
		gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx, 0.0, 0.0);

		gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

		//no forcing
		// gridBuilder->setPressureBoundaryCondition(SideType::PY, 0.0);
		// gridBuilder->setPressureBoundaryCondition(SideType::MY, 0.0);

		// gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
		// gridBuilder->setPressureBoundaryCondition(SideType::MX, 0.0);

		bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
		bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
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

		gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemManager, communicator);
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


	Simulation sim;
	SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
	SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
	SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
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
				std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (const std::bad_alloc e)
			{
				std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (...)
			{
				std::cout << "unknown exeption" << std::endl;
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