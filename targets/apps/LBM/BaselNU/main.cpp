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


void multipleLevel(const std::string& configPath)
{
	std::ofstream logFile("F:/Basel2019NU/grid/gridGeneratorLog.txt");
	logging::Logger::addStream(&logFile);

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
	SPtr<GridProvider> gridGenerator;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool useGridGenerator = false;

	if (useGridGenerator) {

		real dx = 4.0;
		real vx = 0.05;

		TriangularMesh* BaselSTL = TriangularMesh::make("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl");

		gridBuilder->addCoarseGrid(-256.0, -256.0, -8.0,
			                        256.0, 256.0, 160.0, dx);

		gridBuilder->addGeometry(BaselSTL);

		//////////////////////////////////////////////////////////////////////////

		Cuboid* refBoxMX = new Cuboid( -300, -300, - 20, 
			                           -254,  300,  200 );
		Cuboid* refBoxPX = new Cuboid(  254, -300, - 20, 
			                            300,  300,  200 );
		Cuboid* refBoxMY = new Cuboid( -300, -300, - 20, 
			                            300, -254,  200 );
		Cuboid* refBoxPY = new Cuboid( -300,  254, - 20, 
			                            300,  300,  200 );

		Conglomerate* refRegion = new Conglomerate();
		
		refRegion->add(refBoxMX);
		refRegion->add(refBoxPX);
		refRegion->add(refBoxMY);
		refRegion->add(refBoxPY);

		refRegion->add(BaselSTL);

		gridBuilder->setNumberOfLayers(12, 8);

		gridBuilder->addGrid(refRegion, 2);

		//////////////////////////////////////////////////////////////////////////

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
		gridBuilder->writeGridsToVtk("F:/Basel2019NU/grids/BaselUniMergedXAllStreets/Basel_Grid");
		SimulationFileWriter::write("F:/Basel2019NU/grids/BaselUniMergedXAllStreets/", gridBuilder, FILEFORMAT::BINARY);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		StreetPointFinder finder;
		uint streetLevel = 2;

		finder.readStreets("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Streets.txt");

		finder.writeVTK("F:/Basel2019NU/results/ExampleStreets.vtk");

		finder.findIndicesLB(gridBuilder->getGrid(streetLevel), 7.0);

		//Merged for Wind in X Direction
		finder.writeConnectionVTK("F:/Basel2019NU/grids/BaselUniMergedXAllStreets/Basel_Grid/ExampleStreetsConnection.vtk", gridBuilder->getGrid(streetLevel));
		finder.writeSimulationFile("F:/Basel2019NU/grids/BaselUniMergedXAllStreets/", 1.0, gridBuilder->getNumberOfLevels(), streetLevel);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		return;

		gridGenerator = GridGenerator::make(gridBuilder, para);

	}
	else
	{
		gridGenerator = GridReader::make(FileFormat::BINARY, para);
		//gridGenerator = GridReader::make(FileFormat::ASCII, para);
	}

	logFile.close();

	//return;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//std::ifstream stream;
	//stream.open(configPath.c_str(), std::ios::in);
	//if (stream.fail())
	//    throw std::runtime_error("can not open config file!");

	//UPtr<input::Input> input = input::Input::makeInput(stream, "config");

	Simulation sim;
	SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
	sim.init(para, gridGenerator, fileWriter);
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
				*logging::out << logging::Logger::ERROR << e.what() << "\n";
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
				multipleLevel("F:/Basel2019NU/configBasel.txt"); //Baumbart 2
				//multipleLevel("F:/Work/Computations/gridGenerator/inp/configTest.txt");
				//multipleLevel("C:/Users/hiwi/Desktop/configBasel.txt"); //Gamling
			}
			catch (const std::exception& e)
			{

				*logging::out << logging::Logger::ERROR << e.what() << "\n";
				//std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (const std::bad_alloc e)
			{

				*logging::out << logging::Logger::ERROR << "Bad Alloc:" << e.what() << "\n";
				//std::cout << e.what() << std::flush;
				//MPI_Abort(MPI_COMM_WORLD, -1);
			}
			catch (...)
			{
				*logging::out << logging::Logger::ERROR << "Unknown exception!\n";
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
