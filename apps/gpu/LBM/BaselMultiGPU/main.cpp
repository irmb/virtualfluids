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
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"

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
    //std::ofstream logFile( "F:/Work/Computations/gridGenerator/grid/gridGeneratorLog.txt" );
    //std::ofstream logFile( "grid/gridGeneratorLog.txt" );
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
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
	SPtr<CudaMemoryManager> cudaMemManager = CudaMemoryManager::make(para);
	SPtr<GridProvider> gridGenerator;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	const uint generatePart = comm->getPID();

    std::ofstream logFile2;
	std::string gridpath = "/work/marschoe/Basel4GPU/";
	logFile2.open(gridpath + std::to_string(generatePart) + "/gridGeneratorLog.txt");//Phoenix
	//logFile2.open(std::string("M:/Basel2019/grids4/") + std::to_string(generatePart) + "/gridGeneratorLog.txt");//Baumbart

	logging::Logger::addStream(&logFile2);

    bool useGridGenerator = false;

    if(useGridGenerator){
        real dx = 1.0;
        real vx = 0.05;

        SPtr<TriangularMesh> BaselSTL;

		if (generatePart == 0)
			BaselSTL = std::make_shared<TriangularMesh>("/work/marschoe/Basel4GPU/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl"); //Phoenix
			//BaselSTL = std::make_shared<TriangularMesh>("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND.stl"); //Baumbart
		if (generatePart == 1)
			BaselSTL = std::make_shared<TriangularMesh>("/work/marschoe/Basel4GPU/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_X.stl"); //Phoenix
			//BaselSTL = std::make_shared<TriangularMesh>("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_X.stl"); //Baumbart
		if (generatePart == 2)
			BaselSTL = std::make_shared<TriangularMesh>("/work/marschoe/Basel4GPU/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_X_Y.stl"); //Phoenix
			//BaselSTL = std::make_shared<TriangularMesh>("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_X_Y.stl"); //Baumbart
		if (generatePart == 3)
			BaselSTL = std::make_shared<TriangularMesh>("/work/marschoe/Basel4GPU/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_Y.stl"); //Phoenix
			//BaselSTL = std::make_shared<TriangularMesh>("M:/Basel2019/stl/BaselUrbanProfile_066_deg_bridge_3_All_CLOSED_WIDE_GROUND_MIRROR_Y.stl"); //Baumbart

		real lengthInXDirection = 512.0;
		real lengthInYDirection = 512.0;
		real lengthInXDirectionOverlap = 520.0;
		real lengthInYDirectionOverlap = 520.0;
		real centerInX;
		real centerInY;

		if (generatePart == 0){	centerInX =    0.0;	centerInY =   0.0;}
		if (generatePart == 1){	centerInX = -512.0;	centerInY =   0.0;}
		if (generatePart == 2){	centerInX = -512.0;	centerInY = 512.0;}
		if (generatePart == 3){	centerInX =    0.0;	centerInY = 512.0;}

		gridBuilder->addCoarseGrid( (centerInX - lengthInXDirectionOverlap * 0.5), (centerInY - lengthInYDirectionOverlap * 0.5), -8.0,
									(centerInX + lengthInXDirectionOverlap * 0.5), (centerInY + lengthInYDirectionOverlap * 0.5),  160.0, dx);

        gridBuilder->addGeometry(BaselSTL);

		gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>( (centerInX - lengthInXDirection * 0.5), (centerInX + lengthInXDirection * 0.5),
																	(centerInY - lengthInYDirection * 0.5), (centerInY + lengthInYDirection * 0.5),
																	-200.0, 200.0));

		gridBuilder->setPeriodicBoundaryCondition(false, false, false);

        gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

		//////////////////////////////////////////////////////////////////////////

		uint xNeighborRank, yNeighborRank;
		if (generatePart == 0) { xNeighborRank = 1; yNeighborRank = 3; }
		if (generatePart == 1) { xNeighborRank = 0; yNeighborRank = 2; }
		if (generatePart == 2) { xNeighborRank = 3; yNeighborRank = 1; }
		if (generatePart == 3) { xNeighborRank = 2; yNeighborRank = 0; }

		gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
		gridBuilder->setCommunicationProcess (CommunicationDirections::PX, xNeighborRank);
		gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
		gridBuilder->setCommunicationProcess (CommunicationDirections::MX, xNeighborRank);
		gridBuilder->findCommunicationIndices(CommunicationDirections::PY);
		gridBuilder->setCommunicationProcess (CommunicationDirections::PY, yNeighborRank);
		gridBuilder->findCommunicationIndices(CommunicationDirections::MY);
		gridBuilder->setCommunicationProcess (CommunicationDirections::MY, yNeighborRank);

		//////////////////////////////////////////////////////////////////////////

        gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vx , 0.0, 0.0);
        gridBuilder->setVelocityBoundaryCondition(SideType::MZ, vx , 0.0, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);

		if (generatePart == 0 || generatePart == 3) gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
		if (generatePart == 1 || generatePart == 2)	gridBuilder->setPressureBoundaryCondition(SideType::MX, 0.0);
		if (generatePart == 2 || generatePart == 3) gridBuilder->setPressureBoundaryCondition(SideType::PY, 0.0);
		if (generatePart == 0 || generatePart == 1) gridBuilder->setPressureBoundaryCondition(SideType::MY, 0.0);

		//////////////////////////////////////////////////////////////////////////

		gridBuilder->writeGridsToVtk(gridpath + std::to_string(generatePart) + "/Basel_Grid");
		SimulationFileWriter::write(gridpath + std::to_string(generatePart )+ "/", gridBuilder, FILEFORMAT::BINARY);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (generatePart == 0) {
			StreetPointFinder finder;
			finder.readStreets("/work/marschoe/Basel4GPU/source/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");//phoenix
			//finder.readStreets("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");//Baumbart
			//finder.writeVTK("M:/Basel2019/results/ExampleStreets.vtk");
			finder.findIndicesLB(gridBuilder->getGrid(0), 7.0);
			//finder.writeConnectionVTK("M:/Basel2019/grids/BaselUniMergedX/Basel_Grid/ExampleStreetsConnection.vtk", gridBuilder->getGrid(0));
			finder.writeSimulationFile(gridpath + std::to_string(generatePart) + "/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		}
		else
		{
			StreetPointFinder finder;
			finder.writeSimulationFile(gridpath + std::to_string(generatePart) + "/", 1.0, gridBuilder->getNumberOfLevels(), 0);
		}
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

    logFile2.close();

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
    sim.init(para, gridGenerator, fileWriter, cudaMemManager);
    sim.run();
	sim.free();
}


int main( int argc, char* argv[])
{
     MPI_Init(&argc, &argv);
    std::string str, str2; 
    if ( argv != NULL )
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
                *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
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
				//multipleLevel("C:/Users/schoen/Desktop/bin/ReleaseBasel/configBasel.txt"); //Baumbart
				multipleLevel("/work/marschoe/Basel4GPU/source/configBasel.txt"); //Phoenix
				//multipleLevel("F:/Work/Computations/gridGenerator/inp/configTest.txt");
            }
            catch (const std::exception& e)
            {
                
                *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
                //std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (const std::bad_alloc e)
            {
                
                *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
                //std::cout << e.what() << std::flush;
                //MPI_Abort(MPI_COMM_WORLD, -1);
            }
            catch (...)
            {
                *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
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
   if ( argv != NULL )
      MPE_Finish_log( argv[0] );
   if ( str != "" )
      MPE_Finish_log( str.c_str() );
   else
      MPE_Finish_log( "TestLog" );
#endif

   MPI_Finalize();
   return 0;
}
