//#define MPI_LOGGING

//Martin Branch

#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <sstream>
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

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"

#include "utilities/math/Math.h"
#include "utilities/communication.h"
#include "utilities/transformator/TransformatorImp.h"

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

	const real L = 1.0;

    const real Re = 1000;

    const uint nx = 64;

	const real vx = 0.05; // LB units
	const real vy = 0.05; // LB units

	const real velocity = sqrt(vx*vx + vy*vy);

    const real viscosity = nx * velocity / Re; // LB units

    *logging::out << logging::Logger::INFO_HIGH << "velocity = " << velocity << " s\n";

    *logging::out << logging::Logger::INFO_HIGH << "viscosity = " << viscosity << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	real dx = L / real(nx);

	gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L,
								0.5 * L,  0.5 * L,  0.5 * L, dx);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(LBM, false); // buildGrids() has to be called before setting the BCs!!!!

	gridBuilder->setVelocityBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
	gridBuilder->setVelocityBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
	gridBuilder->setVelocityBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
	gridBuilder->setVelocityBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
	gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  vx,  vy, 0.0);
	gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

	//////////////////////////////////////////////////////////////////////////

	SPtr<Parameter> para = Parameter::make(configData, comm);
	//SPtr<Parameter> para = Parameter::make();

	SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

	SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setVelocity( velocity );

    para->setViscosity( viscosity );

    para->setVelocityRatio( 1.0 / velocity );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);
    sim.run();
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
