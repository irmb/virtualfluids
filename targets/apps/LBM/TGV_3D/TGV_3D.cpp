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

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/Kernel/Utilities/Mapper/KernelMapper/KernelMapper.h"

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

#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/STLReaderWriter/STLReader.h"
#include "io/STLReaderWriter/STLWriter.h"

#include "utilities/math/Math.h"
#include "utilities/communication.h"
#include "utilities/transformator/TransformatorImp.h"

//////////////////////////////////////////////////////////////////////////
const real Re =  160000.0;

const uint dtPerL = 250*4;

const uint nx = 64;
const uint gpuIndex = 0;

std::string path("F:/Work/Computations/out/TaylorGreen3DNew/"); //LEGOLAS
//std::string path("E:/DrivenCavity/results/"); //TESLA03

std::string simulationName("TaylorGreen3D");
//////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath, uint nx, uint gpuIndex)
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

    std::shared_ptr<KernelMapper> kernelMapper = KernelMapper::getInstance();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real PI = 3.141592653589793238462643383279;

    real L = nx / ( 2.0 * PI );

    const real velocity = 64.0 / ( dtPerL * 2.0 * PI );

    const real viscosity = nx / ( 2.0 * PI ) * velocity / Re;

    *logging::out << logging::Logger::INFO_HIGH << "velocity = " << velocity << " s\n";

    *logging::out << logging::Logger::INFO_HIGH << "viscosity = " << viscosity << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	real dx = 2.0 * PI / real(nx);

	gridBuilder->addCoarseGrid(-PI, -PI, -PI,
								PI,  PI,  PI, dx);

	gridBuilder->setPeriodicBoundaryCondition(true, true, true);

	gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
	//std::stringstream _path;
 //   std::stringstream _prefix;

 //   //_path << "F:/Work/Computations/TaylorGreenVortex_3D/TGV_LBM/" << nx << "_Re_1.6e4";
 //   //_path << "F:/Work/Computations/TaylorGreenVortex_3D/TGV_LBM/" << nx << "_neqInit";
 //   _path << "F:/Work/Computations/TaylorGreenVortex_3D/TGV_LBM/Re_1600/AA2016/" << nx << "_FD_O8";

 //   //_path << "./results/AA2016/" << nx;
 //   //_path << "./results/CumOne/" << nx;
 //   //_path << "./results/F3_2018/" << nx;

 //   _prefix << "TGV_3D_" << nx << "_" ;

 //   para->setOutputPath(_path.str());
 //   para->setOutputPrefix(_prefix.str());
 //   para->setFName(_path.str() + "/" + _prefix.str());

 //   para->setDevices(std::vector<uint>{gpuIndex});

    //////////////////////////////////////////////////////////////////////////
    
    para->setOutputPath( path );
    para->setOutputPrefix( simulationName );

    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

    para->setPrintFiles(true);

    para->setTEnd( 40 * lround(L/velocity) );	
	para->setTOut(  5 * lround(L/velocity) );

    para->setVelocity( velocity );

    para->setViscosity( viscosity );

    para->setVelocityRatio( 1.0 / velocity );

    para->setInitialCondition( [&]( real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz){

        real a = 1.0;
        real b = 1.0;
        real c = 1.0;

        rho = 3.0 * ((velocity * velocity) / 16.0 * ( cos( 2.0 * a * coordX ) + cos( 2.0 * b * coordY ) ) * ( cos( 2.0 * c * coordZ ) + 2.0 ) );
        vx  =  velocity * sin( a * coordX ) * cos( b * coordY ) * cos( c * coordZ );
        vy  = -velocity * cos( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );
        vz  = 0.0;

    } );

    para->setMainKernel(kernelMapper->getEnum("CumulantAA2016CompSP27"));

    //para->setQuadricLimiters( 1000000.0, 0.01, 0.01 );
    //para->setQuadricLimiters( 0.01, 1000000.0, 1000000.0 );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);
    //SPtr<GridProvider> gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);

    Simulation sim;
    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    SPtr<KernelFactoryImp> kernelFactory = KernelFactoryImp::getInstance();
    SPtr<PreProcessorFactoryImp> preProcessorFactory = PreProcessorFactoryImp::getInstance();
    sim.setFactories(kernelFactory, preProcessorFactory);
    sim.init(para, gridGenerator, fileWriter, cudaMemoryManager);
    
    sim.run();
    sim.free();
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

            //if( argc > 1 ) gpuIndex = atoi( argv[1] );

            //if( argc > 2 ) nx = atoi( argv[2] );
			std::string targetPath;

			targetPath = __FILE__;

#ifdef _WIN32
			targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
			targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

			multipleLevel(targetPath + "config.txt", nx, gpuIndex);

            //////////////////////////////////////////////////////////////////////////
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

        //std::cout << "\nConfiguration file must be set!: lbmgm <config file>" << std::endl << std::flush;
        //MPI_Abort(MPI_COMM_WORLD, -1);
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
