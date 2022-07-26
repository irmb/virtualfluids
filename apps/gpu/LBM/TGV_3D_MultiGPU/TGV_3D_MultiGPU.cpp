//#define MPI_LOGGING

//Martin Branch

#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

//#include "metis.h"

#include "Core/LbmOrGks.h"
#include "Core/StringUtilities/StringUtil.h"
#include "basics/config/ConfigurationFile.h"

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// from https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
real Re =  1600.0;

uint dtPerL = 250;

uint nx = 64;
uint gpuIndex = 0;

bool useLimiter = false;
bool useWale = false;

int mpirank;
int mpiWorldSize;

std::string kernel( "CumulantK20Comp" );

//std::string path("F:/Work/Computations/out/TaylorGreen3DNew/"); //LEGOLAS
//std::string path("results/"); //PHOENIX
//std::string path("E:/DrivenCavity/results/"); //TESLA03
std::string path("E:/TaylorGreen3D/"); //AzultecPC

std::string simulationName("TGV_3D");
//////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int sideLengthX, sideLengthY, sideLengthZ, rankX, rankY, rankZ;

    
    if      (mpiWorldSize == 1 ) { sideLengthX = 1; sideLengthY = 1; sideLengthZ = 1; }
    else if (mpiWorldSize == 2 ) { sideLengthX = 2; sideLengthY = 1; sideLengthZ = 1; }
    else if (mpiWorldSize == 4 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 1; }
    else if (mpiWorldSize == 8 ) { sideLengthX = 2; sideLengthY = 2; sideLengthZ = 2; }
    else if (mpiWorldSize == 16) { sideLengthX = 4; sideLengthY = 2; sideLengthZ = 2; }
    else if (mpiWorldSize == 32) { sideLengthX = 4; sideLengthY = 4; sideLengthZ = 2; }

    rankX =   mpirank %   sideLengthX;
    rankY = ( mpirank % ( sideLengthX * sideLengthY ) ) /   sideLengthX;
    rankZ =   mpirank                                   / ( sideLengthY * sideLengthX );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    logging::Logger::addStream(&std::cout);
    
    std::ofstream logFile( path + simulationName + "_rank_" + std::to_string(mpirank) + ".log" );
    logging::Logger::addStream(&logFile);

    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    //UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG5");

    auto gridFactory = GridFactory::make();
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::RAYCASTING);
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    //gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_UNDER_TRIANGLE);

    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);
    
    vf::basics::ConfigurationFile config;
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(config, communicator.getNummberOfProcess(), communicator.getPID());
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    *logging::out << logging::Logger::INFO_HIGH << "SideLength = " << sideLengthX << " " << sideLengthY << " " << sideLengthZ << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "rank       = " << rankX << " " << rankY << " " << rankZ << "\n";

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

    real LX = 2.0 * PI / double(sideLengthX);
    real LY = 2.0 * PI / double(sideLengthY);
    real LZ = 2.0 * PI / double(sideLengthZ);

    *logging::out << logging::Logger::INFO_HIGH << "L       = " << LX << " " << LY << " " << LZ << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	real dx = 2.0 * PI / real(nx);

    real xOverlap = ( sideLengthX == 1 ) ? 0.0 : 5.0*dx;
    real yOverlap = ( sideLengthY == 1 ) ? 0.0 : 5.0*dx;
    real zOverlap = ( sideLengthZ == 1 ) ? 0.0 : 5.0*dx;

    *logging::out << logging::Logger::INFO_HIGH << "Domain       = " <<  rankX*LX    - PI - xOverlap << " " <<  rankY*LY    - PI - yOverlap << " " <<  rankZ*LZ    - PI - zOverlap << "\n";
    *logging::out << logging::Logger::INFO_HIGH << "Domain       = " << (rankX*LX+1) - PI + xOverlap << " " << (rankY*LY+1) - PI + yOverlap << " " << (rankZ*LZ+1) - PI + zOverlap << "\n";

    gridBuilder->addCoarseGrid(  rankX   *LX - PI - xOverlap,      rankY   *LY - PI - yOverlap,      rankZ   *LZ - PI - zOverlap,
                                (rankX+1)*LX - PI + xOverlap,     (rankY+1)*LY - PI + yOverlap,     (rankZ+1)*LZ - PI + zOverlap, dx);

    gridBuilder->setSubDomainBox( std::make_shared<BoundingBox>( rankX*LX - PI, (rankX+1)*LX - PI, 
                                                                 rankY*LY - PI, (rankY+1)*LY - PI,
                                                                 rankZ*LZ - PI, (rankZ+1)*LZ - PI  ) );

    gridBuilder->setPeriodicBoundaryCondition(sideLengthX == 1, sideLengthY == 1, sideLengthZ == 1);

	gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

    if( mpiWorldSize > 1 )
    {
        int rankPX = ( (rankX + 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMX = ( (rankX - 1 + sideLengthX) % sideLengthX ) +    rankY                                    * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPY =    rankX                                    + ( (rankY + 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankMY =    rankX                                    + ( (rankY - 1 + sideLengthY) % sideLengthY ) * sideLengthX +    rankZ                                    * sideLengthX * sideLengthY;
        int rankPZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ + 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;
        int rankMZ =    rankX                                    +    rankY                                    * sideLengthX + ( (rankZ - 1 + sideLengthZ) % sideLengthZ ) * sideLengthX * sideLengthY;

        if( sideLengthX > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PX, GKS );
        if( sideLengthX > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MX, GKS );
        if( sideLengthY > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PY, GKS );
        if( sideLengthY > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MY, GKS );
        if( sideLengthZ > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::PZ, GKS );
        if( sideLengthZ > 1 ) gridBuilder->findCommunicationIndices( CommunicationDirections::MZ, GKS );

        if( sideLengthX > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MX, rankMX);
        if( sideLengthY > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MY, rankMY);
        if( sideLengthZ > 1 ) gridBuilder->setCommunicationProcess ( CommunicationDirections::MZ, rankMZ);

        if( sideLengthX > 1 && rankMX != rankPX ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PX, rankPX);
        if( sideLengthY > 1 && rankMY != rankPY ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PY, rankPY);
        if( sideLengthZ > 1 && rankMZ != rankPZ ) gridBuilder->setCommunicationProcess ( CommunicationDirections::PZ, rankPZ);

        if( rankMX == rankPX ) gridBuilder->getGrid(0)->repairCommunicationIndices(CommunicationDirections::MX);
        if( rankMY == rankPY ) gridBuilder->getGrid(0)->repairCommunicationIndices(CommunicationDirections::MY);
        if( rankMZ == rankPZ ) gridBuilder->getGrid(0)->repairCommunicationIndices(CommunicationDirections::MZ);

        *logging::out << logging::Logger::INFO_HIGH << "neighborRanks = " << rankPX << " " << rankMX << " " << rankPY << " " << rankMY << " " << rankPZ << " " << rankMZ << "\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////

    SimulationFileWriter::write(path + "grid/" + std::to_string(mpirank) + "/", gridBuilder, FILEFORMAT::BINARY);

    //////////////////////////////////////////////////////////////////////////

    std::vector<uint> devices( mpiWorldSize );

    std::iota(devices.begin(), devices.end(), 0);

    //para->setDevices(std::vector<uint>{0,1});
    para->setDevices(devices);
	
	para->setMaxDev(mpiWorldSize);

    //////////////////////////////////////////////////////////////////////////
    
    para->setOutputPath( path );
    para->setOutputPrefix( simulationName );

    para->setPathAndFilename(para->getOutputPath() + "/" + para->getOutputPrefix());

    para->setPrintFiles(true);

 //   para->setTimestepEnd( 40 * lround(L/velocity) );	
	//para->setTimestepOut(  5 * lround(L/velocity) );
	para->setTimestepOut(  100  );

    para->setTimestepEnd( 1000 );	
	//para->setTimestepOut(    1 );

    para->setVelocityLB( velocity );

    para->setViscosityLB( viscosity );

    para->setVelocityRatio( 1.0 / velocity );

    para->setInitialCondition( [&]( real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz){

        real a = 1.0;
        real b = 1.0;
        real c = 1.0;

        rho = 3.0 * ((velocity * velocity) / 16.0 * ( cos( 2.0 * a * coordX ) + cos( 2.0 * b * coordY ) ) * ( cos( 2.0 * c * coordZ ) + 2.0 ) );
        vx  =  velocity * sin( a * coordX ) * cos( b * coordY ) * cos( c * coordZ );
        vy  = -velocity * cos( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );
        vz  = 0.0;

        //rho = mpirank;
        //vx  = 0.0;
        //vy  = 0.0;
        //vz  = 0.0;

    } );

    para->setMainKernel(kernel);

    if( !useLimiter )
        para->setQuadricLimiters( 1000000.0, 1000000.0, 1000000.0 );

    if( useWale )
        para->setUseWale( true );

    para->setUseInitNeq( true );

	if (kernel == "CumulantK18Comp" || kernel == "CumulantK20Comp")
		para->setIsF3(true);
	else
		para->setIsF3(false);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
    //SPtr<GridProvider> gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);

    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
    sim.run();
    
    sim.addKineticEnergyAnalyzer( 10 );
    sim.addEnstrophyAnalyzer( 10 );

    sim.run();

    logFile.close();
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
            MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

            //////////////////////////////////////////////////////////////////////////
			std::string targetPath( __FILE__ );

#ifdef _WIN32
			targetPath = targetPath.substr(0, targetPath.find_last_of('\\') + 1);
#else
			targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
#endif

            //////////////////////////////////////////////////////////////////////////

            if( cmdOptionExists( argv, argv+argc, "--Re" ) )
                Re = atof( getCmdOption( argv, argv+argc, "--Re" ) );

            if( cmdOptionExists( argv, argv+argc, "--nx" ) )
                nx = atoi( getCmdOption( argv, argv+argc, "--nx" ) );

            if( cmdOptionExists( argv, argv+argc, "--dtPerL" ) )
                dtPerL = atoi( getCmdOption( argv, argv+argc, "--dtPerL" ) );

            if( cmdOptionExists( argv, argv+argc, "--kernel" ) )
                kernel = getCmdOption( argv, argv+argc, "--kernel" );

            if( cmdOptionExists( argv, argv+argc, "--gpu" ) )
                gpuIndex = atoi( getCmdOption( argv, argv+argc, "--gpu" ) );

            if( cmdOptionExists( argv, argv+argc, "--useLimiter" ) )
                useLimiter = true;

            if( cmdOptionExists( argv, argv+argc, "--useWale" ) )
                useWale = true;

            //////////////////////////////////////////////////////////////////////////

            {
                std::stringstream _path;

                _path << path;
                _path << kernel;
                _path << "MultiGPU/";

                if (useLimiter) _path << "_Limiter";

                path = _path.str();
            }

            //////////////////////////////////////////////////////////////////////////

            {
                std::stringstream _simulationName;

                _simulationName << simulationName;
                _simulationName << "_nx_" << nx;
                _simulationName << "_dtPerL_" << dtPerL << "_";

                simulationName = _simulationName.str();
            }

            //////////////////////////////////////////////////////////////////////////

			multipleLevel(targetPath + "config.txt");

            //////////////////////////////////////////////////////////////////////////
		}
        catch (const std::bad_alloc& e)
        {
                
            *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
            //std::cout << e.what() << std::flush;
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }
        catch (const std::exception& e)
        {
                
            *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
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
