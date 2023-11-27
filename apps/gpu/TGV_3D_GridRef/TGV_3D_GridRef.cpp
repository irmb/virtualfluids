//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TGV_3D.cpp
//! \ingroup Applications
//! \author Martin Schoenherr
//=======================================================================================
#define _USE_MATH_DEFINES
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

//////////////////////////////////////////////////////////////////////////

#include "parallel/MpiCommunicator.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Kernel/KernelTypes.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

uint dtPerL = 500;

uint nx = 64;
uint gpuIndex = 0;

bool useLimiter = false;

std::string kernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

std::string path("D:/out/TGV_3D/"); //MOLLOK

std::string simulationName("TGV_3D_Gridref_noSqPress");
//////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    //vf::gpu::Communicator& communicator = vf::gpu::MpiCommunicator::getInstance();
    vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();

    vf::basics::ConfigurationFile config;
    config.load(configPath);
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real PI = 3.141592653589793238462643383279;

    real L = nx / ( 2.0 * PI );

    const real velocity = 64.0 / ( dtPerL * 2.0 * PI );

    const real viscosity = nx / ( 2.0 * PI ) * velocity / Re;

    VF_LOG_INFO("velocity = {}", velocity);
    VF_LOG_INFO("viscosity = {}", viscosity);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    real dx = 2.0 * PI / real(nx);

    gridBuilder->addCoarseGrid(-PI, -PI, -PI,
                                PI,  PI,  PI, dx);

    gridBuilder->setNumberOfLayers(0, 0);

    auto fineGrid = std::make_shared<Cuboid>(-PI * 0.5, -PI * 0.5, -PI * 0.5, 
                                     0.0,  PI * 0.5,       0.0);

    gridBuilder->addGrid(fineGrid, 1);

    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

    gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

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
 //   para->setPathAndFilename(_path.str() + "/" + _prefix.str());

    //////////////////////////////////////////////////////////////////////////

    {
        std::stringstream _path;

        _path << path;
        _path << kernel;
        _path << "SingleGPU";

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

    para->setDevices(std::vector<uint>{gpuIndex});

    //////////////////////////////////////////////////////////////////////////

    para->setOutputPath( path );
    para->setOutputPrefix( simulationName );

    para->setPrintFiles(true);

    para->setTimestepEnd(40 * lround(L / velocity));
    para->setTimestepOut(5 * lround(L / velocity));
    //para->setTimestepOut(lround(L / velocity));
 //   para->setTimestepEnd(2048);
    //para->setTimestepOut(512);
 //   para->setTimestepStartOut(500);

    para->setVelocityLB( velocity );

    para->setViscosityLB( viscosity );

    para->setVelocityRatio( 1.0 / velocity );

    para->setDensityRatio(1.0);

    para->setInitialCondition( [&]( real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz){

        real a = 1.0;
        real b = 1.0;
        real c = 1.0;

        rho = 3.0 * ((velocity * velocity) / 16.0 * ( cos( 2.0 * a * coordX ) + cos( 2.0 * b * coordY ) ) * ( cos( 2.0 * c * coordZ ) + 2.0 ) );
        vx  =  velocity * sin( a * coordX ) * cos( b * coordY ) * cos( c * coordZ );
        vy  = -velocity * cos( a * coordX ) * sin( b * coordY ) * cos( c * coordZ );
        vz  = 0.0;

    } );

    para->configureMainKernel( kernel );

    if( !useLimiter )
        para->setQuadricLimiters( 1000000.0, 1000000.0, 1000000.0 );

    para->setUseInitNeq( true );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);
    //SPtr<GridProvider> gridGenerator = GridProvider::makeGridReader(FILEFORMAT::BINARY, para, cudaMemoryManager);

    SPtr<FileWriter> fileWriter = SPtr<FileWriter>(new FileWriter());
    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
    sim.run();

    //sim.addKineticEnergyAnalyzer( 10 );
    //sim.addEnstrophyAnalyzer( 10 );

    //sim.run();
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

            multipleLevel(targetPath + "config.txt");

            //////////////////////////////////////////////////////////////////////////
        }
        catch (const std::bad_alloc& e)
        {
            std::cout << "Bad alloc: " << e.what() << std::flush;
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << std::flush;
            //MPI_Abort(MPI_COMM_WORLD, -1);
        }
        catch (...)
        {
            std::cout << "unknown exeption" << std::endl;
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
