
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
#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"
#include "Core/Logger/Logger.h"

#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>


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
#include "VirtualFluids_GPU/Visitor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/Visitor/Probes/PlaneProbe.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

LbmOrGks lbmOrGks = LBM;

const real D = 126.0; // diameter in m

const real L_x = 10*D;
const real L_y = 6*D;
const real L_z = 6*D;

const real viscosity = 1.56e-5;

const real velocity  = 9.0;

const real mach = 0.1;

const uint nodes_per_D = 32;

std::string path(".");

std::string simulationName("ActuatorLine");

const uint timeStepOut = 500;
const uint timeStepEnd = 1000;
const float tEnd = 200;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    auto gridBuilder = MultipleGridBuilder::makeShared();

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
        vf::gpu::Communicator* comm = vf::gpu::Communicator::getInstanz();

        vf::basics::ConfigurationFile config;
        config.load(configPath);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        SPtr<Parameter> para = std::make_shared<Parameter>(config, comm->getNummberOfProcess(), comm->getPID());

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        const real dt = dx * mach / (sqrt(3) * velocity);

        const real velocityLB = velocity * dt / dx; // LB units

        const real viscosityLB = viscosity * dt / (dx * dx); // LB units

        VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
        VF_LOG_INFO("viscosity [10^8 dx^2/dt] = {}", viscosityLB*1e8);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		para->setDevices(std::vector<uint>{(uint)0});

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
        para->setTEnd( uint(tEnd/dt) );

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

        PointProbe* pointProbe = new PointProbe("pointProbe", 100, 500, 100);
        std::vector<real> probeCoordsX = {D,2*D,5*D};
        std::vector<real> probeCoordsY = {3*D,3*D,3*D};
        std::vector<real> probeCoordsZ = {3*D,3*D,3*D};
        pointProbe->setProbePointsFromList(probeCoordsX,probeCoordsY,probeCoordsZ);
        // pointProbe->setProbePointsFromXNormalPlane(2*D, 0.0, 0.0, L_y, L_z, dx, dx);
        pointProbe->addPostProcessingVariable(PostProcessingVariable::Means);
        pointProbe->addPostProcessingVariable(PostProcessingVariable::Variances);
        para->addProbe( pointProbe );

        PlaneProbe* planeProbe = new PlaneProbe("planeProbe", 100, 500, 100);
        planeProbe->setProbePlane(5*D, 0, 0, dx, L_y, L_z);
        planeProbe->addPostProcessingVariable(PostProcessingVariable::Means);
        para->addProbe( planeProbe );




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

            vf::logging::Logger::initalizeLogger();

            if( argc > 1){ path = argv[1]; }

			multipleLevel(path + "/configActuatorLine.txt");

            //////////////////////////////////////////////////////////////////////////
		}
        catch (const spdlog::spdlog_ex &ex) {
            std::cout << "Log initialization failed: " << ex.what() << std::endl;
        }

        catch (const std::bad_alloc& e)
        { 
            VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
        }
        catch (const std::exception& e)
        {   
            VF_LOG_CRITICAL("exception: {}", e.what());
        }
        catch (...)
        {
            VF_LOG_CRITICAL("Unknown exception!");
        }
    }

    MPI_Finalize();
    return 0;
}
