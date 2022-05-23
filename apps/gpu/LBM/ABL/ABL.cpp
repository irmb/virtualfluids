
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
#include "PointerDefinitions.h"

#include "Core/StringUtilities/StringUtil.h"

#include "Core/VectorTypes.h"

#include <basics/config/ConfigurationFile.h>
#include "lbm/constants/NumericConstants.h"

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
#include "VirtualFluids_GPU/PreCollisionInteractor/ActuatorLine.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"

#include "VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h"
#include "VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/PrecursorWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/VelocitySetter.h"

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

const real reference_height = 2000.f; // boundary layer height in m

const real L_x = 6*reference_height;
const real L_y = 4*reference_height;
const real L_z = reference_height;

const real viscosity = 1.56e-5;

const real z0 = 0.1; // roughness length in m
const real u_star = 0.4; //friction velocity in m/s
const real kappa = 0.4; // von Karman constant 


const real velocity  = 0.5*u_star/kappa*log(L_z/z0);

const real mach = 0.1;

const uint nodes_per_height = 64;

std::string path(".");

std::string simulationName("ABL");

real turnOverTime = L_z/velocity;
const float tOut = 10*turnOverTime;
const float tEnd = 1000*turnOverTime; // total time of simulation in s


const float tStartPrecursor = 1*turnOverTime;
const uint nTReadPrecursor = 10;
const uint nTWritePrecursor = 10;
bool readPrecursor = true;
bool writePrecursor = false;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    auto gridFactory = GridFactory::make();
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real dx = reference_height/real(nodes_per_height);

	gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
							   L_x, L_y, L_z, dx);

    gridBuilder->setPeriodicBoundaryCondition(false, true, false);
    
	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    vf::basics::ConfigurationFile config;
    config.load(configPath);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^
    SPtr<Parameter> para = std::make_shared<Parameter>(config, communicator.getNummberOfProcess(), communicator.getPID());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real dt = dx * mach / (sqrt(3) * velocity);

    const real velocityLB = velocity * dt / dx; // LB units

    const real viscosityLB = viscosity * dt / (dx * dx); // LB units

    const real pressureGradient = u_star * u_star / reference_height ;
    const real pressureGradientLB = pressureGradient * (dt*dt)/dx; // LB units

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("viscosity [10^8 dx^2/dt] = {}", viscosityLB*1e8);



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setDevices(std::vector<uint>{(uint)0});

    para->setOutputPrefix( simulationName );

    para->setFName(para->getOutputPath() + "/" + para->getOutputPrefix());

    para->setPrintFiles(true);


    para->setVelocity(velocityLB);
    para->setViscosity(viscosityLB);
    para->setVelocityRatio( dx / dt );
    para->setViscosityRatio( dx*dx/dt );
    para->setMainKernel("TurbulentViscosityCumulantK17CompChim");

    para->setForcing(pressureGradientLB, 0, 0);

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = (u_star/c4o10 * log(coordZ/z0) + c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/L_z)/(pow(coordZ/L_z,c2o1)+c1o1))  * dt / dx; 
        vy  = c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/L_z)/(pow(coordZ/L_z,c2o1)+c1o1)  * dt / dx; 
        vz  = c8o1*u_star/c4o10*(sin(cPi*c8o1*coordY/L_z)*sin(cPi*c8o1*coordZ/L_z)+sin(cPi*c8o1*coordX/L_x))/(pow(L_z*c1o2-coordZ, c2o1)+c1o1) * dt / dx;
    });

    para->setTOut( uint(tOut/dt) );
    para->setTEnd( uint(tEnd/dt) );

    para->setIsBodyForce( false );
    para->setUseAMD( true );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This order works!
    // gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);
    // gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.f,0.f,1.f);
    // gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);

    gridBuilder->setNoSlipBoundaryCondition(SideType::PZ);
    gridBuilder->setSlipBoundaryCondition(SideType::MZ, 0.f, 0.f, 0.f);
    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    // gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.f, 0.f, 0.f);

    if(false)
    {
        auto precursor = SPtr<VTKFileCollection>( new VTKFileCollection("precursor/Precursor") );
        auto velocitySetter = SPtr<VelocitySetter>( new VelocitySetter(precursor, nTReadPrecursor) );

        for(int level=0; level<para->getMaxLevel()+1; level++)
        {
            auto velBC = gridBuilder->getBoundaryCondition(SideType::MX, level);
    
            velocitySetter->setBCArrays(gridBuilder->getGrid(level), velBC);
        }
        para->addActuator(velocitySetter);
    }
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.f);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    SPtr<CudaMemoryManager> cudaMemoryManager = CudaMemoryManager::make(para);

    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager);

    real x_pos = L_x/2;

    if(writePrecursor)
    {
        auto precursor_writer = SPtr<PrecursorWriter>( new PrecursorWriter("Precursor", "precursor", x_pos, 0, L_y, 0, L_z, uint(tStartPrecursor/dt), nTWritePrecursor, 10000) );
        para->addProbe( precursor_writer );
    }

    Simulation sim(communicator);
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
    if ( argv != NULL )
    {
        try
        {
            vf::logging::Logger::initalizeLogger();

            if( argc > 1){ path = argv[1]; }

            multipleLevel(path + "/configABL.txt");
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
    return 0;
}
