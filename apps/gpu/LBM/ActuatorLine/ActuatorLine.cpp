
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

#include <logger/Logger.h>


//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"

#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/TransientBCSetter/TransientBCSetter.h"


//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/ActuatorFarm.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"
#include "VirtualFluids_GPU/Factories/GridScalingFactory.h"

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

std::string path(".");

std::string simulationName("ActuatorLine");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{

    logging::Logger::addStream(&std::cout);
    logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
    logging::Logger::timeStamp(logging::Logger::ENABLE);
    logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

    vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();

    auto gridFactory = GridFactory::make();
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    vf::basics::ConfigurationFile config;
    config.load(configPath);

    const real reference_diameter = config.getValue<real>("ReferenceDiameter");
    const uint nodes_per_diameter = config.getValue<uint>("NodesPerDiameter");
    const real velocity = config.getValue<real>("Velocity");


    const real L_x = 24*reference_diameter;
    const real L_y = 6*reference_diameter;
    const real L_z = 6*reference_diameter;

    const real viscosity = 1.56e-5;

    const real mach = 0.1;


    const float tStartOut   = config.getValue<real>("tStartOut");
    const float tOut        = config.getValue<real>("tOut");
    const float tEnd        = config.getValue<real>("tEnd"); // total time of simulation

    const float tStartAveraging     =  config.getValue<real>("tStartAveraging");
    const float tStartTmpAveraging  =  config.getValue<real>("tStartTmpAveraging");
    const float tAveraging          =  config.getValue<real>("tAveraging");
    const float tStartOutProbe      =  config.getValue<real>("tStartOutProbe");
    const float tOutProbe           =  config.getValue<real>("tOutProbe");
        
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory  = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real dx = reference_diameter/real(nodes_per_diameter);

    real turbPos[3] = {3*reference_diameter, 3*reference_diameter, 3*reference_diameter};

	gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
							   L_x,  L_y,  L_z, dx);

    gridBuilder->setNumberOfLayers(4,0);
    gridBuilder->addGrid( new Cuboid(   turbPos[0]-1.5*reference_diameter,  turbPos[1]-1.5*reference_diameter,  turbPos[2]-1.5*reference_diameter, 
                                        turbPos[0]+10.0*reference_diameter, turbPos[1]+1.5*reference_diameter,  turbPos[2]+1.5*reference_diameter) , 1 );
    para->setMaxLevel(2);
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real dt = dx * mach / (sqrt(3) * velocity);

    const real velocityLB = velocity * dt / dx; // LB units

    const real viscosityLB = viscosity * dt / (dx * dx); // LB units

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("viscosity [10^8 dx^2/dt] = {}", viscosityLB*1e8);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    para->setDevices(std::vector<uint>{(uint)0});

    para->setOutputPrefix( simulationName );

    para->setPrintFiles(true);

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio( dx / dt );
    para->setViscosityRatio( dx*dx/dt );
    para->setMainKernel("CumulantK17");

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = velocityLB;
        vy  = (real)0.0;
        vz  = (real)0.0;
    });

    para->setTimestepStartOut( uint(tStartOut/dt) );
    para->setTimestepOut( uint(tOut/dt) );
    para->setTimestepEnd( uint(tEnd/dt) );

    para->setIsBodyForce( true );
    para->setUseStreams( true );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gridBuilder->setVelocityBoundaryCondition(SideType::MX,  velocityLB,  0.0, 0.0);

    gridBuilder->setVelocityBoundaryCondition(SideType::MY,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  velocityLB,  0.0, 0.0);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

    SPtr<TurbulenceModelFactory> tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int level = 1; // grid level at which the turbine samples velocities and distributes forces
    const real epsilon = dx*exp2(-level)*1.5; // width of gaussian smearing
    const real density = 1.225f;
    const uint nBlades = 3;
    const uint nBladeNodes = 32;
    const real tipspeed_ratio = 7.5f; // tipspeed ratio = angular vel * radius / inflow vel
    const real omega = 2*tipspeed_ratio*velocity/reference_diameter;
    

    SPtr<ActuatorFarm> actuator_farm = std::make_shared<ActuatorFarm>(nBlades, density, nBladeNodes, epsilon, level, dt, dx, true);
    std::vector<real> bladeRadii;
    real dr = reference_diameter/(nBladeNodes*2);
    for(uint node=0; node<nBladeNodes; node++){ bladeRadii.emplace_back(dr*(node+1)); }
    actuator_farm->addTurbine(turbPos[0], turbPos[1], turbPos[2], reference_diameter, omega, 0, 0, bladeRadii);
    para->addActuator( actuator_farm );


    // SPtr<PointProbe> pointProbe = std::make_shared<PointProbe>("pointProbe", para->getOutputPath(), 100, 1, 500, 100);
    // std::vector<real> probeCoordsX = {reference_diameter,2*reference_diameter,5*reference_diameter};
    // std::vector<real> probeCoordsY = {3*reference_diameter,3*reference_diameter,3*reference_diameter};
    // std::vector<real> probeCoordsZ = {3*reference_diameter,3*reference_diameter,3*reference_diameter};
    // pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);
    // // pointProbe->addProbePointsFromXNormalPlane(2*D, 0.0, 0.0, L_y, L_z, (uint)L_y/dx, (uint)L_z/dx);

    // pointProbe->addStatistic(Statistic::Means);
    // pointProbe->addStatistic(Statistic::Variances);
    // para->addProbe( pointProbe );

    // SPtr<PlaneProbe> planeProbe = std::make_shared<PlaneProbe>("planeProbe", para->getOutputPath(), 100, 500, 100, 100);
    // planeProbe->setProbePlane(5*reference_diameter, 0, 0, dx, L_y, L_z);
    // planeProbe->addStatistic(Statistic::Means);
    // para->addProbe( planeProbe );


    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, tmFactory, &scalingFactory);
    sim.run();
}

int main( int argc, char* argv[])
{
    if ( argv != NULL )
    {
        try
        {
            vf::logging::Logger::initalizeLogger();

            if( argc > 1){ path = argv[1]; }

            multipleLevel(path + "/configActuatorLine.txt");
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
