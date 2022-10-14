
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
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"

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

const real reference_diameter = 126.0; // diameter in m

const real L_x = 10*reference_diameter;
const real L_y = 6*reference_diameter;
const real L_z = 6*reference_diameter;

const real viscosity = 1.56e-5;

const real velocity  = 9.0;

const real mach = 0.1;

const uint nodes_per_diameter = 16;

std::string path(".");

std::string simulationName("ActuatorLine");

const float tOut = 100;
const float tEnd = 280; // total time of simulation in s
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
    gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
    auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real dx = reference_diameter/real(nodes_per_diameter);

	gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
							   L_x,  L_y,  L_z, dx);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(lbmOrGks, false); // buildGrids() has to be called before setting the BCs!!!!

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vf::basics::ConfigurationFile config;
    config.load(configPath);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
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

    para->setMaxLevel(1);


    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio( dx / dt );
    para->setViscosityRatio( dx*dx/dt );
    para->setMainKernel("CumulantK17Almighty");

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        rho = (real)0.0;
        vx  = velocityLB;
        vy  = (real)0.0;
        vz  = (real)0.0;
    });

    para->setTimestepOut( uint(tOut/dt) );
    para->setTimestepEnd( uint(tEnd/dt) );

    para->setIsBodyForce( true );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gridBuilder->setVelocityBoundaryCondition(SideType::MX,  velocityLB,  0.0, 0.0);

    gridBuilder->setVelocityBoundaryCondition(SideType::MY,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ,  velocityLB,  0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ,  velocityLB,  0.0, 0.0);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real turbPos[3] = {3*reference_diameter, 3*reference_diameter, 3*reference_diameter};
    real epsilon = 5.f; // width of gaussian smearing
    real density = 1.225f;
    int level = 0;
    uint nBlades = 3;
    uint nBladeNodes = 32;

    SPtr<ActuatorLine> actuator_line =SPtr<ActuatorLine>( new ActuatorLine(nBlades, density, nBladeNodes, epsilon, turbPos[0], turbPos[1], turbPos[2], reference_diameter, level, dt, dx) );
    para->addActuator( actuator_line );

    SPtr<PointProbe> pointProbe = SPtr<PointProbe>( new PointProbe("pointProbe", para->getOutputPath(), 100, 1, 500, 100) );
    std::vector<real> probeCoordsX = {reference_diameter,2*reference_diameter,5*reference_diameter};
    std::vector<real> probeCoordsY = {3*reference_diameter,3*reference_diameter,3*reference_diameter};
    std::vector<real> probeCoordsZ = {3*reference_diameter,3*reference_diameter,3*reference_diameter};
    pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);
    // pointProbe->addProbePointsFromXNormalPlane(2*D, 0.0, 0.0, L_y, L_z, (uint)L_y/dx, (uint)L_z/dx);

    pointProbe->addStatistic(Statistic::Means);
    pointProbe->addStatistic(Statistic::Variances);
    para->addProbe( pointProbe );

    SPtr<PlaneProbe> planeProbe = SPtr<PlaneProbe>( new PlaneProbe("planeProbe", para->getOutputPath(), 100, 500, 100, 100) );
    planeProbe->setProbePlane(5*reference_diameter, 0, 0, dx, L_y, L_z);
    planeProbe->addStatistic(Statistic::Means);
    para->addProbe( planeProbe );


    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    auto gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
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
