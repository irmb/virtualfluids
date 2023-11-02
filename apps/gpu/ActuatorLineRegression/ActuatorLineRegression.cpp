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
//! \file ActuatorLine.cpp
//! \ingroup ActuatorLine
//! \author Henry Korb, Henrik Asmuth
//=======================================================================================
#define _USE_MATH_DEFINES
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/StringUtilities/StringUtil.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"

#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/TransientBCSetter/TransientBCSetter.h"


//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/ActuatorFarm.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/Probe.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"
#include "VirtualFluids_GPU/Factories/GridScalingFactory.h"
#include "VirtualFluids_GPU/Kernel/Utilities/KernelTypes.h"

#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Actuator Line app for regression tests
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string path(".");

std::string simulationName("ActuatorLine");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void multipleLevel(const std::string& configPath)
{
    vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();

    vf::basics::ConfigurationFile config;
    config.load(configPath);

    const real reference_diameter = config.getValue<real>("ReferenceDiameter");
    const uint nodes_per_diameter = config.getValue<uint>("NodesPerDiameter");
    const real velocity = config.getValue<real>("Velocity");

    const real L_x = 10 * reference_diameter;
    const real L_y = 4 * reference_diameter;
    const real L_z = 4 * reference_diameter;

    const real viscosity = 1.56e-5;

    const real mach = 0.1;


    const float tStartOut   = config.getValue<real>("tStartOut");
    const float tOut        = config.getValue<real>("tOut");
    const float tEnd        = config.getValue<real>("tEnd"); // total time of simulation

    const float tStartAveraging     = config.getValue<real>("tStartAveraging");
    const float tStartTmpAveraging  = config.getValue<real>("tStartTmpAveraging");
    const float tAveraging          = config.getValue<real>("tAveraging");
    const float tStartOutProbe      = config.getValue<real>("tStartOutProbe");
    const float tOutProbe           = config.getValue<real>("tOutProbe");
        
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory  = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real dx = reference_diameter/real(nodes_per_diameter);

    real turbPos[3] = {3.0f * reference_diameter, 0.0, 0.0};

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    gridBuilder->addCoarseGrid(0.0, -0.5*L_y, -0.5*L_z,
                               L_x,  0.5*L_y,  0.5*L_z, dx);

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(false); // buildGrids() has to be called before setting the BCs!!!!

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const real dt = dx * mach / (sqrt(3) * velocity);

    const real velocityLB = velocity * dt / dx; // LB units

    const real viscosityLB = viscosity * dt / (dx * dx); // LB units

    VF_LOG_INFO("dx = {}m", dx);
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
    para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

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

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MY, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, velocityLB, 0.0, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, velocityLB, 0.0, 0.0);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

    SPtr<TurbulenceModelFactory> tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int level = 0; // grid level at which the turbine samples velocities and distributes forces
    const real smearing_width = dx*exp2(-level)*2; // width of gaussian smearing
    VF_LOG_INFO("smearing_width = {}m", smearing_width);
    const real density = 1.225f;
    const uint nBlades = 3;
    const uint nBladeNodes = 32;
    const real tipspeed_ratio = 7.5f; // tipspeed ratio = angular vel * radius / inflow vel
    const real rotor_speed = 2*tipspeed_ratio*velocity/reference_diameter;


    SPtr<ActuatorFarm> actuator_farm = std::make_shared<ActuatorFarm>(nBlades, density, nBladeNodes, smearing_width, level, dt, dx, true);
    std::vector<real> bladeRadii;
    real dr = reference_diameter/(nBladeNodes*2);
    for(uint node=0; node<nBladeNodes; node++){ bladeRadii.emplace_back(dr*(node+1)); }
    actuator_farm->addTurbine(turbPos[0], turbPos[1], turbPos[2], reference_diameter, rotor_speed, 0, 0, bladeRadii);
    para->addActuator( actuator_farm );

    std::vector<real> planePositions = {-1*reference_diameter, 1*reference_diameter, 3*reference_diameter};

    for(int i=0; i < planePositions.size(); i++)
    {
        SPtr<PlaneProbe> planeProbe = std::make_shared<PlaneProbe>("planeProbe_" + std::to_string(i), para->getOutputPath(), tStartTmpAveraging/dt, tAveraging/dt, tStartOutProbe/dt, tOutProbe/dt);
        planeProbe->setProbePlane(turbPos[0]+planePositions[i], -0.5 * L_y, -0.5 * L_z, dx, L_y, L_z);
        planeProbe->addStatistic(Statistic::Means);
        planeProbe->addStatistic(Statistic::Variances);
        planeProbe->addStatistic(Statistic::Instantaneous);
        para->addProbe( planeProbe );
    }
    SPtr<PlaneProbe> planeProbeVert = std::make_shared<PlaneProbe>("planeProbeVertical", para->getOutputPath(), tStartTmpAveraging/dt, tAveraging/dt, tStartOutProbe/dt, tOutProbe/dt);
    planeProbeVert->setProbePlane(0, turbPos[1], -0.5 * L_z, L_x, dx, L_z);
    planeProbeVert->addStatistic(Statistic::Means);
    planeProbeVert->addStatistic(Statistic::Variances);
    planeProbeVert->addStatistic(Statistic::Instantaneous);
    para->addProbe( planeProbeVert );

    SPtr<PlaneProbe> planeProbeHorz = std::make_shared<PlaneProbe>("planeProbeHorizontal", para->getOutputPath(), tStartTmpAveraging/dt, tAveraging/dt, tStartOutProbe/dt, tOutProbe/dt);
    planeProbeHorz->setProbePlane(0, -0.5 * L_y, turbPos[2], L_x, L_y, dx);
    planeProbeHorz->addStatistic(Statistic::Means);
    planeProbeHorz->addStatistic(Statistic::Variances);
    planeProbeHorz->addStatistic(Statistic::Instantaneous);
    para->addProbe( planeProbeHorz );


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
            vf::logging::Logger::initializeLogger();

            if( argc > 1){ path = argv[1]; }

            multipleLevel(path + "/apps/gpu/ActuatorLineRegression/configActuatorLine.txt");
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
