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

#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/DataStructureInitializer/GridReaderFiles/GridReader.h"
#include "gpu/core/Factories/BoundaryConditionFactory.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandalone.h"
#include "gpu/core/PreCollisionInteractor/Probes/PointProbe.h"
#include "gpu/core/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          U s e r    s e t t i n g s
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

    const real referenceDiameter = config.getValue<real>("ReferenceDiameter");
    const uint nodesPerDiameter = config.getValue<uint>("NodesPerDiameter");
    const real velocity = config.getValue<real>("Velocity");


    const real lengthX = 24*referenceDiameter;
    const real lengthY = 6*referenceDiameter;
    const real lengthZ = 6*referenceDiameter;

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
        
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(), &config);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    GridScalingFactory scalingFactory  = GridScalingFactory();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const real dx = referenceDiameter/real(nodesPerDiameter);

    std::vector<real>turbinePositionsX{3.f*referenceDiameter};
    std::vector<real>turbinePositionsY{0.5f*lengthY};
    std::vector<real>turbinePositionsZ{0.5f*lengthY};

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

	gridBuilder->addCoarseGrid(0.0, 0.0, 0.0,
							   lengthX,  lengthY,  lengthZ, dx);

    gridBuilder->setNumberOfLayers(4,0);
    gridBuilder->addGrid( std::make_shared<Cuboid>( turbinePositionsX[0]-1.5*referenceDiameter,  turbinePositionsY[0]-1.5*referenceDiameter,  turbinePositionsZ[0]-1.5*referenceDiameter, 
                                                    turbinePositionsX[0]+10.0*referenceDiameter, turbinePositionsY[0]+1.5*referenceDiameter,  turbinePositionsZ[0]+1.5*referenceDiameter) , 1 );
    para->setMaxLevel(2);
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

	gridBuilder->setPeriodicBoundaryCondition(false, false, false);

	gridBuilder->buildGrids(false); // buildGrids() has to be called before setting the BCs!!!!

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
    const real smearingWidth = dx*exp2(-level)*1.5; // width of gaussian smearing
    const real density = 1.225f;
    const uint nBladeNodes = 32;
    const real tipspeedRatio = 7.5f; // tipspeed ratio = angular vel * radius / inflow vel
    const std::vector<real> rotorSpeeds = {2*tipspeedRatio*velocity/referenceDiameter};    

    SPtr<ActuatorFarmStandalone> actuatorFarm = std::make_shared<ActuatorFarmStandalone>(referenceDiameter, nBladeNodes, turbinePositionsX, turbinePositionsZ, turbinePositionsZ, rotorSpeeds, density, smearingWidth, level, dt, dx);
    actuatorFarm->enableOutput("ALM", uint(tStartOutProbe/dt), uint(tOutProbe/dt));
    para->addActuator( actuatorFarm );


    SPtr<PointProbe> pointProbe = std::make_shared<PointProbe>("pointProbe", para->getOutputPath(), 100, 1, 500, 100, false);
    std::vector<real> probeCoordsX = {referenceDiameter,2*referenceDiameter,5*referenceDiameter};
    std::vector<real> probeCoordsY = {3*referenceDiameter,3*referenceDiameter,3*referenceDiameter};
    std::vector<real> probeCoordsZ = {3*referenceDiameter,3*referenceDiameter,3*referenceDiameter};

    pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);
    // pointProbe->addProbePointsFromXNormalPlane(2*D, 0.0, 0.0, lengthY, lengthZ, (uint)lengthY/dx, (uint)lengthZ/dx);

    pointProbe->addStatistic(Statistic::Means);
    pointProbe->addStatistic(Statistic::Variances);
    para->addProbe( pointProbe );

    SPtr<PointProbe> timeseriesProbe = std::make_shared<PointProbe>("timeProbe", para->getOutputPath(), 100, 1, 500, 100, true);
    timeseriesProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);
    timeseriesProbe->addStatistic(Statistic::Instantaneous);
    timeseriesProbe->addStatistic(Statistic::Means);
    para->addProbe( timeseriesProbe );

    SPtr<PlaneProbe> planeProbe = std::make_shared<PlaneProbe>("planeProbe", para->getOutputPath(), 100, 500, 100, 100);
    planeProbe->setProbePlane(5*referenceDiameter, 0, 0, dx, lengthY, lengthZ);
    planeProbe->addStatistic(Statistic::Means);
    para->addProbe( planeProbe );


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
