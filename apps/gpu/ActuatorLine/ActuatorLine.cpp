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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup ActuatorLine
//! \ingroup gpu_apps
//! \{
//! \author Henry Korb, Henrik Asmuth, Anna Wellmann
//=======================================================================================

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>

#include <logger/Logger.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

//////////////////////////////////////////////////////////////////////////

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandalone.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"

//////////////////////////////////////////////////////////////////////////
const std::string defaultConfigFile = "actuatorline.cfg";

void run(const vf::basics::ConfigurationFile& config)
{
    //////////////////////////////////////////////////////////////////////////
    // Simulation parameters
    //////////////////////////////////////////////////////////////////////////

    const std::string simulationName("ActuatorLine");

    const real viscosity = 1.56e-5F;
    const real machNumber = c1o10;
    const uint timeStepAverageTimeSeriesProbe = 1;

    const real rotorDiameter = config.getValue<real>("RotorDiameter");
    const uint nodesPerDiameter = config.getValue<uint>("NodesPerDiameter");
    const real velocity = config.getValue<real>("Velocity");

    const real timeStartOut = config.getValue<real>("tStartOut");
    const real timeOut = config.getValue<real>("tOut");
    const real timeEnd = config.getValue<real>("tEnd");

    // const real timeStartAveraging = config.getValue<real>("tStartAveraging");
    const real timeStartTemporalAveraging = config.getValue<real>("tStartTmpAveraging");
    const real timeAveraging = config.getValue<real>("tAveraging");
    const real timeStartOutProbe = config.getValue<real>("tStartOutProbe");
    const real timeOutProbe = config.getValue<real>("tOutProbe");

    const real lengthX = config.getValue<real>("LengthXinDiameter") * rotorDiameter;
    const real lengthY = config.getValue<real>("LengthYinDiameter") * rotorDiameter;
    const real lengthZ = config.getValue<real>("LengthZinDiameter") * rotorDiameter;

    const std::vector<real> turbinePositionsX = config.getVector<real>("TurbinePositionsX");
    const std::vector<real> turbinePositionsY = config.getVector<real>("TurbinePositionsY");
    const std::vector<real> turbinePositionsZ = config.getVector<real>("TurbinePositionsZ");

    std::vector<real> probePositionsX {}, probePositionsY {}, probePositionsZ {};
    if (config.contains("probePositionsX")) {
        probePositionsX = config.getVector<real>("ProbePositionsX");
        probePositionsY = config.getVector<real>("ProbePositionsY");
        probePositionsZ = config.getVector<real>("ProbePositionsZ");
    }

    //////////////////////////////////////////////////////////////////////////
    // compute parameters in lattice units
    //////////////////////////////////////////////////////////////////////////

    const real deltaX = rotorDiameter / real(nodesPerDiameter);
    const real deltaT = deltaX * machNumber / (sqrt(3) * velocity);
    const real velocityLB = velocity * deltaT / deltaX;              // LB units
    const real viscosityLB = viscosity * deltaT / (deltaX * deltaX); // LB units

    const uint timeStepStartOut = timeStartOut / deltaT;
    const uint timeStepOut = timeOut / deltaT;
    const uint timeStepEnd = timeEnd / deltaT;

    const uint timeStepStartOutProbe = timeStartOutProbe / deltaT;
    const uint timeStepStartTemporalAveraging = timeStartTemporalAveraging / deltaT;
    const uint numberOfAveragingTimeSteps = timeAveraging / deltaT;
    const uint timeStepOutProbe = timeOutProbe / deltaT;

    //////////////////////////////////////////////////////////////////////////
    // create grid
    //////////////////////////////////////////////////////////////////////////

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    gridBuilder->addCoarseGrid(c0o1, -c1o2 * lengthY, -c1o2 * lengthZ, lengthX, c1o2 * lengthY, c1o2 * lengthZ, deltaX);
    gridBuilder->setPeriodicBoundaryCondition(false, false, false);
    gridBuilder->buildGrids(false);

    //////////////////////////////////////////////////////////////////////////
    // set parameters
    //////////////////////////////////////////////////////////////////////////

    auto para = std::make_shared<Parameter>(&config);

    para->worldLength = rotorDiameter;

    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(deltaX / deltaT);
    para->setViscosityRatio(deltaX * deltaX / deltaT);
    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);

    para->setInitialCondition([&](real, real, real, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx = velocityLB;
        vy = c0o1;
        vz = c0o1;
    });

    para->setTimestepStartOut(timeStepStartOut);
    para->setTimestepOut(timeStepOut);
    para->setTimestepEnd(timeStepEnd);

    para->setIsBodyForce(true);
    para->setUseStreams(true);

    //////////////////////////////////////////////////////////////////////////
    // set boundary conditions
    //////////////////////////////////////////////////////////////////////////

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::MY, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, velocityLB, c0o1, c0o1);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, c0o1);

    BoundaryConditionFactory bcFactory;
    bcFactory.setVelocityBoundaryCondition(
        BoundaryConditionFactory::VelocityBC::VelocityWithPressureInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    //////////////////////////////////////////////////////////////////////////
    // add turbine
    //////////////////////////////////////////////////////////////////////////

    const int level = 0; // grid level at which the turbine samples velocities and distributes forces
    const real smearingWidth = deltaX * std::exp2(-level) * c2o1; // width of gaussian smearing
    const real density = 1.225F;
    const uint actuatorNodesPerBlade = 32;
    const real tipSpeedRatio = 7.5F; // tipspeed ratio = angular vel * radius / inflow vel
    const std::vector<real> rotorSpeeds { c2o1 * tipSpeedRatio * velocity / rotorDiameter };

    auto actuatorFarm = std::make_shared<ActuatorFarmStandalone>(
        para, cudaMemoryManager, rotorDiameter, actuatorNodesPerBlade, turbinePositionsX, turbinePositionsY,
        turbinePositionsZ, rotorSpeeds, density, smearingWidth, level, deltaT, deltaX);
    para->addInteractor(actuatorFarm);

    actuatorFarm->enableOutput("ActuatorLineForcesAndVelocities", timeStepStartOutProbe, timeStepOutProbe);

    //////////////////////////////////////////////////////////////////////////
    // add probes
    //////////////////////////////////////////////////////////////////////////

    std::vector<real> planePositions = { -c1o1 * rotorDiameter, c1o1 * rotorDiameter, c3o1 * rotorDiameter };

    for (size_t i = 0; i < planePositions.size(); i++) {
        const std::string name = "planeProbe_" + std::to_string(i);
        auto planeProbe =
            std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), name, timeStepStartTemporalAveraging,
                                    numberOfAveragingTimeSteps, timeStepStartOutProbe, timeStepOutProbe, false, false);
        planeProbe->addProbePlane(turbinePositionsX[0] + planePositions[i], -c1o2 * lengthY, -c1o2 * lengthZ, deltaX,
                                  lengthY, lengthZ);
        planeProbe->addStatistic(Probe::Statistic::Means);
        planeProbe->addStatistic(Probe::Statistic::Variances);
        planeProbe->addStatistic(Probe::Statistic::Instantaneous);
        para->addSampler(planeProbe);
    }

    auto planeProbeVertical = std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "planeProbeVertical",
                                                      timeStepStartTemporalAveraging, numberOfAveragingTimeSteps,
                                                      timeStepStartOutProbe, timeStepOutProbe, false, false);
    planeProbeVertical->addProbePlane(c0o1, turbinePositionsY[0], -c1o2 * lengthZ, lengthX, deltaX, lengthZ);
    planeProbeVertical->addStatistic(Probe::Statistic::Means);
    planeProbeVertical->addStatistic(Probe::Statistic::Variances);
    planeProbeVertical->addStatistic(Probe::Statistic::Instantaneous);
    para->addSampler(planeProbeVertical);

    auto planeProbeHorizontal = std::make_shared<Probe>(
        para, cudaMemoryManager, para->getOutputPath(), "planeProbeHorizontal", timeStepStartTemporalAveraging,
        numberOfAveragingTimeSteps, timeStepStartOutProbe, timeStepOutProbe, false, false);
    planeProbeHorizontal->addProbePlane(c0o1, -c1o2 * lengthY, turbinePositionsZ[0], lengthX, lengthY, deltaX);
    planeProbeHorizontal->addStatistic(Probe::Statistic::Means);
    planeProbeHorizontal->addStatistic(Probe::Statistic::Variances);
    planeProbeHorizontal->addStatistic(Probe::Statistic::Instantaneous);
    para->addSampler(planeProbeHorizontal);

    if (probePositionsX.size() > 0) {
        auto timeseriesProbe = std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "timeProbe",
                                                       timeStepStartTemporalAveraging, timeStepAverageTimeSeriesProbe,
                                                       timeStepStartOutProbe, timeStepOutProbe, true, false);
        timeseriesProbe->addProbePointsFromList(probePositionsX, probePositionsY, probePositionsZ);
        timeseriesProbe->addStatistic(Probe::Statistic::Instantaneous);
        timeseriesProbe->addStatistic(Probe::Statistic::Means);
        para->addSampler(timeseriesProbe);
    }

    //////////////////////////////////////////////////////////////////////////
    // run simulation
    //////////////////////////////////////////////////////////////////////////

    VF_LOG_INFO("Start Running ActuatorLine Showcase...\n");

    VF_LOG_INFO("turbine parameters:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("rotorDiameter [m]      = {}", rotorDiameter);
    VF_LOG_INFO("nodesPerDiameter       = {}", nodesPerDiameter);
    VF_LOG_INFO("actuatorNodesPerBlade  = {}", actuatorNodesPerBlade);
    VF_LOG_INFO("smearingWidth [m]      = {}", smearingWidth);
    VF_LOG_INFO("tipSpeedRatio          = {}", tipSpeedRatio);

    Simulation simulation(para, cudaMemoryManager, gridBuilder, &bcFactory, tmFactory);
    simulation.run();
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        const auto config = vf::basics::loadConfig(argc, argv, defaultConfigFile);
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}

//! \}
