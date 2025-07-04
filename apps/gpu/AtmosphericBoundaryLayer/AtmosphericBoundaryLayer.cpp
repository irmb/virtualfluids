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
//! \addtogroup AtmosphericBoundaryLayer
//! \ingroup gpu_apps
//! \{
//! \author Henry Korb, Henrik Asmuth
//=======================================================================================
#include <memory>
#include <numeric>

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>
#include <basics/geometry3d/Axis.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/TransientBCSetter/TransientBCSetter.h"
#include "GridGenerator/geometries/BoundingBox/BoundingBox.h"
#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/MultipleGridBuilderFacade.h"
#include "GridGenerator/utilities/communication.h"

//////////////////////////////////////////////////////////////////////////

#include "PreCollisionInteractor/CoriolisForce.h"
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Samplers/PlanarAverageProbe.h"
#include "gpu/core/Samplers/PrecursorWriter.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/Samplers/WallModelProbe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"
#include "grid/Grid.h"
#include "grid/GridDimensions.h"

using namespace vf::basics::constant;

const std::string defaultConfigFile = "abl.cfg";

void run(const vf::basics::ConfigurationFile& config)
{
    const std::string simulationName("AtmosphericBoundaryLayer");

    const uint samplingOffsetWallModel = 2;
    const uint averagingTimestepsPlaneProbes = 10;
    const uint maximumNumberOfTimestepsPerPrecursorFile = 1000;
    const real viscosity = 1.56e-5;
    const real machNumber = 0.1;

    const real boundaryLayerHeight = config.getValue("BoundaryLayerHeight", 1000.0);

    const real lengthX = 6 * boundaryLayerHeight;
    const real lengthY = 4 * boundaryLayerHeight;
    const real lengthZ = boundaryLayerHeight;

    const real roughnessLength = config.getValue("RoughnessLength", c1o10);
    real frictionVelocity = config.getValue("FrictionVelocity", c4o10);
    const real vonKarmanConstant = config.getValue("vonKarmanConstant", c4o10);

    const uint nodesPerBoundaryLyerHeight = config.getValue<uint>("NodesPerBoundaryLayerHeight", 64);

    const float periodicShift = config.getValue<float>("PeriodicShift", c0o1);

    const bool writePrecursor = config.getValue("WritePrecursor", false);
    const bool useCoriolisForce = config.getValue("UseCoriolisForce", false);
    const real geostrophicWindSpeed = config.getValue("GeostrophicWindSpeed", c8o1);
    const real geostrophicWindDirection = config.getValue("GeostrophicWindDirection", c0o1);
    const real coriolisParameter = config.getValue("CoriolisParameter", 1e-4F);

    const bool useDistributionsForPrecursor = config.getValue<bool>("UseDistributions", false);
    std::string precursorDirectory = config.getValue<std::string>("PrecursorDirectory", "precursor/");
    if (precursorDirectory.back() != '/')
        precursorDirectory += '/';
    const int timeStepsWritePrecursor = config.getValue<int>("nTimestepsWritePrecursor", 10);
    const real timeStartPrecursor = config.getValue<real>("tStartPrecursor", 36000.);
    const real positionXPrecursorSamplingPlane = config.getValue<real>("posXPrecursor", c1o2 * lengthX);

    const bool usePrecursorInflow = config.getValue("ReadPrecursor", false);
    const int timestepsBetweenReadsPrecursor = config.getValue<int>("nTimestepsReadPrecursor", 10);

    const bool useRefinement = config.getValue<bool>("Refinement", false);

    // all in s
    const float timeStartOut = config.getValue<real>("tStartOut");
    const float timeOut = config.getValue<real>("tOut");
    const float timeEnd = config.getValue<real>("tEnd");

    const float timeStartAveraging = config.getValue<real>("tStartAveraging");
    const float timeStartTemporalAveraging = config.getValue<real>("tStartTmpAveraging");
    const float timeAveraging = config.getValue<real>("tAveraging");
    const float timeStartOutProbe = config.getValue<real>("tStartOutProbe");
    const float timeOutProbe = config.getValue<real>("tOutProbe");

    //////////////////////////////////////////////////////////////////////////
    // compute parameters in lattice units
    //////////////////////////////////////////////////////////////////////////
    if (useCoriolisForce)
        frictionVelocity = geostrophicWindSpeed * vonKarmanConstant / std::log(boundaryLayerHeight / roughnessLength + c1o1);
    const auto velocityProfile = [&](real coordZ) {
        return frictionVelocity / vonKarmanConstant * std::log(coordZ / roughnessLength + c1o1);
    };
    const real velocity = c1o2 * velocityProfile(boundaryLayerHeight);

    const real deltaX = boundaryLayerHeight / real(nodesPerBoundaryLyerHeight);

    const real deltaT = c1oSqrt3 * deltaX * machNumber / velocity;

    const real velocityLB = velocity * deltaT / deltaX;

    const real viscosityLB = viscosity * deltaT / (deltaX * deltaX);

    const real pressureGradient = frictionVelocity * frictionVelocity / boundaryLayerHeight;
    const real pressureGradientLB = pressureGradient * (deltaT * deltaT) / deltaX;

    const uint timeStepStartAveraging = uint(timeStartAveraging / deltaT);
    const uint timeStepStartTemporalAveraging = uint(timeStartTemporalAveraging / deltaT);
    const uint timeStepAveraging = uint(timeAveraging / deltaT);
    const uint timeStepStartOutProbe = uint(timeStartOutProbe / deltaT);
    const uint timeStepOutProbe = uint(timeOutProbe / deltaT);

    const uint timeStepStartPrecursor = uint(timeStartPrecursor / deltaT);

    VF_LOG_INFO("Friction velocity [m/s] {}", frictionVelocity);

    //////////////////////////////////////////////////////////////////////////
    // create grid
    //////////////////////////////////////////////////////////////////////////

    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

    const int numberOfProcesses = communicator.getNumberOfProcesses();
    const int processID = communicator.getProcessID();

    const real overlap = c8o1 * deltaX;

    auto dimensions = std::make_shared<GridDimensions>(c0o1, lengthX, c0o1, lengthY, c0o1, lengthZ, deltaX);

    auto gridBuilder = std::make_unique<MultipleGridBuilderFacade>(dimensions, overlap);
    auto scalingFactory = GridScalingFactory();

    for (int iProcess = 1; iProcess < numberOfProcesses; iProcess++)
        gridBuilder->addDomainSplit(lengthX / numberOfProcesses * iProcess, Axis::x);

    if (useRefinement) {
        gridBuilder->setNumberOfLayersForRefinement(4, 0);
        const real endRefinement = usePrecursorInflow ? lengthX - boundaryLayerHeight : lengthX;
        gridBuilder->addFineGrid(std::make_shared<Cuboid>(c0o1, c0o1, c0o1, endRefinement, lengthY, c1o2 * lengthZ), 1);
        scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);
    }

    gridBuilder->setPeriodicBoundaryCondition(!usePrecursorInflow, true, false);

    if (!usePrecursorInflow)
        gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(periodicShift);

    gridBuilder->createGrids(processID);

    //////////////////////////////////////////////////////////////////////////
    // set parameters
    //////////////////////////////////////////////////////////////////////////
    auto para = std::make_shared<Parameter>(numberOfProcesses, processID, &config);

    para->setOutputPrefix(simulationName);
    para->setPrintFiles(true);

    if (!usePrecursorInflow && !useCoriolisForce)
        para->setForcing(pressureGradientLB, 0, 0);
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(deltaX / deltaT);
    para->setViscosityRatio(deltaX * deltaX / deltaT);
    para->setDensityRatio(c1o1);

    para->setUseStreams(numberOfProcesses > 1);
    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);

    para->setTimestepStartOut(uint(timeStartOut / deltaT));
    para->setTimestepOut(uint(timeOut / deltaT));
    para->setTimestepEnd(uint(timeEnd / deltaT));

    std::vector<uint> devices(10);
    std::iota(devices.begin(), devices.end(), 0);
    para->setDevices(devices);
    para->setMaxDev(numberOfProcesses);
    if (usePrecursorInflow) {
        para->setInitialCondition([&](real, real, real coordZ, real& rho, real& vx, real& vy, real& vz) {
            rho = c0o1;
            vx = velocityProfile(coordZ) * deltaT / deltaX;
            vy = c0o1;
            vz = c0o1;
        });
    } else {
        para->setInitialCondition([&](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) {
            const real relativeX = coordX / lengthX;
            const real relativeY = coordY / lengthY;
            const real relativeZ = coordZ / lengthZ;

            const real horizontalPerturbation = std::sin(cPi * c16o1 * relativeX);
            const real verticalPerturbation = std::sin(cPi * c8o1 * relativeZ) / (std::pow(relativeZ, c2o1) + c1o1);
            const real perturbation = c2o1 * horizontalPerturbation * verticalPerturbation;

            rho = c0o1;
            vx = (velocityProfile(coordZ) + perturbation) * (c1o1 - c1o10 * std::abs(relativeY - c1o2)) * deltaT / deltaX;
            vy = perturbation * deltaT / deltaX;
            vz = c8o1 * frictionVelocity / vonKarmanConstant *
                 (std::sin(cPi * c8o1 * coordY / boundaryLayerHeight) * std::sin(cPi * c8o1 * relativeZ) +
                  std::sin(cPi * c8o1 * relativeX)) /
                 (std::pow(c1o2 * lengthZ - coordZ, c2o1) + c1o1) * deltaT / deltaX;
        });
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (usePrecursorInflow) {
        auto precursor = createFileCollection(precursorDirectory + "precursor", TransientBCFileType::VTK);
        gridBuilder->setPrecursorBoundaryCondition(SideType::MX, precursor, timestepsBetweenReadsPrecursor);
        gridBuilder->setPressureBoundaryCondition(SideType::PX, c0o1);
    }

    gridBuilder->setStressBoundaryCondition(SideType::MZ, c0o1, c0o1, c1o1, samplingOffsetWallModel, roughnessLength,
                                            deltaX);

    gridBuilder->setSlipBoundaryCondition(SideType::PZ, c0o1, c0o1, -c1o1);

    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackPressureCompressible);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipTurbulentViscosityCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    if (useDistributionsForPrecursor) {
        bcFactory.setPrecursorBoundaryCondition(BoundaryConditionFactory::PrecursorBC::PrecursorDistributions);
    } else {
        bcFactory.setPrecursorBoundaryCondition(BoundaryConditionFactory::PrecursorBC::PrecursorNonReflectiveCompressible);
    }

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    if (useCoriolisForce) {
        para->setIsBodyForce(true);
        para->setAllNodesAllFeatures(true);
        auto coriolisForce = std::make_shared<CoriolisForce>(para, cudaMemoryManager, geostrophicWindSpeed*std::cos(geostrophicWindDirection),
                                                             geostrophicWindDirection*std::sin(geostrophicWindDirection), coriolisParameter);
        para->addInteractor(coriolisForce);
    }

    //////////////////////////////////////////////////////////////////////////
    // add probes
    //////////////////////////////////////////////////////////////////////////

    if (!usePrecursorInflow && processID == 0) {
        const auto planarAverageProbe =
            std::make_shared<PlanarAverageProbe>(para, cudaMemoryManager, para->getOutputPath(), "planarAverageProbe",
                                                 timeStepStartAveraging, timeStepStartTemporalAveraging, timeStepAveraging,
                                                 timeStepStartOutProbe, timeStepOutProbe, Axis::z, true, false);
        planarAverageProbe->addAllAvailableStatistics();
        planarAverageProbe->setFileNameToNOut();
        para->addSampler(planarAverageProbe);

        const auto wallModelProbe = std::make_shared<WallModelProbe>(
            para, cudaMemoryManager, para->getOutputPath(), "wallModelProbe", timeStepStartAveraging,
            timeStepStartTemporalAveraging, timeStepAveraging / 4, timeStepStartOutProbe, timeStepOutProbe, false, true,
            true, para->getIsBodyForce());

        para->addSampler(wallModelProbe);

        para->setHasWallModelMonitor(true);
    }

    for (int iPlane = 0; iPlane < 3; iPlane++) {
        const std::string name = "planeProbe" + std::to_string(iPlane);
        const auto horizontalProbe =
            std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), name, timeStepStartAveraging,
                                    averagingTimestepsPlaneProbes, timeStepStartOutProbe, timeStepOutProbe, false, false);
        horizontalProbe->addProbePlane(c0o1, c0o1, iPlane * lengthZ / c4o1, lengthX, lengthY, deltaX);
        horizontalProbe->addAllAvailableStatistics();
        para->addSampler(horizontalProbe);
    }

    auto crossStreamPlane =
        std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "crossStreamPlane", timeStepStartAveraging,
                                averagingTimestepsPlaneProbes, timeStepStartOutProbe, timeStepOutProbe, false, false);
    crossStreamPlane->addProbePlane(c1o2 * lengthX, c0o1, c0o1, deltaX, lengthY, lengthZ);
    crossStreamPlane->addAllAvailableStatistics();
    para->addSampler(crossStreamPlane);

    if (usePrecursorInflow) {
        auto streamwisePlane = std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "streamwisePlane",
                                                       timeStepStartAveraging, averagingTimestepsPlaneProbes,
                                                       timeStepStartOutProbe, timeStepOutProbe, false, false);
        streamwisePlane->addProbePlane(c0o1, c1o2 * lengthY, c0o1, lengthX, deltaX, lengthZ);
        streamwisePlane->addAllAvailableStatistics();
        para->addSampler(streamwisePlane);
    }

    if (writePrecursor) {
        const std::string fullPrecursorDirectory = para->getOutputPath() + precursorDirectory;
        const auto outputVariable = useDistributionsForPrecursor ? PrecursorWriter::OutputVariable::Distributions
                                                                 : PrecursorWriter::OutputVariable::Velocities;
        auto precursorWriter = std::make_shared<PrecursorWriter>(
            para, cudaMemoryManager, fullPrecursorDirectory, "precursor", positionXPrecursorSamplingPlane, c0o1, lengthY,
            c0o1, lengthZ, timeStepStartPrecursor, timeStepsWritePrecursor, outputVariable,
            maximumNumberOfTimestepsPerPrecursorFile);
        para->addSampler(precursorWriter);
    }

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    VF_LOG_INFO("process parameter:");
    VF_LOG_INFO("Number of Processes {} process ID {}", numberOfProcesses, processID);
    printf("\n");
    
    para->worldLength = lengthX + c2o1*deltaX;
    Simulation simulation(para, cudaMemoryManager, gridBuilder->getGridBuilder(), &bcFactory, tmFactory, &scalingFactory);
    simulation.run();
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        auto config = vf::basics::loadConfig(argc, argv, defaultConfigFile);
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}

//! \}
