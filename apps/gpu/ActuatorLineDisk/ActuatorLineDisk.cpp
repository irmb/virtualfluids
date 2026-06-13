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
//! \addtogroup ActuatorLineDisk
//! \ingroup gpu_apps
//! \{
//! \author Nils Horneff
//=======================================================================================

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>

#include <logger/Logger.h>

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PreCollisionInteractor/Actuator/ActuatorFarm.h"
#include "gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandalone.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"

namespace
{
const std::string defaultConfigFile = "actuatorlinedisk.cfg";
using namespace vf::basics::constant;
using namespace vf::gpu;

std::vector<real> parseNumberList(std::string text)
{
    for (char& ch : text) {
        if (ch == ',')
            ch = ' ';
    }

    std::stringstream stream(text);
    std::vector<real> values;
    real value = c0o1;
    while (stream >> value)
        values.push_back(value);

    return values;
}

std::vector<std::vector<real>> parseNestedNumberLists(const std::string& text)
{
    std::vector<std::vector<real>> lists;

    size_t start = text.find('[');
    while (start != std::string::npos) {
        const size_t end = text.find(']', start);
        lists.push_back(parseNumberList(text.substr(start + 1, end - start - 1)));
        start = text.find('[', end);
    }

    return lists;
}

void scaleNestedNumberLists(std::vector<std::vector<real>>& lists, const real factor)
{
    for (auto& list : lists) {
        for (real& value : list)
            value *= factor;
    }
}

std::vector<real> computeBladeRadii(const real diameter, const uint numberOfPointsPerBlade)
{
    std::vector<real> bladeRadii(numberOfPointsPerBlade);
    const real spacing = c1o2 * diameter / static_cast<real>(numberOfPointsPerBlade);

    for (uint i = 0; i < numberOfPointsPerBlade; ++i)
        bladeRadii[i] = spacing * (static_cast<real>(i) + c1o2);

    return bladeRadii;
}

real computeEquivalentDiskNormalCoefficient(const real thrustCoefficient)
{
    const real induction = c1o2 * (c1o1 - std::sqrt(c1o1 - thrustCoefficient));
    return c4o1 * induction / (c1o1 - induction);
}

real getRealConfigValueWithFallback(
    const vf::basics::ConfigurationFile& config,
    const std::string& key,
    const std::string& fallbackKey)
{
    if (config.contains(key))
        return config.getValue<real>(key);
    return config.getValue<real>(fallbackKey);
}

BoundaryConditionFactory::VelocityBC getVelocityBoundaryCondition(
    const vf::basics::ConfigurationFile& config)
{
    if (!config.contains("VELOCITY_BOUNDARY_CONDITION"))
        return BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible;

    const auto velocityBoundaryCondition = config.getValue<std::string>("VELOCITY_BOUNDARY_CONDITION");
    if (velocityBoundaryCondition == "VelocityBounceBack")
        return BoundaryConditionFactory::VelocityBC::VelocityBounceBack;
    if (velocityBoundaryCondition == "VelocityInterpolatedIncompressible")
        return BoundaryConditionFactory::VelocityBC::VelocityInterpolatedIncompressible;
    if (velocityBoundaryCondition == "VelocityInterpolatedCompressible")
        return BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible;
    if (velocityBoundaryCondition == "VelocityWithPressureInterpolatedCompressible")
        return BoundaryConditionFactory::VelocityBC::VelocityWithPressureInterpolatedCompressible;

    throw std::runtime_error(
        "VELOCITY_BOUNDARY_CONDITION must be one of: VelocityBounceBack, "
        "VelocityInterpolatedIncompressible, VelocityInterpolatedCompressible, "
        "VelocityWithPressureInterpolatedCompressible.");
}

real computeDiskHubRadius(const real diameter, const uint numberOfPointsPerBlade)
{
    return c1o4 * diameter / static_cast<real>(numberOfPointsPerBlade);
}

std::vector<real> computeDiskBladeAnnulusAreas(
    const real diameter,
    const uint numberOfBlades,
    const uint numberOfPointsPerBlade)
{
    const auto bladeRadii = computeBladeRadii(diameter, numberOfPointsPerBlade);
    std::vector<real> annulusAreas(numberOfPointsPerBlade);

    const real pi = std::acos(-c1o1);
    for (uint i = 0; i < numberOfPointsPerBlade; ++i) {
        const real innerRadius = (i == 0)
                                     ? computeDiskHubRadius(diameter, numberOfPointsPerBlade)
                                     : c1o2 * (bladeRadii[i] + bladeRadii[i - 1]);
        const real outerRadius = (i + 1 == numberOfPointsPerBlade)
                                     ? c1o2 * diameter
                                     : c1o2 * (bladeRadii[i] + bladeRadii[i + 1]);
        annulusAreas[i] = pi * (outerRadius * outerRadius - innerRadius * innerRadius) /
                          static_cast<real>(numberOfBlades);
    }

    return annulusAreas;
}

std::vector<real> computeDiskBladeNormalCoefficients(
    const real diameter,
    const uint numberOfBlades,
    const uint numberOfPointsPerBlade,
    const real thrustCoefficient)
{
    const auto bladeRadii = computeBladeRadii(diameter, numberOfPointsPerBlade);
    const auto annulusAreas = computeDiskBladeAnnulusAreas(diameter, numberOfBlades, numberOfPointsPerBlade);
    const real diskNormalCoefficient = computeEquivalentDiskNormalCoefficient(thrustCoefficient);

    std::vector<real> bladeNormalCoefficients(numberOfPointsPerBlade);
    for (uint i = 0; i < numberOfPointsPerBlade; ++i) {
        const real previousRadius = (i == 0) ? c0o1 : bladeRadii[i - 1];
        const real nextRadius = (i + 1 == numberOfPointsPerBlade) ? c1o2 * diameter : bladeRadii[i + 1];
        const real bladePointWidth = c1o2 * (nextRadius - previousRadius);
        const real chordArgument = c1o1 - std::pow(c4o1 * bladeRadii[i] / diameter - c1o1, c2o1);
        const real bladeChord = c2o1 * std::sqrt(std::max(chordArgument, c0o1));
        bladeNormalCoefficients[i] = diskNormalCoefficient * annulusAreas[i] / (bladeChord * bladePointWidth);
    }

    return bladeNormalCoefficients;
}

struct DiskActuatorConfig
{
    ActuatorFarm::HubConfig hubConfig;
    real hubDragCoefficient;
    real hubSkinFrictionCoefficient;
    std::vector<real> bladeNormalCoefficients;
};

DiskActuatorConfig buildDiskActuatorConfig(
    const real diameter,
    const uint numberOfBlades,
    const uint numberOfPointsPerBlade,
    const real thrustCoefficient)
{
    const real diskNormalCoefficient = computeEquivalentDiskNormalCoefficient(thrustCoefficient);
    return {
        ActuatorFarm::HubConfig { c0o1, computeDiskHubRadius(diameter, numberOfPointsPerBlade), c0o1, 1 },
        diskNormalCoefficient,
        c0o1,
        computeDiskBladeNormalCoefficients(diameter, numberOfBlades, numberOfPointsPerBlade, thrustCoefficient)
    };
}

uint toPositiveTimestep(const real flowThroughCount, const real timeFlowThroughDomain, const real deltaT)
{
    return std::max(static_cast<uint>(1), static_cast<uint>(timeFlowThroughDomain * flowThroughCount / deltaT));
}

void run(const vf::basics::ConfigurationFile& config)
{
    const real viscosity = config.getValue<real>("VISCOSITY_KINEMATIC");
    const uint numberOfBlades = config.getValue<uint>("NUMBER_BLADES");
    const real diskDiameter = config.getValue<real>("DISK_DIAMETER");
    const real velocityInlet = config.getValue<real>("VELOCITY_INLET");
    const real thrustCoefficient = config.getValue<real>("THRUST_COEFFICIENT");
    const uint numberOfCellsPerDiameter = config.getValue<uint>("NUMBER_CELLS_PER_DIAMETER");
    const uint numberOfPointsPerBlade = config.getValue<uint>("NUMBER_ALM_NODES_PER_SPAN");
    const real smearingWidthPerDx = config.getValue<real>("SMEARING_WIDTH_PER_DX");
    auto positionDomain = parseNestedNumberLists(config.getValue<std::string>("POSITION_DOMAIN"));
    auto positionDiskCenter = parseNestedNumberLists(config.getValue<std::string>("POSITION_DISK_CENTER"));
    auto positionProbePlaneZ = parseNestedNumberLists(config.getValue<std::string>("POSITION_PROBE_PLANE_Z"));
    const real quadricDelimiter = config.getValue<real>("QUADRIC_DELIMITER");
    const real machNumber = config.getValue<real>("LBM_MACH_NUMBER");
    const bool printFiles = config.getValue<bool>("FLAG_PRINT_FILES");
    const real nStartWriteLoads = config.getValue<real>("N_START_WRITE_LOADS");
    const real nStartWriteProbe = config.getValue<real>("N_START_WRITE_PROBE");
    const real nStartWriteField = config.getValue<real>("N_START_WRITE_FIELD");
    const real nStartProbeAveraging = config.getValue<real>("N_START_PROBE_AVERAGING");
    const real nIntervallAveraging = config.getValue<real>("N_INTERVALL_AVERAGING");
    const real nIntervallWriteLoad = config.getValue<real>("N_INTERVALL_WRITE_LOAD");
    const real nIntervallWriteProbe = config.getValue<real>("N_INTERVALL_WRITE_PROBE");
    const real nIntervallWriteField =
        getRealConfigValueWithFallback(config, "N_INTERVALL_WRITE_FIELD", "N_INTERVALL_WRTIE_FIELD");
    const real nEnd = config.getValue<real>("N_END");

    scaleNestedNumberLists(positionDomain, diskDiameter);
    scaleNestedNumberLists(positionDiskCenter, diskDiameter);
    scaleNestedNumberLists(positionProbePlaneZ, diskDiameter);

    if (positionDomain.size() != 3 || positionDomain[0].size() != 2 || positionDomain[1].size() != 2 ||
        positionDomain[2].size() != 2) {
        throw std::runtime_error("POSITION_DOMAIN must be [x_min,x_max],[y_min,y_max],[z_min,z_max].");
    }
    if (positionDiskCenter.size() != 3 || positionDiskCenter[0].size() != 1 || positionDiskCenter[1].size() != 1 ||
        positionDiskCenter[2].size() != 1) {
        throw std::runtime_error("POSITION_DISK_CENTER must be [x],[y],[z].");
    }

    const real xcm = positionDomain[0][0];
    const real xcp = positionDomain[0][1];
    const real ycm = positionDomain[1][0];
    const real ycp = positionDomain[1][1];
    const real zcm = positionDomain[2][0];
    const real zcp = positionDomain[2][1];

    const real diskCenterX = positionDiskCenter[0][0];
    const real diskCenterY = positionDiskCenter[1][0];
    const real diskCenterZ = positionDiskCenter[2][0];

    const real timeFlowThroughDomain = (xcp - xcm) / velocityInlet;
    const real deltaX = diskDiameter / static_cast<real>(numberOfCellsPerDiameter);
    const real deltaT = machNumber * deltaX / (std::sqrt(c3o1) * velocityInlet);
    const real velocityRatio = deltaX / deltaT;
    const real velocityLB = velocityInlet / velocityRatio;
    const real viscosityLB = viscosity / (velocityRatio * deltaX);
    const real smearingWidth = smearingWidthPerDx * deltaX;
    const auto velocityBoundaryCondition = getVelocityBoundaryCondition(config);

    const uint timeStepStartWriteLoads = toPositiveTimestep(nStartWriteLoads, timeFlowThroughDomain, deltaT);
    const uint timeStepStartWriteProbe = toPositiveTimestep(nStartWriteProbe, timeFlowThroughDomain, deltaT);
    const uint timeStepStartWriteField = toPositiveTimestep(nStartWriteField, timeFlowThroughDomain, deltaT);
    const uint timeStepStartProbeAveraging = toPositiveTimestep(nStartProbeAveraging, timeFlowThroughDomain, deltaT);
    const uint timeStepIntervalAveraging = toPositiveTimestep(nIntervallAveraging, timeFlowThroughDomain, deltaT);
    const uint timeStepIntervalWriteLoads = toPositiveTimestep(nIntervallWriteLoad, timeFlowThroughDomain, deltaT);
    const uint timeStepIntervalWriteProbe = toPositiveTimestep(nIntervallWriteProbe, timeFlowThroughDomain, deltaT);
    const uint timeStepIntervalWriteField = toPositiveTimestep(nIntervallWriteField, timeFlowThroughDomain, deltaT);
    const uint timeStepEnd = static_cast<uint>(timeFlowThroughDomain * nEnd / deltaT);

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();
    gridBuilder->addCoarseGrid(xcm, ycm, zcm, xcp, ycp, zcp, deltaX);
    gridBuilder->setPeriodicBoundaryCondition(false, false, false);
    gridBuilder->buildGrids(false);

    auto para = std::make_shared<Parameter>(&config);
    para->worldLength = xcp - xcm;
    para->setOutputPath(config.getValue<std::string>("Path"));
    para->setOutputPrefix("fields/disk");
    para->setPrintFiles(printFiles);
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(velocityRatio);
    para->setViscosityRatio(deltaX * deltaX / deltaT);
    para->setDensityRatio(c1o1);
    para->setOutflowPressureCorrectionFactor(c0o1);
    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);
    para->setInitialCondition([velocityLB](real, real, real, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx = velocityLB;
        vy = c0o1;
        vz = c0o1;
    });
    para->setTimestepStartOut(timeStepStartWriteField);
    para->setTimestepOut(timeStepIntervalWriteField);
    para->setTimestepEnd(timeStepEnd);
    para->setIsBodyForce(true);
    para->setQuadricLimiters(quadricDelimiter, quadricDelimiter, quadricDelimiter);

    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, c0o1, c0o1);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, c1o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::MY, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::PY, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::MZ, velocityLB, c0o1, c0o1);
    gridBuilder->setVelocityBoundaryCondition(SideType::PZ, velocityLB, c0o1, c0o1);

    BoundaryConditionFactory bcFactory;
    bcFactory.setVelocityBoundaryCondition(velocityBoundaryCondition);
    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipTurbulentViscosityCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    GridScalingFactory gridScalingFactory;
    gridScalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    const int gridLevelForAlm = 0;
    const auto diskActuatorConfig =
        buildDiskActuatorConfig(diskDiameter, numberOfBlades, numberOfPointsPerBlade, thrustCoefficient);
    const std::vector<real> turbinePosX { diskCenterX };
    const std::vector<real> turbinePosY { diskCenterY };
    const std::vector<real> turbinePosZ { diskCenterZ };
    const std::vector<real> rotorSpeeds { c0o1 };
    auto actuatorFarm = std::make_shared<ActuatorFarmStandalone>(
        para,
        cudaMemoryManager,
        diskDiameter,
        numberOfPointsPerBlade,
        turbinePosX,
        turbinePosY,
        turbinePosZ,
        rotorSpeeds,
        smearingWidth,
        gridLevelForAlm,
        diskActuatorConfig.hubConfig,
        std::nullopt,
        diskActuatorConfig.hubDragCoefficient,
        diskActuatorConfig.hubSkinFrictionCoefficient,
        std::nullopt,
        numberOfBlades,
        diskActuatorConfig.bladeNormalCoefficients);
    actuatorFarm->enableOutput("loads/ALM", timeStepStartWriteLoads, timeStepIntervalWriteLoads);
    para->addInteractor(actuatorFarm);

    if (positionProbePlaneZ.size() != 3 || positionProbePlaneZ[0].size() != 2 || positionProbePlaneZ[1].size() != 2 ||
        positionProbePlaneZ[2].empty()) {
        throw std::runtime_error("POSITION_PROBE_PLANE_Z must be [x_range],[y_range],[z_positions].");
    }

    const auto& xRange = positionProbePlaneZ[0];
    const auto& yRange = positionProbePlaneZ[1];
    const auto& zPositions = positionProbePlaneZ[2];

    const real x0 = std::min(xRange[0], xRange[1]);
    const real x1 = std::max(xRange[0], xRange[1]);
    const real y0 = std::min(yRange[0], yRange[1]);
    const real y1 = std::max(yRange[0], yRange[1]);
    const real requestedProbePlaneZ = zPositions[0];
    const real level0MinZ = zcm - c1o2 * deltaX;
    const real lowerProbePlaneZ =
        level0MinZ + std::floor((requestedProbePlaneZ - level0MinZ) / deltaX) * deltaX;
    const real upperProbePlaneZ = lowerProbePlaneZ + deltaX;
    const real probePlaneZ = (requestedProbePlaneZ - lowerProbePlaneZ <= upperProbePlaneZ - requestedProbePlaneZ)
                                 ? lowerProbePlaneZ
                                 : upperProbePlaneZ;

    auto probe = std::make_shared<Probe>(
        para,
        cudaMemoryManager,
        para->getOutputPath(),
        "zplane/zPlane",
        timeStepStartProbeAveraging,
        timeStepIntervalAveraging,
        timeStepStartWriteProbe,
        timeStepIntervalWriteProbe,
        false,
        false);
    probe->addProbePlane(x0, y0, probePlaneZ, x1 - x0, y1 - y0, deltaX);
    probe->addAllAvailableStatistics();
    para->addSampler(probe);

    VF_LOG_INFO("Start running ActuatorLineDisk...");
    VF_LOG_INFO("thrust coefficient = {}", thrustCoefficient);
    VF_LOG_INFO("time for one flow-through = {}", timeFlowThroughDomain);
    VF_LOG_INFO("dx = {}", deltaX);
    VF_LOG_INFO("dt = {}", deltaT);
    VF_LOG_INFO("lbm velocity = {}", velocityLB);
    VF_LOG_INFO("smearing width = {}", smearingWidth);
    VF_LOG_INFO("timestep end = {}", timeStepEnd);
    if (!config.contains("VELOCITY_BOUNDARY_CONDITION"))
        VF_LOG_INFO("Using default inlet velocity BC: VelocityInterpolatedCompressible");

    Simulation simulation(para, cudaMemoryManager, gridBuilder, &bcFactory, tmFactory, &gridScalingFactory);
    simulation.run();
}
} // namespace

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
