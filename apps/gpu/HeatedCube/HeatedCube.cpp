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
//! \author Henrik Asmuth, Henry Korb
//=======================================================================================
#define _USE_MATH_DEFINES
#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const std::string simulationName("HeatedCube");
const std::string defaultConfigFile("heated_cube.cfg");

using namespace vf::basics::constant;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run(const vf::basics::ConfigurationFile& config)
{
    const real meanTemperature = c0o1;
    const real temperatureDifference = c1o1;
    const real temperatureColdSide = meanTemperature - c1o2 * temperatureDifference;
    const real temperatureHotSide = meanTemperature + c1o2 * temperatureDifference;

    const real sideLength = config.getValue("domainLength", c1o1); // Side length of the cube
    const uint numberOfNodes = config.getValue<uint>("nNodes", 16);

    const real viscosity = config.getValue("viscosity", c1o10);
    const real gravity = config.getValue("gravity", 9.81f);

    const real machNumber = config.getValue<real>("Ma", c1o10);
    const real prandtlNumber = config.getValue<real>("Pr", 0.7f);
    const real rayleighNumber = config.getValue<real>("Ra", 1000.f);

    const real diffusivity = viscosity / prandtlNumber;
    const real thermalExpansion =
        (rayleighNumber * diffusivity * viscosity) / (gravity * std::pow(sideLength, c3o1) * temperatureDifference);
    // const real Thot = (Ra*diffusivity*viscosity)/(gravity*thermalExpansion*pow(L,c3o1))+Tcold;

    const real velocity = std::sqrt(gravity * thermalExpansion * temperatureDifference * sideLength);

    // all in s
    const real tStartOut = config.getValue<real>("tStartOut");
    const real tOut = config.getValue<real>("tOut");
    const real tEnd = config.getValue<real>("tEnd"); // total time of simulation

    const real deltaX = sideLength / real(numberOfNodes);
    const real deltaT = deltaX * machNumber / (std::sqrt(3) * velocity);
    const real velocityLB = velocity * deltaT / deltaX;
    const real viscosityLB = viscosity * deltaT / (deltaX * deltaX);
    const real diffusivityLB = viscosityLB / prandtlNumber;
    const real thermalDiffusionVelocity = diffusivity / sideLength;

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("dt   = {}", deltaT);
    VF_LOG_INFO("dx   = {}", deltaX);
    VF_LOG_INFO("viscosity [dx^2/dt] = {}", viscosityLB);
    VF_LOG_INFO("Ma = {}", machNumber);
    VF_LOG_INFO("Pr = {}", prandtlNumber);
    VF_LOG_INFO("Ra = {}", rayleighNumber);
    VF_LOG_INFO("T_cold = {}", temperatureColdSide);
    VF_LOG_INFO("T_hot = {}", temperatureHotSide);
    VF_LOG_INFO("velocity = {}", velocity);
    VF_LOG_INFO("thermal_diffusion_velocity = {}", thermalDiffusionVelocity);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Grid and boundary conditions
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
    gridBuilder->addCoarseGrid(-c1o2 * sideLength, -c1o2 * sideLength, -c1o2 * sideLength, c1o2 * sideLength,
                               c1o2 * sideLength, c1o2 * sideLength, deltaX);
    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    gridBuilder->buildGrids(false); // buildGrids() has to be called before setting the BCs!!!!

    const real vxADBC = c0o1;
    const real vyADBC = c0o1;
    const real vzADBC = c0o1;

    gridBuilder->setNoSlipBoundaryCondition(SideType::PX);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MX);
    gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
    gridBuilder->setNoSlipBoundaryCondition(SideType::PZ);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);

    gridBuilder->setADDirichletBoundaryCondition(SideType::MX, temperatureHotSide, vxADBC, vyADBC, vzADBC);
    gridBuilder->setADDirichletBoundaryCondition(SideType::PX, temperatureColdSide, vxADBC, vyADBC, vzADBC);

    gridBuilder->setADNoSlipBoundaryCondition(SideType::MY);
    gridBuilder->setADNoSlipBoundaryCondition(SideType::PY);
    gridBuilder->setADNoSlipBoundaryCondition(SideType::MZ);
    gridBuilder->setADNoSlipBoundaryCondition(SideType::PZ);

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityBounceBack);
    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipDelayBounceBack);
    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip);
    bcFactory.setAdvectionDiffusionNoSlipBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNoSlipBC::NoSlipBounceBack);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^

    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();
    auto para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(), &config);
    auto memoryManager = std::make_shared<CudaMemoryManager>(para);
    para->setInitialCondition([&](real, real, real, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx = c0o1;
        vy = c0o1;
        vz = c0o1;
    });
    para->setInitialConditionAD([&](real, real, real) { return meanTemperature; });
    para->setInitialLocalReferenceTemperature([&](real, real, real) -> real { return meanTemperature; });

    para->setOutputPrefix(simulationName);
    para->setPrintFiles(true);

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(deltaX / deltaT);
    para->setViscosityRatio(deltaX * deltaX / deltaT);
    para->setDensityRatio(1.0);

    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);

    para->setTimestepStartOut(uint(tStartOut / deltaT));
    para->setTimestepOut(uint(tOut / deltaT));
    para->setTimestepEnd(uint(tEnd / deltaT));
    para->worldLength = sideLength;

    // Advection Diffusion
    para->setDiffOn(true);
    para->setBuoyancyEnabled(true);
    para->setTurbulentPrandtlNumber(prandtlNumber);
    para->setADKernel(vf::advectionDiffusionKernel::compressible::F16);
    para->setDiffusivity(diffusivityLB);
    para->setBuoyancyFactor(thermalExpansion * gravity * (deltaT * deltaT / deltaX));

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto midPlane =
        std::make_shared<Probe>(para, memoryManager, para->getOutputPath(), "midPlane", 0, uint(tOut / deltaT), 0, uint(tOut / deltaT), false, false, true);
    midPlane->addProbePlane(-c1o2 * sideLength, 0.0, -c1o2 * sideLength, sideLength, deltaX, sideLength);
    midPlane->addAllAvailableStatistics();
    para->addSampler(midPlane);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    Simulation sim(para, memoryManager, gridBuilder, &bcFactory, tmFactory);
    sim.run();
}

int main(int argc, char* argv[])
{
    if (argv == NULL)
        return 1;

    try {
        vf::logging::Logger::initializeLogger();

        auto config = vf::basics::loadConfig(argc, argv, defaultConfigFile);

        run(config);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    } catch (const std::bad_alloc& e) {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    } catch (const std::exception& e) {
        VF_LOG_CRITICAL("exception: {}", e.what());
    } catch (...) {
        VF_LOG_CRITICAL("Unknown exception!");
    }

    return 0;
}
