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
#include <basics/geometry3d/Axis.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/Side.h"
#include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "gpu/GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "gpu/GridGenerator/grid/GridDimensions.h"
#include "gpu/GridGenerator/grid/MultipleGridBuilderFacade.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PreCollisionInteractor/BuoyancyProvider/BuoyancyProvider.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"

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
    const real thermalDiffusionVelocity = c1o1;
    const real sideLength = c1o1;
    const real gravity = 9.81F;
    const real diffusivity = thermalDiffusionVelocity * sideLength;
    const real prandtlNumber = 0.71F;

    const uint numberOfNodes = config.getValue<uint>("nNodes", 16);
    const real diffusivityLB = config.getValue("diffusivityLB", c1o100 * c1o10);
    const real rayleighNumber = config.getValue<real>("Ra", 1000.F);

    const real deltaX = sideLength / real(numberOfNodes);
    const real deltaT = diffusivityLB * deltaX * deltaX / diffusivity;
    const real velocityLB = thermalDiffusionVelocity * deltaT / deltaX;

    const real viscosity = prandtlNumber * diffusivity;
    const real viscosityLB = prandtlNumber * diffusivityLB;
    const real thermalExpansion =
        (rayleighNumber * diffusivity * viscosity) / (gravity * std::pow(sideLength, c3o1) * temperatureDifference);
    const real machNumber = std::sqrt(c3o1) * velocityLB;

    // all in s
    const real tStartOut = config.getValue<real>("tStartOut");
    const real tOut = config.getValue<real>("tOut");
    const real tEnd = config.getValue<real>("tEnd"); // total time of simulation

    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("dt   = {}", deltaT);
    VF_LOG_INFO("dx   = {}", deltaX);
    VF_LOG_INFO("viscosity [dx^2/dt] = {}", viscosityLB);
    VF_LOG_INFO("Ma = {}", machNumber);
    VF_LOG_INFO("Pr = {}", prandtlNumber);
    VF_LOG_INFO("Ra = {}", rayleighNumber);
    VF_LOG_INFO("T_cold = {}", temperatureColdSide);
    VF_LOG_INFO("T_hot = {}", temperatureHotSide);
    VF_LOG_INFO("thermal_diffusion_velocity = {}", thermalDiffusionVelocity);

    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Grid and boundary conditions
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto gridBuilder = std::make_unique<MultipleGridBuilderFacade>(
        std::make_shared<GridDimensions>(-c1o2 * sideLength, c1o2 * sideLength, -c1o2 * sideLength, c1o2 * sideLength,
                                         -c1o2 * sideLength, c1o2 * sideLength, deltaX),
        c2o1 * deltaX);
    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    for (int iProcess = 1; iProcess < communicator.getNumberOfProcesses(); iProcess++)
        gridBuilder->addDomainSplit(c0o1, Axis::x);

    gridBuilder->createGrids(communicator.getProcessID()); // buildGrids() has to be called before setting the BCs!!!!

    const real vxADBC = c0o1;
    const real vyADBC = c0o1;
    const real vzADBC = c0o1;

    gridBuilder->setNoSlipBoundaryCondition(SideType::PX);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MX);
    gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
    gridBuilder->setNoSlipBoundaryCondition(SideType::PZ);
    gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);

    // only using all kinds of Bc as showcase
    gridBuilder->setADDirichletBoundaryCondition(SideType::MX, temperatureHotSide, vxADBC, vyADBC, vzADBC);
    gridBuilder->setADDirichletBoundaryCondition(SideType::PX, temperatureColdSide, vxADBC, vyADBC, vzADBC);

    gridBuilder->setADNoFluxBoundaryCondition(SideType::MY);
    gridBuilder->setADNeumannBoundaryCondition(SideType::PY, c0o1, c0o1, c0o1, c0o1, deltaX);
    gridBuilder->setADFluxBoundaryCondition(SideType::MZ, c0o1, c0o1, c1o1, c0o1, deltaX);
    gridBuilder->setADFluxBoundaryCondition(SideType::PZ, c0o1, c0o1, -c1o1, c0o1, deltaX);
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityBounceBack);
    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipDelayBounceBack);
    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(
        BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip);
    bcFactory.setAdvectionDiffusionNeumannBoundaryCondition(
        BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip);
    bcFactory.setAdvectionDiffusionNoFluxBoundaryCondition(
        BoundaryConditionFactory::AdvectionDiffusionNoFluxBC::NoFluxBounceBack);
    bcFactory.setAdvectionDiffusionFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxBounceBack);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    para->setIsBodyForce(true);
    para->setBuoyancyEnabled(true);
    para->setTurbulentPrandtlNumber(prandtlNumber);
    para->setADKernel(vf::advectionDiffusionKernel::compressible::F16);
    para->setDiffusivity(diffusivityLB);
    para->setBuoyancyFactor(thermalExpansion * gravity * (deltaT * deltaT / deltaX));
    auto buyoancyProvider = std::make_shared<BuoyancyProviderConstantValue>(para, memoryManager);
    para->addInteractor(buyoancyProvider);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto midPlane = std::make_shared<Probe>(para, memoryManager, para->getOutputPath(), "midPlane", 0, uint(tOut / deltaT),
                                            0, uint(tOut / deltaT), false, false, true);
    midPlane->addProbePlane(-c1o2 * sideLength, 0.0, -c1o2 * sideLength, sideLength, deltaX, sideLength);
    midPlane->addAllAvailableStatistics();
    para->addSampler(midPlane);

    auto sidePlane = std::make_shared<Probe>(para, memoryManager, para->getOutputPath(), "sidePlane", 0, uint(tOut / deltaT),
                                             0, uint(tOut / deltaT), false, false, true);
    sidePlane->addProbePlane(c1o2 * sideLength - deltaX, -c1o2 * sideLength, -c1o2 * sideLength, deltaX, sideLength,
                             sideLength);
    sidePlane->addAllAvailableStatistics();
    para->addSampler(sidePlane);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    Simulation sim(para, memoryManager, gridBuilder->getGridBuilder(), &bcFactory, tmFactory);
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
