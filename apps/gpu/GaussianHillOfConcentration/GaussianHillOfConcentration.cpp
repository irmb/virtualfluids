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

//////////////////////////////////////////////////////////////////////////

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>
#include <basics/constants/NumericConstants.h>

#include <logger/Logger.h>
#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/MultipleGridBuilderFacade.h"
#include "GridGenerator/grid/GridDimensions.h"
//////////////////////////////////////////////////////////////////////////

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"

#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/TurbulenceModels/TurbulenceModelFactory.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace vf::basics::constant;

const std::string defaultConfigPath("gaussian_hill.cfg");
const std::string simulationName("GaussianHillOfConcentration");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run(const vf::basics::ConfigurationFile& config)
{
    const real prandtlNumber = c1o1;

    const real deltaT = c1o1;
    const real sigma0 = c1o1;
    const real pecletNumber = config.getValue<real>("PecletNumber");
    const real nodesPerSigma0 = config.getValue<real>("NodesPerSigma0");
    const real diffusivityLB = config.getValue<real>("DiffusivityLB");
    const real transports = config.getValue<real>("Transports");

    const bool useDiffusionVelocity = pecletNumber < c1o1;

    const real deltaX = sigma0 / nodesPerSigma0;
    const real domainSize = c18o1 * sigma0;
    const real diffusivity = diffusivityLB * deltaX * deltaX / deltaT;
    const real advectionVelocityLB = pecletNumber * diffusivityLB / nodesPerSigma0;
    const real advectionVelocity = advectionVelocityLB * deltaX / deltaT;
    const real diffusionVelocity = diffusivity / sigma0;
    const real referenceVelocity = useDiffusionVelocity ? diffusionVelocity : advectionVelocity;
    const real velocityLB = referenceVelocity * deltaT / deltaX; // LB units
    const real viscosityLB = prandtlNumber * diffusivityLB;      // LB units

    const real transportLength = transports * sigma0;

    const real time = transportLength / referenceVelocity;

    const uint numberOfTimesteps = time / deltaT;
    const real machNumber = referenceVelocity * std::sqrt(c3o1) * deltaT / deltaX;

    VF_LOG_INFO("reference Velocity = {}", referenceVelocity);
    VF_LOG_INFO("velocity  [dx/dt] = {}", velocityLB);
    VF_LOG_INFO("dt   = {}", deltaT);
    VF_LOG_INFO("dx   = {}", deltaX);
    VF_LOG_INFO("viscosity [dx^2/dt] = {}", viscosityLB);
    VF_LOG_INFO("Peclet number = {}", pecletNumber);
    VF_LOG_INFO("Prandtl number = {}", prandtlNumber);
    VF_LOG_INFO("Mach number = {}", machNumber);
    VF_LOG_INFO("advection velocity = {}", advectionVelocity);
    VF_LOG_INFO("diffusion velocity = {}", diffusionVelocity);
    VF_LOG_INFO("tEnd = {}", time);
    VF_LOG_INFO("Number of Timesteps = {}", numberOfTimesteps);

    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

    auto grid = std::make_shared<GridDimensions>(-c1o2 * domainSize, c1o2 * domainSize, -c1o2 * domainSize,
                                                 c1o2 * domainSize, -c1o2 * domainSize, c1o2 * domainSize, deltaX);
    auto gridBuilder = std::make_shared<MultipleGridBuilderFacade>(grid, c2o1 * deltaX);

    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    for (int iProcess = 1; iProcess < communicator.getNumberOfProcesses(); iProcess++)
        gridBuilder->addDomainSplit(c0o1, Axis::x);

    gridBuilder->createGrids(communicator.getProcessID());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNumberOfProcesses(), communicator.getProcessID(),&config);

    para->setInitialCondition([&](real, real, real, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx = advectionVelocityLB;
        vy = advectionVelocityLB;
        vz = advectionVelocityLB;
    });

    para->setInitialConditionAD([&](real coordX, real coordY, real coordZ) -> real {
        const real shiftedX = coordX + time * advectionVelocity * c1o2;
        const real shiftedY = coordY + time * advectionVelocity * c1o2;
        const real shiftedZ = coordZ + time * advectionVelocity * c1o2;
        const real distSquared = shiftedX * shiftedX + shiftedY * shiftedY + shiftedZ * shiftedZ;
        return std::exp(-c1o2 * distSquared / (sigma0 * sigma0));
    });

    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);
    para->worldLength = domainSize + 2 * deltaX;
    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);
    para->setVelocityRatio(deltaX / deltaT);
    para->setViscosityRatio(deltaX * deltaX / deltaT);

    para->configureMainKernel(vf::collision_kernel::compressible::K17CompressibleNavierStokes);

    para->setTimestepStartOut(0);
    para->setTimestepOut(numberOfTimesteps);
    para->setTimestepEnd(numberOfTimesteps);

    para->setDiffOn(true);
    para->setUseStreams(true);
    para->setAllNodesAllFeatures(true);
    para->setADKernel(vf::advectionDiffusionKernel::compressible::F16);
    para->setDiffusivity(diffusivityLB);

    SPtr<TurbulenceModelFactory> tmFactory = std::make_shared<TurbulenceModelFactory>(para);
    tmFactory->readConfigFile(config);

    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    Simulation sim(para, gridBuilder->getGridBuilder(), &bcFactory, tmFactory);
    sim.run();
}

int main(int argc, char* argv[])
{
    if (argv == NULL)
        return 0;

    try {
        vf::logging::Logger::initializeLogger();
        auto config = vf::basics::loadConfig(argc, argv, defaultConfigPath);
        run(config);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    }

    catch (const std::bad_alloc& e) {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    } catch (const std::exception& e) {
        VF_LOG_CRITICAL("exception: {}", e.what());
    } catch (...) {
        VF_LOG_CRITICAL("Unknown exception!");
    }

    return 0;
}
