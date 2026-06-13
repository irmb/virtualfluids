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
//! \addtogroup TaylorGreenVortex
//! \ingroup gpu_apps
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

#include "GridGenerator/geometries/Conglomerate/Conglomerate.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/io/GridVTKWriter/GridVTKWriter.h"
#include "GridGenerator/io/STLReaderWriter/STLReader.h"
#include "GridGenerator/io/STLReaderWriter/STLWriter.h"
#include "GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h"

#include "constants/NumericConstants.h"
#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/Parameter/Parameter.h"

using namespace vf::gpu;

void run(const vf::basics::ConfigurationFile& config)
{
    //////////////////////////////////////////////////////////////////////////
    // Simulation parameters
    //////////////////////////////////////////////////////////////////////////
    std::string path("./output/TaylorGreenVortex");
    std::string simulationName("TaylorGreenVortex");

    uint numberOfNodesX = 64;
    real reynoldsNumber = 1600.0;

    bool useLimiter = false;

    uint dtPerL = 250;

    uint gpuIndex = 0;

    std::string kernel("K17CompressibleNavierStokes");

    if (config.contains("output_path"))
        path = config.getValue<std::string>("output_path");

    if (config.contains("numberOfNodesX"))
        numberOfNodesX = config.getValue<uint>("numberOfNodesX");

    if (config.contains("reynoldsNumber"))
        numberOfNodesX = config.getValue<uint>("reynoldsNumber");

    if (config.contains("useLimiter"))
        useLimiter = config.getValue<bool>("useLimiter");

    if (config.contains("kernel"))
        kernel = config.getValue<std::string>("kernel");

    if (config.contains("dtPerL"))
        dtPerL = config.getValue<uint>("dtPerL");

    if (config.contains("gpuIndex"))
        gpuIndex = config.getValue<uint>("gpuIndex");

    real length = numberOfNodesX / (2.0 * vf::basics::constant::cPi);

    const real velocityLB = 64.0 / (dtPerL * 2.0 * vf::basics::constant::cPi);
    const real viscosityLB = numberOfNodesX / (2.0 * vf::basics::constant::cPi) * velocityLB / reynoldsNumber;

    //////////////////////////////////////////////////////////////////////////
    // create grid
    //////////////////////////////////////////////////////////////////////////

    real deltaX = 2.0 * vf::basics::constant::cPi / real(numberOfNodesX);
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    gridBuilder->addCoarseGrid(-vf::basics::constant::cPi, -vf::basics::constant::cPi, -vf::basics::constant::cPi,
                               vf::basics::constant::cPi, vf::basics::constant::cPi, vf::basics::constant::cPi, deltaX);

    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

    gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

    //////////////////////////////////////////////////////////////////////////
    // set parameters
    //////////////////////////////////////////////////////////////////////////

    auto para = std::make_shared<Parameter>();

    para->setDevices(std::vector<uint> { gpuIndex });

    para->setOutputPath(path);
    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);

    para->setTimestepEnd(40 * lround(length / velocityLB));
    para->setTimestepOut(5 * lround(length / velocityLB));

    para->setVelocityLB(velocityLB);

    para->setViscosityLB(viscosityLB);

    para->setVelocityRatio(1.0 / velocityLB);

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) {
        real a = 1.0;
        real b = 1.0;
        real c = 1.0;

        rho = 3.0 * ((velocityLB * velocityLB) / 16.0 * (cos(2.0 * a * coordX) + cos(2.0 * b * coordY)) *
                     (cos(2.0 * c * coordZ) + 2.0));
        vx = velocityLB * sin(a * coordX) * cos(b * coordY) * cos(c * coordZ);
        vy = -velocityLB * cos(a * coordX) * sin(b * coordY) * cos(c * coordZ);
        vz = 0.0;
    });

    para->configureMainKernel(kernel);

    if (!useLimiter)
        para->setQuadricLimiters(1000000.0, 1000000.0, 1000000.0);

    para->setUseInitNeq(true);

    //////////////////////////////////////////////////////////////////////////

    BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

    //////////////////////////////////////////////////////////////////////////
    // run simulation
    //////////////////////////////////////////////////////////////////////////

    Simulation simulation(para, gridBuilder, &bcFactory);
    simulation.run();

    simulation.addKineticEnergyAnalyzer(10);
    simulation.addEnstrophyAnalyzer(10);

    simulation.run();
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        vf::basics::ConfigurationFile config =
            vf::basics::loadConfig(argc, argv, "./apps/gpu/TaylorGreenVortex/taylorGreenVortex.cfg");

        run(config);

    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }

    return 0;
}

//! \}
