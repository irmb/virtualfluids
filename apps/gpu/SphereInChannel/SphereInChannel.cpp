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
//! \addtogroup SphereInChannel
//! \ingroup gpu_apps
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
#include <string>

#include <basics/DataTypes.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/Calculation/Simulation.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Samplers/Probe.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"

void run(const vf::basics::ConfigurationFile& config)
{
    //////////////////////////////////////////////////////////////////////////
    // Simulation parameters
    //////////////////////////////////////////////////////////////////////////
    std::string path("output/Sphere");
    std::string simulationName("SphereInChannel");

    const real reynoldsNumber = 300.0;
    const real velocity = 1.0;
    const real velocityLB = (real)0.5e-2; // LB units
    const uint numberOfNodesSphere = 10;
    const real dSphere =
        0.2; // Caution, this App uses an STL for the sphere, so this setting does not change the size of the sphere

    const uint timeStepOut = 10000;
    const uint timeStepEnd = 10000;

    bool refine = false;
    if (config.contains("refine"))
        refine = config.getValue<bool>("refine");

    if (config.contains("output_path"))
        path = config.getValue<std::string>("output_path");

    //////////////////////////////////////////////////////////////////////////
    // compute parameters in lattice units
    //////////////////////////////////////////////////////////////////////////

    const real deltaX = dSphere / real(numberOfNodesSphere);
    const real viscosityLB = velocityLB / reynoldsNumber; // LB units

    //////////////////////////////////////////////////////////////////////////
    // create grid
    //////////////////////////////////////////////////////////////////////////

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    gridBuilder->addCoarseGrid(-3 * dSphere, -3 * dSphere, -3 * dSphere, 10 * dSphere, 3 * dSphere, 3 * dSphere, deltaX);

    // add geometry: use primitive
    // auto sphere = std::make_shared<Sphere>(0.0, 0.0, 0.0, dSphere / 2.0);

    // add geometry: use stl
    std::string stlPath = "./apps/gpu/SphereInChannel/sphere02.stl";
    std::cout << "Reading stl from " << stlPath << "." << std::endl;
    auto sphere = std::make_shared<TriangularMesh>(stlPath);

    gridBuilder->addGeometry(sphere);
    gridBuilder->setPeriodicBoundaryCondition(false, false, false);

    // add two fine grids
    if (refine) {
        gridBuilder->addGrid(std::make_shared<Sphere>(0., 0., 0., dSphere), 2);
    }

    GridScalingFactory scalingFactory = GridScalingFactory();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    gridBuilder->buildGrids(false);

    //////////////////////////////////////////////////////////////////////////
    // set parameters
    //////////////////////////////////////////////////////////////////////////
    SPtr<Parameter> para = std::make_shared<Parameter>();

    para->setOutputPath(path);
    para->setOutputPrefix(simulationName);

    para->setPrintFiles(true);

    para->setVelocityLB(velocityLB);
    para->setViscosityLB(viscosityLB);

    para->setVelocityRatio(velocity / velocityLB);
    para->setDensityRatio(1.0);

    para->setTimestepOut(timeStepOut);
    para->setTimestepEnd(timeStepEnd);

    para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

    //////////////////////////////////////////////////////////////////////////
    // set boundary conditions
    //////////////////////////////////////////////////////////////////////////

    gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0);
    gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
    gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);

    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);

    BoundaryConditionFactory bcFactory;

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

    //////////////////////////////////////////////////////////////////////////
    // setup probe(s)
    //////////////////////////////////////////////////////////////////////////

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);

    const uint tStartAveraging = 0;
    const uint tAveraging = 100;
    const uint tStartOutProbe = 0;
    const uint tOutProbe = para->getTimestepOut();
    SPtr<Probe> pointProbe = std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "pointProbe",
                                                     tStartAveraging, tAveraging, tStartOutProbe, tOutProbe, false, false);
    std::vector<real> probeCoordsX = { 0.3, 0.5 };
    std::vector<real> probeCoordsY = { 0.0, 0.0 };
    std::vector<real> probeCoordsZ = { 0.0, 0.0 };
    pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);

    pointProbe->addStatistic(Probe::Statistic::Instantaneous);
    pointProbe->addStatistic(Probe::Statistic::Means);
    pointProbe->addStatistic(Probe::Statistic::Variances);
    para->addSampler(pointProbe);

    SPtr<Probe> planeProbe = std::make_shared<Probe>(para, cudaMemoryManager, para->getOutputPath(), "planeProbe",
                                                     tStartAveraging, tAveraging, tStartOutProbe, tOutProbe, false, false);
    planeProbe->addProbePlane(0.4, 0, 0, 0.3, 0.01, 0.1);
    planeProbe->addStatistic(Probe::Statistic::Instantaneous);
    para->addSampler(planeProbe);

    //////////////////////////////////////////////////////////////////////////
    // initial state of the flow field
    //////////////////////////////////////////////////////////////////////////

    para->setInitialCondition([&](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) {
        rho = c0o1;
        vx = velocityLB;
        vy = c0o1;
        vz = c0o1;
    });

    //////////////////////////////////////////////////////////////////////////
    // run simulation
    //////////////////////////////////////////////////////////////////////////

    Simulation simulation(para, cudaMemoryManager, gridBuilder, &bcFactory, &scalingFactory);
    simulation.run();
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        vf::basics::ConfigurationFile config =
            vf::basics::loadConfig(argc, argv, "/workspaces/VirtualFluids/apps/gpu/SphereInChannel/sphere_1level.cfg");
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}

//! \}
