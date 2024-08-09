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
//! \addtogroup LaminarPipeFlow
//! \ingroup gpu_apps
//! \{
//! \author Anna Wellmann
//=======================================================================================
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/config/ConfigurationFile.h>

#include <logger/Logger.h>

#include <parallel/MPICommunicator.h>

#include <GridGenerator/geometries/Conglomerate/Conglomerate.h>
#include <GridGenerator/geometries/Cuboid/Cuboid.h>
#include <GridGenerator/geometries/Cylinder/Cylinder.h>
#include <GridGenerator/grid/BoundaryConditions/BoundaryCondition.h>
#include <GridGenerator/grid/BoundaryConditions/Side.h>
#include <GridGenerator/grid/GridBuilder/LevelGridBuilder.h>
#include <GridGenerator/grid/GridBuilder/MultipleGridBuilder.h>

#include <gpu/core/BoundaryConditions/BoundaryConditionFactory.h>
#include <gpu/core/Calculation/Simulation.h>
#include <gpu/core/Cuda/CudaMemoryManager.h>
#include <gpu/core/DataStructureInitializer/GridProvider.h>
#include <gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h>
#include <gpu/core/GridScaling/GridScalingFactory.h>
#include <gpu/core/Kernel/KernelTypes.h>
#include <gpu/core/Output/FileWriter.h>
#include <gpu/core/Parameter/Parameter.h>

void run(const vf::basics::ConfigurationFile& config)
{
    // simulation parameters

    std::string path("./output/LaminarPipeFlow");
    std::string simulationName("LaminarPipeFlow");

    const std::array<real, 3> length = { 128, 64, 64 };
    const real radius = length[1] / 2.0;
    const real deltaX = config.getValue<real>("deltaX");
    const real diameterLB = length[1] / deltaX;
    const real lengthLB = length[0] / deltaX;
    const real radiusLB = diameterLB / 2.;
    const real density1 = 1e-5;
    const real density2 = 0.0;
    const real viscosityLB = 0.0064 / deltaX; // same max velocity and Re an all resolutions (diffusive scaling)

    const real pressureDifferenceLB = (density1 / 3. - density2 / 3.);
    const real velocityMax = radiusLB * radiusLB / (4. * viscosityLB) * (pressureDifferenceLB / lengthLB);
    const real reynoldsNumber = velocityMax * 2 * radiusLB / viscosityLB;

    const uint timeStepOut = config.getValue<uint>("timeStepOutput");
    const uint timeStepEnd = config.getValue<uint>("timeStepEnd");

    // create grid

    auto gridBuilder = std::make_shared<MultipleGridBuilder>();

    // create grid: bounding box
    const real boundingBoxMinX1 = 0.0;
    const real boundingBoxMinX2 = -length[1] / 2.0;
    const real boundingBoxMinX3 = -length[2] / 2.0;
    const real boundingBoxMaxX1 = length[0];
    const real boundingBoxMaxX2 = length[1] / 2.0;
    const real boundingBoxMaxX3 = length[2] / 2.0;
    gridBuilder->addCoarseGrid(boundingBoxMinX1, boundingBoxMinX2, boundingBoxMinX3, boundingBoxMaxX1, boundingBoxMaxX2,
                               boundingBoxMaxX3, deltaX);

    // create grid: pipe
    const real padding = 2.0 * deltaX;
    std::array<double, 3> cylinderAxisCoordinate1 { boundingBoxMinX1 - padding, 0., 0. };
    std::array<double, 3> cylinderAxisCoordinate2 { boundingBoxMaxX1 + padding, 0., 0. };
    auto cylinder = std::make_shared<Cylinder>(cylinderAxisCoordinate1, cylinderAxisCoordinate2, radius);
    auto pipe = std::make_shared<Conglomerate>();
    pipe->add(std::make_shared<Cuboid>(boundingBoxMinX1 - padding, boundingBoxMinX2 - padding, boundingBoxMinX3 - padding,
                                       boundingBoxMaxX1 + padding, boundingBoxMaxX2 + padding, boundingBoxMaxX3 + padding));
    pipe->subtract(cylinder);
    gridBuilder->addGeometry(pipe);

    // create grid: cuboidal pipe as refinement
    auto cuboid = std::make_shared<Cuboid>(boundingBoxMinX1, 0.5 * boundingBoxMinX2, 0.5 * boundingBoxMinX3,
                                           boundingBoxMaxX1, 0.5 * boundingBoxMaxX2, 0.5 * boundingBoxMaxX3);
    auto pipeRefinement = std::make_shared<Conglomerate>();
    pipeRefinement->add(std::make_shared<Cuboid>(boundingBoxMinX1 - padding, boundingBoxMinX2 - padding,
                                                 boundingBoxMinX3 - padding, boundingBoxMaxX1 + padding,
                                                 boundingBoxMaxX2 + padding, boundingBoxMaxX3 + padding));
    pipeRefinement->subtract(cuboid);
    gridBuilder->addGrid(pipeRefinement);

    // create grid: round pipe as refinement
    // auto cylinderRefinement =
    //     std::make_shared<Cylinder>(cylinderAxisCoordinate1, cylinderAxisCoordinate2, 0.5 * radius);
    // auto pipeRefinement = std::make_shared<Conglomerate>();
    // pipeRefinement->add(std::make_shared<Cuboid>(boundingBoxMinX1 - padding, boundingBoxMinX2 - padding,
    //                                              boundingBoxMinX3 - padding, boundingBoxMaxX1 + padding,
    //                                              boundingBoxMaxX2 + padding, boundingBoxMaxX3 + padding));
    // pipeRefinement->subtract(cylinderRefinement);
    // gridBuilder->addGrid(pipeRefinement);

    GridScalingFactory scalingFactory = GridScalingFactory();
    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

    // create grid: build grids
    gridBuilder->buildGrids(false);

    // set parameters

    auto para = std::make_shared<Parameter>();
    para->setOutputPath(path);
    para->setOutputPrefix(simulationName);
    para->setPrintFiles(true);
    para->setViscosityLB(viscosityLB);
    para->setDensityRatio(1.0);
    para->setVelocityRatio(1.0);
    para->setTimestepOut(timeStepOut);
    para->setTimestepEnd(timeStepEnd);

    para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

    // set boundary conditions

    gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0, 0, 0);
    gridBuilder->setPressureBoundaryCondition(SideType::MX, density1);
    gridBuilder->setPressureBoundaryCondition(SideType::PX, density2);

    // set velocity profile as inflow
    // gridBuilder->setVelocityBoundaryCondition(SideType::MX, u_max, 0.0, 0.0);
    // const uint level = 0;
    // auto inflowBC =
    //     std::dynamic_pointer_cast<VelocityBoundaryCondition>(gridBuilder->getBoundaryCondition(SideType::MX, level));
    // inflowBC->setVelocityProfile(gridBuilder->getGrid(level), [&](real x, real y, real z, real& vx, real& vy, real& vz) {
    //     auto distanceToCenter = sqrt(y * y + z * z);
    //     vx = -u_max / (radiusLB * radiusLB) * distanceToCenter * distanceToCenter + u_max; // parabolic profile
    //     vy = 0.0;
    //     vz = 0.0;
    // });

    BoundaryConditionFactory bcFactory;
    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

    // prepare mesh copy from the grid builder to the simulation

    vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

    auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
    SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

    // log parameters

    VF_LOG_INFO("parameters:");
    VF_LOG_INFO("--------------");
    VF_LOG_INFO("pressure1                     = {}", density1 / 3.);
    VF_LOG_INFO("pressure2                     = {}", density2 / 3.);
    VF_LOG_INFO("maximum velocity              = {}", velocityMax);
    VF_LOG_INFO("Re (reference: diameter pipe) = {}", reynoldsNumber);
    VF_LOG_INFO("deltaX                        = {}", deltaX);
    VF_LOG_INFO("number of nodes in diameter   = {}", diameterLB);
    VF_LOG_INFO("write_nth_timestep            = {}", timeStepOut);
    VF_LOG_INFO("last timestep                 = {}", timeStepEnd);
    VF_LOG_INFO("output_path                   = {}", path);

    // run simulation

    Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
    VF_LOG_INFO("Starting simulation.");
    sim.run();
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        vf::basics::ConfigurationFile config =
            vf::basics::loadConfig(argc, argv, "./apps/gpu/LaminarPipeFlowGPU/laminarpipeflow.cfg");
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}