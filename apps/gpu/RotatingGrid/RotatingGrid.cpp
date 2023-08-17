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
//! \file RotatingGrid.cpp
//! \ingroup Applications
//! \author Anna Wellmann
//=======================================================================================
#define _USE_MATH_DEFINES
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
#include <logger/Logger.h>
#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/VerticalCylinder/VerticalCylinder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Factories/GridScalingFactory.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/Kernel/Utilities/KernelTypes.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Output/NeighborDebugWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

//////////////////////////////////////////////////////////////////////////

int main()
{
    try {
         vf::logging::Logger::initializeLogger();
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////
        enum RotationOrInterpolation {Rot, Int};
        const RotationOrInterpolation rotOrInt = Int;

        const std::string path("./output/RotatingGrid");
        const std::string simulationName = rotOrInt == Int ? "RotatingGridInterpolationTest" : "RotatingGrid";

        const real L = 1.0;
        const real Re = 2000.0;
        const real velocity = 1.0;
        const real velocityLB = 0.05; // LB units
        const uint nx = 64;

        const uint timeStepOut = 10000;
        const uint timeStepEnd = 10000;

        //////////////////////////////////////////////////////////////////////////
        // compute parameters in lattice units
        //////////////////////////////////////////////////////////////////////////

        const real dx = L / real(nx);
        const real dt  = velocityLB / velocity * dx;

        const real viscosityLB = nx * velocityLB / Re; // LB units

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////
        auto gridBuilder = std::make_shared<MultipleGridBuilder>();

        gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -1.5 * L, 0.5 * L, 0.5 * L, 1.5 * L, dx);

        if (rotOrInt == Rot) gridBuilder->addGridRotatingGrid(std::make_shared<VerticalCylinder>(0.0, 0.0, 0.0, 0.3 * L, 1. * L));
        if (rotOrInt == Int) gridBuilder->addGrid(std::make_shared<VerticalCylinder>(0.0, 0.0, 0.0, 0.3 * L, 1. * L), 1);

        GridScalingFactory scalingFactory = GridScalingFactory();
        scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

        gridBuilder->setPeriodicBoundaryCondition(false, false, false);

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

        para->setVelocityRatio(1.0); // velocity / velocityLB);
        para->setDensityRatio(1.0);

        para->setTimestepOut(timeStepOut);
        para->setTimestepEnd(timeStepEnd);

        para->configureMainKernel(vf::CollisionKernel::Compressible::K17CompressibleNavierStokes);

        //////////////////////////////////////////////////////////////////////////
        // set boundary conditions
        //////////////////////////////////////////////////////////////////////////

        gridBuilder->setSlipBoundaryCondition(SideType::MX, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::PX, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::MZ, 0.0, 0.0, velocityLB);
        gridBuilder->setPressureBoundaryCondition(SideType::PZ, 0.0);

        BoundaryConditionFactory bcFactory;

        bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
        // bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
        //////////////////////////////////////////////////////////////////////////
        // set copy mesh to simulation
        //////////////////////////////////////////////////////////////////////////

        vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();

        auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
        SPtr<GridProvider> gridGenerator =
            GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

        //////////////////////////////////////////////////////////////////////////
        // run simulation
        //////////////////////////////////////////////////////////////////////////

        VF_LOG_INFO("Start Running DrivenCavity Showcase...");
        printf("\n");
        VF_LOG_INFO("world parameter:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("dt [s]                 = {}", dt);
        VF_LOG_INFO("world_length   [m]     = {}", L);
        VF_LOG_INFO("world_velocity [m/s]   = {}", velocity);
        VF_LOG_INFO("dx [m]                 = {}", dx);
        printf("\n");
        VF_LOG_INFO("LB parameter:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("Re                     = {}", Re);
        VF_LOG_INFO("lb_velocity [dx/dt]    = {}", velocityLB);
        VF_LOG_INFO("lb_viscosity [dx^2/dt] = {}", viscosityLB);
        printf("\n");
        VF_LOG_INFO("simulation parameter:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("nx                     = {}", nx);
        VF_LOG_INFO("ny                     = {}", nx);
        VF_LOG_INFO("nz                     = {}", nx);
        VF_LOG_INFO("number of nodes        = {}", nx * nx * nx);
        VF_LOG_INFO("n timesteps            = {}", timeStepOut);
        VF_LOG_INFO("write_nth_timestep     = {}", timeStepEnd);
        VF_LOG_INFO("output_path            = {}", path);

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);

        const std::string gridName = rotOrInt == Rot ? "rot_grid" : "grid";
        gridBuilder->writeGridsToVtk(para->getOutputPath() + gridName);
        // NeighborDebugWriter::writeNeighborLinkLinesDebug(para.get());

        // sim.run();

    } catch (const spdlog::spdlog_ex &ex) {
        std::cout << "Log initialization failed: " << ex.what() << std::endl;
    } catch (const std::bad_alloc &e) {
        VF_LOG_CRITICAL("Bad Alloc: {}", e.what());
    } catch (const std::exception &e) {
        VF_LOG_CRITICAL("exception: {}", e.what());
    } catch (...) {
        VF_LOG_CRITICAL("Unknown exception!");
    }

    return 0;
}
