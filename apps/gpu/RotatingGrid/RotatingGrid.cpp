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
#include "GridGenerator/geometries/Cylinder/Cylinder.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

//////////////////////////////////////////////////////////////////////////

#include "core/DataStructureInitializer/GridProvider.h"
#include "core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "core/Factories/BoundaryConditionFactory.h"
#include "core/Factories/GridScalingFactory.h"
#include "core/GPU/CudaMemoryManager.h"
#include "core/Kernel/KernelTypes.h"
#include "core/LBM/Simulation.h"
#include "core/Output/FileWriter.h"
#include "core/Output/NeighborDebugWriter.h"
#include "core/Parameter/Parameter.h"
#include "core/Parameter/ParameterRotatingGrid.h"

//////////////////////////////////////////////////////////////////////////

int main()
{
    try {
         vf::logging::Logger::initializeLogger();
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////
        enum RotationOrInterpolation {Rot, Int};
        const RotationOrInterpolation rotOrInt = Rot;
        if (rotOrInt == Int) VF_LOG_INFO("Use interpolation.");
        if (rotOrInt == Rot) VF_LOG_INFO("Use rotation.");

        const std::string path("./output/RotatingGrid");
        const std::string simulationName = rotOrInt == Int ? "RotatingGridInterpolationTest" : "RotatingGrid"; //PiHalbe";

        const real L = 1.0;
        const real Re = 2000.0;
        const real velocity = 1.0;
        const real velocityLB = 0.05; // LB units
        const uint nx = 64;

        const uint timeStepOut = 1;
        const uint timeStepEnd = 30;
        const uint timeStepStartOutput = 10;

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

        gridBuilder->addCoarseGrid(-0.5 * L, -0.5 * L, -0.5 * L, 0.5 * L, 0.5 * L, 0.5 * L, dx);

        if (rotOrInt == Rot) gridBuilder->addGridRotatingGrid(std::make_shared<Cylinder>(0.0, 0.0, 0.0, 0.25 * L, 0.75 * L, Axis::z));
        if (rotOrInt == Int) gridBuilder->addGrid(std::make_shared<Cylinder>(0.2, 0.1, 0.1, 0.25 * L, 0.8 * L, Axis::x), 1);

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
        para->setTimestepStartOut(timeStepStartOutput);

        para->configureMainKernel(vf::collisionKernel::compressible::K15CompressibleNavierStokes);

        //////////////////////////////////////////////////////////////////////////
        // set boundary conditions
        //////////////////////////////////////////////////////////////////////////

        // gridBuilder->setSlipBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
        // gridBuilder->setSlipBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
        // gridBuilder->setSlipBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);
        // gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);

        gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
        gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
        gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);
        gridBuilder->setNoSlipBoundaryCondition(SideType::PZ);

        gridBuilder->setVelocityBoundaryCondition(SideType::PX, -velocityLB, 0.0, 0.0);
        gridBuilder->setPressureBoundaryCondition(SideType::MX, 0.0);

        BoundaryConditionFactory bcFactory;

        bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
        // bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);

        //////////////////////////////////////////////////////////////////////////
        // set initial condition
        //////////////////////////////////////////////////////////////////////////

        auto setPressPoint = [](std::array<real, 3> coordPressPoint, real rhoPressPoint, real dx,
                                std::array<real, 3> coordinates) -> real {
            auto isCoordinateInPressPoint = [&](uint dimension) -> bool {
                return coordinates[dimension] > coordPressPoint[dimension] - dx &&
                       coordinates[dimension] < coordPressPoint[dimension] + dx;
            };
            if (isCoordinateInPressPoint(0) && isCoordinateInPressPoint(1) && isCoordinateInPressPoint(2))
                return rhoPressPoint;
            else
                return (real)0.0;
        };

        // para->setInitialCondition([&](real coordX, real coordY, real coordZ, real &rho, real &vx, real &vy, real &vz) {
        //     // rho = (real) 0.0;
        //     // if (coordX > -0.52 && coordY > 0.27 && coordZ > -0.1 && coordX < -0.49 && coordY < 0.3 && coordZ < 0.11) rho = 1e-5;
        //     rho =  setPressPoint({-0.2, 0.24, 0.0}, 1e-5, dx, {coordX, coordY, coordZ});
        //     vx = 0.0;
        //     vy = 0.0;
        //     vz = 0.0;
        // });

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
        // gridBuilder->writeGridsToVtk(para->getOutputPath() + gridName);
        // NeighborDebugWriter::writeNeighborLinkLinesDebug(para.get());

        SPtr<ParameterRotatingGrid> paraRot = para->getRotatingGridParameter();
        sim.run();

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
