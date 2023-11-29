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
//! \author Martin Schoenherr
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

#include <GridGenerator/geometries/Cuboid/Cuboid.h>
#include <GridGenerator/grid/BoundaryConditions/Side.h>
#include <GridGenerator/grid/GridBuilder/LevelGridBuilder.h>
#include <GridGenerator/grid/GridBuilder/MultipleGridBuilder.h>

#include <gpu/core/BoundaryConditions/BoundaryConditionFactory.h>
#include <gpu/core/DataStructureInitializer/GridProvider.h>
#include <gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h>
#include <gpu/core/GPU/CudaMemoryManager.h>
#include <gpu/core/GridScaling/GridScalingFactory.h>
#include <gpu/core/Kernel/KernelTypes.h>
#include <gpu/core/LBM/Simulation.h>
#include <gpu/core/Output/FileWriter.h>
#include <gpu/core/Parameter/Parameter.h>

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();

        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////

        std::string path("./output/DrivenCavity");
        std::string simulationName("LidDrivenCavity");

        const real length = 1.0;
        const real reynoldsNumber = 1000.0;
        const real velocity = 1.0;
        const uint numberOfNodesX = 64;

        const uint timeStepOut = 1000;
        const uint timeStepEnd = 10000;

        real velocityLB = 0.05; // LB units

        vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv);

        bool refine = false;
        if (config.contains("refine"))
            refine = config.getValue<bool>("refine");

        if (config.contains("output_path"))
            path = config.getValue<std::string>("output_path");

        if (config.contains("velocityLB"))
            velocityLB = config.getValue<real>("velocityLB");

        //////////////////////////////////////////////////////////////////////////
        // compute parameters in lattice units
        //////////////////////////////////////////////////////////////////////////

        const real deltaX = length / real(numberOfNodesX);
        const real deltaT = velocityLB / velocity * deltaX;

        const real vxLB = velocityLB / sqrt(2.0); // LB units
        const real vyLB = velocityLB / sqrt(2.0); // LB units

        const real viscosityLB = numberOfNodesX * velocityLB / reynoldsNumber; // LB units

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////

        auto gridBuilder = std::make_shared<MultipleGridBuilder>();

        gridBuilder->addCoarseGrid(-0.5 * length, -0.5 * length, -0.5 * length, 0.5 * length, 0.5 * length, 0.5 * length, deltaX);
        if (refine)
            gridBuilder->addGrid(std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25), 1);

        GridScalingFactory scalingFactory = GridScalingFactory();
        scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleCompressible);

        gridBuilder->setPeriodicBoundaryCondition(false, false, false);

        gridBuilder->buildGrids(false);

        //////////////////////////////////////////////////////////////////////////
        // set parameters
        //////////////////////////////////////////////////////////////////////////
        auto para = std::make_shared<Parameter>();

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

        gridBuilder->setNoSlipBoundaryCondition(SideType::PX);
        gridBuilder->setNoSlipBoundaryCondition(SideType::MX);
        gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
        gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
        gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);
        gridBuilder->setVelocityBoundaryCondition(SideType::PZ, vxLB, vyLB, 0.0);

        BoundaryConditionFactory bcFactory;

        bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocitySimpleBounceBackCompressible);

        //////////////////////////////////////////////////////////////////////////
        // set copy mesh to simulation
        //////////////////////////////////////////////////////////////////////////

        vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

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
        VF_LOG_INFO("delta t [s]            = {}", deltaT);
        VF_LOG_INFO("world_length   [m]     = {}", length);
        VF_LOG_INFO("world_velocity [m/s]   = {}", velocity);
        VF_LOG_INFO("delta x [m]            = {}", deltaX);
        printf("\n");
        VF_LOG_INFO("LB parameter:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("Re                     = {}", reynoldsNumber);
        VF_LOG_INFO("lb_velocity [dx/dt]    = {}", velocityLB);
        VF_LOG_INFO("lb_viscosity [dx^2/dt] = {}", viscosityLB);
        VF_LOG_INFO("lb_vx [dx/dt] (lb_velocity/sqrt(2)) = {}", vxLB);
        VF_LOG_INFO("lb_vy [dx/dt] (lb_velocity/sqrt(2)) = {}", vyLB);
        printf("\n");
        VF_LOG_INFO("simulation parameter:");
        VF_LOG_INFO("--------------");
        VF_LOG_INFO("number of nodes in x   = {}", numberOfNodesX);
        VF_LOG_INFO("number of nodes        = {}", numberOfNodesX * numberOfNodesX * numberOfNodesX);
        VF_LOG_INFO("write_nth_timestep     = {}", timeStepOut);
        VF_LOG_INFO("last timestep          = {}", timeStepEnd);
        VF_LOG_INFO("output_path            = {}", path);

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
        sim.run();

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
