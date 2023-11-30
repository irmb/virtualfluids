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
//! \author Martin Schoenherr, Stephan Lenz, Anna Wellmann
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
#include <basics/config/ConfigurationFile.h>
#include <logger/Logger.h>
#include <parallel/MPICommunicator.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/geometries/Cuboid/Cuboid.h"
#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

//////////////////////////////////////////////////////////////////////////

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/GridScaling/GridScalingFactory.h"
#include "gpu/core/Kernel/KernelTypes.h"
#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/PreCollisionInteractor/Probes/PlaneProbe.h"
#include "gpu/core/PreCollisionInteractor/Probes/PointProbe.h"

//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        vf::parallel::Communicator& communicator = *vf::parallel::MPICommunicator::getInstance();

        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////
        std::string path("output/Sphere");
        std::string simulationName("SphereInChannel");

        const real length = 1.0;
        const real reynoldsNumber = 300.0;
        const real velocity = 1.0;
        const real velocityLB = (real)0.5e-2; // LB units
        const uint numberOfNodesX = 50;
        const real dSphere = 0.2;

        const uint timeStepOut = 10000;
        const uint timeStepEnd = 10000;

        vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv);

        bool refine = false;
        if (config.contains("refine"))
            refine = config.getValue<bool>("refine");

        if (config.contains("output_path"))
            path = config.getValue<std::string>("output_path");

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

        gridBuilder->addCoarseGrid(-0.6 * length, -0.6 * length, -0.6 * length,
                                    2.0 * length,  0.6 * length,  0.6 * length, deltaX);

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

        const uint tStartAveraging = 0;
        const uint tAveraging = 100;
        const uint tStartOutProbe = 0;
        const uint tOutProbe = para->getTimestepOut();
        SPtr<PointProbe> pointProbe = std::make_shared<PointProbe>("pointProbe", para->getOutputPath(), tStartAveraging,
                                                                   tAveraging, tStartOutProbe, tOutProbe);
        std::vector<real> probeCoordsX = { 0.3, 0.5 };
        std::vector<real> probeCoordsY = { 0.0, 0.0 };
        std::vector<real> probeCoordsZ = { 0.0, 0.0 };
        pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);

        pointProbe->addStatistic(Statistic::Instantaneous);
        pointProbe->addStatistic(Statistic::Means);
        pointProbe->addStatistic(Statistic::Variances);
        para->addProbe(pointProbe);

        SPtr<PlaneProbe> planeProbe = std::make_shared<PlaneProbe>("planeProbe", para->getOutputPath(), tStartAveraging,
                                                                   tAveraging, tStartOutProbe, tOutProbe);
        planeProbe->setProbePlane(0.4, 0, 0, 0.3, 0.01, 0.1);
        planeProbe->addStatistic(Statistic::Instantaneous);
        para->addProbe(planeProbe);

        //////////////////////////////////////////////////////////////////////////
        // set copy mesh to simulation
        //////////////////////////////////////////////////////////////////////////

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
