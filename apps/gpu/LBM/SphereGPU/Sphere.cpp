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
//! \file LidDrivenCavity.cpp
//! \ingroup Applications
//! \author Martin Schoenherr, Stephan Lenz, Anna Wellmann
//=======================================================================================
#define _USE_MATH_DEFINES
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

//////////////////////////////////////////////////////////////////////////
#include <basics/PointerDefinitions.h>
#include <basics/DataTypes.h>
#include <logger/Logger.h>
#include <basics/PointerDefinitions.h>
#include <basics/config/ConfigurationFile.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/Communication/MpiCommunicator.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Factories/GridScalingFactory.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PointProbe.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/Probes/PlaneProbe.h"

//////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    try {
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////

        const bool useConfigFile = true;

        const real L = 1.0;
        const real dSphere = 0.2;
        const real Re = 300.0; // related to the sphere's diameter
        const real velocity = 1.0;
        const real dt = (real)0.5e-3;
        const uint nx = 50;

        const uint timeStepOut = 10000;
        const uint timeStepEnd = 10000;

        //////////////////////////////////////////////////////////////////////////
        // setup simulation parameters (with or without config file)
        //////////////////////////////////////////////////////////////////////////

        SPtr<Parameter> para;
        BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
        GridScalingFactory scalingFactory = GridScalingFactory();
        vf::basics::ConfigurationFile config;
        if (useConfigFile) {
            VF_LOG_TRACE("For the default config path to work, execute the app from the project root.");
            vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv, "./apps/gpu/LBM/SphereGPU/config.txt");
            para = std::make_shared<Parameter>(&config);
        } else {
            para = std::make_shared<Parameter>();
        }


        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////
        auto gridBuilder = std::make_shared<MultipleGridBuilder>();

        real dx = L / real(nx);
        gridBuilder->addCoarseGrid(-1.0 * L, -0.6 * L, -0.6 * L,
                                    3.0 * L,  0.6 * L,  0.6 * L, dx);

        // use primitive
        // auto sphere = std::make_shared<Sphere>(0.0, 0.0, 0.0, dSphere / 2.0);

        // use stl
        std::string stlPath = "./apps/gpu/LBM/SphereGPU/sphere02.stl";
        if (useConfigFile && config.contains("STLPath")) {
            stlPath = config.getValue<std::string>("STLPath");
        }
        std::cout << "Reading stl from " << stlPath << "." << std::endl;
        auto sphere = std::make_shared<TriangularMesh>(stlPath);

        gridBuilder->addGeometry(sphere);
        gridBuilder->setPeriodicBoundaryCondition(false, false, false);

        //////////////////////////////////////////////////////////////////////////
        // add grid refinement
        //////////////////////////////////////////////////////////////////////////

        // gridBuilder->setNumberOfLayers(10, 8);
        // gridBuilder->addGrid(std::make_shared<Sphere>(0.0, 0.0, 0.0, 2.0 * dSphere), 1);
        // para->setMaxLevel(2);
        // scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleK17);

        //////////////////////////////////////////////////////////////////////////
        // build grid
        //////////////////////////////////////////////////////////////////////////

        gridBuilder->buildGrids(false);  // buildGrids() has to be called before setting the BCs!!!!

        //////////////////////////////////////////////////////////////////////////
        // compute parameters in lattice units
        //////////////////////////////////////////////////////////////////////////

        const real velocityLB = velocity * dt / dx; // LB units
        const real viscosityLB =  (dSphere / dx) * velocityLB / Re; // LB units

        VF_LOG_INFO("LB parameters:");
        VF_LOG_INFO("velocity LB [dx/dt]              = {}", velocityLB);
        VF_LOG_INFO("viscosity LB [dx/dt]             = {}", viscosityLB);

        //////////////////////////////////////////////////////////////////////////
        // set parameters
        //////////////////////////////////////////////////////////////////////////

        para->setPrintFiles(true);

        para->setVelocityLB(velocityLB);
        para->setViscosityLB(viscosityLB);

        para->setVelocityRatio(velocity / velocityLB);
        para->setDensityRatio((real)1.0);

        para->setTimestepOut(timeStepOut);
        para->setTimestepEnd(timeStepEnd);

        //////////////////////////////////////////////////////////////////////////
        // set boundary conditions
        //////////////////////////////////////////////////////////////////////////

        gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);

        gridBuilder->setSlipBoundaryCondition(SideType::PY, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::MY, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::PZ, 0.0, 0.0, 0.0);
        gridBuilder->setSlipBoundaryCondition(SideType::MZ, 0.0, 0.0, 0.0);

        gridBuilder->setVelocityBoundaryCondition(SideType::GEOMETRY, 0.0, 0.0, 0.0);
        gridBuilder->setPressureBoundaryCondition(SideType::PX, 0.0); // set pressure boundary condition last

        bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
        bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
        bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
        bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);

        //////////////////////////////////////////////////////////////////////////
        // setup probe(s)
        //////////////////////////////////////////////////////////////////////////

        const uint tStartAveraging = 0;
        const uint tAveraging      = 100;
        const uint tStartOutProbe  = 0;
        const uint tOutProbe       = para->getTimestepOut();
        SPtr<PointProbe> pointProbe = std::make_shared<PointProbe>("pointProbe", para->getOutputPath(), tStartAveraging, tAveraging, tStartOutProbe, tOutProbe);
        std::vector<real> probeCoordsX = {0.3, 0.5};
        std::vector<real> probeCoordsY = {0.0, 0.0};
        std::vector<real> probeCoordsZ = {0.0, 0.0};
        pointProbe->addProbePointsFromList(probeCoordsX, probeCoordsY, probeCoordsZ);

        pointProbe->addStatistic(Statistic::Instantaneous);
        pointProbe->addStatistic(Statistic::Means);
        pointProbe->addStatistic(Statistic::Variances);
        para->addProbe( pointProbe );

        SPtr<PlaneProbe> planeProbe = std::make_shared<PlaneProbe>("planeProbe", para->getOutputPath(), tStartAveraging, tAveraging, tStartOutProbe, tOutProbe);
        planeProbe->setProbePlane(dSphere, 0, 0, 0.5, 0.1, 0.1);
        planeProbe->addStatistic(Statistic::Means);
        para->addProbe( planeProbe );

        //////////////////////////////////////////////////////////////////////////
        // setup to copy mesh to simulation
        //////////////////////////////////////////////////////////////////////////
        vf::gpu::Communicator& communicator = vf::gpu::MpiCommunicator::getInstance();
        auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
        SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

        //////////////////////////////////////////////////////////////////////////
        // run simulation
        //////////////////////////////////////////////////////////////////////////

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory, &scalingFactory);
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
