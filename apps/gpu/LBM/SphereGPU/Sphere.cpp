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

#include "Core/DataTypes.h"
#include "Core/LbmOrGks.h"
#include "Core/Logger/Logger.h"
#include "Core/VectorTypes.h"
#include "PointerDefinitions.h"
#include "config/ConfigurationFile.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
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
        const real Re = 1000.0; // related to the sphere's diameter
        const real velocity = 1.0;
        const real dt = (real)0.5e-3;
        const uint nx = 64;

        const uint timeStepOut = 1000;
        const uint timeStepEnd = 10000;

        //////////////////////////////////////////////////////////////////////////
        // setup logger
        //////////////////////////////////////////////////////////////////////////

        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
        logging::Logger::timeStamp(logging::Logger::ENABLE);
        logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

        //////////////////////////////////////////////////////////////////////////
        // setup simulation parameters (with or without config file)
        //////////////////////////////////////////////////////////////////////////

        vf::gpu::Communicator& communicator = vf::gpu::Communicator::getInstance();;
        SPtr<Parameter> para;
        BoundaryConditionFactory bcFactory = BoundaryConditionFactory();
        vf::basics::ConfigurationFile config;
        if (useConfigFile) {
            //////////////////////////////////////////////////////////////////////////
            // read simulation parameters from config file
            //////////////////////////////////////////////////////////////////////////

            // assuming that a config files is stored parallel to this file.
            std::filesystem::path configPath = __FILE__;

            // the config file's default name can be replaced by passing a command line argument
            std::string configName("config.txt");
            if (argc == 2) {
                configName = argv[1];
                std::cout << "Using configFile command line argument: " << configName << std::endl;
            }

            configPath.replace_filename(configName);
            config.load(configPath.string());

            para = std::make_shared<Parameter>(config);
        } else {
            para = std::make_shared<Parameter>();
        }

        //////////////////////////////////////////////////////////////////////////
        // setup gridGenerator
        //////////////////////////////////////////////////////////////////////////

        auto gridFactory = GridFactory::make();
        gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
        auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////

        real dx = L / real(nx);
        gridBuilder->addCoarseGrid(-1.0 * L, -0.8 * L, -0.8 * L,
                                    6.0 * L,  0.8 * L,  0.8 * L, dx);

        // use primitive
        Object *sphere = new Sphere(0.0, 0.0, 0.0, dSphere / 2.0);

        // use stl
        // std::string stlPath = "stl/sphere02.stl";
        // if (useConfigFile && config.contains("STLPath")) {
        //     stlPath = config.getValue<std::string>("STLPath");
        // }
        // std::cout << "Reading stl from " << stlPath << "." << std::endl;
        // Object *sphere = TriangularMesh::make(stlPath);

        gridBuilder->addGeometry(sphere);
        gridBuilder->setPeriodicBoundaryCondition(false, false, false);
        gridBuilder->buildGrids(LBM, false);  // buildGrids() has to be called before setting the BCs!!!!

        //////////////////////////////////////////////////////////////////////////
        // compute parameters in lattice units
        //////////////////////////////////////////////////////////////////////////

        const real velocityLB = velocity * dt / dx; // LB units
        const real viscosityLB =  (dSphere / dx) * velocityLB / Re; // LB units

        *logging::out << logging::Logger::INFO_HIGH << "velocity  [dx/dt] = " << velocityLB << " \n";
        *logging::out << logging::Logger::INFO_HIGH << "viscosity [dx^2/dt] = " << viscosityLB << "\n";

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
        SPtr<PointProbe> pointProbe = std::make_shared<PointProbe>( "pointProbe", para->getOutputPath(), tStartAveraging, tAveraging, tStartOutProbe, tOutProbe);
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

        auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
        SPtr<GridProvider> gridGenerator = GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

        //////////////////////////////////////////////////////////////////////////
        // run simulation
        //////////////////////////////////////////////////////////////////////////

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
        sim.run();

    } catch (const std::bad_alloc &e) {
        *logging::out << logging::Logger::LOGGER_ERROR << "Bad Alloc:" << e.what() << "\n";
    } catch (const std::exception &e) {
        *logging::out << logging::Logger::LOGGER_ERROR << e.what() << "\n";
    } catch (std::string &s) {
        *logging::out << logging::Logger::LOGGER_ERROR << s << "\n";
    } catch (...) {
        *logging::out << logging::Logger::LOGGER_ERROR << "Unknown exception!\n";
    }

    return 0;
}
