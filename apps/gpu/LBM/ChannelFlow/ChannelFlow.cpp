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
//! \file ChannelFlow.cpp
//! \ingroup Applications
//! \author Anna Wellmann
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

#include "mpi.h"

//////////////////////////////////////////////////////////////////////////

#include "Core/DataTypes.h"
#include "Core/LbmOrGks.h"
#include "Core/Logger/Logger.h"
#include "Core/VectorTypes.h"
#include "PointerDefinitions.h"
#include "config/ConfigurationFile.h"
#include "logger/Logger.h"

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/GridFactory.h"

#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/utilities/communication.h"

//////////////////////////////////////////////////////////////////////////

#include "VirtualFluids_GPU/BoundaryConditions/BoundaryConditionFactory.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "VirtualFluids_GPU/GPU/CudaMemoryManager.h"
#include "VirtualFluids_GPU/LBM/Simulation.h"
#include "VirtualFluids_GPU/Output/FileWriter.h"
#include "VirtualFluids_GPU/Parameter/Parameter.h"

//////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    try {
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////

        const real channelWidth = 1.0;
        const real Re = 10000.0;
        const uint nx = 700; // ~60 GB on A100 single precision
        const real velocityLB = 0.05; // LB units

        const uint timeStepOut = 10000;
        const uint timeStepEnd = 100000;

        //////////////////////////////////////////////////////////////////////////
        // setup simulation parameters (without config file)
        //////////////////////////////////////////////////////////////////////////

        vf::gpu::Communicator &communicator = vf::gpu::Communicator::getInstance();
        SPtr<Parameter> para = std::make_shared<Parameter>(communicator.getNummberOfProcess(), communicator.getPID());
        para->setDevices({ 0, 1, 2, 3, 4 });
        para->setMaxDev(5);
        BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

        //////////////////////////////////////////////////////////////////////////
        // setup logger
        //////////////////////////////////////////////////////////////////////////

        std::ofstream logFile("output/log_process" + std::to_string(vf::gpu::Communicator::getInstance().getPID()) +
                              ".txt");
        logging::Logger::addStream(&logFile);
        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
        logging::Logger::timeStamp(logging::Logger::ENABLE);
        logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

        vf::logging::Logger::changeLogPath("output/vflog_process" +
                                           std::to_string(vf::gpu::Communicator::getInstance().getPID()) + ".txt");
        vf::logging::Logger::initalizeLogger();

        //////////////////////////////////////////////////////////////////////////
        // setup gridGenerator
        //////////////////////////////////////////////////////////////////////////

        auto gridFactory = GridFactory::make();
        gridFactory->setTriangularMeshDiscretizationMethod(TriangularMeshDiscretizationMethod::POINT_IN_OBJECT);
        auto gridBuilder = MultipleGridBuilder::makeShared(gridFactory);

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////

        const real xGridMin = 0.0 * channelWidth;
        // const real xGridMax = 5.0 * channelWidth;
        const real yGridMin = 0.0 * channelWidth;
        const real yGridMax = 1.0 * channelWidth;
        const real zGridMin = 0.0 * channelWidth;
        const real zGridMax = 1.0 * channelWidth;

        real dx = channelWidth / real(nx);

        //////////////////////////////////////////////////////////////////////////
        // compute parameters in lattice units
        //////////////////////////////////////////////////////////////////////////

        const real viscosityLB = (channelWidth / dx) * velocityLB / Re; // LB units

        VF_LOG_INFO("LB parameters:");
        VF_LOG_INFO("velocity LB [dx/dt]              = {}", velocityLB);
        VF_LOG_INFO("viscosity LB [dx/dt]             = {}", viscosityLB);

        //////////////////////////////////////////////////////////////////////////
        // set parameters
        //////////////////////////////////////////////////////////////////////////

        para->setPrintFiles(true);

        para->setVelocityLB(velocityLB);
        para->setViscosityLB(viscosityLB);

        para->setVelocityRatio((real)1.0);
        para->setDensityRatio((real)1.0);

        para->setTimestepOut(timeStepOut);
        para->setTimestepEnd(timeStepEnd);

        para->setOutputPrefix("ChannelFlow");
        para->setMainKernel("CumulantK17CompChimStream");

        const uint generatePart = vf::gpu::Communicator::getInstance().getPID();
        real overlap = (real)8.0 * dx;

        if (communicator.getNummberOfProcess() == 5) {

            //////////////////////////////////////////////////////////////////////////
            // add coarse grids
            //////////////////////////////////////////////////////////////////////////

            if (generatePart == 0) {
                gridBuilder->addCoarseGrid(xGridMin, yGridMin, zGridMin, 1.0 * channelWidth + overlap, yGridMax,
                                           zGridMax, dx);
            }
            if (generatePart == 1) {
                gridBuilder->addCoarseGrid(1.0 * channelWidth - overlap, yGridMin, zGridMin,
                                           2.0 * channelWidth + overlap, yGridMax, zGridMax, dx);
            }
            if (generatePart == 2) {
                gridBuilder->addCoarseGrid(2.0 * channelWidth - overlap, yGridMin, zGridMin,
                                           3.0 * channelWidth + overlap, yGridMax, zGridMax, dx);
            }
            if (generatePart == 3) {
                gridBuilder->addCoarseGrid(3.0 * channelWidth - overlap, yGridMin, zGridMin,
                                           4.0 * channelWidth + overlap, yGridMax, zGridMax, dx);
            }
            if (generatePart == 4) {
                gridBuilder->addCoarseGrid(4.0 * channelWidth - overlap, yGridMin, zGridMin, 5.0 * channelWidth,
                                           yGridMax, zGridMax, dx);
            }

            //////////////////////////////////////////////////////////////////////////
            // set subdomain dimensions
            //////////////////////////////////////////////////////////////////////////

            if (generatePart == 0) {
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(xGridMin, 1.0 * channelWidth, yGridMin,
                                                                           yGridMax, zGridMin, zGridMax));
            }
            if (generatePart == 1) {
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(1.0 * channelWidth, 2.0 * channelWidth,
                                                                           yGridMin, yGridMax, zGridMin, zGridMax));
            }
            if (generatePart == 2) {
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(2.0 * channelWidth, 3.0 * channelWidth,
                                                                           yGridMin, yGridMax, zGridMin, zGridMax));
            }
            if (generatePart == 3) {
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(3.0 * channelWidth, 4.0 * channelWidth,
                                                                           yGridMin, yGridMax, zGridMin, zGridMax));
            }
            if (generatePart == 4) {
                gridBuilder->setSubDomainBox(std::make_shared<BoundingBox>(4.0 * channelWidth, 5.0 * channelWidth,
                                                                           yGridMin, yGridMax, zGridMin, zGridMax));
            }

            //////////////////////////////////////////////////////////////////////////
            // build grids
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->buildGrids(LBM, true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////
            // configure communication neighbors
            //////////////////////////////////////////////////////////////////////////

            if (generatePart == 0) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 1);
            }
            if (generatePart == 1) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 0);
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 2);
            }
            if (generatePart == 2) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 1);
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 3);
            }
            if (generatePart == 3) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 2);
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, 4);
            }
            if (generatePart == 4) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX, LBM);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, 3);
            }

            //////////////////////////////////////////////////////////////////////////
            // set boundary conditions
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            if (generatePart == 0) {
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
            }
            if (generatePart == 4) {
                gridBuilder->setPressureBoundaryCondition(SideType::PX,
                                                          0.0); // set pressure boundary condition last
                bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
            }
            gridBuilder->setNoSlipBoundaryCondition(SideType::MY);
            gridBuilder->setNoSlipBoundaryCondition(SideType::PY);
            gridBuilder->setNoSlipBoundaryCondition(SideType::MZ);
            gridBuilder->setNoSlipBoundaryCondition(SideType::PZ);
            bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
            bcFactory.setVelocityBoundaryCondition(
                BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
        } else {
            VF_LOG_CRITICAL("This app has no setup for {} GPUs.", para->getNumprocs());
        }

        //////////////////////////////////////////////////////////////////////////
        // setup to copy mesh to simulation
        //////////////////////////////////////////////////////////////////////////

        auto cudaMemoryManager = std::make_shared<CudaMemoryManager>(para);
        SPtr<GridProvider> gridGenerator =
            GridProvider::makeGridGenerator(gridBuilder, para, cudaMemoryManager, communicator);

        //////////////////////////////////////////////////////////////////////////
        // run simulation
        //////////////////////////////////////////////////////////////////////////

        Simulation sim(para, cudaMemoryManager, communicator, *gridGenerator, &bcFactory);
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
