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
#include <numeric>
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

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <basics/StringUtilities/StringUtil.h>
#include <basics/config/ConfigurationFile.h>

//////////////////////////////////////////////////////////////////////////

#include "GridGenerator/grid/BoundaryConditions/Side.h"
#include "GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

#include "GridGenerator/geometries/Sphere/Sphere.h"
#include "GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "GridGenerator/utilities/communication.h"

//////////////////////////////////////////////////////////////////////////

#include "gpu/core/BoundaryConditions/BoundaryConditionFactory.h"
#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h"
#include "gpu/core/GPU/CudaMemoryManager.h"
#include "gpu/core/LBM/Simulation.h"
#include "gpu/core/Output/FileWriter.h"
#include "gpu/core/Parameter/Parameter.h"
#include "gpu/core/Kernel/Utilities/KernelTypes.h"

//////////////////////////////////////////////////////////////////////////

#include <parallel/MPICommunicator.h>

int main(int argc, char *argv[])
{
    try {
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////

        const real channelWidth = 1.0;
        const real Re = 10000.0;
        const uint nx = 700;          // 700 nodes need ~60 GB on A100 (single precision)
        const real velocityLB = 0.05; // LB units

        const uint timeStepOut = 10000;
        const uint timeStepEnd = 100000;

        //////////////////////////////////////////////////////////////////////////
        // setup simulation parameters (without config file)
        //////////////////////////////////////////////////////////////////////////

        vf::parallel::Communicator &communicator = *vf::parallel::MPICommunicator::getInstance();
        const int numberOfProcesses = communicator.getNumberOfProcesses();
        const auto processID = communicator.getProcessID();
        SPtr<Parameter> para = std::make_shared<Parameter>(numberOfProcesses, processID);
        std::vector<uint> devices(10);
        std::iota(devices.begin(), devices.end(), 0);
        para->setDevices(devices);
        para->setMaxDev(communicator.getNumberOfProcesses());
        BoundaryConditionFactory bcFactory = BoundaryConditionFactory();

        //////////////////////////////////////////////////////////////////////////
        // setup logger
        //////////////////////////////////////////////////////////////////////////
        vf::logging::Logger::changeLogPath("output/vflog_process" +
                                           std::to_string(processID) + ".txt");
        vf::logging::Logger::initializeLogger();

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////

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
        para->configureMainKernel(vf::collisionKernel::compressible::K17CompressibleNavierStokes);

        real overlap = (real)8.0 * dx;

        if (numberOfProcesses > 1) {

            //////////////////////////////////////////////////////////////////////////
            // add coarse grids
            //////////////////////////////////////////////////////////////////////////

            real subdomainMinX = channelWidth * processID;
            real subdomainMinXoverlap = subdomainMinX;
            real subdomainMaxX = subdomainMinX + channelWidth;
            real subdomainMaxXoverlap = subdomainMaxX;

            if (processID != 0)
                subdomainMinXoverlap -= overlap;

            if (processID != numberOfProcesses - 1)
                subdomainMaxXoverlap += overlap;

            auto gridBuilder = std::make_shared<MultipleGridBuilder>();

            gridBuilder->addCoarseGrid(subdomainMinXoverlap, yGridMin, zGridMin, subdomainMaxXoverlap, yGridMax,
                                       zGridMax, dx);

            //////////////////////////////////////////////////////////////////////////
            // set subdomain dimensions
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setSubDomainBox(
                std::make_shared<BoundingBox>(subdomainMinX, subdomainMaxX, yGridMin, yGridMax, zGridMin, zGridMax));

            //////////////////////////////////////////////////////////////////////////
            // build grids
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->buildGrids(true); // buildGrids() has to be called before setting the BCs!!!!

            //////////////////////////////////////////////////////////////////////////
            // configure communication neighbors
            //////////////////////////////////////////////////////////////////////////

            if (processID != 0) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::MX);
                gridBuilder->setCommunicationProcess(CommunicationDirections::MX, processID - 1);
            }

            if (processID != numberOfProcesses - 1) {
                gridBuilder->findCommunicationIndices(CommunicationDirections::PX);
                gridBuilder->setCommunicationProcess(CommunicationDirections::PX, processID + 1);
            }

            //////////////////////////////////////////////////////////////////////////
            // set boundary conditions
            //////////////////////////////////////////////////////////////////////////

            gridBuilder->setPeriodicBoundaryCondition(false, false, false);

            if (processID == 0) {
                gridBuilder->setVelocityBoundaryCondition(SideType::MX, velocityLB, 0.0, 0.0);
            }
            if (processID == numberOfProcesses - 1) {
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
            VF_LOG_CRITICAL("This app has no setup for a single GPU");
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
