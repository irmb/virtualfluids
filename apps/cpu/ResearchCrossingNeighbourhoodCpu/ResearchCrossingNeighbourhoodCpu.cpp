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
//! \addtogroup City
//! \ingroup cpu_apps
//! \{
//! \author HAH
//=======================================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <basics/utilities/UbTuple.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "VirtualFluids.h"

using namespace std;

int main(int argc, char* argv[])
{
    using namespace vf::lbm::dir;
    try {
        vf::logging::Logger::initializeLogger();
        //////////////////////////////////////////////////////////////////////////
        // Simulation parameters
        //////////////////////////////////////////////////////////////////////////
        vf::basics::ConfigurationFile config = vf::basics::loadConfig(argc, argv);

        // set your output path here
        string path = config.getValue<string>("path");
        //const real L = config.getValue<real>("L");
        const real Re = config.getValue<real>("Re");
        const real velocity = config.getValue<real>("velocity");
        //const real dt = config.getValue<real>("dt");
        //const unsigned int nx = config.getValue<unsigned int>("nx");

        const real timeStepOut = config.getValue<real>("timeStepOut");
        const real timeStepEnd = config.getValue<real>("timeStepEnd");
        string CityFilename     = config.getValue<string>("CityFilename");
        string BodenFilename    = config.getValue<string>("BodenFilename");
        const real dx          = config.getValue<real>("deltaX");
        const real fastWindingThreshold = config.getValue<real>("FastWindingThreshold");

        // Number of OpenMP threads
        int numOfThreads = config.getValue<int>("numOfThreads");

         //////////////////////////////////////////////////////////////////////////
        //Geometry
        SPtr<GbTriFaceMesh3D> city     = std::make_shared<GbTriFaceMesh3D>();
        SPtr<GbTriFaceMesh3D> boden    = std::make_shared<GbTriFaceMesh3D>();

        bool removeRedundant = true;
        city->readMeshFromSTLFileBinary(CityFilename, removeRedundant);
        boden->readMeshFromSTLFileBinary(BodenFilename, removeRedundant);

        real initialX = city->getX1Minimum();
        real initialY = city->getX2Minimum();
        real initialZ = city->getX3Minimum();
        city->translate(-initialX, -initialY, -initialZ);
        boden->translate(-initialX, -initialY, -initialZ);
        gb_system_3d::writeGeoObject(city.get(), path + "/geo/city", WbWriterVtkXmlBinary::getInstance());
        gb_system_3d::writeGeoObject(boden.get(), path + "/geo/boden", WbWriterVtkXmlBinary::getInstance());

        //////////////////////////////////////////////////////////////////////////
        // create grid
        //////////////////////////////////////////////////////////////////////////
        const std::vector<real> gridMin = config.getVector<real>("GridCubeMin");
        const std::vector<real> gridMax = config.getVector<real>("GridCubeMax");
        if (gridMin.size() != 3 || gridMax.size() != 3) {
            throw std::runtime_error("GridCubeMin and GridCubeMax must have 3 entries each.");
        }
        const std::vector<int> blockNx = config.getVector<int>("blocknx");
        if (blockNx.size() != 3) {
            throw std::runtime_error("blocknx must have 3 entries.");
        }

        real g_minX1 = gridMin[0];
        real g_minX2 = gridMin[1];
        real g_minX3 = gridMin[2];

        real g_maxX1 = gridMax[0];
        real g_maxX2 = gridMax[1];
        real g_maxX3 = gridMax[2];

        real boundingBoxSizeX = g_maxX1 - g_minX1;
        real boundingBoxSizeY = g_maxX2 - g_minX2;
        real boundingBoxSizeZ = g_maxX3 - g_minX3;

        //////////////////////////////////////////////////////////////////////////
        real L  = boundingBoxSizeX;
        real nx = L/dx;
        const real velocityLB = velocity;    // LB units
        const real viscosityLB = nx * velocityLB / Re; // LB unit

        SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
        int myid = comm->getProcessID();

        // new grid object
        SPtr<Grid3D> grid(new Grid3D(comm));
        // set grid spacing
        grid->setDeltaX(dx);
        // match GPU setup: periodic in X and Y, non-periodic in Z
        grid->setPeriodicX1(true);
        grid->setPeriodicX2(true);
        grid->setPeriodicX3(false);
        // set block size for three dimensions (ensure integer number of blocks)
        const int nodesX = static_cast<int>(std::lround(boundingBoxSizeX / dx));
        const int nodesY = static_cast<int>(std::lround(boundingBoxSizeY / dx));
        const int nodesZ = static_cast<int>(std::lround(boundingBoxSizeZ / dx));

        auto validateBlock = [](int nodes, int block, const char* axis) {
            if (nodes <= 0) {
                throw std::runtime_error("Invalid grid extent along axis " + std::string(axis) + ".");
            }
            if (block <= 0) {
                throw std::runtime_error("blocknx values must be positive.");
            }
            if (nodes % block != 0) {
                throw std::runtime_error("blocknx does not evenly divide grid nodes along axis " +
                                         std::string(axis) + ".");
            }
        };

        const int blockSizeX = blockNx[0];
        const int blockSizeY = blockNx[1];
        const int blockSizeZ = blockNx[2];
        validateBlock(nodesX, blockSizeX, "X");
        validateBlock(nodesY, blockSizeY, "Y");
        validateBlock(nodesZ, blockSizeZ, "Z");
        grid->setBlockNX(blockSizeX, blockSizeY, blockSizeZ);

        // Create simulation bounding box
        SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
        gb_system_3d::writeGeoObject(gridCube.get(), path + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

        if (myid == comm->getRoot()) {
            UBLOG(logINFO, "Research Crossing Neighbourhood:");
            UBLOG(logINFO, "Domain size = " << boundingBoxSizeX / dx << " x " << boundingBoxSizeY / dx << " x "
                                             << boundingBoxSizeZ / dx);
            UBLOG(logINFO, "Block size = " << blockSizeX << " x " << blockSizeY << " x " << blockSizeZ);
            UBLOG(logINFO, "velocity    = " << velocity << " m/s");
            UBLOG(logINFO, "velocityLB  = " << velocityLB);
            UBLOG(logINFO, "viscosityLB = " << viscosityLB);
            UBLOG(logINFO, "u  = " << velocityLB);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << dx);
            //UBLOG(logINFO, "dt = " << dt);
            UBLOG(logINFO, "City Braunschweig:"
                               << " X:[" << city->getX1Minimum() << ", " << city->getX1Maximum() << "]"
                               << " Y:[" << city->getX2Minimum() << ", " << city->getX2Maximum() << "]"
                               << " Z:[" << city->getX3Minimum() << ", " << city->getX3Maximum() << "]");
            UBLOG(logINFO, "Preprocess - start");
        }

        // Generate block grid
        GenBlocksGridVisitor genBlocks(gridCube);
        grid->accept(genBlocks);



        // Create LBM kernel
        auto kernel = std::make_shared<K17CompressibleNavierStokes>();
       
        //////////////////////////////////////////////////////////////////////////
        // Create boundary conditions (BC)
        //////////////////////////////////////////////////////////////////////////
        // Create no-slip BC
        auto noSlipBC = std::make_shared<NoSlipBC>();
        noSlipBC->setBCStrategy(std::make_shared<NoSlipInterpolated>());

        // Velocity BC
        mu::Parser fct;
        fct.SetExpr("u");
        fct.DefineConst("u", velocityLB);
        // Set the same velocity in x and y-direction
        auto velBC = std::make_shared<VelocityBC>(true, true, false, fct, 0, BCFunction::INFCONST);
        velBC->setBCStrategy(std::make_shared<VelocityInterpolated>());

        // Add velocity boundary condition to visitor. No-slip boundary
        BoundaryConditionsBlockVisitor bcVisitor;

        // Create boundary conditions
        SPtr<BCSet> bcProc;
        bcProc = std::make_shared<BCSet>();
        kernel->setBCSet(bcProc);

        // Create boundary conditions geometry (only top Z boundary is a wall; bottom is given by geometry)
        GbCuboid3DPtr wallZmax(
            new GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx));
        gb_system_3d::writeGeoObject(wallZmax.get(), path + "/geo/wallZmax", WbWriterVtkXmlASCII::getInstance());

        // Add boundary conditions to grid generator
        SPtr<D3Q27Interactor> wallZmaxInt(new D3Q27Interactor(wallZmax, grid, velBC, Interactor3D::SOLID));
        SPtr<D3Q27GridWindingInteractor> CityInt(new D3Q27GridWindingInteractor(
            city, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::POINTS));
        // Ground (boden) obstacle, classified the same way as the city
        SPtr<D3Q27GridWindingInteractor> BodenInt(new D3Q27GridWindingInteractor(
            boden, grid, noSlipBC, Interactor3D::SOLID, Interactor3D::POINTS));

        SPtr<Grid3DVisitor> metisVisitor(
            new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, d00M));
        InteractorsHelper intHelper(grid, metisVisitor);
        intHelper.addInteractor(wallZmaxInt);
        intHelper.addInteractor(CityInt);
        intHelper.addInteractor(BodenInt);

        intHelper.selectBlocks();

        if (myid == 0)
            VF_LOG_INFO("{}", utilities::toString(grid, comm->getNumberOfProcesses()));

        // Generate grid
        SetKernelBlockVisitor kernelVisitor(kernel, viscosityLB);
        grid->accept(kernelVisitor);

        // Apply interactors so BC/solid data is present before BC output.
        intHelper.setBC();

        // Write block grid to VTK-file
        auto ppblocks = std::make_shared<WriteBlocksSimulationObserver>(grid, SPtr<UbScheduler>(new UbScheduler(1)),
                                                                        path, WbWriterVtkXmlBinary::getInstance(), comm);
        ppblocks->update(0);
        ppblocks.reset();

        // Initialization of distributions
        InitDistributionsBlockVisitor initVisitor;
        initVisitor.setRho(vf::basics::constant::c0o1);
        initVisitor.setVx1(static_cast<real>(0.1));
        initVisitor.setVx2(static_cast<real>(0.1));
        initVisitor.setVx3(static_cast<real>(0.0));
        grid->accept(initVisitor);

        // Set connectivity between blocks
        OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
        grid->accept(setConnsVisitor);

        // Create lists of boundary nodes
        grid->accept(bcVisitor);

        // Write grid with boundary conditions information to VTK-file
        SPtr<UbScheduler> geoSch(new UbScheduler(1));
        WriteBoundaryConditionsSimulationObserver ppgeo(grid, geoSch, path, WbWriterVtkXmlBinary::getInstance(), comm);
        ppgeo.update(0);

        // Create coprocessor object for writing macroscopic quantities to VTK-file
        SPtr<UbScheduler> visSch(new UbScheduler(timeStepOut));
        SPtr<SimulationObserver> mqSimulationObserver(new WriteMacroscopicQuantitiesSimulationObserver(
            grid, visSch, path, WbWriterVtkXmlBinary::getInstance(),
            SPtr<LBMUnitConverter>(new LBMUnitConverter(L, velocity, 1.0, nx, velocityLB)), comm));
        mqSimulationObserver->update(0);

        // Grid-winding based boundary diagnostics for the city geometry
        // (subgrid Q-values, missing links).
#if defined(VF_HAS_FAST_WINDING)
        vf::grid_winding::cpu::BoundaryProcessingConfig boundaryConfig;
        boundaryConfig.fastWindingThreshold = static_cast<float>(fastWindingThreshold);

        std::vector<SPtr<BC>> windingAdapters{ noSlipBC };

        SPtr<UbScheduler> gwSch(new UbScheduler(1));
        auto gwObserver = std::make_shared<GridWindingDiagnosticsSimulationObserver>(
            grid,
            gwSch,
            std::vector<SPtr<GbTriFaceMesh3D>>{ city, boden },
            bcProc,
            windingAdapters,
            nullptr,
            path,
            comm,
            boundaryConfig);

        gwObserver->update(0);
        gwObserver.reset();
#endif

        // Create coprocessor object for writing NUPS
        SPtr<UbScheduler> nupsSch(new UbScheduler(100, 100));
        SPtr<SimulationObserver> nupsSimulationObserver(
            new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

        // OpenMP threads control
#ifdef _OPENMP
        omp_set_num_threads(numOfThreads);
#endif
        // Create simulation
        SPtr<Simulation> simulation(new Simulation(grid, visSch, timeStepEnd));
        simulation->addSimulationObserver(nupsSimulationObserver);
        simulation->addSimulationObserver(mqSimulationObserver);

        //////////////////////////////////////////////////////////////////////////
        // Run simulation
        //////////////////////////////////////////////////////////////////////////
        if (myid == comm->getRoot()) {
            UBLOG(logINFO, "Preprocess - end");
            UBLOG(logINFO, "Total Physical Memory (RAM): " << utilities::getTotalPhysMem() / 1e9 << " GB");
            UBLOG(logINFO, "Physical Memory currently used: " << utilities::getPhysMemUsed() / 1e9 << " GB");
            //UBLOG(logINFO, "Physical Memory currently used by current process: " << utilities::getPhysMemUsedByMe() / 1e9 << " GB");
            UBLOG(logINFO, "Simulation-start");
        }
        simulation->run();
        if (myid == comm->getRoot()) {
            UBLOG(logINFO, "Simulation-end");
        }
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}

//! \}
