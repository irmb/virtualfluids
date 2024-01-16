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
//! \addtogroup LaminarPlaneFlow
//! \ingroup cpu_apps
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include <iostream>
#include <string>

#include <VirtualFluids.h>

using namespace std;

void run(const vf::basics::ConfigurationFile& config)
{
    using namespace vf::lbm::dir;

    string pathname = config.getValue<string>("pathname");
    int numOfThreads = config.getValue<int>("numOfThreads");
    vector<int> blocknx = config.getVector<int>("blocknx");
    vector<real> boundingBox = config.getVector<real>("boundingBox");
    real endTime = config.getValue<real>("endTime");
    real outTime = config.getValue<real>("outTime");
    int refineLevel = config.getValue<int>("refineLevel");
    real restartStep = config.getValue<real>("restartStep");
    real deltax = config.getValue<real>("deltax");
    real cpStep = config.getValue<real>("cpStep");
    real cpStart = config.getValue<real>("cpStart");
    bool newStart = config.getValue<bool>("newStart");

    SPtr<vf::parallel::Communicator> comm = vf::parallel::MPICommunicator::getInstance();
    int myid = comm->getProcessID();

    real rhoLB1 = 0.00001;
    real rhoLB2 = 0.0;
    real nu = 0.0064;

    real h = boundingBox[1] / 2.0;
    real dp = (rhoLB1 / 3. - rhoLB2 / 3.);
    real L = boundingBox[0];
    real u_max = h * h / (2. * nu) * (dp / L);
    real Re = u_max * 2 * h / nu;

    SPtr<LBMUnitConverter> conv = SPtr<LBMUnitConverter>(new LBMUnitConverter());

    // bounding box
    real g_minX1 = 0;
    real g_minX2 = 0;
    real g_minX3 = 0;

    real g_maxX1 = boundingBox[0];
    real g_maxX2 = boundingBox[1];
    real g_maxX3 = boundingBox[2];

    real blockLength = 3.0 * deltax;

    // bc
    SPtr<BC> noSlipBC(new NoSlipBC());
    noSlipBC->setBCStrategy(SPtr<BCStrategy>(new NoSlipInterpolated()));

    SPtr<BC> pressureBC1(new PressureBC(rhoLB1));
    pressureBC1->setBCStrategy(SPtr<BCStrategy>(new PressureNonEquilibrium()));

    SPtr<BC> pressureBC2(new PressureBC(rhoLB2));
    pressureBC2->setBCStrategy(SPtr<BCStrategy>(new PressureNonEquilibrium()));

    // BS visitor
    BoundaryConditionsBlockVisitor bcVisitor;

    SPtr<Grid3D> grid(new Grid3D(comm));
    grid->setPeriodicX1(false);
    grid->setPeriodicX2(false);
    grid->setPeriodicX3(true);
    grid->setDeltaX(deltax);
    grid->setBlockNX(blocknx[0], blocknx[1], blocknx[2]);

    SPtr<GbObject3D> gridCube(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2, g_maxX3));
    if (myid == 0)
        GbSystem3D::writeGeoObject(gridCube.get(), pathname + "/geo/gridCube", WbWriterVtkXmlBinary::getInstance());

    real k1 = 4;

    SPtr<GbObject3D> refineCube1_1(new GbCuboid3D(g_minX1, g_minX2, g_minX3, g_maxX1, g_maxX2 / k1 - 1.0, g_maxX3));
    if (myid == 0)
        GbSystem3D::writeGeoObject(refineCube1_1.get(), pathname + "/geo/refineCube1_1",
                                   WbWriterVtkXmlBinary::getInstance());

    SPtr<GbObject3D> refineCube1_2(
        new GbCuboid3D(g_minX1, g_maxX2 - g_maxX2 / k1 + 1.0, g_minX3, g_maxX1, g_maxX2, g_maxX3));
    if (myid == 0)
        GbSystem3D::writeGeoObject(refineCube1_2.get(), pathname + "/geo/refineCube1_2",
                                   WbWriterVtkXmlBinary::getInstance());

    SPtr<LBMKernel> kernel;
    kernel = SPtr<LBMKernel>(new K17CompressibleNavierStokes());

    SPtr<BCSet> bcProc;
    bcProc = std::make_shared<BCSet>();
    kernel->setBCSet(bcProc);

    SPtr<Grid3DVisitor> metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, d00M));

    //////////////////////////////////////////////////////////////////////////
    // restart
    SPtr<UbScheduler> mSch(new UbScheduler(cpStep, cpStart));
    SPtr<MPIIOMigrationSimulationObserver> restart(
        new MPIIOMigrationSimulationObserver(grid, mSch, metisVisitor, pathname, comm));
    restart->setLBMKernel(kernel);
    restart->setBCSet(bcProc);
    //////////////////////////////////////////////////////////////////////////

    if (newStart) {
        GenBlocksGridVisitor genBlocks(gridCube);
        grid->accept(genBlocks);

        if (myid == 0) {
            UBLOG(logINFO, "Parameters:");
            UBLOG(logINFO, "h = " << h);
            UBLOG(logINFO, "nue = " << nu);
            UBLOG(logINFO, "Re = " << Re);
            UBLOG(logINFO, "dx = " << deltax);
            UBLOG(logINFO, "dpLB = " << dp);
            UBLOG(logINFO, "Umax = " << u_max);
            UBLOG(logINFO, "number of levels = " << refineLevel + 1);
            UBLOG(logINFO, "numOfThreads = " << numOfThreads);
            UBLOG(logINFO, "path = " << pathname);
            UBLOG(logINFO, "Preprozess - start");
        }

        // walls
        GbCuboid3DPtr addWallYmin(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength,
                                                 g_maxX1 + blockLength, g_minX2, g_maxX3 + blockLength));
        if (myid == 0)
            GbSystem3D::writeGeoObject(addWallYmin.get(), pathname + "/geo/addWallYmin", WbWriterVtkXmlASCII::getInstance());

        GbCuboid3DPtr addWallYmax(new GbCuboid3D(g_minX1 - blockLength, g_maxX2, g_minX3 - blockLength,
                                                 g_maxX1 + blockLength, g_maxX2 + blockLength, g_maxX3 + blockLength));
        if (myid == 0)
            GbSystem3D::writeGeoObject(addWallYmax.get(), pathname + "/geo/addWallYmax", WbWriterVtkXmlASCII::getInstance());

        // inflow
        GbCuboid3DPtr geoInflow(new GbCuboid3D(g_minX1 - blockLength, g_minX2 - blockLength, g_minX3 - blockLength, g_minX1,
                                               g_maxX2 + blockLength, g_maxX3 + blockLength));
        if (myid == 0)
            GbSystem3D::writeGeoObject(geoInflow.get(), pathname + "/geo/geoInflow", WbWriterVtkXmlASCII::getInstance());

        // outflow
        GbCuboid3DPtr geoOutflow(new GbCuboid3D(g_maxX1, g_minX2 - blockLength, g_minX3 - blockLength, g_maxX1 + blockLength,
                                                g_maxX2 + blockLength, g_maxX3 + blockLength));
        if (myid == 0)
            GbSystem3D::writeGeoObject(geoOutflow.get(), pathname + "/geo/geoOutflow", WbWriterVtkXmlASCII::getInstance());

        SPtr<SimulationObserver> ppblocks(new WriteBlocksSimulationObserver(
            grid, SPtr<UbScheduler>(new UbScheduler(1)), pathname, WbWriterVtkXmlBinary::getInstance(), comm));

        if (refineLevel > 0) {
            if (myid == 0)
                UBLOG(logINFO, "Refinement - start");
            RefineCrossAndInsideGbObjectHelper refineHelper(grid, refineLevel, comm);
            refineHelper.addGbObject(refineCube1_1, 1);
            refineHelper.addGbObject(refineCube1_2, 1);
            refineHelper.refine();
            if (myid == 0)
                UBLOG(logINFO, "Refinement - end");
        }

        // walls
        SPtr<D3Q27Interactor> addWallYminInt(new D3Q27Interactor(addWallYmin, grid, noSlipBC, Interactor3D::SOLID));
        SPtr<D3Q27Interactor> addWallYmaxInt(new D3Q27Interactor(addWallYmax, grid, noSlipBC, Interactor3D::SOLID));

        // inflow
        SPtr<D3Q27Interactor> inflowInt =
            SPtr<D3Q27Interactor>(new D3Q27Interactor(geoInflow, grid, pressureBC1, Interactor3D::SOLID));

        // outflow
        SPtr<D3Q27Interactor> outflowInt =
            SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, pressureBC2, Interactor3D::SOLID));
        // SPtr<D3Q27Interactor> outflowSolidInt = SPtr<D3Q27Interactor>(new D3Q27Interactor(geoOutflow, grid, noSlipBC,
        // Interactor3D::SOLID));

        ////////////////////////////////////////////
        // METIS
        SPtr<Grid3DVisitor> metisVisitor(
            new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, d00M));
        ////////////////////////////////////////////
        /////delete solid blocks
        if (myid == 0)
            UBLOG(logINFO, "deleteSolidBlocks - start");
        InteractorsHelper intHelper(grid, metisVisitor);

        intHelper.addInteractor(addWallYminInt);
        intHelper.addInteractor(addWallYmaxInt);

        intHelper.addInteractor(inflowInt);

        intHelper.addInteractor(outflowInt);

        intHelper.selectBlocks();
        if (myid == 0)
            UBLOG(logINFO, "deleteSolidBlocks - end");
        //////////////////////////////////////

        ppblocks->update(0);
        ppblocks.reset();

        if (myid == 0)
            VF_LOG_INFO("{}", Utilities::toString(grid, comm->getNumberOfProcesses()));

        SetKernelBlockVisitor kernelVisitor(kernel, nu);
        grid->accept(kernelVisitor);

        if (refineLevel > 0) {
            SetUndefinedNodesBlockVisitor undefNodesVisitor;
            grid->accept(undefNodesVisitor);
        }

        // walls
        intHelper.setBC();

        grid->accept(bcVisitor);

        InitDistributionsBlockVisitor initVisitor;
        grid->accept(initVisitor);

        // Postprocess
        SPtr<UbScheduler> geoSch(new UbScheduler(1));
        SPtr<SimulationObserver> ppgeo(new WriteBoundaryConditionsSimulationObserver(
            grid, geoSch, pathname, WbWriterVtkXmlBinary::getInstance(), comm));
        ppgeo->update(0);
        ppgeo.reset();

        if (myid == 0)
            UBLOG(logINFO, "Preprocess - end");
    } else {
        restart->restart((int)restartStep);
        grid->setTimeStep(restartStep);

        if (myid == 0)
            UBLOG(logINFO, "Restart - end");
    }

    grid->accept(bcVisitor);

    // set connectors
    OneDistributionSetConnectorsBlockVisitor setConnsVisitor(comm);
    grid->accept(setConnsVisitor);

    SPtr<Interpolator> iProcessor(new CompressibleOffsetMomentsInterpolator());
    SetInterpolationConnectorsBlockVisitor setInterConnsVisitor(comm, nu, iProcessor);
    grid->accept(setInterConnsVisitor);

    SPtr<UbScheduler> nupsSch(new UbScheduler(10, 30, 100));
    SPtr<SimulationObserver> npr(new NUPSCounterSimulationObserver(grid, nupsSch, numOfThreads, comm));

    // write data for visualization of macroscopic quantities
    SPtr<UbScheduler> visSch(new UbScheduler(outTime));
    SPtr<WriteMacroscopicQuantitiesSimulationObserver> writeMQSimulationObserver(
        new WriteMacroscopicQuantitiesSimulationObserver(grid, visSch, pathname, WbWriterVtkXmlASCII::getInstance(),
                                                         SPtr<LBMUnitConverter>(new LBMUnitConverter()), comm));

    // OpenMP threads control
#ifdef _OPENMP
    omp_set_num_threads(numOfThreads);
#endif
    // Create simulation
    SPtr<Simulation> simulation(new Simulation(grid, visSch, endTime));
    simulation->addSimulationObserver(npr);
    simulation->addSimulationObserver(writeMQSimulationObserver);

    //////////////////////////////////////////////////////////////////////////
    // Run simulation
    //////////////////////////////////////////////////////////////////////////

    UBLOG(logINFO, "Simulation-start");
    simulation->run();
    UBLOG(logINFO, "Simulation-end");
}

int main(int argc, char* argv[])
{
    try {
        vf::logging::Logger::initializeLogger();
        auto config = vf::basics::loadConfig(argc, argv);
        run(config);
    } catch (const std::exception& e) {
        VF_LOG_WARNING("{}", e.what());
        return 1;
    }
    return 0;
}

//! \}
