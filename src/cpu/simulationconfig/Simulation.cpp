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
#include <simulationconfig/Simulation.h>

#include <memory>
#include <string>
#include <set>
#include <utility>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <mpi.h>

#include <basics/utilities/UbScheduler.h>
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbSystem3D.h>

#include <BoundaryConditions/BCStrategy.h>
#include <BoundaryConditions/BCSet.h>
#include <SimulationObservers/NUPSCounterSimulationObserver.h>
#include <SimulationObservers/WriteBlocksSimulationObserver.h>
#include <SimulationObservers/WriteBoundaryConditionsSimulationObserver.h>
#include <SimulationObservers/WriteMacroscopicQuantitiesSimulationObserver.h>

#include <Simulation/Simulation.h>
#include <Simulation/Grid3D.h>
#include <Interactors/InteractorsHelper.h>
#include <LBM/Interpolation/CompressibleOffsetMomentsInterpolator.h>
#include <LBM/LBMKernel.h>
#include <LBM/LBMUnitConverter.h>
#include <parallel/MPICommunicator.h>
#include <Visitors/GenBlocksGridVisitor.h>
#include <Visitors/InitDistributionsBlockVisitor.h>
#include <Visitors/MetisPartitioningGridVisitor.h>
#include <Visitors/SetConnectorsBlockVisitor.h>
#include <Visitors/SetKernelBlockVisitor.h>

#include <simulationconfig/SimulationParameters.h>

#include <lbm/constants/D3Q27.h>

#include <logger/Logger.h>


CPUSimulation::CPUSimulation()
{
    this->communicator = vf::parallel::MPICommunicator::getInstance();
    this->grid = std::make_shared<Grid3D>(communicator);
}

void CPUSimulation::setGridParameters(std::shared_ptr<GridParameters> parameters)
{
    this->gridParameters = std::move(parameters);
}

void CPUSimulation::setPhysicalParameters(std::shared_ptr<PhysicalParameters> parameters)
{
    this->physicalParameters = std::move(parameters);
}

void CPUSimulation::setRuntimeParameters(std::shared_ptr<RuntimeParameters> parameters)
{
    this->simulationParameters = std::move(parameters);
}

void CPUSimulation::addObject(const std::shared_ptr<GbObject3D> &object, const std::shared_ptr<BC> &bcAdapter, int state, const std::string &folderPath)
{
    const bool is_in = registeredAdapters.find(bcAdapter) != registeredAdapters.end();
    if (!is_in) addBCAdapter(bcAdapter);
    this->interactors.push_back(lbmSystem->makeInteractor(object, this->grid, bcAdapter, state));
    if (communicator->getProcessID() != 0) return;

    GbSystem3D::writeGeoObject(object, writerConfig.outputPath + folderPath, writerConfig.getWriter());
}

void CPUSimulation::addBCAdapter(const std::shared_ptr<BC> &bcAdapter)
{
    registeredAdapters.insert(bcAdapter);
}

void CPUSimulation::setKernelConfiguration(const std::shared_ptr<LBMKernelConfiguration> &kernel)
{
    this->kernelConfig = kernel;
    this->lbmKernel = kernelFactory.makeKernel(kernel->kernelType);
    this->lbmSystem = kernelFactory.makeLBMSystem(kernel->kernelType);
}

void CPUSimulation::setWriterConfiguration(WriterConfiguration config)
{
    this->writerConfig = config;
}

void CPUSimulation::run()
{
    VF_LOG_INFO("Beginning simulation setup for MPI rank {}", communicator->getProcessID());
    grid->setDeltaX(gridParameters->nodeDistance);
    grid->setPeriodicX1(gridParameters->periodicBoundaryInX1);
    grid->setPeriodicX2(gridParameters->periodicBoundaryInX2);
    grid->setPeriodicX3(gridParameters->periodicBoundaryInX3);

    std::shared_ptr<LBMUnitConverter> converter = makeLBMUnitConverter();

    int &nodesInX1 = gridParameters->numberOfNodesPerDirection[0];
    int &nodesInX2 = gridParameters->numberOfNodesPerDirection[1];
    int &nodesInX3 = gridParameters->numberOfNodesPerDirection[2];

    if (isMainProcess())
        logSimulationData(nodesInX1, nodesInX2, nodesInX3);

    setBlockSize(nodesInX1, nodesInX2, nodesInX3);
    auto gridCube = makeSimulationBoundingBox();

    generateBlockGrid(gridCube);
    setKernelForcing(lbmKernel, converter);
    setBoundaryConditionProcessor(lbmKernel);

    auto metisVisitor = std::make_shared<MetisPartitioningGridVisitor>(communicator,
                                                                       MetisPartitioningGridVisitor::LevelBased,
                                                                       vf::lbm::dir::d00M, MetisPartitioner::RECURSIVE);

    InteractorsHelper intHelper(grid, metisVisitor);
    for (auto const &interactor : interactors)
        intHelper.addInteractor(interactor);

    intHelper.selectBlocks();

    int numberOfProcesses = communicator->getNumberOfProcesses();
    SetKernelBlockVisitor kernelVisitor(lbmKernel, physicalParameters->latticeViscosity,
                                        numberOfProcesses);
    grid->accept(kernelVisitor);
    intHelper.setBC();


    writeBlocksToFile(); // important: run this after metis & intHelper.selectBlocks()
    setConnectors();
    initializeDistributions();
    grid->accept(bcVisitor);
    writeBoundaryConditions();

#ifdef _OPENMP
    omp_set_num_threads(simulationParameters->numberOfThreads);
    if (isMainProcess())
        VF_LOG_INFO("OpenMP is set to run with {} threads", simulationParameters->numberOfThreads);
#endif

    auto visualizationScheduler = std::make_shared<UbScheduler>(simulationParameters->timeStepLogInterval);
    auto mqCoProcessor = makeMacroscopicQuantitiesCoProcessor(converter,
                                                              visualizationScheduler);

    auto nupsCoProcessor = makeNupsCoProcessor();

    auto simulation = std::make_shared<Simulation>(grid, visualizationScheduler,
                                                        simulationParameters->numberOfTimeSteps);
    simulation->addSimulationObserver(nupsCoProcessor);
    simulation->addSimulationObserver(mqCoProcessor);

    if (isMainProcess()) VF_LOG_TRACE("Simulation start");
    simulation->run();
    if (isMainProcess()) VF_LOG_TRACE("Simulation end");
}

bool CPUSimulation::isMainProcess()
{
    return communicator->getProcessID() == 0;
}

void CPUSimulation::setKernelForcing(const std::shared_ptr<LBMKernel> &kernel,
                             std::shared_ptr<LBMUnitConverter> &converter) const
{
    kernel->setWithForcing(kernelConfig->useForcing);
    kernel->setForcingX1(kernelConfig->forcingX1 * converter->getFactorForceWToLb());
    kernel->setForcingX2(kernelConfig->forcingX2 * converter->getFactorForceWToLb());
    kernel->setForcingX3(kernelConfig->forcingX3 * converter->getFactorForceWToLb());
}

void CPUSimulation::logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
{
    VF_LOG_INFO("Domain size = {} x {} x {}", nodesInX1, nodesInX2, nodesInX3);
    VF_LOG_INFO("dx          = {} m", gridParameters->nodeDistance);
    VF_LOG_INFO("latticeViscosity    = {}", physicalParameters->latticeViscosity);
}

void CPUSimulation::generateBlockGrid(const std::shared_ptr<GbObject3D> &gridCube) const
{
    VF_LOG_TRACE("Generate block grid");
    GenBlocksGridVisitor genBlocks(gridCube);
    grid->accept(genBlocks);
}

void CPUSimulation::setBoundaryConditionProcessor(const std::shared_ptr<LBMKernel> &kernel)
{
    VF_LOG_TRACE("Create boundary conditions processor");
    auto bcProc = std::make_shared<BCSet>();
    kernel->setBCSet(bcProc);
}

void CPUSimulation::setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const
{
    int blockSizeX1 = nodesInX1 / gridParameters->blocksPerDirection[0];
    int blockSizeX2 = nodesInX2 / gridParameters->blocksPerDirection[1];
    int blockSizeX3 = nodesInX3 / gridParameters->blocksPerDirection[2];
    VF_LOG_INFO("Block size  = {} x {} x {}", blockSizeX1, blockSizeX2, blockSizeX3);
    grid->setBlockNX(blockSizeX1,
                     blockSizeX2,
                     blockSizeX3);
}

std::shared_ptr<LBMUnitConverter>
CPUSimulation::makeLBMUnitConverter()
{
    return std::make_shared<LBMUnitConverter>();
}



void CPUSimulation::writeBoundaryConditions() const
{
    auto geoSch = std::make_shared<UbScheduler>(1);
    WriteBoundaryConditionsSimulationObserver ppgeo(grid, geoSch, writerConfig.outputPath, writerConfig.getWriter(),
                                             communicator);
    ppgeo.update(0);
}

void CPUSimulation::writeBlocksToFile() const
{
    VF_LOG_TRACE("Write block grid to VTK-file");
    auto scheduler = std::make_shared<UbScheduler>(1);
    auto ppblocks = std::make_shared<WriteBlocksSimulationObserver>(grid,
                                                             scheduler,
                                                             writerConfig.outputPath,
                                                             writerConfig.getWriter(),
                                                             communicator);
    ppblocks->update(0);
    ppblocks.reset();
}

std::shared_ptr<GbObject3D> CPUSimulation::makeSimulationBoundingBox()
{
    auto box = gridParameters->boundingBox();
    auto gridCube = std::make_shared<GbCuboid3D>(box->minX1, box->minX2, box->minX3, box->maxX1, box->maxX2, box->maxX3);

    if (isMainProcess()) {
        VF_LOG_INFO("Bounding box dimensions = [({}},{}},{}}); ({}}, {}}, {}})]", box->minX1, box->minX2, box->minX3, box->maxX1, box->maxX2, box->maxX3);

        GbSystem3D::writeGeoObject(gridCube.get(), writerConfig.outputPath + "/geo/gridCube", writerConfig.getWriter());
    }

    return gridCube;
}

void CPUSimulation::setConnectors()
{
    OneDistributionSetConnectorsBlockVisitor setConnsVisitor(communicator);
    grid->accept(setConnsVisitor);
}

void CPUSimulation::initializeDistributions()
{
    InitDistributionsBlockVisitor initVisitor;
    grid->accept(initVisitor);
}

std::shared_ptr<SimulationObserver> CPUSimulation::makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter, const std::shared_ptr<UbScheduler> &visualizationScheduler) const
{
    auto mqCoProcessor = std::make_shared<WriteMacroscopicQuantitiesSimulationObserver>(grid, visualizationScheduler, writerConfig.outputPath, writerConfig.getWriter(), converter, communicator);
    mqCoProcessor->update(0);
    return mqCoProcessor;
}

std::shared_ptr<SimulationObserver> CPUSimulation::makeNupsCoProcessor() const
{
    auto scheduler = std::make_shared<UbScheduler>(100, 100);
    return std::make_shared<NUPSCounterSimulationObserver>(grid, scheduler,
                                                    simulationParameters->numberOfThreads,
                                                    communicator);
}
