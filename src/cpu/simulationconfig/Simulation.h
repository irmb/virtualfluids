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
#ifndef VIRTUALFLUIDSBUILD_H
#define VIRTUALFLUIDSBUILD_H

#include <string>
#include <memory>
#include <set>

#include <parallel/Communicator.h>

#include <geometry3d/GbPoint3D.h>
#include <Interactors/Interactor3D.h>
#include <BoundaryConditions/BC.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <SimulationObservers/SimulationObserver.h>
#include <LBM/LBMUnitConverter.h>
#include "KernelFactory.h"
#include "AbstractLBMSystem.h"
#include "KernelConfigStructs.h"
#include "SimulationParameters.h"
#include "WriterConfiguration.h"



class CPUSimulation
{
public:
    CPUSimulation();

    void setWriterConfiguration(WriterConfiguration config);

    void setGridParameters(std::shared_ptr<GridParameters> parameters);

    void setPhysicalParameters(std::shared_ptr<PhysicalParameters> parameters);

    void setRuntimeParameters(std::shared_ptr<RuntimeParameters> parameters);

    void setKernelConfiguration(const std::shared_ptr<LBMKernelConfiguration> &kernel);

    void addObject(const std::shared_ptr<GbObject3D> &object, const std::shared_ptr<BC> &bcAdapter, int state,
                   const std::string &folderPath);

    void addBCAdapter(const std::shared_ptr<BC> &bcAdapter);

    void run();

private:
    bool isMainProcess();

    std::shared_ptr<GbObject3D> makeSimulationBoundingBox();

    void writeBlocksToFile() const;

    void writeBoundaryConditions() const;

    std::shared_ptr<SimulationObserver> makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                           const std::shared_ptr<UbScheduler> &visualizationScheduler) const;

    static std::shared_ptr<LBMUnitConverter> makeLBMUnitConverter();

    void setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    static void setBoundaryConditionProcessor(const std::shared_ptr<LBMKernel> &kernel);

    void generateBlockGrid(const std::shared_ptr<GbObject3D> &gridCube) const;

    void logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    void setKernelForcing(const std::shared_ptr<LBMKernel> &kernel, std::shared_ptr<LBMUnitConverter> &converter) const;

    void setConnectors();

    void initializeDistributions();

    std::shared_ptr<SimulationObserver> makeNupsCoProcessor() const;

private:
    KernelFactory kernelFactory = KernelFactory();

    std::shared_ptr<LBMKernel> lbmKernel;
    std::shared_ptr<AbstractLBMSystem> lbmSystem;
    std::shared_ptr<vf::parallel::Communicator> communicator;

    std::shared_ptr<Grid3D> grid;
    std::vector<std::shared_ptr<Interactor3D>> interactors;
    BoundaryConditionsBlockVisitor bcVisitor {};
    std::set<std::shared_ptr<BC>> registeredAdapters;

    std::shared_ptr<LBMKernelConfiguration> kernelConfig;
    std::shared_ptr<RuntimeParameters> simulationParameters;
    std::shared_ptr<GridParameters> gridParameters;
    std::shared_ptr<PhysicalParameters> physicalParameters;

    WriterConfiguration writerConfig;
};


#endif