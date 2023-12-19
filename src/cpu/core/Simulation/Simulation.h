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
//! \addtogroup cpu_Simulation Simulation
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef Simulation_H
#define Simulation_H

#include <PointerDefinitions.h>
#include <vector>
#include "lbm/constants/D3Q27.h"

class Grid3D;
class UbScheduler;
class Block3D;
class Block3DConnector;
class SimulationObserver;

//! \class Simulation
//! \brief A base class for main simulation loop

class Simulation
{
public:
    Simulation(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps);
    virtual ~Simulation();
    //! control of coProcessors
    void addSimulationObserver(SPtr<SimulationObserver> coProcessor);
    void notifyObservers(real step);

    virtual void run();

protected:
    virtual void initLocalConnectors();
    virtual void initRemoteConnectors();
    void initConnectors(std::vector<SPtr<Block3DConnector>> &connectors);
    void deleteBlocks();
    void deleteConnectors();
    void deleteConnectors(std::vector<std::vector<SPtr<Block3DConnector>>> &conns);

    void calculateBlocks(int startLevel, int maxInitLevel, int calcStep);
    void swapDistributions(int startLevel, int maxInitLevel);
    void exchangeBlockData(int startLevel, int maxInitLevel);
    void connectorsPrepareLocal(std::vector<SPtr<Block3DConnector>> &connectors);
    void connectorsSendLocal(std::vector<SPtr<Block3DConnector>> &connectors);
    void connectorsReceiveLocal(std::vector<SPtr<Block3DConnector>> &connectors);
    void connectorsPrepareRemote(std::vector<SPtr<Block3DConnector>> &connectors);
    void connectorsSendRemote(std::vector<SPtr<Block3DConnector>> &connectors);
    void connectorsReceiveRemote(std::vector<SPtr<Block3DConnector>> &connectors);
    void interpolation(int startLevel, int maxInitLevel);
    void applyPreCollisionBC(int startLevel, int maxInitLevel);
    void applyPostCollisionBC(int startLevel, int maxInitLevel);

    int minLevel, maxLevel;
    int startTimeStep;
    int numberOfTimeSteps;
    std::vector<std::vector<SPtr<Block3DConnector>>> localConns;
    std::vector<std::vector<SPtr<Block3DConnector>>> remoteConns;

    bool refinement;
    SPtr<Grid3D> grid;
    SPtr<UbScheduler> additionalGhostLayerUpdateScheduler;
    std::vector<std::vector<SPtr<Block3D>>> blocks;

    // localInterConns and remoteInterConns save interpolation connectors
    // every element save CF connectors for current level and FC connectors for next level
    // e.g.
    // localInterConns[0] = CF(0), FC(1)
    // localInterConns[1] = CF(1), FC(2)
    // localInterConns[2]
    std::vector<std::vector<SPtr<Block3DConnector>>> localInterConns;
    std::vector<std::vector<SPtr<Block3DConnector>>> remoteInterConns;

    std::vector<SPtr<SimulationObserver>> simulationObserver;
};

#endif

//! \}
