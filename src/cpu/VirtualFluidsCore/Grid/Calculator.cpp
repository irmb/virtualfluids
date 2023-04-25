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
//! \file Calculator.cpp
//! \ingroup Grid
//! \author Konstantin Kutscher
//=======================================================================================

#include "Calculator.h"

#include "Block3D.h"
#include "Block3DConnector.h"
#include "CoProcessor.h"
#include "Grid3D.h"
#include "UbScheduler.h"

#include <basics/utilities/UbException.h>

Calculator::Calculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps)
    : grid(grid), additionalGhostLayerUpdateScheduler(additionalGhostLayerUpdateScheduler),
      numberOfTimeSteps(numberOfTimeSteps)
{
    this->grid    = grid;
    startTimeStep = int(grid->getTimeStep()) + 1;
    minLevel      = grid->getCoarsestInitializedLevel();
    maxLevel      = grid->getFinestInitializedLevel();
    if (maxLevel > 0)
        refinement = true;
    else
        refinement = false;
    blocks.resize(maxLevel + 1);
    localConns.resize(maxLevel + 1);
    remoteConns.resize(maxLevel + 1);
    localInterConns.resize(maxLevel);
    remoteInterConns.resize(maxLevel);

    int gridRank = grid->getRank();

    for (int level = minLevel; level <= maxLevel; level++) {
        std::vector<SPtr<Block3D>> blockVector;
        grid->getBlocks(level, gridRank, true, blockVector);
        for (const auto &block : blockVector)
            if (block)
                blocks[block->getLevel()].push_back(block);
    }

    initLocalConnectors();
    initRemoteConnectors();
}
//////////////////////////////////////////////////////////////////////////
Calculator::~Calculator() = default;
//////////////////////////////////////////////////////////////////////////
void Calculator::addCoProcessor(SPtr<CoProcessor> coProcessor) { coProcessors.push_back(coProcessor); }
//////////////////////////////////////////////////////////////////////////
void Calculator::coProcess(real step)
{
    for (SPtr<CoProcessor> cp : coProcessors) {
        cp->process(step);
    }
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initLocalConnectors()
{
    UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - start");

    for (int l = minLevel; l <= maxLevel; l++) {
        for (SPtr<Block3D> block : blocks[l]) {
            block->pushBackLocalSameLevelConnectors(localConns[l]);

            if (l != maxLevel)
                block->pushBackLocalInterpolationConnectorsCF(localInterConns[l]);
        }
        if (l != maxLevel) {
            for (SPtr<Block3D> block : blocks[l + 1]) {
                block->pushBackLocalInterpolationConnectorsFC(localInterConns[l]);
            }
        }
        UBLOG(logDEBUG5, "Calculator::initConnectors()-initConnectors(localConns[" << l << "])");
        initConnectors(localConns[l]);

        if (l != maxLevel) {
            UBLOG(logDEBUG5, "Calculator::initConnectors()-initConnectors(localInterConns[" << l << "])");
            initConnectors(localInterConns[l]);
        }
    }

    UBLOG(logDEBUG1, "Calculator::initLocalConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initRemoteConnectors()
{
    std::vector<std::vector<SPtr<Block3DConnector>>> remoteInterConnsCF;
    std::vector<std::vector<SPtr<Block3DConnector>>> remoteInterConnsFC;
    remoteInterConnsCF.resize(maxLevel + 1);
    remoteInterConnsFC.resize(maxLevel + 1);

    for (int l = minLevel; l <= maxLevel; l++) {
        std::vector<SPtr<Block3D>> blockVector;
        // grid->getBlocks(level, gridRank, true, blockVector);
        grid->getBlocks(l, blockVector);
        for (SPtr<Block3D> block : blockVector) {
            int block_level = block->getLevel();
            block->pushBackRemoteSameLevelConnectors(remoteConns[block_level]);

            block->pushBackRemoteInterpolationConnectorsCF(remoteInterConnsCF[block_level]);
            block->pushBackRemoteInterpolationConnectorsFC(remoteInterConnsFC[block_level]);
        }
    }

    for (int l = minLevel; l <= maxLevel; l++) {
        UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteConns[" << l << "])");
        initConnectors(remoteConns[l]);
        if (l != maxLevel) {
            for (size_t i = 0; i < remoteInterConnsCF[l].size(); i++)
                remoteInterConns[l].push_back(remoteInterConnsCF[l][i]);
            for (size_t i = 0; i < remoteInterConnsFC[l + 1].size(); i++)
                remoteInterConns[l].push_back(remoteInterConnsFC[l + 1][i]);
        }
    }
    //////////////////////////////////////////////////////////////////////////
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5, "Calculator::initRemoteConnectors()-initConnectors(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->init();
        }
    }
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - end");
    //////////////////////////////////////////////////////////////////////////
    // sendTransmitterDataSize
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5,
                  "Calculator::initRemoteConnectors()-sendTransmitterDataSize(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->sendTransmitterDataSize();
        }
    }
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
    // receiveTransmitterDataSize
    // wenn er hier bei verteilten berechnungen stopped, dann ist vermutlich auf einer seite ein nicht aktiver block!!!
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5,
                  "Calculator::initRemoteConnectors()-receiveTransmitterDataSize(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->receiveTransmitterDataSize();
        }
    }
    // UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////
void Calculator::initConnectors(std::vector<SPtr<Block3DConnector>> &connectors)
{
    UBLOG(logDEBUG1, "Calculator::initConnectors() - start");

    // initialization
    //////////////////////////////////////////////////////////////////////////
    // initialize connectors
    UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->init();
    UBLOG(logDEBUG5, "Calculator::initConnectors() - connectoren initialisieren - end");
    //////////////////////////////////////////////////////////////////////////
    // sendTransmitterDataSize
    UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->sendTransmitterDataSize();
    UBLOG(logDEBUG5, "Calculator::initConnectors() - sendTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
    // receiveTransmitterDataSize
    // wenn er hier bei verteilten berechnungen stopped, dann ist vermutlich auf einer seite ein nicht aktiver block!!!
    UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->receiveTransmitterDataSize();
    UBLOG(logDEBUG5, "Calculator::initConnectors() - receiveTransmitterDataSize - end");

    UBLOG(logDEBUG1, "Calculator::initConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteBlocks()
{
    for (std::vector<SPtr<Block3D>> &bs : blocks)
        bs.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors()
{
    deleteConnectors(localConns);
    deleteConnectors(remoteConns);

    deleteConnectors(localInterConns);
    deleteConnectors(remoteInterConns);
}
//////////////////////////////////////////////////////////////////////////
void Calculator::deleteConnectors(std::vector<std::vector<SPtr<Block3DConnector>>> &conns)
{
    for (std::vector<SPtr<Block3DConnector>> &c : conns)
        c.resize(0);
}
//////////////////////////////////////////////////////////////////////////
