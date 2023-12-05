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
//! \file Simulation.cpp
//! \ingroup Simulation
//! \author Konstantin Kutscher
//=======================================================================================

#include "Simulation.h"

#include "Block3D.h"
#include "Block3DConnector.h"
#include "SimulationObserver.h"
#include "Grid3D.h"
#include "BCSet.h"
#include "Block3DConnector.h"
#include "LBMKernel.h"
#include "UbLogger.h"
#include "UbScheduler.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#define OMP_SCHEDULE guided

// #define TIMING
#ifdef TIMING
#include <basics/Timer/Timer.h>
using namespace vf::basics;
#endif

#include <basics/utilities/UbException.h>


Simulation::Simulation(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps)
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
Simulation::~Simulation() = default;
//////////////////////////////////////////////////////////////////////////
void Simulation::addSimulationObserver(SPtr<SimulationObserver> coProcessor) { simulationObserver.push_back(coProcessor); }
//////////////////////////////////////////////////////////////////////////
void Simulation::notifyObservers(real step)
{
    for (SPtr<SimulationObserver> cp : simulationObserver) {
        cp->update(step);
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::initLocalConnectors()
{
    UBLOG(logDEBUG1, "Simulation::initLocalConnectors() - start");

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
        UBLOG(logDEBUG5, "Simulation::initConnectors()-initConnectors(localConns[" << l << "])");
        initConnectors(localConns[l]);

        if (l != maxLevel) {
            UBLOG(logDEBUG5, "Simulation::initConnectors()-initConnectors(localInterConns[" << l << "])");
            initConnectors(localInterConns[l]);
        }
    }

    UBLOG(logDEBUG1, "Simulation::initLocalConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Simulation::initRemoteConnectors()
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
        UBLOG(logDEBUG5, "Simulation::initRemoteConnectors()-initConnectors(remoteConns[" << l << "])");
        initConnectors(remoteConns[l]);
        if (l != maxLevel) {
            for (size_t i = 0; i < remoteInterConnsCF[l].size(); i++)
                remoteInterConns[l].push_back(remoteInterConnsCF[l][i]);
            for (size_t i = 0; i < remoteInterConnsFC[l + 1].size(); i++)
                remoteInterConns[l].push_back(remoteInterConnsFC[l + 1][i]);
        }
    }
    //////////////////////////////////////////////////////////////////////////
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - connectoren initialisieren - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5, "Simulation::initRemoteConnectors()-initConnectors(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->init();
        }
    }
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - connectoren initialisieren - end");
    //////////////////////////////////////////////////////////////////////////
    // sendTransmitterDataSize
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - sendTransmitterDataSize - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5,
                  "Simulation::initRemoteConnectors()-sendTransmitterDataSize(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->sendTransmitterDataSize();
        }
    }
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - sendTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
    // receiveTransmitterDataSize
    // if it stops here during distributed calculations, then there is probably an inactive block on one side!!!
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - receiveTransmitterDataSize - start");
    for (int l = minLevel; l <= maxLevel; l++) {
        if (l != maxLevel) {
            UBLOG(logDEBUG5,
                  "Simulation::initRemoteConnectors()-receiveTransmitterDataSize(remoteInterConns[" << l << "])");
            for (SPtr<Block3DConnector> c : remoteInterConns[l])
                c->receiveTransmitterDataSize();
        }
    }
    // UBLOG(logDEBUG5, "Simulation::initConnectors() - receiveTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////
void Simulation::initConnectors(std::vector<SPtr<Block3DConnector>> &connectors)
{
    UBLOG(logDEBUG1, "Simulation::initConnectors() - start");

    // initialization
    //////////////////////////////////////////////////////////////////////////
    // initialize connectors
    UBLOG(logDEBUG5, "Simulation::initConnectors() - connectoren initialisieren - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->init();
    UBLOG(logDEBUG5, "Simulation::initConnectors() - connectoren initialisieren - end");
    //////////////////////////////////////////////////////////////////////////
    // sendTransmitterDataSize
    UBLOG(logDEBUG5, "Simulation::initConnectors() - sendTransmitterDataSize - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->sendTransmitterDataSize();
    UBLOG(logDEBUG5, "Simulation::initConnectors() - sendTransmitterDataSize - end");
    //////////////////////////////////////////////////////////////////////////
    // receiveTransmitterDataSize
    // if it stops here during distributed calculations, then there is probably an inactive block on one side!!!
    UBLOG(logDEBUG5, "Simulation::initConnectors() - receiveTransmitterDataSize - start");
    for (SPtr<Block3DConnector> c : connectors)
        c->receiveTransmitterDataSize();
    UBLOG(logDEBUG5, "Simulation::initConnectors() - receiveTransmitterDataSize - end");

    UBLOG(logDEBUG1, "Simulation::initConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void Simulation::deleteBlocks()
{
    for (std::vector<SPtr<Block3D>> &bs : blocks)
        bs.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Simulation::deleteConnectors()
{
    deleteConnectors(localConns);
    deleteConnectors(remoteConns);

    deleteConnectors(localInterConns);
    deleteConnectors(remoteInterConns);
}
//////////////////////////////////////////////////////////////////////////
void Simulation::deleteConnectors(std::vector<std::vector<SPtr<Block3DConnector>>> &conns)
{
    for (std::vector<SPtr<Block3DConnector>> &c : conns)
        c.resize(0);
}
//////////////////////////////////////////////////////////////////////////
void Simulation::run()
{
    UBLOG(logDEBUG1, "OMPSimulation::calculate() - started");
    int calcStep = 0;
    try {
        int minInitLevel = minLevel;
        int maxInitLevel = maxLevel - minLevel;
        int straightStartLevel = minInitLevel;
        int internalIterations = 1 << (maxInitLevel - minInitLevel);
        int threshold;

#ifdef TIMING
        Timer timer;
        real time[6];
#endif

        for (calcStep = startTimeStep; calcStep <= numberOfTimeSteps; calcStep++) {
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            UBLOG(logINFO, "calcStep = " << calcStep);
#endif
            //////////////////////////////////////////////////////////////////////////

            for (int staggeredStep = 1; staggeredStep <= internalIterations; staggeredStep++) {
                if (staggeredStep == internalIterations)
                    straightStartLevel = minInitLevel;
                else {
                    for (straightStartLevel = maxInitLevel, threshold = 1; (staggeredStep & threshold) != threshold; straightStartLevel--, threshold <<= 1)
                        ;
                }
#ifdef TIMING   
                timer.start();
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPreCollisionBC(straightStartLevel, maxInitLevel);

                // do collision for all blocks
                calculateBlocks(straightStartLevel, maxInitLevel, calcStep);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[0] = timer.getCurrentRuntimeInSeconds();
                UBLOG(logINFO, "calculateBlocks time = " << time[0]);
                timer.start();
#endif
                //////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////
                // exchange data between blocks
                exchangeBlockData(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[1] = timer.getCurrentRuntimeInSeconds();
                UBLOG(logINFO, "exchangeBlockData time = " << time[1]);
                timer.start();
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPostCollisionBC(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[2] = timer.getCurrentRuntimeInSeconds();
                UBLOG(logINFO, "applyBCs time = " << time[2]);
                timer.start();
#endif
                //////////////////////////////////////////////////////////////////////////
                // swap distributions in kernel
                swapDistributions(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[3] = timer.getCurrentRuntimeInSeconds();
                UBLOG(logINFO, "swapDistributions time = " << time[3]);
                timer.start();
#endif
                //////////////////////////////////////////////////////////////////////////
                if (refinement) {
                    if (straightStartLevel < maxInitLevel) exchangeBlockData(straightStartLevel, maxInitLevel);
                        //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[4] = timer.getCurrentRuntimeInSeconds();
                    UBLOG(logINFO, "refinement exchangeBlockData time = " << time[4]);
                    timer.start();
#endif
                    //////////////////////////////////////////////////////////////////////////
                    // now ghost nodes have actual values
                    // interpolation of interface nodes between grid levels
                    interpolation(straightStartLevel, maxInitLevel);
                    //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[5] = timer.getCurrentRuntimeInSeconds();
                    UBLOG(logINFO, "refinement interpolation time = " << time[5]);
                    timer.start();
#endif
                    //////////////////////////////////////////////////////////////////////////
                }
            }
            // exchange data between blocks for visualization
            if (additionalGhostLayerUpdateScheduler->isDue(calcStep)) {
                exchangeBlockData(straightStartLevel, maxInitLevel);
            }
            notifyObservers((real)(calcStep));
            // now ghost nodes have actual values
        }
        UBLOG(logDEBUG1, "OMPSimulation::calculate() - stoped");
    } catch (std::exception &e) {
        UBLOG(logERROR, e.what());
        UBLOG(logERROR, " step = " << calcStep);
        // throw e;
        // exit(EXIT_FAILURE);
    } catch (std::string &s) {
        UBLOG(logERROR, s);
        // exit(EXIT_FAILURE);
        // throw s;
    } catch (...) {
        UBLOG(logERROR, "unknown exception");
        // exit(EXIT_FAILURE);
        // throw;
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::calculateBlocks(int startLevel, int maxInitLevel, int calcStep)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        SPtr<Block3D> blockTemp;
        // startLevel bis maxInitLevel
        for (int level = startLevel; level <= maxInitLevel; level++) {
#ifdef TIMING
            Timer timer;
            timer.start();
#endif
            // call LBM kernel
            int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp for schedule(OMP_SCHEDULE)
#endif
            for (int i = 0; i < size; i++) {
                try {
                    blockTemp = blocks[level][i];
                    blockTemp->getKernel()->calculate(calcStep);
                } catch (std::exception &e) {
                    UBLOG(logERROR, e.what());
                    UBLOG(logERROR, blockTemp->toString() << " step = " << calcStep);
                    std::exit(EXIT_FAILURE);
                }
            }
#ifdef TIMING
            timer.end();
            UBLOG(logINFO, "level = " << level << " blocks = " << blocks[level].size()
                                      << " collision time = " << timer.getTimeInSeconds());
#endif
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::exchangeBlockData(int startLevel, int maxInitLevel)
{
    // startLevel bis maxInitLevel
    for (int level = startLevel; level <= maxInitLevel; level++) {
        // connectorsPrepareLocal(localConns[level]);
        connectorsSendLocal(localConns[level]);
        // connectorsReceiveLocal(localConns[level]);

        connectorsPrepareRemote(remoteConns[level]);
        connectorsSendRemote(remoteConns[level]);
        connectorsReceiveRemote(remoteConns[level]);
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::swapDistributions(int startLevel, int maxInitLevel)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        // startLevel bis maxInitLevel
        for (int level = startLevel; level <= maxInitLevel; level++) {
            int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp for schedule(OMP_SCHEDULE)
#endif
            for (int i = 0; i < size; i++) {
                blocks[level][i]->getKernel()->swapDistributions();
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::connectorsPrepareLocal(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
    for (int i = 0; i < size; i++) {
        try {
            connectors[i]->prepareForReceive();
            connectors[i]->prepareForSend();
        } catch (std::exception &e) {
            UBLOG(logERROR, e.what());
            std::exit(EXIT_FAILURE);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::connectorsSendLocal(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
    for (int i = 0; i < size; i++) {
        try {
            connectors[i]->fillSendVectors();
            connectors[i]->sendVectors();
        } catch (std::exception &e) {
            UBLOG(logERROR, e.what());
            std::exit(EXIT_FAILURE);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::connectorsReceiveLocal(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
    for (int i = 0; i < size; i++) {
        connectors[i]->receiveVectors();
        connectors[i]->distributeReceiveVectors();
    }
}
void Simulation::connectorsPrepareRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->prepareForReceive();
        connectors[i]->prepareForSend();
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::connectorsSendRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->fillSendVectors();
        connectors[i]->sendVectors();
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::connectorsReceiveRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->receiveVectors();
        connectors[i]->distributeReceiveVectors();
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::interpolation(int startLevel, int maxInitLevel)
{
    for (int level = startLevel; level < maxInitLevel; level++) {
        connectorsPrepareLocal(localInterConns[level]);
        connectorsPrepareRemote(remoteInterConns[level]);
    }

    for (int level = startLevel; level < maxInitLevel; level++) {
        connectorsSendLocal(localInterConns[level]);
        connectorsSendRemote(remoteInterConns[level]);
    }

    for (int level = startLevel; level < maxInitLevel; level++) {
        connectorsReceiveLocal(localInterConns[level]);
        connectorsReceiveRemote(remoteInterConns[level]);
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
    // from startLevel to maxInitLevel
    for (int level = startLevel; level <= maxInitLevel; level++) {
        int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
        for (int i = 0; i < size; i++) {
            try {
                blocks[level][i]->getKernel()->getBCSet()->applyPreCollisionBC();
            } catch (std::exception &e) {
                UBLOG(logERROR, e.what());
                exit(EXIT_FAILURE);
            } catch (std::string &s) {
                UBLOG(logERROR, s);
                exit(EXIT_FAILURE);
            } catch (...) {
                UBLOG(logERROR, "unknown exception");
                exit(EXIT_FAILURE);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Simulation::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
    //  from startLevel to maxInitLevel
    for (int level = startLevel; level <= maxInitLevel; level++) {
        int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
        for (int i = 0; i < size; i++) {
            try {
                blocks[level][i]->getKernel()->getBCSet()->applyPostCollisionBC();
            } catch (std::exception &e) {
                UBLOG(logERROR, e.what());
                exit(EXIT_FAILURE);
            } catch (std::string &s) {
                UBLOG(logERROR, s);
                exit(EXIT_FAILURE);
            } catch (...) {
                UBLOG(logERROR, "unknown exception");
                exit(EXIT_FAILURE);
            }
        }
    }
}