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
//! \file BasicCalculator.cpp
//! \ingroup Grid
//! \author Konstantin Kutscher
//=======================================================================================

#include "BasicCalculator.h"

#include "BCProcessor.h"
#include "Block3D.h"
#include "Block3DConnector.h"
#include "LBMKernel.h"
#include "UbLogger.h"
#include "UbScheduler.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#define OMP_SCHEDULE guided

//#define TIMING
//#include "UbTiming.h"

BasicCalculator::BasicCalculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler,
                                 int numberOfTimeSteps)
    : Calculator(grid, additionalGhostLayerUpdateScheduler, numberOfTimeSteps)
{
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::calculate()
{
    UBLOG(logDEBUG1, "OMPCalculator::calculate() - started");
    int calcStep = 0;
    try {
        int minInitLevel       = minLevel;
        int maxInitLevel       = maxLevel - minLevel;
        int straightStartLevel = minInitLevel;
        int internalIterations = 1 << (maxInitLevel - minInitLevel);
        int threshold;

#ifdef TIMING
        UbTimer timer;
        real time[6];
#endif

        for (calcStep = startTimeStep; calcStep <= numberOfTimeSteps; calcStep++) {
            //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
            UBLOG(logINFO, "calcStep = " << calcStep);
#endif
         //////////////////////////////////////////////////////////////////////////

         for (int staggeredStep = 1; staggeredStep <= internalIterations; staggeredStep++)
         {
            if (staggeredStep == internalIterations) straightStartLevel = minInitLevel;
            else
            {
               for (straightStartLevel = maxInitLevel, threshold = 1;
                  (staggeredStep & threshold) != threshold; straightStartLevel--, threshold <<= 1);
            }
#ifdef TIMING
                timer.resetAndStart();
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPreCollisionBC(straightStartLevel, maxInitLevel);

                // do collision for all blocks
                calculateBlocks(straightStartLevel, maxInitLevel, calcStep);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[0] = timer.stop();
                UBLOG(logINFO, "calculateBlocks time = " << time[0]);
#endif
                //////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////
                // exchange data between blocks
                exchangeBlockData(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[1] = timer.stop();
                UBLOG(logINFO, "exchangeBlockData time = " << time[1]);
#endif
                //////////////////////////////////////////////////////////////////////////
                applyPostCollisionBC(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[2] = timer.stop();
                UBLOG(logINFO, "applyBCs time = " << time[2]);
#endif
                //////////////////////////////////////////////////////////////////////////
                // swap distributions in kernel
                swapDistributions(straightStartLevel, maxInitLevel);
                //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                time[3] = timer.stop();
                UBLOG(logINFO, "swapDistributions time = " << time[3]);
#endif
                //////////////////////////////////////////////////////////////////////////
                if (refinement) {
                    if (straightStartLevel < maxInitLevel)
                        exchangeBlockData(straightStartLevel, maxInitLevel);
                        //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[4] = timer.stop();
                    UBLOG(logINFO, "refinement exchangeBlockData time = " << time[4]);
#endif
                    //////////////////////////////////////////////////////////////////////////
                    // now ghost nodes have actual values
                    // interpolation of interface nodes between grid levels
                    interpolation(straightStartLevel, maxInitLevel);
                    //////////////////////////////////////////////////////////////////////////
#ifdef TIMING
                    time[5] = timer.stop();
                    UBLOG(logINFO, "refinement interpolation time = " << time[5]);
#endif
                    //////////////////////////////////////////////////////////////////////////
                }
            }
            // exchange data between blocks for visualization
            if (additionalGhostLayerUpdateScheduler->isDue(calcStep)) {
                exchangeBlockData(straightStartLevel, maxInitLevel);
            }
            coProcess((real)(calcStep));
            // now ghost nodes have actual values
        }
        UBLOG(logDEBUG1, "OMPCalculator::calculate() - stoped");
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
void BasicCalculator::calculateBlocks(int startLevel, int maxInitLevel, int calcStep)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        SPtr<Block3D> blockTemp;
        // startLevel bis maxInitLevel
        for (int level = startLevel; level <= maxInitLevel; level++) {
            // timer.resetAndStart();
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
            // timer.stop();
            // UBLOG(logINFO, "level = " << level << " blocks = " << blocks[level].size() << " collision time = " <<
            // timer.getTotalTime());
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::exchangeBlockData(int startLevel, int maxInitLevel)
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
void BasicCalculator::swapDistributions(int startLevel, int maxInitLevel)
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
void BasicCalculator::connectorsPrepareLocal(std::vector<SPtr<Block3DConnector>> &connectors)
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
void BasicCalculator::connectorsSendLocal(std::vector<SPtr<Block3DConnector>> &connectors)
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
void BasicCalculator::connectorsReceiveLocal(std::vector<SPtr<Block3DConnector>> &connectors)
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
void BasicCalculator::connectorsPrepareRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->prepareForReceive();
        connectors[i]->prepareForSend();
    }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsSendRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->fillSendVectors();
        connectors[i]->sendVectors();
    }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::connectorsReceiveRemote(std::vector<SPtr<Block3DConnector>> &connectors)
{
    int size = (int)connectors.size();
    for (int i = 0; i < size; i++) {
        connectors[i]->receiveVectors();
        connectors[i]->distributeReceiveVectors();
    }
}
//////////////////////////////////////////////////////////////////////////
void BasicCalculator::interpolation(int startLevel, int maxInitLevel)
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
void BasicCalculator::applyPreCollisionBC(int startLevel, int maxInitLevel)
{
    // from startLevel to maxInitLevel
    for (int level = startLevel; level <= maxInitLevel; level++) {
        int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
        for (int i = 0; i < size; i++) {
            try {
                blocks[level][i]->getKernel()->getBCProcessor()->applyPreCollisionBC();
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
void BasicCalculator::applyPostCollisionBC(int startLevel, int maxInitLevel)
{
    //  from startLevel to maxInitLevel
    for (int level = startLevel; level <= maxInitLevel; level++) {
        int size = (int)blocks[level].size();
#ifdef _OPENMP
#pragma omp parallel for schedule(OMP_SCHEDULE)
#endif
        for (int i = 0; i < size; i++) {
            try {
                blocks[level][i]->getKernel()->getBCProcessor()->applyPostCollisionBC();
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
