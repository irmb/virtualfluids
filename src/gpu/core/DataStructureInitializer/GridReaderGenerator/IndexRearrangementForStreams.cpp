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
//! \addtogroup gpu_DataStructureInitializer DataStructureInitializer
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann
//=======================================================================================
#include "IndexRearrangementForStreams.h"

#include <algorithm>
#include <iostream>

#include <logger/Logger.h>

#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <parallel/Communicator.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "lbm/constants/D3Q27.h"

bool indexInArray(const uint* array, uint numberOfElements, uint index)
{
    return std::find(array, array + numberOfElements, index) != array + numberOfElements;
}

bool indexInVector(const std::vector<uint>& vector, uint index)
{
    return std::find(vector.begin(), vector.end(), index) != vector.end();
}

IndexRearrangementForStreams::IndexRearrangementForStreams(std::shared_ptr<Parameter> para,
                                                           std::shared_ptr<GridBuilder> builder,
                                                           vf::parallel::Communicator& communicator)
    : para(std::move(para)), builder(std::move(builder)), communicator(communicator)
{
}

std::array<ProcessNeighbor27, 4>
IndexRearrangementForStreams::initCommunicationArraysForCommAfterFinetoCoarse(const ProcessNeighbor27& sendNeighborHost,
                                                                              const ProcessNeighbor27& sendNeighborDevice,
                                                                              const ProcessNeighbor27& recvNeighborHost,
                                                                              const ProcessNeighbor27& recvNeighborDevice,
                                                                              const int level, const int direction) const
{
    VF_LOG_INFO("Communication: reorder send indices in direction {}", direction);

    const std::vector<uint> sendIndexPositions =
        reorderSendIndicesForCommAfterFtoC(sendNeighborHost.index, direction, level);

    const std::vector<uint> recvIndexPositions =
        exchangeIndicesForCommAfterFtoC(sendNeighborHost, recvNeighborHost, sendIndexPositions);

    reorderRecvIndicesForCommAfterFtoC(recvNeighborHost.index, direction, level, recvIndexPositions);

    const ProcessNeighbor27 sendNeighborAfterFtoCHost = makeProcessNeighborToCommAfterFtoC(sendNeighborHost, uint(sendIndexPositions.size()));
    const ProcessNeighbor27 sendNeighborAfterFtoCDevice = makeProcessNeighborToCommAfterFtoC(sendNeighborDevice, uint(sendIndexPositions.size()));
    const ProcessNeighbor27 recvNeighborAfterFtoCHost = makeProcessNeighborToCommAfterFtoC(recvNeighborHost, uint(recvIndexPositions.size()));
    const ProcessNeighbor27 recvNeighborAfterFtoCDevice = makeProcessNeighborToCommAfterFtoC(recvNeighborDevice, uint(recvIndexPositions.size()));

    return { sendNeighborAfterFtoCHost, sendNeighborAfterFtoCDevice, recvNeighborAfterFtoCHost, 
             recvNeighborAfterFtoCDevice };
}

std::vector<uint> IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoC(uint* indices,
                                                                                   int direction, int level) const
{
    VF_LOG_INFO("Reorder send indices for communication after fine to coarse: level: {} direction: {}", level, direction);
    
    const ICells& fineToCoarse = para->getParH(level)->fineToCoarse;
    if (para->getParH(level)->coarseToFine.numberOfCells == 0 || fineToCoarse.numberOfCells == 0)
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): para->getParH(level)->intCF needs to be initialized "
                        "before calling this function");

    std::vector<uint> sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions;

    const uint numberOfSendIndices = builder->getNumberOfSendIndices(direction, level);

    const std::vector<uint> aggregatedCoarseNodesForCtoF = aggregateCoarseNodesForCtoF(level);

    // coarse cells of interpolation fine to coarse (iCellFCC)
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        const uint sparseIndexSend = indices[posInSendIndices];
        if (sparseIndexSend == 0 || indexInVector(sendIndicesAfterFtoC, sparseIndexSend))
            continue;
        if (indexInArray(fineToCoarse.coarseCellIndices, fineToCoarse.numberOfCells, sparseIndexSend) ||
            indexInVector(aggregatedCoarseNodesForCtoF, sparseIndexSend)) {
            sendIndicesAfterFtoC.push_back(sparseIndexSend);
            sendIndicesForCommAfterFtoCPositions.push_back(posInSendIndices);
        }
    }

    const std::vector<uint> sendIndicesOther =
        findIndicesNotInCommAfterFtoC(numberOfSendIndices, indices, sendIndicesAfterFtoC);

    std::copy(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), indices);
    std::copy(sendIndicesOther.begin(), sendIndicesOther.end(), indices + sendIndicesAfterFtoC.size());

    VF_LOG_INFO("Reorder send indices: process {}, numberOfSendNodesAfterFtoC {}", communicator.getProcessID(),
                sendIndicesAfterFtoC.size());

    if (sendIndicesAfterFtoC.size() + sendIndicesOther.size() != numberOfSendIndices) {
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfSendNodesAfterFtoC = {}, sendIndicesOther.size() = {}, numberOfSendIndices = {}",
                        sendIndicesAfterFtoC.size(), sendIndicesOther.size(), numberOfSendIndices);
    }
    return sendIndicesForCommAfterFtoCPositions;
}

std::vector<uint>
IndexRearrangementForStreams::exchangeIndicesForCommAfterFtoC(const ProcessNeighbor27& sendNeighbor,
                                                              const ProcessNeighbor27& recvNeighbor,
                                                              const std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
{
    // fill the receive vector with zeros as placeholders
    // give vector an arbitrary size (larger than needed) // TODO: Find a better way
    std::vector<uint> recvIndicesForCommAfterFtoCPositions((size_t)sendNeighbor.numberOfNodes * 2, 0);

    communicator.receiveSend(recvIndicesForCommAfterFtoCPositions.data(), (int)recvIndicesForCommAfterFtoCPositions.size(),
                             recvNeighbor.rankNeighbor, sendIndicesForCommAfterFtoCPositions.data(),
                             (int)sendIndicesForCommAfterFtoCPositions.size(), sendNeighbor.rankNeighbor);

    // resize receiving vector to correct size
    if (!recvIndicesForCommAfterFtoCPositions.empty()) {
        auto it = std::unique(
            recvIndicesForCommAfterFtoCPositions.begin(),
            recvIndicesForCommAfterFtoCPositions.end()); // finds the second zero when there are multiple zeros in a row
        recvIndicesForCommAfterFtoCPositions.erase(std::prev(it, 1), // begin erasing at the first zero
                                                   recvIndicesForCommAfterFtoCPositions.end()); // TODO: Find a better way
    }

    return recvIndicesForCommAfterFtoCPositions;
}


void IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoC(
    uint* indices, const int direction, const int level,
    const std::vector<uint>& recvIndicesForCommAfterFtoCPositions) const
{
    VF_LOG_INFO("Reorder recv indices for communication after fine to coarse: level: {} direction: {}", level, direction);

    if (recvIndicesForCommAfterFtoCPositions.empty())
        VF_LOG_WARNING("ReorderRecvIndicesForCommAfterFtoC(): sendIndicesForCommAfterFtoCPositions is empty.");

    const uint numberOfRecvIndices = builder->getNumberOfReceiveIndices(direction, level);
    
    std::vector<uint> recvIndicesAfterFtoC;
    recvIndicesAfterFtoC.reserve(recvIndicesForCommAfterFtoCPositions.size());
   
    // find recvIndices for Communication after fine to coarse
    for (uint vectorPos : recvIndicesForCommAfterFtoCPositions)
        recvIndicesAfterFtoC.push_back(indices[vectorPos]);

    const std::vector<uint> recvIndicesOther = findIndicesNotInCommAfterFtoC(numberOfRecvIndices, indices, recvIndicesAfterFtoC);

    // copy new vectors back to recvIndices array
    std::copy(recvIndicesAfterFtoC.begin(), recvIndicesAfterFtoC.end(), indices);
    std::copy(recvIndicesOther.begin(), recvIndicesOther.end(), indices + recvIndicesAfterFtoC.size());

    VF_LOG_INFO("Reorder recv indices: process {}, numberOfRecvNodesAfterFtoC {}", communicator.getProcessID(),
                recvIndicesAfterFtoC.size());

    if (recvIndicesAfterFtoC.size() + recvIndicesOther.size() != numberOfRecvIndices) {
        VF_LOG_CRITICAL("reorderRecvIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfRecvNodesAfterFtoC = {}, recvIndicesOther.size() = {}, numberOfRecvIndices = {}",
                        recvIndicesAfterFtoC.size(), recvIndicesOther.size(), numberOfRecvIndices);
    }
}

ProcessNeighbor27 IndexRearrangementForStreams::makeProcessNeighborToCommAfterFtoC(const ProcessNeighbor27& neighbor,
                                                                      const uint numberOfNodes)
{
    ProcessNeighbor27 neighborAfterFtoC(neighbor);
    neighborAfterFtoC.numberOfNodes = numberOfNodes;
    neighborAfterFtoC.numberOfFs = vf::lbm::dir::NUMBER_Of_DIRECTIONS * numberOfNodes;
    neighborAfterFtoC.memsizeIndex = sizeof(uint) * numberOfNodes;
    neighborAfterFtoC.memsizeFs = sizeof(real) * vf::lbm::dir::NUMBER_Of_DIRECTIONS * numberOfNodes;
    return neighborAfterFtoC;
}

std::vector<uint> IndexRearrangementForStreams::aggregateCoarseNodesForCtoF(const int level) const
{
    std::vector<uint> nodesCFC;
    const uint* neighborX = para->getParH(level)->neighborX;
    const uint* neighborY = para->getParH(level)->neighborY;
    const uint* neighborZ = para->getParH(level)->neighborZ;

    for (uint x = 0; x < para->getParH(level)->coarseToFine.numberOfCells; x++) {
        const uint sparseIndex = para->getParH(level)->coarseToFine.coarseCellIndices[x];
        nodesCFC.push_back(sparseIndex);
        nodesCFC.push_back(neighborX[sparseIndex]);
        nodesCFC.push_back(neighborY[sparseIndex]);
        nodesCFC.push_back(neighborZ[sparseIndex]);
        nodesCFC.push_back(neighborY[neighborX[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborX[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborY[sparseIndex]]);
        nodesCFC.push_back(neighborZ[neighborY[neighborX[sparseIndex]]]);
    }

    // remove duplicate nodes
    std::sort(nodesCFC.begin(), nodesCFC.end());
    auto iterator = std::unique(nodesCFC.begin(), nodesCFC.end());
    nodesCFC.erase(iterator, nodesCFC.end());
    return nodesCFC;
}

std::vector<uint> IndexRearrangementForStreams::findIndicesNotInCommAfterFtoC(const uint numberOfIndices,
                                                                              const uint* indices,
                                                                              const std::vector<uint>& indicesAfterFtoC)
{
    std::vector<uint> otherIndices;
    for (uint posInSendIndices = 0; posInSendIndices < numberOfIndices; posInSendIndices++) {
        const uint sparseIndexSend = indices[posInSendIndices];
        if (!indexInVector(indicesAfterFtoC, sparseIndexSend))
            otherIndices.push_back(sparseIndexSend);
    }
    return otherIndices;
}


//! \}
