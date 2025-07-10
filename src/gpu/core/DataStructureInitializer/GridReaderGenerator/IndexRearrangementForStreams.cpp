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

bool indexInArray(uint* array, uint numberOfElements, uint index)
{
    return std::find(array, array + numberOfElements, index) != array + numberOfElements;
}

bool indexInVector(std::vector<uint>& vector, uint index)
{
    return std::find(vector.begin(), vector.end(), index) != vector.end();
}

IndexRearrangementForStreams::IndexRearrangementForStreams(std::shared_ptr<Parameter> para,
                                                           std::shared_ptr<GridBuilder> builder,
                                                           vf::parallel::Communicator& communicator)
    : para(para), builder(builder), communicator(communicator)
{
}

void IndexRearrangementForStreams::initCommunicationArraysForCommAfterFinetoCoarse(
    ProcessNeighbor27& sendNeighborHost, ProcessNeighbor27& sendNeighborDevice, ProcessNeighbor27& sendNeighborAfterFtoCHost,
    ProcessNeighbor27& sendNeighborAfterFtoCDevice, ProcessNeighbor27& recvNeighborHost,
    ProcessNeighbor27& recvNeighborDevice, ProcessNeighbor27& recvNeighborAfterFtoCHost,
    ProcessNeighbor27& recvNeighborAfterFtoCDevice, const int level, const int direction) const
{
    VF_LOG_INFO("Communication: reorder send indices in direction {}", direction);

    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    const uint numberOfNodesSend =
        reorderSendIndicesForCommAfterFtoC(sendNeighborHost, direction, level, sendIndicesForCommAfterFtoCPositions);
    setNumberOfNodes(sendNeighborAfterFtoCHost, sendNeighborAfterFtoCDevice, numberOfNodesSend);

    const std::vector<uint> recvIndicesForCommAfterFtoCPositions = exchangeIndicesForCommAfterFtoC(
        sendNeighborAfterFtoCHost, recvNeighborAfterFtoCHost, sendIndicesForCommAfterFtoCPositions);

    const uint numberOfNodesRecv =
        reorderRecvIndicesForCommAfterFtoC(recvNeighborHost, direction, level, recvIndicesForCommAfterFtoCPositions);
    setNumberOfNodes(recvNeighborAfterFtoCHost, recvNeighborAfterFtoCDevice, numberOfNodesRecv);

    copyProcessNeighborToCommAfterFtoC(recvNeighborHost, recvNeighborAfterFtoCHost);
    copyProcessNeighborToCommAfterFtoC(sendNeighborHost, sendNeighborAfterFtoCHost);
    copyProcessNeighborToCommAfterFtoC(recvNeighborDevice, recvNeighborAfterFtoCDevice);
    copyProcessNeighborToCommAfterFtoC(sendNeighborDevice, sendNeighborAfterFtoCDevice);
}

std::vector<uint>
IndexRearrangementForStreams::exchangeIndicesForCommAfterFtoC(ProcessNeighbor27& sendNeighbor,
                                                              ProcessNeighbor27& recvNeighbor,
                                                              std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
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

void IndexRearrangementForStreams::copyProcessNeighborToCommAfterFtoC(ProcessNeighbor27& neighbor,
                                                                      ProcessNeighbor27& neighborAfterFtoC) const
{
    neighborAfterFtoC.populations[0] = neighbor.populations[0];
    neighborAfterFtoC.index = neighbor.index;
    neighborAfterFtoC.rankNeighbor = neighbor.rankNeighbor;
}

void IndexRearrangementForStreams::setNumberOfNodes(ProcessNeighbor27& neighborAfterFtoCHost,
                                                    ProcessNeighbor27& neighborAfterFtoCDevice,
                                                    const uint numberOfNodes) const
{
    neighborAfterFtoCHost.numberOfNodes = numberOfNodes;
    neighborAfterFtoCHost.numberOfFs = 27 * numberOfNodes;
    neighborAfterFtoCHost.memsizeFs = sizeof(real) * 27 * numberOfNodes;

    neighborAfterFtoCDevice.numberOfNodes = numberOfNodes;
    neighborAfterFtoCDevice.numberOfFs = 27 * numberOfNodes;
    neighborAfterFtoCDevice.memsizeFs = sizeof(real) * 27 * numberOfNodes;
}

uint IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoC(
    ProcessNeighbor27& sendNeighbor, int direction, int level, std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
{
    ICells& fineToCoarse = para->getParH(level)->fineToCoarse;
    VF_LOG_INFO("Reorder send indices for communication after fine to coarse: level: {} direction: {}", level, direction);
    if (para->getParH(level)->coarseToFine.numberOfCells == 0 || fineToCoarse.numberOfCells == 0)
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): para->getParH(level)->intCF needs to be initialized "
                        "before calling this function");
    std::vector<uint> sendIndicesAfterFtoC;
    std::vector<uint> sendIndicesOther;
    const uint numberOfSendIndices = builder->getNumberOfSendIndices(direction, level);

    // coarse cells of interpolation fine to coarse (iCellFCC)
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        const uint sparseIndexSend = sendNeighbor.index[posInSendIndices];
        if (indexInArray(fineToCoarse.coarseCellIndices, fineToCoarse.numberOfCells, sparseIndexSend)) {
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndexSend, sendIndicesForCommAfterFtoCPositions,
                                                 posInSendIndices);
        }
    }

    // coarse cells of interpolation coarse to fine (iCellCFC)
    std::vector<uint> coarseCellsForCtoF;
    aggregateCoarseNodesForCtoF(level, coarseCellsForCtoF);
    for (auto sparseIndex : coarseCellsForCtoF)
        findIfSparseIndexIsInSendIndicesAndAddToCommVectors(sparseIndex, sendNeighbor.index, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);

    const uint numberOfNodes = (uint)sendIndicesAfterFtoC.size();

    findIndicesNotInCommAfterFtoC(numberOfSendIndices, sendNeighbor.index, sendIndicesAfterFtoC, sendIndicesOther);

    std::copy_n(sendIndicesAfterFtoC.begin(), numberOfNodes, sendNeighbor.index);
    std::copy_n(sendIndicesOther.begin(), sendIndicesOther.size(), sendNeighbor.index + numberOfNodes);

    VF_LOG_INFO("Reorder send indices: process {}, numberOfSendNodesAfterFtoC {}", communicator.getProcessID(),
                numberOfNodes);

    if (numberOfNodes + sendIndicesOther.size() != numberOfSendIndices) {
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfSendNodesAfterFtoC = {}, sendIndicesOther.size() = {}, numberOfSendIndices = {}",
                        numberOfNodes, sendIndicesOther.size(), numberOfSendIndices);
    }
    return numberOfNodes;
}

void IndexRearrangementForStreams::aggregateCoarseNodesForCtoF(int level, std::vector<uint>& nodesCFC) const
{
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
}

void IndexRearrangementForStreams::addUniqueIndexToCommunicationVectors(
    std::vector<uint>& sendIndicesAfterFtoC, uint sparseIndexSend, std::vector<uint>& sendIndicesForCommAfterFtoCPositions,
    uint posInSendIndices) const
{
    // add index to corresponding vectors, but omit indices which are already in sendIndicesAfterFtoC
    if (!indexInVector(sendIndicesAfterFtoC, sparseIndexSend)) {
        sendIndicesAfterFtoC.push_back(sparseIndexSend);
        sendIndicesForCommAfterFtoCPositions.push_back(posInSendIndices);
    }
}

void IndexRearrangementForStreams::findIfSparseIndexIsInSendIndicesAndAddToCommVectors(
    uint sparseIndex, const uint* sendIndices, uint numberOfSendIndices, std::vector<uint>& sendIndicesAfterFtoC,
    std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
{
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        if (sparseIndex == sendIndices[posInSendIndices]) {
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndex, sendIndicesForCommAfterFtoCPositions,
                                                 posInSendIndices);
            break;
        }
    }
}

void IndexRearrangementForStreams::findIndicesNotInCommAfterFtoC(const uint& numberOfSendOrRecvIndices,
                                                                 const uint* sendOrReceiveIndices,
                                                                 std::vector<uint>& sendOrReceiveIndicesAfterFtoC,
                                                                 std::vector<uint>& sendOrIndicesOther)
{
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendOrRecvIndices; posInSendIndices++) {
        const uint sparseIndexSend = sendOrReceiveIndices[posInSendIndices];
        if (!indexInVector(sendOrReceiveIndicesAfterFtoC, sparseIndexSend))
            sendOrIndicesOther.push_back(sparseIndexSend);
    }
}

uint IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoC(
    ProcessNeighbor27& recvNeighbor, int direction, int level,
    const std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
{
    VF_LOG_INFO("Reorder recv indices for communication after fine to coarse: level: {} direction: {}", level, direction);

    if (sendIndicesForCommAfterFtoCPositions.empty())
        VF_LOG_WARNING("ReorderRecvIndicesForCommAfterFtoC(): sendIndicesForCommAfterFtoCPositions is empty.");

    uint numberOfRecvIndices = builder->getNumberOfReceiveIndices(direction, level);
    std::vector<uint> recvIndicesAfterFtoC;
    std::vector<uint> recvIndicesOther;

    // find recvIndices for Communication after fine to coarse
    for (uint vectorPos : sendIndicesForCommAfterFtoCPositions)
        recvIndicesAfterFtoC.push_back(recvNeighbor.index[vectorPos]);

    findIndicesNotInCommAfterFtoC(numberOfRecvIndices, recvNeighbor.index, recvIndicesAfterFtoC, recvIndicesOther);

    const uint numberOfNodes = static_cast<uint>(recvIndicesAfterFtoC.size());

    // copy new vectors back to recvIndices array
    std::copy_n(recvIndicesAfterFtoC.begin(), numberOfNodes, recvNeighbor.index);
    std::copy_n(recvIndicesOther.begin(), recvIndicesOther.size(), recvNeighbor.index + numberOfNodes);

    VF_LOG_INFO("Reorder recv indices: process {}, numberOfRecvNodesAfterFtoC {}", communicator.getProcessID(),
                numberOfNodes);

    if (numberOfNodes + recvIndicesOther.size() != numberOfRecvIndices) {
        VF_LOG_CRITICAL("reorderRecvIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfRecvNodesAfterFtoC = {}, recvIndicesOther.size() = {}, numberOfRecvIndices = {}",
                        numberOfNodes, recvIndicesOther.size(), numberOfRecvIndices);
    }
    return numberOfNodes;
}

//! \}
