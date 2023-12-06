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
//! \author Anna Wellmann
//! \details See [master thesis of Anna Wellmann]
//=======================================================================================
#include "IndexRearrangementForStreams.h"

#include <algorithm>
#include <iostream>

#include <logger/Logger.h>

#include <GridGenerator/grid/Grid.h>
#include <GridGenerator/grid/GridBuilder/GridBuilder.h>

#include <parallel/Communicator.h>

#include "Parameter/Parameter.h"

IndexRearrangementForStreams::IndexRearrangementForStreams(std::shared_ptr<Parameter> para,
                                                           std::shared_ptr<GridBuilder> builder,
                                                           vf::parallel::Communicator &communicator)
    : para(para), builder(builder), communicator(communicator)
{
}

void IndexRearrangementForStreams::initCommunicationArraysForCommAfterFinetoCoarseX(uint level,
                                                                                    int indexOfProcessNeighbor,
                                                                                    int direction) const
{
    VF_LOG_INFO("Communication: reorder send indices in x direction");
    std::vector<uint> sendIndicesForCommAfterFtoCPositions =
        initSendIndicesForCommAfterFToCX(level, indexOfProcessNeighbor, direction);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions =
        exchangeIndicesForCommAfterFtoCX(level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);

    initRecvIndicesForCommAfterFToCX(level, indexOfProcessNeighbor, direction, recvIndicesForCommAfterFtoCPositions);

    copyProcessNeighborToCommAfterFtoCX(level, indexOfProcessNeighbor);
}

void IndexRearrangementForStreams::initCommunicationArraysForCommAfterFinetoCoarseY(uint level,
                                                                                    int indexOfProcessNeighbor,
                                                                                    int direction) const
{
    VF_LOG_INFO("Communication: reorder send indices in x direction");
    std::vector<uint> sendIndicesForCommAfterFtoCPositions =
        initSendIndicesForCommAfterFToCY(level, indexOfProcessNeighbor, direction);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions =
        exchangeIndicesForCommAfterFtoCY(level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);

    initRecvIndicesForCommAfterFToCY(level, indexOfProcessNeighbor, direction, recvIndicesForCommAfterFtoCPositions);

    copyProcessNeighborToCommAfterFtoCY(level, indexOfProcessNeighbor);
}

void IndexRearrangementForStreams::initCommunicationArraysForCommAfterFinetoCoarseZ(uint level,
                                                                                    int indexOfProcessNeighbor,
                                                                                    int direction) const
{
    VF_LOG_INFO("Communication: reorder send indices in z direction");
    std::vector<uint> sendIndicesForCommAfterFtoCPositions =
        initSendIndicesForCommAfterFToCZ(level, indexOfProcessNeighbor, direction);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions =
        exchangeIndicesForCommAfterFtoCZ(level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);

    initRecvIndicesForCommAfterFToCZ(level, indexOfProcessNeighbor, direction, recvIndicesForCommAfterFtoCPositions);

    copyProcessNeighborToCommAfterFtoCZ(level, indexOfProcessNeighbor);
}

std::vector<uint> IndexRearrangementForStreams::initSendIndicesForCommAfterFToCX(uint level, int indexOfProcessNeighbor,
                                                                                 int direction) const
{
    para->initProcessNeighborsAfterFtoCX(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCX(direction, level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCX(
        para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
    return sendIndicesForCommAfterFtoCPositions;
}

std::vector<uint> IndexRearrangementForStreams::initSendIndicesForCommAfterFToCY(uint level, int indexOfProcessNeighbor,
                                                                                 int direction) const
{
    para->initProcessNeighborsAfterFtoCY(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCY(direction, level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCY(
        para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
    return sendIndicesForCommAfterFtoCPositions;
}

std::vector<uint> IndexRearrangementForStreams::initSendIndicesForCommAfterFToCZ(uint level, int indexOfProcessNeighbor,
                                                                                 int direction) const
{
    para->initProcessNeighborsAfterFtoCZ(level);
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    reorderSendIndicesForCommAfterFtoCZ(direction, level, indexOfProcessNeighbor, sendIndicesForCommAfterFtoCPositions);
    para->setSendProcessNeighborsAfterFtoCZ(
        para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
    return sendIndicesForCommAfterFtoCPositions;
}

std::vector<uint> IndexRearrangementForStreams::exchangeIndicesForCommAfterFtoCX(
    uint level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    // fill the receive vector with zeros as placeholders
    // give vector an arbitrary size (larger than needed) // TODO: Find a better way
    std::vector<uint> recvIndicesForCommAfterFtoCPositions(
        (size_t)para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes * 2, 0);

    communicator.receiveSend(
        recvIndicesForCommAfterFtoCPositions.data(), (int)recvIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].rankNeighbor,
        sendIndicesForCommAfterFtoCPositions.data(), (int)sendIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->sendProcessNeighborX[indexOfProcessNeighbor].rankNeighbor);

    // resize receiving vector to correct size
    if ((uint)recvIndicesForCommAfterFtoCPositions.size() > 0) {
        auto it = std::unique(
            recvIndicesForCommAfterFtoCPositions.begin(),
            recvIndicesForCommAfterFtoCPositions.end()); // finds the second zero when there are multiple zeros in a row
        recvIndicesForCommAfterFtoCPositions.erase(
            std::prev(it, 1),                            // begin erasing at the first zero
            recvIndicesForCommAfterFtoCPositions.end()); // TODO: Find a better way
    }

    return recvIndicesForCommAfterFtoCPositions;
}

std::vector<uint> IndexRearrangementForStreams::exchangeIndicesForCommAfterFtoCY(
    uint level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    // fill the receive vector with zeros as placeholders
    // give vector an arbitrary size (larger than needed) // TODO: Find a better way
    std::vector<uint> recvIndicesForCommAfterFtoCPositions(
        (size_t)para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes * 2, 0);

    communicator.receiveSend(
        recvIndicesForCommAfterFtoCPositions.data(), (int)recvIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].rankNeighbor,
        sendIndicesForCommAfterFtoCPositions.data(), (int)sendIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->sendProcessNeighborY[indexOfProcessNeighbor].rankNeighbor);

    // resize receiving vector to correct size
    if ((uint)recvIndicesForCommAfterFtoCPositions.size() > 0) {
        auto it = std::unique(
            recvIndicesForCommAfterFtoCPositions.begin(),
            recvIndicesForCommAfterFtoCPositions.end()); // finds the second zero when there are multiple zeros in a row
        recvIndicesForCommAfterFtoCPositions.erase(
            std::prev(it, 1),                            // begin erasing at the first zero
            recvIndicesForCommAfterFtoCPositions.end()); // TODO: Find a better way
    }

    return recvIndicesForCommAfterFtoCPositions;
}

std::vector<uint> IndexRearrangementForStreams::exchangeIndicesForCommAfterFtoCZ(
    uint level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    // fill the receive vector with zeros as placeholders
    // give vector an arbitrary size (larger than needed) // TODO: Find a better way
    std::vector<uint> recvIndicesForCommAfterFtoCPositions(
        (size_t)para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes * 2, 0);

    communicator.receiveSend(
        recvIndicesForCommAfterFtoCPositions.data(), (int)recvIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].rankNeighbor,
        sendIndicesForCommAfterFtoCPositions.data(), (int)sendIndicesForCommAfterFtoCPositions.size(),
        para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighbor].rankNeighbor);

    // resize receiving vector to correct size
    if ((uint)recvIndicesForCommAfterFtoCPositions.size() > 0) {
        auto it = std::unique(
            recvIndicesForCommAfterFtoCPositions.begin(),
            recvIndicesForCommAfterFtoCPositions.end()); // finds the second zero when there are multiple zeros in a row
        recvIndicesForCommAfterFtoCPositions.erase(
            std::prev(it, 1),                            // begin erasing at the first zero
            recvIndicesForCommAfterFtoCPositions.end()); // TODO: Find a better way
    }

    return recvIndicesForCommAfterFtoCPositions;
}

void IndexRearrangementForStreams::initRecvIndicesForCommAfterFToCX(
    uint level, int indexOfProcessNeighbor, int direction,
    std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const
{
    reorderRecvIndicesForCommAfterFtoCX(direction, level, indexOfProcessNeighbor, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCX(
        para->getParH(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
}

void IndexRearrangementForStreams::initRecvIndicesForCommAfterFToCY(
    uint level, int indexOfProcessNeighbor, int direction,
    std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const
{
    reorderRecvIndicesForCommAfterFtoCY(direction, level, indexOfProcessNeighbor, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCY(
        para->getParH(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
}
void IndexRearrangementForStreams::initRecvIndicesForCommAfterFToCZ(
    uint level, int indexOfProcessNeighbor, int direction,
    std::vector<uint> &recvIndicesForCommAfterFtoCPositions) const
{
    reorderRecvIndicesForCommAfterFtoCZ(direction, level, indexOfProcessNeighbor, recvIndicesForCommAfterFtoCPositions);
    para->setRecvProcessNeighborsAfterFtoCZ(
        para->getParH(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes, level,
        indexOfProcessNeighbor);
}

void IndexRearrangementForStreams::copyProcessNeighborToCommAfterFtoCX(uint level, int indexOfProcessNeighbor) const
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->sendProcessNeighborX[indexOfProcessNeighbor].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->sendProcessNeighborX[indexOfProcessNeighbor].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->recvProcessNeighborX[indexOfProcessNeighbor].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].index =
        para->getParD(level)->sendProcessNeighborX[indexOfProcessNeighbor].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].index =
        para->getParH(level)->sendProcessNeighborX[indexOfProcessNeighbor].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].index =
        para->getParD(level)->recvProcessNeighborX[indexOfProcessNeighbor].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].index =
        para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->sendProcessNeighborX[indexOfProcessNeighbor].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].rankNeighbor;
}

void IndexRearrangementForStreams::copyProcessNeighborToCommAfterFtoCY(uint level, int indexOfProcessNeighbor) const
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->sendProcessNeighborY[indexOfProcessNeighbor].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->sendProcessNeighborY[indexOfProcessNeighbor].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->recvProcessNeighborY[indexOfProcessNeighbor].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].index =
        para->getParD(level)->sendProcessNeighborY[indexOfProcessNeighbor].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].index =
        para->getParH(level)->sendProcessNeighborY[indexOfProcessNeighbor].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].index =
        para->getParD(level)->recvProcessNeighborY[indexOfProcessNeighbor].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].index =
        para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->sendProcessNeighborY[indexOfProcessNeighbor].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].rankNeighbor;
}

void IndexRearrangementForStreams::copyProcessNeighborToCommAfterFtoCZ(uint level, int indexOfProcessNeighbor) const
{
    // init f[0]*
    para->getParD(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->sendProcessNeighborZ[indexOfProcessNeighbor].f[0];
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighbor].f[0];
    para->getParD(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].f[0] =
        para->getParD(level)->recvProcessNeighborZ[indexOfProcessNeighbor].f[0];
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].f[0] =
        para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].f[0];

    // init index*
    para->getParD(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].index =
        para->getParD(level)->sendProcessNeighborZ[indexOfProcessNeighbor].index;
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].index =
        para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighbor].index;
    para->getParD(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].index =
        para->getParD(level)->recvProcessNeighborZ[indexOfProcessNeighbor].index;
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].index =
        para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].index;

    // rank neighbor
    para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighbor].rankNeighbor;
    para->getParH(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].rankNeighbor =
        para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].rankNeighbor;
}

void IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoCX(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *sendIndices = para->getParH(level)->sendProcessNeighborX[indexOfProcessNeighbor].index;
    int &numberOfSendNodesAfterFtoC =
        para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoCY(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *sendIndices = para->getParH(level)->sendProcessNeighborY[indexOfProcessNeighbor].index;
    int &numberOfSendNodesAfterFtoC =
        para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoCZ(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *sendIndices = para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighbor].index;
    int &numberOfSendNodesAfterFtoC =
        para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes;
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderSendIndicesForCommAfterFtoC(
    int *sendIndices, int &numberOfSendNodesAfterFtoC, int direction, int level,
    std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    VF_LOG_INFO("Reorder send indices for communication after fine to coarse: level: {} direction: {}", level,
                direction);
    if (para->getParH(level)->coarseToFine.numberOfCells == 0 || para->getParH(level)->fineToCoarse.numberOfCells == 0)
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): para->getParH(level)->intCF needs to be initialized "
                        "before calling this function");

    int sparseIndexSend;
    std::vector<int> sendIndicesAfterFtoC;
    std::vector<int> sendIndicesOther;
    uint numberOfSendIndices = builder->getNumberOfSendIndices(direction, level);

    // coarse cells of interpolation fine to coarse (iCellFCC)
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        sparseIndexSend = sendIndices[posInSendIndices];
        if (isSparseIndexInCoarseIndexForFtoC(para->getParH(level)->fineToCoarse.numberOfCells, sparseIndexSend, level)) {
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndexSend,
                                                 sendIndicesForCommAfterFtoCPositions, posInSendIndices);
        }
    }

    // coarse cells of interpolation coarse to fine (iCellCFC)
    std::vector<uint> coarseCellsForCtoF;
    aggregateCoarseNodesForCtoF(level, coarseCellsForCtoF);
    for (auto sparseIndex : coarseCellsForCtoF)
        findIfSparseIndexIsInSendIndicesAndAddToCommVectors(sparseIndex, sendIndices, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);

    numberOfSendNodesAfterFtoC = (int)sendIndicesAfterFtoC.size();

    findIndicesNotInCommAfterFtoC(numberOfSendIndices, sendIndices, sendIndicesAfterFtoC, sendIndicesOther);

    // copy new vectors back to sendIndices array
    for (int i = 0; i < numberOfSendNodesAfterFtoC; i++)
        sendIndices[i] = sendIndicesAfterFtoC[i];
    for (uint i = 0; i < (uint)sendIndicesOther.size(); i++)
        sendIndices[i + numberOfSendNodesAfterFtoC] = sendIndicesOther[i];

    VF_LOG_INFO("Reorder send indices: process {}, numberOfSendNodesAfterFtoC {}", communicator.getProcessID(),
                numberOfSendNodesAfterFtoC);

    if (numberOfSendNodesAfterFtoC + sendIndicesOther.size() != numberOfSendIndices) {
        VF_LOG_CRITICAL("reorderSendIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfSendNodesAfterFtoC = {}, sendIndicesOther.size() = {}, numberOfSendIndices = {}",
                        numberOfSendNodesAfterFtoC, sendIndicesOther.size(), numberOfSendIndices);
    }
}

bool IndexRearrangementForStreams::isSparseIndexInCoarseIndexForFtoC(uint numberOfCoarseNodesForFtoC, int sparseIndex, int level) const
{
    for (uint j = 0; j < numberOfCoarseNodesForFtoC; j++) {
        if (sparseIndex < 0)
            return false;
        if (para->getParH(level)->fineToCoarse.coarseCellIndices[j] == (uint)sparseIndex) {
            return true;
        }
    }
    return false;
}

void IndexRearrangementForStreams::aggregateCoarseNodesForCtoF(int level, std::vector<uint> &nodesCFC) const
{
    uint sparseIndex;
    uint *neighborX = para->getParH(level)->neighborX;
    uint *neighborY = para->getParH(level)->neighborY;
    uint *neighborZ = para->getParH(level)->neighborZ;

    for (uint x = 0; x < para->getParH(level)->coarseToFine.numberOfCells; x++) {
        sparseIndex = para->getParH(level)->coarseToFine.coarseCellIndices[x];
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
    std::vector<int> &sendIndicesAfterFtoC, int &sparseIndexSend,
    std::vector<unsigned int> &sendIndicesForCommAfterFtoCPositions, uint &posInSendIndices) const
{
    // add index to corresponding vectors, but omit indices which are already in sendIndicesAfterFtoC
    if (std::find(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), sparseIndexSend) ==
        sendIndicesAfterFtoC.end()) {
        sendIndicesAfterFtoC.push_back(sparseIndexSend);
        sendIndicesForCommAfterFtoCPositions.push_back(posInSendIndices);
    }
}

void IndexRearrangementForStreams::findIfSparseIndexIsInSendIndicesAndAddToCommVectors(
    int sparseIndex, int *sendIndices, uint numberOfSendIndices, std::vector<int> &sendIndicesAfterFtoC,
    std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int sparseIndexSend;
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        sparseIndexSend = sendIndices[posInSendIndices];
        if (sparseIndex == sparseIndexSend) {
            addUniqueIndexToCommunicationVectors(sendIndicesAfterFtoC, sparseIndex,
                                                 sendIndicesForCommAfterFtoCPositions, posInSendIndices);
            break;
        }
    }
}

void IndexRearrangementForStreams::findIndicesNotInCommAfterFtoC(const uint &numberOfSendOrRecvIndices,
                                                                 int *sendOrReceiveIndices,
                                                                 std::vector<int> &sendOrReceiveIndicesAfterFtoC,
                                                                 std::vector<int> &sendOrIndicesOther)
{
    int sparseIndexSend;
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendOrRecvIndices; posInSendIndices++) {
        sparseIndexSend = sendOrReceiveIndices[posInSendIndices];
        if (std::find(sendOrReceiveIndicesAfterFtoC.begin(), sendOrReceiveIndicesAfterFtoC.end(), sparseIndexSend) ==
            sendOrReceiveIndicesAfterFtoC.end())
            sendOrIndicesOther.push_back(sparseIndexSend);
    }
}

void IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoCX(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *recvIndices = para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].index;
    int &numberOfRecvNodesAfterFtoC =
        para->getParH(level)->recvProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoCY(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *recvIndices = para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].index;
    int &numberOfRecvNodesAfterFtoC =
        para->getParH(level)->recvProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoCZ(
    int direction, int level, int indexOfProcessNeighbor, std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    int *recvIndices = para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].index;
    int &numberOfRecvNodesAfterFtoC =
        para->getParH(level)->recvProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes;
    reorderRecvIndicesForCommAfterFtoC(recvIndices, numberOfRecvNodesAfterFtoC, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
}

void IndexRearrangementForStreams::reorderRecvIndicesForCommAfterFtoC(
    int *recvIndices, int &numberOfRecvNodesAfterFtoC, int direction, int level,
    std::vector<uint> &sendIndicesForCommAfterFtoCPositions) const
{
    VF_LOG_INFO("Reorder recv indices for communication after fine to coarse: level: {} direction: {}", level,
                direction);

    if (sendIndicesForCommAfterFtoCPositions.size() == 0)
        VF_LOG_WARNING("ReorderRecvIndicesForCommAfterFtoC(): sendIndicesForCommAfterFtoCPositions is empty.");

    uint numberOfRecvIndices = builder->getNumberOfReceiveIndices(direction, level);
    std::vector<int> recvIndicesAfterFtoC;
    std::vector<int> recvIndicesOther;

    // find recvIndices for Communication after fine to coarse
    for (uint vectorPos : sendIndicesForCommAfterFtoCPositions)
        recvIndicesAfterFtoC.push_back(recvIndices[vectorPos]);

    findIndicesNotInCommAfterFtoC(numberOfRecvIndices, recvIndices, recvIndicesAfterFtoC, recvIndicesOther);

    numberOfRecvNodesAfterFtoC = (int)recvIndicesAfterFtoC.size();

    // copy new vectors back to recvIndices array
    for (int i = 0; i < numberOfRecvNodesAfterFtoC; i++)
        recvIndices[i] = recvIndicesAfterFtoC[i];
    for (uint i = 0; i < (uint)recvIndicesOther.size(); i++)
        recvIndices[i + numberOfRecvNodesAfterFtoC] = recvIndicesOther[i];

    VF_LOG_INFO("Reorder send indices: process {}, numberOfRecvNodesAfterFtoC {}", communicator.getProcessID(),
                numberOfRecvNodesAfterFtoC);

    if (numberOfRecvNodesAfterFtoC + recvIndicesOther.size() != numberOfRecvIndices) {
        VF_LOG_CRITICAL("reorderRecvIndicesForCommAfterFtoC(): incorrect number of nodes");
        VF_LOG_CRITICAL("numberOfRecvNodesAfterFtoC = {}, recvIndicesOther.size() = {}, numberOfRecvIndices = {}",
                        numberOfRecvNodesAfterFtoC, recvIndicesOther.size(), numberOfRecvIndices);
    }
}
