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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef EDGENODEDEBUG_HPP
#define EDGENODEDEBUG_HPP

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <basics/StringUtilities/StringUtil.h>
#include <basics/utilities/UbSystem.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include <lbm/constants/D3Q27.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

namespace edge_node_debug_writer
{

void addCoordinatesToNodeVector(LBMSimulationParameter* parH, std::vector<UbTupleFloat3>& nodesVec, uint indexInNodesVector,
                                uint sparseIndexOfNode)
{
    nodesVec[indexInNodesVector] =
        makeUbTuple(float(parH->coordinateX[sparseIndexOfNode]), float(parH->coordinateY[sparseIndexOfNode]),
                    float(parH->coordinateZ[sparseIndexOfNode]));
}

void writeEdgeNodesXZ_Send(SPtr<Parameter>& para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<std::string> datanames = { "SparseIndex", "ProcessNeighbor", "IndexInSendVector", "AfterFtoC" };
    std::vector<std::vector<double>> nodedata;

    uint numberOfNodes = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        numberOfNodes += (uint)para->getParH(level)->edgeNodesXtoZ.size();
    }
    nodesVec.resize(numberOfNodes);
    nodedata.resize(datanames.size(), std::vector<double>(numberOfNodes));

    uint nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (uint u = 0; u < numberOfNodes; u++) {
            const uint indexOfProcessNeighborSend = para->getParH(level)->edgeNodesXtoZ[u].indexOfProcessNeighborSend;
            const uint indexInSendBuffer = para->getParH(level)->edgeNodesXtoZ[u].indexInSendBuffer;
            const uint sparseIndex =
                para->getParH(level)->sendProcessNeighborsZ[indexOfProcessNeighborSend].index[indexInSendBuffer];
            nodedata[0][nodeCount] = sparseIndex;
            nodedata[1][nodeCount] = indexOfProcessNeighborSend;
            nodedata[2][nodeCount] = indexInSendBuffer;
            nodedata[3][nodeCount] =
                double(indexInSendBuffer <
                       para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighborSend].numberOfNodes);

            addCoordinatesToNodeVector(para->getParH(level).get(), nodesVec, nodeCount, sparseIndex);

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeEdgeNodesXZ_Send_PID_" + std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

void writeEdgeNodesXZ_Recv(SPtr<Parameter>& para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<std::string> datanames = { "SparseIndex", "ProcessNeighbor", "IndexInRecvVector", "AfterFtoC" };
    std::vector<std::vector<double>> nodedata;

    uint numberOfNodes = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        numberOfNodes += (uint)para->getParH(level)->edgeNodesXtoZ.size();
    }
    nodesVec.resize(numberOfNodes);
    nodedata.resize(datanames.size(), std::vector<double>(numberOfNodes));

    uint nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (uint u = 0; u < numberOfNodes; u++) {
            uint indexOfProcessNeighborRecv = para->getParH(level)->edgeNodesXtoZ[u].indexOfProcessNeighborRecv;
            uint indexInRecvBuffer = para->getParH(level)->edgeNodesXtoZ[u].indexInRecvBuffer;
            uint sparseIndex =
                para->getParH(level)->recvProcessNeighborsX[indexOfProcessNeighborRecv].index[indexInRecvBuffer];
            nodedata[0][nodeCount] = sparseIndex;
            nodedata[1][nodeCount] = indexOfProcessNeighborRecv;
            nodedata[2][nodeCount] = indexInRecvBuffer;
            nodedata[3][nodeCount] = double(
                indexInRecvBuffer < para->getParH(level)->recvProcessNeighborsX[indexOfProcessNeighborRecv].numberOfNodes);

            addCoordinatesToNodeVector(para->getParH(level).get(), nodesVec, nodeCount, sparseIndex);

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeEdgeNodesXZ_Recv_PID_" + std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}
} // namespace edge_node_debug_writer

#endif

//! \}
