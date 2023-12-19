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
//! \addtogroup gpu_Parameter_tests Parameter
//! \ingroup gpu_core_tests core
//! \{
//! \author Anna Wellmann
//! \brief Functions for finding edge nodes in the multi-gpu implementation
//! \details Edge nodes are nodes, which are part of the communication in multiple directions
//! \ref master thesis of Anna Wellmann (p. 54-57)
//=======================================================================================
#include <gmock/gmock.h>

#include <filesystem>

#include "EdgeNodeFinder.h"
#include "Parameter.h"
#include "basics/config/ConfigurationFile.h"

bool compareEdgeNodesRecv(const std::vector<LBMSimulationParameter::EdgeNodePositions> &actual,
                          const std::vector<std::pair<int, int>> &expected)
{
    for (int i = 0; i < (int)expected.size(); i++) {
        if (actual[i].indexOfProcessNeighborRecv != expected[i].first) {
            return false;
        }
        if (actual[i].indexInRecvBuffer != expected[i].second) {
            return false;
        }
    }
    return true;
}

bool compareEdgeNodesSend(const std::vector<LBMSimulationParameter::EdgeNodePositions> &actual,
                          const std::vector<std::pair<int, int>> &expected)
{
    for (int i = 0; i < (int)expected.size(); i++) {
        if (actual[i].indexOfProcessNeighborSend != expected[i].first) {
            return false;
        }
        if (actual[i].indexInSendBuffer != expected[i].second) {
            return false;
        }
    }
    return true;
}

class EdgeNodeFinderTest_findEdgeNodes : public testing::Test
{
protected:
    std::shared_ptr<Parameter> para;
    const int level = 0;

private:
    void SetUp() override
    {
        para = std::make_shared<Parameter>();
        para->initLBMSimulationParameter();
    }
};

TEST_F(EdgeNodeFinderTest_findEdgeNodes, shouldReturnCorrectVectorForXY)
{
    para->parH[level]->recvProcessNeighborX.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborY.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborY.push_back(ProcessNeighbor27());

    const int numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborX.size() - 1;
    const int numSendNeighbor = (int)para->parH[level]->sendProcessNeighborY.size() - 1;

    const int sizeRecv = 6;
    const int sizeSend = 10;
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborY[numSendNeighbor].numberOfNodes = sizeSend;

    int recvNeighbors[sizeRecv] = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].index = recvNeighbors;

    int sendNeighbors[sizeSend] = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->sendProcessNeighborY[numSendNeighbor].index = sendNeighbors;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<int, int>> expectedEdgeNodesXtoYRecv = { std::pair(numRecvNeighbor, 0),
                                                                         std::pair(numRecvNeighbor, 4),
                                                                         std::pair(numRecvNeighbor, 5) };

    const std::vector<std::pair<int, int>> expectedEdgeNodesXtoYSend = { std::pair(numSendNeighbor, 1),
                                                                         std::pair(numSendNeighbor, 6),
                                                                         std::pair(numSendNeighbor, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoY.size(), testing::Eq(expectedEdgeNodesXtoYRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesXtoY, expectedEdgeNodesXtoYRecv))
        << "the edgeNodesXtoY for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesXtoY, expectedEdgeNodesXtoYSend))
        << "the edgeNodesXtoY for the send process do not match the expected nodes";
}

TEST_F(EdgeNodeFinderTest_findEdgeNodes, shouldReturnCorrectVectorForXZ)
{
    para->parH[level]->recvProcessNeighborX.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());

    const int numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborX.size() - 1;
    const int numSendNeighbor = (int)para->parH[level]->sendProcessNeighborZ.size() - 1;

    const int sizeRecv = 10;
    const int sizeSend = 6;
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborZ[numSendNeighbor].numberOfNodes = sizeSend;

    int recvNeighbors[sizeRecv] = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].index = recvNeighbors;

    int sendNeighbors[sizeSend] = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->sendProcessNeighborZ[numSendNeighbor].index = sendNeighbors;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<int, int>> expectedEdgeNodesXtoZRecv = { std::pair(numRecvNeighbor, 1),
                                                                         std::pair(numRecvNeighbor, 4),
                                                                         std::pair(numRecvNeighbor, 6) };
    const std::vector<std::pair<int, int>> expectedEdgeNodesXtoZSend = { std::pair(numSendNeighbor, 0),
                                                                         std::pair(numSendNeighbor, 5),
                                                                         std::pair(numSendNeighbor, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoZ.size(), testing::Eq(expectedEdgeNodesXtoZRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesXtoZ, expectedEdgeNodesXtoZRecv))
        << "the edgeNodesXtoZ for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesXtoZ, expectedEdgeNodesXtoZSend))
        << "the edgeNodesXtoZ for the send process do not match the expected nodes";
}

TEST_F(EdgeNodeFinderTest_findEdgeNodes, shouldReturnCorrectVectorForYZ)
{
    para->parH[level]->recvProcessNeighborY.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());

    const int sizeRecv = 10;
    const int sizeSend1 = 6;
    const int sizeSend2 = 5;

    para->parH[level]->recvProcessNeighborY[0].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborZ[0].numberOfNodes = sizeSend1;
    para->parH[level]->sendProcessNeighborZ[1].numberOfNodes = sizeSend2;

    int recvNeighbors[sizeRecv] = { 20, 1, 9, 22, 6, 23, 5, 24, 11, 26 };
    para->parH[level]->recvProcessNeighborY[0].index = recvNeighbors;

    int sendNeighbors1[sizeSend1] = { 1, 2, 3, 4, 5, 6 };
    int sendNeighbors2[sizeSend2] = { 7, 8, 9, 10, 11 };
    para->parH[level]->sendProcessNeighborZ[0].index = sendNeighbors1;
    para->parH[level]->sendProcessNeighborZ[1].index = sendNeighbors2;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<int, int>> expectedEdgeNodesYtoZRecv = { std::pair(0, 1), std::pair(0, 2),
                                                                         std::pair(0, 4), std::pair(0, 6),
                                                                         std::pair(0, 8) };
    const std::vector<std::pair<int, int>> expectedEdgeNodesYtoZSend = { std::pair(0, 0), std::pair(1, 2),
                                                                         std::pair(0, 5), std::pair(0, 4),
                                                                         std::pair(1, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesYtoZ.size(), testing::Eq(expectedEdgeNodesYtoZRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesYtoZ, expectedEdgeNodesYtoZRecv))
        << "the edgeNodesYtoZ for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesYtoZ, expectedEdgeNodesYtoZSend))
        << "the edgeNodesYtoZ for the send process do not match the expected nodes";
}

//! \}
