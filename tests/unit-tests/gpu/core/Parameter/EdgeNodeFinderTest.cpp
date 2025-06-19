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
//=======================================================================================
#include <gmock/gmock.h>

#include <filesystem>

#include <basics/config/ConfigurationFile.h>

#include <gpu/core/Parameter/EdgeNodeFinder.h>
#include <gpu/core/Parameter/Parameter.h>

bool compareEdgeNodesRecv(const std::vector<LBMSimulationParameter::EdgeNodePositions> &actual,
                          const std::vector<std::pair<uint, uint>> &expected)
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
                          const std::vector<std::pair<uint, uint>> &expected)
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
    para->parH[level]->recvProcessNeighborsX.emplace_back();
    para->parH[level]->sendProcessNeighborsY.emplace_back();
    para->parH[level]->sendProcessNeighborsY.emplace_back();
    const uint numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborsX.size() - 1;
    const uint numSendNeighbor = (int)para->parH[level]->sendProcessNeighborsY.size() - 1;

    const uint sizeRecv = 6;
    const uint sizeSend = 10;
    para->parH[level]->recvProcessNeighborsX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborsY[numSendNeighbor].numberOfNodes = sizeSend;

    uint recvNeighbors[sizeRecv] = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->recvProcessNeighborsX[numRecvNeighbor].index = recvNeighbors;

    uint sendNeighbors[sizeSend] = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->sendProcessNeighborsY[numSendNeighbor].index = sendNeighbors;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<uint, uint>> expectedEdgeNodesXtoYRecv = { std::pair(numRecvNeighbor, 0U),
                                                                         std::pair(numRecvNeighbor, 4U),
                                                                         std::pair(numRecvNeighbor, 5U) };

    const std::vector<std::pair<uint, uint>> expectedEdgeNodesXtoYSend = { std::pair(numSendNeighbor, 1U),
                                                                         std::pair(numSendNeighbor, 6U),
                                                                         std::pair(numSendNeighbor, 4U) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoY.size(), testing::Eq(expectedEdgeNodesXtoYRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesXtoY, expectedEdgeNodesXtoYRecv))
        << "the edgeNodesXtoY for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesXtoY, expectedEdgeNodesXtoYSend))
        << "the edgeNodesXtoY for the send process do not match the expected nodes";
}

TEST_F(EdgeNodeFinderTest_findEdgeNodes, shouldReturnCorrectVectorForXZ)
{
    para->parH[level]->recvProcessNeighborsX.emplace_back();
    para->parH[level]->sendProcessNeighborsZ.emplace_back();
    para->parH[level]->sendProcessNeighborsZ.emplace_back();

    const uint numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborsX.size() - 1;
    const uint numSendNeighbor = (int)para->parH[level]->sendProcessNeighborsZ.size() - 1;

    const uint sizeRecv = 10;
    const uint sizeSend = 6;
    para->parH[level]->recvProcessNeighborsX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborsZ[numSendNeighbor].numberOfNodes = sizeSend;

    uint recvNeighbors[sizeRecv] = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->recvProcessNeighborsX[numRecvNeighbor].index = recvNeighbors;

    uint sendNeighbors[sizeSend] = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->sendProcessNeighborsZ[numSendNeighbor].index = sendNeighbors;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<uint, uint>> expectedEdgeNodesXtoZRecv = { std::pair(numRecvNeighbor, 1U),
                                                                         std::pair(numRecvNeighbor, 4U),
                                                                         std::pair(numRecvNeighbor, 6U) };
    const std::vector<std::pair<uint, uint>> expectedEdgeNodesXtoZSend = { std::pair(numSendNeighbor, 0U),
                                                                         std::pair(numSendNeighbor, 5U),
                                                                         std::pair(numSendNeighbor, 4U) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoZ.size(), testing::Eq(expectedEdgeNodesXtoZRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesXtoZ, expectedEdgeNodesXtoZRecv))
        << "the edgeNodesXtoZ for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesXtoZ, expectedEdgeNodesXtoZSend))
        << "the edgeNodesXtoZ for the send process do not match the expected nodes";
}

TEST_F(EdgeNodeFinderTest_findEdgeNodes, shouldReturnCorrectVectorForYZ)
{
    para->parH[level]->recvProcessNeighborsY.emplace_back();
    para->parH[level]->sendProcessNeighborsZ.emplace_back();
    para->parH[level]->sendProcessNeighborsZ.emplace_back();

    const uint sizeRecv = 10;
    const uint sizeSend1 = 6;
    const uint sizeSend2 = 5;

    para->parH[level]->recvProcessNeighborsY[0].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborsZ[0].numberOfNodes = sizeSend1;
    para->parH[level]->sendProcessNeighborsZ[1].numberOfNodes = sizeSend2;

    uint recvNeighbors[sizeRecv] = { 20, 1, 9, 22, 6, 23, 5, 24, 11, 26 };
    para->parH[level]->recvProcessNeighborsY[0].index = recvNeighbors;

    uint sendNeighbors1[sizeSend1] = { 1, 2, 3, 4, 5, 6 };
    uint sendNeighbors2[sizeSend2] = { 7, 8, 9, 10, 11 };
    para->parH[level]->sendProcessNeighborsZ[0].index = sendNeighbors1;
    para->parH[level]->sendProcessNeighborsZ[1].index = sendNeighbors2;

    vf::gpu::findEdgeNodesCommMultiGPU(*para);

    const std::vector<std::pair<uint, uint>> expectedEdgeNodesYtoZRecv = { std::pair(0U, 1U), std::pair(0U, 2U),
                                                                         std::pair(0U, 4U), std::pair(0U, 6U),
                                                                         std::pair(0U, 8U) };
    const std::vector<std::pair<uint, uint>> expectedEdgeNodesYtoZSend = { std::pair(0U, 0U), std::pair(1U, 2U),
                                                                         std::pair(0U, 5U), std::pair(0U, 4U),
                                                                         std::pair(1U, 4U) };

    EXPECT_THAT(para->parH[level]->edgeNodesYtoZ.size(), testing::Eq(expectedEdgeNodesYtoZRecv.size()));
    EXPECT_TRUE(compareEdgeNodesRecv(para->parH[level]->edgeNodesYtoZ, expectedEdgeNodesYtoZRecv))
        << "the edgeNodesYtoZ for the receive process do not match the expected nodes";
    EXPECT_TRUE(compareEdgeNodesSend(para->parH[level]->edgeNodesYtoZ, expectedEdgeNodesYtoZSend))
        << "the edgeNodesYtoZ for the send process do not match the expected nodes";
}

//! \}
