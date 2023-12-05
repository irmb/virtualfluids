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
//! \author Martin Schoenherr
//======================================================================================

#include <gmock/gmock.h>

#include <filesystem>

#include "ExchangeData27.h"
#include "gpu/core/Calculation/Calculation.h"

#include <basics/config/ConfigurationFile.h>

void setUpFsByCopyingF0(std::vector<real> &distributionVector, int numberOfNodes)
{
    for (uint direction = 0; direction < vf::lbm::dir::ENDDIR; direction++) {
        distributionVector.insert(distributionVector.end(), distributionVector.begin(),
                                  distributionVector.begin() + numberOfNodes);
    }
}

class ExchangeData27Test_CopyEdgeNodesXZTest : public testing::Test
{
protected:
    SPtr<Parameter> para;
    int level    = 0;
    int numNodes = 10;
    std::vector<real> recvFs;
    std::vector<real> sendFs;
    std::vector<ProcessNeighbor27> sendProcessNeighborHost;
    std::vector<ProcessNeighbor27> recvProcessNeighborHost;

    void SetUp() override
    {
        para = std::make_shared<Parameter>();
        para->setMaxLevel(level + 1);       // setMaxLevel resizes parH
        para->initLBMSimulationParameter(); // init parH

        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 1, 0, 1);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 6, 0, 6);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 2, 0, 3);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
    }

    void setUpRecvProcessNeighbors(int numberOfNodesInRecv)
    {
        recvFs.resize(numberOfNodesInRecv);
        std::fill(recvFs.begin(), recvFs.end(), 0.5); // 0.5s should not be copied
        for (LBMSimulationParameter::EdgeNodePositions edgeNode : para->getParH(level)->edgeNodesXtoZ) {
            if (edgeNode.indexInRecvBuffer > numberOfNodesInRecv) {
                continue;
            }
            recvFs[edgeNode.indexInRecvBuffer] = 0.1; // 0.1s should be copied
        }
        setUpFsByCopyingF0(recvFs, numberOfNodesInRecv);

        recvProcessNeighborHost.resize(1);
        recvProcessNeighborHost[0].f[0]          = recvFs.data();
        recvProcessNeighborHost[0].numberOfNodes = numberOfNodesInRecv;
    }

    void setUpSendProcessNeighbors(int numberOfNodesInSend)
    {
        sendFs.resize(27 * numberOfNodesInSend);
        std::fill(sendFs.begin(), sendFs.end(), 0.0);

        sendProcessNeighborHost.resize(1);
        sendProcessNeighborHost[0].f[0]          = sendFs.data();
        sendProcessNeighborHost[0].numberOfNodes = numberOfNodesInSend;
    }
};

TEST_F(ExchangeData27Test_CopyEdgeNodesXZTest, copyEdgeNodes_XZ_CommunicationAfterFtoC_recvVectorFullSize)
{
    int numNodesAfterFtoC = 5; // indexInSend < 5 --> mode is in AfterFToC
    setUpRecvProcessNeighbors(numNodes);
    setUpSendProcessNeighbors(numNodesAfterFtoC);

    // expected
    std::vector<real> expectedFs(numNodesAfterFtoC, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    setUpFsByCopyingF0(expectedFs, numNodesAfterFtoC);

    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, recvProcessNeighborHost, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHost[0].f[0], sendProcessNeighborHost[0].f[0] + 27 * numNodesAfterFtoC);

    EXPECT_THAT(result, testing::Eq(expectedFs));
}

TEST_F(ExchangeData27Test_CopyEdgeNodesXZTest, copyEdgeNodes_XZ_CommunicationAfterFtoC_recvVectorShort)
{
    int numNodesAfterFtoC = 5; // indexInSend < 5 --> mode is in AfterFToC
    setUpRecvProcessNeighbors(numNodesAfterFtoC);
    setUpSendProcessNeighbors(numNodesAfterFtoC);

    // expected
    std::vector<real> expectedFs(numNodesAfterFtoC, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    setUpFsByCopyingF0(expectedFs, numNodesAfterFtoC);

    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, recvProcessNeighborHost, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHost[0].f[0], sendProcessNeighborHost[0].f[0] + 27 * numNodesAfterFtoC);

    EXPECT_THAT(result, testing::Eq(expectedFs));
}

TEST_F(ExchangeData27Test_CopyEdgeNodesXZTest, copyEdgeNodes_XZ_CommunicateAll)
{
    setUpRecvProcessNeighbors(numNodes);
    setUpSendProcessNeighbors(numNodes);

    // expected
    std::vector<real> expectedFs(numNodes, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    expectedFs[6] = 0.1;
    expectedFs[8] = 0.1;
    setUpFsByCopyingF0(expectedFs, numNodes);

    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, recvProcessNeighborHost, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHost[0].f[0], sendProcessNeighborHost[0].f[0] + 27 * numNodes);

    EXPECT_THAT(result, testing::Eq(expectedFs));
}