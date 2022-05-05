#include <gmock/gmock.h>

#include <filesystem>

#include "ExchangeData27.h"

#include <basics/config/ConfigurationFile.h>

SPtr<Parameter> initParameterClass()
{
    std::filesystem::path filePath = __FILE__; //  assuming that the config file is stored parallel to this file.
    filePath.replace_filename("ExchangeData27Test.cfg");
    vf::basics::ConfigurationFile config;
    config.load(filePath.string());
    return std::make_shared<Parameter>(config, 1, 0);
}

void setUpFsByCopyingF0(std::vector<real> &distributionVector, int numberOfNodes)
{
    for (uint direction = 0; direction < dirEND; direction++) {
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
    
    void SetUp() override
    {
        para = initParameterClass();
        para->setMaxLevel(level + 1);       // setMaxLevel resizes parH
        para->initLBMSimulationParameter(); // init parH

        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 1, 0, 1);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 6, 0, 6);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 2, 0, 3);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
        para->getParH(level)->edgeNodesXtoZ.emplace_back(0, 7, 0, 8);
    }

    SPtr<std::vector<ProcessNeighbor27>> setUpRecvProcessNeighbors(int numberOfNodesInRecv){
        SPtr<std::vector<ProcessNeighbor27>>recvProcessNeighborHost=std::make_shared<std::vector<ProcessNeighbor27>>(1);
        recvFs.resize(numberOfNodesInRecv);
        std::fill(recvFs.begin(), recvFs.end(), 0.5); // 0.5s should not be copied
        for (LBMSimulationParameter::EdgeNodePositions edgeNode : para->getParH(level)->edgeNodesXtoZ) {
            recvFs[edgeNode.indexInRecvBuffer] = 0.1; // 0.1s should be copied
        }
        setUpFsByCopyingF0(recvFs, numberOfNodesInRecv);
        (*recvProcessNeighborHost)[0].f[0]          = recvFs.data();
        (*recvProcessNeighborHost)[0].numberOfNodes = numberOfNodesInRecv;
        return recvProcessNeighborHost;
    }   

};

TEST_F(ExchangeData27Test_CopyEdgeNodesXZTest, copyEdgeNodes_XZ_CommunicationAfterFtoC_recvVectorFullSize)
{
    int numNodesAfterFtoC = 5; // indexInSend < 5 --> mode is in AfterFToC

    // recvProcessNeighborHost
    SPtr<std::vector<ProcessNeighbor27>>recvProcessNeighborHost = setUpRecvProcessNeighbors(numNodes);

    // sendProcessNeighborHost
    std::vector<ProcessNeighbor27> sendProcessNeighborHost(1);
    std::vector<real> sendFs(27 * numNodesAfterFtoC, 0.0);
    sendProcessNeighborHost[0].f[0]          = sendFs.data();
    sendProcessNeighborHost[0].numberOfNodes = numNodesAfterFtoC;

    // expected
    std::vector<real> expectedFs(numNodesAfterFtoC, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    setUpFsByCopyingF0(expectedFs, numNodesAfterFtoC);

    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, *recvProcessNeighborHost, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHost[0].f[0], sendProcessNeighborHost[0].f[0] + 27 * numNodesAfterFtoC);

    EXPECT_THAT(result, testing::Eq(expectedFs));
}

TEST_F(ExchangeData27Test_CopyEdgeNodesXZTest, copyEdgeNodes_XZ_CommunicateAll)
{
    // recvProcessNeighborHost
    SPtr<std::vector<ProcessNeighbor27>>recvProcessNeighborHost = setUpRecvProcessNeighbors(numNodes);

    // sendProcessNeighborHost
    std::vector<ProcessNeighbor27> sendProcessNeighborHost(1);
    std::vector<real> sendFs(27 * numNodes, 0.0);
    sendProcessNeighborHost[0].f[0]          = sendFs.data();
    sendProcessNeighborHost[0].numberOfNodes = numNodes;

    // expected
    std::vector<real> expectedFs(numNodes, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    expectedFs[6] = 0.1;
    expectedFs[8] = 0.1;
    setUpFsByCopyingF0(expectedFs, numNodes);

    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, *recvProcessNeighborHost, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHost[0].f[0], sendProcessNeighborHost[0].f[0] + 27 * numNodes);

    EXPECT_THAT(result, testing::Eq(expectedFs));
}