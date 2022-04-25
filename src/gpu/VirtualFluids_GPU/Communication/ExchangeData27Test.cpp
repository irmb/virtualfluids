#include <gmock/gmock.h>

#include <filesystem>

#include "ExchangeData27.h"

#include <basics/config/ConfigurationFile.h>

auto RealEq = [](auto value) { 
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value); 
#else 
    return testing::FloatEq(value);
#endif
};


SPtr<Parameter> initParameterClass(std::shared_ptr<Parameter> &para)
{
    std::filesystem::path filePath = __FILE__; //  assuming that the config file is stored parallel to this file.
    filePath.replace_filename("ExchangeData27Test.cfg");
    vf::basics::ConfigurationFile config;
    config.load(filePath.string());
    return std::make_shared<Parameter>(config, 1, 0);
}

TEST(ExchangeData27Test, copyEdgeNodes_XZ_CommunicationAfterFtoC)
{
    int level = 0;
    SPtr<Parameter> para = initParameterClass(para);
    para->setMaxLevel(level + 1); // setMaxLevel resizes parH
    para->initLBMSimulationParameter(); // init parH

    // indexInSend < 5 --> in AfterFToC
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,1,0,1);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,2,0,3);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,6,0,6);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,7,0,8);

    int numNodes = 10;

    std::vector<ProcessNeighbor27> recvProcessNeighborHostAllNodes(1);
    std::vector<real> recvFs(numNodes*27, 0.1);
    recvProcessNeighborHostAllNodes[0].f[0] = recvFs.data();
    recvProcessNeighborHostAllNodes[0].numberOfNodes = numNodes;
    
    std::vector<ProcessNeighbor27> sendProcessNeighborHostAllNodes(1);
    std::vector<real> sendFs(numNodes*27, 0.0);
    sendProcessNeighborHostAllNodes[0].f[0] = sendFs.data();
    sendProcessNeighborHostAllNodes[0].numberOfNodes = numNodes;

    std::vector<ProcessNeighbor27> sendProcessNeighborHost(1);
    sendProcessNeighborHost[0].numberOfNodes = 5;

    // expected
    std::vector<real> expectedFs(numNodes, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    std::vector<real> expectedFsAllDirections;
    for (uint direction = 0; direction <= dirEND; direction ++){
        expectedFsAllDirections.insert(expectedFsAllDirections.end(), expectedFs.begin(), expectedFs.end());
    }
    
    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, recvProcessNeighborHostAllNodes, sendProcessNeighborHostAllNodes, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHostAllNodes[0].f[0],  sendProcessNeighborHostAllNodes[0].f[0] + 27*numNodes);

    EXPECT_THAT(result, testing::Eq(expectedFsAllDirections));
}

TEST(ExchangeData27Test, copyEdgeNodes_XZ_CommunicateAll)
{
    int level = 0;
    SPtr<Parameter> para = initParameterClass(para);
    para->setMaxLevel(level + 1); // setMaxLevel resizes parH
    para->initLBMSimulationParameter(); // init parH

    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,1,0,1);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,2,0,3);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,6,0,6);
    para->getParH(level)->edgeNodesXtoZ.emplace_back(0,7,0,8);

    int numNodes = 10;

    std::vector<ProcessNeighbor27> recvProcessNeighborHostAllNodes(1);
    std::vector<real> recvFs(numNodes*27, 0.1);
    recvProcessNeighborHostAllNodes[0].f[0] = recvFs.data();
    recvProcessNeighborHostAllNodes[0].numberOfNodes = numNodes;
    
    std::vector<ProcessNeighbor27> sendProcessNeighborHostAllNodes(1);
    std::vector<real> sendFs(numNodes*27, 0.0);
    sendProcessNeighborHostAllNodes[0].f[0] = sendFs.data();
    sendProcessNeighborHostAllNodes[0].numberOfNodes = numNodes;

    std::vector<ProcessNeighbor27> sendProcessNeighborHost(1);
    sendProcessNeighborHost[0].numberOfNodes = numNodes;

    // expected
    std::vector<real> expectedFs(numNodes, 0.0);
    expectedFs[1] = 0.1;
    expectedFs[3] = 0.1;
    expectedFs[6] = 0.1;
    expectedFs[8] = 0.1;
    std::vector<real> expectedFsAllDirections;
    for (uint direction = 0; direction <= dirEND; direction ++){
        expectedFsAllDirections.insert(expectedFsAllDirections.end(), expectedFs.begin(), expectedFs.end());
    }
    
    // act
    copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, recvProcessNeighborHostAllNodes, sendProcessNeighborHostAllNodes, sendProcessNeighborHost);

    // convert result to std::vector
    std::vector<real> result;
    result.assign(sendProcessNeighborHostAllNodes[0].f[0],  sendProcessNeighborHostAllNodes[0].f[0] + 27*numNodes);

    EXPECT_THAT(result, testing::Eq(expectedFsAllDirections));
}