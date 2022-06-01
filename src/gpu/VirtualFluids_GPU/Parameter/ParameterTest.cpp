#include <gmock/gmock.h>

#include <filesystem>
#include <iostream>
#include <string>

#include "Parameter.h"
#include <basics/config/ConfigurationFile.h>

auto RealEq = [](auto value) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value);
#else
    return testing::FloatEq(value);
#endif
};

TEST(ParameterTest, passingEmptyFileWithoutPath_ShouldThrow)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("parameterTest_emptyfile.cfg");

    vf::basics::ConfigurationFile config;
    config.load(filePath.string());

    EXPECT_THROW(Parameter para(config, 1, 0), std::runtime_error);
}

// TODO: test setPossNeighborFilesX
// TODO: test default values

TEST(ParameterTest, check_all_Parameter_CanBePassedToConstructor)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("parameterTest.cfg");

    vf::basics::ConfigurationFile config;
    config.load(filePath.string());

    Parameter para(config, 1, 0);

    // this two parameters need to be defined in each config file
    EXPECT_THAT(para.getOutputPath(), testing::Eq("/output/path"));
    EXPECT_THAT(para.getgeoVec(), testing::Eq("/path/to/grid/geoVec.dat"));
    // ... all grid files could be tested as well

    // test optional parameter
    EXPECT_THAT(para.getMaxDev(), testing::Eq(2));
    EXPECT_THAT(para.getDevices(), testing::ElementsAreArray({ 2, 3 }));
    EXPECT_THAT(para.getOutputPrefix(), testing::Eq("MyPrefix"));
    EXPECT_THAT(para.getPrintFiles(), testing::Eq(true));
    EXPECT_THAT(para.getIsGeometryValues(), testing::Eq(true));
    EXPECT_THAT(para.getCalc2ndOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalc3rdOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcHighOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcMedian(), testing::Eq(true));
    EXPECT_THAT(para.getCalcCp(), testing::Eq(true));
    EXPECT_THAT(para.getCalcDragLift(), testing::Eq(true));
    EXPECT_THAT(para.getWriteVeloASCIIfiles(), testing::Eq(true));
    EXPECT_THAT(para.getCalcPlaneConc(), testing::Eq(true));
    EXPECT_THAT(para.getConcFile(), testing::Eq(true));
    EXPECT_THAT(para.isStreetVelocityFile(), testing::Eq(true));
    EXPECT_THAT(para.getUseMeasurePoints(), testing::Eq(true));
    EXPECT_THAT(para.getUseWale(), testing::Eq(true));
    EXPECT_THAT(para.getUseInitNeq(), testing::Eq(true));
    EXPECT_THAT(para.getSimulatePorousMedia(), testing::Eq(true));

    EXPECT_THAT(para.getD3Qxx(), testing::Eq(99));
    EXPECT_THAT(para.getTEnd(), testing::Eq(33));
    EXPECT_THAT(para.getTOut(), testing::Eq(22));
    EXPECT_THAT(para.getTStartOut(), testing::Eq(11));
    EXPECT_THAT(para.getTimeCalcMedStart(), testing::Eq(22));
    EXPECT_THAT(para.getTimeCalcMedEnd(), testing::Eq(44));
    EXPECT_THAT(para.getPressInID(), testing::Eq(25));
    EXPECT_THAT(para.getPressOutID(), testing::Eq(26));
    EXPECT_THAT(para.getPressInZ(), testing::Eq(27));
    EXPECT_THAT(para.getPressOutZ(), testing::Eq(28));

    EXPECT_THAT(para.getDiffOn(), testing::Eq(true));
    EXPECT_THAT(para.getDiffMod(), testing::Eq(99));
    EXPECT_THAT(para.getDiffusivity(), RealEq(1.11));
    EXPECT_THAT(para.getTemperatureInit(), RealEq(2.22));
    EXPECT_THAT(para.getTemperatureBC(), RealEq(3.33));

    EXPECT_THAT(para.getViscosity(), RealEq(4.44));
    EXPECT_THAT(para.getVelocity(), RealEq(5.55));
    EXPECT_THAT(para.getViscosityRatio(), RealEq(6.66));
    EXPECT_THAT(para.getVelocityRatio(), RealEq(7.77));
    EXPECT_THAT(para.getDensityRatio(), RealEq(8.88));
    EXPECT_THAT(para.getPressRatio(), RealEq(9.99));

    EXPECT_THAT(para.getRealX(), RealEq(0.1));
    EXPECT_THAT(para.getRealY(), RealEq(0.2));
    EXPECT_THAT(para.getFactorPressBC(), RealEq(0.3));

    EXPECT_THAT(para.getReadGeo(), testing::Eq(true));
    EXPECT_THAT(para.getGeometryFileC(), testing::Eq("/pass/to/c"));
    EXPECT_THAT(para.getGeometryFileM(), testing::Eq("/pass/to/m"));
    EXPECT_THAT(para.getGeometryFileF(), testing::Eq("/pass/to/f"));

    EXPECT_THAT(para.getclockCycleForMP(), RealEq(0.4));
    EXPECT_THAT(para.getTimestepForMP(), testing::Eq(4));

    std::vector<real> forces{ 2.0, 2.1, 2.2 };
    double *forces_actual = para.getForcesDouble();
    for (size_t i = 0; i < forces.size(); ++i) {
        EXPECT_THAT((real)forces_actual[i], RealEq(forces[i]));
    }

    std::vector<real> limiters{ 3.0, 3.1, 3.2 };
    double *limiters_actual = para.getQuadricLimitersDouble();
    for (size_t i = 0; i < limiters.size(); ++i) {
        EXPECT_THAT((real)limiters_actual[i], RealEq(limiters[i]));
    }

    EXPECT_THAT(para.getCalcParticle(), testing::Eq(true));
    EXPECT_THAT(para.getParticleBasicLevel(), testing::Eq(1));
    EXPECT_THAT(para.getParticleInitLevel(), testing::Eq(2));
    EXPECT_THAT(para.getNumberOfParticles(), testing::Eq(1111));
    EXPECT_THAT(para.getStartXHotWall(), RealEq(4.1));
    EXPECT_THAT(para.getEndXHotWall(), RealEq(4.2));

    EXPECT_THAT(para.getTimeDoCheckPoint(), testing::Eq(33));
    EXPECT_THAT(para.getTimeDoRestart(), testing::Eq(44));
    EXPECT_THAT(para.getDoCheckPoint(), testing::Eq(true));
    EXPECT_THAT(para.getDoRestart(), testing::Eq(true));
    EXPECT_THAT(para.getMaxLevel(), testing::Eq(1)); // NOGL - 1

    EXPECT_THAT(para.getGridX(), testing::ElementsAreArray({ 100, 101 }));
    EXPECT_THAT(para.getGridY(), testing::ElementsAreArray({ 200, 201 }));
    EXPECT_THAT(para.getGridZ(), testing::ElementsAreArray({ 300, 301 }));
    EXPECT_THAT(para.getDistX(), testing::ElementsAreArray({ 400, 401 }));
    EXPECT_THAT(para.getDistY(), testing::ElementsAreArray({ 500, 501 }));
    EXPECT_THAT(para.getDistZ(), testing::ElementsAreArray({ 600, 601 }));

    EXPECT_THAT(para.getMainKernel(), testing::Eq("KernelName"));
    EXPECT_THAT(para.getMultiKernelOn(), testing::Eq(true));
    EXPECT_THAT(para.getMultiKernelLevel(), testing::ElementsAreArray({ 3, 2, 1 }));

    std::vector<std::string> kernel{ "Kernel1", "Kernel2", "Kernel3" };
    auto kernel_actual = para.getMultiKernel();
    for (size_t i = 0; i < kernel.size(); ++i) {
        EXPECT_THAT(kernel_actual[i], testing::Eq(kernel[i]));
    }

    EXPECT_THAT(para.getCoarse(), testing::Eq(0));
    EXPECT_THAT(para.getFine(), testing::Eq(1)); // NOGL - 1
    EXPECT_THAT(para.parH.size(), testing::Eq(2));
    EXPECT_THAT(para.parD.size(), testing::Eq(2));
}

static std::shared_ptr<Parameter> initParameterClass()
{
    std::filesystem::path filePath = __FILE__; //  assuming that the config file is stored parallel to this file.
    filePath.replace_filename("parameterTest.cfg");
    vf::basics::ConfigurationFile config;
    config.load(filePath.string());
    return std::make_shared<Parameter>(config, 1, 0);
}

TEST(ParameterTest, findEdgeNodesXY_shouldReturnCorrectVector)
{

    auto para = initParameterClass();
    para->initLBMSimulationParameter();

    int level = 0;
    para->parH[level]->recvProcessNeighborX.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborY.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborY.push_back(ProcessNeighbor27());

    int numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborX.size() - 1;
    int numSendNeighbor = (int)para->parH[level]->sendProcessNeighborY.size() - 1;

    const int sizeRecv                                                    = 6;
    const int sizeSend                                                    = 10;
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborY[numSendNeighbor].numberOfNodes = sizeSend;

    int recvNeighbors[sizeRecv]                                   = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].index = recvNeighbors;

    int sendNeighbors[sizeSend]                                   = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->sendProcessNeighborY[numSendNeighbor].index = sendNeighbors;

    para->findEdgeNodesCommMultiGPU();

    std::vector<std::pair<int, int>> expectedEdgeNodesXtoYRecv = { std::pair(numRecvNeighbor, 0),
                                                                   std::pair(numRecvNeighbor, 4),
                                                                   std::pair(numRecvNeighbor, 5) };

    std::vector<std::pair<int, int>> expectedEdgeNodesXtoYSend = { std::pair(numSendNeighbor, 1),
                                                                   std::pair(numSendNeighbor, 6),
                                                                   std::pair(numSendNeighbor, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoY.size(), testing::Eq(expectedEdgeNodesXtoYRecv.size()));

    bool vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoYRecv.size(); i++) {
        if (para->parH[level]->edgeNodesXtoY[i].indexOfProcessNeighborRecv != expectedEdgeNodesXtoYRecv[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesXtoY[i].indexInRecvBuffer != expectedEdgeNodesXtoYRecv[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);

    vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoYSend.size(); i++) {
        if (para->parH[level]->edgeNodesXtoY[i].indexOfProcessNeighborSend != expectedEdgeNodesXtoYSend[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesXtoY[i].indexInSendBuffer != expectedEdgeNodesXtoYSend[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);
}

TEST(ParameterTest, findEdgeNodesXZ_shouldReturnCorrectVector)
{

    auto para = initParameterClass();
    para->initLBMSimulationParameter();

    int level = 0;
    para->parH[level]->recvProcessNeighborX.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());

    int numRecvNeighbor = (int)para->parH[level]->recvProcessNeighborX.size() - 1;
    int numSendNeighbor = (int)para->parH[level]->sendProcessNeighborZ.size() - 1;

    const int sizeRecv = 10;
    const int sizeSend = 6;

    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborZ[numSendNeighbor].numberOfNodes = sizeSend;

    int recvNeighbors[sizeRecv]                                   = { 20, 1, 21, 22, 6, 23, 5, 24, 25, 26 };
    para->parH[level]->recvProcessNeighborX[numRecvNeighbor].index = recvNeighbors;

    int sendNeighbors[sizeSend]                                   = { 1, 2, 3, 4, 5, 6 };
    para->parH[level]->sendProcessNeighborZ[numSendNeighbor].index = sendNeighbors;

    para->findEdgeNodesCommMultiGPU();

    std::vector<std::pair<int, int>> expectedEdgeNodesXtoZRecv = { std::pair(numRecvNeighbor, 1),
                                                                   std::pair(numRecvNeighbor, 4),
                                                                   std::pair(numRecvNeighbor, 6) };
    std::vector<std::pair<int, int>> expectedEdgeNodesXtoZSend = { std::pair(numSendNeighbor, 0),
                                                                   std::pair(numSendNeighbor, 5),
                                                                   std::pair(numSendNeighbor, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesXtoZ.size(), testing::Eq(expectedEdgeNodesXtoZRecv.size()));

    bool vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoZRecv.size(); i++) {
        if (para->parH[level]->edgeNodesXtoZ[i].indexOfProcessNeighborRecv != expectedEdgeNodesXtoZRecv[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesXtoZ[i].indexInRecvBuffer != expectedEdgeNodesXtoZRecv[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);

    vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoZRecv.size(); i++) {
        if (para->parH[level]->edgeNodesXtoZ[i].indexOfProcessNeighborSend != expectedEdgeNodesXtoZSend[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesXtoZ[i].indexInSendBuffer != expectedEdgeNodesXtoZSend[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);
}

TEST(ParameterTest, findEdgeNodesYZ_shouldReturnCorrectVector)
{

    auto para = initParameterClass();
    para->initLBMSimulationParameter();

    int level = 0;

    para->parH[level]->recvProcessNeighborY.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());
    para->parH[level]->sendProcessNeighborZ.push_back(ProcessNeighbor27());

    const int sizeRecv  = 10;
    const int sizeSend1 = 6;
    const int sizeSend2 = 5;

    para->parH[level]->recvProcessNeighborY[0].numberOfNodes = sizeRecv;
    para->parH[level]->sendProcessNeighborZ[0].numberOfNodes = sizeSend1;
    para->parH[level]->sendProcessNeighborZ[1].numberOfNodes = sizeSend2;

    int recvNeighbors[sizeRecv]                     = { 20, 1, 9, 22, 6, 23, 5, 24, 11, 26 };
    para->parH[level]->recvProcessNeighborY[0].index = recvNeighbors;

    int sendNeighbors1[sizeSend1]                   = { 1, 2, 3, 4, 5, 6 };
    int sendNeighbors2[sizeSend2]                   = { 7, 8, 9, 10, 11 };
    para->parH[level]->sendProcessNeighborZ[0].index = sendNeighbors1;
    para->parH[level]->sendProcessNeighborZ[1].index = sendNeighbors2;

    para->findEdgeNodesCommMultiGPU();

    std::vector<std::pair<int, int>> expectedEdgeNodesXtoZRecv = { std::pair(0, 1), std::pair(0, 2), std::pair(0, 4),
                                                                   std::pair(0, 6), std::pair(0, 8) };
    std::vector<std::pair<int, int>> expectedEdgeNodesXtoZSend = { std::pair(0, 0), std::pair(1, 2), std::pair(0, 5),
                                                                   std::pair(0, 4), std::pair(1, 4) };

    EXPECT_THAT(para->parH[level]->edgeNodesYtoZ.size(), testing::Eq(expectedEdgeNodesXtoZRecv.size()));

    bool vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoZRecv.size(); i++) {
        if (para->parH[level]->edgeNodesYtoZ[i].indexOfProcessNeighborRecv != expectedEdgeNodesXtoZRecv[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesYtoZ[i].indexInRecvBuffer != expectedEdgeNodesXtoZRecv[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);

    vectorsAreIdentical = true;
    for (int i = 0; i < (int)expectedEdgeNodesXtoZRecv.size(); i++) {
        if (para->parH[level]->edgeNodesYtoZ[i].indexOfProcessNeighborSend != expectedEdgeNodesXtoZSend[i].first) {
            vectorsAreIdentical = false;
            break;
        }
        if (para->parH[level]->edgeNodesYtoZ[i].indexInSendBuffer != expectedEdgeNodesXtoZSend[i].second) {
            vectorsAreIdentical = false;
            break;
        }
    }

    EXPECT_TRUE(vectorsAreIdentical);
}