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
//! \addtogroup gpu_DataStructureInitializer_tests DataStructureInitializer
//! \ingroup gpu_core_tests core
//! \{
//! \author Anna Wellmann
//=======================================================================================
#include <gmock/gmock.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../../testUtilitiesGPU.h"

#include "Calculation/Calculation.h"
#include "DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "basics/config/ConfigurationFile.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"

#include <parallel/NullCommunicator.h>

namespace index_rearrangement_tests
{
template <typename T>
bool vectorsAreEqual(const T *vector1, const std::vector<T>& vectorExpected)
{
    for (uint i = 0; i < vectorExpected.size(); i++) {
        if (vector1[i] != vectorExpected[i])
            return false;
    }
    return true;
}

class LevelGridBuilderDouble : public LevelGridBuilder
{
private:
    SPtr<Grid> grid;
    LevelGridBuilderDouble() = default;

    uint numberOfSendIndices;
    uint numberOfRecvIndices;

public:
    LevelGridBuilderDouble(SPtr<Grid> grid) : LevelGridBuilder(), grid(grid){};
    SPtr<Grid> getGrid(uint level) override
    {
        return grid;
    };
    void setNumberOfSendIndices(uint numberOfSendIndices)
    {
        this->numberOfSendIndices = numberOfSendIndices;
    };
    uint getNumberOfSendIndices(int direction, uint level) override
    {
        return numberOfSendIndices;
    };
    uint getNumberOfReceiveIndices(int direction, uint level) override
    {
        return numberOfRecvIndices;
    };
    void setNumberOfRecvIndices(uint numberOfRecvIndices)
    {
        this->numberOfRecvIndices = numberOfRecvIndices;
    };
};

class GridImpDouble : public GridImp
{
private:
    std::vector<uint> fluidNodeIndicesBorder;

public:
    GridImpDouble(SPtr<Object> object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                  Distribution d, uint level)
        : GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, d, level)
    {
    }

    static SPtr<GridImpDouble> makeShared(SPtr<Object> object, real startX, real startY, real startZ, real endX, real endY,
                                          real endZ, real delta, Distribution d, uint level)
    {
        SPtr<GridImpDouble> grid(std::make_shared<GridImpDouble>(object, startX, startY, startZ, endX, endY, endZ, delta, d, level));
        return grid;
    }
};

} // namespace index_rearrangement_tests

//////////////////////////////////////////////////////////////////////////
// Test reorderSendIndices
//////////////////////////////////////////////////////////////////////////
using namespace index_rearrangement_tests;

struct SendIndicesForCommAfterFtoCX {
    // data to work on
    std::vector<uint> sendIndices = { 10, 11, 12, 13, 14, 15, 16 };
    const int level = 0;
    const int direction = communication_directions::MX;
    const int numberOfProcessNeighbors = 1;
    const int indexOfProcessNeighbor = 0;

    std::vector<uint> interpolationCellCoarseToFineCoarse = { 8, 10, 12 };
    std::vector<uint> interpolationCellFineToCoarseCoarse = { 14, 16, 18 };
    const uint numNodesCtoF = (uint)interpolationCellCoarseToFineCoarse.size();
    const uint numNodesFtoC = (uint)interpolationCellFineToCoarseCoarse.size();
    uint neighborX[18] = { 0u };
    uint neighborY[18] = { 0u };
    uint neighborZ[18] = { 0u };

    // output data
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;

    // expected data
    const std::vector<uint> sendIndicesForCommAfterFtoCPositions_expected = { 4, 6, 0, 2 };
    const std::vector<uint> sendProcessNeighborX_expected = { 14, 16, 10, 12, 11, 13, 15 };
    const int numberOfSendNodesAfterFtoC_expected = (int)sendIndicesForCommAfterFtoCPositions_expected.size();
};

class IndexRearrangementForStreamsTest_reorderSendIndices : public testing::Test
{
protected:
    SendIndicesForCommAfterFtoCX sendIndices;
    SPtr<Parameter> para;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

    void act()
    {
        const uint numberOfNodes = testSubject->reorderSendIndicesForCommAfterFtoC(para->getParH(sendIndices.level)->sendProcessNeighborsX[sendIndices.indexOfProcessNeighbor], sendIndices.direction, sendIndices.level,
                                                         sendIndices.sendIndicesForCommAfterFtoCPositions);
        testSubject->setNumberOfNodes(para->getParH(sendIndices.level)->sendProcessNeighborsAfterFtoCX[sendIndices.indexOfProcessNeighbor], para->getParH(sendIndices.level)->sendProcessNeighborsAfterFtoCX[sendIndices.indexOfProcessNeighbor], numberOfNodes);
    };

private:
    void SetUp() override
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);
        builder->setNumberOfSendIndices((uint)sendIndices.sendIndices.size());

        para = testing::vf::createParameterForLevel(sendIndices.level);

        para->getParH(sendIndices.level)->fineToCoarse.numberOfCells = sendIndices.numNodesFtoC;
        para->getParH(sendIndices.level)->fineToCoarse.coarseCellIndices = &(sendIndices.interpolationCellFineToCoarseCoarse.front());
        para->getParH(sendIndices.level)->coarseToFine.coarseCellIndices = &(sendIndices.interpolationCellCoarseToFineCoarse.front());
        para->getParH(sendIndices.level)->coarseToFine.numberOfCells = sendIndices.numNodesCtoF;
        para->getParH(sendIndices.level)->neighborX = sendIndices.neighborX;
        para->getParH(sendIndices.level)->neighborY = sendIndices.neighborY;
        para->getParH(sendIndices.level)->neighborZ = sendIndices.neighborZ;

        para->setNumberOfProcessNeighborsX(sendIndices.numberOfProcessNeighbors, sendIndices.level, "send");
        para->getParH(sendIndices.level)->sendProcessNeighborsX[sendIndices.indexOfProcessNeighbor].index = sendIndices.sendIndices.data();
        para->initProcessNeighborsAfterFtoCX(sendIndices.level);

        testSubject = std::make_unique<IndexRearrangementForStreams>(
            IndexRearrangementForStreams(para, builder, communicator));
    };

    vf::parallel::NullCommunicator communicator;
};

TEST_F(IndexRearrangementForStreamsTest_reorderSendIndices, reorderSendIndicesForCommAfterFtoCX)
{
    act();

    EXPECT_THAT(sendIndices.sendIndicesForCommAfterFtoCPositions.size(),
                testing::Eq(sendIndices.sendIndicesForCommAfterFtoCPositions_expected.size()));
    EXPECT_THAT(sendIndices.sendIndicesForCommAfterFtoCPositions, testing::Eq(sendIndices.sendIndicesForCommAfterFtoCPositions_expected));

    EXPECT_THAT(para->getParH(sendIndices.level)->sendProcessNeighborsAfterFtoCX[sendIndices.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(sendIndices.numberOfSendNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(sendIndices.level)->sendProcessNeighborsX[sendIndices.indexOfProcessNeighbor].index,
                                sendIndices.sendProcessNeighborX_expected))
        << "sendProcessNeighborX[].index does not match the expected vector";
}

//////////////////////////////////////////////////////////////////////////
// Test exchangeIndicesForCommAfterFtoC
//////////////////////////////////////////////////////////////////////////

class CommunicatorDouble : public vf::parallel::NullCommunicator
{
public:
    void receiveSend(uint *buffer_receive, int, int, uint *, int, int) const override
    {
        for (int i = 0; i < (int)receivedIndices.size(); ++i) {
            *(buffer_receive + i) = receivedIndices[i];
        }
    }

    void receiveSend(real *buffer_send, int size_buffer_send, real *buffer_receive, int size_buffer_recv,
                     int neighbor_rank) const override
    {
    }

    void setReceivedIndices(const std::vector<uint>& receivedIndices)
    {
        this->receivedIndices = receivedIndices;
    }

private:
    std::vector<uint> receivedIndices;
};

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC : public testing::Test
{

public:
    void createTestSubject(vf::parallel::Communicator &Communicator)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, Communicator);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoC(sendProcess, recvProcess, sendIndicesForCommAfterFtoCPositions);
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderDouble> builder;
    std::unique_ptr<IndexRearrangementForStreams> sut;
    ProcessNeighbor27 sendProcess, recvProcess;
    const uint level = 1;
    const int indexOfProcessNeighbor = 0;
    const uint numberOfProcessNeighbors = 2;
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = { 1, 2, 3 };

private:
    void SetUp() override
    {
        para = testing::vf::createParameterForLevel(level);

        sendProcess.numberOfNodes = 3;
        recvProcess.rankNeighbor = 0;

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        builder = std::make_shared<LevelGridBuilderDouble>(grid);
    };
};

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, emptyRecvInX)
{
    CommunicatorDouble communicator;
    communicator.setReceivedIndices(std::vector<uint>());
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, zeroRecvIndexX)
{
    CommunicatorDouble communicator;
    communicator.setReceivedIndices({ 0 });
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();

    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, oneRecvIndexX)
{
    CommunicatorDouble communicator;
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesByComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, threeRecvIndicesX)
{
    CommunicatorDouble communicator;
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesByComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, sixRecvIndicesX)
{
    // this test shows the limits of the current approach. The last index is always deleted
    CommunicatorDouble communicator;
    std::vector<uint> expected = { 10, 20, 30, 40, 50 };
    std::vector<uint> receivedIndicesByComm = { 10, 20, 30, 40, 50, 60 };
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(5));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoC, recvIndicesXContainZero)
{
    CommunicatorDouble communicator;
    std::vector<uint> expected = { 0, 20, 30, 40 };
    std::vector<uint> receivedIndicesByComm(6, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(4));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

//////////////////////////////////////////////////////////////////////////
// Test reorderReceiveIndices
//////////////////////////////////////////////////////////////////////////

struct RecvIndicesForCommAfterFtoC {
    // data to work on
    std::vector<uint> recvIndices = { 10, 11, 12, 13, 14, 15, 16 };
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = {};
    const int level = 0;
    const int direction = communication_directions::MX;
    const int numberOfProcessNeighbors = 1;
    const int indexOfProcessNeighbor = 0;

    // output data
    int numberOfRecvNodesAfterFtoC;
    // and reordered recvIndices
};

class IndexRearrangementForStreamsTest_reorderRecvIndices : public testing::Test
{
protected:
    RecvIndicesForCommAfterFtoC ri;
    SPtr<Parameter> para;

    ProcessNeighbor27 neighbor;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

    void act()
    {
        ri.numberOfRecvNodesAfterFtoC = testSubject->reorderRecvIndicesForCommAfterFtoC(neighbor, 
                                                        ri.direction, ri.level,
                                                        ri.sendIndicesForCommAfterFtoCPositions);
    };

private:
    void SetUp() override
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);
        builder->setNumberOfRecvIndices((uint)ri.recvIndices.size());

        para = testing::vf::createParameterForLevel(ri.level);

        testSubject = std::make_unique<IndexRearrangementForStreams>(
            IndexRearrangementForStreams(para, builder, communicator));
        neighbor.index = ri.recvIndices.data();
    };

    vf::parallel::NullCommunicator communicator;
};

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndices, noSendIndicesForCommunicationAfterScalingFineToCoarse_receiveIndicesAreUnchanged)
{
    ri.sendIndicesForCommAfterFtoCPositions = {};
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<uint> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndices, someSendIndicesForCommunicationAfterScalingFineToCoarse_receiveIndicesAreReorderedCorrectly)
{
    ri.sendIndicesForCommAfterFtoCPositions = { 0, 2, 4, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<uint> recvIndices_expected = { 10, 12, 14, 16, 11, 13, 15 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndices, allIndicesAreSendIndicesForCommunicationAfterScalingFineToCoarse_receiveIndicesAreReorderedCorrectly)
{
    ri.sendIndicesForCommAfterFtoCPositions = { 0, 1, 2, 3, 4, 5, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<uint> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

//! \}
