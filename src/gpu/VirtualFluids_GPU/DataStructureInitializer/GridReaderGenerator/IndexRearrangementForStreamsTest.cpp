#include <gmock/gmock.h>

#include <algorithm>
#include <iostream>
#include <sys/types.h>
#include <vector>

#include "Logger/Logger.h"
#include "Utilities/testUtilitiesGPU.h"

#include "Communication/Communicator.h"
#include "DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "basics/config/ConfigurationFile.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"
#include "gpu/VirtualFluids_GPU/Communication/Communicator.cpp"

namespace indexRearrangementTests
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

} // namespace indexRearrangementTests

//////////////////////////////////////////////////////////////////////////
// Test reorderSendIndices
//////////////////////////////////////////////////////////////////////////
using namespace indexRearrangementTests;

struct SendIndicesForCommAfterFtoCX {
    // data to work on
    std::vector<int> sendIndices = { 10, 11, 12, 13, 14, 15, 16 };
    const int level = 0;
    const int direction = CommunicationDirections::MX;
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
    const std::vector<int> sendProcessNeighborX_expected = { 14, 16, 10, 12, 11, 13, 15 };
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
        testSubject->reorderSendIndicesForCommAfterFtoCX(sendIndices.direction, sendIndices.level, sendIndices.indexOfProcessNeighbor,
                                                         sendIndices.sendIndicesForCommAfterFtoCPositions);
    };

private:
    void SetUp() override
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);
        builder->setNumberOfSendIndices((uint)sendIndices.sendIndices.size());

        para = testingVF::createParameterForLevel(sendIndices.level);

        para->getParH(sendIndices.level)->fineToCoarse.numberOfCells = sendIndices.numNodesFtoC;
        para->getParH(sendIndices.level)->fineToCoarse.coarseCellIndices = &(sendIndices.interpolationCellFineToCoarseCoarse.front());
        para->getParH(sendIndices.level)->coarseToFine.coarseCellIndices = &(sendIndices.interpolationCellCoarseToFineCoarse.front());
        para->getParH(sendIndices.level)->coarseToFine.numberOfCells = sendIndices.numNodesCtoF;
        para->getParH(sendIndices.level)->neighborX = sendIndices.neighborX;
        para->getParH(sendIndices.level)->neighborY = sendIndices.neighborY;
        para->getParH(sendIndices.level)->neighborZ = sendIndices.neighborZ;

        para->setNumberOfProcessNeighborsX(sendIndices.numberOfProcessNeighbors, sendIndices.level, "send");
        para->getParH(sendIndices.level)->sendProcessNeighborX[sendIndices.indexOfProcessNeighbor].index = sendIndices.sendIndices.data();
        para->initProcessNeighborsAfterFtoCX(sendIndices.level);

        testSubject = std::make_unique<IndexRearrangementForStreams>(
            IndexRearrangementForStreams(para, builder, vf::gpu::Communicator::getInstance()));
    };
};

TEST_F(IndexRearrangementForStreamsTest_reorderSendIndices, reorderSendIndicesForCommAfterFtoCX)
{
    act();

    EXPECT_THAT(sendIndices.sendIndicesForCommAfterFtoCPositions.size(),
                testing::Eq(sendIndices.sendIndicesForCommAfterFtoCPositions_expected.size()));
    EXPECT_THAT(sendIndices.sendIndicesForCommAfterFtoCPositions, testing::Eq(sendIndices.sendIndicesForCommAfterFtoCPositions_expected));

    EXPECT_THAT(para->getParH(sendIndices.level)->sendProcessNeighborsAfterFtoCX[sendIndices.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(sendIndices.numberOfSendNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(sendIndices.level)->sendProcessNeighborX[sendIndices.indexOfProcessNeighbor].index,
                                sendIndices.sendProcessNeighborX_expected))
        << "sendProcessNeighborX[].index does not match the expected vector";
}

//////////////////////////////////////////////////////////////////////////
// Test exchangeIndicesForCommAfterFtoC
//////////////////////////////////////////////////////////////////////////

class CommunicationRoutineDouble : public vf::gpu::CommunicationRoutine
{
public:
    void receive_send(uint *buffer_receive, int, int, uint *, int, int) const override
    {
        for (int i = 0; i < (int)receivedIndices.size(); ++i) {
            *(buffer_receive + i) = receivedIndices[i];
        }
    }

    int getPID() const override
    {
        return 0;
    }

    void setReceivedIndices(const std::vector<uint>& receivedIndices)
    {
        this->receivedIndices = receivedIndices;
    }

private:
    std::vector<uint> receivedIndices;
};

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX : public testing::Test
{

public:
    void createTestSubject(vf::gpu::CommunicationRoutine &CommunicationRoutine)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, CommunicationRoutine);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCX(level, indexOfProcessNeighbor,
                                                     sendIndicesForCommAfterFtoCPositions);
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderDouble> builder;
    std::unique_ptr<IndexRearrangementForStreams> sut;
    const uint level = 1;
    const int indexOfProcessNeighbor = 0;
    const uint numberOfProcessNeighbors = 2;
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = { 1, 2, 3 };

private:
    void SetUp() override
    {
        para = testingVF::createParameterForLevel(level);

        para->setNumberOfProcessNeighborsX(numberOfProcessNeighbors, level, "send");
        para->initProcessNeighborsAfterFtoCX(level);
        para->getParH(level)->sendProcessNeighborsAfterFtoCX[indexOfProcessNeighbor].numberOfNodes = 3;

        para->setNumberOfProcessNeighborsX(numberOfProcessNeighbors, level, "recv");
        para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighbor].rankNeighbor = 0;

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        builder = std::make_shared<LevelGridBuilderDouble>(grid);
    };
};

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, emptyRecvInX)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices(std::vector<uint>());
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, zeroRecvIndexX)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices({ 0 });
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, oneRecvIndexX)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesByComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, threeRecvIndicesX)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesByComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, sixRecvIndicesX)
{
    // this test shows the limits of the current approach. The last index is always deleted
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30, 40, 50 };
    std::vector<uint> receivedIndicesByComm = { 10, 20, 30, 40, 50, 60 };
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(5));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, recvIndicesXContainZero)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 0, 20, 30, 40 };
    std::vector<uint> receivedIndicesByComm(6, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(4));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY : public testing::Test
{

public:
    void createTestSubject(vf::gpu::CommunicationRoutine &CommunicationRoutine)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, CommunicationRoutine);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCY(level, indexOfProcessNeighbor,
                                                     sendIndicesForCommAfterFtoCPositions);
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderDouble> builder;
    std::unique_ptr<IndexRearrangementForStreams> sut;
    const uint level = 1;
    const int indexOfProcessNeighbor = 0;
    const uint numberOfProcessNeighbors = 2;
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = { 1, 2, 3 };

private:
    void SetUp() override
    {
        para = testingVF::createParameterForLevel(level);

        para->setNumberOfProcessNeighborsY(numberOfProcessNeighbors, level, "send");
        para->initProcessNeighborsAfterFtoCY(level);
        para->getParH(level)->sendProcessNeighborsAfterFtoCY[indexOfProcessNeighbor].numberOfNodes = 3;

        para->setNumberOfProcessNeighborsY(numberOfProcessNeighbors, level, "recv");
        para->getParH(level)->recvProcessNeighborY[indexOfProcessNeighbor].rankNeighbor = 0;

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        builder = std::make_shared<LevelGridBuilderDouble>(grid);
    };
};

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, emptyRecvInY)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices(std::vector<uint>());
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, zeroRecvIndexY)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices({ 0 });
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, oneRecvIndexY)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesByComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, threeRecvIndicesY)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesByComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, sixRecvIndicesY)
{
    // this test shows the limits of the current approach. The last index is always deleted
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30, 40, 50 };
    std::vector<uint> receivedIndicesByComm = { 10, 20, 30, 40, 50, 60 };
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(5));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, recvIndicesYContainZero)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 0, 20, 30, 40 };
    std::vector<uint> receivedIndicesByComm(6, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(4));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ : public testing::Test
{

public:
    void createTestSubject(vf::gpu::CommunicationRoutine &CommunicationRoutine)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, CommunicationRoutine);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCZ(level, indexOfProcessNeighbor,
                                                     sendIndicesForCommAfterFtoCPositions);
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderDouble> builder;
    std::unique_ptr<IndexRearrangementForStreams> sut;
    const uint level = 1;
    const int indexOfProcessNeighbor = 0;
    const uint numberOfProcessNeighbors = 2;
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = { 1, 2, 3 };

private:
    void SetUp() override
    {
        para = testingVF::createParameterForLevel(level);

        para->setNumberOfProcessNeighborsZ(numberOfProcessNeighbors, level, "send");
        para->initProcessNeighborsAfterFtoCZ(level);
        para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighbor].numberOfNodes = 3;

        para->setNumberOfProcessNeighborsZ(numberOfProcessNeighbors, level, "recv");
        para->getParH(level)->recvProcessNeighborZ[indexOfProcessNeighbor].rankNeighbor = 0;

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        builder = std::make_shared<LevelGridBuilderDouble>(grid);
    };
};

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, emptyRecvInZ)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices(std::vector<uint>());
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, zeroRecvIndexZ)
{
    CommunicationRoutineDouble communicator;
    communicator.setReceivedIndices({ 0 });
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, oneRecvIndexZ)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesBZComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesBZComm.begin());
    communicator.setReceivedIndices(receivedIndicesBZComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, threeRecvIndicesZ)
{
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesBZComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesBZComm.begin());
    communicator.setReceivedIndices(receivedIndicesBZComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, sixRecvIndicesYZ)
{
    // this test shows the limits of the current approach. The last index is always deleted
    CommunicationRoutineDouble communicator;
    std::vector<uint> expected = { 10, 20, 30, 40, 50 };
    std::vector<uint> receivedIndicesByComm = { 10, 20, 30, 40, 50, 60 };
    communicator.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(communicator);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(5));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, recvIndicesZContainZero)
{
    CommunicationRoutineDouble communicator;
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
    std::vector<int> recvIndices = { 10, 11, 12, 13, 14, 15, 16 };
    std::vector<uint> sendIndicesForCommAfterFtoCPositions = {};
    const int level = 0;
    const int direction = CommunicationDirections::MX;
    const int numberOfProcessNeighbors = 1;
    const int indexOfProcessNeighbor = 0;

    // output data
    int numberOfRecvNodesAfterFtoC;
    // and reordered recvIndices
};

class IndexRearrangementForStreamsTest_reorderRecvIndicesX : public testing::Test
{
protected:
    RecvIndicesForCommAfterFtoC ri;
    SPtr<Parameter> para;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

    void setUpParaInX()
    {
        para->setNumberOfProcessNeighborsX(ri.numberOfProcessNeighbors, ri.level, "recv");
        para->initProcessNeighborsAfterFtoCX(ri.level);
        para->getParH(ri.level)->recvProcessNeighborX[ri.indexOfProcessNeighbor].index = ri.recvIndices.data();
    }

    void act()
    {
        testSubject->reorderRecvIndicesForCommAfterFtoC(ri.recvIndices.data(), ri.numberOfRecvNodesAfterFtoC,
                                                        ri.direction, ri.level,
                                                        ri.sendIndicesForCommAfterFtoCPositions);
    };

    void actWithX()
    {
        testSubject->reorderRecvIndicesForCommAfterFtoCX(ri.direction, ri.level, ri.indexOfProcessNeighbor,
                                                         ri.sendIndicesForCommAfterFtoCPositions);
    };

private:
    void SetUp() override
    {
        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::WARNING);

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);
        builder->setNumberOfRecvIndices((uint)ri.recvIndices.size());

        para = testingVF::createParameterForLevel(ri.level);

        testSubject = std::make_unique<IndexRearrangementForStreams>(
            IndexRearrangementForStreams(para, builder, vf::gpu::Communicator::getInstance()));
    };
};

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, noRecvIndices)
{
    ri.sendIndicesForCommAfterFtoCPositions = {};
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, someRecvIndices)
{
    ri.sendIndicesForCommAfterFtoCPositions = { 0, 2, 4, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 12, 14, 16, 11, 13, 15 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, allRecvIndices)
{
    ri.sendIndicesForCommAfterFtoCPositions = { 0, 1, 2, 3, 4, 5, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    act();
    EXPECT_THAT(ri.numberOfRecvNodesAfterFtoC, testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_THAT(ri.recvIndices, testing::Eq(recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, noRecvIndices_callWithX)
{
    setUpParaInX();

    ri.sendIndicesForCommAfterFtoCPositions = {};
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    this->actWithX();

    EXPECT_THAT(para->getParH(ri.level)->recvProcessNeighborsAfterFtoCX[ri.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(ri.level)->recvProcessNeighborX[ri.indexOfProcessNeighbor].index,
                                recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, someRecvIndices_callWithX)
{
    setUpParaInX();

    ri.sendIndicesForCommAfterFtoCPositions = { 0, 2, 4, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 12, 14, 16, 11, 13, 15 };

    actWithX();

    EXPECT_THAT(para->getParH(ri.level)->recvProcessNeighborsAfterFtoCX[ri.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(ri.level)->recvProcessNeighborX[ri.indexOfProcessNeighbor].index,
                                recvIndices_expected));
}

TEST_F(IndexRearrangementForStreamsTest_reorderRecvIndicesX, allRecvIndices_callWithX)
{
    setUpParaInX();

    ri.sendIndicesForCommAfterFtoCPositions = { 0, 1, 2, 3, 4, 5, 6 };
    auto numberOfRecvNodesAfterFtoC_expected = ri.sendIndicesForCommAfterFtoCPositions.size();
    std::vector<int> recvIndices_expected = { 10, 11, 12, 13, 14, 15, 16 };

    actWithX();

    EXPECT_THAT(para->getParH(ri.level)->recvProcessNeighborsAfterFtoCX[ri.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(numberOfRecvNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(ri.level)->recvProcessNeighborX[ri.indexOfProcessNeighbor].index,
                                recvIndices_expected));
}
