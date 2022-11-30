#include <gmock/gmock.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <vector>

#include "../utilities/testUtilitiesGPU.h"

#include "Communication/Communicator.h"
#include "DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "basics/config/ConfigurationFile.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"
#include "gpu/VirtualFluids_GPU/Communication/Communicator.cpp"

template <typename T>
bool vectorsAreEqual(const T * vector1, const std::vector<T> vectorExpected)
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

public:
    LevelGridBuilderDouble(SPtr<Grid> grid) : LevelGridBuilder(), grid(grid){};
    SPtr<Grid> getGrid(uint level) override
    {
        return grid;
    };
    std::shared_ptr<Grid> getGrid(int level, int box) override
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
};

class GridImpDouble : public GridImp
{
private:
    std::vector<uint> fluidNodeIndicesBorder;

public:
    GridImpDouble(Object *object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                  Distribution d, uint level)
        : GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, d, level)
    {
    }

    static SPtr<GridImpDouble> makeShared(Object *object, real startX, real startY, real startZ, real endX, real endY,
                                          real endZ, real delta, Distribution d, uint level)
    {
        SPtr<GridImpDouble> grid(new GridImpDouble(object, startX, startY, startZ, endX, endY, endZ, delta, d, level));
        return grid;
    }

    void setFluidNodeIndicesBorder(std::vector<uint> fluidNodeIndicesBorder)
    {
        this->fluidNodeIndicesBorder = fluidNodeIndicesBorder;
    }

    bool isSparseIndexInFluidNodeIndicesBorder(uint &sparseIndex) const override
    {
        return std::find(this->fluidNodeIndicesBorder.begin(), this->fluidNodeIndicesBorder.end(), sparseIndex) !=
               this->fluidNodeIndicesBorder.end();
    }
};

//////////////////////////////////////////////////////////////////////////
// Test reorderSendIndices
//////////////////////////////////////////////////////////////////////////

struct SendIndicesForCommAfterFtoCX {
    // data to work on
    std::vector<int> sendIndices = { 10, 11, 12, 13, 14, 15, 16 };
    const int level = 0;
    const int direction = CommunicationDirections::MX;
    const int numberOfProcessNeighbors = 1;
    const int indexOfProcessNeighbor = 0;

    std::vector<uint> iCellCFC = { 8, 10, 12 };
    std::vector<uint> iCellFCC = { 14, 16, 18 };
    const uint kCF = (uint)iCellCFC.size();
    const uint kFC = (uint)iCellFCC.size();
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
    SendIndicesForCommAfterFtoCX si;
    SPtr<Parameter> para;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

    void act()
    {
        testSubject->reorderSendIndicesForCommAfterFtoCX(si.direction, si.level, si.indexOfProcessNeighbor,
                                                         si.sendIndicesForCommAfterFtoCPositions);
    };

private:
    std::unique_ptr<IndexRearrangementForStreams> createTestSubjectForReorderSendIndices()
    {
        logging::Logger::addStream(&std::cout);

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);
        builder->setNumberOfSendIndices((uint)si.sendIndices.size());

        para = testingVF::createParameterForLevel(si.level);


        para->getParH(si.level)->intFC.kFC = si.kFC;
        para->getParH(si.level)->intFC.ICellFCC = &(si.iCellFCC.front());
        para->getParH(si.level)->intCF.ICellCFC = &(si.iCellCFC.front());
        para->getParH(si.level)->intCF.kCF = si.kCF;
        para->getParH(si.level)->neighborX = si.neighborX;
        para->getParH(si.level)->neighborY = si.neighborY;
        para->getParH(si.level)->neighborZ = si.neighborZ;

        para->setNumberOfProcessNeighborsX(si.numberOfProcessNeighbors, si.level, "send");
        para->getParH(si.level)->sendProcessNeighborX[si.indexOfProcessNeighbor].index = si.sendIndices.data();
        para->initProcessNeighborsAfterFtoCX(si.level);

        return std::make_unique<IndexRearrangementForStreams>(
            IndexRearrangementForStreams(para, builder, vf::gpu::Communicator::getInstance()));
    };

    void SetUp() override
    {
        para = std::make_shared<Parameter>();
        testSubject = createTestSubjectForReorderSendIndices();
    };
};

TEST_F(IndexRearrangementForStreamsTest_reorderSendIndices, reorderSendIndicesForCommAfterFtoCX)
{
    act();

    EXPECT_THAT(si.sendIndicesForCommAfterFtoCPositions.size(),
                testing::Eq(si.sendIndicesForCommAfterFtoCPositions_expected.size()));
    EXPECT_THAT(si.sendIndicesForCommAfterFtoCPositions, testing::Eq(si.sendIndicesForCommAfterFtoCPositions_expected));

    EXPECT_THAT(para->getParH(si.level)->sendProcessNeighborsAfterFtoCX[si.indexOfProcessNeighbor].numberOfNodes,
                testing::Eq(si.numberOfSendNodesAfterFtoC_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(si.level)->sendProcessNeighborX[si.indexOfProcessNeighbor].index,
                                si.sendProcessNeighborX_expected))
        << "sendProcessNeighborX[].index does not match the expected vector";
}

//////////////////////////////////////////////////////////////////////////
// Test exchangeIndicesForCommAfterFtoC
//////////////////////////////////////////////////////////////////////////

class CommunicatorDouble : public vf::gpu::Communicator
{
public:
    static CommunicatorDouble &getInstance()
    {
        static CommunicatorDouble comm;
        return comm;
    }

    void exchangeIndices(uint *rbuf, int count_r, int nb_rank_r, uint *sbuf, int count_s, int nb_rank_s) override
    {
        for (int i = 0; i < (int)receivedIndices.size(); ++i) {
            *(rbuf + i) = receivedIndices[i];
        }
    }

    void setReceivedIndices(std::vector<uint> receivedIndices)
    {
        this->receivedIndices = receivedIndices;
    }

private:
    std::vector<uint> receivedIndices;
};

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX : public testing::Test
{

public:
    void createTestSubject(vf::gpu::Communicator &communicator)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, communicator);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCX(level, indexOfProcessNeighbor, CommunicationDirections::PX,
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
        logging::Logger::addStream(&std::cout);

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
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices(std::vector<uint>());
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, zeroRecvIndexX)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices({0});
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, oneRecvIndexX)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesByComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, threeRecvIndicesX)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesByComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCX, sixRecvIndicesX)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30, 40, 50, 60 };
    std::vector<uint> receivedIndicesByComm(expected.size(), 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(6));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY : public testing::Test
{

public:
    void createTestSubject(vf::gpu::Communicator &communicator)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, communicator);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCY(level, indexOfProcessNeighbor, CommunicationDirections::PY,
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
        logging::Logger::addStream(&std::cout);

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
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices(std::vector<uint>());
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, zeroRecvIndexY)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices({0});
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, oneRecvIndexY)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesByComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, threeRecvIndicesY)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesByComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCY, sixRecvIndicesY)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30, 40, 50, 60 };
    std::vector<uint> receivedIndicesByComm(expected.size(), 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesByComm.begin());
    comm.setReceivedIndices(receivedIndicesByComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(6));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

class IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ : public testing::Test
{

public:
    void createTestSubject(vf::gpu::Communicator &communicator)
    {
        sut = std::make_unique<IndexRearrangementForStreams>(para, builder, communicator);
    }

protected:
    std::vector<uint> act()
    {
        return sut->exchangeIndicesForCommAfterFtoCZ(level, indexOfProcessNeighbor, CommunicationDirections::PZ,
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
        logging::Logger::addStream(&std::cout);

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
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices(std::vector<uint>());
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, zeroRecvIndexZ)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    comm.setReceivedIndices({0});
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(0));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, oneRecvIndexZ)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10 };
    std::vector<uint> receivedIndicesBZComm(4, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesBZComm.begin());
    comm.setReceivedIndices(receivedIndicesBZComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(1));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, threeRecvIndicesZ)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30 };
    std::vector<uint> receivedIndicesBZComm(5, 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesBZComm.begin());
    comm.setReceivedIndices(receivedIndicesBZComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(3));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}

TEST_F(IndexRearrangementForStreamsTest_exchangeIndicesForCommAfterFtoCZ, sixRecvIndicesZ)
{
    CommunicatorDouble &comm = CommunicatorDouble::getInstance();
    std::vector<uint> expected = { 10, 20, 30, 40, 50, 60 };
    std::vector<uint> receivedIndicesBZComm(expected.size(), 0);
    std::copy(expected.begin(), expected.end(), receivedIndicesBZComm.begin());
    comm.setReceivedIndices(receivedIndicesBZComm);
    createTestSubject(comm);

    std::vector<uint> recvIndicesForCommAfterFtoCPositions = act();
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions.size(), testing::Eq(6));
    EXPECT_THAT(recvIndicesForCommAfterFtoCPositions, testing::Eq(expected));
}
