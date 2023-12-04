#include "GridGenerator.h"
#include <gmock/gmock.h>

#include <basics/DataTypes.h>
#include "GPU/CudaMemoryManager.h"
#include "IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"

#include <parallel/NullCommunicator.h>

namespace GridGeneratorTest
{

class LevelGridBuilderStub : public LevelGridBuilder
{
private:
    const SPtr<Grid> grid;
    LevelGridBuilderStub() = default;

public:
    uint numberOfSendIndices = 0;

    explicit LevelGridBuilderStub(SPtr<Grid> grid) : LevelGridBuilder(), grid(grid){};

    uint getCommunicationProcess(int direction) override
    {
        uint process = 0;
        if (direction != CommunicationDirections::MX)
            process = (uint)INVALID_INDEX;
        return process;
    }

    uint getNumberOfGridLevels() const override
    {
        return 2;
    }

    uint getNumberOfSendIndices(int direction, uint level) override
    {
        return numberOfSendIndices;
    }

    uint getNumberOfReceiveIndices(int direction, uint level) override
    {
        return 0;
    }

    void getSendIndices(int *sendIndices, int direction, int level) override
    {
    }

    void getReceiveIndices(int *sendIndices, int direction, int level) override
    {
    }
};

class CudaMemoryManagerDouble : public CudaMemoryManager
{
public:
    explicit CudaMemoryManagerDouble(std::shared_ptr<Parameter> parameter) : CudaMemoryManager(parameter){};

    void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor) override{};
    void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor) override{};
};

class IndexRearrangementForStreamsDouble : public IndexRearrangementForStreams
{
public:
    IndexRearrangementForStreamsDouble(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder,
                                       vf::parallel::Communicator &communicator)
        : IndexRearrangementForStreams(para, builder, communicator){};

    void initCommunicationArraysForCommAfterFinetoCoarseX(uint level, int indexOfProcessNeighbor,
                                                          int direction) const override{};
    void initCommunicationArraysForCommAfterFinetoCoarseY(uint level, int indexOfProcessNeighbor,
                                                          int direction) const override{};
    void initCommunicationArraysForCommAfterFinetoCoarseZ(uint level, int indexOfProcessNeighbor,
                                                          int direction) const override{};
};

} // namespace GridGeneratorTest

using namespace GridGeneratorTest;

class GridGeneratorTests_initalValuesDomainDecompostion : public testing::Test
{
public:
    void act() const
    {
        gridGenerator->initalValuesDomainDecompostion();
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderStub> builder;

    const uint level = 1;
    const uint direction = CommunicationDirections::MX;

    SPtr<GridGenerator> gridGenerator;

private:
    void SetUp() override
    {
        para = std::make_shared<Parameter>();
        para->setMaxLevel(level + 1); // setMaxLevel resizes parH and parD
        for (uint i = 0; i <= level; i++) {
            para->parH[i] = std::make_shared<LBMSimulationParameter>();
            para->parD[i] = std::make_shared<LBMSimulationParameter>();
        }
        para->setNumprocs(2);

        builder = std::make_shared<LevelGridBuilderStub>(nullptr);
        vf::parallel::NullCommunicator communicator;

        gridGenerator = std::make_shared<GridGenerator>(builder, para, std::make_shared<CudaMemoryManagerDouble>(para),
                                                        communicator);
        gridGenerator->setIndexRearrangementForStreams(
            std::make_unique<IndexRearrangementForStreamsDouble>(para, builder, communicator));
    }
};

TEST_F(GridGeneratorTests_initalValuesDomainDecompostion, whenNoCommunication_sendProcessNeighborShouldNotExist)
{
    act();
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborX.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborY.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborZ.size(), testing::Eq(0));
}

TEST_F(GridGeneratorTests_initalValuesDomainDecompostion, whenCommunicationInX_sendProcessNeighborShouldExistInX)
{
    builder->numberOfSendIndices = 1;
    act();
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborX.size(),
                testing::Eq(1)); // one entry for CommunicationDirections::MX
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborY.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborZ.size(), testing::Eq(0));
}
