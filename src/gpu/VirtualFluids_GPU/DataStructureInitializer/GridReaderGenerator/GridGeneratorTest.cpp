#include "GridGenerator.h"
#include <gmock/gmock.h>

#include "Communication/Communicator.h"
#include "GPU/CudaMemoryManager.h"
#include "IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"

namespace GridGeneratorTest
{

class LevelGridBuilderStub : public LevelGridBuilder
{
private:
    SPtr<Grid> grid;
    LevelGridBuilderStub() = default;

public:
    LevelGridBuilderStub(SPtr<Grid> grid) : LevelGridBuilder(), grid(grid){};

    uint getCommunicationProcess(int direction) override
    {
        return 0;
    }

    uint getNumberOfGridLevels() const override
    {
        return 2;
    }

    uint getNumberOfSendIndices(int direction, uint level) override
    {
        return 0;
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
    CudaMemoryManagerDouble(std::shared_ptr<Parameter> parameter) : CudaMemoryManager(parameter){};

    void cudaAllocProcessNeighborX(int lev, unsigned int processNeighbor) override{};
    void cudaCopyProcessNeighborXIndex(int lev, unsigned int processNeighbor) override{};
};

class IndexRearrangementForStreamsDouble : public IndexRearrangementForStreams
{
public:
    IndexRearrangementForStreamsDouble(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder,
                                       vf::gpu::Communicator &communicator)
        : IndexRearrangementForStreams(para, builder, communicator){};

    void initCommunicationArraysForCommAfterFinetoCoarseX(uint level, int j, int direction) override {};
    void initCommunicationArraysForCommAfterFinetoCoarseY(uint level, int j, int direction) override {};
    void initCommunicationArraysForCommAfterFinetoCoarseZ(uint level, int j, int direction) override {};
};

} // namespace GridGeneratorTest

using namespace GridGeneratorTest;

class GridGeneratorTests_initalValuesDomainDecompostion : public testing::Test
{
public:
    void act()
    {
        gridGenerator->initalValuesDomainDecompostion();
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderStub> builder;

    uint level = 1;
    uint direction = CommunicationDirections::MX;

    SPtr<GridGenerator> gridGenerator;

private:
    void SetUp() override
    {
        logging::Logger::addStream(&std::cout);
        logging::Logger::setDebugLevel(logging::Logger::WARNING);

        para = std::make_shared<Parameter>();
        para->setMaxLevel(level + 1); // setMaxLevel resizes parH and parD
        for (uint i = 0; i <= level; i++) {
            para->parH[i] = std::make_shared<LBMSimulationParameter>();
            para->parD[i] = std::make_shared<LBMSimulationParameter>();
        }
        para->setNumprocs(2);

        builder = std::make_shared<LevelGridBuilderStub>(nullptr);
        vf::gpu::Communicator &communicator = vf::gpu::Communicator::getInstance();

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
