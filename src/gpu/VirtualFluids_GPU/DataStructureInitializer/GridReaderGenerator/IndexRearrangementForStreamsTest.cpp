#include <gmock/gmock.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <mpi.h>

#include "DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h"
#include "Parameter/Parameter.h"
#include "basics/config/ConfigurationFile.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"
#include "gpu/VirtualFluids_GPU/Communication/Communicator.cpp"

template <typename T>
bool vectorsAreEqual(T *vector1, std::vector<T> vectorExpected)
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

struct SendIndicesForCommAfterFtoCX {
    // data to work on
    std::vector<int> sendIndices = { 10, 11, 12, 13, 14, 15, 16 };
    int level = 0;
    int direction = CommunicationDirections::MX;
    int numberOfProcessNeighbors = 1;
    int indexOfProcessNeighbor = 0;

    std::vector<uint> iCellCFC = { 8, 10, 12 };
    std::vector<uint> iCellFCC = { 14, 16, 18 };
    uint kCF = (uint)iCellCFC.size();
    uint kFC = (uint)iCellFCC.size();
    uint neighborX[18] = { 0u };
    uint neighborY[18] = { 0u };
    uint neighborZ[18] = { 0u };

    // output data
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;

    // expected data
    std::vector<uint> sendIndicesForCommAfterFtoCPositions_expected = { 4, 6, 0, 2 };
    std::vector<int> sendProcessNeighborX_expected = { 14, 16, 10, 12, 11, 13, 15 };
    int numberOfSendNodesAfterFtoC_expected = (int)sendIndicesForCommAfterFtoCPositions_expected.size();
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
    std::unique_ptr<IndexRearrangementForStreams> createTestSubjectReorderSendIndices()
    {
        logging::Logger::addStream(&std::cout);

        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);

        builder->setNumberOfSendIndices((uint)si.sendIndices.size());
        para->setMaxLevel(si.level + 1); // setMaxLevel resizes parH and parD
        para->parH[si.level] = std::make_shared<LBMSimulationParameter>();
        para->parD[si.level] = std::make_shared<LBMSimulationParameter>();

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
        testSubject = createTestSubjectReorderSendIndices();
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
