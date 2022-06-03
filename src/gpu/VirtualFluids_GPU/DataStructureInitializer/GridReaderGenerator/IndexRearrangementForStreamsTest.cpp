#include <gmock/gmock.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <mpi.h>

#include "Parameter/Parameter.h"
#include "basics/config/ConfigurationFile.h"
#include "DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h"
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
    SPtr<Grid> getGrid(uint level) override { return grid; };
    std::shared_ptr<Grid> getGrid(int level, int box) override { return grid; };
    void setNumberOfSendIndices(uint numberOfSendIndices) { this->numberOfSendIndices = numberOfSendIndices; };
    uint getNumberOfSendIndices(int direction, uint level) override { return numberOfSendIndices; };
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
                                          real endZ, real delta, Distribution d,
                                          uint level)
    {
        SPtr<GridImpDouble> grid(
            new GridImpDouble(object, startX, startY, startZ, endX, endY, endZ, delta, d, level));
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

struct CFBorderBulk {
    // data to work on
    std::vector<uint> fluidNodeIndicesBorder = { 10, 11, 12, 13, 14, 15, 16 };
    std::vector<uint> iCellCFC               = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<uint> iCellCFF               = { 2, 12, 4, 14, 6, 16, 8 };
    uint sizeOfICellCf                       = (uint)iCellCFC.size();
    uint neighborX_SP[17]                    = { 0u };
    uint neighborY_SP[17]                    = { 0u };
    uint neighborZ_SP[17]                    = { 0u };
    int level                                = 0;
    std::vector<real> offsetCFx              = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<real> offsetCFy              = { 101, 111, 103, 113, 105, 115, 107 };
    std::vector<real> offsetCFz              = { 1001, 1011, 1003, 1013, 1005, 1015, 1007 };

    // expected data
    std::vector<uint> iCellCfcBorder_expected   = { 11, 13, 15 };
    std::vector<uint> iCellCfcBulk_expected     = { 1, 3, 5, 7 };
    std::vector<uint> iCellCffBorder_expected   = { 12, 14, 16 };
    std::vector<uint> iCellCffBulk_expected     = { 2, 4, 6, 8 };
    std::vector<real> offsetCFx_Border_expected = { 11, 13, 15 };
    std::vector<real> offsetCFx_Bulk_expected   = { 1, 3, 5, 7 };
    std::vector<real> offsetCFy_Border_expected = { 111, 113, 115 };
    std::vector<real> offsetCFy_Bulk_expected   = { 101, 103, 105, 107 };
    std::vector<real> offsetCFz_Border_expected = { 1011, 1013, 1015 };
    std::vector<real> offsetCFz_Bulk_expected   = { 1001, 1003, 1005, 1007 };
};

static SPtr<Parameter> initParameterClass()
{
    std::filesystem::path filePath = __FILE__; //  assuming that the config file is stored parallel to this file.
    filePath.replace_filename("IndexRearrangementForStreamsTest.cfg");
    vf::basics::ConfigurationFile config;
    config.load(filePath.string());
    return std::make_shared<Parameter>(config, 1, 0);
}

class IndexRearrangementForStreamsTest_IndicesCFBorderBulkTest : public testing::Test
{
protected:
    CFBorderBulk cf;
    SPtr<Parameter> para;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

private:
    std::unique_ptr<IndexRearrangementForStreams> createTestSubjectCFBorderBulk()
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        grid->setFluidNodeIndicesBorder(cf.fluidNodeIndicesBorder);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);

        para->setMaxLevel(cf.level + 1); // setMaxLevel resizes parH and parD
        para->parH[cf.level]                    = std::make_shared<LBMSimulationParameter>();
        para->parD[cf.level]                    = std::make_shared<LBMSimulationParameter>();
        para->getParH(cf.level)->intCF.ICellCFC = &(cf.iCellCFC.front());
        para->getParH(cf.level)->intCF.ICellCFF = &(cf.iCellCFF.front());
        para->getParH(cf.level)->neighborX_SP   = cf.neighborX_SP;
        para->getParH(cf.level)->neighborY_SP   = cf.neighborY_SP;
        para->getParH(cf.level)->neighborZ_SP   = cf.neighborZ_SP;
        para->getParH(cf.level)->intCF.kCF      = cf.sizeOfICellCf;
        para->getParH(cf.level)->offCF.xOffCF   = &(cf.offsetCFx.front());
        para->getParH(cf.level)->offCF.yOffCF   = &(cf.offsetCFy.front());
        para->getParH(cf.level)->offCF.zOffCF   = &(cf.offsetCFz.front());

        return std::make_unique<IndexRearrangementForStreams>(para, builder, vf::gpu::Communicator::getInstance());
    };

    void SetUp() override
    {
        para        = initParameterClass();
        testSubject = createTestSubjectCFBorderBulk();
    }
};

TEST_F(IndexRearrangementForStreamsTest_IndicesCFBorderBulkTest, splitCoarseToFineIntoBorderAndBulk)
{
    testSubject->splitCoarseToFineIntoBorderAndBulk(cf.level);

    EXPECT_THAT(para->getParH(cf.level)->intCFBorder.kCF + para->getParH(cf.level)->intCFBulk.kCF,
                testing::Eq(cf.sizeOfICellCf))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->intCFBorder.kCF, testing::Eq((uint)cf.iCellCfcBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->intCFBorder.ICellCFC, cf.iCellCfcBorder_expected))
        << "intCFBorder.ICellCFC does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->intCFBorder.kCF, testing::Eq((uint)cf.iCellCffBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->intCFBorder.ICellCFF, cf.iCellCffBorder_expected))
        << "intCFBorder.ICellCFF does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->intCFBulk.kCF, testing::Eq((uint)cf.iCellCfcBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->intCFBulk.ICellCFC, cf.iCellCfcBulk_expected))
        << "intCFBulk.ICellCFC does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->intCFBulk.kCF, testing::Eq((uint)cf.iCellCffBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->intCFBulk.ICellCFF, cf.iCellCffBulk_expected))
        << "intCFBulk.ICellCFF does not match the expected bulk vector";

    // check offset cells
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCF.xOffCF, cf.offsetCFx_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCFBulk.xOffCF, cf.offsetCFx_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCF.yOffCF, cf.offsetCFy_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCFBulk.yOffCF, cf.offsetCFy_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCF.zOffCF, cf.offsetCFz_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->offCFBulk.zOffCF, cf.offsetCFz_Bulk_expected));
}

struct FCBorderBulk {
    // data to work on
    std::vector<uint> fluidNodeIndicesBorder = { 110, 111, 112, 113, 114, 115, 116 };
    std::vector<uint> iCellFCC               = { 11, 111, 13, 113, 15, 115, 17 };
    std::vector<uint> iCellFCF               = { 12, 112, 14, 114, 16, 116, 18 };
    uint sizeOfICellFC                       = (uint)iCellFCC.size();
    int level                                = 1;

    // expected data
    std::vector<uint> iCellFccBorder_expected = { 111, 113, 115 };
    std::vector<uint> iCellFccBulk_expected   = { 11, 13, 15, 17 };
    std::vector<uint> iCellFcfBorder_expected = { 112, 114, 116 };
    std::vector<uint> iCellFcfBulk_expected   = { 12, 14, 16, 18 };
};

class IndexRearrangementForStreamsTest_IndicesFCBorderBulkTest : public testing::Test
{
protected:
    FCBorderBulk fc;
    SPtr<Parameter> para;
    std::unique_ptr<IndexRearrangementForStreams> testSubject;

private:
    std::unique_ptr<IndexRearrangementForStreams> createTestSubjectFCBorderBulk()
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1);
        grid->setFluidNodeIndicesBorder(fc.fluidNodeIndicesBorder);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);

        para->setMaxLevel(fc.level + 1); // setMaxLevel resizes parH and parD
        para->parH[fc.level]                    = std::make_shared<LBMSimulationParameter>();
        para->parD[fc.level]                    = std::make_shared<LBMSimulationParameter>();
        para->getParH(fc.level)->intFC.ICellFCC = &(fc.iCellFCC.front());
        para->getParH(fc.level)->intFC.ICellFCF = &(fc.iCellFCF.front());
        para->getParH(fc.level)->intFC.kFC      = fc.sizeOfICellFC;

        return std::make_unique<IndexRearrangementForStreams>(para, builder, vf::gpu::Communicator::getInstance());
    };

    void SetUp() override
    {
        para        = initParameterClass();
        testSubject = createTestSubjectFCBorderBulk();
    }
};

TEST_F(IndexRearrangementForStreamsTest_IndicesFCBorderBulkTest, splitFineToCoarseIntoBorderAndBulk)
{
    testSubject->splitFineToCoarseIntoBorderAndBulk(fc.level);

    EXPECT_THAT(para->getParH(fc.level)->intFCBorder.kFC + para->getParH(fc.level)->intFCBulk.kFC,
                testing::Eq(fc.sizeOfICellFC))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->intFCBorder.kFC, testing::Eq((uint)fc.iCellFccBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->intFCBorder.ICellFCC, fc.iCellFccBorder_expected))
        << "intFCBorder.ICellFCC does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->intFCBorder.kFC, testing::Eq((uint)fc.iCellFcfBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->intFCBorder.ICellFCF, fc.iCellFcfBorder_expected))
        << "intFCBorder.ICellFCF does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->intFCBulk.kFC, testing::Eq((uint)fc.iCellFccBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->intFCBulk.ICellFCC, fc.iCellFccBulk_expected))
        << "intFCBulk.ICellFCC does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->intFCBulk.kFC, testing::Eq((uint)fc.iCellFcfBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->intFCBulk.ICellFCF, fc.iCellFcfBulk_expected))
        << "intFCBulk.ICellFCF does not match the expected bulk vector";
}
struct SendIndicesForCommAfterFtoCX {
    // data to work on
    std::vector<int> sendIndices = { 10, 11, 12, 13, 14, 15, 16 };
    int level                    = 0;
    int direction                = CommunicationDirections::MX;
    int numberOfProcessNeighbors = 1;
    int indexOfProcessNeighbor   = 0;

    std::vector<uint> iCellCFC = { 8, 10, 12 };
    std::vector<uint> iCellFCC = { 14, 16, 18 };
    uint kCF                   = (uint)iCellCFC.size();
    uint kFC                   = (uint)iCellFCC.size();
    uint neighborX_SP[18]      = { 0u };
    uint neighborY_SP[18]      = { 0u };
    uint neighborZ_SP[18]      = { 0u };

    // output data
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;

    // expected data
    std::vector<uint> sendIndicesForCommAfterFtoCPositions_expected = { 4, 6, 0, 2 };
    std::vector<int> sendProcessNeighborX_expected                  = { 14, 16, 10, 12, 11, 13, 15 };
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

        para->getParH(si.level)->intFC.kFC      = si.kFC;
        para->getParH(si.level)->intFC.ICellFCC = &(si.iCellFCC.front());
        para->getParH(si.level)->intCF.ICellCFC = &(si.iCellCFC.front());
        para->getParH(si.level)->intCF.kCF      = si.kCF;
        para->getParH(si.level)->neighborX_SP   = si.neighborX_SP;
        para->getParH(si.level)->neighborY_SP   = si.neighborY_SP;
        para->getParH(si.level)->neighborZ_SP   = si.neighborZ_SP;

        para->setNumberOfProcessNeighborsX(si.numberOfProcessNeighbors, si.level, "send");
        para->getParH(si.level)->sendProcessNeighborX[si.indexOfProcessNeighbor].index = si.sendIndices.data();
        para->initProcessNeighborsAfterFtoCX(si.level);

        return std::make_unique<IndexRearrangementForStreams>(IndexRearrangementForStreams(para, builder, vf::gpu::Communicator::getInstance()));
    };

    void SetUp() override
    {
        para        = initParameterClass();
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