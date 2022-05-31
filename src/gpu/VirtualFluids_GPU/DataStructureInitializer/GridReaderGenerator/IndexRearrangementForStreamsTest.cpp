#include <gmock/gmock.h>

#include <algorithm>
#include <filesystem>
#include <iostream>

#include <Parameter/Parameter.h>
#include <basics/config/ConfigurationFile.h>

#include <DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h>
#include <gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h>
#include <gpu/GridGenerator/grid/GridImp.h>

class LevelGridBuilderDouble : public LevelGridBuilder
{
private:
    SPtr<Grid> grid;
    LevelGridBuilderDouble() = default;

public:
    LevelGridBuilderDouble(SPtr<Grid> grid) : LevelGridBuilder(Device(), ""), grid(grid){};
    SPtr<Grid> getGrid(uint level) override { return grid; };
    std::shared_ptr<Grid> getGrid(int level, int box) override { return grid; };
};

class GridImpDouble : public GridImp
{
private:
    std::vector<uint> fluidNodeIndicesBorder;

public:
    GridImpDouble(Object *object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                  SPtr<GridStrategy> gridStrategy, Distribution d, uint level)
        : GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, d, level)
    {
    }

    static SPtr<GridImpDouble> makeShared(Object *object, real startX, real startY, real startZ, real endX, real endY,
                                          real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d,
                                          uint level)
    {
        SPtr<GridImpDouble> grid(
            new GridImpDouble(object, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, d, level));
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

template <typename T>
bool vectorsAreEqual(T *vector1, std::vector<T> vectorExpected)
{
    for (uint i = 0; i < vectorExpected.size(); i++) {
        if (vector1[i] != vectorExpected[i])
            return false;
    }
    return true;
}

class IndexRearrangementForStreamsTest_IndicesCFBorderBulkTest : public testing::Test
{
public:
    CFBorderBulk cf;
    SPtr<Parameter> para;

private:
    static std::unique_ptr<IndexRearrangementForStreams> createTestSubjectCFBorderBulk(CFBorderBulk &cf,
                                                                                       std::shared_ptr<Parameter> para)
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, nullptr, Distribution(), 1);
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

        return std::make_unique<IndexRearrangementForStreams>(para, builder);
    };

    void SetUp() override
    {
        para             = initParameterClass();
        auto testSubject = createTestSubjectCFBorderBulk(cf, para);
        testSubject->splitCoarseToFineIntoBorderAndBulk(cf.level);
    }
};

TEST_F(IndexRearrangementForStreamsTest_IndicesCFBorderBulkTest, splitCoarseToFineIntoBorderAndBulk)
{
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
