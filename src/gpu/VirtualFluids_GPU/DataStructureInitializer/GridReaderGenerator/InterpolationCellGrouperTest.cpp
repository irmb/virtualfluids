#include <gmock/gmock.h>
#include "Utilities/testUtilitiesGPU.h"

#include <iostream>

#include "DataStructureInitializer/GridReaderGenerator/InterpolationCellGrouper.h"
#include "Parameter/Parameter.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"

template <typename T>
bool vectorsAreEqual(const T * vector1, const std::vector<T>& vectorExpected)
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

public:
    explicit LevelGridBuilderDouble(SPtr<Grid> grid) : LevelGridBuilder(), grid(grid){};
    SPtr<Grid> getGrid(uint) override
    {
        return grid;
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

    static SPtr<GridImpDouble> makeShared()
    {
        SPtr<GridImpDouble> grid(new GridImpDouble(nullptr, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, Distribution(), 1));
        return grid;
    }

    void setFluidNodeIndicesBorder(const std::vector<uint>& fluidNodeIndicesBorder)
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
    std::vector<uint> iCellCFC = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<uint> iCellCFF = { 2, 12, 4, 14, 6, 16, 8 };
    const uint sizeOfICellCf = (uint)iCellCFC.size();
    uint neighborX[17] = { 0u };
    uint neighborY[17] = { 0u };
    uint neighborZ[17] = { 0u };
    const int level = 0;
    std::vector<real> offsetCFx = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<real> offsetCFy = { 101, 111, 103, 113, 105, 115, 107 };
    std::vector<real> offsetCFz = { 1001, 1011, 1003, 1013, 1005, 1015, 1007 };

    // expected data
    std::vector<uint> iCellCfcBorder_expected = { 11, 13, 15 };
    std::vector<uint> iCellCfcBulk_expected = { 1, 3, 5, 7 };
    std::vector<uint> iCellCffBorder_expected = { 12, 14, 16 };
    std::vector<uint> iCellCffBulk_expected = { 2, 4, 6, 8 };
    std::vector<real> offsetCFx_Border_expected = { 11, 13, 15 };
    std::vector<real> offsetCFx_Bulk_expected = { 1, 3, 5, 7 };
    std::vector<real> offsetCFy_Border_expected = { 111, 113, 115 };
    std::vector<real> offsetCFy_Bulk_expected = { 101, 103, 105, 107 };
    std::vector<real> offsetCFz_Border_expected = { 1011, 1013, 1015 };
    std::vector<real> offsetCFz_Bulk_expected = { 1001, 1003, 1005, 1007 };
};

class InterpolationCellGrouperTest_IndicesCFBorderBulkTest : public testing::Test
{
protected:
    CFBorderBulk cf;
    SPtr<Parameter> para;
    std::unique_ptr<InterpolationCellGrouper> testSubject;

private:
    std::unique_ptr<InterpolationCellGrouper> createTestSubjectCFBorderBulk()
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared();
        grid->setFluidNodeIndicesBorder(cf.fluidNodeIndicesBorder);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);

        para = testingVF::createParameterForLevel(cf.level);
        para->getParH(cf.level)->coarseToFine.coarseCellIndices = &(cf.iCellCFC.front());
        para->getParH(cf.level)->coarseToFine.fineCellIndices = &(cf.iCellCFF.front());
        para->getParH(cf.level)->neighborX = cf.neighborX;
        para->getParH(cf.level)->neighborY = cf.neighborY;
        para->getParH(cf.level)->neighborZ = cf.neighborZ;
        para->getParH(cf.level)->coarseToFine.numberOfCells = cf.sizeOfICellCf;
        para->getParH(cf.level)->neighborCF.x = &(cf.offsetCFx.front());
        para->getParH(cf.level)->neighborCF.y = &(cf.offsetCFy.front());
        para->getParH(cf.level)->neighborCF.z = &(cf.offsetCFz.front());

        return std::make_unique<InterpolationCellGrouper>(para->getParHallLevels(), para->getParDallLevels(), builder);
    };

    void SetUp() override
    {
        testSubject = createTestSubjectCFBorderBulk();
    }
};

TEST_F(InterpolationCellGrouperTest_IndicesCFBorderBulkTest, splitCoarseToFineIntoBorderAndBulk)
{
    testSubject->splitCoarseToFineIntoBorderAndBulk(cf.level);

    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBorder.numberOfCells + para->getParH(cf.level)->coarseToFineBulk.numberOfCells,
                testing::Eq(cf.sizeOfICellCf))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBorder.numberOfCells, testing::Eq((uint)cf.iCellCfcBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBorder.coarseCellIndices, cf.iCellCfcBorder_expected))
        << "coarseToFineBorder.ICellCFC does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBorder.numberOfCells, testing::Eq((uint)cf.iCellCffBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBorder.fineCellIndices, cf.iCellCffBorder_expected))
        << "coarseToFineBorder.ICellCFF does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBulk.numberOfCells, testing::Eq((uint)cf.iCellCfcBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBulk.coarseCellIndices, cf.iCellCfcBulk_expected))
        << "coarseToFineBulk.ICellCFC does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBulk.numberOfCells, testing::Eq((uint)cf.iCellCffBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBulk.fineCellIndices, cf.iCellCffBulk_expected))
        << "coarseToFineBulk.ICellCFF does not match the expected bulk vector";

    // check offset cells
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCF.x, cf.offsetCFx_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCFBulk.x, cf.offsetCFx_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCF.y, cf.offsetCFy_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCFBulk.y, cf.offsetCFy_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCF.z, cf.offsetCFz_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCFBulk.z, cf.offsetCFz_Bulk_expected));
}

struct FCBorderBulk {
    // data to work on
    std::vector<uint> fluidNodeIndicesBorder = { 110, 111, 112, 113, 114, 115, 116 };
    std::vector<uint> iCellFCC = { 11, 111, 13, 113, 15, 115, 17 };
    std::vector<uint> iCellFCF = { 12, 112, 14, 114, 16, 116, 18 };
    const uint sizeOfICellFC = (uint)iCellFCC.size();
    const int level = 1;
    std::vector<real> offsetFCx = { 11, 111, 13, 113, 15, 115, 17 };
    std::vector<real> offsetFCy = { 1101, 1111, 1103, 1113, 1105, 1115, 1107 };
    std::vector<real> offsetFCz = { 11001, 11011, 11003, 11013, 11005, 11015, 11007 };

    // expected data
    std::vector<uint> iCellFccBorder_expected = { 111, 113, 115 };
    std::vector<uint> iCellFccBulk_expected = { 11, 13, 15, 17 };
    std::vector<uint> iCellFcfBorder_expected = { 112, 114, 116 };
    std::vector<uint> iCellFcfBulk_expected = { 12, 14, 16, 18 };
    std::vector<real> offsetFCx_Border_expected = { 111, 113, 115 };
    std::vector<real> offsetFCx_Bulk_expected = { 11, 13, 15, 17 };
    std::vector<real> offsetFCy_Border_expected = { 1111, 1113, 1115 };
    std::vector<real> offsetFCy_Bulk_expected = { 1101, 1103, 1105, 1107 };
    std::vector<real> offsetFCz_Border_expected = { 11011, 11013, 11015 };
    std::vector<real> offsetFCz_Bulk_expected = { 11001, 11003, 11005, 11007 };
};

class InterpolationCellGrouperTest_IndicesFCBorderBulkTest : public testing::Test
{
protected:
    FCBorderBulk fc;
    SPtr<Parameter> para;
    std::unique_ptr<InterpolationCellGrouper> testSubject;

private:
    std::unique_ptr<InterpolationCellGrouper> createTestSubjectFCBorderBulk()
    {
        SPtr<GridImpDouble> grid =
            GridImpDouble::makeShared();
        grid->setFluidNodeIndicesBorder(fc.fluidNodeIndicesBorder);
        std::shared_ptr<LevelGridBuilderDouble> builder = std::make_shared<LevelGridBuilderDouble>(grid);

        para = testingVF::createParameterForLevel(fc.level);
        para->getParH(fc.level)->fineToCoarse.coarseCellIndices = &(fc.iCellFCC.front());
        para->getParH(fc.level)->fineToCoarse.fineCellIndices = &(fc.iCellFCF.front());
        para->getParH(fc.level)->fineToCoarse.numberOfCells = fc.sizeOfICellFC;
        para->getParH(fc.level)->neighborFC.x = &(fc.offsetFCx.front());
        para->getParH(fc.level)->neighborFC.y = &(fc.offsetFCy.front());
        para->getParH(fc.level)->neighborFC.z = &(fc.offsetFCz.front());

        return std::make_unique<InterpolationCellGrouper>(para->getParHallLevels(), para->getParDallLevels(), builder);
    };

    void SetUp() override
    {
        testSubject = createTestSubjectFCBorderBulk();
    }
};

TEST_F(InterpolationCellGrouperTest_IndicesFCBorderBulkTest, splitFineToCoarseIntoBorderAndBulk)
{
    testSubject->splitFineToCoarseIntoBorderAndBulk(fc.level);

    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBorder.numberOfCells + para->getParH(fc.level)->fineToCoarseBulk.numberOfCells,
                testing::Eq(fc.sizeOfICellFC))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBorder.numberOfCells, testing::Eq((uint)fc.iCellFccBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBorder.coarseCellIndices, fc.iCellFccBorder_expected))
        << "fineToCoarseBorder.ICellFCC does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBorder.numberOfCells, testing::Eq((uint)fc.iCellFcfBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBorder.fineCellIndices, fc.iCellFcfBorder_expected))
        << "fineToCoarseBorder.ICellFCF does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBulk.numberOfCells, testing::Eq((uint)fc.iCellFccBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBulk.coarseCellIndices, fc.iCellFccBulk_expected))
        << "fineToCoarseBulk.ICellFCC does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBulk.numberOfCells, testing::Eq((uint)fc.iCellFcfBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBulk.fineCellIndices, fc.iCellFcfBulk_expected))
        << "fineToCoarseBulk.ICellFCF does not match the expected bulk vector";

    // check offset cells
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFC.x, fc.offsetFCx_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFCBulk.x, fc.offsetFCx_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFC.y, fc.offsetFCy_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFCBulk.y, fc.offsetFCy_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFC.z, fc.offsetFCz_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFCBulk.z, fc.offsetFCz_Bulk_expected));
}
