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
    GridImpDouble(SPtr<Object> object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
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

struct CoarseToFineBorderBulk {
    // data to work on
    std::vector<uint> fluidNodeIndicesBorder = { 10, 11, 12, 13, 14, 15, 16 };
    std::vector<uint> intCtoFcoarse = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<uint> fineCellIndices = { 2, 12, 4, 14, 6, 16, 8 };
    const uint sizeOfInterpolationCoarseToFine = (uint)intCtoFcoarse.size();
    uint neighborX[17] = { 0u };
    uint neighborY[17] = { 0u };
    uint neighborZ[17] = { 0u };
    const int level = 0;
    std::vector<real> neighborCFx = { 1, 11, 3, 13, 5, 15, 7 };
    std::vector<real> neighborCFy = { 101, 111, 103, 113, 105, 115, 107 };
    std::vector<real> neighborCFz = { 1001, 1011, 1003, 1013, 1005, 1015, 1007 };

    // expected data
    std::vector<uint> intCtoFcoarseBorder_expected = { 11, 13, 15 };
    std::vector<uint> intCtoFcoarseBulk_expected = { 1, 3, 5, 7 };
    std::vector<uint> fineCellIndicesBorder_expected = { 12, 14, 16 };
    std::vector<uint> fineCellIndicesBulk_expected = { 2, 4, 6, 8 };
    std::vector<real> neighborCFx_Border_expected = { 11, 13, 15 };
    std::vector<real> neighborCFx_Bulk_expected = { 1, 3, 5, 7 };
    std::vector<real> neighborCFy_Border_expected = { 111, 113, 115 };
    std::vector<real> neighborCFy_Bulk_expected = { 101, 103, 105, 107 };
    std::vector<real> neighborCFz_Border_expected = { 1011, 1013, 1015 };
    std::vector<real> neighborCFz_Bulk_expected = { 1001, 1003, 1005, 1007 };
};

class InterpolationCellGrouperTest_IndicesCFBorderBulkTest : public testing::Test
{
protected:
    CoarseToFineBorderBulk cf;
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
        para->getParH(cf.level)->coarseToFine.coarseCellIndices = &(cf.intCtoFcoarse.front());
        para->getParH(cf.level)->coarseToFine.fineCellIndices = &(cf.fineCellIndices.front());
        para->getParH(cf.level)->neighborX = cf.neighborX;
        para->getParH(cf.level)->neighborY = cf.neighborY;
        para->getParH(cf.level)->neighborZ = cf.neighborZ;
        para->getParH(cf.level)->coarseToFine.numberOfCells = cf.sizeOfInterpolationCoarseToFine;
        para->getParH(cf.level)->neighborCoarseToFine.x = &(cf.neighborCFx.front());
        para->getParH(cf.level)->neighborCoarseToFine.y = &(cf.neighborCFy.front());
        para->getParH(cf.level)->neighborCoarseToFine.z = &(cf.neighborCFz.front());

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
                testing::Eq(cf.sizeOfInterpolationCoarseToFine))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBorder.numberOfCells, testing::Eq((uint)cf.intCtoFcoarseBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBorder.coarseCellIndices, cf.intCtoFcoarseBorder_expected))
        << "coarseToFineBorder.intCtoFcoarse does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBorder.numberOfCells, testing::Eq((uint)cf.fineCellIndicesBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBorder.fineCellIndices, cf.fineCellIndicesBorder_expected))
        << "coarseToFineBorder.fineCellIndices does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBulk.numberOfCells, testing::Eq((uint)cf.intCtoFcoarseBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBulk.coarseCellIndices, cf.intCtoFcoarseBulk_expected))
        << "coarseToFineBulk.intCtoFcoarse does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(cf.level)->coarseToFineBulk.numberOfCells, testing::Eq((uint)cf.fineCellIndicesBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->coarseToFineBulk.fineCellIndices, cf.fineCellIndicesBulk_expected))
        << "coarseToFineBulk.fineCellIndices does not match the expected bulk vector";

    // check neighbor cells
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFine.x, cf.neighborCFx_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFineBulk.x, cf.neighborCFx_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFine.y, cf.neighborCFy_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFineBulk.y, cf.neighborCFy_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFine.z, cf.neighborCFz_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(cf.level)->neighborCoarseToFineBulk.z, cf.neighborCFz_Bulk_expected));
}

struct FineToCoarseBorderBulk {
    // data to work on
    std::vector<uint> fluidNodeIndicesBorder = { 110, 111, 112, 113, 114, 115, 116 };
    std::vector<uint> coarseCellIndices = { 11, 111, 13, 113, 15, 115, 17 };
    std::vector<uint> fineCellIndices = { 12, 112, 14, 114, 16, 116, 18 };
    const uint sizeOfIntFineToCoarse = (uint)coarseCellIndices.size();
    const int level = 1;
    std::vector<real> neighborx = { 11, 111, 13, 113, 15, 115, 17 };
    std::vector<real> neighbory = { 1101, 1111, 1103, 1113, 1105, 1115, 1107 };
    std::vector<real> neighborz = { 11001, 11011, 11003, 11013, 11005, 11015, 11007 };

    // expected data
    std::vector<uint> coarseCellIndicesBorder_expected = { 111, 113, 115 };
    std::vector<uint> coarseCellIndicesBulk_expected = { 11, 13, 15, 17 };
    std::vector<uint> fineCellIndicesBorder_expected = { 112, 114, 116 };
    std::vector<uint> fineCellIndicesBulk_expected = { 12, 14, 16, 18 };
    std::vector<real> neighborx_Border_expected = { 111, 113, 115 };
    std::vector<real> neighborx_Bulk_expected = { 11, 13, 15, 17 };
    std::vector<real> neighbory_Border_expected = { 1111, 1113, 1115 };
    std::vector<real> neighbory_Bulk_expected = { 1101, 1103, 1105, 1107 };
    std::vector<real> neighborz_Border_expected = { 11011, 11013, 11015 };
    std::vector<real> neighborz_Bulk_expected = { 11001, 11003, 11005, 11007 };
};

class InterpolationCellGrouperTest_IndicesFCBorderBulkTest : public testing::Test
{
protected:
    FineToCoarseBorderBulk fc;
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
        para->getParH(fc.level)->fineToCoarse.coarseCellIndices = &(fc.coarseCellIndices.front());
        para->getParH(fc.level)->fineToCoarse.fineCellIndices = &(fc.fineCellIndices.front());
        para->getParH(fc.level)->fineToCoarse.numberOfCells = fc.sizeOfIntFineToCoarse;
        para->getParH(fc.level)->neighborFineToCoarse.x = &(fc.neighborx.front());
        para->getParH(fc.level)->neighborFineToCoarse.y = &(fc.neighbory.front());
        para->getParH(fc.level)->neighborFineToCoarse.z = &(fc.neighborz.front());

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
                testing::Eq(fc.sizeOfIntFineToCoarse))
        << "The number of interpolation cells from coarse to fine changed during reordering.";

    // check coarse to fine border (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBorder.numberOfCells, testing::Eq((uint)fc.coarseCellIndicesBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBorder.coarseCellIndices, fc.coarseCellIndicesBorder_expected))
        << "fineToCoarseBorder.coarseCellIndices does not match the expected border vector";
    // check coarse to fine border (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBorder.numberOfCells, testing::Eq((uint)fc.fineCellIndicesBorder_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBorder.fineCellIndices, fc.fineCellIndicesBorder_expected))
        << "fineToCoarseBorder.fineCellIndices does not match the expected border vector";

    // check coarse to fine bulk (coarse nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBulk.numberOfCells, testing::Eq((uint)fc.coarseCellIndicesBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBulk.coarseCellIndices, fc.coarseCellIndicesBulk_expected))
        << "fineToCoarseBulk.coarseCellIndices does not match the expected bulk vector";
    // check coarse to fine bulk (fine nodes)
    EXPECT_THAT(para->getParH(fc.level)->fineToCoarseBulk.numberOfCells, testing::Eq((uint)fc.fineCellIndicesBulk_expected.size()));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->fineToCoarseBulk.fineCellIndices, fc.fineCellIndicesBulk_expected))
        << "fineToCoarseBulk.fineCellIndices does not match the expected bulk vector";

    // check neighbor cells
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarse.x, fc.neighborx_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarseBulk.x, fc.neighborx_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarse.y, fc.neighbory_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarseBulk.y, fc.neighbory_Bulk_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarse.z, fc.neighborz_Border_expected));
    EXPECT_TRUE(vectorsAreEqual(para->getParH(fc.level)->neighborFineToCoarseBulk.z, fc.neighborz_Bulk_expected));
}
