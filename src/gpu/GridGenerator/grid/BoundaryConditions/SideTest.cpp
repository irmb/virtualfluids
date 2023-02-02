#include "Side.h"
#include "PointerDefinitions.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"
#include "lbm/constants/D3Q27.h"
#include "gmock/gmock.h"
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

using namespace vf::gpu;
using namespace vf::lbm::dir;

class SideForTest : public Side
{

public:
    void setQs(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, uint index)
    {
        Side::setQs(grid, boundaryCondition, index);
    };
    int sideDirection = POSITIVE_DIR;
    int coordinateDirection = X_INDEX;
    SideType mySide = SideType::PX;

private:
    void correctNeighborForPeriodicBoundaries(Grid *grid, real x, real y, real z, real *coords, real neighborX,
                                              real neighborY, real neighborZ) const override
    {
    }

    int getDirection() const override
    {
        return sideDirection;
    }

    void addIndices(std::vector<SPtr<Grid>> grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override
    {
    }

    int getCoordinate() const override
    {
        return coordinateDirection;
    }

    SideType whoAmI() const override
    {
        return mySide;
    }
};

class GridDouble : public GridImp
{

public:
    int endDirection = -1;

    GridDouble()
    {
        this->distribution = DistributionHelper::getDistribution27();
    }

    void transIndexToCoords(uint index, real &x, real &y, real &z) const override
    {
        x = 0;
        y = 0;
        z = 0;
    }

    real getDelta() const override
    {
        return 1.0;
    }

    uint transCoordToIndex(const real &x, const real &y, const real &z) const override
    {
        return 0;
    }

    char getFieldEntry(uint /*matrixIndex*/) const override
    {
        return STOPPER_OUT_OF_GRID_BOUNDARY;
    }

    int getEndDirection() const override
    {
        return endDirection;
    }
};

class BoundaryConditionSpy : public gg::BoundaryCondition
{
public:
    char getType() const override
    {
        return 't';
    };
    const std::vector<std::vector<real>> &getQs()
    {
        return this->qs;
    }
    void resetQVector()
    {
        this->qs.clear();
    }
};

TEST(SideTest, setQs2D_DiagonalEdgeNodesAreSetOnlyOnce)
{
    SideForTest side;
    SPtr<GridDouble> grid = std::make_shared<GridDouble>();
    SPtr<BoundaryConditionSpy> bc = std::make_shared<BoundaryConditionSpy>();
    uint index = 0;
    grid->endDirection = 10;

    // first bc for this node

    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(11, -1);
    expectedFirstBC[DIR_P00] = 0.5;
    expectedFirstBC[DIR_PP0] = 0.5;
    expectedFirstBC[DIR_PM0] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in PY direction

    bc->resetQVector();
    grid->getBCAlreadySet().push_back(SideType::PY);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ(11, -1);
    expectedPreviousQ[DIR_P00] = 0.5;
    expectedPreviousQ[DIR_PM0] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));
}

class SideTestBC : public testing::Test
{
protected:
    SideForTest side;
    SPtr<GridDouble> grid = std::make_shared<GridDouble>();
    SPtr<BoundaryConditionSpy> bc = std::make_shared<BoundaryConditionSpy>();
    uint index = 0;

    std::vector<real> noBC;

    void SetUp() override
    {
        grid->endDirection = 26;
    }
};

TEST_F(SideTestBC, setQs3D_DiagonalEdgeNodesAreNotInfluencedByOppositeNodes)
{
    // first bc for this node

    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(27, -1);
    expectedFirstBC[DIR_P00] = 0.5;
    expectedFirstBC[DIR_PP0] = 0.5;
    expectedFirstBC[DIR_PM0] = 0.5;
    expectedFirstBC[DIR_P0P] = 0.5;
    expectedFirstBC[DIR_P0M] = 0.5;
    expectedFirstBC[DIR_PPP] = 0.5;
    expectedFirstBC[DIR_PMP] = 0.5;
    expectedFirstBC[DIR_PPM] = 0.5;
    expectedFirstBC[DIR_PMM] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in MX direction, does not change anything

    bc->resetQVector();
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedFirstBC));
}

TEST_F(SideTestBC, setQs3D_DiagonalEdgeNodesAreSetOnlyOnce_setBcPX)
{
    // first bc for this node
    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(27, -1);
    expectedFirstBC[DIR_P00] = 0.5;
    expectedFirstBC[DIR_PP0] = 0.5;
    expectedFirstBC[DIR_PM0] = 0.5;
    expectedFirstBC[DIR_P0P] = 0.5;
    expectedFirstBC[DIR_P0M] = 0.5;
    expectedFirstBC[DIR_PPP] = 0.5;
    expectedFirstBC[DIR_PMP] = 0.5;
    expectedFirstBC[DIR_PPM] = 0.5;
    expectedFirstBC[DIR_PMM] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in PY direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ(27, -1);
    expectedPreviousQ[DIR_P00] = 0.5;
    expectedPreviousQ[DIR_PM0] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_PMP] = 0.5;
    expectedPreviousQ[DIR_PMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MY direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_P00] = 0.5;
    expectedPreviousQ[DIR_PP0] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_PPP] = 0.5;
    expectedPreviousQ[DIR_PPM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));
    
    // node already has BC in PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_P00] = 0.5;
    expectedPreviousQ[DIR_PP0] = 0.5;
    expectedPreviousQ[DIR_PM0] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_PPM] = 0.5;
    expectedPreviousQ[DIR_PMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_P00] = 0.5;
    expectedPreviousQ[DIR_PP0] = 0.5;
    expectedPreviousQ[DIR_PM0] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_PPP] = 0.5;
    expectedPreviousQ[DIR_PMP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PY and MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    auto qsTwoPreviousQs = bc->getQs()[0];

    std::vector<real> expectedTwoPreviousQs(27, -1);
    expectedTwoPreviousQs[DIR_P00] = 0.5;
    expectedTwoPreviousQs[DIR_PM0] = 0.5;
    expectedTwoPreviousQs[DIR_P0P] = 0.5;
    expectedTwoPreviousQs[DIR_PMP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_P00] = 0.5;
    expectedTwoPreviousQs[DIR_PM0] = 0.5;
    expectedTwoPreviousQs[DIR_P0M] = 0.5;
    expectedTwoPreviousQs[DIR_PMM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_P00] = 0.5;
    expectedTwoPreviousQs[DIR_PP0] = 0.5;
    expectedTwoPreviousQs[DIR_P0M] = 0.5;
    expectedTwoPreviousQs[DIR_PPM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_P00] = 0.5;
    expectedTwoPreviousQs[DIR_PP0] = 0.5;
    expectedTwoPreviousQs[DIR_P0P] = 0.5;
    expectedTwoPreviousQs[DIR_PPP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));
}

TEST_F(SideTestBC, setQs3D_DiagonalEdgeNodesAreSetOnlyOnce_setBcMX)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    // first bc for this node
    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(27, -1);
    expectedFirstBC[DIR_M00] = 0.5;
    expectedFirstBC[DIR_MP0] = 0.5;
    expectedFirstBC[DIR_MM0] = 0.5;
    expectedFirstBC[DIR_M0P] = 0.5;
    expectedFirstBC[DIR_M0M] = 0.5;
    expectedFirstBC[DIR_MPP] = 0.5;
    expectedFirstBC[DIR_MMP] = 0.5;
    expectedFirstBC[DIR_MPM] = 0.5;
    expectedFirstBC[DIR_MMM] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in PY direction

    bc->resetQVector();
    grid->getBCAlreadySet().push_back(SideType::PY);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ(27, -1);
    expectedPreviousQ[DIR_M00] = 0.5;
    expectedPreviousQ[DIR_MM0] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_MMP] = 0.5;
    expectedPreviousQ[DIR_MMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MY direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_M00] = 0.5;
    expectedPreviousQ[DIR_MP0] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_MPP] = 0.5;
    expectedPreviousQ[DIR_MPM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_M00] = 0.5;
    expectedPreviousQ[DIR_MP0] = 0.5;
    expectedPreviousQ[DIR_MM0] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_MPM] = 0.5;
    expectedPreviousQ[DIR_MMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_M00] = 0.5;
    expectedPreviousQ[DIR_MP0] = 0.5;
    expectedPreviousQ[DIR_MM0] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_MPP] = 0.5;
    expectedPreviousQ[DIR_MMP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PY and MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    auto qsTwoPreviousQs = bc->getQs()[0];

    std::vector<real> expectedTwoPreviousQs(27, -1);
    expectedTwoPreviousQs[DIR_M00] = 0.5;
    expectedTwoPreviousQs[DIR_MM0] = 0.5;
    expectedTwoPreviousQs[DIR_M0P] = 0.5;
    expectedTwoPreviousQs[DIR_MMP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_M00] = 0.5;
    expectedTwoPreviousQs[DIR_MM0] = 0.5;
    expectedTwoPreviousQs[DIR_M0M] = 0.5;
    expectedTwoPreviousQs[DIR_MMM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and PZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::PZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_M00] = 0.5;
    expectedTwoPreviousQs[DIR_MP0] = 0.5;
    expectedTwoPreviousQs[DIR_M0M] = 0.5;
    expectedTwoPreviousQs[DIR_MPM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and MZ direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::MZ);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_M00] = 0.5;
    expectedTwoPreviousQs[DIR_MP0] = 0.5;
    expectedTwoPreviousQs[DIR_M0P] = 0.5;
    expectedTwoPreviousQs[DIR_MPP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));
}

TEST_F(SideTestBC, setQs3D_DiagonalEdgeNodesAreSetOnlyOnce_setBcMZ)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    // first bc for this node
    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(27, -1);
    expectedFirstBC[DIR_00M] = 0.5;
    expectedFirstBC[DIR_P0M] = 0.5;
    expectedFirstBC[DIR_M0M] = 0.5;
    expectedFirstBC[DIR_0PM] = 0.5;
    expectedFirstBC[DIR_0MM] = 0.5;
    expectedFirstBC[DIR_PPM] = 0.5;
    expectedFirstBC[DIR_MPM] = 0.5;
    expectedFirstBC[DIR_PMM] = 0.5;
    expectedFirstBC[DIR_MMM] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in MY direction

    bc->resetQVector();
    grid->getBCAlreadySet().push_back(SideType::MY);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ(27, -1);
    expectedPreviousQ[DIR_00M] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_0PM] = 0.5;
    expectedPreviousQ[DIR_PPM] = 0.5;
    expectedPreviousQ[DIR_MPM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PY direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00M] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_0MM] = 0.5;
    expectedPreviousQ[DIR_PMM] = 0.5;
    expectedPreviousQ[DIR_MMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00M] = 0.5;
    expectedPreviousQ[DIR_M0M] = 0.5;
    expectedPreviousQ[DIR_0PM] = 0.5;
    expectedPreviousQ[DIR_0MM] = 0.5;
    expectedPreviousQ[DIR_MPM] = 0.5;
    expectedPreviousQ[DIR_MMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00M] = 0.5;
    expectedPreviousQ[DIR_P0M] = 0.5;
    expectedPreviousQ[DIR_0PM] = 0.5;
    expectedPreviousQ[DIR_0MM] = 0.5;
    expectedPreviousQ[DIR_PPM] = 0.5;
    expectedPreviousQ[DIR_PMM] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MY and PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    auto qsTwoPreviousQs = bc->getQs()[0];

    std::vector<real> expectedTwoPreviousQs(27, -1);
    expectedTwoPreviousQs[DIR_00M] = 0.5;
    expectedTwoPreviousQs[DIR_M0M] = 0.5;
    expectedTwoPreviousQs[DIR_0PM] = 0.5;
    expectedTwoPreviousQs[DIR_MPM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00M] = 0.5;
    expectedTwoPreviousQs[DIR_P0M] = 0.5;
    expectedTwoPreviousQs[DIR_0PM] = 0.5;
    expectedTwoPreviousQs[DIR_PPM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00M] = 0.5;
    expectedTwoPreviousQs[DIR_M0M] = 0.5;
    expectedTwoPreviousQs[DIR_0MM] = 0.5;
    expectedTwoPreviousQs[DIR_MMM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00M] = 0.5;
    expectedTwoPreviousQs[DIR_P0M] = 0.5;
    expectedTwoPreviousQs[DIR_0MM] = 0.5;
    expectedTwoPreviousQs[DIR_PMM] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));
}

TEST_F(SideTestBC, setQs3D_DiagonalEdgeNodesAreSetOnlyOnce_setBcPZ)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;

    // first bc for this node
    side.setQs(grid, bc, index);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC(27, -1);
    expectedFirstBC[DIR_00P] = 0.5;
    expectedFirstBC[DIR_P0P] = 0.5;
    expectedFirstBC[DIR_M0P] = 0.5;
    expectedFirstBC[DIR_0PP] = 0.5;
    expectedFirstBC[DIR_0MP] = 0.5;
    expectedFirstBC[DIR_PPP] = 0.5;
    expectedFirstBC[DIR_MPP] = 0.5;
    expectedFirstBC[DIR_PMP] = 0.5;
    expectedFirstBC[DIR_MMP] = 0.5;
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has BC in MY direction

    bc->resetQVector();
    grid->getBCAlreadySet().push_back(SideType::MY);

    side.setQs(grid, bc, index);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ(27, -1);
    expectedPreviousQ[DIR_00P] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_0PP] = 0.5;
    expectedPreviousQ[DIR_PPP] = 0.5;
    expectedPreviousQ[DIR_MPP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PY direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00P] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_0MP] = 0.5;
    expectedPreviousQ[DIR_PMP] = 0.5;
    expectedPreviousQ[DIR_MMP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00P] = 0.5;
    expectedPreviousQ[DIR_M0P] = 0.5;
    expectedPreviousQ[DIR_0PP] = 0.5;
    expectedPreviousQ[DIR_0MP] = 0.5;
    expectedPreviousQ[DIR_MPP] = 0.5;
    expectedPreviousQ[DIR_MMP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsPreviousQ = bc->getQs()[0];

    expectedPreviousQ.clear();
    expectedPreviousQ.resize(27, -1);
    expectedPreviousQ[DIR_00P] = 0.5;
    expectedPreviousQ[DIR_P0P] = 0.5;
    expectedPreviousQ[DIR_0PP] = 0.5;
    expectedPreviousQ[DIR_0MP] = 0.5;
    expectedPreviousQ[DIR_PPP] = 0.5;
    expectedPreviousQ[DIR_PMP] = 0.5;
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // node already has BC in MY and PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    auto qsTwoPreviousQs = bc->getQs()[0];

    std::vector<real> expectedTwoPreviousQs(27, -1);
    expectedTwoPreviousQs[DIR_00P] = 0.5;
    expectedTwoPreviousQs[DIR_M0P] = 0.5;
    expectedTwoPreviousQs[DIR_0PP] = 0.5;
    expectedTwoPreviousQs[DIR_MPP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in MY and MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::MY);
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00P] = 0.5;
    expectedTwoPreviousQs[DIR_P0P] = 0.5;
    expectedTwoPreviousQs[DIR_0PP] = 0.5;
    expectedTwoPreviousQs[DIR_PPP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and PX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::PX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00P] = 0.5;
    expectedTwoPreviousQs[DIR_M0P] = 0.5;
    expectedTwoPreviousQs[DIR_0MP] = 0.5;
    expectedTwoPreviousQs[DIR_MMP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));

    // node already has BC in PY and MX direction

    bc->resetQVector();
    grid->getBCAlreadySet().clear();
    grid->getBCAlreadySet().push_back(SideType::PY);
    grid->getBCAlreadySet().push_back(SideType::MX);

    side.setQs(grid, bc, index);
    qsTwoPreviousQs = bc->getQs()[0];

    expectedTwoPreviousQs.clear();
    expectedTwoPreviousQs.resize(27, -1);
    expectedTwoPreviousQs[DIR_00P] = 0.5;
    expectedTwoPreviousQs[DIR_P0P] = 0.5;
    expectedTwoPreviousQs[DIR_0MP] = 0.5;
    expectedTwoPreviousQs[DIR_PMP] = 0.5;
    EXPECT_THAT(qsTwoPreviousQs, testing::Eq(expectedTwoPreviousQs));
}
