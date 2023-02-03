#include "Side.h"
#include "PointerDefinitions.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"
#include "lbm/constants/D3Q27.h"
#include "gmock/gmock.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace vf::gpu;
using namespace vf::lbm::dir;

class SideTestSpecificSubclass : public Side
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
    void correctNeighborForPeriodicBoundaries(const Grid *grid, real x, real y, real z, real *coords, real& neighborX,
                                              real& neighborY, real& neighborZ) const override
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

class SideTestBC : public testing::Test
{
protected:
    SideTestSpecificSubclass side;
    SPtr<GridDouble> grid = std::make_shared<GridDouble>();
    SPtr<BoundaryConditionSpy> bc = std::make_shared<BoundaryConditionSpy>();
    uint index = 0;

    std::vector<real> noBC;

    void SetUp() override
    {
        grid->endDirection = 26;
    }
};

TEST_F(SideTestBC, setQs2D_whenSettingPX_setAllQsNormalToBC)
{
    grid->endDirection = 10;
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(11, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs2D_givenPYhasBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    grid->endDirection = 10;
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(11, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMXhasBeenSet_thenSetPX_setAllQsNormalToPX)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;

    // no previous BC on this node

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));

    // node already has BC in MX direction, but this does not change anything

    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    actualQs = bc->getQs()[0];

    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenGeometryBCInVector_thenSetPX_throws)
{
    // do not add Geometry BC to this vector, as it has an invalid normal
    grid->addBCalreadySet(SideType::GEOMETRY);

    EXPECT_THROW(side.setQs(grid, bc, index), std::out_of_range);
}

TEST_F(SideTestBC, setQs3D_whenSettingPX_setAllQsNormalToBC)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYhasBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYhasBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPZhasBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMZhasBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_P00] = 0.5;
    expectedQs[DIR_PP0] = 0.5;
    expectedQs[DIR_PM0] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYandMZhaveBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_P00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PM0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandPZhaveBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_P00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PM0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandPZhaveBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_P00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PP0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandMZhaveBeenSet_thenSetPX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_P00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PP0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingMX_setAllQsNormalToBC)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_M00] = 0.5;
    expectedQs[DIR_MP0] = 0.5;
    expectedQs[DIR_MM0] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYhasBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_M00] = 0.5;
    expectedQs[DIR_MM0] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYhasBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_M00] = 0.5;
    expectedQs[DIR_MP0] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPZhasBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_M00] = 0.5;
    expectedQs[DIR_MP0] = 0.5;
    expectedQs[DIR_MM0] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMZhasBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_M00] = 0.5;
    expectedQs[DIR_MP0] = 0.5;
    expectedQs[DIR_MM0] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYandMZhaveBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_M00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MM0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandPZhaveBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_M00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MM0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandPZhaveBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::PZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_M00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MP0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandMZhaveBeenSet_thenSetMX_doNotSetSameQsAgain)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::MZ);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_M00] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MP0] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingMZ_setAllQsNormalToBC)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00M] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_0PM] = 0.5;
    expectedQs[DIR_0MM] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYhasBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00M] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_0PM] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYhasBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00M] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_0MM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPXhasBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00M] = 0.5;
    expectedQs[DIR_M0M] = 0.5;
    expectedQs[DIR_0PM] = 0.5;
    expectedQs[DIR_0MM] = 0.5;
    expectedQs[DIR_MPM] = 0.5;
    expectedQs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMXhasBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00M] = 0.5;
    expectedQs[DIR_P0M] = 0.5;
    expectedQs[DIR_0PM] = 0.5;
    expectedQs[DIR_0MM] = 0.5;
    expectedQs[DIR_PPM] = 0.5;
    expectedQs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYandPXhaveBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0PM] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandMXhaveBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0PM] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PPM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandPXhaveBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0MM] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandMXhaveBeenSet_thenSetMZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0M] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0MM] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingPZ_setAllQsNormalToBC)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00P] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_0PP] = 0.5;
    expectedQs[DIR_0MP] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYhasBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00P] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_0PP] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPYhasBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00P] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_0MP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenPXhasBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00P] = 0.5;
    expectedQs[DIR_M0P] = 0.5;
    expectedQs[DIR_0PP] = 0.5;
    expectedQs[DIR_0MP] = 0.5;
    expectedQs[DIR_MPP] = 0.5;
    expectedQs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMXhasBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[DIR_00P] = 0.5;
    expectedQs[DIR_P0P] = 0.5;
    expectedQs[DIR_0PP] = 0.5;
    expectedQs[DIR_0MP] = 0.5;
    expectedQs[DIR_PPP] = 0.5;
    expectedQs[DIR_PMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQs));
}

TEST_F(SideTestBC, setQs3D_givenMYandPXhaveBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0PP] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenMYandMXhaveBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::MY);
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0PP] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandPXhaveBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::PX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_M0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0MP] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_MMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_givenPYandMXhaveBeenSet_thenSetPZ_doNotSetSameQsAgain)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;
    grid->addBCalreadySet(SideType::PY);
    grid->addBCalreadySet(SideType::MX);

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQsForTwoPreviousBCs(27, -1);
    expectedQsForTwoPreviousBCs[DIR_00P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_P0P] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_0MP] = 0.5;
    expectedQsForTwoPreviousBCs[DIR_PMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}
