#include "Side.h"
#include "PointerDefinitions.h"

#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/GridImp.h"
#include "grid/NodeValues.h"
#include "lbm/constants/D3Q27.h"
#include "gmock/gmock.h"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <ostream>
#include <vector>

using namespace vf::gpu;
using namespace vf::lbm::dir;

class SideForTest : public Side
{

public:
    void setQs(SPtr<Grid> grid, SPtr<gg::BoundaryCondition> boundaryCondition, uint index, bool nodeIsDifferentBC)
    {
        Side::setQs(grid, boundaryCondition, index, nodeIsDifferentBC);
    };

private:
    void correctNeighborForPeriodicBoundaries(Grid *grid, real x, real y, real z, real *coords, real neighborX,
                                              real neighborY, real neighborZ) const override
    {
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }

    // all directions with a positive x component are aligned with the normal (! depends on order in D3Q27.h !)
    
    std::vector<bool> alignedWithPX = {false, true, false, false, false, false, false, true, false, true, false, true, false, true, false, false, false, false, false, true, false, true, false, true, false, true, false};
    bool isAlignedWithNormal(Grid * /*grid*/, int dir) const override
    {
        return alignedWithPX[dir];
    }

    void addIndices(std::vector<SPtr<Grid>> grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override
    {
    }

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    SideType whoAmI() const override
    {
        return SideType::PX;
    }
};

class GridDouble : public GridImp
{
private:
    int dir[1] = { -1 };

    int *directions = dir;

public:
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
        return 10;
    }

    int *getDirection() const override
    {
        return directions;
    };
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
};

TEST(SideTest, setQs)
{
    SideForTest side;
    SPtr<GridDouble> grid = std::make_shared<GridDouble>();
    SPtr<BoundaryConditionSpy> bc = std::make_shared<BoundaryConditionSpy>();
    uint index = 0;

    // first bc for this node

    bool nodeIsDifferentBC = false;
    side.setQs(grid, bc, index, nodeIsDifferentBC);
    auto qsFirstBC = bc->getQs()[0];

    std::vector<real> expectedFirstBC = { -1, 0.5, -1, -1, -1, -1, -1, 0.5, -1, 0.5, -1 };
    EXPECT_THAT(qsFirstBC, testing::Eq(expectedFirstBC));

    // node already has different BC

    nodeIsDifferentBC = true;
    side.setQs(grid, bc, index, nodeIsDifferentBC);
    auto qsPreviousQ = bc->getQs()[0];

    std::vector<real> expectedPreviousQ = { -1, 0.5, -1, -1, -1, -1, -1, -1, -1, 0.5, -1 };
    EXPECT_THAT(qsPreviousQ, testing::Eq(expectedPreviousQ));

    // the 7th entry (DIR_PP0) is different
    EXPECT_TRUE(qsFirstBC[7] != qsPreviousQ[7]);
}
