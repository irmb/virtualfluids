//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_grid_tests grid
//! \ingroup gpu_GridGenerator_tests GridGenerator
//! \{
//=======================================================================================
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
    void correctNeighborForPeriodicBoundaries(const Grid *grid, std::array<real, 3>& coords, std::array<real, 3>& neighbors) const override
    {
    }

    int getDirection() const override
    {
        return sideDirection;
    }

    void addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<gg::BoundaryCondition> boundaryCondition) override
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dPM0] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPM0] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dPM0] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dPMP] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dPMM] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dPM0] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dPMP] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dPMM] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPM0] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dPMP] = 0.5;
    expectedQs[dPMM] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dPPM] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dPM0] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dPMM] = 0.5;
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
    expectedQs[dP00] = 0.5;
    expectedQs[dPP0] = 0.5;
    expectedQs[dPM0] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dPMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[dP00] = 0.5;
    expectedQsForTwoPreviousBCs[dPM0] = 0.5;
    expectedQsForTwoPreviousBCs[dP0P] = 0.5;
    expectedQsForTwoPreviousBCs[dPMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[dP00] = 0.5;
    expectedQsForTwoPreviousBCs[dPM0] = 0.5;
    expectedQsForTwoPreviousBCs[dP0M] = 0.5;
    expectedQsForTwoPreviousBCs[dPMM] = 0.5;
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
    expectedQsForTwoPreviousBCs[dP00] = 0.5;
    expectedQsForTwoPreviousBCs[dPP0] = 0.5;
    expectedQsForTwoPreviousBCs[dP0M] = 0.5;
    expectedQsForTwoPreviousBCs[dPPM] = 0.5;
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
    expectedQsForTwoPreviousBCs[dP00] = 0.5;
    expectedQsForTwoPreviousBCs[dPP0] = 0.5;
    expectedQsForTwoPreviousBCs[dP0P] = 0.5;
    expectedQsForTwoPreviousBCs[dPPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingMX_setAllQsNormalToBC)
{
    side.coordinateDirection = X_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[dM00] = 0.5;
    expectedQs[dMP0] = 0.5;
    expectedQs[dMM0] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[dMPP] = 0.5;
    expectedQs[dMMP] = 0.5;
    expectedQs[dMPM] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[dM00] = 0.5;
    expectedQs[dMM0] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[dMMP] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[dM00] = 0.5;
    expectedQs[dMP0] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[dMPP] = 0.5;
    expectedQs[dMPM] = 0.5;
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
    expectedQs[dM00] = 0.5;
    expectedQs[dMP0] = 0.5;
    expectedQs[dMM0] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[dMPM] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[dM00] = 0.5;
    expectedQs[dMP0] = 0.5;
    expectedQs[dMM0] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[dMPP] = 0.5;
    expectedQs[dMMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[dM00] = 0.5;
    expectedQsForTwoPreviousBCs[dMM0] = 0.5;
    expectedQsForTwoPreviousBCs[dM0P] = 0.5;
    expectedQsForTwoPreviousBCs[dMMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[dM00] = 0.5;
    expectedQsForTwoPreviousBCs[dMM0] = 0.5;
    expectedQsForTwoPreviousBCs[dM0M] = 0.5;
    expectedQsForTwoPreviousBCs[dMMM] = 0.5;
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
    expectedQsForTwoPreviousBCs[dM00] = 0.5;
    expectedQsForTwoPreviousBCs[dMP0] = 0.5;
    expectedQsForTwoPreviousBCs[dM0M] = 0.5;
    expectedQsForTwoPreviousBCs[dMPM] = 0.5;
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
    expectedQsForTwoPreviousBCs[dM00] = 0.5;
    expectedQsForTwoPreviousBCs[dMP0] = 0.5;
    expectedQsForTwoPreviousBCs[dM0P] = 0.5;
    expectedQsForTwoPreviousBCs[dMPP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingMZ_setAllQsNormalToBC)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = NEGATIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[d00M] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[d0PM] = 0.5;
    expectedQs[d0MM] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dMPM] = 0.5;
    expectedQs[dPMM] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[d00M] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[d0PM] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dMPM] = 0.5;
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
    expectedQs[d00M] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[d0MM] = 0.5;
    expectedQs[dPMM] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[d00M] = 0.5;
    expectedQs[dM0M] = 0.5;
    expectedQs[d0PM] = 0.5;
    expectedQs[d0MM] = 0.5;
    expectedQs[dMPM] = 0.5;
    expectedQs[dMMM] = 0.5;
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
    expectedQs[d00M] = 0.5;
    expectedQs[dP0M] = 0.5;
    expectedQs[d0PM] = 0.5;
    expectedQs[d0MM] = 0.5;
    expectedQs[dPPM] = 0.5;
    expectedQs[dPMM] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00M] = 0.5;
    expectedQsForTwoPreviousBCs[dM0M] = 0.5;
    expectedQsForTwoPreviousBCs[d0PM] = 0.5;
    expectedQsForTwoPreviousBCs[dMPM] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00M] = 0.5;
    expectedQsForTwoPreviousBCs[dP0M] = 0.5;
    expectedQsForTwoPreviousBCs[d0PM] = 0.5;
    expectedQsForTwoPreviousBCs[dPPM] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00M] = 0.5;
    expectedQsForTwoPreviousBCs[dM0M] = 0.5;
    expectedQsForTwoPreviousBCs[d0MM] = 0.5;
    expectedQsForTwoPreviousBCs[dMMM] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00M] = 0.5;
    expectedQsForTwoPreviousBCs[dP0M] = 0.5;
    expectedQsForTwoPreviousBCs[d0MM] = 0.5;
    expectedQsForTwoPreviousBCs[dPMM] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

TEST_F(SideTestBC, setQs3D_whenSettingPZ_setAllQsNormalToBC)
{
    side.coordinateDirection = Z_INDEX;
    side.sideDirection = POSITIVE_DIR;

    side.setQs(grid, bc, index);
    auto actualQs = bc->getQs()[0];

    std::vector<real> expectedQs(27, -1);
    expectedQs[d00P] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[d0PP] = 0.5;
    expectedQs[d0MP] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dMPP] = 0.5;
    expectedQs[dPMP] = 0.5;
    expectedQs[dMMP] = 0.5;
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
    expectedQs[d00P] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[d0PP] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dMPP] = 0.5;
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
    expectedQs[d00P] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[d0MP] = 0.5;
    expectedQs[dPMP] = 0.5;
    expectedQs[dMMP] = 0.5;
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
    expectedQs[d00P] = 0.5;
    expectedQs[dM0P] = 0.5;
    expectedQs[d0PP] = 0.5;
    expectedQs[d0MP] = 0.5;
    expectedQs[dMPP] = 0.5;
    expectedQs[dMMP] = 0.5;
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
    expectedQs[d00P] = 0.5;
    expectedQs[dP0P] = 0.5;
    expectedQs[d0PP] = 0.5;
    expectedQs[d0MP] = 0.5;
    expectedQs[dPPP] = 0.5;
    expectedQs[dPMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00P] = 0.5;
    expectedQsForTwoPreviousBCs[dM0P] = 0.5;
    expectedQsForTwoPreviousBCs[d0PP] = 0.5;
    expectedQsForTwoPreviousBCs[dMPP] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00P] = 0.5;
    expectedQsForTwoPreviousBCs[dP0P] = 0.5;
    expectedQsForTwoPreviousBCs[d0PP] = 0.5;
    expectedQsForTwoPreviousBCs[dPPP] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00P] = 0.5;
    expectedQsForTwoPreviousBCs[dM0P] = 0.5;
    expectedQsForTwoPreviousBCs[d0MP] = 0.5;
    expectedQsForTwoPreviousBCs[dMMP] = 0.5;
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
    expectedQsForTwoPreviousBCs[d00P] = 0.5;
    expectedQsForTwoPreviousBCs[dP0P] = 0.5;
    expectedQsForTwoPreviousBCs[d0MP] = 0.5;
    expectedQsForTwoPreviousBCs[dPMP] = 0.5;
    EXPECT_THAT(actualQs, testing::Eq(expectedQsForTwoPreviousBCs));
}

//! \}
