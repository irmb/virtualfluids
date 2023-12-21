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
#include <array>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <memory>
#include <ostream>

#include <gpu/GridGenerator/grid/GridImp.h>

#include "PointerDefinitions.h"
#include "grid/Field.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"
#include "grid/distributions/Distribution.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "utilities/communication.h"

// This test is commented out because it causes a compiler error in Clang 10 --> The bug is fixed in Clang 14 (https://github.com/google/googletest/issues/2271)

// class FieldDouble : public Field
// {
// public:
//     FieldDouble() : Field(1)
//     {
//         this->allocateMemory();
//     };

//     void setToStopper(uint index)
//     {
//         this->field[index] = vf::gpu::STOPPER_SOLID;
//     }
// };

// class GridImpDouble : public GridImp
// {
// public:
//     std::array<real, 3> coordsOfTestedNode;
//     GridImpDouble(Object *object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
//                   Distribution d, uint level)
//         : GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, d, level)
//     {
//         this->neighborIndexX = new int[5];
//         this->neighborIndexY = new int[5];
//         this->neighborIndexZ = new int[5];
//     }

//     static SPtr<GridImpDouble> makeShared(Object *object, real startX, real startY, real startZ, real endX, real endY,
//                                           real endZ, real delta, Distribution d, uint level)
//     {
//         SPtr<GridImpDouble> grid(new GridImpDouble(object, startX, startY, startZ, endX, endY, endZ, delta, d, level));
//         return grid;
//     }

//     void transIndexToCoords(uint, real &x, real &y, real &z) const override
//     {
//         x = coordsOfTestedNode[0];
//         y = coordsOfTestedNode[1];
//         z = coordsOfTestedNode[2];
//     }

//     uint transCoordToIndex(const real &, const real &, const real &) const override
//     {
//         return 0;
//     }

//     void setStopperNeighborCoords(uint index) override
//     {
//         GridImp::setStopperNeighborCoords(index);
//     }

//     void setField(Field &field)
//     {
//         this->field = field;
//     }

//     MOCK_METHOD(int, getSparseIndex, (const real &x, const real &y, const real &z), (const, override));
// };

// // This is test is highly dependent on the implementation. Maybe it should be removed :(
// TEST(GridImp, setStopperNeighborCoords)
// {
//     real end = 1.0;
//     real delta = 0.1;

//     SPtr<GridImpDouble> gridImp =
//         GridImpDouble::makeShared(nullptr, 0.0, 0.0, 0.0, end, end, end, delta, Distribution(), 0);
//     FieldDouble field;
//     field.setToStopper(0);
//     gridImp->setField(field);

//     gridImp->coordsOfTestedNode = { end - ((real)0.5 * delta), end - ((real)0.5 * delta), end - ((real)0.5 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)0.51 * delta), end - ((real)0.51 * delta),
//                                     end - ((real)0.51 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);
//     gridImp->coordsOfTestedNode = { end - ((real)0.99 * delta), end - ((real)0.99 * delta),
//                                     end - ((real)0.99 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - delta, end - delta, end - delta };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)1.01 * delta), end - ((real)1.01 * delta),
//                                     end - ((real)1.01 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(3);
//     gridImp->setStopperNeighborCoords(0);

//     // The grid should not be like this, so this should be fine...
//     gridImp->coordsOfTestedNode = { end, end, end };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(0);
//     gridImp->setStopperNeighborCoords(0);

//     gridImp->coordsOfTestedNode = { end - ((real)0.25 * delta), end - ((real)0.25 * delta),
//                                     end - ((real)0.25 * delta) };
//     EXPECT_CALL(*gridImp, getSparseIndex).Times(0);
//     gridImp->setStopperNeighborCoords(0);
// }

std::array<int, 3> countInvalidNeighbors(SPtr<Grid> grid)
{
    auto countInvalidX = 0;
    auto countInvalidY = 0;
    auto countInvalidZ = 0;
    for (uint index = 0; index < grid->getSize(); index++) {
        if (grid->getNeighborsX()[index] == -1)
            countInvalidX++;
        if (grid->getNeighborsY()[index] == -1)
            countInvalidY++;
        if (grid->getNeighborsZ()[index] == -1)
            countInvalidZ++;
    }
    return { countInvalidX, countInvalidY, countInvalidZ };
}

std::array<int, 3> testFluidNodeNeighbors(SPtr<Grid> grid)
{
    auto countInvalidX = 0;
    auto countInvalidXY = 0;
    auto countInvalidXYZ = 0;
    for (uint index = 0; index < grid->getSize(); index++) {
        if (grid->getFieldEntry(index) != vf::gpu::FLUID) {
            continue;
        }

        auto neighX = grid->getNeighborsX()[index];
        if (neighX == -1) {
            countInvalidX++;
            continue;
        }

        auto neighXY = grid->getNeighborsY()[neighX];
        if (neighXY == -1) {
            countInvalidXY++;
            continue;
        }

        auto neighXYZ = grid->getNeighborsZ()[neighXY];
        if (neighXYZ == -1) {
            countInvalidXYZ++;
            continue;
        }
    }

    return { countInvalidX, countInvalidXY, countInvalidXYZ };
}

class findNeighborsIntegrationTest : public ::testing::Test
{
protected:
    SPtr<MultipleGridBuilder> gridBuilder;

    void SetUp() override
    {
        gridBuilder = std::make_shared<MultipleGridBuilder>();
    }
};

TEST_F(findNeighborsIntegrationTest, grid1)
{
    const real dx = 0.15;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    // Only the last layer of nodes should have invalid neighbors. The grid is a cube with a side length of 9 nodes
    // -> 9 * 9 = 81 invalid nodes are expected
    auto numberOfInvalidNeighbors = countInvalidNeighbors(grid);
    auto expected = 9 * 9;
    EXPECT_THAT(numberOfInvalidNeighbors[0], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[1], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[2], testing::Eq(expected));

    // additional test: all fluid nodes should have valid neighbors
    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, grid2)
{
    const real dx = 1.0 / 64;
    gridBuilder->addCoarseGrid(-0.6, -0.6, -0.6, 0.6, 0.6, 0.6, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    // Only the last layer of nodes should have invalid neighbors. The grid is a cube with a side length of 79 nodes
    // -> 79 * 79 invalid nodes are expected
    auto numberOfInvalidNeighbors = countInvalidNeighbors(grid);
    auto expected = 79 * 79;
    EXPECT_THAT(numberOfInvalidNeighbors[0], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[1], testing::Eq(expected));
    EXPECT_THAT(numberOfInvalidNeighbors[2], testing::Eq(expected));

    // additional test: all fluid nodes should have valid neighbors
    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, validFluidNeighbors1)
{
    real dx = 0.17;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}

TEST_F(findNeighborsIntegrationTest, validFluidNeighbors2)
{
    real dx = 0.18;
    gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, dx);

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);

    auto numberInvalidFluidNeighbors = testFluidNodeNeighbors(grid);
    EXPECT_THAT(numberInvalidFluidNeighbors[0], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[1], testing::Eq(0));
    EXPECT_THAT(numberInvalidFluidNeighbors[2], testing::Eq(0));
}


class PeriodicBoundaryShiftIntegrationTest : public testing::TestWithParam<std::tuple<int, int>>
{
protected:
    SPtr<MultipleGridBuilder> gridBuilder;
    const real dx{1.0};
    const int nx{5}, ny{std::get<0>(GetParam())}, nz{5};
    const int direction{std::get<1>(GetParam())};

    void SetUp() override
    {
        gridBuilder = std::make_shared<MultipleGridBuilder>();
        gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, nx*dx, ny*dx, nz*dx, dx);
    }
};

int getIndex(std::shared_ptr<Grid> grid, uint ix, uint iy, uint iz)
{
    return ix + grid->getNumberOfNodesX() * (iy + grid->getNumberOfNodesY() * iz);
}

void compareNodeToCoordinates(SPtr<Grid> grid, int index, real xExpected, real yExpected, real zExpected , const std::string& message)
{
    real xNode, yNode, zNode;
    grid->transIndexToCoords(index, xNode, yNode, zNode);
    EXPECT_FLOAT_EQ(xNode, xExpected)  << message;
    EXPECT_FLOAT_EQ(yNode, yExpected) << message;
    EXPECT_FLOAT_EQ(zNode, zExpected) << message;
}

void iterateOverXLine(std::shared_ptr<Grid> grid, uint iy, uint iz, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at X Line with iy: " + std::to_string(iy) + " and iz: " + std::to_string(iz);
    for(uint ix = 2; ix < grid->getNumberOfNodesX()-2; ix++){
        func(getIndex(grid, ix, iy, iz), msg);
    }
}

void iterateOverYLine(std::shared_ptr<Grid> grid, uint ix, uint iz, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at Y Line with ix: " + std::to_string(ix) + " and iz: " + std::to_string(iz);
    for(uint iy = 2; iy < grid->getNumberOfNodesY()-2; iy++){
        func(getIndex(grid, ix, iy, iz), msg);
    }
}

void iterateOverZLine(std::shared_ptr<Grid> grid, uint ix, uint iy, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at Z Line with ix: " + std::to_string(ix) + " and iy: " + std::to_string(iy);
    for(uint iz = 2; iz < grid->getNumberOfNodesZ()-2; iz++){
        func(getIndex(grid, ix, iy, iz), msg);
    }
}
void iterateOverYZPlane(std::shared_ptr<Grid> grid, uint ix, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at YZ Plane with ix: " + std::to_string(ix);
    for(uint iz = 2; iz < grid->getNumberOfNodesZ()-2; iz++){
    for(uint iy = 2; iy < grid->getNumberOfNodesY()-2; iy++){
        func(getIndex(grid, ix, iy, iz), msg);
    }}
}

void iterateOverXZPlane(std::shared_ptr<Grid> grid, uint iy, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at XZ Plane with iy: " + std::to_string(iy);
    for(uint iz = 2; iz < grid->getNumberOfNodesZ()-2; iz++){
    for(uint ix = 2; ix < grid->getNumberOfNodesX()-2; ix++){
        func(getIndex(grid, ix, iy, iz), msg);
    }}
}
void iterateOverXYPlane(std::shared_ptr<Grid> grid, uint iz, std::function<void(int, const std::string&)> func)
{    
    const std::string msg = "Failed at XY Plane with iz: " + std::to_string(iz);
    for(uint iy = 2; iy < grid->getNumberOfNodesY()-2; iy++){
    for(uint ix = 2; ix < grid->getNumberOfNodesX()-2; ix++){
        func(getIndex(grid, ix, iy, iz), msg);
    }}
}
void iterateOverInnerDomain(std::shared_ptr<Grid> grid, std::function<void(int, const std::string&)> func)
{
    const std::string msg = "Failed at Inner Domain";
    for(uint iz = 2; iz < grid->getNumberOfNodesZ()-2; iz++){
    for(uint iy = 2; iy < grid->getNumberOfNodesY()-2; iy++){
    for(uint ix = 2; ix < grid->getNumberOfNodesX()-2; ix++){
        func(getIndex(grid, ix, iy, iz), msg);
    }}}
}



TEST_P(PeriodicBoundaryShiftIntegrationTest, NoPeriodicity)
{
    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);
    auto func = [&](int index, const std::string msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    };
    iterateOverInnerDomain(grid, func);

    iterateOverYZPlane(grid, 1, func);
    iterateOverYZPlane(grid, nx, func);

    iterateOverXZPlane(grid, 1, func);
    iterateOverXZPlane(grid, ny, func);

    iterateOverXYPlane(grid, 1, func);
    iterateOverXYPlane(grid, nz, func);

    iterateOverXLine(grid, 1, 1, func);
    iterateOverXLine(grid, 1, nz, func);
    iterateOverXLine(grid, ny, 1, func);
    iterateOverXLine(grid, ny, nz, func);

    iterateOverYLine(grid, 1, 1, func);
    iterateOverYLine(grid, 1, nz, func);
    iterateOverYLine(grid, nx, 1, func);
    iterateOverYLine(grid, nx, nz, func);

    iterateOverZLine(grid, 1, 1, func);
    iterateOverZLine(grid, 1, ny, func);
    iterateOverZLine(grid, nx, 1, func);
    iterateOverZLine(grid, nx, ny, func);
}

TEST_P(PeriodicBoundaryShiftIntegrationTest, PeriodicNoShift)
{
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);
    iterateOverInnerDomain(grid, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });
    iterateOverYZPlane(grid, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], (nx-0.5)*dx, y-dx, z-dx, msg);
    });
    iterateOverYZPlane(grid, nx, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });
    iterateOverXLine(grid, 1, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, (ny-0.5)*dx, (nz-0.5)*dx, msg);
    });
    iterateOverXLine(grid, ny, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, 0.5*dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });
}

real wrap(real coord, real length){
    if(coord<0.0)
        return coord + length;
    if(coord>length)
        return coord - length;
    return coord;
}

real shiftCoordForward(real coord, real shift, real delta, int n, int direction, int caseNumber){
    if(caseNumber == direction)
        return wrap(coord + shift, n*delta);
    return coord;
}

real shiftCoordForward(real coord, real shift, real delta, int n, int direction, int caseNumber1, int caseNumber2)
{   
    if(direction==caseNumber1 || direction==caseNumber2)
        return shiftCoordForward(coord, shift, delta, n, 0, 0);
    return shiftCoordForward(coord, shift, delta, n, 0, 1);

}

real shiftCoordBackward(real coord, real shift, real delta, int n, int direction, int caseNumber){
    if(caseNumber == direction)
        return wrap(coord - (shift+delta), n*delta);
    return wrap(coord - delta, n*delta);
}

real shiftCoordBackward(real coord, real shift, real delta, int n, int direction, int caseNumber1, int caseNumber2){
    if(direction==caseNumber1 || direction==caseNumber2)
        return shiftCoordBackward(coord, shift, delta, n, 0, 0);
    return shiftCoordBackward(coord, shift, delta, n, 0, 1);
}


TEST_P(PeriodicBoundaryShiftIntegrationTest, PeriodicWithShift){
    const real shift = 3*dx; // must be integer multiple of dx because gridBuilder shortens shift to multiple of dx
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
    }

    gridBuilder->buildGrids(false);
    auto grid = gridBuilder->getGrid(0);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test inner domain
    iterateOverInnerDomain(grid, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test faces of the domain
    iterateOverYZPlane(grid, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 0);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 1);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], (nx-0.5)*dx, y_shifted, z_shifted, msg);
    });

    iterateOverXZPlane(grid, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 2);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted, (ny-0.5)*dx , z_shifted, msg);
    });

    iterateOverXYPlane(grid, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 4);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 5);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted, y_shifted, (nz-0.5)*dx, msg);
    });

    iterateOverYZPlane(grid, nx, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real y_shifted = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real z_shifted = shiftCoordForward(z, shift, dx, nz, direction, 1);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y_shifted, z_shifted, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });

    iterateOverXZPlane(grid, ny, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real z_shifted = shiftCoordForward(z, shift, dx, nz, direction, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x_shifted, 0.5*dx, z_shifted, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx , z-dx, msg);
    });

    iterateOverXYPlane(grid, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted = shiftCoordForward(y, shift, dx, ny, direction, 5);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x_shifted, y_shifted, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx , z-dx, msg);
    });


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test edges of x-z faces
    iterateOverYLine(grid, 1, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 4);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 0, 5);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 1);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted, y_shifted, z_shifted, msg);
    });

    iterateOverYLine(grid, 1, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real y_shifted_on_x = shiftCoordBackward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordBackward(z, shift, dx, nz, direction, 1);
        const real x_shifted_on_z = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordForward(y, shift, dx, ny, direction, 5);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x_shifted_on_z, y_shifted_on_z, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], (nx-0.5)*dx, y_shifted_on_x, z_shifted_on_x, msg);
    });

    iterateOverYLine(grid, nx, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real y_shifted_on_x = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordForward(z, shift, dx, nz, direction, 1);
        const real x_shifted_on_z = shiftCoordBackward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordBackward(y, shift, dx, ny, direction, 5);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y_shifted_on_x, z_shifted_on_x, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted_on_z, y_shifted_on_z, (nz-0.5)*dx, msg);
    });

    iterateOverYLine(grid, nx, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real y_shifted_on_x = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordForward(z, shift, dx, nz, direction, 1);
        const real x_shifted_on_z = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordForward(y, shift, dx, ny, direction, 5);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y_shifted_on_x, z_shifted_on_x, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x_shifted_on_z, y_shifted_on_z, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test edges of x-y faces
    iterateOverZLine(grid, 1, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 2);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 0);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 1, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted, y_shifted, z_shifted, msg);
    });

    iterateOverZLine(grid, 1, ny, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real y_shifted_on_x = shiftCoordBackward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordBackward(z, shift, dx, nz, direction, 1);
        const real z_shifted_on_y = shiftCoordForward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x_shifted_on_y, 0.5*dx, z_shifted_on_y, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], (nx-0.5)*dx, y_shifted_on_x, z_shifted_on_x, msg);
    });

    iterateOverZLine(grid, nx, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordBackward(x, shift, dx, nx, direction, 2);
        const real y_shifted_on_x = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordForward(z, shift, dx, nz, direction, 1);
        const real z_shifted_on_y = shiftCoordBackward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y_shifted_on_x, z_shifted_on_x, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted_on_y, (ny-0.5)*dx, z_shifted_on_y, msg);
    });

    iterateOverZLine(grid, nx, ny, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real y_shifted_on_x = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real z_shifted_on_x = shiftCoordForward(z, shift, dx, nz, direction, 1);
        const real z_shifted_on_y = shiftCoordForward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], 0.5*dx, y_shifted_on_x, z_shifted_on_x, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x_shifted_on_y, 0.5*dx, z_shifted_on_y, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test edges of y-z faces
    iterateOverXLine(grid, 1, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 2, 4);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 5);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted, y_shifted, z_shifted, msg);
    });

    iterateOverXLine(grid, 1, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordBackward(x, shift, dx, nx, direction, 2);
        const real x_shifted_on_z = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordForward(y, shift, dx, ny, direction, 5);
        const real z_shifted_on_y = shiftCoordBackward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x_shifted_on_z, y_shifted_on_z, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted_on_y, (ny-0.5)*dx, z_shifted_on_y, msg);
    });

    iterateOverXLine(grid, ny, 1, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real x_shifted_on_z = shiftCoordBackward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordBackward(y, shift, dx, ny, direction, 5);
        const real z_shifted_on_y = shiftCoordForward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x_shifted_on_y, 0.5*dx, z_shifted_on_y, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x_shifted_on_z, y_shifted_on_z, (nz-0.5)*dx, msg);
    });

    iterateOverXLine(grid, ny, nz, [&](int index, const std::string& msg){
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);
        const real x_shifted_on_y = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real x_shifted_on_z = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_z = shiftCoordForward(y, shift, dx, ny, direction, 5);
        const real z_shifted_on_y = shiftCoordForward(z, shift, dx, nz, direction, 3);

        compareNodeToCoordinates(grid, grid->getNeighborsX()[index], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[index], x_shifted_on_y, 0.5*dx, z_shifted_on_y, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[index], x_shifted_on_z, y_shifted_on_z, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[index], x-dx, y-dx, z-dx, msg);
    });

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!Test corners

    {
        const int ix{1}, iy{1}, iz{1};
        const std::string msg = "Failed at Corner 0";
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, iz), x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 2, 4);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 0, 5);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 1, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[getIndex(grid, ix, iy, iz)], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[getIndex(grid, ix, iy, iz)], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[getIndex(grid, ix, iy, iz)], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[getIndex(grid, ix, iy, iz)], x_shifted, y_shifted, z_shifted, msg);
    }

    {
        const int ix{nx}, iy{1}, iz{1};
        const std::string msg = "Failed at Corner 1";
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, iz), x, y, z);
        const real x_shifted = shiftCoordBackward(x, shift, dx, nx, direction, 2, 4);
        const real y_shifted_on_x = shiftCoordForward(y, shift, dx, ny, direction, 0);
        const real y_shifted_on_z = shiftCoordBackward(y, shift, dx, ny, direction, 5);
        const real z_shifted_on_x = shiftCoordForward(z, shift, dx, nz, direction, 1);
        const real z_shifted_on_y = shiftCoordBackward(z, shift, dx, nz, direction, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[getIndex(grid, ix, iy, iz)], 0.5*dx, y_shifted_on_x, z_shifted_on_x, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[getIndex(grid, ix, iy, iz)], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[getIndex(grid, ix, iy, iz)], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[getIndex(grid, ix, iy, iz)], x_shifted, y_shifted_on_z, z_shifted_on_y, msg);
    }
    {
        const int ix{1}, iy{ny}, iz{1};
        const std::string msg = "Failed at Corner 2";
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, iz), x, y, z);
        const real x_shifted_on_y = shiftCoordForward(x, shift, dx, nx, direction, 2);
        const real x_shifted_on_z = shiftCoordBackward(x, shift, dx, nx, direction, 4);
        const real y_shifted = shiftCoordBackward(y, shift, dx, ny, direction, 0, 5);
        const real z_shifted_on_x = shiftCoordBackward(z, shift, dx, nz, direction, 1);
        const real z_shifted_on_y = shiftCoordForward(z, shift, dx, nz, direction, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[getIndex(grid, ix, iy, iz)], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[getIndex(grid, ix, iy, iz)], x_shifted_on_y, 0.5*dx, z_shifted_on_y, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[getIndex(grid, ix, iy, iz)], x, y, z+dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[getIndex(grid, ix, iy, iz)], x_shifted_on_z, y_shifted, z_shifted_on_x, msg);
    }
    {
        const int ix{1}, iy{1}, iz{nz};
        const std::string msg = "Failed at Corner 3";
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, iz), x, y, z);
        const real x_shifted_on_y = shiftCoordBackward(x, shift, dx, nx, direction, 2);
        const real x_shifted_on_z = shiftCoordForward(x, shift, dx, nx, direction, 4);
        const real y_shifted_on_x = shiftCoordBackward(y, shift, dx, ny, direction, 0);
        const real y_shifted_on_z = shiftCoordForward(y, shift, dx, ny, direction, 5);
        const real z_shifted = shiftCoordBackward(z, shift, dx, nz, direction, 1, 3);
        compareNodeToCoordinates(grid, grid->getNeighborsX()[getIndex(grid, ix, iy, iz)], x+dx, y, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsY()[getIndex(grid, ix, iy, iz)], x, y+dx, z, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsZ()[getIndex(grid, ix, iy, iz)], x_shifted_on_z, y_shifted_on_z, 0.5*dx, msg);
        compareNodeToCoordinates(grid, grid->getNeighborsNegative()[getIndex(grid, ix, iy, iz)], x_shifted_on_y, y_shifted_on_x, z_shifted, msg);
    }
}

INSTANTIATE_TEST_SUITE_P(PeriodicBoundaryShiftSingleGPU, PeriodicBoundaryShiftIntegrationTest,
    testing::Combine(testing::Values(7, 100), testing::Range(0,6)));

class PeriodicBoundaryShiftMultiGPUIntegrationTest : public ::testing::TestWithParam<std::tuple<int, int>>
{
protected:
    SPtr<MultipleGridBuilder> gridBuilder;
    SPtr<BoundingBox> subdomain;
    const real dx{1.0};
    const int nx{5}, ny{std::get<0>(GetParam())}, nz{5};
    const int direction{std::get<1>(GetParam())};

    void SetUp() override
    {
        gridBuilder = std::make_shared<MultipleGridBuilder>();
        gridBuilder->addCoarseGrid(0.0, 0.0, 0.0, (nx+2)*dx, (ny+2)*dx, (nz+2)*dx, dx);
        gridBuilder->setPeriodicBoundaryCondition(true, true, true);
        subdomain = std::make_shared<BoundingBox>(dx, (nx+1)*dx, dx, (ny+1)*dx, dx,  (nz+1)*dx);
    }
};


TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, PX)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }

    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::PX, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iz=1; iz<grid->getNumberOfNodesZ()-1; iz++){
    for(uint iy=1; iy<grid->getNumberOfNodesY()-1; iy++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, nx+1, iy, iz), x, y, z);
        if(direction==0)
            y = wrap(y+shift, (ny+2)*dx);
        if(direction==1)
            z = wrap(z+shift, (nz+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::PX, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::PX, nodeIndex), x+dx, y, z, "Receive");
        nodeIndex++;
    }}; 
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::PX));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::PX));
}

TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, MX)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }
    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::MX, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iz=1; iz<grid->getNumberOfNodesZ()-1; iz++){
    for(uint iy=1; iy<grid->getNumberOfNodesY()-1; iy++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, 2, iy, iz), x, y, z);
        if(direction==0)
            y = wrap(y-shift, (ny+2)*dx);
        if(direction==1)
            z = wrap(z-shift, (nz+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::MX, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::MX, nodeIndex), x-dx, y, z, "Receive");
        nodeIndex++;
    }};
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::MX));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::MX));
}

TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, PY)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }

    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::PY, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iz=1; iz<grid->getNumberOfNodesZ()-1; iz++){
    for(uint ix=1; ix<grid->getNumberOfNodesX()-1; ix++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, ny+1, iz), x, y, z);
        if(direction==2)
            x = wrap(x+shift, (nx+2)*dx);
        if(direction==3)
            z = wrap(z+shift, (nz+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::PY, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::PY, nodeIndex), x, y+dx, z, "Receive");
        nodeIndex++;
    }}; 
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::PY));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::PY));
}

TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, MY)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }
    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::MY, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iz=1; iz<grid->getNumberOfNodesZ()-1; iz++){
    for(uint ix=1; ix<grid->getNumberOfNodesX()-1; ix++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, 2, iz), x, y, z);
        if(direction==2)
            x = wrap(x-shift, (nx+2)*dx);
        if(direction==3)
            z = wrap(z-shift, (nz+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::MY, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::MY, nodeIndex), x, y-dx, z, "Receive");
        nodeIndex++;
    }};
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::MY));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::MY));
}

TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, PZ)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);

    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }

    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::PZ, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iy=1; iy<grid->getNumberOfNodesY()-1; iy++){
    for(uint ix=1; ix<grid->getNumberOfNodesX()-1; ix++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, nx+1), x, y, z);
        if(direction==4)
            x = wrap(x+shift, (nx+2)*dx);
        if(direction==5)
            y = wrap(y+shift, (ny+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::PZ, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::PZ, nodeIndex), x, y, z+dx, "Receive");
        nodeIndex++;
    }}; 
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::PZ));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::PZ));
}

TEST_P(PeriodicBoundaryShiftMultiGPUIntegrationTest, MZ)
{
    const real shift = 3.0f*dx;
    gridBuilder->setPeriodicBoundaryCondition(true, true, true);
    switch(direction){
        case 0: gridBuilder->setPeriodicShiftOnXBoundaryInYDirection(shift); break;
        case 1: gridBuilder->setPeriodicShiftOnXBoundaryInZDirection(shift); break;
        case 2: gridBuilder->setPeriodicShiftOnYBoundaryInXDirection(shift); break;
        case 3: gridBuilder->setPeriodicShiftOnYBoundaryInZDirection(shift); break;
        case 4: gridBuilder->setPeriodicShiftOnZBoundaryInXDirection(shift); break;
        case 5: gridBuilder->setPeriodicShiftOnZBoundaryInYDirection(shift); break;
        case 6: break;
    }
    gridBuilder->buildGrids();
    auto grid = gridBuilder->getGrid(0);
    grid->findCommunicationIndices(CommunicationDirections::MZ, subdomain, direction!=6);
    uint nodeIndex=0;
    for(uint iy=1; iy<grid->getNumberOfNodesY()-1; iy++){
    for(uint ix=1; ix<grid->getNumberOfNodesX()-1; ix++){
        real x, y, z;
        grid->transIndexToCoords(getIndex(grid, ix, iy, 2), x, y, z);
        if(direction==4)
            x = wrap(x-shift, (nx+2)*dx);
        if(direction==5)
            y = wrap(y-shift, (ny+2)*dx);
        compareNodeToCoordinates(grid, grid->getSendIndex(CommunicationDirections::MZ, nodeIndex), x, y, z, "Send");
        compareNodeToCoordinates(grid, grid->getReceiveIndex(CommunicationDirections::MZ, nodeIndex), x, y, z-dx, "Receive");
        nodeIndex++;
    }};
    EXPECT_EQ(nodeIndex, grid->getNumberOfSendNodes(CommunicationDirections::MZ));
    EXPECT_EQ(nodeIndex, grid->getNumberOfReceiveNodes(CommunicationDirections::MZ));
}
INSTANTIATE_TEST_SUITE_P(PeriodicBoundaryShiftMultiGPUIntegration, PeriodicBoundaryShiftMultiGPUIntegrationTest, testing::Combine(testing::Values(5), testing::Range(0,7)));

//! \}
