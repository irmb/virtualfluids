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
#include <gmock/gmock.h>
#include <memory>
#include <stdexcept>

#include "MultipleGridBuilderFacade.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Sphere/Sphere.h"
#include "core/BoundaryConditions/BoundaryConditionKernelManager.h"
#include "grid/BoundaryConditions/Side.h"
#include "basics/constants/NumericConstants.h"
#include "utilities/communication.h"
#include "grid/GridBuilder/MultipleGridBuilder.h"

class MultipleGridBuilderFacadeTest : public MultipleGridBuilderFacade
{public:
    MultipleGridBuilderFacadeTest(SPtr<MultipleGridBuilder> gridBuilder, UPtr<GridDimensions> gridDimensions, std::optional<real> overlapOfSubdomains =  std::nullopt)
        : MultipleGridBuilderFacade(std::move(gridBuilder), std::move(gridDimensions), overlapOfSubdomains)
    {
    }
    MultipleGridBuilderFacadeTest(SPtr<MultipleGridBuilder> gridBuilder, real minX, real maxX, real minY, real maxY,
                                  real minZ, real maxZ, real dxGrid, std::optional<real> overlapOfSubdomains = std::nullopt)
        : MultipleGridBuilderFacade(std::move(gridBuilder),
                                    std::make_unique<GridDimensions>(minX, maxX, minY, maxY, minZ, maxZ, dxGrid), overlapOfSubdomains)
    {
    }

    void setNumberOfGrids(uint x, uint y, uint z)
    {
        numberOfSubdomains = { x, y, z };
    };

    uint getX3D(uint index1D)
    {
        return MultipleGridBuilderFacade::getX3D(index1D);
    }
    uint getY3D(uint index1D)
    {
        return MultipleGridBuilderFacade::getY3D(index1D);
    }
    uint getZ3D(uint index1D)
    {
        return MultipleGridBuilderFacade::getZ3D(index1D);
    }
    uint getIndex1D(uint xIndex, uint yIndex, uint zIndex)
    {
        return MultipleGridBuilderFacade::getIndex1D(xIndex, yIndex, zIndex);
    }
    uint getIndex1D(const std::array<uint, 3>& index3D)
    {
        return MultipleGridBuilderFacade::getIndex1D(index3D);
    }
};

class MockMultipleGridBuilder : public MultipleGridBuilder
{
public:
    MOCK_METHOD(void, addCoarseGrid, (real startX, real startY, real startZ, real endX, real endY, real endZ, real delta),
                (override));
    void addCoarseGrid(const GridDimensions& gridDimensions) override {};
    MOCK_METHOD(void, setSubDomainBox, (SPtr<BoundingBox> subDomainBox), (override));
    MOCK_METHOD(void, findCommunicationIndices, (int direction, bool doShift), (override));
    MOCK_METHOD(void, setCommunicationProcess, (int direction, uint process), (override));

    MOCK_METHOD(void, setNoSlipBoundaryCondition, (SideType side), (override));
    MOCK_METHOD(void, setStressBoundaryCondition,
                (SideType sideType, real normalX, real normalY, real normalZ, uint samplingOffset, real z0, real dx),
                (override));
    MOCK_METHOD(void, setVelocityBoundaryCondition, (SideType sideType, real vx, real vy, real vz), (override));
    MOCK_METHOD(void, setPressureBoundaryCondition, (SideType sideType, real rho), (override));
    MOCK_METHOD(void, setSlipBoundaryCondition, (SideType sideType, real normalX, real normalY, real normalZ),
                (override));
    MOCK_METHOD(void, setPrecursorBoundaryCondition,
                (SideType sideType, SPtr<FileCollection> fileCollection, int timeStepsBetweenReads, real velocityX,
                 real velocityY, real velocityZ, std::vector<uint> fileLevelToGridLevelMap),
                (override));
    MOCK_METHOD(void, setPeriodicBoundaryCondition, (bool periodic_X, bool periodic_Y, bool periodic_Z), (override));
    MOCK_METHOD(void, setNumberOfLayers, (uint numberOfLayersFine, uint numberOfLayersBetweenLevels), (override));
    MOCK_METHOD(void, addGrid, (std::shared_ptr<Object> gridShape, uint levelFine), (override));
    MOCK_METHOD(void, addGrid, (std::shared_ptr<Object> gridShape), (override));
    MOCK_METHOD(void, addGeometry, (std::shared_ptr<Object> gridShape), (override));

    void buildGrids(bool enableThinWalls) override {};
};

TEST(MultipleGridBuilderFacadeTest, transform1dCoordinateToComponents)
{
    //            y ^
    //              |
    //              |---->
    //             /      x
    //         z  /

    //    zIndex = 0:    -----------
    //                   |  4 |  5 |
    //                   -----------
    //                   |  2 |  3 |
    //                   -----------
    //                   |  0 |  1 |
    //                   -----------

    //    zIndex = 3:    -----------
    //                   | 22 | 23 |
    //                   -----------
    //                   | 20 | 21 |
    //                   -----------
    //                   | 18 | 19 |
    //                   -----------

    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    sut.setNumberOfGrids(2, 3, 4);

    // x
    EXPECT_THAT(sut.getX3D(0), testing::Eq(0));
    EXPECT_THAT(sut.getX3D(11), testing::Eq(1));
    EXPECT_THAT(sut.getX3D(12), testing::Eq(0));
    EXPECT_THAT(sut.getX3D(23), testing::Eq(1));

    // y
    EXPECT_THAT(sut.getY3D(0), testing::Eq(0));
    EXPECT_THAT(sut.getY3D(9), testing::Eq(1));
    EXPECT_THAT(sut.getY3D(10), testing::Eq(2));
    EXPECT_THAT(sut.getY3D(13), testing::Eq(0));
    EXPECT_THAT(sut.getY3D(23), testing::Eq(2));

    // z
    EXPECT_THAT(sut.getZ3D(0), testing::Eq(0));
    EXPECT_THAT(sut.getZ3D(6), testing::Eq(1));
    EXPECT_THAT(sut.getZ3D(17), testing::Eq(2));
    EXPECT_THAT(sut.getZ3D(23), testing::Eq(3));

    sut.setNumberOfGrids(1, 1, 1);
    EXPECT_THAT(sut.getX3D(0), testing::Eq(0));
    EXPECT_THAT(sut.getY3D(0), testing::Eq(0));
    EXPECT_THAT(sut.getZ3D(0), testing::Eq(0));
}

TEST(MultipleGridBuilderFacadeTest, transformComponentsTo1DCoordinate)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    sut.setNumberOfGrids(2, 3, 4);

    EXPECT_THAT(sut.getIndex1D(0, 0, 0), testing::Eq(0));
    EXPECT_THAT(sut.getIndex1D({0, 0, 0}), testing::Eq(0));
    EXPECT_THAT(sut.getIndex1D(1, 2, 1), testing::Eq(11));
    EXPECT_THAT(sut.getIndex1D({1, 2, 1}), testing::Eq(11));
    EXPECT_THAT(sut.getIndex1D(0, 0, 2), testing::Eq(12));
    EXPECT_THAT(sut.getIndex1D({0, 0, 2}), testing::Eq(12));
    EXPECT_THAT(sut.getIndex1D(1, 2, 3), testing::Eq(23));
    EXPECT_THAT(sut.getIndex1D({1, 2, 3}), testing::Eq(23));
}

TEST(MultipleGridBuilderFacadeTest, noOverlapOnMultiGpu_Throws)
{
    auto gridBuilder = std::make_shared<MultipleGridBuilder>();
    auto sut = MultipleGridBuilderFacadeTest(gridBuilder, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    sut.addDomainSplit(1.0, Axis::x);
    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, createGrids1Gpu)
{
    auto mockGridBuilder = std::make_shared<MockMultipleGridBuilder>();
    auto sut = MultipleGridBuilderFacadeTest(mockGridBuilder, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.1);

    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(-1.0, -2.0, -3.0, 1.0, 2.0, 3.0, 0.1)).Times(1);
    EXPECT_CALL(*mockGridBuilder, setSubDomainBox).Times(0);
    sut.createGrids(0);
}

TEST(MultipleGridBuilderFacadeTest, createGrid_ProcessIdIsLargerThanNumberOfSubdomains_Throws)
{
    auto mockGridBuilder = std::make_shared<MockMultipleGridBuilder>();
    auto sut = MultipleGridBuilderFacadeTest(mockGridBuilder, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.1);

    EXPECT_THROW(sut.createGrids(1), std::runtime_error); // generatePart >= number of subdomains --> should throw
}

class MultipleGridBuilderFacadeTest_24subdomains : public testing::Test
{
protected:
    //            y ^
    //              |
    //              |---->
    //             /      x
    //         z  /

    //    zIndex = 0:    -----------
    //                   |  4 |  5 |
    //                   -----------
    //                   |  2 |  3 |
    //                   -----------
    //                   |  0 |  1 |
    //                   -----------

    //    zIndex = 3:    -----------
    //                   | 22 | 23 |
    //                   -----------
    //                   | 20 | 21 |
    //                   -----------
    //                   | 18 | 19 |
    //                   -----------

    std::shared_ptr<testing::NiceMock<MockMultipleGridBuilder>> mockGridBuilder =
        std::make_shared<testing::NiceMock<MockMultipleGridBuilder>>();
    std::unique_ptr<MultipleGridBuilderFacadeTest> sut;

    void SetUp() override
    {
        createNewSut();
    }

    void createNewSut()
    {
        sut = std::make_unique<MultipleGridBuilderFacadeTest>(mockGridBuilder, std::make_unique<GridDimensions>(0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1), 0.1);
        sut->addDomainSplit(1.0, Axis::x);
        sut->addDomainSplit(1.0, Axis::y);
        sut->addDomainSplit(2.0, Axis::y);
        sut->addDomainSplit(1.0, Axis::z);
        sut->addDomainSplit(2.0, Axis::z);
        sut->addDomainSplit(3.0, Axis::z);
    }
};

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, createGridsMultiGpu)
{
    // addCoarseGrid works with overlaps
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(0.0, 0.0, 0.0, 1.1, 1.1, 1.1, 0.1));
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(0.9, 0.9, 0.9, 2.0, 2.1, 2.1, 0.1));
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(0.9, 1.9, 0.9, 2.0, 3.0, 2.1, 0.1));
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(0.0, 0.0, 1.9, 1.1, 1.1, 3.1, 0.1));
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid(0.9, 1.9, 2.9, 2.0, 3.0, 4.0, 0.1));

    // setSubDomainBox has to be called once per subdomain
    EXPECT_CALL(*mockGridBuilder, setSubDomainBox(testing::_)).Times(5);

    sut->createGrids(0);
    this->createNewSut();
    sut->createGrids(9);
    this->createNewSut();
    sut->createGrids(11);
    this->createNewSut();
    sut->createGrids(12);
    this->createNewSut();
    sut->createGrids(23);
}

TEST(MultipleGridBuilderFacadeTest, xSplitToLarge)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(10.0, Axis::x); // xSplit > maxX

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, xSplitToSmall)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(-1.0, Axis::x); // xSplit < minX

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, ySplitToLarge)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::y);  // valid ySplit
    sut.addDomainSplit(10.0, Axis::y); // ySplit > maxY

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, ySplitToSmall)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::y);
    sut.addDomainSplit(-1.0, Axis::y); // ySplit < minY

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, zSplitToLarge)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::z);
    sut.addDomainSplit(10.0, Axis::z); // zSplit > maxZ
    sut.addDomainSplit(2.0, Axis::z);

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, zSplitToSmall)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::z);
    sut.addDomainSplit(-1.0, Axis::z); // zSplit < minZ
    sut.addDomainSplit(2.0, Axis::z);

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, sameSplitTwiceY)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::y);
    sut.addDomainSplit(2.0, Axis::y);
    sut.addDomainSplit(1.0, Axis::y);

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, sameSplitTwiceZ)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1);
    sut.addDomainSplit(0.9, Axis::z);
    sut.addDomainSplit(0.9, Axis::z);

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, sameSplitTwiceX)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1, 0.1);
    sut.addDomainSplit(1.0, Axis::x);
    sut.addDomainSplit(1.0, Axis::x);

    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
    EXPECT_THROW(sut.createGrids(1), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, setUpCommunicationNeighbors1Gpu)
{
    // use nice mock to suppress uninteresting call warning in createGrids()
    auto niceMockGridBuilder = std::make_shared<testing::NiceMock<MockMultipleGridBuilder>>();
    auto sut = MultipleGridBuilderFacadeTest(niceMockGridBuilder, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.1);

    EXPECT_CALL(*niceMockGridBuilder, findCommunicationIndices).Times(0);
    EXPECT_CALL(*niceMockGridBuilder, setCommunicationProcess).Times(0);

    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, setUpCommunicationNeighborsMultiGpu)
{
    // process index 0
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(1, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(3, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(5, false));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(1, 1));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(3, 2));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(5, 6));
    sut->createGrids(0);

    // process index 9
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(0, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(2, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(3, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(4, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(5, false));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(0, 8));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(2, 7));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(3, 11));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(4, 3));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(5, 15));
    sut->createGrids(9);

    // process index 23
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(0, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(2, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(4, false));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(0, 22));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(2, 21));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(4, 17));
    sut->createGrids(23);
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, noSlipPX)
{
    // process index 0
    EXPECT_CALL(*mockGridBuilder, setNoSlipBoundaryCondition).Times(0);
    sut->createGrids(0);
    sut->setNoSlipBoundaryCondition(SideType::PX);

    // process index 9
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setNoSlipBoundaryCondition(SideType::PX));
    sut->createGrids(9);
    sut->setNoSlipBoundaryCondition(SideType::PX);

    // process index 18
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setNoSlipBoundaryCondition).Times(0);
    sut->createGrids(18);
    sut->setNoSlipBoundaryCondition(SideType::PX);

    // process index 23
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setNoSlipBoundaryCondition(SideType::PX));
    sut->createGrids(23);
    sut->setNoSlipBoundaryCondition(SideType::PX);

    EXPECT_CALL(*mockGridBuilder, setNoSlipBoundaryCondition(SideType::PX)).Times(12);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setNoSlipBoundaryCondition(SideType::PX);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, stressMX)
{
    SideType sideType = SideType::MX;
    real normalX = 0.1;
    real normalY = 0.2;
    real normalZ = 0.3;
    uint samplingOffset = 1;
    real z0 = 0.5;
    real dx = 0.7;

    // process index
    EXPECT_CALL(*mockGridBuilder,
                setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx));
    sut->createGrids(0);
    sut->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);

    // process index 9
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setStressBoundaryCondition).Times(0);
    sut->createGrids(9);
    sut->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);

    // process index 18
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder,
                setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx));
    sut->createGrids(18);
    sut->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);

    // process index 23
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setStressBoundaryCondition).Times(0);
    sut->createGrids(23);
    sut->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);

    EXPECT_CALL(*mockGridBuilder,
                setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx))
        .Times(12);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setStressBoundaryCondition(sideType, normalX, normalY, normalZ, samplingOffset, z0, dx);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, velocityPY)
{
    SideType sideType = SideType::PY;
    real vx = 0.1;
    real vy = 0.2;
    real vz = 0.3;

    // process index 0
    EXPECT_CALL(*mockGridBuilder, setVelocityBoundaryCondition).Times(0);
    sut->createGrids(0);
    sut->setVelocityBoundaryCondition(sideType, vx, vy, vz);

    // process index 9
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setVelocityBoundaryCondition).Times(0);
    sut->createGrids(9);
    sut->setVelocityBoundaryCondition(sideType, vx, vy, vz);

    // process index 16
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setVelocityBoundaryCondition(sideType, vx, vy, vz));
    sut->createGrids(16);
    sut->setVelocityBoundaryCondition(sideType, vx, vy, vz);

    // process index 23
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setVelocityBoundaryCondition(sideType, vx, vy, vz));
    sut->createGrids(23);
    sut->setVelocityBoundaryCondition(sideType, vx, vy, vz);

    EXPECT_CALL(*mockGridBuilder, setVelocityBoundaryCondition(sideType, vx, vy, vz)).Times(8);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setVelocityBoundaryCondition(sideType, vx, vy, vz);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, pressureMY)
{
    SideType sideType = SideType::MY;
    real rho = 0.1;

    // process index 0
    EXPECT_CALL(*mockGridBuilder, setPressureBoundaryCondition(sideType, rho));
    sut->createGrids(0);
    sut->setPressureBoundaryCondition(sideType, rho);

    // process index 9
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setPressureBoundaryCondition).Times(0);
    sut->createGrids(9);
    sut->setPressureBoundaryCondition(sideType, rho);

    // process index 16
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setPressureBoundaryCondition).Times(0);
    sut->createGrids(16);
    sut->setPressureBoundaryCondition(sideType, rho);

    // process index 23
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setPressureBoundaryCondition).Times(0);
    sut->createGrids(23);
    sut->setPressureBoundaryCondition(sideType, rho);

    EXPECT_CALL(*mockGridBuilder, setPressureBoundaryCondition).Times(8);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setPressureBoundaryCondition(sideType, rho);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, slipPZ)
{
    SideType sideType = SideType::PZ;
    real normalX = 0.1;
    real normalY = 0.2;
    real normalZ = 0.3;

    // process index 0
    EXPECT_CALL(*mockGridBuilder, setSlipBoundaryCondition).Times(0);
    this->createNewSut();
    sut->createGrids(0);
    sut->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ);

    // process index 17
    EXPECT_CALL(*mockGridBuilder, setSlipBoundaryCondition).Times(0);
    this->createNewSut();
    sut->createGrids(17);
    sut->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ);

    // process index 18
    EXPECT_CALL(*mockGridBuilder, setSlipBoundaryCondition(sideType, normalX, normalY, normalZ));
    this->createNewSut();
    sut->createGrids(18);
    sut->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ);

    // process index 23
    EXPECT_CALL(*mockGridBuilder, setSlipBoundaryCondition(sideType, normalX, normalY, normalZ));
    this->createNewSut();
    sut->createGrids(23);
    sut->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ);

    EXPECT_CALL(*mockGridBuilder, setSlipBoundaryCondition).Times(6);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setSlipBoundaryCondition(sideType, normalX, normalY, normalZ);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, precursorMZ)
{
    SideType sideType = SideType::MZ;
    SPtr<FileCollection> fileCollection = nullptr;
    int timeStepsBetweenReads = 10;
    real velocityX = c1o2;
    real velocityY = c1o3;
    real velocityZ = c1o4;
    std::vector<uint> fileLevelToGridLevelMap = { 0 };

    EXPECT_CALL(*mockGridBuilder,
                setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY,
                                              velocityZ, fileLevelToGridLevelMap));
    sut->createGrids(0);
    sut->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY, velocityZ,
                                       fileLevelToGridLevelMap);

    // process index 17
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setPrecursorBoundaryCondition).Times(0);
    sut->createGrids(17);
    sut->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY, velocityZ,
                                       fileLevelToGridLevelMap);

    // process index 18
    this->createNewSut();
    EXPECT_CALL(*mockGridBuilder, setPrecursorBoundaryCondition).Times(0);
    sut->createGrids(18);
    sut->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY, velocityZ,
                                       fileLevelToGridLevelMap);

    EXPECT_CALL(*mockGridBuilder, setPrecursorBoundaryCondition).Times(6);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setPrecursorBoundaryCondition(sideType, fileCollection, timeStepsBetweenReads, velocityX, velocityY,
                                           velocityZ, fileLevelToGridLevelMap);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, notPeriodic)
{
    bool periodic_X = false;
    bool periodic_Y = false;
    bool periodic_Z = false;

    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(periodic_X, periodic_Y, periodic_Z)).Times(24);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setPeriodicBoundaryCondition(periodic_X, periodic_Y, periodic_Z);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, periodicAllDirectionsMultiGPU_allProcesses)
{
    bool periodic_X = true;
    bool periodic_Y = true;
    bool periodic_Z = true;

    // no local periodicity, periodicity is realized by setting up inter-gpu communication
    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(periodic_X, periodic_Y, periodic_Z)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(false, false, false)).Times(24);
    for (int i = 0; i < 24; i++) {
        this->createNewSut();
        sut->createGrids(i);
        sut->setPeriodicBoundaryCondition(periodic_X, periodic_Y, periodic_Z);
    }
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, periodicAllDirectionsMultiGPU_process0)
{
    sut->createGrids(0);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MX, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MY, false));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MZ, false));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MX, 1));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MY, 4));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MZ, 18));

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, false)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PY, false)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PZ, false)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PZ, testing::_)).Times(0);

    sut->setPeriodicBoundaryCondition(true, true, true);
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, periodicAllDirectionsMultiGPU_process9)
{
    sut->createGrids(9);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, false));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, 8));

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PY, false)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PZ, false)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PZ, testing::_)).Times(0);

    sut->setPeriodicBoundaryCondition(true, true, true);
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, periodicAllDirectionsMultiGPU_process17)
{
    sut->createGrids(17);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, testing::_));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PY, testing::_));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, 16));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PY, 13));

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PZ, testing::_)).Times(0);

    sut->setPeriodicBoundaryCondition(true, true, true);
}

TEST_F(MultipleGridBuilderFacadeTest_24subdomains, periodicAllDirectionsMultiGPU_process23)
{
    sut->createGrids(23);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, testing::_));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PY, testing::_));
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PZ, testing::_));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, 22));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PY, 19));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PZ, 5));

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MZ, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MX, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MY, testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MZ, testing::_)).Times(0);

    sut->setPeriodicBoundaryCondition(true, true, true);
}

class MultipleGridBuilderFacadeTest_CreateMockAndSut : public testing::Test
{
protected:
    std::shared_ptr<MockMultipleGridBuilder> mockGridBuilder =
        std::make_shared<testing::NiceMock<MockMultipleGridBuilder>>();
    MultipleGridBuilderFacade sut =
        MultipleGridBuilderFacade(mockGridBuilder, std::make_unique<GridDimensions>(0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.1));
};

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, periodicAllDirectionsSingleGPU)
{
    sut.setOverlapOfSubdomains(0.1);
    sut.createGrids(0);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices).Times(0);
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess).Times(0);
    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(true, true, true));

    sut.setPeriodicBoundaryCondition(true, true, true);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, periodicXY2GPUs)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, testing::_));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, 1));
    sut.createGrids(0);
    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::MX, testing::_));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::MX, 1));
    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(false, true, false));
    sut.setPeriodicBoundaryCondition(true, true, false);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, periodicYZ2GPUs)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);

    EXPECT_CALL(*mockGridBuilder, findCommunicationIndices(CommunicationDirections::PX, testing::_));
    EXPECT_CALL(*mockGridBuilder, setCommunicationProcess(CommunicationDirections::PX, 1));
    sut.createGrids(0);
    EXPECT_CALL(*mockGridBuilder, setPeriodicBoundaryCondition(false, true, true));
    sut.setPeriodicBoundaryCondition(false, true, true);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, setNumberOfLayersForRefinement_callsSetNumberOfLayersInGridBuilder)
{
    EXPECT_CALL(*mockGridBuilder, setNumberOfLayers(10, 8));
    sut.setNumberOfLayersForRefinement(10, 8);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addFineGrid_createGrid_callsAddGridAfterAddCoarseGrid)
{
    std::shared_ptr<Object> fineGrid = std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
    sut.addFineGrid(fineGrid, 1);

    testing::Expectation addCoarseGrid = EXPECT_CALL(*mockGridBuilder, addCoarseGrid);
    EXPECT_CALL(*mockGridBuilder, addGrid(fineGrid, 1)).After(addCoarseGrid);

    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, noFineGrid_createGrid_addGridNotCalled)
{
    EXPECT_CALL(*mockGridBuilder, addCoarseGrid);
    EXPECT_CALL(*mockGridBuilder, addGrid(testing::_)).Times(0);
    EXPECT_CALL(*mockGridBuilder, addGrid(testing::_, testing::_)).Times(0);

    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addFineGrid_createGridCallsMultipleFineGridsInOriginalOrder)
{
    std::shared_ptr<Object> fineGrid = nullptr;
    sut.addFineGrid(fineGrid, 2);
    sut.addFineGrid(fineGrid, 1);

    testing::Expectation addCoarseGrid = EXPECT_CALL(*mockGridBuilder, addCoarseGrid);
    testing::Expectation addFine2 = EXPECT_CALL(*mockGridBuilder, addGrid(fineGrid, 2)).After(addCoarseGrid);
    EXPECT_CALL(*mockGridBuilder, addGrid(fineGrid, 1)).After(addFine2);

    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addFineGrid_createGridCallsFineGridsInProcess0)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);
    std::shared_ptr<Object> fineGrid = std::make_shared<Sphere>(0.0, 0.0, 0.0, 10.0);
    sut.addFineGrid(fineGrid, 1);

    testing::Expectation addCoarseGrid = EXPECT_CALL(*mockGridBuilder, addCoarseGrid);
    EXPECT_CALL(*mockGridBuilder, addGrid(fineGrid, 1)).After(addCoarseGrid);
    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addFineGrid_createGridCallsFineGridsInProcess1)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);
    std::shared_ptr<Object> fineGrid = std::make_shared<Sphere>(0.0, 0.0, 0.0, 10.0);
    sut.addFineGrid(fineGrid, 1);

    testing::Expectation addCoarseGrid = EXPECT_CALL(*mockGridBuilder, addCoarseGrid);
    EXPECT_CALL(*mockGridBuilder, addGrid(fineGrid, 1)).After(addCoarseGrid);
    sut.createGrids(1);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addGeometry_createGrid_callsAddGeometryFunctionOfGridBuilder)
{
    std::shared_ptr<Object> geometry = std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
    sut.addGeometry(geometry);

    EXPECT_CALL(*mockGridBuilder, addGeometry);

    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, noFineGrid_createGrid_doesNotCallAddGeometryFunctionOfGridBuilder)
{
    EXPECT_CALL(*mockGridBuilder, addGeometry).Times(0);
    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addGeometry_createGridCallsAddGeometryFunctionOfGridBuilderProcess0)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);
    std::shared_ptr<Object> geometry = std::make_shared<Sphere>(0.0, 0.0, 0.0, 10.0);
    sut.addGeometry(geometry);

    EXPECT_CALL(*mockGridBuilder, addGeometry).Times(1);
    sut.createGrids(0);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addGeometry_createGridCallsAddGeometryFunctionOfGridBuilderProcess1)
{
    sut.addDomainSplit(1.0, Axis::x);
    sut.setOverlapOfSubdomains(0.1);
    std::shared_ptr<Object> geometry = std::make_shared<Sphere>(0.0, 0.0, 0.0, 10.0);
    sut.addGeometry(geometry);

    EXPECT_CALL(*mockGridBuilder, addGeometry).Times(1);
    sut.createGrids(1);
}

TEST(MultipleGridBuilderFacadeTest, setBoundaryCondition_createGridsNotCalledBeforeBC_throws)
{
    auto sut = MultipleGridBuilderFacadeTest(nullptr, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.1);
    EXPECT_THROW(sut.setNoSlipBoundaryCondition(SideType::PX), std::runtime_error);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addFineGrid_calledAfterCreateGrids_throws)
{
    sut.createGrids(0);
    std::shared_ptr<Object> fineGrid = std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
    EXPECT_THROW(sut.addFineGrid(fineGrid, 1), std::runtime_error);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addGeometry_calledAfterCreateGrids_throws)
{
    sut.createGrids(0);
    std::shared_ptr<Object> geometry = std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
    EXPECT_THROW(sut.addGeometry(geometry), std::runtime_error);
}

TEST_F(MultipleGridBuilderFacadeTest_CreateMockAndSut, addSplit_calledAfterCreateGrids_throws)
{
    sut.createGrids(0);
    std::shared_ptr<Object> fineGrid = std::make_shared<Cuboid>(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
    EXPECT_THROW(sut.addDomainSplit(1.0, Axis::x), std::runtime_error);
}

TEST(MultipleGridBuilderFacadeTest, createGrids_createGridsMoreThanOnceForSamePart_throws)
{
    auto niceMockGridBuilder = std::make_shared<testing::NiceMock<MockMultipleGridBuilder>>();
    auto sut = MultipleGridBuilderFacadeTest(niceMockGridBuilder, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, 0.1);
    sut.createGrids(0);
    EXPECT_THROW(sut.createGrids(0), std::runtime_error);
}
//! \}
