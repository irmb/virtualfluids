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
//! \addtogroup gpu_Output_tests Output
//! \ingroup gpu_core_tests core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "WriterUtilities.h"
#include "gpu/core/Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>
#include <numeric>

class WriterUtilitiesPeriodicCellTest : public testing::Test
{
protected:
    LBMSimulationParameter parH = LBMSimulationParameter();
    const uint level = 0;
    const uint baseNodeIndex = 0;
    const uint otherNodeIndex = 1;
    std::array<real, 2> coordinates = {0.0, 1.0};

    void SetUp() override
    {
        // create a domain with only three layers of nodes
        // nodes are at the coordinates 0.0, 1.0 and 2.0

        parH.gridSpacing = 1.0;

        parH.coordinateX = new real[2];
        parH.coordinateY = new real[2];
        parH.coordinateZ = new real[2];

        parH.coordinateX[baseNodeIndex] = coordinates[baseNodeIndex];
        parH.coordinateY[baseNodeIndex] = coordinates[baseNodeIndex];
        parH.coordinateZ[baseNodeIndex] = coordinates[baseNodeIndex];
        parH.coordinateX[otherNodeIndex] = coordinates[otherNodeIndex];
        parH.coordinateY[otherNodeIndex] = coordinates[otherNodeIndex];
        parH.coordinateZ[otherNodeIndex] = coordinates[otherNodeIndex];
    }

    void TearDown() override
    {
        delete[] parH.coordinateX;
        delete[] parH.coordinateY;
        delete[] parH.coordinateZ;
    }
};

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsNotPeriodic)
{
    EXPECT_FALSE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_FALSE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInX)
{
    parH.coordinateX[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInY)
{
    parH.coordinateY[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInZ)
{
    parH.coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}
TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXY)
{
    parH.coordinateX[1] = 2.0;
    parH.coordinateY[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXZ)
{
    parH.coordinateX[1] = 2.0;
    parH.coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInYZ)
{
    parH.coordinateY[1] = 2.0;
    parH.coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

TEST_F(WriterUtilitiesPeriodicCellTest, cellIsPeriodicInXYZ)
{
    parH.coordinateX[1] = 2.0;
    parH.coordinateY[1] = 2.0;
    parH.coordinateZ[1] = 2.0;
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, baseNodeIndex, otherNodeIndex));
    EXPECT_TRUE(WriterUtilities::isPeriodicCell(parH, otherNodeIndex, baseNodeIndex));
}

class WriterUtilitiesNeighborOctTest : public testing::Test
{
    static void setUpNeighborsNeighborsForOct(LBMSimulationParameter& parH, const std::array<uint, 8>& nodeIndices)
    {
        // node indices: MMM, PMM, PPM, MPM,
        //               MMP, PMP, PPP, MPP

        for (uint i = 0; i < (uint)nodeIndices.size(); i++) {
            const uint currentNodeIndex = nodeIndices[i];
            if (i < 4)
                parH.neighborZ[currentNodeIndex] = nodeIndices[i + 4];
            else
                parH.neighborZ[currentNodeIndex] = 99;

            if (i == 0 || i == 4)
                parH.neighborY[currentNodeIndex] = nodeIndices[i + 3];
            else if (i == 1 || i == 5)
                parH.neighborY[currentNodeIndex] = nodeIndices[i + 1];
            else
                parH.neighborY[currentNodeIndex] = 999;

            if (i == 0 || i == 4)
                parH.neighborX[currentNodeIndex] = nodeIndices[i + 1];
            else if (i == 3 || i == 7)
                parH.neighborX[currentNodeIndex] = nodeIndices[i - 1];
            else
                parH.neighborX[currentNodeIndex] = 9999;
        }
    }

public:
    LBMSimulationParameter parH = LBMSimulationParameter();
    std::array<uint, 8> nodeIndices;

    void SetUp() override
    {
        // set up some node indices from 0 to 7
        std::iota(nodeIndices.begin(), nodeIndices.end(), 0);
        std::reverse(nodeIndices.begin(), nodeIndices.end());

        parH.neighborX = new uint[8];
        parH.neighborY = new uint[8];
        parH.neighborZ = new uint[8];
        setUpNeighborsNeighborsForOct(parH, nodeIndices);
    }

    void TearDown() override
    {
        delete[] parH.neighborX;
        delete[] parH.neighborY;
        delete[] parH.neighborZ;
    }
};

TEST_F(WriterUtilitiesNeighborOctTest, getIndicesOfAllNodesInOct)
{
    std::array<uint, 8> resultingNodeIndices;
    WriterUtilities::getIndicesOfAllNodesInOct(resultingNodeIndices, nodeIndices[0], parH);
    for (uint i = 0; i < 8; i++)
        EXPECT_THAT(resultingNodeIndices[i], testing::Eq(nodeIndices[i])) << "for index i = " << i << " in nodeIndices";
}

TEST(WriterUtilitiesTest, calculateRelativeNodeIndexInPart)
{
    std::array<uint, 8> indicesOfOct = { 10, 12, 14, 16, 20, 22, 24, 26 };
    uint startPositionOfPart = 10;
    std::array<uint, 8> expected = { 0, 2, 4, 6, 10, 12, 14, 16 };

    std::array<uint, 8> result;
    WriterUtilities::calculateRelativeNodeIndexInPart(result, indicesOfOct, startPositionOfPart);
    EXPECT_THAT(result, testing::Eq(expected));
}

class WriterUtilitiesTestNodeValidity : public testing::Test
{
protected:
    LBMSimulationParameter parH = LBMSimulationParameter();
    std::array<uint, 8> nodeIndices;
    std::array<uint, 8> typeOfGridNode;

    void SetUp() override
    {
        // set up node indices from 0 to 7
        std::iota(nodeIndices.begin(), nodeIndices.end(), 0);

        std::fill(typeOfGridNode.begin(), typeOfGridNode.end(), GEO_FLUID);
        parH.typeOfGridNode = typeOfGridNode.data();
    }
};

TEST_F(WriterUtilitiesTestNodeValidity, allNodesInOctValidForWriting)
{
    uint endPositionOfPart = 7;
    EXPECT_TRUE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH, endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NodeOutOfPart)
{
    uint endPositionOfPart = 6;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH, endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NonFluidNode)
{
    uint endPositionOfPart = 7;
    typeOfGridNode[0] = GEO_SOLID;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH, endPositionOfPart));
}

TEST_F(WriterUtilitiesTestNodeValidity, areAllNodesInOctValidForWriting_NonFluidNodeAtEnd)
{
    uint endPositionOfPart = 7;
    typeOfGridNode[7] = GEO_SOLID;
    EXPECT_FALSE(WriterUtilities::areAllNodesInOctValidForWriting(nodeIndices, parH, endPositionOfPart));
}

//! \}
