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
//! \author Anna Wellmann
//=======================================================================================
#include <gmock/gmock.h>

#include "NeighborDebugWriter.hpp"
#include "gpu/core/Utilities/testUtilitiesGPU.h"

class WbWriterVtkXmlBinarySpy : public WbWriter
{
public:
    std::string writeLines(const std::string & /*filename*/, std::vector<UbTupleFloat3> &nodes,
                           std::vector<UbTupleInt2> &lines) override
    {
        this->nodes = nodes;
        this->lines = lines;
        return "";
    }
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> lines;

    std::string getFileExtension() override { return ""; }
};

class NeighborDebugWriterTest : public testing::Test
{
protected:
    void SetUp() override
    {
        typeOfGridNode = std::vector<uint>(numberOfNodes, GEO_FLUID);
        neighbors = std::vector<uint>(numberOfNodes, 2);
        coordinates = std::vector<real>(numberOfNodes, 1.0);
        coordinates[2] = 3.0;

        parH->numberOfNodes = numberOfNodes;
        parH->coordinateX = coordinates.data();
        parH->coordinateY = coordinates.data();
        parH->coordinateZ = coordinates.data();
        parH->neighborX = neighbors.data();
        parH->typeOfGridNode = typeOfGridNode.data();
    }

    const int level = 0;
    const unsigned long long numberOfNodes = 3;
    const uint direction = vf::lbm::dir::dP00; // x
    std::unique_ptr<LBMSimulationParameter> parH = std::make_unique<LBMSimulationParameter>();
    WbWriterVtkXmlBinarySpy writerSpy;
    std::vector<uint> typeOfGridNode;
    std::vector<uint> neighbors;
    std::vector<real> coordinates;
};

TEST_F(NeighborDebugWriterTest, writeNeighborLinkLines_onlyFLuidNodes_writesAllNodes)
{
    UbTupleFloat3 oneCoord(1.0, 1.0, 1.0);
    UbTupleFloat3 threeCoord(3.0, 3.0, 3.0);
    std::vector<UbTupleFloat3> expectedNodes = { oneCoord, threeCoord, oneCoord, threeCoord, threeCoord, threeCoord };
    std::vector<UbTupleInt2> expectedLines = { UbTupleInt2(0, 1), UbTupleInt2(2, 3), UbTupleInt2(4, 5) };

    NeighborDebugWriter::writeNeighborLinkLines(parH.get(), direction, "name", &writerSpy);

    EXPECT_THAT(writerSpy.nodes.size(), testing::Eq(numberOfNodes * 2));
    EXPECT_THAT(writerSpy.lines.size(), testing::Eq(numberOfNodes));
    EXPECT_THAT(writerSpy.nodes, testing::Eq(expectedNodes));
    EXPECT_THAT(writerSpy.lines, testing::Eq(expectedLines));
}

TEST_F(NeighborDebugWriterTest, writeNeighborLinkLines_fluidAndSolidNodes_writesOnlyFluidNodes)
{
    typeOfGridNode[2] = GEO_SOLID;
    
    UbTupleFloat3 oneCoord(1.0, 1.0, 1.0);
    UbTupleFloat3 threeCoord(3.0, 3.0, 3.0);
    std::vector<UbTupleFloat3> expectedNodes = { oneCoord, threeCoord, oneCoord, threeCoord};
    std::vector<UbTupleInt2> expectedLines = { UbTupleInt2(0, 1), UbTupleInt2(2, 3)};

    NeighborDebugWriter::writeNeighborLinkLines(parH.get(), direction, "name", &writerSpy);

    EXPECT_THAT(writerSpy.nodes.size(), testing::Eq((numberOfNodes-1) * 2));
    EXPECT_THAT(writerSpy.lines.size(), testing::Eq(numberOfNodes-1));
    EXPECT_THAT(writerSpy.nodes, testing::Eq(expectedNodes));
    EXPECT_THAT(writerSpy.lines, testing::Eq(expectedLines));
}

//! \}
