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
#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

#include <gpu/core/Output/QDebugVtkWriter.hpp>

MATCHER(DoubleNear5, "") {
    return abs(std::get<0>(arg) - std::get<1>(arg)) < 0.00001;
}

using namespace QDebugVtkWriter;

double calcVectorLength(const std::array<double, 3> coords, const std::array<double, 3> neighborCoords)
{
    return std::sqrt(std::pow((neighborCoords[0] - coords[0]), 2) + std::pow((neighborCoords[1] - coords[1]), 2) +
                     std::pow((neighborCoords[2] - coords[2]), 2));
}

TEST(QDebugVtkWriterTest, modifyLineLengthsForQsSameCoords3)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 1, 1 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.3, 0.3 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}

TEST(QDebugVtkWriterTest, modifyLineLengthDifferentCoords)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 2, 3 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.6, 0.9 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}

TEST(QDebugVtkWriterTest, modifyLineLengthNegativeCoord)
{
    const std::array<double, 3> coords = { 0, 0, 0 };
    std::array<double, 3> neighborCoords = { 1, 2, -3 };
    const real q = 0.3;
    const real initialLength = calcVectorLength(coords, neighborCoords);

    modifyLineLengthsForQs(coords, neighborCoords, q);

    std::array<double, 3> expectedNeighborCoords = { 0.3, 0.6, -0.9 };
    EXPECT_THAT(neighborCoords,testing::Pointwise(DoubleNear5(), expectedNeighborCoords));
    EXPECT_THAT(calcVectorLength(coords, neighborCoords), testing::DoubleNear(q*initialLength, 0.00001));
}

//! \}
