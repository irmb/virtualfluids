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
//! \addtogroup gpu_geometries_tests geometries
//! \ingroup gpu_GridGenerator_tests GridGenerator
//! \{
#include <gmock/gmock.h>
#include <array>

#include "VerticalCylinder.h"

class VerticalCylinderTest : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    VerticalCylinder cylinder = VerticalCylinder(center[0], center[1], center[2], radius, height);
};

TEST_F(VerticalCylinderTest, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(VerticalCylinderTest, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-0.5 * height + center[2]));
}

TEST_F(VerticalCylinderTest, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(0.5 * height + center[2]));
}

TEST_F(VerticalCylinderTest, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - radius + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + radius + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - radius + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + radius + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - 0.5 * height + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - 0.5 * height - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + 0.5 * height - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + 0.5 * height + epsilon, 0.0, 0.0));
}
//! \}
