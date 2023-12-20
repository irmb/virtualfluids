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
//=======================================================================================
#include <gmock/gmock.h>

#include "Cylinder.h"


void testChangeSizeByDelta(Cylinder &cylinder)
{
    const double delta = 0.2;
    const double tolerance = 0.0000000001;
    const double oldRadius = cylinder.getRadius();
    const double oldHeight = cylinder.getHeight();
    const double oldMinimumX1 = cylinder.getX1Minimum();
    const double oldMinimumX2 = cylinder.getX2Minimum();
    const double oldMinimumX3 = cylinder.getX3Minimum();
    const double oldMaximumX1 = cylinder.getX1Maximum();
    const double oldMaximumX2 = cylinder.getX2Maximum();
    const double oldMaximumX3 = cylinder.getX3Maximum();

    cylinder.changeSizeByDelta(delta);

    EXPECT_THAT(cylinder.getRadius(), testing::Eq(oldRadius + delta));
    EXPECT_THAT(cylinder.getHeight(), testing::Eq(oldHeight + 2 * delta));
    EXPECT_THAT(cylinder.getX1Minimum(), testing::DoubleNear(oldMinimumX1 - delta, tolerance));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::DoubleNear(oldMinimumX2 - delta, tolerance));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::DoubleNear(oldMinimumX3 - delta, tolerance));
    EXPECT_THAT(cylinder.getX1Maximum(), testing::DoubleNear(oldMaximumX1 + delta, tolerance));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::DoubleNear(oldMaximumX2 + delta, tolerance));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::DoubleNear(oldMaximumX3 + delta, tolerance));
}

// CylinderTestAxisNormalToX ////////////////////////////////////////////////

class CylinderTestAxisNormalToX : public testing::Test
{
public:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder(center[0], center[1], center[2], radius, height, Axis::x);
};

TEST_F(CylinderTestAxisNormalToX, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToX, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-0.5 * height + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToX, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(0.5 * height + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToX, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - 0.5 * height + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - 0.5 * height - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + 0.5 * height - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + 0.5 * height + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - radius + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + radius - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + radius + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - radius + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - radius - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + radius - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + radius + epsilon, 0.0, 0.0));
}

TEST_F(CylinderTestAxisNormalToX, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}

// CylinderTestAxisNormalToY ////////////////////////////////////////////////

class CylinderTestAxisNormalToY : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Axis::y);
};

TEST_F(CylinderTestAxisNormalToY, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToY, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-0.5 * height + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToY, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(0.5 * height + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(radius + center[2]));
}

TEST_F(CylinderTestAxisNormalToY, isPointInObject)
{
    double epsilon = 0.0001;
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2], 0.0, 0.0));

    // x
    EXPECT_TRUE(cylinder.isPointInObject(center[0] - radius + epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] - radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0] + radius - epsilon, center[1], center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0] + radius + epsilon, center[1], center[2], 0.0, 0.0));

    // y
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] - 0.5 * height + epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] - 0.5 * height - epsilon, center[2], 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1] + 0.5 * height - epsilon, center[2], 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1] + 0.5 * height + epsilon, center[2], 0.0, 0.0));

    // z
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] - radius + epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] - radius - epsilon, 0.0, 0.0));
    EXPECT_TRUE(cylinder.isPointInObject(center[0], center[1], center[2] + radius - epsilon, 0.0, 0.0));
    EXPECT_FALSE(cylinder.isPointInObject(center[0], center[1], center[2] + radius + epsilon, 0.0, 0.0));
}

TEST_F(CylinderTestAxisNormalToY, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}

// CylinderTestAxisNormalToZ ////////////////////////////////////////////////

class CylinderTestAxisNormalToZ : public testing::Test
{
protected:
    std::array<double, 3> center = { 0.1, 0.2, 0.3 };
    double radius = 2.0;
    double height = 8.0;
    Cylinder cylinder = Cylinder({ center[0], center[1], center[2] }, radius, height, Axis::z);
};

TEST_F(CylinderTestAxisNormalToZ, getCentroid)
{
    EXPECT_THAT(cylinder.getX1Centroid(), testing::Eq(center[0]));
    EXPECT_THAT(cylinder.getX2Centroid(), testing::Eq(center[1]));
    EXPECT_THAT(cylinder.getX3Centroid(), testing::Eq(center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, getMinimum)
{
    EXPECT_THAT(cylinder.getX1Minimum(), testing::Eq(-radius + center[0]));
    EXPECT_THAT(cylinder.getX2Minimum(), testing::Eq(-radius + center[1]));
    EXPECT_THAT(cylinder.getX3Minimum(), testing::Eq(-0.5 * height + center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, getMaximum)
{
    EXPECT_THAT(cylinder.getX1Maximum(), testing::Eq(radius + center[0]));
    EXPECT_THAT(cylinder.getX2Maximum(), testing::Eq(radius + center[1]));
    EXPECT_THAT(cylinder.getX3Maximum(), testing::Eq(0.5 * height + center[2]));
}

TEST_F(CylinderTestAxisNormalToZ, isPointInObject)
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

TEST_F(CylinderTestAxisNormalToZ, changeSizeByDelta)
{
    testChangeSizeByDelta(cylinder);
}

//! \}
