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
#include "gmock/gmock.h"
#include "basics/tests/testUtilities.h"

#include "geometries/Triangle/Triangle.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Vertex/Vertex.h"


using namespace testing;

TEST(BoundingBoxExactTest, findMinMaxFromTriangle)
{
    BoundingBox box = BoundingBox::makeInvalidMinMaxBox();

    real minX = 1.0f;
    real minY = 23.0f;
    real minZ = 1.1222f;

    real maxX = 110.0f;
    real maxY = 50.0f;
    real maxZ = 12122.23f;
    Vertex v1 = Vertex(maxX, maxY - 10, minZ + 2);
    Vertex v2 = Vertex(minX, maxY, maxZ);
    Vertex v3 = Vertex(minX + 3, minY, minZ);
    Vertex normal = Vertex(0.0f, 0.0f, 0.0f);
    Triangle t = Triangle(v1, v2, v3, normal);

    box.setMinMax(t);

    EXPECT_THAT(box.minX, RealEq(minX));
    EXPECT_THAT(box.minY, RealEq(minY));
    EXPECT_THAT(box.minZ, RealEq(minZ));

    EXPECT_THAT(box.maxX, RealEq(maxX));
    EXPECT_THAT(box.maxY, RealEq(maxY));
    EXPECT_THAT(box.maxZ, RealEq(maxZ));
}

TEST(BoundingBoxTest, isInside_true)
{
    BoundingBox box = BoundingBox();

    box.minX = 0.0f;
    box.minY = 0.0f;
    box.minZ = 0.0f;

    box.maxX = 10.0f;
    box.maxY = 10.0f;
    box.maxZ = 10.0f;

    Vertex v1 = Vertex(1, 1, 1);
    Vertex v2 = Vertex(2, 2, 2);
    Vertex v3 = Vertex(3, 3, 3);
    Vertex normal = Vertex(0.0f, 0.0f, 0.0f);
    Triangle t = Triangle(v1, v2, v3, normal);

    EXPECT_TRUE(box.isInside(t));
}

TEST(BoundingBoxTest, isInside_false)
{
    BoundingBox box = BoundingBox();

    box.minX = 0.0f;
    box.minY = 0.0f;
    box.minZ = 0.0f;

    box.maxX = 10.0f;
    box.maxY = 10.0f;
    box.maxZ = 10.0f;

    Vertex v1 = Vertex(1, 1, 1);
    Vertex v2 = Vertex(2, 2, 2);
    Vertex v3 = Vertex(3, 3, 11);
    Vertex normal = Vertex(0.0f, 0.0f, 0.0f);
    Triangle t = Triangle(v1, v2, v3, normal);

    EXPECT_FALSE(box.isInside(t));
}



//! \}
