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
//! \addtogroup utilities_test utilities
//! \ingroup basics_test
//! \{
//! \author Hussein Alihussein
//=======================================================================================

#include <gmock/gmock.h>

#include <cmath>

#include <basics/geometry3d/GbTriFaceMesh3D.h>

TEST(GbTriFaceMesh3DDedupTest, RemovesDuplicateVerticesAndDegenerateTriangles)
{
    GbTriFaceMesh3D mesh;
    auto *nodes = mesh.getNodes();
    auto *tris  = mesh.getTriangles();
    nodes->clear();
    tris->clear();

    auto addVertex = [&](float x, float y, float z) {
        GbTriFaceMesh3D::Vertex v;
        v.x = x;
        v.y = y;
        v.z = z;
        nodes->push_back(v);
    };

    addVertex(0.0f, 0.0f, 0.0f); // 0
    addVertex(1.0f, 0.0f, 0.0f); // 1
    addVertex(0.0f, 1.0f, 0.0f); // 2
    addVertex(0.0f, 0.0f, 0.0f); // 3 (duplicate)
    addVertex(1.0f, 0.0f, 0.0f); // 4 (duplicate)
    addVertex(2.0f, 0.0f, 0.0f); // 5 (colinear)

    tris->emplace_back(0, 1, 2); // valid
    tris->emplace_back(3, 4, 2); // duplicate vertices -> should map to (0,1,2)
    tris->emplace_back(0, 0, 1); // duplicate index
    tris->emplace_back(0, 1, 4); // becomes duplicate index after dedup
    tris->emplace_back(0, 1, 5); // colinear

    mesh.deleteRedundantNodes();

    ASSERT_EQ(nodes->size(), 4u);
    ASSERT_EQ(tris->size(), 2u);

    for (auto &tri : *tris) {
        EXPECT_NE(tri.getIndexVertex1(), tri.getIndexVertex2());
        EXPECT_NE(tri.getIndexVertex1(), tri.getIndexVertex3());
        EXPECT_NE(tri.getIndexVertex2(), tri.getIndexVertex3());
        EXPECT_LT(tri.getIndexVertex1(), static_cast<int>(nodes->size()));
        EXPECT_LT(tri.getIndexVertex2(), static_cast<int>(nodes->size()));
        EXPECT_LT(tri.getIndexVertex3(), static_cast<int>(nodes->size()));

        const double area = tri.getArea(*nodes);
        EXPECT_GT(area, 1.0e-6);
    }

    {
        GbTriFaceMesh3D mesh2;
        auto *nodes2 = mesh2.getNodes();
        auto *tris2  = mesh2.getTriangles();
        nodes2->clear();
        tris2->clear();

        auto addVertex2 = [&](float x, float y, float z) {
            GbTriFaceMesh3D::Vertex v;
            v.x = x;
            v.y = y;
            v.z = z;
            nodes2->push_back(v);
        };

        const float eps = 1.0e-7f;
        addVertex2(0.0f, 0.0f, 0.0f);         // 0
        addVertex2(0.49f * eps, 0.0f, 0.0f);  // 1 (merge with 0)
        addVertex2(0.51f * eps, 0.0f, 0.0f);  // 2 (distinct)
        addVertex2(1.0f, 0.0f, 0.0f);         // 3
        addVertex2(0.0f, 1.0f, 0.0f);         // 4

        tris2->emplace_back(0, 3, 4);
        tris2->emplace_back(1, 3, 4);
        tris2->emplace_back(2, 3, 4);

        mesh2.deleteRedundantNodes();

        ASSERT_EQ(nodes2->size(), 4u);
        ASSERT_EQ(tris2->size(), 3u);

        const double expectedX = 0.51e-7;
        bool foundDistinct = false;
        for (const auto &node : *nodes2) {
            if (std::abs(static_cast<double>(node.x) - expectedX) < 1.0e-10 &&
                std::abs(static_cast<double>(node.y)) < 1.0e-12 &&
                std::abs(static_cast<double>(node.z)) < 1.0e-12) {
                foundDistinct = true;
                break;
            }
        }
        EXPECT_TRUE(foundDistinct);
    }

    {
        GbTriFaceMesh3D mesh3;
        auto *nodes3 = mesh3.getNodes();
        auto *tris3  = mesh3.getTriangles();
        nodes3->clear();
        tris3->clear();

        auto addVertex3 = [&](float x, float y, float z) {
            GbTriFaceMesh3D::Vertex v;
            v.x = x;
            v.y = y;
            v.z = z;
            nodes3->push_back(v);
        };

        const float height = 5.1e-8f;
        addVertex3(0.0f, 0.0f, 0.0f); // 0
        addVertex3(1.0e-4f, 0.0f, 0.0f); // 1 (below area threshold)
        addVertex3(3.0e-2f, 0.0f, 0.0f); // 2 (above area threshold)
        addVertex3(0.0f, height, 0.0f); // 3

        tris3->emplace_back(0, 1, 3);
        tris3->emplace_back(0, 2, 3);

        mesh3.deleteRedundantNodes();

        ASSERT_EQ(nodes3->size(), 4u);
        ASSERT_EQ(tris3->size(), 1u);
        const double area = (*tris3)[0].getArea(*nodes3);
        EXPECT_GT(area, 1.0e-11);
    }
}

//! \}
