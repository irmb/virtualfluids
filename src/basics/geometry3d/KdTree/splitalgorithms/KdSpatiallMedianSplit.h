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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef SPATIALLMEDIANSPLIT_H
#define SPATIALLMEDIANSPLIT_H

#include <basics/utilities/UbMath.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

namespace Kd
{
template <typename T>
class SpatialMedianSplit : public SplitAlgorithm<T>
{
    /* ======================================================================================= */
    SplitCandidate<T> findBestSplitCandidate(const int &level, const int &maxLevel, Node<T> &node) const override
    {
        if (node.getTriFaces()->size() <= 24 // max triangles in node
            || level >= maxLevel) {
            return SplitCandidate<T>();
        }

        T dx = std::fabs(node.x[1] - node.x[0]);
        T dy = std::fabs(node.y[1] - node.y[0]);
        T dz = std::fabs(node.z[1] - node.z[0]);

        if (UbMath::equal(dx, UbMath::max(dx, dy, dz)))
            return SplitCandidate<T>(Axis::X, node.x[0] + 0.5 * dx, 0, 0, 0);
        else if (UbMath::equal(dy, UbMath::max(dy, dz)))
            return SplitCandidate<T>(Axis::Y, node.y[0] + 0.5 * dy, 0, 0, 0);

        return SplitCandidate<T>(Axis::Z, node.z[0] + 0.5 * dz, 0, 0, 0);
    }
    /* ======================================================================================= */
    void distributeTriFaces(const SplitCandidate<T> &candidate,
                            std::vector<GbTriFaceMesh3D::TriFace> &primitives_child1,
                            std::vector<GbTriFaceMesh3D::TriFace> &primitives_child2, Node<T> &node) const override
    {
        if (!node.getTriFaces())
            throw UbException(UB_EXARGS, "null pointer");

        std::vector<GbTriFaceMesh3D::TriFace> &srcTriFaces = *node.getTriFaces();
        std::vector<GbTriFaceMesh3D::Vertex> &srcNodes     = node.getNodes();
        std::vector<T> projection;

        for (std::size_t i = 0; i < srcTriFaces.size(); i++) {
            GbTriFaceMesh3D::TriFace &triFace = srcTriFaces[i];
            Kd::project2Axis(triFace, srcNodes, candidate.axis, projection);

            T &min = projection[0];
            T &max = projection[2];

            // case 1 : object inside plane
            if (UbMath::equal(min, max)) {
                if (UbMath::equal(min, candidate.position)) {
                    primitives_child1.push_back(triFace);
                    primitives_child2.push_back(triFace);
                } else if (UbMath::less(min, candidate.position)) {
                    primitives_child1.push_back(triFace);
                } else if (UbMath::greater(min, candidate.position)) {
                    primitives_child2.push_back(triFace);
                }
            }
            // case 2 : object on left side of plane
            else if (UbMath::lessEqual(max, candidate.position)) {
                primitives_child1.push_back(triFace);
            }
            // case 3 : object on right side of plane
            else if (UbMath::greaterEqual(min, candidate.position)) {
                primitives_child2.push_back(triFace);
            }
            // case 4 : object in both nodes
            else {
                primitives_child1.push_back(triFace);
                primitives_child2.push_back(triFace);
            }
        }

        node.deleteTriFaces();
    }
};
} // namespace Kd

#endif // SPATIALLMEDIANSPLIT_H

//! \}
