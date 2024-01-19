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
#ifndef KDCOUNTLINEINTERSECTIONHANDLER_H
#define KDCOUNTLINEINTERSECTIONHANDLER_H

#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbTuple.h>

#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/KdNode.h>
#include <geometry3d/KdTree/KdUtilities.h>
#include <geometry3d/KdTree/intersectionhandler/KdLineIntersectionHandler.h>

#include <set>

namespace Kd
{
template <typename T>
class CountLineIntersectionHandler : public LineIntersectionHandler<T>
{
public:
    bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2, Node<T> &parent, Node<T> *&child1,
                       Node<T> *&child2) const override
    {
        if (parent.intersectLineBoundingBox(n1, n2) == Intersection::INTERSECTION) {
            if (parent.isLeaf()) {
                std::vector<GbTriFaceMesh3D::TriFace> &triFaces = *parent.getTriFaces();
                std::vector<GbTriFaceMesh3D::Vertex> &nodes     = parent.getNodes();

                for (std::size_t i = 0; i < triFaces.size(); i++) {
                    GbTriFaceMesh3D::TriFace &triFace = triFaces[i];

                    if (Kd::intersectLine(n1, n2, triFace, nodes))
                        return true;
                }
                return false;
            } else {
                if (child1) {
                    if (child1->intersectLine(n1, n2, *this))
                        return true;
                }
                if (child2) {
                    if (child2->intersectLine(n1, n2, *this))
                        return true;
                }
            }
        }
        return false;
    }
    /* ======================================================================================= */
};
} // namespace Kd

#endif // KDCOUNTLINEINTERSECTIONHANDLER_H

//! \}
