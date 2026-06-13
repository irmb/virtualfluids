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
//! \author Henry Korb
//=======================================================================================
#include "GbSpatialData3D.h"
#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <queue>
#include <vector>

std::unique_ptr<Node> insertNodes(const std::vector<real3>& coordinates)
{
    using iter = std::vector<uint>::iterator;
    std::function<std::unique_ptr<Node>(iter, iter, const std::vector<real3>&, uint)> recurse =
        [&](iter startIndices, iter endIndices, const std::vector<real3>& coordinates, uint level) {
            if (startIndices == endIndices)
                return std::unique_ptr<Node>();

            auto midPointer = std::next(startIndices, std::distance(startIndices, endIndices) / 2);

            std::nth_element(startIndices, midPointer, endIndices, [&](const auto& a, const auto& b) {
                if (level % 3 == 0)
                    return coordinates[a].x < coordinates[b].x;
                if (level % 3 == 1)
                    return coordinates[a].y < coordinates[b].y;
                return coordinates[a].z < coordinates[b].z;
            });

            return std::make_unique<Node>(*midPointer, recurse(startIndices, midPointer, coordinates, level + 1),
                                          recurse(std::next(midPointer), endIndices, coordinates, level + 1));
        };
    std::vector<uint> indices(coordinates.size());
    std::iota(indices.begin(), indices.end(), 0);
    return recurse(indices.begin(), indices.end(), coordinates, 0);
}

std::vector<uint> findNearestPoints(const Node* root, real3 point, uint nPoints, const std::vector<real3>& coordinates)
{
    std::vector<uint> result;
    nPoints = std::min<uint>(nPoints, (uint)coordinates.size());
    using DistPoint = std::pair<real, uint>;
    struct Compare
    {
        bool operator()(const DistPoint& a, const DistPoint& b) const
        {
            return a.first < b.first; // max-heap by distance
        }
    };

    std::priority_queue<DistPoint, std::vector<DistPoint>, Compare> pq;

    std::function<void(const Node*, uint)> recurse = [&](const Node* node, uint depth) {
        if (node == nullptr)
            return;
        const real3 coordinate = coordinates[node->index];
        const real dist = std::sqrt(square(coordinate - point));
        if (pq.size() < nPoints) {
            pq.emplace(dist, node->index);
        } else if (dist < pq.top().first) {
            pq.pop();
            pq.emplace(dist, node->index);
        }

        const uint axis = depth % 3;

        const real diff = (axis == 0 ? point.x - coordinate.x : axis == 1 ? point.y - coordinate.y : point.z - coordinate.z);

        const Node* first = (diff <= 0) ? node->left.get() : node->right.get();
        const Node* second = (diff <= 0) ? node->right.get() : node->left.get();

        recurse(first, depth + 1);

        const real worstDist = pq.size() < nPoints ? std::numeric_limits<real>::infinity() : pq.top().first;
        if (std::abs(diff) < worstDist) {
            recurse(second, depth + 1);
        }
    };

    recurse(root, 0);

    while (!pq.empty()) {
        result.push_back(pq.top().second);
        pq.pop();
    }
    std::reverse(result.begin(), result.end());
    return result;
}

void applyToTree(Node* root, std::function<void(uint, uint)> func)
{
    std::function<void(Node*, uint)> recurse = [&](Node* node, uint depth) {
        if (node == nullptr)
            return;
        func(node->index, depth);
        recurse(node->left.get(), depth + 1);
        recurse(node->right.get(), depth + 1);
    };
    recurse(root, 0);
}

