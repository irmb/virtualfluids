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
#ifndef KDSAHSPLIT_H
#define KDSAHSPLIT_H

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbMath.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/KdNode.h>
#include <geometry3d/KdTree/KdSplitCandidateManager.h>
#include <geometry3d/KdTree/KdUtilities.h>
#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

#include <cmath>
#include <vector>

namespace Kd
{
template <typename T>
class SAHSplit : public SplitAlgorithm<T>
{
public:
    /* ======================================================================================= */
    SplitCandidate<T> findBestSplitCandidate(const int &level, const int &maxLevel, Node<T> &node) const override
    {
        if (!node.getTriFaces())
            throw UbException(UB_EXARGS, "triFace NULL pointer");

        if (node.getTriFaces()->size() <= 1 // max triangles in node
            || level >= maxLevel) {
            return SplitCandidate<T>();
        }

        SplitCandidate<T> bestSplitCandidate;
        T minCN = Ub::inf;

        for (int splitAxis = 0; splitAxis < 3; splitAxis++) {
            SplitCandidateManager<T> sc;
            findPossibleSplitCandidates(splitAxis, node, sc);

            // incremental sweep to find best split position
            for (std::size_t i = 0; i < sc.size(); i++) {
                if (i == 0) {
                    sc[i].nl = sc.objects_starting_outside_left + sc.objects_fully_outside_node;
                    sc[i].nr = int(node.getTriFaces()->size()) - sc[0].np - sc[0].ending;
                } else {
                    sc[i].nl = sc[i - 1].nl + sc[i - 1].starting + sc[i - 1].np;
                    sc[i].nr = sc[i - 1].nr - sc[i].ending - sc[i].np;
                }

                this->calcSAH(sc[i], node);

                if (sc[i].Cn < minCN) {
                    minCN              = sc[i].Cn;
                    bestSplitCandidate = sc[i];
                }
            }
        }

        // automatic termination criterion (SAH)
        if (bestSplitCandidate.isValid && bestSplitCandidate.Cn >= node.getTriFaces()->size() * Ci) {
            return SplitCandidate<T>();
        }

        return bestSplitCandidate;
    }
    /* ======================================================================================= */
    void distributeTriFaces(const SplitCandidate<T> &candidate,
                            std::vector<GbTriFaceMesh3D::TriFace> &triFacesForChild1,
                            std::vector<GbTriFaceMesh3D::TriFace> &triFacesForChild2, Node<T> &node) const override
    {
        if (!node.getTriFaces())
            throw UbException(UB_EXARGS, "null pointer at triface list");

        std::vector<GbTriFaceMesh3D::TriFace> &srcTriFaces = *node.getTriFaces();
        std::vector<GbTriFaceMesh3D::Vertex> &srcNodes     = node.getNodes();
        std::vector<T> projection;

        for (std::size_t i = 0; i < srcTriFaces.size(); i++) {
            GbTriFaceMesh3D::TriFace &triFace = srcTriFaces[i];
            Kd::project2Axis(triFace, srcNodes, candidate.axis, projection);

            T &min = projection[0];
            T &max = projection[2];
            // --------------------------------------------------- //
            // case 1 : object inside plane
            if (UbMath::equal(min, max)) {
                if (UbMath::equal(min, candidate.position)) {
                    if (candidate.np_left) {
                        triFacesForChild1.push_back(triFace);
                    } else if (candidate.np_right) {
                        triFacesForChild2.push_back(triFace);
                    }
                } else if (UbMath::less(min, candidate.position)) {
                    triFacesForChild1.push_back(triFace);
                } else // if( UbMath::greater(min, candidate.position)
                {
                    triFacesForChild2.push_back(triFace);
                }
            } //
            // --------------------------------------------------- //
            // case 2 : object on left side of plane
            else if (UbMath::lessEqual(max, candidate.position)) {
                triFacesForChild1.push_back(triFace);
            } // --------------------------------------------------- //
            // case 3 : object on right side of plane
            else if (UbMath::greaterEqual(min, candidate.position)) {
                triFacesForChild2.push_back(triFace);
            } //
            // --------------------------------------------------- //
            // case 4 : object in both nodes
            else {
                triFacesForChild1.push_back(triFace);
                triFacesForChild2.push_back(triFace);
            } //
              // --------------------------------------------------- //
        }

        node.deleteTriFaces();
    }

private:
    /* ======================================================================================= */
    // cost function
    inline T calcCosts(const int &nl, const int &nr, const T &SA_VL, const T &SA_VR, const T &SA_V) const
    {
        return Ct + Ci * (nl * SA_VL / SA_V + nr * SA_VR / SA_V);
    }
    /* ======================================================================================= */
    void findPossibleSplitCandidates(const int &splitAxis, Node<T> &node,
                                     SplitCandidateManager<T> &splitCandidateManager) const
    {
        T p1_node = (splitAxis == Axis::X ? node.x[0] : splitAxis == Axis::Y ? node.y[0] : node.z[0]);
        T p2_node = (splitAxis == Axis::X ? node.x[1] : splitAxis == Axis::Y ? node.y[1] : node.z[1]);

        if (!node.getTriFaces())
            throw UbException(UB_EXARGS, "null pointer");

        std::vector<GbTriFaceMesh3D::TriFace> &srcTriFaces = *node.getTriFaces();
        std::vector<GbTriFaceMesh3D::Vertex> &srcNodes     = node.getNodes();
        std::vector<T> projection;

        for (std::size_t i = 0; i < srcTriFaces.size(); i++) {
            GbTriFaceMesh3D::TriFace &triFace = srcTriFaces[i];

            // project object to axis
            Kd::project2Axis(triFace, srcNodes, splitAxis, projection);
            // left point
            T &p1 = projection[0];
            // right point
            T &p2 = projection[2];

            // --------------------------------------------------- //
            // --------------------------------------------------- //
            // case 1 : object is fully inside the current node
            if (UbMath::greaterEqual(p1, p1_node) && UbMath::lessEqual(p2, p2_node)) {
                if (UbMath::equal(p1, p2)) {
                    // object is inside the plane
                    splitCandidateManager.add(p1, splitAxis, 0, 0, 1);
                } else {
                    splitCandidateManager.add(p1, splitAxis, 1, 0, 0);
                    splitCandidateManager.add(p2, splitAxis, 0, 1, 0);
                }
            } //
            // --------------------------------------------------- //
            // --------------------------------------------------- //
            // case 2 : just the right point (p2) is inside the current node
            else if (UbMath::less(p1, p1_node) && UbMath::lessEqual(p2, p2_node) && UbMath::greaterEqual(p2, p1_node)) {
                splitCandidateManager.add(p2, splitAxis, 0, 1, 0);
                splitCandidateManager.objects_starting_outside_left++;
            } //
            // --------------------------------------------------- //
            // --------------------------------------------------- //
            // case 3 : just the left point (p1) is inside the current node
            else if (UbMath::greaterEqual(p1, p1_node) && UbMath::greater(p2, p2_node) &&
                     UbMath::lessEqual(p1, p2_node)) {
                splitCandidateManager.add(p1, splitAxis, 1, 0, 0);
            } //
            // --------------------------------------------------- //
            // --------------------------------------------------- //
            // case 4 : left and right point are outside the current node
            else if (UbMath::less(p1, p1_node) && UbMath::greater(p2, p2_node)) {
                splitCandidateManager.objects_fully_outside_node++;
            } //
              // --------------------------------------------------- //
              // --------------------------------------------------- //
        }

        splitCandidateManager.createSortedArray();
    }

    /* ======================================================================================= */
    // calculates the costs for a given splitCandidate based on the Surface Area Heuristic (SAH)
    void calcSAH(SplitCandidate<T> &candidate, Node<T> &node) const
    {
        T p1_node = (candidate.axis == Axis::X ? node.x[0] : candidate.axis == Axis::Y ? node.y[0] : node.z[0]);

        // edges of (root) voxel
        T dx = std::fabs(node.x[1] - node.x[0]);
        T dy = std::fabs(node.y[1] - node.y[0]);
        T dz = std::fabs(node.z[1] - node.z[0]);

        // surface area (root) voxel
        T SA_V = T((2.0 * dx * dy) + (2.0 * dx * dz) + (2.0 * dy * dz));

        T delta  = (candidate.axis == Axis::X ? dx : candidate.axis == Axis::Y ? dy : dz);
        T deltaL = std::fabs(candidate.position - p1_node);
        T deltaR = std::fabs(delta - deltaL);

        // edges of sub voxel left
        T dx_l = (candidate.axis == Axis::X ? deltaL : dx), dy_l = (candidate.axis == Axis::Y ? deltaL : dy),
          dz_l = (candidate.axis == Axis::Z ? deltaL : dz);

        // surface area sub voxel left
        T SA_VL = T((2.0 * dx_l * dy_l) + (2.0 * dx_l * dz_l) + (2.0 * dy_l * dz_l));

        // edges of sub voxel right
        T dx_r = (candidate.axis == Axis::X ? deltaR : dx), dy_r = (candidate.axis == Axis::Y ? deltaR : dy),
          dz_r = (candidate.axis == Axis::Z ? deltaR : dz);

        // surface area sub voxel right
        T SA_VR = T((2.0 * dx_r * dy_r) + (2.0 * dx_r * dz_r) + (2.0 * dy_r * dz_r));

        if (candidate.np == 0) {
            candidate.Cn = calcCosts(candidate.nl, candidate.nr, SA_VL, SA_VR, SA_V);
            return;
        }

        // once putting np with nl, and once with nr - and select the one with lowest cost
        // see: Wald, Havran: "On building fast kd-Trees for Ray Tracing, and doing that in O(N log N)", 2006
        T CP_L = calcCosts(candidate.nl + candidate.np, candidate.nr, SA_VL, SA_VR, SA_V);
        T CP_R = calcCosts(candidate.nl, candidate.nr + candidate.np, SA_VL, SA_VR, SA_V);

        if (CP_L < CP_R) {
            candidate.Cn       = CP_L;
            candidate.np_right = true;
        } else {
            candidate.Cn      = CP_R;
            candidate.np_left = true;
        }
    }
    /* ======================================================================================= */

protected:
    static const T Ct; // = 3.0; traversal cost
    static const T Ci; // = 4.0; ray-patch-intersection-cost
};

template <typename T>
const T SAHSplit<T>::Ct = 3.0; // traversal cost
template <typename T>
const T SAHSplit<T>::Ci = 4.0; // ray-patch-intersection-cost
} // namespace Kd

#endif // KDSAHSPLIT_H

//! \}
