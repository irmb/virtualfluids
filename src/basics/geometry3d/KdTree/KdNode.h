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
#ifndef KDNODE_H
#define KDNODE_H

#include <basics/memory/MbSmartPtr.h>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbTuple.h>

#include <geometry3d/GbTriFaceMesh3D.h>
#include <geometry3d/KdTree/KdRay.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>
#include <geometry3d/KdTree/KdUtilities.h>
#include <geometry3d/KdTree/intersectionhandler/KdLineIntersectionHandler.h>
#include <geometry3d/KdTree/intersectionhandler/KdRayIntersectionHandler.h>
#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

#include <string>
#include <vector>

namespace Kd
{
template <typename T>
class Node
{
public:
    Node(const T &x1, const T &y1, const T &z1, const T &x2, const T &y2, const T &z2,
         const MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFaces,
         std::vector<GbTriFaceMesh3D::Vertex> *ptrNodes)
        : child1(NULL), child2(NULL), triFaces(triFaces), ptrNodes(ptrNodes)
    {
        if (x1 < x2) {
            this->x[0] = x1;
            this->x[1] = x2;
        } else {
            this->x[0] = x2;
            this->x[1] = x1;
        }

        if (y1 < y2) {
            this->y[0] = y1;
            this->y[1] = y2;
        } else {
            this->y[0] = y2;
            this->y[1] = y1;
        }

        if (z1 < z2) {
            this->z[0] = z1;
            this->z[1] = z2;
        } else {
            this->z[0] = z2;
            this->z[1] = z1;
        }
    }
    /* ======================================================================================= */
    ~Node()
    {
        if (child1) {
            delete child1;
            child1 = NULL;
        }
        if (child2) {
            delete child2;
            child2 = NULL;
        }
    }
    /* ======================================================================================= */
    bool isLeaf() { return child1 == NULL && child2 == NULL; }
    /* ======================================================================================= */
    void deleteTriFaces() { triFaces = MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>>(); }
    /* ======================================================================================= */
    const MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> &getTriFaces() { return triFaces; }
    /* ======================================================================================= */
    std::vector<GbTriFaceMesh3D::Vertex> &getNodes()
    {
        if (!ptrNodes)
            throw UbException(UB_EXARGS, "ups,no nodes");
        return *ptrNodes;
    }

    /* ======================================================================================= */
    void buildTree(const int &level, const int &maxLevel, const SplitAlgorithm<T> &splitAlg)
    {
        SplitCandidate<T> splitCandidate = splitAlg.findBestSplitCandidate(level, maxLevel, *this);

        if (splitCandidate.isValid) {

            MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFacesForChild1(
                new std::vector<GbTriFaceMesh3D::TriFace>);
            MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFacesForChild2(
                new std::vector<GbTriFaceMesh3D::TriFace>);

            splitAlg.distributeTriFaces(splitCandidate, *triFacesForChild1, *triFacesForChild2, *this);

            //////////////////////////////////////////////////////////////////////////
            // calculate center points and edges of new child nodes
            T x1_l = x[0], y1_l = y[0], z1_l = z[0];
            T x2_l = x[1], y2_l = y[1], z2_l = z[1];
            T x1_r = x[0], y1_r = y[0], z1_r = z[0];
            T x2_r = x[1], y2_r = y[1], z2_r = z[1];

            if (splitCandidate.axis == Axis::X) {
                x2_l = splitCandidate.position;
                x1_r = splitCandidate.position;
            } else if (splitCandidate.axis == Axis::Y) {
                y2_l = splitCandidate.position;
                y1_r = splitCandidate.position;
            } else {
                z2_l = splitCandidate.position;
                z1_r = splitCandidate.position;
            }
            // ----------------------------------------------------------------------
            // ----------------------------------------------------------------------

            if (triFacesForChild1->size() > 0) {
                if (this->child1)
                    delete this->child1;
                this->child1 = new Node(x1_l, y1_l, z1_l, x2_l, y2_l, z2_l, triFacesForChild1, ptrNodes);
                this->child1->buildTree(level + 1, maxLevel, splitAlg);
            }

            if (triFacesForChild2->size() > 0) {
                if (this->child2)
                    delete this->child2;
                this->child2 = new Node(x1_r, y1_r, z1_r, x2_r, y2_r, z2_r, triFacesForChild2, ptrNodes);
                this->child2->buildTree(level + 1, maxLevel, splitAlg);
            }
        }
    }
    /* ======================================================================================= */
    int intersectLineBoundingBox(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2)
    {
        const T &n1X = val<1>(n1);
        const T &n1Y = val<2>(n1);
        const T &n1Z = val<3>(n1);

        const T &n2X = val<1>(n2);
        const T &n2Y = val<2>(n2);
        const T &n2Z = val<3>(n2);

        if (UbMath::greater(UbMath::max(((n1X <= n2X ? x[0] : x[1]) - n1X) / (n2X - n1X),
                                        ((n1Y <= n2Y ? y[0] : y[1]) - n1Y) / (n2Y - n1Y),
                                        ((n1Z <= n2Z ? z[0] : z[1]) - n1Z) / (n2Z - n1Z)),
                            UbMath::min(((n1X > n2X ? x[0] : x[1]) - n1X) / (n2X - n1X),
                                        ((n1Y > n2Y ? y[0] : y[1]) - n1Y) / (n2Y - n1Y),
                                        ((n1Z > n2Z ? z[0] : z[1]) - n1Z) / (n2Z - n1Z)))) {
            return Intersection::NO_INTERSECTION;
        } else {
            return Intersection::INTERSECTION;
        }
    }
    /* ======================================================================================= */
    int intersectRayBoundingBox(const Ray<T> &ray)
    {
        T tmin = (x[ray.signX] - ray.originX) * ray.inv_directionX;
        T tmax = (x[1 - ray.signX] - ray.originX) * ray.inv_directionX;

        T tymin = (y[ray.signY] - ray.originY) * ray.inv_directionY;
        T tymax = (y[1 - ray.signY] - ray.originY) * ray.inv_directionY;

        if ((tmin > tymax) || (tymin > tmax)) {
            return false;
        }
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        T tzmin = (z[ray.signZ] - ray.originZ) * ray.inv_directionZ;
        T tzmax = (z[1 - ray.signZ] - ray.originZ) * ray.inv_directionZ;

        // if( (UbMath::greater( tmin, tzmax) ) || ( UbMath::greater( tzmin, tmax) ) )
        if ((tmin > tzmax) || (tzmin > tmax)) {
            return false;
        }
        // if(tzmin > tmin) tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        // return ( (tmin =< t1) && (tmax >= t0) );
        if (UbMath::greaterEqual(tmax, T(0.0))) {
            return Intersection::INTERSECTION;
        } else {
            return Intersection::NO_INTERSECTION;
        }
    }
    /* ======================================================================================= */
    bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2,
                       const LineIntersectionHandler<T> &iHandler)
    {
        return iHandler.intersectLine(n1, n2, *this, child1, child2);
    }
    /* ======================================================================================= */
    int intersectRay(const Ray<T> &ray, const RayIntersectionHandler<T> &iHandler, std::set<UbKeys::Key3<int>> &mailbox)
    {
        return iHandler.intersectRay(ray, *this, child1, child2, mailbox);
    }
    /* ======================================================================================= */
    int getNumOfTriFaces()
    {
        if (!child1 && !child2) {
            if (triFaces)
                return (int)triFaces->size();
            else
                return 0;
        } else {
            int sum = 0;

            if (child1)
                sum += child1->getNumOfTriFaces();
            if (child2)
                sum += child2->getNumOfTriFaces();

            return sum;
        }
    }
    /* ======================================================================================= */
    int getNumOfNodes()
    {
        if (!child1 && !child2) {
            return 1;
        } else {
            int sum = 0;
            if (child1)
                sum += child1->getNumOfNodes();
            if (child2)
                sum += child2->getNumOfNodes();

            return 1 + sum;
        }
    }
    /* ======================================================================================= */
    std::string toString()
    {
        return ""; //"[" + x1 + "," + y1 + "," + z1 + "]  -" + "  [" + x2 + "," + y2 + "," + z2 + "]";
    }
    /* ======================================================================================= */
    void addCubeInfo(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt8> &cells,
                     std::vector<std::string> &datanames, std::vector<std::vector<double>> &celldata)
    {
        nodes.push_back(makeUbTuple(float(x[0]), float(y[0]), float(z[0])));
        nodes.push_back(makeUbTuple(float(x[1]), float(y[0]), float(z[0])));
        nodes.push_back(makeUbTuple(float(x[1]), float(y[1]), float(z[0])));
        nodes.push_back(makeUbTuple(float(x[0]), float(y[1]), float(z[0])));

        nodes.push_back(makeUbTuple(float(x[0]), float(y[0]), float(z[1])));
        nodes.push_back(makeUbTuple(float(x[1]), float(y[0]), float(z[1])));
        nodes.push_back(makeUbTuple(float(x[1]), float(y[1]), float(z[1])));
        nodes.push_back(makeUbTuple(float(x[0]), float(y[1]), float(z[1])));

        cells.push_back(makeUbTuple(int(nodes.size() - 8), int(nodes.size() - 7), int(nodes.size() - 6),
                                    int(nodes.size() - 5), int(nodes.size() - 4), int(nodes.size() - 3),
                                    int(nodes.size() - 2), int(nodes.size() - 1)));
        datanames.resize(1);
        datanames[0] = "childs";
        celldata.resize(datanames.size());
        if (child1 && child2)
            celldata[0].push_back(2);
        else if (child1 || child2)
            celldata[0].push_back(1);
        else
            celldata[0].push_back(0);

        if (child1)
            child1->addCubeInfo(nodes, cells, datanames, celldata);
        if (child2)
            child2->addCubeInfo(nodes, cells, datanames, celldata);
    }

public:
    T x[2], y[2], z[2];

private:
    Node *child1;
    Node *child2;

    MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFaces;
    std::vector<GbTriFaceMesh3D::Vertex> *ptrNodes; // lediglich f�r Zugriff auf die Knoten!!!
};
} // namespace Kd
#endif // KDNODE_H

//! \}
