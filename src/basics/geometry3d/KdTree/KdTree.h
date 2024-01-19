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
#ifndef KDTREE_H
#define KDTREE_H

#include <basics/utilities/UbKeys.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/KdNode.h>
#include <geometry3d/KdTree/KdRay.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>
#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

#include <cmath>
#include <string>

namespace Kd
{
template <typename T>
class Tree
{
public:
    /* ======================================================================================= */
    Tree(GbTriFaceMesh3D &mesh, const SplitAlgorithm<T> &splitAlg) : rootNode(NULL) { this->buildTree(mesh, splitAlg); }
    /* ======================================================================================= */
    ~Tree()
    {
        if (rootNode) {
            delete rootNode;
            rootNode = NULL;
        }
    }
    /* ======================================================================================= */
    // the IntersectionHandler specifies how to handle the intersection
    bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2,
                       const LineIntersectionHandler<T> &iHandler)
    {
        return rootNode->intersectLine(n1, n2, iHandler);
    }
    /* ======================================================================================= */
    // the IntersectionHandler specifies how to handle the intersection
    int intersectRay(const Ray<T> &ray, const RayIntersectionHandler<T> &iHandler)
    {
        std::set<UbKeys::Key3<int>> mailbox;
        return rootNode->intersectRay(ray, iHandler, mailbox);
    }
    /* ======================================================================================= */
    int getNumOfNodes()
    {
        if (rootNode)
            return rootNode->getNumOfNodes();
        return 0;
    }
    /* ======================================================================================= */
    int getNumOfTriFaces()
    {
        if (rootNode)
            return rootNode->getNumOfTriFaces();
        return 0;
    }
    /* ======================================================================================= */
    std::string toString()
    {
        return ""; // Tree:: num of nodes: " + rootNode.getNumOfNodes() + ", primitives:" +
                   // rootNode.getNumOfPrimitives() + ", root_primitives:" + getNumOfPrimitives() + ", max_level:" +
                   // max_level;
    }
    /* ======================================================================================= */
    void buildTree(GbTriFaceMesh3D &mesh, const SplitAlgorithm<T> &splitAlg)
    {
        if (rootNode)
            delete rootNode;

        // create a copy of triangles
        MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFaces(
            new std::vector<GbTriFaceMesh3D::TriFace>(*mesh.getTriangles()));

        const int maxLevel =
            static_cast<int>(std::lround(8.0 + 1.3 * std::log((double)triFaces->size()))); // TODO: remove magic numbers

        rootNode =
            new Node<T>(T(mesh.getX1Minimum()), T(mesh.getX2Minimum()), T(mesh.getX3Minimum()), T(mesh.getX1Maximum()),
                        T(mesh.getX2Maximum()), T(mesh.getX3Maximum()), triFaces, mesh.getNodes());

        rootNode->buildTree(0, maxLevel, splitAlg);
    }
    void writeTree(const std::string &filename, WbWriter *writer = WbWriterVtkXmlBinary::getInstance())
    {
        if (rootNode) {
            std::vector<UbTupleFloat3> nodes;
            std::vector<UbTupleInt8> cubes;
            std::vector<std::string> datanames;
            std::vector<std::vector<double>> cubesdata;
            rootNode->addCubeInfo(nodes, cubes, datanames, cubesdata);
            writer->writeOctsWithCellData(filename, nodes, cubes, datanames, cubesdata);
        }
    }

private:
    Node<T> *rootNode;
};
} // namespace Kd

#endif // KDTREE_H

//! \}
