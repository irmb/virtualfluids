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
//! \addtogroup gpu_geometries geometries
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "TriangleRefinement.h"

#include <GridGenerator/geometries/Triangle/Triangle.h>

TriangleRefinement::TriangleRefinement(std::vector<Triangle> *triangles)
{
    this->triangles = triangles;
}

TriangleRefinement::~TriangleRefinement()
{

}

void TriangleRefinement::redoubleTriangles()
{
    int counter = 0;
    for (size_t i = 0; i < triangles->size(); i++)
    {
        refine((int)i);
        counter++;
        if (counter % 50 == 0)
            printf("triangle refine: %d\n", (counter / 2));
    }
}

void TriangleRefinement::refineUntilMinDistance(double d_min)
{
    double d = 10e9;
    int counter = 0;

    while (d > d_min) {
        int triangleToRefine = findIndexFromTriangleWithLongestEdge(&d);
        refine(triangleToRefine);

        counter++;
        if (counter % 50 == 0)
            printf("triangle refine: %d, actual dMAX = %2.6f, d_min = %2.6f\n", counter, d, d_min);

    }
}

void TriangleRefinement::refineUntilcountTriangle(int countTri)
{
    double d = 10e9;
    int counter = 0;

    while (counter < countTri) {
        int triangleToRefine = findIndexFromTriangleWithLongestEdge(&d);
        refine(triangleToRefine);

        counter++;
        if (counter % 100 == 0)
            printf("triangle refine: %d, countTri = %d\n", counter, countTri);
    }
}

int TriangleRefinement::findIndexFromTriangleWithLongestEdge(double *d)
{
    *d = 0;
    double dTemp;
    int triangleToRefine = 0;
    for (size_t i = 0; i < triangles->size(); i++)
    {
        dTemp = getLongestEdgeDistance((*triangles)[i]);
        if (*d < dTemp) {
            *d = dTemp;
            triangleToRefine = (int)i;
        }
    }
    return triangleToRefine;
}


void TriangleRefinement::refine(int iTriangle)
{
    sortNeighborIndices();

    Triangle t = (*triangles)[iTriangle];
    int edge = getEdgeWithLongestDistance(t);

    std::vector<Vertex> v = getVertexArray(iTriangle);
    Vertex newVertex = getNewhalfVertexFromTrianglesEdge(v, edge);

    createTwoTriangles(v, edge, newVertex);

    int indexNeighbor = this->indices.ptr[iTriangle * indices.DIM + edge];
    v = getVertexArray(indexNeighbor);
    int commonEdgeNeigbor = findCommonEdgeFromTriangles(indexNeighbor, iTriangle);
    createTwoTriangles(v, commonEdgeNeigbor, newVertex);

    eraseOldTrianglesFromVector(iTriangle, indexNeighbor);
}

void TriangleRefinement::sortNeighborIndices()
{
    indices.DIM = 3;
    indices.size = (int)triangles->size();
    indices.ptr = new int[indices.DIM * indices.size];
    TriangleNeighborFinder neighbors((*triangles).data(), indices.size);
    neighbors.fillWithNeighborIndices(&indices, (*triangles).data());
}

std::vector<Vertex> TriangleRefinement::getVertexArray(int iTriangle)
{
    Triangle t = (*triangles)[iTriangle];

    std::vector<Vertex> v;
    v.resize(4);
    v[0] = t.v1;
    v[1] = t.v2;
    v[2] = t.v3;
    v[3] = t.normal;
    return v;
}

Vertex TriangleRefinement::getNewhalfVertexFromTrianglesEdge(std::vector<Vertex>  v, int edge)
{
    return getHalfVertex(v[edge], v[edge == 2 ? 0 : edge + 1]);
}

void TriangleRefinement::createTwoTriangles(std::vector<Vertex> v, int edge, Vertex newEdge)
{
    Vertex againstNewEdge = v[edge - 1 < 0 ? 2 : edge - 1];
    Vertex firstOldVertex = v[edge];
    Vertex secondOldVertex = v[edge + 1 > 2 ? 0 : edge + 1];

    Triangle firstNewTriangle(newEdge, againstNewEdge, firstOldVertex, v[3]);
    Triangle secondNewTriangle(newEdge, secondOldVertex, againstNewEdge, v[3]);

    (*triangles).push_back(firstNewTriangle);
    (*triangles).push_back(secondNewTriangle);
}

int TriangleRefinement::findCommonEdgeFromTriangles(int indexNeighbor, int iTriangle)
{
    int commonEdgeNeigbor = -1;
    for (int i = 0; i < indices.DIM; i++) {
        if (indices.ptr[indexNeighbor * indices.DIM + i] == iTriangle)
            commonEdgeNeigbor = i;
    }
    return commonEdgeNeigbor;
}

void TriangleRefinement::eraseOldTrianglesFromVector(int iTriangle, int indexNeighbor)
{
    (*triangles).erase((*triangles).begin() + iTriangle);
    if (iTriangle < indexNeighbor)
        indexNeighbor--;
    (*triangles).erase((*triangles).begin() + indexNeighbor);
}

void TriangleRefinement::refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle) 
{
    int edge = getEdgeWithLongestDistance(t);

    std::vector<Vertex> v;
    v.resize(4);
    v[0] = t.v1;
    v[1] = t.v2;
    v[2] = t.v3;
    v[3] = t.normal;
    Vertex newEdge = getHalfVertex(v[edge], v[edge == 2 ? 0 : edge + 1]);

    Vertex againstNewEdge = v[edge - 1 < 0 ? 2 : edge - 1];
    Vertex firstOldVertex = v[edge];
    Vertex secondOldVertex = v[edge + 1 > 2 ? 0 : edge + 1];

    firstNewTriangle = Triangle(newEdge, againstNewEdge, firstOldVertex, v[3]);
    secondNewTriangle = Triangle(newEdge, secondOldVertex, againstNewEdge, v[3]);
}


int TriangleRefinement::getEdgeWithLongestDistance(Triangle &t)
{
    real d1 = t.v2.getEuclideanDistanceTo(t.v1);
    real d2 = t.v3.getEuclideanDistanceTo(t.v2);
    real d3 = t.v1.getEuclideanDistanceTo(t.v3);

    real max = d1;
    int edge = 0;

    if (d2 > d1) {
        edge = 1;
        max = d2;
    }

    if (d3 > max){
        edge = 2;
        max = d3;
    }

    return edge;
}

real TriangleRefinement::getLongestEdgeDistance(Triangle &t) {

    int edge = getEdgeWithLongestDistance(t);
    Vertex v[3];
    v[0] = t.v1;
    v[1] = t.v2;
    v[2] = t.v3;

    if (edge == 2)
        return v[0].getEuclideanDistanceTo(v[2]);

    return v[edge + 1].getEuclideanDistanceTo(v[edge]);
}

Vertex TriangleRefinement::getHalfVertex(const Vertex &v, const Vertex &w)
{
    Vertex r;
    r.x = (v.x + w.x) / 2.0f;
    r.y = (v.y + w.y) / 2.0f;
    r.z = (v.z + w.z) / 2.0f;
    return r;
}

//! \}
