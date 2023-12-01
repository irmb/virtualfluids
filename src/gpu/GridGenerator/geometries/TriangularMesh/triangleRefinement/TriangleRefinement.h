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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file TriangleRefinement.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef TriangleRefinement_h
#define TriangleRefinement_h

#include <vector>

#include "global.h"

#include "geometries/TriangularMesh/triangleNeighborFinder/TriangleNeighborFinder.h"

struct Triangle;
struct Vertex;
struct IntegerPtr2D;

class  TriangleRefinement
{
public:
    TriangleRefinement(std::vector<Triangle> *triangles);
    ~TriangleRefinement();

    void refine(int iTriangle);
    static void refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle);
    
    void refineUntilMinDistance(double d_min);
    void refineUntilcountTriangle(int countTri);
    void redoubleTriangles();

    static Vertex getHalfVertex(const Vertex &v, const Vertex &w);

private:
    std::vector<Triangle> *triangles;
    IntegerPtr2D indices;

    int findIndexFromTriangleWithLongestEdge(double *d);
    std::vector<Vertex> getVertexArray(int iTriangle);
    Vertex getNewhalfVertexFromTrianglesEdge(std::vector<Vertex> v, int edge);
    void sortNeighborIndices();
    void createTwoTriangles(std::vector<Vertex> v, int edge, Vertex newEdge);
    void eraseOldTrianglesFromVector(int iTriangle, int indexNeighbor);
    int findCommonEdgeFromTriangles(int indexNeighbor, int iTriangle);

public:
    static int getEdgeWithLongestDistance(Triangle &t);
    static real getLongestEdgeDistance(Triangle &t);
};


#endif
