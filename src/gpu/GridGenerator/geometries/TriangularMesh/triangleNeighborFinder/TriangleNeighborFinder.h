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
//! \file TriangleNeighborFinder.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef TriangleNeighborFinder_h
#define TriangleNeighborFinder_h

#include "global.h"
#include <vector>

struct IDS {
    enum IDs { vertexID = 0, firstVertexID = 1, coordinateID = 2, uniqueCoordID = 3, x = 3, y = 4, z = 5 };
};

struct Triangle;
struct Vertex;
class TriangularMesh;

struct IntegerPtr2D {
    int *ptr;
    int size;
    int DIM;
};

class TriangleNeighborFinder
{
public:
    TriangleNeighborFinder(Triangle *triangles, int size);
    ~TriangleNeighborFinder();
    
    std::vector<int> getTriangleIDsWithCommonVertex(int vertexID) const;
    std::vector< std::vector<Triangle> > getTrianglesPerVertex() const;

    void fillWithNeighborIndices(IntegerPtr2D *indices, Triangle *triangles);
    void fillWithNeighborAngles(TriangularMesh *geom) const;

    void printSortedToTriangles() const;
    void printSortedInSpace() const;

private:
    void initalSortedInSpaceWithCoords(Triangle *triangles_ptr, int size);
    void fillSortedInSpaceWithFirstVertexAndCoordinateIDs(int numberOfRows);
    void fillVectorWithIndicesOfTriangleNeighbors();
    unsigned int findTriangleID(unsigned int uniqueCoordID);
    Vertex getCoordinatesIDfromTriangle(int triangleID);
    bool isNeighborNotAlreadyInside(unsigned int iTriangle, unsigned int jTriangle);

    bool isTriangleNeighborOfParentTriangle(Vertex, Vertex);

    int numberOfRows;
    Triangle *triangles;

public:
    std::vector<std::vector<uint> > indicesOfTriangleNeighbors;
    real **sortedToTriangles;
    real **sortedInSpace;
};

#endif
