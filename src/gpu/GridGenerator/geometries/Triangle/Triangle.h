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
//! \file Triangle.h
//! \ingroup geometries
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef Triangle_h
#define Triangle_h

#include <memory>

#include "global.h"
#include "geometries/Vertex/Vertex.h"

class TriangleMemento;

struct Triangle
{
    Vertex v1, v2, v3, normal;
    real alphaAngles[3];
    real layerThickness;
    
    uint patchIndex;

    Triangle(Vertex &v1, Vertex &v2, Vertex &v3, Vertex &normal);
    Triangle(Vertex &v1, Vertex &v2, Vertex &v3);
    Triangle();

    void set(const Vertex &v1, const Vertex &v2, const Vertex &v3);
    void set(int index, Vertex value);
    Vertex get(int index);
    void calcNormal();

    void initalLayerThickness(real delta);


    Vertex getCenterOfMass() const;
    real getHalfAngleBetweenToAdjacentTriangle(const Triangle &t2) const;
    int isEqual(const Triangle &t2) const;
    bool doesNormalsShowToEachOther(const  Triangle &t2) const;
    int getCommonEdge(const Triangle &t2) const;

    bool contains(const Vertex& v)const;
    int getNumberOfCommonEdge(const Triangle &t2) const;
    int getTriangleIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnTri, real &qVal) const;
    void print() const;

    char isUnderFace(const Vertex &point) const;

    bool isUnterExtendedFace(const Vertex & point, real &s) const;
    bool isNotNextToFace(const Vertex &point) const;
    bool isUnderAngleToNeighbors(const Vertex &point) const;
    void getClosestPointsOnEdges(Vertex arr[], const Vertex &P) const;
    real getPerpedicularDistanceFrom(const Vertex &P) const;
    Vertex getPerpedicularPointFrom(const Vertex &P) const;
    bool isQNode(const Vertex & point, const real &s) const;
    bool isNegativeDirectionBorder(const Vertex & point) const;

    bool operator==(const Triangle &t) const;

    TriangleMemento getState() const;
    void setState(const TriangleMemento &memento);


    void setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ) const;
};

#endif
